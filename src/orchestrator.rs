use crate::bam::BamReader;
use crate::config::{Config, StartStage};
use crate::log::{self, Stage};
use crate::{count, dedup, filter, gtf};
use anyhow::{bail, Context, Result};
use noodles::sam::alignment::RecordBuf;
use std::collections::HashMap;
use std::fs::File;
use std::os::unix::process::CommandExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::{Duration, Instant};

// STAR instances = min(threads/32, mem/(index_GB+20)), at least 1; 32 threads each.
pub fn star_instances(cfg: &Config) -> Result<usize> {
    let idx_gb = dir_size_gb(&cfg.reference.star_index)?;
    let by_cpu = cfg.num_threads / 32;
    let by_mem = (cfg.mem_limit as f64 / (idx_gb + 20.0)).floor() as usize;
    Ok(by_cpu.min(by_mem).max(1))
}

fn dir_size_gb(path: &str) -> Result<f64> {
    let mut bytes = 0u64;
    for entry in std::fs::read_dir(path).with_context(|| format!("reading {path}"))? {
        let m = entry?.metadata()?;
        if m.is_file() {
            bytes += m.len();
        }
    }
    Ok(bytes as f64 / 1e9)
}

// Run a child process with stdout+stderr redirected to a file (NOT piped — a pipe can
// deadlock if a grandchild inherits the write-end and blocks). On non-zero exit, print
// ONE FAILURE line + the last 15 log lines and bail with a short message; the full child
// output stays in `log_path` for inspection.
fn checked(cmd: &mut Command, what: &str, log_path: &Path) -> Result<()> {
    let f = File::create(log_path).with_context(|| format!("creating {}", log_path.display()))?;
    let ferr = f.try_clone().with_context(|| format!("cloning {}", log_path.display()))?;
    cmd.stdin(Stdio::null()).stdout(Stdio::from(f)).stderr(Stdio::from(ferr));
    let status = cmd.status().with_context(|| format!("failed to launch {what}"))?;
    if !status.success() {
        let code = status.code().map(|c| c.to_string()).unwrap_or_else(|| "signal".into());
        let tail = tail_lines(log_path, 15);
        log::failure(what, &format!("exit {code}"), Some(log_path));
        eprintln!("--- last 15 lines of {} ---\n{tail}\n---", log_path.display());
        bail!("{what} failed");
    }
    Ok(())
}

fn tail_lines(p: &Path, n: usize) -> String {
    let s = std::fs::read_to_string(p).unwrap_or_default();
    let lines: Vec<&str> = s.lines().collect();
    lines[lines.len().saturating_sub(n)..].join("\n")
}

pub fn run(cfg: &Config, start: Instant) -> Result<()> {
    // Honor num_threads for every rayon par_iter in the pipeline (build_global may only be
    // called once per process; a second call is a harmless no-op for our single-run CLI).
    let _ = rayon::ThreadPoolBuilder::new().num_threads(cfg.num_threads).build_global();

    let out = Path::new(&cfg.out_dir);
    for sub in ["", "filtered", "star", "qc", "expression", "logs"] {
        std::fs::create_dir_all(out.join(sub)).with_context(|| format!("creating {}/{sub}", cfg.out_dir))?;
    }
    let logs = out.join("logs");

    let stage = cfg.start_stage;
    let final_annot = build_final_annot(cfg)?;
    let n = if stage <= StartStage::Filtering {
        star_instances(cfg)?
    } else {
        count_shards(&out.join("filtered"))?
    };

    if stage <= StartStage::Filtering {
        let st = Stage::begin("FILTERING", format!("{n} shard(s)"));
        let fs = filter::filter(cfg, n, &st)?;
        run_fastqc(cfg, &st)?;
        st.done(format!("{} reads, {} passed", fs.total, fs.passed));
    }

    if stage <= StartStage::Mapping {
        let st = Stage::begin("MAPPING", format!("{n} STAR instance(s)"));
        let cdna_len = detect_cdna_len(out)?;
        run_star(cfg, n, &final_annot, cdna_len, &st)?;
        merge_bams(cfg, n, &st)?;
        st.done(format!("{n} instance(s) mapped + merged"));
    }

    if stage <= StartStage::Counting {
        let st = Stage::begin("COUNTING", "annotate, count, dedup");
        let ann = gtf::build(final_annot.to_str().unwrap(), cfg.counting_opts.introns)?;
        st.step(format!("annotation loaded ({} genes)", ann.genes.len()));
        let rt = count::count(cfg, &ann, &st)?;
        dedup::collapse(&rt, cfg, &ann, &st)?;
        let final_bam = out.join(format!("{}.bam", cfg.project));
        st.step("indexing final BAM");
        checked(
            Command::new("samtools").args(["index", "-@", &cfg.num_threads.to_string()]).arg(&final_bam),
            "samtools index",
            &logs.join("samtools_index.log"),
        )?;
        st.done(format!("{} cells, {} count rows", rt.barcodes.len(), rt.rows.len()));
    }

    if stage <= StartStage::Summarising {
        let st = Stage::begin("SUMMARISING", "dgecounts + QC report");
        run_summarise(cfg, &st)?;
        st.done("reports written");
    }

    log::success(&cfg.out_dir, start.elapsed());
    Ok(())
}

fn count_shards(dir: &Path) -> Result<usize> {
    let n = (0..).take_while(|i| dir.join(format!("shard_{i}.bam")).exists()).count();
    if n == 0 {
        bail!("no filtered/shard_*.bam found; cannot resume past Filtering");
    }
    Ok(n)
}

// GTF plus a one-exon "User" gene per additional_files contig (for STAR --genomeFastaFiles spike-ins).
fn build_final_annot(cfg: &Config) -> Result<PathBuf> {
    let path = Path::new(&cfg.out_dir).join("star").join(format!("{}.final_annot.gtf", cfg.project));
    std::fs::copy(&cfg.reference.gtf_file, &path)
        .with_context(|| format!("copying GTF {} -> {}", cfg.reference.gtf_file, path.display()))?;
    if let Some(fa) = &cfg.reference.additional_files {
        let mut extra = String::new();
        for (name, len) in fasta_lengths(fa)? {
            extra.push_str(&format!(
                "{name}\tUser\texon\t1\t{len}\t.\t+\t.\tgene_id \"{name}\"; transcript_id \"{name}\"; exon_number \"1\"; gene_name \"{name}\"; gene_biotype \"User\"; transcript_name \"{name}\"; exon_id \"{name}\"\n"
            ));
        }
        use std::io::Write;
        let mut f = std::fs::OpenOptions::new().append(true).open(&path).context("appending additional annotation")?;
        f.write_all(extra.as_bytes())?;
    }
    Ok(path)
}

fn fasta_lengths(path: &str) -> Result<Vec<(String, usize)>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading fasta {path}"))?;
    let mut out = Vec::new();
    let mut name: Option<String> = None;
    let mut len = 0usize;
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                out.push((n, len));
            }
            name = Some(h.split_whitespace().next().unwrap_or("").to_string());
            len = 0;
        } else {
            len += line.trim().len();
        }
    }
    if let Some(n) = name.take() {
        out.push((n, len));
    }
    Ok(out)
}

// mode of the first ~1000 shard read1 lengths -> sjdbOverhang basis
fn detect_cdna_len(out: &Path) -> Result<usize> {
    let mut reader = BamReader::open(&out.join("filtered").join("shard_0.bam"))?;
    let mut rec = RecordBuf::default();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    let mut seen = 0;
    while seen < 1000 && reader.next(&mut rec)? {
        if crate::bam::flags(&rec) & 0x40 != 0 {
            *counts.entry(crate::bam::seq_len(&rec)).or_default() += 1;
            seen += 1;
        }
    }
    Ok(counts.into_iter().max_by_key(|&(_, c)| c).map(|(l, _)| l).unwrap_or(100))
}

// Send SIGTERM (then, after a grace period, SIGKILL) to every STAR instance's process GROUP so
// the inherited `samtools view` grandchild that holds the FIFO open is reaped too — no orphan can
// be left holding a FIFO and wedging the run. Finally remove the star_tmp FIFO dirs.
fn teardown(kids: &mut [(usize, Child, PathBuf, PathBuf)]) {
    for (_, child, _, _) in kids.iter_mut() {
        unsafe { libc::killpg(child.id() as i32, libc::SIGTERM) };
    }
    std::thread::sleep(Duration::from_millis(200));
    for (_, child, _, _) in kids.iter_mut() {
        unsafe { libc::killpg(child.id() as i32, libc::SIGKILL) };
        let _ = child.wait();
    }
    for (_, _, _, tmp) in kids.iter() {
        let _ = std::fs::remove_dir_all(tmp);
    }
}

fn run_star(cfg: &Config, n: usize, final_annot: &Path, cdna_len: usize, stage: &Stage) -> Result<()> {
    let out = Path::new(&cfg.out_dir);
    let star_dir = out.join("star");
    let logs = out.join("logs");
    let threads = (cfg.num_threads / n).max(1);
    let overhang = cdna_len.saturating_sub(1).max(1);

    // Spawn all n STAR instances up front (n is small — bounded by the cpu/mem budget — so
    // spawning them all at once IS the intended concurrency). Each runs in its OWN process group
    // (process_group(0)) with stdout/stderr redirected to a file (never a pipe → no deadlock).
    let mut kids: Vec<(usize, Child, PathBuf, PathBuf)> = Vec::with_capacity(n);
    for i in 0..n {
        let prefix = star_dir.join(format!("{}.{i}.", cfg.project));
        let tmp = Path::new(&cfg.star_tmp).join(format!("{}.star.{i}", cfg.project));
        let _ = std::fs::remove_dir_all(&tmp); // STAR requires outTmpDir to not pre-exist
        let shard = out.join("filtered").join(format!("shard_{i}.bam"));
        let logf = logs.join(format!("star.{i}.log"));
        // Any failure here (log open, fd clone, fork/exec) must NOT orphan the STAR instances
        // already spawned (0..i) — each holds a samtools-view FIFO grandchild that would wedge the
        // run. Do the fallible work in an inner closure and tear everything down on error.
        let spawn_one = || -> Result<Child> {
            let f = File::create(&logf).with_context(|| format!("creating {}", logf.display()))?;
            let ferr = f.try_clone().with_context(|| format!("cloning {}", logf.display()))?;

            let mut cmd = Command::new("STAR");
            cmd.args(["--genomeDir", &cfg.reference.star_index])
                .arg("--readFilesIn")
                .arg(&shard)
                .args(["--readFilesType", "SAM", "PE"])
                .args(["--readFilesCommand", "samtools", "view"])
                .args(["--outSAMtype", "BAM", "SortedByCoordinate"])
                .arg("--outTmpDir")
                .arg(&tmp)
                .arg("--outFileNamePrefix")
                .arg(&prefix)
                .arg("--sjdbGTFfile")
                .arg(final_annot)
                .args(["--sjdbOverhang", &overhang.to_string()])
                .args(["--genomeLoad", "NoSharedMemory"])
                .args(["--runThreadN", &threads.to_string()]);
            if let Some(fa) = &cfg.reference.additional_files {
                cmd.args(["--genomeFastaFiles", fa]);
            }
            for a in cfg.reference.additional_star_params.split_whitespace() {
                cmd.arg(a);
            }
            cmd.process_group(0).stdin(Stdio::null()).stdout(Stdio::from(f)).stderr(Stdio::from(ferr));
            cmd.spawn().with_context(|| format!("launching STAR instance {i}"))
        };
        match spawn_one() {
            Ok(child) => kids.push((i, child, prefix, tmp)),
            Err(e) => {
                teardown(&mut kids);
                return Err(e);
            }
        }
    }
    stage.step(format!("{n} STAR instance(s) launched ({threads} threads each)"));

    // Poll: on the FIRST non-zero exit, tear down every peer immediately (prompt failure, seconds —
    // no waiting for the slowest survivor) and bail. Otherwise heartbeat every 10s until all exit.
    loop {
        let mut running = 0usize;
        let mut failure: Option<(usize, Option<i32>, PathBuf)> = None;
        let mut wait_err: Option<(usize, anyhow::Error)> = None;
        for (i, child, prefix, _) in kids.iter_mut() {
            match child.try_wait() {
                Ok(Some(st)) if !st.success() => {
                    failure = Some((*i, st.code(), prefix.clone()));
                    break;
                }
                Ok(Some(_)) => {}     // succeeded (status cached; subsequent try_wait is a no-op)
                Ok(None) => running += 1, // still running
                Err(e) => {
                    // A waitpid error must not orphan the remaining running instances (each holds a
                    // samtools-view FIFO grandchild) — tear them all down before propagating.
                    wait_err = Some((*i, e.into()));
                    break;
                }
            }
        }
        if let Some((i, e)) = wait_err {
            teardown(&mut kids);
            return Err(e).with_context(|| format!("waiting on STAR instance {i}"));
        }
        if let Some((i, code, prefix)) = failure {
            teardown(&mut kids);
            let logf = logs.join(format!("star.{i}.log"));
            let code_s = code.map(|c| c.to_string()).unwrap_or_else(|| "signal".into());
            log::failure(
                &format!("STAR instance {i}"),
                &format!("exit {code_s}; see {}Log.final.out", prefix.display()),
                Some(&logf),
            );
            eprintln!("--- last 15 lines of {} ---\n{}\n---", logf.display(), tail_lines(&logf, 15));
            bail!("STAR instance {i} failed");
        }
        if running == 0 {
            break;
        }
        stage.beat(Duration::from_secs(10), || {
            format!("{running}/{n} STAR instance(s) still running (elapsed {})", log::fmt_dur(stage.elapsed()))
        });
        std::thread::sleep(Duration::from_millis(500));
    }
    Ok(())
}

fn merge_bams(cfg: &Config, n: usize, stage: &Stage) -> Result<()> {
    let star_dir = Path::new(&cfg.out_dir).join("star");
    let logs = Path::new(&cfg.out_dir).join("logs");
    let merged = star_dir.join(format!("{}.merged.bam", cfg.project));
    stage.step(format!("merging {n} BAM(s)"));
    let mut cmd = Command::new("samtools");
    cmd.args(["merge", "-f", "-@", &cfg.num_threads.to_string()]).arg(&merged);
    for i in 0..n {
        cmd.arg(star_dir.join(format!("{}.{i}.Aligned.sortedByCoord.out.bam", cfg.project)));
    }
    checked(&mut cmd, "samtools merge", &logs.join("samtools_merge.log"))
}

// FastQC on the sampled filtered reads -> flat qc/.
fn run_fastqc(cfg: &Config, stage: &Stage) -> Result<()> {
    let qc = Path::new(&cfg.out_dir).join("qc");
    let logs = Path::new(&cfg.out_dir).join("logs");
    for label in ["R1", "R2"] {
        let fq = Path::new(&cfg.star_tmp).join(format!("qc_{label}.fastq.gz"));
        if !fq.exists() {
            continue;
        }
        stage.step(format!("FastQC {label}"));
        // -D/-R/-S take paths relative to CWD (the -o flag is ignored), so pass absolute qc/ paths
        let named = |suffix: &str| qc.join(format!("{}_{label}_{suffix}", cfg.project)).to_string_lossy().into_owned();
        checked(
            Command::new("falco")
                .arg("-D")
                .arg(named("fastqc_data.txt"))
                .arg("-R")
                .arg(named("fastqc.html"))
                .arg("-S")
                .arg(named("summary.txt"))
                .arg(&fq),
            &format!("FastQC {label}"),
            &logs.join(format!("falco.{label}.log")),
        )?;
        let _ = std::fs::remove_file(&fq);
    }
    Ok(())
}

fn run_summarise(cfg: &Config, stage: &Stage) -> Result<()> {
    let repo = repo_root()?;
    let logs = Path::new(&cfg.out_dir).join("logs");
    stage.step("finalize.R (dgecounts)");
    checked(
        Command::new("Rscript").arg(repo.join("R").join("finalize.R")).arg(&cfg.out_dir).arg(&cfg.project),
        "finalize.R",
        &logs.join("finalize.R.log"),
    )?;
    // qc_report.R also gets the config path so the report can show the run setup
    stage.step("qc_report.R");
    checked(
        Command::new("Rscript")
            .arg(repo.join("R").join("qc_report.R"))
            .arg(&cfg.out_dir)
            .arg(&cfg.project)
            .arg(&cfg.config_path),
        "qc_report.R",
        &logs.join("qc_report.R.log"),
    )?;
    Ok(())
}

// repo root = two levels above target/release/quartsx
fn repo_root() -> Result<PathBuf> {
    let exe = std::env::current_exe().context("locating quartsx executable")?;
    let repo = exe
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .context("resolving repo root from executable path")?;
    Ok(repo.to_path_buf())
}
