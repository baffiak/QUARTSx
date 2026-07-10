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

    // §6 lazy folders: create ONLY the out root + logs/ up front (logging is needed immediately).
    // filtered/ star/ qc/ expression/ are each created at their first write per stage (see
    // build_final_annot, the filtering block, run_fastqc, and count/dedup's own expression/ creation)
    // so a resume-from-a-later-stage run never litters unused empty folders.
    let out = Path::new(&cfg.out_dir);
    std::fs::create_dir_all(out).with_context(|| format!("creating {}", cfg.out_dir))?;
    let logs = out.join("logs");
    std::fs::create_dir_all(&logs).with_context(|| format!("creating {}/logs", cfg.out_dir))?;

    let stage = cfg.start_stage;

    // Preflight STAR BEFORE the (hours-long) filtering stage so a bad additional_STAR_params,
    // an un-clippable per-mate adapter, or a corrupt index fails in SECONDS, not after filtering.
    // Only when we will actually run STAR (stage <= Mapping). Returns the validated + per-mate-
    // expanded extra STAR argument tokens to hand to run_star.
    let star_params: Option<Vec<String>> =
        if stage <= StartStage::Mapping { Some(preflight_star(cfg)?) } else { None };

    // final_annot (GTF + spike-in exon lines) is only needed by MAPPING and COUNTING; skip it (and
    // the star/ folder + GTF copy it entails) on a Summarising-only resume (§6 lazy work).
    let final_annot: Option<PathBuf> =
        if stage <= StartStage::Counting { Some(build_final_annot(cfg)?) } else { None };

    let n = if stage <= StartStage::Filtering {
        star_instances(cfg)?
    } else {
        count_shards(&out.join("filtered"), &cfg.project)?
    };

    if stage <= StartStage::Filtering {
        // filtered/ is written here first (the checkpoint-1 shard(s)); create it lazily.
        std::fs::create_dir_all(out.join("filtered"))
            .with_context(|| format!("creating {}/filtered", cfg.out_dir))?;
        let st = Stage::begin("FILTERING", format!("{n} shard(s)"));
        let fs = filter::filter(cfg, n, &st)?;
        run_fastqc(cfg, &st)?;
        st.done(format!("{} reads, {} passed", fs.total, fs.passed));
    }

    if stage <= StartStage::Mapping {
        let st = Stage::begin("MAPPING", format!("{n} STAR instance(s)"));
        let annot = final_annot.as_ref().expect("final_annot built for stage <= Mapping");
        let cdna_len = detect_cdna_len(&out.join("filtered"), &cfg.project, n)?;
        run_star(cfg, n, annot, cdna_len, star_params.as_deref().unwrap_or(&[]), &st)?;
        // §2: single instance writes the coordinate-sorted checkpoint DIRECTLY (no merge/rename).
        // n>1: samtools-merge the per-instance sorted BAMs into the one no-number checkpoint.
        if n > 1 {
            merge_bams(cfg, n, &st)?;
            st.done(format!("{n} instance(s) mapped + merged"));
        } else {
            st.done("1 instance mapped (checkpoint written directly)");
        }
    }

    if stage <= StartStage::Counting {
        let st = Stage::begin("COUNTING", "annotate, count, dedup");
        let annot = final_annot.as_ref().expect("final_annot built for stage <= Counting");
        let ann = gtf::build(annot.to_str().unwrap(), cfg.counting_opts.introns)?;
        st.step(format!("annotation loaded ({} genes)", ann.genes.len()));
        let rt = count::count(cfg, &ann, &st)?;
        dedup::collapse(&rt, cfg, &ann, &st)?;
        // Final indexable checkpoint-3 target is star/{project}.tagged.bam (was {project}.bam).
        let final_bam = out.join("star").join(format!("{}.tagged.bam", cfg.project));
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

    // Remove the now-empty star_tmp parent left behind after each stage cleaned
    // its own contents. NON-RECURSIVE on purpose: if star_tmp is a shared/populated
    // location (e.g. /tmp), remove_dir fails harmlessly instead of deleting data.
    let _ = std::fs::remove_dir(&cfg.star_tmp);

    log::success(&cfg.out_dir, start.elapsed());
    Ok(())
}

// Count the filtered shard(s) present on disk when resuming past FILTERING, via the single naming
// authority bam::filtered_shard_path (no shard_{i}.bam literals). The multi-instance scheme is
// {project}.filtered.{i}.bam (probe with n=2 to select the numbered branch); the common single
// path is the no-suffix {project}.filtered.bam.
fn count_shards(dir: &Path, project: &str) -> Result<usize> {
    let multi = (0..)
        .take_while(|&i| crate::bam::filtered_shard_path(dir, project, i, 2).exists())
        .count();
    if multi > 0 {
        return Ok(multi);
    }
    if crate::bam::filtered_shard_path(dir, project, 0, 1).exists() {
        return Ok(1);
    }
    bail!("no {}/{project}.filtered*.bam found; cannot resume past Filtering", dir.display());
}

// GTF plus a one-exon "User" gene per additional_files contig (for STAR --genomeFastaFiles spike-ins).
fn build_final_annot(cfg: &Config) -> Result<PathBuf> {
    // star/ is created lazily here — the first write into it (§6).
    let star_dir = Path::new(&cfg.out_dir).join("star");
    std::fs::create_dir_all(&star_dir).with_context(|| format!("creating {}", star_dir.display()))?;
    let path = star_dir.join(format!("{}.final_annot.gtf", cfg.project));
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

// Mode of the first ~1000 shard read1 lengths -> STAR --sjdbOverhang basis. STAR manual (§2.2.3)
// recommends sjdbOverhang = max(mate length) - 1 for splice-junction database construction; here we
// derive that max from the actual filtered read1 length distribution (SS3xpress R1 is fixed-length
// after tag/UMI stripping). The tie-break (below) makes the pick deterministic (§1).
fn detect_cdna_len(filtered_dir: &Path, project: &str, n: usize) -> Result<usize> {
    let shard0 = crate::bam::filtered_shard_path(filtered_dir, project, 0, n);
    let mut reader = BamReader::open(&shard0, 0)?; // single-threaded: reading only ~1000 records
    let mut rec = RecordBuf::default();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    let mut seen = 0;
    while seen < 1000 && reader.next(&mut rec)? {
        if crate::bam::flags(&rec) & 0x40 != 0 {
            *counts.entry(crate::bam::seq_len(&rec)).or_default() += 1;
            seen += 1;
        }
    }
    // §1 determinism fix: HashMap iteration order is nondeterministic, so a plain max_by_key on the
    // count alone would break ties arbitrarily run-to-run. The key length `l` is unique across map
    // entries, so ordering by (count, Reverse(length)) is a TOTAL order -> the pick is fully
    // deterministic regardless of iteration order (ties broken toward the shorter length).
    Ok(counts
        .into_iter()
        .max_by_key(|&(l, c)| (c, std::cmp::Reverse(l)))
        .map(|(l, _)| l)
        .unwrap_or(100))
}

// STAR flags QUARTSx sets itself; a user override for any of these in additional_STAR_params would
// duplicate the flag and trigger STAR's fatal parameter error. Reject them in preflight.
const RESERVED_STAR_FLAGS: &[&str] = &[
    "--genomeDir",
    "--readFilesIn",
    "--readFilesType",
    "--readFilesCommand",
    "--outSAMtype",
    "--outFileNamePrefix",
    "--outTmpDir",
    "--sjdbGTFfile",
    "--sjdbOverhang",
    "--genomeLoad",
    "--runThreadN",
    "--outFilterMultimapNmax",
    "--outSAMmultNmax",
    "--winAnchorMultimapNmax",
    "--genomeFastaFiles",
];

// §2 preflight: fail in seconds (before filtering) on a broken STAR setup. Validates the STAR index
// integrity, dry-checks additional_STAR_params, and auto-expands the per-mate clip parameters for our
// 2-mate PE shard. Returns the validated + expanded extra STAR argument tokens.
fn preflight_star(cfg: &Config) -> Result<Vec<String>> {
    // 1) STAR index integrity. A STAR genome index (STAR --runMode genomeGenerate) always writes these
    //    core files (read/mmapped at load). Missing any of them means a truncated/incompatible index
    //    -> catch it now, not after filtering.
    let idx = Path::new(&cfg.reference.star_index);
    for f in ["Genome", "SA", "SAindex", "genomeParameters.txt", "chrNameLength.txt"] {
        if !idx.join(f).exists() {
            bail!(
                "STAR index at {} is missing '{f}'; it is incomplete or was not built with \
                 STAR --runMode genomeGenerate — regenerate the index",
                cfg.reference.star_index
            );
        }
    }

    // 2) dry-check additional_STAR_params tokens.
    let mut toks: Vec<String> =
        cfg.reference.additional_star_params.split_whitespace().map(str::to_string).collect();
    if let Some(first) = toks.first() {
        if !first.starts_with("--") {
            bail!("additional_STAR_params must begin with a --flag, found '{first}'");
        }
    }
    for t in &toks {
        if t.starts_with("--") && RESERVED_STAR_FLAGS.contains(&t.as_str()) {
            bail!(
                "additional_STAR_params sets '{t}', which QUARTSx controls internally; remove it \
                 (multimapper caps go under counting_opts, spike-ins under reference.additional_files)"
            );
        }
    }

    // 3) auto-expand per-mate clip params. Our shard is a 2-mate PE SAM (R1 tagged + R2 internal), so
    //    STAR 2.7.10b requires ONE value per mate for --clip3pAdapterSeq / --clip3pAdapterMMp; a lone
    //    value aborts with "--clip3pAdapterSeq has to contain 2 values to match the number of mates".
    //    Duplicate a single adapter sequence to both mates and pair it with the STAR-default mismatch
    //    prop 0.1/0.1.
    expand_clip_params(&mut toks);
    Ok(toks)
}

// Duplicate a lone --clip3pAdapterSeq value across both PE mates and ensure a matching 2-value
// --clip3pAdapterMMp (STAR default 0.1 per mate). No-op when the user gave no adapter or already
// supplied per-mate values.
fn expand_clip_params(toks: &mut Vec<String>) {
    // (index of flag, number of value tokens that follow before the next --flag), if present.
    fn locate(toks: &[String], flag: &str) -> Option<(usize, usize)> {
        let i = toks.iter().position(|t| t == flag)?;
        let vals = toks[i + 1..].iter().take_while(|t| !t.starts_with("--")).count();
        Some((i, vals))
    }
    if let Some((i, vals)) = locate(toks, "--clip3pAdapterSeq") {
        if vals == 1 {
            let v = toks[i + 1].clone();
            toks.insert(i + 2, v); // now SEQ SEQ
        }
        match locate(toks, "--clip3pAdapterMMp") {
            None => {
                toks.push("--clip3pAdapterMMp".into());
                toks.push("0.1".into());
                toks.push("0.1".into());
            }
            Some((j, 1)) => {
                let v = toks[j + 1].clone();
                toks.insert(j + 2, v); // 0.1 -> 0.1 0.1
            }
            _ => {}
        }
    }
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

// Surface STAR's OWN standard progress lines (the "..... started mapping" markers STAR prints to
// stderr, which we redirect to `logf`) as they are appended — WITHOUT re-reading already-shown bytes.
// The verbose per-instance STAR Log.out stays file-only; the noisy Log.progress.out is never shown
// (§2 CLOSED). `offset` advances only over complete (newline-terminated) lines so a line mid-write is
// left for the next poll. Non-UTF8 bytes are lossily decoded (STAR progress is ASCII).
fn surface_star_progress(logf: &Path, offset: &mut u64, i: usize, n: usize, stage: &Stage) {
    use std::io::{Read, Seek, SeekFrom};
    let mut f = match File::open(logf) {
        Ok(f) => f,
        Err(_) => return, // log not created yet
    };
    if f.seek(SeekFrom::Start(*offset)).is_err() {
        return;
    }
    let mut bytes = Vec::new();
    if f.read_to_end(&mut bytes).is_err() {
        return;
    }
    let end = match bytes.iter().rposition(|&b| b == b'\n') {
        Some(p) => p + 1, // consume up to and including the last newline
        None => return,   // no complete line yet
    };
    let text = String::from_utf8_lossy(&bytes[..end]);
    for line in text.lines() {
        // STAR standard progress lines look like "Jul 08 12:34:56 ..... started mapping". Show only
        // those; drop STAR's own leading timestamp (Stage::step re-stamps) by taking the text after
        // the "..... " marker.
        if let Some(idx) = line.find("..... ") {
            let msg = line[idx + 6..].trim_end();
            if n > 1 {
                stage.step(format!("[STAR {i}] {msg}"));
            } else {
                stage.step(msg);
            }
        }
    }
    *offset += end as u64;
}

fn run_star(
    cfg: &Config,
    n: usize,
    final_annot: &Path,
    cdna_len: usize,
    star_params: &[String],
    stage: &Stage,
) -> Result<()> {
    let out = Path::new(&cfg.out_dir);
    let star_dir = out.join("star");
    let logs = out.join("logs");
    let threads = (cfg.num_threads / n).max(1);
    let overhang = cdna_len.saturating_sub(1).max(1);

    // §2/§4: Unique mode caps multimappers to 1 locus (1/1) so STAR marks them UNMAPPED at MAPPING;
    // the resolver modes inject no cap flags and use STAR's own defaults.
    let unique_mode = cfg.counting_opts.multi_mappers == crate::config::MultiMapperMode::Unique;

    // Spawn all n STAR instances up front (n is small — bounded by the cpu/mem budget — so
    // spawning them all at once IS the intended concurrency). Each runs in its OWN process group
    // (process_group(0)) with stdout/stderr redirected to a file (never a pipe → no deadlock).
    let mut kids: Vec<(usize, Child, PathBuf, PathBuf)> = Vec::with_capacity(n);
    for i in 0..n {
        // §2 naming: single instance writes STAR's own checkpoint name directly (prefix "{project}.",
        // no number, no merge). Multiple instances use "{project}_{i}." and are merged afterwards.
        let prefix = if n == 1 {
            star_dir.join(format!("{}.", cfg.project))
        } else {
            star_dir.join(format!("{}_{i}.", cfg.project))
        };
        let tmp = Path::new(&cfg.star_tmp).join(format!("{}.star.{i}", cfg.project));
        let _ = std::fs::remove_dir_all(&tmp); // STAR requires outTmpDir to not pre-exist
        let shard = crate::bam::filtered_shard_path(&out.join("filtered"), &cfg.project, i, n);
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
            if unique_mode {
                cmd.args(["--outFilterMultimapNmax", "1"]).args(["--outSAMmultNmax", "1"]);
            } else {
                let cap = cfg.counting_opts.multimapper_cap.to_string();
                cmd.args(["--outFilterMultimapNmax", &cap]).args(["--outSAMmultNmax", &cap]);
            }
            if let Some(fa) = &cfg.reference.additional_files {
                cmd.args(["--genomeFastaFiles", fa]);
            }
            for a in star_params {
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

    // Poll: surface each instance's own STAR progress lines as they appear (no generic
    // "N/N still running" heartbeat — §2 CLOSED). On the FIRST non-zero exit, tear down every peer
    // immediately (prompt failure, seconds — no waiting for the slowest survivor) and bail.
    let mut offsets = vec![0u64; n];
    loop {
        for i in 0..n {
            surface_star_progress(&logs.join(format!("star.{i}.log")), &mut offsets[i], i, n, stage);
        }
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
        std::thread::sleep(Duration::from_millis(500));
    }
    // Final flush so the closing "..... finished successfully" line is surfaced.
    for i in 0..n {
        surface_star_progress(&logs.join(format!("star.{i}.log")), &mut offsets[i], i, n, stage);
    }
    Ok(())
}

// §2 (n>1 only): samtools-merge the per-instance coordinate-sorted BAMs into the ONE no-number
// checkpoint, then remove the per-instance files (keep exactly one ~90 GB checkpoint).
fn merge_bams(cfg: &Config, n: usize, stage: &Stage) -> Result<()> {
    let star_dir = Path::new(&cfg.out_dir).join("star");
    let logs = Path::new(&cfg.out_dir).join("logs");
    // STAR's ACTUAL emitted sorted-BAM name is "<prefix>Aligned.sortedByCoord.out.bam" (lowercase
    // 'sorted', capital 'Coord'); the merged checkpoint uses the no-number "{project}." prefix.
    let checkpoint = star_dir.join(format!("{}.Aligned.sortedByCoord.out.bam", cfg.project));
    stage.step(format!("merging {n} coordinate-sorted BAM(s) -> checkpoint"));
    let mut cmd = Command::new("samtools");
    cmd.args(["merge", "-f", "-@", &cfg.num_threads.to_string()]).arg(&checkpoint);
    let mut parts = Vec::with_capacity(n);
    for i in 0..n {
        let p = star_dir.join(format!("{}_{i}.Aligned.sortedByCoord.out.bam", cfg.project));
        cmd.arg(&p);
        parts.push(p);
    }
    checked(&mut cmd, "samtools merge", &logs.join("samtools_merge.log"))?;
    // IMPL NOTE (verified): `samtools merge` PRESERVES coordinate sort order from sorted inputs and
    // re-sort is needed ONLY if the inputs' @SQ header orders conflict (samtools-merge(1) manpage).
    // All instances map against the SAME STAR index -> identical @SQ order -> the merged output is
    // already coordinate-sorted; no separate `samtools sort` pass required.
    // Remove the per-instance BAMs now folded into the checkpoint (keep exactly ONE ~90 GB file).
    for p in parts {
        let _ = std::fs::remove_file(&p);
    }
    Ok(())
}

// FastQC (falco engine) on the sampled filtered reads -> flat qc/. QUARTSx's OWN outputs stay neutral
// (qc/, fastqc.*.log); only the invoked binary and the tool's in-report name are falco's (§6, §0).
fn run_fastqc(cfg: &Config, stage: &Stage) -> Result<()> {
    let qc = Path::new(&cfg.out_dir).join("qc");
    let logs = Path::new(&cfg.out_dir).join("logs");
    for label in ["R1", "R2"] {
        let fq = Path::new(&cfg.star_tmp).join(format!("qc_{label}.fastq.gz"));
        if !fq.exists() {
            continue;
        }
        // qc/ is created lazily at the first QC write (§6).
        std::fs::create_dir_all(&qc).with_context(|| format!("creating {}", qc.display()))?;
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
            &logs.join(format!("fastqc.{label}.log")),
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
