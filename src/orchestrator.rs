use crate::bam::BamReader;
use crate::config::{Config, StartStage};
use crate::{count, dedup, filter, gtf};
use anyhow::{bail, Context, Result};
use noodles::sam::alignment::RecordBuf;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

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

fn checked(cmd: &mut Command, what: &str) -> Result<()> {
    let out = cmd.output().with_context(|| format!("failed to launch {what}"))?;
    if !out.status.success() {
        bail!("{what} failed ({}):\n{}", out.status, String::from_utf8_lossy(&out.stderr));
    }
    Ok(())
}

pub fn run(cfg: &Config) -> Result<()> {
    let out = Path::new(&cfg.out_dir);
    for sub in ["", "filtered", "star", "qc", "expression"] {
        std::fs::create_dir_all(out.join(sub)).with_context(|| format!("creating {}/{sub}", cfg.out_dir))?;
    }

    let stage = cfg.start_stage;
    let final_annot = build_final_annot(cfg)?;
    let n = if stage <= StartStage::Filtering {
        star_instances(cfg)?
    } else {
        count_shards(&out.join("filtered"))?
    };

    if stage <= StartStage::Filtering {
        eprintln!("[quartsx] FILTERING: {n} shard(s)");
        filter::filter(cfg, n)?;
        run_fastqc(cfg)?;
    }

    if stage <= StartStage::Mapping {
        eprintln!("[quartsx] MAPPING: {n} STAR instance(s)");
        let cdna_len = detect_cdna_len(out)?;
        run_star(cfg, n, &final_annot, cdna_len)?;
        merge_bams(cfg, n)?;
    }

    if stage <= StartStage::Counting {
        eprintln!("[quartsx] COUNTING: annotate, count, dedup");
        let ann = gtf::build(final_annot.to_str().unwrap(), cfg.counting_opts.introns)?;
        let rt = count::count(cfg, &ann)?;
        dedup::collapse(&rt, cfg, &ann)?;
        let final_bam = out.join(format!("{}.bam", cfg.project));
        checked(
            Command::new("samtools").args(["index", "-@", &cfg.num_threads.to_string()]).arg(&final_bam),
            "samtools index",
        )?;
    }

    if stage <= StartStage::Summarising {
        eprintln!("[quartsx] SUMMARISING: dgecounts + QC report");
        run_summarise(cfg)?;
    }

    eprintln!("[quartsx] done: {}", cfg.out_dir);
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

fn run_star(cfg: &Config, n: usize, final_annot: &Path, cdna_len: usize) -> Result<()> {
    let out = Path::new(&cfg.out_dir);
    let star_dir = out.join("star");
    let threads = (cfg.num_threads / n).max(1);
    let overhang = cdna_len.saturating_sub(1).max(1);

    (0..n).into_par_iter().try_for_each(|i| -> Result<()> {
        let prefix = star_dir.join(format!("{}.{i}.", cfg.project));
        let tmp = Path::new(&cfg.star_tmp).join(format!("{}.star.{i}", cfg.project));
        let _ = std::fs::remove_dir_all(&tmp); // STAR requires outTmpDir to not pre-exist
        let shard = out.join("filtered").join(format!("shard_{i}.bam"));

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
        checked(&mut cmd, &format!("STAR instance {i}"))
    })
}

fn merge_bams(cfg: &Config, n: usize) -> Result<()> {
    let star_dir = Path::new(&cfg.out_dir).join("star");
    let merged = star_dir.join(format!("{}.merged.bam", cfg.project));
    let mut cmd = Command::new("samtools");
    cmd.args(["merge", "-f", "-@", &cfg.num_threads.to_string()]).arg(&merged);
    for i in 0..n {
        cmd.arg(star_dir.join(format!("{}.{i}.Aligned.sortedByCoord.out.bam", cfg.project)));
    }
    checked(&mut cmd, "samtools merge")
}

// FastQC on the sampled filtered reads -> flat qc/.
fn run_fastqc(cfg: &Config) -> Result<()> {
    let qc = Path::new(&cfg.out_dir).join("qc");
    for label in ["R1", "R2"] {
        let fq = Path::new(&cfg.star_tmp).join(format!("qc_{label}.fastq.gz"));
        if !fq.exists() {
            continue;
        }
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
        )?;
        let _ = std::fs::remove_file(&fq);
    }
    Ok(())
}

fn run_summarise(cfg: &Config) -> Result<()> {
    let repo = repo_root()?;
    checked(
        Command::new("Rscript").arg(repo.join("R").join("finalize.R")).arg(&cfg.out_dir).arg(&cfg.project),
        "finalize.R",
    )?;
    // qc_report.R also gets the config path so the report can show the run setup
    checked(
        Command::new("Rscript")
            .arg(repo.join("R").join("qc_report.R"))
            .arg(&cfg.out_dir)
            .arg(&cfg.project)
            .arg(&cfg.config_path),
        "qc_report.R",
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
