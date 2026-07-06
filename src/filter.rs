use crate::bam::ShardSet;
use crate::barcode::{is_tagged, EditBudget, IndexTable};
use crate::config::{Config, Geometry};
use crate::trim::{trim, TrimParams};
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Duration;

const QC_CAP: usize = 1_000_000; // reservoir size of filtered reads sampled for FastQC

#[derive(Serialize, Default)]
pub struct PerBc {
    pub tagged: u64,
    pub internal: u64,
}

#[derive(Serialize, Default)]
pub struct FilterStats {
    pub total: u64,
    pub passed: u64,
    pub dropped_short_r2: u64,
    pub dropped_no_barcode: u64,
    pub dropped_bc_quality: u64,
    pub dropped_umi_quality: u64,
    pub tagged: u64,
    pub internal: u64,
    pub i7_orient: String,
    pub i5_orient: String,
    pub bc_exact: u64,      // both indexes matched exactly
    pub bc_corrected: u64,  // recovered via indel/substitution sphere (net-positive gain)
    pub per_barcode: BTreeMap<String, PerBc>,
}

// Reservoir sample of the aligned cDNA reads, written to star_tmp for FastQC after filtering.
struct QcSampler {
    seen: u64,
    r1: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>, // (name, seq, qual)
    r2: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>,
    rng: StdRng,
}

impl QcSampler {
    fn new() -> QcSampler {
        QcSampler { seen: 0, r1: Vec::new(), r2: Vec::new(), rng: StdRng::seed_from_u64(1) }
    }

    fn offer(&mut self, name: &[u8], r1: (&[u8], &[u8]), r2: (&[u8], &[u8])) {
        let i = self.seen as usize;
        self.seen += 1;
        let mk = |s: &[u8], q: &[u8]| (name.to_vec(), s.to_vec(), q.to_vec());
        if i < QC_CAP {
            self.r1.push(mk(r1.0, r1.1));
            self.r2.push(mk(r2.0, r2.1));
        } else {
            let j = self.rng.gen_range(0..=i);
            if j < QC_CAP {
                self.r1[j] = mk(r1.0, r1.1);
                self.r2[j] = mk(r2.0, r2.1);
            }
        }
    }

    fn write(&self, path: &Path) -> Result<()> {
        for (reads, name) in [(&self.r1, "qc_R1.fastq.gz"), (&self.r2, "qc_R2.fastq.gz")] {
            let f = File::create(path.join(name)).with_context(|| format!("creating {name}"))?;
            let mut w = GzEncoder::new(BufWriter::new(f), Compression::default());
            for (n, s, q) in reads {
                w.write_all(b"@")?;
                w.write_all(n)?;
                w.write_all(b"\n")?;
                w.write_all(s)?;
                w.write_all(b"\n+\n")?;
                w.write_all(q)?;
                w.write_all(b"\n")?;
            }
            w.finish().with_context(|| format!("closing {name}"))?;
        }
        Ok(())
    }
}

fn qname(id: &[u8]) -> &[u8] {
    // first whitespace-delimited token, without a trailing /1 or /2
    let end = id.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(id.len());
    let tok = &id[..end];
    if tok.len() >= 2 && tok[tok.len() - 2] == b'/' {
        &tok[..tok.len() - 2]
    } else {
        tok
    }
}

// count bases below `phred` in a quality slice (ASCII phred+33)
fn low_bases(qual: &[u8], phred: u8) -> usize {
    qual.iter().filter(|&&q| q.saturating_sub(33) < phred).count()
}

const BATCH: usize = 32_768;

// One input quadruplet (all 4 FASTQ mates), bytes owned so the batch can be processed off-reader.
struct Quad {
    tagged_id: Vec<u8>,      // id of the tagged read (only field needed for the qname)
    seqs: Vec<Vec<u8>>,      // 4 sequences, in reader order
    quals: Vec<Vec<u8>>,     // 4 qualities, in reader order
}

enum DropKind {
    ShortR2,
    NoBarcode,
    BcQuality,
    UmiQuality,
}

// Per-read decision produced by the pure, parallel `process`. Drops are split at the barcode
// quality gate: `DropAfterBc` reads already passed it, so the serial drain must still count them
// into bc_exact/bc_corrected before recording the drop — exactly matching the original ordering.
enum Outcome {
    Drop(DropKind),
    DropAfterBc { corrected: bool, kind: DropKind },
    Pass {
        corrected: bool,
        tagged: bool,
        name: Vec<u8>,
        r1_seq: Vec<u8>,
        r1_qual: Vec<u8>,
        r2_seq: Vec<u8>,
        r2_qual: Vec<u8>,
        bc: String,
        umi: String,
    },
}

// Pure per-read classification over immutable, shared inputs. Called in parallel; performs NO
// counter mutation, sharding, RNG draw, or I/O — those all happen in the serial drain in input order
// so the shards, QC reservoir and stats are byte-identical regardless of thread count.
fn process(q: &Quad, table: &IndexTable, params: &TrimParams, g: &Geometry, cfg: &Config) -> Outcome {
    let seqs = &q.seqs;
    let quals = &q.quals;

    // ---- R2 read-filtering FIRST: adapter clip + quality trim + min_length ----
    let r2_seq = &seqs[g.internal_file];
    let r2_qual = &quals[g.internal_file];
    let (s, e) = match trim(r2_seq, r2_qual, params) {
        Some(v) => v,
        None => return Outcome::Drop(DropKind::ShortR2),
    };
    let r2_cdna = &r2_seq[s..e];
    let r2_cdna_q = &r2_qual[s..e];

    // ---- cell barcode: i7 (bc[0], I1) and i5 (bc[1], I2) corrected separately ----
    let s0 = &g.bc[0];
    let s1 = &g.bc[1];
    if seqs[s0.file].len() < s0.end || seqs[s1.file].len() < s1.end {
        return Outcome::Drop(DropKind::NoBarcode);
    }
    let raw_i7 = &seqs[s0.file][s0.start..s0.end];
    let raw_i5 = &seqs[s1.file][s1.start..s1.end];
    // concatenated quality preserves the existing BC quality-gate semantics
    let mut raw_bc_q = Vec::with_capacity((s0.end - s0.start) + (s1.end - s1.start));
    raw_bc_q.extend_from_slice(&quals[s0.file][s0.start..s0.end]);
    raw_bc_q.extend_from_slice(&quals[s1.file][s1.start..s1.end]);

    let (bc, corrected) = match table.assign_pair(raw_i7, raw_i5) {
        Some(v) => v,
        None => return Outcome::Drop(DropKind::NoBarcode),
    };
    if low_bases(&raw_bc_q, cfg.filter_cutoffs.bc_filter.phred) > cfg.filter_cutoffs.bc_filter.num_bases {
        return Outcome::Drop(DropKind::BcQuality);
    }

    // ---- tagged/internal split + UMI from the tagged read (R1) ----
    let t_seq = &seqs[g.tagged_file];
    let t_qual = &quals[g.tagged_file];
    let tagged = is_tagged(t_seq, &g.tag, g.tag_mismatch);
    let (r1_cdna, r1_cdna_q, umi): (Vec<u8>, Vec<u8>, String) = if tagged {
        if t_seq.len() < g.umi.end || t_seq.len() < g.cdna_start {
            return Outcome::DropAfterBc { corrected, kind: DropKind::ShortR2 };
        }
        let umi_q = &t_qual[g.umi.start..g.umi.end];
        if low_bases(umi_q, cfg.filter_cutoffs.umi_filter.phred) > cfg.filter_cutoffs.umi_filter.num_bases {
            return Outcome::DropAfterBc { corrected, kind: DropKind::UmiQuality };
        }
        let umi = String::from_utf8_lossy(&t_seq[g.umi.start..g.umi.end]).into_owned();
        (t_seq[g.cdna_start..].to_vec(), t_qual[g.cdna_start..].to_vec(), umi)
    } else {
        (t_seq.to_vec(), t_qual.to_vec(), String::new())
    };

    Outcome::Pass {
        corrected,
        tagged,
        name: qname(&q.tagged_id).to_vec(),
        r1_seq: r1_cdna,
        r1_qual: r1_cdna_q,
        r2_seq: r2_cdna.to_vec(),
        r2_qual: r2_cdna_q.to_vec(),
        bc,
        umi,
    }
}

pub fn filter(cfg: &Config, n_shards: usize, stage: &crate::log::Stage) -> Result<FilterStats> {
    let out = Path::new(&cfg.out_dir);
    let filtered_dir = out.join("filtered");
    std::fs::create_dir_all(&filtered_dir).context("creating filtered/")?;

    let budget = EditBudget { max_total: cfg.barcodes.barcode_binning as u32, max_indel: 1 };
    let mut table = IndexTable::load(&cfg.barcodes.index_table, budget)?;
    for w in &table.warnings {
        stage.step(format!("panel-guard: {w}"));
    }
    let params = TrimParams::load(&cfg.read_filtering)?;
    let g: Geometry = cfg.geometry()?;

    // ---- orientation pre-pass: sample the first N records of I1 (i7) and I2 (i5) ----
    const ORIENT_SAMPLE: usize = 50_000;
    let (s7, e7) = (g.bc[0].start, g.bc[0].end);
    let (s5, e5) = (g.bc[1].start, g.bc[1].end);
    let mut i7_samples: Vec<Vec<u8>> = Vec::new();
    let mut i5_samples: Vec<Vec<u8>> = Vec::new();
    {
        let mut r7 = needletail::parse_fastx_file(&g.files[g.bc[0].file])
            .with_context(|| format!("opening {}", g.files[g.bc[0].file]))?;
        while i7_samples.len() < ORIENT_SAMPLE {
            match r7.next() {
                Some(rec) => {
                    let rec = rec?;
                    let s = rec.seq();
                    if s.len() >= e7 {
                        i7_samples.push(s[s7..e7].to_vec());
                    }
                }
                None => break,
            }
        }
        let mut r5 = needletail::parse_fastx_file(&g.files[g.bc[1].file])
            .with_context(|| format!("opening {}", g.files[g.bc[1].file]))?;
        while i5_samples.len() < ORIENT_SAMPLE {
            match r5.next() {
                Some(rec) => {
                    let rec = rec?;
                    let s = rec.seq();
                    if s.len() >= e5 {
                        i5_samples.push(s[s5..e5].to_vec());
                    }
                }
                None => break,
            }
        }
    }
    table.detect_and_set_orientation(&i7_samples, &i5_samples)?;
    stage.step(format!("orientation i7={:?} i5={:?}", table.i7_orient, table.i5_orient));

    let mut shards = ShardSet::create(&filtered_dir, n_shards, "quartsx filter")?;
    let mut qc = QcSampler::new();
    let mut stats = FilterStats::default();
    stats.i7_orient = format!("{:?}", table.i7_orient);
    stats.i5_orient = format!("{:?}", table.i5_orient);

    let mut readers: Vec<_> = g
        .files
        .iter()
        .map(|f| needletail::parse_fastx_file(f).with_context(|| format!("opening {f}")))
        .collect::<Result<_>>()?;

    stage.step(format!("reading FASTQs -> {n_shards} shard(s)"));
    // Read a batch of quadruplets serially (owning the bytes), classify the batch in parallel with
    // rayon (order-preserving par_iter), then drain the ordered outcomes serially. The serial drain
    // keeps every order-sensitive action — passed-rank shard routing, the seed-1 QC reservoir RNG,
    // and all counters — in exact input order, so shards + filter_stats.json are byte-identical to
    // the single-threaded pipeline at any thread count. One in-flight batch bounds memory.
    let mut batch: Vec<Quad> = Vec::with_capacity(BATCH);
    let mut eof = false;
    loop {
        batch.clear();
        while batch.len() < BATCH {
            let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(4);
            let mut quals: Vec<Vec<u8>> = Vec::with_capacity(4);
            let mut tagged_id: Vec<u8> = Vec::new();
            let mut ended = 0;
            for (fi, r) in readers.iter_mut().enumerate() {
                match r.next() {
                    Some(rec) => {
                        let rec = rec?;
                        if fi == g.tagged_file {
                            tagged_id = rec.id().to_vec();
                        }
                        seqs.push(rec.seq().to_vec());
                        quals.push(rec.qual().context("fastq missing qualities")?.to_vec());
                    }
                    None => {
                        ended += 1;
                        seqs.push(Vec::new());
                        quals.push(Vec::new());
                    }
                }
            }
            if ended == 4 {
                eof = true;
                break;
            }
            if ended != 0 {
                bail!("input fastqs have different lengths");
            }
            batch.push(Quad { tagged_id, seqs, quals });
        }
        if batch.is_empty() {
            if eof {
                break;
            } else {
                continue;
            }
        }

        let outcomes: Vec<Outcome> = batch.par_iter().map(|q| process(q, &table, &params, &g, cfg)).collect();

        for oc in outcomes {
            stats.total += 1;
            match oc {
                Outcome::Drop(kind) => match kind {
                    DropKind::ShortR2 => stats.dropped_short_r2 += 1,
                    DropKind::NoBarcode => stats.dropped_no_barcode += 1,
                    DropKind::BcQuality => stats.dropped_bc_quality += 1,
                    DropKind::UmiQuality => stats.dropped_umi_quality += 1,
                },
                Outcome::DropAfterBc { corrected, kind } => {
                    if corrected {
                        stats.bc_corrected += 1;
                    } else {
                        stats.bc_exact += 1;
                    }
                    match kind {
                        DropKind::ShortR2 => stats.dropped_short_r2 += 1,
                        DropKind::UmiQuality => stats.dropped_umi_quality += 1,
                        _ => {}
                    }
                }
                Outcome::Pass { corrected, tagged, name, r1_seq, r1_qual, r2_seq, r2_qual, bc, umi } => {
                    if corrected {
                        stats.bc_corrected += 1;
                    } else {
                        stats.bc_exact += 1;
                    }
                    let shard = (stats.passed as usize) % n_shards;
                    shards.write_pair(shard, &name, &r1_seq, &r1_qual, &r2_seq, &r2_qual, &bc, &umi)?;
                    qc.offer(&name, (&r1_seq, &r1_qual), (&r2_seq, &r2_qual));
                    stats.passed += 1;
                    let e = stats.per_barcode.entry(bc).or_default();
                    if tagged {
                        stats.tagged += 1;
                        e.tagged += 1;
                    } else {
                        stats.internal += 1;
                        e.internal += 1;
                    }
                }
            }
        }

        let (total, passed) = (stats.total, stats.passed);
        stage.beat(Duration::from_secs(5), || {
            format!("{:.1}M reads, {:.1}M passed", total as f64 / 1e6, passed as f64 / 1e6)
        });
        if eof {
            break;
        }
    }

    shards.finish()?;
    qc.write(Path::new(&cfg.star_tmp))?;

    let json = serde_json::to_string_pretty(&stats).context("serializing filter stats")?;
    std::fs::write(filtered_dir.join("filter_stats.json"), json).context("writing filter_stats.json")?;
    Ok(stats)
}
