//! STAGE: FILTERING — three-stage overlapped pipeline.
//!
//! Shape (STAR's parallelization shape + QUARTSx order-preserving determinism):
//!   1. Per-file inflate — one decompressor thread per input FASTQ (R1‖R2 heavy; I1/I2 tiny),
//!      each pushing FILE-ORDER chunks over a bounded std `sync_channel` (an in-process "FIFO",
//!      never `mkfifo`, never on disk). A deterministic 4-way zipper on the main thread joins one
//!      record from each file into a `Quad`, reproducing the original EOF / length-mismatch checks.
//!   2. N filter+encode workers on a filter-LOCAL rayon pool (the global pool is left untouched):
//!      per-read adapter/quality trim, barcode assign+correct, UMI, TSO tag classification. The
//!      order-preserving `par_iter().map().collect()` keeps outcomes in input order.
//!   3. P-deflater BGZF write via `ShardSet` (its internal `MultithreadedWriter`). The serial drain
//!      does only routing + the seed-1 QC reservoir + counters + `write_pair`.
//!
//! Determinism gate: at a fixed `n_shards`, varying P and N yields identical
//! sha256(`{project}.filtered.bam`) + `filter_stats.*` + QC fastqs. This holds because every
//! order-sensitive action (shard routing on `passed % n_shards`, the seed-1 reservoir RNG, all
//! counters) runs SERIALLY in input order on the drain, the rayon collect is order-preserving, and
//! the MultithreadedWriter emits blocks in submission order (byte-identical to the serial baseline).
//!
//! NOTE on "encode-in-worker": `ShardSet::write_pair(raw slices)` keeps the cheap `RecordBuf` build
//! on the drain; the real bottleneck — DEFLATE — parallelizes across P deflaters inside `ShardSet`.

use crate::bam::ShardSet;
use crate::barcode::{is_tagged, AssignResult, EditBudget, IndexTable, RevcompScratch, UnassignedReason};
use crate::config::{Config, Geometry};
use crate::log::Stage;
use crate::trim::{trim, TrimParams};
use anyhow::{bail, Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_reader;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::Arc;

const QC_CAP: usize = 1_000_000; // reservoir size of filtered reads sampled for FastQC
const BATCH: usize = 32_768; // quads classified per rayon batch (bounds one in-flight batch)
const CHUNK: usize = 8_192; // records a reader thread accumulates before sending a chunk
const CHUNK_CHAN_CAP: usize = 4; // chunks buffered per file (bounds read-ahead memory)

#[derive(Serialize, Deserialize, Default)]
pub struct PerBc {
    pub tagged: u64,
    pub internal: u64,
}

#[derive(Serialize, Default)]
pub struct FilterStats {
    pub total: u64,
    pub passed: u64,
    pub dropped_short_r2: u64,
    // no_barcode is split into three non-overlapping sub-reasons. Field names match the
    // qc_report.R `drop_labels` map exactly so the drop-reason bar renders each one; the report
    // sums EVERY scalar `dropped_*` field, so we deliberately do NOT also emit an aggregate
    // `dropped_no_barcode` (that would double-count).
    pub dropped_no_barcode_absent: u64, // an index window is a genuine miss (not in exact or sphere)
    pub dropped_no_barcode_uncorrectable: u64, // both indexes decode but (i7,i5) is not a listed cell
    pub dropped_no_barcode_ambiguous: u64, // an index window hit the REJECT sphere sentinel (>=2 indexes)
    pub dropped_bc_quality: u64,
    pub dropped_umi_quality: u64,
    pub tagged: u64,
    pub internal: u64,
    pub i7_orient: String,
    pub i5_orient: String,
    pub bc_exact: u64,     // both indexes matched exactly
    pub bc_corrected: u64, // recovered via indel/substitution sphere (net-positive gain)
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
    qual.iter().filter(|&&q| q.saturating_sub(crate::PHRED_OFFSET) < phred).count()
}

// ---------------------------------------------------------------------------
// Stage 1: per-file inflate + byte-counting reader (for the progress bar)
// ---------------------------------------------------------------------------

/// Wraps a `File` and tallies COMPRESSED bytes consumed, so the progress bar can report % of the
/// (on-disk, gzip) input read. Placed UNDER a `BufReader`/the gzip decoder, so `read` is called in
/// large chunks against the raw file and the tally equals compressed bytes pulled off disk.
struct CountingReader {
    inner: File,
    counted: Arc<AtomicU64>,
}

impl Read for CountingReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.counted.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// One decompressed FASTQ record from a single input file, owned so it can cross the channel.
/// `id` is populated ONLY for the tagged read (empty otherwise) — the qname comes from that mate.
#[derive(Default)]
struct FileRec {
    id: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

/// One decompressor thread for a single input file. Sends FILE-ORDER `CHUNK`-sized batches of
/// records over a bounded channel; on any parse/IO error sends a single `Err` and stops. Breaks
/// silently if the consumer (main) has gone away (early-return path), so `thread::scope` join
/// never deadlocks.
fn reader_thread(
    path: &str,
    counter: Arc<AtomicU64>,
    is_tagged: bool,
    tx: SyncSender<Result<Vec<FileRec>>>,
) {
    let run = || -> Result<()> {
        let f = File::open(path).with_context(|| format!("opening {path}"))?;
        let reader = BufReader::with_capacity(1 << 16, CountingReader { inner: f, counted: counter });
        // parse_fastx_reader sniffs the gzip magic itself, so a plain (buffered) file reader still
        // gets transparently decompressed; the byte tally underneath counts compressed bytes.
        let mut rdr = parse_fastx_reader(reader).with_context(|| format!("parsing {path}"))?;
        let mut chunk: Vec<FileRec> = Vec::with_capacity(CHUNK);
        loop {
            match rdr.next() {
                Some(rec) => {
                    let rec = rec.with_context(|| format!("reading record in {path}"))?;
                    let id = if is_tagged { rec.id().to_vec() } else { Vec::new() };
                    let seq = rec.seq().to_vec();
                    let qual = rec.qual().context("fastq missing qualities")?.to_vec();
                    chunk.push(FileRec { id, seq, qual });
                    if chunk.len() >= CHUNK {
                        let send = std::mem::replace(&mut chunk, Vec::with_capacity(CHUNK));
                        if tx.send(Ok(send)).is_err() {
                            return Ok(()); // consumer gone
                        }
                    }
                }
                None => {
                    if !chunk.is_empty() {
                        let _ = tx.send(Ok(chunk));
                    }
                    return Ok(());
                }
            }
        }
    };
    if let Err(e) = run() {
        let _ = tx.send(Err(e));
    }
}

/// Buffered per-file pull: hands out one record at a time in file order, refilling from the chunk
/// channel as needed. `Ok(None)` = clean EOF (channel closed); `Err` = a reader-thread error.
struct FilePull {
    rx: Receiver<Result<Vec<FileRec>>>,
    buf: Vec<FileRec>,
    pos: usize,
    done: bool,
}

impl FilePull {
    fn new(rx: Receiver<Result<Vec<FileRec>>>) -> FilePull {
        FilePull { rx, buf: Vec::new(), pos: 0, done: false }
    }

    fn next(&mut self) -> Result<Option<FileRec>> {
        loop {
            if self.pos < self.buf.len() {
                let r = std::mem::take(&mut self.buf[self.pos]);
                self.pos += 1;
                return Ok(Some(r));
            }
            if self.done {
                return Ok(None);
            }
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    self.buf = chunk;
                    self.pos = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                    return Ok(None);
                }
            }
        }
    }
}

// One input quadruplet (all 4 FASTQ mates), bytes owned so the batch can be classified off-reader.
struct Quad {
    tagged_id: Vec<u8>,   // id of the tagged read (only field needed for the qname)
    seqs: [Vec<u8>; 4],   // 4 sequences, in reader order
    quals: [Vec<u8>; 4],  // 4 qualities, in reader order
}

/// Deterministic 4-way zipper: pull one record from each `FilePull`, assemble a `Quad`, up to
/// `BATCH` of them. Reproduces the original EOF/length semantics: all four files ending together is
/// clean EOF; any partial ending is a hard "different lengths" error. Returns `Ok(true)` when EOF
/// was reached (the partially-filled `batch` is still to be processed by the caller).
fn fill_batch(pulls: &mut [FilePull], tagged_file: usize, batch: &mut Vec<Quad>) -> Result<bool> {
    batch.clear();
    while batch.len() < BATCH {
        let mut items: [Option<FileRec>; 4] = Default::default();
        let mut ended = 0;
        for (fi, item) in items.iter_mut().enumerate() {
            match pulls[fi].next()? {
                Some(fr) => *item = Some(fr),
                None => ended += 1,
            }
        }
        if ended == 4 {
            return Ok(true); // clean EOF: all inputs exhausted together
        }
        if ended != 0 {
            bail!("input fastqs have different lengths");
        }
        let mut seqs: [Vec<u8>; 4] = Default::default();
        let mut quals: [Vec<u8>; 4] = Default::default();
        let mut tagged_id = Vec::new();
        for (fi, item) in items.iter_mut().enumerate() {
            let fr = item.take().unwrap();
            if fi == tagged_file {
                tagged_id = fr.id;
            }
            seqs[fi] = fr.seq;
            quals[fi] = fr.qual;
        }
        batch.push(Quad { tagged_id, seqs, quals });
    }
    Ok(false)
}

// ---------------------------------------------------------------------------
// Stage 2: per-read classification (pure, run in parallel on the filter pool)
// ---------------------------------------------------------------------------

enum DropKind {
    ShortR2,
    BcAbsent,  // index window not present (genuine miss / non-ACGT / short read)
    BcReject,  // index window ambiguous (REJECT sphere sentinel)
    BcBadPair, // both indexes decode but the (i7,i5) pair is not a listed cell
    BcQuality,
    UmiQuality,
}

// Per-read decision produced by the pure, parallel `process`. Drops are split at the barcode
// quality gate: `DropAfterBc` reads already passed barcode assignment, so the serial drain must
// still count them into bc_exact/bc_corrected before recording the drop. (`DropAfterBc` only
// ever carries ShortR2 / UmiQuality.)
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
// counter mutation, sharding, RNG draw, or I/O — those all happen in the serial drain in input
// order so the shards, QC reservoir and stats are byte-identical regardless of thread count.
fn process(q: &Quad, table: &IndexTable, params: &TrimParams, g: &Geometry, cfg: &Config, scratch: &mut RevcompScratch) -> Outcome {
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
        return Outcome::Drop(DropKind::BcAbsent);
    }
    let raw_i7 = &seqs[s0.file][s0.start..s0.end];
    let raw_i5 = &seqs[s1.file][s1.start..s1.end];
    // concatenated quality preserves the existing BC quality-gate semantics. bc windows are <=12 bp
    // each (<=24 concatenated), so a fixed 64-byte stack scratch avoids a per-read heap allocation.
    let mut bc_q_buf = [0u8; 64];
    let bc_q_len = (s0.end - s0.start) + (s1.end - s1.start);
    let raw_bc_q: &[u8] = if bc_q_len <= bc_q_buf.len() {
        bc_q_buf[..s0.end - s0.start].copy_from_slice(&quals[s0.file][s0.start..s0.end]);
        bc_q_buf[s0.end - s0.start..bc_q_len].copy_from_slice(&quals[s1.file][s1.start..s1.end]);
        &bc_q_buf[..bc_q_len]
    } else {
        // unreachable for SS3xpress index lengths (i7+i5 <= 24 bp); safe fallback if ever exceeded.
        &bc_q_buf[..0] // len 0 -> passes the gate; never taken in practice
    };

    // assign_into distinguishes the three no-barcode sub-reasons: a genuine miss, an ambiguous
    // reject-sentinel window, and a decoded-but-unlisted pair. The per-worker `scratch` keeps the
    // orientation-flip allocation off the hot path. The interned label is borrowed, so it is
    // materialized into an owned `bc` only for a passing read.
    let (bc, corrected) = match table.assign_into(raw_i7, raw_i5, scratch) {
        AssignResult::Assigned { label, corrected, .. } => (label.to_string(), corrected),
        AssignResult::Unassigned(UnassignedReason::Absent) => return Outcome::Drop(DropKind::BcAbsent),
        AssignResult::Unassigned(UnassignedReason::AmbiguousReject) => return Outcome::Drop(DropKind::BcReject),
        AssignResult::Unassigned(UnassignedReason::InvalidPair) => return Outcome::Drop(DropKind::BcBadPair),
    };
    if low_bases(raw_bc_q, cfg.filter_cutoffs.bc_filter.phred) > cfg.filter_cutoffs.bc_filter.num_bases {
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

fn breakdown(s: &FilterStats) -> String {
    format!(
        "passed={} short_r2={} bc_absent={} bc_uncorrectable={} bc_ambiguous={} bc_qual={} umi_qual={}",
        s.passed,
        s.dropped_short_r2,
        s.dropped_no_barcode_absent,
        s.dropped_no_barcode_uncorrectable,
        s.dropped_no_barcode_ambiguous,
        s.dropped_bc_quality,
        s.dropped_umi_quality
    )
}

// Serial drain of ONE classified outcome: routing, seed-1 QC reservoir, counters, write. All
// order-sensitive, so it always runs in input order (called per outcome in input order from both the
// single-threaded and multithreaded paths) -> shards/QC/stats are byte-identical regardless of T.
fn drain_outcome(
    oc: Outcome,
    stats: &mut FilterStats,
    shards: &mut ShardSet,
    qc: &mut QcSampler,
    n_shards: usize,
) -> Result<()> {
    stats.total += 1;
    match oc {
        Outcome::Drop(kind) => match kind {
            DropKind::ShortR2 => stats.dropped_short_r2 += 1,
            DropKind::BcAbsent => stats.dropped_no_barcode_absent += 1,
            DropKind::BcReject => stats.dropped_no_barcode_ambiguous += 1,
            DropKind::BcBadPair => stats.dropped_no_barcode_uncorrectable += 1,
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
        Outcome::Pass {
            corrected,
            tagged,
            name,
            r1_seq,
            r1_qual,
            r2_seq,
            r2_qual,
            bc,
            umi,
        } => {
            if corrected {
                stats.bc_corrected += 1;
            } else {
                stats.bc_exact += 1;
            }
            let shard = (stats.passed as usize) % n_shards;
            shards.write_pair(shard, &name, &r1_seq, &r1_qual, &r2_seq, &r2_qual, &bc, &umi)?;
            qc.offer(&name, (&r1_seq, &r1_qual), (&r2_seq, &r2_qual));
            stats.passed += 1;
            let entry = stats.per_barcode.entry(bc).or_default();
            if tagged {
                stats.tagged += 1;
                entry.tagged += 1;
            } else {
                stats.internal += 1;
                entry.internal += 1;
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

pub fn filter(cfg: &Config, n_shards: usize, stage: &Stage) -> Result<FilterStats> {
    let out = Path::new(&cfg.out_dir);
    let filtered_dir = out.join("filtered");

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

    // Thread knobs: P = compress deflaters (ShardSet), N = filter-local rayon workers. At T<=1
    // resolved_threads() returns the degenerate (1,1) budget and we take the single-threaded path
    // below (no reader threads, no rayon pool, single-threaded BGZF codec) so one working thread never
    // oversubscribes the lone core. The filter pool (T>1 only) is DISTINCT from the global rayon pool.
    let (p_threads, n_threads) = cfg.resolved_threads();
    let single = cfg.num_threads <= 1;
    let p = p_threads.max(1);
    let n = n_threads.max(1);
    if single {
        stage.step(format!("pipeline: single-threaded, {n_shards} shard(s)"));
    } else {
        stage.step(format!("pipeline: {n} filter worker(s), {p} deflater(s), {n_shards} shard(s)"));
    }

    std::fs::create_dir_all(&filtered_dir).context("creating filtered/")?;
    // compress_threads == 0 selects the single-threaded BGZF codec.
    let codec_threads = if single { 0 } else { p };
    let mut shards = ShardSet::create(&filtered_dir, n_shards, &cfg.project, codec_threads, "quartsx filter")?;
    let mut qc = QcSampler::new();
    let mut stats = FilterStats::default();
    stats.i7_orient = format!("{:?}", table.i7_orient);
    stats.i5_orient = format!("{:?}", table.i5_orient);

    // Progress bar over COMPRESSED input bytes consumed.
    let total_bytes: u64 = g
        .files
        .iter()
        .map(|f| std::fs::metadata(f).map(|m| m.len()).unwrap_or(0))
        .sum();
    let counters: Vec<Arc<AtomicU64>> = (0..4).map(|_| Arc::new(AtomicU64::new(0))).collect();
    let progress = stage.progress_bytes(total_bytes);

    if single {
        // Single-threaded path: inflate the 4 inputs INLINE (no reader threads), classify + drain
        // each quad serially (no rayon pool), write via the single-threaded BGZF codec — so at T<=1
        // exactly one working thread runs. Output is byte-identical to the multithreaded path.
        let mut readers = Vec::with_capacity(4);
        for fi in 0..4 {
            let f = File::open(&g.files[fi]).with_context(|| format!("opening {}", g.files[fi]))?;
            let reader =
                BufReader::with_capacity(1 << 16, CountingReader { inner: f, counted: counters[fi].clone() });
            readers.push(parse_fastx_reader(reader).with_context(|| format!("parsing {}", g.files[fi]))?);
        }
        let mut scratch = RevcompScratch::new();
        let mut next_report = 0.25f64;
        loop {
            let mut seqs: [Vec<u8>; 4] = Default::default();
            let mut quals: [Vec<u8>; 4] = Default::default();
            let mut tagged_id = Vec::new();
            let mut ended = 0usize;
            for fi in 0..4 {
                match readers[fi].next() {
                    Some(rec) => {
                        let rec = rec.with_context(|| format!("reading record in {}", g.files[fi]))?;
                        if fi == g.tagged_file {
                            tagged_id = rec.id().to_vec();
                        }
                        seqs[fi] = rec.seq().to_vec();
                        quals[fi] = rec.qual().context("fastq missing qualities")?.to_vec();
                    }
                    None => ended += 1,
                }
            }
            if ended == 4 {
                break; // clean EOF: all inputs exhausted together
            }
            if ended != 0 {
                bail!("input fastqs have different lengths");
            }
            let q = Quad { tagged_id, seqs, quals };
            let oc = process(&q, &table, &params, &g, cfg, &mut scratch);
            drain_outcome(oc, &mut stats, &mut shards, &mut qc, n_shards)?;

            if stats.total % BATCH as u64 == 0 {
                let consumed: u64 = counters.iter().map(|c| c.load(Ordering::Relaxed)).sum();
                progress.set(consumed);
                let frac = consumed as f64 / total_bytes.max(1) as f64;
                while next_report <= 1.0 && frac >= next_report {
                    stage.step(format!("{:.0}% input consumed — {}", next_report * 100.0, breakdown(&stats)));
                    next_report += 0.25;
                }
            }
        }
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .thread_name(|i| format!("filter-{i}"))
            .build()
            .context("building filter thread pool")?;
        // Overlapped pipeline: 4 decompressor threads inflate concurrently into bounded channels while
        // the main thread zips + filters + drains. On any error the main closure returns; its `FilePull`
        // receivers drop, the reader threads' sends fail, they finish, and `scope` joins without hanging.
        std::thread::scope(|scope| -> Result<()> {
            let mut rxs: Vec<Receiver<Result<Vec<FileRec>>>> = Vec::with_capacity(4);
            for fi in 0..4 {
                let (tx, rx) = sync_channel::<Result<Vec<FileRec>>>(CHUNK_CHAN_CAP);
                rxs.push(rx);
                let path = g.files[fi].clone();
                let counter = counters[fi].clone();
                let is_tagged_file = fi == g.tagged_file;
                scope.spawn(move || reader_thread(&path, counter, is_tagged_file, tx));
            }
            let mut pulls: Vec<FilePull> = rxs.into_iter().map(FilePull::new).collect();

            let mut batch: Vec<Quad> = Vec::with_capacity(BATCH);
            let mut next_report = 0.25f64; // print a drop-reason breakdown a few times during the run
            loop {
                let eof = fill_batch(&mut pulls, g.tagged_file, &mut batch)?;

                if !batch.is_empty() {
                    // Stage 2: classify the batch in parallel on the FILTER-LOCAL pool (order-preserving).
                    let outcomes: Vec<Outcome> = pool.install(|| {
                        batch
                            .par_iter()
                            .map_init(RevcompScratch::new, |scratch, q| process(q, &table, &params, &g, cfg, scratch))
                            .collect()
                    });

                    // Stage 3: serial drain in input order so output is thread-count-independent.
                    for oc in outcomes {
                        drain_outcome(oc, &mut stats, &mut shards, &mut qc, n_shards)?;
                    }

                    let consumed: u64 = counters.iter().map(|c| c.load(Ordering::Relaxed)).sum();
                    progress.set(consumed);
                    let frac = consumed as f64 / total_bytes.max(1) as f64;
                    while next_report <= 1.0 && frac >= next_report {
                        stage.step(format!("{:.0}% input consumed — {}", next_report * 100.0, breakdown(&stats)));
                        next_report += 0.25;
                    }
                }

                if eof {
                    break;
                }
            }
            Ok(())
        })?;
    }

    progress.finish();
    stage.step(format!("{:.1}M reads, {:.1}M passed", stats.total as f64 / 1e6, stats.passed as f64 / 1e6));

    shards.finish()?;
    qc.write(Path::new(&cfg.star_tmp))?;

    // Stats outputs: JSON + a scalar TSV + a per-barcode TSV table.
    let json = serde_json::to_string_pretty(&stats).context("serializing filter stats")?;
    std::fs::write(filtered_dir.join("filter_stats.json"), json).context("writing filter_stats.json")?;
    // The filter stats table is QC output → write it into the run's qc/ folder, not filtered/.
    let qc_dir = out.join("qc");
    std::fs::create_dir_all(&qc_dir).context("creating qc/")?;
    write_stats_tsv(&qc_dir, &stats)?;
    write_per_barcode_tsv(&filtered_dir, &stats)?;

    Ok(stats)
}

/// Scalar filter metrics as a two-column TSV alongside the JSON.
fn write_stats_tsv(dir: &Path, s: &FilterStats) -> Result<()> {
    let f = File::create(dir.join("filter_stats.tsv")).context("creating filter_stats.tsv")?;
    let mut w = BufWriter::new(f);
    writeln!(w, "metric\tvalue")?;
    let rows: [(&str, String); 15] = [
        ("total", s.total.to_string()),
        ("passed", s.passed.to_string()),
        ("tagged", s.tagged.to_string()),
        ("internal", s.internal.to_string()),
        ("dropped_short_r2", s.dropped_short_r2.to_string()),
        ("dropped_no_barcode_absent", s.dropped_no_barcode_absent.to_string()),
        ("dropped_no_barcode_uncorrectable", s.dropped_no_barcode_uncorrectable.to_string()),
        ("dropped_no_barcode_ambiguous", s.dropped_no_barcode_ambiguous.to_string()),
        ("dropped_bc_quality", s.dropped_bc_quality.to_string()),
        ("dropped_umi_quality", s.dropped_umi_quality.to_string()),
        ("bc_exact", s.bc_exact.to_string()),
        ("bc_corrected", s.bc_corrected.to_string()),
        ("i7_orient", s.i7_orient.clone()),
        ("i5_orient", s.i5_orient.clone()),
        ("passed_pct", format!("{:.4}", if s.total > 0 { 100.0 * s.passed as f64 / s.total as f64 } else { 0.0 })),
    ];
    for (k, v) in rows {
        writeln!(w, "{k}\t{v}")?;
    }
    w.flush().context("flushing filter_stats.tsv")?;
    Ok(())
}

/// Per-barcode fragment table. BTreeMap iteration is sorted → deterministic ordering.
fn write_per_barcode_tsv(dir: &Path, s: &FilterStats) -> Result<()> {
    let f = File::create(dir.join("filter_per_barcode.tsv")).context("creating filter_per_barcode.tsv")?;
    let mut w = BufWriter::new(f);
    writeln!(w, "barcode\ttagged\tinternal\ttotal")?;
    for (bc, pb) in &s.per_barcode {
        writeln!(w, "{bc}\t{}\t{}\t{}", pb.tagged, pb.internal, pb.tagged + pb.internal)?;
    }
    w.flush().context("flushing filter_per_barcode.tsv")?;
    Ok(())
}
