// All noodles BAM/SAM I/O lives here so the (fast-moving) noodles API surface is in one place.
//
// BGZF codec: the on-disk checkpoints (`{project}.filtered.bam`,
// `{project}.tagged.bam`) are written through `BgzfW` and checkpoints are read back through `BgzfR`.
// Each is multithreaded for T>1 and single-threaded for T<=1 (a `threads` argument of 0 selects the
// single-threaded branch — no worker/write/reader threads). Every branch keeps noodles' DEFAULT
// codec (level-6), so the compressed bytes are IDENTICAL regardless of branch/worker count: BGZF
// block boundaries depend only on `MAX_BUF_SIZE` and the byte stream, and the multithreaded writer
// re-orders deflated blocks back into submission order. Result: byte-identical checkpoints for any T.
use anyhow::{Context, Result};
use noodles::bam;
use noodles::bgzf;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::Header;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

/// Clamp a thread count to a valid worker count (>= 1) for the multithreaded BGZF codecs.
fn worker_count(n: usize) -> NonZeroUsize {
    // Safe: max(1) guarantees a non-zero value.
    NonZeroUsize::new(n.max(1)).expect("n.max(1) >= 1")
}

/// BGZF write layer: multithreaded (T>1) or single-threaded (T<=1). `create` picks the branch from
/// the deflate-worker budget — `threads == 0` => single-threaded (`bgzf::Writer`, no threads).
enum BgzfW {
    Mt(bgzf::MultithreadedWriter<BufWriter<File>>),
    St(bgzf::Writer<BufWriter<File>>),
}

impl BgzfW {
    fn create(file: File, threads: usize) -> BgzfW {
        let inner = BufWriter::new(file);
        if threads == 0 {
            BgzfW::St(bgzf::Writer::new(inner))
        } else {
            BgzfW::Mt(bgzf::MultithreadedWriter::with_worker_count(worker_count(threads), inner))
        }
    }

    /// Append the BGZF EOF block (joining the deflater/writer threads for the multithreaded codec)
    /// and flush the underlying file buffer.
    fn finish(self) -> Result<()> {
        let mut inner = match self {
            BgzfW::Mt(mut w) => w.finish().context("finalizing bgzf stream")?,
            BgzfW::St(w) => w.finish().context("finalizing bgzf stream")?,
        };
        inner.flush().context("flushing bgzf file buffer")?;
        Ok(())
    }
}

impl Write for BgzfW {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            BgzfW::Mt(w) => w.write(buf),
            BgzfW::St(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            BgzfW::Mt(w) => w.flush(),
            BgzfW::St(w) => w.flush(),
        }
    }
}

/// BGZF read layer: multithreaded (T>1) or single-threaded (T<=1). `threads == 0` => single-threaded
/// (`bgzf::Reader`, no reader thread). Both serve decoded bytes in FILE ORDER, so the record stream
/// is identical regardless of branch/worker count.
enum BgzfR {
    Mt(bgzf::MultithreadedReader<File>),
    St(bgzf::Reader<File>),
}

impl BgzfR {
    fn open(file: File, threads: usize) -> BgzfR {
        if threads == 0 {
            BgzfR::St(bgzf::Reader::new(file))
        } else {
            BgzfR::Mt(bgzf::MultithreadedReader::with_worker_count(worker_count(threads), file))
        }
    }
}

impl Read for BgzfR {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        match self {
            BgzfR::Mt(r) => r.read(buf),
            BgzfR::St(r) => r.read(buf),
        }
    }
}

/// SINGLE authority for the filtered-shard checkpoint name. The orchestrator and the filter
/// stage both call this — NEVER `format!` a shard path inline.
///   n == 1  -> dir/{project}.filtered.bam        (COMMON path: no numeric suffix, no `shard_` prefix)
///   n  > 1  -> dir/{project}.filtered.{i}.bam     (internal-only, when N parallel STAR jobs are used)
pub fn filtered_shard_path(dir: &Path, project: &str, i: usize, n: usize) -> PathBuf {
    if n == 1 {
        dir.join(format!("{project}.filtered.bam"))
    } else {
        dir.join(format!("{project}.filtered.{i}.bam"))
    }
}

fn ubam_header(cl: &str) -> Result<Header> {
    let text = format!(
        "@HD\tVN:1.6\tSO:unsorted\n@RG\tID:quartsx\n@PG\tID:quartsx\tPN:quartsx\tCL:{cl}\n"
    );
    text.parse::<Header>().context("building uBAM header")
}

fn build_record(name: &[u8], flags: u16, seq: &[u8], qual_ascii: &[u8]) -> RecordBuf {
    let mut rec = RecordBuf::default();
    *rec.name_mut() = Some(name.to_vec().into());
    *rec.flags_mut() = Flags::from(flags);
    *rec.sequence_mut() = Sequence::from(seq.to_vec());
    let phred: Vec<u8> = qual_ascii.iter().map(|&q| q.saturating_sub(crate::PHRED_OFFSET)).collect();
    *rec.quality_scores_mut() = QualityScores::from(phred);
    rec
}

fn add_bc_ub(rec: &mut RecordBuf, bc: &str, ub: &str) {
    let data = rec.data_mut();
    data.insert(Tag::from([b'B', b'C']), Value::String(bc.as_bytes().into()));
    data.insert(Tag::from([b'U', b'B']), Value::String(ub.as_bytes().into()));
}

// The BAM writer layered ON TOP of the BGZF codec (`BgzfW`: multithreaded for T>1, single-threaded
// for T<=1). Attached with `bam::io::Writer::from` (NOT `::new`, which would add a SECOND bgzf layer
// and double-compress); `BgzfW: std::io::Write` satisfies the `bam::io::Writer<W: Write>` bound.
type ShardW = bam::io::Writer<BgzfW>;

pub struct ShardSet {
    writers: Vec<ShardW>,
    header: Header,
}

impl ShardSet {
    /// Create the N filtered-shard writers. `compress_threads` (= P from the spec's thread knobs) is
    /// the TOTAL deflate budget; it is split evenly across the N shard writers (min 1 each) so N>1
    /// parallel STAR feeds do not oversubscribe. Worker count does NOT affect the output bytes
    /// (byte-identical checkpoint), only throughput.
    pub fn create(
        dir: &Path,
        n: usize,
        project: &str,
        compress_threads: usize,
        cl: &str,
    ) -> Result<ShardSet> {
        let header = ubam_header(cl)?;
        // compress_threads == 0 => single-threaded codec (T<=1); else split P evenly across the N
        // shard writers (min 1 each). Worker count never affects the output bytes (byte-identical).
        let per_shard = if compress_threads == 0 { 0 } else { (compress_threads / n.max(1)).max(1) };
        let mut writers = Vec::with_capacity(n);
        for i in 0..n {
            let path = filtered_shard_path(dir, project, i, n);
            let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
            let mut w = bam::io::Writer::from(BgzfW::create(file, per_shard));
            w.write_header(&header).context("writing shard header")?;
            writers.push(w);
        }
        Ok(ShardSet { writers, header })
    }

    pub fn write_pair(
        &mut self,
        shard: usize,
        name: &[u8],
        r1_seq: &[u8],
        r1_qual: &[u8],
        r2_seq: &[u8],
        r2_qual: &[u8],
        bc: &str,
        ub: &str,
    ) -> Result<()> {
        // Barcode+UMI ride as native SAM tags (BC:Z / UB:Z, UB empty on internal reads).
        // 77  = paired|unmapped|mate-unmapped|first  (0x1|0x4|0x8|0x40)
        // 141 = paired|unmapped|mate-unmapped|second (0x1|0x4|0x8|0x80)
        let mut r1 = build_record(name, 77, r1_seq, r1_qual);
        let mut r2 = build_record(name, 141, r2_seq, r2_qual);
        add_bc_ub(&mut r1, bc, ub);
        add_bc_ub(&mut r2, bc, ub);
        let w = &mut self.writers[shard];
        w.write_alignment_record(&self.header, &r1).context("writing shard read1")?;
        w.write_alignment_record(&self.header, &r2).context("writing shard read2")?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        for w in self.writers.drain(..) {
            // Pull the BGZF layer out of the bam writer and finish it: appends the BGZF EOF block
            // (joining the deflater/writer threads for the multithreaded codec) and flushes the file.
            w.into_inner().finish()?;
        }
        Ok(())
    }
}

pub struct BamReader {
    inner: bam::io::Reader<BgzfR>,
    pub header: Header,
    pub ref_names: Vec<String>,
}

impl BamReader {
    /// Open a BAM for reading with `decode_threads` BGZF inflate workers (`decode_threads == 0` =>
    /// single-threaded codec, T<=1). Both codecs serve decompressed bytes in FILE ORDER, so records
    /// come out deterministically regardless of worker count. Attached with `bam::io::Reader::from`
    /// (NOT `::new`, which would add a second `bgzf::Reader` layer); `read_header`/`read_record_buf`
    /// need only `R: std::io::Read`, which `BgzfR` satisfies.
    pub fn open(path: &Path, decode_threads: usize) -> Result<BamReader> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut inner = bam::io::Reader::from(BgzfR::open(file, decode_threads));
        let header = inner.read_header().context("reading BAM header")?;
        let ref_names = header
            .reference_sequences()
            .keys()
            .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
            .collect();
        Ok(BamReader { inner, header, ref_names })
    }

    /// Reads the next record into `rec`; returns false at EOF.
    pub fn next(&mut self, rec: &mut RecordBuf) -> Result<bool> {
        let n = self
            .inner
            .read_record_buf(&self.header, rec)
            .context("reading BAM record")?;
        Ok(n != 0)
    }
}

// The tagged-BAM writer uses the same `BgzfW` codec and `::from` attachment / `into_inner()` finish
// pattern as ShardSet.
type TaggedW = bam::io::Writer<BgzfW>;

pub struct TaggedWriter {
    inner: TaggedW,
    header: Header,
}

impl TaggedWriter {
    pub fn create(path: &Path, header: &Header, compress_threads: usize) -> Result<TaggedWriter> {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut inner = bam::io::Writer::from(BgzfW::create(file, compress_threads));
        inner.write_header(header).context("writing tagged BAM header")?;
        Ok(TaggedWriter { inner, header: header.clone() })
    }

    /// Adds GE/GI string tags to a record (empty = no exon/intron assignment).
    pub fn write(&mut self, rec: &mut RecordBuf, ge: &str, gi: &str) -> Result<()> {
        let data = rec.data_mut();
        data.insert(Tag::from([b'G', b'E']), Value::String(ge.as_bytes().into()));
        data.insert(Tag::from([b'G', b'I']), Value::String(gi.as_bytes().into()));
        self.inner.write_alignment_record(&self.header, rec).context("writing tagged record")?;
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        self.inner.into_inner().finish()
    }
}

// ---- record field accessors (kept here so count.rs never touches noodles directly) ----

pub fn flags(rec: &RecordBuf) -> u16 {
    u16::from(rec.flags())
}

pub fn name(rec: &RecordBuf) -> &[u8] {
    rec.name().map(|n| n.as_ref()).unwrap_or(b"")
}

pub fn ref_id(rec: &RecordBuf) -> Option<usize> {
    rec.reference_sequence_id()
}

pub fn start(rec: &RecordBuf) -> Option<i64> {
    rec.alignment_start().map(|p| usize::from(p) as i64)
}

pub fn seq_len(rec: &RecordBuf) -> usize {
    rec.sequence().len()
}

/// Reads a BAM string (Z) tag (e.g. BC/UB).
pub fn tag_string(rec: &RecordBuf, tag: [u8; 2]) -> Option<String> {
    match rec.data().get(&Tag::from(tag)) {
        Some(Value::String(s)) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
        _ => None,
    }
}

/// Reads a BAM integer tag (e.g. NH) regardless of its stored width.
pub fn tag_int(rec: &RecordBuf, tag: [u8; 2]) -> Option<i64> {
    match rec.data().get(&Tag::from(tag)) {
        Some(Value::Int8(v)) => Some(*v as i64),
        Some(Value::UInt8(v)) => Some(*v as i64),
        Some(Value::Int16(v)) => Some(*v as i64),
        Some(Value::UInt16(v)) => Some(*v as i64),
        Some(Value::Int32(v)) => Some(*v as i64),
        Some(Value::UInt32(v)) => Some(*v as i64),
        _ => None,
    }
}

/// Reference intervals (1-based inclusive) covered by aligned read bases (M/=/X); D and N advance
/// the reference cursor but are not counted as covered.
///
/// This is a pure GENOMIC-block accessor. Transcript-space projection of these blocks lives in
/// `gtf.rs` (`Gene::project_transcript_span`) and the length calc in `count.rs`.
pub fn covered_blocks(rec: &RecordBuf, start_pos: i64) -> Vec<(i64, i64)> {
    use noodles::sam::alignment::record::cigar::op::Kind;
    let mut blocks = Vec::new();
    let mut pos = start_pos;
    for op in rec.cigar().as_ref().iter() {
        let len = op.len() as i64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                blocks.push((pos, pos + len - 1));
                pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                pos += len;
            }
            _ => {}
        }
    }
    blocks
}
