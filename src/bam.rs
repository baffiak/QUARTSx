// All noodles BAM/SAM I/O lives here so the (fast-moving) noodles API surface is in one place.
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
use std::io::BufWriter;

const PHRED_OFFSET: u8 = 33;

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
    let phred: Vec<u8> = qual_ascii.iter().map(|&q| q.saturating_sub(PHRED_OFFSET)).collect();
    *rec.quality_scores_mut() = QualityScores::from(phred);
    rec
}

fn add_bc_ub(rec: &mut RecordBuf, bc: &str, ub: &str) {
    let data = rec.data_mut();
    data.insert(Tag::from([b'B', b'C']), Value::String(bc.as_bytes().into()));
    data.insert(Tag::from([b'U', b'B']), Value::String(ub.as_bytes().into()));
}

type ShardW = bam::io::Writer<bgzf::Writer<BufWriter<File>>>;

pub struct ShardSet {
    writers: Vec<ShardW>,
    header: Header,
}

impl ShardSet {
    pub fn create(dir: &std::path::Path, n: usize, cl: &str) -> Result<ShardSet> {
        let header = ubam_header(cl)?;
        let mut writers = Vec::with_capacity(n);
        for i in 0..n {
            let path = dir.join(format!("shard_{i}.bam"));
            let file = File::create(&path).with_context(|| format!("creating {}", path.display()))?;
            let mut w = bam::io::Writer::new(BufWriter::new(file));
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
        // Barcode+UMI ride as native SAM tags (BC:Z / UB:Z, UB empty on internal reads); STAR copies
        // input tags through to the aligned BAM, so both mates keep the identical read name.
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
        for mut w in self.writers.drain(..) {
            w.try_finish().context("finalizing shard bgzf stream")?;
        }
        Ok(())
    }
}

pub struct BamReader {
    inner: bam::io::Reader<bgzf::Reader<File>>,
    pub header: Header,
    pub ref_names: Vec<String>,
}

impl BamReader {
    pub fn open(path: &std::path::Path) -> Result<BamReader> {
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut inner = bam::io::Reader::new(file);
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

pub struct TaggedWriter {
    inner: bam::io::Writer<bgzf::Writer<BufWriter<File>>>,
    header: Header,
}

impl TaggedWriter {
    pub fn create(path: &std::path::Path, header: &Header) -> Result<TaggedWriter> {
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut inner = bam::io::Writer::new(BufWriter::new(file));
        inner.write_header(header).context("writing tagged BAM header")?;
        Ok(TaggedWriter { inner, header: header.clone() })
    }

    /// Adds GE/GI string tags to a record (empty = no exon/intron assignment). BC/UB were carried
    /// through the mapping by STAR, so the record already keeps them.
    pub fn write(&mut self, rec: &mut RecordBuf, ge: &str, gi: &str) -> Result<()> {
        let data = rec.data_mut();
        data.insert(Tag::from([b'G', b'E']), Value::String(ge.as_bytes().into()));
        data.insert(Tag::from([b'G', b'I']), Value::String(gi.as_bytes().into()));
        self.inner.write_alignment_record(&self.header, rec).context("writing tagged record")?;
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        self.inner.try_finish().context("finalizing tagged bgzf stream")?;
        Ok(())
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

/// Reads a BAM string (Z) tag, e.g. BC/UB carried through from the uBAM by STAR.
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
