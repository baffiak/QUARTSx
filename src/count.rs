use crate::bam::{self, BamReader, TaggedWriter};
use crate::config::{Config, MultiMapperMode};
use crate::gtf::{Annotation, Gene, Tree};
use anyhow::{bail, Context, Result};
use coitrees::IntervalTree;
use noodles::sam::alignment::RecordBuf;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::fs::File;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Duration;

pub const NONE: u32 = u32::MAX;

// gene-body 5'->3' coverage resolution
const NBINS: usize = 100;

// records processed between two heartbeat checks (internal gate only — NOT a printed record count;
// the old "N M records annotated" counter that exceeded the read count is removed per §3).
const BEAT_EVERY: u64 = 262_144;

const BATCH: usize = 32_768;

pub struct Row {
    pub bc: u32,
    pub umi: Vec<u8>, // empty on internal reads
    pub ge: u32,      // exon gene-set id or NONE
    pub gi: u32,      // intron gene-set id or NONE
    pub uniq: bool,   // uniquely mapped (NH==1)? drives unique-evidence-first molecule base (§4a)
}

// read-type taxonomy per cell (mapped fragments; unmapped is derived from filter totals later)
#[derive(Default, Clone)]
pub struct CellReads {
    pub exon: u64,
    pub intron: u64,
    pub intergenic: u64,
    pub ambiguity: u64,
    pub user: u64,
}

impl CellReads {
    pub fn mapped(&self) -> u64 {
        self.exon + self.intron + self.intergenic + self.ambiguity + self.user
    }
}

pub struct ReadTable {
    pub rows: Vec<Row>,
    pub gene_sets: Vec<Vec<u32>>,
    pub barcodes: Vec<String>,
    pub percell: Vec<CellReads>,
}

impl ReadTable {
    fn intern_bc(&mut self, map: &mut HashMap<String, u32>, bc: &str) -> u32 {
        if let Some(&i) = map.get(bc) {
            return i;
        }
        let i = self.barcodes.len() as u32;
        self.barcodes.push(bc.to_string());
        self.percell.push(CellReads::default());
        map.insert(bc.to_string(), i);
        i
    }
    pub fn set_genes(&self, set: u32) -> &[u32] {
        self.gene_sets[set as usize].as_slice()
    }
}

// A fragment being assembled from its mate(s) / multimapping loci. Holds ONLY currently-open
// fragments — for a unique (NH==1) concordant pair this is just mate1 until mate2 arrives, so the
// resident set is bounded by local coverage, NOT by the genome-wide read count. This is what makes
// the old 40 GB global name-keyed `assign` map unnecessary (§3): tagging is fused into the counting
// pass, so no name->gene map has to survive to a second pass.
struct OpenFrag {
    name: Vec<u8>,
    exon_ov: Vec<(u32, i64)>,
    intron_ov: Vec<(u32, i64)>,
    blocks: Vec<(usize, i64, i64)>, // (ref_id, start, end) covered ref intervals, all mates/loci
    mapped_len: i64,
    r1_rev: Option<bool>,
    r2_rev: Option<bool>,
    bc: String,
    ub: String,
    is_multi: bool,          // NH>1 at any record: gene set is the UNION across loci (§4)
    have_r1_primary: bool,   // completion tracking for the unique concordant fast path
    have_r2_primary: bool,
    recs: Vec<RecordBuf>,    // the actual records to write to the tagged BAM (both mates / all loci)
}

impl OpenFrag {
    fn new(a: Annotated) -> OpenFrag {
        let (r1_rev, r2_rev) = if a.is_r1 { (Some(a.rev), None) } else { (None, Some(a.rev)) };
        OpenFrag {
            name: a.name,
            exon_ov: a.exon_ov,
            intron_ov: a.intron_ov,
            blocks: a.frag_blocks,
            mapped_len: a.mapped_len,
            r1_rev,
            r2_rev,
            bc: a.bc,
            ub: a.ub,
            is_multi: a.nh > 1,
            have_r1_primary: a.is_r1 && a.is_primary,
            have_r2_primary: !a.is_r1 && a.is_primary,
            recs: vec![a.rec],
        }
    }

    // fold a second record of the SAME fragment (the mate, or another multimapping locus)
    fn merge(&mut self, a: Annotated) {
        for (g, l) in a.exon_ov {
            add(&mut self.exon_ov, g, l);
        }
        for (g, l) in a.intron_ov {
            add(&mut self.intron_ov, g, l);
        }
        self.blocks.extend(a.frag_blocks);
        self.mapped_len += a.mapped_len;
        if a.is_r1 {
            self.r1_rev = Some(a.rev);
            self.have_r1_primary |= a.is_primary;
        } else {
            self.r2_rev = Some(a.rev);
            self.have_r2_primary |= a.is_primary;
        }
        self.is_multi |= a.nh > 1;
        self.recs.push(a.rec);
    }

    // absorb another OpenFrag with the same read name (a cross-contig mate rejoined in the leftover
    // pass). Keep this (earlier-seen) frag's bc/ub — the first-seen-mate convention is unchanged.
    fn absorb(&mut self, other: OpenFrag) {
        for (g, l) in other.exon_ov {
            add(&mut self.exon_ov, g, l);
        }
        for (g, l) in other.intron_ov {
            add(&mut self.intron_ov, g, l);
        }
        self.blocks.extend(other.blocks);
        self.mapped_len += other.mapped_len;
        if other.r1_rev.is_some() {
            self.r1_rev = other.r1_rev;
        }
        if other.r2_rev.is_some() {
            self.r2_rev = other.r2_rev;
        }
        self.have_r1_primary |= other.have_r1_primary;
        self.have_r2_primary |= other.have_r2_primary;
        self.is_multi |= other.is_multi;
        self.recs.extend(other.recs);
    }

    fn complete(&self) -> bool {
        self.have_r1_primary && self.have_r2_primary
    }
}

// Cheap per-record data extracted serially in BAM (coordinate) order, then handed to the parallel
// annotator. Owns the RecordBuf so the fragment can carry it to the tagged-BAM write.
struct RawRec {
    rec: RecordBuf,
    name: Vec<u8>,
    is_r1: bool,
    rev: bool,
    is_primary: bool,
    is_supplementary: bool,
    nh: i64,
    bc: String,
    ub: String,
    rid: Option<usize>,
    blocks: Vec<(i64, i64)>, // covered ref intervals for THIS record (empty for supplementary)
    mapped_len: i64,
}

// Result of the pure, parallel interval-tree annotation of one `RawRec`.
struct Annotated {
    rec: RecordBuf,
    name: Vec<u8>,
    is_r1: bool,
    rev: bool,
    is_primary: bool,
    is_supplementary: bool,
    nh: i64,
    bc: String,
    ub: String,
    rid: Option<usize>,
    exon_ov: Vec<(u32, i64)>,
    intron_ov: Vec<(u32, i64)>,
    frag_blocks: Vec<(usize, i64, i64)>,
    mapped_len: i64,
}

// Pure interval-tree annotation of one `RawRec` (two immutable queries). Shared by the serial (T<=1)
// and rayon (T>1) paths so both produce identical `Annotated` values in batch order.
fn annotate(r: RawRec, ann: &Annotation, ref_names: &[String]) -> Annotated {
    let mut exon_ov = Vec::new();
    let mut intron_ov = Vec::new();
    if let Some(rid) = r.rid {
        if !r.is_supplementary {
            let chrom = &ref_names[rid];
            accumulate(&ann.exon_tree, chrom, &r.blocks, &mut exon_ov);
            accumulate(&ann.intron_tree, chrom, &r.blocks, &mut intron_ov);
        }
    }
    let frag_blocks: Vec<(usize, i64, i64)> = match r.rid {
        Some(rid) => r.blocks.iter().map(|&(s, e)| (rid, s, e)).collect(),
        None => Vec::new(),
    };
    Annotated {
        rec: r.rec,
        name: r.name,
        is_r1: r.is_r1,
        rev: r.rev,
        is_primary: r.is_primary,
        is_supplementary: r.is_supplementary,
        nh: r.nh,
        bc: r.bc,
        ub: r.ub,
        rid: r.rid,
        exon_ov,
        intron_ov,
        frag_blocks,
        mapped_len: r.mapped_len,
    }
}

pub fn count(cfg: &Config, ann: &Annotation, stage: &crate::log::Stage) -> Result<ReadTable> {
    let out = Path::new(&cfg.out_dir);
    let star_dir = out.join("star");
    std::fs::create_dir_all(&star_dir).context("creating star/")?;
    // Checkpoint 2 — STAR's own coordinate-sorted output name (lowercase `sorted`, cap `Coord`).
    let input = star_dir.join(format!("{}.Aligned.sortedByCoord.out.bam", cfg.project));
    // Checkpoint 3 is written unsorted (fragments finalize out of coordinate order because a mate's
    // gene can only be resolved once BOTH mates are seen), then samtools-sorted — "STAR's way" (§3).
    let unsorted = star_dir.join(format!("{}.tagged.unsorted.bam", cfg.project));
    let tagged = star_dir.join(format!("{}.tagged.bam", cfg.project));

    // P deflate/inflate workers for the BGZF reader+writer (§3 parallel I/O). At T<=1 take the
    // single-threaded path: io_threads == 0 selects the single-threaded BGZF codec and the annotation
    // runs serially (no rayon fan-out), so one working thread never oversubscribes the lone core. For
    // T>1 the global rayon pool (num_threads) drives the annotation and these are its codec pools.
    let (p, _n) = cfg.resolved_threads();
    let single = cfg.num_threads <= 1;
    let io_threads = if single { 0 } else { p.max(1) };

    stage.step("counting: fused tag + count over coordinate-sorted BAM");
    let mut reader = BamReader::open(&input, io_threads)?;
    let ref_names: Vec<String> = reader.ref_names.clone();

    // genome-position progress denominator (§3 progress bar). reference_sequences() is an IndexMap in
    // ref_id order, so values() line up with ref_names / bam::ref_id indices.
    let ref_lengths: Vec<u64> = reader
        .header
        .reference_sequences()
        .values()
        .map(|rs| usize::from(rs.length()) as u64)
        .collect();
    let mut genome_off: Vec<u64> = Vec::with_capacity(ref_lengths.len());
    let mut gacc = 0u64;
    for &l in &ref_lengths {
        genome_off.push(gacc);
        gacc += l;
    }
    let genome_len = gacc.max(1);
    let n_refs = ref_names.len();

    let mut state = CountState::new();
    state.ref_names = ref_names.clone();
    state.additional_chroms =
        ann.genes.iter().filter(|g| g.additional).map(|g| g.chrom.clone()).collect();

    let mut writer = TaggedWriter::create(&unsorted, &reader.header, io_threads)?;

    // Per-CONTIG pending (invariant 5): unique concordant pairs whose second mate has not yet arrived.
    // Freed at every contig boundary — a >100 kb intron keeps a pair open within a contig (fine), only
    // genuine cross-contig/singleton mates survive to `leftover`.
    let mut pending: HashMap<Vec<u8>, OpenFrag> = HashMap::new();
    // cross-contig mates + singletons, rejoined in a deterministic name-sorted leftover pass.
    let mut leftover: HashMap<Vec<u8>, OpenFrag> = HashMap::new();
    // multimapping fragments (NH>1) — their loci are scattered across the genome, so they cannot be
    // paired within one contig; accumulate the full UNION here and finalize in the leftover pass.
    // In the default Unique mode STAR emits no NH>1 records, so this stays empty.
    let mut multi: HashMap<Vec<u8>, OpenFrag> = HashMap::new();
    let mut current_rid: Option<usize> = None;

    let mut n_seen: u64 = 0;
    let mut raws: Vec<RawRec> = Vec::with_capacity(BATCH);
    let mut eof = false;
    loop {
        while raws.len() < BATCH {
            let mut rec = RecordBuf::default();
            if !reader.next(&mut rec)? {
                eof = true;
                break;
            }
            n_seen += 1;
            let fl = bam::flags(&rec);
            let is_supplementary = fl & 0x800 != 0;
            let is_primary = fl & 0x100 == 0 && !is_supplementary;
            let rid = bam::ref_id(&rec);
            if n_seen % BEAT_EVERY == 0 {
                if let (Some(r), Some(pos)) = (rid, bam::start(&rec)) {
                    let done = genome_off.get(r).copied().unwrap_or(genome_len) + pos as u64;
                    let pct = (done as f64 / genome_len as f64 * 100.0).min(100.0);
                    stage.beat(Duration::from_secs(10), || {
                        format!("counting: {pct:.1}% of genome (contig {}/{n_refs})", r + 1)
                    });
                }
            }
            let (blocks, mapped_len) = match (is_supplementary, rid, bam::start(&rec)) {
                (false, Some(_), Some(pos)) => {
                    let b = bam::covered_blocks(&rec, pos);
                    let ml: i64 = b.iter().map(|&(s, e)| e - s + 1).sum();
                    (b, ml)
                }
                _ => (Vec::new(), 0),
            };
            // BC/UB were passed through the mapping by STAR (empty UB => internal read); only the
            // first-seen mate's copy is used at insert time, exactly as before.
            let bc = bam::tag_string(&rec, [b'B', b'C']).unwrap_or_default();
            let ub = bam::tag_string(&rec, [b'U', b'B']).unwrap_or_default();
            let nh = bam::tag_int(&rec, [b'N', b'H']).unwrap_or(1);
            let name = bam::name(&rec).to_vec();
            raws.push(RawRec {
                rec,
                name,
                is_r1: fl & 0x40 != 0,
                rev: fl & 0x10 != 0,
                is_primary,
                is_supplementary,
                nh,
                bc,
                ub,
                rid,
                blocks,
                mapped_len,
            });
        }
        if raws.is_empty() {
            if eof {
                break;
            } else {
                continue;
            }
        }
        let batch = std::mem::replace(&mut raws, Vec::with_capacity(BATCH));
        // Heavy, embarrassingly-parallel step: two immutable interval-tree queries per record. Both
        // the serial (T<=1) and rayon (T>1) collects preserve batch order, so the serial drain below
        // sees records in exact coordinate (file) order regardless of thread count -> deterministic
        // rt.rows.
        let annotated: Vec<Annotated> = if single {
            batch.into_iter().map(|r| annotate(r, ann, &ref_names)).collect()
        } else {
            batch.into_par_iter().map(|r| annotate(r, ann, &ref_names)).collect()
        };

        for mut a in annotated {
            // contig boundary: free this contig's pending into the leftover pass (invariant 5).
            if let Some(rid) = a.rid {
                if current_rid != Some(rid) {
                    flush_pending(&mut pending, &mut leftover);
                    current_rid = Some(rid);
                }
            }
            // supplementary (chimeric) records: not part of the primary fragment span; write through
            // untagged in place (they're rare and off by default) rather than corrupt the mate-join.
            if a.is_supplementary {
                writer.write(&mut a.rec, "", "")?;
                continue;
            }
            // multimappers: route to the global union accumulator (finalized in the leftover pass).
            if a.nh > 1 {
                match multi.entry(a.name.clone()) {
                    std::collections::hash_map::Entry::Occupied(mut e) => e.get_mut().merge(a),
                    std::collections::hash_map::Entry::Vacant(e) => {
                        e.insert(OpenFrag::new(a));
                    }
                }
                continue;
            }
            // unique (NH==1) concordant fast path: pair the two mates within the contig.
            let name = a.name.clone();
            match pending.remove(&name) {
                None => {
                    pending.insert(name, OpenFrag::new(a));
                }
                Some(mut of) => {
                    of.merge(a);
                    if of.complete() {
                        finalize_and_write(&mut state, of, cfg, ann, &mut writer)?;
                    } else {
                        pending.insert(name, of); // pathological (e.g. two read1s) — keep open
                    }
                }
            }
        }
        if eof {
            break;
        }
    }

    // Finalize everything still open, in a deterministic (read-name) order so rt.rows push order — and
    // hence seeded downsampling — is reproducible regardless of thread count or HashMap seed.
    flush_pending(&mut pending, &mut leftover);
    let mut lo: Vec<OpenFrag> = leftover.into_values().collect();
    lo.sort_by(|a, b| a.name.cmp(&b.name));
    for of in lo {
        finalize_and_write(&mut state, of, cfg, ann, &mut writer)?;
    }
    let mut mp: Vec<OpenFrag> = multi.into_values().collect();
    mp.sort_by(|a, b| a.name.cmp(&b.name));
    for of in mp {
        finalize_and_write(&mut state, of, cfg, ann, &mut writer)?;
    }

    writer.finish()?;

    // "Merge = STAR's way": the unsorted tagged BAM -> samtools sort -@T -> coordinate-sorted, indexable
    // checkpoint 3; then drop the intermediate so exactly ONE ~90 GB checkpoint remains (§3 IMPL NOTE).
    stage.step("counting: samtools sort tagged BAM");
    samtools_sort(&unsorted, &tagged, cfg.num_threads, &out.join("logs"))?;

    state.write_qc(out, ann, &cfg.project)?;
    write_gene_names(out, ann, &cfg.project)?;
    Ok(state.rt)
}

// move a contig's still-open fragments into the leftover map, merging a same-named entry that was
// already parked there (the other half of a cross-contig pair seen on an earlier contig).
fn flush_pending(pending: &mut HashMap<Vec<u8>, OpenFrag>, leftover: &mut HashMap<Vec<u8>, OpenFrag>) {
    for (name, of) in pending.drain() {
        match leftover.entry(name) {
            std::collections::hash_map::Entry::Occupied(mut e) => e.get_mut().absorb(of),
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(of);
            }
        }
    }
}

fn finalize_and_write(
    state: &mut CountState,
    of: OpenFrag,
    cfg: &Config,
    ann: &Annotation,
    writer: &mut TaggedWriter,
) -> Result<()> {
    let (ge_set, gi_set) = state.finalize(&of, cfg, ann);
    let ge_str = set_to_str(ge_set, &state.rt, ann);
    let gi_str = set_to_str(gi_set, &state.rt, ann);
    let mut recs = of.recs;
    for rec in recs.iter_mut() {
        writer.write(rec, &ge_str, &gi_str)?;
    }
    Ok(())
}

fn samtools_sort(unsorted: &Path, sorted: &Path, threads: usize, logs: &Path) -> Result<()> {
    std::fs::create_dir_all(logs).ok();
    let logf = logs.join("samtools_sort.tagged.log");
    let f = File::create(&logf).with_context(|| format!("creating {}", logf.display()))?;
    let ferr = f.try_clone().with_context(|| format!("cloning {}", logf.display()))?;
    let status = Command::new("samtools")
        .args(["sort", "-@", &threads.to_string(), "-o"])
        .arg(sorted)
        .arg(unsorted)
        .stdin(Stdio::null())
        .stdout(Stdio::from(f))
        .stderr(Stdio::from(ferr))
        .status()
        .context("launching samtools sort")?;
    if !status.success() {
        let code = status.code().map(|c| c.to_string()).unwrap_or_else(|| "signal".into());
        bail!("samtools sort of tagged BAM failed (exit {code}); see {}", logf.display());
    }
    std::fs::remove_file(unsorted).ok();
    Ok(())
}

struct CountState {
    rt: ReadTable,
    bc_map: HashMap<String, u32>,
    set_map: HashMap<Vec<u32>, u32>,
    ref_names: Vec<String>,
    additional_chroms: HashSet<String>,      // contigs from additional_fasta (spike-ins)
    gb_tagged: HashMap<u32, [f64; NBINS]>,   // per-gene 5'->3' profile, UMI-tagged reads
    gb_internal: HashMap<u32, [f64; NBINS]>, // per-gene 5'->3' profile, internal reads
    ins_tagged: HashMap<(bool, i64), u64>,   // insert-size histogram, UMI-tagged, keyed (additional, size)
    ins_internal: HashMap<(bool, i64), u64>, // insert-size histogram, internal, keyed (additional, size)
}

impl CountState {
    fn new() -> CountState {
        CountState {
            rt: ReadTable { rows: Vec::new(), gene_sets: Vec::new(), barcodes: Vec::new(), percell: Vec::new() },
            bc_map: HashMap::new(),
            set_map: HashMap::new(),
            ref_names: Vec::new(),
            additional_chroms: HashSet::new(),
            gb_tagged: HashMap::new(),
            gb_internal: HashMap::new(),
            ins_tagged: HashMap::new(),
            ins_internal: HashMap::new(),
        }
    }

    // Resolve a fragment's gene set(s), fold its QC, push its per-fragment row, and RETURN the interned
    // (exon, intron) gene-set ids so the caller can tag the fragment's records. No global name-keyed
    // map is produced (§3): the assignment is consumed immediately.
    fn finalize(&mut self, acc: &OpenFrag, cfg: &Config, ann: &Annotation) -> (u32, u32) {
        let tagged = !acc.ub.is_empty();
        let strand = if let Some(r) = acc.r1_rev {
            if r { '-' } else { '+' }
        } else if let Some(r) = acc.r2_rev {
            if r { '+' } else { '-' } // read2 is antisense to the fragment
        } else {
            '+'
        };

        let (ge_genes, ge_amb) = resolve(&acc.exon_ov, tagged, strand, acc.mapped_len, acc.is_multi, cfg, ann);
        let (gi_genes, gi_amb) = resolve(&acc.intron_ov, tagged, strand, acc.mapped_len, acc.is_multi, cfg, ann);

        // gene-body 5'->3' coverage: fold each exon-assigned gene's own transcript profile (gene strand)
        {
            let ref_names = &self.ref_names;
            let map = if tagged { &mut self.gb_tagged } else { &mut self.gb_internal };
            for &g in &ge_genes {
                let prof = map.entry(g).or_insert([0.0; NBINS]);
                add_coverage(prof, &ann.genes[g as usize], &acc.blocks, ref_names);
            }
        }

        // insert size (§5): TRANSCRIPT-SPACE span, not genomic. Requires a single resolved gene
        // (ge_genes.len()==1) and both mapped mates; each mate's aligned blocks are projected onto the
        // gene's exon-union transcript coordinate (introns collapse, a distant mate projects onto ~no
        // exon so it cannot inflate), and the span is tmax-tmin+1. Category from gene.additional.
        if acc.r1_rev.is_some() && acc.r2_rev.is_some() && ge_genes.len() == 1 {
            let g = &ann.genes[ge_genes[0] as usize];
            if let Some((tmin, tmax)) = g.project_transcript_span(&acc.blocks, &self.ref_names) {
                let isize = tmax - tmin + 1;
                if isize > 0 {
                    let hist = if tagged { &mut self.ins_tagged } else { &mut self.ins_internal };
                    *hist.entry((g.additional, isize)).or_insert(0) += 1;
                }
            }
        }

        let bc_id = self.rt.intern_bc(&mut self.bc_map, &acc.bc);
        {
            let cell = &mut self.rt.percell[bc_id as usize];
            if !ge_genes.is_empty() {
                if ann.genes[ge_genes[0] as usize].additional {
                    cell.user += 1;
                } else {
                    cell.exon += 1;
                }
            } else if ge_amb {
                cell.ambiguity += 1;
            } else if !gi_genes.is_empty() {
                cell.intron += 1;
            } else if gi_amb {
                cell.ambiguity += 1;
            } else {
                cell.intergenic += 1;
            }
        }

        let ge_set = intern_set(&mut self.rt, &mut self.set_map, ge_genes);
        let gi_set = intern_set(&mut self.rt, &mut self.set_map, gi_genes);
        if ge_set != NONE || gi_set != NONE {
            self.rt.rows.push(Row { bc: bc_id, umi: acc.ub.clone().into_bytes(), ge: ge_set, gi: gi_set, uniq: !acc.is_multi });
        }
        (ge_set, gi_set)
    }

    fn write_qc(&self, out: &Path, ann: &Annotation, project: &str) -> Result<()> {
        write_genebody(out, ann, &self.gb_tagged, &self.gb_internal, project)?;
        write_insertsize(out, &self.ins_tagged, &self.ins_internal, project)?;
        Ok(())
    }
}

// Fold a fragment's exonic coverage into a gene's 100-bin 5'->3' profile. Reference bases are mapped
// into the gene's exon-union transcript coordinate; the 5' end is the gene start on +, gene end on -.
fn add_coverage(prof: &mut [f64; NBINS], gene: &Gene, blocks: &[(usize, i64, i64)], ref_names: &[String]) {
    let len = gene.length;
    if len <= 0 {
        return;
    }
    for &(rid, bs, be) in blocks {
        if ref_names[rid] != gene.chrom {
            continue;
        }
        let mut cum = 0i64; // transcript bases before the current exon
        for &(es, ee) in &gene.exons {
            let a = bs.max(es);
            let b = be.min(ee);
            if a <= b {
                for pos in a..=b {
                    let fwd = cum + (pos - es);
                    let off = if gene.strand == '-' { len - 1 - fwd } else { fwd };
                    let bin = ((off * NBINS as i64 / len) as usize).min(NBINS - 1);
                    prof[bin] += 1.0;
                }
            }
            cum += ee - es + 1;
        }
    }
}

// Mean of per-gene max-normalized profiles over the genes of one class (equal weight).
fn aggregate(map: &HashMap<u32, [f64; NBINS]>, ann: &Annotation, additional: bool) -> Option<[f64; NBINS]> {
    let mut agg = [0.0f64; NBINS];
    let mut n = 0u32;
    // Sort gene ids so the per-gene FP reduction runs in a fixed order (HashMap order is unstable).
    let mut genes: Vec<u32> = map.keys().copied().filter(|&g| ann.genes[g as usize].additional == additional).collect();
    genes.sort_unstable();
    for g in genes {
        let prof = &map[&g];
        let max = prof.iter().copied().fold(0.0f64, f64::max);
        if max <= 0.0 {
            continue;
        }
        for i in 0..NBINS {
            agg[i] += prof[i] / max;
        }
        n += 1;
    }
    if n == 0 {
        return None;
    }
    for v in agg.iter_mut() {
        *v /= n as f64;
    }
    Some(agg)
}

fn write_genebody(
    out: &Path,
    ann: &Annotation,
    tagged: &HashMap<u32, [f64; NBINS]>,
    internal: &HashMap<u32, [f64; NBINS]>,
    project: &str,
) -> Result<()> {
    let qc = out.join("qc");
    std::fs::create_dir_all(&qc).context("creating qc/")?;
    let mut s = String::from("category\tread_type\tpercentile\tcoverage\n");
    for (read_type, map) in [("UMI-tagged", tagged), ("Internal", internal)] {
        for (category, additional) in [("Transcriptome", false), ("Spike-in", true)] {
            if let Some(prof) = aggregate(map, ann, additional) {
                for (i, &cov) in prof.iter().enumerate() {
                    let _ = write!(s, "{category}\t{read_type}\t{}\t{cov}\n", i + 1);
                }
            }
        }
    }
    std::fs::write(qc.join(format!("{project}.genebody_coverage.txt")), s).context("writing genebody_coverage.txt")?;
    Ok(())
}

fn write_insertsize(
    out: &Path,
    tagged: &HashMap<(bool, i64), u64>,
    internal: &HashMap<(bool, i64), u64>,
    project: &str,
) -> Result<()> {
    let qc = out.join("qc");
    std::fs::create_dir_all(&qc).context("creating qc/")?;
    let mut s = String::from("category\tread_type\tinsert_size\tcount\n");
    for (read_type, hist) in [("UMI-tagged", tagged), ("Internal", internal)] {
        for (category, additional) in [("Transcriptome", false), ("Spike-in", true)] {
            let mut sizes: Vec<i64> = hist.keys().filter(|(a, _)| *a == additional).map(|(_, sz)| *sz).collect();
            sizes.sort_unstable();
            for sz in sizes {
                let _ = write!(s, "{category}\t{read_type}\t{sz}\t{}\n", hist[&(additional, sz)]);
            }
        }
    }
    std::fs::write(qc.join(format!("{project}.insertsize.txt")), s).context("writing insertsize.txt")?;
    Ok(())
}

// Returns (winning genes, was_ambiguous).
// - `multimapper` (NH>1): the gene set is the UNION across the read's loci (§4) — every distinct gene
//   with a real overlap (strand-filtered), never "ambiguous". Distribution across the set is handled
//   downstream by `write_multimapper_matrices` (normalize over |ug|, never over NH).
// - unique: best-overlap gene(s) passing the fraction gate; a >1-gene tie is ambiguous unless
//   multi_overlap is set (then the tie becomes a multi-gene set).
fn resolve(
    ov: &[(u32, i64)],
    tagged: bool,
    strand: char,
    mapped_len: i64,
    multimapper: bool,
    cfg: &Config,
    ann: &Annotation,
) -> (Vec<u32>, bool) {
    if ov.is_empty() {
        return (Vec::new(), false);
    }
    let stranded = cfg.counting_opts.strand == 1 && tagged;

    if multimapper {
        let mut genes: Vec<u32> = ov
            .iter()
            .filter(|&&(g, l)| l >= 1 && !(stranded && ann.genes[g as usize].strand != strand))
            .map(|&(g, _)| g)
            .collect();
        genes.sort_unstable();
        genes.dedup();
        return (genes, false);
    }

    let mut best = 0i64;
    for &(g, l) in ov {
        if stranded && ann.genes[g as usize].strand != strand {
            continue;
        }
        if l < 1 || (mapped_len > 0 && (l as f64) / (mapped_len as f64) < cfg.counting_opts.fraction_overlap) {
            continue;
        }
        if l > best {
            best = l;
        }
    }
    if best < 1 {
        return (Vec::new(), false);
    }
    let mut winners: Vec<u32> = Vec::new();
    for &(g, l) in ov {
        if stranded && ann.genes[g as usize].strand != strand {
            continue;
        }
        if l == best {
            winners.push(g);
        }
    }
    winners.sort_unstable();
    if cfg.counting_opts.multi_overlap || winners.len() == 1 {
        (winners, false)
    } else {
        (Vec::new(), true) // ambiguous multi-gene tie
    }
}

fn intern_set(rt: &mut ReadTable, set_map: &mut HashMap<Vec<u32>, u32>, genes: Vec<u32>) -> u32 {
    if genes.is_empty() {
        return NONE;
    }
    if let Some(&i) = set_map.get(&genes) {
        return i;
    }
    let i = rt.gene_sets.len() as u32;
    rt.gene_sets.push(genes.clone());
    set_map.insert(genes, i);
    i
}

fn set_to_str(set: u32, rt: &ReadTable, ann: &Annotation) -> String {
    if set == NONE {
        return String::new();
    }
    rt.set_genes(set)
        .iter()
        .map(|&g| ann.genes[g as usize].id.as_str())
        .collect::<Vec<_>>()
        .join(",")
}

fn accumulate(trees: &HashMap<String, Tree>, chrom: &str, blocks: &[(i64, i64)], ov: &mut Vec<(u32, i64)>) {
    let tree = match trees.get(chrom) {
        Some(t) => t,
        None => return,
    };
    for &(s, e) in blocks {
        tree.query(s as i32, e as i32, |node| {
            let os = s.max(node.first as i64);
            let oe = e.min(node.last as i64);
            let olen = oe - os + 1;
            if olen > 0 {
                add(ov, node.metadata, olen);
            }
        });
    }
}

fn add(v: &mut Vec<(u32, i64)>, g: u32, len: i64) {
    for e in v.iter_mut() {
        if e.0 == g {
            e.1 += len;
            return;
        }
    }
    v.push((g, len));
}

// gene_names.txt CARRIES a header row now (§6): `gene_id\tgene_name`.
fn write_gene_names(out: &Path, ann: &Annotation, project: &str) -> Result<()> {
    let dir = out.join("expression");
    std::fs::create_dir_all(&dir).context("creating expression/")?;
    let mut s = String::from("gene_id\tgene_name\n");
    for g in &ann.genes {
        let name = g.name.as_deref().unwrap_or(&g.id);
        let _ = write!(s, "{}\t{}\n", g.id, name);
    }
    std::fs::write(dir.join(format!("{project}.gene_names.txt")), s).context("writing gene_names.txt")?;
    Ok(())
}

fn mode_label(m: MultiMapperMode) -> &'static str {
    match m {
        MultiMapperMode::Unique => "unique",
        MultiMapperMode::Uniform => "uniform",
        MultiMapperMode::Rescue => "rescue",
        MultiMapperMode::PropUnique => "propunique",
        MultiMapperMode::EM => "em",
    }
}

/// Distribute multi-gene units across their genes per the selected `--soloMultiMappers` mode and
/// write the resulting MULTI-ONLY mass into the counts.tsv.gz handoff (quant = `<family>_mult_<mode>`)
/// for the given feature/level.
///
/// Called from dedup::collapse per LAYER and LEVEL, once per FAMILY (`umicount`=tagged molecules,
/// `readcount_internal`=internal fragments), so the resolver mirrors the integer matrix on every
/// axis (§4). `counts[(bc, gene set)]` is the collapsed count for a (cell, gene-set) — molecules on
/// the INTERSECTED set for tagged, internal fragments on the per-fragment set for internal. The
/// always-emitted integer unique-gene matrix `gEu` is dedup's standard integer block; `gEu +
/// <family>_mult_<mode>` reconstructs the `UniqueAndMult-<mode>` matrix.
///
/// Per-cell terms: `gEu[g]` = count of units whose gene set is exactly `{g}`; `gEuniform[g]` = Σ over
/// multi-gene sets S∋g of `count(S)/|S|`. Modes:
///   Uniform    — M[g] = gEuniform[g]
///   PropUnique — per set S, norm=Σ_{h∈S} gEu[h]; norm=0 -> uniform fallback, else g∈S gets
///                count(S)·gEu[g]/norm
///   Rescue     — per set S, w[g]=gEuniform[g]+gEu[g]; Z=Σ_{h∈S} w[h]; g gets count(S)·w[g]/Z
///   EM         — θ init = gEuniform+gEu; each iteration θ_new=gEu, then per set S add
///                count(S)·θ_old[g]/Σ_{h∈S}θ_old[h]; prune θ_new<0.01 -> 0; stop at maxAbsChange<0.01
///                or after 100 iters
/// The published mass is the distributed MULTI portion only (EM = θ_final − gEu). Normalization is
/// ALWAYS over the number of distinct written genes |ug|, NEVER over NH (§4). Determinism: cells in
/// bc-id order, per-set genes ascending, fixed EM iters/tolerance -> byte-identical run to run.
pub fn write_multimapper_matrices(
    w: &mut impl std::io::Write,
    rt: &ReadTable,
    ann: &Annotation,
    mode: MultiMapperMode,
    counts: &HashMap<(u32, Vec<u32>), f64>,
    family: &str,
    feature: &str,
    level: &str,
) -> Result<()> {
    if mode == MultiMapperMode::Unique {
        return Ok(()); // Unique publishes only the integer gEu matrix; no multi mass to distribute.
    }
    // group the (cell, gene-set) counts by cell in a FIXED (bc, gene set) order for reproducible FP sums.
    let mut keys: Vec<(&(u32, Vec<u32>), &f64)> = counts.iter().collect();
    keys.sort_unstable_by(|a, b| a.0.cmp(b.0));
    let mut by_cell: BTreeMap<u32, Vec<(&[u32], f64)>> = BTreeMap::new();
    for (k, &c) in keys {
        by_cell.entry(k.0).or_default().push((k.1.as_slice(), c));
    }

    let mut out: BTreeMap<(u32, u32), f64> = BTreeMap::new(); // (gene, bc) -> multi-only mass
    for (&bc, sets) in &by_cell {
        let mut geu: BTreeMap<u32, f64> = BTreeMap::new();
        let mut guni: BTreeMap<u32, f64> = BTreeMap::new();
        let mut multis: Vec<(&[u32], f64)> = Vec::new();
        for &(genes, c) in sets {
            if genes.len() == 1 {
                *geu.entry(genes[0]).or_default() += c;
            } else if genes.len() > 1 {
                let inv = c / genes.len() as f64;
                for &g in genes {
                    *guni.entry(g).or_default() += inv;
                }
                multis.push((genes, c));
            }
        }

        let mut mass: BTreeMap<u32, f64> = BTreeMap::new();
        match mode {
            MultiMapperMode::Unique => {}
            MultiMapperMode::Uniform => {
                // M[g] = gEuniform[g]
                for (&g, &v) in &guni {
                    *mass.entry(g).or_default() += v;
                }
            }
            MultiMapperMode::PropUnique => {
                for (genes, c) in &multis {
                    let norm: f64 = genes.iter().map(|&g| *geu.get(&g).unwrap_or(&0.0)).sum();
                    if norm > 0.0 {
                        for &g in *genes {
                            let u = *geu.get(&g).unwrap_or(&0.0);
                            if u > 0.0 {
                                *mass.entry(g).or_default() += c * u / norm;
                            }
                        }
                    } else {
                        let inv = c / genes.len() as f64;
                        for &g in *genes {
                            *mass.entry(g).or_default() += inv;
                        }
                    }
                }
            }
            MultiMapperMode::Rescue => {
                for (genes, c) in &multis {
                    let z: f64 = genes
                        .iter()
                        .map(|&g| geu.get(&g).unwrap_or(&0.0) + guni.get(&g).unwrap_or(&0.0))
                        .sum();
                    if z > 0.0 {
                        for &g in *genes {
                            let wgt = geu.get(&g).unwrap_or(&0.0) + guni.get(&g).unwrap_or(&0.0);
                            if wgt > 0.0 {
                                *mass.entry(g).or_default() += c * wgt / z;
                            }
                        }
                    } else {
                        let inv = c / genes.len() as f64;
                        for &g in *genes {
                            *mass.entry(g).or_default() += inv;
                        }
                    }
                }
            }
            MultiMapperMode::EM => {
                // θ init = gEuniform + gEu
                let mut theta: BTreeMap<u32, f64> = BTreeMap::new();
                for (&g, &v) in &geu {
                    *theta.entry(g).or_default() += v;
                }
                for (&g, &v) in &guni {
                    *theta.entry(g).or_default() += v;
                }
                for _iter in 0..100 {
                    let mut theta_new: BTreeMap<u32, f64> = geu.clone(); // start from the unique anchor
                    for (genes, c) in &multis {
                        let z: f64 = genes.iter().map(|&g| *theta.get(&g).unwrap_or(&0.0)).sum();
                        if z > 0.0 {
                            for &g in *genes {
                                let t = *theta.get(&g).unwrap_or(&0.0);
                                if t > 0.0 {
                                    *theta_new.entry(g).or_default() += c * t / z;
                                }
                            }
                        } else {
                            let inv = c / genes.len() as f64;
                            for &g in *genes {
                                *theta_new.entry(g).or_default() += inv;
                            }
                        }
                    }
                    for v in theta_new.values_mut() {
                        if *v < 0.01 {
                            *v = 0.0;
                        }
                    }
                    let mut both: BTreeSet<u32> = theta.keys().copied().collect();
                    both.extend(theta_new.keys().copied());
                    let mut maxchg = 0.0f64;
                    for g in both {
                        let a = *theta.get(&g).unwrap_or(&0.0);
                        let b = *theta_new.get(&g).unwrap_or(&0.0);
                        maxchg = maxchg.max((a - b).abs());
                    }
                    theta = theta_new;
                    if maxchg < 0.01 {
                        break;
                    }
                }
                // published mass = θ_final − gEu (the distributed multi portion)
                for (&g, &t) in &theta {
                    let m = t - geu.get(&g).unwrap_or(&0.0);
                    if m > 0.0 {
                        *mass.entry(g).or_default() += m;
                    }
                }
            }
        }

        for (g, m) in mass {
            if m != 0.0 {
                *out.entry((g, bc)).or_default() += m;
            }
        }
    }

    // Emit in canonical (gene_id string, cell string) order so the block is byte-identical run to run
    // and thread-count independent, matching dedup's `emit`.
    let quant = format!("{family}_mult_{}", mode_label(mode));
    let mut ks: Vec<&(u32, u32)> = out.keys().collect();
    ks.sort_unstable_by(|a, b| {
        ann.genes[a.0 as usize]
            .id
            .cmp(&ann.genes[b.0 as usize].id)
            .then_with(|| rt.barcodes[a.1 as usize].cmp(&rt.barcodes[b.1 as usize]))
    });
    for k in ks {
        let v = out[k];
        writeln!(
            w,
            "{quant}\t{feature}\t{level}\t{}\t{}\t{v}",
            ann.genes[k.0 as usize].id, rt.barcodes[k.1 as usize]
        )?;
    }
    Ok(())
}
