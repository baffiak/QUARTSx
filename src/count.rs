use crate::bam::{self, BamReader, TaggedWriter};
use crate::config::Config;
use crate::gtf::{Annotation, Gene, Tree};
use anyhow::{Context, Result};
use coitrees::IntervalTree;
use noodles::sam::alignment::RecordBuf;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;
use std::path::Path;

pub const NONE: u32 = u32::MAX;

// gene-body 5'->3' coverage resolution (RSeQC-style)
const NBINS: usize = 100;

pub struct Row {
    pub bc: u32,
    pub umi: Vec<u8>, // empty on internal reads
    pub ge: u32,      // exon gene-set id or NONE
    pub gi: u32,      // intron gene-set id or NONE
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

struct MateAcc {
    name: Vec<u8>,
    exon_ov: Vec<(u32, i64)>,
    intron_ov: Vec<(u32, i64)>,
    blocks: Vec<(usize, i64, i64)>, // (ref_id, start, end) covered ref intervals, both mates
    mapped_len: i64,
    r1_rev: Option<bool>,
    r2_rev: Option<bool>,
    bc: String,
    ub: String,
}

// Cheap per-record data extracted serially in BAM order (survives the primary/multimapper filters),
// then handed to the parallel annotator. `blocks` are covered reference intervals for THIS mate.
struct RawRec {
    name: Vec<u8>,
    is_r1: bool,
    rev: bool,
    bc: String,
    ub: String,
    rid: usize,
    blocks: Vec<(i64, i64)>,
    mapped_len: i64,
}

// Result of the pure, parallel interval-tree annotation of one `RawRec`.
struct Annotated {
    name: Vec<u8>,
    is_r1: bool,
    rev: bool,
    bc: String,
    ub: String,
    exon_ov: Vec<(u32, i64)>,
    intron_ov: Vec<(u32, i64)>,
    frag_blocks: Vec<(usize, i64, i64)>,
    mapped_len: i64,
}

const BATCH: usize = 32_768;

pub fn count(cfg: &Config, ann: &Annotation, stage: &crate::log::Stage) -> Result<ReadTable> {
    use std::time::Duration;
    let out = Path::new(&cfg.out_dir);
    let merged = out.join("star").join(format!("{}.merged.bam", cfg.project));

    let mut state = CountState::new();

    // ---- pass 1: pair mates, assign genes ----
    // The interval-tree annotation (covered_blocks feeding two immutable-tree queries) is the heavy
    // per-record cost and is embarrassingly parallel. We read a batch serially in BAM order, run the
    // pure annotation in parallel (into_par_iter preserves batch order), then drain the ordered batch
    // serially so the mate-join, interning and rt.rows push order are byte-identical to single-thread.
    stage.step("pass 1: annotate + pair mates");
    let mut reader = BamReader::open(&merged)?;
    let ref_names: Vec<String> = reader.ref_names.clone();
    state.ref_names = ref_names.clone();
    state.additional_chroms = ann.genes.iter().filter(|g| g.additional).map(|g| g.chrom.clone()).collect();
    let mut pending: HashMap<Vec<u8>, MateAcc> = HashMap::new();
    let mut rec = RecordBuf::default();
    let mut n_rec = 0u64;
    let mut raws: Vec<RawRec> = Vec::with_capacity(BATCH);
    let mut eof = false;
    loop {
        while raws.len() < BATCH {
            if !reader.next(&mut rec)? {
                eof = true;
                break;
            }
            n_rec += 1;
            if n_rec % 1_000_000 == 0 {
                stage.beat(Duration::from_secs(5), || format!("{:.1}M records annotated", n_rec as f64 / 1e6));
            }
            let fl = bam::flags(&rec);
            if fl & 0x100 != 0 || fl & 0x800 != 0 {
                continue; // primary alignments only for assignment
            }
            // primary_hit=false == zUMIs countMultiMappingReads=FALSE: drop multimappers entirely
            if !cfg.counting_opts.primary_hit && bam::tag_int(&rec, [b'N', b'H']).unwrap_or(1) > 1 {
                continue;
            }
            let rid = match bam::ref_id(&rec) {
                Some(r) => r,
                None => continue,
            };
            let pos = match bam::start(&rec) {
                Some(p) => p,
                None => continue,
            };
            let blocks = bam::covered_blocks(&rec, pos);
            let mapped_len: i64 = blocks.iter().map(|&(s, e)| e - s + 1).sum();
            // BC/UB were passed through the mapping by STAR (empty UB => internal read); only the
            // first-seen mate's copy is used at insert time, exactly as before.
            let bc = bam::tag_string(&rec, [b'B', b'C']).unwrap_or_default();
            let ub = bam::tag_string(&rec, [b'U', b'B']).unwrap_or_default();
            raws.push(RawRec {
                name: bam::name(&rec).to_vec(),
                is_r1: fl & 0x40 != 0,
                rev: fl & 0x10 != 0,
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
        let annotated: Vec<Annotated> = batch
            .into_par_iter()
            .map(|r| {
                let chrom = &ref_names[r.rid];
                let mut exon_ov = Vec::new();
                let mut intron_ov = Vec::new();
                accumulate(&ann.exon_tree, chrom, &r.blocks, &mut exon_ov);
                accumulate(&ann.intron_tree, chrom, &r.blocks, &mut intron_ov);
                let frag_blocks: Vec<(usize, i64, i64)> = r.blocks.iter().map(|&(s, e)| (r.rid, s, e)).collect();
                Annotated {
                    name: r.name,
                    is_r1: r.is_r1,
                    rev: r.rev,
                    bc: r.bc,
                    ub: r.ub,
                    exon_ov,
                    intron_ov,
                    frag_blocks,
                    mapped_len: r.mapped_len,
                }
            })
            .collect();
        for a in annotated {
            match pending.remove(&a.name) {
                None => {
                    pending.insert(
                        a.name.clone(),
                        MateAcc {
                            name: a.name,
                            bc: a.bc,
                            ub: a.ub,
                            exon_ov: a.exon_ov,
                            intron_ov: a.intron_ov,
                            blocks: a.frag_blocks,
                            mapped_len: a.mapped_len,
                            r1_rev: if a.is_r1 { Some(a.rev) } else { None },
                            r2_rev: if a.is_r1 { None } else { Some(a.rev) },
                        },
                    );
                }
                Some(mut acc) => {
                    for (gene, l) in a.exon_ov {
                        add(&mut acc.exon_ov, gene, l);
                    }
                    for (gene, l) in a.intron_ov {
                        add(&mut acc.intron_ov, gene, l);
                    }
                    acc.blocks.extend(a.frag_blocks);
                    acc.mapped_len += a.mapped_len;
                    if a.is_r1 {
                        acc.r1_rev = Some(a.rev);
                    } else {
                        acc.r2_rev = Some(a.rev);
                    }
                    state.finalize(&acc, cfg, ann);
                }
            }
        }
        if eof {
            break;
        }
    }
    // leftover singletons (one mate mapped) — drained in a deterministic (read-name) order so that
    // rt.rows push order (and hence downsampling) is reproducible regardless of HashMap seed.
    let mut leftovers: Vec<MateAcc> = pending.into_values().collect();
    leftovers.sort_by(|a, b| a.name.cmp(&b.name));
    for acc in &leftovers {
        state.finalize(acc, cfg, ann);
    }

    state.write_qc(out, ann, &cfg.project)?;

    // ---- pass 2: tag every record with GE/GI (BC/UB rode through from the uBAM) ----
    // GE/GI string building is parallel per batch; records are written serially in BAM order so the
    // final BAM byte layout is identical to single-thread.
    stage.step("pass 2: writing tagged BAM");
    let mut reader = BamReader::open(&merged)?;
    let final_bam = out.join(format!("{}.bam", cfg.project));
    let mut writer = TaggedWriter::create(&final_bam, &reader.header)?;
    let mut batch: Vec<RecordBuf> = Vec::with_capacity(BATCH);
    loop {
        batch.clear();
        let mut eof = false;
        while batch.len() < BATCH {
            let mut r = RecordBuf::default();
            if !reader.next(&mut r)? {
                eof = true;
                break;
            }
            batch.push(r);
        }
        if batch.is_empty() {
            break;
        }
        let tags: Vec<(String, String)> = batch
            .par_iter()
            .map(|rec| {
                let name = bam::name(rec).to_vec();
                let (ge, gi) = state.assign.get(&name).copied().unwrap_or((NONE, NONE));
                (set_to_str(ge, &state.rt, ann), set_to_str(gi, &state.rt, ann))
            })
            .collect();
        for (rec, (ge_str, gi_str)) in batch.iter_mut().zip(tags) {
            writer.write(rec, &ge_str, &gi_str)?;
        }
        if eof {
            break;
        }
    }
    writer.finish()?;

    write_gene_names(out, ann, &cfg.project)?;
    Ok(state.rt)
}

struct CountState {
    rt: ReadTable,
    bc_map: HashMap<String, u32>,
    set_map: HashMap<Vec<u32>, u32>,
    assign: HashMap<Vec<u8>, (u32, u32)>, // read name -> (ge_set, gi_set)
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
            assign: HashMap::new(),
            ref_names: Vec::new(),
            additional_chroms: HashSet::new(),
            gb_tagged: HashMap::new(),
            gb_internal: HashMap::new(),
            ins_tagged: HashMap::new(),
            ins_internal: HashMap::new(),
        }
    }

    fn finalize(&mut self, acc: &MateAcc, cfg: &Config, ann: &Annotation) {
        let tagged = !acc.ub.is_empty();
        let strand = if let Some(r) = acc.r1_rev {
            if r { '-' } else { '+' }
        } else if let Some(r) = acc.r2_rev {
            if r { '+' } else { '-' } // read2 is antisense to the fragment
        } else {
            '+'
        };

        let (ge_genes, ge_amb) = resolve(&acc.exon_ov, tagged, strand, acc.mapped_len, cfg, ann);
        let (gi_genes, gi_amb) = resolve(&acc.intron_ov, tagged, strand, acc.mapped_len, cfg, ann);

        // gene-body 5'->3' coverage: fold each exon-assigned gene's own transcript profile (gene strand)
        {
            let ref_names = &self.ref_names;
            let map = if tagged { &mut self.gb_tagged } else { &mut self.gb_internal };
            for &g in &ge_genes {
                let prof = map.entry(g).or_insert([0.0; NBINS]);
                add_coverage(prof, &ann.genes[g as usize], &acc.blocks, ref_names);
            }
        }

        // insert size: outer reference span of both mapped mates on the same contig,
        // split by category (spike-in contig or main transcriptome) like the coverage panel
        if acc.r1_rev.is_some() && acc.r2_rev.is_some() && !acc.blocks.is_empty() {
            let rid0 = acc.blocks[0].0;
            if acc.blocks.iter().all(|b| b.0 == rid0) {
                let min_s = acc.blocks.iter().map(|b| b.1).min().unwrap();
                let max_e = acc.blocks.iter().map(|b| b.2).max().unwrap();
                let isize = max_e - min_s + 1;
                if isize > 0 {
                    let additional = self.additional_chroms.contains(&self.ref_names[rid0]);
                    let hist = if tagged { &mut self.ins_tagged } else { &mut self.ins_internal };
                    *hist.entry((additional, isize)).or_insert(0) += 1;
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
            self.rt.rows.push(Row { bc: bc_id, umi: acc.ub.clone().into_bytes(), ge: ge_set, gi: gi_set });
            self.assign.insert(acc.name.clone(), (ge_set, gi_set));
        }
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

// Mean of per-gene max-normalized profiles over the genes of one class (equal weight), RSeQC-style.
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

// Returns (winning genes, was_ambiguous). Ambiguous = a real overlap that a >1-gene tie discarded.
fn resolve(ov: &[(u32, i64)], tagged: bool, strand: char, mapped_len: i64, cfg: &Config, ann: &Annotation) -> (Vec<u32>, bool) {
    if ov.is_empty() {
        return (Vec::new(), false);
    }
    let stranded = cfg.counting_opts.strand == 1 && tagged;
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

fn write_gene_names(out: &Path, ann: &Annotation, project: &str) -> Result<()> {
    let dir = out.join("expression");
    std::fs::create_dir_all(&dir).context("creating expression/")?;
    let mut s = String::new();
    for g in &ann.genes {
        let name = g.name.as_deref().unwrap_or(&g.id);
        let _ = write!(s, "{}\t{}\n", g.id, name);
    }
    std::fs::write(dir.join(format!("{project}.gene_names.txt")), s).context("writing gene_names.txt")?;
    Ok(())
}
