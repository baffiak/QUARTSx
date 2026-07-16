use crate::config::{Config, MultiMapperMode};
use crate::count::{ReadTable, Row, NONE};
use crate::gtf::Annotation;
use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_distr::{Binomial, Distribution};
use rayon::prelude::*;
use crate::filter::PerBc;
use serde::Deserialize;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Deserialize)]
struct FilterStatsFile {
    per_barcode: HashMap<String, PerBc>,
}

#[derive(Clone, Copy, PartialEq)]
enum Layer {
    Exon,
    Intron,
    Inex,
}

impl Layer {
    fn feature(self) -> &'static str {
        match self {
            Layer::Exon => "exon",
            Layer::Intron => "intron",
            Layer::Inex => "inex",
        }
    }
    // gene set for a row under this layer, or NONE if the row doesn't belong to the layer
    fn set_of(self, r: &Row) -> u32 {
        match self {
            Layer::Exon => r.ge,
            Layer::Intron => {
                if r.ge == NONE {
                    r.gi
                } else {
                    NONE
                }
            }
            Layer::Inex => {
                if r.ge != NONE {
                    r.ge
                } else {
                    r.gi
                }
            }
        }
    }
}

// Per-UMI accumulator within one cell: read count plus the running INTERSECTION of the per-fragment
// gene sets of the reads carrying this UMI, kept SEPARATELY by mapping uniqueness so unique evidence
// can dominate. `uniq` = intersection over this UMI's uniquely-mapped (NH==1) reads; `multi` =
// intersection over its multimapping (NH>1) reads. Each is `None` until the first such read; each list
// stays sorted ascending. The molecule's gene set is taken from `uniq` when present, else `multi`.
struct UmiAcc {
    count: u32,
    uniq: Option<Vec<u32>>,
    multi: Option<Vec<u32>>,
}

// Per-(layer, level) collapsed tallies. `mol` is the TAGGED family keyed by the molecule's INTERSECTED
// gene set; `read`/`read_internal` are per-fragment-set fragment tallies (union, no collapse).
#[derive(Default)]
struct LayerCounts {
    mol: HashMap<(u32, Vec<u32>), u64>,      // umicount: (bc, intersected gene set) -> molecules
    read: HashMap<(u32, u32), u64>,          // readcount: (bc, per-fragment set id) -> fragments
    read_internal: HashMap<(u32, u32), u64>, // readcount_internal: (bc, per-fragment set id) -> internal fragments
}

pub fn collapse(rt: &ReadTable, cfg: &Config, ann: &Annotation, stage: &crate::log::Stage) -> Result<()> {
    stage.step("dedup: collapsing UMIs + writing counts");
    let out = Path::new(&cfg.out_dir);
    let expr = out.join("expression");
    std::fs::create_dir_all(&expr).context("creating expression/")?;

    let stats_text = std::fs::read_to_string(out.join("filtered").join("filter_stats.json"))
        .context("reading filter_stats.json")?;
    let stats: FilterStatsFile = serde_json::from_str(&stats_text).context("parsing filter_stats.json")?;

    write_bcumistats(&expr, &stats, &cfg.project)?;

    // per-barcode filtered read totals drive the cell floor (zUMIs nReadsperCell)
    let totals: Vec<u64> = rt
        .barcodes
        .iter()
        .map(|b| stats.per_barcode.get(b).map(|p| p.tagged + p.internal).unwrap_or(0))
        .collect();
    let min_reads = cfg.barcodes.n_reads_per_cell;
    let kept: Vec<bool> = totals.iter().map(|&n| n >= min_reads).collect();

    write_readspercell(&expr, rt, &stats, min_reads, &cfg.project)?;

    // downsampling depth = mapped reads that reached counting per cell
    let depth: Vec<u64> = rt.percell.iter().map(|c| c.mapped()).collect();

    let layers: Vec<Layer> = if cfg.counting_opts.introns {
        vec![Layer::Exon, Layer::Inex, Layer::Intron]
    } else {
        vec![Layer::Exon]
    };

    // pre-split rows per layer to kept cells (also used for downsampling)
    let mut layer_rows: Vec<Vec<usize>> = Vec::with_capacity(layers.len());
    for &layer in &layers {
        let mut v = Vec::new();
        for (i, r) in rt.rows.iter().enumerate() {
            if kept[r.bc as usize] && layer.set_of(r) != NONE {
                v.push(i);
            }
        }
        layer_rows.push(v);
    }

    let f = File::create(expr.join(format!("{}.counts.tsv.gz", cfg.project))).context("creating counts.tsv.gz")?;
    let mut w = GzEncoder::new(BufWriter::new(f), Compression::default());
    writeln!(w, "quant\tfeature\tlevel\tgene\tcell\tvalue")?;

    let mode = cfg.counting_opts.multi_mappers;

    // ---- "all" level ----
    // Each layer emits its integer blocks (umicount/readcount/readcount_internal) AND the STARsolo-mode
    // resolver matrices for BOTH families (umicount=tagged, readcount_internal=internal), mirroring the
    // integer matrix on every axis. The resolver redistributes multi-gene sets by the mode formula;
    // formulas live in count::write_multimapper_matrices.
    let mut exon_internal: Option<BTreeMap<(u32, u32), f64>> = None;
    for (li, &layer) in layers.iter().enumerate() {
        let lc = collapse_rows(rt, &layer, layer_rows[li].iter().copied());
        emit_level(&mut w, rt, ann, layer, "all", &lc)?;
        emit_resolver(&mut w, rt, ann, mode, layer, "all", &lc)?;
        if layer == Layer::Exon {
            exon_internal = Some(set_value_map(&lc.read_internal, rt));
        }
    }

    // ---- rpkm (exon/all, from internal reads only) ----
    if let Some(vals) = &exon_internal {
        let fpkm = fpkm_from_internal(vals, ann);
        emit(&mut w, rt, ann, "rpkm", "exon", "all", &fpkm)?;
    }

    // ---- downsampling ----
    stage.step("dedup: downsampling levels");
    let mut rng = StdRng::seed_from_u64(42);
    let max_n = depth.iter().copied().max().unwrap_or(0);
    for (nmin, nmax, label) in parse_splits(&cfg.counting_opts.downsampling, max_n) {
        let level = format!("downsampled_{label}");
        for (li, &layer) in layers.iter().enumerate() {
            let mut by_cell: HashMap<u32, Vec<usize>> = HashMap::new();
            for &ri in &layer_rows[li] {
                by_cell.entry(rt.rows[ri].bc).or_default().push(ri);
            }
            // deterministic cell order so seeded draws are reproducible
            let mut cells: Vec<u32> = by_cell.keys().copied().collect();
            cells.sort_unstable();
            let mut sampled: Vec<usize> = Vec::new();
            for bc in cells {
                let mut entries = by_cell.remove(&bc).unwrap();
                let n = depth[bc as usize];
                if n < nmin {
                    continue;
                }
                let exn = entries.len() as u64;
                let target = if n <= nmax {
                    exn
                } else {
                    let p = exn as f64 / n as f64;
                    Binomial::new(nmax, p.min(1.0)).unwrap().sample(&mut rng).min(exn)
                };
                entries.shuffle(&mut rng);
                sampled.extend(entries.into_iter().take(target as usize));
            }
            let lc = collapse_rows(rt, &layer, sampled.into_iter());
            emit_level(&mut w, rt, ann, layer, &level, &lc)?;
            emit_resolver(&mut w, rt, ann, mode, layer, &level, &lc)?;
        }
    }

    w.finish().context("closing counts.tsv.gz")?;
    Ok(())
}

/// Collapse the rows of one layer into per-(cell, gene-set) tallies.
///  - `read`/`read_internal`: every fragment contributes ONE unit to its per-fragment (UNION) set,
///    grouped (bc, set); internal fragments (no UMI) additionally feed `read_internal`. No collapse.
///  - `mol` (TAGGED family): each cell's UMI-carrying reads are directional-collapsed into molecules.
///    A molecule's gene set is the INTERSECTION of its member reads' per-fragment sets, but UNIQUE
///    (NH==1) reads dominate: if a molecule has any uniquely-mapped read its set is the
///    intersection of its unique reads only (multimapper reads never alter the integer base and never
///    reintroduce ambiguity); only molecules with NO unique read fall back to their multimapper
///    intersection. An empty intersection drops the molecule from the layer. Cells are independent, so
///    molecule resolution runs in parallel; the integer per-set counts are order-independent to merge.
fn collapse_rows<I: Iterator<Item = usize>>(rt: &ReadTable, layer: &Layer, rows: I) -> LayerCounts {
    let mut read: HashMap<(u32, u32), u64> = HashMap::new();
    let mut read_internal: HashMap<(u32, u32), u64> = HashMap::new();
    let mut tagged: HashMap<u32, HashMap<Vec<u8>, UmiAcc>> = HashMap::new();
    for ri in rows {
        let r = &rt.rows[ri];
        let set = layer.set_of(r);
        if set == NONE {
            continue;
        }
        *read.entry((r.bc, set)).or_default() += 1;
        if r.umi.is_empty() {
            *read_internal.entry((r.bc, set)).or_default() += 1;
        } else {
            let genes = rt.set_genes(set);
            let a = tagged
                .entry(r.bc)
                .or_default()
                .entry(r.umi.clone())
                .or_insert_with(|| UmiAcc { count: 0, uniq: None, multi: None });
            a.count += 1;
            // Fold this read into the intersection of its own uniqueness class, keeping unique and
            // multimapper evidence apart so unique reads can dominate at resolution.
            let slot = if r.uniq { &mut a.uniq } else { &mut a.multi };
            *slot = Some(match slot.take() {
                None => genes.to_vec(),
                Some(cur) => intersect_sorted(&cur, genes),
            });
        }
    }

    let per_cell: Vec<HashMap<(u32, Vec<u32>), u64>> =
        tagged.par_iter().map(|(&bc, umis)| resolve_molecules(bc, umis)).collect();
    let mut mol: HashMap<(u32, Vec<u32>), u64> = HashMap::new();
    for m in per_cell {
        for (k, v) in m {
            *mol.entry(k).or_default() += v;
        }
    }

    LayerCounts { mol, read, read_internal }
}

/// Directional-collapse one cell's UMIs into molecules. A molecule's gene set is the INTERSECTION of its
/// member reads' per-fragment sets, with UNIQUE (NH==1) reads dominating: when the network has any
/// uniquely-mapped read the set is the intersection of its unique reads ONLY (multimapper reads are
/// ignored, so a unique read on gene X is never annihilated by a co-UMI multimapper on {Y,Z}); a network
/// with no unique read falls back to its multimapper intersection. Empty intersections are dropped.
/// Returns (bc, set) tallies.
fn resolve_molecules(bc: u32, umis: &HashMap<Vec<u8>, UmiAcc>) -> HashMap<(u32, Vec<u32>), u64> {
    // Deterministic node order (count desc, UMI asc) pins network formation regardless of HashMap order.
    let mut items: Vec<(&[u8], &UmiAcc)> = umis.iter().map(|(k, v)| (k.as_slice(), v)).collect();
    items.sort_unstable_by(|a, b| b.1.count.cmp(&a.1.count).then_with(|| a.0.cmp(b.0)));
    let nodes: Vec<(&[u8], u32)> = items.iter().map(|&(u, a)| (u, a.count)).collect();

    let mut out: HashMap<(u32, Vec<u32>), u64> = HashMap::new();
    for net in directional_networks(&nodes) {
        // Fold the network's unique and multimapper intersections independently.
        let mut uniq: Option<Vec<u32>> = None;
        let mut multi: Option<Vec<u32>> = None;
        for &i in &net {
            let a = items[i].1;
            if let Some(u) = &a.uniq {
                uniq = Some(match uniq {
                    None => u.clone(),
                    Some(cur) => intersect_sorted(&cur, u),
                });
            }
            if let Some(m) = &a.multi {
                multi = Some(match multi {
                    None => m.clone(),
                    Some(cur) => intersect_sorted(&cur, m),
                });
            }
        }
        // Unique evidence dominates: commit to the unique intersection whenever the network has ANY
        // uniquely-mapped read (even if it annihilates to empty -> molecule dropped, exactly as Unique
        // mode would since it never sees the multimapper reads); only a network with no unique read at
        // all falls back to its multimapper intersection.
        if let Some(gv) = uniq.or(multi) {
            if !gv.is_empty() {
                *out.entry((bc, gv)).or_default() += 1;
            }
        }
    }
    out
}

/// Intersection of two ASCENDING-sorted gene-id slices (result stays sorted).
fn intersect_sorted(a: &[u32], b: &[u32]) -> Vec<u32> {
    let mut out = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

fn umi_value_map(mol: &HashMap<(u32, Vec<u32>), u64>) -> BTreeMap<(u32, u32), f64> {
    // gEu = UNIQUE-ONLY INTEGER base: a molecule contributes ONLY when its intersected gene set
    // is exactly one gene. Multi-gene molecules go solely to the `umicount_mult_<mode>` resolver, never
    // here, so `base + umicount_mult_<mode>` reconstructs UniqueAndMult with no double-count. Each
    // single-gene set {g} is one HashMap key, so each (gene, bc) gets one integer add -> order-safe.
    let mut vals: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    for (k, &c) in mol {
        if k.1.len() == 1 {
            *vals.entry((k.1[0], k.0)).or_default() += c as f64;
        }
    }
    vals
}

fn set_value_map(m: &HashMap<(u32, u32), u64>, rt: &ReadTable) -> BTreeMap<(u32, u32), f64> {
    // gEu = UNIQUE-ONLY INTEGER base: a fragment contributes ONLY when its per-fragment gene set
    // is exactly one gene. Multi-gene (multi-overlap / multimapper) fragments go solely to the
    // `<family>_mult_<mode>` resolver, never the base. Each single-gene set interns to one id, so each
    // (gene, bc) gets one integer add -> order-safe.
    let mut vals: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    for (&(bc, set), &c) in m {
        let genes = rt.set_genes(set);
        if genes.len() == 1 {
            *vals.entry((genes[0], bc)).or_default() += c as f64;
        }
    }
    vals
}

fn emit_level<W: Write>(
    w: &mut W,
    rt: &ReadTable,
    ann: &Annotation,
    layer: Layer,
    level: &str,
    lc: &LayerCounts,
) -> Result<()> {
    emit(w, rt, ann, "umicount", layer.feature(), level, &umi_value_map(&lc.mol))?;
    emit(w, rt, ann, "readcount", layer.feature(), level, &set_value_map(&lc.read, rt))?;
    emit(w, rt, ann, "readcount_internal", layer.feature(), level, &set_value_map(&lc.read_internal, rt))?;
    Ok(())
}

/// Emit the STARsolo-mode resolver matrices for one (layer, level) for BOTH families, mirroring the
/// integer matrix: `umicount` from tagged molecules (INTERSECTED sets), `readcount_internal` from
/// internal fragments (per-fragment sets). Unique mode emits nothing (resolver early-returns).
fn emit_resolver<W: Write>(
    w: &mut W,
    rt: &ReadTable,
    ann: &Annotation,
    mode: MultiMapperMode,
    layer: Layer,
    level: &str,
    lc: &LayerCounts,
) -> Result<()> {
    if mode == MultiMapperMode::Unique {
        return Ok(());
    }
    let tagged: HashMap<(u32, Vec<u32>), f64> =
        lc.mol.iter().map(|(k, &c)| (k.clone(), c as f64)).collect();
    crate::count::write_multimapper_matrices(w, rt, ann, mode, &tagged, "umicount", layer.feature(), level)?;
    let internal: HashMap<(u32, Vec<u32>), f64> = lc
        .read_internal
        .iter()
        .map(|(&(bc, set), &c)| ((bc, rt.set_genes(set).to_vec()), c as f64))
        .collect();
    crate::count::write_multimapper_matrices(w, rt, ann, mode, &internal, "readcount_internal", layer.feature(), level)?;
    Ok(())
}

// One long-format block: (quant, feature, level, gene_id, cell, value) per non-zero entry. Lines are
// emitted in canonical (gene_id string, cell string) order so counts.tsv.gz is byte-identical run to
// run and thread-count-independent (Rust HashMap iteration order is otherwise nondeterministic).
fn emit<W: Write>(
    w: &mut W,
    rt: &ReadTable,
    ann: &Annotation,
    quant: &str,
    feature: &str,
    level: &str,
    vals: &BTreeMap<(u32, u32), f64>,
) -> Result<()> {
    let mut keys: Vec<&(u32, u32)> = vals.keys().collect();
    keys.sort_unstable_by(|a, b| {
        ann.genes[a.0 as usize]
            .id
            .cmp(&ann.genes[b.0 as usize].id)
            .then_with(|| rt.barcodes[a.1 as usize].cmp(&rt.barcodes[b.1 as usize]))
    });
    for k in keys {
        let v = vals[k];
        writeln!(w, "{quant}\t{feature}\t{level}\t{}\t{}\t{v}", ann.genes[k.0 as usize].id, rt.barcodes[k.1 as usize])?;
    }
    Ok(())
}

fn fpkm_from_internal(internal: &BTreeMap<(u32, u32), f64>, ann: &Annotation) -> BTreeMap<(u32, u32), f64> {
    // BTreeMap iteration is sorted, so the libsize FP sum is reduced in a fixed order.
    let mut libsize: BTreeMap<u32, f64> = BTreeMap::new();
    for (&(_, bc), &v) in internal {
        *libsize.entry(bc).or_default() += v;
    }
    let mut fpkm: BTreeMap<(u32, u32), f64> = BTreeMap::new();
    for (&(gene, bc), &v) in internal {
        let ls = 1e-6 * libsize[&bc];
        let kb = ann.genes[gene as usize].length as f64 / 1000.0;
        if ls > 0.0 && kb > 0.0 {
            fpkm.insert((gene, bc), v / ls / kb);
        }
    }
    fpkm
}

/// UMI-tools directional method, edit distance 1. Distinct UMIs are nodes with
/// read counts; a directed edge a->b exists when hamming(a,b)==1 and count[a] >= 2*count[b] - 1.
/// Molecules = directed networks: visiting nodes high->low count, each unvisited node seeds a network
/// that absorbs everything reachable along out-edges. Returns each network's member node indices
/// (into `nodes`, which the caller pre-sorts count-desc/UMI-asc for deterministic membership).
fn directional_networks(nodes: &[(&[u8], u32)]) -> Vec<Vec<usize>> {
    let n = nodes.len();
    if n == 0 {
        return Vec::new();
    }
    let index: HashMap<&[u8], usize> = nodes.iter().enumerate().map(|(i, &(k, _))| (k, i)).collect();

    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (i, &(umi, ci)) in nodes.iter().enumerate() {
        for nb in neighbours(umi) {
            if let Some(&j) = index.get(nb.as_slice()) {
                if ci as u64 >= 2 * nodes[j].1 as u64 - 1 {
                    adj[i].push(j);
                }
            }
        }
    }

    let mut visited = vec![false; n];
    let mut networks: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut members = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(u) = stack.pop() {
            members.push(u);
            for &v in &adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    stack.push(v);
                }
            }
        }
        networks.push(members);
    }
    networks
}

/// All single-substitution variants of `umi` (Hamming distance 1).
fn neighbours(umi: &[u8]) -> Vec<Vec<u8>> {
    const BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];
    let mut out = Vec::with_capacity(umi.len() * 4);
    let mut buf = umi.to_vec();
    for p in 0..umi.len() {
        let old = buf[p];
        for &b in &BASES {
            if b != old {
                buf[p] = b;
                out.push(buf.clone());
            }
        }
        buf[p] = old;
    }
    out
}

fn write_bcumistats(expr: &Path, stats: &FilterStatsFile, project: &str) -> Result<()> {
    let mut bcs: Vec<&String> = stats.per_barcode.keys().collect();
    bcs.sort();
    let mut s = String::from("XC\tnNontagged\tnUMItag\n");
    for xc in bcs {
        let p = &stats.per_barcode[xc];
        let _ = write!(s, "{xc}\t{}\t{}\n", p.internal, p.tagged);
    }
    std::fs::write(expr.join(format!("{project}.BCUMIstats.txt")), s).context("writing BCUMIstats.txt")?;
    Ok(())
}

fn write_readspercell(expr: &Path, rt: &ReadTable, stats: &FilterStatsFile, min_reads: u64, project: &str) -> Result<()> {
    let percell: HashMap<&str, &crate::count::CellReads> =
        rt.barcodes.iter().enumerate().map(|(i, b)| (b.as_str(), &rt.percell[i])).collect();

    let types = ["Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped", "User"];
    let mut bad = [0u64; 6];
    let mut s = String::from("RG\tN\ttype\n");

    let mut bcs: Vec<&String> = stats.per_barcode.keys().collect();
    bcs.sort();
    for xc in bcs {
        let f = &stats.per_barcode[xc];
        let total = f.tagged + f.internal;
        let c = percell.get(xc.as_str()).copied().cloned().unwrap_or_default();
        let unmapped = total.saturating_sub(c.mapped());
        let counts = [c.exon, c.intron, c.intergenic, c.ambiguity, unmapped, c.user];
        if total >= min_reads {
            for (n, ty) in counts.iter().zip(types) {
                let _ = write!(s, "{xc}\t{n}\t{ty}\n");
            }
        } else {
            for (i, n) in counts.iter().enumerate() {
                bad[i] += n;
            }
        }
    }
    for (n, ty) in bad.iter().zip(types) {
        let _ = write!(s, "bad\t{n}\t{ty}\n");
    }
    std::fs::write(expr.join(format!("{project}.readspercell.txt")), s).context("writing readspercell.txt")?;
    Ok(())
}

/// Parse downsampling into (min, max, label) splits. "0" => none; "x" => [x,x]; "a-b" => [a,b].
/// Drop any split whose min exceeds the deepest cell.
fn parse_splits(down: &str, max_n: u64) -> Vec<(u64, u64, String)> {
    let mut out = Vec::new();
    for tok in down.split(',') {
        let tok = tok.trim();
        if tok.is_empty() || tok == "0" {
            continue;
        }
        let (a, b) = if let Some((x, y)) = tok.split_once('-') {
            (x.parse().unwrap_or(0), y.parse().unwrap_or(0))
        } else {
            let x = tok.parse().unwrap_or(0);
            (x, x)
        };
        if a > max_n {
            continue;
        }
        out.push((a, b, tok.to_string()));
    }
    out
}
