use crate::config::Config;
use crate::count::{ReadTable, Row, NONE};
use crate::gtf::Annotation;
use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_distr::{Binomial, Distribution};
use serde::Deserialize;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Deserialize)]
struct FilterStatsFile {
    per_barcode: HashMap<String, PerBc>,
}
#[derive(Deserialize)]
struct PerBc {
    tagged: u64,
    internal: u64,
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

struct GroupAcc {
    umis: HashMap<Vec<u8>, u32>,
    readcount: u64,
    readcount_internal: u64,
}

pub fn collapse(rt: &ReadTable, cfg: &Config, ann: &Annotation) -> Result<()> {
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

    // ---- "all" level ----
    let mut exon_internal: Option<HashMap<(u32, u32), f64>> = None;
    for (li, &layer) in layers.iter().enumerate() {
        let groups = collapse_rows(rt, &layer, layer_rows[li].iter().copied());
        emit_level(&mut w, rt, ann, layer, "all", &groups)?;
        if layer == Layer::Exon {
            exon_internal = Some(value_map(&groups, Quant::ReadInternal, rt));
        }
    }

    // ---- rpkm (exon/all, from internal reads only) ----
    if let Some(vals) = &exon_internal {
        let fpkm = fpkm_from_internal(vals, ann);
        emit(&mut w, rt, ann, "rpkm", "exon", "all", &fpkm)?;
    }

    // ---- downsampling ----
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
            let groups = collapse_rows(rt, &layer, sampled.into_iter());
            emit_level(&mut w, rt, ann, layer, &level, &groups)?;
        }
    }

    w.finish().context("closing counts.tsv.gz")?;
    Ok(())
}

fn collapse_rows<I: Iterator<Item = usize>>(rt: &ReadTable, layer: &Layer, rows: I) -> HashMap<(u32, u32), GroupAcc> {
    let mut groups: HashMap<(u32, u32), GroupAcc> = HashMap::new();
    for ri in rows {
        let r = &rt.rows[ri];
        let set = layer.set_of(r);
        if set == NONE {
            continue;
        }
        let g = groups
            .entry((r.bc, set))
            .or_insert_with(|| GroupAcc { umis: HashMap::new(), readcount: 0, readcount_internal: 0 });
        g.readcount += 1;
        if r.umi.is_empty() {
            g.readcount_internal += 1;
        } else {
            *g.umis.entry(r.umi.clone()).or_insert(0) += 1;
        }
    }
    groups
}

#[derive(Clone, Copy)]
enum Quant {
    Umi,
    Read,
    ReadInternal,
}

fn value_map(groups: &HashMap<(u32, u32), GroupAcc>, quant: Quant, rt: &ReadTable) -> HashMap<(u32, u32), f64> {
    let mut vals: HashMap<(u32, u32), f64> = HashMap::new();
    for (&(bc, set), g) in groups {
        let v = match quant {
            Quant::Umi => directional(&g.umis) as f64,
            Quant::Read => g.readcount as f64,
            Quant::ReadInternal => g.readcount_internal as f64,
        };
        if v == 0.0 {
            continue;
        }
        // multi-overlap sets share the value equally across their genes
        let genes = rt.set_genes(set);
        let share = v / genes.len() as f64;
        for &gene in genes {
            *vals.entry((gene, bc)).or_default() += share;
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
    groups: &HashMap<(u32, u32), GroupAcc>,
) -> Result<()> {
    for (quant, name) in [(Quant::Umi, "umicount"), (Quant::Read, "readcount"), (Quant::ReadInternal, "readcount_internal")] {
        let vals = value_map(groups, quant, rt);
        emit(w, rt, ann, name, layer.feature(), level, &vals)?;
    }
    Ok(())
}

// One long-format block: (quant, feature, level, gene_id, cell, value) per non-zero entry.
fn emit<W: Write>(
    w: &mut W,
    rt: &ReadTable,
    ann: &Annotation,
    quant: &str,
    feature: &str,
    level: &str,
    vals: &HashMap<(u32, u32), f64>,
) -> Result<()> {
    for (&(gene, bc), &v) in vals {
        writeln!(w, "{quant}\t{feature}\t{level}\t{}\t{}\t{v}", ann.genes[gene as usize].id, rt.barcodes[bc as usize])?;
    }
    Ok(())
}

fn fpkm_from_internal(internal: &HashMap<(u32, u32), f64>, ann: &Annotation) -> HashMap<(u32, u32), f64> {
    let mut libsize: HashMap<u32, f64> = HashMap::new();
    for (&(_, bc), &v) in internal {
        *libsize.entry(bc).or_default() += v;
    }
    let mut fpkm: HashMap<(u32, u32), f64> = HashMap::new();
    for (&(gene, bc), &v) in internal {
        let ls = 1e-6 * libsize[&bc];
        let kb = ann.genes[gene as usize].length as f64 / 1000.0;
        if ls > 0.0 && kb > 0.0 {
            fpkm.insert((gene, bc), v / ls / kb);
        }
    }
    fpkm
}

/// UMI-tools directional (Smith/Heger/Sudbery 2017), edit distance 1. Distinct UMIs are nodes with
/// read counts; a directed edge a->b exists when hamming(a,b)==1 and count[a] >= 2*count[b] - 1.
/// Molecules = number of directed networks: visiting nodes high->low count, each unvisited node seeds
/// a network that absorbs everything reachable along out-edges.
fn directional(umis: &HashMap<Vec<u8>, u32>) -> u64 {
    let n = umis.len();
    if n <= 1 {
        return n as u64;
    }
    let nodes: Vec<(&[u8], u32)> = umis.iter().map(|(k, &v)| (k.as_slice(), v)).collect();
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

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| nodes[b].1.cmp(&nodes[a].1)); // descending count
    let mut visited = vec![false; n];
    let mut networks = 0u64;
    for &start in &order {
        if visited[start] {
            continue;
        }
        networks += 1;
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(u) = stack.pop() {
            for &v in &adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    stack.push(v);
                }
            }
        }
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
