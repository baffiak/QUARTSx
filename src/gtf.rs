use anyhow::{Context, Result};
use coitrees::{COITree, Interval, IntervalTree};
use std::collections::{HashMap, HashSet};

pub type Tree = COITree<u32, u32>;

pub struct Gene {
    pub id: String,
    pub name: Option<String>,
    pub strand: char,
    pub chrom: String,
    pub exons: Vec<(i64, i64)>, // 1-based inclusive, strand-reduced union, sorted ascending
    pub length: i64,            // sum of exon widths (FPKM gene length)
    pub additional: bool,       // contig from additional_fasta (GTF source "User")
}

pub struct Annotation {
    pub genes: Vec<Gene>,
    pub exon_tree: HashMap<String, Tree>,
    pub intron_tree: HashMap<String, Tree>,
}

struct GeneBuild {
    id: String,
    chrom: String,
    strand: char,
    name: Option<String>,
    additional: bool,
    exons: Vec<(i64, i64)>,
}

fn attr(attrs: &str, key: &str) -> Option<String> {
    let pat = format!("{key} \"");
    let start = attrs.find(&pat)? + pat.len();
    let rest = &attrs[start..];
    let end = rest.find('"')?;
    Some(rest[..end].to_string())
}

pub fn build(gtf: &str, introns: bool) -> Result<Annotation> {
    let text = std::fs::read_to_string(gtf).with_context(|| format!("reading GTF {gtf}"))?;
    let mut index: HashMap<String, usize> = HashMap::new();
    let mut builds: Vec<GeneBuild> = Vec::new();

    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let gid = match attr(f[8], "gene_id") {
            Some(g) => g,
            None => continue,
        };
        let start: i64 = f[3].parse().with_context(|| format!("bad exon start: {line}"))?;
        let end: i64 = f[4].parse().with_context(|| format!("bad exon end: {line}"))?;
        let strand = f[6].chars().next().unwrap_or('+');
        let idx = *index.entry(gid.clone()).or_insert_with(|| {
            builds.push(GeneBuild {
                id: gid.clone(),
                chrom: f[0].to_string(),
                strand,
                name: attr(f[8], "gene_name"),
                additional: f[1] == "User",
                exons: Vec::new(),
            });
            builds.len() - 1
        });
        builds[idx].exons.push((start, end));
    }

    // builds are in first-appearance order, so gene indices stay stable
    let mut genes = Vec::with_capacity(builds.len());
    for mut b in builds {
        b.exons.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::new();
        for &(s, e) in &b.exons {
            if let Some(last) = merged.last_mut() {
                if s <= last.1 + 1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            merged.push((s, e));
        }
        let length = merged.iter().map(|&(s, e)| e - s + 1).sum();
        genes.push(Gene {
            id: b.id,
            name: b.name,
            strand: b.strand,
            chrom: b.chrom,
            exons: merged,
            length,
            additional: b.additional,
        });
    }

    let exon_tree = build_exon_tree(&genes);
    let intron_tree = if introns {
        build_intron_tree(&genes)
    } else {
        HashMap::new()
    };

    Ok(Annotation { genes, exon_tree, intron_tree })
}

fn build_exon_tree(genes: &[Gene]) -> HashMap<String, Tree> {
    let mut per_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();
    for (gi, g) in genes.iter().enumerate() {
        let v = per_chrom.entry(g.chrom.clone()).or_default();
        for &(s, e) in &g.exons {
            v.push(Interval::new(s as i32, e as i32, gi as u32));
        }
    }
    per_chrom.into_iter().map(|(c, iv)| (c, COITree::new(&iv))).collect()
}

fn build_intron_tree(genes: &[Gene]) -> HashMap<String, Tree> {
    // group exon intervals and gene spans by chrom
    let mut exons_by_chrom: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    let mut spans_by_chrom: HashMap<String, Vec<(i64, i64, u32)>> = HashMap::new();
    for (gi, g) in genes.iter().enumerate() {
        if g.exons.is_empty() {
            continue;
        }
        let ev = exons_by_chrom.entry(g.chrom.clone()).or_default();
        for &(s, e) in &g.exons {
            ev.push((s, e));
        }
        let gs = g.exons.first().unwrap().0;
        let ge = g.exons.last().unwrap().1;
        spans_by_chrom.entry(g.chrom.clone()).or_default().push((gs, ge, gi as u32));
    }

    let mut per_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();
    for (chrom, mut exs) in exons_by_chrom {
        // union of all exons on the chrom -> gaps between them are candidate introns
        exs.sort_unstable();
        let mut union: Vec<(i64, i64)> = Vec::new();
        for (s, e) in exs {
            if let Some(last) = union.last_mut() {
                if s <= last.1 + 1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            union.push((s, e));
        }
        let mut gaps: Vec<(i64, i64)> = Vec::new();
        for w in union.windows(2) {
            let (gs, ge) = (w[0].1 + 1, w[1].0 - 1);
            if gs <= ge {
                gaps.push((gs, ge));
            }
        }

        let spans = spans_by_chrom.get(&chrom).cloned().unwrap_or_default();
        let regions = single_gene_regions(&spans); // (start, end, gene) disjoint, sorted by start

        let out = per_chrom.entry(chrom).or_default();
        for &(gs, ge) in &gaps {
            let width = ge - gs + 1;
            if width <= 10 || width >= 100_000 {
                continue;
            }
            if let Some(gi) = containing_gene(&regions, gs, ge) {
                out.push(Interval::new(gs as i32, ge as i32, gi));
            }
        }
    }
    per_chrom.into_iter().map(|(c, iv)| (c, COITree::new(&iv))).collect()
}

/// Maximal intervals covered by exactly one gene span (disjoint, sorted). Sweep over span edges.
fn single_gene_regions(spans: &[(i64, i64, u32)]) -> Vec<(i64, i64, u32)> {
    if spans.is_empty() {
        return Vec::new();
    }
    let mut events: Vec<(i64, i32, u32)> = Vec::new();
    for &(s, e, g) in spans {
        events.push((s, 1, g));
        events.push((e + 1, -1, g));
    }
    events.sort_unstable_by_key(|&(p, _, _)| p);

    let mut regions: Vec<(i64, i64, u32)> = Vec::new();
    let mut active: HashSet<u32> = HashSet::new();
    let mut i = 0;
    while i < events.len() {
        let pos = events[i].0;
        while i < events.len() && events[i].0 == pos {
            let (_, d, g) = events[i];
            if d == 1 {
                active.insert(g);
            } else {
                active.remove(&g);
            }
            i += 1;
        }
        let next = if i < events.len() { events[i].0 } else { pos };
        if next > pos && active.len() == 1 {
            let g = *active.iter().next().unwrap();
            let (lo, hi) = (pos, next - 1);
            if let Some(last) = regions.last_mut() {
                if last.2 == g && last.1 + 1 == lo {
                    last.1 = hi;
                    continue;
                }
            }
            regions.push((lo, hi, g));
        }
    }
    regions
}

fn containing_gene(regions: &[(i64, i64, u32)], gs: i64, ge: i64) -> Option<u32> {
    let mut lo = 0usize;
    let mut hi = regions.len();
    while lo < hi {
        let mid = (lo + hi) / 2;
        if regions[mid].1 < gs {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo < regions.len() && regions[lo].0 <= gs && regions[lo].1 >= ge {
        Some(regions[lo].2)
    } else {
        None
    }
}
