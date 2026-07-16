use crate::config::ReadFiltering;
use anyhow::{Context, Result};

pub struct TrimParams {
    pub adapters: Vec<Vec<u8>>,
    pub quality: u8,
    pub min_length: usize,
}

impl TrimParams {
    pub fn load(rf: &ReadFiltering) -> Result<TrimParams> {
        let adapters = match &rf.adapter_fasta {
            Some(path) => read_fasta_seqs(path)?,
            None => Vec::new(),
        };
        Ok(TrimParams { adapters, quality: rf.quality, min_length: rf.min_length })
    }
}

fn read_fasta_seqs(path: &str) -> Result<Vec<Vec<u8>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading adapter_fasta {path}"))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('>') {
            continue;
        }
        out.push(l.as_bytes().to_vec());
    }
    Ok(out)
}

#[inline]
fn phred(q: u8) -> u8 {
    q.saturating_sub(crate::PHRED_OFFSET)
}

/// Adapter clip -> both-ends quality trim -> 4 bp sliding window, all at `quality`.
/// Returns the kept slice [start, end) of the original read, or None if shorter than min_length.
pub fn trim(seq: &[u8], qual: &[u8], p: &TrimParams) -> Option<(usize, usize)> {
    let mut start = 0usize;
    let mut end = seq.len();
    for a in &p.adapters {
        end = end.min(clip_adapter(&seq[..end], a));
    }

    while start < end && phred(qual[start]) < p.quality {
        start += 1;
    }
    while end > start && phred(qual[end - 1]) < p.quality {
        end -= 1;
    }
    let win = 4usize;
    let mut i = start;
    while i + win <= end {
        let sum: u32 = qual[i..i + win].iter().map(|&q| phred(q) as u32).sum();
        if sum < win as u32 * p.quality as u32 {
            end = i;
            break;
        }
        i += 1;
    }

    if end - start >= p.min_length {
        Some((start, end))
    } else {
        None
    }
}

/// Leftmost 3' adapter alignment (<=10% mismatches over the overlap, min 5 bp seed).
fn clip_adapter(seq: &[u8], adapter: &[u8]) -> usize {
    let n = seq.len();
    let m = adapter.len();
    if m == 0 {
        return n;
    }
    for i in 0..n {
        let overlap = m.min(n - i);
        if overlap < 5 {
            break;
        }
        let mm = (0..overlap).filter(|&k| seq[i + k] != adapter[k]).count();
        if (mm as f64) <= 0.1 * overlap as f64 {
            return i;
        }
    }
    n
}
