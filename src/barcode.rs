use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};

/// Sentinel id stored in the sphere map for windows reachable from >=2 distinct indexes.
/// A read landing on such a window is rejected (unassigned), never mis-corrected.
pub const REJECT: u32 = u32::MAX;

/// Edit ceiling for per-index correction. `max_total` is the yaml `BarcodeBinning`
/// (default/expected 2); `max_indel` is hard-fixed to 1 by the caller.
#[derive(Clone, Copy)]
pub struct EditBudget {
    pub max_total: u32,
    pub max_indel: u32,
}

/// Read-vs-table strand for a single index column, auto-detected per index.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orient {
    Forward,
    RevComp,
}

// ----------------------------------------------------------------------------
// Packing + sequence helpers
// ----------------------------------------------------------------------------

#[inline]
fn base2(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Pack an ACGT sequence into 2 bits/base, MSB-first. `None` if any base is non-ACGT.
/// Supports L<=32; the index windows here are L<=12 (<=24 bits).
#[inline]
fn pack(seq: &[u8]) -> Option<u64> {
    let mut k = 0u64;
    for &b in seq {
        k = (k << 2) | base2(b)?;
    }
    Some(k)
}

#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement `seq` into `dst` (cleared first); any non-ACGT base becomes 'N'.
/// Allocation-free when `dst` already has capacity — the hot path hands in per-worker scratch so
/// the orientation flip costs no per-read heap traffic (§1 alloc-churn fix).
#[inline]
pub fn revcomp_into(seq: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.reserve(seq.len());
    dst.extend(seq.iter().rev().map(|&b| comp(b)));
}

/// Reverse complement over ACGT; any non-ACGT base becomes 'N'. Convenience allocator over
/// `revcomp_into` (kept for tests and the cold orientation-vote path).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    let mut v = Vec::with_capacity(seq.len());
    revcomp_into(seq, &mut v);
    v
}

// ----------------------------------------------------------------------------
// Sequence-Levenshtein distance
// ----------------------------------------------------------------------------

/// Levenshtein DP with a FREE (unpenalized) 3' / right end on BOTH sequences: the returned
/// distance is the minimum over the final DP row `min_j D[n][j]` AND the final DP column
/// `min_i D[i][m]`. Because the index is 5'-anchored and the instrument emits a fixed-length
/// read, trailing bases of either sequence beyond a common 5' prefix match cost nothing, so a
/// single boundary indel scores exactly 1 REGARDLESS of argument order (the metric is symmetric).
/// Used only for the panel guard (min pairwise distance), never on the hot path.
pub fn seqlev(a: &[u8], b: &[u8]) -> u32 {
    let n = a.len();
    let m = b.len();
    let mut prev: Vec<u32> = (0..=m as u32).collect();
    let mut cur = vec![0u32; m + 1];
    // Track the last-column minimum (free right end of `a`): start with D[0][m] = m.
    let mut col_min = prev[m];
    for i in 1..=n {
        cur[0] = i as u32;
        for j in 1..=m {
            let sub = prev[j - 1] + if a[i - 1] == b[j - 1] { 0 } else { 1 };
            let del = prev[j] + 1; // delete from a
            let ins = cur[j - 1] + 1; // insert into a
            cur[j] = sub.min(del).min(ins);
        }
        col_min = col_min.min(cur[m]); // D[i][m]
        std::mem::swap(&mut prev, &mut cur);
    }
    // free right end of BOTH: min over the final row (now in `prev` after the swap) and the
    // final column tracked in `col_min`.
    let row_min = *prev.iter().min().unwrap_or(&0);
    row_min.min(col_min)
}

// ----------------------------------------------------------------------------
// IndexList — one per index column (i7, i5)
// ----------------------------------------------------------------------------

pub struct IndexList {
    /// Table-orientation index strings; id = position in this vec.
    pub codes: Vec<String>,
    /// Packed exact L-mer -> id. Authoritative fast path, never holds REJECT.
    exact: HashMap<u64, u32>,
    /// Packed L-window -> id OR REJECT. Indel/substitution recovery.
    sphere: HashMap<u64, u32>,
    /// Minimum pairwise Sequence-Levenshtein distance over `codes`.
    pub min_dist: u32,
    /// Safely correctable radius: floor((min_dist - 1) / 2).
    pub safe_radius: u32,
    /// (id_a, id_b, seqlev) for pairs with seqlev <= 2*max_total (panel report).
    pub collisions: Vec<(usize, usize, u32)>,
}

impl IndexList {
    pub fn build(codes: Vec<String>, budget: EditBudget) -> Result<IndexList> {
        if codes.is_empty() {
            bail!("empty index list");
        }
        // Uniform-length guard (local only; the length is not stored — decode() packs windows of
        // any L<=32 and never validates against a stored column length, so the field was dead).
        let len = codes[0].len();
        for c in &codes {
            if c.len() != len {
                bail!("non-uniform index length: {c} is {} bp, expected {len}", c.len());
            }
            if pack(c.as_bytes()).is_none() {
                bail!("non-ACGT index: {c}");
            }
        }

        // exact map + duplicate detection
        let mut exact: HashMap<u64, u32> = HashMap::new();
        for (id, c) in codes.iter().enumerate() {
            let k = pack(c.as_bytes()).unwrap();
            if let Some(&other) = exact.get(&k) {
                bail!("duplicate/zero-distance index: {} == {}", codes[other as usize], c);
            }
            exact.insert(k, id as u32);
        }

        // pairwise min distance + collision report
        let two_t = 2 * budget.max_total;
        let mut min_dist = u32::MAX;
        let mut collisions = Vec::new();
        for i in 0..codes.len() {
            for j in (i + 1)..codes.len() {
                let d = seqlev(codes[i].as_bytes(), codes[j].as_bytes());
                if d == 0 {
                    bail!("duplicate/zero-distance index: {} == {}", codes[i], codes[j]);
                }
                if d < min_dist {
                    min_dist = d;
                }
                if d <= two_t {
                    collisions.push((i, j, d));
                }
            }
        }
        if codes.len() == 1 {
            min_dist = u32::MAX;
        }
        let safe_radius = if min_dist == u32::MAX { u32::MAX } else { (min_dist - 1) / 2 };

        // sphere enumeration
        let mut sphere: HashMap<u64, u32> = HashMap::new();
        for (id, c) in codes.iter().enumerate() {
            let id = id as u32;
            enumerate_neighbors(c.as_bytes(), budget, &mut |window: &[u8]| {
                if let Some(k) = pack(window) {
                    register(&mut sphere, k, id);
                }
            });
        }

        Ok(IndexList { codes, exact, sphere, min_dist, safe_radius, collisions })
    }

    /// Thin wrapper over `decode_full` collapsing reject and miss to `None`; used by the sphere
    /// unit tests. `obs` must already be in table orientation.
    #[cfg(test)]
    #[inline]
    pub fn decode(&self, obs: &[u8]) -> Option<u32> {
        match self.decode_full(obs) {
            Decode::Hit { id, .. } => Some(id),
            Decode::Reject | Decode::Miss => None,
        }
    }

    /// Full decode outcome — distinguishes an exact hit, a sphere-corrected hit, an ambiguous
    /// REJECT-sentinel window, and a plain miss, and reports the exact-vs-corrected flag inline so
    /// `assign_into` needs no separate `exact_hit` re-pack (§1) and can split drop reasons (§6).
    #[inline]
    fn decode_full(&self, obs: &[u8]) -> Decode {
        let k = match pack(obs) {
            Some(k) => k,
            None => return Decode::Miss, // non-ACGT window
        };
        if let Some(&id) = self.exact.get(&k) {
            return Decode::Hit { id, exact: true };
        }
        match self.sphere.get(&k) {
            Some(&REJECT) => Decode::Reject,
            Some(&id) => Decode::Hit { id, exact: false },
            None => Decode::Miss,
        }
    }
}

/// Internal decode outcome for a single index column (see `IndexList::decode_full`).
#[derive(Clone, Copy)]
enum Decode {
    /// Resolved to index `id`; `exact` = matched the authoritative exact map (no correction).
    Hit { id: u32, exact: bool },
    /// Window sits on a REJECT sentinel — reachable from >=2 indexes; refuse to guess.
    Reject,
    /// Non-ACGT, or outside the edit budget of every index.
    Miss,
}

/// Registration rule for the sphere map — unambiguous by construction:
/// absent -> insert id; same id -> keep; different non-REJECT id -> overwrite with REJECT;
/// already REJECT -> keep REJECT. Never last-writer-wins.
#[inline]
fn register(sphere: &mut HashMap<u64, u32>, k: u64, id: u32) {
    match sphere.get(&k) {
        None => {
            sphere.insert(k, id);
        }
        Some(&existing) => {
            if existing != id && existing != REJECT {
                sphere.insert(k, REJECT);
            }
        }
    }
}

/// Enumerate every fixed-length L-window reachable from `code` within the budget and hand each
/// to `emit`. Reachable set (max_total=2, max_indel=1): exact; <=2 subs; 1 indel + <=1 sub.
/// Deletion: drop a base, fill the last window position with each of A/C/G/T (the free base a
/// fixed-length sequencer emits). Insertion: insert a base, truncate the right overhang.
fn enumerate_neighbors(code: &[u8], budget: EditBudget, emit: &mut dyn FnMut(&[u8])) {
    let l = code.len();
    let s_free = budget.max_total; // subs allowed when no indel
    let s_after_indel = budget.max_total.saturating_sub(1); // subs alongside the one indel

    // 0 indels: <= s_free substitutions on the L-mer (includes d=0 exact).
    sub_variants(code, s_free, emit);

    if budget.max_indel < 1 {
        return;
    }

    // 1 deletion: drop position p -> L-1 bases, append each free fill base at position L.
    // Layer <= s_after_indel subs on the first L-1 retained positions only (never the fill base).
    for p in 0..l {
        let mut base: Vec<u8> = Vec::with_capacity(l);
        base.extend_from_slice(&code[..p]);
        base.extend_from_slice(&code[p + 1..]); // L-1 retained bases
        for &x in b"ACGT" {
            let mut w = base.clone();
            w.push(x); // free fill base at position L
            // subs only on the first L-1 positions
            sub_variants_prefix(&w, l - 1, s_after_indel, emit);
        }
    }

    // 1 insertion: insert base b at position q (length L+1), truncate to first L bases.
    // Layer <= s_after_indel subs on the resulting L window.
    for q in 0..=l {
        for &b in b"ACGT" {
            let mut w: Vec<u8> = Vec::with_capacity(l + 1);
            w.extend_from_slice(&code[..q]);
            w.push(b);
            w.extend_from_slice(&code[q..]);
            w.truncate(l);
            sub_variants(&w, s_after_indel, emit);
        }
    }
}

/// Emit every window at Hamming distance 0..=max_subs from `seq` (subs over all positions).
fn sub_variants(seq: &[u8], max_subs: u32, emit: &mut dyn FnMut(&[u8])) {
    sub_variants_prefix(seq, seq.len(), max_subs, emit);
}

/// Emit every window formed by substituting <= max_subs of the first `sub_span` positions of
/// `seq` (positions >= sub_span are left untouched — used to protect the free deletion fill base).
fn sub_variants_prefix(seq: &[u8], sub_span: usize, max_subs: u32, emit: &mut dyn FnMut(&[u8])) {
    let mut w = seq.to_vec();
    emit(&w); // d = 0
    if max_subs == 0 || sub_span == 0 {
        return;
    }
    recurse_subs(&mut w, 0, sub_span, max_subs, emit);
}

fn recurse_subs(w: &mut Vec<u8>, start: usize, span: usize, remaining: u32, emit: &mut dyn FnMut(&[u8])) {
    if remaining == 0 {
        return;
    }
    for pos in start..span {
        let orig = w[pos];
        for &b in b"ACGT" {
            if b == orig {
                continue;
            }
            w[pos] = b;
            emit(w);
            recurse_subs(w, pos + 1, span, remaining - 1, emit);
        }
        w[pos] = orig;
    }
}

// ----------------------------------------------------------------------------
// CSV parser (shared) + dims probe
// ----------------------------------------------------------------------------

struct ParsedTable {
    i7: Vec<String>,
    i5: Vec<String>,
    pairs: Vec<(u32, u32)>,
}

/// Table dimensions for config-time validation (charset/length/columns checked here).
/// Only the two column lengths are surfaced — config::validate cross-checks them against the BC
/// slice geometry. The row/index counts (`n_i7`/`n_i5`/`n_pairs`) were never read by any consumer
/// and were removed (§6 dead-code cleanup); tests that need those counts call `parse_table`.
pub struct TableDims {
    pub i7_len: usize,
    pub i5_len: usize,
}

/// Cheap dims/charset/columns probe (config::validate). Calls the full parser.
pub fn probe_table(path: &str) -> Result<TableDims> {
    let t = parse_table(path)?;
    Ok(TableDims {
        i7_len: t.i7.first().map(|s| s.len()).unwrap_or(0),
        i5_len: t.i5.first().map(|s| s.len()).unwrap_or(0),
    })
}

fn detect_delim(header: &str) -> Result<char> {
    let mut best = (' ', 0usize);
    for d in [',', '\t', ';'] {
        let c = header.chars().filter(|&ch| ch == d).count();
        if c > best.1 {
            best = (d, c);
        }
    }
    if best.1 == 0 {
        bail!("cannot detect delimiter (no comma/tab/semicolon in header)");
    }
    Ok(best.0)
}

fn parse_table(path: &str) -> Result<ParsedTable> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading index table {path}"))?;
    // Strip a leading UTF-8 BOM (Excel on Windows/European locale writes one). str::trim() does
    // NOT remove U+FEFF (it has White_Space=No), so it would otherwise glue onto header[0].
    let text = text.strip_prefix('\u{feff}').map(str::to_string).unwrap_or(text);

    // Defensive line splitting: split on either \n or \r (handles \n, \r\n, and bare-\r files),
    // then trim; drop blank lines. str::lines() would NOT split a bare-\r file.
    let lines: Vec<String> = text
        .split(|c| c == '\n' || c == '\r')
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect();
    if lines.is_empty() {
        bail!("index table {path} is empty");
    }

    let delim = detect_delim(&lines[0])?;

    // header: locate i7_index / i5_index columns case-insensitively, at any position.
    let header: Vec<String> = lines[0].split(delim).map(|c| c.trim().to_lowercase()).collect();
    let i7_col = header
        .iter()
        .position(|c| c == "i7_index")
        .context("index table missing an 'i7_index' column")?;
    let i5_col = header
        .iter()
        .position(|c| c == "i5_index")
        .context("index table missing an 'i5_index' column")?;

    let mut i7_ids: HashMap<String, u32> = HashMap::new();
    let mut i5_ids: HashMap<String, u32> = HashMap::new();
    let mut i7: Vec<String> = Vec::new();
    let mut i5: Vec<String> = Vec::new();
    let mut pair_set: HashSet<(u32, u32)> = HashSet::new();
    let mut pairs: Vec<(u32, u32)> = Vec::new();
    let mut i7_len: Option<usize> = None;
    let mut i5_len: Option<usize> = None;

    let is_acgt = |s: &str| s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'));

    for (n, line) in lines.iter().enumerate().skip(1) {
        let cells: Vec<&str> = line.split(delim).collect();
        if cells.len() <= i7_col.max(i5_col) {
            // short/ragged row (e.g. trailing blank) — skip if it's effectively empty
            if line.chars().all(|c| c == delim || c.is_whitespace()) {
                continue;
            }
            bail!("index table row {n} has too few columns: {line}");
        }
        let c7 = cells[i7_col].trim().to_uppercase();
        let c5 = cells[i5_col].trim().to_uppercase();
        if c7.is_empty() && c5.is_empty() {
            continue;
        }
        if !is_acgt(&c7) {
            bail!("non-ACGT base in i7_index row {n}: {c7}");
        }
        if !is_acgt(&c5) {
            bail!("non-ACGT base in i5_index row {n}: {c5}");
        }
        match i7_len {
            None => i7_len = Some(c7.len()),
            Some(l) if l != c7.len() => {
                bail!("i7_index length mismatch: {c7} is {} bp, expected {l}", c7.len())
            }
            _ => {}
        }
        match i5_len {
            None => i5_len = Some(c5.len()),
            Some(l) if l != c5.len() => {
                bail!("i5_index length mismatch: {c5} is {} bp, expected {l}", c5.len())
            }
            _ => {}
        }

        let id7 = *i7_ids.entry(c7.clone()).or_insert_with(|| {
            i7.push(c7.clone());
            (i7.len() - 1) as u32
        });
        let id5 = *i5_ids.entry(c5.clone()).or_insert_with(|| {
            i5.push(c5.clone());
            (i5.len() - 1) as u32
        });
        if pair_set.insert((id7, id5)) {
            pairs.push((id7, id5));
        }
    }

    if i7.is_empty() || i5.is_empty() {
        bail!("index table {path} has no data rows");
    }
    Ok(ParsedTable { i7, i5, pairs })
}

// ----------------------------------------------------------------------------
// IndexTable — top-level object
// ----------------------------------------------------------------------------

pub struct IndexTable {
    pub i7: IndexList,
    pub i5: IndexList,
    /// Valid (i7_id, i5_id) pairs → pair ordinal. Key packed ((i7_id as u64) << 32) | i5_id.
    /// The ordinal indexes `labels`; it is also returned to callers as a compact integer barcode id
    /// so the hot path can carry a `u32` instead of re-hashing a string (§1 alloc-churn fix).
    pairs: HashMap<u64, u32>,
    /// Interned `i7 ++ i5` label per valid pair, built ONCE at load. `assign_into` hands back a
    /// borrow of the ready string — no per-read `format!` concatenation.
    labels: Vec<String>,
    pub i7_orient: Orient,
    pub i5_orient: Orient,
    /// Panel-guard messages (bcl2fastq-style "distance too small").
    pub warnings: Vec<String>,
}

/// Per-worker revcomp scratch: one owned buffer per index column, reused across reads so the
/// orientation flip in `assign_into` allocates zero heap on the hot path (Forward columns don't
/// even touch it — they borrow the raw window directly). Each filter/encode worker owns one.
#[derive(Default)]
pub struct RevcompScratch {
    o7: Vec<u8>,
    o5: Vec<u8>,
}

impl RevcompScratch {
    pub fn new() -> RevcompScratch {
        RevcompScratch::default()
    }
}

/// Why a read could not be assigned to a cell — lets `filter.rs` split its `no_barcode` bucket into
/// absent/uncorrectable vs ambiguous-reject for drop-reason reporting (§6). `AmbiguousReject` is the
/// distinct, scientifically-meaningful case: the read window is correctable-distance from >=2 listed
/// indexes, so the sphere holds a REJECT sentinel and we refuse to guess (never mis-assign).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum UnassignedReason {
    /// An index window landed on a REJECT sentinel (reachable from >=2 indexes) — ambiguous.
    AmbiguousReject,
    /// An index window is non-ACGT, or not within the edit budget of ANY listed index (a plain
    /// miss: absent / uncorrectable).
    Absent,
    /// Both indexes decoded to real ids, but the (i7, i5) COMBINATION is not a listed sample cell.
    InvalidPair,
}

/// Result of assigning one read's index pair. `Assigned` borrows the interned label from the table.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AssignResult<'a> {
    Assigned {
        /// Compact pair id (ordinal into the valid-pair list) — a ready integer barcode key.
        pair_id: u32,
        /// Interned `i7 ++ i5` table-orientation label (borrowed; no per-read allocation).
        label: &'a str,
        /// True iff either index needed the sphere (i.e. not both an exact hit).
        corrected: bool,
    },
    Unassigned(UnassignedReason),
}

#[inline]
fn pair_key(id7: u32, id5: u32) -> u64 {
    ((id7 as u64) << 32) | id5 as u64
}

impl IndexTable {
    pub fn load(path: &str, budget: EditBudget) -> Result<IndexTable> {
        let parsed = parse_table(path)?;
        let pair_list = parsed.pairs; // move out before the partial moves of i7/i5 below
        let i7 = IndexList::build(parsed.i7, budget)?;
        let i5 = IndexList::build(parsed.i5, budget)?;
        // Intern one label per valid pair ONCE (ordinal = insertion order over the parsed rows), and
        // map its packed key -> ordinal. Kills the per-read `format!` on the assign hot path.
        let mut pairs: HashMap<u64, u32> = HashMap::with_capacity(pair_list.len());
        let mut labels: Vec<String> = Vec::with_capacity(pair_list.len());
        for &(a, b) in &pair_list {
            let key = pair_key(a, b);
            if let std::collections::hash_map::Entry::Vacant(e) = pairs.entry(key) {
                e.insert(labels.len() as u32);
                let mut label = String::with_capacity(i7.codes[a as usize].len() + i5.codes[b as usize].len());
                label.push_str(&i7.codes[a as usize]);
                label.push_str(&i5.codes[b as usize]);
                labels.push(label);
            }
        }

        let mut warnings = Vec::new();
        for (name, list) in [("i7", &i7), ("i5", &i5)] {
            if budget.max_total > list.safe_radius {
                let n = list.collisions.len();
                let example = list
                    .collisions
                    .first()
                    .map(|&(a, b, d)| format!("{}<->{} (d={d})", list.codes[a], list.codes[b]))
                    .unwrap_or_else(|| "n/a".to_string());
                warnings.push(format!(
                    "{name} panel min Sequence-Levenshtein distance = {}; requested budget {} exceeds safe radius {}. {} colliding index pair(s), e.g. {}. Reads landing in overlaps are rejected (reject sentinels), not misassigned — reduce BarcodeBinning to {} to correct fully.",
                    list.min_dist, budget.max_total, list.safe_radius, n, example, list.safe_radius
                ));
            }
        }

        Ok(IndexTable {
            i7,
            i5,
            pairs,
            labels,
            i7_orient: Orient::Forward,
            i5_orient: Orient::Forward,
            warnings,
        })
    }

    /// Sample-driven per-index orientation detection; sets `i7_orient`/`i5_orient`.
    /// Bails on low-confidence or mixed orientation.
    pub fn detect_and_set_orientation(
        &mut self,
        i7_samples: &[Vec<u8>],
        i5_samples: &[Vec<u8>],
    ) -> Result<()> {
        let v7 = vote_orientation(i7_samples, &self.i7);
        self.i7_orient = resolve_orientation(&v7, "i7")?;
        let v5 = vote_orientation(i5_samples, &self.i5);
        self.i5_orient = resolve_orientation(&v5, "i5")?;
        Ok(())
    }

    /// Hot-path assign: correct i7 and i5 separately (applying detected orientation via per-worker
    /// `scratch`, so a Forward column borrows the raw window with ZERO copies and a RevComp column
    /// reuses the scratch buffer), validate the pair, and hand back the interned label + pair id.
    ///
    /// Zero per-read heap allocation on the assign path (§1): no revcomp `Vec` churn, no `format!`.
    /// The decode reason is preserved so `filter.rs` can split drop reasons (§6):
    /// a REJECT-sentinel hit (ambiguous, reachable from >=2 indexes in the Sequence-Levenshtein
    /// sphere) is `AmbiguousReject`; anything else non-assigned is `Absent`
    /// (plain miss) or `InvalidPair` (both indexes real, combination not a listed cell).
    #[inline]
    pub fn assign_into<'a>(
        &'a self,
        raw_i7: &[u8],
        raw_i5: &[u8],
        scratch: &mut RevcompScratch,
    ) -> AssignResult<'a> {
        // Borrow the raw window on Forward; fill (and borrow) reusable scratch on RevComp. The two
        // scratch fields are disjoint so both borrows coexist.
        let o7: &[u8] = match self.i7_orient {
            Orient::Forward => raw_i7,
            Orient::RevComp => {
                revcomp_into(raw_i7, &mut scratch.o7);
                &scratch.o7
            }
        };
        let o5: &[u8] = match self.i5_orient {
            Orient::Forward => raw_i5,
            Orient::RevComp => {
                revcomp_into(raw_i5, &mut scratch.o5);
                &scratch.o5
            }
        };

        let d7 = self.i7.decode_full(o7);
        let d5 = self.i5.decode_full(o5);

        // Ambiguous-reject takes priority over a plain miss on the other column: if EITHER index
        // hit a REJECT sentinel we refuse to guess, and that is the more informative drop reason.
        if matches!(d7, Decode::Reject) || matches!(d5, Decode::Reject) {
            return AssignResult::Unassigned(UnassignedReason::AmbiguousReject);
        }
        let (id7, ex7) = match d7 {
            Decode::Hit { id, exact } => (id, exact),
            _ => return AssignResult::Unassigned(UnassignedReason::Absent),
        };
        let (id5, ex5) = match d5 {
            Decode::Hit { id, exact } => (id, exact),
            _ => return AssignResult::Unassigned(UnassignedReason::Absent),
        };

        match self.pairs.get(&pair_key(id7, id5)) {
            Some(&pair_id) => AssignResult::Assigned {
                pair_id,
                label: &self.labels[pair_id as usize],
                corrected: !(ex7 && ex5),
            },
            None => AssignResult::Unassigned(UnassignedReason::InvalidPair),
        }
    }

    /// Interned label for a pair ordinal returned by `assign_into`; used by the tests.
    #[cfg(test)]
    #[inline]
    pub fn label(&self, pair_id: u32) -> &str {
        &self.labels[pair_id as usize]
    }

    /// Convenience allocator over `assign_into` — allocates a throwaway scratch and owns the label.
    /// Returns `(label, corrected)` or `None` for any unassigned reason; used by the tests (the hot
    /// filter loop uses `assign_into` with per-worker scratch).
    #[cfg(test)]
    pub fn assign_pair(&self, raw_i7: &[u8], raw_i5: &[u8]) -> Option<(String, bool)> {
        let mut scratch = RevcompScratch::new();
        match self.assign_into(raw_i7, raw_i5, &mut scratch) {
            AssignResult::Assigned { label, corrected, .. } => Some((label.to_string(), corrected)),
            AssignResult::Unassigned(_) => None,
        }
    }
}

impl IndexList {
    /// True iff `obs` is an exact (packable) member of the index list.
    #[inline]
    fn exact_hit(&self, obs: &[u8]) -> bool {
        match pack(obs) {
            Some(k) => self.exact.contains_key(&k),
            None => false,
        }
    }
}

// ----------------------------------------------------------------------------
// Orientation detection (palindrome / shared-index tolerant)
// ----------------------------------------------------------------------------

pub struct Vote {
    pub fwd: u64,
    pub rc: u64,
    pub ambig: u64,
    pub unmatched: u64,
}

/// Tally, over sampled raw index windows, exact forward vs revcomp membership.
/// A sample matching BOTH orientations (palindromic, or its revcomp equals a different listed
/// index) is counted `ambig` and EXCLUDED from the decision — never treated as a conflict.
pub fn vote_orientation(samples: &[Vec<u8>], list: &IndexList) -> Vote {
    let mut v = Vote { fwd: 0, rc: 0, ambig: 0, unmatched: 0 };
    for s in samples {
        let fwd_hit = list.exact_hit(s);
        let rc_hit = list.exact_hit(&revcomp(s));
        match (fwd_hit, rc_hit) {
            (true, true) => v.ambig += 1,
            (true, false) => v.fwd += 1,
            (false, true) => v.rc += 1,
            (false, false) => v.unmatched += 1,
        }
    }
    v
}

const MIN_DISCRIM: u64 = 200;
const DOMINANCE: f64 = 0.90;

/// Decide orientation from the discriminating votes. Bails on low confidence (~0 discriminating
/// reads) or strong two-sided support ("mixed").
pub fn resolve_orientation(v: &Vote, which: &str) -> Result<Orient> {
    let discrim = v.fwd + v.rc;
    if discrim < MIN_DISCRIM {
        bail!(
            "low-confidence orientation for {which}: only {discrim} discriminating reads in {} sampled — likely wrong index table or wrong read files (I1/I2 swapped or wrong lengths).",
            v.fwd + v.rc + v.ambig + v.unmatched
        );
    }
    // Decide purely on dominance once we have enough discriminating reads. A near-even or exact
    // tie fails dominance and is reported as "mixed" rather than silently resolved to Forward.
    let major = v.fwd.max(v.rc);
    if (major as f64) / (discrim as f64) < DOMINANCE {
        bail!("mixed orientation for {which} (fwd={}, rc={})", v.fwd, v.rc);
    }
    Ok(if v.fwd >= v.rc { Orient::Forward } else { Orient::RevComp })
}

// ----------------------------------------------------------------------------
// Kept unchanged: fixed-position Hamming of R1's first bases against the TSO tag.
// ----------------------------------------------------------------------------

/// Fixed-position Hamming (NOT edit distance) of R1's first bases against the TSO tag.
pub fn is_tagged(r1: &[u8], tag: &[u8], max_mm: u8) -> bool {
    if r1.len() < tag.len() {
        return false;
    }
    let mut mm = 0u8;
    for i in 0..tag.len() {
        if r1[i] != tag[i] {
            mm += 1;
            if mm > max_mm {
                return false;
            }
        }
    }
    true
}

// ----------------------------------------------------------------------------
// Tests
// ----------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    const BUDGET: EditBudget = EditBudget { max_total: 2, max_indel: 1 };

    fn tmp(name: &str, content: &[u8]) -> String {
        let mut p = std::env::temp_dir();
        p.push(format!("quartsx_bc_test_{}_{}", std::process::id(), name));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(content).unwrap();
        p.to_string_lossy().into_owned()
    }

    // ---- seqlev ----
    #[test]
    fn seqlev_basic() {
        assert_eq!(seqlev(b"ACGTACGTAC", b"ACGTACGTAC"), 0);
        // strict prefix: with a two-sided free right end the trailing base of the longer sequence
        // is unpenalized, so a prefix is indistinguishable at fixed read length -> distance 0.
        assert_eq!(seqlev(b"ACGTACGTAC", b"ACGTACGTA"), 0);
        // one internal deletion (not a prefix; b shorter by 1) = 1
        assert_eq!(seqlev(b"ACGTACGTAC", b"AGTACGTAC"), 1);
        // two substitutions = 2
        assert_eq!(seqlev(b"ACGTACGTAC", b"AGGTACGTAG"), 2);
        // one substitution = 1
        assert_eq!(seqlev(b"ACGTACGTAC", b"AGGTACGTAC"), 1);
    }

    #[test]
    fn seqlev_symmetric_front_indel() {
        // A single left/front indel apart: the true Sequence-Levenshtein distance is 1 and the
        // metric must be symmetric in argument order (two-sided free end).
        assert_eq!(seqlev(b"CGTACGTACA", b"ACGTACGTAC"), 1);
        assert_eq!(seqlev(b"ACGTACGTAC", b"CGTACGTACA"), 1);
    }

    // ---- sphere + reject sentinel ----
    #[test]
    fn sphere_reject_sentinel() {
        // Two codes at Hamming distance 2. A single-sub window between them (AAAAAAAAAC) is
        // reachable from both -> REJECT. Each code still decodes to itself via exact.
        let a = "AAAAAAAAAA".to_string();
        let b = "AAAAAAAACC".to_string();
        let list = IndexList::build(vec![a.clone(), b.clone()], BUDGET).unwrap();
        assert_eq!(list.decode(b"AAAAAAAAAA"), Some(0));
        assert_eq!(list.decode(b"AAAAAAAACC"), Some(1));
        // shared 1-sub window -> reject
        assert_eq!(list.decode(b"AAAAAAAAAC"), None);
    }

    #[test]
    fn sphere_recovers_unique_neighbor() {
        // Well-separated codes: a 1-sub neighbor recovers uniquely.
        let list = IndexList::build(vec!["ACGTACGTAC".into(), "TGCATGCATG".into()], BUDGET).unwrap();
        assert_eq!(list.decode(b"ACGTACGTAC"), Some(0)); // exact
        assert_eq!(list.decode(b"ACGTACGTAG"), Some(0)); // 1 sub of code 0
        assert_eq!(list.decode(b"TGCATGCATA"), Some(1)); // 1 sub of code 1
    }

    // ---- budget enforcement ----
    #[test]
    fn budget_enforcement() {
        // Single isolated code; sphere built around it, no collisions.
        let code = "ACGTACGTAC";
        let list = IndexList::build(vec![code.into()], BUDGET).unwrap();

        // exact
        assert_eq!(list.decode(b"ACGTACGTAC"), Some(0));
        // 2 substitutions (pos 2 and 8): C->A, A->C
        assert_eq!(list.decode(b"ACATACGTCC"), Some(0));
        // 1 deletion + 1 sub: delete pos0 'A' -> "CGTACGTAC" + free fill 'X', then sub one base.
        // observed with a boundary deletion then one mismatch. "CGTACGTACA" is del+fill; sub pos1.
        assert_eq!(list.decode(b"CGTACGTACA"), Some(0)); // 1 indel, 0 sub
        assert_eq!(list.decode(b"CTTACGTACA"), Some(0)); // 1 indel + 1 sub (pos1 G->T)
        // 1 indel + 2 subs must NOT decode
        assert_eq!(list.decode(b"CTTTCGTACA"), None);
    }

    #[test]
    fn budget_three_subs_rejected() {
        let list = IndexList::build(vec!["ACGTACGTAC".into()], BUDGET).unwrap();
        // exactly 3 substitutions from ACGTACGTAC -> AAATACGTAG (pos1 C->A? let's build cleanly)
        // base ACGTACGTAC; change pos0 A->C, pos1 C->A, pos2 G->A => "CAATACGTAC" (3 subs)
        assert_eq!(list.decode(b"CAATACGTAC"), None);
        // 2 subs of that (pos0,pos1) should decode
        assert_eq!(list.decode(b"CAGTACGTAC"), Some(0));
    }

    #[test]
    fn budget_two_indels_rejected() {
        let list = IndexList::build(vec!["ACGTACGTAC".into()], BUDGET).unwrap();
        // two deletions can't be represented within a single-indel sphere -> None for a
        // window two boundary-deletions away, e.g. drop first two bases + two free fills.
        assert_eq!(list.decode(b"GTACGTACAA"), None);
    }

    // ---- orientation vote ----
    #[test]
    fn orient_vote_counts_and_resolve() {
        // rc-only reads against a single-code list: read whose fwd is NOT in the list but whose
        // revcomp IS -> counted `rc`.
        let rc_list = IndexList::build(vec!["ACGTAAGGCC".into()], BUDGET).unwrap();
        let rc_read = revcomp(b"ACGTAAGGCC");
        assert!(!rc_list.exact_hit(&rc_read));
        assert!(rc_list.exact_hit(&revcomp(&rc_read)));

        // A list containing both a code and its revcomp -> a read equal to the code hits BOTH
        // orientations -> ambiguous (excluded from the vote).
        let both_list =
            IndexList::build(vec!["ACGTAAGGCC".into(), revcomp_str("ACGTAAGGCC")], BUDGET).unwrap();
        let both = b"ACGTAAGGCC".to_vec();
        assert!(both_list.exact_hit(&both) && both_list.exact_hit(&revcomp(&both)));
        let vboth = vote_orientation(&[both], &both_list);
        assert_eq!(vboth.ambig, 1);

        // 800 rc-only reads -> RevComp.
        let samples: Vec<Vec<u8>> = (0..800).map(|_| rc_read.clone()).collect();
        let v = vote_orientation(&samples, &rc_list);
        assert_eq!(v.rc, 800);
        assert_eq!(v.fwd, 0);
        assert_eq!(resolve_orientation(&v, "i7").unwrap(), Orient::RevComp);
    }

    fn revcomp_str(s: &str) -> String {
        String::from_utf8(revcomp(s.as_bytes())).unwrap()
    }

    #[test]
    fn orient_low_confidence_errors() {
        let list = IndexList::build(vec!["ACGTAAGGCC".into()], BUDGET).unwrap();
        // 50 unmatched reads -> discrim = 0 -> low confidence error
        let samples: Vec<Vec<u8>> = (0..50).map(|_| b"TTTTTTTTTT".to_vec()).collect();
        let v = vote_orientation(&samples, &list);
        assert_eq!(v.unmatched, 50);
        assert!(resolve_orientation(&v, "i7").is_err());
    }

    #[test]
    fn orient_mixed_errors() {
        let v = Vote { fwd: 500, rc: 500, ambig: 0, unmatched: 0 };
        assert!(resolve_orientation(&v, "i7").is_err());
    }

    #[test]
    fn orient_near_even_tie_errors() {
        // discrim=200 (>= MIN_DISCRIM) but an exact 100/100 split: must error "mixed", never
        // silently resolve to Forward via the old minor<MIN_DISCRIM hole.
        let v = Vote { fwd: 100, rc: 100, ambig: 0, unmatched: 0 };
        assert!(resolve_orientation(&v, "i7").is_err());
        // Slight majority still under the dominance floor -> still mixed.
        let v2 = Vote { fwd: 120, rc: 100, ambig: 0, unmatched: 0 };
        assert!(resolve_orientation(&v2, "i7").is_err());
    }

    #[test]
    fn orient_palindrome_excluded_resolves_revcomp() {
        // 20% palindromic (ambig) + 80% rc -> RevComp, never mixed.
        let v = Vote { fwd: 0, rc: 800, ambig: 200, unmatched: 0 };
        assert_eq!(resolve_orientation(&v, "i7").unwrap(), Orient::RevComp);
    }

    // ---- CSV parsing ----
    #[test]
    fn csv_crlf_semicolon() {
        let content = b"CellID;i7_index;i5_index\r\n1;ACGTACGTAC;TTGGTACGCG\r\n2;TACAACCTCA;TTGGTACGCG\r\n";
        let path = tmp("crlf.csv", content);
        let d = probe_table(&path).unwrap();
        assert_eq!(d.i7_len, 10);
        assert_eq!(d.i5_len, 10);
        // row/index counts now come from the full parser (TableDims dropped its count fields).
        let t = parse_table(&path).unwrap();
        assert_eq!(t.i7.len(), 2);
        assert_eq!(t.i5.len(), 1);
        assert_eq!(t.pairs.len(), 2);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn csv_leading_bom() {
        // Excel-style UTF-8 BOM before the header, with i7_index as the FIRST column: the BOM must
        // be stripped so the column lookup still finds "i7_index".
        let content = "\u{feff}i7_index;i5_index\r\nACGTACGTAC;TTGGTACGCG\r\n".as_bytes();
        let path = tmp("bom.csv", content);
        let d = probe_table(&path).unwrap();
        assert_eq!(d.i7_len, 10);
        assert_eq!(parse_table(&path).unwrap().i7.len(), 1);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn csv_three_delimiters() {
        let comma = b"a,i7_index,i5_index\n1,ACGTACGTAC,TTGGTACGCG\n";
        let tab = b"a\ti7_index\ti5_index\n1\tACGTACGTAC\tTTGGTACGCG\n";
        let semi = b"a;i7_index;i5_index\n1;ACGTACGTAC;TTGGTACGCG\n";
        for (name, c) in [("c.csv", &comma[..]), ("t.csv", &tab[..]), ("s.csv", &semi[..])] {
            let p = tmp(name, c);
            let d = probe_table(&p).unwrap();
            assert_eq!(d.i7_len, 10);
            assert_eq!(parse_table(&p).unwrap().i7.len(), 1);
            std::fs::remove_file(&p).ok();
        }
    }

    #[test]
    fn csv_reordered_columns() {
        // i5 before i7
        let content = b"CellID,i5_index,i7_index\n1,TTGGTACGCG,ACGTACGTAC\n";
        let p = tmp("reorder.csv", content);
        let t = parse_table(&p).unwrap();
        assert_eq!(t.i7[0], "ACGTACGTAC");
        assert_eq!(t.i5[0], "TTGGTACGCG");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn csv_variable_length() {
        for (name, l7, l5) in [("l8", 8, 8), ("l10", 10, 10), ("l12", 12, 12)] {
            let seven: String = "A".repeat(l7);
            let five: String = "C".repeat(l5);
            let content = format!("i7_index,i5_index\n{seven},{five}\n");
            let p = tmp(name, content.as_bytes());
            let d = probe_table(&p).unwrap();
            assert_eq!(d.i7_len, l7);
            assert_eq!(d.i5_len, l5);
            std::fs::remove_file(&p).ok();
        }
    }

    #[test]
    fn csv_non_acgt_errors() {
        let content = b"i7_index,i5_index\nACGTNCGTAC,TTGGTACGCG\n";
        let p = tmp("nonacgt.csv", content);
        assert!(probe_table(&p).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn csv_missing_column_errors() {
        let content = b"CellID,i5_index\n1,TTGGTACGCG\n";
        let p = tmp("missingcol.csv", content);
        assert!(probe_table(&p).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn csv_length_mismatch_errors() {
        let content = b"i7_index,i5_index\nACGTACGTAC,TTGGTACGCG\nACGTACG,TTGGTACGCG\n";
        let p = tmp("lenmismatch.csv", content);
        assert!(probe_table(&p).is_err());
        std::fs::remove_file(&p).ok();
    }

    // ---- pair validity + full IndexTable ----
    #[test]
    fn pair_validity_and_label() {
        // Two i7 x two i5, but only 3 of the 4 pairs present.
        let content = b"i7_index,i5_index\n\
            AAAAAAAAAA,CCCCCCCCCC\n\
            AAAAAAAAAA,GGGGGGGGGG\n\
            TTTTTTTTTT,CCCCCCCCCC\n";
        let p = tmp("pairs.csv", content);
        let table = IndexTable::load(&p, BUDGET).unwrap();
        // valid pair -> label = i7 ++ i5
        let (label, corrected) =
            table.assign_pair(b"AAAAAAAAAA", b"CCCCCCCCCC").unwrap();
        assert_eq!(label, "AAAAAAAAAACCCCCCCCCC");
        assert!(!corrected);
        // each index individually valid but the PAIR (TTTTTTTTTT, GGGGGGGGGG) is not a row
        assert!(table.assign_pair(b"TTTTTTTTTT", b"GGGGGGGGGG").is_none());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn assign_pair_corrected_flag() {
        let content = b"i7_index,i5_index\nACGTACGTAC,TGCATGCATG\n";
        let p = tmp("corr.csv", content);
        let table = IndexTable::load(&p, BUDGET).unwrap();
        // one substitution on i7 -> corrected = true
        let (label, corrected) = table.assign_pair(b"ACGTACGTAA", b"TGCATGCATG").unwrap();
        assert_eq!(label, "ACGTACGTACTGCATGCATG");
        assert!(corrected);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn duplicate_index_is_hard_error() {
        assert!(IndexList::build(vec!["ACGTACGTAC".into(), "ACGTACGTAC".into()], BUDGET).is_err());
    }

    // ---- assign_into: low-alloc path + drop-reason split (§1/§6) ----
    #[test]
    fn assign_into_reasons_split() {
        // i7 has two Hamming-distance-2 codes so "AAAAAAAAAC" is a shared 1-sub window -> REJECT.
        let content = b"i7_index,i5_index\n\
            AAAAAAAAAA,CCCCCCCCCC\n\
            AAAAAAAACC,CCCCCCCCCC\n";
        let p = tmp("reasons.csv", content);
        let table = IndexTable::load(&p, BUDGET).unwrap();
        let mut sc = RevcompScratch::new();

        // exact valid pair -> Assigned, not corrected, interned label reachable via label().
        match table.assign_into(b"AAAAAAAAAA", b"CCCCCCCCCC", &mut sc) {
            AssignResult::Assigned { pair_id, label, corrected } => {
                assert_eq!(label, "AAAAAAAAAACCCCCCCCCC");
                assert_eq!(table.label(pair_id), "AAAAAAAAAACCCCCCCCCC");
                assert!(!corrected);
            }
            other => panic!("expected Assigned, got {other:?}"),
        }
        // shared 1-sub i7 window -> ambiguous reject (not a plain miss).
        assert_eq!(
            table.assign_into(b"AAAAAAAAAC", b"CCCCCCCCCC", &mut sc),
            AssignResult::Unassigned(UnassignedReason::AmbiguousReject)
        );
        // i7 far from any listed index -> plain miss / absent.
        assert_eq!(
            table.assign_into(b"TTTTTTTTTT", b"CCCCCCCCCC", &mut sc),
            AssignResult::Unassigned(UnassignedReason::Absent)
        );
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn assign_into_invalid_pair_reason() {
        // both indexes decode, but the (TTTTTTTTTT, GGGGGGGGGG) combination is not a listed row.
        let content = b"i7_index,i5_index\n\
            AAAAAAAAAA,CCCCCCCCCC\n\
            AAAAAAAAAA,GGGGGGGGGG\n\
            TTTTTTTTTT,CCCCCCCCCC\n";
        let p = tmp("invalidpair.csv", content);
        let table = IndexTable::load(&p, BUDGET).unwrap();
        let mut sc = RevcompScratch::new();
        assert_eq!(
            table.assign_into(b"TTTTTTTTTT", b"GGGGGGGGGG", &mut sc),
            AssignResult::Unassigned(UnassignedReason::InvalidPair)
        );
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn assign_into_revcomp_scratch_path() {
        // Orientation RevComp: the raw read is the revcomp of the listed code; scratch is reused
        // across calls (no per-read revcomp Vec). One sub on i7 -> corrected = true.
        let content = b"i7_index,i5_index\nACGTACGTAC,TGCATGCATG\n";
        let p = tmp("rc.csv", content);
        let mut table = IndexTable::load(&p, BUDGET).unwrap();
        table.i7_orient = Orient::RevComp;
        table.i5_orient = Orient::RevComp;
        let raw_i7 = revcomp(b"ACGTACGTAC");
        let raw_i5 = revcomp(b"TGCATGCATG");
        let mut sc = RevcompScratch::new();
        match table.assign_into(&raw_i7, &raw_i5, &mut sc) {
            AssignResult::Assigned { label, corrected, .. } => {
                assert_eq!(label, "ACGTACGTACTGCATGCATG");
                assert!(!corrected);
            }
            other => panic!("expected Assigned, got {other:?}"),
        }
        // reuse the SAME scratch on a second, sphere-corrected call (one sub on the fwd i7 window).
        let raw_i7_1sub = revcomp(b"ACGTACGTAA");
        match table.assign_into(&raw_i7_1sub, &raw_i5, &mut sc) {
            AssignResult::Assigned { label, corrected, .. } => {
                assert_eq!(label, "ACGTACGTACTGCATGCATG");
                assert!(corrected);
            }
            other => panic!("expected Assigned, got {other:?}"),
        }
        std::fs::remove_file(&p).ok();
    }
}
