use anyhow::{bail, Context, Result};
use serde::de::Deserializer;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Deserialize, Serialize, Clone)]
pub struct Config {
    pub project: String,
    pub sequence_files: SequenceFiles,
    pub reference: Reference,
    pub out_dir: String,
    pub star_tmp: String,
    pub num_threads: usize,
    pub mem_limit: usize,
    // Thread knobs. T = num_threads; the FILTERING pipeline overlaps N filter/encode workers with
    // P BGZF deflaters, so the two concurrent pools must satisfy N + P <= T to avoid oversubscription.
    // Both default (None) => P = T/2, N = T - P. See resolved_threads().
    #[serde(default)]
    pub compress_threads: Option<usize>, // P: BGZF deflate workers (default T/2)
    #[serde(default)]
    pub filter_threads: Option<usize>, // N: filter/encode workers (default T - P)
    pub read_filtering: ReadFiltering,
    pub filter_cutoffs: FilterCutoffs,
    pub barcodes: Barcodes,
    pub counting_opts: CountingOpts,
    #[serde(default)]
    pub start_stage: StartStage,
    #[serde(skip)]
    pub config_path: String,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct SequenceFiles {
    pub file1: SeqFile,
    pub file2: SeqFile,
    pub file3: SeqFile,
    pub file4: SeqFile,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct SeqFile {
    pub name: String,
    pub base_definition: Vec<String>,
    #[serde(default)]
    pub find_pattern: Option<String>,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct Reference {
    #[serde(rename = "STAR_index")]
    pub star_index: String,
    #[serde(rename = "GTF_file")]
    pub gtf_file: String,
    #[serde(default, rename = "additional_files")]
    pub additional_files: Option<String>,
    #[serde(default, rename = "additional_STAR_params")]
    pub additional_star_params: String,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct ReadFiltering {
    #[serde(default)]
    pub adapter_fasta: Option<String>,
    pub quality: u8,
    pub min_length: usize,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct FilterCutoffs {
    #[serde(rename = "BC_filter")]
    pub bc_filter: BaseQual,
    #[serde(rename = "UMI_filter")]
    pub umi_filter: BaseQual,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct BaseQual {
    pub num_bases: usize,
    pub phred: u8,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct Barcodes {
    pub index_table: String,
    #[serde(rename = "BarcodeBinning")]
    pub barcode_binning: u8,
    #[serde(rename = "nReadsperCell")]
    pub n_reads_per_cell: u64,
}

#[derive(Deserialize, Serialize, Clone)]
pub struct CountingOpts {
    pub introns: bool,
    pub strand: u8,
    #[serde(deserialize_with = "de_downsampling")]
    pub downsampling: String,
    pub multi_overlap: bool,
    pub fraction_overlap: f64,
    // Multimapper handling. Mode selects the STARsolo --soloMultiMappers resolution strategy.
    #[serde(default)]
    pub multi_mappers: MultiMapperMode,
    // Multimapper cap for the resolver modes (Uniform/Rescue/PropUnique/EM). Ignored when
    // multi_mappers is Unique (which hardcodes 1/1). Defaults to 20.
    #[serde(default = "default_multimap_nmax")]
    pub multimapper_cap: u32,
}

fn default_multimap_nmax() -> u32 {
    20
}

/// STARsolo `--soloMultiMappers` multi-gene resolution mode. Serde parses the bare mode name
/// straight from YAML (e.g. `multi_mappers: EM`). The five modes are Unique, Uniform, Rescue,
/// PropUnique, EM; distribution formulas live in count.rs.
#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Eq, Debug)]
pub enum MultiMapperMode {
    Unique,
    Uniform,
    Rescue,
    PropUnique,
    EM,
}

impl Default for MultiMapperMode {
    fn default() -> Self {
        MultiMapperMode::Unique
    }
}

#[derive(Deserialize, Serialize, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum StartStage {
    Filtering,
    Mapping,
    Counting,
    Summarising,
}

impl Default for StartStage {
    fn default() -> Self {
        StartStage::Filtering
    }
}

// downsampling accepts a bare int (0, 10000) or a string ("0", "5000,10000", "5000-10000")
fn de_downsampling<'de, D: Deserializer<'de>>(d: D) -> Result<String, D::Error> {
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum IntOrStr {
        Int(i64),
        Str(String),
    }
    Ok(match IntOrStr::deserialize(d)? {
        IntOrStr::Int(i) => i.to_string(),
        IntOrStr::Str(s) => s,
    })
}

// One barcode/UMI/cDNA span from a base_definition, 0-based [start, end).
#[derive(Clone, Copy)]
pub struct Segment {
    pub file: usize,
    pub start: usize,
    pub end: usize,
}

// Resolved SS3xpress read layout: where the cell barcode, UMI, TSO tag and cDNA live.
pub struct Geometry {
    pub files: Vec<String>, // file1..file4 paths, in reader order
    pub bc: Vec<Segment>,   // concatenated (I1 then I2) to form the cell barcode
    pub umi: Segment,       // on the tagged read
    pub tagged_file: usize, // read carrying the TSO tag + UMI (R1)
    pub cdna_start: usize,  // 0-based cDNA start on the tagged read
    pub internal_file: usize, // the internal cDNA read that gets R2 filtering (R2)
    pub tag: Vec<u8>,       // find_pattern sequence
    pub tag_mismatch: u8,   // find_pattern mismatch budget
}

fn parse_base(def: &str, file: usize) -> Result<(String, Segment)> {
    let open = def.find('(').with_context(|| format!("base_definition missing '(': {def}"))?;
    let close = def.find(')').with_context(|| format!("base_definition missing ')': {def}"))?;
    let kind = def[..open].trim().to_string();
    let (a, b) = def[open + 1..close]
        .split_once('-')
        .with_context(|| format!("base_definition range needs a-b: {def}"))?;
    let start: usize = a.trim().parse().with_context(|| format!("bad range start: {def}"))?;
    let end: usize = b.trim().parse().with_context(|| format!("bad range end: {def}"))?;
    Ok((kind, Segment { file, start: start - 1, end })) // 1-based inclusive -> 0-based [start, end)
}

impl Config {
    pub fn geometry(&self) -> Result<Geometry> {
        let sf = &self.sequence_files;
        let files = [&sf.file1, &sf.file2, &sf.file3, &sf.file4];

        let mut bc = Vec::new();
        let mut umi = None;
        let mut tagged_file = None;
        let mut cdna_start = None;
        let mut internal_file = None;

        for (i, f) in files.iter().enumerate() {
            let mut has_umi = false;
            let mut has_cdna = false;
            let mut cstart = 0usize;
            for def in &f.base_definition {
                let (kind, seg) = parse_base(def, i)?;
                match kind.as_str() {
                    "BC" => bc.push(seg),
                    "UMI" => {
                        umi = Some(seg);
                        has_umi = true;
                    }
                    "cDNA" => {
                        has_cdna = true;
                        cstart = seg.start;
                    }
                    other => bail!("unknown base_definition kind '{other}' in {def}"),
                }
            }
            if f.find_pattern.is_some() || has_umi {
                tagged_file = Some(i);
                if has_cdna {
                    cdna_start = Some(cstart);
                }
            } else if has_cdna {
                internal_file = Some(i);
            }
        }

        let tagged_file = tagged_file.context("no read defines a UMI / find_pattern (tagged read)")?;
        let internal_file = internal_file.context("no internal cDNA read (R2)")?;
        let umi = umi.context("no UMI base_definition")?;
        let cdna_start = cdna_start.context("tagged read has no cDNA base_definition")?;
        if bc.is_empty() {
            bail!("no BC base_definition on any read");
        }
        bc.sort_by_key(|s| (s.file, s.start));

        let pat = files[tagged_file]
            .find_pattern
            .as_deref()
            .context("tagged read needs a find_pattern")?;
        let (seq, mm) = pat.split_once(';').unwrap_or((pat, "0"));
        let tag_mismatch = mm.trim().parse().with_context(|| format!("bad find_pattern mismatch: {pat}"))?;

        Ok(Geometry {
            files: files.iter().map(|f| f.name.clone()).collect(),
            bc,
            umi,
            tagged_file,
            cdna_start,
            internal_file,
            tag: seq.trim().as_bytes().to_vec(),
            tag_mismatch,
        })
    }

    /// Resolve the FILTERING thread split into `(P, N)` = (BGZF deflaters, filter/encode workers).
    /// Defaults: `P = compress_threads.unwrap_or(T/2)`, `N = filter_threads.unwrap_or(T-P)`.
    /// Both pools run CONCURRENTLY, so the invariant `N + P <= T` is enforced here (never oversubscribe
    /// the allocated cores); each pool is kept >= 1. If an explicit split would exceed T, P is trimmed
    /// first (BGZF deflate overlaps its ordered write thread and tolerates fewer cores better than the
    /// CPU-bound filter workers).
    ///
    /// T == 1 has NO valid concurrent split (any two positive pools already exceed the one core), so it
    /// is the SINGLE-THREADED case: `(1, 1)` is returned only as the degenerate budget, and callers MUST
    /// take their single-threaded path (gated on `num_threads <= 1`) — running inline with no separate
    /// rayon pool and a single BGZF codec worker — instead of overlapping two CPU-bound pools on the
    /// lone core.
    pub fn resolved_threads(&self) -> (usize, usize) {
        let t = self.num_threads.max(1);
        if t == 1 {
            return (1, 1);
        }
        let mut p = self.compress_threads.unwrap_or(t / 2).max(1);
        let mut n = self.filter_threads.unwrap_or(t.saturating_sub(p)).max(1);
        if p + n > t {
            // Keep at least one core for the filter workers, then give N the remainder.
            p = p.min(t.saturating_sub(1)).max(1);
            n = t.saturating_sub(p).max(1);
        }
        (p, n)
    }
}

pub fn load(path: &str) -> Result<Config> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading config {path}"))?;
    let mut cfg: Config = serde_yaml::from_str(&text).context("parsing config yaml")?;
    cfg.config_path = path.to_string();
    // Expand `~` / `$HOME` in ALL path fields BEFORE validate() so existence checks,
    // FIFO probing and every downstream consumer (orchestrator, filter, count) see absolute paths.
    expand_paths(&mut cfg);
    validate(&cfg)?;
    Ok(cfg)
}

/// Expand a leading `~` / `~/` and `$HOME` / `${HOME}` in a single path string. Non-path YAML
/// (e.g. STAR params) is not passed through here. Empty strings and paths with no marker pass through
/// unchanged. Shell-only forms (`~user`, arbitrary `$VARS`) are intentionally NOT expanded — only the
/// two documented markers, so behavior is identical across macOS / WSL2 / Linux.
fn expand_path(s: &str) -> String {
    if s.is_empty() {
        return String::new();
    }
    let home = match std::env::var("HOME") {
        Ok(h) if !h.is_empty() => h,
        _ => return s.to_string(), // no HOME => leave path untouched rather than corrupt it
    };
    let home_trimmed = home.trim_end_matches('/');
    let mut out = if s == "~" {
        home.clone()
    } else if let Some(rest) = s.strip_prefix("~/") {
        format!("{home_trimmed}/{rest}")
    } else {
        s.to_string()
    };
    // Expand $HOME / ${HOME} anywhere (braced form first — "$HOME" is not a substring of "${HOME}").
    out = out.replace("${HOME}", &home).replace("$HOME", &home);
    out
}

// Apply expand_path to every path-valued field. Kept in one place so the set of expanded fields
// is auditable; orchestrator/filter/count consume the already-expanded paths and never re-expand.
fn expand_opt(field: &mut Option<String>) {
    if let Some(s) = field {
        *s = expand_path(s);
    }
}

fn expand_paths(cfg: &mut Config) {
    cfg.out_dir = expand_path(&cfg.out_dir);
    cfg.star_tmp = expand_path(&cfg.star_tmp);
    cfg.sequence_files.file1.name = expand_path(&cfg.sequence_files.file1.name);
    cfg.sequence_files.file2.name = expand_path(&cfg.sequence_files.file2.name);
    cfg.sequence_files.file3.name = expand_path(&cfg.sequence_files.file3.name);
    cfg.sequence_files.file4.name = expand_path(&cfg.sequence_files.file4.name);
    cfg.reference.star_index = expand_path(&cfg.reference.star_index);
    cfg.reference.gtf_file = expand_path(&cfg.reference.gtf_file);
    expand_opt(&mut cfg.reference.additional_files);
    cfg.barcodes.index_table = expand_path(&cfg.barcodes.index_table);
    expand_opt(&mut cfg.read_filtering.adapter_fasta);
}

fn validate(cfg: &Config) -> Result<()> {
    let sf = &cfg.sequence_files;
    for (name, p) in [
        ("sequence_files.file1", &sf.file1.name),
        ("sequence_files.file2", &sf.file2.name),
        ("sequence_files.file3", &sf.file3.name),
        ("sequence_files.file4", &sf.file4.name),
        ("reference.GTF_file", &cfg.reference.gtf_file),
        ("reference.STAR_index", &cfg.reference.star_index),
        ("barcodes.index_table", &cfg.barcodes.index_table),
    ] {
        if !Path::new(p).exists() {
            bail!("{name} path does not exist: {p}");
        }
    }
    if let Some(fa) = &cfg.reference.additional_files {
        if !Path::new(fa).exists() {
            bail!("reference.additional_files does not exist: {fa}");
        }
    }
    if let Some(fa) = &cfg.read_filtering.adapter_fasta {
        if !Path::new(fa).exists() {
            bail!("read_filtering.adapter_fasta does not exist: {fa}");
        }
    }

    let g = cfg.geometry()?; // fails fast on a malformed base_definition / find_pattern

    // Index table: charset/length/columns/delimiter validated in probe_table; cross-check that the
    // yaml BC() slice lengths on I1/I2 match the table's i7/i5 column lengths.
    let dims = crate::barcode::probe_table(&cfg.barcodes.index_table)?;
    if g.bc.len() != 2 {
        bail!("SS3xpress expects exactly two BC segments (I1 i7 + I2 i5), found {}", g.bc.len());
    }
    let i7_seg_len = g.bc[0].end - g.bc[0].start; // I1 / file3
    let i5_seg_len = g.bc[1].end - g.bc[1].start; // I2 / file4
    if i7_seg_len != dims.i7_len {
        bail!("BC slice on I1 is {i7_seg_len} bp but i7_index column is {} bp", dims.i7_len);
    }
    if i5_seg_len != dims.i5_len {
        bail!("BC slice on I2 is {i5_seg_len} bp but i5_index column is {} bp", dims.i5_len);
    }

    // An EXPLICIT split that oversubscribes the cores is a config error, caught here
    // (resolved_threads() also clamps defensively, but an explicit N+P>T is worth a clear message).
    if cfg.num_threads == 0 {
        bail!("num_threads must be >= 1");
    }
    if let (Some(p), Some(n)) = (cfg.compress_threads, cfg.filter_threads) {
        if p == 0 || n == 0 {
            bail!("compress_threads and filter_threads must each be >= 1 (got P={p}, N={n})");
        }
        if p + n > cfg.num_threads {
            bail!(
                "compress_threads (P={p}) + filter_threads (N={n}) = {} exceeds num_threads (T={}); \
                 the filter/encode and BGZF-deflate pools run concurrently and must satisfy N+P<=T",
                p + n,
                cfg.num_threads
            );
        }
    }

    let out = Path::new(cfg.out_dir.trim_end_matches('/'));
    if Path::new(&cfg.star_tmp).starts_with(out) {
        bail!(
            "star_tmp ({}) must not live under out_dir ({}); it must be on a native FS",
            cfg.star_tmp, cfg.out_dir
        );
    }
    probe_fifo(&cfg.star_tmp)?;

    Ok(())
}

// STAR streams the shard through a FIFO in star_tmp; overlay/network filesystems reject mkfifo, which
// would otherwise fail deep inside STAR. Probe up front.
fn probe_fifo(star_tmp: &str) -> Result<()> {
    use std::os::unix::ffi::OsStrExt;
    std::fs::create_dir_all(star_tmp).with_context(|| format!("creating star_tmp {star_tmp}"))?;
    let probe = Path::new(star_tmp).join(".quartsx_fifo_probe");
    let _ = std::fs::remove_file(&probe);
    let cpath = std::ffi::CString::new(probe.as_os_str().as_bytes()).context("star_tmp probe path")?;
    let rc = unsafe { libc::mkfifo(cpath.as_ptr(), 0o600) };
    if rc != 0 {
        let err = std::io::Error::last_os_error();
        match err.raw_os_error() {
            Some(libc::EPERM) | Some(libc::ENOSYS) | Some(libc::EOPNOTSUPP) => bail!(
                "star_tmp ({star_tmp}) is on a filesystem that does not support FIFOs (mkfifo failed: {err}).\n\
                 QUARTSx streams the decompressed shard into STAR through a named pipe here, so star_tmp \
                 MUST be on a native POSIX filesystem (ext4 on Linux/WSL2, apfs on macOS).\n\
                 Common cause: under WSL2 a Windows drive mounted via drvfs/9p (e.g. /mnt/c, /mnt/e) \
                 cannot host FIFOs. Move star_tmp onto the Linux filesystem, e.g. /home/<user>/quartsx_tmp \
                 or /tmp, and keep it OFF /mnt/*."
            ),
            _ => bail!("probing star_tmp ({star_tmp}) for FIFO support failed (mkfifo: {err})"),
        }
    }
    let _ = std::fs::remove_file(&probe);
    Ok(())
}
