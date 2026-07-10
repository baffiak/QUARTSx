# QUARTSx

**QUantification And Reporting Tool for Smart-seq3(x)press.** A from-scratch
reimplementation of [zUMIs](https://github.com/sdparekh/zUMIs), specialised to
SmartSeq3xpress — output-compatible but not numerically identical, and built as a single-pass,
low-disk pipeline. One Rust binary runs the whole thing
(FILTERING → MAPPING → COUNTING → SUMMARISING); STAR does the alignment, everything else is
native Rust in a single one-pass, low-disk reader. R is used only for the final QC report.

## What it is

QUARTSx reprocesses SmartSeq3xpress FASTQs the way zUMIs does — same gene-assignment model,
same output structure — but rewrites the three steps that make zUMIs slow and disk-hungry on this
chemistry: read filtering, cell-barcode demultiplexing, and molecule counting. The design
goal is that a QUARTSx run stays directly comparable to a published zUMIs analysis (Parekh
*et al.* 2018), while the I/O and error-correction steps are reworked for speed and lower disk use.

Concretely, QUARTSx keeps from zUMIs:

- **STAR** for alignment, unchanged (pinned to `star=2.7.10b` in `env.yaml`; Dobin *et al.* 2013).
  The barcode/UMI tags ride through STAR untouched.
- The **exon / intron / inex** feature model and largest-overlap gene assignment, with a
  fraction-overlap gate, `multi_overlap` sharing, and the exon-preference
  inex union (a read counts toward its exon gene if it has any exonic overlap, otherwise toward
  its intron gene).
- The per-cell read taxonomy (Exon / Intron / Intergenic / Ambiguity / Unmapped / User),
  binomial downsampling, and the nested `dgecounts.rds` layout
  (`umicount` / `readcount` / `readcount_internal` × `exon` / `inex` / `intron`).

The fixed SS3xpress read layout is assumed throughout: **R1** = TSO tag + UMI + 5′ cDNA,
**R2** = internal cDNA, **I1 + I2** = dual sample index → cell barcode. Reads whose
R1 matches the TSO tag (fixed-position, within a mismatch budget) are *tagged* and carry an
UMI; the rest are *internal* (untagged, tiled across the transcript by tagmentation).

## What is different

Two things are deliberately different from zUMIs — UMI deduplication and barcode demultiplexing. Everything else (gene assignment, intron/exon
counting, the count layers) reproduces zUMIs.

### Molecule counting: UMI-tools directional deduplication

UMIs are deduplicated per (cell, gene-set) with the **UMI-tools *directional*** network method
(Smith, Heger & Sudbery 2017). Each distinct UMI is a node carrying its read count; a directed
edge A → B is drawn when the two UMIs differ by exactly one base **and** `count(A) ≥ 2·count(B) − 1`. Nodes
are visited from highest count down; each unvisited node seeds a network that absorbs everything
reachable along its out-edges, and the number of networks is the molecule count.

### Cell-barcode demultiplexing: indel-aware Sequence-Levenshtein

The dual sample index (I1 = i7, I2 = i5) is decoded against the index table with an **indel-aware**
scheme instead of zUMIs' fixed-whitelist Hamming binning. Each index is corrected within a budget
of **≤ 2 total edits, of which ≤ 1 is an indel**, using the **Sequence-Levenshtein** distance
(Buschmann & Bystrykh 2013) realised as a precomputed FREE-style edit-sphere lookup (Hawkins
*et al.* 2018 — only their FREE edit-sphere *decode* model is used).

At build time every reachable read-window of every index is enumerated into an
exact-then-sphere hash map, so decoding is an O(1) lookup. A window reachable from two different
indices is marked with a reject sentinel, so an ambiguous read is dropped, never misassigned.

A per-index guard warns when the configured budget exceeds the panel's safe correction radius
(⌊(min pairwise distance − 1)/2⌋), and each index's orientation (forward vs reverse-complement)
is auto-detected from the first reads, because i5 is read in either orientation depending on the
instrument.

Note this is a structurally different barcode model from zUMIs (dual i7/i5 index table, not a
single whitelist).

### Counting model: reproduces zUMIs

SmartSeq3xpress yields two read populations:

- **Tagged** reads mark the molecule's 5′ end and carry a UMI → deduplicated to UMI molecule counts
  (`umicount`, directional dedup). With `strand: 1` these reads must match the gene's strand.
- **Internal** reads tile the transcript body and have no UMI → counted as reads. They populate
  `readcount_internal` (per cell × gene `readcount_internal`) and, for the exon/all
  layer, a length- and library-size-normalised `rpkm`. Internal reads are unstranded.

Gene assignment, introns and layers all follow zUMIs: a read goes to the gene it overlaps most
(featureCounts `largestOverlap`, Liao *et al.* 2014 — the mode zUMIs configures), introns are the
single-gene gaps of the chromosome-wide exon union, and counts are reported over
exon / intron / exon∪intron (`inex`).

## Install

Recommended: create the environment with conda from `env.yaml` (channels conda-forge and bioconda),
then activate it and work inside it.

```bash
git clone https://github.com/baffiak/QUARTSx.git
cd QUARTSx
conda env create -n quartsx -f env.yaml
conda activate quartsx
```

`env.yaml` records the exact combination QUARTSx was tested with. Treat it as a worked example of a
setup that runs, not a promise that future package versions will keep working — bioinformatics tools
change behaviour between releases.

QUARTSx is compiled from source: `quartsx.sh` runs `cargo build --release` once on the first run and
reuses the binary afterwards. Building needs the Rust toolchain, which `env.yaml` provides.

### Dependencies (installing without conda)

conda only resolves dependencies if you use conda. Install by hand instead and you install and
resolve every package yourself. Below is the short list of what the pipeline calls directly, with
the versions it was tested with, other versions may or may not work:

| Package | Tested with |
| --- | --- |
| rust | 1.96.1 |
| star | 2.7.10b |
| samtools | 1.23.1 |
| falco | 1.3.2 |
| r-base | 4.5.3 |
| r-matrix | 1.7-5 |
| r-data.table | 1.17.8 |
| r-jsonlite | 2.0.0 |
| r-ggplot2 | 4.0.3 |
| r-svglite | 2.2.2 |
| r-yaml | 2.3.12 |

Use STAR 2.7.10b: the bundled `testdata/` index was built with it.

The table lists only the packages the pipeline calls directly. The complete resolved environment —
every package and exact version of a validated install — is in `env.tested.txt` (a `conda list`
dump).

## Test

A small self-contained dataset ships in `testdata/`: a mini reference (6 multi-exon genes on both
strands, with a prebuilt STAR index) and synthetic R1/R2/I1/I2 FASTQs covering exact/substitution/
indel barcodes, tagged and internal reads, and dedup groups. `test_config.yaml` in the repo root
points at it with repo-relative paths, so run it from the repo root:

```bash
./quartsx.sh -y test_config.yaml
```

Outputs land in `testdata/run/`.

## Usage

```bash
./quartsx.sh -y config.yaml
```

`-y` is the config path. Activate your environment first so `cargo`, STAR, samtools, fastqc and
`Rscript` are on `PATH`; everything runs there. All configuration is in the YAML — copy
`test_config.yaml` and edit. The fields:

- **`project`** — output basename.
- **`out_dir`** — flat output directory.
- **`star_tmp`** — parent for STAR's temp dir; keep it on a native filesystem, not under `out_dir`.
- **`num_threads`**, **`mem_limit`** (GB) — resources.

**`sequence_files`** — the four SS3xpress FASTQs and their layout:
- `file1` (R1): `base_definition` `cDNA(25-150)` + `UMI(12-19)`, and `find_pattern:
  "ATTGCGCAATG;1"` — the TSO tag and its mismatch budget (`;N`) that splits tagged from internal
  reads.
- `file2` (R2): internal `cDNA(1-150)` — the target of the read-filtering step.
- `file3` / `file4` (I1 / I2): `BC(1-10)` — the i7 and i5 indices; lengths must match the
  `i7_index` / `i5_index` columns of the index table.

**`reference`** —
- `STAR_index`,
- `GTF_file`,
- optional `additional_files` (extra FASTA such as spike-ins; one User gene per contig) and
- `additional_STAR_params`.

**`read_filtering`** — the native R2 pre-filter:
- `adapter_fasta` (shipped SS3x adapters; clipped from R2, `null` to skip),
- `quality` (phred for both-ends + 4 bp sliding-window trim),
- `min_length` (drop if trimmed R2 is shorter).

**`filter_cutoffs`** — drop a read if too many

- barcode (`BC_filter`) or
- UMI (`UMI_filter`)

bases fall below `phred` (`num_bases` = how many low-quality bases are tolerated).

**`barcodes`** —
- `index_table` (CSV/TSV with `i7_index` and `i5_index` columns, delimiter auto-detected),
- `BarcodeBinning` (max total edits per index, ≤ 2, of which ≤ 1 indel; Sequence-Levenshtein),
- `nReadsperCell` (minimum filtered reads for a cell to be kept).

**`counting_opts`**:
- `introns` — build the intron/inex layers
- `strand` — `1` = tagged reads stranded, internal unstranded; `0` = unstranded
- `downsampling` — `0` = none; a number; `"5000,10000"`; or a `"5000-10000"` range
- `multi_overlap` — share a largest-overlap tie across the tied genes, value split 1/n
- `fraction_overlap` — minimum overlap/mapped-length per gene
- `multi_mappers` — how reads mapping to multiple genes are counted: `Unique` (default), `Uniform`, `Rescue`, `PropUnique`, `EM`. With `Unique` only uniquely-mapped reads are counted. The resolver modes (`Uniform`, `Rescue`, `PropUnique`, `EM`) distribute multi-gene reads.
- `multimapper_cap` — (default `20`) caps how many loci a multimapping read may map to under the resolver modes; it is ignored when `multi_mappers` is `Unique` (which drops multimappers).

**`start_stage`** — where a rerun begins, using the kept `filtered/` and `star/` intermediates:

| Stage | Does |
| --- | --- |
| `Filtering` | one-pass read of R1/R2/I1/I2: adapter clip + quality trim + length filter on R2, barcode/UMI quality filter, index demux, write per-cell uBAM shards |
| `Mapping` | STAR alignment of the filtered reads, BC/UB tags carried through |
| `Counting` | gene assignment, exon/intron/inex counting, directional UMI dedup, internal FPKM, downsampling |
| `Summarising` | write `dgecounts.rds`, per-cell stats, and the QC report |

## Outputs (in `out_dir`)

- `<project>.bam` (+ `.bai`) — final coordinate-sorted BAM; every read tagged `BC`/`UB`/`GE`/`GI`.
- `filtered/` — uBAM shards + `filter_stats.json` (resume from Mapping).
- `star/` — STAR outputs, logs and merged BAM (resume from Counting).
- `expression/` — `<project>.dgecounts.rds` (nested count list), `<project>.gene_names.txt`,
  `<project>.readspercell.txt`, `<project>.BCUMIstats.txt`.
- `qc/` — fastqc report on the filtered reads + the features QC report.

## References

Main resources used to build this and thanks to them:
- Parekh S, Ziegenhain C, Vieth B, Enard W, Hellmann I. **zUMIs — A fast and flexible pipeline to
  process RNA sequencing data with UMIs.** *GigaScience* 2018;7(6):giy059.
  https://doi.org/10.1093/gigascience/giy059
- Andrews S. **FastQC: a quality control tool for high throughput sequence data.** Babraham
  Bioinformatics, 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Dobin A, Davis CA, Schlesinger F, *et al.* **STAR: ultrafast universal RNA-seq aligner.**
  *Bioinformatics* 2013;29(1):15–21. https://doi.org/10.1093/bioinformatics/bts635
- Kaminow B, Yunusov D, Dobin A. **STARsolo: accurate, fast and versatile mapping/quantification
  of single-cell and single-nucleus RNA-seq data.** *bioRxiv* 2021.
  https://doi.org/10.1101/2021.05.05.442755

Additional resources used as inspiration.
- Smith T, Heger A, Sudbery I. **UMI-tools: modeling sequencing errors in Unique Molecular
  Identifiers to improve quantification accuracy.** *Genome Research* 2017;27(3):491–499.
  https://doi.org/10.1101/gr.209601.116
- Buschmann T, Bystrykh LV. **Levenshtein error-correcting barcodes for multiplexed DNA
  sequencing.** *BMC Bioinformatics* 2013;14:272. https://doi.org/10.1186/1471-2105-14-272
- Hawkins JA, Jones SK Jr, Finkelstein IJ, Press WH. **Indel-correcting DNA barcodes for
  high-throughput sequencing.** *PNAS* 2018;115(27):E6217–E6226.
  https://doi.org/10.1073/pnas.1802640115
- Hagemann-Jensen M, Ziegenhain C, Chen P, *et al.* **Single-cell RNA counting at allele and
  isoform resolution using Smart-seq3.** *Nature Biotechnology* 2020;38(6):708–714.
  https://doi.org/10.1038/s41587-020-0497-0
- Hagemann-Jensen M, Ziegenhain C, Sandberg R. **Scalable single-cell RNA sequencing from full
  transcripts with Smart-seq3xpress.** *Nature Biotechnology* 2022;40(10):1452–1457.
  https://doi.org/10.1038/s41587-022-01311-4
- Liao Y, Smyth GK, Shi W. **featureCounts: an efficient general purpose program for assigning
  sequence reads to genomic features.** *Bioinformatics* 2014;30(7):923–930.
  https://doi.org/10.1093/bioinformatics/btt656
