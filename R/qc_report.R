#!/usr/bin/env Rscript
# Build one self-contained HTML QC report: qc/<project>.qc_report.html.
# Each ggplot is rendered to vector SVG and embedded inline (svg+xml data URI, so it
# stays crisp/zoomable); tables are plain HTML. Laid out to fit a single A4 page.
# No rmarkdown/knitr, no PDF.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(jsonlite)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: qc_report.R <out_dir> <project>")
out_dir <- args[1]
project <- args[2]
qc_dir  <- file.path(out_dir, "qc")

featColors <- c(Exon = "#1A5084", `Intron+Exon` = "#914614", Intron = "#118730",
                Unmapped = "grey33", Ambiguity = "tan1", MultiMapping = "#631879FF",
                Intergenic = "gold1", `Unused BC` = "grey73", User = "firebrick3")
feat_lab <- c(exon = "Exon", intron = "Intron", inex = "Intron+Exon")
cat_levels <- c("Transcriptome", "Spike-in")
rt_levels  <- c("UMI-tagged", "Internal")

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
comma <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")

# ggh4x provides facet_grid2(independent="y"): a true grid (aligned category x read_type
# strips) with a genuinely independent y scale per panel, which base facet_grid cannot do.
# Load it, installing on demand so the report builds unattended on mac/WSL2/Ubuntu; a clear
# error is raised only if both the load and the install fail.
ensure_ggh4x <- function() {
  if (!requireNamespace("ggh4x", quietly = TRUE))
    try(install.packages("ggh4x", repos = "https://cloud.r-project.org"), silent = TRUE)
  if (!requireNamespace("ggh4x", quietly = TRUE))
    stop("Package 'ggh4x' is required for per-facet free-y QC panels but could not be loaded or installed. Install it with install.packages('ggh4x') and re-run.")
}

# render a ggplot to vector SVG, embed as an svg+xml data URI (each SVG is an isolated
# document, so clip-path ids never collide across panels); w/h are inches (aspect ratio)
img <- function(plot, w, h) {
  f <- tempfile(fileext = ".svg")
  svglite::svglite(f, width = w, height = h)
  print(plot)
  invisible(dev.off())
  b64 <- jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
  file.remove(f)
  sprintf('<img src="data:image/svg+xml;base64,%s" alt="plot">', gsub("\n", "", b64, fixed = TRUE))
}

# ONE builder for every horizontal 100%-stacked percent bar in the report (read-type
# composition AND discarded reads), so both are identical by construction: same geom,
# coord_flip, theme, proportions and legend styling. Given a table with a category column
# and a raw-count column, it computes each category's % share and renders the bar. Category
# order: a pre-set factor keeps its (semantic) level order; otherwise segments stack
# smallest -> largest for a tidy bar. Colours use one shared mechanism: scale_fill_manual on
# a named palette. `palette` supplies fixed colours where it names a level (composition ->
# featColors); any level it does not cover is filled from the colourblind-safe Okabe-Ito
# qualitative palette (recycled), so no bar ever needs a bespoke, hand-picked colour set.
pct_stacked_bar <- function(d, category_col, value_col, axis_label, palette = NULL) {
  cat_v <- as.character(d[[category_col]])
  val_v <- as.numeric(d[[value_col]])
  lv <- if (is.factor(d[[category_col]])) levels(d[[category_col]]) else cat_v[order(val_v)]
  b  <- data.table(category = factor(cat_v, levels = lv), pct = 100 * val_v / sum(val_v))
  okabe <- unname(grDevices::palette.colors(palette = "Okabe-Ito"))
  okabe <- okabe[okabe != "#000000"]                       # drop black for readable fills
  pal <- setNames(okabe[((seq_along(lv) - 1L) %% length(okabe)) + 1L], lv)
  if (!is.null(palette)) { have <- intersect(lv, names(palette)); pal[have] <- palette[have] }
  ggplot(b, aes(x = 1, y = pct, fill = category)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = pal, guide = guide_legend(nrow = 1)) +
    labs(x = NULL, y = axis_label) +
    theme_classic(base_size = 9) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), legend.title = element_blank(),
          legend.position = "bottom", legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 7))
}

## data shared by several panels --------------------------------------------
dge     <- readRDS(file.path(out_dir, "expression", paste0(project, ".dgecounts.rds")))
percell <- fread(file.path(out_dir, "expression", paste0(project, ".readspercell.txt")))
config_path <- if (length(args) >= 3) args[3] else ""

kv <- function(k, v) sprintf("<tr><td>%s</td><td>%s</td></tr>", k, esc(as.character(v)))
name_or <- function(x, alt = "none") {
  if (is.null(x) || length(x) == 0) return(alt)
  x <- as.character(x)
  if (is.na(x) || !nzchar(x)) alt else basename(x)
}

## (1) run summary: reference inputs (top), filtered stats (mid), setup (bottom)
# top: reference file names
ref_rows <- kv("Genome (STAR index)", "n/a")
# mid: filter stats — read the scalar metrics TSV from the QC folder (relocated there from
# filtered/, since the filter stats table is QC output). Two columns: metric<TAB>value.
fs_path <- file.path(qc_dir, "filter_stats.tsv")
fs <- list()
if (file.exists(fs_path)) {
  ft <- fread(fs_path, sep = "\t", header = TRUE, colClasses = "character")
  fs <- as.list(setNames(ft$value, ft$metric))
}
num <- function(k) { v <- suppressWarnings(as.numeric(fs[[k]])); if (length(v) == 0 || is.na(v)) 0 else v }

# Friendly labels for drop reasons; unknown `dropped_*` fields are humanised generically, so
# extra sub-reasons emitted by filter.rs appear automatically without changes here.
drop_labels <- c(
  dropped_short_r2                 = "R2 too short",
  dropped_no_barcode               = "No barcode",
  dropped_no_barcode_absent        = "Barcode absent",
  dropped_no_barcode_uncorrectable = "Barcode uncorrectable",
  dropped_no_barcode_ambiguous     = "Barcode ambiguous",
  dropped_bc_quality               = "BC low quality",
  dropped_umi_quality              = "UMI low quality")
humanize_drop <- function(k) {
  if (k %in% names(drop_labels)) return(unname(drop_labels[k]))
  s <- sub("^dropped_", "", k); s <- gsub("_", " ", s)
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}
# dropped reads split by reason (every scalar `dropped_*` metric), count > 0, largest first.
drop_dt <- NULL
if (length(fs)) {
  dk <- names(fs)[grepl("^dropped_", names(fs))]
  dk <- dk[vapply(dk, function(k) !is.na(suppressWarnings(as.numeric(fs[[k]]))), logical(1))]
  drop_dt <- data.table(reason = vapply(dk, humanize_drop, ""), count = vapply(dk, num, 0))
  drop_dt <- drop_dt[count > 0][order(-count)]
}

# (3) reads-per-category stats table: hierarchical parent->child over the ACTUAL filter
# categories, with count, % of total input, and % of the immediate parent's count.
stats_table_html <- ""
if (length(fs)) {
  total <- num("total"); passed <- num("passed"); dropped <- total - passed
  rows <- list(
    list("Total input reads", total,          NA_real_, 0L),
    list("Passed filter",     passed,          total,   1L),
    list("UMI-tagged",        num("tagged"),   passed,  2L),
    list("Internal",          num("internal"), passed,  2L),
    list("Dropped",           dropped,         total,   1L))
  if (!is.null(drop_dt) && nrow(drop_dt))
    for (i in seq_len(nrow(drop_dt)))
      rows <- c(rows, list(list(drop_dt$reason[i], drop_dt$count[i], dropped, 2L)))
  render_row <- function(r) {
    pct_tot <- if (total > 0) sprintf("%.1f%%", 100 * r[[2]] / total) else "&mdash;"
    pct_par <- if (is.na(r[[3]])) "&mdash;" else if (r[[3]] > 0) sprintf("%.1f%%", 100 * r[[2]] / r[[3]]) else "&mdash;"
    sprintf("<tr class='lvl%d'><td class='cat'>%s</td><td class='num'>%s</td><td class='num'>%s</td><td class='num'>%s</td></tr>",
            r[[4]], esc(r[[1]]), comma(r[[2]]), pct_tot, pct_par)
  }
  stats_table_html <- sprintf(
    "<table class='stats hier'><tr><th>Category</th><th class='num'>Count</th><th class='num'>%% total</th><th class='num'>%% parent</th></tr>%s</table>",
    paste(vapply(rows, render_row, ""), collapse = ""))
}
# run setup section (goes at the BOTTOM of the report): command + the actual run
# configuration (identity + IMPORTANT, result-affecting parameters; no rds/downsampling
# noise). Always renders: project comes from args, config-derived fields are added only
# when the config is available (a NULL field is simply dropped, never fabricated).
setup_items <- c(Project = project)
cmd_html <- ""
mm_html  <- ""   # top Run-summary metric (multimapper reads); filled below when a resolver mode is set
if (nzchar(config_path) && file.exists(config_path)) {
  y  <- yaml::read_yaml(config_path)
  rf <- y$read_filtering; fc <- y$filter_cutoffs; bc <- y$barcodes; co <- y$counting_opts
  ref_rows <- paste0(
    kv("Genome (STAR index)", name_or(y$reference$STAR_index)),
    kv("Annotation (GTF)",    name_or(y$reference$GTF_file)),
    kv("Additional files",    name_or(y$reference$additional_files)))
  # quartsx.sh accepts ONLY -y <config>; the previously shown -c flag never existed (spec §6).
  cmd_html <- sprintf("<p class='cmd'><code>quartsx.sh -y %s</code></p>", esc(config_path))
  setup_items <- c(setup_items,
    `Start stage`      = name_or(y$start_stage),
    `Multimapper mode` = co$multi_mappers,
    introns            = co$introns,
    strand             = co$strand,
    `R2 adapter`       = name_or(rf$adapter_fasta),
    `R2 quality`       = rf$quality,
    `R2 min length`    = rf$min_length,
    `BC filter`        = sprintf("%s bp < Q%s", fc$BC_filter$num_bases,  fc$BC_filter$phred),
    `UMI filter`       = sprintf("%s bp < Q%s", fc$UMI_filter$num_bases, fc$UMI_filter$phred),
    `Barcode binning`  = bc$BarcodeBinning,
    `Min reads/cell`   = bc$nReadsperCell)
  # When a multimapper RESOLVER is used (mode != Unique), surface (as a top Run-summary metric)
  # how many reads STAR mapped to multiple loci, read straight from the STAR final log. Omitted
  # entirely for Unique mode.
  mm_mode <- co$multi_mappers
  if (!is.null(mm_mode) && nzchar(as.character(mm_mode)) &&
      !identical(tolower(as.character(mm_mode)), "unique")) {
    star_log <- file.path(out_dir, "star", paste0(project, ".Log.final.out"))
    if (file.exists(star_log)) {
      ll <- readLines(star_log, warn = FALSE)
      grab <- function(pat) {
        hit <- grep(pat, ll, value = TRUE, fixed = TRUE)
        if (!length(hit)) return(NA_character_)
        trimws(sub(".*\\|", "", hit[1]))
      }
      n_multi <- grab("Number of reads mapped to multiple loci")
      p_multi <- grab("% of reads mapped to multiple loci")
      if (!is.na(n_multi) && nzchar(n_multi)) {
        val <- comma(n_multi)
        if (!is.na(p_multi) && nzchar(p_multi)) val <- sprintf("%s (%s of input)", val, p_multi)
        mm_html <- sprintf("<table class='meta'>%s</table>", kv("Multimapper reads", val))
      }
    }
  }
}
dl <- paste(sprintf("<div><dt>%s</dt><dd>%s</dd></div>", names(setup_items), esc(as.character(setup_items))), collapse = "")
setup_html <- sprintf(
  "<section class='panel setup'><h2>Run setup</h2>
     %s
     <dl class='kv'>%s</dl></section>",
  cmd_html, dl)
summary_html <- sprintf(
  "<section class='panel summary'><h2>Run summary</h2>
     <div class='sumcols'>
       <div class='sumcol refcol'><h3>Reference inputs</h3><table class='meta'>%s</table></div>
       <div class='sumcol statcol'><h3>Reads per category</h3>%s%s</div>
     </div></section>",
  ref_rows, stats_table_html, mm_html)

## (2) discarded-reads breakdown --------------------------------------------
# One compact 100%-stacked horizontal bar of the DROPPED reads, split into one segment per
# reason; segment size = share of the dropped reads (% of filtered), read off the axis. Built
# by the SAME pct_stacked_bar() as Read-type composition, so the two bars are identical by
# construction. Source-agnostic: drop_dt already holds every scalar `dropped_*` metric with
# count > 0, so if filter.rs adds sub-reasons they appear as extra segments with no change here.
drops_html <- ""
if (!is.null(drop_dt) && nrow(drop_dt)) {
  dbar <- pct_stacked_bar(drop_dt, "reason", "count", "% of filtered")
  drops_html <- sprintf(
    "<section class='panel drops'><h2>Discarded reads</h2>%s</section>",
    img(dbar, 7.4, 1.15))
}

per_cell <- function(quant, detected) {
  rbindlist(lapply(names(dge[[quant]]), function(f) {
    m <- dge[[quant]][[f]]$all
    v <- if (detected) Matrix::colSums(m >= 1) else Matrix::colSums(m)
    data.table(type = feat_lab[[f]], Count = as.numeric(v))
  }))
}

gc <- per_cell("readcount", TRUE);  gc$metric <- "Genes"
uc <- per_cell("umicount",  FALSE); uc$metric <- "UMIs"
cnt <- rbind(gc, uc)
cnt[, type   := factor(type, levels = c("Exon", "Intron", "Intron+Exon"))]
cnt[, metric := factor(metric, levels = c("Genes", "UMIs"))]
med <- cnt[, .(n = round(median(Count))), by = .(metric, type)]
p_counts <- ggplot(cnt, aes(type, Count, fill = type)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  geom_text(data = med, aes(type, n, label = n), colour = "orange", size = 2.6, vjust = -0.4) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = featColors) +
  labs(x = NULL, y = "Per cell") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

bar_levels <- rev(c("Exon", "Intron", "Intergenic", "Ambiguity", "MultiMapping", "Unmapped", "User", "Unused BC"))
sumBar <- percell[RG != "bad", .(tot = sum(N)), by = type]
sumBar <- rbind(sumBar, data.table(type = "Unused BC", tot = percell[RG == "bad", sum(N)]))
sumBar[, type := factor(type, levels = bar_levels)]
bar <- pct_stacked_bar(sumBar, "type", "tot", "% of total reads", featColors)

# per-cell % of reads in each feature type (distribution across kept cells)
box_levels <- c("Exon", "Intron", "Intergenic", "Ambiguity", "Unmapped", "User")
pc <- percell[RG != "bad"]
pc[, perc := 100 * N / sum(N), by = RG]
pc[, type := factor(type, levels = box_levels)]
pctbox <- ggplot(pc, aes(type, perc, fill = type)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_manual(values = featColors) +
  labs(x = NULL, y = "% reads/cell") +
  theme_bw(base_size = 9) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

counts_html <- sprintf(
  "<section class='panel counts'><h2>Per-cell counts</h2>%s</section>", img(p_counts, 3.6, 1.7))
compos_html <- sprintf(
  "<section class='panel compos'><h2>Read-type composition</h2>%s</section>", img(bar, 7.4, 1.15))
pctbox_html <- sprintf(
  "<section class='panel pctbox'><h2>Per-cell read-type %%</h2>%s</section>", img(pctbox, 3.6, 1.7))

## helper: 2x2 facet by category (rows) x read_type (cols), used by both panels
facet_cats <- function(d) intersect(cat_levels, unique(d$category))

## (3) gene-body 5'->3' coverage ---------------------------------------------
# The upstream .genebody_coverage.txt holds per-percentile NORMALISED coverage (0-100 gene-body
# percentile, 5'->3'); here we only plot it. The y-axis is FREED per facet via
# ggh4x::facet_grid2(independent="y"): SS3xpress UMI-tagged vs internal reads sit in very different
# normalised ranges (~0.02-0.72), so a shared 0..1 axis flattened the real shape into a near-flat
# line. Per-panel free y lets each facet scale to its own data (spec §6).
cover_html <- ""
gb_path <- file.path(qc_dir, paste0(project, ".genebody_coverage.txt"))
if (file.exists(gb_path)) {
  gb <- fread(gb_path)
  if (nrow(gb)) {
    rts  <- intersect(rt_levels, unique(gb$read_type))
    # drive the category factor from the levels PRESENT in the data (droplevels), so a
    # facet ROW is only rendered for categories that actually occur (no empty Spike-in row
    # when there are no spike/additional genes) rather than a hardcoded 2-level factor.
    gb[, category  := droplevels(factor(category, levels = facet_cats(gb)))]
    gb[, read_type := factor(read_type, levels = rts)]
    cats <- levels(gb$category)
    ensure_ggh4x()
    # rows=category, cols=read_type; independent="y" gives each panel its OWN y range
    # (base facet_grid would share y across a whole row of read_types).
    gcov <- ggplot(gb, aes(percentile, coverage)) +
      geom_line(colour = "#118730", linewidth = 0.6) +
      ggh4x::facet_grid2(vars(category), vars(read_type),
                         scales = "free_y", independent = "y") +
      labs(x = "Gene body percentile (5'->3')", y = "Norm. coverage") +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.5))
    cover_html <- sprintf(
      "<section class='panel cover'><h2>Gene-body coverage (5'&rarr;3')</h2>%s</section>",
      img(gcov, 7.4, 0.7 * length(cats) + 0.7))
  }
}

## (4) insert-size distribution, same category x read_type split as coverage --
# The upstream insert size is computed in TRANSCRIPT space (exon-union projection, spec §5), so
# introns are already collapsed and only genuine fragment length remains; discordant/edge cases can
# still fall outside the plotted window. Plot (spec §5): window to a plausible 100 bp-10 kb range
# (drop <100 and >10 kb BEFORE plotting), scale_x_log10() so the ~250 bp mode is legible instead of
# squashed against a long tail, geom_area/line instead of geom_col over thousands of 1 bp bins, and
# y FREED per facet via ggh4x::facet_grid2(independent="y") (a shared y would crush Transcriptome
# under Spike-in).
insert_html <- ""
is_path <- file.path(qc_dir, paste0(project, ".insertsize.txt"))
if (file.exists(is_path)) {
  ins <- fread(is_path)
  # window to plausible fragment lengths BEFORE any fraction is computed, so each facet's
  # fraction integrates over the plotted 100 bp-10 kb window (spec §5).
  ins <- ins[insert_size >= 100 & insert_size <= 10000]
  if (nrow(ins)) {
    rts  <- intersect(rt_levels, unique(ins$read_type))
    # drive the category factor from the levels PRESENT in the data (droplevels), so a
    # facet ROW is only rendered for categories that actually occur (no empty Spike-in row
    # when there are no spike/additional genes) rather than a hardcoded 2-level factor.
    ins[, category  := droplevels(factor(category, levels = facet_cats(ins)))]
    ins[, read_type := factor(read_type, levels = rts)]
    cats <- levels(ins$category)
    ins[, frac := count / sum(count), by = .(category, read_type)]
    ensure_ggh4x()
    isize <- ggplot(ins, aes(insert_size, frac)) +
      geom_area(fill = "#1A5084", alpha = 0.18) +
      geom_line(colour = "#1A5084", linewidth = 0.4) +
      ggh4x::facet_grid2(vars(category), vars(read_type),
                         scales = "free_y", independent = "y") +
      scale_x_log10(limits = c(100, 10000),
                    breaks = c(100, 250, 500, 1000, 2500, 5000, 10000),
                    labels = c("100", "250", "500", "1k", "2.5k", "5k", "10k")) +
      labs(x = "Insert size (bp, log scale)", y = "Fraction of fragments",
           caption = "Transcript-space fragment length (exon-union projection); window 100 bp-10 kb, log x, free y per facet.") +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.5),
            plot.caption = element_text(size = 6.5, colour = "#666", hjust = 0))
    insert_html <- sprintf(
      "<section class='panel insert'><h2>Insert-size distribution</h2>%s</section>",
      img(isize, 7.4, 0.7 * length(cats) + 0.9))
  }
}

## assemble the single A4 HTML file ------------------------------------------
css <- "
@page{size:A4 portrait;margin:8mm}
*{box-sizing:border-box}
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;color:#222;margin:0;
  -webkit-print-color-adjust:exact;print-color-adjust:exact}
.report{max-width:194mm;margin:4px auto;padding:0 8px}
h1{font-size:15px;margin:0 0 5px;border-bottom:2px solid #1A5084;padding-bottom:3px}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:6px;
  grid-template-areas:'summary summary' 'drops drops' 'counts pctbox' 'compos compos' 'cover cover' 'insert insert' 'setup setup'}
.summary{grid-area:summary}.drops{grid-area:drops}.counts{grid-area:counts}.pctbox{grid-area:pctbox}
.compos{grid-area:compos}.cover{grid-area:cover}.insert{grid-area:insert}.setup{grid-area:setup}
.panel h2{font-size:12px;color:#1A5084;margin:0 0 3px}
.panel h3{font-size:11px;margin:4px 0 1px;color:#555}
.panel img{width:100%;height:auto;display:block;border:1px solid #eee;border-radius:3px}
table{border-collapse:collapse;width:100%}
td{padding:2px 6px;border-bottom:1px solid #eee;font-size:11px}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.meta td:first-child,.stats td:first-child{color:#555}
th{padding:2px 6px;border-bottom:1px solid #ccc;font-size:11px;text-align:left;color:#555;font-weight:600}
th.num{text-align:right}
.hier td.cat{color:#333}
.hier tr.lvl0 .cat{font-weight:600;color:#222}
.hier tr.lvl1 .cat{padding-left:16px}
.hier tr.lvl2 .cat{padding-left:32px;color:#666}
.sumcols{display:flex;gap:16px}.sumcol{flex:1 1 0;min-width:0}
.refcol{flex:0 1 34%}.statcol{flex:1 1 0;min-width:0}
.cmd{margin:2px 0;font-size:10px}
.cmd code{background:#f4f4f4;padding:1px 4px;border-radius:3px;word-break:break-all}
.kv{column-count:2;column-gap:14px;margin:3px 0}
.kv>div{break-inside:avoid;display:flex;justify-content:space-between;gap:8px;
  border-bottom:1px solid #eee;padding:1px 0;font-size:10px}
.kv dt{color:#555;white-space:nowrap}.kv dd{margin:0;text-align:right;overflow:hidden;text-overflow:ellipsis}
.dge{font-size:10px;color:#333;margin:4px 0 0;word-break:break-word}
@media print{.report{max-width:none;margin:0;padding:0}}
"
body <- paste0(summary_html, drops_html, counts_html, compos_html, pctbox_html, cover_html, insert_html, setup_html)
html <- sprintf(
  "<!doctype html><html><head><meta charset='utf-8'>
   <title>QUARTSx QC report &mdash; %s</title><style>%s</style></head>
   <body><div class='report'><h1>QUARTSx QC report &mdash; %s</h1>
   <div class='grid'>%s</div></div></body></html>",
  esc(project), css, esc(project), body)

out_html <- file.path(qc_dir, paste0(project, ".qc_report.html"))
writeLines(html, out_html)
cat("wrote", out_html, "\n")
