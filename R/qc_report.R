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
# mid: filtered-read stats
fs_path <- file.path(out_dir, "filtered", "filter_stats.json")
stat_rows <- ""
if (file.exists(fs_path)) {
  fs <- jsonlite::read_json(fs_path, simplifyVector = TRUE)
  rows <- list(c("Total input reads", fs$total), c("Passed filter", fs$passed),
               c("UMI-tagged", fs$tagged), c("Internal", fs$internal),
               c("Dropped", fs$total - fs$passed))
  stat_rows <- paste(vapply(rows, function(r)
    sprintf("<tr><td>%s</td><td class='num'>%s</td></tr>", r[1], comma(r[2])), ""), collapse = "")
}
# run setup section (goes at the BOTTOM of the report): command + only the
# IMPORTANT, result-affecting parameters (no compute/rds/downsampling noise)
setup_html <- ""
if (nzchar(config_path) && file.exists(config_path)) {
  y  <- yaml::read_yaml(config_path)
  rf <- y$read_filtering; fc <- y$filter_cutoffs; bc <- y$barcodes; co <- y$counting_opts
  ref_rows <- paste0(
    kv("Genome (STAR index)", name_or(y$reference$STAR_index)),
    kv("Annotation (GTF)",    name_or(y$reference$GTF_file)),
    kv("Additional files",    name_or(y$reference$additional_files)))
  items <- c(
    `R2 adapter`      = name_or(rf$adapter_fasta),
    `R2 quality`      = rf$quality,
    `R2 min length`   = rf$min_length,
    `BC filter`       = sprintf("%s bp < Q%s", fc$BC_filter$num_bases,  fc$BC_filter$phred),
    `UMI filter`      = sprintf("%s bp < Q%s", fc$UMI_filter$num_bases, fc$UMI_filter$phred),
    `Barcode binning` = bc$BarcodeBinning,
    `Min reads/cell`  = bc$nReadsperCell,
    introns = co$introns, strand = co$strand)
  dl <- paste(sprintf("<div><dt>%s</dt><dd>%s</dd></div>", names(items), esc(as.character(items))), collapse = "")
  setup_html <- sprintf(
    "<section class='panel setup'><h2>Run setup</h2>
       <p class='cmd'><code>quartsx.sh -c -y %s</code></p>
       <dl class='kv'>%s</dl></section>",
    esc(config_path), dl)
}
summary_html <- sprintf(
  "<section class='panel summary'><h2>Run summary</h2>
     <div class='sumcols'>
       <div class='sumcol'><h3>Reference inputs</h3><table class='meta'>%s</table></div>
       <div class='sumcol'><h3>Filtered reads</h3><table class='stats'>%s</table></div>
     </div></section>",
  ref_rows, stat_rows)

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
sumBar[, perc := 100 * tot / sum(tot)]
sumBar[, type := factor(type, levels = bar_levels)]
bar <- ggplot(sumBar, aes(x = 1, y = perc, fill = type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = featColors, guide = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "% of total reads") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), legend.title = element_blank(),
        legend.position = "bottom", legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

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
cover_html <- ""
gb_path <- file.path(qc_dir, paste0(project, ".genebody_coverage.txt"))
if (file.exists(gb_path)) {
  gb <- fread(gb_path)
  if (nrow(gb)) {
    cats <- facet_cats(gb)
    gb[, category  := factor(category,  levels = cats)]
    gb[, read_type := factor(read_type, levels = intersect(rt_levels, unique(read_type)))]
    gcov <- ggplot(gb, aes(percentile, coverage)) +
      geom_line(colour = "#118730", linewidth = 0.6) +
      facet_grid(category ~ read_type) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "Gene body percentile (5'->3')", y = "Norm. Coverage") +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.5))
    cover_html <- sprintf(
      "<section class='panel cover'><h2>Gene-body coverage (5'&rarr;3')</h2>%s</section>",
      img(gcov, 7.4, 0.7 * length(cats) + 0.7))
  }
}

## (4) insert-size distribution, same 2x2 split as coverage ------------------
insert_html <- ""
is_path <- file.path(qc_dir, paste0(project, ".insertsize.txt"))
if (file.exists(is_path)) {
  ins <- fread(is_path)
  if (nrow(ins)) {
    cats <- facet_cats(ins)
    ins[, category  := factor(category,  levels = cats)]
    ins[, read_type := factor(read_type, levels = intersect(rt_levels, unique(read_type)))]
    ins[, frac := count / sum(count), by = .(category, read_type)]
    xmax <- ins[, { o <- order(insert_size); insert_size[o][which(cumsum(count[o]) >= 0.99 * sum(count))[1]] }]
    isize <- ggplot(ins, aes(insert_size, frac)) +
      geom_col(fill = "#1A5084", width = 1) +
      facet_grid(category ~ read_type) +
      coord_cartesian(xlim = c(0, xmax + 1)) +
      labs(x = "Insert size (bp)", y = "Fraction of fragments") +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 8.5))
    insert_html <- sprintf(
      "<section class='panel insert'><h2>Insert-size distribution</h2>%s</section>",
      img(isize, 7.4, 0.7 * length(cats) + 0.7))
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
  grid-template-areas:'summary summary' 'counts pctbox' 'compos compos' 'cover cover' 'insert insert' 'setup setup'}
.summary{grid-area:summary}.counts{grid-area:counts}.pctbox{grid-area:pctbox}
.compos{grid-area:compos}.cover{grid-area:cover}.insert{grid-area:insert}.setup{grid-area:setup}
.panel h2{font-size:12px;color:#1A5084;margin:0 0 3px}
.panel h3{font-size:11px;margin:4px 0 1px;color:#555}
.panel img{width:100%;height:auto;display:block;border:1px solid #eee;border-radius:3px}
table{border-collapse:collapse;width:100%}
td{padding:2px 6px;border-bottom:1px solid #eee;font-size:11px}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.meta td:first-child,.stats td:first-child{color:#555}
.sumcols{display:flex;gap:16px}.sumcol{flex:1 1 0;min-width:0}
.cmd{margin:2px 0;font-size:10px}
.cmd code{background:#f4f4f4;padding:1px 4px;border-radius:3px;word-break:break-all}
.kv{column-count:2;column-gap:14px;margin:3px 0}
.kv>div{break-inside:avoid;display:flex;justify-content:space-between;gap:8px;
  border-bottom:1px solid #eee;padding:1px 0;font-size:10px}
.kv dt{color:#555;white-space:nowrap}.kv dd{margin:0;text-align:right;overflow:hidden;text-overflow:ellipsis}
.dge{font-size:10px;color:#333;margin:4px 0 0;word-break:break-word}
@media print{.report{max-width:none;margin:0;padding:0}}
"
body <- paste0(summary_html, counts_html, compos_html, pctbox_html, cover_html, insert_html, setup_html)
html <- sprintf(
  "<!doctype html><html><head><meta charset='utf-8'>
   <title>QUARTSx QC report &mdash; %s</title><style>%s</style></head>
   <body><div class='report'><h1>QUARTSx QC report &mdash; %s</h1>
   <div class='grid'>%s</div></div></body></html>",
  esc(project), css, esc(project), body)

out_html <- file.path(qc_dir, paste0(project, ".qc_report.html"))
writeLines(html, out_html)
cat("wrote", out_html, "\n")
