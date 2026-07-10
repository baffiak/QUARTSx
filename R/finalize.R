#!/usr/bin/env Rscript
# Build the nested dgecounts list from the Rust counts handoff, saveRDS, drop the handoff.
#
# Input handoff (out/expression/{project}.counts.tsv.gz) is a single long-format table written by
# src/dedup.rs::emit  ->  columns: quant  feature  level  gene  cell  value
# Rows are byte-deterministic (canonical gene/cell order); one row per non-zero matrix entry.
#
# Nesting keys are `quant` -> `feature` -> `level`, taken verbatim from the handoff, so the loop
# below is fully generic and every block lands where its labels put it -- no per-quant branch.
#
# Integer blocks:  umicount / readcount / readcount_internal / rpkm.
# Resolver blocks (multimapper DECISION, build-spec §4): for the ONE selected mode
# (uniform|rescue|propunique|em) count.rs emits two families that MIRROR the integer matrix on all
# three axes -- quant `umicount_mult_<mode>` and `readcount_internal_mult_<mode>`, feature in
# {exon,intron,inex} (whichever layers the integer matrix used), level `all` plus every
# `downsampled_<label>`. Their `value` is FRACTIONAL (a multi-gene unit split across its |ug| genes);
# downstream reconstructs STARsolo `UniqueAndMult-<mode>` as `umicount + umicount_mult_<mode>` and
# `readcount_internal + readcount_internal_mult_<mode>`. The generic loop nests these exactly like the
# integer families with NO special case; the only hazard is truncating the fractional values, which is
# why `value` is pinned numeric on read below.
#
# NOTE: gene_id<->gene_name mapping lives in out/expression/{project}.gene_names.txt (written with a
# `gene_id\tgene_name` header by count.rs). finalize.R does NOT read that file -- dgecounts is keyed
# by gene_id straight from the handoff -- so the added header row is a no-op here.
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: finalize.R <out_dir> <project>")
out_dir <- args[1]
project <- args[2]

handoff <- file.path(out_dir, "expression", paste0(project, ".counts.tsv.gz"))
counts  <- fread(cmd = paste("gzip -dc", shQuote(handoff)),
                 colClasses = list(character = c("quant", "feature", "level", "gene", "cell"),
                                   numeric   = "value"))

# one sparse gene x cell matrix from the rows of a single (quant,feature,level) group
# (dgCMatrix stores x as double, so fractional *_mult_<mode> resolver values are preserved as-is)
make_matrix <- function(d) {
  genes <- unique(d$gene)
  cells <- unique(d$cell)
  sparseMatrix(i = match(d$gene, genes), j = match(d$cell, cells), x = d$value,
               dims = c(length(genes), length(cells)), dimnames = list(genes, cells))
}

dgecounts <- list()
for (grp in split(counts, by = c("quant", "feature", "level"), drop = TRUE)) {
  quant <- grp$quant[1]; feature <- grp$feature[1]; level <- grp$level[1]
  m <- make_matrix(grp)
  # rpkm is the only dense block; every other quant (umicount/readcount/readcount_internal and the
  # umicount_mult_<mode>/readcount_internal_mult_<mode> UniqueAndMult resolver matrices) stays sparse.
  if (quant == "rpkm") m <- as.matrix(m)
  if (level == "all") {
    dgecounts[[quant]][[feature]]$all <- m
  } else {
    dgecounts[[quant]][[feature]]$downsampling[[level]] <- m
  }
}

saveRDS(dgecounts, file.path(out_dir, "expression", paste0(project, ".dgecounts.rds")))
invisible(file.remove(handoff))
