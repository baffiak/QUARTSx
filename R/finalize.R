#!/usr/bin/env Rscript
# Build the nested dgecounts list from the Rust counts handoff, saveRDS, drop the handoff.
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
                 colClasses = list(character = c("quant", "feature", "level", "gene", "cell")))

# one sparse gene x cell matrix from the rows of a single (quant,feature,level) group
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
  if (quant == "rpkm") m <- as.matrix(m)
  if (level == "all") {
    dgecounts[[quant]][[feature]]$all <- m
  } else {
    dgecounts[[quant]][[feature]]$downsampling[[level]] <- m
  }
}

saveRDS(dgecounts, file.path(out_dir, "expression", paste0(project, ".dgecounts.rds")))
invisible(file.remove(handoff))
