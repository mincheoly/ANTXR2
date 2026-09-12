# Export raw counts + metadata + tSNE from the IBD colon atlas's discovery-
# cohort Seurat objects (train.{Epi,Fib,Imm}.seur.rds), for two purposes:
#
# 1. Epi and Imm: their full-cohort gene_sorted-*.matrix.mtx files downloaded
#    from HCA are CORRUPTED AT THE SOURCE (verified via independent S3
#    byte-range probing, not a download artifact -- see ANALYSIS_SUMMARY.md
#    and config.py's comment above IBD_COLON_RDS_* for the full story: Epi is
#    truncated to ~53% of its declared line count, Imm to ~4.6%, both ending
#    mid-line with no trailing newline). These RDS objects are this project's
#    only complete alternative raw-count source for those two compartments --
#    but they only cover the "discovery cohort" (17 of 30 donors), not the
#    full cohort Fib's intact mtx provides.
# 2. All 3 compartments: pre-computed tSNE coordinates, as the "processed with
#    visualization coordinates" artifact the user asked for (also
#    discovery-cohort-only, since that's what these objects contain).
#
# Deliberately does NOT require the Seurat package (not installed in r-env,
# and installing it is a heavy, slow build) -- a Seurat object's S4 slots are
# readable directly via attr()/readRDS() without the defining package loaded,
# since dgCMatrix (the counts slot's class) only needs the much lighter
# `Matrix` package, already installed. Confirmed working in interactive
# testing before writing this as a script.
#
# Usage: conda run -n r-env Rscript export_ibd_colon_rds.R
suppressMessages(library(Matrix))

raw_dir <- "/data/ANTXR2/raw/ibd_colon_atlas"
out_dir <- file.path(raw_dir, "rds_export")
dir.create(out_dir, showWarnings = FALSE)

compartments <- c("Epi", "Fib", "Imm")

for (comp in compartments) {
  cat("===", comp, "===\n")
  obj <- readRDS(file.path(raw_dir, paste0("train.", comp, ".seur.rds")))

  rna <- attr(attr(obj, "assays"), "names")
  rna_assay <- attr(obj, "assays")[["RNA"]]
  counts <- attr(rna_assay, "counts")
  meta <- attr(obj, "meta.data")
  reductions <- attr(obj, "reductions")
  tsne <- reductions[["tsne"]]
  emb <- if (!is.null(tsne)) attr(tsne, "cell.embeddings") else NULL

  cat("  counts:", nrow(counts), "genes x", ncol(counts), "cells\n")
  cat("  meta.data:", nrow(meta), "x", ncol(meta), "-- cols:", paste(colnames(meta), collapse=","), "\n")
  cat("  tsne:", if (!is.null(emb)) paste(dim(emb), collapse=" x ") else "NONE", "\n")

  comp_dir <- file.path(out_dir, comp)
  dir.create(comp_dir, showWarnings = FALSE)

  writeMM(counts, file.path(comp_dir, "counts.mtx"))  # genes x cells, same orientation as the HCA mtx files
  writeLines(rownames(counts), file.path(comp_dir, "genes.tsv"))
  writeLines(colnames(counts), file.path(comp_dir, "barcodes.tsv"))
  write.csv(meta, file.path(comp_dir, "meta.csv"), row.names = TRUE)
  if (!is.null(emb)) {
    write.csv(emb, file.path(comp_dir, "tsne.csv"), row.names = TRUE)
  }
  cat("  wrote to", comp_dir, "\n")
}
cat("done\n")
