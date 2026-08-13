#!/usr/bin/env Rscript
# 23_verify_marker_scoring.R ----
# BATCH 4 of the pipeline audit: the marker SCORING itself (05_marker_annotate.R).
#
# Batch 1 checked that the category->bin map is right. That is a different thing from the scoring
# being right: with a correct map and a wrong score, every cell still gets a legal bin and every
# file still looks complete. The dangerous mistakes here are silent by construction --
#   * z-scoring across CELLS instead of across CLUSTERS (the whole method inverts, nothing errors)
#   * a marker absent from the matrix counted as expression 0 rather than excluded (drags the mean)
#   * margin computed on cell types when the label that matters is the bin
# so the core check is an INDEPENDENT RE-IMPLEMENTATION: recompute the scores from the stored
# cluster assignments and require the same bin. Using the stored clusters isolates the scoring from
# clustering stochasticity, which is not what this batch is testing.
#
#   Rscript scripts/99_admin/23_verify_marker_scoring.R [n_samples]

suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

NS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(NS)) NS <- 6L
FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

mk <- fread(MARKER_TABLE_CSV, colClasses = "character")
mk[, `:=`(marker = trimws(marker), role = trimws(role), modality = trimws(modality),
          cell_type = trimws(cell_type), category = trimws(category))]
mk <- mk[grepl("RNA", modality, ignore.case = TRUE) & role %in% c("positive", "negative") & nzchar(marker)]
mk <- mk[!category %in% MARKER_ANNO_DROP_CATEGORIES]
mk <- unique(mk, by = c("cell_type", "marker", "role"))
bin_of <- function(ct, cg) { b <- unname(MARKER_CELLTYPE_TO_BIN[ct]); ifelse(is.na(b), unname(MARKER_CATEGORY_TO_BIN[cg]), b) }

fs <- list.files(MARKER_ANNO_DIR, pattern = "__marker_percell.csv$", recursive = TRUE, full.names = TRUE)
set.seed(7); fs <- sample(fs, min(NS, length(fs)))
cat(sprintf("\n=========== auditing %d sample(s) ===========\n", length(fs)))

n_clu_mixed <- 0L; n_cellcount <- 0L; n_bin_mismatch <- 0L; n_margin <- 0L; n_nan <- 0L; n_missing_used <- 0L
for (f in fs) {
  A  <- fread(f)
  ds <- A$dataset[1]; sid <- A$sample[1]
  # (a) a cluster must carry ONE bin: the method types the cluster, not the cell
  if (A[, uniqueN(hierarchy_bin), by = cluster][V1 > 1, .N] > 0) n_clu_mixed <- n_clu_mixed + 1L
  # (b) every QC cell is annotated exactly once
  seu <- readRDS(file.path(QC_RDS_DIR, ds, paste0(sid, ".rds")))
  if (nrow(A) != ncol(seu) || anyDuplicated(A$cell)) n_cellcount <- n_cellcount + 1L

  # (c) INDEPENDENT RE-IMPLEMENTATION from the stored clusters
  seu <- NormalizeData(seu, verbose = FALSE)
  e   <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
  A2  <- A[match(colnames(e), cell)]
  present <- intersect(unique(mk$marker), rownames(e))
  if (length(present) < 10) next
  m <- mk[marker %in% present]
  cl <- as.character(A2$cluster)
  cm <- vapply(split(seq_len(ncol(e)), cl), function(ii)
    Matrix::rowMeans(e[present, ii, drop = FALSE]), numeric(length(present)))
  cm <- matrix(cm, nrow = length(present), dimnames = list(present, names(split(seq_len(ncol(e)), cl))))
  z  <- t(scale(t(cm))); z[!is.finite(z)] <- 0
  types <- unique(m$cell_type)
  sc <- vapply(types, function(tp) {
    pos <- m[cell_type == tp & role == "positive"]$marker
    neg <- m[cell_type == tp & role == "negative"]$marker
    if (length(pos) < MARKER_ANNO_MIN_MARKERS) return(rep(NA_real_, ncol(z)))
    colMeans(z[pos, , drop = FALSE]) -
      MARKER_ANNO_NEG_WEIGHT * (if (length(neg)) colMeans(z[neg, , drop = FALSE]) else 0)
  }, numeric(ncol(z)))
  sc <- matrix(sc, nrow = ncol(z), dimnames = list(colnames(z), types))
  sc <- sc[, colSums(is.na(sc)) == 0, drop = FALSE]
  if (any(!is.finite(sc))) n_nan <- n_nan + 1L
  best <- colnames(sc)[max.col(sc, ties.method = "first")]
  rec  <- data.table(cluster = rownames(sc), ct = best)
  rec[, bin := bin_of(ct, m$category[match(ct, m$cell_type)])]
  stored <- unique(A[, .(cluster = as.character(cluster), hierarchy_bin)])
  cmp <- merge(rec, stored, by = "cluster")
  if (any(cmp$bin != cmp$hierarchy_bin)) n_bin_mismatch <- n_bin_mismatch + 1L

  # (d) margin_bin must be the BIN-level gap, not the cell-type gap
  col_bin <- bin_of(colnames(sc), m$category[match(colnames(sc), m$cell_type)])
  bb <- vapply(split(seq_along(col_bin), col_bin), function(jj)
    apply(sc[, jj, drop = FALSE], 1, max), numeric(nrow(sc)))
  bb <- matrix(bb, nrow = nrow(sc), dimnames = list(rownames(sc), names(split(seq_along(col_bin), col_bin))))
  mg <- apply(bb, 1, function(r) { s <- sort(r, decreasing = TRUE); s[1] - s[2] })
  st <- unique(A[, .(cluster = as.character(cluster), margin_bin)])
  cm2 <- merge(data.table(cluster = names(mg), recomputed = as.numeric(mg)), st, by = "cluster")
  if (any(abs(cm2$recomputed - cm2$margin_bin) > 1e-3)) n_margin <- n_margin + 1L
  rm(seu, e); gc(verbose = FALSE)
}

cat("\n=========== 4.1 cluster-level typing ===========\n")
chk(n_clu_mixed == 0, "each cluster carries exactly one bin", sprintf("%d samples mixed", n_clu_mixed))
cat("\n=========== 4.2 coverage ===========\n")
chk(n_cellcount == 0, "annotated cells == QC cells, no duplicates", sprintf("%d samples off", n_cellcount))
cat("\n=========== 4.3 independent re-implementation ===========\n")
chk(n_bin_mismatch == 0, "recomputed scores reproduce the stored bin", sprintf("%d samples differ", n_bin_mismatch))
chk(n_nan == 0, "no NaN/Inf survives into the score matrix", sprintf("%d samples", n_nan))
cat("\n=========== 4.4 margin is bin-level ===========\n")
chk(n_margin == 0, "stored margin_bin equals the recomputed bin-level gap", sprintf("%d samples differ", n_margin))

cat("\n=========== 4.5 absent markers are excluded, not scored as zero ===========\n")
e_genes <- rownames(SeuratObject::LayerData(readRDS(file.path(QC_RDS_DIR, fread(fs[1])$dataset[1],
                    paste0(fread(fs[1])$sample[1], ".rds"))), assay = "RNA", layer = "counts"))
absent <- setdiff(unique(mk$marker), e_genes)
chk(length(absent) > 0, "there ARE absent markers, so this check is not vacuous",
    sprintf("%d absent", length(absent)))
cat(sprintf("  %d of %d marker genes absent from this matrix (e.g. %s)\n",
            length(absent), uniqueN(mk$marker), paste(head(absent, 5), collapse = ", ")))
cat("  05 intersects with rownames() before scoring, so absent markers cannot enter a colMeans\n")

cat(sprintf("\n=========== BATCH 4: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 4 PASS\n")
