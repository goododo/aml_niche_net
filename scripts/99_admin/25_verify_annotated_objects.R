#!/usr/bin/env Rscript
# 25_verify_annotated_objects.R ----
# BATCH 6 of the pipeline audit: the annotated Seurat objects written by 07.
#
# Two things are being checked, and the second matters as much as the first.
#
# (1) The metadata is joined by cell id, not by position. A position-based join produces an object
#     where every cell carries its neighbour's label: nothing is NA, every count is right, and the
#     annotation is silently scrambled. So the check is cell-by-cell equality against the source
#     table, not a row count.
#
# (2) The QC objects were NOT modified. The entire justification for writing to 04_annotated
#     instead of back into 01_per_sample_qc is that touching those files would mark all 179
#     finished inferCNV runs stale via the mtime guard in 44_infercnv_run_one.R. That claim is
#     only worth anything if it is verified rather than asserted.
#
#   Rscript scripts/99_admin/25_verify_annotated_objects.R [n_samples]

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

NS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(NS)) NS <- 10L
FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

R <- qc_rds_roster(on_extra = "ignore")
objs <- list.files(ANNO_OBJ_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)

cat("\n=========== 6.1 coverage ===========\n")
chk(length(objs) == nrow(R), "one annotated object per roster sample",
    sprintf("objects %d vs roster %d", length(objs), nrow(R)))

cat("\n=========== 6.2 QC objects were NOT modified ===========\n")
# The annotated objects were all written after the QC objects were last touched; if any QC object
# is NEWER than the annotation run, something wrote to it.
run_start <- min(file.mtime(objs))
qc_newer  <- R[file.mtime(rds) > run_start]
chk(nrow(qc_newer) == 0, "no QC object was written during/after the annotation run",
    sprintf("%d newer", nrow(qc_newer)))
# and the inferCNV freshness guard must still consider every burden CSV current.
#
# THE ROSTER IS R, NOT THE TASK FILES. This used to glob ^infercnv_tasks(_healthy)?\.tsv$, which
# matches 2 of the 6 task files in the tree and covers 179 of the 212 samples that have burden
# CSVs -- so 33 samples were a blind spot while the check reported PASS. That is the same failure
# 22_verify_malignancy_percell.R already hit and documents: a roster extended by hand every time
# work is queued silently shrinks the audit exactly when new work needs auditing. R comes from
# qc_rds_roster(), the roster every producing script uses, so it cannot drift from them.
tasks <- R[, .(dataset, sample, q = rds)]
tasks[, b := file.path(INFERCNV_BURDEN_ROOT, dataset, paste0(sample, "_infercnv_burden.csv"))]
tk <- tasks[file.exists(b) & file.exists(q)]
# NON-VACUITY: a check that reaches 0 samples passes for the wrong reason.
chk(nrow(tk) >= 200, "the burden-freshness check actually reaches the cohort",
    sprintf("%d of %d roster samples have both a burden CSV and a QC object", nrow(tk), nrow(tasks)))
stale <- tk[file.mtime(b) < file.mtime(q)]
chk(nrow(stale) == 0, "no burden CSV became stale (the whole reason for a separate directory)",
    sprintf("%d would recompute", nrow(stale)))

cat("\n=========== 6.3 metadata is joined by cell id, cell by cell ===========\n")
set.seed(5); pick <- sample(objs, min(NS, length(objs)))
n_align <- 0L; n_mismatch <- 0L; n_na <- 0L; n_assay <- 0L; n_missing_col <- 0L
KEY <- c("hierarchy_bin", "anno_source", "bin_bmm", "bin_marker", "in_ccc_graph")
for (p in pick) {
  seu <- readRDS(p)
  ds  <- basename(dirname(p)); sid <- sub("\\.rds$", "", basename(p))
  A   <- fread(file.path(ANNO_RECONCILED_DIR, ds, paste0(sid, "__anno_percell.csv")))
  # 07 normalises "" -> NA when it loads this CSV (fwrite writes NA as an empty field), so compare
  # against the same normalisation or every blank bin reads as a mismatch that is not one.
  for (cc in names(A)) if (is.character(A[[cc]])) set(A, which(!nzchar(trimws(A[[cc]]))), cc, NA_character_)
  if (!identical(rownames(seu@meta.data), colnames(seu))) n_align <- n_align + 1L
  if (!all(KEY %in% names(seu@meta.data))) { n_missing_col <- n_missing_col + 1L; next }
  A <- A[match(colnames(seu), cell)]
  for (k in KEY) {
    a <- as.character(A[[k]]); b <- as.character(seu@meta.data[[k]])
    if (!identical(a, b)) { n_mismatch <- n_mismatch + 1L; break }
  }
  if (anyNA(seu$hierarchy_bin) || any(!seu$hierarchy_bin %in% c(HIERARCHY_BINS, "unassigned")))
    n_na <- n_na + 1L
  # the expression data must survive the metadata edit untouched
  cnt <- tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(cnt) || ncol(cnt) != ncol(seu)) n_assay <- n_assay + 1L
  rm(seu); gc(verbose = FALSE)
}
chk(n_align == 0, "meta.data rownames align with colnames", sprintf("%d objects", n_align))
chk(n_missing_col == 0, "all key annotation columns present", sprintf("%d objects", n_missing_col))
chk(n_mismatch == 0, "every cell's annotation equals the reconciled table", sprintf("%d objects differ", n_mismatch))
chk(n_na == 0, "hierarchy_bin has no NA and only legal values", sprintf("%d objects", n_na))
chk(n_assay == 0, "counts assay intact and same width", sprintf("%d objects", n_assay))

cat("\n=========== 6.4 the join is not vacuously correct ===========\n")
# If a sample's cells were all one bin, a scrambled join would still compare equal. Confirm the
# audited objects actually carry several bins, so 6.3 has something to catch.
seu <- readRDS(pick[1])
nb <- length(unique(seu$hierarchy_bin))
chk(nb > 1, "the audited object carries more than one bin (so a scrambled join WOULD differ)",
    sprintf("%d bins", nb))
cat(sprintf("  %s: %d cells, %d bins, sources: %s\n", basename(pick[1]), ncol(seu), nb,
            paste(sprintf("%s=%d", names(table(seu$anno_source)), table(seu$anno_source)), collapse = " ")))

cat(sprintf("\n=========== BATCH 6: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 6 PASS\n")
