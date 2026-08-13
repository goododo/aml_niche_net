#!/usr/bin/env Rscript
# 22_verify_malignancy_percell.R ----
# BATCH 3 of the pipeline audit: burden -> per-cell call (41) -> consensus (50) -> rollup (60).
#
# The failure this batch is built around did not crash anything. A driver script parsed a
# tab-separated job list with `IFS=$'\t' read`, and because tab is whitespace bash collapsed the
# empty fields -- so the inferCNV per-cell path was handed to 50 as --numbat. Every sample still
# produced a complete, well-formed consensus record. The only visible symptom was that 128 samples
# claimed a numbat arm when just 8 had usable allele data: the OUTPUT disagreed with the INPUTS
# on disk. That class of error is what 3.3 checks for, and it is why "it ran and produced files"
# is not evidence of anything.
#
#   Rscript scripts/99_admin/22_verify_malignancy_percell.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

tasks <- rbindlist(lapply(list.files(".", pattern = "^infercnv_tasks(_healthy)?\\.tsv$"),
                          function(p) fread(p, header = FALSE, col.names = c("dataset", "sample"))))
tasks[, burden   := file.path(INFERCNV_BURDEN_ROOT, dataset, paste0(sample, "_infercnv_burden.csv"))]
tasks[, percell  := file.path(INFERCNV_ROOT, dataset, sample, paste0(sample, "__infercnv_percell.csv"))]
tasks[, cons     := file.path(DIR_MALIGNANCY, dataset, paste0(sample, "__consensus_summary.csv"))]

cat("\n=========== 3.0 which samples are mid-recompute (excluded from consistency) ===========\n")
tasks[, has_b := file.exists(burden)][, has_p := file.exists(percell)]
tasks[, pending := has_b & has_p & (file.mtime(burden) > file.mtime(percell))]
cat(sprintf("  samples: %d | burden present: %d | percell present: %d | burden newer than percell: %d\n",
            nrow(tasks), sum(tasks$has_b), sum(tasks$has_p), sum(tasks$pending)))
if (any(tasks$pending))
  cat("  (these are being re-scored right now; 41 has not caught up, so they are reported, not failed)\n")
stable <- tasks[has_b & has_p & !pending]
chk(nrow(stable) > 0, "there is a stable set to check")

cat("\n=========== 3.1 per-cell tables contain ONLY the sample's own cells ===========\n")
set.seed(11); samp <- stable[sample(.N, min(25L, .N))]
viol_foreign <- 0L; viol_count <- 0L; viol_thr <- 0L
for (i in seq_len(nrow(samp))) {
  b <- fread(samp$burden[i], select = c("cell", "infercnv_burden", "group"))
  p <- fread(samp$percell[i], select = c("cell", "malignant", "score"))
  own <- b[group %in% c("observation", "reference_normal")]
  if (nrow(p) != nrow(own)) viol_count <- viol_count + 1L
  if (any(p$cell %in% b[!group %in% c("observation", "reference_normal")]$cell)) viol_foreign <- viol_foreign + 1L
  is_ref <- grepl("reference", b$group, ignore.case = TRUE)
  thr <- as.numeric(quantile(b$infercnv_burden[is_ref], INFERCNV_SCORE_Q, na.rm = TRUE))
  # reference_normal cells are forced 0 regardless of score, so only observations are re-derivable
  o <- merge(p, b[group == "observation", .(cell, infercnv_burden)], by = "cell")
  if (nrow(o) && any(o$malignant != as.integer(o$infercnv_burden > thr))) viol_thr <- viol_thr + 1L
}
chk(viol_foreign == 0, "no foreign reference cell leaks into the per-cell table",
    sprintf("%d/%d samples", viol_foreign, nrow(samp)))
chk(viol_count == 0, "per-cell row count == observation + autologous reference",
    sprintf("%d/%d samples", viol_count, nrow(samp)))
chk(viol_thr == 0, "malignant flag reproduces the reference-P95 threshold",
    sprintf("%d/%d samples", viol_thr, nrow(samp)))

cat("\n=========== 3.2 reference cells are never called malignant ===========\n")
bad_ref <- 0L
for (i in seq_len(nrow(samp))) {
  b <- fread(samp$burden[i], select = c("cell", "group"))
  p <- fread(samp$percell[i], select = c("cell", "malignant"))
  refc <- b[group == "reference_normal"]$cell
  if (any(p[cell %in% refc]$malignant == 1L)) bad_ref <- bad_ref + 1L
}
chk(bad_ref == 0, "autologous reference cells are all malignant=0", sprintf("%d samples", bad_ref))

cat("\n=========== 3.3 consensus ARMS match the evidence actually on disk ===========\n")
usable_numbat <- function(ds, sid) {
  f <- file.path(LARGE1_DIR, "03_cnv_snv", "numbat", ds, sid, "numbat", paste0(sid, "__numbat_percell.csv"))
  if (!file.exists(f)) return(FALSE)
  d <- tryCatch(fread(f, nrows = 5), error = function(e) NULL)
  !is.null(d) && "compartment_opt" %in% names(d)      # degraded files lack the clone call
}
has_author <- function(ds, sid)
  file.exists(file.path(LARGE1_DIR, "03_cnv_snv", "author", ds, sid, paste0(sid, "__author_percell.csv")))

S <- tasks[file.exists(cons)]
S[, recorded := vapply(seq_len(.N), function(i) as.character(fread(cons[i])$arms[1]), character(1))]
S[, exp_numbat := mapply(usable_numbat, dataset, sample)]
S[, exp_author := mapply(has_author, dataset, sample)]
S[, rec_numbat := grepl("numbat", recorded)]
S[, rec_author := grepl("author", recorded)]
S[, rec_icnv   := grepl("infercnv", recorded)]
chk(all(S$rec_icnv), "every consensus record claims the inferCNV arm",
    sprintf("%d missing", sum(!S$rec_icnv)))
chk(all(S$rec_numbat == S$exp_numbat), "recorded numbat arm == usable numbat file on disk",
    sprintf("%d mismatched", sum(S$rec_numbat != S$exp_numbat)))
chk(all(S$rec_author == S$exp_author), "recorded author arm == author file on disk",
    sprintf("%d mismatched", sum(S$rec_author != S$exp_author)))
cat(sprintf("  arms observed: %s\n",
            paste(sprintf("%s=%d", names(table(S$recorded)), table(S$recorded)), collapse = "  ")))

cat("\n=========== 3.4 degraded numbat abstains (does not vote normal) ===========\n")
deg <- tasks[, .(dataset, sample)][, f := file.path(LARGE1_DIR, "03_cnv_snv", "numbat", dataset, sample,
                                                   "numbat", paste0(sample, "__numbat_percell.csv"))]
deg <- deg[file.exists(f)]
deg[, usable := mapply(usable_numbat, dataset, sample)]
degraded <- deg[usable == FALSE]
if (nrow(degraded)) {
  bad_vote <- merge(degraded[, .(dataset, sample)], S[, .(dataset, sample, rec_numbat)],
                    by = c("dataset", "sample"))[rec_numbat == TRUE]
  chk(nrow(bad_vote) == 0, "no degraded numbat file is counted as an arm",
      sprintf("%d counted", nrow(bad_vote)))
  cat(sprintf("  degraded numbat files: %d (correctly abstaining)\n", nrow(degraded)))
} else cat("  [skip] no degraded numbat files found\n")

cat("\n=========== 3.5 rollup covers exactly the per-sample records ===========\n")
if (file.exists(file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv"))) {
  A <- fread(file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv"))
  chk(nrow(A) == nrow(S), "rollup rows == per-sample summary files",
      sprintf("rollup %d vs files %d", nrow(A), nrow(S)))
  chk(!anyDuplicated(A[, .(dataset, sample)]), "no duplicated sample in the rollup")
  chk(all(A$malignant_frac >= 0 & A$malignant_frac <= 1, na.rm = TRUE), "malignant_frac in [0,1]")
  chk(all(A$n_malignant <= A$n_called, na.rm = TRUE), "n_malignant <= n_called")
} else cat("  [skip] rollup not built\n")

cat(sprintf("\n=========== BATCH 3: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 3 PASS\n")
