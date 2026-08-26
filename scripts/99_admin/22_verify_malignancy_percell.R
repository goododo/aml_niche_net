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

# THE ROSTER COMES FROM THE COHORT, NOT FROM THE TASK FILES. This used to glob
# ^infercnv_tasks(_healthy)?\.tsv$, which matches 2 of the 6 task files now in the tree. The
# 33-sample external-reference backfill was queued through infercnv_tasks_sorted.tsv, so it was
# INVISIBLE here: 3.5 failed loudly (rollup 212 vs files 179) but 3.1-3.4 had been quietly
# validating a 179-sample subset and reporting PASS. A roster that has to be extended by hand
# every time work is queued will silently shrink the audit exactly when new work needs auditing.
# qc_rds_roster() is the roster every producing script already uses, so it cannot drift from them.
tasks <- qc_rds_roster(on_extra = "ignore")[, .(dataset, sample)]
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

cat("\n=========== 3.4b where the CNV caller is blind, and whether anything rescues it ===========\n")
# THIS CHECK WAS WRONG WHEN FIRST WRITTEN (2026-08-25) AND IS INVERTED HERE.
#
# It originally FAILED on 3853_Dg and 6323_Dg because inferCNV called 0.000 malignant there while
# the consensus called 0.765 / 0.497, and read that as numbat overriding the primary caller. The
# karyotypes say otherwise:
#
#   3853_Dg   46,XX,t(6;9)(p22;q34)/46,XX   balanced translocation, COPY-NUMBER NEUTRAL   blast 85%
#   6323_Dg   46,XX                          normal karyotype                              blast 50%
#
# Both marrows are full of leukaemia and neither carries a copy-number change. inferCNV can only
# see copy number, so 0.000 is the CORRECT inferCNV answer -- it is blind to this disease, not
# wrong about it. numbat uses allelic imbalance and can see it. The consensus following numbat is
# the right outcome, and the old check punished it.
#
# The real failure mode is the opposite one: a sample where the caller sees nothing and NOTHING
# rescues it. Measured over the 59 AML samples with a recorded clinical blast percentage:
# inferCNV malignant_frac has NO relationship to clinical blast burden (Spearman rho = +0.093,
# p = 0.52 restricted to unsorted diagnosis samples; median clinical blast 65.5% vs median called
# 0.072, a ~9x under-call), and 23 of 59 have < 5% called against >= 20% clinical blasts. None of
# the 5 samples carrying a second arm is in that group.
#
# The baseline below counts the FULL roster, not the CCC-eligible subset: 39 of the 88 samples that
# carry a clinical blast percentage. The 23/59 figure above is the same measurement restricted to
# the 138-sample analysis cohort, which is the number that matters for the results.
#
# Under-detection on this scale does NOT explain the stemness finding, and that was tested without
# reusing inferCNV: if elevated stemness in the non-malignant pool were uncalled blasts, it would
# rise with how much leukaemia the caller missed. It does not -- Spearman rho = -0.036, p = 0.79
# over 59 samples, flat in every dataset separately, and flat across under-call tertiles (median
# stemness 0.0145 / 0.0181 / 0.0165 while median under-call goes 0.006 -> 0.465 -> 0.743). That
# test is worth more than the earlier CNV-burden purity sweep, which purified using the very
# signal inferCNV is blind to and was therefore circular.
BLIND_BASELINE <- 39L   # recorded debt, not an approval -- this must not grow silently
if (nrow(S)) {
  bl <- merge(S[, .(dataset, sample, exp_numbat, exp_author, cons)],
              fread("results/tables/01_preprocess/00_curated_manifest.csv")[, .(sample, blast_pct_clinical)],
              by = "sample")
  bl <- bl[!is.na(blast_pct_clinical)]
  bl[, mfrac := vapply(cons, function(p) { d <- fread(p); as.numeric(d$malignant_frac[1]) }, numeric(1))]
  bl[, rescued := exp_numbat | exp_author]
  bl[, blind := is.finite(mfrac) & mfrac < 0.05 & blast_pct_clinical >= 20]
  cat(sprintf("  samples with a clinical blast %%: %d | blind: %d | of those, rescued by a 2nd arm: %d\n",
              nrow(bl), sum(bl$blind), bl[blind == TRUE & rescued == TRUE, .N]))
  chk(nrow(bl) >= 40, "the blind-caller check reaches enough samples to mean anything",
      sprintf("only %d samples carry a clinical blast %%", nrow(bl)))
  chk(sum(bl$blind) <= BLIND_BASELINE,
      sprintf("the blind-caller count has not grown past its recorded %d", BLIND_BASELINE),
      sprintf("now %d -- a new sample lost its malignancy call; find out which caller stopped seeing it",
              sum(bl$blind)))
  chk(bl[blind == TRUE & rescued == TRUE, .N] == 0,
      "no sample with a second evidence arm is left blind",
      sprintf("%d rescued sample(s) still blind", bl[blind == TRUE & rescued == TRUE, .N]))
} else cat("  [skip] no consensus summaries\n")

cat("\n=========== 3.5 rollup covers exactly the per-sample records ===========\n")
# Counted straight off the filesystem as well as through the roster. If the two ever disagree the
# roster is wrong, and the roster is what the rest of this batch trusts -- so a roster-only count
# would report PASS on a subset, which is the bug this section was written to catch.
disk <- list.files(DIR_MALIGNANCY, pattern = "__consensus_summary\\.csv$", recursive = TRUE)
cat(sprintf("  roster samples: %d | with a consensus summary: %d | summaries on disk: %d\n",
            nrow(tasks), nrow(S), length(disk)))
chk(nrow(S) == length(disk), "roster sees every consensus summary on disk",
    sprintf("roster %d vs disk %d -- an orphan output or a roster gap", nrow(S), length(disk)))
if (file.exists(file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv"))) {
  A <- fread(file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv"))
  chk(nrow(A) == nrow(S), "rollup rows == per-sample summary files",
      sprintf("rollup %d vs files %d", nrow(A), nrow(S)))
  chk(!anyDuplicated(A[, .(dataset, sample)]), "no duplicated sample in the rollup")
  chk(all(A$malignant_frac >= 0 & A$malignant_frac <= 1, na.rm = TRUE), "malignant_frac in [0,1]")
  chk(all(A$n_malignant <= A$n_called, na.rm = TRUE), "n_malignant <= n_called")

  ## -- 3.6 the coverage gate on malignant_frac ------------------------------------------------
  # malignant_frac divides by n_called, not by n_qc_cells. When a sample's only two arms are both
  # EXPRESSION arms, strict_majority() returns NA on every cell they disagree about, so n_called
  # collapses and the published fraction can describe a fraction of a percent of the sample --
  # measured at 8 of 3,910 cells for GSE239721/PT15A. These assert the collapse stays visible.
  if (all(c("called_frac", "expr_conflict_frac", "evidence_tier") %in% names(A))) {
    chk(all(abs(A$called_frac - A$n_called / A$n_qc_cells) < 1e-4, na.rm = TRUE),
        "called_frac equals n_called / n_qc_cells",
        sprintf("%d rows disagree", sum(abs(A$called_frac - A$n_called / A$n_qc_cells) >= 1e-4, na.rm = TRUE)))
    chk(all((A$evidence_tier == "D_low_coverage") == (A$called_frac < CONSENSUS_MIN_CALLED_FRAC)),
        sprintf("D_low_coverage iff called_frac < %.2f", CONSENSUS_MIN_CALLED_FRAC),
        sprintf("%d rows violate", sum((A$evidence_tier == "D_low_coverage") != (A$called_frac < CONSENSUS_MIN_CALLED_FRAC))))
    # NON-VACUITY. Every assertion above is satisfied by an empty set. If no sample is ever tiered
    # D, the gate is untested and indistinguishable from one that never fires.
    chk(A[evidence_tier == "D_low_coverage", .N] > 0,
        "at least one sample actually trips the coverage gate (the check is not vacuous)",
        sprintf("%d tiered D", A[evidence_tier == "D_low_coverage", .N]))
    # expr_conflict_frac must be populated exactly where two expression arms exist, and NA elsewhere
    two_expr <- A[, grepl("infercnv", arms) & grepl("author", arms)]
    chk(all(!is.na(A$expr_conflict_frac[two_expr])) && all(is.na(A$expr_conflict_frac[!two_expr])),
        "expr_conflict_frac is present iff the sample has two expression arms",
        sprintf("%d with two expr arms, %d non-NA", sum(two_expr), sum(!is.na(A$expr_conflict_frac))))
    cat(sprintf("  coverage: %d samples below the gate, worst called_frac %.4f, max arm disagreement %.4f\n",
                A[evidence_tier == "D_low_coverage", .N], min(A$called_frac),
                max(A$expr_conflict_frac, na.rm = TRUE)))
  } else cat("  [SKIP] summary predates called_frac/expr_conflict_frac -- rerun 50 and 60\n")
} else cat("  [skip] rollup not built\n")

cat(sprintf("\n=========== BATCH 3: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 3 PASS\n")
