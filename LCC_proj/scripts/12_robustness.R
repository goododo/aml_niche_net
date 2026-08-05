#!/usr/bin/env Rscript
# 12_robustness.R ----
# Two robustness questions about the 11 result, answered together because they share machinery.
#
# Q1. AML328 contributes 4 of the 9 TP53-mut samples, so 11's "9 pairs" are really 6 patients and its
#     p-values are optimistic. -> DESIGN B collapses every patient to one value and re-tests at n = 6.
#
# Q2. "Try some other TP53-WT controls and see whether the pathway result comes alive."
#     Re-drawing control sets and keeping the draw that turns a null significant is p-hacking -- the
#     p-value stops meaning anything the moment the control set is chosen for its outcome. The
#     question BEHIND the request is legitimate and answerable: does the conclusion depend on WHICH
#     valid WT samples happened to be picked? So instead of trying a few sets and keeping the best,
#     this script enumerates R VALID matchings at random and reports the whole DISTRIBUTION of
#     p-values per feature. Robust findings are significant in most draws; artefacts of one lucky
#     match are significant in few. Nothing is selected on its result.
#
# READ THE OUTPUT THIS WAY: frac_p05 is the fraction of valid control sets in which the feature
# reached p < 0.05. It is the answer to "would another analyst, matching equally correctly, have
# seen this?" -- not a p-value, and it must not be reported as one.
#
# INPUT  : LCC_TAB_DIR/09_tp53_groups.csv, 04_detection_by_sample.csv, 04_pathway_sample.csv
# OUTPUT : LCC_TAB_DIR/12_robustness_genes.csv     per gene, both designs, over R draws
#          LCC_TAB_DIR/12_robustness_pathways.csv  per pathway, both designs, over R draws
#          LCC_TAB_DIR/12_design_B_pairs.csv       the patient-collapsed pairing actually used
# Usage  : Rscript LCC_proj/scripts/12_robustness.R [--draws 500]

suppressPackageStartupMessages({ library(data.table); library(here); library(optparse) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
set.seed(SEED)
opt <- parse_args(OptionParser(option_list = list(
  make_option("--draws", type = "integer", default = 500L, help = "random valid matchings [500]"))))
R <- opt$draws
MIN_CELLS <- 30L

## -- Step 1. arms and the per-sample value matrices ----
message("[1] loading")
grp <- fread(file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))
det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))
aml0 <- grp[timepoint != "Healthy" & !is.na(uid_patient) & uid_patient != "",
            .(dataset, sample, timepoint, uid_patient, tp53_tier)]
# Two mut definitions, run side by side so the cost of the looser one is visible rather than argued:
#   strict = tiers A+B only  (genotype or allele-LOH; 9 samples / 6 patients)
#   broad  = A+B+C          (adds 5 expression-CNV-only samples -> 14 samples / 11 patients)
# Tier C was measured at 0/2 sensitivity against genotype truth, so most of the added samples are
# expected to be mislabelled. Broad buys patients and spends specificity; the tables show the trade.
MUT_DEFS <- list(strict = c("A_genotype", "B_allele_loh"),
                 broad  = c("A_genotype", "B_allele_loh", "C_expr_only"))
mk_aml <- function(tiers) {
  a <- copy(aml0)
  a[, arm := fifelse(tp53_tier %in% tiers, "mut",
              fifelse(tp53_tier %in% c("WT_genotyped", "WT_presumed"), "WT", NA_character_))]
  a[!is.na(arm), .(dataset, sample, timepoint, uid_patient, arm)]
}
aml <- mk_aml(MUT_DEFS$strict)

alias <- c(PVRL4 = "NECTIN4", PVRL1 = "NECTIN1", PVRL2 = "NECTIN2", PVRL3 = "NECTIN3")
gv <- det[stratum == "all" & n_cells >= MIN_CELLS & gene_present == TRUE & !is.na(mean_lognorm)]
gv[gene %in% names(alias), gene := alias[gene]]
gv <- gv[, .(v = max(mean_lognorm)), by = .(dataset, sample, feature = gene)][, kind := "gene"]

pw <- fread(file.path(LCC_TAB_DIR, "04_pathway_sample.csv"))[unit == "sample" & group == "all"]
sc <- grep("_UCell$", names(pw), value = TRUE)
pv <- melt(pw, id.vars = c("dataset", "sample"), measure.vars = sc,
           variable.name = "feature", value.name = "v")[!is.na(v)]
pv[, `:=`(feature = as.character(feature), kind = "pathway")]
val <- rbind(gv, pv)
message("    ", uniqueN(val$feature), " features over ", uniqueN(val$sample), " samples")

## -- Step 2. the exact one-sided paired signed-rank, cached by k ----
# For tie-free |differences| the null depends only on k, so the DP is computed once per k.
null_cache <- new.env(parent = emptyenv())
sr_null <- function(r) {
  key <- paste(r, collapse = "_")
  if (!is.null(null_cache[[key]])) return(null_cache[[key]])
  k <- length(r); mx <- sum(r); w <- numeric(mx + 1L); w[1L] <- 1
  for (i in seq_len(k)) w <- w + c(numeric(r[i]), w[seq_len(mx + 1L - r[i])])
  cw <- rev(cumsum(rev(w))) / sum(w)         # cw[s+1] = P(W+ >= s)
  null_cache[[key]] <- cw; cw
}
p_paired <- function(dm, dw) {
  ok <- !is.na(dm) & !is.na(dw); d <- dm[ok] - dw[ok]; d <- d[d != 0]
  k <- length(d); if (k < 4L) return(NA_real_)
  r <- rank(abs(d)); if (any(r != round(r))) r <- rank(abs(d), ties.method = "first")
  cw <- sr_null(as.integer(r)); s <- as.integer(round(sum(r[d > 0])))
  cw[min(s + 1L, length(cw))]
}

## -- Step 3. matching machinery ----
# A draw picks, for each mut unit, a random WT unit from the same dataset (and timepoint, for the
# sample-level design), without reusing a WT sample or a WT patient-timepoint. Identical constraints
# to 11 -- only the tie-break changes from "nearest cell count" to random, which is what makes the
# draws a sample from the space of equally valid matchings.
draw_pairs <- function(mut_u, wt_u, strat_cols) {
  used_s <- character(0); used_pt <- character(0); out <- integer(nrow(mut_u))
  ord <- sample(seq_len(nrow(mut_u)))
  for (i in ord) {
    cond <- rep(TRUE, nrow(wt_u))
    for (cc in strat_cols) cond <- cond & wt_u[[cc]] == mut_u[[cc]][i]
    cond <- cond & !(wt_u$unit %in% used_s) & !(paste(wt_u$uid_patient, wt_u$timepoint) %in% used_pt)
    cand <- which(cond)
    if (!length(cand)) { out[i] <- NA_integer_; next }
    j <- cand[sample.int(length(cand), 1L)]
    used_s <- c(used_s, wt_u$unit[j]); used_pt <- c(used_pt, paste(wt_u$uid_patient[j], wt_u$timepoint[j]))
    out[i] <- j
  }
  out
}

run_design <- function(label, unit_col, strat_cols, mut_def = "strict") {
  message("[", label, "] building units")
  aml <- mk_aml(MUT_DEFS[[mut_def]])
  u <- unique(aml[, .(unit = get(unit_col), dataset, timepoint, uid_patient, arm)])
  # patient-level: a patient spanning timepoints gets its modal timepoint, and its value is the mean
  # over its samples -- one patient contributes exactly one number, which is the whole point.
  if (unit_col == "uid_patient")
    u <- u[, .(dataset = dataset[1], timepoint = names(sort(table(timepoint), decreasing = TRUE))[1],
               uid_patient = unit[1], arm = arm[1]), by = unit]
  v <- merge(val, aml[, .(dataset, sample, unit = get(unit_col))], by = c("dataset", "sample"))
  v <- v[, .(v = mean(v)), by = .(unit, feature, kind)]
  mut_u <- u[arm == "mut"]; wt_u <- u[arm == "WT"]
  message("    ", nrow(mut_u), " mut units vs ", nrow(wt_u), " WT candidates")

  feats <- unique(v[, .(feature, kind)])
  vw <- dcast(v, unit ~ feature, value.var = "v")
  setkey(vw, unit)
  mv <- as.matrix(vw[J(mut_u$unit), -1L])
  wv_all <- as.matrix(vw[J(wt_u$unit), -1L])
  fn <- colnames(mv)

  acc <- matrix(NA_real_, nrow = R, ncol = length(fn), dimnames = list(NULL, fn))
  first_pairs <- NULL
  for (r in seq_len(R)) {
    idx <- draw_pairs(mut_u, wt_u, strat_cols)
    if (all(is.na(idx))) next
    if (r == 1L) first_pairs <- data.table(mut_unit = mut_u$unit, wt_unit = wt_u$unit[idx],
                                           dataset = mut_u$dataset, timepoint = mut_u$timepoint)
    wv <- wv_all[idx, , drop = FALSE]
    acc[r, ] <- vapply(seq_along(fn), function(j) p_paired(mv[, j], wv[, j]), 0)
    if (r %% 100 == 0) message("    draw ", r, "/", R)
  }
  res <- data.table(feature = fn,
                    median_p = apply(acc, 2, median, na.rm = TRUE),
                    min_p    = apply(acc, 2, min,  na.rm = TRUE),
                    frac_p05 = colMeans(acc < 0.05, na.rm = TRUE),
                    frac_p01 = colMeans(acc < 0.01, na.rm = TRUE),
                    n_draws_evaluable = colSums(!is.na(acc)))
  res[feats, kind := i.kind, on = "feature"]
  res[, `:=`(design = label, mut_def = mut_def, n_mut_units = nrow(mut_u),
             n_wt_candidates = nrow(wt_u))]
  list(res = res, pairs = first_pairs)
}

## -- Step 4. run all four designs: {sample, patient} x {strict, broad} ----
RUNS <- list(
  list(l = "A_sample_strict",  u = "sample",      s = c("dataset", "timepoint"), m = "strict"),
  list(l = "B_patient_strict", u = "uid_patient", s = "dataset",                 m = "strict"),
  list(l = "C_sample_broad",   u = "sample",      s = c("dataset", "timepoint"), m = "broad"),
  list(l = "D_patient_broad",  u = "uid_patient", s = "dataset",                 m = "broad"))
got <- lapply(RUNS, function(r) run_design(r$l, r$u, r$s, r$m))
B <- got[[2]]
if (!is.null(B$pairs)) fwrite_safe(B$pairs, file.path(LCC_TAB_DIR, "12_design_B_pairs.csv"))
out <- rbindlist(lapply(got, `[[`, "res"), fill = TRUE)
fwrite_safe(out[kind == "gene"],    file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))
fwrite_safe(out[kind == "pathway"], file.path(LCC_TAB_DIR, "12_robustness_pathways.csv"))

message("[4] summary")
cat("\n== design sizes ==\n")
print(unique(out[, .(design, mut_def, n_mut_units, n_wt_candidates)]))
HITS <- c("SLC2A1", "GP9", "PGK1", "C1QB", "MAF", "LYVE1", "C1QC", "SPP1", "OSM", "PDGFB",
          "IL11", "PDGFRB", "C1QA", "MRC1", "TIGIT", "CD226", "NECTIN4")
cat("\n== the 11 hits across all four designs: frac of valid control sets reaching p<0.05 ==\n")
print(dcast(out[feature %in% HITS], feature ~ design, value.var = "frac_p05")[
      order(-A_sample_strict)])
cat("\n== same, median p ==\n")
print(dcast(out[feature %in% HITS], feature ~ design, value.var = "median_p")[
      order(A_sample_strict)])
for (dz in unique(out$design)) {
  cat("\n== PATHWAYS under", dz, "-- does ANY control set rescue them? ==\n")
  print(out[design == dz & kind == "pathway"
            ][order(-frac_p05), .(pathway = sub("_UCell$", "", feature),
                                  median_p = signif(median_p, 3), best_p = signif(min_p, 3),
                                  frac_p05 = round(frac_p05, 3))])
}
cat("\n== genes most robust to the control-set choice, per design ==\n")
for (dz in unique(out$design))
  print(out[design == dz & kind == "gene"][order(-frac_p05)][1:12,
        .(design, feature, median_p = signif(median_p, 3), frac_p05 = round(frac_p05, 3))])
message("[done] 12_robustness")
