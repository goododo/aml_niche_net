# 05_undercall_contamination.R ----
# INPUT  : DIR_CCC/ccc_node_features.csv                      (mean_stemness_normal, frac_malignant, n_cells)
#          DIR_CCC/ccc_sample_manifest.csv                    (healthy label, CCC eligibility)
#          DIR_MALIGNANCY/<ds>/<sample>__consensus_summary.csv (arms, malignant_frac)
#          results/tables/01_preprocess/00_curated_manifest.csv (blast_pct_clinical, karyotype)
# OUTPUT : DIR_CCC/undercall_contamination.csv                (per sample: under-call and stemness)
#          DIR_CCC/undercall_by_dataset.csv                   (the dataset-level table)
# WHAT   : Asks whether the elevated stemness of the NON-MALIGNANT pool in AML is really uncalled
#          blasts, using clinical blast percentage -- evidence the CNV caller never saw -- as the
#          external truth.
#
# WHY THIS EXISTS, AND WHY 04's ANSWER WAS NOT ENOUGH. 04_stemness_purity_sweep.R already asked
# whether the stemness signal survives removing likely-malignant cells from the "normal" pool. It
# purified by CNV BURDEN, which is derived from inferCNV -- so it can only remove cells inferCNV
# already half-sees. It is structurally unable to remove a blast that carries no copy-number
# change, and that is exactly the population at issue:
#
#   3853_Dg   46,XX,t(6;9)(p22;q34)/46,XX   balanced translocation, CN-neutral   clinical blast 85%
#   6323_Dg   46,XX                          normal karyotype                     clinical blast 50%
#
# inferCNV calls 0.000 malignant in both. That is the correct inferCNV answer -- it sees only copy
# number -- and numbat, which uses allelic imbalance, calls 0.765 and 0.497 there. So 04's sweep is
# circular for the case that matters, and a non-circular test needs a measurement that does not come
# from the CNV caller. Clinical blast percentage is that measurement.
#
# THE TEST. If elevated stemness in the non-malignant pool is uncalled leukaemia, it must RISE with
# how much leukaemia the caller missed. Under-call is (clinical blast fraction - called malignant
# fraction) per sample. A flat relationship is evidence the signal is not contamination; a positive
# one would retract the finding.
#
# READING IT:
#   rho ~ 0 pooled AND within every dataset  -> the finding is not uncalled blasts
#   rho > 0                                  -> the "normal" pool is contaminated; retract
#   rho < 0 strongly                         -> something else is confounded; do not celebrate it
#
# Usage : Rscript scripts/05_ccc/05_undercall_contamination.R [--min_blast 20] [--blind_frac 0.05]
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--min_blast",  type = "double", default = 20,
              help = "clinical blast %% at or above which a sample is unambiguously leukaemic"),
  make_option("--blind_frac", type = "double", default = 0.05,
              help = "called malignant fraction below which the caller is treated as blind")
)))
set.seed(SEED)

## -- Step 1. assemble one row per AML sample ----
CUR <- fread(file.path(DIR_PREPROCESS, "00_curated_manifest.csv"),
             select = c("sample", "blast_pct_clinical", "karyotype"))
man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))[ccc_eligible == TRUE,
             .(dataset, sample, Timepoint, sorted_sublib)]
man[, healthy := Timepoint == "Healthy"]

cons_files <- list.files(DIR_MALIGNANCY, pattern = "__consensus_summary\\.csv$",
                         recursive = TRUE, full.names = TRUE)
CON <- rbindlist(lapply(cons_files, function(p) fread(p)[1, .(sample, dataset, arms, malignant_frac)]),
                 fill = TRUE)

# Node features are per (sample, hierarchy_bin). Collapse to one number per sample, weighting each
# bin by its cell count -- an unweighted mean would let a 12-cell Stromal node count as much as a
# 2,000-cell Mono_DC node.
NF <- fread(file.path(DIR_CCC, "ccc_node_features.csv"))[has_graph == TRUE]
P <- NF[, .(stem_normal = weighted.mean(mean_stemness_normal, n_cells, na.rm = TRUE),
            frac_mal_nf = weighted.mean(frac_malignant,       n_cells, na.rm = TRUE),
            cells       = sum(n_cells, na.rm = TRUE)), by = .(dataset, sample)]

D <- Reduce(function(a, b) merge(a, b, by = intersect(names(a), names(b))),
            list(man, P, CON[, .(dataset, sample, arms, malignant_frac)]))
D <- merge(D, CUR, by = "sample", all.x = TRUE)
message(sprintf("[1] CCC-eligible samples with node features and a consensus record: %d (healthy %d / AML %d)",
                nrow(D), sum(D$healthy), sum(!D$healthy)))

A <- D[healthy == FALSE & is.finite(malignant_frac)]
A[, has_second_arm := grepl("numbat|author", arms)]
A[, blast_frac := blast_pct_clinical / 100]
A[, undercall  := blast_frac - malignant_frac]
B <- A[is.finite(undercall) & is.finite(stem_normal)]
stopifnot(nrow(B) > 0)
message(sprintf("[1] AML samples carrying a clinical blast %%: %d of %d", nrow(B), nrow(A)))
if (nrow(B) < 30)
  stop("fewer than 30 AML samples carry a clinical blast percentage -- every correlation below ",
       "would be noise. Do not report this script's verdict; fix the metadata first.")

## -- Step 2. does the CNV caller track leukaemia burden at all? ----
message("\n[2] does the called malignant fraction track CLINICAL blast burden?")
rho_of <- function(X, lbl) {
  X <- X[is.finite(malignant_frac) & is.finite(blast_frac)]
  if (nrow(X) < 8) { message(sprintf("    %-38s n=%-3d too few", lbl, nrow(X))); return(NULL) }
  ct <- suppressWarnings(cor.test(X$malignant_frac, X$blast_frac, method = "spearman"))
  message(sprintf("    %-38s n=%-3d rho=%+.3f p=%.3f | med blast %.0f%% vs med called %.3f",
                  lbl, nrow(X), ct$estimate, ct$p.value, 100 * median(X$blast_frac), median(X$malignant_frac)))
  data.table(stratum = lbl, n = nrow(X), rho = as.numeric(ct$estimate), p = ct$p.value)
}
strata <- rbindlist(list(
  rho_of(B,                                             "all AML"),
  rho_of(B[has_second_arm == FALSE],                    "  inferCNV arm only"),
  rho_of(B[has_second_arm == FALSE & sorted_sublib == FALSE], "  ...unsorted")
), fill = TRUE)

message("\n    per dataset (the between-dataset spread is the point):")
byds <- B[, .(n = .N, med_blast = round(100 * median(blast_frac), 1),
              med_called = round(median(malignant_frac), 3)), by = dataset][order(-n)]
for (i in seq_len(nrow(byds)))
  message(sprintf("      %-14s n=%-3d clinical blast %4.1f%%  ->  called %.3f",
                  byds$dataset[i], byds$n[i], byds$med_blast[i], byds$med_called[i]))

## -- Step 3. how many samples is the caller simply blind on ----
B[, blind := malignant_frac < opt$blind_frac & blast_pct_clinical >= opt$min_blast]
message(sprintf("\n[3] blind samples (called <%.2f while clinical blast >=%.0f%%): %d of %d (%.0f%%)",
                opt$blind_frac, opt$min_blast, sum(B$blind), nrow(B), 100 * mean(B$blind)))
message(sprintf("    of those, rescued by a second evidence arm: %d",
                B[blind == TRUE & has_second_arm == TRUE, .N]))
# NON-VACUITY: if nothing is blind the test below has no range to work with, and a flat correlation
# would mean nothing at all.
if (sum(B$blind) == 0)
  stop("no sample is blind by these thresholds -- the contamination test has no dynamic range. ",
       "Either the thresholds are wrong or the caller improved; check before reporting a flat result.")

## -- Step 4. THE TEST ----
message("\n[4] CONTAMINATION TEST -- does normal-pool stemness rise with what the caller missed?")
ct <- suppressWarnings(cor.test(B$stem_normal, B$undercall, method = "spearman"))
message(sprintf("    pooled        rho = %+.3f   p = %.4f   (n=%d)", ct$estimate, ct$p.value, nrow(B)))
per_ds <- rbindlist(lapply(B[, .N, by = dataset][N >= 8]$dataset, function(d) {
  X <- B[dataset == d]
  c2 <- suppressWarnings(cor.test(X$stem_normal, X$undercall, method = "spearman"))
  message(sprintf("    %-13s rho = %+.3f   p = %.4f   (n=%d)", d, c2$estimate, c2$p.value, nrow(X)))
  data.table(dataset = d, n = nrow(X), rho = as.numeric(c2$estimate), p = c2$p.value)
}), fill = TRUE)

# Tertiles say the same thing without assuming monotonicity.
B[, tert := cut(undercall, quantile(undercall, 0:3/3), include.lowest = TRUE,
                labels = c("least missed", "middle", "most missed"))]
message("\n    by under-call tertile:")
TT <- B[, .(n = .N, med_undercall = round(median(undercall), 3),
            med_stem_normal = round(median(stem_normal), 4)), by = tert][order(tert)]
for (i in seq_len(nrow(TT)))
  message(sprintf("      %-13s n=%-3d  under-call %.3f  ->  stemness %.4f",
                  as.character(TT$tert[i]), TT$n[i], TT$med_undercall[i], TT$med_stem_normal[i]))

## -- Step 5. verdict ----
pos_any <- ct$p.value < 0.05 && ct$estimate > 0
pos_ds  <- nrow(per_ds) && any(per_ds$p < 0.05 & per_ds$rho > 0)
message("\n[5] VERDICT")
if (pos_any || pos_ds) {
  message("    normal-pool stemness RISES with what the caller missed.")
  message("    -> the finding is consistent with uncalled blasts in the normal pool. RETRACT it")
  message("       until the malignancy call is fixed; do not report mean_stemness_normal.")
} else {
  message(sprintf("    no relationship, pooled or within any dataset (max |rho| = %.3f).",
                  max(abs(c(ct$estimate, per_ds$rho)), na.rm = TRUE)))
  message("    -> the elevated stemness of the non-malignant pool is NOT explained by uncalled")
  message("       blasts, over a range of under-call this cohort actually spans.")
  message("    This supersedes 04's CNV-burden purity sweep as the argument: that sweep purified")
  message("    using the very signal inferCNV cannot see, and was circular for CN-neutral disease.")
}
message("\n    NOTE: this says nothing about whether the malignant/normal SPLIT is usable as a")
message("    feature. It is not -- see the dataset table in [2]. It says only that the stemness")
message("    contrast survives the split being wrong.")

fwrite_safe(B[, .(dataset, sample, Timepoint, arms, has_second_arm, cells,
                  blast_pct_clinical, malignant_frac, undercall, stem_normal, blind, karyotype)],
            file.path(DIR_CCC, "undercall_contamination.csv"))
fwrite_safe(byds, file.path(DIR_CCC, "undercall_by_dataset.csv"))
message("\n[done] wrote ", file.path(DIR_CCC, "undercall_contamination.csv"))
