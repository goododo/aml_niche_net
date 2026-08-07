#!/usr/bin/env Rscript
# =============================================================================
# 97_calibration_confound_check.R
#
# Recompute the calibration gate's two arms with the autologous reference cells
# taken OUT of the denominator, and report the gate both ways.
#
# WHY. 41_infercnv_to_percell sets malignant=0 on every reference cell, so a
# reference cell cannot be a false positive however it scores. Per sample that is
# correct -- those cells are the definition of normal for that sample. Across arms
# it is a confound, because the reference is drawn from mature T/B cells and the
# arms differ systematically in how many of those they contain:
#
#     Healthy   median frac_ref 0.31
#     Diagnosis median frac_ref 0.13   (blast-packed marrow)
#
# So a healthy sample has ~31% of its cells defined normal and a diagnosis sample
# ~13%. That pushes healthy frac_malignant DOWN relative to diagnosis -- i.e. it
# pushes the gate toward PASSING. The v1 gate failed despite this tailwind.
#
# This script does not change the gate. 70_residual_stratum keeps its own
# definition so its verdict stays comparable to v1. This is the control that says
# whether a pass is real or is the asymmetry.
#
# INPUT : DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv
#         REFNORM_REF_CELL_DIR/<ds>/<sample>_ref_norm_cells.txt
#         01_preprocess/02_sample_split.csv        (Timepoint)
# OUTPUT: DIR_MALIGNANCY/calibration_confound_check.csv
#
# Usage: Rscript scripts/02_malignancy/97_calibration_confound_check.R
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(here)})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

R <- qc_rds_roster(on_extra = "error")
sp <- fread(file.path(DIR_PREPROCESS, "02_sample_split.csv"))
if (!"sample" %in% names(sp) && "Sample" %in% names(sp)) setnames(sp, "Sample", "sample")
R <- merge(R, sp[, .(dataset, sample, Timepoint)], by = c("dataset", "sample"), all.x = TRUE)

per <- rbindlist(lapply(seq_len(nrow(R)), function(i) {
  ds <- R$dataset[i]; sid <- R$sample[i]
  cf <- file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
  if (!file.exists(cf)) return(NULL)
  d  <- fread(cf, select = c("cell", "malignant"))
  d  <- d[!is.na(malignant)]
  if (!nrow(d)) return(NULL)
  rf <- file.path(REFNORM_REF_CELL_DIR, ds, paste0(sid, "_ref_norm_cells.txt"))
  rc <- if (file.exists(rf)) readLines(rf) else character(0)
  d[, is_ref := cell %in% rc]
  nr <- d[is_ref == FALSE]
  data.table(dataset = ds, sample = sid, Timepoint = R$Timepoint[i],
             n_cells = nrow(d), n_ref = sum(d$is_ref), frac_ref = round(mean(d$is_ref), 4),
             n_mal = sum(d$malignant == 1),
             malignant_frac          = round(mean(d$malignant == 1), 4),
             malignant_frac_excl_ref = if (nrow(nr)) round(mean(nr$malignant == 1), 4) else NA_real_)
}), fill = TRUE)

if (!nrow(per)) stop("[97] no consensus per-cell files yet -- run the consensus stage first")
fwrite_safe(per, file.path(DIR_MALIGNANCY, "calibration_confound_check.csv"))

auc_mw <- function(pos, neg) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

cat("\n============ reference-cell share by arm ============\n")
print(per[, .(n = .N, med_frac_ref = round(median(frac_ref), 3)), by = Timepoint][order(-n)])

cat("\n============ the gate, computed both ways ============\n")
for (col in c("malignant_frac", "malignant_frac_excl_ref")) {
  h <- per[Timepoint == "Healthy"][[col]]; d <- per[Timepoint == "Diagnosis"][[col]]
  cat(sprintf("\n  %s\n", col))
  cat(sprintf("    healthy   n=%3d median %.4f\n", length(h), median(h, na.rm = TRUE)))
  cat(sprintf("    diagnosis n=%3d median %.4f\n", length(d), median(d, na.rm = TRUE)))
  cat(sprintf("    AUC = %.3f  (gate needs >= %.2f)\n", auc_mw(d, h), RESIDUAL_CALIB_MIN_AUC))
  cat(sprintf("    healthy median below the deep cut (%.2f)? %s\n",
              RESIDUAL_CALIB_MAX_HEALTHY_MEDIAN,
              ifelse(median(h, na.rm = TRUE) <= RESIDUAL_CALIB_MAX_HEALTHY_MEDIAN, "YES", "NO")))
}

cat("\n[reading it] The first block is what 70_residual_stratum gates on and what v1 is\n")
cat("  comparable to. The second removes the cells that were DEFINED normal, so it is the\n")
cat("  arm-comparable number. If the gate passes on the first and fails on the second, the\n")
cat("  pass is the reference-fraction asymmetry, not better CNV calling.\n")
cat("\n[done] ", file.path(DIR_MALIGNANCY, "calibration_confound_check.csv"), "\n")
