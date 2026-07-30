# =============================================================================
# 70_residual_stratum.R   (stage 03.5, runs after 50_consensus_malignancy.R)
#
# Derives the H3 treatment-pressure stratum from RESIDUAL DISEASE BURDEN rather
# than from the depositing author's timepoint label.
#
# Why not the label: "Complete_remission" (GSE201966) is clinically an MRD draw,
# and GSE116256 maps every D<n> to Post_treatment whether n is 14 or 171. Only
# GSE227903 ever reaches the MRD class, because its authors happened to use the
# word. See the note above timepoint_days() in config_qc.R.
#
# Scope: ONLY samples whose nominal Timepoint is in RESIDUAL_NOMINAL_SCOPE.
# Diagnosis / Relapse / Healthy keep their nominal label untouched -- a low-blast
# diagnosis marrow is not an MRD sample.
#
# BLOCKING CALIBRATION GATE: frac_malignant is only a residual-disease proxy if
# it separates healthy marrow from AML marrow. This script REFUSES to write the
# stratum when it does not, and writes a diagnostic instead. That is the whole
# point -- a stratum built on a non-discriminating score would hand H3 a result
# manufactured from caller noise.
#
# CIRCULARITY WARNING (carried into every downstream H3 report):
#   frac_malignant is simultaneously (a) the variable defining this stratum and
#   (b) the covariate R5 requires be regressed out of the topology scores. H3
#   must therefore be reported on BOTH axes -- nominal label and residual
#   stratum -- and its primary criterion is that the topology signal survives
#   continuous frac_malignant as a covariate WITHIN stratum. Otherwise H3
#   degrades to "samples with more blasts have more blasts".
#
# INPUT  : DIR_MALIGNANCY/ALL_consensus_summary.csv   (sample, dataset, malignant_frac, evidence_tier)
#          DIR_PREPROCESS/02_sample_split.csv         (dataset, Sample, Timepoint, Patient_ID)
#          DIR_MALIGNANCY/malignancy_fpr_healthy.csv  (optional; per-sample healthy FPR)
# OUTPUT : DIR_MALIGNANCY/residual_stratum.csv        (one row per sample)
#          DIR_MALIGNANCY/residual_stratum_sweep.csv  (threshold sensitivity)
#          DIR_MALIGNANCY/residual_calibration.csv    (the gate's evidence, always written)
# Usage  : Rscript scripts/02_malignancy/70_residual_stratum.R [--force]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--force", action = "store_true", default = FALSE,
              help = "write the stratum even if the calibration gate fails (diagnostic use only)")
)))

f_cons  <- file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv")
f_split <- file.path(DIR_PREPROCESS, "02_sample_split.csv")
f_fpr   <- file.path(DIR_MALIGNANCY, "malignancy_fpr_healthy.csv")
for (f in c(f_cons, f_split)) if (!file.exists(f)) stop("[70] missing input: ", f)

cons  <- fread(f_cons)
split <- fread(f_split)
if (!"sample" %in% names(split) && "Sample" %in% names(split)) setnames(split, "Sample", "sample")

keep_cols <- intersect(c("dataset","sample","Patient_ID","uid_patient","Timepoint",
                         "Timepoint_days","subtype_stratum","study_role"), names(split))
d <- merge(cons[, .(dataset, sample, malignant_frac, evidence_tier, conflict_frac,
                    n_called, n_qc_cells)],
           split[, ..keep_cols], by = c("dataset","sample"), all.x = TRUE)

n_unmatched <- d[is.na(Timepoint), .N]
if (n_unmatched > 0)
  message("[70] WARNING: ", n_unmatched, " consensus rows have no split-table match (Timepoint NA)")

# ---------------------------------------------------------------------------
# [1] Calibration gate -- does frac_malignant separate healthy from AML at all?
# ---------------------------------------------------------------------------
# AUC via the Mann-Whitney U identity (no pROC dependency): the probability that
# a random diagnosis sample scores above a random healthy one.
auc_mw <- function(pos, neg) {
  pos <- pos[is.finite(pos)]; neg <- neg[is.finite(neg)]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

healthy_v <- d[Timepoint == "Healthy",   malignant_frac]
dx_v      <- d[Timepoint == "Diagnosis", malignant_frac]

healthy_median <- if (length(healthy_v)) median(healthy_v, na.rm = TRUE) else NA_real_
dx_median      <- if (length(dx_v))      median(dx_v,      na.rm = TRUE) else NA_real_
auc            <- auc_mw(dx_v, healthy_v)

healthy_fpr_median <- NA_real_
if (file.exists(f_fpr)) {
  fpr <- fread(f_fpr)
  fpr_col <- if ("fpr_ccc_bins" %in% names(fpr)) "fpr_ccc_bins" else "fpr"
  healthy_fpr_median <- median(fpr[[fpr_col]], na.rm = TRUE)
}

# Fraction of samples resting on a single evidence arm -- reported, not gated on,
# because gating there would currently exclude ~96% of the cohort.
tier_c_frac <- d[, mean(grepl("^C", evidence_tier), na.rm = TRUE)]

checks <- data.table(
  check = c("healthy_median_below_deep_cut", "healthy_vs_diagnosis_auc", "healthy_fpr_median"),
  value = c(healthy_median, auc, healthy_fpr_median),
  bound = c(RESIDUAL_CALIB_MAX_HEALTHY_MEDIAN, RESIDUAL_CALIB_MIN_AUC,
            RESIDUAL_CALIB_MAX_HEALTHY_FPR),
  rule  = c("value <= bound", "value >= bound", "value <= bound"))
checks[, pass := c(!is.na(value[1]) && value[1] <= bound[1],
                   !is.na(value[2]) && value[2] >= bound[2],
                   is.na(value[3])  || value[3] <= bound[3])]

calib <- rbindlist(list(
  checks,
  data.table(check = c("n_healthy","n_diagnosis","diagnosis_median","frac_single_arm_evidence"),
             value = c(length(healthy_v), length(dx_v), dx_median, tier_c_frac),
             bound = NA_real_, rule = "reported only", pass = NA)),
  use.names = TRUE, fill = TRUE)

fwrite_safe(calib, file.path(DIR_MALIGNANCY, "residual_calibration.csv"))

message("[1] calibration:")
message(sprintf("      healthy   n=%3d median frac_malignant = %s", length(healthy_v), format(healthy_median, digits = 3)))
message(sprintf("      diagnosis n=%3d median frac_malignant = %s", length(dx_v),      format(dx_median,      digits = 3)))
message(sprintf("      healthy-vs-diagnosis AUC = %s (need >= %.2f)", format(auc, digits = 3), RESIDUAL_CALIB_MIN_AUC))
message(sprintf("      median healthy FPR       = %s (need <= %.2f)", format(healthy_fpr_median, digits = 3), RESIDUAL_CALIB_MAX_HEALTHY_FPR))
message(sprintf("      single-arm evidence      = %.1f%% of samples", 100 * tier_c_frac))

gate_ok <- all(checks$pass)
if (!gate_ok) {
  message("[1] CALIBRATION FAILED on: ", paste(checks[pass == FALSE, check], collapse = ", "))
  if (!opt$force) {
    message("[70] REFUSING to write residual_stratum.csv. frac_malignant does not separate ",
            "healthy from AML marrow, so any stratum built on it would encode caller noise, ",
            "not treatment pressure. Fix the malignancy consensus (or add a second evidence ",
            "arm) and re-run. Evidence written to residual_calibration.csv.")
    quit(status = 0)
  }
  message("[70] --force given: writing anyway. Results are NOT usable for H3.")
}

# ---------------------------------------------------------------------------
# [2] Assign the stratum (post-treatment draws only)
# ---------------------------------------------------------------------------
cut_stratum <- function(x, deep, resid) {
  fifelse(is.na(x), NA_character_,
  fifelse(x <  deep,  RESIDUAL_LABELS[1],
  fifelse(x <  resid, RESIDUAL_LABELS[2], RESIDUAL_LABELS[3])))
}

d[, in_scope := Timepoint %in% RESIDUAL_NOMINAL_SCOPE]
d[, residual_stratum := NA_character_]
d[in_scope == TRUE,
  residual_stratum := cut_stratum(malignant_frac, RESIDUAL_CUTS[["deep"]], RESIDUAL_CUTS[["residual"]])]

# tx_stage: one column H3 can use directly -- nominal label everywhere except the
# post-treatment draws, which carry their burden-derived class instead.
d[, tx_stage := fifelse(in_scope == TRUE & !is.na(residual_stratum), residual_stratum, Timepoint)]

# Confidence is carried, never used to silently drop a sample.
d[, stratum_confidence := fifelse(grepl("^A", evidence_tier), "high",
                          fifelse(grepl("^B", evidence_tier), "medium", "low"))]
d[, calibration_passed := gate_ok]

out_cols <- intersect(c("dataset","sample","Patient_ID","uid_patient","Timepoint","Timepoint_days",
                        "malignant_frac","evidence_tier","conflict_frac","n_called","n_qc_cells",
                        "in_scope","residual_stratum","tx_stage","stratum_confidence",
                        "calibration_passed","subtype_stratum","study_role"), names(d))
setorder(d, dataset, sample)
fwrite_safe(d[, ..out_cols], file.path(DIR_MALIGNANCY, "residual_stratum.csv"))

message("[2] stratum assigned for ", d[in_scope == TRUE, .N], " post-treatment samples:")
print(d[in_scope == TRUE, .N, by = .(Timepoint, residual_stratum)][order(Timepoint, residual_stratum)])

# ---------------------------------------------------------------------------
# [3] Threshold sensitivity sweep
# ---------------------------------------------------------------------------
# The 5%/20% cuts are a clinical convention, not a property of this cohort. Any
# H3 claim has to be shown to survive moving them, so the full grid is written
# out rather than left as a footnote.
sweep <- rbindlist(lapply(RESIDUAL_SWEEP_DEEP, function(dp) {
  rbindlist(lapply(RESIDUAL_SWEEP_RESIDUAL, function(rs) {
    if (rs <= dp) return(NULL)
    s <- cut_stratum(d[in_scope == TRUE, malignant_frac], dp, rs)
    data.table(cut_deep = dp, cut_residual = rs,
               n_deep     = sum(s == RESIDUAL_LABELS[1], na.rm = TRUE),
               n_residual = sum(s == RESIDUAL_LABELS[2], na.rm = TRUE),
               n_active   = sum(s == RESIDUAL_LABELS[3], na.rm = TRUE),
               n_na       = sum(is.na(s)))
  }))
}))
fwrite_safe(sweep, file.path(DIR_MALIGNANCY, "residual_stratum_sweep.csv"))

message("[3] wrote residual_stratum.csv / residual_stratum_sweep.csv / residual_calibration.csv")
message("[done] 70_residual_stratum")
