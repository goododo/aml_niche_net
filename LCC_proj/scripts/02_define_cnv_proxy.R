# 02_define_cnv_proxy.R ----
# Define, and then VALIDATE OR REJECT, the CNV-based TP53-aberrant proxy groups.
# Two independent tracks, each validated against the same negative control before being used:
#   CK track  : >= LCC_ARM_EVENT_MIN arms carrying an arm-level CNV event  (complex-karyotype-like)
#   17p track : >= LCC_ARM_FRAC_MIN of 17p genes in a loss state in a predominantly malignant subclone
#
# INPUT  : LCC_TAB_DIR/01_subclone_cnv.csv, 01_arm_events_long.csv, 01_parse_qc.csv
#          LCC_TAB_DIR/03_sample_manifest.csv via load_sample_meta()  -- timepoint-CORRECTED, and
#            the dataset selection read from the config constants, NOT from 06's output file
#            (06 reads THIS script's anchors, so depending on 06's file here would be circular)
#          TAB_DIR/05_ccc/ccc_node_features.csv          (downstream feasibility)
# OUTPUT : LCC_TAB_DIR/02_sample_cnv_proxy.csv    per-sample calls + group assignment
#          LCC_TAB_DIR/02_sensitivity_grid.csv    calls across the full threshold grid
#          LCC_TAB_DIR/02_arm_specificity.csv     per-arm Diagnosis-vs-Healthy hit rates
#          LCC_TAB_DIR/02_validation.csv          the pass/fail verdict for each track
# Usage  : Rscript LCC_proj/scripts/02_define_cnv_proxy.R
#
# THE NEGATIVE CONTROL IS THE POINT OF THIS SCRIPT. Healthy bone-marrow donors must not be called
# positive. The main line already documents that expression-CNV labelling misfires badly on some
# healthy samples (results/tables/02_malignancy/malignancy_fpr_healthy.csv: BM5-34p38n FPR 0.756,
# HSC_MPP bin FPR 0.398), so a proxy that is not explicitly held to a healthy FPR is worthless here.
#
# NO POSITIVE CONTROL EXISTS FOR THE 17p TRACK. The cohort's only genotype-confirmed TP53 sample
# (Petti2019 809653) carries TP53 E286G (missense) with del(17q) and no 17p loss, so it can validate
# the CK track only. Any 17p-track result is therefore specificity-validated but NOT sensitivity-
# validated -- state this wherever the group is used.

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--all_datasets", action = "store_true", default = FALSE,
              help = "ignore the dataset selection and score all 13 datasets (for the methods appendix)")
)))
# --all_datasets writes to *_all_datasets.csv so the selection-restricted tables are never clobbered.
# 07 reads the all-datasets version, because the genotype truth set lives mostly in GSE116256, which
# the selection excludes for its unusable healthy arm.
SFX <- if (opt$all_datasets) "_all_datasets" else ""
out_path <- function(stem) file.path(LCC_TAB_DIR, paste0(stem, SFX, ".csv"))

## -- Step 1. assemble ----
message("[1] loading 01 output + sample metadata")
sub <- fread(file.path(LCC_TAB_DIR, "01_subclone_cnv.csv"))
arm <- fread(file.path(LCC_TAB_DIR, "01_arm_events_long.csv"))
# load_sample_meta() applies the config_qc.R healthy rule, so Chen2023's 8 NBM normal-marrow
# samples are Healthy here even though the main-line manifest on disk stores them as Diagnosis.
# Never read the raw Timepoint column -- see config_lcc.R::load_sample_meta.
man <- load_sample_meta()
sub[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient, study = i.study,
              in_selection = i.in_selection), on = c("dataset", "sample")]
if (!opt$all_datasets) {
  n0 <- uniqueN(sub[, paste(dataset, sample)])
  sub <- sub[in_selection == TRUE]
  message("    restricted to the selected datasets: ",
          uniqueN(sub[, paste(dataset, sample)]), " of ", n0, " samples")
}

# subclone-level gates. n_cells_labelled (not n_cells): inferCNV ran on a larger pre-QC-refresh cell
# set, so only labelled cells carry malignancy information [see 01 header].
sub_ok <- sub[n_cells_labelled >= LCC_MIN_SUBCLUST_CELLS]
arm    <- merge(arm, sub_ok[, .(dataset, sample, subclone, is_reference, frac_malignant,
                                clone_frac_of_malignant, timepoint, uid_patient, study)],
                by = c("dataset", "sample", "subclone"))

samples <- unique(sub_ok[, .(dataset, sample, timepoint, uid_patient, study)])
message("    ", nrow(samples), " samples with at least one evaluable subclone")

## -- Step 2. per-sample internal noise floor from the inferCNV reference cells ----
message("[2] internal reference baseline (the per-sample noise floor)")
base <- sub_ok[is_reference == TRUE,
               .(n_ref_subclone = .N, n_ref_cells = sum(n_cells_labelled),
                 ref_arm_mean = sum(n_arm_events * n_cells_labelled) / sum(n_cells_labelled),
                 ref_arm_max  = max(n_arm_events)), by = .(dataset, sample)]
samples <- merge(samples, base, by = c("dataset", "sample"), all.x = TRUE)
samples[is.na(n_ref_subclone), `:=`(n_ref_subclone = 0L, n_ref_cells = 0L)]
samples[, has_baseline := n_ref_subclone >= LCC_MIN_REF_SUBCLONES & n_ref_cells >= LCC_MIN_REF_CELLS]

## -- Step 3. call both tracks over the full threshold grid ----
message("[3] sweeping the threshold grid (both tracks, with the healthy FPR at every cell)")
call_grid <- function(af, ae, fm) {
  ev  <- arm[is_event_at(frac_arm, af) & is_reference == FALSE & frac_malignant >= fm]
  ck  <- ev[, .(n_arms = uniqueN(arm)), by = .(dataset, sample, subclone)
            ][, .(ck_max_arms = max(n_arms)), by = .(dataset, sample)]
  p17 <- unique(ev[arm == LCC_TP53_ARM & direction == "loss", .(dataset, sample)])[, hit17p := TRUE]
  x <- merge(samples[, .(dataset, sample, timepoint)], ck, by = c("dataset", "sample"), all.x = TRUE)
  x <- merge(x, p17, by = c("dataset", "sample"), all.x = TRUE)
  x[is.na(ck_max_arms), ck_max_arms := 0L][is.na(hit17p), hit17p := FALSE]
  x[, ck_pos := ck_max_arms >= ae]
  nD <- sum(x$timepoint == "Diagnosis"); nH <- sum(x$timepoint == "Healthy")
  data.table(arm_frac_min = af, arm_event_min = ae, subclone_mal_min = fm,
             n_Dg = nD, n_Healthy = nH,
             ck_Dg = sum(x$ck_pos & x$timepoint == "Diagnosis"),
             ck_Healthy = sum(x$ck_pos & x$timepoint == "Healthy"),
             ck_healthy_fpr = if (nH) sum(x$ck_pos & x$timepoint == "Healthy") / nH else NA_real_,
             p17_Dg = sum(x$hit17p & x$timepoint == "Diagnosis"),
             p17_Healthy = sum(x$hit17p & x$timepoint == "Healthy"),
             p17_healthy_fpr = if (nH) sum(x$hit17p & x$timepoint == "Healthy") / nH else NA_real_,
             ck_auc_Dg_vs_H = auc_(x[timepoint == "Diagnosis"]$ck_max_arms,
                                   x[timepoint == "Healthy"]$ck_max_arms))
}
# helpers kept local: an "event" is a >= af gene-fraction call; AUC is the Mann-Whitney statistic.
is_event_at <- function(frac, af) frac >= af
auc_ <- function(a, b) if (!length(a) || !length(b)) NA_real_ else
  as.numeric(wilcox.test(a, b)$statistic) / (length(a) * length(b))

grid <- rbindlist(lapply(LCC_SENS_ARM_FRAC_MIN, function(af)
  rbindlist(lapply(LCC_SENS_ARM_EVENT_MIN, function(ae)
    rbindlist(lapply(LCC_SENS_SUBCLONE_MAL, function(fm) call_grid(af, ae, fm)))))))
fwrite_safe(grid, out_path("02_sensitivity_grid"))

## -- Step 4. lesion specificity: are calls concentrated on real AML lesions? ----
message("[4] per-arm Diagnosis-vs-Healthy specificity")
ev  <- arm[frac_arm >= LCC_ARM_FRAC_MIN & is_reference == FALSE]
hit <- unique(ev[, .(dataset, sample, arm, direction)])
hit <- merge(hit, samples[, .(dataset, sample, timepoint)], by = c("dataset", "sample"))
nD  <- sum(samples$timepoint == "Diagnosis"); nH <- sum(samples$timepoint == "Healthy")
spec <- dcast(hit[timepoint %in% c("Diagnosis", "Healthy"),
                  .(n = uniqueN(paste(dataset, sample))), by = .(arm, direction, timepoint)],
              arm + direction ~ timepoint, value.var = "n", fill = 0)
spec[, `:=`(pct_Diagnosis = round(100 * Diagnosis / nD, 1), pct_Healthy = round(100 * Healthy / nH, 1))]
spec[, `:=`(enrichment_pp = pct_Diagnosis - pct_Healthy,
            is_canonical  = paste0(arm, ":", direction) %in% LCC_CANONICAL_LESIONS)]
fwrite_safe(spec[order(-enrichment_pp)], out_path("02_arm_specificity"))

## -- Step 5. primary calls at the locked thresholds ----
message("[5] primary calls")
evp <- arm[frac_arm >= LCC_ARM_FRAC_MIN & is_reference == FALSE &
           frac_malignant >= LCC_SUBCLONE_MAL_MIN]
ck  <- evp[, .(n_arms = uniqueN(arm)), by = .(dataset, sample, subclone)
           ][, .(ck_max_arms = max(n_arms)), by = .(dataset, sample)]
p17 <- evp[arm == LCC_TP53_ARM & direction == "loss",
           .(n_subclone_17p = uniqueN(subclone), max_17p_frac = max(frac_arm),
             max_17p_clone_frac = max(clone_frac_of_malignant, na.rm = TRUE)), by = .(dataset, sample)]
out <- merge(samples, ck,  by = c("dataset", "sample"), all.x = TRUE)
out <- merge(out,     p17, by = c("dataset", "sample"), all.x = TRUE)
out[is.na(ck_max_arms), ck_max_arms := 0L][is.na(n_subclone_17p), n_subclone_17p := 0L]
out[, `:=`(ck_pos  = ck_max_arms >= LCC_ARM_EVENT_MIN,
           p17_pos = n_subclone_17p > 0)]
out[, cnv_proxy_class := fifelse(p17_pos & ck_pos, "CK_and_17p",
                          fifelse(p17_pos, "17p_only",
                           fifelse(ck_pos, "CK_only", "negative")))]
out[, proxy_confidence := fifelse(cnv_proxy_class %in% LCC_PROXY_HIGH, "high",
                           fifelse(cnv_proxy_class %in% LCC_PROXY_MID, "mid", "negative"))]
# primary analysis set: Diagnosis only, one sample per patient (largest evaluable sample wins)
out[, in_primary_set := timepoint == LCC_TIMEPOINT_PRIMARY]
setorder(out, uid_patient, -ck_max_arms)
out[in_primary_set == TRUE, dup_patient := duplicated(uid_patient)]
out[in_primary_set == TRUE & dup_patient == TRUE, in_primary_set := FALSE]
out[, dup_patient := NULL]
fwrite_safe(out, out_path("02_sample_cnv_proxy"))

## -- Step 6. verdict: does each track survive its own controls? ----
message("[6] validation verdicts")
g0 <- grid[arm_frac_min == LCC_ARM_FRAC_MIN & arm_event_min == LCC_ARM_EVENT_MIN &
           subclone_mal_min == LCC_SUBCLONE_MAL_MIN]
s809 <- out[sample == "809653"]
ver <- rbindlist(list(
  data.table(track = "CK (complex-karyotype-like)", metric = "healthy FPR",
             value = sprintf("%.3f (%d/%d)", g0$ck_healthy_fpr, g0$ck_Healthy, g0$n_Healthy),
             verdict = fifelse(g0$ck_healthy_fpr <= 0.05, "PASS", "FAIL")),
  data.table(track = "CK (complex-karyotype-like)", metric = "AUC Diagnosis vs Healthy (arm burden)",
             value = sprintf("%.3f", g0$ck_auc_Dg_vs_H),
             verdict = fifelse(g0$ck_auc_Dg_vs_H >= 0.70, "PASS", "FAIL")),
  data.table(track = "CK (complex-karyotype-like)", metric = "positive control 809653 called CK+",
             value = as.character(s809$ck_pos), verdict = fifelse(isTRUE(s809$ck_pos), "PASS", "FAIL")),
  data.table(track = "17p loss", metric = "healthy FPR",
             value = sprintf("%.3f (%d/%d)", g0$p17_healthy_fpr, g0$p17_Healthy, g0$n_Healthy),
             verdict = fifelse(g0$p17_healthy_fpr <= 0.05, "PASS", "FAIL")),
  data.table(track = "17p loss", metric = "positive control",
             value = "none exists (809653 is TP53 E286G missense + del(17q), no 17p loss)",
             verdict = "NOT TESTABLE"),
  data.table(track = "17p loss", metric = "primary-set positives (Diagnosis, 1/patient)",
             value = sprintf("%d positive vs %d negative",
                             out[in_primary_set == TRUE & p17_pos == TRUE, .N],
                             out[in_primary_set == TRUE & p17_pos == FALSE, .N]),
             verdict = fifelse(out[in_primary_set == TRUE & p17_pos == TRUE, .N] >= 10,
                               "PASS", "UNDERPOWERED"))))
fwrite_safe(ver, out_path("02_validation"))
print(ver)
message("[done] 02_define_cnv_proxy")
