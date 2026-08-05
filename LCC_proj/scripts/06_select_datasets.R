# 06_select_datasets.R ----
# Build the per-dataset quality scorecard, apply the locked dataset selection, and emit the analysis
# manifest that every downstream LCC script reads. This is the EVIDENCE behind the dataset list in
# config_lcc.R (LCC_CORE_DATASETS / LCC_STROMA_REF_DATASETS / LCC_EXCLUDED_DATASETS).
#
# SELECTION PRINCIPLE: keep only datasets that contain BOTH AML-diagnosis and healthy-donor samples,
# so every contrast is WITHIN study. That removes the cross-study confounder by design instead of
# modelling it away -- and cross-study variance is what killed the complex-karyotype proxy (see the
# within-study re-test in Step 4 below, which shows the CK failure is NOT a cross-study artefact).
#
# NUMBERING NOTE: this script is numbered 06 because it was written after 01-05, but it sits UPSTREAM
# of them in the data flow -- 02/04/05 read its manifest. Run order is 01 -> 03 -> 04 -> 06 -> 02 -> 05.
# (The main line has the same pattern; CLAUDE.md: "08_scoring numbering is narrative order, not a
# data dependency".)
#
# INPUT  : LCC_TAB_DIR/03_sample_manifest.csv          (via load_sample_meta(), timepoint-corrected)
#          LCC_TAB_DIR/01_parse_qc.csv, 01_subclone_cnv.csv, 01_arm_events_long.csv
#          LCC_TAB_DIR/04_nectin_gate.csv, 04_detection_by_sample.csv
#          LCC_TAB_DIR/02_sample_cnv_proxy.csv         (17p anchors; optional, skipped if absent)
#          DIR_MALIGNANCY/malignancy_fpr_healthy.csv   (main-line negative control)
#          TAB_DIR/05_ccc/ccc_node_features.csv
# OUTPUT : LCC_TAB_DIR/06_dataset_scorecard.csv        every metric, every dataset
#          LCC_TAB_DIR/06_selection_rationale.csv      keep/drop + the disqualifying metric
#          LCC_TAB_DIR/06_analysis_manifest.csv        THE analysis sample set, timepoint-corrected
#          LCC_TAB_DIR/06_ck_within_study.csv          within-study CK re-test (methods evidence)
# Usage  : Rscript LCC_proj/scripts/06_select_datasets.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

## -- Step 1. corrected sample metadata ----
message("[1] loading timepoint-corrected sample metadata")
man <- load_sample_meta(include_stroma_ref = TRUE)
message("    ", nrow(man), " samples; timepoint corrected for ", sum(man$timepoint_was_fixed),
        " (Chen2023 NBM normal-marrow donors stored as Diagnosis on disk)")

## -- Step 2. per-dataset metrics ----
message("[2] assembling the scorecard")
sc <- man[, .(n_samples = .N,
              aml_dg_samples  = sum(timepoint == "Diagnosis"),
              aml_dg_patients = uniqueN(uid_patient[timepoint == "Diagnosis"]),
              healthy_samples = sum(timepoint == "Healthy"),
              healthy_patients= uniqueN(uid_patient[timepoint == "Healthy"]),
              other_samples   = sum(!timepoint %in% c("Diagnosis", "Healthy"))), by = dataset]
sc[, within_study_contrast := aml_dg_patients > 0 & healthy_patients > 0]

join <- function(x, y) merge(x, y, by = "dataset", all.x = TRUE)

# 2a. inferCNV integrity: are the malignancy labels aligned to the cells inferCNV actually saw,
#     and does each sample have an internal reference_normal baseline to be judged against?
qc <- fread(file.path(LCC_TAB_DIR, "01_parse_qc.csv"))
sc <- join(sc, qc[status == "ok", .(icnv_labelled_med = round(median(icnv_labelled_frac, na.rm = TRUE), 2),
                                    n_icnv_lowlabel  = sum(icnv_labelled_frac < 0.90, na.rm = TRUE),
                                    n_icnv_samples   = .N), by = dataset])
sub <- fread(file.path(LCC_TAB_DIR, "01_subclone_cnv.csv"))[n_cells_labelled >= LCC_MIN_SUBCLUST_CELLS]
base <- sub[is_reference == TRUE, .(n_ref_sub = .N, n_ref_cells = sum(n_cells_labelled)),
            by = .(dataset, sample)]
sc <- join(sc, base[, .(n_usable_ref_baseline = sum(n_ref_sub >= LCC_MIN_REF_SUBCLONES &
                                                    n_ref_cells >= LCC_MIN_REF_CELLS)), by = dataset])

# 2b. main-line healthy false-positive rate: how trustworthy is this dataset's healthy arm?
f <- file.path(DIR_MALIGNANCY, "malignancy_fpr_healthy.csv")
if (file.exists(f)) sc <- join(sc, fread(f)[, .(healthy_fpr_med = round(median(fpr), 3),
                                                healthy_fpr_max = round(max(fpr), 3)), by = dataset])

# 2c. cell-type availability for the planned analyses
nf <- fread(LCC_CCC_NODE_FEAT)
wf <- dcast(nf, sample ~ hierarchy_bin, value.var = "n_cells", fill = 0)
sc <- join(sc, merge(man[, .(dataset, sample)], wf, by = "sample")[
             , .(n_ccc_ready = sum(HSC_MPP >= 30 & Mono_DC >= 30, na.rm = TRUE)), by = dataset])
det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))
sc <- join(sc, det[gene == "COL1A1" & stratum == "Stromal", .(stromal_cells = sum(n_cells)), by = dataset])

# 2d. the two target readouts
sc <- join(sc, fread(file.path(LCC_TAB_DIR, "04_nectin_gate.csv"))[
             , .(nectin4_max_pct = round(max(pct_nonzero), 2),
                 n_nectin4_gt1pct = sum(pct_nonzero > 1)), by = dataset])
# Prefer the 11-dataset --all_datasets table. The 4-dataset one it replaced was retired: keeping
# both meant this line could silently score datasets against the smaller of the two.
f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy_all_datasets.csv")
if (!file.exists(f)) f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy.csv")
if (file.exists(f)) sc <- join(sc, fread(f)[p17_pos == TRUE, .(n_17p_anchors = .N), by = dataset])
sc[, genotype_anchor := fifelse(dataset == "Petti2019", LCC_GENOTYPE_ANCHOR, "")]

for (j in c("stromal_cells", "n_17p_anchors", "n_ccc_ready", "n_usable_ref_baseline"))
  if (j %in% names(sc)) sc[is.na(get(j)), (j) := 0L]
sc[, dataset_role := fifelse(dataset %in% LCC_CORE_DATASETS, "core",
                      fifelse(dataset %in% LCC_STROMA_REF_DATASETS, "stroma_reference", "excluded"))]
setorder(sc, -aml_dg_patients)
fwrite_safe(sc, file.path(LCC_TAB_DIR, "06_dataset_scorecard.csv"))
print(sc[, .(dataset, dataset_role, aml_dg_patients, healthy_patients, icnv_labelled_med,
             n_usable_ref_baseline, healthy_fpr_max, n_ccc_ready, stromal_cells, n_17p_anchors)])

## -- Step 3. selection rationale, explicit ----
rat <- merge(sc[, .(dataset, dataset_role, aml_dg_patients, healthy_patients)],
             LCC_EXCLUDED_DATASETS, by = "dataset", all.x = TRUE)
rat[dataset_role == "core", reason :=
      sprintf("KEEP: within-study contrast %d AML-Dg vs %d healthy patients",
              aml_dg_patients, healthy_patients)]
rat[dataset_role == "stroma_reference", reason := "KEEP as MSC reference only (healthy donors, no AML)"]
setorder(rat, dataset_role, -aml_dg_patients)
fwrite_safe(rat, file.path(LCC_TAB_DIR, "06_selection_rationale.csv"))

## -- Step 4. within-study CK re-test: is the CK failure a cross-study artefact? NO. ----
message("[3] within-study complex-karyotype re-test (methods evidence)")
# If cross-study technical variance were the cause of the cohort-wide CK failure, the arm-event
# burden should separate AML from healthy INSIDE a single study run on one protocol. It does not --
# and in Petti2019 healthy donors carry significantly MORE arm events than AML, which is only
# possible if the metric tracks cell-type composition rather than somatic CNV.
arm <- fread(file.path(LCC_TAB_DIR, "01_arm_events_long.csv"))
sub[man, timepoint := i.timepoint, on = c("dataset", "sample")]
A <- merge(arm, sub[, .(dataset, sample, subclone, is_reference)], by = c("dataset", "sample", "subclone"))
samp <- unique(sub[, .(dataset, sample, timepoint)])
ck_test <- function(ds, th) {
  s <- if (ds == "ALL") samp else samp[dataset == ds]
  a <- if (ds == "ALL") A else A[dataset == ds]
  cnt <- a[frac_arm >= th & is_reference == FALSE, .(na = uniqueN(arm)), by = .(dataset, sample, subclone)
           ][, .(mx = max(na)), by = .(dataset, sample)]
  m <- merge(s, cnt, by = c("dataset", "sample"), all.x = TRUE); m[is.na(mx), mx := 0L]
  d <- m[timepoint == "Diagnosis"]$mx; h <- m[timepoint == "Healthy"]$mx
  if (length(d) < 3 || length(h) < 3) return(NULL)
  w <- suppressWarnings(wilcox.test(d, h))
  data.table(dataset = ds, arm_frac_min = th, n_dg = length(d), n_healthy = length(h),
             median_dg = as.numeric(median(d)), median_healthy = as.numeric(median(h)),
             auc = round(as.numeric(w$statistic) / (length(d) * length(h)), 3),
             p = signif(w$p.value, 3))
}
ck <- rbindlist(lapply(c(LCC_CORE_DATASETS, "GSE116256", "ALL"),
                       function(d) rbindlist(lapply(c(0.5, 0.8), function(t) ck_test(d, t)))))
ck[, verdict := fifelse(auc >= 0.70, "separates", fifelse(auc <= 0.30, "INVERTED (healthy > AML)", "no separation"))]
fwrite_safe(ck, file.path(LCC_TAB_DIR, "06_ck_within_study.csv"))
print(ck)

## -- Step 5. the analysis manifest ----
message("[4] writing the analysis manifest")
out <- man[in_selection == TRUE]
out[, in_primary_contrast := timepoint %in% c("Diagnosis", "Healthy") & dataset_role == "core"]
# one sample per patient for the primary contrast; prefer the library with more cells
out[nf[, .(n_ccc_cells = sum(n_cells)), by = sample], n_ccc_cells := i.n_ccc_cells, on = "sample"]
setorder(out, uid_patient, -n_ccc_cells, sample)
out[in_primary_contrast == TRUE, is_patient_representative := !duplicated(uid_patient)]
fwrite_safe(out, file.path(LCC_TAB_DIR, "06_analysis_manifest.csv"))
message(sprintf("    %d samples selected; primary contrast = %d AML-Dg vs %d healthy patients",
        nrow(out),
        out[in_primary_contrast & is_patient_representative & timepoint == "Diagnosis", .N],
        out[in_primary_contrast & is_patient_representative & timepoint == "Healthy", .N]))
message("[done] 06_select_datasets")
