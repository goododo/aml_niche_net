#!/usr/bin/env Rscript
# f05_platform_deviation_control.R ----
# Normal-cell CONTROL for the cross-dataset projection-deviation gradient.
#
# WHY: the headline "high_error %" per dataset (healthy 3.8% < GSE227903 9.1% < GSE239721 14.3%)
# CONFLATES two things: (a) real malignant deviation from the healthy reference, and (b) technical
# platform/chemistry/sampling distance from the BMM reference. A single per-dataset scalar cannot
# separate them. This control uses NON-MALIGNANT bins (T_NK, B_Plasma -- normal bystanders in AML)
# as a per-dataset TECHNICAL BASELINE: normal T/NK cells should carry NO disease deviation, so their
# high_error rate is (almost) pure platform effect. If T_NK shows the same dataset gradient, the
# gradient is largely platform, not biology.
#
# OUTPUT (results/tables/04_hierarchy/qc_platform_control/):
#   deviation_by_bin.csv   dataset x hierarchy_bin: n_cells, high_error_frac (+ mean_prob)
#   platform_control.csv   per dataset: T_NK baseline vs myeloid-malignant bins, and the DIFFERENCE
#                          (myeloid high_error - T_NK baseline) = platform-corrected "suspected biology"
#
# REUSABLE QC STEP: run this for every dataset as it enters the study to quantify its technical
# deviation baseline. Reads only percell_binned/*/*_binned.csv (no heavy recompute).
#
# Usage:
#   Rscript f05_platform_deviation_control.R --datasets=E-MTAB-11536,GSE227903,GSE239721
#   Rscript f05_platform_deviation_control.R --datasets=GSE227903 --min_bin=30

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")
OUT_DIR            <- file.path(FAST_DIR, "results/tables/04_hierarchy/qc_platform_control")

# non-malignant reference bins (normal bystanders) vs myeloid bins where AML blasts concentrate
NORMAL_REF_BINS <- c("T_NK", "B_Plasma")
MYELOID_BINS    <- c("HSC_MPP", "LMPP_GMP")
BASELINE_BIN    <- "T_NK"   # primary technical baseline (most universally non-malignant)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", help = "comma-separated dataset ids"),
  make_option("--min_bin", type = "integer", default = 30L,
              help = "min cells for a dataset x bin cell to be reported (default 30)")
)))
if (is.null(opt$datasets)) stop("--datasets is required")
datasets <- strsplit(opt$datasets, ",")[[1]]
MIN_BIN <- opt$min_bin

## -- load all binned cells for the requested datasets ----
load_ds <- function(ds) {
  fs <- list.files(file.path(PERCELL_BINNED_DIR, ds), pattern = "_binned\\.csv$", full.names = TRUE)
  if (!length(fs)) { message("[!] no binned csv for ", ds); return(NULL) }
  d <- rbindlist(lapply(fs, function(f)
    fread(f, select = c("Dataset", "Sample", "hierarchy_bin",
                        "bmm_celltype_fine_prob", "mapping_qc_flag", "platform"))), fill = TRUE)
  d
}
all <- rbindlist(lapply(datasets, load_ds), fill = TRUE)
if (!nrow(all)) stop("[f05] nothing loaded")
all[, is_high_error := as.integer(mapping_qc_flag == "high_error")]
message("[1] ", nrow(all), " cells across ", uniqueN(all$Dataset), " dataset(s)")

## -- main table: dataset x bin high_error fraction ----
by_bin <- all[, .(n_cells = .N,
                  high_error_frac = round(mean(is_high_error), 4),
                  mean_prob = round(mean(bmm_celltype_fine_prob, na.rm = TRUE), 3)),
              by = .(Dataset, hierarchy_bin, platform)]
by_bin <- by_bin[n_cells >= MIN_BIN]
setorder(by_bin, Dataset, hierarchy_bin)
fwrite_safe(by_bin, file.path(OUT_DIR, "deviation_by_bin.csv"))

## -- platform control: T_NK baseline vs myeloid, per dataset ----
# baseline = high_error in the normal reference bin(s); myeloid = high_error in HSC/LMPP bins
base_dt <- all[hierarchy_bin == BASELINE_BIN, .(baseline_n = .N,
              baseline_high_error = mean(is_high_error)), by = Dataset]
mye_dt  <- all[hierarchy_bin %in% MYELOID_BINS, .(myeloid_n = .N,
              myeloid_high_error = mean(is_high_error)), by = Dataset]
refall_dt <- all[hierarchy_bin %in% NORMAL_REF_BINS, .(normalref_n = .N,
              normalref_high_error = mean(is_high_error)), by = Dataset]

ctrl <- Reduce(function(a, b) merge(a, b, by = "Dataset", all = TRUE),
               list(base_dt, refall_dt, mye_dt))
ctrl[, `:=`(
  baseline_high_error   = round(baseline_high_error, 4),
  normalref_high_error  = round(normalref_high_error, 4),
  myeloid_high_error    = round(myeloid_high_error, 4),
  # platform-corrected suspected biology = myeloid deviation MINUS normal-cell technical baseline
  suspected_biology_vs_TNK      = round(myeloid_high_error - baseline_high_error, 4),
  suspected_biology_vs_normalref= round(myeloid_high_error - normalref_high_error, 4)
)]
setcolorder(ctrl, c("Dataset", "baseline_n", "baseline_high_error",
                    "normalref_n", "normalref_high_error",
                    "myeloid_n", "myeloid_high_error",
                    "suspected_biology_vs_TNK", "suspected_biology_vs_normalref"))
setorder(ctrl, Dataset)
fwrite_safe(ctrl, file.path(OUT_DIR, "platform_control.csv"))

## -- console readout ----
message("\n[done] f05 -> ", OUT_DIR)
message("\n== high_error fraction by dataset x bin (n>=", MIN_BIN, ") ==")
print(dcast(by_bin, hierarchy_bin ~ Dataset, value.var = "high_error_frac"))
message("\n== platform control: does the normal-cell baseline ALSO show the gradient? ==")
print(ctrl[, .(Dataset, TNK_baseline = baseline_high_error, myeloid = myeloid_high_error,
               suspected_biology = suspected_biology_vs_TNK)])
message("\nRead: if TNK_baseline rises across datasets in the SAME order as myeloid, the raw gradient",
        "\n      is largely PLATFORM. The 'suspected_biology' column is the platform-corrected residual.")
