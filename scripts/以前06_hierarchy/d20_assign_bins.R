#!/usr/bin/env Rscript
# d20_assign_bins.R ----
# Phase 2 step 2: collapse the broad (24) cell type to the 8 hierarchy bins, and compute the
# per-cell mapping_qc_flag (pass/high_error) using the per-platform threshold from d15.
# This is the MALIGNANCY-INDEPENDENT backbone -- no Phase-1b labels needed here. Per-sample,
# serial, resume-skip. Writes a sibling 'percell_binned' tree, leaving d10's output immutable.
#
# Usage:
#   Rscript d20_assign_bins.R --dataset=GSE227903
#   Rscript d20_assign_bins.R --dataset=GSE227903 --sample=5143_Dg

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")  # <ds>/<sample>_binned.csv

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", help = "dataset id"),
  make_option("--sample",  type = "character", default = NULL, help = "optional single sample")
)))
if (is.null(opt$dataset)) stop("--dataset is required")

## -- bin map + per-platform threshold ----
bm <- load_bin_map()                                              # celltype_broad, hierarchy_bin, in_ccc_graph
broad2bin <- setNames(bm$hierarchy_bin, bm$celltype_broad)
bin2ccc   <- setNames(as.logical(bm$in_ccc_graph), bm$hierarchy_bin)  # bin-level flag (Stromal = FALSE)

plat    <- fread(DATASET_PLAT_TSV, sep = "\t")
ds_plat <- plat[dataset == opt$dataset, platform]
if (!length(ds_plat)) stop("[d20] dataset ", opt$dataset, " not in ", DATASET_PLAT_TSV)
ds_plat <- ds_plat[1]

if (!file.exists(THRESH_TSV)) stop("[d20] thresholds not found; run d15 first: ", THRESH_TSV)
thr  <- fread(THRESH_TSV, sep = "\t")
trow <- thr[platform == ds_plat]
if (!nrow(trow)) stop("[d20] platform ", ds_plat, " not in thresholds table")
THR_VAL <- as.numeric(trow$threshold[1]); THR_SRC <- trow$threshold_source[1]
message("[1] dataset=", opt$dataset, " platform=", ds_plat,
        " threshold=", ifelse(is.na(THR_VAL), "<bmm_default>", round(THR_VAL, 3)), " (", THR_SRC, ")")

## -- samples ----
pc_dir <- file.path(PERCELL_DIR, opt$dataset)
f <- list.files(pc_dir, pattern = "_percell\\.csv$", full.names = TRUE)
if (!is.null(opt$sample)) f <- f[basename(f) == paste0(opt$sample, "_percell.csv")]
if (!length(f)) stop("[d20] no per-cell csv for ", opt$dataset, " (run d10 first)")
message("[2] ", length(f), " sample(s)")

assign_one <- function(csv) {
  samp <- sub("_percell\\.csv$", "", basename(csv))
  out  <- file.path(PERCELL_BINNED_DIR, opt$dataset, paste0(samp, "_binned.csv"))
  if (file.exists(out)) { message("    [skip] ", samp); return(invisible()) }
  dt <- fread(csv)

  # broad -> hierarchy_bin (validate vocabulary: refuse silent mismaps)
  miss <- setdiff(unique(dt$bmm_celltype_broad), names(broad2bin))
  if (length(miss)) stop("[d20] ", samp, " broad labels not in bin map: ", paste(miss, collapse = ", "))
  dt[, hierarchy_bin := broad2bin[bmm_celltype_broad]]
  dt[, in_ccc_graph  := bin2ccc[hierarchy_bin]]

  # mapping_qc_flag (PRIMARY): scale-free KNN projection confidence. high_error = prob below
  # MIN_PROJ_PROB. Comparable across platform/depth/dataset (prob is a 0-1 proportion).
  if (!"bmm_celltype_fine_prob" %in% names(dt)) {
    stop("[d20] ", samp, " lacks bmm_celltype_fine_prob (re-run d10 to add it)")
  }
  dt[, mapping_qc_flag := fifelse(bmm_celltype_fine_prob < MIN_PROJ_PROB, "high_error", "pass")]

  # Absolute-error flag kept as a REFERENCE column only (not the primary flag).
  if (is.na(THR_VAL)) {
    dt[, mapping_qc_flag_abs := fifelse(mapping_error_QC_bmm == "Fail", "high_error", "pass")]
  } else {
    dt[, mapping_qc_flag_abs := fifelse(mapping_error > THR_VAL, "high_error", "pass")]
  }
  dt[, `:=`(mapping_threshold = THR_VAL, threshold_source = THR_SRC, platform = ds_plat,
            min_proj_prob = MIN_PROJ_PROB)]

  keep <- c("cell", "Dataset", "Sample", "bmm_celltype_fine", "bmm_celltype_broad",
            "hierarchy_bin", "in_ccc_graph", "bmm_celltype_fine_prob", "predicted_pseudotime",
            "mapping_error", "mapping_qc_flag", "min_proj_prob",
            "mapping_qc_flag_abs", "mapping_threshold", "threshold_source", "platform")
  keep <- intersect(keep, names(dt))
  fwrite_safe(dt[, ..keep], out)
  message("    [write] ", out, "  (", nrow(dt), " cells; high_error=",
          dt[mapping_qc_flag == "high_error", .N], ")")
  invisible()
}

for (x in f) assign_one(x)
message("[done] d20 for ", opt$dataset)