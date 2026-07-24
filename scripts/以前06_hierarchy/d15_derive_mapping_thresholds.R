#!/usr/bin/env Rscript
# d15_derive_mapping_thresholds.R ----
# Phase 2 step 1.5: derive the per-PLATFORM high_error threshold from HEALTHY cells.
#
# Rationale: the BMM reference is 10x V2 only, so raw mapping error carries a platform component.
# We stratify by platform and set threshold = MAPPING_QC_QUANTILE (default P95) of the mapping
# error of HEALTHY cells on that platform. Then "high_error" (applied in d20) means "worse than
# same-platform healthy cells" -- platform absorbed, biology isolated. Platforms without a healthy
# anchor fall back to BMM's own MAD Pass/Fail (threshold = NA, source = bmm_default).
#
# A manual override file (THRESH_OVERRIDE_TSV: platform<TAB>threshold) supersedes auto-derived
# values for listed platforms -- edit one line to re-flag without re-projecting healthy samples.
#
# PREREQ: d10 has been run on the HEALTHY_DATASETS (their per-cell csvs must exist).
# Usage:  Rscript d15_derive_mapping_thresholds.R

suppressPackageStartupMessages({
  library(data.table)
})
source("d00_config.R")

## -- Platform map ----
if (!file.exists(DATASET_PLAT_TSV)) stop("[d15] platform map not found: ", DATASET_PLAT_TSV)
plat <- fread(DATASET_PLAT_TSV, sep = "\t")                 # dataset, platform, note
stopifnot(all(c("dataset", "platform") %in% names(plat)))
ds2plat <- setNames(plat$platform, plat$dataset)
message("[1] platform map: ", nrow(plat), " datasets across ",
        length(unique(plat$platform)), " platforms")

## -- Gather healthy per-cell mapping error ----
read_percell_me <- function(ds) {
  d <- file.path(PERCELL_DIR, ds)
  f <- list.files(d, pattern = "_percell\\.csv$", full.names = TRUE)
  if (!length(f)) { message("    [warn] no per-cell csv for healthy dataset ", ds,
                            " (run d10 on it first)"); return(NULL) }
  rbindlist(lapply(f, function(x) fread(x, select = c("cell", "Dataset", "mapping_error"))))
}

message("[2] reading healthy datasets: ", paste(HEALTHY_DATASETS, collapse = ", "))
he <- rbindlist(lapply(HEALTHY_DATASETS, read_percell_me), use.names = TRUE)
if (is.null(he) || nrow(he) == 0) {
  stop("[d15] no healthy per-cell data found. Run d10 on HEALTHY_DATASETS first:\n  ",
       paste(HEALTHY_DATASETS, collapse = ", "))
}
he[, platform := ds2plat[Dataset]]
he <- he[!is.na(platform) & is.finite(mapping_error)]
message("    healthy cells with platform + finite error: ", nrow(he))

## -- Per-platform quantile threshold from healthy cells ----
thr_healthy <- he[, .(
  n_healthy_cells    = .N,
  n_healthy_datasets = uniqueN(Dataset),
  threshold          = as.numeric(quantile(mapping_error, MAPPING_QC_QUANTILE, names = FALSE)),
  threshold_source   = paste0("healthy_p", round(MAPPING_QC_QUANTILE * 100))
), by = platform]

## -- All platforms in the map; mark those without a healthy anchor as bmm_default ----
all_plats <- unique(plat$platform)
out <- merge(data.table(platform = all_plats), thr_healthy, by = "platform", all.x = TRUE)
out[is.na(threshold_source), `:=`(
  n_healthy_cells    = 0L,
  n_healthy_datasets = 0L,
  threshold          = NA_real_,                # NA => d20 uses BMM MAD Pass/Fail instead
  threshold_source   = "bmm_default"
)]

## -- Apply manual override if present ----
if (file.exists(THRESH_OVERRIDE_TSV)) {
  ov <- fread(THRESH_OVERRIDE_TSV, sep = "\t")
  stopifnot(all(c("platform", "threshold") %in% names(ov)))
  for (i in seq_len(nrow(ov))) {
    p <- ov$platform[i]; v <- as.numeric(ov$threshold[i])
    out[platform == p, `:=`(threshold = v, threshold_source = "manual_override")]
    message("    [override] platform ", p, " -> threshold ", v)
  }
}

setorder(out, -n_healthy_cells, platform)
fwrite_safe(out, THRESH_TSV, sep = "\t")
message("[3] wrote ", THRESH_TSV)
print(out)

## -- Coverage report ----
n_anch <- out[threshold_source != "bmm_default", .N]
message("[done] ", n_anch, "/", nrow(out), " platforms healthy-anchored; ",
        out[threshold_source == "bmm_default", paste(platform, collapse = ", ")],
        " -> BMM MAD fallback")
