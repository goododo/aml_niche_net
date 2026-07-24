# ============================================================================
# ingest_GSE253355.R   --   Bandyopadhyay 2024 Cell (healthy BM niche, stroma)
# Format : prefixed 10x mtx  "<GSM>_<donor>_{barcodes,features,matrix}"
#          donor ids like "H14_MACS", "H21", ...
# Role   : healthy reference + auxiliary stroma. All baseline, no timepoints.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE253355"
STUDY   <- "Bandyopadhyay 2024 Cell"
RAW     <- file.path(GEO_RAW_DIR, "GSE253355_RAW")

# ----------------------------------------------------------------------------
# [1] Discover library prefixes ----
# ----------------------------------------------------------------------------
message("[1] listing library prefixes")
bc_files <- list.files(RAW, pattern = "_barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("_barcodes\\.tsv\\.gz$", "", bc_files)   # e.g. "GSM8019212_H14_MACS"
message("    ", length(prefixes), " libraries found")

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per donor (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-donor Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  donor <- sub("^GSM\\d+_", "", prefix)   # "H14_MACS" / "H21" / ...
  message("    - ", prefix, "  (donor=", donor, ")")

  counts <- read_mtx_prefixed(RAW, prefix, sep = "_", features_name = "features")
  obj_list[[prefix]] <- make_seurat(
    counts    = counts,
    sample    = prefix,
    dataset   = DATASET,
    patient   = donor,        # healthy donor as Patient_ID
    timepoint = "Baseline",
    study     = STUDY
  )
}

# ----------------------------------------------------------------------------
# [3] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[3] merging ", length(obj_list), " donors")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE253355 finished")
