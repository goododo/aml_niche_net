# ============================================================================
# ingest_GSE289435.R   --   Zeng 2025 BCD (AML BMMC, cross-sectional)
# Format : prefixed 10x mtx with DOT separator
#          "<GSM>_<MLL_id>.{barcodes,features,matrix}.tsv.gz"
# Excluded: MLL_29512_PDX (patient-derived xenograft).
# Timepoint: cross-sectional -> "Diagnosis".
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE289435"
STUDY   <- "Zeng 2025 BCD"
RAW     <- file.path(GEO_RAW_DIR, "GSE289435_RAW")

# ----------------------------------------------------------------------------
# [1] Discover prefixes and drop PDX ----
# ----------------------------------------------------------------------------
message("[1] listing library prefixes")
bc_files <- list.files(RAW, pattern = "\\.barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("\\.barcodes\\.tsv\\.gz$", "", bc_files)   # e.g. "GSM8791432_MLL_14666"
prefixes <- prefixes[!grepl("PDX", prefixes)]              # drop xenograft
message("    ", length(prefixes), " libraries to process")

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per patient (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-patient Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  patient <- sub("^GSM\\d+_", "", prefix)   # "MLL_14666"
  message("    - ", prefix, "  (patient=", patient, ")")

  # NOTE the dot separator before barcodes/features/matrix.
  counts <- read_mtx_prefixed(RAW, prefix, sep = ".", features_name = "features")
  obj_list[[prefix]] <- make_seurat(
    counts    = counts,
    sample    = prefix,
    dataset   = DATASET,
    patient   = patient,
    timepoint = "Diagnosis",
    study     = STUDY
  )
}

# ----------------------------------------------------------------------------
# [3] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[3] merging ", length(obj_list), " patients")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE289435 finished")
