# ============================================================================
# ingest_GSE207356.R   --   Nicosia 2023 Cancer Cell (CCS1477 time course)
# Format : prefixed 10x mtx  "<GSM>_TS173_##_Sample_X_{barcodes,features,matrix}"
# Scope  : exploratory, SINGLE patient (TS173), 3 on-treatment timepoints
#          (Sample_D / E / F) merged under orig.ident = "TS173".
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE207356"
STUDY   <- "Nicosia 2023 Cancer Cell"
RAW     <- file.path(GEO_RAW_DIR, "GSE207356_RAW")

# ----------------------------------------------------------------------------
# [1] Discover library prefixes ----
# ----------------------------------------------------------------------------
message("[1] listing library prefixes")
bc_files <- list.files(RAW, pattern = "_barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("_barcodes\\.tsv\\.gz$", "", bc_files)   # "GSM6284912_TS173_04_Sample_D"
message("    ", length(prefixes), " libraries found")

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per timepoint (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-timepoint Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  sid     <- sub("^GSM\\d+_", "", prefix)              # "TS173_04_Sample_D"
  patient <- stringr::str_extract(sid, "^TS\\d+")      # "TS173" (single patient)
  tp      <- stringr::str_extract(sid, "Sample_[A-Z]") # "Sample_D" as timepoint
  message("    - ", prefix, "  (patient=", patient, ", timepoint=", tp, ")")

  counts <- read_mtx_prefixed(RAW, prefix, sep = "_", features_name = "features")
  obj_list[[prefix]] <- make_seurat(
    counts    = counts,
    sample    = prefix,
    dataset   = DATASET,
    patient   = patient,
    timepoint = tp,
    study     = STUDY
  )
}

# ----------------------------------------------------------------------------
# [3] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[3] merging ", length(obj_list), " timepoints")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE207356 finished")
