# ============================================================================
# ingest_GSE201966.R   --   Relapse fingerprint, acute monocytic AML (post-HSCT)
# Format : prefixed 10x mtx with genes.tsv (v2)
#          "<GSM>_<AMLxxMx>_<TP>_{barcodes,genes,matrix}"
# Patient: AMLxxMx (e.g. AML01M4) -> Patient_ID; TP in {P, RE, CR}:
#          P = Primary, RE = Relapse, CR = Complete remission -> Timepoint.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE201966"
STUDY   <- "GSE201966 (acute monocytic AML relapse)"
RAW     <- file.path(GEO_RAW_DIR, "GSE201966_RAW")

# ----------------------------------------------------------------------------
# [1] Discover library prefixes ----
# ----------------------------------------------------------------------------
message("[1] listing library prefixes")
bc_files <- list.files(RAW, pattern = "_barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("_barcodes\\.tsv\\.gz$", "", bc_files)   # "GSM6085176_AML01M4_CR"
message("    ", length(prefixes), " libraries found")

# Map short timepoint token to a readable label.
TP_LABEL <- c(P = "Primary", RE = "Relapse", CR = "Complete_remission")

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-library Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  sid     <- sub("^GSM\\d+_", "", prefix)          # "AML01M4_CR"
  patient <- sub("_(P|RE|CR)$", "", sid)           # "AML01M4"
  tp_tok  <- stringr::str_extract(sid, "(P|RE|CR)$")
  tp      <- ifelse(is.na(tp_tok), NA_character_, TP_LABEL[[tp_tok]])
  message("    - ", prefix, "  (patient=", patient, ", timepoint=", tp, ")")

  # genes.tsv instead of features.tsv.
  counts <- read_mtx_prefixed(RAW, prefix, sep = "_", features_name = "genes")
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
message("[3] merging ", length(obj_list), " libraries")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE201966 finished")
