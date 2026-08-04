# ============================================================================
# ingest_Petti2019.R   --   Petti 2019 Nat Commun (mutation-anchored AML + ND)
# Format : flat dir, DOT-prefixed 10x mtx with genes.tsv
#          "<ID>.{barcodes,genes,matrix}.tsv.gz"
# Sample types:
#   numeric ID (508084, 548327, ...)        -> AML, Timepoint = "Diagnosis"
#   ND_*  /  Normal_sorted_*                -> EXCLUDED (see below)
# Each ID is its own patient (cross-sectional; no longitudinal merge here).
#
# SCOPE: the curated metadata is AML-ONLY. The 4 normal donors this deposit ships alongside the
# AML samples are excluded via INGEST_EXCLUDE, so `classify()` below now only ever returns AML.
# They are kept in the code path (not deleted) because the deposit still contains them and a
# future decision to re-admit should be a one-line table edit, not a rewrite.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "Petti2019"
STUDY   <- "Petti 2019 Nat Commun"
RAW     <- file.path(ZENODO_DIR, "10.5281_zenodo.3345981")

# ----------------------------------------------------------------------------
# [1] Discover sample prefixes ----
# ----------------------------------------------------------------------------
message("[1] listing sample prefixes")
bc_files <- list.files(RAW, pattern = "\\.barcodes\\.tsv\\.gz$", full.names = FALSE)
prefixes <- sub("\\.barcodes\\.tsv\\.gz$", "", bc_files)   # "508084", "ND_083017", ...
message("    ", length(prefixes), " samples found")
prefixes <- ingest_keep(DATASET, prefixes, "sample")
message("    ", length(prefixes), " samples retained")

# Classify a sample id into disease state + timepoint.
classify <- function(id) {
  if (grepl("^ND|^Normal", id)) list(state = "Healthy", tp = "Baseline")
  else                          list(state = "AML",     tp = "Diagnosis")
}

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per sample (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-sample Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  cl <- classify(prefix)
  message("    - ", prefix, "  (state=", cl$state, ", timepoint=", cl$tp, ")")

  # DOT separator, genes.tsv.
  counts <- read_mtx_prefixed(RAW, prefix, sep = ".", features_name = "genes")
  s <- make_seurat(
    counts    = counts,
    sample    = prefix,
    dataset   = DATASET,
    patient   = prefix,        # each ID is its own patient
    timepoint = cl$tp,         # raw ("Baseline"/"Diagnosis"); canonicalised inside make_seurat
    study     = STUDY,
    disease   = cl$state       # AML/Healthy -> drives canonical Timepoint (Petti is disease-driven)
  )
  # (make_seurat now sets Disease_state and canonical Timepoint from `disease`)
  obj_list[[prefix]] <- s
}

# ----------------------------------------------------------------------------
# [3] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[3] merging ", length(obj_list), " samples")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] Petti2019 finished")
