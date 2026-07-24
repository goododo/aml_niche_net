# ============================================================================
# ingest_Chen2023.R   --   Chen 2023 Blood Cancer Discovery (paired BM niche)
# Format : nested 10x dirs
#          "<root>/filtered_feature_bc_matrix.7z_<SAMPLE>/filtered_feature_bc_matrix/"
#          SAMPLE like AML103_CD34, AML162_niche_immune, NBM1_CD34, ...
# Patient: strip "_CD34"/"_Niche_Immune" (case-insensitive) -> Patient_ID;
#          keep the fraction in a Compartment column. orig.ident = Patient_ID,
#          so the two sorted libraries of a patient group together but stay
#          distinguishable by Compartment.
# Timepoint: NBM* -> healthy "Baseline"; AML* -> "Diagnosis".
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "Chen2023"
STUDY   <- "Chen 2023 Blood Cancer Discovery"
ROOT    <- file.path(ZENODO_DIR, "10.17632_gwjh3w6ztm.2")

# ----------------------------------------------------------------------------
# [1] Discover per-sample 10x directories ----
# ----------------------------------------------------------------------------
message("[1] listing sample directories")
sample_dirs <- list.dirs(ROOT, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("filtered_feature_bc_matrix\\.7z_", basename(sample_dirs))]
message("    ", length(sample_dirs), " sample folders found")

# Derive Compartment (normalised) and Patient_ID from a sample name.
compartment_of <- function(s) if (grepl("CD34", s, ignore.case = TRUE)) "CD34" else "Niche_Immune"
patient_of     <- function(s) sub("_(CD34|niche.*|NICHE.*|Niche.*)$", "", s)

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-library Seurat objects")
obj_list <- list()

for (d in sample_dirs) {
  sample      <- sub("^filtered_feature_bc_matrix\\.7z_", "", basename(d))
  inner       <- file.path(d, "filtered_feature_bc_matrix")
  patient     <- patient_of(sample)
  compartment <- compartment_of(sample)
  tp          <- if (grepl("^NBM", patient)) "Baseline" else "Diagnosis"
  message("    - ", sample, "  (patient=", patient, ", compartment=", compartment, ")")

  counts <- read_10x_dir(inner)
  s <- make_seurat(
    counts    = counts,
    sample    = sample,
    dataset   = DATASET,
    patient   = patient,
    timepoint = tp,
    study     = STUDY
  )
  s$Compartment <- compartment   # CD34 (malignant/HSPC) vs Niche_Immune (receiver)
  obj_list[[sample]] <- s
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

message("[ok] Chen2023 finished")
