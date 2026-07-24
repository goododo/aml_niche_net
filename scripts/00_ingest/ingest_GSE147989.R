# ============================================================================
# ingest_GSE147989.R   --   Riether 2020 Nat Med (cusatuzumab, sorted CD34 LSC)
# Format : 10x HDF5  "<GSM>_<ID>_filtered_feature_bc_matrix.h5"
#          ID like 009SCRBL / 007C1D1BL
# Patient: leading digits (007 / 009). Timepoint: SCR (screening) or C1D1
#          (cycle 1 day 1). Trailing "BL" = peripheral blood -> Compartment.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE147989"
STUDY   <- "Riether 2020 Nat Med"
RAW     <- file.path(GEO_RAW_DIR, "GSE147989_RAW")

# ----------------------------------------------------------------------------
# [1] Discover h5 files ----
# ----------------------------------------------------------------------------
message("[1] listing h5 files")
h5_files <- list.files(RAW, pattern = "\\.h5$", full.names = TRUE)
message("    ", length(h5_files), " h5 files found")

# "GSM4451499_009SCRBL_filtered_feature_bc_matrix.h5" -> "009SCRBL"
sample_of <- function(f) {
  b <- sub("^GSM\\d+_", "", basename(f))
  sub("_filtered_feature_bc_matrix\\.h5$", "", b)
}

parse_pt_tp <- function(sid) {
  patient    <- stringr::str_extract(sid, "^\\d+")            # "009"
  rest       <- sub("^\\d+", "", sid)                        # "SCRBL" / "C1D1BL"
  compartment <- if (grepl("BL$", rest)) "Blood" else NA_character_
  tp         <- sub("BL$", "", rest)                         # "SCR" / "C1D1"
  list(patient = patient, timepoint = tp, compartment = compartment)
}

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-library Seurat objects")
obj_list <- list()

for (f in h5_files) {
  sid <- sample_of(f)
  pt  <- parse_pt_tp(sid)
  message("    - ", sid, "  (patient=", pt$patient, ", timepoint=", pt$timepoint, ")")

  counts <- read_10x_h5(f)
  s <- make_seurat(
    counts    = counts,
    sample    = sid,
    dataset   = DATASET,
    patient   = pt$patient,
    timepoint = pt$timepoint,
    study     = STUDY
  )
  s$Compartment <- pt$compartment
  obj_list[[sid]] <- s
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

message("[ok] GSE147989 finished")
