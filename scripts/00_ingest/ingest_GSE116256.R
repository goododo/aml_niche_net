# ============================================================================
# ingest_GSE116256.R   --   van Galen 2019 Cell (Seq-Well)
# Format : per-sample dense matrix  "<GSM>_<SAMPLE>.dem.txt.gz"
#          paired annotation        "<GSM>_<SAMPLE>.anno.txt.gz"
# Special: merge author annotations (CellType / PredictionRefined / ...) into
#          meta.data as "anno_<col>".
# Excluded in this dataset: *.nanopore.txt.gz, cell lines MUTZ3 / OCI-AML3.
# Patient merge: AML<id>-D<day> -> Patient_ID = AML<id>, Timepoint = D<day>;
#                BM<n> / BM5-34p... -> healthy donor, Timepoint = "Baseline".
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE116256"
STUDY   <- "van Galen 2019 Cell"
RAW     <- file.path(GEO_RAW_DIR, "GSE116256_RAW")

# ----------------------------------------------------------------------------
# [1] Discover .dem files and drop excluded samples ----
# ----------------------------------------------------------------------------
message("[1] listing .dem matrices")
dem_files <- list.files(RAW, pattern = "\\.dem\\.txt\\.gz$", full.names = TRUE)
dem_files <- dem_files[!grepl("nanopore", dem_files)]          # drop long-read
dem_files <- dem_files[!grepl("MUTZ3|OCI-AML3", dem_files)]    # drop cell lines
message("    ", length(dem_files), " samples to process")

# ----------------------------------------------------------------------------
# [2] Helpers to derive sample name and patient / timepoint ----
# ----------------------------------------------------------------------------
# "GSM3587923_AML1012-D0.dem.txt.gz" -> "AML1012-D0"
sample_of <- function(f) sub("\\.dem\\.txt\\.gz$", "", sub("^GSM\\d+_", "", basename(f)))

parse_pt_tp <- function(s) {
  pid <- sub("-.*$", "", s)            # AML1012-D0 -> AML1012 ; BM5-34p -> BM5
  if (grepl("^BM", s)) {
    tp <- "Baseline"                   # healthy donor reference
  } else {
    tp <- sub("^[^-]*-", "", s)        # AML1012-D0 -> D0
    if (identical(tp, s)) tp <- "D0"   # safety: no dash -> assume diagnosis
  }
  list(patient = pid, timepoint = tp)
}

anno_path_of <- function(dem_file) sub("\\.dem\\.txt\\.gz$", ".anno.txt.gz", dem_file)

# ----------------------------------------------------------------------------
# [3] Build one Seurat object per sample (no filtering) ----
# ----------------------------------------------------------------------------
message("[3] building per-sample Seurat objects")
obj_list <- list()

for (f in dem_files) {
  smp <- sample_of(f)
  pt  <- parse_pt_tp(smp)
  message("    - ", smp, "  (patient=", pt$patient, ", timepoint=", pt$timepoint, ")")

  counts <- read_dense_txt(f, gene_col = 1)

  # Read paired annotation if present; first column is the cell id.
  anno_file <- anno_path_of(f)
  extra <- if (file.exists(anno_file)) read_meta_txt(anno_file, cell_col = 1) else NULL

  obj_list[[smp]] <- make_seurat(
    counts     = counts,
    sample     = smp,
    dataset    = DATASET,
    patient    = pt$patient,
    timepoint  = pt$timepoint,
    study      = STUDY,
    extra_meta = extra
  )
}

# ----------------------------------------------------------------------------
# [4] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[4] merging ", length(obj_list), " samples")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [5] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[5] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE116256 finished")
