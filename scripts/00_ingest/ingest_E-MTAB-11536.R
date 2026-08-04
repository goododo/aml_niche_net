# ============================================================================
# ingest_E-MTAB-11536.R   --   Dominguez Conde 2022 Science (cross-tissue immune)
# Format : per-donor HDF5  "<donor>-<COMP>.h5"  (COMP = BMA bone marrow / BLD blood)
#          Type unknown a priori -> read_h5_auto() probes 10x vs AnnData.
# Role   : healthy control (immune-focused). All baseline, no timepoints.
# Patient: donor before "-" (e.g. 621B). Compartment = BMA only (see below).
#
# SCOPE (curated metadata, 10 rows, all tissue=BM):
#   - BONE MARROW ONLY. This is a cross-TISSUE atlas; the deposit also ships blood (-BLD).
#     v1 ingested all 12 h5 files, so 4 blood libraries were sitting in the healthy-marrow arm
#     and their composition was being averaged into the B_healthy barycenter. Dropped via
#     INGEST_EXCLUDE.
#   - EXPECTED_BM_DONORS is the curated roster. Two of the ten (D496, D503) are NOT in this
#     deposit's local copy, and the script now says so loudly instead of silently shipping 8.
#     They are the Columbia arm: 3' v3.1 (the other 8 are 5' v1/v2), CD66b+CD235ab+CD326
#     triple-depleted, and 72% of the atlas's marrow cells. Admitting them needs BOTH the
#     erythroid-node mask (curation C10) and donor-equal barycenter weighting (C11) -- without
#     those two they would dominate B_healthy with a depleted composition.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

suppressPackageStartupMessages(library(hdf5r))

DATASET <- "E-MTAB-11536"
STUDY   <- "Dominguez Conde 2022 Science"
# NOTE actual location on disk is under Zenodo, not ArrayExpress.
RAW     <- file.path(ZENODO_DIR, "E-MTAB-11536")

# ----------------------------------------------------------------------------
# [1] H5 auto-reader: 10x CellRanger or AnnData ----
# ----------------------------------------------------------------------------
# Returns a genes x cells sparse matrix regardless of on-disk h5 flavour.
read_h5_auto <- function(path) {
  # First try the standard 10x reader (CellRanger v2/v3).
  out <- tryCatch(read_10x_h5(path), error = function(e) NULL)
  if (!is.null(out)) { message("      [h5] 10x CellRanger format"); return(out) }

  # Otherwise inspect structure with hdf5r.
  h5  <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  top <- names(h5)

  # AnnData (.h5ad): expects /X (CSR), /var, /obs.
  if ("X" %in% top) {
    message("      [h5] AnnData format")
    x      <- h5[["X"]]
    data_v <- x[["data"]]$read()
    indices<- x[["indices"]]$read()
    indptr <- x[["indptr"]]$read()
    shape  <- tryCatch(x[["shape"]]$read(),
                       error = function(e) h5attr(x, "shape"))
    n_obs  <- shape[1]; n_var <- shape[2]

    # AnnData X is cells x genes in CSR; build then transpose to genes x cells.
    m <- Matrix::sparseMatrix(
      j = indices + 1,
      p = indptr,
      x = data_v,
      dims = c(n_obs, n_var)
    )
    genes <- h5[["var"]][["_index"]]$read()
    cells <- h5[["obs"]][["_index"]]$read()
    rownames(m) <- cells; colnames(m) <- genes
    return(Matrix::t(m))   # -> genes x cells
  }

  stop("Unrecognised h5 structure (top-level: ", paste(top, collapse = ", "), ")")
}

# Curated roster of bone-marrow donors (meta_E-MTAB-11536_v2.3.csv, 10 rows).
EXPECTED_BM_DONORS <- c("A29", "A31", "A35", "A36", "A37", "621B", "637C", "640C", "D496", "D503")

# ----------------------------------------------------------------------------
# [2] Discover h5 files, then restrict to the bone-marrow compartment ----
# ----------------------------------------------------------------------------
message("[2] listing h5 files")
h5_files <- list.files(RAW, pattern = "\\.h5$", full.names = TRUE)
message("    ", length(h5_files), " h5 files found")

# INGEST_EXCLUDE drops "-BLD" (blood). Match on the basename, not the full path.
keep_base <- ingest_keep(DATASET, sub("\\.h5$", "", basename(h5_files)), "library")
h5_files  <- h5_files[sub("\\.h5$", "", basename(h5_files)) %in% keep_base]
message("    ", length(h5_files), " bone-marrow libraries retained")

# Roster check: a curated donor with no file on disk is a DOWNLOAD gap, not a design decision,
# and must never be discovered later as "the healthy arm looks small".
present <- sub("-.*$", "", keep_base)
missing <- setdiff(EXPECTED_BM_DONORS, present)
if (length(missing))
  warning("E-MTAB-11536: curated BM donor(s) absent from the local deposit: ",
          paste(missing, collapse = ", "),
          ". The healthy arm is incomplete until these are downloaded ",
          "(processed object at tissueimmunecellatlas.org / CELLxGENE; this accession ",
          "hosts FASTQ only, matrix_level=none). See curation Q3.", call. = FALSE)
if (length(setdiff(present, EXPECTED_BM_DONORS)))
  warning("E-MTAB-11536: donor(s) on disk but NOT in the curated roster: ",
          paste(setdiff(present, EXPECTED_BM_DONORS), collapse = ", "), call. = FALSE)

# "621B-BMA.h5" -> donor "621B", compartment "BMA"
parse_donor_comp <- function(f) {
  b     <- sub("\\.h5$", "", basename(f))   # "621B-BMA"
  donor <- sub("-.*$", "", b)               # "621B"
  comp  <- sub("^.*-", "", b)               # "BMA" / "BLD"
  list(donor = donor, compartment = comp)
}

# ----------------------------------------------------------------------------
# [3] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[3] building per-library Seurat objects")
obj_list <- list()

for (f in h5_files) {
  dc  <- parse_donor_comp(f)
  sid <- sub("\\.h5$", "", basename(f))      # "621B-BMA"
  message("    - ", sid, "  (donor=", dc$donor, ", compartment=", dc$compartment, ")")

  counts <- read_h5_auto(f)
  s <- make_seurat(
    counts    = counts,
    sample    = sid,
    dataset   = DATASET,
    patient   = dc$donor,      # healthy donor as Patient_ID
    timepoint = "Baseline",
    study     = STUDY
  )
  s$Compartment <- dc$compartment   # BMA (bone marrow) vs BLD (blood)
  obj_list[[sid]] <- s
}

# ----------------------------------------------------------------------------
# [4] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[4] merging ", length(obj_list), " libraries")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [5] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[5] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] E-MTAB-11536 finished")
