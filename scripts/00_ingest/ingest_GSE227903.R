# ============================================================================
# ingest_GSE227903.R   --   Ennis 2023 iScience (BM niche interactome; Dx/MRD/R)
# Format : per-sample dense matrix  "<GSM>_<patient>_<TP>_count_mat.txt.gz"
#          paired metadata          "<GSM>_<patient>_<TP>_meta.txt.gz"
# IMPORTANT: count_mat is CELLS x GENES (transposed) with an off-by-one header:
#            line1 = 36601 gene names (no corner label), data rows = barcode +
#            36601 counts. Gene names are make.names()-sanitised (dots). We read
#            it explicitly, transpose to genes x cells, and normalise symbols
#            against the GRCh38-2020-A reference.
# Special: merge author metadata into meta.data as "anno_<col>".
# Patient: numeric patient id + timepoint token Dg / MRD / R / R2.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE227903"
STUDY   <- "Ennis 2023 iScience"
RAW     <- file.path(GEO_RAW_DIR, "GSE227903_RAW")

# ----------------------------------------------------------------------------
# [1] Reader for the transposed dense matrix ----
# ----------------------------------------------------------------------------
# line1  : 36601 gene names (header, no barcode-column label)
# line2+ : barcode + 36601 counts  (36602 fields)
# Returns a genes x cells sparse matrix with normalised gene symbols.
read_count_mat_t <- function(path) {
  con   <- gzfile(path, "r")
  genes <- strsplit(readLines(con, n = 1), "\t", fixed = TRUE)[[1]]   # 36601 dotted names
  close(con)

  body  <- data.table::fread(path, sep = "\t", header = FALSE, skip = 1, data.table = FALSE)
  bc    <- as.character(body[[1]])                 # cell barcodes
  mat   <- as.matrix(body[, -1, drop = FALSE])     # cells x genes
  storage.mode(mat) <- "double"

  stopifnot(ncol(mat) == length(genes))            # guard the off-by-one assumption
  colnames(mat) <- make.unique(normalize_symbols(genes))   # dotted -> standard symbols
  rownames(mat) <- bc

  Matrix::Matrix(t(mat), sparse = TRUE)            # -> genes x cells
}

# ----------------------------------------------------------------------------
# [2] Discover count matrices + parse patient / timepoint ----
# ----------------------------------------------------------------------------
message("[2] listing count matrices")
mat_files <- list.files(RAW, pattern = "_count_mat\\.txt\\.gz$", full.names = TRUE)
message("    ", length(mat_files), " samples to process")

# "GSM7110543_1216_Dg_count_mat.txt.gz" -> "1216_Dg"
sample_of <- function(f) sub("_count_mat\\.txt\\.gz$", "", sub("^GSM\\d+_", "", basename(f)))

TP_LABEL <- c(Dg = "Diagnosis", MRD = "MRD", R = "Relapse", R2 = "Relapse2")

parse_pt_tp <- function(sid) {
  tp_tok  <- stringr::str_extract(sid, "(Dg|MRD|R2|R)$")            # R2 before R
  patient <- sub("_(Dg|MRD|R2|R)$", "", sid)                       # "1216"
  tp      <- ifelse(is.na(tp_tok), NA_character_, TP_LABEL[[tp_tok]])
  list(patient = patient, timepoint = tp)
}

meta_path_of <- function(mat_file) sub("_count_mat\\.txt\\.gz$", "_meta.txt.gz", mat_file)

# ----------------------------------------------------------------------------
# [3] Build one Seurat object per sample (no filtering) ----
# ----------------------------------------------------------------------------
message("[3] building per-sample Seurat objects")
obj_list <- list()

for (f in mat_files) {
  sid <- sample_of(f)
  pt  <- parse_pt_tp(sid)
  message("    - ", sid, "  (patient=", pt$patient, ", timepoint=", pt$timepoint, ")")

  counts <- read_count_mat_t(f)

  # Paired metadata; first column is the cell id.
  meta_file <- meta_path_of(f)
  extra <- if (file.exists(meta_file)) read_meta_txt(meta_file, cell_col = 1) else NULL
  if (!is.null(extra)) {
    n_overlap <- length(intersect(colnames(counts), rownames(extra)))
    if (n_overlap == 0) warning("Sample ", sid, ": meta barcodes do not overlap the matrix")
  }

  obj_list[[sid]] <- make_seurat(
    counts     = counts,
    sample     = sid,
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

message("[ok] GSE227903 finished")
