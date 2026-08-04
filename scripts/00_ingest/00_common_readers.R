# 00_common_readers.R ----
# Shared readers + Seurat builder for the 00_ingest stage (formerly 00_common_qc.R).
# Every ingest_<dataset>.R sources this (via here) after the config layer.
#
# Design rules (unchanged):
#   - NO cell filtering. Every raw cell is kept in the RDS (min.cells=0, min.features=0).
#   - QC metrics (percent.mt/ribo/hb) computed PER SAMPLE before merge.
#   - One merged Seurat object saved per dataset.
# Batch-1 changes vs legacy 00_common_qc.R:
#   - paths/constants now come from config_paths.R + config_qc.R (no local path defs)
#   - make_seurat: adds uid_patient (C3), Disease_state default (C4), canonical Timepoint (C1)
#   - add_qc_metrics: warns if a gene-pattern set is empty (A1)
#   - normalize_symbols: guards against duplicate cleaned symbols (A2)
#   - summarise_qc: adds min/max alongside median (feedback #4)

suppressPackageStartupMessages({
  library(here)
  library(Seurat); library(Matrix); library(data.table)
  library(readr);  library(dplyr);  library(stringr); library(ggplot2)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

## -- QC metrics on a single-sample object ----
# hb pattern "^HB[^P]" per project convention; MITO_SHORT covers prefix-stripped datasets.
add_qc_metrics <- function(seu) {
  genes <- rownames(seu)
  mito_genes <- unique(c(grep(PAT_MITO, genes, value = TRUE), intersect(genes, MITO_SHORT)))
  ribo_genes <- grep(PAT_RIBO, genes, value = TRUE)
  hb_genes   <- grep(PAT_HB,   genes, value = TRUE)

  # A1: warn (don't silently assign 0) when a pattern matches nothing -- otherwise a
  # prefix-stripped dataset could carry percent.mt==0 for every cell and go unnoticed.
  if (!length(mito_genes)) warning("add_qc_metrics: NO mitochondrial genes matched (", PAT_MITO,
                                   " / MITO_SHORT). percent.mt will be 0 for all cells.")
  if (!length(ribo_genes)) warning("add_qc_metrics: NO ribosomal genes matched (", PAT_RIBO, ").")
  if (!length(hb_genes))   warning("add_qc_metrics: NO hemoglobin genes matched (", PAT_HB, ").")

  seu[["percent.mt"]]   <- if (length(mito_genes)) PercentageFeatureSet(seu, features = mito_genes) else 0
  seu[["percent.ribo"]] <- if (length(ribo_genes)) PercentageFeatureSet(seu, features = ribo_genes) else 0
  seu[["percent.hb"]]   <- if (length(hb_genes))   PercentageFeatureSet(seu, features = hb_genes)   else 0
  seu
}

## -- Seurat builder: one sample -> one object, no filtering ----
# counts     : dgCMatrix (genes x cells)
# sample     : library-level unique id (e.g. a GSM or "AML1012-D0")
# dataset    : dataset id (e.g. "GSE116256")
# patient    : Patient_ID; orig.ident is set to this (patient-level identity)
# timepoint  : RAW timepoint label (kept in Timepoint_detail; canonicalised into Timepoint)
# disease    : Disease_state (used by canonical timepoint for Petti/Lasry; default "Unknown")
# extra_meta : optional df, rownames = barcodes, author annotations merged as "anno_<col>"
make_seurat <- function(counts, sample, dataset, patient, timepoint,
                        study = NA_character_, disease = NA_character_, extra_meta = NULL) {
  seu <- CreateSeuratObject(counts = counts, project = patient,
                            min.cells = 0, min.features = 0)

  disease <- if (is.na(disease) || !nzchar(disease)) "Unknown" else as.character(disease)  # C4

  seu$Dataset          <- as.character(dataset)
  seu$Study            <- as.character(study)
  seu$Sample           <- as.character(sample)
  seu$Patient_ID       <- as.character(patient)     # never let bare numeric ids sneak in
  seu$uid_patient      <- paste0(as.character(dataset), ":", as.character(patient))  # C3 namespaced key
  seu$Disease_state    <- disease                   # C4 always present
  seu$Timepoint_detail <- as.character(timepoint)   # RAW value preserved
  # sample/patient are passed so canonicalize_timepoint() rule 0 can catch healthy donors that sit
  # inside an AML dataset (Chen2023 NBM*, Petti2019 ND_*, GSE116256 BM<n>) -- dataset-level rules
  # label those "Diagnosis", which put normal marrow into every AML comparison.
  seu$Timepoint        <- canonicalize_timepoint(dataset, timepoint, disease,
                                                 sample_id = sample, patient_id = patient)  # C1 canonical
  # NOMINAL, not a residual-disease statement -- see the note above timepoint_days() in
  # config_qc.R. The H3 stratum is derived in 03.5_residual_stratum.R from frac_malignant.
  seu$Timepoint_days   <- timepoint_days(timepoint) # D<n> where the deposit encodes it, else NA
  .plat                <- platform_of(dataset)      # library-level override applied by the ingest
  seu$Platform         <- .plat$platform            # if the dataset is not platform-homogeneous
  seu$Chemistry        <- .plat$chemistry
  # BM unless the sample is one of the 6 known peripheral-blood draws. A PB sample structurally
  # cannot contain stroma / HSPC / erythroid, so this must be readable downstream as a TECHNICAL
  # absence, never as a biological one.
  seu$Tissue           <- tissue_of(dataset, sample)
  # Defaults: one donor per library, cells attributable to that donor. The two multiplexed
  # deposits (GSE185381, GSE185991) override these after the object is built.
  seu$N_donors_in_library <- 1L
  seu$Patient_resolved    <- TRUE
  seu$orig.ident       <- as.character(patient)     # merge unit is the patient

  if (!is.null(extra_meta)) {
    common <- intersect(colnames(seu), rownames(extra_meta))
    for (cn in colnames(extra_meta)) {
      v <- setNames(rep(NA, ncol(seu)), colnames(seu))
      if (length(common)) v[common] <- extra_meta[common, cn]
      seu[[paste0("anno_", cn)]] <- v
    }
  }

  seu <- add_qc_metrics(seu)
  seu
}

## -- Readers: one helper per on-disk format ----

# Map make.names()-sanitised symbols back to standard symbols via a reference features.tsv (col2).
normalize_symbols <- function(genes, ref_features_path = REF_FEATURES_2020A) {
  ref       <- data.table::fread(ref_features_path, header = FALSE, data.table = FALSE)[[2]]
  ref_clean <- make.names(ref, unique = FALSE)
  if (length(ref) == length(genes) && all(ref_clean == genes)) return(ref)  # exact, same order
  # A2: duplicate cleaned symbols make setNames() keep only the last -> wrong lookup. Warn.
  dup <- anyDuplicated(ref_clean)
  if (dup > 0L) warning("normalize_symbols: reference has duplicate cleaned symbols (e.g. '",
                        ref_clean[dup], "'); lookup keeps the last occurrence for those.")
  out <- unname(setNames(ref, ref_clean)[genes])
  na  <- is.na(out)
  out[na] <- genes[na]
  if (any(na)) message("        normalize_symbols: ", sum(na), " gene(s) kept original (no ref match)")
  out
}

# 10x mtx triplet with a per-sample prefix in a flat GEO_RAW folder.
read_mtx_prefixed <- function(dir, prefix, sep = "_", features_name = "features") {
  build <- function(suffix) file.path(dir, paste0(prefix, sep, suffix))
  mtx <- build("matrix.mtx.gz"); bc <- build("barcodes.tsv.gz")
  ft  <- build(paste0(features_name, ".tsv.gz"))
  stopifnot(file.exists(mtx), file.exists(bc), file.exists(ft))
  tryCatch(
    ReadMtx(mtx = mtx, cells = bc, features = ft, feature.column = 2),
    error = function(e) {
      message("    Fallback: features column 1 (orig error: ", e$message, ")")
      ReadMtx(mtx = mtx, cells = bc, features = ft, feature.column = 1)
    })
}

read_10x_dir <- function(dir) Read10X(data.dir = dir)
read_10x_h5  <- function(path) Read10X_h5(filename = path)

# Dense text expression matrix (genes x cells); first column = gene names.
read_dense_txt <- function(path, gene_col = 1) {
  dt    <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  genes <- as.character(dt[[gene_col]])
  mat   <- as.matrix(dt[, -gene_col, drop = FALSE])
  rownames(mat) <- make.unique(genes)
  Matrix::Matrix(mat, sparse = TRUE)
}

# Per-cell metadata / annotation table; first column = cell id.
read_meta_txt <- function(path, cell_col = 1, sep = "\t") {
  dt <- data.table::fread(path, sep = sep, header = TRUE, data.table = FALSE)
  rn <- as.character(dt[[cell_col]]); dt[[cell_col]] <- NULL
  rownames(dt) <- make.unique(rn)
  dt
}

## -- Merge a named list of per-sample objects (unique barcodes via add.cell.ids) ----
merge_samples <- function(obj_list) {
  stopifnot(length(obj_list) >= 1)
  if (length(obj_list) == 1) return(obj_list[[1]])
  merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = names(obj_list))
}

## -- Outputs: per-cell table + per-sample summary + RDS ----
# feedback #4: min/max alongside median so a bad sample (e.g. max percent.mt ~90%) is visible.
summarise_qc <- function(seu) {
  seu@meta.data %>%
    dplyr::group_by(Dataset, Sample, Patient_ID, Timepoint) %>%
    dplyr::summarise(
      n_cells         = dplyr::n(),
      median_nCount   = median(nCount_RNA), min_nCount = min(nCount_RNA), max_nCount = max(nCount_RNA),
      median_nFeature = median(nFeature_RNA), min_nFeature = min(nFeature_RNA), max_nFeature = max(nFeature_RNA),
      median_pct_mt   = median(percent.mt),   max_pct_mt   = max(percent.mt),
      median_pct_ribo = median(percent.ribo), max_pct_ribo = max(percent.ribo),
      median_pct_hb   = median(percent.hb),   max_pct_hb   = max(percent.hb),
      .groups = "drop")
}

write_qc_outputs <- function(seu, dataset) {
  percell <- seu@meta.data; percell$cell_barcode <- rownames(percell)
  fwrite_safe(percell, file.path(DIR_INGEST, paste0(dataset, "_qc_percell.csv.gz")))
  smry <- summarise_qc(seu)
  fwrite_safe(smry, file.path(DIR_INGEST, paste0(dataset, "_qc_summary.csv")))
  saveRDS(seu, file.path(RDS_INGEST_DIR, paste0(dataset, ".rds")))
  message("[done] ", dataset, ": ", ncol(seu), " cells across ", nrow(smry), " samples")
  invisible(smry)
}

## -- reusable QC plot theme (for the later distribution-review step) ----
vln_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

message("[init] 00_common_readers loaded.")