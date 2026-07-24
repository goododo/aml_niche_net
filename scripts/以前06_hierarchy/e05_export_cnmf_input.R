#!/usr/bin/env Rscript
# e05_export_cnmf_input.R ----
# Phase 3 step 0 (R): export per-sample RAW COUNTS of the cells cNMF will factorize, as an mtx
# triple + meta.json. No HVG here (cNMF picks HVGs itself); we only (a) select cells via the
# three-stage gate and (b) drop technical genes (mito/ribo/HB). Keeps R as the single place that
# touches Seurat v5 objects + the d35 malignant label + d20 bins (reusing get_counts/core16).
#
# THREE-STAGE INPUT GATE (DECISION Q2a):
#   union-malignant cells >= MIN_MALIGNANT_CNMF        -> input_mode="malignant"
#   else fallback-bin cells >= MIN_CELLS_CNMF          -> input_mode="bin_fallback"
#   else                                                -> SKIP (meta.json note="skip_too_few", no mtx)
# Cell-cycle genes are KEPT (proliferation may be a real malignant program; flagged later, not dropped).
#
# OUTPUT: CNMF_INPUT_DIR/<ds>/<sample>/{counts.mtx, barcodes.tsv, genes.tsv, meta.json}
#   counts.mtx is genes x cells (Matrix Market), integer raw counts.
# Per-sample, serial, resume-skip.
# Usage:
#   Rscript e05_export_cnmf_input.R --dataset=GSE227903 --malignancy_dir=/FAST/.../03_malignancy_infercnv/GSE227903
#   Rscript e05_export_cnmf_input.R --dataset=GSE239721 --malignancy_dir=.../GSE239721 --sample=PT10A

suppressPackageStartupMessages({ library(optparse); library(data.table); library(Matrix) })
source("d00_config.R")

PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", help = "dataset id"),
  make_option("--sample",  type = "character", default = NULL, help = "optional single sample"),
  make_option("--malignancy_dir", type = "character", default = NULL,
              help = "dir with <sample>__consensus_percell.csv (default: MALIGNANCY_INFERCNV_DIR/<ds>)")
)))
if (is.null(opt$dataset)) stop("--dataset is required")
mal_dir <- if (!is.null(opt$malignancy_dir)) opt$malignancy_dir else file.path(MALIGNANCY_INFERCNV_DIR, opt$dataset)

## -- helpers ----
core16 <- function(x) sub("-\\d+$", "", sub("^.*_", "", x))
is_tech_gene <- function(g) Reduce(`|`, lapply(TECH_GENE_PATTERNS, function(p) grepl(p, g)))

write_meta <- function(path, lst) {
  # tiny JSON writer (no jsonlite dependency): flat list of scalars/short vectors
  enc <- function(v) {
    if (is.character(v)) {
      if (length(v) == 1) paste0("\"", v, "\"") else
        paste0("[", paste0("\"", v, "\"", collapse = ", "), "]")
    } else if (is.logical(v)) {
      tolower(as.character(v))
    } else {
      if (length(v) == 1) as.character(v) else paste0("[", paste(v, collapse = ", "), "]")
    }
  }
  body <- paste0("  \"", names(lst), "\": ", vapply(lst, enc, character(1)), collapse = ",\n")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0("{\n", body, "\n}"), path)
}

samps <- list_qc_samples(opt$dataset)
if (!is.null(opt$sample)) samps <- samps[sample == opt$sample]
if (!nrow(samps)) stop("[e05] no QC samples for ", opt$dataset)
message("[1] ", nrow(samps), " sample(s); malignancy dir = ", mal_dir)

process_one <- function(ds, samp, rds) {
  out_dir <- file.path(CNMF_INPUT_DIR, ds, samp)
  meta_f  <- file.path(out_dir, "meta.json")
  if (file.exists(meta_f) && file.exists(file.path(out_dir, "counts.mtx"))) {
    message("    [skip] ", samp); return(invisible())
  }

  ## -- malignant label (union, d35) ----
  cons <- file.path(mal_dir, paste0(samp, "__consensus_percell.csv"))
  mal_core <- character(0)
  if (file.exists(cons)) {
    cp <- fread(cons)                          # cell, malignant, confidence
    mal_core <- core16(cp$cell[cp$malignant == 1L])
  } else {
    message("    [warn] ", samp, " no consensus_percell; malignant set empty")
  }

  ## -- counts from QC object (reuse v5/v4-safe get_counts) ----
  obj <- readRDS(rds)
  cnt <- get_counts(obj)                       # genes x cells, raw
  rm(obj)
  cell_core <- core16(colnames(cnt))

  ## -- three-stage gate ----
  in_mal <- cell_core %in% mal_core
  n_mal  <- sum(in_mal)
  bin_csv <- file.path(PERCELL_BINNED_DIR, ds, paste0(samp, "_binned.csv"))
  if (n_mal >= MIN_MALIGNANT_CNMF) {
    sel <- in_mal; mode <- "malignant"; bins_used <- NA_character_
  } else if (file.exists(bin_csv)) {
    bn <- fread(bin_csv, select = c("cell", "hierarchy_bin"))
    bn[, core := core16(cell)]
    fb_core <- bn$core[bn$hierarchy_bin %in% CNMF_FALLBACK_BINS]
    sel <- cell_core %in% fb_core
    if (sum(sel) >= MIN_CELLS_CNMF) {
      mode <- "bin_fallback"; bins_used <- paste(CNMF_FALLBACK_BINS, collapse = "|")
    } else {
      mode <- "skip_too_few"
    }
  } else {
    message("    [warn] ", samp, " no binned csv; cannot fallback")
    mode <- "skip_too_few"
  }

  if (mode == "skip_too_few") {
    write_meta(meta_f, list(dataset = ds, sample = samp, input_mode = "skip_too_few",
                            n_malignant = n_mal, n_cells_selected = 0L,
                            min_malignant = MIN_MALIGNANT_CNMF, min_cells = MIN_CELLS_CNMF,
                            note = "too few cells for cNMF", seed = SEED))
    message(sprintf("    [skip-too-few] %s  n_malignant=%d", samp, n_mal)); return(invisible())
  }

  ## -- subset cells, drop technical genes, keep genes expressed in >=1 selected cell ----
  m <- cnt[, sel, drop = FALSE]
  keep_g <- !is_tech_gene(rownames(m))
  m <- m[keep_g, , drop = FALSE]
  gtot <- Matrix::rowSums(m)
  m <- m[gtot > 0, , drop = FALSE]             # drop all-zero genes (cNMF dislikes them)

  ## -- write mtx triple (genes x cells) + meta ----
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(as(m, "CsparseMatrix"), file.path(out_dir, "counts.mtx"))
  writeLines(colnames(m), file.path(out_dir, "barcodes.tsv"))
  writeLines(rownames(m), file.path(out_dir, "genes.tsv"))
  write_meta(meta_f, list(
    dataset = ds, sample = samp, input_mode = mode,
    n_malignant = n_mal, n_cells_selected = ncol(m), n_genes = nrow(m),
    fallback_bins = if (is.na(bins_used)) "NA" else bins_used,
    tech_patterns = TECH_GENE_PATTERNS,
    min_malignant = MIN_MALIGNANT_CNMF, min_cells = MIN_CELLS_CNMF,
    cnmf_numhvg = CNMF_NUM_HVG, seed = SEED
  ))
  message(sprintf("    [write] %-10s mode=%-12s cells=%d genes=%d (n_mal=%d)",
                  samp, mode, ncol(m), nrow(m), n_mal))
  invisible()
}

for (i in seq_len(nrow(samps))) process_one(samps$dataset[i], samps$sample[i], samps$rds[i])
message("[done] e05 for ", opt$dataset, " -> ", file.path(CNMF_INPUT_DIR, opt$dataset))
