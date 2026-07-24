#!/usr/bin/env Rscript
# d10_project_bmm.R ----
# Phase 2 step 1: project each per-sample QC Seurat object onto BoneMarrowMap (Symphony),
# transfer FINE cell-type (55) + pseudotime for ALL cells (no dropping), and record the
# continuous mapping error. We do NOT compute the high_error flag here -- that needs the
# per-platform threshold from d15. We keep BMM's own MAD Pass/Fail for reference only.
#
# Per-sample, integration-free (vars = NULL), serial, resume-skip on existing per-cell csv.
#
# Usage:
#   Rscript d10_project_bmm.R --dataset GSE227903                    # all samples in a dataset
#   Rscript d10_project_bmm.R --dataset GSE227903 --sample 5143_Dg   # one sample
#
# FIRST-RUN NOTE: d10 logs EVERY new column each BMM function adds ([SCHEMA] lines), so we can
# pin the exact output schema before building d15/d20/d30. Run it on one healthy + one AML
# sample first and send back the [SCHEMA] lines + head() of the per-cell csv.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratObject)
  library(symphony)
  library(BoneMarrowMap)
  library(data.table)
})

source("d00_config.R")
set.seed(SEED)

## -- Args ----
opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset",  type = "character", help = "dataset id (subdir under QC_RDS_DIR)"),
  make_option("--sample",   type = "character", default = NULL, help = "optional single sample id"),
  make_option("--save_obj", type = "logical",   default = TRUE, help = "save projected Seurat obj to LARGE1")
)))
if (is.null(opt$dataset)) stop("--dataset is required")

## -- Load reference once ----
message("[1] loading BoneMarrowMap Symphony reference ...")
ref <- readRDS(BMM_REF_RDS)
ref$save_uwot_path <- BMM_UWOT_MODEL                 # REQUIRED for query UMAP projection
stopifnot(all(c("Z_corr", "R", "centroids") %in% names(ref)))
if (!file.exists(BMM_UWOT_MODEL)) stop("[1] uwot model not found: ", BMM_UWOT_MODEL)
message("    ref ok; uwot model = ", BMM_UWOT_MODEL)

## -- Sample list ----
samps <- list_qc_samples(opt$dataset)
if (!is.null(opt$sample)) samps <- samps[sample == opt$sample]
if (nrow(samps) == 0) {
  stop("no matching samples for dataset=", opt$dataset,
       if (!is.null(opt$sample)) paste0(" sample=", opt$sample) else "")
}
message("[2] ", nrow(samps), " sample(s) to project for ", opt$dataset)

## -- Per-sample projection ----
project_one <- function(ds, samp, rds_path, save_obj = TRUE) {
  out_csv <- file.path(PERCELL_DIR, ds, paste0(samp, "_percell.csv"))
  if (file.exists(out_csv)) { message("    [skip] ", ds, "/", samp, " (exists)"); return(invisible()) }
  message("    [proj] ", ds, "/", samp)

  obj <- readRDS(rds_path)
  cts <- get_counts(obj)

  # 1) Symphony mapping. Per-sample => no within-query batch => vars = NULL.
  q <- map_Query(exp_query = cts, metadata_query = obj@meta.data, ref_obj = ref, vars = NULL)

  # 2) Mapping error: continuous score + BMM's own MAD Pass/Fail. Capture added cols.
  cb <- colnames(q@meta.data)
  q  <- calculate_MappingError(q, reference = ref, MAD_threshold = MAD_THRESHOLD_BMM)
  add_me <- setdiff(colnames(q@meta.data), cb)
  message("      [SCHEMA] calculate_MappingError added: ", paste(add_me, collapse = ", "))

  # 3) FINE cell type (55). We need labels for ALL cells (high-error malignant blasts are SIGNAL,
  #    not garbage), but predict_CellTypes sets NA where mapping_error_QC == 'Fail'. So predict
  #    twice: (a) QC forced to 'Pass' -> canonical all-cell label; (b) real QC -> NA-on-fail ref.
  stopifnot("mapping_error_QC" %in% colnames(q@meta.data))
  qc_real <- q$mapping_error_QC
  q$mapping_error_QC <- "Pass"
  cb <- colnames(q@meta.data)
  q  <- predict_CellTypes(query_obj = q, ref_obj = ref, final_label = "bmm_fine_all")
  add_ct <- setdiff(colnames(q@meta.data), cb)
  message("      [SCHEMA] predict_CellTypes added: ", paste(add_ct, collapse = ", "))
  q$mapping_error_QC <- qc_real
  q  <- predict_CellTypes(query_obj = q, ref_obj = ref, final_label = "bmm_fine_bmmqc")

  # 4) Pseudotime (cheap; feeds LSC stemness ordering later). Non-fatal if it fails.
  cb <- colnames(q@meta.data)
  q  <- tryCatch(
    predict_Pseudotime(query_obj = q, ref_obj = ref,
                       initial_label = "bmm_pt_init", final_label = "bmm_pt_all"),
    error = function(e) { message("      [warn] predict_Pseudotime failed: ", conditionMessage(e)); q }
  )
  add_pt <- setdiff(colnames(q@meta.data), cb)
  if (length(add_pt)) message("      [SCHEMA] predict_Pseudotime added: ", paste(add_pt, collapse = ", "))

  md <- q@meta.data
  # Identify the continuous mapping-error column (numeric, among calculate_MappingError's adds).
  me_num       <- add_me[vapply(add_me, function(c) is.numeric(md[[c]]), logical(1))]
  me_score_col <- if (length(me_num)) me_num[1] else NA_character_
  me_qc_col    <- if ("mapping_error_QC" %in% colnames(md)) "mapping_error_QC" else NA_character_
  message("      [SCHEMA] mapping_error score col = ", me_score_col, " ; QC col = ", me_qc_col)

  # Build per-cell table. bmm_celltype_fine = ALL cells (canonical, never dropped).
  dt <- data.table(
    cell                    = rownames(md),
    Dataset                 = ds,
    Sample                  = samp,
    bmm_celltype_fine       = as.character(md[["bmm_fine_all"]]),       # fine 55, ALL cells (canonical)
    bmm_celltype_broad      = as.character(md[["bmm_fine_all_Broad"]]), # broad 24, direct from predict_CellTypes
    bmm_celltype_fine_prob  = if ("bmm_fine_all_prob" %in% colnames(md)) as.numeric(md[["bmm_fine_all_prob"]]) else NA_real_,  # KNN confidence
    bmm_celltype_fine_bmmqc = as.character(md[["bmm_fine_bmmqc"]]),     # fine, NA on BMM MAD fail (reference)
    mapping_error           = if (!is.na(me_score_col)) as.numeric(md[[me_score_col]]) else NA_real_,
    mapping_error_QC_bmm    = if (!is.na(me_qc_col)) as.character(md[[me_qc_col]]) else NA_character_,
    predicted_pseudotime    = if ("bmm_pt_all" %in% colnames(md)) as.numeric(md[["bmm_pt_all"]]) else NA_real_
  )

  # Projected UMAP coords (for QC plots), if present.
  umap <- tryCatch(as.data.table(Embeddings(q, "umap_projected"), keep.rownames = "cell"),
                   error = function(e) NULL)
  if (!is.null(umap) && ncol(umap) >= 3) {
    setnames(umap, 2:3, c("umap1", "umap2"))
    dt <- merge(dt, umap[, .(cell, umap1, umap2)], by = "cell", all.x = TRUE, sort = FALSE)
  }

  fwrite_safe(dt, out_csv)
  message("      [write] ", out_csv, "  (", nrow(dt), " cells)")

  if (save_obj) {
    obj_out <- file.path(PROJ_OBJ_DIR, ds, paste0(samp, ".rds"))
    dir.create(dirname(obj_out), recursive = TRUE, showWarnings = FALSE)
    saveRDS(q, obj_out)
    message("      [save] ", obj_out)
  }
  invisible()
}

for (i in seq_len(nrow(samps))) {
  project_one(samps$dataset[i], samps$sample[i], samps$rds[i], save_obj = opt$save_obj)
}
message("[done] d10 projection finished for ", opt$dataset)
