#!/usr/bin/env Rscript
# 40_infercnv_run.R ----
# Phase 1.2 (inferCNV arm of the 2/3 malignancy consensus) -- LOOPED multi-sample runner.
# Runs inferCNV per sample with a HYBRID reference strategy driven by the 20 ref_norm_summary.csv:
#   * autologous_ok            -> reference = this sample's ref_norm cells
#   * method == "sorted_guard" -> SKIP inferCNV (Numbat+VarTrix carry it)
#   * singleR + fallback       -> EXTERNAL reference (BoneMarrowMap subsample)
#
# Terminal safety: for autologous samples we measure REFERENCE FLATNESS (reference cells should be
# CNV-flat). If a reference is not flat (contaminated by blasts), the sample is FLAGGED so it can
# be re-routed to the external reference.
#
# INPUT  : REFNORM_SUMMARY_CSV (routing) ; per-sample QC rds under QC_RDS_DIR
# OUTPUT : INFERCNV_BURDEN_ROOT/<ds>/<sample>_infercnv_burden.csv (per-cell CNV burden) +
#          INFERCNV_SUMMARY_CSV (one row per sample + ref flatness flag) + full inferCNV dirs
#
# inferCNV is heavy; runs serially (disk-constrained). Re-running skips samples that already
# produced a burden file (resume-friendly). For parallel runs use 44 + 45 (array).
#
# 3a-2 rewire (no logic change): source -> here 3-line config; all inferCNV params + external-ref
#   constants now come from config_malignancy.R (removed the USER-FILL block); output roots use the
#   canonical INFERCNV_ROOT / INFERCNV_BURDEN_ROOT, QC dir is QC_RDS_DIR, external-ref seed is
#   INFERCNV_EXT_REF_SEED (== 20260613), INFERCNV_THREADS now env-overridable via config (default 4,
#   unchanged); bare fwrite -> fwrite_safe. Inline .get_counts kept (get_counts adoption deferred).
#
# Run (from PROJECT ROOT): Rscript scripts/02_malignancy/40_infercnv_run.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(infercnv)
  library(data.table)
})

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

stopifnot(file.exists(INFERCNV_GENE_ORDER))
dir.create(INFERCNV_ROOT,        recursive = TRUE, showWarnings = FALSE)
dir.create(INFERCNV_BURDEN_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(INFERCNV_SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Helpers ----
# ---------------------------------------------------------------------
.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

# Build (or load) the external healthy reference: a stratified BMM subsample (raw counts).
get_external_ref <- function() {
  if (file.exists(INFERCNV_EXT_REF_CACHE)) {
    message(sprintf("[ext-ref] loading cached external reference: %s", INFERCNV_EXT_REF_CACHE))
    return(readRDS(INFERCNV_EXT_REF_CACHE))
  }
  message("[ext-ref] building external reference from BoneMarrowMap subsample (one-time) ...")
  stopifnot(file.exists(INFERCNV_EXT_REF_RDS))
  bbm <- readRDS(INFERCNV_EXT_REF_RDS)
  lab <- as.character(bbm@meta.data[[INFERCNV_EXT_REF_LABELCOL]])
  set.seed(INFERCNV_EXT_REF_SEED)
  idx <- unlist(lapply(split(seq_along(lab), lab), function(ii)
    if (length(ii) <= INFERCNV_EXT_REF_PER_TYPE) ii else sample(ii, INFERCNV_EXT_REF_PER_TYPE)))
  cnt <- .get_counts(bbm)[, idx, drop = FALSE]
  colnames(cnt) <- paste0("EXTREF_", colnames(cnt))   # namespace to avoid id collisions
  ext <- list(counts = cnt, cells = colnames(cnt))
  dir.create(dirname(INFERCNV_EXT_REF_CACHE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ext, INFERCNV_EXT_REF_CACHE)
  rm(bbm); gc(verbose = FALSE)
  message(sprintf("[ext-ref] %d external reference cells cached.", ncol(cnt)))
  ext
}

# Run inferCNV for one sample given a counts matrix + annotation + ref group name(s).
run_one <- function(counts, anno_dt, ref_groups, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  anno_file <- file.path(out_dir, "annotations.tsv")
  fwrite_safe(anno_dt, anno_file, sep = "\t", col.names = FALSE)
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file  = anno_file,
    gene_order_file   = INFERCNV_GENE_ORDER,
    ref_group_names   = ref_groups
  )
  infercnv::run(
    obj, cutoff = INFERCNV_CUTOFF, out_dir = out_dir,
    cluster_by_groups = TRUE, analysis_mode = INFERCNV_ANALYSIS,
    denoise = INFERCNV_DENOISE, HMM = INFERCNV_HMM,
    num_threads = INFERCNV_THREADS, no_plot = FALSE, save_rds = TRUE
  )
}

# Per-cell CNV burden from the final inferCNV residual matrix (centered ~1).
burden_from_obj <- function(obj) {
  ex <- obj@expr.data                      # genes x cells
  colMeans((ex - 1)^2)
}

# ---------------------------------------------------------------------
# Route table from the 20 summary ----
# ---------------------------------------------------------------------
stopifnot(file.exists(REFNORM_SUMMARY_CSV))
summ <- fread(REFNORM_SUMMARY_CSV)
# routing: autologous | external | skip
summ[, route := fifelse(decision == "autologous_ok", "autologous",
                 fifelse(method == "sorted_guard",   "skip", "external"))]
message(sprintf("[0] routes: autologous=%d, external=%d, skip=%d",
                sum(summ$route == "autologous"), sum(summ$route == "external"),
                sum(summ$route == "skip")))

ext_ref <- NULL  # lazily built only if any external sample is processed
out_rows <- list()

for (i in seq_len(nrow(summ))) {
  sid <- summ$sample_id[i]; ds <- summ$dataset[i]; route <- summ$route[i]
  burden_csv <- file.path(INFERCNV_BURDEN_ROOT, ds, paste0(sid, "_infercnv_burden.csv"))
  out_dir    <- file.path(INFERCNV_ROOT, ds, sid)

  row <- data.table(dataset = ds, sample_id = sid, route = route,
                    status = NA_character_, n_ref = 0L, n_obs = 0L,
                    ref_flatness = NA_real_, ref_flag = "", obs_mean_burden = NA_real_)

  if (route == "skip") {
    row[, status := "skipped_sorted"]; out_rows[[i]] <- row
    message(sprintf("[%d/%d] %s::%s  SKIP (sorted)", i, nrow(summ), ds, sid)); next
  }
  if (file.exists(burden_csv)) {  # resume
    row[, status := "already_done"]; out_rows[[i]] <- row
    message(sprintf("[%d/%d] %s::%s  already done", i, nrow(summ), ds, sid)); next
  }

  message(sprintf("[%d/%d] %s::%s  route=%s", i, nrow(summ), ds, sid, route))
  rds_in <- file.path(QC_RDS_DIR, ds, paste0(sid, ".rds"))
  if (!file.exists(rds_in)) { row[, status := "missing_rds"]; out_rows[[i]] <- row; next }

  res <- tryCatch({
    seu  <- readRDS(rds_in)
    s_cnt <- .get_counts(seu)

    if (route == "autologous") {
      ref_txt <- file.path(REFNORM_REF_CELL_DIR, ds, paste0(sid, "_ref_norm_cells.txt"))
      if (!file.exists(ref_txt)) stop("ref_norm cell list missing for autologous sample")
      ref_cells <- readLines(ref_txt)
      ref_cells <- intersect(ref_cells, colnames(s_cnt))
      grp <- ifelse(colnames(s_cnt) %in% ref_cells, "reference_normal", "observation")
      counts <- s_cnt
      anno   <- data.table(cell = colnames(s_cnt), group = grp)
      ref_groups <- "reference_normal"
    } else {                       # external
      if (is.null(ext_ref)) ext_ref <<- get_external_ref()
      common <- intersect(rownames(s_cnt), rownames(ext_ref$counts))
      stopifnot(length(common) > 1000)
      counts <- cbind(s_cnt[common, , drop = FALSE], ext_ref$counts[common, , drop = FALSE])
      anno   <- data.table(cell = colnames(counts),
                           group = c(rep("observation", ncol(s_cnt)),
                                     rep("reference_external", ncol(ext_ref$counts))))
      ref_groups <- "reference_external"
    }

    n_obs <- sum(anno$group == "observation")
    if (n_obs < INFERCNV_MIN_OBS) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(sprintf("skipped: %d observation cell(s) (< INFERCNV_MIN_OBS=%d); near-all-normal / degenerate",
                         n_obs, INFERCNV_MIN_OBS), file.path(out_dir, "lowobs.skip"))
      message(sprintf("  [skip] %s: only %d observation cell(s) -> inferCNV not applicable", sid, n_obs))
      row[, `:=`(status = "skipped_low_obs", n_ref = sum(anno$group != "observation"), n_obs = n_obs)]
      rm(seu, s_cnt, counts); gc(verbose = FALSE)
      row
    } else {
      obj <- run_one(counts, anno, ref_groups, out_dir)
      burden <- burden_from_obj(obj)

      bdt <- data.table(cell = names(burden), infercnv_burden = as.numeric(burden))
      bdt[, group := anno$group[match(cell, anno$cell)]]
      dir.create(dirname(burden_csv), recursive = TRUE, showWarnings = FALSE)
      fwrite_safe(bdt, burden_csv)

      ref_ids  <- anno$cell[anno$group %in% ref_groups]
      flat     <- mean(bdt$infercnv_burden[bdt$cell %in% ref_ids], na.rm = TRUE)
      obs_mean <- mean(bdt$infercnv_burden[!(bdt$cell %in% ref_ids)], na.rm = TRUE)

      row[, `:=`(
        status = "ok",
        n_ref  = length(ref_ids),
        n_obs  = sum(!(anno$cell %in% ref_ids)),
        ref_flatness = round(flat, 5),
        obs_mean_burden = round(obs_mean, 5),
        ref_flag = if (route == "autologous" && !is.na(flat) && flat > REF_FLATNESS_MAX)
                     "CHECK_reference_not_flat->consider_external" else ""
      )]
      rm(seu, s_cnt, counts, obj); gc(verbose = FALSE)
      row
    }
  }, error = function(e) { row[, status := paste0("error: ", conditionMessage(e))]; row })

  out_rows[[i]] <- res
  fwrite_safe(rbindlist(out_rows, fill = TRUE), INFERCNV_SUMMARY_CSV)  # checkpoint each sample
}

summary_dt <- rbindlist(out_rows, fill = TRUE)
fwrite_safe(summary_dt, INFERCNV_SUMMARY_CSV)
message(sprintf("[done] summary -> %s", INFERCNV_SUMMARY_CSV))
print(summary_dt[, .(dataset, sample_id, route, status, n_ref, ref_flatness, ref_flag)])

# ---------------------------------------------------------------------
# WHAT TO CHECK:
#   * ref_flag == "CHECK_reference_not_flat..." -> autologous reference is contaminated (the low
#     val_lineage outliers). Re-route those to external: set their decision to fallback in the 20
#     output (or add to a manual external list) and re-run.
#   * Calibrate REF_FLATNESS_MAX from the ref_flatness distribution of clean samples before trusting
#     the flag.
#   * Per-cell `infercnv_burden` feeds the 2/3 consensus (next: 41 -> 50).
# ---------------------------------------------------------------------
