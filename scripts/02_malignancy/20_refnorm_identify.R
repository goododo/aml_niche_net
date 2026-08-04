#!/usr/bin/env Rscript
# 20_refnorm_identify.R ----
# Phase 1.1 -- Autologous normal-reference via SingleR vs BoneMarrowMap.
# Per QC'd sample, classify cells against the BoneMarrowMap broad labels with SingleR, keep
# high-confidence mature T/B cells (blast-gated) as the autologous normal reference for inferCNV,
# and export the cell-ids (reusable as Numbat's expression reference).
#
# NOTE: this only anchors the CNV baseline. Malignancy is still called on ALL cells by
# inferCNV+Numbat+VarTrix -> hierarchy (Phase 2) stays orthogonal to malignant identity
# (blueprint R2). SingleR here is a label-transfer to FIND normal lymphocytes, NOT the Phase-2
# projection.
#
# INPUT  : QC objects under QC_RDS_DIR/<Dataset>/<sample>.rds ; BMM_ANNOTATED_RDS (SingleR ref)
# OUTPUT : REFNORM_SUMMARY_CSV (routing table read by 40/44) ;
#          REFNORM_REF_CELL_DIR/<Dataset>/<sample>_ref_norm_cells.txt (autologous ref cell ids)
#
# 3a-2 rewire (no logic change): source -> here 3-line config; the annotated BMM ref is now
#   BMM_ANNOTATED_RDS (was BMM_REF_RDS, collision-fixed); the adequacy threshold is now
#   REFNORM_MIN_NORMAL_CELLS (was MIN_REF_CELLS; value 30->50 comes from config); QC dir is the
#   canonical QC_RDS_DIR (was REFNORM_QC_OBJ_DIR); bare fwrite -> fwrite_safe. The inline
#   .get_counts/.get_data extraction is kept as-is (get_counts adoption deferred to the utils pass).
#   SingleR classification, blast gate, and routing logic are byte-preserved.
#
# Run (from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/20_refnorm_identify.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleR)
  library(SummarizedExperiment)
  library(data.table)
})

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

dir.create(dirname(REFNORM_SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Helpers ----
# ---------------------------------------------------------------------
.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))
.get_data <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "data"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "data"))

.cell_pos_any <- function(cnt, markers) {
  p <- intersect(markers, rownames(cnt))
  if (!length(p)) rep(FALSE, ncol(cnt)) else Matrix::colSums(cnt[p, , drop = FALSE] > 0) >= 1L
}
.cell_n_markers <- function(cnt, markers) {
  p <- intersect(markers, rownames(cnt))
  if (!length(p)) rep(0L, ncol(cnt)) else as.integer(Matrix::colSums(cnt[p, , drop = FALSE] > 0))
}
.decide <- function(n_ref) {
  if (n_ref >= REFNORM_MIN_NORMAL_CELLS) return("autologous_ok")
  switch(REF_FALLBACK_MODE,
         external = "fallback_external", skip_infercnv = "fallback_skip_infercnv",
         halt = "fallback_halt", "fallback_external")
}

# ---------------------------------------------------------------------
# Build (or load) the aggregated BoneMarrowMap reference ONCE ----
#   We cache the pseudobulk-aggregated reference and call the high-level SingleR() per sample:
#   SingleR() auto-intersects genes (robust to per-dataset gene-set differences that break
#   train/classify), and the aggregated ref keeps per-sample marker training cheap.
# ---------------------------------------------------------------------
if (file.exists(SINGLER_AGG_RDS)) {
  message(sprintf("[ref] loading cached aggregated reference: %s", SINGLER_AGG_RDS))
  agg <- readRDS(SINGLER_AGG_RDS)
} else {
  message(sprintf("[ref] aggregating reference from %s (one-time) ...", basename(BMM_ANNOTATED_RDS)))
  stopifnot(file.exists(BMM_ANNOTATED_RDS))
  bbm <- readRDS(BMM_ANNOTATED_RDS)
  ref_lab <- as.character(bbm@meta.data[[BMM_LABEL_COL]])
  stopifnot(!all(is.na(ref_lab)))
  ref_log <- tryCatch({ bbm <- NormalizeData(bbm, verbose = FALSE); .get_data(bbm) },
                      error = function(e) .get_data(bbm))
  keep_ref <- !is.na(ref_lab)
  set.seed(20260613)
  agg <- aggregateReference(ref_log[, keep_ref], labels = ref_lab[keep_ref])
  dir.create(dirname(SINGLER_AGG_RDS), recursive = TRUE, showWarnings = FALSE)
  saveRDS(agg, SINGLER_AGG_RDS)
  rm(bbm, ref_log); gc(verbose = FALSE)
  message("[ref] aggregated reference cached.")
}

# ---------------------------------------------------------------------
# Core: select reference cells for one sample ----
# ---------------------------------------------------------------------
select_ref <- function(seu, sample_id, dataset, agg) {
  n_total <- ncol(seu)
  summ <- data.table(
    sample_id = sample_id, dataset = dataset, n_total = n_total,
    n_singler_TB = 0L, n_ref = 0L, frac_ref = 0,
    val_lineage_pos = NA_real_, val_blast_drop = 0L,
    method = NA_character_, decision = NA_character_,
    fallback = TRUE, fallback_mode = REF_FALLBACK_MODE
  )
  seu$ref_norm <- FALSE

  # --- sorted-library guard: no autologous attempt ----
  if (dataset %in% SORTED_LIBRARY_DATASETS ||
      sample_id %in% SORTED_LIBRARY_SAMPLES ||
      (nzchar(SORTED_LIBRARY_REGEX) && grepl(SORTED_LIBRARY_REGEX, sample_id, ignore.case = TRUE))) {
    message(sprintf("    [guard] sorted/enriched library -> %s", .decide(0L)))
    summ[, `:=`(method = "sorted_guard", decision = .decide(0L))]
    return(list(seu = seu, ref_cells = character(0), summary = summ))
  }

  # --- SingleR classification (whole-profile, ambient-robust) ----
  seu <- NormalizeData(seu, verbose = FALSE)
  pred <- SingleR(test = .get_data(seu), ref = agg, labels = agg$label,
                  fine.tune = SINGLER_FINE_TUNE)
  lab  <- pred$pruned.labels                       # NA where low-confidence (pruned)
  is_TB <- !is.na(lab) & lab %in% REF_KEEP_LABELS
  n_TB  <- sum(is_TB)

  # --- lenient blast-negative gate (SingleR is the primary filter) ----
  cnt     <- .get_counts(seu)
  myelo_n <- .cell_n_markers(cnt, MYELOID_BLAST)
  keep      <- is_TB & (myelo_n <= BLAST_GATE_MAX_MARKERS)
  ref_cells <- colnames(seu)[keep]

  # --- self-validation: do kept cells actually express a lineage marker? ----
  lin_pos     <- .cell_pos_any(cnt, LINEAGE_VALID)
  val_lineage <- if (length(ref_cells)) mean(lin_pos[keep]) else NA_real_

  dec <- .decide(length(ref_cells))
  if (dec == "autologous_ok") {
    seu$ref_norm[keep] <- TRUE
  } else {
    if (dec == "fallback_halt")
      warning(sprintf("    [HALT-FLAG] %s::%s only %d ref cells (< %d)",
                      dataset, sample_id, length(ref_cells), REFNORM_MIN_NORMAL_CELLS))
    ref_cells <- character(0)
  }

  message(sprintf("    SingleR T/B=%d -> ref=%d (lineage+ %s, blast-dropped %d) -> %s",
                  n_TB, length(ref_cells),
                  ifelse(is.na(val_lineage), "NA", sprintf("%.0f%%", 100 * val_lineage)),
                  n_TB - length(ref_cells), dec))

  summ[, `:=`(
    n_singler_TB    = n_TB,
    n_ref           = length(ref_cells),
    frac_ref        = ifelse(n_total > 0, length(ref_cells) / n_total, 0),
    val_lineage_pos = round(val_lineage, 3),
    val_blast_drop  = n_TB - length(ref_cells),
    method          = "singleR",
    decision        = dec,
    fallback        = dec != "autologous_ok"
  )]
  list(seu = seu, ref_cells = ref_cells, summary = summ)
}

# ---------------------------------------------------------------------
# Build the sample list ----
# ---------------------------------------------------------------------
if (nzchar(REFNORM_MANIFEST_CSV) && file.exists(REFNORM_MANIFEST_CSV)) {
  man <- fread(REFNORM_MANIFEST_CSV)
  stopifnot(all(c("sample_id", "dataset", "rds_in") %in% names(man)))
} else {
  rds_in <- list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  rds_in <- rds_in[!grepl("/_", rds_in)]            # drop scratch dirs (e.g. /_dryrun/)
  stopifnot(length(rds_in) > 0)
  man <- data.table(sample_id = sub("\\.rds$", "", basename(rds_in)),
                    dataset   = basename(dirname(rds_in)),
                    rds_in    = rds_in)
}
if (length(REFNORM_SKIP_DATASETS)) {
  n0 <- nrow(man); man <- man[!(dataset %in% REFNORM_SKIP_DATASETS)]
  message(sprintf("[0] skipped %d sample(s) from healthy-only datasets: %s",
                  n0 - nrow(man), paste(REFNORM_SKIP_DATASETS, collapse = ", ")))
}
# Samples whose cells belong to TWO patients cannot get a per-sample CNV profile. inferCNV's
# reference subtraction assumes one genome; on a 2-donor chimera it does not error -- it returns
# a smooth, plausible-looking profile that is the average of two people. That is the worst
# possible failure mode here, because nothing downstream can detect it.
.rman <- file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv")
if (file.exists(.rman)) {
  .m <- fread(.rman)
  if ("patient_resolved" %in% names(.m)) {
    .un <- .m[patient_resolved == FALSE, .(dataset, sample_id = Sample)]
    if (nrow(.un)) {
      n0 <- nrow(man); man <- man[!.un, on = .(dataset, sample_id)]
      message(sprintf("[0] skipped %d unresolved multi-patient sample(s): %s",
                      n0 - nrow(man), paste(.un$sample_id, collapse = ", ")))
    }
  }
}
man[, ref_txt := file.path(REFNORM_REF_CELL_DIR, dataset, paste0(sample_id, "_ref_norm_cells.txt"))]
man[, rds_out := file.path(REFNORM_OUT_OBJ_DIR,  dataset, paste0(sample_id, "_refnorm.rds"))]
message(sprintf("[0] %d sample(s) across %d dataset(s) queued", nrow(man), uniqueN(man$dataset)))

# ---------------------------------------------------------------------
# Main loop -- ONE sample at a time (disk-constrained) ----
# ---------------------------------------------------------------------
all_summ <- vector("list", nrow(man))
for (i in seq_len(nrow(man))) {
  sid <- man$sample_id[i]; ds <- man$dataset[i]
  message(sprintf("[%d/%d] %s :: %s", i, nrow(man), ds, sid))

  seu <- readRDS(man$rds_in[i])
  res <- select_ref(seu, sample_id = sid, dataset = ds, agg = agg)

  if (length(res$ref_cells) > 0L) {
    dir.create(dirname(man$ref_txt[i]), recursive = TRUE, showWarnings = FALSE)
    writeLines(res$ref_cells, man$ref_txt[i])
  }
  if (isTRUE(REFNORM_WRITE_RDS)) {
    dir.create(dirname(man$rds_out[i]), recursive = TRUE, showWarnings = FALSE)
    saveRDS(res$seu, man$rds_out[i])
  }
  all_summ[[i]] <- res$summary
  rm(seu, res); gc(verbose = FALSE)
}

summary_dt <- rbindlist(all_summ, fill = TRUE)
fwrite_safe(summary_dt, REFNORM_SUMMARY_CSV)
message(sprintf("[done] summary -> %s", REFNORM_SUMMARY_CSV))
print(summary_dt[, .(dataset, sample_id, n_total, n_ref, frac_ref, val_lineage_pos, decision)])

# ---------------------------------------------------------------------
# WHAT TO CHECK in the summary:
#   * val_lineage_pos  -> precision audit. Among kept reference cells, the fraction expressing a
#     canonical T/B marker. Expect HIGH (~>0.9). Low val_lineage_pos -> SingleR mislabelling there.
#   * decision == "autologous_ok"  -> has a usable autologous reference.
#   * fallback_external  -> blast-packed / sorted; inferCNV uses external ref, Numbat+VarTrix still
#     provide allele evidence (2/3 consensus intact).
# DOWNSTREAM (40/44 inferCNV driver): verify the ref_norm group is CNV-flat.
# ---------------------------------------------------------------------
