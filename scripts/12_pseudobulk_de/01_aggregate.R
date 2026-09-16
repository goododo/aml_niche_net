#!/usr/bin/env Rscript
# 01_aggregate.R ----
# Stage A of the pseudobulk AML-vs-healthy DE. Sums RAW COUNTS per (sample, hierarchy_bin) over the
# production cell set, and records the per-bin sequencing depth that Stage B needs as a covariate.
# Does NOT normalize, filter genes, or fit anything -- voom needs raw counts and its own library
# sizes, so every decision that shrinks the data belongs downstream.
#
# INPUT  : results/tables/07_fgw/fgw_input_index.csv          (arm: the boolean `healthy`)
#          results/tables/01_preprocess/02_study_split.csv     (Discovery / Validation, dataset-level)
#          results/tables/01_preprocess/00_curated_manifest.csv(library_id, timepoint, sorting)
#          QC_RDS_DIR/<ds>/<sample>.rds                        (Seurat v5, counts-only layer)
#          ANNO_RECONCILED_DIR/<ds>/<sample>__anno_percell.csv  (hierarchy_bin + QC flags + bmm_broad)
# OUTPUT : DIR_PSEUDOBULK/pb_sample_manifest.csv               (the 72-row roster; --row indexes it)
#          PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds            (list; see BUILD below)
#          DIR_PSEUDOBULK/depth/<ds>__<sample>__pb_depth.csv   (per-bin depth + conservation check)
# Usage  :
#   Rscript scripts/12_pseudobulk_de/01_aggregate.R --list                     # write manifest, print roster
#   Rscript scripts/12_pseudobulk_de/01_aggregate.R --dataset=Chen2023 --sample=NBM1   # one sample
#   Rscript scripts/12_pseudobulk_de/01_aggregate.R --row=$SLURM_ARRAY_TASK_ID  # array element
#
# WHY THIS EXISTS AND WHAT IS BINDING. The design is pre-registered in
# PREREGISTRATION_pseudobulk_de.md, committed before this script was written. Three of its clauses
# are implemented here and cannot be changed without amending that document:
#
#   PREREG 4.2  The cell set is FOUR conditions, not the three CellChat uses:
#               in_ccc_graph & !high_error & hierarchy_bin %in% CCC_NODES & !(bmm_broad %in% CCC_EXCLUDE_FINE)
#               Bins come from ANNO_RECONCILED_DIR, not the raw projection: on GSE116256 alone 6.8%
#               of cells change hierarchy_bin after reconciliation. The fourth condition drops
#               "Early Lymphoid", which is 0.95% of the healthy non-malignant T_NK pool against
#               10.64% of the AML one (config_ccc.R:119-123) -- an 11x arm-asymmetric contaminant
#               sitting inside one of the three primary bins. It is the filter behind
#               ccc_node_features.csv, which is what the prereg's binding inclusion table counts.
#
#   PREREG 4.3  Per-(sample, bin) depth is written HERE because it is free at aggregation time and
#               unrecoverable later without re-reading every Seurat object. No per-bin depth column
#               exists anywhere in results/tables/. The prereg's mandatory depth-sensitivity arm
#               (7.1) uses the per-bin median_umi, NOT the sample-level med_ncount_final that every
#               earlier stage used.
#
#   PREREG 8.1  Two EXACT conservation equalities, not an order-of-magnitude check. Summing over
#               bins must reproduce the filtered totals to the count. An order-of-magnitude test
#               passes on dropped cells, double-counted barcodes and bad joins; this does not.
#
# THE ROSTER IS THE 72 BOTH-ARM SAMPLES, AND THAT IS A STATISTICAL DECISION, NOT A CONVENIENCE.
# 7 of 10 datasets carry only one arm, so for those 66 samples dataset and disease status are
# perfectly collinear -- no covariate and no batch correction separates perfectly collinear factors.
# They are excluded in full (PREREG 3.1). The roster is cross-checked against qc_rds_roster() so a
# sample that is in the FGW index but has no QC object fails loudly here rather than silently
# shrinking the cohort downstream.

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here); library(Seurat); library(Matrix)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))   # ANNO_RECONCILED_DIR lives here
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--list",    action = "store_true", default = FALSE, help = "rebuild + print the roster, then exit"),
  make_option("--row",     type = "integer",   default = NA_integer_,   help = "1-based row of pb_sample_manifest.csv"),
  make_option("--dataset", type = "character", default = NA_character_),
  make_option("--sample",  type = "character", default = NA_character_),
  make_option("--force",   action = "store_true", default = FALSE, help = "recompute even if output exists")
)))

BOTH_ARM_DATASETS <- c("Chen2023", "GSE185381", "GSE116256")   # PREREG 3.1, the only both-arm datasets

## ---------------------------------------------------------------- roster ----
build_roster <- function() {
  fgw <- fread(file.path(TAB_DIR, "07_fgw", "fgw_input_index.csv"))
  spl <- fread(file.path(DIR_PREPROCESS, "02_study_split.csv"))[, .(dataset, split)]
  man <- fread(file.path(DIR_PREPROCESS, "00_curated_manifest.csv"),
               select = c("dataset", "sample", "library_id", "sorting_r", "cell_prep_r", "sample_state"))

  R <- fgw[dataset %in% BOTH_ARM_DATASETS,
           .(dataset, sample, timepoint, arm = fifelse(healthy %in% TRUE, "healthy", "AML"))]
  R <- spl[R, on = "dataset"]
  R <- man[R, on = c("dataset", "sample")]

  # PREREG 3.4: Validation drops BM5-34p (CD34+ sorted + fresh against an otherwise
  # viability-only / ficoll / cryopreserved cohort) and keeps AML only at Diagnosis (4 of its 16
  # samples are post-treatment marrows from patients already in the arm). Discovery is untouched --
  # it is already 100% Diagnosis with zero repeated patients.
  R[, in_analysis := TRUE]
  R[, drop_reason := NA_character_]
  R[split == "Validation" & sample == "BM5-34p",
    `:=`(in_analysis = FALSE, drop_reason = "PREREG 3.4a: CD34pos-sorted, fresh; arm-confounded")]
  R[split == "Validation" & arm == "AML" & timepoint != "Diagnosis",
    `:=`(in_analysis = FALSE, drop_reason = "PREREG 3.4b: post-treatment marrow, patient repeated at Diagnosis")]

  # The roster must exist on disk as a QC object. qc_rds_roster() is the project's only sanctioned
  # cohort source (utils.R:74); building a roster from ls() or from a task .tsv has silently shrunk
  # audits in this repo before.
  QC <- qc_rds_roster(on_extra = "ignore")[, .(dataset, sample, rds)]
  R  <- QC[R, on = c("dataset", "sample")]
  if (anyNA(R$rds))
    stop(sprintf("%d roster sample(s) have no QC object: %s",
                 sum(is.na(R$rds)), paste(R[is.na(rds), paste0(dataset, "/", sample)], collapse = ", ")))

  setorder(R, dataset, sample)
  R[, row := .I]
  setcolorder(R, c("row", "dataset", "sample", "split", "arm", "timepoint", "in_analysis",
                   "drop_reason", "library_id", "sorting_r", "cell_prep_r", "sample_state", "rds"))
  R[]
}

manifest_path <- file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv")

if (opt$list || !file.exists(manifest_path)) {
  R <- build_roster()
  fwrite(R, manifest_path)
  message("[roster] wrote ", manifest_path, " : ", nrow(R), " rows (", sum(R$in_analysis), " in analysis)")
  print(R[, .N, by = .(split, arm, in_analysis)][order(split, arm)])
  if (opt$list) quit(save = "no", status = 0L)
}
MAN <- fread(manifest_path)

## ------------------------------------------------------------- one sample ----
aggregate_one <- function(ds, smp) {
  out_rds <- file.path(PB_RDS_DIR, ds, paste0(smp, "__pseudobulk.rds"))
  out_csv <- file.path(DIR_PSEUDOBULK, "depth", paste0(ds, "__", smp, "__pb_depth.csv"))
  if (!opt$force && file.exists(out_rds) && file.exists(out_csv)) {
    message("  [skip] ", ds, "/", smp, " exists"); return(invisible(NULL))
  }
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

  ann_f <- file.path(ANNO_RECONCILED_DIR, ds, paste0(smp, "__anno_percell.csv"))
  if (!file.exists(ann_f)) stop("no reconciled annotation for ", ds, "/", smp, " at ", ann_f)

  d <- fread(ann_f, select = c("cell", "hierarchy_bin", "in_ccc_graph", "high_error", "bmm_broad"))
  # fwrite writes NA as "" and fread reads it back as "" -- normalise or every is.na() is FALSE.
  # Copied from 03_node_features.R:117, where the same trap was already paid for.
  d[!nzchar(trimws(hierarchy_bin)), hierarchy_bin := NA_character_]
  d[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]

  n_ann <- nrow(d)
  d <- d[in_ccc_graph == TRUE & high_error == FALSE & hierarchy_bin %in% CCC_NODES]
  n_after3 <- nrow(d)
  if (length(CCC_EXCLUDE_FINE)) d <- d[!(bmm_broad %in% CCC_EXCLUDE_FINE)]   # PREREG 4.2, 4th condition
  n_after4 <- nrow(d)
  if (!nrow(d)) stop("0 cells survive the PREREG 4.2 filter for ", ds, "/", smp)

  obj <- readRDS(file.path(QC_RDS_DIR, ds, paste0(smp, ".rds")))
  cnt <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "counts"),
                  error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts"))
  # Exact barcode join, in the object's own column order. An inexact join here is the failure mode
  # that produces a correctly shaped matrix of the wrong cells, so it is checked, not assumed.
  j <- match(colnames(cnt), d$cell)
  hit <- !is.na(j)
  frac_matched <- mean(d$cell %in% colnames(cnt))
  if (frac_matched < 0.999)
    stop(sprintf("barcode join %.4f < 0.999 for %s/%s -- barcode convention has drifted, do not report this run",
                 frac_matched, ds, smp))
  cnt <- cnt[, hit, drop = FALSE]
  bin <- factor(d$hierarchy_bin[j[hit]], levels = CCC_NODES)
  stopifnot(ncol(cnt) == length(bin), !anyNA(bin))

  # Per-cell totals BEFORE aggregation: they are the depth covariate (PREREG 4.3) and one side of
  # the conservation check (PREREG 8.1).
  umi_cell  <- Matrix::colSums(cnt)
  gene_cell <- Matrix::colSums(cnt > 0)

  present <- levels(bin)[table(bin) > 0L]
  pb <- vapply(present, function(b) as.numeric(Matrix::rowSums(cnt[, bin == b, drop = FALSE])),
               numeric(nrow(cnt)))
  dimnames(pb) <- list(rownames(cnt), present)

  stats <- data.table(
    dataset = ds, sample = smp, hierarchy_bin = present,
    n_cells      = as.integer(table(bin)[present]),
    total_counts = as.numeric(colSums(pb)),
    median_umi   = vapply(present, function(b) as.numeric(median(umi_cell[bin == b])),  numeric(1)),
    mean_umi     = vapply(present, function(b) as.numeric(mean(umi_cell[bin == b])),    numeric(1)),
    median_genes = vapply(present, function(b) as.numeric(median(gene_cell[bin == b])), numeric(1)),
    n_genes_detected = as.integer(colSums(pb > 0))
  )
  stats[, passes_30 := n_cells >= CCC_MIN_CELLS_PER_OCCUPIED_BIN]   # PREREG 5.1, unit gate

  ## ---- PREREG 8.1 : exact conservation. Not "within 10x" -- equality. ----
  chk_counts <- sum(pb) == sum(umi_cell)
  chk_cells  <- sum(stats$n_cells) == ncol(cnt)
  if (!chk_counts)
    stop(sprintf("PREREG 8.1 FAILED (counts) %s/%s: sum(pseudobulk)=%.0f vs sum(cell UMI)=%.0f",
                 ds, smp, sum(pb), sum(umi_cell)))
  if (!chk_cells)
    stop(sprintf("PREREG 8.1 FAILED (cells) %s/%s: sum(n_cells)=%d vs filtered cells=%d",
                 ds, smp, sum(stats$n_cells), ncol(cnt)))
  stats[, `:=`(chk_counts_exact = chk_counts, chk_cells_exact = chk_cells,
               n_cells_annotated = n_ann, n_cells_after_3cond = n_after3,
               n_cells_after_4cond = n_after4,
               n_excluded_fine = n_after3 - n_after4, barcode_join_frac = frac_matched)]

  saveRDS(list(dataset = ds, sample = smp,
               genes = rownames(pb), bins = present,
               counts = pb,                      # genes x bins, RAW SUMMED COUNTS
               stats = stats,
               prereg = "PREREGISTRATION_pseudobulk_de.md sec 4.2/4.3/8.1",
               cell_filter = "in_ccc_graph & !high_error & hierarchy_bin %in% CCC_NODES & !(bmm_broad %in% CCC_EXCLUDE_FINE)",
               bin_source = "ANNO_RECONCILED_DIR"),
          out_rds)
  fwrite(stats, out_csv)

  message(sprintf("  [%s/%s] cells %d -> %d (3cond) -> %d (4cond, -%d Early Lymphoid) | bins %s | conservation OK",
                  ds, smp, n_ann, n_after3, n_after4, n_after3 - n_after4,
                  paste0(present, "(", stats$n_cells, ")", collapse = ",")))
  invisible(stats)
}

## ------------------------------------------------------------------ main ----
if (!is.na(opt$row)) {
  if (opt$row < 1L || opt$row > nrow(MAN)) stop("--row out of range 1..", nrow(MAN))
  r <- MAN[opt$row]
  if (!r$in_analysis) { message("[skip] row ", opt$row, " ", r$dataset, "/", r$sample,
                                " excluded: ", r$drop_reason); quit(save = "no", status = 0L) }
  aggregate_one(r$dataset, r$sample)
} else if (!is.na(opt$dataset) && !is.na(opt$sample)) {
  aggregate_one(opt$dataset, opt$sample)
} else {
  # Serial whole-cohort pass. Usable directly -- the array exists for wall-clock, not necessity.
  todo <- MAN[in_analysis == TRUE]
  message("[serial] ", nrow(todo), " samples")
  for (i in seq_len(nrow(todo))) aggregate_one(todo$dataset[i], todo$sample[i])
}
