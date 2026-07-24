#!/usr/bin/env Rscript
# 44_infercnv_run_one.R ----
# SINGLE-SAMPLE inferCNV, array-friendly. Logic is IDENTICAL to 40 (same routing from
# ref_norm_summary, same gene_order, same inferCNV params, same burden output) so results are
# comparable to the already-finished runs. The only change vs 40 is that the per-sample work is
# parameterized (--dataset/--sample) instead of looped. Submitted by 45 as a SLURM array.
#
# INPUT  : REFNORM_SUMMARY_CSV (routing) ; per-sample QC rds under QC_RDS_DIR
# OUTPUT : INFERCNV_BURDEN_ROOT/<ds>/<sample>_infercnv_burden.csv (+ full inferCNV dir)
#
# 3a-2 rewire (no logic change): source -> here 3-line config; all inferCNV params + external-ref
#   constants now come from config_malignancy.R (removed the script-local block); output roots use
#   the canonical INFERCNV_ROOT / INFERCNV_BURDEN_ROOT, QC dir is QC_RDS_DIR, external-ref seed is
#   INFERCNV_EXT_REF_SEED (== 20260613), BMM ref is INFERCNV_EXT_REF_RDS(=BMM_ANNOTATED_RDS); bare
#   fwrite -> fwrite_safe. Inline .get_counts kept (get_counts adoption deferred).
#
# Run the external-reference cache build ONCE before launching the array (external samples share
# it; otherwise parallel tasks race to build it):
#   Rscript scripts/02_malignancy/44_infercnv_run_one.R --prebuild_ext_ref
# Then per sample (via 45 array):
#   Rscript scripts/02_malignancy/44_infercnv_run_one.R --dataset GSE227903 --sample 1216_Dg

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(infercnv); library(data.table); library(optparse)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = NA),
  make_option("--sample",  type = "character", default = NA),
  make_option("--prebuild_ext_ref", action = "store_true", default = FALSE)
)))

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
stopifnot(file.exists(INFERCNV_GENE_ORDER))

.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

get_external_ref <- function() {
  if (file.exists(INFERCNV_EXT_REF_CACHE)) return(readRDS(INFERCNV_EXT_REF_CACHE))
  message("[ext-ref] building BoneMarrowMap subsample (one-time) ...")
  stopifnot(file.exists(INFERCNV_EXT_REF_RDS))
  bbm <- readRDS(INFERCNV_EXT_REF_RDS)
  lab <- as.character(bbm@meta.data[[INFERCNV_EXT_REF_LABELCOL]])
  set.seed(INFERCNV_EXT_REF_SEED)
  idx <- unlist(lapply(split(seq_along(lab), lab), function(ii)
    if (length(ii) <= INFERCNV_EXT_REF_PER_TYPE) ii else sample(ii, INFERCNV_EXT_REF_PER_TYPE)))
  cnt <- .get_counts(bbm)[, idx, drop = FALSE]
  colnames(cnt) <- paste0("EXTREF_", colnames(cnt))
  ext <- list(counts = cnt, cells = colnames(cnt))
  dir.create(dirname(INFERCNV_EXT_REF_CACHE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ext, INFERCNV_EXT_REF_CACHE)
  rm(bbm); gc(verbose = FALSE)
  message(sprintf("[ext-ref] %d external reference cells cached.", ncol(cnt)))
  ext
}

run_one <- function(counts, anno_dt, ref_groups, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  anno_file <- file.path(out_dir, "annotations.tsv")
  fwrite_safe(anno_dt, anno_file, sep = "\t", col.names = FALSE)
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts, annotations_file = anno_file,
    gene_order_file = INFERCNV_GENE_ORDER, ref_group_names = ref_groups)
  infercnv::run(obj, cutoff = INFERCNV_CUTOFF, out_dir = out_dir,
    cluster_by_groups = TRUE, analysis_mode = INFERCNV_ANALYSIS,
    denoise = INFERCNV_DENOISE, HMM = INFERCNV_HMM,
    num_threads = INFERCNV_THREADS, no_plot = FALSE, save_rds = TRUE)
}
burden_from_obj <- function(obj) colMeans((obj@expr.data - 1)^2)

## -- prebuild mode: build the external ref cache and exit ----
if (opt$prebuild_ext_ref) { invisible(get_external_ref()); message("[prebuild] external ref ready."); quit(save = "no") }

stopifnot(!is.na(opt$dataset), !is.na(opt$sample))

## -- route this sample from the summary (same derivation as 40) ----
stopifnot(file.exists(REFNORM_SUMMARY_CSV))
summ <- fread(REFNORM_SUMMARY_CSV)
summ[, route := fifelse(decision == "autologous_ok", "autologous",
              fifelse(method == "sorted_guard", "skip", "external"))]
row <- summ[dataset == opt$dataset & sample_id == opt$sample]
if (!nrow(row)) stop("sample not in summary: ", opt$dataset, "::", opt$sample)
route <- row$route[1]

out_dir    <- file.path(INFERCNV_ROOT, opt$dataset, opt$sample)
burden_csv <- file.path(INFERCNV_BURDEN_ROOT, opt$dataset, paste0(opt$sample, "_infercnv_burden.csv"))
if (file.exists(burden_csv)) { message("[skip] already done: ", burden_csv); quit(save = "no") }
if (route == "skip") { message("[skip] sorted_guard: ", opt$sample); quit(save = "no") }

rds_in <- file.path(QC_RDS_DIR, opt$dataset, paste0(opt$sample, ".rds"))
stopifnot(file.exists(rds_in))
message(sprintf("[run] %s::%s route=%s", opt$dataset, opt$sample, route))
seu <- readRDS(rds_in); s_cnt <- .get_counts(seu)

if (route == "autologous") {
  ref_txt <- file.path(REFNORM_REF_CELL_DIR, opt$dataset, paste0(opt$sample, "_ref_norm_cells.txt"))
  stopifnot(file.exists(ref_txt))
  ref_cells <- intersect(readLines(ref_txt), colnames(s_cnt))
  grp <- ifelse(colnames(s_cnt) %in% ref_cells, "reference_normal", "observation")
  counts <- s_cnt
  anno   <- data.table(cell = colnames(s_cnt), group = grp)
  ref_groups <- "reference_normal"
} else {
  ext <- get_external_ref()
  common <- intersect(rownames(s_cnt), rownames(ext$counts)); stopifnot(length(common) > 1000)
  counts <- cbind(s_cnt[common, , drop = FALSE], ext$counts[common, , drop = FALSE])
  anno   <- data.table(cell = colnames(counts),
                       group = c(rep("observation", ncol(s_cnt)),
                                 rep("reference_external", ncol(ext$counts))))
  ref_groups <- "reference_external"
}

obj <- {
  n_obs <- sum(anno$group == "observation")
  if (n_obs < INFERCNV_MIN_OBS) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(sprintf("skipped: %d observation cell(s) (< INFERCNV_MIN_OBS=%d); near-all-normal / degenerate for inferCNV",
                       n_obs, INFERCNV_MIN_OBS), file.path(out_dir, "lowobs.skip"))
    message(sprintf("[skip] %s: only %d observation cell(s) -> inferCNV not applicable", opt$sample, n_obs))
    quit(save = "no")
  }
  run_one(counts, anno, ref_groups, out_dir)
}
burden <- burden_from_obj(obj)
bdt <- data.table(cell = names(burden), infercnv_burden = as.numeric(burden))
bdt[, group := anno$group[match(cell, anno$cell)]]
dir.create(dirname(burden_csv), recursive = TRUE, showWarnings = FALSE)
fwrite_safe(bdt, burden_csv)
message("[done] ", opt$sample, " -> ", burden_csv)
