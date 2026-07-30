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
# Helpers (.get_counts / get_external_ref / .ext_ref_lineage_block / build_infercnv_input /
# run_one / burden_from_obj / infercnv_routes) are SHARED with 40_infercnv_run.R. They used to
# be duplicated here, which meant a fix applied to 40 silently left this array runner -- the one
# that actually produces the cohort -- on the old behaviour.
source(here::here("scripts", "02_malignancy", "00_infercnv_common.R"))
stopifnot(file.exists(INFERCNV_GENE_ORDER))

## -- prebuild mode: build the external ref cache and exit ----
if (opt$prebuild_ext_ref) { invisible(get_external_ref()); message("[prebuild] external ref ready."); quit(save = "no") }

stopifnot(!is.na(opt$dataset), !is.na(opt$sample))

## -- route this sample from the summary (shared derivation with 40) ----
summ <- infercnv_routes()
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

inp <- build_infercnv_input(s_cnt, route, opt$dataset, opt$sample)
counts <- inp$counts; anno <- inp$anno; ref_groups <- inp$ref_groups

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
