#!/usr/bin/env Rscript
# 41_infercnv_to_percell.R ----
# Turn a finished inferCNV run into a standardized per-cell malignancy table for 50 (consensus).
# PREFERRED input is the 40/44 burden CSV (<sample>_infercnv_burden.csv: cell, infercnv_burden,
# group) -- it already carries the reference cells, so we threshold the burden off the (CNV-flat)
# reference distribution. Falls back to the run.final.infercnv_obj two-metric method if no burden
# CSV is given.
#
# INPUT  : --burden_csv (preferred) OR --infercnv_dir with run.final.infercnv_obj (fallback)
# OUTPUT : <outdir>/<sample>__infercnv_percell.csv
#          schema shared by all arms:  cell, malignant(0/1), score, method, sample
#
# 3a-1 rewire (no logic change): here-anchored config source; --score_q default now reads
#   INFERCNV_SCORE_Q (== 0.95, unchanged) from config_malignancy.R; bare fwrite -> fwrite_safe.
#
# Usage (burden-CSV, preferred; run from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/41_infercnv_to_percell.R --sample 1216_Dg \
#     --burden_csv /LARGE1/.../05_cnv_snv/infercnv_burden/GSE227903/1216_Dg_infercnv_burden.csv \
#     --outdir /LARGE1/.../05_cnv_snv/infercnv/GSE227903/1216_Dg [--score_q 0.95]
# Usage (object fallback):
#   Rscript scripts/02_malignancy/41_infercnv_to_percell.R --sample S \
#     --infercnv_dir <dir-with-run.final.infercnv_obj> --outdir <dir>

suppressPackageStartupMessages({ library(optparse); library(data.table) })

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample",       type = "character"),
  make_option("--burden_csv",   type = "character", default = "", help = "40/44 burden CSV (preferred)"),
  make_option("--infercnv_dir", type = "character", default = "", help = "dir with run.final.infercnv_obj (fallback)"),
  make_option("--outdir",       type = "character"),
  make_option("--score_q",      type = "double", default = INFERCNV_SCORE_Q,
              help = "reference-cell burden quantile for threshold [default from config: %default]"),
  make_option("--cor_min",      type = "double", default = 0.20),
  make_option("--top_frac",     type = "double", default = 0.10)
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(opt$outdir, paste0(opt$sample, "__infercnv_percell.csv"))

## ---- PATH A: burden CSV (preferred) ----
if (nzchar(opt$burden_csv) && file.exists(opt$burden_csv)) {
  message("[A] burden CSV: ", opt$burden_csv)
  d <- fread(opt$burden_csv)                       # cell, infercnv_burden, group
  stopifnot(all(c("cell", "infercnv_burden", "group") %in% names(d)))
  is_ref <- grepl("reference", d$group, ignore.case = TRUE)
  if (sum(is_ref) < 10)
    warning("few reference cells (", sum(is_ref), "); threshold may be unstable")
  thr <- as.numeric(quantile(d$infercnv_burden[is_ref], opt$score_q, na.rm = TRUE))
  message("    burden threshold (ref q", opt$score_q, ") = ", signif(thr, 3))
  d[, malignant := as.integer(infercnv_burden > thr)]
  d[is_ref, malignant := 0L]                       # reference cells normal by design
  # FOREIGN reference cells are not this sample's cells -> drop them from output. Match on "any
  # reference group that is not reference_normal" rather than on the reference_external prefix:
  # 00_infercnv_common now also emits reference_matched__<bin> (this dataset's healthy donors),
  # and a prefix test would have silently passed those foreign cells through as sample cells.
  d <- d[!(is_ref & group != "reference_normal")]
  res <- d[, .(cell, malignant, score = infercnv_burden, method = "infercnv", sample = opt$sample)]
  fwrite_safe(res, out)
  message("[done] ", sum(res$malignant), "/", nrow(res), " malignant -> ", out); quit(save = "no")
}

## ---- PATH B: object fallback (two-metric: burden + correlation) ----
message("[B] no burden CSV; falling back to run.final.infercnv_obj in ", opt$infercnv_dir)
obj_path <- file.path(opt$infercnv_dir, "run.final.infercnv_obj")
if (!file.exists(obj_path)) {
  cands <- list.files(opt$infercnv_dir, pattern = "infercnv_obj$", recursive = TRUE, full.names = TRUE)
  if (length(cands)) obj_path <- cands[1] else stop("no burden CSV and no infercnv_obj found")
}
icnv <- readRDS(obj_path)
expr <- icnv@expr.data; cells <- colnames(expr)
ref_idx <- tryCatch(unlist(icnv@reference_grouped_cell_indices), error = function(e) integer(0))
obs_idx <- tryCatch(unlist(icnv@observation_grouped_cell_indices), error = function(e) integer(0))
if (!length(obs_idx)) obs_idx <- setdiff(seq_along(cells), ref_idx)
cnv_score <- colMeans((expr - 1)^2)
ord <- obs_idx[order(cnv_score[obs_idx], decreasing = TRUE)]
top <- head(ord, max(20L, ceiling(opt$top_frac * length(obs_idx))))
mp <- rowMeans(expr[, top, drop = FALSE])
xc <- expr - rowMeans(expr); pc <- mp - mean(mp)
cnv_cor <- as.numeric(crossprod(pc, xc)) / sqrt(sum(pc^2) * colSums(xc^2)); cnv_cor[is.na(cnv_cor)] <- 0
if (length(ref_idx) >= 20) {
  thr_s <- as.numeric(quantile(cnv_score[ref_idx], opt$score_q, na.rm = TRUE))
  thr_c <- max(opt$cor_min, as.numeric(quantile(cnv_cor[ref_idx], opt$score_q, na.rm = TRUE)))
} else { thr_s <- as.numeric(quantile(cnv_score[obs_idx], 0.75, na.rm = TRUE)); thr_c <- opt$cor_min }
mal <- as.integer(cnv_score > thr_s & cnv_cor > thr_c); mal[ref_idx] <- 0L
res <- data.table(cell = cells, malignant = mal, score = cnv_cor, method = "infercnv", sample = opt$sample)
res <- res[setdiff(seq_along(cells), ref_idx)]    # drop reference cells if external
fwrite_safe(res, out)
message("[done] ", sum(res$malignant), "/", nrow(res), " malignant -> ", out)
