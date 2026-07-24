# =====================================================================
# f02_qc_percell_cohort.R  --  per-cell before/after QC violins, all datasets
# ---------------------------------------------------------------------
# RAW (before) objects are ONE merged Seurat object per dataset:
#   <RAW_OBJ_DIR>/<Dataset>.rds  with a `Sample` column (= filtered basename)
#   and nCount_RNA/nFeature_RNA/percent.mt/percent.hb already computed.
# FILTERED (after) objects are per-sample: <FILTERED_OBJ_DIR>/<Dataset>/<sample>.rds
#
# Strategy: read each dataset's merged raw ONCE (split by Sample) for the
# `before` side; read per-sample filtered for the `after` side. Restrict
# `before` to PASS samples (those with a filtered file) for apples-to-apples.
# Downsample per (dataset,stage) for plotting only. Pilot via ONLY_DATASETS.
# =====================================================================
.here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                  error = function(e) ".")
if (length(.here) == 0 || .here == "") .here <- "."
source(file.path(.here, "f00_fig_config.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

ONLY_DATASETS  <- NULL   # NULL = all, test sets = c("GSE227903", "Chen2023")
RAW_SAMPLE_COL <- "Sample"                      # per-cell sample id in the merged raw object

## -- filtered (PASS) per-sample index: dataset = parent dir, sample = basename ----
filt <- data.table(path = list.files(FILTERED_OBJ_DIR, "\\.rds$", recursive = TRUE, full.names = TRUE))
stopifnot(nrow(filt) > 0)
filt[, sample  := sub("\\.rds$", "", basename(path))]
filt[, dataset := basename(dirname(path))]
if (!is.null(ONLY_DATASETS)) filt <- filt[dataset %in% ONLY_DATASETS]
stopifnot(nrow(filt) > 0)

acc <- list(); k <- 0L

## -- AFTER: per-sample filtered objects ----
for (i in seq_len(nrow(filt))) {
  ds <- filt$dataset[i]; sm <- filt$sample[i]
  message(sprintf("[f02][after] (%d/%d) %s / %s", i, nrow(filt), ds, sm))
  seu <- tryCatch(readRDS(filt$path[i]), error = function(e) NULL)
  if (!is.null(seu)) { k <- k + 1L; acc[[k]] <- cell_metrics(seu, ds, sm, "after") }
  rm(seu); gc(FALSE)
}

## -- BEFORE: one merged raw object per dataset, restricted to PASS samples ----
pass_by_ds <- split(filt$sample, filt$dataset)
for (ds in names(pass_by_ds)) {
  rp <- file.path(RAW_OBJ_DIR, paste0(ds, ".rds"))
  if (!file.exists(rp)) { message(sprintf("[f02][before] ! no raw object for %s -> skip", ds)); next }
  message(sprintf("[f02][before] %s (merged raw)", ds))
  seu <- readRDS(rp)
  md  <- seu@meta.data
  if (!RAW_SAMPLE_COL %in% names(md)) { message(sprintf("   ! '%s' col missing in %s -> skip", RAW_SAMPLE_COL, ds)); rm(seu); gc(FALSE); next }
  mt <- if ("percent.mt" %in% names(md)) md[["percent.mt"]] else pct_feature(seu, MT_PATTERN)
  hb <- if ("percent.hb" %in% names(md)) md[["percent.hb"]] else pct_feature(seu, HB_PATTERN)
  bdt <- data.table(dataset = ds, sample = as.character(md[[RAW_SAMPLE_COL]]), stage = "before",
                    cell = colnames(seu),
                    nCount = as.numeric(md[["nCount_RNA"]]), nFeature = as.numeric(md[["nFeature_RNA"]]),
                    pct_mt = as.numeric(mt), pct_hb = as.numeric(hb))
  bdt <- bdt[sample %in% pass_by_ds[[ds]]]    # parity with the `after` (PASS) set
  k <- k + 1L; acc[[k]] <- bdt
  rm(seu, bdt); gc(FALSE)
}

dt <- rbindlist(acc[seq_len(k)], use.names = TRUE)
dt[, stage := factor(stage, levels = c("before", "after"))]
dt[, Dataset := factor(dataset, levels = sort(unique(dataset)))]

## -- downsample per (dataset, stage) for plotting ----
set.seed(20260605)
plt <- dt[, .SD[if (.N > DOWNSAMPLE_PER_GROUP) sample(.N, DOWNSAMPLE_PER_GROUP) else seq_len(.N)],
          by = .(Dataset, stage)]

## -- long form ----
long <- melt(plt, id.vars = c("Dataset", "stage"),
             measure.vars = c("nCount", "nFeature", "pct_mt", "pct_hb"),
             variable.name = "metric", value.name = "value")
long[, metric := factor(metric, levels = c("nCount", "nFeature", "pct_mt", "pct_hb"),
                        labels = c("nCount (UMI)", "nFeature (genes)", "% mito", "% hemoglobin"))]

mk <- function(m, logy) {
  g <- ggplot(long[metric == m], aes(Dataset, value, fill = stage)) +
    geom_violin(scale = "width", position = position_dodge(0.8), linewidth = 0.2, alpha = 0.9) +
    scale_fill_manual(values = STAGE_COLS) +
    labs(title = m, x = NULL, y = NULL, fill = NULL) + theme_qc()
  if (logy) g <- g + scale_y_log10(labels = scales::comma)
  g
}
ps <- list(mk("nCount (UMI)", TRUE), mk("nFeature (genes)", TRUE),
           mk("% mito", FALSE), mk("% hemoglobin", FALSE))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- Reduce(`/`, ps) +
    plot_annotation(title = "AML niche_net  --  per-cell QC, before vs after (cohort)")
  save_fig(combo, "f02_qc_percell_cohort", w = 12, h = 16, subdir = "01_qc")
} else {
  for (j in seq_along(ps)) save_fig(ps[[j]], sprintf("f02_qc_percell_cohort_%d", j), w = 12, h = 4, subdir = "01_qc")
}
message("[f02] done.")