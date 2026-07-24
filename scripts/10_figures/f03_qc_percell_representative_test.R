# =====================================================================
# f03_qc_percell_representative.R  --  detailed before/after for REP_SAMPLES
# ---------------------------------------------------------------------
# RAW (before) = merged per-dataset object <RAW_OBJ_DIR>/<Dataset>.rds,
# subset to one sample via the `Sample` column. FILTERED (after) =
# per-sample <FILTERED_OBJ_DIR>/<Dataset>/<sample>.rds.
# Per representative sample: (1) before/after violins of the 4 metrics;
# (2) nCount-vs-nFeature scatter colored by %mt, REMOVED cells (raw cells
# absent from the filtered object) highlighted.
# =====================================================================
.here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                  error = function(e) ".")
if (length(.here) == 0 || .here == "") .here <- "."
source(file.path(.here, "f00_fig_config.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

RAW_SAMPLE_COL <- "Sample"

find_filt <- function(ds, samp) {
  p <- file.path(FILTERED_OBJ_DIR, ds, paste0(samp, ".rds"))
  if (file.exists(p)) return(p)
  h <- list.files(FILTERED_OBJ_DIR, paste0("^", samp, "\\.rds$"), recursive = TRUE, full.names = TRUE)
  if (!length(h)) NA_character_ else h[1]
}

# per-cell metrics for a subset of a merged raw object (uses precomputed %mt/%hb if present)
raw_sample_metrics <- function(seu_r, cells, ds, sm) {
  md <- seu_r@meta.data[cells, , drop = FALSE]
  mt <- if ("percent.mt" %in% names(md)) md[["percent.mt"]] else pct_feature(seu_r, MT_PATTERN)[cells]
  hb <- if ("percent.hb" %in% names(md)) md[["percent.hb"]] else pct_feature(seu_r, HB_PATTERN)[cells]
  data.table(dataset = ds, sample = sm, stage = "before", cell = cells,
             nCount = as.numeric(md[["nCount_RNA"]]), nFeature = as.numeric(md[["nFeature_RNA"]]),
             pct_mt = as.numeric(mt), pct_hb = as.numeric(hb))
}

for (r in seq_len(nrow(REP_SAMPLES))) {
  ds <- REP_SAMPLES$dataset[r]; sm <- REP_SAMPLES$sample[r]
  tag <- paste0(ds, "__", sm)
  message(sprintf("[f03] %s", tag))
  rp <- file.path(RAW_OBJ_DIR, paste0(ds, ".rds")); fp <- find_filt(ds, sm)
  if (!file.exists(rp) || is.na(fp)) { message("   ! missing raw or filtered object -> skip"); next }

  seu_r <- readRDS(rp)
  if (!RAW_SAMPLE_COL %in% names(seu_r@meta.data)) { message("   ! Sample col missing -> skip"); rm(seu_r); gc(FALSE); next }
  cells_this <- colnames(seu_r)[as.character(seu_r@meta.data[[RAW_SAMPLE_COL]]) == sm]
  if (!length(cells_this)) { message(sprintf("   ! no raw cells with Sample==%s -> skip", sm)); rm(seu_r); gc(FALSE); next }
  mr <- raw_sample_metrics(seu_r, cells_this, ds, sm)
  rm(seu_r); gc(FALSE)

  seu_f <- readRDS(fp)
  mf   <- cell_metrics(seu_f, ds, sm, "after")
  kept <- colnames(seu_f)
  rm(seu_f); gc(FALSE)

  n_overlap <- length(intersect(cells_this, kept))
  if (n_overlap == 0L)
    message(sprintf("   ! WARNING: 0 barcode overlap raw<->filtered for %s (naming mismatch?) -> all shown as Removed", tag))

  ## -- (1) before/after violins ----
  long <- melt(rbind(mr[, .(stage, nCount, nFeature, pct_mt, pct_hb)],
                     mf[, .(stage, nCount, nFeature, pct_mt, pct_hb)]),
               id.vars = "stage", variable.name = "metric", value.name = "value")
  long[, stage := factor(stage, levels = c("before", "after"))]
  long[, metric := factor(metric, levels = c("nCount", "nFeature", "pct_mt", "pct_hb"),
                          labels = c("nCount", "nFeature", "% mito", "% hb"))]
  p_vio <- ggplot(long, aes(stage, value, fill = stage)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = STAGE_COLS, guide = "none") +
    labs(title = sprintf("%s  --  per-cell metrics before/after QC (n_before=%d, n_after=%d)",
                         tag, nrow(mr), length(kept)), x = NULL, y = NULL) +
    theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  ## -- (2) scatter colored by %mt, removed cells highlighted ----
  mr[, fate := fifelse(cell %in% kept, "Kept", "Removed")]
  p_sc <- ggplot(mr[order(fate)], aes(nCount, nFeature)) +
    geom_point(aes(color = pct_mt, shape = fate, alpha = fate), size = 0.5) +
    scale_x_log10(labels = scales::comma) + scale_y_log10(labels = scales::comma) +
    scale_color_viridis_c(option = "plasma", name = "% mito") +
    scale_shape_manual(values = c(Kept = 16, Removed = 4), name = NULL) +
    scale_alpha_manual(values = c(Kept = 0.35, Removed = 0.9), guide = "none") +
    labs(title = sprintf("%s  --  nCount vs nFeature (x = removed by QC)", tag),
         x = "nCount (log10)", y = "nFeature (log10)") +
    theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    save_fig(p_vio / p_sc, paste0("f03_rep_", tag), w = 11, h = 9, subdir = "01_qc")
  } else {
    save_fig(p_vio, paste0("f03_rep_", tag, "_violin"),  w = 11, h = 4, subdir = "01_qc")
    save_fig(p_sc,  paste0("f03_rep_", tag, "_scatter"), w = 8,  h = 6, subdir = "01_qc")
  }
  rm(mr, mf, long); gc(FALSE)
}
message("[f03] done.")