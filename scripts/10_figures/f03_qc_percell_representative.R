# =====================================================================
# f03_qc_percell_representative.R  --  per-sample detailed before/after QC
# ---------------------------------------------------------------------
# Runs over ALL PASS samples (one 2-panel figure each), grouped by dataset
# so each merged raw object <RAW_OBJ_DIR>/<Dataset>.rds is read only ONCE.
#   (1) before/after violins of the 4 metrics;
#   (2) nCount-vs-nFeature scatter colored by %mt, REMOVED cells (raw cells
#       absent from the filtered object) highlighted.
# Output: <FIG_DIR>/01_qc/persample/f03_persample_<Dataset>__<sample>.{pdf,png}
# Pilot with ONLY_DATASETS; cap scatter points with SCATTER_MAX.
# =====================================================================
.here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                  error = function(e) ".")
if (length(.here) == 0 || .here == "") .here <- "."
source(file.path(.here, "f00_fig_config.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

ONLY_DATASETS  <- NULL          # NULL = all PASS samples; or c("GSE227903") to pilot
RAW_SAMPLE_COL <- "Sample"
SCATTER_MAX    <- 40000L        # cap points drawn in the scatter (per sample)
OUT_SUBDIR     <- "01_qc/persample"

# per-cell metrics for a subset of a merged raw object (uses precomputed %mt/%hb if present)
raw_sample_metrics <- function(seu_r, cells, ds, sm) {
  md <- seu_r@meta.data[cells, , drop = FALSE]
  mt <- if ("percent.mt" %in% names(md)) md[["percent.mt"]] else pct_feature(seu_r, MT_PATTERN)[cells]
  hb <- if ("percent.hb" %in% names(md)) md[["percent.hb"]] else pct_feature(seu_r, HB_PATTERN)[cells]
  data.table(dataset = ds, sample = sm, stage = "before", cell = cells,
             nCount = as.numeric(md[["nCount_RNA"]]), nFeature = as.numeric(md[["nFeature_RNA"]]),
             pct_mt = as.numeric(mt), pct_hb = as.numeric(hb))
}

## -- filtered (PASS) per-sample index ----
filt <- data.table(path = list.files(FILTERED_OBJ_DIR, "\\.rds$", recursive = TRUE, full.names = TRUE))
stopifnot(nrow(filt) > 0)
filt[, sample  := sub("\\.rds$", "", basename(path))]
filt[, dataset := basename(dirname(path))]
if (!is.null(ONLY_DATASETS)) filt <- filt[dataset %in% ONLY_DATASETS]
stopifnot(nrow(filt) > 0)

by_ds <- split(seq_len(nrow(filt)), filt$dataset)
message(sprintf("[f03] %d sample(s) across %d dataset(s)", nrow(filt), length(by_ds)))
done <- 0L; n_warn <- 0L

for (ds in names(by_ds)) {
  rp <- file.path(RAW_OBJ_DIR, paste0(ds, ".rds"))
  if (!file.exists(rp)) { message(sprintf("[f03] ! no raw object for %s -> skip dataset", ds)); next }
  message(sprintf("[f03] %s  (loading merged raw once)", ds))
  seu_r <- readRDS(rp)
  if (!RAW_SAMPLE_COL %in% names(seu_r@meta.data)) {
    message(sprintf("   ! '%s' col missing in %s -> skip dataset", RAW_SAMPLE_COL, ds)); rm(seu_r); gc(FALSE); next
  }
  samp_vec <- as.character(seu_r@meta.data[[RAW_SAMPLE_COL]])

  for (i in by_ds[[ds]]) {
    sm <- filt$sample[i]; fp <- filt$path[i]; tag <- paste0(ds, "__", sm)
    cells_this <- colnames(seu_r)[samp_vec == sm]
    if (!length(cells_this)) { message(sprintf("   ! %s: no raw cells (Sample==%s) -> skip", tag, sm)); n_warn <- n_warn + 1L; next }

    mr   <- raw_sample_metrics(seu_r, cells_this, ds, sm)
    seu_f <- tryCatch(readRDS(fp), error = function(e) NULL)
    if (is.null(seu_f)) { message(sprintf("   ! %s: filtered read failed -> skip", tag)); n_warn <- n_warn + 1L; next }
    mf   <- cell_metrics(seu_f, ds, sm, "after")
    kept <- colnames(seu_f); rm(seu_f); gc(FALSE)

    if (length(intersect(cells_this, kept)) == 0L) {
      message(sprintf("   ! %s: 0 barcode overlap raw<->filtered (naming mismatch?) -> all Removed", tag)); n_warn <- n_warn + 1L
    }

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
      labs(title = sprintf("%s  --  before/after QC (n_before=%d, n_after=%d)",
                           tag, nrow(mr), length(kept)), x = NULL, y = NULL) +
      theme_qc() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

    ## -- (2) scatter colored by %mt, removed cells highlighted (downsampled) ----
    mr[, fate := fifelse(cell %in% kept, "Kept", "Removed")]
    mr_sc <- if (nrow(mr) > SCATTER_MAX) mr[sample(.N, SCATTER_MAX)] else mr
    p_sc <- ggplot(mr_sc[order(fate)], aes(nCount, nFeature)) +
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
      save_fig(p_vio / p_sc, paste0("f03_persample_", tag), w = 11, h = 9, subdir = OUT_SUBDIR)
    } else {
      save_fig(p_vio, paste0("f03_persample_", tag, "_violin"),  w = 11, h = 4, subdir = OUT_SUBDIR)
      save_fig(p_sc,  paste0("f03_persample_", tag, "_scatter"), w = 8,  h = 6, subdir = OUT_SUBDIR)
    }
    done <- done + 1L
    rm(mr, mr_sc, mf, long); gc(FALSE)
  }
  rm(seu_r); gc(FALSE)
}
message(sprintf("[f03] done. %d figure(s) written, %d sample(s) skipped/flagged.", done, n_warn))
