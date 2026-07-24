# =====================================================================
# f00_fig_config.R  --  shared config + helpers for QC / ref_norm figures
# ---------------------------------------------------------------------
# Sourced by f01..f04. Outputs PDF (vector, paper) + PNG (tracker) to FIG_DIR.
# =====================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## -- roots (override here if your global config names differ) ----
if (!exists("LARGE1_DIR")) LARGE1_DIR <- "/LARGE1/gr10634/gaozy/aml_niche_net"
if (!exists("FAST_DIR"))   FAST_DIR   <- "/FAST/gr10634/gaozy/aml_niche_net"

## -- inputs (USER-FILL where noted) ----
RAW_OBJ_DIR      <- file.path(LARGE1_DIR, "01_processed_counts", "rds")          # before QC
FILTERED_OBJ_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "01_per_sample_qc") # after QC
QC_REPORT_CSV    <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/01_qc/03_qc_report__ALL.csv"        # <-- POINT to your CSV
REFNORM_SUMMARY_CSV <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/03_malignancy/ref_norm_summary.csv"

## -- outputs ----
FIG_DIR <- file.path(FAST_DIR, "results", "figures")
dir.create(file.path(FIG_DIR, "01_qc"),      recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIG_DIR, "04_malignancy", "refnorm"), recursive = TRUE, showWarnings = FALSE)

## -- representative samples for the detailed per-cell figure (f03) ----
# dataset + sample (sample = RDS basename without .rds). Change freely.
REP_SAMPLES <- data.table(
  dataset = c("GSE227903",          "Chen2023"),
  sample  = c("5750_Dg",            "AML320_NICHE_IMMUNE")
)

## -- knobs ----
DOWNSAMPLE_PER_GROUP <- 20000L     # cells per (dataset, stage) for violin plotting only
MT_PATTERN <- "^MT[-.]"            # matches your QC pipeline
HB_PATTERN <- "^HB[^P]"
FIG_DPI    <- 300

## -- palette / theme ----
STAGE_COLS <- c(before = "#9aa7b1", after = "#2c7fb8")
ROUTE_COLS <- c(autologous = "#2c7fb8", external = "#f4a259", skip = "#9aa7b1")
theme_qc <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
}

## -- save helper: one PDF + one PNG ----
save_fig <- function(p, name, w = 10, h = 6, subdir = "01_qc") {
  outdir <- file.path(FIG_DIR, subdir)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)   # 自建,避免名字/scratch 问题
  base <- file.path(outdir, name)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = cairo_pdf, limitsize = FALSE)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = FIG_DPI, limitsize = FALSE)
  message(sprintf("[fig] %s.{pdf,png}", base))
}

## -- version-robust counts accessor + per-cell QC metrics from an object ----
.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

pct_feature <- function(seu, pattern) {
  cnt <- .get_counts(seu)
  g <- grep(pattern, rownames(cnt), value = TRUE)
  if (!length(g)) return(rep(0, ncol(cnt)))
  100 * Matrix::colSums(cnt[g, , drop = FALSE]) / pmax(Matrix::colSums(cnt), 1)
}

# Returns a data.table of per-cell QC metrics; computes %mt/%hb if absent.
cell_metrics <- function(seu, dataset, sample, stage) {
  md <- seu@meta.data
  nc <- md[["nCount_RNA"]];  nf <- md[["nFeature_RNA"]]
  mt <- if ("percent.mt" %in% names(md)) md[["percent.mt"]] else pct_feature(seu, MT_PATTERN)
  hb <- if ("percent.hb" %in% names(md)) md[["percent.hb"]] else pct_feature(seu, HB_PATTERN)
  data.table(dataset = dataset, sample = sample, stage = stage,
             cell = colnames(seu), nCount = as.numeric(nc), nFeature = as.numeric(nf),
             pct_mt = as.numeric(mt), pct_hb = as.numeric(hb))
}
