#!/usr/bin/env Rscript
# 42_exprcnv_run.R ----
# Run an EXPRESSION-CNV malignancy caller (CopyKAT default, or SCEVAN) directly on a per-sample
# QC Seurat RDS -- no FASTQ needed. Use this only as the FALLBACK second expression arm, i.e.
# when a sample has no FASTQ/allele/SNV path AND no author CNV call, so inferCNV alone would be
# a single-type (Tier C) consensus. Adding this lifts it to expression-multi (Tier B).
#
# NOTE: this is still EXPRESSION-CNV (correlated with inferCNV). It does NOT make the consensus
# orthogonal (Tier A); 50 collapses all expression callers into one type-level vote.
#
# INPUT  : --qc_rds (per-sample QC Seurat .rds), --method {copykat|scevan}
# OUTPUT : <outdir>/<sample>__<method>_percell.csv
#          schema shared by all arms:  cell, malignant(0/1/NA), score, method, sample
#
# 3a-1 rewire (no logic change): here-anchored config source; bare fwrite -> fwrite_safe.
#   NOTE: inline counts extraction (GetAssayData slot="counts") deliberately kept as-is this
#   batch; get_counts() adoption (v5/v4 layer/slot) is deferred to the later utils pass because
#   it changes the extraction path (and copykat needs a dense as.matrix).
#
# Usage (run from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/42_exprcnv_run.R --sample PT10A --method copykat \
#     --qc_rds .../PT10A.rds --outdir .../03_cnv_snv/copykat/GSE239721/PT10A \
#     [--normal_col anno_SingleR] [--normal_values T,NK,CD4,CD8] [--ncores 8]

suppressPackageStartupMessages({ library(optparse); library(Matrix); library(data.table) })

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample",        type = "character"),
  make_option("--qc_rds",        type = "character"),
  make_option("--outdir",        type = "character"),
  make_option("--method",        type = "character", default = "copykat",
              help = "copykat | scevan [copykat]"),
  make_option("--normal_col",    type = "character", default = "",
              help = "metadata column naming a normal reference cell type (optional)"),
  make_option("--normal_values", type = "character", default = "",
              help = "comma-separated normal cell-type values to anchor the reference"),
  make_option("--ncores",        type = "integer", default = 8L)
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
owd <- setwd(opt$outdir)   # CopyKAT writes outputs to the working dir

message("[1] Loading QC object: ", opt$sample)
obj <- readRDS(opt$qc_rds)
counts <- tryCatch(as.matrix(SeuratObject::GetAssayData(obj, slot = "counts")),
                   error = function(e) as.matrix(obj@assays$RNA@counts))
md <- tryCatch(obj@meta.data, error = function(e) as.data.frame(obj[[]]))
message("    ", nrow(counts), " genes x ", ncol(counts), " cells")

# optional normal-cell anchor (helps CopyKAT define the diploid baseline)
norm_cells <- character(0)
if (nzchar(opt$normal_col) && opt$normal_col %in% colnames(md) && nzchar(opt$normal_values)) {
  nv <- trimws(strsplit(opt$normal_values, ",")[[1]])
  norm_cells <- rownames(md)[as.character(md[[opt$normal_col]]) %in% nv]
  message("    normal anchor: ", length(norm_cells), " cells from ", opt$normal_col)
}

if (tolower(opt$method) == "copykat") {
  suppressPackageStartupMessages(library(copykat))
  message("[2] Running CopyKAT")
  res <- copykat(rawmat = counts, id.type = "S", ngene.chr = 5, win.size = 25,
                 KS.cut = 0.1, sam.name = opt$sample, distance = "euclidean",
                 norm.cell.names = norm_cells, output.seg = FALSE,
                 plot.genes = FALSE, genome = "hg20", n.cores = opt$ncores)
  pred <- as.data.table(res$prediction)            # cell.names, copykat.pred
  setnames(pred, c("cell.names", "copykat.pred"), c("cell", "pred"), skip_absent = TRUE)
  pred[, malignant := fifelse(grepl("aneuploid", pred, ignore.case = TRUE), 1L,
                       fifelse(grepl("diploid",  pred, ignore.case = TRUE), 0L, NA_integer_))]
  dt <- pred[, .(cell, malignant, score = NA_real_, method = "copykat", sample = opt$sample)]

} else if (tolower(opt$method) == "scevan") {
  suppressPackageStartupMessages(library(SCEVAN))
  message("[2] Running SCEVAN")
  res <- SCEVAN::pipelineCNA(count_mtx = counts, sample = opt$sample,
                             par_cores = opt$ncores, SUBCLONES = FALSE)
  pred <- as.data.table(res, keep.rownames = "cell")   # has 'class' = tumor/normal/filtered
  cls <- intersect(c("class", "pred"), names(pred))[1]
  pred[, malignant := fifelse(grepl("tumor",  get(cls), ignore.case = TRUE), 1L,
                       fifelse(grepl("normal", get(cls), ignore.case = TRUE), 0L, NA_integer_))]
  dt <- pred[, .(cell, malignant, score = NA_real_, method = "scevan", sample = opt$sample)]
} else stop("unknown --method: ", opt$method)

setwd(owd)
out <- file.path(opt$outdir, paste0(opt$sample, "__", tolower(opt$method), "_percell.csv"))
fwrite_safe(dt, out)
message("[done] ", opt$method, ": malignant=", sum(dt$malignant == 1, na.rm = TRUE),
        " normal=", sum(dt$malignant == 0, na.rm = TRUE),
        " NA=", sum(is.na(dt$malignant)), " -> ", out)
