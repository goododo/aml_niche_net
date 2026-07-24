#!/usr/bin/env Rscript
# 43_author_to_percell.R ----
# Extract an AUTHOR-provided malignancy/CNV call (e.g. CopyKAT aneuploid/diploid, SCEVAN,
# or any per-cell label already in the processed object) from a per-sample QC Seurat RDS
# into the standardized per-cell schema. No raw FASTQ needed; barcodes are exactly the QC
# object's, so they align trivially in 50.
#
# INPUT  : --qc_rds (per-sample QC Seurat .rds), --col (author metadata column)
# OUTPUT : <outdir>/<sample>__author_percell.csv
#          schema shared by all arms:  cell, malignant(0/1/NA), score, method, sample
#
# 3a-1 rewire (no logic change): here-anchored config source; bare fwrite -> fwrite_safe.
#
# Usage (run from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/43_author_to_percell.R --sample PT10A \
#     --qc_rds /LARGE1/.../01_per_sample_qc/GSE239721/PT10A.rds \
#     --outdir /LARGE1/.../03_cnv_snv/author/GSE239721/PT10A \
#     [--col copykat] [--malignant_values aneuploid] [--normal_values diploid]
#
# Defaults target CopyKAT. For other callers pass --col / --malignant_values / --normal_values.

suppressPackageStartupMessages({ library(optparse); library(data.table) })

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample", type = "character"),
  make_option("--qc_rds", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--col",    type = "character", default = "copykat",
              help = "metadata column holding the author call [copykat]"),
  make_option("--malignant_values", type = "character", default = "aneuploid",
              help = "comma-separated values meaning malignant [aneuploid]"),
  make_option("--normal_values",    type = "character", default = "diploid",
              help = "comma-separated values meaning normal [diploid]")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

message("[1] Loading QC object: ", opt$qc_rds)
obj <- readRDS(opt$qc_rds)
md  <- tryCatch(obj@meta.data, error = function(e) as.data.frame(obj[[]]))
cells <- rownames(md)

# resolve the column (case-insensitive; also try an anno_ prefixed copy)
cands <- c(opt$col, paste0("anno_", opt$col))
hit <- intersect(tolower(cands), tolower(colnames(md)))
if (!length(hit)) stop("column '", opt$col, "' (or anno_", opt$col, ") not found. Available: ",
                       paste(head(colnames(md), 60), collapse = ", "))
col <- colnames(md)[match(hit[1], tolower(colnames(md)))]
message("[2] Using author column: ", col)

vals <- tolower(as.character(md[[col]]))
mal  <- tolower(trimws(strsplit(opt$malignant_values, ",")[[1]]))
nor  <- tolower(trimws(strsplit(opt$normal_values,    ",")[[1]]))

malignant <- rep(NA_integer_, length(vals))
malignant[vals %in% mal] <- 1L
malignant[vals %in% nor] <- 0L

dt <- data.table(cell = cells, malignant = malignant, score = NA_real_,
                 method = paste0("author_", opt$col), sample = opt$sample)
out <- file.path(opt$outdir, paste0(opt$sample, "__author_percell.csv"))
fwrite_safe(dt, out)
message("[done] malignant=", sum(malignant == 1, na.rm = TRUE),
        " normal=", sum(malignant == 0, na.rm = TRUE),
        " NA=", sum(is.na(malignant)), " -> ", out)
