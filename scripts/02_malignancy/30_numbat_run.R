#!/usr/bin/env Rscript
# 30_numbat_run.R ----
# Numbat CNV/clone inference for one sample (cluster). Uses a custom expression reference RDS if
# provided (--ref_rds), else the built-in ref_hca. Handles the "no CNV detected" case without
# crashing (writes a degraded stub whose note == NUMBAT_DEGRADED_NOTE, which 50 reads as ABSTAIN).
#
# INPUT  : --matrix (STARsolo Gene/filtered dir), --allele (allele_counts.tsv.gz)
# OUTPUT : <outdir>/<sample>__numbat_percell.csv
#            valid run  : clone_post columns (incl. compartment_opt) + sample
#            degraded    : cell, sample, p_cnv=NA, note=NUMBAT_DEGRADED_NOTE
#
# 3a-1 rewire (value-preserving): here-anchored config source; --max_entropy / --t defaults now
#   read NUMBAT_MAX_ENTROPY (0.5, [D-A7]) / NUMBAT_T (1e-5) from config; the degraded-stub note
#   literal is now the shared NUMBAT_DEGRADED_NOTE constant (couples this writer to the 50 reader
#   so they cannot drift); bare fwrite -> fwrite_safe. Run logic unchanged.

suppressPackageStartupMessages({
  library(optparse); library(numbat); library(Matrix); library(data.table)
})

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample",  type = "character"),
  make_option("--matrix",  type = "character", help = "STARsolo Gene/filtered dir"),
  make_option("--allele",  type = "character", help = "allele_counts.tsv.gz"),
  make_option("--outdir",  type = "character"),
  make_option("--genome",  type = "character", default = "hg38"),
  make_option("--ncores",  type = "integer",   default = 16L),
  make_option("--ref_rds", type = "character", default = ""),  # optional custom ref
  make_option("--max_entropy", type = "double", default = NUMBAT_MAX_ENTROPY,
              help = "Numbat single-cell CNV entropy filter [default from config: %default; D-A7]"),
  make_option("--t", type = "double", default = NUMBAT_T,
              help = "Numbat HMM transition prob [default from config: %default]")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Read raw counts (genes x cells) from STARsolo filtered output ----
message("[1] Reading STARsolo filtered matrix: ", opt$matrix)
read_star_filtered <- function(dir) {
  f <- function(b) { p <- file.path(dir, b); pg <- paste0(p, ".gz")
    if (file.exists(p)) p else if (file.exists(pg)) pg else stop("missing ", b, " in ", dir) }
  m  <- as(Matrix::readMM(f("matrix.mtx")), "CsparseMatrix")
  bc <- fread(f("barcodes.tsv"), header = FALSE)[[1]]
  ft <- fread(f("features.tsv"), header = FALSE)
  gene <- if (ncol(ft) >= 2) ft[[2]] else ft[[1]]
  rownames(m) <- make.unique(gene); colnames(m) <- bc; m
}
count_mat <- read_star_filtered(opt$matrix)
message("    matrix: ", nrow(count_mat), " genes x ", ncol(count_mat), " cells")

## ---- 2. Allele counts ----
message("[2] Reading allele counts: ", opt$allele)
df_allele <- fread(opt$allele)

## ---- 3. Expression reference: custom RDS if given, else ref_hca ----
lambdas_ref <- NULL
if (nzchar(opt$ref_rds) && file.exists(opt$ref_rds)) {
  lambdas_ref <- tryCatch(readRDS(opt$ref_rds), error = function(e) NULL)
  if (!is.null(lambdas_ref)) message("[3] Using custom reference: ", opt$ref_rds)
}
if (is.null(lambdas_ref)) {
  message("[3] Using built-in ref_hca")
  data(ref_hca, package = "numbat"); lambdas_ref <- ref_hca
}

## ---- 4. Run Numbat ----
message("[4] run_numbat() ...  (t=", opt$t, ", max_entropy=", opt$max_entropy, ")")
out <- run_numbat(
  count_mat = count_mat, lambdas_ref = lambdas_ref, df_allele = df_allele,
  genome = opt$genome, t = opt$t, max_entropy = opt$max_entropy,
  ncores = opt$ncores, plot = TRUE, out_dir = opt$outdir
)

## ---- 5. Tidy per-cell table (graceful on no-CNV) ----
message("[5] Writing per-cell CNV/clone table")
out_csv <- file.path(opt$outdir, paste0(opt$sample, "__numbat_percell.csv"))
nb    <- tryCatch(Numbat$new(out_dir = opt$outdir), error = function(e) NULL)
clone <- if (!is.null(nb)) tryCatch(nb$clone_post, error = function(e) NULL) else NULL

if (!is.null(clone) && is.data.frame(clone) && nrow(clone) > 0) {
  dt <- as.data.table(clone); dt[, sample := opt$sample]
  fwrite_safe(dt, out_csv)
  message("    wrote per-cell clone table: ", nrow(dt), " cells -> ", out_csv)
} else {
  stub <- data.table(cell = colnames(count_mat), sample = opt$sample,
                     p_cnv = NA_real_, note = NUMBAT_DEGRADED_NOTE)
  fwrite_safe(stub, out_csv)
  message("    NOTE: no CNV detected; wrote stub -> ", out_csv,
          " (sample contributes no Numbat CNV evidence).")
}
message("[done] Numbat outputs -> ", opt$outdir)
