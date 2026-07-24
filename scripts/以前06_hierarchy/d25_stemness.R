#!/usr/bin/env Rscript
# d25_stemness.R ----
# Phase 2 step 2.5: per-cell stemness scores from versioned signatures (stemness_signatures.tsv).
# Weighted signatures (all genes have weights, e.g. LSC17) -> z-score each gene across the sample's
# cells, then weighted sum (the LSC17 convention). Unweighted signatures -> Seurat AddModuleScore.
# Sidecar output (does NOT touch d10/d20); d40 joins it for per-bin mean stemness.
#
# NOTE: scores are scaled WITHIN each sample, so they rank cells within a sample; cross-sample
# comparison should use rank-normalisation downstream (consistent with the rank-based philosophy).
#
# Per-sample, serial, resume-skip.
# Usage:
#   Rscript d25_stemness.R --dataset=GSE227903
#   Rscript d25_stemness.R --dataset=GSE239721 --sample=PT10A

suppressPackageStartupMessages({ library(optparse); library(Seurat); library(SeuratObject); library(data.table) })
source("d00_config.R")

PERCELL_STEM_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_stem")
SIG_TSV          <- file.path(SCRIPT_DIR, "stemness_signatures.tsv")

# Old<->new symbol aliases (try both against the object).
ALIASES <- c(KIAA0125 = "FAM30A", FAM30A = "KIAA0125", NGFRAP1 = "BEX3", BEX3 = "NGFRAP1",
             GPR56 = "ADGRG1", ADGRG1 = "GPR56", SMIM24 = "C19orf77", C19orf77 = "SMIM24")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character"),
  make_option("--sample",  type = "character", default = NULL)
)))
if (is.null(opt$dataset)) stop("--dataset is required")

if (!file.exists(SIG_TSV)) stop("[d25] signature file not found: ", SIG_TSV)
sig_lines <- readLines(SIG_TSV)
sig_lines <- sig_lines[!grepl("^\\s*#", sig_lines) & nzchar(trimws(sig_lines))]   # drop comments/blanks
sig <- fread(text = paste(sig_lines, collapse = "\n"), sep = "\t", header = TRUE)   # signature, gene, weight
sig <- sig[!is.na(gene) & gene != "" & !startsWith(gene, "#")]
sigs <- split(sig, sig$signature)
message("[1] signatures: ", paste(names(sigs), collapse = ", "))

resolve_gene <- function(g, present) {
  if (g %in% present) return(g)
  if (g %in% names(ALIASES)) {
    a <- ALIASES[[g]]
    if (a %in% present) return(a)
  }
  NA_character_
}

score_one <- function(rds, ds, samp) {
  out <- file.path(PERCELL_STEM_DIR, ds, paste0(samp, "_stem.csv"))
  if (file.exists(out)) { message("    [skip] ", samp); return(invisible()) }
  obj <- readRDS(rds)
  # The QC RDS often has an EMPTY 'data' layer that NormalizeData() does not reliably backfill
  # under some Seurat v5 layer structures (observed: GSE227903 -> empty data -> garbage z-scores).
  # So compute log-CP10k directly from counts; independent of any data-layer state.
  cnt <- get_counts(obj)
  libsize <- Matrix::colSums(cnt); libsize[libsize == 0] <- 1
  expr <- cnt
  expr@x <- log1p(expr@x / rep(libsize, diff(expr@p)) * 1e4)   # log1p(CP10k), per the LogNormalize convention
  present <- rownames(expr)
  res <- data.table(cell = colnames(obj))

  for (sn in names(sigs)) {
    s  <- sigs[[sn]]
    gg <- vapply(s$gene, resolve_gene, character(1), present = present)
    ok <- !is.na(gg)
    message("    [", samp, "] ", sn, ": ", sum(ok), "/", nrow(s), " genes found")
    if (sum(ok) < 3) { res[[paste0(sn, "_score")]] <- NA_real_; next }
    m  <- as.matrix(expr[gg[ok], , drop = FALSE])
    mz <- t(scale(t(m)))                                   # z-score each gene across cells
    mz[is.na(mz)] <- 0                                     # zero-variance genes -> 0
    weighted <- all(!is.na(s$weight))
    if (weighted) {
      res[[paste0(sn, "_score")]] <- as.numeric(t(mz) %*% s$weight[ok])      # weighted sum (LSC17)
    } else {
      res[[paste0(sn, "_score")]] <- colMeans(mz)                            # unweighted = mean z-score
    }
  }
  res[, `:=`(Dataset = ds, Sample = samp)]
  fwrite_safe(res, out)
  message("      [write] ", out, "  (", nrow(res), " cells)")
  invisible()
}

samps <- list_qc_samples(opt$dataset)
if (!is.null(opt$sample)) samps <- samps[sample == opt$sample]
if (!nrow(samps)) stop("[d25] no samples for ", opt$dataset)
message("[2] ", nrow(samps), " sample(s)")
for (i in seq_len(nrow(samps))) score_one(samps$rds[i], samps$dataset[i], samps$sample[i])
message("[done] d25 for ", opt$dataset)