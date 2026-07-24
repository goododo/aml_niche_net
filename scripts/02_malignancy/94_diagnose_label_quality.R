#!/usr/bin/env Rscript
# 94_diagnose_label_quality.R ----
# Two diagnostics that decide whether upstream QC is confounding the malignancy labels and the
# meta-programs built on top of them. READ-ONLY.
#
#  (1) QUALITY-DRIVEN MALIGNANCY: is per-sample malignant_frac correlated with per-sample data
#      quality (median nFeature / nCount)? If yes, "how much tumour a sample has" is partly "how
#      deeply it was sequenced", which confounds 02/03/04 AND the CCC frac_malignant node feature.
#      Motivation: GSE239721 median malignant_frac = 0.011 at DIAGNOSIS (implausible; real AML
#      aspirates are 20-80% blasts) and GSE201966 = 0.058 with ~182 genes/cell.
#
#  (2) MP8 DOUBLET TEST: are MP8-dominant cells (Plasma/B Ig program, diffuse across bins) enriched
#      for doublet-like features (higher nFeature_RNA, higher Doublet_Score)? Motivation: the
#      doublet consensus is "intersection", recovering only ~37% of expected doublets; a
#      blast+plasma doublet carries Ig into any bin, which is exactly MP8's diffuse signature.
#
#   Rscript scripts/02_malignancy/94_diagnose_label_quality.R [--datasets GSE227903,Chen2023]

suppressPackageStartupMessages({ library(Seurat); library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", default = "",
              help = "dataset filter for the (slower) per-cell test [default: all]"),
  make_option("--max_samples", type = "integer", default = 40L,
              help = "cap QC objects read for the per-cell test [%default]")
)))

## ================= (1) malignant_frac vs data quality =================
qc <- data.table(rds = list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE))
qc <- qc[!grepl("/_", rds)][, sample := sub("\\.rds$", "", basename(rds))][, dataset := basename(dirname(rds))]
summ <- fread(file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv"))
man <- merge(summ[, .(dataset, sample, malignant_frac, n_qc_cells, evidence_tier)], qc, by = c("dataset", "sample"))
message(sprintf("[1] reading QC metadata for %d samples ...", nrow(man)))

qstat <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  m <- tryCatch(readRDS(man$rds[i])@meta.data, error = function(e) NULL)
  if (is.null(m) || !"nFeature_RNA" %in% names(m)) return(NULL)
  data.table(dataset = man$dataset[i], sample = man$sample[i],
             med_nFeature = median(m$nFeature_RNA, na.rm = TRUE),
             med_nCount   = median(m$nCount_RNA,  na.rm = TRUE))
}), fill = TRUE)
d <- merge(man[, .(dataset, sample, malignant_frac, n_qc_cells, evidence_tier)], qstat, by = c("dataset", "sample"))

cat("\n================ (1) malignant_frac vs data quality ================\n")
cat(sprintf("n = %d samples\n", nrow(d)))
cat(sprintf("  Spearman(malignant_frac, med_nFeature) = %.3f  (p = %.3g)\n",
            cor(d$malignant_frac, d$med_nFeature, method = "spearman", use = "complete.obs"),
            cor.test(d$malignant_frac, d$med_nFeature, method = "spearman")$p.value))
cat(sprintf("  Spearman(malignant_frac, med_nCount)   = %.3f\n",
            cor(d$malignant_frac, d$med_nCount, method = "spearman", use = "complete.obs")))
cat("\n-- per dataset: quality vs called malignancy --\n")
print(d[, .(n = .N,
            med_nFeature = as.integer(median(med_nFeature)),
            med_nCount   = as.integer(median(med_nCount)),
            med_mal_frac = round(median(malignant_frac), 3)), by = dataset][order(med_nFeature)])
cat("\n-- samples with implausibly low malignancy despite decent depth (candidate under-calls) --\n")
print(d[malignant_frac < 0.05 & med_nFeature > 800][order(malignant_frac),
        .(dataset, sample, n_qc_cells, med_nFeature = as.integer(med_nFeature),
          malignant_frac, evidence_tier)][1:min(20, .N)])
cat("\n-- lowest-complexity samples (CCC-unusable: cannot co-detect ligand+receptor) --\n")
print(d[order(med_nFeature)][1:min(15, .N),
        .(dataset, sample, n_qc_cells, med_nFeature = as.integer(med_nFeature), malignant_frac)])
fwrite(d, file.path(DIR_MALIGNANCY, "sample_quality_vs_malignancy.csv"))
cat("\n[interpretation] a strong positive Spearman means malignancy calling tracks sequencing depth,\n")
cat("  i.e. frac_malignant is partly a QC readout -- a confounder for 02/03/04 and for the CCC\n")
cat("  frac_malignant node feature. Near-zero means the labels are quality-robust.\n")

## ================= (2) MP8-dominant cells: doublet-like? =================
usage_dir <- file.path(LARGE1_DIR, "04_cnmf", "malignant", "mp_usage")
uf <- list.files(usage_dir, pattern = "__mp_usage\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(uf)) {
  u <- rbindlist(lapply(uf, fread), fill = TRUE)[, .(cell, sample, dataset, hierarchy_bin, top_MP)]
  sel <- unique(u[, .(dataset, sample)])
  if (nzchar(opt$datasets)) sel <- sel[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
  sel <- sel[seq_len(min(opt$max_samples, .N))]
  message(sprintf("\n[2] per-cell QC metadata for %d samples ...", nrow(sel)))
  meta <- rbindlist(lapply(seq_len(nrow(sel)), function(i) {
    p <- man[dataset == sel$dataset[i] & sample == sel$sample[i]]$rds
    if (!length(p)) return(NULL)
    m <- tryCatch(readRDS(p[1])@meta.data, error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.table(cell = rownames(m),
               nFeature_RNA = m$nFeature_RNA,
               Doublet_Score = if ("Doublet_Score" %in% names(m)) m$Doublet_Score else NA_real_)
  }), fill = TRUE)
  j <- merge(u, meta, by = "cell")
  cat("\n================ (2) MP8 doublet test ================\n")
  cat(sprintf("n = %d cells with QC metadata\n", nrow(j)))
  st <- j[, .(n = .N,
              med_nFeature = as.integer(median(nFeature_RNA, na.rm = TRUE)),
              med_doublet  = round(median(Doublet_Score, na.rm = TRUE), 3)), by = top_MP][order(-med_nFeature)]
  print(st)
  if (nrow(j[top_MP == "MP8"]) > 10) {
    w <- wilcox.test(j[top_MP == "MP8"]$nFeature_RNA, j[top_MP != "MP8"]$nFeature_RNA)
    cat(sprintf("\n  MP8 vs rest, nFeature_RNA: median %d vs %d (Wilcoxon p = %.3g)\n",
                as.integer(median(j[top_MP == "MP8"]$nFeature_RNA, na.rm = TRUE)),
                as.integer(median(j[top_MP != "MP8"]$nFeature_RNA, na.rm = TRUE)), w$p.value))
  }
  cat("\n[interpretation] MP8 sitting at the TOP of med_nFeature (and/or med_doublet) supports a\n")
  cat("  residual-doublet origin -> keep the likely_contaminant flag and exclude MP8 downstream.\n")
  cat("  If MP8 is unremarkable, the Ig signal is more likely ambient RNA or true plasma cells\n")
  cat("  miscalled as malignant -- a different fix (malignancy calling, not doublet removal).\n")
} else message("[2] no mp_usage tables found; skipping")
