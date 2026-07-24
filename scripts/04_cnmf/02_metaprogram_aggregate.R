#!/usr/bin/env Rscript
# 02_metaprogram_aggregate.R ----
# Phase L1.2 -- cluster the per-sample NMF programs (from 01) into robust cross-sample META-PROGRAMS
# (GeneNMF::getMetaPrograms), per compartment. Then, when BOTH compartments are aggregated, compare
# malignant MPs vs normal MPs by gene-overlap and flag each malignant MP as tumor-specific or
# normal-like (co-opted / lineage-contaminant) -- the tumor-specific set is what feeds downstream.
#
# INPUT : CNMF_RES_DIR/<compartment>/nmf_res.rds        (from 01)
# OUTPUT: CNMF_TAB_DIR/<compartment>_metaprograms.tsv   (MP, rank, gene, weight)
#         CNMF_TAB_DIR/<compartment>_mp_metrics.csv     (sampleCoverage, silhouette, ... per MP)
#         CNMF_FIG_DIR/<compartment>_metaprograms.png   (program tree / heatmap)
#   and (once both compartments exist):
#         CNMF_TAB_DIR/malignant_vs_normal_similarity.csv   (Jaccard matrix, long)
#         CNMF_TAB_DIR/malignant_mp_annotation.csv          (per malignant MP: tumor_specific | normal_like)
#
# Run per compartment; the malignant run auto-does the comparison if normal MPs already exist.
#   Rscript scripts/04_cnmf/02_metaprogram_aggregate.R --compartment malignant [--nMP 12]

suppressPackageStartupMessages({ library(GeneNMF); library(data.table); library(optparse); library(ggplot2) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--compartment", type = "character", default = "malignant", help = "malignant | normal"),
  make_option("--nMP",         type = "integer",   default = CNMF_NMP, help = "target number of meta-programs [%default]")
)))
stopifnot(opt$compartment %in% c("malignant", "normal"))

res_rds <- file.path(CNMF_RES_DIR, opt$compartment, "nmf_res.rds")
stopifnot(file.exists(res_rds))
nmf <- readRDS(res_rds)
message(sprintf("[0] %s: %d per-sample program sets", opt$compartment, length(nmf)))

## ---- aggregate into meta-programs ----
mp <- getMetaPrograms(nmf, nMP = opt$nMP, metric = CNMF_METRIC, specificity.weight = CNMF_SPEC_WT,
                      max.genes = CNMF_MAX_GENES, min.confidence = CNMF_MIN_CONF)
genes <- mp$metaprograms.genes
wts   <- mp$metaprograms.genes.weights
met   <- as.data.table(mp$metaprograms.metrics, keep.rownames = "MP")

## signatures (long: MP, rank, gene, weight)
sig <- rbindlist(lapply(names(genes), function(m) {
  g <- genes[[m]]; w <- wts[[m]][g]
  data.table(MP = m, rank = seq_along(g), gene = g, weight = round(as.numeric(w), 5))
}), fill = TRUE)
fwrite_safe(sig, file.path(CNMF_TAB_DIR, paste0(opt$compartment, "_metaprograms.tsv")))
fwrite_safe(met, file.path(CNMF_TAB_DIR, paste0(opt$compartment, "_mp_metrics.csv")))

cat(sprintf("\n================ %s meta-programs (n=%d) ================\n", opt$compartment, length(genes)))
met_print <- copy(met); met_print[, top_genes := sapply(genes[MP], function(g) paste(head(g, 8), collapse = ","))]
print(met_print[, .(MP, sampleCoverage, silhouette, numberGenes, numberPrograms, top_genes)])
cat("\n[guide] trust MPs with high sampleCoverage (>~0.3) AND silhouette (>0): those are robust,\n")
cat("        recurrent across samples. Low-coverage MPs are patient-private or noise.\n")

## tree/heatmap plot (best-effort; GeneNMF's plotMetaPrograms)
png_out <- file.path(CNMF_FIG_DIR, paste0(opt$compartment, "_metaprograms.png"))
tryCatch({ png(png_out, width = 1400, height = 1000, res = 130); print(plotMetaPrograms(mp)); dev.off()
           message("[fig] ", png_out) }, error = function(e) message("[fig] plotMetaPrograms skipped: ", conditionMessage(e)))

saveRDS(mp, file.path(CNMF_RES_DIR, opt$compartment, "metaprograms.rds"))

## ---- malignant vs normal comparison (only meaningful from the malignant run, once normal exists) ----
norm_sig_f <- file.path(CNMF_TAB_DIR, "normal_metaprograms.tsv")
mal_sig_f  <- file.path(CNMF_TAB_DIR, "malignant_metaprograms.tsv")
if (opt$compartment == "malignant" && file.exists(norm_sig_f)) {
  jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  mal_g  <- split(sig$gene, sig$MP)
  norm_sig <- fread(norm_sig_f); norm_g <- split(norm_sig$gene, norm_sig$MP)
  simdt <- rbindlist(lapply(names(mal_g), function(m)
    rbindlist(lapply(names(norm_g), function(n)
      data.table(malignant_MP = m, normal_MP = n, jaccard = round(jacc(mal_g[[m]], norm_g[[n]]), 4))))))
  fwrite_safe(simdt, file.path(CNMF_TAB_DIR, "malignant_vs_normal_similarity.csv"))

  ann <- simdt[, .(max_jaccard = max(jaccard), nearest_normal_MP = normal_MP[which.max(jaccard)]), by = malignant_MP]
  ann[met, `:=`(sampleCoverage = i.sampleCoverage, silhouette = i.silhouette), on = c(malignant_MP = "MP")]
  ann[, class := fifelse(max_jaccard > CNMF_MP_SIM_MAX, "normal_like", "tumor_specific")]
  setorder(ann, -max_jaccard)
  fwrite_safe(ann, file.path(CNMF_TAB_DIR, "malignant_mp_annotation.csv"))

  cat("\n================ malignant MP tumor-specificity (vs normal MPs) ================\n")
  cat(sprintf("   normal_like = shares > %.2f Jaccard with some normal MP (co-opted / lineage contaminant)\n", CNMF_MP_SIM_MAX))
  print(ann)
  cat(sprintf("\n   -> %d tumor-specific, %d normal-like of %d malignant MPs\n",
              ann[class == "tumor_specific", .N], ann[class == "normal_like", .N], nrow(ann)))
} else if (opt$compartment == "malignant") {
  cat("\n[comparison] normal_metaprograms.tsv not found yet -- run the normal compartment, then re-run\n")
  cat("             this (malignant) to get the tumor-specific vs normal-like annotation.\n")
}
