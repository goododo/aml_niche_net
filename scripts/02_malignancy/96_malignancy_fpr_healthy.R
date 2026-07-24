#!/usr/bin/env Rscript
# 96_malignancy_fpr_healthy.R ----
# Use the healthy-donor samples as a BUILT-IN NEGATIVE CONTROL for the malignancy calls. A normal
# marrow contains no malignant cells, so every cell the consensus flags there is a false positive.
# This turns "do we trust the malignancy labels?" from an argument into a measurement, and gives a
# per-dataset / per-bin false-positive rate (FPR) to report, calibrate against, or subtract.
#
# Motivation: healthy samples in this cohort carry called malignant fractions from 0.08 up to 0.76.
# Everything built on the labels inherits that error -- the CCC frac_malignant node feature, the
# per-bin malignant fractions in 03_hierarchy, and the malignant arm of 04_cnmf.
#
# INPUT : QC objects (Timepoint) ; consensus per-cell ; BMM projection per-cell (for the per-bin FPR)
# OUTPUT: DIR_MALIGNANCY/malignancy_fpr_healthy.csv        (per healthy sample)
#         DIR_MALIGNANCY/malignancy_fpr_by_bin.csv         (per hierarchy bin, pooled)
#
#   Rscript scripts/02_malignancy/96_malignancy_fpr_healthy.R

suppressPackageStartupMessages({ library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
HIER_PROJ_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")

## ---- which samples are healthy (label first, name pattern as a safety net) ----
rds <- list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
rds <- rds[!grepl("/_", rds)]
info <- rbindlist(lapply(rds, function(f) {
  m <- tryCatch(readRDS(f)@meta.data, error = function(e) NULL); if (is.null(m)) return(NULL)
  s <- sub("\\.rds$", "", basename(f))
  data.table(dataset = basename(dirname(f)), sample = s,
             timepoint = if ("Timepoint" %in% names(m)) as.character(m$Timepoint[1]) else NA_character_,
             n_cells_qc = nrow(m))
}), fill = TRUE)
info[, healthy := timepoint == "Healthy" | mapply(function(s) isTRUE(is_healthy_sample(s)), sample)]
message(sprintf("[0] %d samples, %d healthy donors", nrow(info), info[healthy == TRUE, .N]))

## ---- per-sample FPR = fraction of cells called malignant in a healthy marrow ----
res <- rbindlist(lapply(which(info$healthy), function(i) {
  ds <- info$dataset[i]; sid <- info$sample[i]
  cf <- file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
  if (!file.exists(cf)) return(NULL)
  c <- fread(cf, select = c("cell", "malignant"))
  pf <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  bin <- if (file.exists(pf)) fread(pf, select = c("cell", "hierarchy_bin", "in_ccc_graph")) else NULL
  if (!is.null(bin)) c <- merge(c, bin, by = "cell", all.x = TRUE)
  data.table(dataset = ds, sample = sid,
             n_evaluable = sum(!is.na(c$malignant)),
             n_false_pos = sum(c$malignant == 1, na.rm = TRUE),
             fpr = round(mean(c$malignant == 1, na.rm = TRUE), 4),
             fpr_ccc_bins = if (!is.null(bin))
               round(mean(c[in_ccc_graph == TRUE]$malignant == 1, na.rm = TRUE), 4) else NA_real_)
}), fill = TRUE)
setorder(res, -fpr)
fwrite(res, file.path(DIR_MALIGNANCY, "malignancy_fpr_healthy.csv"))

cat("\n================ per-sample false-positive rate (healthy donors) ================\n")
print(res)
cat("\n-- pooled by dataset --\n")
print(res[, .(n_healthy = .N, cells = sum(n_evaluable), false_pos = sum(n_false_pos),
              FPR = round(sum(n_false_pos) / sum(n_evaluable), 4)), by = dataset][order(-FPR)])
cat(sprintf("\n  OVERALL FPR = %.4f  (%d of %d cells in healthy marrow called malignant)\n",
            sum(res$n_false_pos) / sum(res$n_evaluable), sum(res$n_false_pos), sum(res$n_evaluable)))

## ---- per-bin FPR: does the caller misfire preferentially in some compartments? ----
bins <- rbindlist(lapply(which(info$healthy), function(i) {
  ds <- info$dataset[i]; sid <- info$sample[i]
  cf <- file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
  pf <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  if (!file.exists(cf) || !file.exists(pf)) return(NULL)
  m <- merge(fread(cf, select = c("cell", "malignant")),
             fread(pf, select = c("cell", "hierarchy_bin", "in_ccc_graph")), by = "cell")
  m[, .(dataset = ds, sample = sid, hierarchy_bin, in_ccc_graph, malignant)]
}), fill = TRUE)

if (nrow(bins)) {
  bt <- bins[, .(n = .N, false_pos = sum(malignant == 1, na.rm = TRUE)), by = .(hierarchy_bin, in_ccc_graph)]
  bt[, FPR := round(false_pos / n, 4)][order(-FPR)]
  fwrite(bt, file.path(DIR_MALIGNANCY, "malignancy_fpr_by_bin.csv"))
  cat("\n================ false-positive rate by hierarchy bin (healthy donors) ================\n")
  print(bt[order(-FPR)])
  cat("\n[reading it] A bin with a high FPR is one where the CNV caller cannot be trusted. Compare\n")
  cat("  these numbers against the AML per-bin malignant fractions from 03_hierarchy: a bin whose\n")
  cat("  AML fraction is not clearly above its healthy FPR carries no usable malignancy signal.\n")
}

cat("\n[caveat] These healthy samples span several platforms/depths, so the FPR is a cohort-level\n")
cat("  estimate rather than a per-sample correction factor. It is a floor for how much of the\n")
cat("  called malignancy is noise -- not a value to subtract cell-by-cell.\n")
