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
# REFNORM_REF_CELL_DIR lives here: the reference-cells-out-of-the-denominator arm reads the
# per-sample ref_norm cell lists, so this file is a hard dependency of that arm, not optional.
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
HIER_PROJ_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")

## ---- which samples are healthy (label first, name pattern as a safety net) ----
# Roster from the QC report, not from ls: this script IS the calibration gate, so
# a cohort silently widened by leftover objects from a previous ingest would move
# the very number the gate is read off. See qc_rds_roster() in utils.R.
R <- qc_rds_roster(on_extra = "error")
info <- rbindlist(lapply(seq_len(nrow(R)), function(i) {
  m <- tryCatch(readRDS(R$rds[i])@meta.data, error = function(e) NULL); if (is.null(m)) return(NULL)
  data.table(dataset = R$dataset[i], sample = R$sample[i],
             timepoint = if ("Timepoint" %in% names(m)) as.character(m$Timepoint[1]) else NA_character_,
             n_cells_qc = nrow(m))
}), fill = TRUE)
info[, healthy := timepoint == "Healthy" | mapply(function(s) isTRUE(is_healthy_sample(s)), sample)]
message(sprintf("[0] %d samples, %d healthy donors", nrow(info), info[healthy == TRUE, .N]))

## ---- the autologous reference cells are forced malignant=0 by construction ----
# 41_infercnv_to_percell sets `malignant := 0L` for every reference cell, so those
# cells cannot be false positives no matter how they score. That is defensible per
# sample, but it makes the FPR depend on how BIG the reference is -- and the
# reference fraction is not comparable across arms: healthy marrow is rich in the
# mature T/B cells the reference is drawn from (median frac_ref 0.31) while
# diagnosis marrow is blast-packed (0.13). The bias runs toward a LOWER healthy
# score, i.e. toward the gate passing. So report the FPR with those cells removed
# from the denominator as well, and read the gate against both.
ref_cells_of <- function(ds, sid) {
  f <- file.path(REFNORM_REF_CELL_DIR, ds, paste0(sid, "_ref_norm_cells.txt"))
  if (file.exists(f)) readLines(f) else character(0)
}

## ---- per-sample FPR = fraction of cells called malignant in a healthy marrow ----
res <- rbindlist(lapply(which(info$healthy), function(i) {
  ds <- info$dataset[i]; sid <- info$sample[i]
  cf <- file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
  if (!file.exists(cf)) return(NULL)
  c <- fread(cf, select = c("cell", "malignant"))
  pf <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  bin <- if (file.exists(pf)) fread(pf, select = c("cell", "hierarchy_bin", "in_ccc_graph")) else NULL
  if (!is.null(bin)) c <- merge(c, bin, by = "cell", all.x = TRUE)
  rc <- ref_cells_of(ds, sid)
  c[, is_ref := cell %in% rc]
  nr <- c[is_ref == FALSE]
  data.table(dataset = ds, sample = sid,
             n_evaluable = sum(!is.na(c$malignant)),
             n_false_pos = sum(c$malignant == 1, na.rm = TRUE),
             fpr = round(mean(c$malignant == 1, na.rm = TRUE), 4),
             fpr_ccc_bins = if (!is.null(bin))
               round(mean(c[in_ccc_graph == TRUE]$malignant == 1, na.rm = TRUE), 4) else NA_real_,
             n_ref_cells = length(rc),
             frac_ref = round(mean(c$is_ref), 4),
             # same numerator, denominator excludes the cells that were defined normal
             n_evaluable_excl_ref = nrow(nr),
             fpr_excl_ref = if (nrow(nr)) round(mean(nr$malignant == 1, na.rm = TRUE), 4) else NA_real_)
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
if ("fpr_excl_ref" %in% names(res)) {
  cat(sprintf("  OVERALL FPR excluding the autologous reference cells = %.4f (%d of %d)\n",
              sum(res$n_false_pos) / sum(res$n_evaluable_excl_ref),
              sum(res$n_false_pos), sum(res$n_evaluable_excl_ref)))
  cat(sprintf("  reference cells are %.1f%% of healthy cells and are defined normal, so the first\n",
              100 * (1 - sum(res$n_evaluable_excl_ref) / sum(res$n_evaluable))))
  cat("  number is the one the pipeline acts on and the second is the one that is comparable\n")
  cat("  across arms. Report both; the gap IS the reference-fraction confound.\n")
}

## ---- per-bin FPR: does the caller misfire preferentially in some compartments? ----
bins <- rbindlist(lapply(which(info$healthy), function(i) {
  ds <- info$dataset[i]; sid <- info$sample[i]
  cf <- file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
  pf <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  if (!file.exists(cf) || !file.exists(pf)) return(NULL)
  hdr <- names(fread(cf, nrows = 0L))
  cc  <- if ("is_ref" %in% hdr) c("cell", "malignant", "is_ref") else c("cell", "malignant")
  cm  <- fread(cf, select = cc)
  if (!"is_ref" %in% names(cm)) cm[, is_ref := NA]
  m <- merge(cm, fread(pf, select = c("cell", "hierarchy_bin", "in_ccc_graph")), by = "cell")
  m[, .(dataset = ds, sample = sid, hierarchy_bin, in_ccc_graph, malignant, is_ref = as.logical(is_ref))]
}), fill = TRUE)

if (nrow(bins)) {
  # BOTH DENOMINATORS, ALWAYS. 41_infercnv_to_percell.R ASSIGNS malignant=0 to every autologous
  # reference cell -- those cells cannot be false positives whatever their burden, so a rate over
  # all cells is arithmetically FPR_evaluable x evaluable_share, not a false-positive rate.
  # Measured effect when this was a single number: T_NK published 0.0557 against a true evaluable
  # rate of 0.2598 (78.6% of that bin is forced zeros), and the BIN RANKING INVERTED -- T_NK was
  # plotted as the cleanest compartment when it is among the dirtiest. The per-sample arm of this
  # same script already reported both via fpr_excl_ref; only the per-bin arm did not.
  if (bins[!is.na(is_ref), .N] == 0)
    warning("no is_ref column in the consensus files -- rerun 41 and 50; reporting the padded rate only")
  bt <- bins[, .(n = .N,
                 false_pos = sum(malignant == 1, na.rm = TRUE),
                 n_excl_ref = sum(!(is_ref %in% TRUE)),
                 false_pos_excl_ref = sum(malignant == 1 & !(is_ref %in% TRUE), na.rm = TRUE)),
             by = .(hierarchy_bin, in_ccc_graph)]
  bt[, FPR := round(false_pos / n, 4)]
  bt[, FPR_excl_ref := round(fifelse(n_excl_ref > 0, false_pos_excl_ref / n_excl_ref, NA_real_), 4)]
  bt[, ref_share := round(1 - n_excl_ref / n, 4)]
  setorder(bt, -FPR_excl_ref)
  fwrite(bt, file.path(DIR_MALIGNANCY, "malignancy_fpr_by_bin.csv"))
  cat("\n================ false-positive rate by hierarchy bin (healthy donors) ================\n")
  # ordered by the COMPARABLE rate, not the padded one -- ordering by FPR is what put T_NK at the
  # bottom of this table (0.0557) and had it read as the cleanest compartment
  print(bt[order(-FPR_excl_ref)])
  cat("\n  FPR is over ALL cells (what the pipeline acts on); FPR_excl_ref drops the autologous\n")
  cat("  reference cells, which are DEFINED normal and cannot be positive. ref_share is how much of\n")
  cat("  the bin that is. Read the second column when comparing bins -- the first one ranks bins by\n")
  cat("  how large their reference is as much as by how often the caller misfires.\n")
  cat("\n[reading it] A bin with a high FPR is one where the CNV caller cannot be trusted. Compare\n")
  cat("  these numbers against the AML per-bin malignant fractions from 03_hierarchy: a bin whose\n")
  cat("  AML fraction is not clearly above its healthy FPR carries no usable malignancy signal.\n")

  ## ---- bin x dataset: is a bin's FPR a property of the caller, or of one dataset? ----
  # The pooled per-bin FPR above is dominated by whichever dataset contributes the most cells
  # to that bin, so a bin can look unusable because ONE cohort misfires there. Splitting by
  # dataset is what separates "the caller cannot read this compartment" from "this dataset is
  # the problem" -- and only the first justifies dropping the bin cohort-wide.
  bd <- bins[, .(n = .N, false_pos = sum(malignant == 1, na.rm = TRUE)),
             by = .(hierarchy_bin, in_ccc_graph, dataset)]
  bd[, FPR := round(false_pos / n, 4)]
  bd[, n_healthy := bins[, uniqueN(sample), by = dataset][match(bd$dataset, dataset)]$V1]
  setorder(bd, hierarchy_bin, -FPR)
  fwrite(bd, file.path(DIR_MALIGNANCY, "malignancy_fpr_by_bin_dataset.csv"))
  cat("\n=========== false-positive rate by hierarchy bin x dataset (healthy donors) ===========\n")
  print(dcast(bd[in_ccc_graph == TRUE], hierarchy_bin ~ dataset, value.var = "FPR"))
  cat("\n  cell counts behind each rate (an FPR over few hundred cells is not a stable estimate):\n")
  print(dcast(bd[in_ccc_graph == TRUE], hierarchy_bin ~ dataset, value.var = "n", fill = 0L))
  cat("\n[reading it] A bin whose FPR is high in EVERY dataset is a caller limitation. A bin that is\n")
  cat("  high in one dataset only is a cohort problem, and dropping the bin everywhere would throw\n")
  cat("  away usable signal from the other datasets.\n")
}

cat("\n[caveat] These healthy samples span several platforms/depths, so the FPR is a cohort-level\n")
cat("  estimate rather than a per-sample correction factor. It is a floor for how much of the\n")
cat("  called malignancy is noise -- not a value to subtract cell-by-cell.\n")
