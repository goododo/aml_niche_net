#!/usr/bin/env Rscript
# =============================================================================
# 06_qc_figures.R   --   cohort-level QC figure set for the v2 re-run
#
# 03_per_sample_qc.R already writes one bar chart per dataset. That is enough to
# see whether a dataset ran; it is not enough to answer the questions this re-run
# was launched to answer. Each panel below is aimed at ONE item from the audit,
# so a reviewer can check the fix instead of taking it on faith:
#
#   Q1 cohort shape    did the sample SET change the way it was supposed to?
#                      (E-MTAB 12->8, Petti 9->5, GSE185381 53->52 donors,
#                       Chen2023 20 libs -> 10 donors, healthy arm 53->40)
#   Q2 filter strength per-sample filter strength varies 10-100x across this
#                      cohort and every score is compared BETWEEN samples, so
#                      "what counts as a cell" differing 10x propagates straight
#                      into edge weights. Shows loss split MAD vs ABS.
#   Q3 thresholds      the ABS floors were defined in config and wired to
#                      nothing until this re-run. Shows where the binding
#                      threshold actually came from, per sample.
#   Q4 doublets        the v1 headline "37% of expected doublets recovered"
#                      pooled author-prefiltered deposits with raw ones, which
#                      makes the number meaningless. Split by that flag.
#   Q5 vocabulary      the v2 timepoint vocabulary and the BM/PB split, as they
#                      actually landed -- not as the config says they should.
#
# Input : results/tables/01_preprocess/03_qc_report__ALL.csv
#         results/tables/01_preprocess/01_sample_role_manifest.csv (optional)
# Output: results/figures/01_preprocess/06_qc_*.{png,pdf}
#
# Usage: Rscript scripts/01_preprocess/06_qc_figures.R
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(data.table); library(ggplot2); library(scales)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

dir.create(FIG_PREPROCESS, showWarnings = FALSE, recursive = TRUE)

RPT <- file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv")
if (!file.exists(RPT)) stop("missing ", RPT, " -- run 03_per_sample_qc.R (all datasets) first")
R <- fread(RPT)
message("[0] ", nrow(R), " samples across ", uniqueN(R$dataset), " datasets")

MAN <- file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv")
M   <- if (file.exists(MAN)) fread(MAN) else NULL

# Consistent theme + saver. PDF alongside PNG because these end up in the deck.
th <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        strip.text        = element_text(size = 7.5),
        legend.position   = "top", legend.title = element_text(size = 8))

save2 <- function(p, name, w, h) {
  for (ext in c("png", "pdf"))
    ggsave(file.path(FIG_PREPROCESS, paste0(name, ".", ext)), p,
           width = w, height = h, dpi = 200, limitsize = FALSE)
  message("    wrote ", name)
}

# ---------------------------------------------------------------------------
# Q1  Cohort shape: raw -> post-MAD -> final, per dataset
# ---------------------------------------------------------------------------
message("[1] cohort shape")
ds <- R[, .(n_samples = .N,
            n_pass    = sum(status == "PASS", na.rm = TRUE),
            cells_raw = sum(n_raw, na.rm = TRUE),
            cells_fin = sum(n_final, na.rm = TRUE)), by = dataset]
setorder(ds, -cells_raw)

long <- melt(ds, id.vars = "dataset", measure.vars = c("cells_raw", "cells_fin"),
             variable.name = "stage", value.name = "cells")
long[, stage := factor(stage, levels = c("cells_raw", "cells_fin"),
                       labels = c("ingested", "after QC + doublets"))]
long[, dataset := factor(dataset, levels = ds$dataset)]

p1 <- ggplot(long, aes(dataset, cells, fill = stage)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(data = ds, aes(x = dataset, y = cells_raw, label = paste0(n_pass, "/", n_samples)),
            inherit.aes = FALSE, vjust = -0.4, size = 2.5, colour = "grey25") +
  scale_y_continuous(labels = label_comma(), expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c("ingested" = "#b8c9e0", "after QC + doublets" = "#3b7dd8")) +
  labs(title = "Q1  Cohort shape after the v2 sample-set fixes",
       subtitle = "label = samples passing / samples ingested",
       x = NULL, y = "cells", fill = NULL) +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1))
save2(p1, "06_qc_Q1_cohort_shape", 8, 4.2)

# ---------------------------------------------------------------------------
# Q2  Filter strength per sample, split by which rule did the cutting
# ---------------------------------------------------------------------------
message("[2] filter strength")
F <- R[!is.na(n_raw) & n_raw > 0]
F[, loss_mad := fifelse(is.na(n_fail_mad), NA_real_, n_fail_mad / n_raw)]
F[, loss_abs := fifelse(is.na(n_fail_abs), NA_real_, n_fail_abs / n_raw)]
F[, loss_tot := fifelse(is.na(n_after_mad), NA_real_, 1 - n_after_mad / n_raw)]

fl <- melt(F[, .(dataset, Sample, loss_mad, loss_abs)],
           id.vars = c("dataset", "Sample"), variable.name = "rule", value.name = "frac")
fl[, rule := factor(rule, levels = c("loss_mad", "loss_abs"),
                    labels = c("MAD outlier", "absolute floor/ceiling"))]
fl <- fl[!is.na(frac)]

if (nrow(fl)) {
  p2 <- ggplot(fl, aes(dataset, frac, colour = rule)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75),
                 width = 0.6, fill = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
               size = 0.7, alpha = 0.5) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_colour_manual(values = c("MAD outlier" = "#3b7dd8",
                                   "absolute floor/ceiling" = "#d1495b")) +
    labs(title = "Q2  Per-sample filter strength, split by which rule cut the cell",
         subtitle = paste("The absolute bounds were defined in config and wired to nothing before this",
                          "re-run.\nA sample where the red points dominate was being kept by MAD alone."),
         x = NULL, y = "fraction of ingested cells failing", colour = NULL) +
    th + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save2(p2, "06_qc_Q2_filter_strength", 8.5, 4.6)

  # Spread of total loss: the 10-100x claim, stated as a number per dataset.
  sp <- F[!is.na(loss_tot), .(min = min(loss_tot), med = median(loss_tot), max = max(loss_tot),
                              ratio = ifelse(min(loss_tot) > 0, max(loss_tot) / min(loss_tot), NA_real_)),
          by = dataset][order(-ratio)]
  fwrite_safe(sp, file.path(DIR_PREPROCESS, "06_qc_filter_strength_spread.csv"))
  message("    filter-strength spread written; worst within-dataset ratio: ",
          round(sp[1, ratio], 1), "x (", sp[1, dataset], ")")
}

# ---------------------------------------------------------------------------
# Q3  Which threshold actually bound, per sample
# ---------------------------------------------------------------------------
message("[3] binding thresholds")
# mad_keep() reports the EFFECTIVE threshold as max(MAD_lo, ABS) / min(MAD_hi, ABS).
# Comparing effective against the MAD-only value says which rule was in control.
T3 <- R[!is.na(nfeat_lo) & !is.na(nfeat_lo_mad),
        .(dataset, Sample, eff = nfeat_lo, mad = nfeat_lo_mad)]
if (nrow(T3)) {
  T3[, bound_by := fifelse(eff > mad + 1e-9, "absolute floor", "MAD")]
  p3 <- ggplot(T3, aes(mad, eff, colour = bound_by)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
    geom_hline(yintercept = ABS_MIN_NFEAT, linetype = 3, colour = "#d1495b") +
    geom_point(size = 1.1, alpha = 0.75) +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c("MAD" = "#3b7dd8", "absolute floor" = "#d1495b")) +
    labs(title = "Q3  nFeature floor: MAD-derived vs the threshold actually applied",
         subtitle = paste0("Points above the diagonal are samples the absolute floor (",
                           ABS_MIN_NFEAT, " genes, dotted) rescued from a\nMAD floor that ",
                           "log-scale had pushed implausibly low."),
         x = "MAD-derived floor (genes)", y = "effective floor (genes)", colour = "bound by") +
    th
  save2(p3, "06_qc_Q3_binding_threshold", 6, 4.6)
  message("    ", T3[bound_by == "absolute floor", .N], "/", nrow(T3),
          " samples bound by the absolute floor")
}

# ---------------------------------------------------------------------------
# Q4  Doublet recovery, split by whether the deposit was already filtered
# ---------------------------------------------------------------------------
message("[4] doublet recovery")
D <- R[!is.na(dbl_rate_exp) & dbl_rate_exp > 0 & !is.na(dbl_rate_obs)]
if (nrow(D)) {
  D[, upstream := fifelse(as.logical(upstream_filtered) %in% TRUE,
                          "author-prefiltered deposit", "raw matrix")]
  D[, recovery := dbl_rate_obs / dbl_rate_exp]
  p4 <- ggplot(D, aes(upstream, recovery, colour = upstream)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
    geom_boxplot(outlier.shape = NA, width = 0.45, fill = NA) +
    geom_point(position = position_jitter(width = 0.13), size = 1.1, alpha = 0.6) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_colour_manual(values = c("raw matrix" = "#3b7dd8",
                                   "author-prefiltered deposit" = "#8f8f8f"), guide = "none") +
    labs(title = "Q4  Doublet recovery: observed rate / expected rate",
         subtitle = paste("A low value on the grey side is EXPECTED -- those authors already removed",
                          "doublets.\nPooling the two groups is what made the v1 headline 37% uninterpretable."),
         x = NULL, y = "observed / expected") +
    th
  save2(p4, "06_qc_Q4_doublet_recovery", 5.5, 4.4)
  message("    median recovery  raw: ",
          round(D[upstream == "raw matrix", median(recovery, na.rm = TRUE)], 3),
          "  |  prefiltered: ",
          round(D[upstream != "raw matrix", median(recovery, na.rm = TRUE)], 3))
}

# ---------------------------------------------------------------------------
# Q5  Cohort composition under the v2 vocabulary
# ---------------------------------------------------------------------------
if (!is.null(M) && "Timepoint" %in% names(M)) {
  message("[5] cohort composition")
  # INNER join, and PASS only. The manifest is built at INGEST level, so it still contains
  # GSE185381's 27 UNDEMUX__* rows -- cells the author metadata could not assign to a donor,
  # removed by the demux prefilter before the per-sample split. An outer join with
  # `status == "PASS" | is.na(status)` admits exactly those (they have no QC row, so status is
  # NA), which inflated the healthy arm from 40 to 47: ten un-demuxed all-control LIBRARIES
  # counted as healthy donors. A sample that produced no QC row is not a sample.
  Mp <- merge(M, R[, .(dataset, Sample, status)], by = c("dataset", "Sample"))
  Mp <- Mp[status == "PASS"]
  Mp[, Timepoint := factor(Timepoint, levels = intersect(CANONICAL_TIMEPOINTS, unique(Timepoint)))]

  p5 <- ggplot(Mp[!is.na(Timepoint)], aes(dataset, fill = Timepoint)) +
    geom_bar(width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Q5  Cohort composition under the v2 timepoint vocabulary",
         subtitle = paste("Post_treatment and MRD no longer exist. Refractory and On_treatment are new;",
                          "without them\nGSE207356's screening draw had nowhere to go but Diagnosis."),
         x = NULL, y = "samples passing QC", fill = NULL) +
    th + theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    guides(fill = guide_legend(nrow = 2))
  save2(p5, "06_qc_Q5_composition", 8, 4.8)

  cmp <- dcast(Mp[!is.na(Timepoint), .N, by = .(dataset, Timepoint)],
               dataset ~ Timepoint, value.var = "N", fill = 0)
  fwrite_safe(cmp, file.path(DIR_PREPROCESS, "06_qc_composition.csv"))
  print(cmp)

  if ("Tissue" %in% names(Mp)) {
    message("\n    tissue: ", paste(sprintf("%s=%d", names(table(Mp$Tissue)), table(Mp$Tissue)),
                                    collapse = "  "))
    if ("PB" %in% Mp$Tissue) print(Mp[Tissue == "PB", .(dataset, Sample, Timepoint)])
  }
  if ("patient_resolved" %in% names(Mp)) {
    ur <- Mp[patient_resolved == FALSE]
    message("    patient_resolved=FALSE: ", nrow(ur),
            if (nrow(ur)) paste0(" (", paste(ur$Sample, collapse = ", "), ")") else "")
  }
  hh <- Mp[Timepoint == "Healthy", .N, by = dataset][order(-N)]
  message("\n    HEALTHY ARM: ", sum(hh$N), " samples")
  print(hh)
} else message("[5] skipped: no manifest (run 01_dataset_roles.R first)")

message("\n[done] figures in ", FIG_PREPROCESS)
