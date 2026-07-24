# =====================================================================
# f01_qc_overview.R  --  cohort-level QC overview from 03_qc_report__ALL.csv
# ---------------------------------------------------------------------
# No RDS loading. Panels: retention waterfall / doublet rate /
# sample pass-exclude / median depth / low-complexity flags.
# Column names are auto-resolved (case-insensitive) and printed for review.
# =====================================================================
.here <- tryCatch(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                  error = function(e) ".")
if (length(.here) == 0 || .here == "") .here <- "."
source(file.path(.here, "f00_fig_config.R"))

stopifnot(file.exists(QC_REPORT_CSV))
qc <- fread(QC_REPORT_CSV)
message("[f01] columns in QC report:")
print(names(qc))

## -- column resolver ----
pick <- function(cands, required = TRUE) {
  hit <- names(qc)[tolower(names(qc)) %in% tolower(cands)]
  if (!length(hit)) {
    if (required) stop(sprintf("missing column; tried: %s", paste(cands, collapse = ", ")))
    return(NA_character_)
  }
  hit[1]
}
C_ds   <- pick(c("Dataset", "dataset"))
C_samp <- pick(c("Sample", "sample"), required = FALSE)
C_raw  <- pick(c("n_raw", "nRaw"))
C_mad  <- pick(c("n_after_mad", "n_postmad", "n_after_MAD"))
C_dbl  <- pick(c("n_doublet", "n_doublets"))
C_fin  <- pick(c("n_final", "n_kept"))
C_stat <- pick(c("status"), required = FALSE)
C_reas <- pick(c("reason", "exclude_reason"), required = FALSE)
C_med  <- pick(c("med_ncount_final", "median_ncount_final", "med_ncount"), required = FALSE)
C_lowc <- pick(c("lowcomplexity_flag", "lowcomplexity", "low_complexity_flag"), required = FALSE)

setnames(qc, C_ds, "Dataset")
qc[, Dataset := factor(Dataset, levels = sort(unique(Dataset)))]
is_pass <- if (!is.na(C_stat)) toupper(qc[[C_stat]]) %chin% c("PASS", "OK", "KEEP", "INCLUDE") else rep(TRUE, nrow(qc))

panels <- list()

## -- A. retention waterfall (PASS samples; cells removed by stage) ----
wf <- qc[is_pass, .(
  mad_removed     = sum(pmax(get(C_raw) - get(C_mad), 0), na.rm = TRUE),
  doublet_removed = sum(get(C_dbl), na.rm = TRUE),
  kept            = sum(get(C_fin), na.rm = TRUE)
), by = Dataset]
wf[, c("mad_removed","doublet_removed","kept") := lapply(.SD, as.numeric),
   .SDcols = c("mad_removed","doublet_removed","kept")]

wfl <- melt(wf, id.vars = "Dataset", variable.name = "fate", value.name = "cells")
wfl[, fate := factor(fate, levels = c("kept", "doublet_removed", "mad_removed"),
                     labels = c("Kept (final)", "Removed: doublet", "Removed: MAD outlier"))]
panels$A <- ggplot(wfl, aes(Dataset, cells, fill = fate)) +
  geom_col() +
  scale_fill_manual(values = c("Kept (final)" = "#2c7fb8",
                               "Removed: doublet" = "#f4a259",
                               "Removed: MAD outlier" = "#cbd5dd")) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "A. Cell retention by QC stage", x = NULL, y = "Cells", fill = NULL) +
  theme_qc()

## -- B. doublet rate per dataset ----
qc[, .dbl_rate := get(C_dbl) / pmax(get(C_mad), 1)]
panels$B <- ggplot(qc[is_pass], aes(Dataset, .dbl_rate)) +
  geom_boxplot(outlier.size = 0.6, fill = "#e9eef2") +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "B. Doublet rate (consensus) per sample", x = NULL, y = "doublet / post-MAD") +
  theme_qc()

## -- C. sample pass / exclude ----
if (!is.na(C_stat)) {
  st <- qc[, .N, by = .(Dataset, status = get(C_stat))]
  panels$C <- ggplot(st, aes(Dataset, N, fill = status)) +
    geom_col() +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "C. Samples passed vs excluded", x = NULL, y = "n samples", fill = NULL) +
    theme_qc()
}

## -- D. median sequencing depth (final) ----
if (!is.na(C_med)) {
  panels$D <- ggplot(qc[is_pass], aes(Dataset, get(C_med))) +
    geom_boxplot(outlier.size = 0.6, fill = "#e9eef2") +
    geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
    scale_y_log10(labels = scales::comma) +
    labs(title = "D. Median nCount (final) per sample", x = NULL, y = "median UMI / cell (log10)") +
    theme_qc()
}

## -- E. low-complexity flags ----
if (!is.na(C_lowc)) {
  lc <- qc[is_pass, .(flagged = sum(as.logical(get(C_lowc)), na.rm = TRUE),
                      total = .N), by = Dataset]
  panels$E <- ggplot(lc, aes(Dataset, flagged)) +
    geom_col(fill = "#d9784e") +
    geom_text(aes(label = sprintf("%d/%d", flagged, total)), vjust = -0.3, size = 2.8) +
    labs(title = "E. Low-complexity-flagged samples", x = NULL, y = "n samples") +
    theme_qc()
}

## -- assemble ----
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combo <- Reduce(`/`, panels)
  combo <- combo + plot_annotation(title = "AML niche_net  --  per-sample QC overview")
  save_fig(combo, "f01_qc_overview", w = 11, h = 4 * length(panels), subdir = "01_qc")
} else {
  message("[f01] patchwork not installed -> saving panels individually. install.packages('patchwork')")
  for (nm in names(panels)) save_fig(panels[[nm]], paste0("f01_qc_overview_", nm), w = 11, h = 4, subdir = "01_qc")
}
message("[f01] done.")
