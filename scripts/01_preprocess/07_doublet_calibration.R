#!/usr/bin/env Rscript
# =============================================================================
# 07_doublet_calibration.R
#
# THE QUESTION: doublet removal currently takes 1.6x to 3.9x the expected number
# of cells, and the multiplier tracks library size (small samples lose the most,
# relative to expectation). Composition is the FGW node marginal, so a
# size-correlated removal rate injects a technical gradient into exactly the
# quantity H2 measures. Before choosing a consensus rule we need to know WHICH
# caller produces the excess and whether either is calibrated.
#
# The per-sample QC report stores only the merged n_doublet, so union vs
# intersection cannot be compared from existing output. This script runs both
# callers separately on samples spanning the size range and records:
#   rate_sc, rate_df, rate_union, rate_inter, jaccard  -- each against expected
#
# It CHANGES NOTHING. It writes a table and a figure; the consensus rule stays
# whatever config says.
#
# WHY IT MATTERS WHICH WAY WE ERR: a doublet is not noise in a CCC graph. An
# HSC-monocyte doublet is a synthetic cell co-expressing both programs, which is
# precisely the signal CellChat reads as communication -- doublets MANUFACTURE
# edges. A discarded real cell only adds variance. The loss function is
# asymmetric, so "be conservative" is not automatically the safe choice here.
#
# Usage: Rscript scripts/01_preprocess/07_doublet_calibration.R [--n_per_bin 2]
# =============================================================================

options(aml_qc.source_only = TRUE)
suppressPackageStartupMessages({
  library(here); library(data.table); library(optparse); library(ggplot2); library(scales)
})
suppressMessages(source(here::here("scripts", "01_preprocess", "03_per_sample_qc.R")))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--n_per_bin", type = "integer", default = 2L,
              help = "samples drawn per size bin [2]"),
  make_option("--dbr_sd", type = "double", default = -1,
              help = "scDblFinder dbr.sd; -1 = package default (40%% of dbr) [-1]")
)))
DBR_SD <- if (opt$dbr_sd < 0) NULL else opt$dbr_sd

set.seed(SEED)

# ---------------------------------------------------------------------------
# [1] Choose samples spanning the size range
# ---------------------------------------------------------------------------
R <- fread(file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv"))[status == "PASS" & n_after_mad > 0]
R[, size_bin := cut(n_after_mad, breaks = c(0, 1e3, 3e3, 1e4, Inf),
                    labels = c("<1k", "1-3k", "3-10k", ">10k"))]
sel <- R[, .SD[sample(.N, min(.N, opt$n_per_bin))], by = size_bin][
        , .(dataset, Sample, n_after_mad, dbl_rate_exp, dbl_rate_obs, size_bin)]
setorder(sel, n_after_mad)
message("[1] ", nrow(sel), " samples selected across ", uniqueN(sel$size_bin), " size bins")
print(sel)

# ---------------------------------------------------------------------------
# [2] Run both callers separately
# ---------------------------------------------------------------------------
# scDblFinder is re-called here rather than reused from 03, because dbr.sd is the
# parameter under test and 03 does not expose it.
call_sc <- function(obj, rate, dbr.sd) {
  sce <- Seurat::as.SingleCellExperiment(obj)
  sce <- if (is.null(dbr.sd)) scDblFinder::scDblFinder(sce, dbr = rate)
         else                 scDblFinder::scDblFinder(sce, dbr = rate, dbr.sd = dbr.sd)
  setNames(as.character(sce$scDblFinder.class) == "doublet", colnames(sce))
}

res <- vector("list", nrow(sel))
for (i in seq_len(nrow(sel))) {
  ds <- sel$dataset[i]; sm <- sel$Sample[i]
  f  <- file.path(QC_RDS_DIR, ds, paste0(sm, ".rds"))
  if (!file.exists(f)) { message("[skip] no QC object: ", f); next }
  message(sprintf("[2] %d/%d  %s :: %s  (%d cells)", i, nrow(sel), ds, sm, sel$n_after_mad[i]))

  obj  <- readRDS(f)
  # The QC object is POST doublet removal, so the callers must be re-run on it
  # only to compare their relative behaviour -- absolute rates are therefore a
  # LOWER bound. Recorded explicitly so the numbers are not over-read.
  rate <- sel$dbl_rate_exp[i]
  sc <- tryCatch(call_sc(obj, rate, DBR_SD), error = function(e) { message("  sc fail: ", conditionMessage(e)); NULL })
  df <- tryCatch(call_doubletfinder(obj, rate), error = function(e) { message("  df fail: ", conditionMessage(e)); NULL })
  if (is.null(sc) || is.null(df)) next

  sc <- sc[colnames(obj)]; df <- df[colnames(obj)]; n <- length(sc)
  res[[i]] <- data.table(
    dataset = ds, Sample = sm, n_cells = n, size_bin = sel$size_bin[i], rate_exp = rate,
    rate_sc    = sum(sc) / n,
    rate_df    = sum(df) / n,
    rate_union = sum(sc | df) / n,
    rate_inter = sum(sc & df) / n,
    jaccard    = if (sum(sc | df) > 0) sum(sc & df) / sum(sc | df) else NA_real_)
  print(res[[i]])
  rm(obj); gc(verbose = FALSE)
}

D <- rbindlist(res)
if (!nrow(D)) stop("no samples produced calls")
for (k in c("sc", "df", "union", "inter"))
  D[[paste0("ratio_", k)]] <- D[[paste0("rate_", k)]] / D$rate_exp

fwrite_safe(D, file.path(DIR_PREPROCESS, "07_doublet_calibration.csv"))

# ---------------------------------------------------------------------------
# [3] Verdict
# ---------------------------------------------------------------------------
message("\n[3] observed / expected, by caller")
print(D[, .(n_cells, rate_exp = round(rate_exp, 4),
            sc = round(ratio_sc, 2), df = round(ratio_df, 2),
            union = round(ratio_union, 2), inter = round(ratio_inter, 2),
            jaccard = round(jaccard, 2)), by = .(dataset, Sample)])

message("\n    median ratio   sc=", round(D[, median(ratio_sc)], 2),
        "  df=", round(D[, median(ratio_df)], 2),
        "  union=", round(D[, median(ratio_union)], 2),
        "  inter=", round(D[, median(ratio_inter)], 2))

# The decisive number: does the multiplier depend on library size? A calibrated
# rule has slope ~0 against log(n_cells). A size-dependent one puts a technical
# gradient into composition, which is the FGW node marginal.
if (nrow(D) >= 4) {
  message("\n    Spearman(ratio, n_cells) -- near 0 means size-independent:")
  for (k in c("sc", "df", "union", "inter"))
    message(sprintf("      %-6s rho = %+.2f", k,
                    suppressWarnings(cor(D[[paste0("ratio_", k)]], D$n_cells, method = "spearman"))))
}

L <- melt(D, id.vars = c("dataset", "Sample", "n_cells"),
          measure.vars = paste0("ratio_", c("sc", "df", "union", "inter")),
          variable.name = "rule", value.name = "ratio")
L[, rule := factor(sub("^ratio_", "", rule), levels = c("df", "inter", "sc", "union"))]

p <- ggplot(L, aes(n_cells, ratio, colour = rule)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, formula = y ~ x) +
  scale_x_log10(labels = label_comma()) + scale_y_log10() +
  labs(title = "Doublet calls relative to expectation, by rule and library size",
       subtitle = paste("Flat = calibrated. A sloped line puts a library-size artefact into",
                        "composition,\nwhich is the FGW node marginal."),
       x = "cells in sample (log)", y = "observed / expected (log)", colour = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "top")
for (ext in c("png", "pdf"))
  ggsave(file.path(FIG_PREPROCESS, paste0("07_doublet_calibration.", ext)), p,
         width = 6.5, height = 4.5, dpi = 200)

message("\n[done] ", file.path(DIR_PREPROCESS, "07_doublet_calibration.csv"))
