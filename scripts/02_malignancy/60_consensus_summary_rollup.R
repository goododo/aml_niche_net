#!/usr/bin/env Rscript
# 60_consensus_summary_rollup.R ----
# Roll every per-sample __consensus_summary.csv (all datasets) into one table + cross-cohort stats.
# Read-only; run from PROJECT ROOT.

suppressPackageStartupMessages({ library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))

DIR <- DIR_MALIGNANCY   # /FAST/.../results/tables/02_malignancy
f <- list.files(DIR, pattern = "__consensus_summary\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(f) > 0)
message(sprintf("[0] %d summary files under %s", length(f), DIR))

all <- rbindlist(lapply(f, fread), fill = TRUE)
setorder(all, dataset, sample)
out_all <- file.path(DIR, "ALL_consensus_summary.csv")
fwrite(all, out_all)

cat(sprintf("\n================ %d samples across %d datasets ================\n", nrow(all), uniqueN(all$dataset)))

cat("\n-- per dataset: n samples, arms, tier mix, median malignant_frac --\n")
per_ds <- all[, .(
  n              = .N,
  with_numbat    = sum(grepl("numbat", arms)),
  A_concordant   = sum(evidence_tier == "A_concordant"),
  B_multi_partial= sum(evidence_tier == "B_multi_partial"),
  C_single       = sum(evidence_tier == "C_single"),
  med_mal_frac   = round(median(malignant_frac, na.rm = TRUE), 3),
  med_cells      = as.integer(median(n_qc_cells, na.rm = TRUE))
), by = dataset][order(dataset)]
print(per_ds)

cat("\n-- overall tier distribution --\n")
print(all[, .N, by = evidence_tier][order(-N)])

cat("\n-- arms coverage --\n")
print(all[, .N, by = arms][order(-N)])

cat("\n-- malignant_frac distribution (all samples) --\n")
print(round(quantile(all$malignant_frac, probs = c(0, .25, .5, .75, .9, 1), na.rm = TRUE), 3))

cat("\n-- samples flagged for attention (very low or extreme malignant_frac, or tiny n) --\n")
flagged <- all[malignant_frac < 0.01 | malignant_frac > 0.9 | n_qc_cells < 200,
               .(dataset, sample, n_qc_cells, arms, evidence_tier, malignant_frac)]
print(flagged[order(dataset, sample)])

cat(sprintf("\n[done] combined table -> %s\n", out_all))
