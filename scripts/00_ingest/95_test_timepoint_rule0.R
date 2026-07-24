#!/usr/bin/env Rscript
# 95_test_timepoint_rule0.R ----
# DRY RUN for the canonicalize_timepoint() rule-0 fix (sample-level healthy donors). Reads the
# CURRENT Timepoint out of every QC object, applies the new rule to the same sample names, and
# reports exactly which samples would change -- so the regex can be audited BEFORE paying for a
# re-ingest. Writes nothing except a report CSV.
#
# What to look for:
#   * WOULD_CHANGE rows: should be exactly the known normal-marrow samples (Chen2023 NBM1-4,
#     Petti2019 ND_*/Normal_sorted_*, GSE116256 BM<n>). Anything else is a FALSE POSITIVE.
#   * ALREADY_HEALTHY rows: samples the dataset-level rules already handled (sanity check).
#   * A malignant_frac column is joined where available: a "healthy" sample with a high called
#     malignant fraction is an inferCNV false positive, which is the whole point of this fix.
#
#   Rscript scripts/00_ingest/95_test_timepoint_rule0.R

suppressPackageStartupMessages({ library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))

stopifnot(exists("is_healthy_sample"))          # fails loudly if the patched config isn't deployed

rds <- list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
rds <- rds[!grepl("/_", rds)]
message(sprintf("[0] inspecting %d QC objects ...", length(rds)))

cur <- rbindlist(lapply(rds, function(f) {
  m <- tryCatch(readRDS(f)@meta.data, error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.table(dataset = basename(dirname(f)),
             sample  = sub("\\.rds$", "", basename(f)),
             patient = if ("Patient_ID" %in% names(m)) as.character(m$Patient_ID[1]) else NA_character_,
             tp_now  = if ("Timepoint" %in% names(m)) as.character(m$Timepoint[1]) else NA_character_,
             tp_raw  = if ("Timepoint_detail" %in% names(m)) as.character(m$Timepoint_detail[1]) else NA_character_,
             disease = if ("Disease_state" %in% names(m)) as.character(m$Disease_state[1]) else NA_character_,
             n_cells = ncol(m))
}), fill = TRUE)

# what rule 0 alone would say (suppress the per-hit message here; we tabulate instead)
cur[, rule0_healthy := mapply(function(s, p) isTRUE(is_healthy_sample(s, p)), sample, patient)]
cur[, tp_new := suppressMessages(mapply(function(d, r, di, s, p)
      canonicalize_timepoint(d, r, di, sample_id = s, patient_id = p),
      dataset, tp_raw, disease, sample, patient))]
cur[, status := fifelse(tp_now == tp_new, "unchanged",
                 fifelse(tp_new == "Healthy" & tp_now != "Healthy", "WOULD_CHANGE -> Healthy", "OTHER_CHANGE"))]

# join the called malignant fraction where we have it (the evidence this matters)
mf <- file.path(DIR_MALIGNANCY, "ALL_consensus_summary.csv")
if (file.exists(mf)) cur[fread(mf)[, .(dataset, sample, malignant_frac)],
                         malignant_frac := i.malignant_frac, on = c("dataset", "sample")]

cat("\n================ status summary ================\n")
print(cur[, .N, by = status][order(-N)])

cat("\n================ WOULD CHANGE to Healthy ================\n")
cat("(these are healthy donors currently labelled as AML timepoints)\n")
print(cur[status == "WOULD_CHANGE -> Healthy"][order(dataset, sample),
      .(dataset, sample, tp_now, tp_raw, n_cells, malignant_frac)])

cat("\n================ already Healthy (dataset-level rules) ================\n")
print(cur[tp_now == "Healthy", .N, by = dataset])

cat("\n================ FALSE-POSITIVE CHECK ================\n")
cat("Every row below matched the healthy name patterns. Confirm each really is normal marrow;\n")
cat("if any is an AML sample, tighten TP_HEALTHY_SAMPLE_PATTERNS before re-ingesting.\n")
print(cur[rule0_healthy == TRUE][order(dataset, sample),
      .(dataset, sample, patient, tp_now, tp_raw, disease, malignant_frac)])

cat("\n================ MISSED-HEALTHY CHECK ================\n")
cat("Samples NOT matched by rule 0 whose name still looks donor-ish -- extend the patterns if real:\n")
print(cur[rule0_healthy == FALSE & tp_now != "Healthy" &
          grepl("normal|healthy|control|donor|^BM|NBM|^ND", sample, ignore.case = TRUE)][
      order(dataset, sample), .(dataset, sample, tp_now, malignant_frac)])

out <- file.path(FAST_DIR, "results", "tables", "00_ingest", "timepoint_rule0_dryrun.csv")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(cur, out)
cat("\n[done] full report -> ", out, "\n", sep = "")
