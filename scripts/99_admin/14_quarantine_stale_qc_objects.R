#!/usr/bin/env Rscript
# =============================================================================
# 14_quarantine_stale_qc_objects.R
#
# Move QC objects that are NOT in the current PASS roster out of QC_RDS_DIR.
#
# WHY THIS EXISTS. Eight scripts build their sample list with
# list.files(QC_RDS_DIR) -- among them the BMM projection, the cNMF stages, and
# BOTH malignancy calibration diagnostics. That makes "present on disk" the
# operative definition of "in the cohort". The v2 re-ingest dropped samples the
# v1 ingest had written (Chen2023 sorted fractions, the v1 GSE185381 roster,
# Petti2019 ND_/Normal_sorted_, E-MTAB -BLD), but their .rds files stayed where
# they were, so the old definition kept returning them. J3 duly processed 74 of
# them and reported success.
#
# NOTHING IS DELETED. Objects move to QC_RDS_DIR/../_stale_qc_objects_v1/ with a
# manifest, so v1 stays available for comparison and the move is one mv away from
# being undone. The selection is made against the ROSTER, never against mtime --
# a fresh timestamp on a sample that left the cohort is still stale.
#
# Usage: Rscript scripts/99_admin/14_quarantine_stale_qc_objects.R [--apply]
#        (default is a dry run)
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(optparse); library(here)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--apply", action = "store_true", default = FALSE,
              help = "actually move the files (default: dry run)")
)))

source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))

QUARANTINE <- file.path(dirname(QC_RDS_DIR), "_stale_qc_objects_v1")

R <- qc_rds_roster(on_extra = "ignore")          # the roster IS the QC report
on_disk <- list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
on_disk <- on_disk[!grepl("/_", on_disk)]

S <- data.table(rds = setdiff(on_disk, R$rds))
message("[1] roster (PASS) = ", nrow(R), " | on disk = ", length(on_disk),
        " | not in roster = ", nrow(S))
if (!nrow(S)) { message("[done] QC_RDS_DIR already matches the roster; nothing to do."); quit(save = "no") }

S[, `:=`(dataset = basename(dirname(rds)),
         sample  = sub("\\.rds$", "", basename(rds)),
         mtime   = file.mtime(rds),
         mb      = round(file.size(rds) / 1e6, 1))]

# Say WHY each one is out, from the QC report -- "absent from the report" and
# "present but P2-exclude" are different facts and get reported as such.
Q <- fread(file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv"))
Q[, k := paste(dataset, Sample)]
S[, k := paste(dataset, sample)]
S[, why := Q$status[match(k, Q$k)]]
S[is.na(why), why := "absent from the QC report (dropped by the v2 ingest)"]

message("\n[2] by dataset and reason:")
print(S[, .(n = .N, oldest = min(mtime), newest = max(mtime), gb = round(sum(mb) / 1e3, 2)),
        by = .(dataset, why)][order(-n)])

if (!opt$apply) {
  message("\n[dry run] re-run with --apply to move ", nrow(S), " object(s) to:\n  ", QUARANTINE)
  quit(save = "no")
}

dir.create(QUARANTINE, recursive = TRUE, showWarnings = FALSE)
S[, dest := file.path(QUARANTINE, dataset, basename(rds))]
for (d in unique(dirname(S$dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
S[, moved := file.rename(rds, dest)]
if (!all(S$moved)) stop("failed to move ", sum(!S$moved), " file(s); nothing else changed")

fwrite_safe(S[, .(dataset, sample, why, mtime, mb, from = rds, to = dest)],
            file.path(QUARANTINE, "quarantine_manifest.csv"))
message("\n[3] moved ", nrow(S), " object(s) -> ", QUARANTINE)

# Re-assert. The point of the exercise is that qc_rds_roster() stops complaining,
# not that some files moved.
invisible(qc_rds_roster(on_extra = "error"))
message("[done] QC_RDS_DIR now matches the PASS roster exactly (", nrow(R), " samples).")
