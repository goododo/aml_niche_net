#!/usr/bin/env Rscript
# 00_share_eligibility.R ----
# Doublet-pairing redesign (share design), Stage 1 gate (prereg_share/PREREGISTRATION_SHARE_F1.md section 4; user 2026-11-12).
# INPUT  : scripts/dbl/dbl_lanes_share.tsv (analysis lanes); DBL_TAB_DIR/dbl_class_a_composition_share_<ds>.csv (stage 02 --share)
# OUTPUT : DBL_TAB_DIR/dbl_share_stage1_eligibility.csv   per patient: lanes, strict-T and LEUK class-A cells summed over the
#                                                         patient's analysis lanes, eligible (T >= DBL_SHARE_MIN_T, LEUK >= DBL_SHARE_MIN_LEUK)
#          DBL_TAB_DIR/dbl_share_stage1_verdict.csv       n eligible by dataset and age group; the verdict (n >= DBL_SHARE_MIN_N)
# What it does: the pre-registered Stage 1 count -- fewer than DBL_SHARE_MIN_N eligible patients closes the line with BN-1.
# Usage  : Rscript scripts/dbl/00_share_eligibility.R   (only after the freeze-1 tag; stage 02 --share must have run on every dataset)

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
stopifnot(system2("git", c("-C", here::here(), "rev-parse", "-q", "--verify", "refs/tags/prereg-dbl-share-f1"), stdout = FALSE) == 0)
LT <- fread(here::here("scripts", "dbl", DBL_SHARE_LANES_FILE))[role == "analysis"]
LT[, donor := fifelse(dataset == "GSE185381", patient, lane)]     # stage 01's donor id: Sample (GSE185381) or the lane itself
LEUK <- c(DBL_LEUK_STATES, "LEUK_other")
CP <- rbindlist(lapply(unique(LT$dataset), function(d) {
  f <- file.path(DBL_TAB_DIR, paste0("dbl_class_a_composition_share_", d, ".csv")); stopifnot(file.exists(f)); fread(f, colClasses = list(character = c("donor", "lane"))) }))
X <- merge(LT[, .(dataset, patient, age_group, lane, donor)], CP, by = c("lane", "donor"), all.x = TRUE)
stopifnot(LT[, all(paste(lane, donor) %in% CP[, paste(lane, donor)])])   # every analysis lane x patient has class-A cells counted
E <- X[, .(lanes = uniqueN(lane), T = sum(N[compartment == "T"], na.rm = TRUE), LEUK = sum(N[compartment %in% LEUK], na.rm = TRUE),
           classA = sum(N, na.rm = TRUE)), by = .(dataset, patient, age_group)]
E[, eligible := T >= DBL_SHARE_MIN_T & LEUK >= DBL_SHARE_MIN_LEUK]
n <- E[, sum(eligible)]
V <- E[, .(candidates = .N, eligible = sum(eligible)), by = .(dataset, adult = grepl("^adult", age_group))][order(dataset, -adult)]
verdict <- if (n >= DBL_SHARE_MIN_N) sprintf("Stage 1 pass: %d eligible patients (>= %d) -> Stage 2 calibration", n, DBL_SHARE_MIN_N) else
  sprintf("BN-1: lineage-level T-LEUK pairing has no validated null in these data (n = %d < %d)", n, DBL_SHARE_MIN_N)
print(E[order(-eligible, dataset, patient)]); print(V); message(verdict)
fwrite(E, file.path(DBL_TAB_DIR, "dbl_share_stage1_eligibility.csv"))
fwrite(rbind(V[, .(dataset, adult = as.character(adult), candidates, eligible)], data.table(dataset = "ALL", adult = "", candidates = nrow(E), eligible = n))[, verdict := verdict],
       file.path(DBL_TAB_DIR, "dbl_share_stage1_verdict.csv"))
