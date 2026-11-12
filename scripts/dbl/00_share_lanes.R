#!/usr/bin/env Rscript
# 00_share_lanes.R ----
# Doublet-pairing redesign (share design), Stage 1 step S1a (reports/.../redesign_2026-11-11/STAGE1_plan.md; user 2026-11-12).
# INPUT  : results/tables/01_preprocess/00_curated_manifest.csv (disease, timepoint, platform per patient)
#          TAB_DIR/00_ingest/<dataset>_qc_percell.csv.gz (lane, author assignment, author age group; no cell labels)
#          GEO_RAW_DIR/GSE185381_RAW/<lane>_features.tsv.gz (are there hashtag rows?)
# OUTPUT : scripts/dbl/dbl_lanes_share.tsv   one row per lane x patient: dataset, gsm, lane, patient, age_group, chemistry,
#          donors_in_lane, n_assigned, hashed, k1_lane, role (analysis / excluded), reason
# What it does: lists every candidate lane of the share test BY METADATA ONLY -- AML diagnosis marrow of GSE185381 (3'),
#   GSE289435 and Petti 2019 (5') -- with the exclusions fixed in the plan. No composition or cell label is read: screening on
#   compositions belongs to the pre-registered Stage 1 gate (after freeze 1).
# Usage  : Rscript scripts/dbl/00_share_lanes.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
OUT <- here::here("scripts", "dbl", "dbl_lanes_share.tsv")
EXCL_GSM <- sprintf("GSM56137%02d", 47:55)                       # GSE185381 T-ALL co-capture pools (memory gse185381-healthy-tall-cocapture)
EXCL_TP  <- sprintf("GSM56137%02d", 57:59)                       # pool hashed D / EoT / Relapse per the authors' code (2019-11-04-count-1-CREATE.Rmd);
                                                                 #   the deposit labels only the patient and has no HTO rows: timepoints not separable
MAN <- fread(here::here("results", "tables", "01_preprocess", "00_curated_manifest.csv"))
qc  <- function(d, cols) fread(file.path(TAB_DIR, "00_ingest", paste0(d, "_qc_percell.csv.gz")), select = cols)
top <- function(x) names(sort(table(x), decreasing = TRUE))[1]   # majority value (author age group per patient)

## -- GSE185381: every lane that holds author-assigned barcodes of an AML diagnosis patient ----
aml <- MAN[dataset == "GSE185381" & disease == "AML" & timepoint == "Diagnosis", unique(patient_id)]
Q <- qc("GSE185381", c("cell_barcode", "Sample", "Demuxed", "GSM", "Library", "anno_age_group"))
dl <- Q[Demuxed == TRUE, .(donors_in_lane = uniqueN(Sample)), by = Library]
hto <- vapply(unique(Q$Library), function(L) { f <- file.path(GEO_RAW_DIR, "GSE185381_RAW", paste0(L, "_features.tsv.gz"))
  if (!file.exists(f)) return(NA); ft <- fread(f, header = FALSE); any(ft$V3 == "Antibody Capture" & grepl(DBL_HTO_PREFIX, ft$V1)) }, TRUE)
G <- Q[Demuxed == TRUE & Sample %in% aml, .(n_assigned = .N, age_group = top(anno_age_group)), by = .(gsm = GSM, lane = Library, patient = Sample)]
G <- merge(G, dl, by.x = "lane", by.y = "Library")
G[, `:=`(dataset = "GSE185381", chemistry = "10x_3p_v3", hashed = unname(hto[lane]))]
G[, k1_lane := hashed %in% TRUE & donors_in_lane >= 2]          # cross-donor barcodes of a hashed pool can serve K1
G[, `:=`(role = fifelse(gsm %in% c(EXCL_GSM, EXCL_TP), "excluded", "analysis"),
         reason = fifelse(gsm %in% EXCL_GSM, "T-ALL co-capture pool (GSM5613747-55)",
                  fifelse(gsm %in% EXCL_TP, "pool hashed Dx/EoT/Relapse, deposit has no HTO rows: diagnosis cells not separable",
                  fifelse(donors_in_lane >= 2 & !(hashed %in% TRUE), "multiplexed without hashtag rows: author-assigned barcodes only", ""))))]
## a patient's lanes must be different libraries, not one library sequenced twice: barcode-set Jaccard < 0.5
bcs <- Q[Demuxed == TRUE & Sample %in% aml, .(bc = sub("^.*_", "", cell_barcode)), by = .(patient = Sample, lane = Library)]
ov <- bcs[, { L <- unique(lane); if (length(L) < 2) NULL else rbindlist(combn(L, 2, function(p) {
  a <- bc[lane == p[1]]; b <- bc[lane == p[2]]; data.table(l1 = p[1], l2 = p[2], jaccard = length(intersect(a, b)) / length(union(a, b))) }, simplify = FALSE)) }, by = patient]
print(ov)
stopifnot(nrow(ov) == 0 || all(ov$jaccard < 0.5))

## -- GSE289435 and Petti 2019: one lane per patient, no hashtags ----
O <- rbindlist(lapply(c("GSE289435", "Petti2019"), function(d) {
  q <- qc(d, c("Sample", "Patient_ID"))[, .(n_assigned = .N), by = .(lane = as.character(Sample), patient = as.character(Patient_ID))]
  m <- MAN[dataset == d, .(patient = as.character(patient_id), gsm = gsm_id, disease, timepoint, chemistry = platform)]
  x <- merge(q, m, by = "patient", all.x = TRUE)
  x[, `:=`(dataset = d, age_group = "adult", donors_in_lane = 1L, hashed = FALSE, k1_lane = FALSE,   # both are adult cohorts
           role = fifelse(disease == "AML" & timepoint == "Diagnosis", "analysis", "excluded"),
           reason = fifelse(disease == "AML" & timepoint == "Diagnosis", "5' chemistry: M1 untested (Stage 2 lane gates decide)",
                            paste0("not AML diagnosis (", disease, ")")))][, !c("disease", "timepoint")] }))

## -- write, with checks ----
TAB <- rbind(G, O, use.names = TRUE)[, .(dataset, gsm, lane, patient, age_group, chemistry, donors_in_lane, n_assigned, hashed, k1_lane, role, reason)]
setorder(TAB, dataset, patient, lane)
stopifnot(!anyDuplicated(TAB[, .(lane, patient)]), !anyNA(TAB$role), all(aml %in% TAB$patient))
print(TAB[, .(lanes = .N, patients = uniqueN(patient)), by = .(dataset, role)])
print(TAB[role == "analysis", .(patients = uniqueN(patient)), by = .(dataset, age_group)])
fwrite(TAB, OUT, sep = "\t"); ts("[write] ", OUT, " (", nrow(TAB), " rows)")
