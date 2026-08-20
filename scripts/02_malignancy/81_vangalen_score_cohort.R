#!/usr/bin/env Rscript
# 81_vangalen_score_cohort.R ----
# Score all 214 samples with the van Galen malignant-state signatures from 80, set the threshold on
# healthy donors the signature has never seen, and report the sensitivity that threshold buys.
#
# THE ORDER OF OPERATIONS IS THE POINT, and it is the opposite of the CNV gate's.
#
#   The inferCNV gate takes the 95th percentile of REFERENCE-CELL burden and calls anything above
#   it malignant. Nothing in that construction looks at a healthy marrow, so the false-positive
#   rate is whatever it turns out to be -- measured afterwards at 0.186 against a nominal 0.05,
#   and inverted on negative controls (healthy median malignant_frac 0.108 vs autologous AML 0.063).
#
#   Here the operating point is DEFINED as the healthy 95th percentile, so FPR 0.05 holds by
#   construction, and the open question becomes sensitivity -- which is reported, not assumed.
#   A gate can be honest about at most one of the two; this picks the one that protects the claim.
#
# THE CALIBRATION DONORS ARE NOT THE TRAINING DONORS. GSE116256's 4 healthy donors ARE contrast B
# in 80 -- they helped select the genes. Calibrating on them would be circular. The 33 healthy
# donors in the other 4 datasets were never seen, and they are 10x while van Galen is Seq-Well, so
# they test cross-platform transfer and calibration in the same measurement.
#
# THE AXIS HAS ONE BIN OF VALIDITY. Only the HSC and Prog signatures cleared 80's held-out AUC gate
# (0.86 and 0.76); Mono (0.68) and cDC (0.49) did not, and GMP/ProMono produced no consistent genes
# at all. That is a real result -- a malignant monocyte transcriptionally resembles a normal
# monocyte -- and it means this axis speaks only for the stem/progenitor compartment. Scores are
# written for every cell because the distribution outside the domain is diagnostic, but the CALL is
# made only in VG_CALL_BINS. A call in a bin where the signature was never shown to work is not
# evidence, and dressing one up as evidence is the failure this whole detour exists to avoid.
#
# WHAT COUNTS AS A TRUE POSITIVE. van Galen genotyped single cells: 939 carry a mutant transcript.
# Those are malignant beyond argument. Cells with only wild-type transcripts are NOT true negatives
# -- allelic dropout means a malignant cell routinely fails to show its mutant read -- so they are
# reported separately and never used as negatives. Sensitivity is measured on mut+ cells alone.
#
# OUTPUT: VG_SCORE_DIR/<dataset>/<sample>__vgmalig_percell.csv
#         results/tables/02_malignancy/vg_calibration.csv
#         results/tables/02_malignancy/vg_sample_summary.csv
#
#   Rscript scripts/02_malignancy/81_vangalen_score_cohort.R [--dataset X] [--limit N] [--force]

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(UCell) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = ""),
  make_option("--limit",   type = "integer",   default = 0L),
  make_option("--force",   action = "store_true", default = FALSE)
)))

## ---------------------------------------------------------------------------------------------
## signatures
## ---------------------------------------------------------------------------------------------
stopifnot(file.exists(VG_SIG_TSV))
SIG <- fread(VG_SIG_TSV)
UPS <- SIG[dir == "up"]
FEATS <- split(UPS$gene, UPS$signature)
stopifnot(length(FEATS) > 0)
message(sprintf("[0] %d signatures: %s", length(FEATS),
                paste(sprintf("%s(%d)", names(FEATS), lengths(FEATS)), collapse = " ")))
# the two that cleared the AUC gate; their mean is the axis
CORE <- intersect(c("malig_HSC", "malig_Prog"), names(FEATS))
stopifnot(length(CORE) > 0)

R_ALL <- qc_rds_roster(on_extra = "ignore")

## -- THE SCORE MUST BE BUILT FROM THE SAME GENES IN EVERY SAMPLE ------------------------------
# UCell imputes an absent gene to expression 0. That does not error and does not fail anything --
# it just lowers the score in whichever datasets lack the gene. GPX1 is absent from 8 of the 13
# datasets (60 samples) because their ingests used different reference annotations, and two of
# those datasets supply healthy calibration donors. A threshold set on one gene universe and
# applied to another is not a threshold. So the signature is restricted to genes present
# EVERYWHERE, before anything is scored, and what that costs is printed.
gene_universe <- function(roster) Reduce(intersect, lapply(roster[, .SD[1], by = dataset]$rds,
                                                           function(p) rownames(readRDS(p))))
SHARED  <- gene_universe(R_ALL)
dropped <- setdiff(unique(UPS$gene), SHARED)
FEATS   <- lapply(FEATS, function(g) intersect(g, SHARED))
FEATS   <- FEATS[lengths(FEATS) > 0]
CORE    <- intersect(CORE, names(FEATS))
stopifnot(length(CORE) > 0)
cat(sprintf("[0b] shared gene universe over %d datasets: %d signature genes kept, %d dropped%s\n",
            uniqueN(R_ALL$dataset), length(intersect(unique(UPS$gene), SHARED)), length(dropped),
            if (length(dropped)) paste0(" (", paste(dropped, collapse = ", "), ")") else ""))
cat(sprintf("     sizes after restriction: %s\n",
            paste(sprintf("%s(%d)", names(FEATS), lengths(FEATS)), collapse = " ")))

R <- R_ALL
if (nzchar(opt$dataset)) R <- R[dataset %in% strsplit(opt$dataset, ",")[[1]]]
if (opt$limit > 0) R <- head(R, opt$limit)
stopifnot(nrow(R) > 0)

dst_of  <- function(ds, sid) file.path(VG_SCORE_DIR, ds, paste0(sid, "__vgmalig_percell.csv"))
anno_of <- function(ds, sid) file.path(ANNO_RECONCILED_DIR, ds, paste0(sid, "__anno_percell.csv"))
cons_of <- function(ds, sid) file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))

score_one <- function(rds, ds, sid) {
  seu <- readRDS(rds)
  # the healthy flag is read HERE, off the object being scored, so the calibration roster cannot
  # drift from the samples it calibrates on -- the same rule 96 uses
  tp  <- if ("Timepoint" %in% names(seu@meta.data)) as.character(seu@meta.data$Timepoint[1]) else NA_character_
  seu <- NormalizeData(seu, verbose = FALSE)
  dat <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
  # the shared universe came from one sample per dataset; verify it holds for THIS one, because a
  # silently-absent gene is exactly the failure the restriction exists to prevent
  miss <- setdiff(unlist(FEATS, use.names = FALSE), rownames(dat))
  if (length(miss))
    stop(sprintf("%d signature gene(s) absent from this sample (%s) -- the shared universe is wrong",
                 length(miss), paste(head(miss, 5), collapse = ", ")))
  U <- UCell::ScoreSignatures_UCell(dat, features = FEATS, maxRank = CELLSTATE_MAXRANK,
                                    ncores = CELLSTATE_NCORES, name = "")
  A <- as.data.table(as.data.frame(U), keep.rownames = "cell")
  A[, vg_malignant := rowMeans(.SD, na.rm = TRUE), .SDcols = CORE]
  af <- anno_of(ds, sid)
  if (file.exists(af)) A <- merge(A, fread(af, select = c("cell", "hierarchy_bin")), by = "cell", all.x = TRUE)
  else A[, hierarchy_bin := NA_character_]
  cf <- cons_of(ds, sid)
  if (file.exists(cf)) A <- merge(A, fread(cf, select = c("cell", "malignant")), by = "cell", all.x = TRUE)
  else A[, malignant := NA_integer_]
  setnames(A, "malignant", "cnv_malignant")
  A[, `:=`(dataset = ds, sample = sid, timepoint = tp,
           healthy = isTRUE(tp == "Healthy") || isTRUE(is_healthy_sample(sid)))]
  rm(seu, dat, U); gc(verbose = FALSE)
  A[]
}

t0 <- Sys.time()
for (i in seq_len(nrow(R))) {
  ds <- R$dataset[i]; sid <- R$sample[i]; dst <- dst_of(ds, sid)
  if (file.exists(dst) && !opt$force) { message(sprintf("[%d/%d] %s::%s done", i, nrow(R), ds, sid)); next }
  message(sprintf("[%d/%d] %s::%s", i, nrow(R), ds, sid))
  x <- tryCatch(score_one(R$rds[i], ds, sid),
                error = function(e) { message("  [FAIL] ", conditionMessage(e)); NULL })
  if (is.null(x)) next
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, dst)
  if (i == 1) message(sprintf("  [timing] first sample took %.1f min",
                              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------------------------------------------------------------------------
## calibration and validation
## ---------------------------------------------------------------------------------------------
done <- R[file.exists(dst_of(dataset, sample))]
if (!nrow(done)) { message("[!] nothing scored"); quit(save = "no", status = 1) }
P <- rbindlist(lapply(seq_len(nrow(done)), function(i)
  fread(dst_of(done$dataset[i], done$sample[i]))), fill = TRUE)
cat(sprintf("\n================ %d cells across %d samples, %d datasets ================\n",
            nrow(P), uniqueN(P$sample), uniqueN(P$dataset)))
FAIL <- 0L

D <- P[hierarchy_bin %in% VG_CALL_BINS]
cat(sprintf("cells inside the axis's domain (%s): %d (%.1f%% of the cohort)\n",
            paste(VG_CALL_BINS, collapse = "/"), nrow(D), 100 * nrow(D) / nrow(P)))

cat("\n================ 1. DOES IT TRANSFER? healthy score by dataset ================\n")
CAL <- D[healthy == TRUE & dataset != VG_CALIB_EXCLUDE_DATASET]
byds <- CAL[, .(n_samples = uniqueN(sample), n_cells = .N,
                med = round(median(vg_malignant), 4),
                p95 = round(quantile(vg_malignant, 1 - VG_TARGET_HEALTHY_FPR), 4)), by = dataset][order(dataset)]
print(byds)
cat(sprintf("  van Galen is Seq-Well; every dataset above is 10x. UCell ranks WITHIN each cell,\n"))
cat(sprintf("  which is the mechanism that is supposed to make that survivable.\n"))
if (nrow(byds) >= 2) {
  spread <- max(byds$p95) - min(byds$p95)
  okt <- spread <= VG_MAX_CALIB_SPREAD
  if (!okt) FAIL <- FAIL + 1L
  cat(sprintf("  [%s] healthy P95 spread across datasets = %.4f (tolerated %.2f)\n",
              if (okt) "PASS" else "FAIL", spread, VG_MAX_CALIB_SPREAD))
  if (!okt) cat("  A single global threshold is the wrong object -- calibrate per dataset.\n")
}
# THE ACTUAL PLATFORM TEST. Everything above is 10x-vs-10x, so agreement there says the score is
# reproducible across labs, not across chemistries. GSE116256 is the only Seq-Well dataset in the
# cohort and it has its own healthy donors, which makes it the one available cross-platform contrast.
IN <- D[healthy == TRUE & dataset == VG_CALIB_EXCLUDE_DATASET]
if (nrow(IN)) {
  p95_in <- as.numeric(quantile(IN$vg_malignant, 1 - VG_TARGET_HEALTHY_FPR))
  shift  <- abs(p95_in - median(byds$p95))
  cat(sprintf("  %s healthy donors (Seq-Well, %d cells): P95 = %.4f\n",
              VG_CALIB_EXCLUDE_DATASET, nrow(IN), p95_in))
  # Deliberately NOT scored as a gate. It is a property of the data that is MITIGATED below by
  # dataset-matched thresholds, so failing the run on it would be wrong -- but printing "[FAIL]"
  # and then concluding "all gates PASS" is worse than either, so it is labelled for what it is.
  cat(sprintf("  [%s] cross-PLATFORM shift = %.4f (10x median P95 %.4f vs Seq-Well %.4f)\n",
              if (shift <= VG_MAX_CALIB_SPREAD) "no shift" else "SHIFT -- mitigated below",
              shift, median(byds$p95), p95_in))
  if (shift > VG_MAX_CALIB_SPREAD) {
    cat("  The score is NOT platform-invariant. UCell's within-cell ranking survives lab-to-lab\n")
    cat("  variation but not the depth/capture difference between Seq-Well and 10x. A single global\n")
    cat("  threshold is therefore wrong, and the calls below are DATASET-MATCHED wherever possible.\n")
  }
}

cat("\n================ 2. THRESHOLD, SET ON THE NEGATIVE CONTROLS ================\n")
if (nrow(CAL) < VG_MIN_CALIB_CELLS)
  stop(sprintf("only %d calibration cells (need >= %d) -- refusing to set a gate on noise",
               nrow(CAL), VG_MIN_CALIB_CELLS))
THR <- as.numeric(quantile(CAL$vg_malignant, 1 - VG_TARGET_HEALTHY_FPR))
cat(sprintf("  calibration set: %d cells from %d healthy donors in %d datasets, none seen in training\n",
            nrow(CAL), uniqueN(CAL$sample), uniqueN(CAL$dataset)))
cat(sprintf("  pooled 10x threshold = healthy P%.0f = %.4f -> healthy FPR = %.4f by construction\n",
            100 * (1 - VG_TARGET_HEALTHY_FPR), THR, mean(CAL$vg_malignant > THR)))

# A dataset that ships its own healthy donors gets its own threshold: same chemistry, same depth,
# same handling as the AML samples it is applied to. Everything else falls back to the pooled 10x
# value. GSE116256 IS included here -- its healthy donors helped select genes in 80, so this is
# mildly circular, but a platform-mismatched threshold is not the safer error: it silently drives
# sensitivity to zero, which is the more misleading failure of the two.
DTHR <- D[healthy == TRUE, .(n = .N, thr = as.numeric(quantile(vg_malignant, 1 - VG_TARGET_HEALTHY_FPR))),
          by = dataset][n >= VG_MIN_DS_CALIB_CELLS]
cat(sprintf("  dataset-matched thresholds (>= %d healthy cells in-domain):\n", VG_MIN_DS_CALIB_CELLS))
if (nrow(DTHR)) print(DTHR[order(thr)]) else cat("    (none qualify)\n")
# WHAT THE MITIGATION DOES NOT COVER. A dataset with no healthy donors of its own inherits the
# pooled 10x threshold, which is only right if it IS 10x. Given a shift of 0.26 between chemistries,
# a non-10x dataset in this group would have its calls driven almost entirely by platform. The
# cohort is 10x apart from GSE116256 (Seq-Well), and GSE116256 has its own donors -- so the
# fallback group is homogeneous today. It is reported because that is an assumption, not a fact
# the script can check.
fb <- P[!dataset %in% DTHR$dataset, .(n_cells = .N), by = dataset][order(-n_cells)]
cat(sprintf("  the remaining %d datasets use the pooled 10x threshold: %s\n",
            nrow(fb), paste(fb$dataset, collapse = ", ")))
cat(sprintf("  [assumption] all %d are 10x. Seq-Well (GSE116256) is calibrated on its own donors.\n", nrow(fb)))
thr_for <- function(ds) { i <- match(ds, DTHR$dataset); fifelse(is.na(i), THR, DTHR$thr[i]) }
P[, vg_thr := thr_for(dataset)]
P[, vg_call := NA_integer_]
P[hierarchy_bin %in% VG_CALL_BINS, vg_call := as.integer(vg_malignant > vg_thr)]
# D was taken before vg_call existed; rebuild it so every section below sees the same table
D <- P[hierarchy_bin %in% VG_CALL_BINS]

## -- PERSIST THE CALL NEXT TO THE SCORE -------------------------------------------------------
# vg_call and vg_thr were computed here and existed only in memory, so the per-cell files on disk
# carried the SCORE but never the CALL. Any consumer wanting "is this cell malignant by the van
# Galen axis" had either to re-derive the threshold itself -- getting the pooled one wrong for
# Seq-Well, a 0.26 shift -- or to give up: 05_ccc/03_node_features.R asked for vg_call, found no
# such column, and emitted frac_malignant_vg as NA for 100% of rows without failing.
# The threshold is a property of the cohort, so it can only be written after every sample is
# scored; hence a second pass rather than writing it in the scoring loop.
for (i in seq_len(nrow(done))) {
  f <- dst_of(done$dataset[i], done$sample[i])
  x <- fread(f)
  x[, c("vg_thr", "vg_call") := NULL]
  x[, vg_thr := thr_for(done$dataset[i])]
  x[, vg_call := NA_integer_]
  x[hierarchy_bin %in% VG_CALL_BINS, vg_call := as.integer(vg_malignant > vg_thr)]
  fwrite_safe(x, f)
}
message(sprintf("[persist] vg_thr + vg_call written back to %d per-cell files", nrow(done)))

cat("\n================ 3. SENSITIVITY, on genotyped cells ================\n")
# van Galen's own single-cell genotyping. mut+ is malignant beyond argument; wt-only is NOT a
# negative, because allelic dropout makes a malignant cell fail to show its mutant read routinely.
RAW <- file.path(LARGE1_DIR, "00_raw", "public", "GEO", "GSE116256_RAW")
.core <- function(x) sub("^.*_", "", x)
G <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(p) {
  d <- fread(cmd = paste("zcat", shQuote(p)), colClasses = "character", showProgress = FALSE)
  if (!"MutTranscripts" %in% names(d)) return(NULL)
  data.table(smp = sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(p))),
             cell = d$Cell, mut = trimws(d$MutTranscripts),
             wt = if ("WtTranscripts" %in% names(d)) trimws(d$WtTranscripts) else NA_character_)
}), fill = TRUE)
G[, has_mut := nzchar(mut) & mut != "NA"]
G[, has_wt  := nzchar(wt)  & wt  != "NA"]
G[, k := paste(smp, .core(cell))]
V <- P[dataset == "GSE116256"][, k := paste(sample, .core(cell))]
# vg_thr must be carried: without it GD$vg_thr[1] is NULL, sprintf returns character(0), and cat
# prints NOTHING -- the sentence naming both thresholds silently disappeared from the report while
# the numbers around it still printed, leaving a dangling "would be 0.118".
GJ <- merge(G[, .(k, has_mut, has_wt)],
            V[, .(k, hierarchy_bin, vg_malignant, vg_thr, vg_call, cnv_malignant)], by = "k")
GD <- GJ[hierarchy_bin %in% VG_CALL_BINS]
cat(sprintf("  genotyped cells joined: %d | inside the domain: %d | mut+ inside: %d\n",
            nrow(GJ), nrow(GD), GD[has_mut == TRUE, .N]))
if (GD[has_mut == TRUE, .N] >= 20) {
  sens <- mean(GD[has_mut == TRUE]$vg_call == 1L, na.rm = TRUE)
  cat(sprintf("  SENSITIVITY at healthy-FPR %.2f = %.3f  (%d of %d mut+ cells called)\n",
              VG_TARGET_HEALTHY_FPR, sens, GD[has_mut == TRUE & vg_call == 1L, .N], GD[has_mut == TRUE, .N]))
  stopifnot("vg_thr" %in% names(GD), length(GD$vg_thr) > 0)
  cat(sprintf("    using the platform-matched threshold %.4f; under the pooled 10x threshold %.4f it\n",
              GD$vg_thr[1], THR))
  cat(sprintf("    would be %.3f -- the difference is the platform shift, not the signature.\n",
              mean(GD[has_mut == TRUE]$vg_malignant > THR)))
  if (GD[has_mut == TRUE & !is.na(cnv_malignant), .N] >= 20)
    cat(sprintf("  the CNV consensus calls %.3f of the same cells -- the two axes on identical ground truth\n",
                mean(GD[has_mut == TRUE & !is.na(cnv_malignant)]$cnv_malignant == 1L)))
  cat(sprintf("  [context] wt-only cells called: %.3f -- NOT a false-positive rate, dropout makes\n",
              mean(GD[has_mut == FALSE & has_wt == TRUE]$vg_call == 1L, na.rm = TRUE)))
  cat("            many of these genuinely malignant. Reported so it is not mistaken for one.\n")
} else cat("  [SKIP] too few genotyped cells inside the domain to measure sensitivity\n")

cat("\n================ 4. NEGATIVE CONTROL: healthy vs AML, sample level ================\n")
# USE vg_call, NOT A FRESH COMPARISON AGAINST THE POOLED THRESHOLD. This line read
# `mean(vg_malignant > THR)`, which quietly undid the dataset-matched calibration two sections
# above: for GSE116256 (Seq-Well, matched threshold 0.398 vs pooled 0.652) 17 of 24 samples came
# out at exactly frac 0, including all 4 healthy donors. The negative control was then computed on
# the very platform shift the script had just announced it was mitigating, and it still printed
# "all gates PASS". vg_call already carries the per-dataset threshold.
SS <- D[!is.na(vg_call), .(n_cells = .N, frac = mean(vg_call == 1L), healthy = healthy[1]),
        by = .(dataset, sample)]
cat(sprintf("  healthy samples: median frac = %.3f (n=%d)\n", median(SS[healthy == TRUE]$frac), SS[healthy == TRUE, .N]))
cat(sprintf("  AML samples:     median frac = %.3f (n=%d)\n", median(SS[healthy == FALSE]$frac), SS[healthy == FALSE, .N]))
auc_of <- function(score, y) { n1 <- sum(y); n0 <- sum(!y)
  if (n1 < 3 || n0 < 3) return(NA_real_); (sum(rank(score)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
a <- auc_of(SS$frac, SS$healthy == FALSE)
cat(sprintf("  healthy-vs-AML AUC at the sample level = %.3f\n", a))
# THE CONFOUND: healthy samples are concentrated in a few datasets, so a pooled AUC can be driven
# by dataset baseline rather than by disease. Recompute inside datasets that contain BOTH.
both <- SS[, .(nh = sum(healthy), na = sum(!healthy)), by = dataset][nh >= 2 & na >= 2]$dataset
if (length(both)) {
  W <- rbindlist(lapply(both, function(d) data.table(dataset = d,
        n_healthy = SS[dataset == d & healthy == TRUE, .N], n_aml = SS[dataset == d & healthy == FALSE, .N],
        auc = round(auc_of(SS[dataset == d]$frac, SS[dataset == d]$healthy == FALSE), 3))))
  cat("  within-dataset (healthy and AML side by side, no dataset baseline effect):\n")
  print(W)
} else cat("  [note] no dataset holds both healthy and AML samples -- the pooled AUC keeps a dataset confound\n")
cat("  The CNV route scores 0.546 here and is INVERTED on the medians; that is the bar to clear.\n")
oki <- !is.na(a) && a > 0.5 && median(SS[healthy == FALSE]$frac) > median(SS[healthy == TRUE]$frac)
if (!oki) FAIL <- FAIL + 1L
cat(sprintf("  [%s] AML scores higher than healthy, and AUC is above chance\n", if (oki) "PASS" else "FAIL"))

cat("\n================ 5. is the axis orthogonal to CNV? ================\n")
CC <- D[!is.na(cnv_malignant) & !is.na(vg_call)]
if (nrow(CC) > 1000) {
  print(table(vg = CC$vg_call, cnv = CC$cnv_malignant))
  cat(sprintf("  agreement %.3f | phi %.3f\n", mean(CC$vg_call == CC$cnv_malignant),
              suppressWarnings(cor(CC$vg_call, CC$cnv_malignant))))
  cat("  Low correlation is the POINT, not a problem: an axis that merely reproduced the CNV call\n")
  cat("  would add nothing to a cohort that is half cytogenetically normal.\n")
}

dir.create(HIER_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(byds[, `:=`(threshold = THR)], file.path(DIR_MALIGNANCY, "vg_calibration.csv"))
fwrite(SS, file.path(DIR_MALIGNANCY, "vg_sample_summary.csv"))
cat(sprintf("\n[done] %s\n[done] %s\n", file.path(DIR_MALIGNANCY, "vg_calibration.csv"),
            file.path(DIR_MALIGNANCY, "vg_sample_summary.csv")))
if (FAIL > 0) { cat(sprintf("\n[!] %d gate(s) FAILED -- do not use the van Galen axis until resolved.\n", FAIL))
                quit(save = "no", status = 1) }
cat("\nall gates PASS\n")
