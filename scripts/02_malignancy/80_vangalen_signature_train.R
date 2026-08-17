#!/usr/bin/env Rscript
# 80_vangalen_signature_train.R ----
# Derive per-bin MALIGNANT-STATE signatures from van Galen 2019, as an axis that does not depend
# on copy number.
#
# WHY A SECOND MALIGNANCY AXIS AT ALL. The inferCNV route has two measured problems that no amount
# of threshold tuning fixes:
#
#   (1) IT IS INVERTED ON NEGATIVE CONTROLS. Healthy donors call a HIGHER median malignant_frac
#       (0.108) than autologous-route AML samples (0.063). Healthy-donor FPR is 0.186 against the
#       0.05 the P95 threshold is supposed to deliver.
#   (2) IT IS STRUCTURALLY BLIND TO ~HALF OF AML. inferCNV infers copy-number change. Roughly half
#       of AML is cytogenetically normal -- GSE185991 (29 samples) is NPM1-mutant, the textbook
#       normal-karyotype subtype, and it correctly returns ~0. "No CNV" is not "not malignant",
#       but the pipeline currently has no way to say so.
#
# van Galen's PredictionRefined column is a malignant/normal/unclear call made WITHOUT copy number,
# so a signature learned from it is orthogonal to inferCNV on both counts.
#
# TWO CONTRASTS, AND A GENE MUST SURVIVE BOTH.
#
#   A  WITHIN SAMPLE   AML malignant vs AML normal, same donor, same library, same platform.
#                      Donor and batch cancel. But most AML samples have almost no NORMAL myeloid
#                      cells left, so at 10 cells a side only 4-10 samples per bin are testable --
#                      low power, and a fluke can clear a consistency filter over few samples.
#   B  VS HEALTHY      AML malignant vs the 4 healthy donors inside GSE116256 (4,994 normal cells,
#                      same dataset, same protocol). Ample power; donor is NOT controlled.
#
# Each contrast fails in a way the other does not: A is underpowered, B is donor-confounded. A gene
# that clears both cannot be a donor effect (B's confounder would not survive A) and is unlikely to
# be a low-power fluke (A's weakness is covered by B). Neither contrast alone would justify a
# signature that then gets applied to 214 samples.
#
# A GENE EARNS ITS PLACE BY CONSISTENCY, NOT EFFECT SIZE. Within contrast A the filter is "up in at
# least VG_MIN_FRAC_SAMPLES of the testable samples", so one donor with an enormous fold change
# cannot install its own genes. Leave-one-sample-out stability is reported for the same reason.
#
# WHAT THIS SCRIPT DOES NOT DO: it does not decide whether the signature transfers. That is 81's
# job, and it is settled on healthy donors from other datasets and other platforms -- data this
# script never sees.
#
# OUTPUT: scripts/02_malignancy/vangalen_malignant_signatures.tsv  (signature, gene, weight, bin)
#         results/tables/02_malignancy/vangalen_signature_qc.csv
#
#   Rscript scripts/02_malignancy/80_vangalen_signature_train.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })

RAW   <- file.path(LARGE1_DIR, "00_raw", "public", "GEO", "GSE116256_RAW")
DS    <- "GSE116256"
OUT   <- file.path(SCRIPTS_DIR, "02_malignancy", "vangalen_malignant_signatures.tsv")
QC_CSV <- file.path(DIR_MALIGNANCY, "vangalen_signature_qc.csv")
.core <- function(x) sub("^.*_", "", x)

## ---------------------------------------------------------------------------------------------
## labels: van Galen's own malignant call, joined to our bin vocabulary
## ---------------------------------------------------------------------------------------------
A <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(p) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(p)))
  d <- fread(cmd = paste("zcat", shQuote(p)), colClasses = "character", showProgress = FALSE)
  if (!"PredictionRefined" %in% names(d)) return(NULL)
  data.table(smp = s, cell = d$Cell, lab = trimws(d$PredictionRefined))
}), fill = TRUE)[lab %in% c("malignant", "normal")]

rec <- rbindlist(lapply(list.files(file.path(ANNO_RECONCILED_DIR, DS), full.names = TRUE), function(p)
  fread(p, select = c("cell", "hierarchy_bin"))[, smp := sub("__anno_percell\\.csv$", "", basename(p))][]),
  fill = TRUE)
A[, k := paste(smp, .core(cell))]; rec[, k := paste(smp, .core(cell))]
L <- merge(A[, .(k, smp, lab)], rec[, .(k, cell, hierarchy_bin)], by = "k")
L <- L[!is.na(hierarchy_bin) & nzchar(hierarchy_bin)]
message(sprintf("[0] %d labelled cells over %d samples, %d bins",
                nrow(L), uniqueN(L$smp), uniqueN(L$hierarchy_bin)))

## which (sample, bin) cells can be tested at all
CT <- dcast(L[, .N, by = .(smp, hierarchy_bin, lab)], smp + hierarchy_bin ~ lab,
            value.var = "N", fill = 0L)
CT <- CT[malignant >= VG_MIN_CELLS_PER_SIDE & normal >= VG_MIN_CELLS_PER_SIDE]
BINS <- CT[, .N, by = hierarchy_bin][N >= VG_MIN_SAMPLES_PER_BIN]$hierarchy_bin
message(sprintf("[1] testable (sample,bin) pairs: %d | bins with >= %d such samples: %s",
                nrow(CT), VG_MIN_SAMPLES_PER_BIN, paste(BINS, collapse = ", ")))
if (!length(BINS)) stop("no bin has enough paired samples")

## ---------------------------------------------------------------------------------------------
## contrast A: per-sample, within-bin log fold change (donor-controlled)
## ---------------------------------------------------------------------------------------------
R <- qc_rds_roster(on_extra = "ignore")[dataset == DS]
HEALTHY <- grep(VG_HEALTHY_REGEX, unique(L$smp), value = TRUE)
message(sprintf("[2] healthy donors inside %s for contrast B: %s", DS, paste(HEALTHY, collapse = ", ")))

load_sample <- function(s) {
  rp <- R[sample == s]$rds
  if (!length(rp)) return(NULL)
  seu <- NormalizeData(readRDS(rp), verbose = FALSE)
  SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
}
bin_stats <- function(dat, cells) {
  list(mu = Matrix::rowMeans(dat[, cells, drop = FALSE]),
       det = Matrix::rowMeans(dat[, cells, drop = FALSE] > 0))
}

FC <- list(); MAL <- list()
for (s in unique(L$smp)) {
  dat <- tryCatch(load_sample(s), error = function(e) NULL)
  if (is.null(dat)) { message("  [skip] no QC object: ", s); next }
  ls <- L[smp == s][cell %in% colnames(dat)]
  for (b in intersect(BINS, unique(ls$hierarchy_bin))) {
    mc <- ls[hierarchy_bin == b & lab == "malignant"]$cell
    nc <- ls[hierarchy_bin == b & lab == "normal"]$cell
    # contrast A needs both sides in the SAME sample
    if (length(mc) >= VG_MIN_CELLS_PER_SIDE && length(nc) >= VG_MIN_CELLS_PER_SIDE) {
      a <- bin_stats(dat, mc); n <- bin_stats(dat, nc)
      FC[[paste(s, b)]] <- data.table(gene = rownames(dat), smp = s, hierarchy_bin = b,
                                      lfc = a$mu - n$mu, ddet = a$det - n$det)
    }
    # keep the per-(sample,bin) profiles that contrast B needs, malignant and healthy-normal alike
    grp <- if (s %in% HEALTHY) "healthy_normal" else NA_character_
    if (!is.na(grp) && length(nc) >= VG_MIN_CELLS_PER_SIDE) {
      n <- bin_stats(dat, nc)
      MAL[[paste(s, b, "H")]] <- data.table(gene = rownames(dat), smp = s, hierarchy_bin = b,
                                            grp = "healthy_normal", mu = n$mu, det = n$det)
    }
    if (is.na(grp) && length(mc) >= VG_MIN_CELLS_PER_SIDE) {
      a <- bin_stats(dat, mc)
      MAL[[paste(s, b, "M")]] <- data.table(gene = rownames(dat), smp = s, hierarchy_bin = b,
                                            grp = "aml_malignant", mu = a$mu, det = a$det)
    }
  }
  rm(dat); gc(verbose = FALSE)
  message(sprintf("  [profile] %s", s))
}
FC  <- rbindlist(FC, fill = TRUE)
MAL <- rbindlist(MAL, fill = TRUE)
stopifnot(nrow(FC) > 0, nrow(MAL) > 0)

## contrast B: AML-malignant vs healthy-donor normal, per bin, sample means so a big sample cannot
## dominate the comparison by cell count alone
B_STAT <- MAL[, .(mu = mean(mu), det = mean(det), n_smp = uniqueN(smp)), by = .(gene, hierarchy_bin, grp)]
B_STAT <- dcast(B_STAT, gene + hierarchy_bin ~ grp, value.var = c("mu", "det", "n_smp"))
B_STAT[, `:=`(lfc_B = mu_aml_malignant - mu_healthy_normal,
              ddet_B = det_aml_malignant - det_healthy_normal)]
message(sprintf("[3] contrast B built over %d bins", uniqueN(B_STAT$hierarchy_bin)))

## ---------------------------------------------------------------------------------------------
## consistency filter -> signature
## ---------------------------------------------------------------------------------------------
# min_samples is a PARAMETER, not the config constant. The leave-one-out check below re-runs this
# on n-1 samples; with the constant hard-coded, every gene in a 4-sample bin was filtered out for
# having only 3 samples and "retention" came back as exactly 0 -- a property of the test, not of
# the signature. Any stability check that gets stricter as it removes data measures itself.
build <- function(D, direction = 1, min_samples = VG_MIN_SAMPLES_PER_BIN) {
  g <- D[, .(n_samples = .N,
             med_lfc  = median(direction * lfc),
             frac_up  = mean(direction * lfc > VG_MIN_LFC),
             med_ddet = median(direction * ddet)), by = .(gene, hierarchy_bin)]
  g[n_samples >= min_samples & frac_up >= VG_MIN_FRAC_SAMPLES &
    med_lfc >= VG_MIN_LFC & med_ddet >= VG_MIN_DDET][order(hierarchy_bin, -med_lfc)]
}
# contrast B gate, applied on top: the gene must move the same way against healthy donors too
gate_B <- function(g, direction = 1) {
  merge(g, B_STAT[, .(gene, hierarchy_bin, lfc_B, ddet_B)], by = c("gene", "hierarchy_bin"))[
    direction * lfc_B >= VG_MIN_LFC_B & direction * ddet_B >= 0]
}
UP <- gate_B(build(FC,  1),  1)[order(hierarchy_bin, -med_lfc)][, head(.SD, VG_TOP_N), by = hierarchy_bin]
DN <- gate_B(build(FC, -1), -1)[order(hierarchy_bin, -med_lfc)][, head(.SD, VG_TOP_N), by = hierarchy_bin]
cat(sprintf("\n[4] contrast A alone -> %d up genes; after the healthy-donor gate -> %d\n",
            nrow(build(FC, 1)), nrow(UP)))

SIG <- rbind(
  UP[, .(signature = paste0("malig_", hierarchy_bin), gene, weight = 1, bin = hierarchy_bin, dir = "up")],
  DN[, .(signature = paste0("norm_",  hierarchy_bin), gene, weight = 1, bin = hierarchy_bin, dir = "down")])
fwrite(SIG, OUT, sep = "\t")

cat("\n================ signature sizes ================\n")
print(dcast(SIG[, .N, by = .(bin, dir)], bin ~ dir, value.var = "N", fill = 0L))

cat("\n================ how bin-specific are they? (Jaccard between malignant sets) ================\n")
ml <- split(UP$gene, UP$hierarchy_bin)
if (length(ml) > 1) {
  J <- CJ(a = names(ml), b = names(ml))[a < b]
  J[, jaccard := mapply(function(x, y) round(length(intersect(ml[[x]], ml[[y]])) /
                                             length(union(ml[[x]], ml[[y]])), 3), a, b)]
  print(J[order(-jaccard)])
  cat("  High overlap is not automatically wrong -- blasts share a programme -- but a Jaccard near 1\n")
  cat("  would mean the per-bin split buys nothing and one signature would do.\n")
}

cat("\n================ LEAVE-ONE-SAMPLE-OUT STABILITY ================\n")
# If dropping one donor changes the signature substantially, the signature is that donor's.
STAB <- rbindlist(lapply(intersect(BINS, unique(UP$hierarchy_bin)), function(b) {
  full <- UP[hierarchy_bin == b]$gene
  ss <- unique(FC[hierarchy_bin == b]$smp)
  rbindlist(lapply(ss, function(drop) {
    # min_samples drops by one along with the data -- otherwise the filter tightens as the sample
    # shrinks and the test reports instability it created itself
    u <- gate_B(build(FC[hierarchy_bin == b & smp != drop], 1,
                      min_samples = max(2L, VG_MIN_SAMPLES_PER_BIN - 1L)), 1)[order(-med_lfc)][seq_len(min(VG_TOP_N, .N))]
    data.table(bin = b, dropped = drop,
               retained = round(length(intersect(full, u$gene)) / max(1, length(full)), 3))
  }))
}), fill = TRUE)
print(STAB[, .(n_loo = .N, min_retained = min(retained), median_retained = median(retained)), by = bin])
cat("  median_retained is the share of the full signature that survives dropping any one donor.\n")

qc <- merge(SIG[, .N, by = .(bin, dir)],
            STAB[, .(min_retained = min(retained), median_retained = median(retained)), by = bin],
            by = "bin", all.x = TRUE)
fwrite(qc, QC_CSV)
cat(sprintf("\n[done] %s\n[done] %s\n", OUT, QC_CSV))

cat("\n================ top malignant genes per bin ================\n")
for (b in BINS) cat(sprintf("  %-14s %s\n", b, paste(head(UP[hierarchy_bin == b]$gene, 12), collapse = ", ")))
