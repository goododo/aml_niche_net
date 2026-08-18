#!/usr/bin/env Rscript
# 80_vangalen_signature_train.R ----
# Derive per-lineage MALIGNANT-STATE signatures from van Galen 2019, as an axis that does not
# depend on copy number.
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
# THE CONTRAST IS LINEAGE-MATCHED, AND THAT IS THE WHOLE DESIGN.
#
#   van Galen names every malignant cell `X-like` and every normal cell `X`. So `HSC-like` vs `HSC`
#   inside one sample differs in malignancy and in nothing else -- not donor, not batch, not
#   platform, and crucially NOT LINEAGE.
#
#   The first version of this script paired malignant vs normal inside OUR projection bins instead,
#   and it was wrong. The projection places malignant HSC-like/Prog-like blasts into the T_NK (299
#   cells), B_Plasma (112) and Erythroid (1,483) bins, so those contrasts were comparing myeloid
#   blasts against T cells. The tell was in the output and I read it as success: the three
#   lineage-MISMATCHED bins each yielded a full 50 genes, while correctly-matched Mono_DC yielded
#   3. Signature size was tracking lineage mismatch, not signal.
#
# TWO CONTRASTS, AND A GENE MUST SURVIVE BOTH.
#
#   A  WITHIN SAMPLE   `X-like` vs `X`, same donor, same library, same platform, same lineage.
#                      But most AML samples have almost no NORMAL myeloid cells left, so at 10
#                      cells a side only 4-10 samples per type are testable -- low power.
#   B  VS HEALTHY      AML `X-like` vs `X` in the 4 healthy donors inside GSE116256, same dataset
#                      and protocol. Ample power; donor is NOT controlled.
#
# Each contrast fails in a way the other does not: A is underpowered, B is donor-confounded. A gene
# that clears both cannot be a donor effect (B's confounder would not survive A) and is unlikely to
# be a low-power fluke (A's weakness is covered by B).
#
# GENES ARE RANKED BY DETECTION-RATE DIFFERENCE, NOT FOLD CHANGE -- see VG_RANK_BY in the config.
# Ranking by mean log fold change selects for abundance, which fills the signature with GAPDH,
# ACTG1, ENO1 and ribosomal genes; those sit at the top of the within-cell ranking in EVERY cell,
# so a UCell score built on them saturates and cannot discriminate.
#
# VALIDATION IS HELD-OUT AUC, NOT GENE RETENTION. Retention under leave-one-out says the signature
# is reproducible; it does not say it works. For each donor we retrain WITHOUT that donor and ask
# whether the resulting signature separates that donor's malignant cells from their normal cells.
# A type that cannot clear VG_MIN_AUC does not get a signature.
#
# WHAT THIS SCRIPT DOES NOT DO: it does not decide whether the signature transfers ACROSS datasets
# and platforms. That is 81's job, settled on healthy donors from data this script never sees.
#
# OUTPUT: scripts/02_malignancy/vangalen_malignant_signatures.tsv  (signature, gene, weight, ...)
#         results/tables/02_malignancy/vangalen_signature_qc.csv
#
#   Rscript scripts/02_malignancy/80_vangalen_signature_train.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(UCell) })

RAW    <- file.path(LARGE1_DIR, "00_raw", "public", "GEO", "GSE116256_RAW")
DS     <- "GSE116256"
OUT    <- file.path(SCRIPTS_DIR, "02_malignancy", "vangalen_malignant_signatures.tsv")
QC_CSV <- file.path(DIR_MALIGNANCY, "vangalen_signature_qc.csv")
AUC_CSV<- file.path(DIR_MALIGNANCY, "vangalen_heldout_auc.csv")
.core  <- function(x) sub("^.*_", "", x)

## ---------------------------------------------------------------------------------------------
## labels: van Galen's own malignant call AND his own cell type
## ---------------------------------------------------------------------------------------------
A <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(p) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(p)))
  d <- fread(cmd = paste("zcat", shQuote(p)), colClasses = "character", showProgress = FALSE)
  if (!all(c("PredictionRefined", "CellType") %in% names(d))) return(NULL)
  data.table(smp = s, cell = d$Cell, lab = trimws(d$PredictionRefined), ctype = trimws(d$CellType))
}), fill = TRUE)[lab %in% c("malignant", "normal")]

# `X-like` is malignant, `X` is its normal counterpart; strip the suffix to get the pairing key.
A[, ctype_base := sub("-like$", "", ctype)]
# INTEGRITY OF THE PAIRING, asserted rather than assumed: every malignant cell must carry a `-like`
# type and every normal cell must not. If this ever fails the pairing key is meaningless.
stopifnot(A[lab == "malignant" & !grepl("-like$", ctype), .N] == 0,
          A[lab == "normal"    &  grepl("-like$", ctype), .N] == 0)

# our bin is carried along for reporting only -- it takes no part in the contrast
rec <- rbindlist(lapply(list.files(file.path(ANNO_RECONCILED_DIR, DS), full.names = TRUE), function(p)
  fread(p, select = c("cell", "hierarchy_bin"))[, smp := sub("__anno_percell\\.csv$", "", basename(p))][]),
  fill = TRUE)
A[, k := paste(smp, .core(cell))]; rec[, k := paste(smp, .core(cell))]
L <- merge(A[, .(k, smp, lab, ctype, ctype_base)], rec[, .(k, cell, hierarchy_bin)], by = "k")
message(sprintf("[0] %d labelled cells over %d samples, %d paired cell types",
                nrow(L), uniqueN(L$smp), uniqueN(L$ctype_base)))

CT <- dcast(L[, .N, by = .(smp, ctype_base, lab)], smp + ctype_base ~ lab, value.var = "N", fill = 0L)
CT <- CT[malignant >= VG_MIN_CELLS_PER_SIDE & normal >= VG_MIN_CELLS_PER_SIDE]
TYPES <- CT[, .N, by = ctype_base][N >= VG_MIN_SAMPLES_PER_TYPE][order(-N)]$ctype_base
message(sprintf("[1] testable (sample,celltype) pairs: %d | types with >= %d such samples: %s",
                nrow(CT), VG_MIN_SAMPLES_PER_TYPE, paste(TYPES, collapse = ", ")))
if (!length(TYPES)) stop("no cell type has enough paired samples")

R <- qc_rds_roster(on_extra = "ignore")[dataset == DS]
HEALTHY <- grep(VG_HEALTHY_REGEX, unique(L$smp), value = TRUE)
message(sprintf("[2] healthy donors inside %s for contrast B: %s", DS, paste(HEALTHY, collapse = ", ")))

load_sample <- function(s) {
  rp <- R[sample == s]$rds
  if (!length(rp)) return(NULL)
  SeuratObject::LayerData(NormalizeData(readRDS(rp), verbose = FALSE), assay = "RNA", layer = "data")
}
bin_stats <- function(dat, cells)
  list(mu = Matrix::rowMeans(dat[, cells, drop = FALSE]), det = Matrix::rowMeans(dat[, cells, drop = FALSE] > 0))

## ---------------------------------------------------------------------------------------------
## pass 1: per-(sample,celltype) profiles for both contrasts
## ---------------------------------------------------------------------------------------------
FC <- list(); MAL <- list()
for (s in unique(L$smp)) {
  dat <- tryCatch(load_sample(s), error = function(e) NULL)
  if (is.null(dat)) { message("  [skip] no QC object: ", s); next }
  ls <- L[smp == s][cell %in% colnames(dat)]
  for (b in intersect(TYPES, unique(ls$ctype_base))) {
    mc <- ls[ctype_base == b & lab == "malignant"]$cell
    nc <- ls[ctype_base == b & lab == "normal"]$cell
    if (length(mc) >= VG_MIN_CELLS_PER_SIDE && length(nc) >= VG_MIN_CELLS_PER_SIDE) {
      a <- bin_stats(dat, mc); n <- bin_stats(dat, nc)
      FC[[paste(s, b)]] <- data.table(gene = rownames(dat), smp = s, ctype_base = b,
                                      lfc = a$mu - n$mu, ddet = a$det - n$det)
    }
    if (s %in% HEALTHY && length(nc) >= VG_MIN_CELLS_PER_SIDE) {
      n <- bin_stats(dat, nc)
      MAL[[paste(s, b, "H")]] <- data.table(gene = rownames(dat), smp = s, ctype_base = b,
                                            grp = "healthy_normal", mu = n$mu, det = n$det)
    }
    if (!(s %in% HEALTHY) && length(mc) >= VG_MIN_CELLS_PER_SIDE) {
      a <- bin_stats(dat, mc)
      MAL[[paste(s, b, "M")]] <- data.table(gene = rownames(dat), smp = s, ctype_base = b,
                                            grp = "aml_malignant", mu = a$mu, det = a$det)
    }
  }
  rm(dat); gc(verbose = FALSE)
  message(sprintf("  [profile] %s", s))
}
FC  <- rbindlist(FC,  fill = TRUE)
MAL <- rbindlist(MAL, fill = TRUE)
stopifnot(nrow(FC) > 0, nrow(MAL) > 0)

# contrast B: sample means first, so a big sample cannot dominate by cell count alone.
# THIS IS A FUNCTION OF THE DONOR SET, NOT A CONSTANT. Computing it once over every donor and then
# reusing it inside the leave-one-out refits leaked the held-out donor's own malignant profile into
# the signature that was about to be scored on that donor. The leak is small -- one sample in a
# mean over ~10 -- which is exactly why it would have gone unnoticed and inflated every AUC.
build_B <- function(M) {
  b <- M[, .(mu = mean(mu), det = mean(det)), by = .(gene, ctype_base, grp)]
  b <- dcast(b, gene + ctype_base ~ grp, value.var = c("mu", "det"))
  if (!all(c("mu_aml_malignant", "mu_healthy_normal") %in% names(b))) return(b[0])
  b[, `:=`(lfc_B = mu_aml_malignant - mu_healthy_normal,
           ddet_B = det_aml_malignant - det_healthy_normal)][]
}
B_STAT <- build_B(MAL)
B_TYPES <- intersect(TYPES, unique(B_STAT[!is.na(lfc_B)]$ctype_base))
message(sprintf("[3] contrast B available for %d of %d types: %s",
                length(B_TYPES), length(TYPES), paste(B_TYPES, collapse = ", ")))

## ---------------------------------------------------------------------------------------------
## consistency filter -> signature
## ---------------------------------------------------------------------------------------------
# min_samples is a PARAMETER, not the config constant. The leave-one-out checks below re-run this
# on n-1 samples; with the constant hard-coded, every gene in a 4-sample type was filtered out for
# having only 3 samples and retention came back as exactly 0 -- a property of the test, not of the
# signature. Any stability check that gets stricter as it removes data measures itself.
build <- function(D, direction = 1, min_samples = VG_MIN_SAMPLES_PER_TYPE) {
  g <- D[, .(n_samples = .N,
             med_lfc  = median(direction * lfc),
             frac_up  = mean(direction * lfc > VG_MIN_LFC),
             med_ddet = median(direction * ddet)), by = .(gene, ctype_base)]
  g[n_samples >= min_samples & frac_up >= VG_MIN_FRAC_SAMPLES &
    med_lfc >= VG_MIN_LFC & med_ddet >= VG_MIN_DDET]
}
gate_B <- function(g, direction = 1, B = B_STAT) {
  if (!nrow(g) || !nrow(B)) return(g[0])
  merge(g, B[, .(gene, ctype_base, lfc_B, ddet_B)], by = c("gene", "ctype_base"))[
    direction * lfc_B >= VG_MIN_LFC_B & direction * ddet_B >= 0]
}
rank_key <- function(g) if (VG_RANK_BY == "ddet") -g$med_ddet else -g$med_lfc
top_n <- function(g, n = VG_TOP_N) { if (!nrow(g)) return(g); g[order(ctype_base, rank_key(g))][, head(.SD, n), by = ctype_base] }

# one call, reused by the full fit and by every leave-one-out refit, so the test cannot drift
# from the thing it is testing
fit <- function(D, direction = 1, min_samples = VG_MIN_SAMPLES_PER_TYPE, n = VG_TOP_N, B = B_STAT)
  top_n(gate_B(build(D, direction, min_samples), direction, B), n)

UP <- fit(FC,  1); DN <- fit(FC, -1)
cat(sprintf("\n[4] contrast A alone -> %d up genes; after the healthy-donor gate -> %d\n",
            nrow(build(FC, 1)), nrow(UP)))

## ---------------------------------------------------------------------------------------------
## held-out AUC -- the gate that decides whether a signature is used at all
## ---------------------------------------------------------------------------------------------
# Mann-Whitney AUC. y is the malignant indicator.
auc_of <- function(score, y) {
  n1 <- sum(y); n0 <- sum(!y)
  if (n1 < 3 || n0 < 3 || !length(unique(score))) return(NA_real_)
  (sum(rank(score)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
# retrain WITHOUT the donor, then score that donor's own cells -- the signature never saw them
LOO_SIG <- rbindlist(lapply(TYPES, function(b) rbindlist(lapply(CT[ctype_base == b]$smp, function(drop) {
  # drop the donor from contrast B as well as contrast A -- otherwise the "held-out" donor is
  # still inside the healthy-donor gate that selected the genes
  u <- fit(FC[ctype_base == b & smp != drop], 1, max(2L, VG_MIN_SAMPLES_PER_TYPE - 1L),
           B = build_B(MAL[smp != drop]))
  if (!nrow(u)) return(NULL)
  data.table(ctype_base = b, dropped = drop, gene = u$gene)
}), fill = TRUE)), fill = TRUE)

AUC <- list()
for (s in intersect(unique(CT$smp), unique(LOO_SIG$dropped))) {
  need <- LOO_SIG[dropped == s]
  dat <- tryCatch(load_sample(s), error = function(e) NULL)
  if (is.null(dat)) next
  ls <- L[smp == s][cell %in% colnames(dat)]
  for (b in unique(need$ctype_base)) {
    g <- intersect(unique(need[ctype_base == b]$gene), rownames(dat))
    cells <- ls[ctype_base == b & lab %in% c("malignant", "normal")]
    if (length(g) < 3 || nrow(cells) < 20) next
    sc <- UCell::ScoreSignatures_UCell(dat[, cells$cell, drop = FALSE], features = list(sig = g),
                                       maxRank = CELLSTATE_MAXRANK, ncores = 1, name = "")
    AUC[[paste(s, b)]] <- data.table(smp = s, ctype_base = b, n_genes = length(g),
                                     n_mal = sum(cells$lab == "malignant"), n_norm = sum(cells$lab == "normal"),
                                     auc = auc_of(as.numeric(sc[, "sig"]), cells$lab == "malignant"))
  }
  rm(dat); gc(verbose = FALSE)
  message(sprintf("  [heldout] %s", s))
}
AUC <- rbindlist(AUC, fill = TRUE)
fwrite(AUC, AUC_CSV)

cat("\n================ HELD-OUT AUC (trained without the donor, scored on that donor) ================\n")
AS <- AUC[!is.na(auc), .(n_donors = .N, median_auc = round(median(auc), 3),
                         min_auc = round(min(auc), 3), n_genes = round(median(n_genes))), by = ctype_base][order(-median_auc)]
print(AS)
cat(sprintf("  A signature is kept only if median held-out AUC >= %.2f. 0.5 is coin-flip.\n", VG_MIN_AUC))
KEEP <- AS[median_auc >= VG_MIN_AUC]$ctype_base
DROP <- setdiff(TYPES, KEEP)
if (length(DROP)) cat(sprintf("  [dropped] %s -- below the AUC gate, not used regardless of stability\n",
                              paste(DROP, collapse = ", ")))
UP <- UP[ctype_base %in% KEEP]; DN <- DN[ctype_base %in% KEEP]
if (!nrow(UP)) stop("no cell type cleared the held-out AUC gate -- there is no usable signature")

## ---------------------------------------------------------------------------------------------
## pan-blast signature: consistent across lineages, so it is about malignancy and not about lineage
## ---------------------------------------------------------------------------------------------
# a threshold above the number of surviving lineages silently returns nothing and reads as
# "no consistent genes" when it actually means "the question was not asked"
PAN_MIN <- min(VG_PANBLAST_MIN_TYPES, length(KEEP))
PAN <- if (length(KEEP) < 2) UP[0][, .(gene, n_types = integer())] else
  UP[, .(n_types = uniqueN(ctype_base)), by = gene][n_types >= PAN_MIN][order(-n_types)]
cat(sprintf("\n================ pan-blast: genes up in >= %d of %d kept lineages ================\n",
            PAN_MIN, length(KEEP)))
if (length(KEEP) < VG_PANBLAST_MIN_TYPES)
  cat(sprintf("  [note] only %d lineages cleared the AUC gate, below the configured %d -- a pan-blast\n         signature built from so few lineages cannot show that a gene is lineage-independent.\n",
              length(KEEP), VG_PANBLAST_MIN_TYPES))
cat(sprintf("  %d genes: %s\n", nrow(PAN), paste(head(PAN$gene, 40), collapse = ", ")))
cat("  This is what gets applied to cells whose bin has no matched counterpart -- which matters,\n")
cat("  because the projection does put blasts in the T_NK / B_Plasma / Erythroid bins.\n")

SIG <- rbind(
  UP[, .(signature = paste0("malig_", ctype_base), gene, weight = 1, ctype_base, dir = "up",
         bin = unname(VG_TYPE_TO_BIN[ctype_base]))],
  DN[, .(signature = paste0("norm_",  ctype_base), gene, weight = 1, ctype_base, dir = "down",
         bin = unname(VG_TYPE_TO_BIN[ctype_base]))],
  if (nrow(PAN)) PAN[, .(signature = "malig_panblast", gene, weight = 1, ctype_base = "ANY",
                         dir = "up", bin = "ANY")] else NULL, fill = TRUE)
fwrite(SIG, OUT, sep = "\t")

cat("\n================ signature sizes ================\n")
print(dcast(SIG[, .N, by = .(signature, dir)], signature ~ dir, value.var = "N", fill = 0L))

cat("\n================ composition check: is this just housekeeping? ================\n")
# The failure mode that killed the previous version. Ribosomal + a short list of the genes that
# topped every abundance-ranked signature. A high share here means the score will saturate.
HK <- c("GAPDH","ACTB","ACTG1","EEF1A1","EIF1","EIF4A1","ENO1","PKM","TPT1","OAZ1","SERF2","PSMA7",
        "COX4I1","HINT1","YBX1","RAN","NACA","BTF3","PFN1","MYL6","H3F3A","HSP90AB1","HSP90AA1")
comp <- SIG[dir == "up", .(n = .N,
                           ribo = sum(grepl("^RP[SL]", gene)),
                           hk   = sum(gene %in% HK)), by = signature]
comp[, pct_hk_ribo := round(100 * (ribo + hk) / n, 1)]
print(comp[order(-pct_hk_ribo)])

cat("\n================ how lineage-specific are they? (Jaccard between malignant sets) ================\n")
ml <- split(UP$gene, UP$ctype_base)
if (length(ml) > 1) {
  J <- CJ(a = names(ml), b = names(ml))[a < b]
  J[, jaccard := mapply(function(x, y) round(length(intersect(ml[[x]], ml[[y]])) /
                                             length(union(ml[[x]], ml[[y]])), 3), a, b)]
  print(J[order(-jaccard)])
}

cat("\n================ LEAVE-ONE-SAMPLE-OUT GENE RETENTION ================\n")
STAB <- rbindlist(lapply(KEEP, function(b) {
  full <- UP[ctype_base == b]$gene
  rbindlist(lapply(CT[ctype_base == b]$smp, function(drop) {
    u <- LOO_SIG[ctype_base == b & dropped == drop]$gene
    data.table(ctype_base = b, dropped = drop,
               retained = round(length(intersect(full, u)) / max(1, length(full)), 3))
  }))
}), fill = TRUE)
print(STAB[, .(n_loo = .N, min_retained = min(retained), median_retained = median(retained)), by = ctype_base])
cat("  Retention says the signature is REPRODUCIBLE; the AUC table above says whether it WORKS.\n")

qc <- merge(merge(SIG[dir == "up", .(n_genes = .N), by = .(ctype_base)],
                  STAB[, .(min_retained = min(retained), median_retained = median(retained)), by = ctype_base],
                  by = "ctype_base", all.x = TRUE),
            AS, by = "ctype_base", all.x = TRUE)
fwrite(qc, QC_CSV)
cat(sprintf("\n[done] %s\n[done] %s\n[done] %s\n", OUT, QC_CSV, AUC_CSV))

cat("\n================ top malignant genes per lineage ================\n")
for (b in KEEP) cat(sprintf("  %-10s %s\n", b, paste(head(UP[ctype_base == b]$gene, 14), collapse = ", ")))
