#!/usr/bin/env Rscript
# 06_mp_srrs.R ----
# INPUT  : CNMF_RES_DIR/malignant/nmf_res.rds        (606 = 101 samples x k=4:9; each list(w = 2000 x k))
#          CNMF_TAB_DIR/malignant_metaprograms.tsv   (MP -> genes; comma-delimited despite the extension)
#          CNMF_TAB_DIR/malignant_mp_annotation.csv  (per MP: class, sampleCoverage, silhouette)
#          DIR_PREPROCESS/02_sample_split.csv        (dataset -> Discovery / Validation arm)
# OUTPUT : CNMF_TAB_DIR/mp_srrs.csv           one row per MP: Global SRRS, its null, the stability grid
#          CNMF_TAB_DIR/mp_srrs_by_study.csv  one row per MP x study: max-Jaccard, and the arm the study is in
#          CNMF_TAB_DIR/mp_srrs_grid.csv      one row per MP x signature length x tau (the sensitivity sweep)
# WHAT   : Computes the Global part of the two-level reproducibility score -- for each meta-program, the
#          fraction of INDEPENDENT STUDIES that contain a per-sample NMF program similar to it.
#
# WHY THIS EXISTS. What currently stands in for reproducibility is GeneNMF's sampleCoverage, and it
# counts SAMPLES. One large cohort can therefore carry a meta-program on its own: GSE185381 alone
# contributes 52 of 212 samples. Counting independent STUDIES is the protection the two-level design
# was specified for, and it is the missing half of P8 in scripts/DECISIONS_pending.md.
#
# The two scores are not variants of each other. Measured here, Spearman(sampleCoverage, SRRS) is
# ~0.0 -- sampleCoverage ranks MP1 first (a one-gene filler cluster the curation labels "drop")
# while SRRS ranks it last. Anything that has been justified by sampleCoverage needs re-reading.
#
# WHAT THIS SCRIPT DOES NOT DO. Conditional SRRS (per genetic subtype) is not implemented and should
# not be: subtype, karyotype and mutations are 0 of 244 populated in the curated metadata, and
# subtype_dataset_level is constant within every dataset, so stratifying on it IS stratifying on
# study and gives denominators of 1, 1, 1, 2 and 4. That is recorded as a refusal, not a deferral.
#
# READING THE OUTPUT. tau is the Jaccard a per-study program must reach to count as a match. The
# retention bar (SRRS >= 0.6) is inherited from the blueprint, but with only ~9 studies it can take
# 4 distinct values in the region that matters, so read the GRID -- how many of the 16 (length, tau)
# cells a program passes -- rather than the single number.
#
# Usage : Rscript scripts/04_cnmf/06_mp_srrs.R [--tau 0.2] [--n_null 2000] [--compartment malignant]

suppressPackageStartupMessages({
  library(data.table); library(optparse); library(here); library(GeneNMF); library(Matrix)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--tau",         type = "double",    default = 0.2,
              help = "Jaccard at which a per-study program counts as reproducing an MP"),
  make_option("--n_null",      type = "integer",   default = 2000L,
              help = "size- and frequency-matched null draws; 0 skips the null"),
  make_option("--compartment", type = "character", default = "malignant",
              help = "malignant | normal -- which nmf_res.rds to read"),
  make_option("--top_n",       type = "integer",   default = 50L,
              help = "genes taken from each per-sample NMF program")
)))
set.seed(SEED)

## -- Step 1. the meta-programs under test ----
mp_tsv <- file.path(CNMF_TAB_DIR, sprintf("%s_metaprograms.tsv", opt$compartment))
stopifnot(file.exists(mp_tsv))
sig <- fread(mp_tsv, sep = ",")            # .tsv by name, comma-delimited by content
stopifnot(all(c("MP", "gene") %in% names(sig)))
mp_genes <- split(sig$gene, sig$MP)
message(sprintf("[1] %d meta-programs, %d-%d genes each",
                length(mp_genes), min(lengths(mp_genes)), max(lengths(mp_genes))))

## -- Step 2. the per-sample programs, and which STUDY each came from ----
nmf_f <- file.path(CNMF_RES_DIR, opt$compartment, "nmf_res.rds")
stopifnot(file.exists(nmf_f))
nmf <- readRDS(nmf_f)
# Names are "<dataset>..<sample>.k<K>". The study is everything before the first "..", which is the
# unit SRRS counts -- not the sample, which is what sampleCoverage counts.
# getNMFgenes returns ONE ENTRY PER NMF FACTOR, not per fitted object: 3924 programs from the 606
# (sample x k) fits, named "<dataset>..<sample>.k<K>.<factor>". It also drops factors whose weighted
# loading selects no gene (3924 vs a nominal 3939), so the study label must come from the returned
# names -- deriving it from names(nmf) would misalign every program by a shifting offset.
# AND IT RETURNS WEIGHTS, NOT NAMES. Each element is a numeric vector of loadings whose NAMES are the
# genes -- so `prog[[1]]` is c(0.0795, 0.0455, ...) and the gene symbols are in names(prog[[1]]).
# Taking the values as a gene set gives an overlap of exactly 0 against every MP, which reads as
# "nothing reproduces" rather than as an error. It has to be lapply(., names).
prog   <- lapply(getNMFgenes(nmf, max.genes = opt$top_n), names)
nm     <- names(prog)
study  <- sub("\\.\\..*$", "", nm)
smp    <- sub("^[^.]*\\.\\.", "", sub("\\.k[0-9]+\\.[0-9]+$", "", nm))
stopifnot(length(prog) == length(nm), all(nzchar(study)), uniqueN(study) > 1)
message(sprintf("[2] %d per-sample programs from %d samples across %d studies: %s",
                length(prog), uniqueN(smp), uniqueN(study), paste(sort(unique(study)), collapse = ", ")))

## -- Step 3. max-Jaccard of each MP against each study ----
# Jaccard on gene sets. Done with a sparse indicator matrix rather than a double loop because the
# comparison is 10 MPs x 606 programs and the intersection is a single crossprod.
universe <- sort(unique(c(unlist(mp_genes), unlist(prog))))
.ind <- function(sets) {
  i <- unlist(lapply(sets, function(g) match(intersect(g, universe), universe)))
  j <- rep(seq_along(sets), vapply(sets, function(g) length(intersect(g, universe)), integer(1)))
  sparseMatrix(i = i, j = j, dims = c(length(universe), length(sets)), x = 1)
}
A <- .ind(mp_genes); B <- .ind(prog)
inter <- as.matrix(Matrix::crossprod(A, B))                       # |MP ∩ program|
sizeA <- Matrix::colSums(A); sizeB <- Matrix::colSums(B)
jac   <- inter / (outer(sizeA, sizeB, "+") - inter)               # |∩| / |∪|
dimnames(jac) <- list(names(mp_genes), nm)

by_study <- rbindlist(lapply(names(mp_genes), function(mp)
  data.table(MP = mp, study = unique(study),
             max_jaccard = vapply(unique(study), function(s) max(jac[mp, study == s]), numeric(1)))))
by_study[, reproduces := max_jaccard >= opt$tau]

# THE ARM EACH STUDY SITS IN. The MP signatures were discovered on all studies, so a meta-program can
# have been defined partly on Validation data -- which means the mp_* node features screened on
# Discovery are not a clean Discovery-arm feature set. Making the arm a column puts that on the audit
# trail instead of leaving it to be discovered later.
sp_f <- file.path(DIR_PREPROCESS, "02_sample_split.csv")
if (file.exists(sp_f)) {
  arm <- unique(fread(sp_f)[, .(study = dataset, split_arm = split)])
  arm <- arm[, .(split_arm = paste(sort(unique(split_arm)), collapse = "/")), by = study]
  by_study <- merge(by_study, arm, by = "study", all.x = TRUE)
  by_study[is.na(split_arm), split_arm := "unassigned"]
} else by_study[, split_arm := NA_character_]

SR <- by_study[, .(n_studies = .N, n_reproducing = sum(reproduces),
                   SRRS = sum(reproduces) / .N,
                   SRRS_discovery  = sum(reproduces & split_arm %like% "Discovery")  /
                                     max(sum(split_arm %like% "Discovery"), 1L),
                   SRRS_validation = sum(reproduces & split_arm %like% "Validation") /
                                     max(sum(split_arm %like% "Validation"), 1L)), by = MP]

## -- Step 4. is tau above chance? size- and frequency-matched null ----
# A large MP matches more by chance than a small one, and a gene present in many programs matches more
# than a rare one. The null therefore draws sets of the SAME SIZE with gene probabilities proportional
# to how often each gene appears across the real programs. Without both matches the null is too easy
# and every MP looks reproducible.
if (opt$n_null > 0L) {
  freq <- Matrix::rowSums(B); freq <- freq / sum(freq)
  null_srrs <- vapply(names(mp_genes), function(mp) {
    n <- length(intersect(mp_genes[[mp]], universe))
    draws <- vapply(seq_len(opt$n_null), function(b) {
      g <- sample.int(length(universe), n, prob = freq)
      v <- rep(0, length(universe)); v[g] <- 1
      it <- as.numeric(v %*% B)
      jj <- it / (n + sizeB - it)
      mean(vapply(unique(study), function(s) max(jj[study == s]) >= opt$tau, logical(1)))
    }, numeric(1))
    quantile(draws, 0.95)
  }, numeric(1))
  SR[, SRRS_null_q95 := null_srrs[MP]]
  SR[, above_null := SRRS > SRRS_null_q95]
}

## -- Step 5. the sensitivity grid -- tau AND signature length ----
# tau is not the only free parameter and it is not the one the score is most sensitive to. Sweeping
# both is the honest presentation: a program that passes in one corner of the grid and fails in the
# rest is not reproducible, whatever its headline number says.
TAUS <- c(0.15, 0.20, 0.25, 0.30)
LENS <- c(30L, 50L, 75L, 100L)
grid <- rbindlist(lapply(LENS, function(L) {
  pL <- lapply(getNMFgenes(nmf, max.genes = L), names)   # names, not values -- see step 2
  BL <- .ind(pL); szL <- Matrix::colSums(BL)
  iL <- as.matrix(Matrix::crossprod(A, BL))
  jL <- iL / (outer(sizeA, szL, "+") - iL)
  dimnames(jL) <- list(names(mp_genes), names(pL))    # indexed by MP name below; without this it fails
  rbindlist(lapply(TAUS, function(t)
    data.table(MP = names(mp_genes), sig_len = L, tau = t,
               SRRS = vapply(names(mp_genes), function(mp)
                 mean(vapply(unique(study), function(s) max(jL[mp, study == s]) >= t, logical(1))),
                 numeric(1)))))
}))
grid[, passes := SRRS >= 0.6]
SR <- merge(SR, grid[, .(cells_ge_0.6 = sum(passes), n_cells = .N), by = MP], by = "MP", all.x = TRUE)

## -- Step 6. next to what it replaces ----
ann_f <- file.path(CNMF_TAB_DIR, sprintf("%s_mp_annotation.csv", opt$compartment))
if (file.exists(ann_f)) {
  ann <- fread(ann_f); setnames(ann, "malignant_MP", "MP", skip_absent = TRUE)
  SR <- merge(SR, ann[, .(MP, sampleCoverage, silhouette, code_class = class)], by = "MP", all.x = TRUE)
}
lab_f <- file.path(CNMF_SCRIPT_DIR, "mp_labels.tsv")
if (file.exists(lab_f))
  SR <- merge(SR, fread(lab_f)[, .(MP, hand_confidence = confidence)], by = "MP", all.x = TRUE)

setorder(SR, -SRRS)
fwrite_safe(SR,       file.path(CNMF_TAB_DIR, "mp_srrs.csv"))
fwrite_safe(by_study, file.path(CNMF_TAB_DIR, "mp_srrs_by_study.csv"))
fwrite_safe(grid,     file.path(CNMF_TAB_DIR, "mp_srrs_grid.csv"))

message(sprintf("\n[3] GLOBAL SRRS at tau = %.2f, top %d genes, over %d studies",
                opt$tau, opt$top_n, uniqueN(study)))
cols <- intersect(c("MP","SRRS","n_reproducing","n_studies","SRRS_null_q95","above_null",
                    "cells_ge_0.6","n_cells","SRRS_discovery","SRRS_validation",
                    "sampleCoverage","hand_confidence"), names(SR))
print(SR[, ..cols])

# THE COMPARISON THAT MATTERS. If these two ranked programs the same way, the switch would be
# cosmetic. Reporting the correlation makes it impossible to assume that without looking.
if ("sampleCoverage" %in% names(SR)) {
  rho <- suppressWarnings(cor(SR$SRRS, SR$sampleCoverage, method = "spearman", use = "complete.obs"))
  # A constant SRRS column makes the correlation NA. That is not a weak relationship, it is a
  # degenerate score, and it must say so rather than fall through an if() on NA.
  if (!is.finite(rho))
    message("\n[4] Spearman(SRRS, sampleCoverage) is undefined -- SRRS has no variance across the ",
            nrow(SR), " meta-programs. Check the gene sets before reading anything else here.")
  else message(sprintf("\n[4] Spearman(SRRS, sampleCoverage) = %+.3f over %d meta-programs.", rho, nrow(SR)))
  if (is.finite(rho) && abs(rho) < 0.3)
    message("    They carry essentially no shared information: every claim resting on sampleCoverage\n",
            "    needs re-reading against SRRS, not treating as already supported.")
}
if ("hand_confidence" %in% names(SR)) {
  rt <- SR[hand_confidence == "robust_tumor"]
  if (nrow(rt))
    message(sprintf("\n[5] The %d hand-labelled robust_tumor MPs pass %s of the %d grid cells (%s).",
                    nrow(rt), paste(rt$cells_ge_0.6, collapse = "/"), rt$n_cells[1],
                    paste(rt$MP, collapse = ", ")))
}
if ("split_arm" %in% names(by_study)) {
  nD <- uniqueN(by_study[split_arm %like% "Discovery"]$study)
  nV <- uniqueN(by_study[split_arm %like% "Validation"]$study)
  message(sprintf("\n[6] LEAK CHECK: the MP signatures were discovered on %d studies, of which %d are\n",
                  uniqueN(study), nV),
          "    Validation-arm. The mp_* node features screened on Discovery are therefore NOT a clean\n",
          "    Discovery-arm feature set. Report it; do not treat the mp family as held out.\n",
          sprintf("    (Discovery studies: %d, Validation: %d)", nD, nV))
}
message("\n[done] wrote mp_srrs.csv, mp_srrs_by_study.csv, mp_srrs_grid.csv to ", CNMF_TAB_DIR)
