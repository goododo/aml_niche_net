#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/<ds>/ws_pb_<ds>.rds for GSE185381, GSE120221, GSE227690, GSE276360 (26: raw T sums per
#            sample, fine-subtype sums, per-cell kept table); results/tables/ws/ws_s4_pseudobulk_meta_<ds>.csv
#          results/tables/ws/ws_s5_scores.csv (the 25 eligible GSE185381 T samples and the S5 covariates), ws_s5_imprint_genes.csv,
#            ws_s8_t_imprint_ambient_prone.csv, ws_s5_de_T.csv, ws_s8_scores_disc.csv, decontx/GSE185381/ws_s1_ambient_blacklist.csv
#          /LARGE1/.../00_raw/metadata/<GSE>/donor_ages_<GSE>_v1.0.csv and meta_<GSE>_v1.0.csv (donor age and sex; K7)
#          /LARGE1/.../00_raw/metadata/GSE185381/meta_GSE185381_v2.3.csv (ages of the 25), BMM_REF_RDS meta_data (K7, Oetjen)
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE227690/*.rds (author CellType labels, K1 report only)
# OUTPUT : results/tables/ws/ws_s9_known_answers.csv, ws_s9_checks.csv, ws_s9_yardstick.csv  (phase 1, always)
#          results/tables/ws/ws_s9_tests.csv, ws_s9_scores_<cohort>.csv, ws_s9_direction.csv, ws_s9_direction_positions.csv, ws_s9_age_genes.csv,
#            ws_s9_robust_imprint.csv, /LARGE1/.../ws_pseudobulk/ws_s9_null_draws.rds; GATES.md block "A1-S9"  (phase 2)
# what it does: PREREGISTRATION_S9_age_direction.md (tag prereg-s9-v1, sections 2-5 and 10), nothing else. Phase 1 runs the
#   known answers in the order of section 5 (K7, K3, K0, K4, r, K6, K2, K5, K1) without reading the age relation of any new
#   cohort; if one of them fails the script writes its tables and stops. Phase 2 reads the age relation: A1 (f_age in the joint
#   space, deflated by r), A2b (wording modifier), leave-one-donor-out, the CI wording, A3 per cohort, then the reported
#   analyses (W space, sensitivities, A2a, section 3.4). The scoring machinery is S8's (32), copied verbatim.
#   usage: Rscript scripts/ws/33_s9_age_direction.R [--smoke]
#   --smoke: 50 random sets, 100 bootstrap draws, 5 cross-fit folds, donor ages of the new cohorts PERMUTED (the real age
#            relation is not read), outputs suffixed _smoke, GATES.md not written.
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(SeuratObject); library(here) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R")); set.seed(SEED)
options(error = function() { calls <- sys.calls(); n <- length(calls); message("--- TRACEBACK ---")
  for (i in seq(max(1, n - 15), n)) message(i, ": ", substr(paste(deparse(calls[[i]]), collapse = " "), 1, 220)); q(status = 1) })
TABWS <- file.path(TAB_DIR, "ws"); GATES <- here::here("scripts", "ws", "GATES.md"); METADIR <- file.path(RAW_DIR, "metadata")
SFX <- if (SMOKE) "_smoke" else ""; out <- function(x) file.path(TABWS, paste0("ws_s9_", x, SFX, ".csv"))
N_RAND <- if (SMOKE) 50L else 1000L; N_BOOT <- if (SMOKE) 100L else 1000L; N_RM <- if (SMOKE) 50L else 1000L; N_FOLD <- if (SMOKE) 5L else 25L
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list(); chk <- function(name, expected, got) { got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-82s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **")) }
KA <- list(); ka <- function(id, statistic, value, threshold, pass, detail = "") {
  KA[[length(KA) + 1]] <<- data.table(id = id, statistic = statistic, value = as.character(value), threshold = threshold, pass = pass, detail = detail)
  message(sprintf("  [%s] %-60s = %-12s (%s) %s", id, statistic, value, threshold, if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "** FAIL **" else "reported")) }
f3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

## -- 0. registered constants (PREREGISTRATION_S9 sections 2.1-2.3, 3, 4) -----------------------------------------------
GAP_IN_REG <- 1.970; GAP_TOL <- 0.01            # section 2.1: the yardstick written before any new-cohort score; K6 tolerance
DELTA_AGE <- 31.5; AGE_REF_HEALTHY <- 32.5       # median ages of the 15 eligible AML (64) and 10 healthy (32.5)
YOUNG_MAX <- 53; OLD_MIN <- 60                   # young = <= 53, old = >= 60 (section 2.3)
LINE <- 1 / 3                                    # the GATE_A1 (4) line
MIN_T <- 20L                                     # kept T cells per DONOR (section 2.3)
K0_MIN_DONORS <- 15L; K0_MIN_OLD <- 4L; K4_MIN_COV <- 0.80; R_MAX <- 2
Y_GENES <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y")   # section 10.2 (GSE227690, no XIST)
NAIVE_GENES <- c("CCR7", "SELL", "LEF1", "TCF7", "IL7R"); CYTO_GENES <- c("GZMH", "GZMB", "PRF1", "NKG7", "KLRG1", "B3GAT1")   # K1b
COHORTS <- c("GSE120221", "GSE227690", "GSE276360"); DECISIVE <- "GSE120221"
T_FINE <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE)[group4 == "T", bmm_fine]
stopifnot(length(T_FINE) == 9L, "CD8 Naive" %in% T_FINE, "CD4 Naive" %in% T_FINE)

## -- 1. gene sets, blacklist, S8 scoring helpers (32, verbatim) ----------------------------------------------------------
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))
GS <- list(primary = IMP[!gene %in% AMB$gene], full = IMP)
chk("primary set is 370 genes (156 up / 214 down)", "370/156/214", GS$primary[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))])
chk("full set is 402 genes (171 up / 231 down)", "402/171/231", GS$full[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))])
UP <- GS$primary[direction == "up", gene]; DN <- GS$primary[direction == "down", gene]; PRIM <- c(UP, DN)
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
black_for <- function(G) unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene], "XIST", GO[chr == "chrY", gene],
                                  grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
logz <- function(M) {   # 32: blacklist removed; CPM > 1 in at least half of the set's samples; TMM; log CPM prior 1; z across the set
  M <- M[!rownames(M) %in% black_for(rownames(M)), , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L <- cpm(y, log = TRUE, prior.count = 1); Z <- t(scale(t(L)))
  U <- data.table(gene = rownames(L), mean_lcpm = rowMeans(L), det = rowMeans(M[rownames(L), , drop = FALSE] > 0))
  U[, bin := paste(cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE), cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf)))]
  list(Z = Z, U = U, L = L) }
score <- function(Z, up, dn) { up <- intersect(up, rownames(Z)); dn <- intersect(dn, rownames(Z))
  (if (length(up)) colMeans(Z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(Z[dn, , drop = FALSE]) else 0) }
matched <- function(U, genes) { tb <- table(U[gene %in% genes, bin])
  unlist(lapply(names(tb), function(b) { pool <- U[bin == b & !gene %in% genes, gene]; if (!length(pool)) pool <- U[!gene %in% genes, gene]; sample(pool, tb[[b]], replace = length(pool) < tb[[b]]) })) }
split_random <- function(U, up, dn) { g <- matched(U, c(up, dn)); list(up = head(g, length(intersect(up, U$gene))), dn = tail(g, length(intersect(dn, U$gene)))) }

## -- 2. GSE185381: the 25 eligible T pseudobulks, ages, the S5 space, the S5 refit (F for lambda) ---------------------------
ts("[2] GSE185381")
PBG <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
S5 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]; stopifnot(nrow(S5) == 25L, all(S5$sample %in% colnames(PBG$raw$T)))
MG <- fread(file.path(METADIR, "GSE185381", "meta_GSE185381_v2.3.csv"), colClasses = "character")[, sid := sub("^GSE185381__", "", sample_id)]
S5[, age := as.numeric(MG$age_years[match(sample, MG$sid)])]; S5[, status := factor(status, c("healthy", "AML"))]
chk("Delta_age = median age of the 15 AML minus the 10 healthy = 31.5", DELTA_AGE, S5[status == "AML", median(age)] - S5[status == "healthy", median(age)])
chk("healthy median age 32.5 (five at 19-26, five at 39-53)", AGE_REF_HEALTHY, S5[status == "healthy", median(age)])
G_RAW <- PBG$raw$T[, S5$sample]
S5Z <- logz(G_RAW)                                                           # the S5 4.1 space on the 25 (the yardstick space)
sc_in <- score(S5Z$Z, UP, DN); st <- S5$status[match(names(sc_in), S5$sample)]
GAP_IN <- mean(sc_in[st == "AML"]) - mean(sc_in[st == "healthy"]); SD_H <- sd(sc_in[st == "healthy"])
ts(sprintf("    yardstick: gap_in %.4f = %.2f healthy-SD (AML %.3f, healthy %.3f, SD_h %.3f, universe %d genes)", GAP_IN, GAP_IN / SD_H, mean(sc_in[st == "AML"]), mean(sc_in[st == "healthy"]), SD_H, nrow(S5Z$Z)))
# the S5 differential model (27's fit_de, edgeR branch, with F kept) -- used for the cross-fit and for lambda
MAIN <- c("status", "hb", "gran", "sex", "subtype_share"); FDR_CUT <- 0.10; LFC_CUT <- 0.5
fit_de <- function(counts, d, terms, coef) {
  keep_g <- !(rownames(counts) %in% black_for(rownames(counts)))
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = d$status)
  keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = as.data.frame(d)); stopifnot(qr(X)$rank == ncol(X))
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
  r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, F, PValue, FDR)]
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]; r }
D5 <- S5[, .(sample, status, hb, gran, sex = factor(sex), subtype_share)]
DE5 <- fread(file.path(TABWS, "ws_s5_de_T.csv"))
ts("    refitting the S5 model on the 25 samples (for the per-gene F)")
REFIT <- fit_de(G_RAW, D5, MAIN, "statusAML")
m5 <- merge(REFIT[, .(gene, logFC_refit = logFC, F)], DE5[, .(gene, logFC)], by = "gene")
chk("S5 refit reproduces ws_s5_de_T log2FC (max |diff| < 1e-6 over the tested genes)", TRUE, nrow(m5) == nrow(DE5) && max(abs(m5$logFC_refit - m5$logFC)) < 1e-6)
REFIT[, SE := fifelse(F > 1e-12, abs(logFC) / sqrt(F), NA_real_)]

## -- 3. the new cohorts: donor pooling, ages, sex, K7, K3, K0 ----------------------------------------------------------------
ts("[3] new cohorts")
age_group_of <- function(a) ifelse(a <= YOUNG_MAX, "young", ifelse(a >= OLD_MIN, "old", "middle"))
load_cohort <- function(ds) {
  PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", ds, paste0("ws_pb_", ds, ".rds")))
  META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", ds, ".csv")))
  don <- fread(file.path(METADIR, ds, paste0("donor_ages_", ds, "_v1.0.csv")), colClasses = "character")
  meta <- fread(file.path(METADIR, ds, paste0("meta_", ds, "_v1.0.csv")), colClasses = "character")
  lib <- unique(META[, .(sample, n_kept)]); lib[, gsm := sub("_.*$", "", sample)]
  lib[, `:=`(donor = don$donor[match(gsm, don$gsm)], age_lib = as.numeric(don$age[match(gsm, don$gsm)]), sex_lib = don$sex[match(gsm, don$gsm)])]
  stopifnot(!anyNA(lib$donor), all(colnames(PB$raw$T) %in% lib$sample))
  # kept T cells and T raw sums per library, then pooled per donor; the donor's age is that of its first-listed (lowest GSM) library
  nT <- META[group4 == "T", .(sample, n_T = n)]; lib <- merge(lib, nT, by = "sample", all.x = TRUE); lib[is.na(n_T), n_T := 0L]
  donors <- lib[order(gsm), .(n_lib = .N, libraries = paste(sample, collapse = ";"), age = age_lib[1], sex = sex_lib[1], n_T = sum(n_T), n_kept = sum(n_kept)), by = donor]
  donors[, `:=`(eligible = n_T >= MIN_T, grp = age_group_of(age), excluded_K3 = FALSE)]
  T_don <- sapply(donors$donor, function(d) rowSums(PB$raw$T[, lib[donor == d, sample], drop = FALSE])); colnames(T_don) <- donors$donor
  list(ds = ds, PB = PB, META = META, don = don, meta = meta, lib = lib, donors = donors, T_don = T_don) }
CO <- lapply(setNames(COHORTS, COHORTS), load_cohort)
chk("Oetjen donor ages from the first-listed library: C = 60, S = 56", "60/56", CO$GSE120221$donors[donor %in% c("C", "S")][order(donor), paste(age, collapse = "/")])
chk("Oetjen: 10 young (<= 53), 5 middle, 5 old (>= 60) donors", "10/5/5", CO$GSE120221$donors[, sprintf("%d/%d/%d", sum(grp == "young"), sum(grp == "middle"), sum(grp == "old"))])
if (SMOKE) {   # the real age relation of the new cohorts is not read in a smoke run
  for (ds in COHORTS) { d <- CO[[ds]]$donors; p <- sample(d$age); CO[[ds]]$donors[, `:=`(age_real = age, age = p, grp = age_group_of(p))] }
  ts("    SMOKE: donor ages of the new cohorts permuted") }

# K7: per library donor table == meta_<GSE> (age_years, sex); per donor, Oetjen vs the BoneMarrowMap reference (ExactAge, DonorSex)
k7 <- rbindlist(lapply(COHORTS, function(ds) { x <- CO[[ds]]; m <- x$meta[disease == "healthy_donor" & gsm_id %in% x$lib$gsm]   # kept donors' libraries only (section 4, K7)
  m[, gsm := gsm_id]; j <- merge(m[, .(gsm, age_meta = as.numeric(age_years), sex_meta = sex, patient_id)], x$don[, .(gsm, age_tab = as.numeric(age), sex_tab = sex, donor)], by = "gsm")
  j[, ok := age_meta == age_tab & sex_meta == sex_tab & patient_id == donor]; j[, ds := ds]; j }))
n_lib_all <- sum(vapply(CO, function(x) nrow(x$lib), 1L))
ka("K7", "libraries whose meta_<GSE> row equals the donor table (age, sex, donor)", sprintf("%d/%d", sum(k7$ok), nrow(k7)), sprintf("all, %d libraries matched", n_lib_all), isTRUE(nrow(k7) == n_lib_all && all(k7$ok)),
   sprintf("%d of %d libraries matched; disagreeing: %s", nrow(k7), n_lib_all, paste(k7[ok == FALSE, gsm], collapse = " ")))
bmm <- as.data.table(readRDS(BMM_REF_RDS)$meta_data)[Study == "Oetjen2018", .(ExactAge = ExactAge[1], DonorSex = DonorSex[1]), by = Donor]
bmm[, donor := sub("^Oetjen2018_", "", Donor)]; bmm <- bmm[donor %in% CO$GSE120221$donors$donor]   # Ck/Sk entries are libraries, not donors
k7b <- merge(CO$GSE120221$donors[, .(donor, age = if (SMOKE) age_real else age, sex)], bmm, by = "donor")
k7b[, ok := as.numeric(as.character(ExactAge)) == age & toupper(substr(as.character(DonorSex), 1, 1)) == toupper(sex)]
ka("K7", "Oetjen donors whose age and sex equal the BoneMarrowMap reference metadata", sprintf("%d/%d", sum(k7b$ok), nrow(k7b)), "all, every donor found", isTRUE(nrow(k7b) == nrow(CO$GSE120221$donors) && all(k7b$ok)),
   sprintf("%d of %d donors found in the reference; disagreeing: %s", nrow(k7b), nrow(CO$GSE120221$donors), paste(k7b[ok == FALSE, donor], collapse = " ")))
# K3: identity from kept-cell sums (the four group matrices are exactly the kept cells); S4 rule, or the Y-gene rule for GSE227690
kept_sum <- function(PB, genes, samples) { gs <- intersect(genes, PB$genes); if (!length(gs)) return(rep(NA_real_, length(samples)))
  Reduce(`+`, lapply(PB$raw, function(M) colSums(M[gs, samples, drop = FALSE]))) }
k3 <- rbindlist(lapply(COHORTS, function(ds) { x <- CO[[ds]]; s <- x$lib$sample
  d <- data.table(ds = ds, sample = s, donor = x$lib$donor, table_sex = x$lib$sex_lib, n_kept = x$lib$n_kept)
  if (ds == "GSE227690") { d[, y_per_cell := kept_sum(x$PB, Y_GENES, s) / n_kept]
    d[, call := fifelse(y_per_cell >= 1.0, "M", fifelse(y_per_cell <= 0.2, "F", "unclear"))]
  } else { d[, xist := kept_sum(x$PB, "XIST", s)]; d[, rps4y1 := kept_sum(x$PB, "RPS4Y1", s)]
    d[, call := fifelse(xist > 3 * pmax(rps4y1, 1), "F", fifelse(rps4y1 > 3 * pmax(xist, 1), "M", "unclear"))] }
  d[, ok := call == table_sex]; d }), fill = TRUE)
k3d <- k3[, .(ok = all(ok)), by = .(ds, donor)]
ka("K3", "donors whose every library's sex call equals the donor table (29 donors)", sprintf("%d/%d", sum(k3d$ok), nrow(k3d)), "all", isTRUE(all(k3d$ok)), paste(k3d[ok == FALSE, paste(ds, donor)], collapse = " "))
k3_bad <- k3d[ok == FALSE]
if (nrow(k3_bad) == 1L) {   # section 4 K3: exactly one disagreeing donor is excluded before any score and recorded; two or more stop (section 9)
  CO[[k3_bad$ds]]$donors[donor == k3_bad$donor, `:=`(eligible = FALSE, excluded_K3 = TRUE)]
  chk(sprintf("K3: donor %s of %s excluded before any score (identity unresolved)", k3_bad$donor, k3_bad$ds), "excluded", "excluded")
  ts(sprintf("    K3: donor %s of %s excluded (identity unresolved); K0 and K4 use the reduced set", k3_bad$donor, k3_bad$ds)) }
# K0: eligibility of the decisive cohort (kept T cells per donor >= 20)
oe <- CO$GSE120221$donors
ka("K0", "eligible Oetjen donors (>= 20 kept T cells after pooling) / of which >= 60 y", sprintf("%d/%d", sum(oe$eligible), oe[eligible == TRUE & grp == "old", .N]),
   sprintf(">= %d incl. >= %d old", K0_MIN_DONORS, K0_MIN_OLD), isTRUE(sum(oe$eligible) >= K0_MIN_DONORS && oe[eligible == TRUE & grp == "old", .N] >= K0_MIN_OLD),
   sprintf("kept T per donor %d-%d", min(oe$n_T), max(oe$n_T)))

## -- 4. the joint space J-E (section 2.4) and K4 -----------------------------------------------------------------------------
jspace <- function(E_raw) {   # E_raw: genes x donors raw T sums of one external cohort; returns z of G and E with G's parameters
  common <- setdiff(intersect(rownames(G_RAW), rownames(E_raw)), union(black_for(rownames(G_RAW)), black_for(rownames(E_raw))))
  yG <- DGEList(G_RAW[!rownames(G_RAW) %in% black_for(rownames(G_RAW)), , drop = FALSE]); yE <- DGEList(E_raw[!rownames(E_raw) %in% black_for(rownames(E_raw)), , drop = FALSE])
  kG <- rownames(yG)[rowMeans(cpm(yG) > 1) >= 0.5]; kE <- rownames(yE)[rowMeans(cpm(yE) > 1) >= 0.5]
  U_J <- intersect(intersect(kG, kE), common)
  yG <- calcNormFactors(yG[U_J, , keep.lib.sizes = FALSE]); yE <- calcNormFactors(yE[U_J, , keep.lib.sizes = FALSE])
  LG <- cpm(yG, log = TRUE, prior.count = 1); LE <- cpm(yE, log = TRUE, prior.count = 1)
  m <- rowMeans(LG); s <- apply(LG, 1, sd); ok <- s > 1e-8; U_J <- U_J[ok]; LG <- LG[ok, , drop = FALSE]; LE <- LE[ok, , drop = FALSE]; m <- m[ok]; s <- s[ok]
  U <- data.table(gene = U_J, mean_lcpm = m, det = rowMeans(G_RAW[U_J, , drop = FALSE] > 0))
  U[, bin := paste(cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE), cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf)))]
  list(U_J = U_J, U = U, LG = LG, LE = LE, m = m, s = s, ZG = (LG - m) / s, ZE = (LE - m) / s, yE_lib = yE$samples$lib.size,
       d_g = rowMeans(LG[, S5$sample[S5$status == "AML"], drop = FALSE]) - rowMeans(LG[, S5$sample[S5$status == "healthy"], drop = FALSE])) }
zE_from_counts <- function(J, E_raw) {   # re-derive E's log-CPM on the fixed U_J with G's fixed parameters (K5 plant, K2 libraries)
  yE <- DGEList(E_raw[!rownames(E_raw) %in% black_for(rownames(E_raw)), , drop = FALSE]); yE <- calcNormFactors(yE[J$U_J, , keep.lib.sizes = FALSE])
  (cpm(yE, log = TRUE, prior.count = 1) - J$m) / J$s }
J <- list(); for (ds in COHORTS) { x <- CO[[ds]]; el <- x$donors[eligible == TRUE, donor]; J[[ds]] <- jspace(x$T_don[, el, drop = FALSE])
  cov <- mean(PRIM %in% J[[ds]]$U_J)
  ka("K4", sprintf("%s: primary genes in U_J", ds), sprintf("%d/370 (%.3f)", sum(PRIM %in% J[[ds]]$U_J), cov), ">= 0.80", isTRUE(cov >= K4_MIN_COV), sprintf("U_J %d genes; eligible donors %d", length(J[[ds]]$U_J), length(el))) }
D_RES <- sapply(COHORTS, function(ds) { s <- score(J[[ds]]$ZG, UP, DN); st <- S5$status[match(names(s), S5$sample)]; mean(s[st == "AML"]) - mean(s[st == "healthy"]) })
ts(sprintf("    D_res_J: %s", paste(sprintf("%s %.3f", COHORTS, D_RES), collapse = "; ")))

## -- 5. the cross-fitted yardstick r (section 3.5) ---------------------------------------------------------------------------
ts(sprintf("[5] leave-one-sample-out cross-fit (%d folds)", N_FOLD))
folds <- if (SMOKE) sample(seq_len(25L), N_FOLD) else seq_len(25L); oof <- rep(NA_real_, 25L); n_fold_genes <- integer(0); n_empty_half <- 0L
for (i in folds) { d24 <- D5[-i]; terms <- if (d24[, uniqueN(sex)] < 2) setdiff(MAIN, "sex") else MAIN
  r <- fit_de(G_RAW, d24, terms, "statusAML"); fu <- setdiff(r[hit == TRUE & logFC > 0, gene], AMB$gene); fd <- setdiff(r[hit == TRUE & logFC < 0, gene], AMB$gene)
  L24 <- S5Z$L[, -i, drop = FALSE]; mu <- rowMeans(L24); sg <- apply(L24, 1, sd); zi <- (S5Z$L[, i] - mu) / sg; zi <- zi[is.finite(zi)]
  oof[i] <- (if (length(intersect(fu, names(zi)))) mean(zi[intersect(fu, names(zi))]) else 0) - (if (length(intersect(fd, names(zi)))) mean(zi[intersect(fd, names(zi))]) else 0)
  n_fold_genes <- c(n_fold_genes, length(fu) + length(fd)); n_empty_half <- n_empty_half + (length(fu) == 0L || length(fd) == 0L); ts(sprintf("    fold %2d (%s): %d genes (%d up / %d down), oof score %.3f", i, D5$sample[i], length(fu) + length(fd), length(fu), length(fd), oof[i])) }
GAP_OOF <- mean(oof[st == "AML"], na.rm = TRUE) - mean(oof[st == "healthy"], na.rm = TRUE); SD_OOF_H <- sd(oof[st == "healthy"], na.rm = TRUE)
R_DEFL <- if (is.finite(GAP_OOF) && GAP_OOF > 0) max(1, GAP_IN / GAP_OOF) else NA_real_
ts(sprintf("    gap_in %.4f, gap_oof %.4f (%.2f oof-healthy-SD), r = %s", GAP_IN, GAP_OOF, GAP_OOF / SD_OOF_H, f3(R_DEFL)))
ka("3.5", "gap_oof > 0 and r <= 2 (yardstick cross-validates and is resolved)", sprintf("gap_oof %.3f, r %s", GAP_OOF, f3(R_DEFL)), "gap_oof > 0, r <= 2", isTRUE(is.finite(R_DEFL) && R_DEFL <= R_MAX),
   sprintf("%sfold imprint sizes %d-%d; folds with an empty up or down half: %d", if (SMOKE) sprintf("SMOKE: %d of 25 folds only; ", N_FOLD) else "", min(n_fold_genes), max(n_fold_genes), n_empty_half))
D_GATE <- D_RES / ifelse(is.finite(R_DEFL), R_DEFL, 1)   # the deflated denominators D_res_J / r per cohort
rates <- function(r) { ff <- approx(c(1, 1.3, 1.6, 2), c(1.2, 4.7, 9.4, 15.9), xout = min(max(r, 1), 2))$y; fp40 <- approx(c(1, 1.6), c(25, 30), xout = min(max(r, 1), 1.6), rule = 2)$y; fp50 <- approx(c(1, 1.6), c(7, 16), xout = min(max(r, 1), 1.6), rule = 2)$y
  sprintf("false-fail under H_AML ~%.0f%%; false pass at a true age share of 40%% ~%.0f%%, of 50%% ~%.0f%% (section 3.1, restated for r)", ff, fp40, fp50) }
YARD <- data.table(quantity = c("gap_in_primary", "healthy_SD_in", "gap_in_healthy_SD_units", "gap_oof", "healthy_SD_oof", "gap_oof_healthy_SD_units", "r", paste0("D_res_J_", COHORTS), paste0("D_res_J_over_r_", COHORTS), "error_rates_for_r", "cross_fit_folds_with_empty_half", "gap_in_registered", "gap_in_registered_healthy_SD_units", "gap_in_registered_note"),
                   value = c(f3(GAP_IN), f3(SD_H), f3(GAP_IN / SD_H), f3(GAP_OOF), f3(SD_OOF_H), f3(GAP_OOF / SD_OOF_H), f3(R_DEFL), f3(D_RES), f3(D_GATE), if (is.finite(R_DEFL)) rates(R_DEFL) else "r undefined", as.character(n_empty_half), sprintf("%.3f", GAP_IN_REG), "3.87",
                             "registered on 2026-10-23 without the blacklist step of logz(); the exact logz() gives gap_in_primary, which feeds r and the A2a ratio so that gap_in and gap_oof share one gene universe (section 10.4)"))
fwrite(YARD, out("yardstick"))

## -- 6. K6 regression, K2 resolution ------------------------------------------------------------------------------------------
D8 <- fread(file.path(TABWS, "ws_s8_scores_disc.csv")); s15 <- score(logz(G_RAW[, D8$sample])$Z, UP, DN); d8 <- D8$score_primary[match(names(s15), D8$sample)]
ka("K6", "S9 scoring vs S8 (15 AML, z over the 15): Spearman / max |diff|", sprintf("%.4f / %.5f", cor(s15, d8, method = "spearman"), max(abs(s15 - d8))), "> 0.99 / < 0.05", isTRUE(cor(s15, d8, method = "spearman") > 0.99 && max(abs(s15 - d8)) < 0.05))
ka("K6", "yardstick on the 25 reproduces the registered 1.970 (|gap_in - 1.970|)", sprintf("%.4f", abs(GAP_IN - GAP_IN_REG)), "< 0.01", isTRUE(abs(GAP_IN - GAP_IN_REG) < GAP_TOL), sprintf("gap_in %.4f = %.2f healthy-SD", GAP_IN, GAP_IN / SD_H))
oeJ <- J[[DECISIVE]]; oeL <- CO$GSE120221$lib[!donor %in% CO$GSE120221$donors[excluded_K3 == TRUE, donor]]; ZLIB <- zE_from_counts(oeJ, CO$GSE120221$PB$raw$T[, oeL$sample]); slib <- score(ZLIB, UP, DN)
libsc <- function(lab) { v <- unname(slib[oeL[sub("^GSM[0-9]+_", "", sample) == lab, sample]]); if (length(v) == 1L) v else NA_real_ }   # NA if the library's donor was excluded by K3
k2 <- data.table(pair = c("C1 vs C2 (technical)", "S1 vs S2 (technical)", "Ck vs mean(C1, C2) (second aspirate)", "Sk1 vs mean(S1, S2) (second aspirate)", "Sk2 vs mean(S1, S2) (second aspirate)"),
                 diff = c(abs(libsc("C1") - libsc("C2")), abs(libsc("S1") - libsc("S2")), abs(libsc("Ck") - mean(c(libsc("C1"), libsc("C2")))), abs(libsc("Sk1") - mean(c(libsc("S1"), libsc("S2")))), abs(libsc("Sk2") - mean(c(libsc("S1"), libsc("S2"))))))
k2[, line := D_GATE[[DECISIVE]] / 3]
ka("K2", "technical replicates C1-C2 and S1-S2: |score difference| in J-Oetjen", sprintf("%.3f / %.3f", k2$diff[1], k2$diff[2]), sprintf("both < (D_res/r)/3 = %.3f", D_GATE[[DECISIVE]] / 3), isTRUE(any(!is.na(k2$diff[1:2])) && all(k2$diff[1:2] < D_GATE[[DECISIVE]] / 3, na.rm = TRUE)),
   sprintf("second aspirates: Ck %.3f, Sk1 %.3f, Sk2 %.3f", k2$diff[3], k2$diff[4], k2$diff[5]))

## -- 7. the age model (A2a engine) and K5, the planted H_age under permuted ages ------------------------------------------------
age_de <- function(E_raw, dn_tab) {   # edgeR QL on donor pseudobulks: ~ age + sex; blacklist removed; CPM > 1 in half; returns per-gene age coefficient (log2 per year), F, p
  y <- DGEList(E_raw[!rownames(E_raw) %in% black_for(rownames(E_raw)), dn_tab$donor, drop = FALSE]); y <- y[rowMeans(cpm(y) > 1) >= 0.5, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  terms <- if (uniqueN(dn_tab$sex) < 2) "age" else c("age", "sex"); X <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), data = as.data.frame(dn_tab))
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = "age")
  as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC_age = logFC, logCPM, F_age = F, PValue, FDR)] }
fit_age <- function(sc, dn_tab) {   # f_slope and f_contrast of a score vector (named by donor) under a donor table with age, sex, grp
  d <- dn_tab[donor %in% names(sc)]; d[, score := sc[donor]]
  fm <- if (uniqueN(d$sex) >= 2) lm(score ~ age + sex, data = d) else lm(score ~ age, data = d)
  c(slope = unname(coef(fm)["age"]), contrast = mean(d[grp == "old", score]) - mean(d[grp == "young", score])) }
alignment <- function(age_tab, genes) {   # A2b: OLS slope of (log2FC_age x 31.5) on log2FC_AML over `genes`, divided by the set's reliability lambda
  g <- intersect(genes, intersect(age_tab$gene, REFIT[!is.na(SE), gene])); x <- DE5$logFC[match(g, DE5$gene)]; y <- age_tab$logFC_age[match(g, age_tab$gene)] * DELTA_AGE; se <- REFIT$SE[match(g, REFIT$gene)]
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; se <- se[ok]; if (length(x) < 10) return(c(b = NA, lambda = NA, b_over_lambda = NA, rho = NA, concord = NA, n = length(x)))
  b <- unname(coef(lm(y ~ x))[2]); lam <- 1 - mean(se^2) / var(x)
  c(b = b, lambda = lam, b_over_lambda = b / lam, rho = cor(x, y, method = "spearman"), concord = mean(sign(x) == sign(y)), n = length(x)) }
oeD <- CO$GSE120221$donors[eligible == TRUE]; oeD[, sex := factor(sex)]
set.seed(SEED); perm <- sample(oeD$age); tries <- 1L
while (abs(cor(perm, oeD$age, method = "spearman")) > 0.3 && tries < 1000L) { perm <- sample(oeD$age); tries <- tries + 1L }
PERM <- copy(oeD)[, `:=`(age = perm, grp = age_group_of(perm))]
chk("K5 permutation: |Spearman(age_perm, age)| <= 0.3", TRUE, abs(cor(PERM$age, oeD$age, method = "spearman")) <= 0.3)
O_RAW <- CO$GSE120221$T_don[, oeD$donor]; d_g <- oeJ$d_g
f_unpl <- fit_age(score(oeJ$ZE, UP, DN), PERM)
plant <- function(frac) {   # multiply each primary gene's counts by 2^(frac * d_g * (age_perm - 32.5)/31.5), donor by donor (section 4, K5)
  P <- O_RAW; gp <- intersect(PRIM, oeJ$U_J); stopifnot(identical(colnames(P), PERM$donor))
  mult <- 2^(frac * outer(d_g[gp], (PERM$age - AGE_REF_HEALTHY) / DELTA_AGE))          # genes x donors
  P[gp, ] <- round(P[gp, , drop = FALSE] * mult); P }
expected_contrast <- (mean(PERM[grp == "old", age]) - mean(PERM[grp == "young", age])) / DELTA_AGE
k5 <- rbindlist(lapply(c(1, 1 / 3), function(frac) { P <- plant(frac); fp <- fit_age(score(zE_from_counts(oeJ, P), UP, DN), PERM)
  rs <- (fp["slope"] - f_unpl["slope"]) * DELTA_AGE / D_RES[[DECISIVE]]; rc <- (fp["contrast"] - f_unpl["contrast"]) / D_RES[[DECISIVE]]
  data.table(plant = frac, slope_recovery = unname(rs), contrast_recovery = unname(rc), expected_contrast = expected_contrast * frac) }))
ka("K5", "full plant: slope recovery (f_slope planted - unplanted, undeflated)", f3(k5$slope_recovery[1]), "[0.95, 1.03]", isTRUE(k5$slope_recovery[1] >= 0.95 && k5$slope_recovery[1] <= 1.03))
ka("K5", "full plant: contrast recovery vs (mean old age - mean young age)/31.5", sprintf("%.3f vs %.3f", k5$contrast_recovery[1], k5$expected_contrast[1]), "within 0.05", isTRUE(abs(k5$contrast_recovery[1] - k5$expected_contrast[1]) <= 0.05))
ka("K5", "third plant: slope recovery", f3(k5$slope_recovery[2]), "[0.30, 0.37]", isTRUE(k5$slope_recovery[2] >= 0.30 && k5$slope_recovery[2] <= 0.37))
ka("K5", "third plant: contrast recovery vs expected", sprintf("%.3f vs %.3f", k5$contrast_recovery[2], k5$expected_contrast[2]), "within 0.03", isTRUE(abs(k5$contrast_recovery[2] - k5$expected_contrast[2]) <= 0.03))
ts("    K5: age model on the planted counts (planted A2b)")
AGE_PL <- age_de(plant(1), PERM); al_pl <- alignment(AGE_PL, intersect(PRIM, oeJ$U_J))
nullJ <- function(age_tab, Ubins, n) { pool <- Ubins[gene %in% age_tab$gene & gene %in% DE5$gene]
  t(replicate(n, { g <- matched(pool, PRIM); alignment(age_tab, g)[c("b", "lambda", "b_over_lambda", "rho")] })) }
set.seed(SEED); nul_pl <- nullJ(AGE_PL, oeJ$U, N_RAND)
ka("K5", "planted A2b: b/lambda of the imprint vs the random-set P95", sprintf("%.3f vs %.3f", al_pl["b_over_lambda"], quantile(nul_pl[, "b_over_lambda"], .95)), "above P95", isTRUE(unname(al_pl["b_over_lambda"] > quantile(nul_pl[, "b_over_lambda"], .95))),
   sprintf("lambda imprint %.2f, random median %.2f", al_pl["lambda"], median(nul_pl[, "lambda"])))
set.seed(SEED); scr <- replicate(100, { lab <- sample(c(rep("up", length(UP)), rep("down", length(DN)))); fit_age(score(oeJ$ZE, PRIM[lab == "up"], PRIM[lab == "down"]), PERM)["slope"] * DELTA_AGE / D_GATE[[DECISIVE]] })
ka("K5", "100 sign-scrambled imprints (permuted ages): f_slope P5 / P95, denominator D_res_J/r (calibration report)", sprintf("%.3f / %.3f", quantile(scr, .05), quantile(scr, .95)), "reported", NA)

## -- 8. K1 (naive CD8 share vs age), K1b, Griffin author labels -----------------------------------------------------------------
naive_share <- function(x) { c <- x$PB$cells[group4 == "T" & high_error == 0]; c[, donor := x$lib$donor[match(sample, x$lib$sample)]]
  c[, .(n_T = .N, cd8_naive = mean(fine == "CD8 Naive"), cd4_naive = mean(fine == "CD4 Naive")), by = donor] }
NS <- lapply(COHORTS, function(ds) merge(naive_share(CO[[ds]]), CO[[ds]]$donors[excluded_K3 == FALSE, .(donor, age, grp)], by = "donor")[, ds := ds]); names(NS) <- COHORTS
k1 <- suppressWarnings(cor.test(NS$GSE120221$cd8_naive, NS$GSE120221$age, method = "spearman", alternative = "less"))   # exact p when there are no ties, else the t approximation
ka("K1", "Oetjen: Spearman of the CD8 Naive share of projected T cells with age (one-sided, rho < 0)", sprintf("rho %.3f, p %.4f", k1$estimate, k1$p.value), "rho < 0, p < 0.05", isTRUE(unname(k1$estimate < 0 && k1$p.value < 0.05)), sprintf("n %d donors", nrow(NS$GSE120221)))
k1c <- suppressWarnings(cor.test(NS$GSE120221$cd4_naive, NS$GSE120221$age, method = "spearman", alternative = "less"))
ka("K1", "Oetjen: the same for the CD4 Naive share (reported)", sprintf("rho %.3f, p %.4f", k1c$estimate, k1c$p.value), "reported", NA)
for (ds in c("GSE227690", "GSE276360")) { n <- NS[[ds]]; below <- all(outer(n[grp == "old", cd8_naive], n[grp == "young", cd8_naive], `<`))
  ka("K1", sprintf("%s: both old donors' CD8 Naive share below every young donor's (descriptive)", ds), paste(sprintf("%s(%s y) %.3f", n$donor, n$age, n$cd8_naive), collapse = "; "), "descriptive", NA, if (below) "old below every young: yes" else "old below every young: NO -> A_E labelled 'age effect not visible in labels'")
  CO[[ds]]$k1_labels_ok <- below }
gr <- rbindlist(lapply(CO$GSE227690$lib$sample, function(s) { md <- as.data.table(readRDS(file.path(LARGE1_DIR, "02_seurat_objects/04_annotated/GSE227690", paste0(s, ".rds")))@meta.data)
  md[high_error == 0 & bmm_fine %in% T_FINE, .(sample = s, n = .N, ours_cd8n = sum(bmm_fine == "CD8 Naive"), auth_cd8n = sum(anno_CellType == "CD8Naive"), agree_cd8n = sum(bmm_fine == "CD8 Naive" & anno_CellType == "CD8Naive"),
                                               ours_cd4n = sum(bmm_fine == "CD4 Naive"), auth_cd4n = sum(anno_CellType == "CD4Naive"), agree_cd4n = sum(bmm_fine == "CD4 Naive" & anno_CellType == "CD4Naive"))] }))
ka("K1", "GSE227690 author labels vs ours on projected T cells: CD8Naive agreement (of author-labelled) / CD4Naive", sprintf("%.3f / %.3f", sum(gr$agree_cd8n) / sum(gr$auth_cd8n), sum(gr$agree_cd4n) / sum(gr$auth_cd4n)), "reported, no statement", NA,
   sprintf("author CD8Naive share of T cells %.3f, ours %.3f", sum(gr$auth_cd8n) / sum(gr$n), sum(gr$ours_cd8n) / sum(gr$n)))
ts("    A2a age model on the real Oetjen donor pseudobulks (K1b)")
AGE <- age_de(O_RAW, oeD)
ka("K1b", "Oetjen T pseudobulk: mean age coefficient (log2 per year) over naive genes / cytotoxic genes", sprintf("%.4f / %.4f", mean(AGE[gene %in% NAIVE_GENES, logFC_age]), mean(AGE[gene %in% CYTO_GENES, logFC_age])), "< 0 / > 0 (reported)", NA,
   sprintf("naive present %d/5, cytotoxic present %d/6", sum(NAIVE_GENES %in% AGE$gene), sum(CYTO_GENES %in% AGE$gene)))

## -- 9. phase-1 verdict ---------------------------------------------------------------------------------------------------------
KAT <- rbindlist(KA); fwrite(KAT, out("known_answers"))
stop_ids <- KAT[pass == FALSE & id %in% c("K0", "K1", "K2", "K4", "3.5"), paste(id, statistic)]; defect_ids <- KAT[pass == FALSE & id %in% c("K3", "K5", "K6", "K7"), paste(id, statistic)]
if (nrow(k3_bad) == 1L) defect_ids <- defect_ids[!grepl("^K3", defect_ids)]   # one disagreeing donor was excluded above (section 4 K3): not a stop
k4_oe <- KAT[id == "K4" & grepl(DECISIVE, statistic), pass]; k4_desc <- KAT[id == "K4" & !grepl(DECISIVE, statistic) & pass == FALSE, sub(":.*$", "", statistic)]
stop_ids <- stop_ids[!(grepl("^K4", stop_ids) & !grepl(DECISIVE, stop_ids))]   # a coverage failure in a descriptive cohort only demotes that cohort
if (length(defect_ids)) { ts("** identity / machinery known answer failed: ", paste(defect_ids, collapse = " | "))
  C <- rbindlist(CHK); fwrite(C, out("checks"))
  if (!SMOKE) { ts("[stop] phase 2 not run (section 4: correct and record before the real age relation is read)"); q(status = 0) } else ts("    SMOKE: phase 2 is exercised anyway (permuted ages; outputs suffixed _smoke)") }
if (length(stop_ids)) { ts("** stop condition: ", paste(stop_ids, collapse = " | "))
  C <- rbindlist(CHK); fwrite(C, out("checks"))
  if (!SMOKE) { ts("[stop] phase 2 not run (section 5: no age relation is read on the new cohorts)"); q(status = 0) } else ts("    SMOKE: phase 2 is exercised anyway (permuted ages; a K1 failure is expected under permutation)") }
if (!length(defect_ids) && !length(stop_ids)) ts("[9] phase 1 passed: the age relation is read now")

## -- 10. A1: the age gradient within Oetjen (section 3.1) --------------------------------------------------------------------------
TESTS <- list(); tst <- function(...) TESTS[[length(TESTS) + 1]] <<- data.table(...)
oeD[, grp := age_group_of(age)]
scJ <- score(oeJ$ZE, UP, DN); fa <- fit_age(scJ, oeD); DG <- D_GATE[[DECISIVE]]
f_slope <- unname(fa["slope"] * DELTA_AGE / DG); f_contrast <- unname(fa["contrast"] / DG); f_age <- max(f_slope, f_contrast)
fa_un <- lm(scJ[oeD$donor] ~ oeD$age)$coefficients[2] * DELTA_AGE / DG
ts(sprintf("    f_slope %.3f, f_contrast %.3f -> f_age %.3f (line %.3f); unadjusted-for-sex slope f %.3f", f_slope, f_contrast, f_age, LINE, fa_un))
set.seed(SEED); boot <- replicate(N_BOOT, { repeat { i <- sample(nrow(oeD), replace = TRUE); b <- oeD[i]; if (b[grp == "old", .N] >= 2 && b[grp == "young", .N] >= 2) break }
  b[, score := scJ[donor]]; fm <- if (uniqueN(b$sex) >= 2) lm(score ~ age + sex, data = b) else lm(score ~ age, data = b)
  max(unname(coef(fm)["age"]) * DELTA_AGE / DG, (mean(b[grp == "old", score]) - mean(b[grp == "young", score])) / DG) })
CI <- quantile(boot, c(.025, .975))
loo <- sapply(seq_len(nrow(oeD)), function(i) { f <- fit_age(scJ[oeD$donor[-i]], oeD[-i]); max(f["slope"] * DELTA_AGE / DG, f["contrast"] / DG) })
pass_A1 <- f_age < LINE; loo_flip <- pass_A1 && any(loo >= LINE)
set.seed(SEED); nulA1 <- t(replicate(N_RAND, { r <- split_random(oeJ$U, UP, DN); f <- fit_age(score(oeJ$ZE, r$up, r$dn), oeD); c(slope = unname(f["slope"] * DELTA_AGE / DG), contrast = unname(f["contrast"] / DG)) }))
tst(test = "A1_f_slope", cohort = DECISIVE, gene_set = "primary", statistic = f_slope, threshold = "component of f_age", null_P5 = quantile(nulA1[, "slope"], .05), null_P95 = quantile(nulA1[, "slope"], .95), n = nrow(oeD), pass = NA, detail = sprintf("beta_age x 31.5 / (D_res_J/r); sex-adjusted; unadjusted %.3f; rank among random sets %.3f", fa_un, mean(nulA1[, "slope"] < f_slope)))
tst(test = "A1_f_contrast", cohort = DECISIVE, gene_set = "primary", statistic = f_contrast, threshold = "component of f_age", null_P5 = quantile(nulA1[, "contrast"], .05), null_P95 = quantile(nulA1[, "contrast"], .95), n = oeD[grp != "middle", .N], pass = NA, detail = sprintf("old (n %d, mean age %.1f) minus young (n %d, mean age %.1f)", oeD[grp == "old", .N], oeD[grp == "old", mean(age)], oeD[grp == "young", .N], oeD[grp == "young", mean(age)]))
tst(test = "A1_f_age", cohort = DECISIVE, gene_set = "primary", statistic = f_age, threshold = "< 1/3 (point estimate)", null_P5 = NA, null_P95 = NA, n = nrow(oeD), pass = pass_A1, detail = sprintf("95%% donor-bootstrap CI %.3f to %.3f (conditional on the yardstick); r %.3f; LOO range %.3f to %.3f; LOO flip %s", CI[1], CI[2], R_DEFL, min(loo), max(loo), loo_flip))
# full 402 set (sensitivity)
UPf <- GS$full[direction == "up", gene]; DNf <- GS$full[direction == "down", gene]; faf <- fit_age(score(oeJ$ZE, UPf, DNf), oeD)
sGf <- score(oeJ$ZG, UPf, DNf); stf <- S5$status[match(names(sGf), S5$sample)]; DGf <- (mean(sGf[stf == "AML"]) - mean(sGf[stf == "healthy"])) / ifelse(is.finite(R_DEFL), R_DEFL, 1)   # the 402 set's own D_res_J / r
tst(test = "A1_f_age_full402", cohort = DECISIVE, gene_set = "full", statistic = max(faf["slope"] * DELTA_AGE / DGf, faf["contrast"] / DGf), threshold = "sensitivity", null_P5 = NA, null_P95 = NA, n = nrow(oeD), pass = NA, detail = sprintf("f_slope %.3f, f_contrast %.3f over the 402 set's own D_res_J/r %.3f (with the primary set's denominator: f_age %.3f)", faf["slope"] * DELTA_AGE / DGf, faf["contrast"] / DGf, DGf, max(faf["slope"] * DELTA_AGE / DG, faf["contrast"] / DG)))

## -- 11. A2b: gene-level age-AML alignment (wording modifier), then A2a (reported) --------------------------------------------------
al <- alignment(AGE, intersect(PRIM, oeJ$U_J)); set.seed(SEED); nulA2 <- nullJ(AGE, oeJ$U, N_RAND); p95_al <- quantile(nulA2[, "b_over_lambda"], .95)
mod_A2b <- al["b_over_lambda"] < p95_al
tst(test = "A2b_alignment_b_over_lambda", cohort = DECISIVE, gene_set = "primary", statistic = unname(al["b_over_lambda"]), threshold = "below random P95 -> wording", null_P5 = quantile(nulA2[, "b_over_lambda"], .05), null_P95 = p95_al, n = al["n"], pass = mod_A2b,
    detail = sprintf("b %.3f, lambda %.2f (random median lambda %.2f, median b/lambda %.3f = genome-wide rate); Spearman %.3f; sign concordance %.3f", al["b"], al["lambda"], median(nulA2[, "lambda"]), median(nulA2[, "b_over_lambda"]), al["rho"], al["concord"]))
age_sens <- AGE[PValue < 0.05, gene]; age_fdr <- AGE[FDR < 0.10, gene]
rob_up <- setdiff(UP, age_sens); rob_dn <- setdiff(DN, age_sens); sc_rob <- score(S5Z$Z, rob_up, rob_dn); gap_rob <- mean(sc_rob[st == "AML"]) - mean(sc_rob[st == "healthy"])
set.seed(SEED); rr <- replicate(N_RM, { rm <- sample(PRIM, length(intersect(PRIM, age_sens))); s <- score(S5Z$Z, setdiff(UP, rm), setdiff(DN, rm)); (mean(s[st == "AML"]) - mean(s[st == "healthy"])) / GAP_IN })
rem <- AGE[gene %in% intersect(PRIM, age_sens)]; rem[, imprint_dir := fifelse(gene %in% UP, 1, -1)]
tst(test = "A2a_removal_ratio", cohort = DECISIVE, gene_set = "primary", statistic = gap_rob / GAP_IN, threshold = "reported", null_P5 = quantile(rr, .05), null_P95 = quantile(rr, .95), n = length(intersect(PRIM, age_sens)), pass = NA,
    detail = sprintf("%d of 370 removed at raw p < 0.05 (%d at FDR < 0.10); removed genes' age sign = AML sign in %.2f; random-removal median ratio %.3f", length(intersect(PRIM, age_sens)), length(intersect(PRIM, age_fdr)), rem[, mean(sign(logFC_age) == imprint_dir)], median(rr)))
fwrite(data.table(gene = c(rob_up, rob_dn), direction = rep(c("up", "down"), c(length(rob_up), length(rob_dn))), note = "A2a record only; never replaces the primary set (section 7)"), out("robust_imprint"))
AG <- merge(data.table(gene = PRIM, direction = fifelse(PRIM %in% UP, "up", "down"), in_UJ_oetjen = PRIM %in% oeJ$U_J), AGE, by = "gene", all.x = TRUE)
AG <- merge(AG, REFIT[, .(gene, logFC_AML = logFC, F_AML = F, SE_AML = SE)], by = "gene", all.x = TRUE)
AG[, `:=`(age_sensitive_p05 = gene %in% age_sens, age_sensitive_fdr10 = gene %in% age_fdr)]; fwrite(AG[order(PValue)], out("age_genes"))

## -- 12. A3: direction against each external cohort (section 3.3) -------------------------------------------------------------------
DIR <- list(); POS <- list()
for (ds in COHORTS) { x <- CO[[ds]]; Jx <- J[[ds]]; dn <- x$donors[eligible == TRUE]; s_E <- score(Jx$ZE, UP, DN); s_G <- score(Jx$ZG, UP, DN); stG <- S5$status[match(names(s_G), S5$sample)]
  mG_h <- mean(s_G[stG == "healthy"]); sdG_h <- sd(s_G[stG == "healthy"]); Dres <- D_RES[[ds]]; Dg <- D_GATE[[ds]]
  yv <- s_E[dn[grp == "young", donor]]; ov <- s_E[dn[grp == "old", donor]]
  S_E <- mean(yv) - mG_h; A_E <- mean(ov) - mean(yv); D_dir <- mean(s_G[stG == "AML"]) - mean(ov); sd_ratio <- if (length(yv) >= 2) sd(yv) / sdG_h else NA_real_
  interpretable <- abs(S_E) < Dres / 2 && is.finite(sd_ratio) && sd_ratio >= 0.5 && sd_ratio <= 2; descriptive <- length(yv) <= 3
  mid <- mG_h + Dres / 2; corr <- s_E - S_E; n_old_above <- sum(corr[dn[grp == "old", donor]] > mid)
  holds <- A_E / Dg < LINE && n_old_above <= 1
  labels_ok <- if (ds == DECISIVE) TRUE else isTRUE(x$k1_labels_ok); k4_ok <- isTRUE(KAT[id == "K4" & grepl(ds, statistic), pass])
  DIR[[ds]] <- data.table(cohort = ds, n_young = length(yv), n_old = length(ov), D_res_J = Dres, D_res_J_over_r = Dg, S_E = S_E, A_E = A_E, A_E_over_Dg = A_E / Dg, D_dir = D_dir, sd_ratio_young_vs_G_healthy = sd_ratio,
                          interpretable = interpretable, descriptive_S_E = descriptive, n_old_above_midpoint = n_old_above, direction_holds = holds, age_effect_visible_in_labels = labels_ok, coverage_K4_pass = k4_ok)
  POS[[ds]] <- data.table(cohort = ds, donor = names(s_E), age = dn$age[match(names(s_E), dn$donor)], grp = dn$grp[match(names(s_E), dn$donor)], score_J = s_E, score_J_corrected = corr, above_midpoint = corr > mid, midpoint = mid, G_healthy_mean = mG_h, G_AML_mean = mean(s_G[stG == "AML"]))
  tst(test = "A3_direction", cohort = ds, gene_set = "primary", statistic = A_E / Dg, threshold = "A_E/(D_res/r) < 1/3 and <= 1 old donor above the midpoint", null_P5 = NA, null_P95 = NA, n = length(yv) + length(ov), pass = if (descriptive) NA else holds,
      detail = sprintf("S_E %.3f (|S_E| < D_res/2 = %.3f: %s), sd ratio %s, interpretable %s, descriptive (<= 3 young) %s; D_dir %.3f; old above midpoint %d of %d%s", S_E, Dres / 2, abs(S_E) < Dres / 2, f3(sd_ratio), interpretable, descriptive, D_dir, n_old_above, length(ov), if (labels_ok) "" else "; age effect not visible in labels (K1)")) }
DIRT <- rbindlist(DIR); fwrite(DIRT, out("direction")); POST <- rbindlist(POS); fwrite(POST, out("direction_positions"))

## -- 13. reported: W space, subtype and composition sensitivities, ambient partials, section 3.4 ---------------------------------------
W <- logz(O_RAW); scW <- score(W$Z, UP, DN); rhoW <- suppressWarnings(cor.test(scW[oeD$donor], oeD$age, method = "spearman"))
set.seed(SEED); nulW <- replicate(N_RAND, { r <- split_random(W$U, UP, DN); cor(score(W$Z, r$up, r$dn)[oeD$donor], oeD$age, method = "spearman") })
slopeW <- unname(coef(lm(scW[oeD$donor] ~ oeD$age))[2]) * 10; sdW_young <- sd(scW[oeD[grp == "young", donor]])
tst(test = "W_rho_score_age", cohort = DECISIVE, gene_set = "primary", statistic = unname(rhoW$estimate), threshold = "reported", null_P5 = quantile(nulW, .05), null_P95 = quantile(nulW, .95), n = nrow(oeD), pass = NA,
    detail = sprintf("two-sided p %.3f; Fisher-z 95%% CI %.3f to %.3f; slope %.3f W-units per decade; young-donor SD %.3f; |rho| rank among random sets %.3f", rhoW$p.value, tanh(atanh(rhoW$estimate) - 1.96 / sqrt(nrow(oeD) - 3)), tanh(atanh(rhoW$estimate) + 1.96 / sqrt(nrow(oeD) - 3)), slopeW, sdW_young, mean(abs(nulW) < abs(rhoW$estimate))))
# naive-share covariate and the subtype pseudobulks (CD4 Naive; CD8 EM 1 + 2), W space
kn <- CO$GSE120221$PB$cells[kept == TRUE & group4 == "T"][, donor := oeL$donor[match(sample, oeL$sample)]][, .(ns = mean(fine %in% c("CD4 Naive", "CD8 Naive"))), by = donor]   # share among the donor's KEPT T cells (section 3.1)
nsh <- kn$ns[match(oeD$donor, kn$donor)]; dd <- copy(oeD)[, `:=`(score = scJ[donor], naive_share = nsh)]
f_nsh <- unname(coef(lm(score ~ age + sex + naive_share, data = dd))["age"]) * DELTA_AGE / DG
tst(test = "A1_f_slope_naive_share_adjusted", cohort = DECISIVE, gene_set = "primary", statistic = f_nsh, threshold = "sensitivity", null_P5 = NA, null_P95 = NA, n = nrow(oeD), pass = NA, detail = "lm(score ~ age + sex + naive share among kept T cells), J space")
sub_pb <- function(fines) {   # subtype pseudobulk per donor; 26's $fine matrices hold only the libraries where the subtype occurs, so align on all libraries first
  PBo <- CO$GSE120221$PB; lib <- CO$GSE120221$lib; fines <- intersect(fines, names(PBo$fine)); if (!length(fines)) return(matrix(nrow = 0, ncol = 0))
  Ms <- matrix(0, nrow = length(PBo$genes), ncol = nrow(lib), dimnames = list(PBo$genes, lib$sample))
  for (f in fines) { M <- as.matrix(PBo$fine[[f]]); cols <- intersect(colnames(M), lib$sample); Ms[rownames(M), cols] <- Ms[rownames(M), cols] + M[, cols, drop = FALSE] }
  keep_n <- PBo$cells[kept == TRUE & fine %in% fines, .N, by = sample]; dn_n <- keep_n[, .(n = sum(N)), by = .(donor = lib$donor[match(sample, lib$sample)])]; el <- dn_n[n >= MIN_T & donor %in% oeD$donor, donor]; if (!length(el)) return(matrix(nrow = 0, ncol = 0))
  M <- sapply(el, function(d) rowSums(Ms[, lib[donor == d, sample], drop = FALSE])); colnames(M) <- el; M }
for (nm in c("CD4 Naive", "CD8 EM")) { fines <- if (nm == "CD8 EM") c("CD8 Effector Memory 1", "CD8 Effector Memory 2") else nm; M <- sub_pb(fines)
  if (ncol(M) >= 6) { Wz <- logz(M); s <- score(Wz$Z, UP, DN); a <- oeD$age[match(colnames(M), oeD$donor)]; ct <- suppressWarnings(cor.test(s, a, method = "spearman"))
    tst(test = paste0("W_rho_subtype_", gsub(" ", "", nm)), cohort = DECISIVE, gene_set = "primary", statistic = unname(ct$estimate), threshold = "sensitivity", null_P5 = NA, null_P95 = NA, n = ncol(M), pass = NA, detail = sprintf("%s pseudobulk, donors with >= 20 kept cells (computed only with >= 6 such donors, an unregistered minimum); two-sided p %.3f", nm, ct$p.value)) } }
donor_amb <- function(x) {   # per-donor cell-weighted ambient shares and doublet score from 26's library rows (section 6; dbl_mean NA where T3a has not run)
  u <- unique(x$META[, .(sample, n_cells, hb_share_nonery, gran_share_lymph, dbl_mean)]); stopifnot(!anyDuplicated(u$sample)); u[, donor := x$lib$donor[match(sample, x$lib$sample)]]
  u[!is.na(donor), .(hb_share_nonery = sum(hb_share_nonery * n_cells) / sum(n_cells), gran_share_lymph = sum(gran_share_lymph * n_cells) / sum(n_cells),
                     dbl_mean = if (all(is.na(dbl_mean))) NA_real_ else sum(dbl_mean * n_cells, na.rm = TRUE) / sum(n_cells[!is.na(dbl_mean)])), by = donor] }
amb <- donor_amb(CO$GSE120221)
pcor <- function(y, x, z) { r <- resid(lm(y ~ z)); cor(r, resid(lm(x ~ z)), method = "spearman") }
a_ <- amb[match(oeD$donor, donor)]; tst(test = "info_score_vs_ambient_partial", cohort = DECISIVE, gene_set = "primary", statistic = NA, threshold = "informational", null_P5 = NA, null_P95 = NA, n = nrow(oeD), pass = NA,
    detail = sprintf("partial Spearman of the J score with hb_share_nonery %.3f and gran_share_lymph %.3f given age", pcor(scJ[oeD$donor], a_$hb_share_nonery, oeD$age), pcor(scJ[oeD$donor], a_$gran_share_lymph, oeD$age)))
# section 3.4 inside GSE185381
y4 <- S5[status == "AML" & age < 60, sample]; h4 <- S5[status == "healthy" & age >= 39, sample]
lfc34 <- rowMeans(S5Z$L[, y4, drop = FALSE]) - rowMeans(S5Z$L[, h4, drop = FALSE]); g34 <- intersect(PRIM, names(lfc34)); conc34 <- mean(sign(lfc34[g34]) == ifelse(g34 %in% UP, 1, -1))
rA <- suppressWarnings(cor.test(sc_in[S5[status == "AML", sample]], S5[status == "AML", age], method = "spearman")); rH <- suppressWarnings(cor.test(sc_in[S5[status == "healthy", sample]], S5[status == "healthy", age], method = "spearman"))
tst(test = "info_GSE185381_internal_age", cohort = "GSE185381", gene_set = "primary", statistic = conc34, threshold = "informational", null_P5 = NA, null_P95 = NA, n = 25L, pass = NA,
    detail = sprintf("sign concordance of the 370 genes between %d AML < 60 and %d healthy 39-53: %.3f; rho(score, age) within AML %.3f (p %.3f), within healthy %.3f (p %.3f); the four AML < 60 score %s vs healthy range %.2f to %.2f", length(y4), length(h4), conc34, rA$estimate, rA$p.value, rH$estimate, rH$p.value, paste(sprintf("%.2f", sc_in[y4]), collapse = "/"), min(sc_in[st == "healthy"]), max(sc_in[st == "healthy"])))

## -- 14. scores tables, tests, reading, GATES.md --------------------------------------------------------------------------------------
for (ds in COHORTS) { x <- CO[[ds]]; dn <- x$donors; Jx <- J[[ds]]; sJ <- score(Jx$ZE, UP, DN); sJf <- score(Jx$ZE, UPf, DNf)
  SC <- dn[, .(cohort = ds, level = "donor", id = donor, age, sex, grp, n_libraries = n_lib, n_kept_T = n_T, eligible, excluded_K3, score_J_primary = sJ[donor], score_J_full = sJf[donor])]
  am <- donor_amb(x); SC[, `:=`(hb_share_nonery = am$hb_share_nonery[match(id, am$donor)], gran_share_lymph = am$gran_share_lymph[match(id, am$donor)], dbl_mean = am$dbl_mean[match(id, am$donor)])]
  if (ds == DECISIVE) { SC[, score_W_primary := scW[id]]; SC <- rbind(SC, data.table(cohort = ds, level = "library", id = oeL$sample, age = if (SMOKE) NA_real_ else oeL$age_lib, sex = if (SMOKE) NA_character_ else oeL$sex_lib, grp = NA, n_libraries = 1L, n_kept_T = oeL$n_T, eligible = NA, score_J_primary = slib[oeL$sample], score_J_full = NA, score_W_primary = NA), fill = TRUE) }
  fwrite(SC, out(paste0("scores_", ds))) }
TT <- rbindlist(TESTS, fill = TRUE); fwrite(TT, out("tests"))
saveRDS(list(A1 = nulA1, A2b = nulA2, W = nulW, boot = boot, loo = loo, scramble = scr, planted_A2b = nul_pl), file.path(LARGE1_DIR, "ws_pseudobulk", paste0("ws_s9_null_draws", SFX, ".rds")))
# reading (section 5): the gate, then the wording clauses
gate <- pass_A1
bound_txt <- if (CI[2] > LINE) sprintf("below one third at the point estimate; upper bound %.3f (worded as a bound)", CI[2]) else "below one third"
verdict <- if (loo_flip) "GATE_S9 not decisive: a leave-one-donor-out refit moves f_age to >= 1/3 (no claim wording issued; discussion)" else if (gate) sprintf(paste(
  "GATE_S9 PASS: In one healthy marrow cohort (Oetjen 2018: %d of 20 donors, %d-%d y, 10x 3' v2, one protocol, part of the BoneMarrowMap reference), the T imprint's age shift over the GSE185381 age gap is %.3f (95%% CI %.3f to %.3f) of the cross-fitted AML-healthy gap, %s.",
  "The age alternative to the direction is bounded, not excluded. Direction versus healthy remains shown in the defining cohort only (no independent AML + healthy cohort exists in this work); the 32 ambient-prone genes stay excluded from every statement.",
  "The S5 label 'AML-or-age' is narrowed to 'AML, with an age share bounded below one third in healthy marrow (one protocol)'."), nrow(oeD), min(oeD$age), max(oeD$age), f_age, CI[1], CI[2], bound_txt) else sprintf("GATE_S9 FAIL: f_age %.3f >= 1/3; the claim stays AML-or-age with the age share quantified (stop, P5)", f_age)
a3_clause <- function(r) {   # section 5: one sentence per A3 cohort; descriptive cohorts are given with their n and never enter the gate
  if (!r$coverage_K4_pass) sprintf("%s: primary-gene coverage of U_J below 0.80 (K4); its A3 is recorded, not worded", r$cohort)
  else if (r$descriptive_S_E) sprintf("%s is descriptive (n young %d, n old %d): A_E/(D_res/r) %.3f, %d of %d older donors above the AML-healthy midpoint after study correction%s", r$cohort, r$n_young, r$n_old, r$A_E_over_Dg, r$n_old_above_midpoint, r$n_old, if (r$age_effect_visible_in_labels) "" else "; age effect not visible in labels (K1)")
  else if (!r$interpretable) sprintf("%s's study offset or spread is outside the interpretable range (S_E %.3f, SD ratio %s); its at-most-one rule is reported, not worded", r$cohort, r$S_E, f3(r$sd_ratio_young_vs_G_healthy))
  else if (r$direction_holds) sprintf("older healthy donors in %s do not score like AML (A_E/(D_res/r) %.3f, n young %d, n old %d)", r$cohort, r$A_E_over_Dg, r$n_young, r$n_old)
  else sprintf("%d of %d older %s donors score above the AML-healthy midpoint after study correction (contradiction recorded)", r$n_old_above_midpoint, r$n_old, r$cohort) }
clauses <- c(sprintf("f_age %.3f (f_slope %.3f, f_contrast %.3f; 95%% donor-bootstrap CI %.3f to %.3f, conditional on the yardstick); r %.2f; LOO range %.3f to %.3f", f_age, f_slope, f_contrast, CI[1], CI[2], R_DEFL, min(loo), max(loo)),   # reported under every outcome (section 3.1)
             if (gate && !loo_flip) (if (mod_A2b) sprintf("A2b: imprint genes' age effects not more aligned with their AML effects than the random P95 (b/lambda %.3f, P95 %.3f)", al["b_over_lambda"], p95_al) else sprintf("A2b: an age alignment beyond the genome-wide rate is visible (b/lambda %.3f, P95 %.3f)", al["b_over_lambda"], p95_al)),
             if (gate && !loo_flip && f_slope > quantile(nulA1[, "slope"], .95)) sprintf("a weak age association exists (f_slope %.3f above the random P95 %.3f)", f_slope, quantile(nulA1[, "slope"], .95)),
             if (gate && !loo_flip) vapply(seq_len(nrow(DIRT)), function(i) a3_clause(DIRT[i]), character(1)),
             if (gate && !loo_flip && f_age < 0) "age opposes the direction (negative f_age)")
ts("    ", verdict); for (cl in clauses) ts("      - ", cl)
if (!SMOKE) { fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x)); esc <- function(x) gsub("|", "\\|", x, fixed = TRUE)   # a pipe inside a cell (|S_E|, |rho|) would split the markdown table
  block <- c("<!-- A1-S9 start -->", sprintf("## T imprint, age/direction control on older healthy donors (S9, prereg-s9-v1) -- %s", format(Sys.Date())), "",
    sprintf("Scripts: scripts/ws/33_s9_age_direction.R. Primary set = 370 genes. Decisive cohort GSE120221 (Oetjen, %d of 20 donors eligible, %d-%d y); GSE227690 and GSE276360 descriptive.", nrow(oeD), min(oeD$age), max(oeD$age)),
    sprintf("Yardstick: gap_in %.4f (%.2f healthy-SD; registered %.3f = 3.87 healthy-SD, computed on 10-23 without the blacklist step; K6 |diff| %.4f), cross-fitted gap_oof %.3f, r %.2f; D_res_J (Oetjen) %.3f; gate denominator D_res_J / r %.3f. Delta_age 31.5 y.", GAP_IN, GAP_IN / SD_H, GAP_IN_REG, abs(GAP_IN - GAP_IN_REG), GAP_OOF, R_DEFL, D_RES[[DECISIVE]], DG), "",
    "| test | cohort | gene set | statistic | threshold | random P5 | random P95 | n | pass | detail |", "|---|---|---|---|---|---|---|---|---|---|",
    TT[, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |", test, cohort, gene_set, fmt(statistic), esc(threshold), fmt(null_P5), fmt(null_P95), n, pass, esc(detail))], "",
    paste0("Known answers: ", paste(KAT[, sprintf("%s %s = %s (%s)", id, statistic, value, ifelse(is.na(pass), "reported", ifelse(pass, "PASS", "FAIL")))], collapse = "; ")), "",
    sprintf("**Verdict: %s**", verdict), paste0("- ", clauses), "<!-- A1-S9 end -->")
  old <- readLines(GATES); s0 <- grep("<!-- A1-S9 start -->", old, fixed = TRUE); e0 <- grep("<!-- A1-S9 end -->", old, fixed = TRUE)
  writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES); ts("    GATES.md A1-S9 block written") }
C <- rbindlist(CHK); fwrite(C, out("checks"))
ts(sprintf("[done] %d checks (%d FAIL); %d known-answer rows (%d FAIL); %d tests", nrow(C), sum(!C$pass), nrow(KAT), sum(KAT$pass == FALSE, na.rm = TRUE), nrow(TT)))
