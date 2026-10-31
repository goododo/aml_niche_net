#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds ($raw$T, $dx$T, $raw$B, $raw$NK, $cells); results/tables/ws/
#          ws_s5_scores.csv (the 25 T samples and the B/NK designs), ws_s5_de_T.csv (main fit), ws_s5_imprint_genes.csv,
#          ws_s8_t_imprint_ambient_prone.csv, ws_s4_pseudobulk_meta_GSE185381.csv, ws_a1_g02_{oof_scores,perm_labels}.csv,
#          decontx/GSE185381/{ws_s1_ambient_blacklist.csv, GSM*__percell.csv.gz}; gencode_GRCh38_gene_order.txt;
#          04_annotated objects of the 25 samples (percent.mt of kept T cells; Control0005 channel counts);
#          GEO RAW feature references (hashing); reference/savage2021_processing_delay/; outputs of scripts/ws/38
# OUTPUT : results/tables/ws/ws_a1_batch_*.csv (covariates, checks, arms, d1, d8, d9, d10, d11, results, summary);
#          /LARGE1/.../ws_a1_batch/ (draws, jobs, task results); GATES.md block A1-BATCH
# what it does: PREREGISTRATION_A1_batch.md (tag prereg-a1-batch-v1), nothing else.
#   --prep      tables, known answers K0-K4, every reference drawn once and written before use (no arm result is printed)
#   --arms      one fit per arm (G1, W1-W4, D1-D11); writes arm values, no verdict
#   --task=N    slice N of the reference fits (SLURM array; jobs in ws_a1_batch/jobs.csv)
#   --combine   percentiles, verdicts, GATE_A1B, GATES.md block
#   --smoke     with any mode: synthetic counts with real metadata, few jobs; every output ends in _smoke
# This script deletes nothing, does not modify its inputs, never runs or sources scripts 02, 22, 26, 27, 33, 35, 36 or
#   12_pseudobulk_de/*, and opens no Chen2023, paediatric, post-treatment or relapse expression.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(here) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODE <- if ("--prep" %in% .a) "prep" else if ("--arms" %in% .a) "arms" else if ("--combine" %in% .a) "combine" else
        if (any(grepl("^--task=", .a))) "task" else stop("give --prep, --arms, --task=N or --combine")
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
TABWS <- file.path(TAB_DIR, "ws"); BDIR <- file.path(LARGE1_DIR, "ws_a1_batch"); TDIR <- file.path(BDIR, "tasks")
dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
## write guard: the only two places this script writes, never overwriting a real (non-smoke) file
wtab <- function(x) { f <- file.path(TABWS, paste0("ws_a1_batch_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
wb   <- function(x) { f <- file.path(BDIR, sub("(\\.[a-z]+(\\.gz)?)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
rb   <- function(x) file.path(BDIR, sub("(\\.[a-z]+(\\.gz)?)$", paste0(SFX, "\\1"), x))   # read path of a file written by wb()
CHK <- list(); chk <- function(name, expected, got, tol = NULL) {
  ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))   # an empty result never passes
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(format(expected), collapse = ";"), got = paste(format(got), collapse = ";"), pass = ok)
  message(sprintf("  %-76s expect %-12s got %-12s %s", name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (ok) "PASS" else "** FAIL **")) }

FDR_CUT <- 0.10; LFC_CUT <- 0.5; PHI <- 0.40; PHI_RANGE <- c(0.29, 0.51)
N_REF <- 1000L; N_RESERVE <- 20L; N_D6 <- 500L; N_TASKS <- 100L
SEEDS <- c(g2 = 6L, w_period = 11L, w_depth = 12L, w_mito = 13L, w_tdbl = 14L, d6_year = 15L, d6_sub = 16L,
           d8_sets = 21L, d9_sets = 22L, d10_B = 23L, d10_NK = 24L, smoke = 99L)          # SEED + offset, one stream per reference
GAP_OOF <- fread(file.path(TABWS, "ws_a1_g02_summary_replay.csv"))[quantity == "R2_gap_oof", value]   # 1.378164 (A1-G02 R2)
T15 <- paste0("Control", 1:5); P02 <- paste0("Control", 1:4); P01 <- c("AML0310", "AML0612", "AML3133", "AML3762"); P15 <- c("AML1371", "AML4340")

## -- 1. inputs, blacklist, DE engine (copied from 35_a1_g02_replay.R; inputs passed explicitly) ----------------------
PB  <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
BL  <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO  <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
fit_de <- function(counts, d, terms, coef = "statusAML", group = d$status) {   # 27's edgeR branch; sex dropped inside any one-sex fit
  if ("sex" %in% terms && uniqueN(d$sex) < 2) terms <- setdiff(terms, "sex")
  keep_g <- !(rownames(counts) %in% BLACK)
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = group)
  keep <- filterByExpr(y, group = group) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(nrow(X) == nrow(d), qr(X)$rank == ncol(X))
  if (nrow(X) - ncol(X) < 2) stop("residual df < 2: not estimable")      # edgeR segfaults (uncatchable) at 1 residual df
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
  r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, PValue, FDR)]
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]; r }
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")
logz <- function(M) {   # 33's logz: blacklist; CPM > 1 in >= half; TMM; log CPM with prior 1 (label-free)
  M <- M[!rownames(M) %in% BLACK, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1) }

S5ALL <- fread(file.path(TABWS, "ws_s5_scores.csv"))
S5  <- S5ALL[compartment == "T"]                                                  # file order = 27's design order (AML first)
D5  <- S5[, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
dsg <- function(cmp) S5ALL[compartment == cmp, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
MAINFIT <- fread(file.path(TABWS, "ws_s5_de_T.csv"))                               # the S5 main fit (frozen logFC)
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]
stopifnot(nrow(D5) == 25L, nrow(F402) == 402L, nrow(P370) == 370L)
M4  <- fread(file.path(TABWS, "ws_s4_pseudobulk_meta_GSE185381.csv"))[group4 == "T" & sample %in% D5$sample]
OOFT <- fread(file.path(TABWS, "ws_a1_g02_oof_scores.csv"))

## R: retention of the frozen direction (section 3); computed only with >= 0.8 of the set tested
R_of <- function(r, set = P370, main = MAINFIT, min_cov = 0.8) {
  x <- merge(set, r[, .(gene, lfc = logFC)], by = "gene"); x <- merge(x, main[, .(gene, lfc0 = logFC)], by = "gene")
  cov <- nrow(x) / nrow(set)
  if (cov < min_cov) return(list(R = NA_real_, R_up = NA_real_, R_down = NA_real_, cov = cov, conc = NA_real_))
  den <- function(z) sum(abs(z$lfc0))
  list(R = sum(x$s * x$lfc) / den(x), R_up = sum(x[s > 0, lfc]) / den(x[s > 0]), R_down = -sum(x[s < 0, lfc]) / den(x[s < 0]),
       cov = cov, conc = mean(sign(x$lfc) == x$s)) }
fhat <- function(R, phi = PHI) (1 - R) / (1 - phi)                                  # share of the gap carried by the factor
pct_drawn <- function(arm, ref) (1 + sum(ref < arm) + 0.5 * sum(ref == arm)) / (length(ref) + 1)          # drawn references (arm not among them)
pct_enum  <- function(arm, ref_with_arm) (sum(ref_with_arm < arm) + 0.5 * (sum(ref_with_arm == arm) - 1) + 1) / length(ref_with_arm)   # enumerated, arm included

## held-out scoring and cross-fitted gap (35's heldout_scores and R2 loop, with counts, log-CPM and design passed in)
heldout <- function(counts, L, D, h, terms = MAIN) {
  r <- fit_de(counts, D[-h], terms)
  fu <- setdiff(r[hit == TRUE & logFC > 0, gene], AMB); fd <- setdiff(r[hit == TRUE & logFC < 0, gene], AMB)
  Ltr <- L[, D$sample[-h], drop = FALSE]; mu <- rowMeans(Ltr); sg <- apply(Ltr, 1, sd)
  sc <- vapply(h, function(j) score_col(L[, D$sample[j]], mu, sg, fu, fd), 1)
  list(scores = sc, n_genes = length(fu) + length(fd)) }
score_col <- function(x, mu, sg, fu, fd) { z <- (x - mu) / sg; z <- z[is.finite(z)]
  (if (length(intersect(fu, names(z)))) mean(z[intersect(fu, names(z))]) else 0) - (if (length(intersect(fd, names(z)))) mean(z[intersect(fd, names(z))]) else 0) }
cf_gap <- function(counts, D, terms = MAIN, folds = seq_len(nrow(D)), L = NULL) {
  if (is.null(L)) L <- logz(counts[, D$sample, drop = FALSE])
  oof <- vapply(folds, function(i) tryCatch(heldout(counts, L, D, i, terms)$scores, error = function(e) NA_real_), 1)
  st <- D$status[folds]
  list(gap = mean(oof[st == "AML"], na.rm = TRUE) - mean(oof[st == "healthy"], na.rm = TRUE), oof = oof, n_failed = sum(is.na(oof))) }
make_L27 <- function(M25, U) {   # D1: the 25 columns exactly as logz (R2); the channel units normalised against R2's TMM reference column
  M <- M25[!rownames(M25) %in% BLACK, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L25 <- cpm(y, log = TRUE, prior.count = 1); Mk <- M[keep, , drop = FALSE]
  f75 <- apply(t(t(Mk) / colSums(Mk)), 2, quantile, 0.75); ref <- which.min(abs(f75 - mean(f75)))   # calcNormFactors' default reference
  eff25 <- y$samples$lib.size * y$samples$norm.factors
  man <- function(X, eff) { pc <- eff / mean(eff25); log2(t((t(X) + pc) / (eff + 2 * pc)) * 1e6) }  # edgeR's log-CPM with R2's mean library size
  Uk <- as.matrix(U[rownames(Mk), , drop = FALSE])
  fu <- vapply(seq_len(ncol(Uk)), function(j) { z <- calcNormFactors(DGEList(cbind(Mk[, ref], Uk[, j])), refColumn = 1)$samples$norm.factors; z[2] / z[1] }, 1)
  LU <- man(Uk, colSums(Uk) * fu * y$samples$norm.factors[ref]); colnames(LU) <- colnames(U)
  list(L = cbind(L25, LU), check = man(Mk, eff25) - L25) }
vfactor <- function(d) { X <- model.matrix(fml(MAIN), as.data.frame(d)); solve(crossprod(X))["statusAML", "statusAML"] }  # design only

## -- 2. counts (real, or synthetic for smoke runs: real metadata, no real expression) ------------------------------------
G_RAW <- PB$raw$T[, D5$sample]; stopifnot(ncol(G_RAW) == 25L)
synth <- function(M, aml, seed, plant = TRUE) {   # NB counts around one pooled profile, each column's real total; S5 logFC planted on P370 in AML columns
  set.seed(SEEDS[["smoke"]] + seed); p <- rowSums(M) / sum(M); mu <- outer(p, colSums(M))
  if (plant && length(aml)) { lf <- setNames(MAINFIT$logFC, MAINFIT$gene)[intersect(P370$gene, rownames(M))]; mu[names(lf), aml] <- mu[names(lf), aml] * 2^lf }
  matrix(rnbinom(length(mu), size = 10, mu = mu), nrow(M), dimnames = dimnames(M)) }
AMLS <- D5[status == "AML", sample]; HEALTHY <- D5[status == "healthy", sample]
if (SMOKE) G_RAW <- synth(G_RAW, AMLS, 1)

## -- 3. per-sample tables: dates, period, hashing, depth, mito, Control0005 channel units, t_dbl ------------------------
KTC <- PB$cells[kept == TRUE & group4 == "T" & sample %in% D5$sample]
KTC[, `:=`(library = sub("_.*$", "", cell), date = regmatches(cell, regexpr("20[12][0-9]-[01][0-9]-[0-3][0-9]", cell)))]
DATES <- KTC[, .(dates = paste(sort(unique(date)), collapse = ";"), libraries = paste(sort(unique(library)), collapse = ";")), by = sample]
RAWD <- file.path(RAW_DIR, "public", "GEO", "GSE185381_RAW")
hashed_lib <- function(gsm) { f <- list.files(RAWD, pattern = paste0("^", gsm, "_.*feature_reference.csv.gz$"), full.names = TRUE)
  if (!length(f)) return(FALSE); x <- fread(f[1]); any(grepl("(?i)(^hto|hashtag|^healthy_|^control_[0-9]|totalseqa|^cntrl|^aml[0-9]+-)", x$id, perl = TRUE)) }
LIBH <- KTC[, .(hashed = hashed_lib(library[1])), by = library]
COV <- merge(D5[, .(sample, status, sex)], DATES, by = "sample")
COV[, hashed := vapply(strsplit(libraries, ";"), function(l) all(LIBH[library %in% l, hashed]), TRUE)]
COV[, period := as.integer(vapply(strsplit(dates, ";"), function(x) any(x >= "2020-07-29"), TRUE))]   # late = any library dated 2020-07-29 or later
COV[, year2019 := as.integer(vapply(strsplit(dates, ";"), function(x) all(substr(x, 1, 4) == "2019"), TRUE))]
COV[, depth := log10(M4$median_umi[match(sample, M4$sample)])]
COV[, sex_inferred := M4$sex_inferred[match(sample, M4$sample)]]
COV[, library_id := M4$library_id[match(sample, M4$sample)]]
AN_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", "GSE185381")
if (MODE %in% c("prep", "arms")) {   # percent.mt of kept T cells; Control0005 channel units (the only per-cell reads of this script)
  suppressPackageStartupMessages(library(Seurat))
  mito <- vapply(D5$sample, function(s) { o <- readRDS(file.path(AN_DIR, paste0(s, ".rds"))); k <- intersect(KTC[sample == s, cell], colnames(o))
    stopifnot(length(k) == KTC[sample == s, .N]); median(o@meta.data[k, "percent.mt"]) }, 1)
  COV[, mito := mito[sample]]
  o5 <- readRDS(file.path(AN_DIR, "Control0005.rds")); k5 <- KTC[sample == "Control0005", cell]
  C5 <- SeuratObject::LayerData(o5[["RNA"]], layer = "counts")[, k5]
  UNITS <- cbind(GSM5613781 = Matrix::rowSums(C5[, startsWith(k5, "GSM5613781")]), GSM5613787 = Matrix::rowSums(C5[, startsWith(k5, "GSM5613787")]))
  rownames(UNITS) <- rownames(C5); UNITS <- UNITS[rownames(G_RAW), ]
  if (SMOKE) { set.seed(SEEDS[["smoke"]] + 2); u1 <- rbinom(nrow(G_RAW), G_RAW[, "Control0005"], 523 / 934)
    UNITS <- cbind(GSM5613781 = u1, GSM5613787 = G_RAW[, "Control0005"] - u1); rownames(UNITS) <- rownames(G_RAW) }
  rm(o5, C5); gc(FALSE)
}
W4F <- file.path(BDIR, "w4_doublet_scores.csv.gz")
if (file.exists(W4F)) { W4 <- fread(W4F); COV[, t_dbl := W4[kept_T == TRUE, mean(score), by = sample][match(COV$sample, sample), V1]]
  COV[, dbl_mean_registered := W4[, mean(score), by = sample][match(COV$sample, sample), V1]]
} else if (SMOKE) { set.seed(SEEDS[["smoke"]] + 3); COV[, t_dbl := runif(.N)]; COV[, dbl_mean_registered := runif(.N)]
} else if (MODE != "combine") stop("W4 scores missing: run scripts/ws/38_a1_batch_cells.R first (", W4F, ")")
COV <- COV[match(D5$sample, sample)]
D <- merge(D5, COV[, .(sample, period, year2019, depth, mito = if ("mito" %in% names(COV)) mito else NA_real_, t_dbl, hashed)], by = "sample", sort = FALSE)
D <- D[match(D5$sample, sample)]; stopifnot(identical(D$sample, D5$sample))
WCOV <- c(W1 = "period", W2 = "depth", W3 = "mito", W4 = "t_dbl")

## -- 4. MODE prep: tables, known answers, every reference drawn once ----------------------------------------------------
if (MODE == "prep") {
  ts("[prep] K0 tables")
  EXP_DATES <- c(AML0310 = "2019-07-01", AML0612 = "2019-07-01", AML3133 = "2019-07-01", AML3762 = "2019-07-01", Control1 = "2019-08-01",
    Control2 = "2019-08-01", Control3 = "2019-08-01", Control4 = "2019-08-01", Control5 = "2019-08-29", AML2451 = "2019-10-25",
    Control0004 = "2019-10-25", Control0058 = "2019-11-04", AML0048 = "2020-03-18", AML3266 = "2020-03-18", Control0082 = "2020-03-18",
    AML2123 = "2020-07-29", Control4003 = "2020-08-26", AML1371 = "2020-09-11", AML4340 = "2020-09-11", AML4897 = "2020-10-06",
    Control0005 = "2020-10-06;2021-03-09", AML0114 = "2021-03-09", AML0693 = "2021-03-09", AML3730 = "2021-03-09", AML3948 = "2021-03-09")
  chk("K0 file-name dates of the 25 samples (section 2.2)", TRUE, identical(COV$dates, unname(EXP_DATES[COV$sample])))
  chk("K0 period: early 7 AML / 8 healthy, late 8 / 2", c(7L, 8L, 8L, 2L),
      c(COV[status == "AML" & period == 0, .N], COV[status == "healthy" & period == 0, .N], COV[status == "AML" & period == 1, .N], COV[status == "healthy" & period == 1, .N]))
  chk("K0 calendar year 2019: 5 AML / 7 healthy", c(5L, 7L), c(COV[status == "AML" & year2019 == 1, .N], COV[status == "healthy" & year2019 == 1, .N]))
  chk("K0 unhashed samples (feature references)", c("AML0693", "AML3266", "AML3730", "AML3948"), sort(COV[hashed == FALSE, sample]))
  nlib <- function(s) COV[library_id == COV[sample == s, library_id], .N]
  chk("K0 library blocks P01, P02, P15, P05: samples of the 25 per block", c(4L, 4L, 2L, 2L), c(nlib("AML0310"), nlib("Control1"), nlib("AML1371"), nlib("AML2451")))
  # block-respecting label space (G2): P01, P02, P15 whole; everything else (incl. P05's two samples) single
  units <- c(list(P01 = P01, P02 = P02, P15 = P15), as.list(setNames(setdiff(D$sample, c(P01, P02, P15)), setdiff(D$sample, c(P01, P02, P15)))))
  usz <- lengths(units); blk <- c("P01", "P02", "P15"); sing <- setdiff(names(units), blk)
  LAB <- rbindlist(lapply(0:7, function(b) { inc <- blk[bitwAnd(b, c(1L, 2L, 4L)) > 0]; need <- 15L - sum(usz[inc])
    if (need < 0 || need > length(sing)) return(NULL)
    cmb <- combn(sing, need, simplify = FALSE)
    data.table(labels = vapply(cmb, function(s) paste(as.integer(D$sample %in% unlist(units[c(inc, s)])), collapse = ""), "")) }))
  chk("K0 block-respecting labellings with 15 AML (G2)", 22284L, nrow(LAB))
  true_lab <- paste(as.integer(D$status == "AML"), collapse = "")
  set.seed(SEED + SEEDS[["g2"]]); G2L <- LAB[labels != true_lab][sample.int(.N, N_REF + N_RESERVE)]
  G2L[, id := .I]
  # design-only drop references (G1, D3 B2c, D7 P01)
  v_full <- vfactor(D)
  drops <- function(pool, k) { cb <- combn(pool, k, simplify = FALSE)
    data.table(id = seq_along(cb), drop = vapply(cb, paste, "", collapse = ";"), v = vapply(cb, function(x) vfactor(D[!sample %in% x]), 1)) }
  DR5 <- drops(HEALTHY, 5L)
  DR5[, real := vapply(strsplit(drop, ";"), function(x) setequal(x, T15), TRUE)]
  DR5[, `:=`(dv = v - v_full, matched = abs(v / v[real] - 1) <= 0.15)]
  DR4 <- drops(HEALTHY, 4L); DR4[, real := vapply(strsplit(drop, ";"), function(x) setequal(x, P02), TRUE)]; DR4[, matched := abs(v / v[real] - 1) <= 0.15]
  DA4 <- drops(AMLS, 4L); DA4[, real := vapply(strsplit(drop, ";"), function(x) setequal(x, P01), TRUE)]; DA4[, matched := abs(v / v[real] - 1) <= 0.15]
  chk("K0 G1 reference: matched drops of 5 healthy (incl. the real one)", 40L, DR5[matched == TRUE, .N])
  chk("K0 D3 reference: matched drops of 4 healthy (incl. P02)", 13L, DR4[matched == TRUE, .N])
  ts(sprintf("    G1 real drop variance factor %.3f = %.2f x full (rank %d of 252); D7 P01 matched drops: %d", DR5[real == TRUE, v], DR5[real == TRUE, v] / v_full,
             DR5[, sum(v >= v[real])], DA4[matched == TRUE, .N]))
  # covariate table with correlation with status and t cap (a gene equal to the pattern; design only)
  tcap <- function(x) { d <- copy(D); d[, y := x]; summary(lm(fml(MAIN) |> update(y ~ .), data = d))$coefficients["statusAML", "t value"] }
  CT <- data.table(covariate = c(WCOV, D6 = "year2019", G1 = "Control1-5"),
                   cor_status = c(vapply(c(WCOV, "year2019"), function(v) cor(D[[v]], as.integer(D$status == "AML")), 1), cor(as.integer(D$sample %in% T15), as.integer(D$status == "AML"))),
                   t_cap = c(vapply(c(WCOV, "year2019"), function(v) tcap(D[[v]]), 1), tcap(as.integer(D$sample %in% T15))))
  print(CT)
  # residual-permutation draws for covariate references (W1-W4, D6 year): fitted part kept, residual permuted
  WD <- lapply(c(WCOV, D6 = "year2019"), function(v) { s <- SEEDS[[c(period = "w_period", depth = "w_depth", mito = "w_mito", t_dbl = "w_tdbl", year2019 = "d6_year")[[v]]]]
    m <- lm(update(fml(MAIN), cv ~ .), data = cbind(D, cv = D[[v]])); f <- fitted(m); r <- resid(m)
    set.seed(SEED + s); Z <- replicate(N_REF + N_RESERVE, { p <- sample.int(25); while (identical(p, 1:25)) p <- sample.int(25); f + r[p] })
    attr(Z, "fitted") <- f; attr(Z, "resid") <- r; Z })
  chk("K3(e) every covariate draw keeps its fitted part and permutes its residual", TRUE,
      all(vapply(WD, function(Z) max(abs(apply(Z - attr(Z, "fitted"), 2, sort) - sort(attr(Z, "resid")))) < 1e-9, TRUE)))
  # D6 early-stratum subsets: 7 of 15 AML and 8 of 10 healthy, P01 and P02 each whole or absent
  amlo <- setdiff(AMLS, P01); heo <- setdiff(HEALTHY, P02)
  SA <- c(lapply(combn(amlo, 3, simplify = FALSE), function(x) c(P01, x)), combn(amlo, 7, simplify = FALSE))
  SH <- lapply(combn(heo, 4, simplify = FALSE), function(x) c(P02, x))
  EARLY <- D[period == 0, sample]
  SUB <- CJ(a = seq_along(SA), h = seq_along(SH))[, subset := mapply(function(i, j) paste(sort(c(SA[[i]], SH[[j]])), collapse = ";"), a, h)]
  SUB <- SUB[subset != paste(sort(EARLY), collapse = ";")]
  ts(sprintf("    D6 subset space %d (+ the real stratum); drawing %d", nrow(SUB), N_D6))
  set.seed(SEED + SEEDS[["d6_sub"]]); SUB <- SUB[sample.int(.N, N_D6)][, id := .I]
  # D2 splits of the 10 healthy donors into 5 + 5 (all 252, incl. the real one)
  SPL <- data.table(id = 1:252, group1 = vapply(combn(HEALTHY, 5, simplify = FALSE), paste, "", collapse = ";"))
  # random gene sets matched on label-free expression (D8 overlap, D9 ratio: T; D10: T and the compartment)
  L <- logz(G_RAW); ex <- rowMeans(L)
  match_sets <- function(target, universe, e1, e2 = NULL, n, seed) {   # for each target gene, a random non-target gene from the same expression bin
    u <- setdiff(universe, target); b1 <- cut(rank(e1[c(target, u)]), 10, labels = FALSE); names(b1) <- c(target, u)
    bin <- if (is.null(e2)) b1 else { b2 <- cut(rank(e2[c(target, u)]), 5, labels = FALSE); paste(b1, b2) }; names(bin) <- c(target, u)
    pools <- split(u, bin[u]); set.seed(SEED + seed)
    replicate(n, vapply(target, function(g) { pool <- pools[[as.character(bin[[g]])]]; if (!length(pool)) NA_character_ else pool[sample.int(length(pool), 1)] }, ""), simplify = FALSE) }
  RS <- list(d8_up = match_sets(P370[s > 0, gene], rownames(L), ex, n = N_REF, seed = SEEDS[["d8_sets"]]),
             d8_down = match_sets(P370[s < 0, gene], rownames(L), ex, n = N_REF, seed = SEEDS[["d8_sets"]] + 100L))
  RS$d9_up <- RS$d8_up; RS$d9_down <- RS$d8_down                         # D9 ratio uses the same expression-matched sets
  for (cmp in c("B", "NK")) { dc <- dsg(cmp); Mc <- PB$raw[[cmp]][, dc$sample]; if (SMOKE) Mc <- synth(Mc, dc[status == "AML", sample], 4 + (cmp == "NK"))
    Lc <- logz(Mc); uni <- intersect(intersect(rownames(L), rownames(Lc)), MAINFIT$gene); tg <- intersect(P370$gene, uni)
    RS[[paste0("d10_", cmp)]] <- match_sets(tg, uni, ex, rowMeans(Lc), n = N_REF, seed = SEEDS[[paste0("d10_", cmp)]]) }
  # known answers that need fits (no arm evaluated)
  ts("[prep] known answers with fits")
  r0 <- fit_de(G_RAW, D, MAIN)
  if (!SMOKE) {
    chk("K1 the MAIN refit reproduces the 402 S5 T imprint genes", TRUE, setequal(r0[hit == TRUE, gene], F402$gene))
    chk("K1 R(MAIN refit) = 1", 1, R_of(r0)$R, tol = 1e-6)
    chk("K1 md5 of ws_a1_g02_perm_labels.csv", "b2f6e70b7e16f4ea374dce14ef6684e7", unname(tools::md5sum(file.path(TABWS, "ws_a1_g02_perm_labels.csv"))))
    PL <- fread(file.path(TABWS, "ws_a1_g02_perm_labels.csv"), colClasses = c(labels = "character"))
    chk("K1 D2 label file: 1,000 rows, 25 characters, 15 ones", c(1000L, 25L, 15L), c(nrow(PL), unique(nchar(PL$labels)), unique(vapply(strsplit(PL$labels, ""), function(x) sum(x == "1"), 1L))))
    chk("K2 the two Control0005 channel units sum to its $raw$T column", 0, max(abs(rowSums(UNITS) - PB$raw$T[rownames(UNITS), "Control0005"])))
    k3a <- cf_gap(G_RAW, D)
    chk("K3(a)/(d) cross-fit with nothing removed and the true labels reproduces gap_oof", GAP_OOF, k3a$gap, tol = 0.001)
    chk("K3(b) the covariate path with no covariate gives R = 1", 1, R_of(fit_de(G_RAW, D, MAIN))$R, tol = 1e-6)
    dm <- merge(D, S5[, .(sample, multiplexed)], by = "sample", sort = FALSE)[match(D$sample, sample)]
    rm_ <- fit_de(G_RAW, dm, c(MAIN, "multiplexed"))
    chk("K3(c) MAIN + multiplexed reproduces S5's arm (405 hits, 308 in the 402)", c(405L, 308L), c(rm_[hit == TRUE, .N], length(intersect(rm_[hit == TRUE, gene], F402$gene))))
  }
  # K4: D1 arithmetic on synthetic scores (no real D1 score is computed here)
  DAYS <- list(list(day = "2019-10-25", aml = "AML2451", h = "Control0004"), list(day = "2020-03-18", aml = c("AML0048", "AML3266"), h = "Control0082"),
               list(day = "2020-10-06", aml = "AML4897", h = "GSM5613781"))
  w_of <- function(s, hsel = vapply(DAYS, `[[`, "", "h")) {   # W for one arrangement: hsel[d] is the unit given the healthy role in day d
    dl <- unlist(lapply(seq_along(DAYS), function(d) { u <- c(DAYS[[d]]$aml, DAYS[[d]]$h); a <- setdiff(u, hsel[d]); s[a] - s[hsel[d]] })); mean(dl) }
  arrangements <- as.matrix(do.call(CJ, lapply(DAYS, function(x) c(x$aml, x$h))))
  set.seed(SEED); s0 <- setNames(rnorm(7), unlist(lapply(DAYS, function(x) c(x$aml, x$h))))
  s1 <- s0; s1[c("AML0048", "AML3266", "Control0082")] <- s1[c("AML0048", "AML3266", "Control0082")] + 3
  s2 <- s0; s2[c("AML2451", "AML0048", "AML3266", "AML4897")] <- s2[c("AML2451", "AML0048", "AML3266", "AML4897")] + 0.5
  chk("K4(a) a constant added to one day leaves W unchanged", 0, w_of(s1) - w_of(s0), tol = 1e-9)
  chk("K4(b) +0.5 on the 4 AML units raises W by 0.5", 0.5, w_of(s2) - w_of(s0), tol = 1e-9)
  chk("K4(c) the 12 arrangement values of W sum to 0", c(12, 0), c(nrow(arrangements), sum(apply(arrangements, 1, function(h) w_of(s0, h)))), tol = 1e-9)
  oh <- setNames(rep(0, 7), names(s0)); oh["Control0082"] <- 1
  chk("K4(d) one-hot on Control0082: W = (0 - 1 - 1 + 0)/4", -0.5, w_of(oh), tol = 1e-9)
  # L27: the 25 columns exactly as R2; the two channel units normalised against R2's reference column
  L27 <- make_L27(G_RAW, UNITS)
  chk("L27: the manual log-CPM formula reproduces logz on the 25 columns", 0, max(abs(L27$check)), tol = 1e-8)
  if (!SMOKE) { e <- heldout(G_RAW, L27$L, D, which(D$sample == "AML0310"))$scores
    chk("K4(e) scoring AML0310 alone with D1's normalisation reproduces its out-of-fold score", OOFT[sample == "AML0310", oof_score], e, tol = 0.001) }
  C <- rbindlist(CHK); ok <- all(C$pass)
  fwrite(C, if (ok) wtab("checks_prep") else file.path(TABWS, sprintf("ws_a1_batch_checks_prep_FAILED_%s%s.csv", format(Sys.time(), "%Y%m%d%H%M"), SFX)))
  if (!ok) { ts("[prep] known answers FAILED: nothing else written"); quit(status = 1) }
  fwrite(COV, wtab("covariates")); fwrite(CT, wtab("covariate_tcaps"))
  fwrite(G2L, wb("g2_labels.csv")); fwrite(DR5, wb("drops5_healthy.csv")); fwrite(DR4, wb("drops4_healthy.csv")); fwrite(DA4, wb("drops4_aml.csv"))
  fwrite(SUB, wb("d6_subsets.csv")); fwrite(SPL, wb("d2_splits.csv")); saveRDS(WD, wb("w_draws.rds")); saveRDS(RS, wb("random_sets.rds"))
  # jobs: every reference fit, one row each; distributed round-robin after sorting by cost
  J <- rbindlist(list(
    data.table(type = "g1drop", kid = DR5$id, cost = 7), data.table(type = "d3drop4", kid = DR4[matched == TRUE & real == FALSE, id], cost = 7),
    data.table(type = "d7drop", kid = DA4[matched == TRUE & real == FALSE, id], cost = 7),
    rbindlist(lapply(names(WD), function(w) data.table(type = paste0("wref_", w), kid = seq_len(N_REF + N_RESERVE), cost = 7))),
    data.table(type = "d6sub", kid = SUB$id, cost = 6), data.table(type = "d2split", kid = SPL$id, cost = 3),
    data.table(type = "g2", kid = G2L$id, cost = 175), data.table(type = "d3cf", kid = DR5[matched == TRUE & real == FALSE, id], cost = 140),
    data.table(type = "lpo", kid = seq_len(150), cost = 7),
    data.table(type = "calib", kid = as.vector(outer(c(0, DR5[matched == TRUE, id]), c(1, 2, 3), function(i, k) k * 1000 + i)), cost = 7)))
  if (SMOKE) J <- J[, .SD[1:min(.N, 2)], by = type]
  J <- J[order(-cost)][, `:=`(job = .I, task = ((.I - 1L) %% N_TASKS) + 1L)]
  fwrite(J, wb("jobs.csv"))
  ts(sprintf("[prep] done: %d jobs (%.1f CPU-h), %d tasks; known answers %d/%d PASS", nrow(J), sum(J$cost) / 3600, uniqueN(J$task), sum(C$pass), nrow(C)))
  quit(status = 0)
}

## -- 5. shared by arms, task and combine: draws written by --prep --------------------------------------------------------
stopifnot(file.exists(rb("jobs.csv")), file.exists(wtab_path <- file.path(TABWS, paste0("ws_a1_batch_checks_prep", SFX, ".csv"))))
stopifnot(all(fread(wtab_path)$pass))                                              # --prep's known answers all passed
JOBS <- fread(rb("jobs.csv")); G2L <- fread(rb("g2_labels.csv"), colClasses = c(labels = "character"))
DR5 <- fread(rb("drops5_healthy.csv")); DR4 <- fread(rb("drops4_healthy.csv")); DA4 <- fread(rb("drops4_aml.csv"))
SUB <- fread(rb("d6_subsets.csv")); SPL <- fread(rb("d2_splits.csv")); WD <- readRDS(rb("w_draws.rds"))
L25 <- logz(G_RAW)
LFC0 <- setNames(MAINFIT$logFC, MAINFIT$gene)
splitv <- function(x) strsplit(x, ";")[[1]]
with_cov <- function(d, v) { d <- copy(d); d[, cov := v]; d }
Df <- copy(D); Df[sample == "Control4", sex := factor("female", levels(D$sex))]       # D5: Control4 coded female (expression)
d2_fit <- function(g1) { dh <- D[status == "healthy"]; dh[, group := factor(as.integer(sample %in% g1), c(0, 1))]
  r <- fit_de(G_RAW, dh, c("group", "sex"), coef = "group1", group = dh$group)
  x <- merge(P370, r[, .(gene, lfc = logFC)], by = "gene"); mean(sign(x$lfc) == -x$s) }   # share of P370 genes with Control-group sign opposite to the imprint
plant_counts <- function(k) { Gp <- G_RAW; g <- intersect(P370$gene, rownames(Gp))     # calibration: healthy-ward shift on Control1-5's P370 counts
  Gp[g, T15] <- round(Gp[g, T15] * 2^(-k * LFC0[g])); Gp }
CAL_K <- c(`1` = 0.25, `2` = 0.5, `3` = 1.0)

## -- 6. MODE arms: one fit per arm ---------------------------------------------------------------------------------------
if (MODE == "arms") {
  A <- list(); add <- function(arm, metric, value) A[[length(A) + 1]] <<- data.table(arm = arm, metric = metric, value = as.numeric(value))
  addR <- function(arm, r) { x <- R_of(r); for (m in names(x)) add(arm, m, x[[m]]); add(arm, "R402", R_of(r, F402)$R); add(arm, "fhat", fhat(x$R))
    add(arm, "fhat_lo", fhat(x$R, PHI_RANGE[1])); add(arm, "fhat_hi", fhat(x$R, PHI_RANGE[2])); ts(sprintf("    %-26s fitted", arm)) }
  ts("[arms] G1, W1-W4")
  addR("G1", fit_de(G_RAW, D[!sample %in% T15], MAIN))
  for (w in names(WCOV)) addR(w, fit_de(G_RAW, with_cov(D, D[[WCOV[[w]]]]), c(MAIN, "cov")))
  ts("[arms] D1 capture-day contrast")
  L27 <- make_L27(G_RAW, UNITS)$L
  day_scores <- function(aml, hunit, drop_train) { tr <- D[!sample %in% drop_train]; r <- fit_de(G_RAW, tr, MAIN)
    fu <- setdiff(r[hit == TRUE & logFC > 0, gene], AMB); fd <- setdiff(r[hit == TRUE & logFC < 0, gene], AMB)
    Ltr <- L27[, tr$sample, drop = FALSE]; mu <- rowMeans(Ltr); sg <- apply(Ltr, 1, sd)
    vapply(c(aml, hunit), function(u) score_col(L27[, u], mu, sg, fu, fd), 1) }
  DAYS <- list(list(day = "2019-10-25", aml = "AML2451", h = "Control0004", drop = c("AML2451", "Control0004")),
               list(day = "2020-03-18", aml = c("AML0048", "AML3266"), h = "Control0082", drop = c("AML0048", "AML3266", "Control0082")),
               list(day = "2020-10-06", aml = "AML4897", h = "GSM5613781", drop = c("AML4897", "Control0005")),
               list(day = "2021-03-09 (date not established)", aml = c("AML0114", "AML0693", "AML3730", "AML3948"), h = "GSM5613787",
                    drop = c("AML0114", "AML0693", "AML3730", "AML3948", "Control0005")))
  SD1 <- lapply(DAYS, function(x) day_scores(x$aml, x$h, x$drop))
  D1T <- rbindlist(lapply(seq_along(DAYS), function(d) { s <- SD1[[d]]
    data.table(day = DAYS[[d]]$day, unit = names(s), role = c(rep("AML", length(DAYS[[d]]$aml)), "healthy"), score = s,
               delta = c(s[DAYS[[d]]$aml] - s[DAYS[[d]]$h], NA)) }))
  s3 <- unlist(SD1[1:3]); w_of <- function(s, hsel) mean(unlist(lapply(1:3, function(d) { u <- c(DAYS[[d]]$aml, DAYS[[d]]$h); a <- setdiff(u, hsel[d]); s[a] - s[hsel[d]] })))
  arr <- as.matrix(do.call(CJ, lapply(DAYS[1:3], function(x) c(x$aml, x$h))))
  W <- w_of(s3, vapply(DAYS[1:3], `[[`, "", "h")); Wp <- apply(arr, 1, function(h) w_of(s3, h))
  add("D1", "W", W); add("D1", "W_over_gap", W / GAP_OOF); add("D1", "p_W_12", mean(Wp >= W - 1e-9)); add("D1", "same_channel_delta", D1T[unit == "AML2451", delta])
  add("D1", "W_2021_not_established", D1T[day == DAYS[[4]]$day & role == "AML", mean(delta)])
  fwrite(D1T, wtab("d1"))
  ts("[arms] D2-D7")
  add("D2", "opposite_share", d2_fit(T15))
  add("D3", "CF_without_Control1-5", cf_gap(G_RAW, D[!sample %in% T15])$gap)
  addR("D3_P02", fit_de(G_RAW, D[!sample %in% P02], MAIN))
  DX <- round(PB$dx$T[, D$sample]); if (SMOKE) DX <- synth(DX, AMLS, 6)
  rdx <- fit_de(DX, D, MAIN); addR("D4", rdx)
  for (g in c("GNA15", "HMGB3")) add("D4", paste0("logFC_ratio_", g), rdx[gene == g, logFC] / LFC0[[g]])
  add("D4", "CF_over_gap", cf_gap(DX, D)$gap / GAP_OOF)
  addR("D5_Control4_female", fit_de(G_RAW, Df, MAIN)); add("D5_Control4_female", "CF_over_gap", cf_gap(G_RAW, Df, L = L25)$gap / GAP_OOF)
  addR("D5_males_by_expression", fit_de(G_RAW, D[COV$sex_inferred[match(D$sample, COV$sample)] == "male"], setdiff(MAIN, "sex")))
  addR("D6_year", fit_de(G_RAW, with_cov(D, D$year2019), c(MAIN, "cov")))
  addR("D6_early", fit_de(G_RAW, D[period == 0], MAIN))
  r6 <- tryCatch(fit_de(G_RAW, D[period == 0 & !sample %in% c(P01, P02)], MAIN), error = function(e) NULL)
  if (is.null(r6)) add("D6_early_without_P01_P02", "R", NA) else addR("D6_early_without_P01_P02", r6)
  addR("D7_P01", fit_de(G_RAW, D[!sample %in% P01], MAIN)); addR("D7_P15", fit_de(G_RAW, D[!sample %in% P15], MAIN))
  addR("D7_P05", fit_de(G_RAW, D[!sample %in% c("AML2451", "Control0004")], MAIN)); addR("D7_Control5", fit_de(G_RAW, D[sample != "Control5"], MAIN))
  ts("[arms] D8 handling-artefact and delay scores")
  RS <- readRDS(rb("random_sets.rds"))
  Z <- (L25 - rowMeans(L25)) / apply(L25, 1, sd)
  zs <- function(up, down = character(0)) { u <- setdiff(intersect(up, rownames(Z)), F402$gene); d <- setdiff(intersect(down, rownames(Z)), F402$gene)
    list(score = (if (length(u)) colMeans(Z[u, , drop = FALSE]) else 0) - (if (length(d)) colMeans(Z[d, , drop = FALSE]) else 0), n = length(u) + length(d)) }
  SAV <- "/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv"
  stopifnot(unname(tools::md5sum(SAV)) == "5ee444a9284ae52b7afab1c8aeeced2c"); SV <- fread(SAV)
  SETS <- list(IEG = zs(WS_B0_IEG), IEG_noKLF = zs(setdiff(WS_B0_IEG, c("KLF2", "KLF4", "KLF6"))), Savage = zs(SV[direction == "up", gene], SV[direction == "down", gene]))
  oof <- OOFT[match(D$sample, sample), oof_score]
  D8 <- rbindlist(lapply(names(SETS), function(k) { sc <- SETS[[k]]$score; a <- D$status == "AML"
    rr <- R_of(fit_de(G_RAW, with_cov(D, sc), c(MAIN, "cov")))
    data.table(set = k, n_genes = SETS[[k]]$n, diff_AML_minus_healthy = mean(sc[a]) - mean(sc[!a]), wilcox_p = wilcox.test(sc[a], sc[!a])$p.value,
               rho_oof_AML = cor(sc[a], oof[a], method = "spearman"), rho_oof_healthy = cor(sc[!a], oof[!a], method = "spearman"), R = rr$R, R_up = rr$R_up, R_down = rr$R_down) }))
  ov <- function(a, b) length(intersect(a, b))
  P_UP <- P370[s > 0, gene]; P_DN <- P370[s < 0, gene]; S_UP <- SV[direction == "up", gene]; S_DN <- SV[direction == "down", gene]
  OVL <- rbindlist(lapply(list(c("up", "up"), c("down", "down"), c("up", "down"), c("down", "up")), function(x) {
    tg <- if (x[1] == "up") P_UP else P_DN; sv <- if (x[2] == "up") S_UP else S_DN; rs <- if (x[1] == "up") RS$d8_up else RS$d8_down
    obs <- ov(tg, sv); rnd <- vapply(rs, function(g) ov(g, sv), 1)
    data.table(imprint = x[1], delay = x[2], overlap = obs, random_median = median(rnd), random_P95 = quantile(rnd, 0.95), p_ge = (1 + sum(rnd >= obs)) / (length(rnd) + 1)) }))
  fwrite(D8, wtab("d8")); fwrite(OVL, wtab("d8_overlap")); print(D8); print(OVL)
  ts("[arms] D9 T-ALL soup bound")
  PRF <- fread(file.path(BDIR, paste0("d9_soup_profile", SFX, ".csv.gz"))); LN <- fread(file.path(BDIR, paste0("d9_lanes", SFX, ".csv")))
  CK <- fread(file.path(BDIR, paste0("cells_checks", SFX, ".csv"))); k6b <- CK[grepl("^K6b", check), pass]
  p <- setNames(PRF$counts, PRF$gene)[rownames(G_RAW)]; p[is.na(p)] <- 0
  cs <- p / sum(p) * 1e6; nt <- rowSums(G_RAW[, HEALTHY]); ct <- nt / sum(nt) * 1e6; lr <- log2((cs + 1) / (ct + 1))   # T-ALL / normal T, 1 CPM pseudocount
  rat <- function(tg, rs) { obs <- mean(lr[tg]); rnd <- vapply(rs, function(g) mean(lr[g], na.rm = TRUE), 1); c(obs = obs, rnd_med = median(rnd), p_ge = (1 + sum(rnd >= obs)) / (length(rnd) + 1)) }
  ru <- rat(P_UP, RS$d9_up); rd <- rat(P_DN, RS$d9_down)
  LN[, gsm := sub("_.*$", "", lane)]
  dose <- rbindlist(lapply(T15, function(s) { kt <- KTC[sample == s]; libs <- kt[, .N, by = library]
    sh <- sum(libs$N * LN$tall_share[match(libs$library, LN$gsm)], na.rm = TRUE) / sum(libs$N)
    ct_ <- rbindlist(lapply(libs$library, function(l) fread(file.path(TABWS, "decontx", "GSE185381", paste0(l, "__percell.csv.gz")))[cell %in% kt$cell, .(contamination)]))
    data.table(sample = s, tall_share = sh, contam_median = median(ct_$contamination), contam_P95 = quantile(ct_$contamination, 0.95)) }))
  dose[, `:=`(dose_lo = tall_share * contam_median, dose_hi = tall_share * contam_P95 / 0.62)]
  mu <- rowMeans(L25); sg <- apply(L25, 1, sd)
  sc_all <- function(Lx) { g <- intersect(rownames(Lx), names(mu)); z <- (Lx[g, ] - mu[g]) / sg[g]
    colMeans(z[intersect(P_UP, g), , drop = FALSE]) - colMeans(z[intersect(P_DN, g), , drop = FALSE]) }
  s0 <- sc_all(L25); gap_in <- mean(s0[AMLS]) - mean(s0[HEALTHY])
  shift <- function(col) { Gp <- G_RAW; for (s in T15) Gp[, s] <- Gp[, s] + round(p / sum(p) * dose[sample == s][[col]] * sum(G_RAW[, s])); sp <- sc_all(logz(Gp)); mean(sp[HEALTHY]) - mean(s0[HEALTHY]) }
  D9 <- data.table(quantity = c("K6b_pass", "log2ratio_up_obs", "log2ratio_up_random_median", "log2ratio_up_p_ge", "log2ratio_down_obs", "log2ratio_down_random_median",
                                "log2ratio_down_p_ge", "gap_in_same_scoring", "healthy_shift_lo_dose", "healthy_shift_hi_dose", "shift_frac_lo", "shift_frac_hi"),
                   value = c(as.numeric(isTRUE(k6b)), ru, rd, gap_in, sl <- shift("dose_lo"), sh <- shift("dose_hi"), sl / gap_in, sh / gap_in))
  fwrite(D9, wtab("d9")); fwrite(dose, wtab("d9_dose")); print(D9); print(dose)
  ts("[arms] D10 negative-control compartments, D11 hashing")
  D10 <- rbindlist(lapply(c("B", "NK"), function(cmp) { dc <- dsg(cmp); Mc <- PB$raw[[cmp]][, dc$sample]; if (SMOKE) Mc <- synth(Mc, dc[status == "AML", sample], 4 + (cmp == "NK"))
    r <- fit_de(Mc, dc, MAIN); x <- R_of(r, min_cov = 0); lf <- setNames(r$logFC, r$gene)
    rnd <- vapply(RS[[paste0("d10_", cmp)]], function(g) { g <- intersect(g[!is.na(g)], names(lf)); sum(sign(LFC0[g]) * lf[g]) / sum(abs(LFC0[g])) }, 1)
    data.table(compartment = cmp, n_AML = dc[status == "AML", .N], n_healthy = dc[status == "healthy", .N], coverage = x$cov, R = x$R, R_up = x$R_up, R_down = x$R_down,
               random_median = median(rnd), random_P95 = quantile(rnd, 0.95), p_ge = (1 + sum(rnd >= x$R)) / (length(rnd) + 1)) }))
  fwrite(D10, wtab("d10")); print(D10)
  addR("D11_hashed", fit_de(G_RAW, with_cov(D, as.integer(D$hashed)), c(MAIN, "cov")))
  QCR <- fread(file.path(TAB_DIR, "01_preprocess", "03_qc_report__GSE185381.csv"))
  KF <- fread(file.path(TAB_DIR, "00_ingest", "GSE185381_qc_percell.csv.gz"), select = c("GSM", "Demuxed"))[, .(kept_fraction = mean(Demuxed)), by = GSM]
  D11 <- copy(COV)[, `:=`(qc_mito_ceiling = QCR$mt_hi[match(sample, QCR$Sample)],
                          author_labelled_fraction = vapply(strsplit(libraries, ";"), function(l) mean(KF[GSM %in% l, kept_fraction]), 1),
                          tall_lane_share = dose$tall_share[match(sample, dose$sample)])]
  fwrite(D11, wtab("d11"))
  AR <- rbindlist(A); fwrite(AR, wtab("arms"))
  ts(sprintf("[arms] done: %d arm values written (no verdict here)", nrow(AR)))
  quit(status = 0)
}

## -- 7. MODE task: one slice of the reference fits -------------------------------------------------------------------
if (MODE == "task") {
  TASK <- as.integer(sub("^--task=", "", grep("^--task=", .a, value = TRUE)))
  JT <- JOBS[task == TASK]; ts(sprintf("[task %d] %d jobs (%s)", TASK, nrow(JT), paste(names(table(JT$type)), table(JT$type), sep = ":", collapse = " ")))
  PAIRS <- CJ(a = which(D$status == "AML"), h = which(D$status == "healthy"))
  rowsR <- function(r) { x <- R_of(r); data.table(metric = c("R", "R_up", "R_down", "cov"), value = c(x$R, x$R_up, x$R_down, x$cov)) }
  run_job <- function(type, key) {
    if (type == "g1drop") return(rowsR(fit_de(G_RAW, D[!sample %in% splitv(DR5[id == key, drop])], MAIN)))
    if (type == "d3drop4") return(rowsR(fit_de(G_RAW, D[!sample %in% splitv(DR4[id == key, drop])], MAIN)))
    if (type == "d7drop") return(rowsR(fit_de(G_RAW, D[!sample %in% splitv(DA4[id == key, drop])], MAIN)))
    if (startsWith(type, "wref_")) return(rowsR(fit_de(G_RAW, with_cov(D, WD[[sub("^wref_", "", type)]][, key]), c(MAIN, "cov"))))
    if (type == "d6sub") return(rowsR(fit_de(G_RAW, D[sample %in% splitv(SUB[id == key, subset])], MAIN)))
    if (type == "d2split") return(data.table(metric = "opposite_share", value = d2_fit(splitv(SPL[id == key, group1]))))
    if (type == "g2") { Dp <- copy(D); Dp[, status := factor(ifelse(strsplit(G2L[id == key, labels], "")[[1]] == "1", "AML", "healthy"), c("healthy", "AML"))]
      x <- cf_gap(G_RAW, Dp, L = L25); return(data.table(metric = c("gap", "n_failed"), value = c(x$gap, x$n_failed))) }
    if (type == "d3cf") { x <- cf_gap(G_RAW, D[!sample %in% splitv(DR5[id == key, drop])]); return(data.table(metric = c("gap", "n_failed"), value = c(x$gap, x$n_failed))) }
    if (type == "lpo") { x <- heldout(G_RAW, L25, Df, c(PAIRS$a[key], PAIRS$h[key])); return(data.table(metric = c("score_aml", "score_healthy"), value = x$scores)) }
    if (type == "calib") { k <- key %/% 1000L; i <- key %% 1000L; Gp <- plant_counts(CAL_K[[as.character(k)]])
      r <- if (i == 0L) fit_de(Gp, D, MAIN) else fit_de(Gp, D[!sample %in% splitv(DR5[id == i, drop])], MAIN)
      x <- r[gene %in% P370$gene]; return(data.table(metric = paste0("lfc:", x$gene), value = x$logFC)) }
    stop("unknown job type ", type) }
  RES <- rbindlist(lapply(seq_len(nrow(JT)), function(k) { j <- JT[k]
    out <- tryCatch(run_job(j$type, j$kid), error = function(e) data.table(metric = "error", value = NA_real_))
    ts(sprintf("    job %5d %-10s key %6d: %s", j$job, j$type, j$kid, if (any(out$metric == "error")) "FAILED" else "ok"))
    cbind(j[, .(job, type, kid)], out) }))
  fn <- file.path(TDIR, sprintf("task_%03d%s.csv", TASK, SFX)); if (file.exists(fn) && !SMOKE) stop("refusing to overwrite ", fn)
  fwrite(RES, fn); ts("[task ", TASK, "] written ", fn); quit(status = 0)
}

## -- 8. MODE combine: references, percentiles, verdicts, GATE_A1B ------------------------------------------------------
if (MODE == "combine") {
  fs <- list.files(TDIR, pattern = paste0("^task_[0-9]+", SFX, "\\.csv$"), full.names = TRUE)
  if (!SMOKE) stopifnot(length(fs) == N_TASKS)
  T_ <- rbindlist(lapply(fs, fread)); AR <- fread(file.path(TABWS, paste0("ws_a1_batch_arms", SFX, ".csv")))
  av <- function(a, m = "R") AR[arm == a & metric == m, value]                     # (function arguments renamed so they differ from the columns)
  jv <- function(tp, m = "R") T_[type == tp & metric == m][order(kid)]
  nfail <- T_[metric == "error", .N, by = type]
  first_ok <- function(x, n = N_REF) x[!is.na(value)][seq_len(min(.N, n))]               # reserves replace failed draws
  RES <- list(); res <- function(readout, value, reference, n_ref, pct, verdict, note = "")
    RES[[length(RES) + 1]] <<- data.table(readout, value, reference, n_ref, percentile = pct, verdict, note)
  ## G1: drop Control1-5 against the 40 matched drops (arm included); secondary z over all 252
  g1 <- jv("g1drop"); g1 <- merge(g1, DR5[, .(kid = id, matched, real, dv)], by = "kid")
  rG1 <- av("G1"); chk("combine: the job fit of the real G1 drop equals the arm fit", rG1, g1[real == TRUE, value], tol = 1e-8)
  ref40 <- c(rG1, g1[matched == TRUE & real == FALSE & !is.na(value), value]); pG1 <- pct_enum(rG1, ref40)
  z <- g1[, (value - 1) / sqrt(dv)]; zG1 <- (rG1 - 1) / sqrt(DR5[real == TRUE, dv]); pz <- pct_enum(zG1, z)
  vG1 <- if (rG1 < 0.6 && pG1 <= 0.05) "explains" else if (pG1 <= 0.05 && rG1 < 0.85) "weakens" else "no effect"
  res("G1 drop Control1-5: R", rG1, "40 design-matched drops", length(ref40), pG1, vG1, sprintf("fhat %.2f (%.2f-%.2f); z-percentile over 252: %.3f", fhat(rG1), fhat(rG1, PHI_RANGE[1]), fhat(rG1, PHI_RANGE[2]), pz))
  ## G2: block-respecting relabellings of the cross-fitted gap
  g2 <- first_ok(jv("g2", "gap")); pG2 <- (1 + sum(g2$value >= GAP_OOF)) / (nrow(g2) + 1)
  vG2 <- if (pG2 > 0.10) "explains" else if (pG2 > 0.05) "weakens" else "no effect"
  res("G2 block-respecting p of the cross-fitted gap", pG2, "block-respecting relabellings", nrow(g2), pG2, vG2, sprintf("null median %.3f, P95 %.3f; free p (A1-G02) 0.008", median(g2$value), quantile(g2$value, 0.95)))
  ## W1-W4: covariates, residual-permutation references
  for (w in names(WCOV)) { r <- av(w); ref <- first_ok(jv(paste0("wref_", w)))$value; p <- pct_drawn(r, ref)
    res(sprintf("%s %s: R", w, WCOV[[w]]), r, "1,000 residual-permutation draws", length(ref), p, if (p <= 0.05 && r < 0.85) "weakens" else "no effect",
        sprintf("fhat %.2f; covariate cannot explain (t cap below 2.78)", fhat(r))) }
  ## descriptive arms
  cf3 <- jv("d3cf", "gap"); cfr <- av("D3", "CF_without_Control1-5") / median(cf3$value, na.rm = TRUE)
  res("D3 CF without Control1-5 / median matched CF", cfr, "39 matched drops", nrow(cf3), NA, "descriptive")
  r <- av("D3_P02"); ref <- c(r, jv("d3drop4")$value); res("D3 drop P02: R", r, "13 matched drops of 4 healthy", length(ref), pct_enum(r, ref), "descriptive")
  r <- av("D7_P01"); ref <- c(r, jv("d7drop")$value); res("D7 drop P01: R", r, "matched drops of 4 AML", length(ref), pct_enum(r, ref), "descriptive")
  r <- av("D6_year"); ref <- first_ok(jv("wref_D6"))$value; res("D6 calendar year: R", r, "1,000 residual-permutation draws", length(ref), pct_drawn(r, ref), "descriptive")
  r <- av("D6_early"); ref <- jv("d6sub")$value; res("D6 early stratum: R", r, "500 subsets (P01, P02 whole or absent)", length(ref), pct_drawn(r, ref), "descriptive")
  real_split <- SPL[vapply(strsplit(group1, ";"), function(x) setequal(x, T15), TRUE), id]
  r <- av("D2", "opposite_share"); ref <- c(r, jv("d2split", "opposite_share")[kid != real_split, value])
  res("D2 within-healthy: share of P370 with Control1-5 sign opposite to the imprint", r, "252 splits of the 10 healthy", length(ref), pct_enum(r, ref), "descriptive",
      sprintf("p_ge %.3f", mean(ref >= r)))
  lpo <- merge(jv("lpo", "score_aml")[, .(kid, a = value)], jv("lpo", "score_healthy")[, .(kid, h = value)], by = "kid")
  res("D5 leave-pair-out AUC with Control4 coded female", mean((lpo$a > lpo$h) + 0.5 * (lpo$a == lpo$h)), "150 pairs (A1-G02 D1 = 0.927)", nrow(lpo), NA, "descriptive")
  ## calibration report: smallest healthy-ward plant on Control1-5 at which G1's percentile <= 0.05
  cal <- T_[type == "calib"][, `:=`(k = kid %/% 1000L, i = kid %% 1000L)]
  CAL <- rbindlist(lapply(1:3, function(kk) { m <- cal[k == kk & i == 0, .(gene = sub("^lfc:", "", metric), lfc0 = value)]
    Rk <- cal[k == kk & i > 0, { x <- merge(data.table(gene = sub("^lfc:", "", metric), lfc = value), merge(P370, m, by = "gene"), by = "gene"); .(R = sum(x$s * x$lfc) / sum(abs(x$lfc0))) }, by = i]
    Rk <- merge(Rk, DR5[, .(i = id, real)], by = "i"); data.table(plant = CAL_K[[as.character(kk)]], R_arm = Rk[real == TRUE, R], pct = pct_enum(Rk[real == TRUE, R], Rk$R)) }))
  RT <- rbindlist(RES)
  gating <- RT[grepl("^(G1|G2|W[1-4]) ", readout)]
  GATE <- if (any(gating$verdict == "explains")) "FAIL" else if (any(gating$verdict == "weakens")) "PARTIAL" else "PASS"
  C <- rbindlist(CHK)
  fwrite(RT, wtab("results")); fwrite(CAL, wtab("calibration")); fwrite(nfail, wtab("failed_jobs")); fwrite(C, wtab("checks_combine"))
  print(RT); print(CAL); print(nfail); ts(sprintf("[combine] GATE_A1B = %s", GATE))
  ## GATES.md block (real run only; smoke writes it next to the draws)
  blk <- c("<!-- A1-BATCH start -->", sprintf("## T imprint, processing factors inside GSE185381 (A1-BATCH, prereg-a1-batch-v1) -- %s", format(Sys.Date())), "",
           "Scripts: scripts/ws/37_a1_batch.R, scripts/ws/38_a1_batch_cells.R. 15 adult AML vs 10 healthy, S5 design, frozen P370 direction.", "",
           "| readout | value | reference | n | percentile | verdict | note |", "|---|---|---|---|---|---|---|",
           RT[, sprintf("| %s | %.3f | %s | %d | %s | %s | %s |", readout, value, reference, n_ref, ifelse(is.na(percentile), "NA", sprintf("%.3f", percentile)), verdict, gsub("\\|", "\\\\|", note))],
           "", sprintf("Calibration (healthy-ward plant on Control1-5, x main logFC): %s", CAL[, paste(sprintf("%.2f: R %.3f, percentile %.3f", plant, R_arm, pct), collapse = "; ")]),
           "", sprintf("**Verdict: GATE_A1B %s.**", GATE), "<!-- A1-BATCH end -->")
  gf <- if (SMOKE) file.path(BDIR, "gates_block_smoke.md") else here::here("scripts", "ws", "GATES.md")
  if (!SMOKE && any(grepl("A1-BATCH start", readLines(gf)))) stop("GATES.md already holds an A1-BATCH block")
  if (SMOKE) writeLines(blk, gf) else write(c("", blk), gf, append = TRUE)
  ts("[combine] done"); quit(status = 0)
}
