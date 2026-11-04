#!/usr/bin/env Rscript
# INPUT  : the inputs of PREREGISTRATION_A1_G04.md sections 2.1-2.3, each with its declared md5 (checked by K0):
#          /LARGE1/.../02_seurat_objects/04_annotated/Chen2023/<donor>.rds (counts + meta.data Pool, Library, bmm_fine, anno_source,
#          high_error, hierarchy_bin), meta_Chen2023_v2.xlsx, the Zenodo features.tsv.gz, GSE120221_gene_remap.csv (--build only);
#          ws_pb_GSE185381.rds (25 adult columns; the 18 children's $raw$T columns and $cells rows), results/tables/ws/ ws_s5_*,
#          ws_s8_*, ws_a1_g02_*, ws_s9_*; the S9 (Oetjen) inputs of K4
# OUTPUT : results/tables/ws/ws_a1_g04_{checks,primary,loo,sens,children,readouts,summary}.csv; /LARGE1/.../ws_a1_g04/chen_pb.rds (--build);
#          GATES.md block A1-G04 (--combine). Smoke runs write only to /LARGE1/.../ws_a1_g04/smoke/ (*_smoke), never under results/.
# what it does: PREREGISTRATION_A1_G04.md (tag prereg-a1-g04-v1), nothing else.
#   --build     Chen2023 pseudobulks on the remapped axis (2.2), decontX per pool-B library (4.7); prints only cell numbers
#   --known     known answers of section 6 (K0-K9); writes and prints no imprint or delay score
#   --primary   3.1-3.2     --sens  4.1-4.7     --readouts  5.1-5.4     --children  3.4, 4.8-4.11, 5.5, 5.6     --combine  section 7
#   --smoke     with any mode: synthetic Chen2023 and child counts (real IDs only), few draws; outputs end in _smoke
#   usage: sbatch scripts/ws/submit_ws_a1_g04.sbatch --known   (the wrapper runs a dated read-only copy of this file)
# This script deletes nothing, modifies no input, and never runs or sources another script: the functions it needs are copied below
#   with their source lines. It opens only the files of the pre-registration's allowlist (section 8).

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(Matrix); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODES <- c("build", "known", "primary", "sens", "readouts", "children", "combine")
MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of ", paste0("--", MODES, collapse = " "), " [--smoke]")
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
TABWS <- file.path(TAB_DIR, "ws"); METADIR <- file.path(RAW_DIR, "metadata")
GDIR <- file.path(LARGE1_DIR, "ws_a1_g04"); SDIR <- file.path(GDIR, "smoke")
WDIR <- if (SMOKE) SDIR else GDIR                     # matrices and draws
ODIR <- if (SMOKE) SDIR else TABWS                    # tables: smoke tables never go under results/, so the 23:30 cron cannot push them
dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])   # path of the script copy being run
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
## write guards (39's pattern): the only places this script writes; a real (non-smoke) file is never overwritten; every table records the md5
wtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g04_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); d[, script_md5 := MD5]; fwrite(d, wtab(x)); invisible(d) }
rtab <- function(x) fread(file.path(ODIR, paste0("ws_a1_g04_", x, SFX, ".csv")))
wg   <- function(x) { f <- file.path(WDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
rg   <- function(x) file.path(WDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x))
ts(sprintf("[%s%s] %d worker processes", MODE, if (SMOKE) " smoke" else "", NC))
NOTED <- c(k3 = "--k3-noted" %in% .a, k3b = "--k3b-noted" %in% .a)   # section 6: a K3 / K3b mismatch needs a dated note first, then this flag
if (!SMOKE) {   # every target of this mode is checked before any work, so a refused write cannot leave a partial result (code review 1)
  tg <- switch(MODE, build = file.path(GDIR, "chen_pb.rds"), known = file.path(TABWS, "ws_a1_g04_checks.csv"),
               primary = file.path(TABWS, paste0("ws_a1_g04_", c("primary", "loo"), ".csv")), sens = file.path(TABWS, "ws_a1_g04_sens.csv"),
               readouts = c(file.path(TABWS, "ws_a1_g04_readouts.csv"), file.path(GDIR, "readouts_C.rds")),
               children = file.path(TABWS, "ws_a1_g04_children.csv"), combine = file.path(TABWS, "ws_a1_g04_summary.csv"))
  if (any(file.exists(tg))) stop("refusing to overwrite ", paste(tg[file.exists(tg)], collapse = ", "), " (section 6: a rerun needs a dated note and the old file renamed)")
  if (MODE == "combine" && any(grepl("A1-G04 start", readLines(here::here("scripts", "ws", "GATES.md"))))) stop("GATES.md already holds an A1-G04 block")
}

## -- 0. constants (section numbers refer to the pre-registration) ----------------------------------------------------------
OFF <- c(boot = 17L, k6b = 18L, smoke = 19L, thin = 20L, k6c = 25L, depth = 26L, k7c = 27L, k6a = 28L)   # seeds = SEED + offset (section 8)
M_PASS <- 6L                  # pass level m/126 (decision 9.2)
MIN_T <- WS_S4_MIN_CELLS      # 20 kept cells per donor or child (2.2, 2.3)
FB_CELLS <- 50L; MIN_GENES <- 5L; COV_MIN <- 0.80; INFL <- 2 / 3; R_YARD <- 1.430
if (SMOKE) { set.seed(SEED + OFF[["smoke"]]); SMS <- sample.int(.Machine$integer.max, 20L) }   # every smoke draw derives from SEED + 19
N_BOOT <- if (SMOKE) 50L else 2000L; N_K6B <- if (SMOKE) 20L else 1000L; N_THIN <- if (SMOKE) 2L else 20L
DRES_S9 <- 1.96770393747843; GAPIN_K4 <- 1.97138    # K4 targets (ws_s9_direction.csv; review 1)
POOLA <- "poolA_HSPC_myeloid"; POOLB <- "poolB_stroma_lymphoid"
MYEL_BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte")   # the pool-B myeloid/progenitor share (section 0, 3.1)
# the declared donor table (2.2); K2 compares it with the xlsx and with --build's counts
DON <- data.table(donor = c("AML103", "AML24", "AML162", "AML304", "AML320", "NBM1", "NBM2", "NBM3", "NBM4", "AML469"),
                  arm = c(rep("AML", 5), rep("donor", 4), "AML"), age = c(59, 59, 51, 56, 39, 47, 52, 50, 53, 50),
                  sex = c("F", "M", "M", "M", "F", "M", "F", "M", "F", "M"),
                  keptT_B = c(211L, 1047L, 890L, 1445L, 808L, 2306L, 1094L, 2795L, 1709L, 2L), eligible = c(rep(TRUE, 9), FALSE))
ELIG <- DON[eligible == TRUE, donor]; ARM <- setNames(DON$arm, DON$donor)[ELIG]   # AML first, then the 4 donors (the order of 3.1)
stopifnot(identical(unname(ARM), c(rep("AML", 5), rep("donor", 4))))
CHILD <- c("AML001", "AML003", "AML005", "AML006", "AML009", "AML012", "AML022", "AML025", "AML026", "AML028", "AML038", "AML043",
           "AML048", "AML051", "AML052", "AML055", "AML056", "AML073")                                          # 2.3
UNDER18 <- setdiff(CHILD, c("AML022", "AML038", "AML043", "AML073"))                                           # 4.8
FLAGFREE <- c("AML001", "AML005", "AML006", "AML009", "AML025", "AML048", "AML055")                               # 4.10
MIXED <- c("AML056", "AML038", "AML026", "AML051")                                                               # K3b
stopifnot(length(CHILD) == 18L, length(UNDER18) == 14L, all(FLAGFREE %in% CHILD))

## known-answer bookkeeping (39's pattern; section 6's three states, plus "reported")
CHK <- list()
chk <- function(id, scope, name, expected, got, tol = NULL) {    # deterministic: PASS or FAIL_deterministic
  ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "deterministic", check = name, expected = paste(format(expected), collapse = ";"),
                                         got = paste(format(got), collapse = ";"), state = if (ok) "PASS" else "FAIL_deterministic")
  message(sprintf("  %-4s %-78s expect %-14s got %-14s %s", id, name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
bar <- function(id, name, value, ok) {                           # statistical bar or coverage floor: PASS or FAIL_statistical (final)
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = "P", kind = "statistical", check = name, expected = "", got = paste(signif(value, 4), collapse = ";"),
                                         state = if (isTRUE(ok)) "PASS" else "FAIL_statistical")
  message(sprintf("  %-4s %-78s got %-14s %s", id, name, paste(signif(value, 4), collapse = ";"), if (isTRUE(ok)) "PASS" else "** FAIL (statistical) **")) }
rep_ <- function(id, scope, name, value, state = "REPORTED") {    # reported rows never block (except K3's swap state, section 6)
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "reported", check = name, expected = "", got = value, state = state)
  message(sprintf("  %-4s %-78s %s [%s]", id, name, value, state)) }
checks <- function() { f <- file.path(ODIR, paste0("ws_a1_g04_checks", SFX, ".csv")); if (!file.exists(f)) stop("run --known first: ", f); fread(f) }
gate <- function(sc) {   # a mode runs only when every deterministic check of its scope passed, and (scope P) K3's swap rule did not stop it
  C <- checks(); stopifnot(nrow(C[scope == "P" & kind == "deterministic"]) > 0)   # sc, not scope: inside C[...] "scope" is the column
  if (!SMOKE && !identical(unique(C$script_md5), MD5)) stop("the checks table was written by another version of this script; rerun --known (section 6)")
  ok <- all(C[scope %in% c(sc, "PB") & kind == "deterministic", state] == "PASS")
  if (sc == "P" && any(C$state == "SWAP_STOP")) { ts("K3: a cross-arm swap pattern stops the run before the primary (section 6)"); ok <- FALSE }
  ok }
syn <- function(M) { attr(M, "synthetic") <- TRUE; M }           # marks a smoke matrix; need_syn() refuses anything else in a smoke run
need_syn <- function(M, what) if (SMOKE && !isTRUE(attr(M, "synthetic"))) stop("smoke run: a real ", what, " column would enter a computation")

## -- 1. gene sets and blacklist (33:53-62, verbatim apart from names) ------------------------------------------------------------
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]
UP <- P370[s > 0, gene]; DN <- P370[s < 0, gene]; UP402 <- F402[s > 0, gene]; DN402 <- F402[s < 0, gene]
stopifnot(nrow(F402) == 402L, nrow(P370) == 370L, length(UP) == 156L, length(DN) == 214L)
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
black_for <- function(G) unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene], "XIST",
                                  GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
MAP <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE)
TL <- MAP[group4 == "T", bmm_fine]; NAIVE <- c("CD4 Naive", "CD8 Naive"); CD8L <- grep("^CD8 ", TL, value = TRUE)
stopifnot(length(TL) == 9L, all(NAIVE %in% TL))

## -- 2. score engines (copied; section 8) ------------------------------------------------------------------------------------
logzZ <- function(M) {   # 33:64-70 without the matched-bin table: blacklist; CPM > 1 in >= half; TMM; log CPM prior 1; z across the set
  M <- M[!rownames(M) %in% black_for(rownames(M)), , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L <- cpm(y, log = TRUE, prior.count = 1); list(Z = t(scale(t(L))), L = L) }
logzL <- function(M) {   # 39:91-94 = 37:65-68 (returns L); black_for() of the matrix's own axis (= BLACK on the GSE185381 axis, review 1)
  M <- M[!rownames(M) %in% black_for(rownames(M)), , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1) }
score <- function(Z, up, dn) { up <- intersect(up, rownames(Z)); dn <- intersect(dn, rownames(Z))   # 33:71-72
  (if (length(up)) colMeans(Z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(Z[dn, , drop = FALSE]) else 0) }
score_col <- function(x, mu, sg, fu, fd) { z <- (x - mu) / sg; z <- z[is.finite(z)]   # 37:101-102 (one column, one fold's lists)
  (if (length(intersect(fu, names(z)))) mean(z[intersect(fu, names(z))]) else 0) - (if (length(intersect(fd, names(z)))) mean(z[intersect(fd, names(z))]) else 0) }
jspace <- function(G_ref, st_ref, E_raw) {   # 33:171-182 with the reference and its status as arguments (T, B and NK references, 5.3)
  common <- setdiff(intersect(rownames(G_ref), rownames(E_raw)), union(black_for(rownames(G_ref)), black_for(rownames(E_raw))))
  yG <- DGEList(G_ref[!rownames(G_ref) %in% black_for(rownames(G_ref)), , drop = FALSE]); yE <- DGEList(E_raw[!rownames(E_raw) %in% black_for(rownames(E_raw)), , drop = FALSE])
  kG <- rownames(yG)[rowMeans(cpm(yG) > 1) >= 0.5]; kE <- rownames(yE)[rowMeans(cpm(yE) > 1) >= 0.5]
  U_J <- intersect(intersect(kG, kE), common)
  yG <- calcNormFactors(yG[U_J, , keep.lib.sizes = FALSE]); yE <- calcNormFactors(yE[U_J, , keep.lib.sizes = FALSE])
  LG <- cpm(yG, log = TRUE, prior.count = 1); LE <- cpm(yE, log = TRUE, prior.count = 1)
  m <- rowMeans(LG); s <- apply(LG, 1, sd); ok <- s > 1e-8; U_J <- U_J[ok]; LG <- LG[ok, , drop = FALSE]; LE <- LE[ok, , drop = FALSE]; m <- m[ok]; s <- s[ok]
  list(U_J = U_J, LG = LG, LE = LE, m = m, s = s, ZG = (LG - m) / s, ZE = (LE - m) / s,
       d_res = { sg <- score((LG - m) / s, UP, DN); mean(sg[st_ref[names(sg)] == "AML"]) - mean(sg[st_ref[names(sg)] == "healthy"]) }) }
make_L27 <- function(M25, U) {   # 37:108-118 verbatim: the 25 columns exactly as logz; new columns TMM'd against the 25's reference column
  M <- M25[!rownames(M25) %in% BLACK, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L25 <- cpm(y, log = TRUE, prior.count = 1); Mk <- M[keep, , drop = FALSE]
  f75 <- apply(t(t(Mk) / colSums(Mk)), 2, quantile, 0.75); ref <- which.min(abs(f75 - mean(f75)))   # calcNormFactors' default reference
  eff25 <- y$samples$lib.size * y$samples$norm.factors
  man <- function(X, eff) { pc <- eff / mean(eff25); log2(t((t(X) + pc) / (eff + 2 * pc)) * 1e6) }  # edgeR's log-CPM with the 25's mean library size
  Uk <- as.matrix(U[rownames(Mk), , drop = FALSE])
  fu <- vapply(seq_len(ncol(Uk)), function(j) { z <- calcNormFactors(DGEList(cbind(Mk[, ref], Uk[, j])), refColumn = 1)$samples$norm.factors; z[2] / z[1] }, 1)
  LU <- man(Uk, colSums(Uk) * fu * y$samples$norm.factors[ref]); colnames(LU) <- colnames(U)
  list(L = cbind(L25, LU), check = man(Mk, eff25) - L25) }
age_group_of <- function(a) ifelse(a <= 53, "young", ifelse(a >= 60, "old", "middle"))   # 33:43-45 constants and 33:112 (K4 only)
load_cohort <- function(ds) {   # 33:113-126 verbatim (K4 only)
  PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", ds, paste0("ws_pb_", ds, ".rds")))
  META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", ds, ".csv")))
  don <- fread(file.path(METADIR, ds, paste0("donor_ages_", ds, "_v1.0.csv")), colClasses = "character")
  meta <- fread(file.path(METADIR, ds, paste0("meta_", ds, "_v1.0.csv")), colClasses = "character")
  lib <- unique(META[, .(sample, n_kept)]); lib[, gsm := sub("_.*$", "", sample)]
  lib[, `:=`(donor = don$donor[match(gsm, don$gsm)], age_lib = as.numeric(don$age[match(gsm, don$gsm)]), sex_lib = don$sex[match(gsm, don$gsm)])]
  stopifnot(!anyNA(lib$donor), all(colnames(PB$raw$T) %in% lib$sample))
  nT <- META[group4 == "T", .(sample, n_T = n)]; lib <- merge(lib, nT, by = "sample", all.x = TRUE); lib[is.na(n_T), n_T := 0L]
  donors <- lib[order(gsm), .(n_lib = .N, libraries = paste(sample, collapse = ";"), age = age_lib[1], sex = sex_lib[1], n_T = sum(n_T), n_kept = sum(n_kept)), by = donor]
  donors[, `:=`(eligible = n_T >= MIN_T, grp = age_group_of(age), excluded_K3 = FALSE)]
  T_don <- sapply(donors$donor, function(d) rowSums(PB$raw$T[, lib[donor == d, sample], drop = FALSE])); colnames(T_don) <- donors$donor
  list(ds = ds, donors = donors, T_don = T_don) }
sexcall <- function(xist, rps) fifelse(is.na(xist) | is.na(rps), NA_character_,       # 26_s4_pseudobulk.R:91
                                      fifelse(xist > 3 * pmax(rps, 1), "F", fifelse(rps > 3 * pmax(xist, 1), "M", "unclear")))

## -- 3. exact enumeration, bounds, donor influence (3.1-3.2) ----------------------------------------------------------------------
enum <- function(na, nh) {   # every labelling of na + nh columns (AML first) into na "AML" and nh others; W gives each labelling's mean difference
  L <- combn(na + nh, na); obs <- which(apply(L, 2, function(s) all(s == seq_len(na))))
  W <- apply(L, 2, function(s) { w <- rep(-1 / nh, na + nh); w[s] <- 1 / na; w })   # column = weights; crossprod(W, x) = all D values
  list(L = L, W = W, obs = obs, k = colSums(L <= na), na = na, nh = nh) }
E54 <- enum(5L, 4L); E53 <- enum(5L, 3L); E44 <- enum(4L, 4L)
Dall <- function(x, E) drop(crossprod(E$W, x))
pexact <- function(x, E, tail = c("upper", "lower")) { tail <- match.arg(tail); D <- Dall(x, E); d0 <- D[E$obs]
  if (tail == "upper") sum(D >= d0 - 1e-12) / ncol(E$L) else sum(D <= d0 + 1e-12) / ncol(E$L) }
Ucount <- function(x, na) sum(outer(x[seq_len(na)], x[-seq_len(na)], ">"))      # rank-sum count: pairs with the AML donor higher
bounds <- function(x, m = M_PASS, E = E54) {   # 3.1 closed form: t_S = (D_obs - D_S) / (1 - c_S); delta_L, delta_U = m-th smallest, m-th largest
  D <- Dall(x, E); cS <- E$k / E$na - (E$na - E$k) / E$nh; o <- setdiff(seq_along(D), E$obs)
  t <- (D[E$obs] - D[o]) / (1 - cS[o]); c(lo = sort(t)[m], hi = sort(t, decreasing = TRUE)[m]) }
bounds_bisect <- function(x, m = M_PASS, E = E54) {   # K5: the same bounds by inverting the exact test numerically
  cS <- E$k / E$na - (E$na - E$k) / E$nh; D <- Dall(x, E)
  pu <- function(d) sum(D - cS * d >= D[E$obs] - d - 1e-12) / ncol(E$L)   # upper-tail p after subtracting d from every AML score
  pl <- function(d) sum(D - cS * d <= D[E$obs] - d + 1e-12) / ncol(E$L)
  find <- function(f, lo, hi) { for (i in 1:200) { mid <- (lo + hi) / 2; if (f(mid)) hi <- mid else lo <- mid }; (lo + hi) / 2 }
  w <- 100 * (max(abs(D)) + 1)
  c(lo = find(function(d) pu(d) > m / ncol(E$L), -w, w), hi = -find(function(d) pl(-d) > m / ncol(E$L), -w, w)) }
influence <- function(x, na = 5L, sign = 1) {   # 3.2: D without each donor (fixed scores); sign = -1 mirrors the rule for "reversed"
  a <- x[seq_len(na)]; h <- x[-seq_len(na)]
  d <- c(vapply(seq_along(a), function(i) mean(a[-i]) - mean(h), 1), vapply(seq_along(h), function(j) mean(a) - mean(h[-j]), 1))
  names(d) <- names(x); sD <- sign * (mean(a) - mean(h)); sd_ <- sign * d    # in the direction of the result
  list(D_minus = d, flag = sD > 0 && min(sd_) < INFL * sD, named = names(d)[which.min(sd_)]) }

## -- 4. (b) machinery: out-of-fold placement and A_b (2.4, 3.4) --------------------------------------------------------------------
Ab_of <- function(S, C) { x <- (S > C) + 0.5 * (S == C); mean(colMeans(x)) }      # S, C: children x folds; ties count 0.5
grid_boot <- function(S, C, extra = NULL, B = N_BOOT) {   # 3.4: rows (children) and columns (folds) resampled as units on the grid; SEED + 17
  set.seed(SEED + OFF[["boot"]]); n <- nrow(S); k <- ncol(S)                # extra(i, j): a second statistic on the same draw (phi)
  vapply(seq_len(B), function(b) { i <- sample.int(n, n, TRUE); j <- sample.int(k, k, TRUE)
    c(Ab_of(S[i, j, drop = FALSE], C[i, j, drop = FALSE]), if (is.null(extra)) NA_real_ else extra(i, j)) }, numeric(2)) }
q7 <- function(v) quantile(v, c(0.025, 0.975), type = 7, names = FALSE)

## -- 4b. the counterfactual of A1-G01 3.1 (39:151-195 without 39's UMI and Poisson branches; logz -> logzL) ------------------------
inst <- function(ids) paste0(ids, "@", seq_along(ids))            # unique column names for repeated samples (DGEList needs them)
halves <- function(X, j, oth, U, D) {   # u = up + dn: X[g, j] minus the mean over the other columns, averaged over U and (negated) over D
  m <- rowMeans(X[c(U, D), oth, drop = FALSE]); x <- X[c(U, D), j] - m
  c(up = mean(x[seq_along(U)]), dn = -mean(x[length(U) + seq_along(D)])) }
mlist <- function(f, unis) { U <- setdiff(f$U, f$D); D <- setdiff(f$D, f$U)   # genes in every universe; NULL below 5 per half
  for (u in unis) { U <- intersect(U, u); D <- intersect(D, u) }
  if (length(U) < MIN_GENES || length(D) < MIN_GENES) NULL else list(U = U, D = D) }
ref_profiles <- function(C, N, w_arm, w_all, labs, mode) {   # 39:158-164; w_* = how often each sample enters the pool
  pool <- function(k, w) { n <- sum(N[k, ] * w); if (n > 0) drop(C[[k]] %*% w) / n else numeric(nrow(C[[k]])) }   # %*% sums the weighted columns
  na <- vapply(labs, function(k) sum(N[k, ] * w_arm), 1)
  fb <- if (mode == "arm_only") na == 0 else na < FB_CELLS
  R <- vapply(labs, function(k) if (fb[[k]]) pool(k, w_all) else pool(k, w_arm), numeric(nrow(C[[1]])))
  dimnames(R) <- list(rownames(C[[1]]), labs)
  list(R = R, fb = fb, none = vapply(labs, function(k) fb[[k]] && sum(N[k, ] * w_all) == 0, TRUE)) }
est_C <- function(ids, st, LST, C, N, W, labs = TL, mode = "fallback", refs = FALSE) {   # 39:165-195
  smp <- colnames(W); cn <- inst(ids)
  mult <- function(sel) tabulate(match(ids[sel], smp), length(smp))   # multiplicity of each sample among the selected instances
  XT <- logzL(`colnames<-`(W[, ids, drop = FALSE], cn))                # `colnames<-`(M, v) returns M with new column names
  NN <- N[labs, ids, drop = FALSE]
  sc <- data.table(inst = cn, sample = ids, status = st, cH = NA_real_, cA = NA_real_, t = NA_real_, nU = 0L, nD = 0L)
  FB <- list(); RF <- list()
  score1 <- function(YH, YA, f, j, oth) { ml <- mlist(f, list(rownames(YH), rownames(YA), rownames(XT)))
    if (is.null(ml) || !length(oth)) return(c(NA, NA, NA, 0, 0))
    c(sum(halves(YH, j, oth, ml$U, ml$D)), sum(halves(YA, j, oth, ml$U, ml$D)), sum(halves(XT, j, oth, ml$U, ml$D)), length(ml$U), length(ml$D)) }
  for (i in unique(ids)) {
    J <- which(ids == i); oth <- which(ids != i); Y <- list()
    for (o in c("H", "A")) {
      rp <- ref_profiles(C, N, mult(st == (if (o == "H") "healthy" else "AML") & ids != i), mult(ids != i), labs, mode)
      FB[[length(FB) + 1L]] <- data.table(fold = i, order = o, label = labs, fallback = rp$fb, empty = rp$none)
      if (refs) RF[[paste(i, o)]] <- rp$R
      Y[[o]] <- `colnames<-`(rp$R %*% NN, cn) }
    v <- score1(logzL(Y$H), logzL(Y$A), LST[[i]], J[1], oth)
    sc[J, `:=`(cH = v[1], cA = v[2], t = v[3], nU = as.integer(v[4]), nD = as.integer(v[5]))] }
  gp <- function(x) { a <- sc$status == "AML"; mean(x[a], na.rm = TRUE) - mean(x[!a], na.rm = TRUE) }
  gT <- gp(sc$t); gH <- gp(sc$cH); gA <- gp(sc$cA); def <- is.finite(gT) && gT > 0
  CH_ <- if (def) gH / gT else NA_real_; CA_ <- if (def) gA / gT else NA_real_
  list(C_H = CH_, C_A = CA_, C = (CH_ + CA_) / 2, gapT = gT, gapC_H = gH, gapC_A = gA, scores = sc, fallback = rbindlist(FB), refs = if (refs) RF else NULL) }
mc_ok <- function(res, what) {   # 39:229-233: a forked worker that dies returns NULL, an uncaught error a try-error: stop on both
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
  if (any(bad)) stop(what, ": ", sum(bad), " of ", length(res), " results missing; rerun")
  res }

## -- 4c. (a) summaries ----------------------------------------------------------------------------------------------------------------
SV <- fread("/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv")
SVU <- setdiff(SV[direction == "up", gene], F402$gene); SVD <- setdiff(SV[direction == "down", gene], F402$gene)   # 5.2: F402 removed
d54 <- function(x) mean(x[1:5]) - mean(x[6:9])
dres <- function(J, up = UP, dn = DN) { sg <- score(J$ZG, up, dn); mean(sg[AMLS]) - mean(sg[HEAL]) }   # D_res_J on the same U_J and genes
p_rank <- function(x, E = E54) { U <- Ucount(x, E$na); Us <- apply(E$L, 2, function(s) sum(outer(x[s], x[-s], ">"))); sum(Us >= U) / ncol(E$L) }
a_row <- function(what, x, J = NULL, up = UP, dn = DN, d_res = NULL) {   # one (a) analysis: D, exact p (both tails), rank-sum, rho and its bounds
  x <- x[ELIG]; D <- d54(x)
  r <- data.table(analysis = what, D = D, p_upper = pexact(x, E54, "upper"), p_lower = pexact(x, E54, "lower"), U = Ucount(x, 5L), p_rank = p_rank(x),
                  D_res_J = NA_real_, rho = NA_real_, rho_L = NA_real_, rho_U = NA_real_, n_up = NA_integer_, n_down = NA_integer_)
  if (!is.null(J)) { nu <- sum(up %in% J$U_J); nd <- sum(dn %in% J$U_J); r[, `:=`(n_up = nu, n_down = nd)]
    if ((length(up) && nu < MIN_GENES) || (length(dn) && nd < MIN_GENES)) return(r[, `:=`(D = NA_real_, p_upper = NA_real_, p_lower = NA_real_, p_rank = NA_real_)]) }
  if (is.null(d_res) && !is.null(J)) d_res <- dres(J, up, dn)
  if (!is.null(d_res)) { b <- bounds(x); k <- d_res / R_YARD
    r[, `:=`(D_res_J = d_res, rho = D / k, rho_L = b[["lo"]] / k, rho_U = b[["hi"]] / k)] }
  r }

## -- 5. GSE185381 (25 adults; the children only where the mode needs them) ---------------------------------------------------------
if (MODE != "combine") {
  PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
  G185 <- PB$genes; BLACK <- black_for(G185)
  S5ALL <- fread(file.path(TABWS, "ws_s5_scores.csv")); S5 <- S5ALL[compartment == "T"]          # file order = S5's design order
  SAMP <- S5$sample; STAT <- setNames(S5$status, SAMP); AMLS <- SAMP[STAT == "AML"]; HEAL <- SAMP[STAT == "healthy"]
  stopifnot(length(SAMP) == 25L, length(AMLS) == 15L, !any(c(CHILD, "PAUMTZ", "PAXMIJ") %in% c(SAMP, S5ALL$sample)))
  G_RAW <- PB$raw$T[, SAMP]; rownames(G_RAW) <- G185
  REF <- lapply(c(B = "B", NK = "NK"), function(cmp) { s <- S5ALL[compartment == cmp, sample]; M <- PB$raw[[cmp]][, s]; rownames(M) <- G185
    list(M = M, st = setNames(S5ALL[compartment == cmp, status], s)) })   # S5's B (23) and NK (21) designs (5.3)
  NEED_KIDS <- MODE %in% c("known", "children")
  if (NEED_KIDS && !SMOKE) {   # 2.3: the 18 children's T columns and $cells rows (cell, sample, fine, group4, kept), subset at once
    KIDS_T <- PB$raw$T[, CHILD]; rownames(KIDS_T) <- G185
    KIDS_CELLS <- PB$cells[sample %in% CHILD, .(cell, sample, fine, group4, kept)]
  }
  rm(PB); invisible(gc(FALSE))
  if (NEED_KIDS && SMOKE) {   # synthetic children: the pooled healthy-adult profile at random depths; real IDs only
    set.seed(SMS[11]); p <- rowSums(G_RAW[, HEAL]) + 0.05; p <- p / sum(p)
    nk <- setNames(round(exp(runif(18, log(20), log(700)))), CHILD)
    KIDS_T <- sapply(CHILD, function(k) rnbinom(length(p), size = 20, mu = p * nk[[k]] * 1500)); rownames(KIDS_T) <- G185; KIDS_T <- syn(KIDS_T)
    KIDS_CELLS <- rbindlist(lapply(CHILD, function(k) data.table(cell = paste0(k, "_", seq_len(nk[[k]] + 20)), sample = k,
      fine = sample(TL, nk[[k]] + 20, TRUE), group4 = "T", kept = c(rep(TRUE, nk[[k]]), rep(FALSE, 20)))))
  }
  OOFT <- fread(file.path(TABWS, "ws_a1_g02_oof_scores.csv"))
  FG <- fread(file.path(TABWS, "ws_a1_g02_fold_genes.csv"))
  FOLD <- lapply(setNames(SAMP, SAMP), function(s) list(U = FG[left_out == s & direction == "up", gene], D = FG[left_out == s & direction == "down", gene]))
  HF <- OOFT[status == "healthy"][order(fold), sample]; stopifnot(length(HF) == 10L, setequal(HF, HEAL))   # folds 16-25
  L25 <- logzL(G_RAW)                                                     # A1-G02 R2's L
  MU <- lapply(setNames(SAMP, SAMP), function(s) rowMeans(L25[, setdiff(SAMP, s), drop = FALSE]))
  SG <- lapply(setNames(SAMP, SAMP), function(s) apply(L25[, setdiff(SAMP, s), drop = FALSE], 1, sd))
  fold_scores <- function(LX, folds) {   # columns of LX (on L25's rows) x folds: z with the fold's other 24, the fold's lists
    LX <- LX[rownames(L25), , drop = FALSE]
    out <- vapply(folds, function(s) apply(LX, 2, function(x) score_col(x, MU[[s]], SG[[s]], FOLD[[s]]$U, FOLD[[s]]$D)), numeric(ncol(LX)))
    matrix(out, ncol(LX), length(folds), dimnames = list(colnames(LX), folds)) }   # matrix(): vapply returns a vector when LX has one column
}

## -- 6. Chen2023 pseudobulks: --build (2.2) ------------------------------------------------------------------------------------------
CHEN_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", "Chen2023")
ZEN <- file.path(RAW_DIR, "public", "Zenodo", "10.17632_gwjh3w6ztm.2")
XLSX <- file.path(METADIR, "Chen2023", "meta_Chen2023_v2.xlsx")
REMAP <- file.path(TAB_DIR, "00_ingest", "GSE120221_gene_remap.csv")
if (MODE == "build") {
  feat <- fread(file.path(ZEN, "filtered_feature_bc_matrix.7z_AML103_CD34", "filtered_feature_bc_matrix", "features.tsv.gz"), header = FALSE,
                col.names = c("id", "symbol", "type"))
  RM <- fread(REMAP); stopifnot(nrow(RM) == nrow(feat))
  ROWN <- make.unique(gsub("_", "-", feat$symbol))                    # the objects' row names (K1)
  MAPPED <- RM$rule == "mapped_by_ensembl"; GN <- RM$final[MAPPED]       # 2.2: only mapped_by_ensembl rows, under their 2020-A names
  stopifnot(!anyDuplicated(GN))
  K1 <- list(feat_eq_remap = identical(feat$id, RM$ens) && identical(feat$symbol, RM$sym_cr2), n_rows = nrow(feat), n_mapped = sum(MAPPED), rownames_ok = list())
  xl <- as.data.table(readxl::read_excel(XLSX, sheet = "meta_Chen2023"))
  XL <- xl[, .(donor = sub("^Chen2023__", "", sample_id), arm = fifelse(disease == "AML", "AML", "donor"), age = as.numeric(age_years), sex = as.character(sex))]
  synth_object <- function(d) {   # smoke: a synthetic object with the real gene axis and meta.data columns, counts from the healthy-adult profile
    set.seed(SMS[match(d, DON$donor)])
    p0 <- rowSums(G_RAW[, HEAL]); p <- rep(0.02, length(ROWN)); hit <- match(RM$final, names(p0)); p[MAPPED & !is.na(hit)] <- p0[hit[MAPPED & !is.na(hit)]] + 0.02
    p <- p / sum(p); nA <- 120L; nB <- if (d == "AML469") 8L else sample(250:450, 1)
    fineB <- sample(c(TL, "Mature B", "NK", "Stromal", "CD14 Mono"), nB, TRUE, prob = c(rep(0.08, 9), 0.08, 0.06, 0.08, 0.06))
    fine <- c(sample(c("CD14 Mono", "GMP", TL[1]), nA, TRUE, prob = c(0.6, 0.35, 0.05)), fineB)
    md <- data.frame(Pool = rep(c(POOLA, POOLB), c(nA, nB)), Library = rep(paste0(d, c("_CD34", "_Niche_Immune")), c(nA, nB)), bmm_fine = fine,
                     anno_source = sample(c("concordant", "bmm_over_marker"), nA + nB, TRUE, prob = c(0.9, 0.1)), high_error = rbinom(nA + nB, 1, 0.03),
                     hierarchy_bin = ifelse(fine %in% TL | fine == "NK", "T_NK", ifelse(fine == "Mature B", "B_Plasma", ifelse(fine == "Stromal", "Stromal", "Mono_DC"))))
    rownames(md) <- paste0(md$Library, "_", seq_len(nA + nB))
    X <- Matrix(rpois(length(p) * (nA + nB), rep(p * 1500, nA + nB)), length(p), nA + nB, sparse = TRUE, dimnames = list(ROWN, rownames(md)))
    list(md = md, X = X) }
  rs <- function(X, j) if (any(j)) Matrix::rowSums(X[, j, drop = FALSE]) else numeric(nrow(X))
  OUT <- list()
  for (d in DON$donor) {
    if (SMOKE) { ob <- synth_object(d); md <- ob$md; X <- ob$X; rm(ob) } else {
      o <- readRDS(file.path(CHEN_DIR, paste0(d, ".rds")))
      md <- o@meta.data[, c("Pool", "Library", "bmm_fine", "anno_source", "high_error", "hierarchy_bin")]   # the allowlisted columns only (section 8)
      X <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts"); rm(o) }
    K1$rownames_ok[[d]] <- identical(rownames(X), ROWN)
    m <- MAP[match(as.character(md$bmm_fine), bmm_fine)]; g4 <- m$group4; fine <- as.character(md$bmm_fine)
    kept <- !is.na(g4) & g4 != "EXCLUDED" & !(WS_S4_DROP_HIGH_ERROR & md$high_error %in% 1L) & (!WS_S4_REQUIRE_CONCORDANT | md$anno_source %in% "concordant")
    pool <- as.character(md$Pool); lib <- as.character(md$Library); inB <- pool == POOLB; inA <- pool == POOLA
    TB <- kept & inB & g4 %in% "T"
    agree <- all(grepl("niche_immune$", lib[kept & inB], ignore.case = TRUE)) && all(grepl("_CD34$", lib[kept & inA])) && all(inA | inB)
    lab <- vapply(TL, function(k) rs(X, TB & fine == k), numeric(nrow(X)))   # genes x 9 labels (zero-filled)
    ncell <- vapply(TL, function(k) sum(TB & fine == k), 1L)
    sx <- function(j) c(xist = sum(X["XIST", j]), rps4y1 = sum(X["RPS4Y1", j]))
    dxT <- NULL; tot_dx <- NA_real_
    if (d %in% ELIG) {   # 4.7: ws/22's call (22:122), one unit = all cells of the pool-B library; summed over kept T, rounded later
      res <- celda::decontX(x = X[, inB, drop = FALSE], z = NULL, seed = SEED, verbose = FALSE)
      DXc <- res$decontXcounts; jj <- colnames(DXc) %in% colnames(X)[TB]; dxT <- rs(DXc, jj); tot_dx <- sum(DXc[MAPPED, jj]); rm(res, DXc) }
    tot <- c(T = sum(X[MAPPED, TB]), A_T = sum(X[MAPPED, kept & inA & g4 %in% "T"]), B = sum(X[MAPPED, kept & inB & g4 %in% "B"]),   # K2: computed
             NK = sum(X[MAPPED, kept & inB & g4 %in% "NK"]), dxT = tot_dx)                                                            # apart from rs()
    OUT[[d]] <- list(T = rs(X, TB)[MAPPED], lab = lab[MAPPED, , drop = FALSE], ncell = ncell, A_T = rs(X, kept & inA & g4 %in% "T")[MAPPED],
                     B = rs(X, kept & inB & g4 %in% "B")[MAPPED], NK = rs(X, kept & inB & g4 %in% "NK")[MAPPED], dxT = if (is.null(dxT)) NULL else dxT[MAPPED],
                     n = c(keptT_B = sum(TB), keptT_A = sum(kept & inA & g4 %in% "T"), keptB_B = sum(kept & inB & g4 %in% "B"), keptNK_B = sum(kept & inB & g4 %in% "NK"),
                           cells_B = sum(inB), myel_B = sum(inB & md$hierarchy_bin %in% MYEL_BINS)),
                     k2 = c(pool_agree = agree), tot = tot,
                     sexB = sx(kept & inB & g4 %in% c("T", "B", "NK")), sexA = sx(kept & inA))
    ts(sprintf("    %-7s pool B: %5d cells, kept T %5d, B %4d, NK %4d | pool A kept T %4d", d, sum(inB), sum(TB), sum(kept & inB & g4 %in% "B"),
               sum(kept & inB & g4 %in% "NK"), sum(kept & inA & g4 %in% "T")))
    rm(X, md); invisible(gc(FALSE))
  }
  mk <- function(f, ds = DON$donor) { M <- sapply(ds, function(d) f(OUT[[d]])); rownames(M) <- GN; if (SMOKE) syn(M) else M }
  CH <- list(genes = GN, xlsx = XL, k1 = K1, n = rbindlist(lapply(DON$donor, function(d) as.data.table(as.list(OUT[[d]]$n))[, donor := d])),
             k2 = rbindlist(lapply(DON$donor, function(d) as.data.table(as.list(OUT[[d]]$k2))[, donor := d])),
             tot = rbindlist(lapply(DON$donor, function(d) as.data.table(as.list(OUT[[d]]$tot))[, donor := d])),
             sex = rbindlist(lapply(DON$donor, function(d) data.table(donor = d, xist_B = OUT[[d]]$sexB[["xist"]], rps_B = OUT[[d]]$sexB[["rps4y1"]],
                                                                      xist_A = OUT[[d]]$sexA[["xist"]], rps_A = OUT[[d]]$sexA[["rps4y1"]]))),
             T = mk(function(o) o$T), A_T = mk(function(o) o$A_T), B = mk(function(o) o$B), NK = mk(function(o) o$NK), dxT = mk(function(o) o$dxT, ELIG),
             fine = lapply(setNames(TL, TL), function(k) mk(function(o) o$lab[, k], ELIG)),
             ncell = sapply(ELIG, function(d) OUT[[d]]$ncell), synthetic = SMOKE)
  saveRDS(CH, wg("chen_pb.rds")); ts("    written ", rg("chen_pb.rds"))
  quit(save = "no")
}
CH <- readRDS(rg("chen_pb.rds")); if (SMOKE != isTRUE(CH$synthetic)) stop("chen_pb.rds and this run disagree on smoke")
CT <- CH$T[, ELIG]                                                     # synthetic in a smoke run (checked on the line above)

## -- 7. known answers (section 6) ---------------------------------------------------------------------------------------------
if (MODE == "known") {
  ts("[K0] md5 of the declared inputs")
  L1 <- LARGE1_DIR
  IN <- data.table(path = c(file.path(TABWS, c("ws_s5_imprint_genes.csv", "ws_s8_t_imprint_ambient_prone.csv", "ws_s5_scores.csv", "ws_a1_g02_fold_genes.csv",
                                                 "ws_a1_g02_oof_scores.csv", "ws_a1_g02_lpo_pairs.csv", "ws_s9_yardstick.csv", "ws_s9_direction.csv",
                                                 "ws_s4_pseudobulk_meta_GSE120221.csv", "decontx/GSE185381/ws_s1_ambient_blacklist.csv", "ws_s4_pseudobulk_meta_GSE185381.csv")),
                            here::here("scripts", "ws", "s4_cell_groups.tsv"), REMAP,
                            file.path(L1, "ws_pseudobulk", c("GSE185381/ws_pb_GSE185381.rds", "GSE120221/ws_pb_GSE120221.rds")),
                            file.path(METADIR, c("GSE120221/donor_ages_GSE120221_v1.0.csv", "GSE120221/meta_GSE120221_v1.0.csv", "GSE185381/meta_GSE185381_v2.3.csv")),
                            file.path(L1, "reference", "gencode_GRCh38_gene_order.txt"),
                            "/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv", XLSX,
                            file.path(CHEN_DIR, paste0(DON$donor, ".rds"))),
                   md5 = c("977f58b06da076cc39a561ab7445b69d", "acb2b043049a4286e4052d04cea15feb", "58374fc60aabc53ba020c3149764fe8f", "ede3280052ab8dfcc547e60c341f334b",
                           "b1da65a5ddd4d050a1b0849a9aa1206d", "166623a1461ae59f0c0095e226788da4", "dcbea817a8eb6c68575486e3576dfce9", "a95d8caf3607b8a0b05083346b7301b6",
                           "11a282f9586be7f5d4090d93bd047f7d", "acebed4dbeed2ef6ecd5a60791e28fb9", "8533d147444e3b70983693dfd77c5beb",
                           "f5a5a0283b243e5efaf17b1ddfdda5f1", "d8d7c37eeb8cb85f5bd4b9d53d1f0acb",
                           "6c2fb685c49d4a23ef20a6842aa7aa67", "b4ad90bd5c587d48cc66e02a537c5079",
                           "178ced1c73e531b67172ed413b009c07", "222ff9a6e85b08e551cd1ad95395c076", "cf0c1ed37db4dabbac440ce26494f586",
                           "c0950d592084a8b713de34ab30a7ddf8", "5ee444a9284ae52b7afab1c8aeeced2c", "c2407360893ccf5261dff4e2144cd996",
                           "5c45be9c3e98789944c3592afd0c2ad8", "d4fe077315c91cd6ca240c3077e0cf88", "44866549e35fb6bbf50a8d1898859574", "ff64ce4ef252a71f3d1a2dbebae80b74",
                           "e31fd77fbab7b33b3ef87c3addbfa7ad", "29bcb07cc25e58869926ef7f0c4324b5", "5364a9c13432d342d36805c400e8ab93", "3876099a82edb40faeca04b7d7f95636",
                           "de8c3cc44f518bb59cc216ac40b4f085", "4b08d3dfb5aa1aece3043fa6c3e823c9"))
  stopifnot(nrow(IN) == 31L)
  if (SMOKE) IN <- IN[!grepl("Chen2023/|\\.rds$|xlsx$", path)]               # smoke: small files only
  got <- unname(tools::md5sum(IN$path)); bad0 <- is.na(got) | got != IN$md5    # a missing file gives NA: counted as a failure
  chk("K0", "PB", sprintf("md5 of %d declared files (2.1-2.3)", nrow(IN)), nrow(IN), sum(!bad0))
  if (any(bad0)) message("      differ or missing: ", paste(basename(IN$path[bad0]), collapse = ", "))
  if (!SMOKE) {
    fz <- Sys.glob(file.path(ZEN, "filtered_feature_bc_matrix.7z_*", "filtered_feature_bc_matrix", "features.tsv.gz"))
    zm <- vapply(fz, function(f) sub(" .*$", "", system(sprintf("zcat '%s' | md5sum", f), intern = TRUE)), "")
    chk("K0", "P", "features.tsv: md5 of the decompressed content, in 20 of 20 libraries", 20L, sum(zm == "16952af81ef735f42952b64c8402f69a"))
    lm <- sub(" .*$", "", system(sprintf("cd '%s' && md5sum GSM*_metadata.csv.gz | sort -k2 | md5sum", file.path(RAW_DIR, "public", "GEO", "GSE185381_RAW")), intern = TRUE))
    chk("K0", "B", "GEO author metadata: md5 of the sorted md5sum listing (2.3 recipe)", "c9bf8d1d53653aad8d65d1da7f9da4e0", lm)
  }

  ts("[K1-K2] Chen2023 axis and donor table")
  chk("K1", "P", "features.tsv (ID, symbol) equals the remap's (ens, sym_cr2), all rows in order", TRUE, CH$k1$feat_eq_remap)
  chk("K1", "P", "rows / mapped_by_ensembl rows", "33694/32351", sprintf("%d/%d", CH$k1$n_rows, CH$k1$n_mapped))
  chk("K1", "P", "every object's row names equal make.unique(gsub('_', '-', symbol))", 10L, sum(unlist(CH$k1$rownames_ok)))
  chk("K1", "P", "P370 on the remapped axis (genes / up / down)", "358/154/204",
      sprintf("%d/%d/%d", sum(P370$gene %in% CH$genes), sum(UP %in% CH$genes), sum(DN %in% CH$genes)))
  chk("K1", "P", "402 set on the remapped axis", 388L, sum(F402$gene %in% CH$genes))
  xl <- merge(DON[, .(donor, arm, age, sex)], CH$xlsx, by = "donor", suffixes = c("", "_xlsx"))
  chk("K2", "P", "donor table equals the xlsx (arm, age, sex)", 10L, xl[, sum(arm == arm_xlsx & age == age_xlsx & sex == sex_xlsx)])
  chk("K2", "P", "Pool agrees with the library name for every kept cell (donors)", 10L, sum(CH$k2$pool_agree))
  if (!SMOKE) chk("K2", "P", "pool-B kept T per donor equals 2.2 (10 donors)", DON$keptT_B, CH$n[match(DON$donor, donor), keptT_B])
  tt <- CH$tot[match(DON$donor, donor)]
  chk("K2", "P", "saved pseudobulks: column sums equal the cells' summed counts over the mapped rows (T, pool-A T, B, NK; 10 donors)", 0,
      max(abs(colSums(CH$T[, DON$donor]) - tt$T), abs(colSums(CH$A_T[, DON$donor]) - tt$A_T), abs(colSums(CH$B[, DON$donor]) - tt$B), abs(colSums(CH$NK[, DON$donor]) - tt$NK)), tol = 0)
  chk("K2", "P", "saved decontX pool-B T: column sums equal the summed decontX counts (9 donors, max |diff|)", 0,
      max(abs(colSums(CH$dxT[, ELIG]) - tt[match(ELIG, donor), dxT])), tol = 1e-6)
  lab_sum <- Reduce(`+`, lapply(CH$fine, function(M) M[, ELIG]))
  chk("K2", "P", "the nine label matrices sum to whole-T (max |diff|)", 0, max(abs(lab_sum - CT)), tol = 0)
  chk("K2", "P", "label cell numbers sum to the pool-B kept T", unname(CH$n[match(ELIG, donor), keptT_B]), unname(colSums(CH$ncell[, ELIG])))

  ts("[K3] expression sex (reported; the swap rule gates P)")
  sx <- merge(DON[, .(donor, arm, sex)], CH$sex, by = "donor")
  sx[, `:=`(call_B = sexcall(xist_B, rps_B), call_A = sexcall(xist_A, rps_A))]
  for (i in seq_len(nrow(sx))) rep_("K3", "P", sprintf("%s (%s, %s): pool B %s; pool A %s", sx$donor[i], sx$arm[i], sx$sex[i], sx$call_B[i], sx$call_A[i]),
                                     if (isTRUE(sx$call_B[i] == sx$sex[i])) "match" else "mismatch")
  opp <- sx[call_A %in% c("F", "M") & call_B %in% c("F", "M") & call_A != call_B, donor]
  bad <- sx[donor %in% ELIG & call_B %in% c("F", "M") & call_B != sex]
  swap <- nrow(bad[arm == "AML"]) > 0 && nrow(bad[arm == "donor"]) > 0 && any(outer(bad[arm == "AML", sex], bad[arm == "donor", sex], "!="))
  rep_("K3", "P", "pool-A vs pool-B opposite confident calls (library swap within a donor)", if (length(opp)) paste(opp, collapse = " ") else "none")
  rep_("K3", "P", "cross-arm swap pattern (an AML and a donor pool-B library each mismatched, exchange restores both)",
       if (swap) paste(bad$donor, collapse = " ") else "none", if (swap) "SWAP_STOP" else "REPORTED")

  ts("[K3b] children's identity (reported)")
  MG <- fread(file.path(METADIR, "GSE185381", "meta_GSE185381_v2.3.csv"), select = c("sample_id", "sex"), colClasses = "character")[, sid := sub("^GSE185381__", "", sample_id)]
  need_syn(KIDS_T, "child")
  kx <- data.table(child = CHILD, sex = MG$sex[match(CHILD, MG$sid)], xist = colSums(KIDS_T[c("XIST"), CHILD, drop = FALSE]), rps = colSums(KIDS_T[c("RPS4Y1"), CHILD, drop = FALSE]))
  kx[, call := sexcall(xist, rps)]; kx[, sexc := substr(toupper(sex), 1, 1)]
  kx[child == "AML003" & !sexc %in% c("F", "M"), sexc := "M"]          # 2.3 identity default: no Table S1 record; sex from author metadata
  if (SMOKE) kx[, sexc := sample(sexc)]                                  # smoke tables carry no real per-child quantity (section 8)
  for (i in seq_len(nrow(kx))) rep_("K3b", "B", sprintf("%s (metadata %s): %s%s", kx$child[i], kx$sexc[i], kx$call[i],   # the other-sex gene for the 4
                                     if (kx$child[i] %in% MIXED) { o <- if (kx$sexc[i] == "F") "RPS4Y1" else "XIST"; sprintf(" (%s %d UMI)", o, as.integer(if (o == "XIST") kx$xist[i] else kx$rps[i])) } else ""),
                                     if (isTRUE(kx$call[i] == kx$sexc[i])) "match" else "mismatch: dated note before (b) is read; child joins the flags")

  ts("[K4] the J-space code reproduces S9")
  CO <- load_cohort("GSE120221"); el <- CO$donors[eligible == TRUE, donor]
  J9 <- jspace(G_RAW, STAT, CO$T_don[, el, drop = FALSE])
  chk("K4", "P", "D_res_J for GSE120221 (S9 donor pooling)", DRES_S9, J9$d_res, tol = 1e-6)
  zs <- score(logzZ(G_RAW)$Z, UP, DN); gin <- mean(zs[AMLS]) - mean(zs[HEAL])
  chk("K4", "P", "gap_in on the 25 in S5 space", GAPIN_K4, gin, tol = 1e-5)
  rm(CO, J9); invisible(gc(FALSE))

  ts("[K5] enumeration, rank-sum, bounds (synthetic vectors only)")
  chk("K5", "P", "126 distinct labellings, the observed one once", "126/1", sprintf("%d/%d", ncol(unique(E54$L, MARGIN = 2)), sum(apply(E54$L, 2, function(s) all(s == 1:5)))))
  set.seed(SEED + OFF[["k6a"]]); xs <- replicate(400, rnorm(9) + c(rep(0.8, 5), rep(0, 4)))   # K5 test vectors (synthetic; SEED + 28, as K6a)
  pw <- vapply(seq_len(ncol(xs)), function(i) { x <- xs[, i]; U <- Ucount(x, 5L)
    Us <- apply(E54$L, 2, function(s) sum(outer(x[s], x[-s], ">"))); c(sum(Us >= U) / 126, pwilcox(U - 1, 5, 4, lower.tail = FALSE)) }, numeric(2))
  chk("K5", "P", "rank-sum exact p by enumeration equals pwilcox (400 synthetic vectors, max |diff|)", 0, max(abs(pw[1, ] - pw[2, ])), tol = 1e-12)
  bb <- vapply(seq_len(50), function(i) c(bounds(xs[, i]), bounds_bisect(xs[, i])), numeric(4))
  chk("K5", "P", "closed-form bounds equal bisection (50 synthetic vectors, max |diff|)", 0, max(abs(bb[1:2, ] - bb[3:4, ])), tol = 1e-8)
  agree <- vapply(seq_len(ncol(xs)), function(i) (bounds(xs[, i])[["lo"]] > 0) == (pexact(xs[, i], E54) <= M_PASS / 126), TRUE)
  chk("K5", "P", "rho_L > 0 exactly when p <= 6/126 (400 synthetic vectors)", 400L, sum(agree))

  ts("[K7-K9] the (b) machinery (adults; children's depth only)")
  st <- OOFT[match(SAMP, sample)]
  rep_oof <- vapply(SAMP, function(s) score_col(L25[, s], MU[[s]], SG[[s]], FOLD[[s]]$U, FOLD[[s]]$D), 1)
  chk("K7", "B", "stored out-of-fold scores of all 25 folds reproduced (max |diff|)", 0, max(abs(rep_oof - st$oof_score)), tol = 1e-6)
  L27 <- make_L27(G_RAW, G_RAW)                                          # every adult also as a new column (K7b)
  chk("K7", "B", "make_L27 keeps L on the 25 (max |diff| of the 25 columns and of the manual formula)", 0,
      max(abs(L27$L[, seq_along(SAMP)] - L25), abs(L27$check)), tol = 1e-10)
  nb <- vapply(seq_along(SAMP), function(i) score_col(L27$L[, 25L + i], MU[[SAMP[i]]], SG[[SAMP[i]]], FOLD[[SAMP[i]]]$U, FOLD[[SAMP[i]]]$D), 1)
  chk("K7b", "B", "each adult placed as a new column reproduces its stored out-of-fold score (25/25)", 25L, sum(abs(nb - st$oof_score) <= 1e-6))
  need_syn(KIDS_T, "child"); lib_k <- colSums(KIDS_T); tgt <- min(lib_k)
  set.seed(SEED + OFF[["k7c"]])
  TH <- sapply(HEAL, function(h) { x <- G_RAW[, h]; rbinom(length(x), x, tgt / sum(x)) }); rownames(TH) <- G185
  LT <- make_L27(G_RAW, cbind(TH, KIDS_T))$L[, -seq_along(SAMP), drop = FALSE]
  FS <- fold_scores(LT, HF)
  chk("K7c", "B", "healthy adults thinned to the smallest child's T library, and every child, give finite scores in all 10 healthy folds", TRUE, all(is.finite(FS)))
  sh <- vapply(HEAL, function(h) mean(FS[h, ] - fold_scores(L25[, h, drop = FALSE], HF)[1, ]), 1)
  rep_("K7c", "B", "mean score shift of the thinned healthy adults over their 10 healthy folds", paste(sprintf("%.3f", sh), collapse = " "))
  if (!SMOKE) {
    M4 <- fread(file.path(TABWS, "ws_s4_pseudobulk_meta_GSE185381.csv"), select = c("sample", "group4", "n", "eligible", "role", "paediatric"))
    m4 <- M4[group4 == "T" & paediatric == TRUE & eligible == TRUE & sample != "PAUMTZ"]
    chk("K8", "B", "the 18 children are S4's eligible paediatric T samples minus PAUMTZ", TRUE, setequal(m4$sample, CHILD))
    kc <- KIDS_CELLS[kept == TRUE & group4 == "T", .N, by = sample]
    chk("K8", "B", "their kept T counts equal the S4 meta table (18/18)", 18L, sum(m4$n[match(CHILD, m4$sample)] == kc$N[match(CHILD, kc$sample)]))
  } else rep_("K8", "B", "skipped in a smoke run (children's real counts)", "skipped")
  LP <- fread(file.path(TABWS, "ws_a1_g02_lpo_pairs.csv"))
  Sm <- as.matrix(dcast(LP, aml ~ healthy, value.var = "score_aml"), rownames = "aml"); Cm <- as.matrix(dcast(LP, aml ~ healthy, value.var = "score_healthy"), rownames = "aml")
  ab9 <- Ab_of(Sm, Cm); lo9 <- q7(grid_boot(Sm, Cm)[1, ])[1]          # the same grid_boot() as (b)
  chk("K9", "B", "the (b) code on A1-G02's leave-pair-out scores: A_b = 139/150", 139 / 150, ab9, tol = 1e-12)
  chk("K9", "B", "... and AML-side (2.5th percentile of A_b* > 0.5)", TRUE, lo9 > 0.5)

  ts("[K2b] est_C ignores the held-out donor's label counts (5.1)")
  ST9 <- setNames(fifelse(ARM == "AML", "AML", "healthy"), ELIG); LST9 <- lapply(setNames(ELIG, ELIG), function(i) list(U = UP, D = DN))
  G0 <- rownames(CT)[!rownames(CT) %in% black_for(rownames(CT)) & rowSums(CT) > 0]   # logzL drops the other rows anyway
  C0 <- lapply(CH$fine, function(M) M[G0, ELIG]); N0 <- CH$ncell[TL, ELIG]; W0 <- CT[G0, ]
  r0 <- est_C(ELIG, ST9, LST9, C0, N0, W0, refs = TRUE)$refs           # computed, never printed or written
  k2b <- unlist(mc_ok(mclapply(ELIG, function(i) { r2 <- est_C(ELIG, ST9, LST9, lapply(C0, function(M) { M[, i] <- 0; M }), N0, W0, refs = TRUE)$refs
    identical(r0[[paste(i, "H")]], r2[[paste(i, "H")]]) && identical(r0[[paste(i, "A")]], r2[[paste(i, "A")]]) }, mc.cores = NC), "K2b"))
  chk("K2b", "P51", "folds whose est_C references ignore the held-out donor's label counts (both orders)", 9L, sum(k2b))

  ts("[K6a] order invariance (real matrix; only the difference is printed) and plant signs (synthetic)")
  s0 <- score(jspace(G_RAW, STAT, CT)$ZE, UP, DN)[ELIG]; w0 <- score(logzZ(CT)$Z, UP, DN)[ELIG]; dif <- 0
  set.seed(SEED + OFF[["k6a"]])
  for (r in 1:10) { o <- sample(ELIG); dif <- max(dif, abs(score(jspace(G_RAW, STAT, CT[, o])$ZE, UP, DN)[ELIG] - s0), abs(score(logzZ(CT[, o])$Z, UP, DN)[ELIG] - w0)) }
  chk("K6a", "P", "J and W scores invariant to the order of the 9 columns (10 orders, max |diff|)", 0, dif, tol = 1e-12)
  rm(s0, w0)
  stopifnot(all(CH$genes %in% G185)); PH <- rowSums(G_RAW[CH$genes, HEAL]); PH <- PH / sum(PH)   # pooled healthy-adult T profile, Chen axis
  LFC <- setNames(IMP$logFC, IMP$gene)[intersect(P370$gene, CH$genes)]                    # S5 log2FC on P370
  plant <- function(cols, f) { mu <- outer(PH, rep(1e6, 9)); dimnames(mu) <- list(CH$genes, ELIG); mu[names(LFC), cols] <- mu[names(LFC), cols] * 2^(LFC * f); round(mu) }
  dA <- d54(score(jspace(G_RAW, STAT, plant(1:5, 1 / R_YARD))$ZE, UP, DN)[ELIG]); dH <- d54(score(jspace(G_RAW, STAT, plant(6:9, 1 / R_YARD))$ZE, UP, DN)[ELIG])
  chk("K6a", "P", "plant signs, synthetic (no noise): S5 log2FC / r into the 5 AML -> D > 0; into the 4 donors -> D < 0", "+/-",
      paste0(if (dA > 0) "+" else "-", "/", if (dH < 0) "-" else "+"))

  LIB <- colSums(CT)                                                       # realised pool-B kept-T library sizes, each kept with its donor
  detok <- all(vapply(CHK, function(d) !(d$kind == "deterministic" && d$scope %in% c("P", "PB") && d$state != "PASS"), TRUE))
  if (!detok) {   # section 6: the statistical bars run only after the deterministic checks pass
    for (b in c("K4b", "K6b(i)", "K6b(ii)")) CHK[[length(CHK) + 1L]] <- data.table(id = b, scope = "P", kind = "statistical",
      check = "not run: a scope-P deterministic check failed (section 6)", expected = "", got = "", state = "NOT_RUN")
    ts("    statistical bars not run: a scope-P deterministic check failed")
  } else {
  ts("[K4b] coverage of P370 in Chen2023's J space (bar)")
  JC <- jspace(G_RAW, STAT, CT)
  cov <- sum(P370$gene %in% JC$U_J) / 358; nu <- sum(UP %in% JC$U_J); nd <- sum(DN %in% JC$U_J)
  bar("K4b", sprintf("|P370 in U_J| / 358 (up %d, down %d); floor 0.80 and >= 5 per half", nu, nd), cov, cov >= COV_MIN && nu >= MIN_GENES && nd >= MIN_GENES)
  rm(JC)
  ts(sprintf("[K6b] machinery bars: %d synthetic data sets per arm (NB size 10, the donors' pool-B library sizes)", N_K6B))
  sim1 <- function(seed, pl) { set.seed(seed); mu <- outer(PH, LIB); if (pl) mu[names(LFC), 1:5] <- mu[names(LFC), 1:5] * 2^LFC
    E <- matrix(rnbinom(length(mu), size = 10, mu = mu), nrow(mu), dimnames = list(CH$genes, ELIG))
    x <- score(jspace(G_RAW, STAT, E)$ZE, UP, DN)[ELIG]; c(D = d54(x), pu = pexact(x, E54, "upper"), pl = pexact(x, E54, "lower")) }
  set.seed(SEED + OFF[["k6b"]]); SEEDS <- sample.int(.Machine$integer.max, 2L * N_K6B)
  R0 <- do.call(rbind, mc_ok(mclapply(seq_len(N_K6B), function(b) sim1(SEEDS[b], FALSE), mc.cores = NC), "K6b (i)"))
  R1 <- do.call(rbind, mc_ok(mclapply(seq_len(N_K6B), function(b) sim1(SEEDS[N_K6B + b], TRUE), mc.cores = NC), "K6b (ii)"))
  rr <- mean(R0[, "D"] > 0 & R0[, "pu"] <= M_PASS / 126); rv <- mean(R0[, "D"] < 0 & R0[, "pl"] <= M_PASS / 126); pp <- mean(R1[, "D"] > 0 & R1[, "pu"] <= M_PASS / 126)
  bar("K6b(i)", "no plant: 'replicated' rate (machinery, one tail; bar <= 0.08)", rr, rr <= 0.08)
  rep_("K6b(i)", "P", "no plant: 'reversed' rate (depth offset dominates; reported only)", sprintf("%.3f", rv))
  bar("K6b(ii)", "S5 log2FC planted on P370 in the 5 synthetic AML: pass rate (bar >= 0.90)", pp, pp >= 0.90)
  }

  ts(sprintf("[K6c] depth calibration curve (3 deepest healthy donors, %d draws per depth)", N_THIN))
  dep <- sort(colSums(G_RAW[, HEAL]), decreasing = TRUE); PRB <- names(dep)[1:3]; GRID <- c(2e6, 1e6, 5e5, 2.5e5, 1e5)
  chk("K6c", "reported", "the 3 deepest healthy T pseudobulks hold at least 3.0M UMI", 3L, sum(dep[1:3] >= 3e6))
  set.seed(SEED + OFF[["k6c"]]); SEEDC <- sample.int(.Machine$integer.max, 3L * N_THIN)
  job <- CJ(probe = PRB, b = seq_len(N_THIN))                               # CJ: every probe x draw combination
  SH <- do.call(rbind, mc_ok(mclapply(seq_len(nrow(job)), function(r) { set.seed(SEEDC[r]); pr <- job$probe[r]; x <- G_RAW[, pr]
    oth <- G_RAW[, setdiff(HEAL, pr)]; one <- function(v) { E <- cbind(v, oth); colnames(E)[1] <- pr; score(jspace(G_RAW, STAT, E)$ZE, UP, DN)[[pr]] }
    x3 <- rbinom(length(x), x, 3e6 / sum(x)); s3 <- one(x3)
    vapply(GRID, function(t) one(rbinom(length(x3), x3, t / sum(x3))) - s3, 1) }, mc.cores = NC), "K6c"))
  CURVE <- colMeans(SH)
  shift_at <- function(u) { gx <- log(c(3e6, GRID)); gy <- c(0, CURVE)   # zero at >= 3.0M; linear in log UMI; the 0.1M value below it
    ifelse(u >= 3e6, 0, ifelse(u <= 1e5, CURVE[length(CURVE)], approx(gx, gy, xout = log(u))$y)) }
  BIAS <- mean(shift_at(LIB[1:5])) - mean(shift_at(LIB[6:9]))
  rep_("K6c", "reported", "J-score shift at 2.0 / 1.0 / 0.5 / 0.25 / 0.1M UMI (mean of 3 probes)", paste(sprintf("%.4f", CURVE), collapse = " / "))
  rep_("K6c", "reported", "implied D bias at Chen2023's realised pool-B UMI totals (J units)", sprintf("%.6f", BIAS))

  put(rbindlist(CHK), "checks")
  ts(sprintf("known answers: %d rows; deterministic FAIL %d; statistical FAIL %d; swap stop %d", length(CHK),
             sum(vapply(CHK, function(d) d$state == "FAIL_deterministic", TRUE)), sum(vapply(CHK, function(d) d$state == "FAIL_statistical", TRUE)),
             sum(vapply(CHK, function(d) d$state == "SWAP_STOP", TRUE))))
  quit(save = "no")
}

## -- 8. (a): primary, sensitivities, readouts --------------------------------------------------------------------------------
stat_status <- function(C) {   # section 6: outputs of scope P are labelled when a bar failed (K6b(i) concerns only a would-be pass)
  f <- C[kind == "statistical" & state != "PASS", unique(id)]
  if (any(c("K4b", "K6b(ii)") %in% f)) "not decisive, descriptive" else if ("K6b(i)" %in% f) "a 'replicated' result would be not decisive" else "ok" }
if (MODE %in% c("primary", "sens", "readouts")) {
  if (!gate("P")) stop("section 6: a scope-P deterministic known answer has not passed, or K3's swap rule stops the run")
  JC <- jspace(G_RAW, STAT, CT); SC <- score(JC$ZE, UP, DN)[ELIG]
  sp <- sqrt((4 * var(SC[1:5]) + 3 * var(SC[6:9])) / 7)                     # pooled within-arm SD (7 df; descriptive)
}
if (MODE == "primary") {
  ts("[3.1-3.2] primary")
  C <- checks(); bias <- as.numeric(C[id == "K6c" & grepl("^implied D bias", check), got])
  if (any(C[id == "K3", state] == "mismatch") && !NOTED[["k3"]]) stop("K3: a sex mismatch needs a dated note in the pre-registration first; then rerun with --k3-noted (section 6)")
  PR <- a_row("primary: P370, J space", SC, JC)
  iu <- influence(SC, 5L, 1); id_ <- influence(SC, 5L, -1)
  PR[, `:=`(n_up = sum(UP %in% JC$U_J), n_down = sum(DN %in% JC$U_J), D_over_pooled_sd = D / sp, depth_bias_J = bias, depth_bias_sd = bias / sp,
            qual_replicated = iu$flag, qual_replicated_donor = iu$named, qual_reversed = id_$flag, qual_reversed_donor = id_$named,
            statistical_failures = { f <- C[kind == "statistical" & state != "PASS", unique(id)]; if (length(f)) paste(f, collapse = ";") else "none" },
            status = stat_status(C))]
  amb <- intersect(c("LYZ", "S100A8", "S100A9"), rownames(CT))
  DONR <- data.table(donor = ELIG, arm = unname(ARM), score = unname(SC), keptT = CH$n[match(ELIG, donor), keptT_B], umi = unname(colSums(CT)),
                     genes_detected = unname(colSums(CT > 0)), myeloid_share_poolB = CH$n[match(ELIG, donor), myel_B / cells_B],
                     ambient_cpm = unname(colSums(CT[amb, , drop = FALSE]) / colSums(CT) * 1e6), D_minus = unname(iu$D_minus[ELIG]))
  LOO <- rbindlist(lapply(seq_along(ELIG), function(i) { x <- SC[-i]; E <- if (i <= 5) E44 else E53
    data.table(dropped = ELIG[i], arm = ARM[[i]], D_minus = mean(x[seq_len(E$na)]) - mean(x[-seq_len(E$na)]),
               p_upper = pexact(x, E, "upper"), p_lower = pexact(x, E, "lower"), n_labellings = ncol(E$L)) }))
  print(PR); print(DONR); print(LOO)
  put(rbind(cbind(row = "summary", PR), cbind(row = "donor", DONR), fill = TRUE), "primary"); put(LOO, "loo")
}
if (MODE == "sens") {
  ts("[4.1-4.7] sensitivities")
  SVALL <- SV$gene; out <- list()
  out[["4.1"]] <- a_row("4.1 the 402 set", score(JC$ZE, UP402, DN402)[ELIG], JC, UP402, DN402)
  out[["4.2"]] <- a_row("4.2 W space", score(logzZ(CT)$Z, UP, DN)[ELIG])
  CA <- CT + CH$A_T[, ELIG]; JA <- jspace(G_RAW, STAT, CA); out[["4.3"]] <- a_row("4.3 all pools (pool-A kept T added)", score(JA$ZE, UP, DN)[ELIG], JA)
  out[["4.4u"]] <- a_row("4.4 up half", score(JC$ZE, UP, character(0))[ELIG], JC, UP, character(0))
  out[["4.4d"]] <- a_row("4.4 down half", score(JC$ZE, character(0), DN)[ELIG], JC, character(0), DN)
  U5 <- setdiff(UP, SVALL); D5_ <- setdiff(DN, SVALL)
  out[["4.5"]] <- a_row(sprintf("4.5 P370 without the Savage genes (%d up / %d down left)", length(U5), length(D5_)), score(JC$ZE, U5, D5_)[ELIG], JC, U5, D5_)
  tgt <- min(colSums(CT)); set.seed(SEED + OFF[["thin"]])                   # 4.6: every draw recomputes U_J, TMM and D_res_J (and W's filter)
  TH <- lapply(seq_len(N_THIN), function(b) { E <- apply(CT, 2, function(x) rbinom(length(x), x, tgt / sum(x))); dimnames(E) <- dimnames(CT)
    J <- jspace(G_RAW, STAT, E); list(J = score(J$ZE, UP, DN)[ELIG], dr = J$d_res, W = score(logzZ(E)$Z, UP, DN)[ELIG]) })
  xJ <- rowMeans(sapply(TH, `[[`, "J")); xW <- rowMeans(sapply(TH, `[[`, "W"))
  out[["4.6J"]] <- a_row(sprintf("4.6 depth-equalised, J (%d draws, %.0f UMI each)", N_THIN, tgt), xJ, d_res = mean(vapply(TH, `[[`, 1, "dr")))
  out[["4.6W"]] <- a_row("4.6 depth-equalised, W", xW)
  DX <- round(CH$dxT[, ELIG]); JD <- jspace(G_RAW, STAT, DX); out[["4.7"]] <- a_row("4.7 decontX per pool-B library", score(JD$ZE, UP, DN)[ELIG], JD)
  S <- rbindlist(out); S[, status := stat_status(checks())]; print(S); put(S, "sens")
}
if (MODE == "readouts") {
  ts("[5.1-5.4] readouts"); R <- list()
  nc <- CH$ncell[TL, ELIG]
  R[["5.1i"]] <- data.table(readout = "5.1 naive share", donor = ELIG, value = colSums(nc[NAIVE, ]) / colSums(nc))
  R[["5.1i_cd8"]] <- data.table(readout = "5.1 CD8 share", donor = ELIG, value = colSums(nc[CD8L, ]) / colSums(nc))
  R[["5.1i_arm"]] <- rbind(R[["5.1i"]], R[["5.1i_cd8"]])[, .(donor = paste("median", ARM[donor][1]), value = median(value)), by = .(readout, arm = ARM[donor])][, arm := NULL]
  C <- checks()
  if (all(C[id == "K2b", state] == "PASS")) {
    ST9 <- setNames(fifelse(ARM == "AML", "AML", "healthy"), ELIG); LST9 <- lapply(setNames(ELIG, ELIG), function(i) list(U = UP, D = DN))
    G0 <- rownames(CT)[!rownames(CT) %in% black_for(rownames(CT)) & rowSums(CT) > 0]
    ec <- est_C(ELIG, ST9, LST9, lapply(CH$fine, function(M) M[G0, ELIG]), nc, CT[G0, ])
    R[["5.1ii"]] <- data.table(readout = "5.1 counterfactual (u-units)", donor = "all",
      quantity = c("gapT", "gapC_H", "gapC_A", "residual_H", "residual_A", "residual_mean"),
      value = with(ec, c(gapT, gapC_H, gapC_A, gapT - gapC_H, gapT - gapC_A, gapT - (gapC_H + gapC_A) / 2)))
    saveRDS(ec[c("C_H", "C_A", "C", "gapT")], wg("readouts_C.rds"))   # 5.1: C_H, C_A and C are reported only under "direction replicated" (--combine)
    fb <- ec$fallback[fallback == TRUE]
    R[["5.1fb"]] <- data.table(readout = "5.1 fallbacks (expected: H CD4 Regulatory in the 4 donor folds; A CD4 Regulatory AML162, CD8 Tissue Resident Memory AML304)",
                               donor = fb$fold, quantity = paste(fb$order, fb$label), value = 1)
  } else R[["5.1ii"]] <- data.table(readout = "5.1 counterfactual not run: K2b did not pass", donor = "all", value = NA_real_)
  xd <- score(JC$ZE, SVU, SVD)[ELIG]                                         # 5.2 delay score
  cen <- function(v) v - ave(v, ARM)                                          # ave(): each arm's mean, repeated per donor
  R[["5.2"]] <- data.table(readout = "5.2 delay score", donor = "all", quantity = c("D_delay", "p_upper (AML higher)", "p_lower", "spearman_centred", "n_genes"),
                           value = c(d54(xd), pexact(xd, E54, "upper"), pexact(xd, E54, "lower"), cor(cen(SC), cen(xd), method = "spearman"),
                                     sum(c(SVU, SVD) %in% JC$U_J)))
  R[["5.2d"]] <- data.table(readout = "5.2 delay score per donor", donor = ELIG, value = unname(xd))
  for (cmp in c("B", "NK")) {                                                 # 5.3: each compartment's own GSE185381 reference
    n_c <- CH$n[match(ELIG, donor), get(paste0("kept", cmp, "_B"))]; el <- ELIG[n_c >= MIN_T]
    na <- sum(ARM[el] == "AML"); nh <- sum(ARM[el] == "donor"); E <- if (na == 5L && nh == 4L) E54 else enum(na, nh)
    J <- jspace(REF[[cmp]]$M, REF[[cmp]]$st, CH[[cmp]][, el, drop = FALSE]); x <- score(J$ZE, UP, DN)[el]
    nu <- sum(UP %in% J$U_J); nd <- sum(DN %in% J$U_J); ok5 <- nu >= MIN_GENES && nd >= MIN_GENES   # the 5-gene rule (5.3)
    dT <- mean(SC[el][seq_len(na)]) - mean(SC[el][-seq_len(na)]); dC <- if (ok5) mean(x[seq_len(na)]) - mean(x[-seq_len(na)]) else NA_real_
    R[[paste0("5.3", cmp)]] <- data.table(readout = paste("5.3", cmp), donor = paste(el, collapse = ";"),
      quantity = c("n_aml", "n_donor", "D", "p_upper (AML higher)", "coverage P370 in U_J", "n_up", "n_down", "D_T same donors", "D / D_T"),
      value = c(na, nh, dC, if (ok5) pexact(x, E, "upper") else NA_real_, sum(P370$gene %in% J$U_J) / 358, nu, nd, dT, dC / dT)) }
  g5 <- intersect(P370$gene, JC$U_J); LE <- JC$LE[g5, ELIG]; sg5 <- P370[match(g5, gene), s]   # 5.4 gene-level sign concordance
  sh <- apply(E54$L, 2, function(s) mean(sign(rowMeans(LE[, s, drop = FALSE]) - rowMeans(LE[, -s, drop = FALSE])) == sg5))
  R[["5.4"]] <- data.table(readout = "5.4 sign concordance (P370 in U_J)", donor = "all", quantity = c("share", "p (labellings with share >= observed)", "n_genes"),
                           value = c(sh[E54$obs], sum(sh >= sh[E54$obs] - 1e-12) / 126, length(g5)))
  RR <- rbindlist(R, fill = TRUE); RR[, status := stat_status(checks())]; print(RR); put(RR, "readouts")
}
## -- 9. (b): the children (3.4, 4.8-4.11, 5.5, 5.6) -------------------------------------------------------------------------------
if (MODE == "children") {
  if (!gate("B")) stop("section 6: a scope-B deterministic known answer has not passed")
  need_syn(KIDS_T, "child"); C <- checks()
  k3b_bad <- sub(" .*$", "", C[id == "K3b" & state != "match", check])        # K3b: a mismatched or "unclear" child joins the 2.3 flags
  if (length(k3b_bad) && !SMOKE && !NOTED[["k3b"]]) stop("K3b: ", paste(k3b_bad, collapse = " "), " need a dated note first; then rerun with --k3b-noted (section 6)")
  FLAGFREE_RUN <- setdiff(FLAGFREE, k3b_bad)
  ts("[3.4] placement, A_b and phi")
  LK <- make_L27(G_RAW, KIDS_T)$L[, CHILD, drop = FALSE]
  S25 <- fold_scores(LK, SAMP); SH <- S25[, HF, drop = FALSE]                 # children x 25 folds; the 10 healthy folds
  OO <- setNames(OOFT$oof_score, OOFT$sample); Am <- mean(OO[AMLS])
  CMPm <- matrix(OO[HF], nrow(SH), length(HF), byrow = TRUE, dimnames = dimnames(SH))   # comparator of each fold, repeated per child
  phi_of <- function(ri, hj) { hm <- mean(OO[HF[hj]])                        # a drawn fold keeps its comparator and its column of child scores;
    (mean(rowMeans(S25[ri, c(AMLS, HF[hj]), drop = FALSE])) - hm) / (Am - hm) } # conditional on the 15 AML folds
  two_way <- function(rows, folds = seq_along(HF)) {   # 3.4: grid_boot() (K9's code), with phi on the same draws
    if (length(rows) < 2L) return(list(Ab = NA_real_, Ab_ci = c(NA_real_, NA_real_), phi = NA_real_, phi_ci = c(NA_real_, NA_real_), n = length(rows), k = length(folds)))
    S <- SH[rows, folds, drop = FALSE]; Cc <- CMPm[rows, folds, drop = FALSE]
    B <- grid_boot(S, Cc, extra = function(i, j) phi_of(rows[i], folds[j]))
    list(Ab = Ab_of(S, Cc), Ab_ci = q7(B[1, ]), phi = phi_of(rows, folds), phi_ci = q7(B[2, ]), n = length(rows), k = length(folds)) }
  kc <- KIDS_CELLS[kept == TRUE & group4 == "T", .N, by = sample]; nT <- setNames(kc$N[match(CHILD, kc$sample)], CHILD)
  row_of <- function(sec, w, r) data.table(section = sec, item = c("A_b", "phi"), value = c(r$Ab, r$phi), lo = c(r$Ab_ci[1], r$phi_ci[1]),
                                           hi = c(r$Ab_ci[2], r$phi_ci[2]), n_children = r$n, n_folds = r$k, which = w)
  ALL <- two_way(seq_along(CHILD)); OUT <- list(row_of("3.4 primary", "18 children x 10 healthy folds", ALL))
  OUT[[length(OUT) + 1L]] <- data.table(section = "3.4 primary", item = c("category input: 2.5th percentile of A_b* > 0.5", "180 * A_b > 131 (exact Wilcoxon 2.5 % point)"),
                                        value = c(ALL$Ab_ci[1] > 0.5, 180 * ALL$Ab > 131))
  big <- which(nT >= 118); OUT[[length(OUT) + 1L]] <- row_of("3.4 phi, children with >= 118 kept T", paste(CHILD[big], collapse = ";"), two_way(big))[item == "phi"]
  for (x in list(list("4.8 under 18 y", match(UNDER18, CHILD)), list("4.9 >= 50 kept T", which(nT >= 50)), list("4.10 flag-free (2.3 and K3b)", match(FLAGFREE_RUN, CHILD))))
    OUT[[length(OUT) + 1L]] <- row_of(x[[1]], paste(CHILD[x[[2]]], collapse = ";"), two_way(x[[2]]))
  OUT5 <- c("Control0004", "Control0005", "Control0058", "Control0082", "Control4003"); OUT7 <- setdiff(HF, c("Control1", "Control3", "Control4"))
  for (x in list(list("4.11 folds outside P02/P03 (no CI)", OUT5), list("4.11 folds without Control1, Control3, Control4 (no CI)", OUT7)))
    OUT[[length(OUT) + 1L]] <- data.table(section = x[[1]], item = "A_b", value = Ab_of(SH[, x[[2]], drop = FALSE], CMPm[, x[[2]], drop = FALSE]),
                                          n_children = 18L, n_folds = length(x[[2]]), which = paste(x[[2]], collapse = ";"))
  sp_ <- c(AML022 = "Control0082", AML052 = "Control0082", AML009 = "Control0005")   # same-pool differences, no test
  OUT[[length(OUT) + 1L]] <- data.table(section = "3.4 same-pool differences (child minus comparator, comparator's fold)", item = paste(names(sp_), "vs", sp_),
                                        value = vapply(names(sp_), function(k) SH[k, sp_[[k]]] - OO[[sp_[[k]]]], 1))
  ts("[3.4] per-child descriptors (metadata and depth)")
  KC <- copy(KIDS_CELLS)[, key := sub("_([ACGT]+)(-[0-9]+)?$", ":\\1", sub("^GSM[0-9]+_", "", cell))]   # the critic's join key (review 1, c02)
  gf <- Sys.glob(file.path(RAW_DIR, "public", "GEO", "GSE185381_RAW", "GSM*_metadata.csv.gz"))
  gf <- gf[vapply(gf, function(f) "inflammation_group" %in% names(fread(f, nrows = 0)), TRUE)]
  AU <- rbindlist(lapply(gf, function(f) fread(f, select = c("cell", "samples", "Broad_cell_identity", "malignant", "age_group", "inflammation_group"))))[samples %in% CHILD]
  setnames(AU, "cell", "key"); AU <- unique(AU, by = "key")
  KT <- merge(KC[kept == TRUE & group4 == "T"], AU[, .(key, Broad_cell_identity, malignant)], by = "key", all.x = TRUE)
  grp <- AU[, .(inflammation = names(which.max(table(inflammation_group))), age_band = names(which.max(table(age_group)))), by = samples]
  if (SMOKE) grp[, `:=`(inflammation = sample(inflammation), age_band = sample(age_band))]   # smoke tables carry no real per-child quantity
  PC <- data.table(section = "3.4 per child", child = CHILD, score_10_healthy_folds = rowMeans(SH)[CHILD], score_25_folds = rowMeans(S25)[CHILD],
                   kept_T = nT[CHILD], T_library = colSums(KIDS_T)[CHILD],
                   kept_over_T_mapped = KIDS_CELLS[group4 == "T", .(f = sum(kept) / .N), by = sample][match(CHILD, sample), f],
                   naive_share = KT[, .(f = mean(fine %in% NAIVE)), by = sample][match(CHILD, sample), f],
                   CD8_share = KT[, .(f = mean(fine %in% CD8L)), by = sample][match(CHILD, sample), f],
                   author_NK_share = KT[, .(f = mean(Broad_cell_identity %in% "NK")), by = sample][match(CHILD, sample), f],
                   author_malignant_share = KT[, .(f = mean(malignant %in% "malignant")), by = sample][match(CHILD, sample), f],
                   author_matched = KT[, .(f = mean(!is.na(Broad_cell_identity))), by = sample][match(CHILD, sample), f],
                   inflammation = grp$inflammation[match(CHILD, grp$samples)], age_band = grp$age_band[match(CHILD, grp$samples)],
                   under18 = CHILD %in% UNDER18, flag_free = CHILD %in% FLAGFREE_RUN)
  if (!SMOKE) {   # the author join was matched 18,322 of 18,322 kept T cells in review 1; two counts disclosed in 1.3 are known answers here
    if (min(PC$author_matched) < 0.95) stop("author metadata join: matched share below 0.95 for ", paste(PC[author_matched < 0.95, child], collapse = " "))
    k28 <- KT[sample == "AML028", sum(Broad_cell_identity %in% "NK")]; k73 <- KT[sample == "AML073", sum(malignant %in% "malignant")]
    if (k28 != 16L || k73 != 210L) stop(sprintf("author metadata join: AML028 author-NK %d (1.3: 16), AML073 author-malignant %d (1.3: 210)", k28, k73)) }
  ts("[5.5] strata")
  ST <- rbind(PC[, .(section = "5.5 strata", item = paste("inflammation:", inflammation), value = median(score_10_healthy_folds), n_children = .N), by = inflammation][, inflammation := NULL],
              PC[, .(section = "5.5 strata", item = paste("age band:", age_band), value = median(score_10_healthy_folds), n_children = .N), by = age_band][, age_band := NULL])
  ts(sprintf("[5.6] depth readout: adults thinned to each child's T library (%d draws, SEED + 26)", N_THIN))
  libA <- colSums(G_RAW); libK <- colSums(KIDS_T)
  set.seed(SEED + OFF[["depth"]]); SEEDD <- sample.int(.Machine$integer.max, length(CHILD) * N_THIN)
  D56 <- rbindlist(mc_ok(mclapply(seq_along(CHILD), function(ci) { k <- CHILD[ci]; el <- SAMP[libA >= libK[[k]]]
    if (length(el) < 5L) return(data.table(child = k, n_adults = length(el), alpha = NA_real_, beta = NA_real_, noise_sd = NA_real_))
    sh <- rbindlist(lapply(seq_len(N_THIN), function(b) { set.seed(SEEDD[(ci - 1L) * N_THIN + b])
      TH <- vapply(el, function(a) { x <- G_RAW[, a]; rbinom(length(x), x, libK[[k]] / sum(x)) }, numeric(nrow(G_RAW)))
      dimnames(TH) <- list(G185, paste0(el, "@thin"))                    # new names: make_L27 binds them next to the 25 originals
      L <- make_L27(G_RAW, TH)$L[, -seq_along(SAMP), drop = FALSE]
      data.table(adult = el, full = OO[el], thin = vapply(seq_along(el), function(j) score_col(L[, j], MU[[el[j]]], SG[[el[j]]], FOLD[[el[j]]]$U, FOLD[[el[j]]]$D), 1)) }))
    f <- lm(I(thin - full) ~ full, data = sh)
    data.table(child = k, n_adults = length(el), alpha = coef(f)[[1]], beta = coef(f)[[2]], noise_sd = sd(resid(f))) }, mc.cores = NC), "5.6"))
  ZS <- rbindlist(lapply(HF, function(h) { g <- intersect(c(FOLD[[h]]$U, FOLD[[h]]$D), rownames(L25))
    data.table(section = "5.6 zero-count share of the fold's genes", item = h, child = CHILD, value = colMeans(KIDS_T[g, CHILD, drop = FALSE] == 0)) }))
  OUT[[length(OUT) + 1L]] <- data.table(section = "4.10 children removed by K3b", item = if (length(k3b_bad)) paste(k3b_bad, collapse = ";") else "none")
  RES <- rbindlist(c(OUT, list(PC, ST, cbind(section = "5.6 depth readout", D56), ZS)), fill = TRUE)
  print(rbindlist(OUT, fill = TRUE)); print(PC); print(D56)
  put(RES, "children")
}

## -- 10. section 7: categories, summary, GATES block ----------------------------------------------------------------------------------
if (MODE == "combine") {
  C <- checks(); has <- function(x) file.exists(file.path(ODIR, paste0("ws_a1_g04_", x, SFX, ".csv")))
  passed <- function(x) all(C[id == x & kind == "statistical", state] == "PASS")   # statistical rows only (K6b(i) also has a reported row)
  detP <- all(C[scope %in% c("P", "PB") & kind == "deterministic", state] == "PASS"); detB <- all(C[scope %in% c("B", "PB") & kind == "deterministic", state] == "PASS")
  swap <- any(C$state == "SWAP_STOP"); m <- M_PASS / 126; kk <- function(p) round(p * 126)
  ps <- if (detP && !swap && has("primary")) rtab("primary")[row == "summary"] else NULL      # a gate that stopped --primary leaves no table
  rep_ok <- !is.null(ps) && isTRUE(ps$D > 0 && ps$p_upper <= m + 1e-12); rev_ok <- !is.null(ps) && isTRUE(ps$D < 0 && ps$p_lower <= m + 1e-12)
  cat_a <- if (is.null(ps) || !passed("K4b") || !passed("K6b(ii)") || (!passed("K6b(i)") && rep_ok)) "not decisive" else
           if (rep_ok) "direction replicated" else if (rev_ok) "reversed" else "not shown"   # section 7(a), first match wins
  qual <- if (cat_a == "direction replicated" && isTRUE(ps$qual_replicated)) sprintf("rests on one donor (%s)", ps$qual_replicated_donor) else
          if (cat_a == "reversed" && isTRUE(ps$qual_reversed)) sprintf("rests on one donor (%s)", ps$qual_reversed_donor) else ""
  S <- if (has("sens")) rtab("sens") else NULL; R <- if (has("readouts")) rtab("readouts") else NULL
  s46 <- if (!is.null(S)) S[grepl("^4.6 depth-equalised, J", analysis)] else NULL
  dly <- !is.null(R) && { dl <- R[readout == "5.2 delay score"]; isTRUE(dl[quantity == "D_delay", value] > 0 && dl[quantity == "p_upper (AML higher)", value] <= m + 1e-12) }
  wording_a <- switch(cat_a,   # section 0's fixed sentences
    "direction replicated" = paste0("\"In a second study whose arms were cryopreserved, FACS-sorted and sequenced under one published protocol and are age-matched ",
      "(Chen2023: 5 NPM1-mutated AML at diagnosis vs 4 transplant donors; pool-B T cells), AML T cells score higher on the frozen imprint than donor T cells ",
      sprintf("(exact p = %d/126; rho %.2f, one-sided 95.2%% bounds %.2f to %.2f).\" ", kk(ps$p_upper), ps$rho, ps$rho_L, ps$rho_U),
      "This is evidence that the direction is not specific to GSE185381's procurement and processing (commercial vs tissue-bank marrow, its capture dates ",
      "and co-capture pools, T cells captured together with blasts). It cannot exclude a difference between donor and diagnostic marrow that both studies ",
      "share (for example blood admixture, A1-BATCH 7), pre-freeze handling and storage time, capture and sort sessions (all unrecorded in Chen2023), or the ",
      "clinical correlates of diagnosis (handoff section 14). A1-BATCH 12.5's may-not list applies unchanged.",
      if (nzchar(qual)) paste0(" The result ", qual, ".") else "",
      if (dly) " The processing-delay program also differs between the arms in Chen2023, so pre-freeze handling is not excluded there either." else ""),
    "not shown" = paste0(sprintf("\"The direction was not shown in Chen2023 (exact p = %d/126); a common shift larger than %.2f of GSE185381's selection-corrected gap ",
      kk(ps$p_upper), ps$rho_U), "is excluded (one-sided 95.2%, common-shift model).\"", if (ps$rho_U >= 1) " \"5 vs 4 cannot exclude an effect of the GSE185381 size.\"" else "",
      " The direction stays \"shown in one study\"."),
    "reversed" = paste0(sprintf("\"AML T cells score lower than donor T cells in Chen2023 (exact p = %d/126).\" ", kk(ps$p_lower)),
      "The direction does not reproduce in Chen2023 (NPM1-mutated AML vs allogeneic-transplant donors, one sorting protocol); study-, genotype- and ",
      "source-specific explanations cannot be told apart.", if (nzchar(qual)) paste0(" The result ", qual, ".") else "",
      if (!is.null(s46)) sprintf(" Depth-equalised (4.6): D %.3f, mirror p %d/126.", s46$D, kk(s46$p_lower)) else ""),
    "not decisive" = "A known answer could not be passed, a statistical bar or the coverage floor failed, or K3's swap rule applies; results are descriptive and the direction stays \"shown in one study\".")
  if (detB && has("children")) { K <- rtab("children"); ab <- K[section == "3.4 primary" & item == "A_b"]; ph <- K[section == "3.4 primary" & item == "phi"] } else ab <- ph <- NULL
  cat_b <- if (is.null(ab)) "not decisive" else if (isTRUE(ab$lo > 0.5)) "AML-side" else "not readable"
  wording_b <- switch(cat_b,
    "AML-side" = paste0(sprintf("\"Within GSE185381, T cells of 18 paediatric-cohort AML patients (14 under 18 y; 0.2-21 y, two without a recorded age) score higher on the cross-fitted imprint lists than out-of-fold adult healthy donors more often than not (A_b %.2f, 95%% CI %.2f to %.2f; placement phi %.2f, 95%% CI %.2f to %.2f).\" ",
      ab$value, ab$lo, ab$hi, ph$value, ph$lo, ph$hi), "The patients are younger than the donors, so an older-age effect cannot explain it. It does not separate disease from sample source, ",
      "from a shift shared by the children's lanes, from blast-derived ambient RNA, from T-subset composition, or from paediatric T-cell biology."),
    "not readable" = "No statement; A_b, phi and their CIs go to G03.", "not decisive" = "A scope-B known answer could not be passed; descriptive only.")
  CC <- if (cat_a == "direction replicated" && file.exists(rg("readouts_C.rds"))) readRDS(rg("readouts_C.rds")) else NULL   # 5.1: C only here
  gate_a3 <- cat_a == "direction replicated" && !nzchar(qual)                  # decision 9.6, option A: only an unflagged pass counts
  md5_of <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g04_", x, SFX, ".csv")); if (!file.exists(f)) "absent" else paste(unique(fread(f, select = "script_md5")$script_md5), collapse = ";") }
  f3 <- function(x) if (is.null(x) || !length(x) || is.na(x)) "not run" else sprintf("%.3f", x)
  SUM <- data.table(item = c("(a) category", "(a) qualifier", "(a) wording", "(a) D", "(a) exact p (upper / mirror)", "(a) rho [rho_L, rho_U]", "(a) depth-equalised 4.6 (J)",
                             "(a) 5.1 C_H / C_A / C (only under 'direction replicated', gapT > 0)",
                             "(b) category", "(b) wording", "(b) A_b [95% CI]", "(b) phi [95% CI]", "GATE_A3 direction half (decision 9.6 A)", "statistical failures",
                             "script md5: checks / primary / sens / readouts / children / summary"),
                    value = c(cat_a, if (nzchar(qual)) qual else "none", wording_a, f3(ps$D),
                              if (is.null(ps)) "not run" else sprintf("%d/126 / %d/126", kk(ps$p_upper), kk(ps$p_lower)),
                              if (is.null(ps)) "not run" else sprintf("%.3f [%.3f, %.3f]", ps$rho, ps$rho_L, ps$rho_U),
                              if (is.null(s46)) "not run" else sprintf("D %.4f, p %d/126", s46$D, kk(s46$p_upper)),
                              if (is.null(CC) || !isTRUE(CC$gapT > 0)) "not reported" else sprintf("%.3f / %.3f / %.3f", CC$C_H, CC$C_A, CC$C),
                              cat_b, wording_b, if (is.null(ab)) "not run" else sprintf("%.3f [%.3f, %.3f]", ab$value, ab$lo, ab$hi),
                              if (is.null(ph)) "not run" else sprintf("%.3f [%.3f, %.3f]", ph$value, ph$lo, ph$hi),
                              if (gate_a3) "counts" else "does not count",
                              { f <- C[kind == "statistical" & state != "PASS", unique(id)]; if (length(f)) paste(f, collapse = ";") else "none" },
                              paste(c(vapply(c("checks", "primary", "sens", "readouts", "children"), md5_of, ""), MD5), collapse = " / ")))
  blk <- c("<!-- A1-G04 start -->", sprintf("## T imprint direction outside GSE185381's design (A1-G04, prereg-a1-g04-v1) -- %s", format(Sys.Date())), "",
           "Script: scripts/ws/41_a1_g04.R. (a) Chen2023 5 NPM1-mutated AML vs 4 transplant donors, pool-B T cells, J space, P370 (358 on the Chen2023 axis); (b) 18 GSE185381 children vs the 10 healthy donors out of fold.", "",
           "| item | value |", "|---|---|", SUM[, sprintf("| %s | %s |", item, gsub("|", "\\|", value, fixed = TRUE))], "",
           sprintf("**Verdict: A1-G04 (a) = %s%s; (b) = %s.** (b) never counts toward GATE_A3.", cat_a, if (nzchar(qual)) paste0(", ", qual) else "", cat_b),
           sprintf("This supersedes the S9 verdict clause \"no independent AML + healthy cohort exists in this work\" (the S9 block is not edited): Chen2023 was tested here (A1-G04 (a) = %s). %s", cat_a, wording_a),
           "<!-- A1-G04 end -->")
  put(SUM, "summary"); print(SUM)
  gf <- if (SMOKE) file.path(SDIR, "gates_block_smoke.md") else here::here("scripts", "ws", "GATES.md")
  if (SMOKE) writeLines(blk, gf) else write(c("", blk), gf, append = TRUE)
  ts("[combine] done")
}
