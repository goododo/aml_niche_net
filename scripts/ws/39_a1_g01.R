#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds ($raw$T/B/NK/Mono, $fine, $cells; 25 adult S5 samples only);
#          /LARGE1/.../ws_pseudobulk/GSE120221/ws_pb_GSE120221.rds ($fine; the 5.5 loading, healthy donors only);
#          results/tables/ws/ ws_s5_scores.csv, ws_s5_de_T.csv, ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv,
#          ws_a1_g02_fold_genes.csv, ws_a1_g02_oof_scores.csv, ws_s4_pseudobulk_meta_GSE185381.csv,
#          decontx/GSE185381/ws_s1_ambient_blacklist.csv; gencode_GRCh38_gene_order.txt; the Savage 2021 delay set;
#          /LARGE1/.../ws_a1_g01/gdmait_counts.rds and ws_a1_g01_checks_cells.csv (written by scripts/ws/40; --gdmait only)
# OUTPUT : results/tables/ws/ws_a1_g01_{checks,k3,labels,primary,boot,sens,share,regress,axis,gdmait,summary}.csv;
#          /LARGE1/.../ws_a1_g01/ (bootstrap draws and K3 splits, written before use; sens_partial.csv, a checkpoint rewritten
#          after each sensitivity); GATES.md block A1-G01 (--combine)
# what it does: PREREGISTRATION_A1_G01.md (tag prereg-a1-g01-v1), nothing else.
#   --known    known answers of section 6 (deterministic first, then K3); draws written; no readout printed
#   --primary  C_H, C_A, C (3.1) and S_w (3.2) with the 2,000-draw bootstrap; 5.6
#   --sens     sensitivities 4.1-4.7          --share  5.1 (B, NK, monocytes, cross-fitted; K7's durable table)
#   --regress  5.3   --axis  5.5   --gdmait  5.2 (after scripts/ws/40)   --combine  section 7 rules, summary, GATES block
#   --smoke    with any mode: synthetic counts AND synthetic label cell numbers (real sample IDs and status only), few
#              draws; real-value known answers skipped; every output ends in _smoke
#   usage: sbatch scripts/ws/submit_ws_a1_g01.sbatch 39 --known   (the wrapper runs a dated copy of this file)
# This script deletes nothing, does not modify its inputs, never runs or sources scripts 26, 27, 32, 33, 35, 36 or 37
#   (their functions are copied below, with source lines), and opens no Chen2023, paediatric, post-treatment or relapse data.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODES <- c("known", "primary", "sens", "share", "regress", "axis", "gdmait", "combine")
MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of ", paste0("--", MODES, collapse = " "), " [--smoke]")
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
TABWS <- file.path(TAB_DIR, "ws"); GDIR <- file.path(LARGE1_DIR, "ws_a1_g01"); dir.create(GDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])   # path of the script copy being run
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
## write guards: the only places this script writes; a real (non-smoke) file is never overwritten; every table records the md5
wtab <- function(x) { f <- file.path(TABWS, paste0("ws_a1_g01_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); d[, script_md5 := MD5]
  if (length(get0("K3FAIL", envir = globalenv(), ifnotfound = character(0)))) d[, k3_status := "not decisive, descriptive"]   # section 6: a failed K3 bar labels every readout
  fwrite(d, wtab(x)); invisible(d) }
rtab <- function(x) fread(file.path(TABWS, paste0("ws_a1_g01_", x, SFX, ".csv")))
failed_file <- function(x) file.path(TABWS, sprintf("ws_a1_g01_%s_FAILED_%s%s.csv", x, format(Sys.time(), "%Y%m%d%H%M"), SFX))
wg   <- function(x) { f <- file.path(GDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
rg   <- function(x) file.path(GDIR, sub("(\\.[a-z]+)$", paste0(SFX, "\\1"), x))
AFF  <- tryCatch(length(parallel::mcaffinity()), error = function(e) NA_integer_)   # CPUs this process may run on; forked workers inherit it
ts(sprintf("[%s%s] %d worker processes; this process may run on %s CPUs%s", MODE, if (SMOKE) " smoke" else "", NC, AFF,
           if (isTRUE(AFF < NC)) " -- fewer than the workers: is OMP_PROC_BIND=false set? (see the wrapper)" else ""))

## -- 0. constants (section numbers refer to the pre-registration) ------------------------------------------------------
FDR_CUT <- 0.10; LFC_CUT <- 0.5; MIN_CELLS <- WS_S4_MIN_CELLS; FB_CELLS <- 50L; MIN_GENES <- 5L
N_BOOT  <- if (SMOKE) 20L else 2000L; N_BOOT_SENS <- if (SMOKE) 8L else 500L; N_BOOT_K3 <- if (SMOKE) 5L else 200L
N_SPLIT <- if (SMOKE) 3L else 20L;    N_POIS <- if (SMOKE) 2L else 20L;      N_RELAB <- if (SMOKE) 20L else 1000L
OFF <- c(k3_split = 31L, boot = 32L, sens = 33L, pois = 34L, share = 35L, axis = 36L, k3_boot = 37L, k3_thin = 38L, smoke = 99L)  # seeds = SEED + offset
GAP_OOF <- 1.378164                                              # A1-G02 R2 (K0)
CHK <- list()
chk <- function(id, scope, name, expected, got, tol = NULL) {    # deterministic known answer: PASS or FAIL_deterministic
  ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "deterministic", check = name, expected = paste(format(expected), collapse = ";"),
                                         got = paste(format(got), collapse = ";"), state = if (ok) "PASS" else "FAIL_deterministic")
  message(sprintf("  %-4s %-74s expect %-14s got %-14s %s", id, name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
bar <- function(id, name, value, ok) {                           # K3 statistical bar: PASS or FAIL_statistical (final; section 6)
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = "P", kind = "statistical", check = name, expected = "", got = paste(signif(value, 4), collapse = ";"),
                                         state = if (isTRUE(ok)) "PASS" else "FAIL_statistical")
  message(sprintf("  %-4s %-74s got %-14s %s", id, name, paste(signif(value, 4), collapse = ";"), if (isTRUE(ok)) "PASS" else "** FAIL (statistical) **")) }
record <- function(id, name, value) CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = "record", kind = "record", check = name, expected = "", got = value, state = "RECORDED")
gate <- function(scopes) {   # a mode or sensitivity runs only when every deterministic check of its scopes passed (section 6)
  f <- file.path(TABWS, paste0("ws_a1_g01_checks", SFX, ".csv")); if (!file.exists(f)) stop("run --known first: ", f)
  C <- fread(f); stopifnot(nrow(C[scope == "P" & kind == "deterministic"]) > 0)
  all(C[scope %in% scopes & kind == "deterministic", state] == "PASS") }
k3_failed <- function() { C <- fread(file.path(TABWS, paste0("ws_a1_g01_checks", SFX, ".csv"))); C[kind == "statistical" & state != "PASS", paste(id, check)] }

## -- 1. inputs, blacklist, engine (copied from 37_a1_batch.R lines 45-89 and 101-102; inputs passed explicitly) -------
PB  <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
BL  <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO  <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
fit_de <- function(counts, d, terms, coef = "statusAML", group = d$status) {   # 37:52-63 (27's edgeR branch)
  if ("sex" %in% terms && uniqueN(d$sex) < 2) terms <- setdiff(terms, "sex")
  keep_g <- !(rownames(counts) %in% BLACK)
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = group)
  keep <- filterByExpr(y, group = group) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(nrow(X) == nrow(d), qr(X)$rank == ncol(X))
  if (nrow(X) - ncol(X) < 2) stop("residual df < 2: not estimable")
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
  r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, PValue, FDR)]
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]; r }
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")
logz <- function(M) {   # 37:65-68: blacklist; CPM > 1 in >= half of the columns; TMM; log-CPM with prior 1
  M <- M[!rownames(M) %in% BLACK, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1) }
score_col <- function(x, mu, sg, fu, fd) { z <- (x - mu) / sg; z <- z[is.finite(z)]   # 37:101-102 (z-score; K0 and K5 only)
  (if (length(intersect(fu, names(z)))) mean(z[intersect(fu, names(z))]) else 0) - (if (length(intersect(fd, names(z)))) mean(z[intersect(fd, names(z))]) else 0) }

S5ALL <- fread(file.path(TABWS, "ws_s5_scores.csv")); S5 <- S5ALL[compartment == "T"]          # file order = S5's design order
D5  <- S5[, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
dsg <- function(cmp) S5ALL[compartment == cmp, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
MAINFIT <- fread(file.path(TABWS, "ws_s5_de_T.csv")); LFC0 <- setNames(MAINFIT$logFC, MAINFIT$gene)
R_of <- function(r, set, main = MAINFIT, min_cov = 0.8) {   # 37:83-89 (retention of the frozen direction)
  x <- merge(set, r[, .(gene, lfc = logFC)], by = "gene"); x <- merge(x, main[, .(gene, lfc0 = logFC)], by = "gene")
  cov <- nrow(x) / nrow(set)
  if (cov < min_cov) return(list(R = NA_real_, cov = cov, n = nrow(x)))
  list(R = sum(x$s * x$lfc) / sum(abs(x$lfc0)), cov = cov, n = nrow(x)) }
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]
SAMP <- D5$sample; STAT <- as.character(D5$status); AMLS <- SAMP[STAT == "AML"]; HEAL <- SAMP[STAT == "healthy"]
stopifnot(length(SAMP) == 25L, length(AMLS) == 15L, nrow(F402) == 402L, nrow(P370) == 370L)
FG <- fread(file.path(TABWS, "ws_a1_g02_fold_genes.csv"))
FOLD <- lapply(setNames(SAMP, SAMP), function(s) list(U = FG[left_out == s & direction == "up", gene], D = FG[left_out == s & direction == "down", gene]))
OOFT <- fread(file.path(TABWS, "ws_a1_g02_oof_scores.csv"))
SAV  <- "/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv"
MAP  <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE)
TL   <- MAP[group4 == "T", bmm_fine]
SW   <- c("CD4 Naive", "CD8 Effector Memory 2", "CD8 Effector Memory 1", "CD4 Central Memory")   # the four S_w labels (2.2)
NAIVE <- c("CD4 Naive", "CD8 Naive"); EM <- c("CD8 Effector Memory 1", "CD8 Effector Memory 2", "CD4 Effector Memory")
stopifnot(length(TL) == 9L, all(c(SW, NAIVE, EM) %in% TL))

## -- 2. counts: the nine T labels (zero-filled, 2.2), kept-cell numbers, whole-T; synthetic in smoke runs --------------
zf <- function(M, smp) { out <- matrix(0, length(G), length(smp), dimnames = list(G, smp))   # absent label columns -> 0
  h <- intersect(colnames(M), smp); if (length(h)) out[, h] <- M[, h, drop = FALSE]; out }
CNT <- lapply(setNames(TL, TL), function(k) zf(if (is.null(PB$fine[[k]])) matrix(0, length(G), 0, dimnames = list(G, NULL)) else PB$fine[[k]], SAMP))
WT  <- PB$raw$T[, SAMP]
NCELL <- matrix(0, length(TL), length(SAMP), dimnames = list(TL, SAMP))
tmp <- PB$cells[kept == TRUE & sample %in% SAMP & fine %in% TL, .N, by = .(fine, sample)]
NCELL[cbind(tmp$fine, tmp$sample)] <- tmp$N                      # matrix indexing by (row name, column name) pairs
CMP <- lapply(c(B = "B", NK = "NK", Mono = "Mono"), function(cmp) PB$raw[[cmp]][, intersect(colnames(PB$raw[[cmp]]), SAMP), drop = FALSE])
KCELLS <- PB$cells[sample %in% SAMP]; rm(PB); invisible(gc(FALSE))   # only the 25 adult samples are kept from here on (2.1)
if (SMOKE) {   # synthetic label numbers and counts; only the real sample IDs and status are kept (no real per-sample number)
  set.seed(SEED + OFF[["smoke"]])
  sh <- c("CD4 Naive" = 0.20, "CD4 Central Memory" = 0.12, "CD4 Effector Memory" = 0.10, "CD4 Regulatory" = 0.05, "CD8 Naive" = 0.07,
          "CD8 Central Memory" = 0.04, "CD8 Effector Memory 1" = 0.12, "CD8 Effector Memory 2" = 0.25, "CD8 Tissue Resident Memory" = 0.05)[TL]
  prof <- vapply(TL, function(k) { p <- rowSums(CNT[[k]]) + 0.05; p / sum(p) }, numeric(length(G)))   # pooled over all 25 (status-free)
  lf <- setNames(numeric(length(G)), G); g370s <- intersect(P370$gene, G); lf[g370s] <- 0.6 * LFC0[g370s]
  for (s in SAMP) { a <- STAT[SAMP == s] == "AML"; p <- sh * (if (a) ifelse(TL %in% NAIVE, 0.5, 1) else 1)
    NCELL[, s] <- as.vector(rmultinom(1, max(60, round(rlnorm(1, log(400), 0.7))), p / sum(p)))
    for (k in TL) CNT[[k]][, s] <- rnbinom(length(G), size = 20, mu = NCELL[k, s] * 1500 * prof[, k] * (if (a) 2^lf else 1)) }
  WT <- Reduce(`+`, CNT)
  CMP <- lapply(CMP, function(M) { p <- (rowSums(M) + 0.05) / sum(rowSums(M) + 0.05); smp <- colnames(M)
    mu <- outer(p, round(rlnorm(length(smp), log(2e5), 0.5))); am <- smp %in% AMLS; mu[g370s, am] <- mu[g370s, am] * 2^lf[g370s]
    matrix(rnbinom(length(mu), size = 20, mu = mu), nrow(mu), dimnames = dimnames(M)) })
}
G0 <- G[!G %in% BLACK & rowSums(WT) > 0]                          # logz drops these rows anyway; restricting changes no number
CNT0 <- lapply(CNT, function(M) M[G0, , drop = FALSE]); WT0 <- WT[G0, , drop = FALSE]
L25 <- logz(WT); PI <- vapply(TL, function(k) sum(CNT[[k]][rownames(L25), ]) / sum(WT[rownames(L25), ]), 1)   # 2.2 weights (fixed)

## -- 3. the estimators: u-score (2.3), matched lists (2.4), C in both orders (3.1), S_w (3.2) ---------------------------
inst <- function(ids) paste0(ids, "@", seq_along(ids))            # unique column names for repeated samples (DGEList needs them)
halves <- function(X, j, oth, U, D) {   # u = up + dn: X[g, j] minus the mean over the other columns, averaged over U and (negated) over D
  m <- rowMeans(X[c(U, D), oth, drop = FALSE]); x <- X[c(U, D), j] - m
  c(up = mean(x[seq_along(U)]), dn = -mean(x[length(U) + seq_along(D)])) }
mlist <- function(f, unis) { U <- setdiff(f$U, f$D); D <- setdiff(f$D, f$U)   # genes in every universe; NULL below 5 per half
  for (u in unis) { U <- intersect(U, u); D <- intersect(D, u) }
  if (length(U) < MIN_GENES || length(D) < MIN_GENES) NULL else list(U = U, D = D) }
ref_profiles <- function(C, N, w_arm, w_all, labs, mode) {   # 3.1 step 1; w_* = how often each sample enters the pool
  pool <- function(k, w) { n <- sum(N[k, ] * w); if (n > 0) drop(C[[k]] %*% w) / n else numeric(nrow(C[[k]])) }   # %*% sums the weighted columns
  na <- vapply(labs, function(k) sum(N[k, ] * w_arm), 1)
  fb <- if (mode == "arm_only") na == 0 else na < FB_CELLS
  R <- vapply(labs, function(k) if (fb[[k]]) pool(k, w_all) else pool(k, w_arm), numeric(nrow(C[[1]])))
  dimnames(R) <- list(rownames(C[[1]]), labs)                       # set explicitly: an all-zero first label would leave no row names
  list(R = R, fb = fb, none = vapply(labs, function(k) fb[[k]] && sum(N[k, ] * w_all) == 0, TRUE)) }
est_C <- function(ids, st, LST, C = CNT0, N = NCELL, W = WT0, labs = TL, mode = "fallback", umi = FALSE, pois_seed = NA, refs = FALSE) {
  smp <- colnames(W); cn <- inst(ids)
  mult <- function(sel) tabulate(match(ids[sel], smp), length(smp))   # multiplicity of each sample among the selected instances
  XT <- logz(`colnames<-`(W[, ids, drop = FALSE], cn))                # `colnames<-`(M, v) returns M with new column names
  NN <- if (umi) t(vapply(labs, function(k) colSums(C[[k]][, ids, drop = FALSE]), numeric(length(ids)))) else N[labs, ids, drop = FALSE]
  if (umi) rownames(NN) <- labs
  tot <- colSums(W[, ids, drop = FALSE])                             # 4.7: the samples' real whole-T totals on the same genes
  if (!is.na(pois_seed)) set.seed(pois_seed)
  sc <- data.table(inst = cn, sample = ids, status = st, cH = NA_real_, cA = NA_real_, t = NA_real_, nU = 0L, nD = 0L)
  FB <- list(); RF <- list()
  score1 <- function(YH, YA, f, j, oth) { ml <- mlist(f, list(rownames(YH), rownames(YA), rownames(XT)))
    if (is.null(ml) || !length(oth)) return(c(NA, NA, NA, 0, 0))
    c(sum(halves(YH, j, oth, ml$U, ml$D)), sum(halves(YA, j, oth, ml$U, ml$D)), sum(halves(XT, j, oth, ml$U, ml$D)), length(ml$U), length(ml$D)) }
  for (i in unique(ids)) {
    J <- which(ids == i); oth <- which(ids != i); Y <- list()
    for (o in c("H", "A")) {
      rp <- ref_profiles(C, N, mult(st == (if (o == "H") "healthy" else "AML") & ids != i), mult(ids != i), labs, mode)
      R <- rp$R; if (umi) R <- sweep(R, 2, pmax(colSums(R), 1e-300), "/")   # 4.6: profiles scaled to sum 1, weighted by UMI totals
      FB[[length(FB) + 1L]] <- data.table(fold = i, order = o, label = labs, fallback = rp$fb, empty = rp$none)
      if (refs) RF[[paste(i, o)]] <- R
      Y[[o]] <- `colnames<-`(R %*% NN, cn) }
    v <- if (is.na(pois_seed)) score1(logz(Y$H), logz(Y$A), LST[[i]], J[1], oth) else {   # 4.7: Poisson draws, u-scores averaged
      V <- vapply(seq_len(N_POIS), function(p) { yp <- lapply(Y, function(Yo) { m <- sweep(Yo, 2, tot / pmax(colSums(Yo), 1e-300), "*")
        matrix(rpois(length(m), m), nrow(m), dimnames = dimnames(m)) }); score1(logz(yp$H), logz(yp$A), LST[[i]], J[1], oth) }, numeric(5))
      c(rowMeans(V[1:3, , drop = FALSE]), V[4:5, 1]) }
    sc[J, `:=`(cH = v[1], cA = v[2], t = v[3], nU = as.integer(v[4]), nD = as.integer(v[5]))] }
  gp <- function(x) { a <- sc$status == "AML"; mean(x[a], na.rm = TRUE) - mean(x[!a], na.rm = TRUE) }
  gT <- gp(sc$t); gH <- gp(sc$cH); gA <- gp(sc$cA); def <- is.finite(gT) && gT > 0
  CH <- if (def) gH / gT else NA_real_; CA <- if (def) gA / gT else NA_real_
  list(C_H = CH, C_A = CA, C = (CH + CA) / 2, gapT = gT, gapC_H = gH, gapC_A = gA, scores = sc, dropped = unique(sc[is.na(t), sample]),
       fallback = rbindlist(FB), refs = if (refs) RF else NULL) }
est_Sw <- function(ids, st, LST, C = CNT0, N = NCELL, W = WT0, labs = SW, min_arm = 1L) {   # per-label gaps for S_w (3.2)
  rbindlist(lapply(labs, function(k) {
    E <- which(N[k, ids] >= MIN_CELLS); e <- ids[E]; es <- st[E]
    na <- uniqueN(e[es == "AML"]); nh <- uniqueN(e[es == "healthy"])
    out <- data.table(label = k, n_aml = na, n_healthy = nh, kept = FALSE, gap = NA_real_, gapT = NA_real_, gap_up = NA_real_, gapT_up = NA_real_,
                      gap_dn = NA_real_, gapT_dn = NA_real_, n_genes = NA_real_, dropped = "")
    if (na < min_arm || nh < min_arm) return(out)                   # the label leaves all three sums of S_w (3.2)
    cn <- inst(e); Xk <- logz(`colnames<-`(C[[k]][, e, drop = FALSE], cn)); XTk <- logz(`colnames<-`(W[, e, drop = FALSE], cn))
    S <- rbindlist(lapply(unique(e), function(i) { J <- which(e == i); oth <- which(e != i)
      ml <- mlist(LST[[i]], list(rownames(Xk), rownames(XTk)))
      if (is.null(ml) || !length(oth)) return(data.table(j = J, sample = i, su = NA_real_, sd = NA_real_, tu = NA_real_, td = NA_real_, ng = 0L))
      h1 <- halves(Xk, J[1], oth, ml$U, ml$D); h2 <- halves(XTk, J[1], oth, ml$U, ml$D)
      data.table(j = J, sample = i, su = h1[["up"]], sd = h1[["dn"]], tu = h2[["up"]], td = h2[["dn"]], ng = length(ml$U) + length(ml$D)) }))
    S[, status := es[j]]
    gp <- function(x) mean(x[S$status == "AML"], na.rm = TRUE) - mean(x[S$status == "healthy"], na.rm = TRUE)
    out[, `:=`(gap = gp(S$su + S$sd), gapT = gp(S$tu + S$td), gap_up = gp(S$su), gapT_up = gp(S$tu), gap_dn = gp(S$sd), gapT_dn = gp(S$td),
               n_genes = median(S$ng), dropped = paste(unique(S[is.na(su), sample]), collapse = ";"))]
    out[, kept := is.finite(gap) & is.finite(gapT)]; out })) }
sw_of <- function(per, pis = PI, g = "gap", gt = "gapT") {   # S_w = (sum pi) x sum(pi gap) / sum(pi gapT) over the kept labels
  k <- per[kept == TRUE]; if (!nrow(k)) return(NA_real_)
  p <- pis[k$label]; d <- sum(p * k[[gt]]); if (!is.finite(d) || d <= 0) NA_real_ else sum(p) * sum(p * k[[g]]) / d }
sw_paired <- function(pl, no, pis = PI) {   # K3: S_k and S_w from paired gaps (plant minus no plant), labels kept in both
  x <- merge(pl[kept == TRUE, .(label, g1 = gap, t1 = gapT)], no[kept == TRUE, .(label, g0 = gap, t0 = gapT)], by = "label")
  if (!nrow(x)) return(list(S_w = NA_real_, S_k = numeric(0), pi_kept = NA_real_))
  p <- pis[x$label]; d <- sum(p * (x$t1 - x$t0))
  list(S_w = if (is.finite(d) && d > 0) sum(p) * sum(p * (x$g1 - x$g0)) / d else NA_real_, S_k = setNames((x$g1 - x$g0) / (x$t1 - x$t0), x$label), pi_kept = sum(p)) }
ci <- function(v) {   # 95% percentile interval (type 1 = an order statistic); undefined draws enter as -Inf (lower) and +Inf (upper)
  if (is.null(v)) v <- NA_real_
  v[!is.finite(v)] <- NA; nu <- sum(is.na(v)); d <- v[!is.na(v)]
  c(lo = quantile(c(rep(-Inf, nu), d), 0.025, type = 1, names = FALSE), hi = quantile(c(d, rep(Inf, nu)), 0.975, type = 1, names = FALSE), n_undefined = nu) }
draws <- function(B, off, aml, heal) {   # stratified draws (AML block first), one row per draw; off NA = continue the current stream
  if (!is.na(off)) set.seed(SEED + off)
  t(vapply(seq_len(B), function(b) c(sample(aml, length(aml), TRUE), sample(heal, length(heal), TRUE)), character(length(aml) + length(heal)))) }
mc_ok <- function(res, what) {   # a forked worker that dies returns NULL with only a warning, an uncaught error a try-error: stop on both
  bad <- vapply(res, function(r) is.null(r) || inherits(r, "try-error"), TRUE)
  if (any(bad)) stop(what, ": ", sum(bad), " of ", length(res), " results missing (",
                     paste(unique(vapply(res[bad], function(r) if (is.null(r)) "worker died" else trimws(as.character(r)), "")), collapse = " | "), "); rerun")
  res }
boot_run <- function(DR, st, fun) {   # fun(ids, st, b) -> named numeric; undefined draws (NA, 3.1) are counted by ci(); a draw that raises an error stops the mode
  res <- mc_ok(mclapply(seq_len(nrow(DR)), function(b) tryCatch(fun(DR[b, ], st, b), error = function(e) structure(c(error = 1), msg = conditionMessage(e))), mc.cores = NC),
               "bootstrap")
  em <- unique(unlist(lapply(res, attr, "msg")))   # an error is not "cannot be computed" in the sense of 3.1 (gapT <= 0 gives NA, not an error)
  if (length(em)) stop(sum(vapply(res, function(r) !is.null(attr(r, "msg")), TRUE)), " bootstrap draws raised an error (rerun; not counted as undefined): ",
                       paste(head(em, 3), collapse = " | "))
  rbindlist(lapply(seq_along(res), function(b) as.data.table(as.list(c(draw = b, res[[b]])))), fill = TRUE) }
category <- function(CH, CA, C_lo, Sw, Sw_lo, fail = FALSE) {   # section 7, first match wins; Sw has no default, so every caller passes it
  if (fail || !is.finite(CH) || !is.finite(CA)) return(c(category = "not decisive", wording = ""))
  band <- function(x) if (x <= 1/2) 1L else if (x >= 2/3) 3L else 2L
  if (CH >= 2/3 && CA >= 2/3 && isTRUE(C_lo > 1/2)) return(c(category = "composition", wording = ""))
  sw_ok <- length(Sw) == 1L && is.finite(Sw) && isTRUE(Sw_lo > 0.10)   # rule 3: an undefined S_w (estimate or bound) does not meet it
  if (CH <= 1/2 && CA <= 1/2) return(c(category = "beyond T-subset composition",
                                       wording = if (sw_ok) "detected inside the main labels" else "not resolved inside the four main labels"))
  c(category = "mixed", wording = if (band(CH) != band(CA) && abs(CH - CA) >= 0.10) "order-dependent" else if (CH >= 2/3 && CA >= 2/3) "composition not established" else "") }
refit_lists <- function(counts, d, terms = MAIN, amb = TRUE, cov_fun = NULL) {   # 35:97-101 selection, one fold per sample of d
  setNames(mc_ok(mclapply(seq_len(nrow(d)), function(i) { d_tr <- copy(d[-i]); tm <- terms
    if (!is.null(cov_fun)) { d_tr[, cov := cov_fun(i, d_tr$sample)]; tm <- c(terms, "cov") }
    r <- fit_de(counts, d_tr, tm); up <- r[hit == TRUE & logFC > 0, gene]; dn <- r[hit == TRUE & logFC < 0, gene]
    if (amb) list(U = setdiff(up, AMB), D = setdiff(dn, AMB)) else list(U = up, D = dn) }, mc.cores = NC), "refit"), d$sample) }
zgap <- function(excl = character(0)) {   # K0/K5: stored lists, z against the other 24 (post-hoc a1_posthoc_d8_shared_genes.R, oof_from)
  o <- vapply(SAMP, function(s) { tr <- setdiff(SAMP, s); mu <- rowMeans(L25[, tr]); sg <- apply(L25[, tr], 1, sd)
    score_col(L25[, s], mu, sg, setdiff(FOLD[[s]]$U, excl), setdiff(FOLD[[s]]$D, excl)) }, 1)
  mean(o[AMLS]) - mean(o[HEAL]) }
delay_score <- function(L, excl) {   # 37:368-373: Savage up minus down, z over the columns of L, genes in excl left out
  SV <- fread(SAV); Z <- (L - rowMeans(L)) / apply(L, 1, sd)
  u <- setdiff(intersect(SV[direction == "up", gene], rownames(Z)), excl); d <- setdiff(intersect(SV[direction == "down", gene], rownames(Z)), excl)
  (if (length(u)) colMeans(Z[u, , drop = FALSE]) else 0) - (if (length(d)) colMeans(Z[d, , drop = FALSE]) else 0) }

## -- 4. MODE known: deterministic known answers, then K3 (section 6); draws written before use ------------------------
if (MODE == "known") {
  stop_failed <- function(what) { fwrite(rbindlist(CHK), failed_file("checks"))
    ts("[known] ", what, " FAILED: fix the cause, record a dated note in the pre-registration, rerun --known"); quit(status = 1) }
  ts("[known] deterministic known answers")
  if (!SMOKE) {
    chk("K0", "P", "whole-T z gap from the stored fold lists", GAP_OOF, zgap(), tol = 1e-6)
    rl <- refit_lists(WT, D5)
    chk("K0b", "refit", "refit reproduces the 25 stored fold lists exactly (folds)", 25L,
        sum(vapply(SAMP, function(s) setequal(rl[[s]]$U, FOLD[[s]]$U) && setequal(rl[[s]]$D, FOLD[[s]]$D), TRUE)))
    el <- function(k, s) sum(NCELL[k, s] >= MIN_CELLS)
    chk("K1", "P", "eligible AML/healthy: CD4N, CD8 EM2, CD8 EM1, CD4 CM", c(12L, 8L, 12L, 9L, 11L, 7L, 10L, 7L), as.integer(unlist(lapply(SW, function(k) c(el(k, AMLS), el(k, HEAL))))))
    record("K1", "pooled eligible samples per label (record; not split by status)", paste(sprintf("%s %d", TL, vapply(TL, function(k) el(k, SAMP), 1L)), collapse = "; "))
    M4 <- fread(file.path(TABWS, "ws_s4_pseudobulk_meta_GSE185381.csv"))[group4 == "T" & sample %in% c("Control1", "Control3", "Control4")]
    chk("K4", "4.4", "Control1, Control3, Control4 sex_inferred (26:91)", "female;female;female", paste(M4[order(sample), sex_inferred], collapse = ";"))
    chk("K5", "4.2b", "Savage file md5", "5ee444a9284ae52b7afab1c8aeeced2c", unname(tools::md5sum(SAV)))
    chk("K5", "4.2b", "z gap with all 792 Savage genes removed from both halves", 1.3184, zgap(fread(SAV)$gene), tol = 1e-4)
    for (x in list(list("B", 0.692, 174L), list("NK", 0.928, 179L), list("Mono", 0.568, 329L))) { dc <- dsg(x[[1]])[sample %in% SAMP]
      r <- R_of(fit_de(CMP[[x[[1]]]][, dc$sample], dc, MAIN), P370, min_cov = 0)
      chk("K7", "5.1", sprintf("%s: R (D10 formula, min_cov 0) and genes tested", x[[1]]), c(x[[2]], x[[3]]), c(r$R, r$n), tol = c(0.001, 0)) }
    d8 <- copy(D5)[, cov := delay_score(L25, F402$gene)[sample]]
    chk("K8", "4.2a", "MAIN + D8 delay score over the 25 (F402 excluded): R", 0.823, R_of(fit_de(WT, d8, c(MAIN, "cov")), P370)$R, tol = 0.001)
  }
  chk("K2", "P", "nine zero-filled label count vectors sum to whole-T (max |diff|)", 0, max(abs(Reduce(`+`, CNT) - WT)))
  r0 <- est_C(SAMP, STAT, FOLD, refs = TRUE)$refs                  # K2b: est_C's own fold references, then again with sample i's label counts set to zero
  k2b <- unlist(mc_ok(mclapply(SAMP, function(i) { r2 <- est_C(SAMP, STAT, FOLD, C = lapply(CNT0, function(M) { M[, i] <- 0; M }), refs = TRUE)$refs
    same <- function(o) !is.null(r0[[paste(i, o)]]) && identical(r0[[paste(i, o)]], r2[[paste(i, o)]])
    same("H") && same("A") }, mc.cores = NC), "K2b"))
  chk("K2b", "P", "folds whose est_C references ignore the held-out sample's label counts (both orders)", 25L, sum(k2b))
  if (any(rbindlist(CHK)$state == "FAIL_deterministic")) stop_failed("a deterministic known answer")

  ts("[known] K3 on the 10 healthy donors")
  set.seed(SEED + OFF[["k3_split"]]); cmb <- combn(HEAL, 5, simplify = FALSE); SPL <- cmb[sample.int(length(cmb), N_SPLIT)]
  C10 <- lapply(CNT0, function(M) M[, HEAL, drop = FALSE]); N10 <- NCELL[, HEAL, drop = FALSE]
  g370 <- intersect(P370$gene, G0); fct <- 2^LFC0[g370]
  plant <- function(C, donors, labs) { for (k in labs) C[[k]][g370, donors] <- C[[k]][g370, donors] * fct; C }
  set.seed(SEED + OFF[["k3_thin"]])                                # K3b thinning: one rbinom pass per split, in split order
  THIN <- lapply(SPL, function(a) { h <- setdiff(HEAL, a); C <- C10; N <- N10
    for (x in list(list(NAIVE, a), list(EM, h))) for (k in x[[1]]) for (d in x[[2]]) { C[[k]][, d] <- rbinom(nrow(C[[k]]), C[[k]][, d], 0.5); N[k, d] <- floor(N[k, d] / 2) }
    list(C = C, N = N) })
  set.seed(SEED + OFF[["k3_boot"]]); K3DR <- lapply(SPL, function(a) draws(N_BOOT_K3, NA, a, setdiff(HEAL, a)))
  LP370 <- list(U = P370[s > 0, gene], D = P370[s < 0, gene]); ld <- function(C) Reduce(`+`, C)   # whole-T rebuilt as the sum of the labels
  st10 <- rep(c("AML", "healthy"), each = 5)
  K3 <- rbindlist(mc_ok(mclapply(seq_along(SPL), function(si) {
    a <- SPL[[si]]; h <- setdiff(HEAL, a); ids <- c(a, h)
    hp <- function(ks) rowSums(Reduce(`+`, lapply(ks, function(k) C10[[k]][, h, drop = FALSE])))
    cpmv <- function(v) v / sum(v) * 1e6; keepg <- rownames(logz(ld(C10)))   # the whole-T logz universe of the 10 (unplanted)
    ld3 <- log2((cpmv(hp(EM)) + 1) / (cpmv(hp(NAIVE)) + 1))[keepg]
    o3 <- order(-ld3, names(ld3)); L3b <- list(U = names(ld3)[head(o3, 100)], D = names(ld3)[tail(o3, 100)])
    L3e <- list(U = union(L3b$U, LP370$U), D = union(L3b$D, LP370$D)); L3e <- list(U = setdiff(L3e$U, L3e$D), D = setdiff(L3e$D, L3e$U))
    ST <- list(a = list(C = plant(C10, a, TL), N = N10, L = LP370), b = list(C = THIN[[si]]$C, N = THIN[[si]]$N, L = L3b),
               c = list(C = plant(C10, a, "CD4 Naive"), N = N10, L = LP370), d = list(C = plant(C10, a, setdiff(TL, SW)), N = N10, L = LP370),
               e = list(C = plant(THIN[[si]]$C, a, "CD4 Naive"), N = THIN[[si]]$N, L = L3e))
    same <- function(L, ii) setNames(rep(list(L), length(unique(ii))), unique(ii))   # one list for every fold
    run <- function(C, N, LL, ii = ids, refs = FALSE) { x <- est_C(ii, st10, LL, C = C, N = N, W = ld(C), refs = refs)
      x$sw <- est_Sw(ii, st10, LL, C = C, N = N, W = ld(C), min_arm = 2L); x }
    NONE <- list()                                                 # the unplanted split, once per list
    none <- function(key, LL, refs) { k <- paste(key, refs); if (is.null(NONE[[k]])) NONE[[k]] <<- run(C10, N10, LL, refs = refs); NONE[[k]] }
    rows <- list(); add <- function(...) rows[[length(rows) + 1L]] <<- data.table(split = si, ...)
    for (p in names(ST)) { S <- ST[[p]]
      for (lt in c("own", "real")) {
        LL <- if (lt == "own") same(S$L, ids) else FOLD[ids]; rf <- lt == "own" && p %in% c("a", "c", "d")
        n0 <- none(if (lt == "own") switch(p, b = "L3b", e = "L3e", "P370") else "real", LL, rf); n1 <- run(S$C, S$N, LL, refs = rf)
        pr <- function(g1, g0) (g1 - g0) / (n1$gapT - n0$gapT); sp <- sw_paired(n1$sw, n0$sw)
        qk <- if (length(sp$S_k)) paste0("paired_S_k:", names(sp$S_k)) else character(0)
        add(plant = p, list = lt, quantity = c("paired_C_H", "paired_C_A", "paired_S_w", "S_w_bar_rhs", "C_H", "C_A", "C", "S_w", qk),
            value = c(pr(n1$gapC_H, n0$gapC_H), pr(n1$gapC_A, n0$gapC_A), sp$S_w, 0.3 * sp$pi_kept, n1$C_H, n1$C_A, n1$C, sw_of(n1$sw), unname(sp$S_k)))
        w1 <- n1$sw; f1 <- n1$fallback[, .(n = sum(fallback)), by = .(order, fold)]   # reported (6): unpaired S_k, eligibility, labels kept in the paired S_w, fallbacks per order and i
        add(plant = p, list = lt, quantity = c(paste0("S_k:", w1$label), paste0("eligible_AML:", w1$label), paste0("eligible_healthy:", w1$label),
                                               paste0("kept_in_paired_S_w:", w1$label), paste0("n_fallback_", f1$order, ":", f1$fold)),
            value = c(w1$gap / w1$gapT, w1$n_aml, w1$n_healthy, as.numeric(w1$label %in% names(sp$S_k)), f1$n))
        if (rf) {   # deterministic reference tests (section 6)
          pl <- switch(p, a = TL, c = "CD4 Naive", d = setdiff(TL, SW)); Cm <- plant(C10, h, pl)
          mir <- est_C(ids, st10, LL, C = Cm, N = N10, W = ld(Cm), refs = TRUE)   # the same plant on the pseudo-healthy donors instead
          fbH <- n0$fallback[order == "H"]; fbA <- n0$fallback[order == "A"]
          same_ref <- function(x, y, i, o, k) identical(x$refs[[paste(i, o)]][, k, drop = FALSE], y$refs[[paste(i, o)]][, k, drop = FALSE])
          tH <- all(vapply(unique(ids), function(i) same_ref(n0, n1, i, "H", fbH[fold == i & fallback == FALSE, label]), TRUE))
          tA <- all(vapply(unique(ids), function(i) { k <- intersect(fbA[fold == i & fallback == FALSE, label], pl)
            all(vapply(k, function(kk) !identical(n0$refs[[paste(i, "A")]][, kk], n1$refs[[paste(i, "A")]][, kk]), TRUE)) }, TRUE))
          tM <- all(vapply(unique(ids), function(i) same_ref(n0, mir, i, "A", fbA[fold == i & fallback == FALSE, label]), TRUE))
          add(plant = p, list = "test", quantity = c("ref_H_unchanged", "ref_A_changed", "ref_A_mirror_unchanged", "n_fallback_H", "n_fallback_A"),
              value = c(tH, tA, tM, fbH[, sum(fallback)], fbA[, sum(fallback)])) }
        if (lt == "own" && p %in% c("b", "d", "e")) {   # operating characteristics, reported only
          bt <- rbindlist(lapply(seq_len(nrow(K3DR[[si]])), function(b) { ib <- K3DR[[si]][b, ]
            x <- tryCatch(run(S$C, S$N, same(S$L, ib), ii = ib), error = function(e) NULL)
            if (is.null(x)) data.table(C = NA_real_, S_w = NA_real_) else data.table(C = x$C, S_w = sw_of(x$sw)) }))
          cc <- ci(bt$C); cw <- ci(bt$S_w); ct <- category(n1$C_H, n1$C_A, cc[["lo"]], sw_of(n1$sw), cw[["lo"]])
          add(plant = p, list = "own", quantity = c("C_lo", "S_w_lo", paste0("category:", ct[["category"]], if (nzchar(ct[["wording"]])) paste0(" (", ct[["wording"]], ")") else "")),
              value = c(cc[["lo"]], cw[["lo"]], 1)) } } }
    ts(sprintf("    K3 split %d done", si)); rbindlist(rows) }, mc.cores = NC), "K3"))
  stopifnot(setequal(unique(K3$split), seq_along(SPL)))            # every split delivered
  tst <- function(p, q) all(K3[plant == p & list == "test" & quantity == q, value] == 1)
  for (p in c("a", "c", "d")) {
    chk("K3", "P", sprintf("K3%s: C_H references of non-fallback labels unchanged by the plant", p), TRUE, tst(p, "ref_H_unchanged"))
    chk("K3", "P", sprintf("K3%s: C_A references of planted non-fallback labels change", p), TRUE, tst(p, "ref_A_changed"))
    chk("K3", "P", sprintf("K3%s: C_A references unchanged by the same plant on the pseudo-healthy donors", p), TRUE, tst(p, "ref_A_mirror_unchanged")) }
  if (any(rbindlist(CHK)$state == "FAIL_deterministic")) { fwrite(K3[list == "test"], failed_file("k3")); stop_failed("a K3 reference test") }   # test rows only: no bar input is shown
  put(K3, "k3")
  med <- function(p, q) median(K3[plant == p & list == "own" & quantity == q, value], na.rm = TRUE)
  sk <- K3[plant == "a" & list == "own" & startsWith(quantity, "paired_S_k:"), .(m = median(value, na.rm = TRUE)), by = quantity]
  ma <- median(abs(K3[plant == "a" & list == "own" & quantity == "paired_C_A", value]), na.rm = TRUE)   # the median of |paired C_A| (review 4's bar)
  bar("K3a", "median |paired C_A| <= 0.10", ma, ma <= 0.10)
  bar("K3a", "median paired S_k of every kept label in [0.3, 1.3]", sk$m, nrow(sk) > 0 && all(sk$m >= 0.3 & sk$m <= 1.3))
  bar("K3a", "median paired S_w >= 0.3 x (sum of pi over kept labels)", c(med("a", "paired_S_w"), med("a", "S_w_bar_rhs")), med("a", "paired_S_w") >= med("a", "S_w_bar_rhs"))
  bar("K3b", "median paired C_H >= 0.8 and paired C_A >= 0.8", c(med("b", "paired_C_H"), med("b", "paired_C_A")), med("b", "paired_C_H") >= 0.8 && med("b", "paired_C_A") >= 0.8)
  bar("K3c", "median paired S_w >= 0.2", med("c", "paired_S_w"), med("c", "paired_S_w") >= 0.2)
  record("K3d", "share of splits whose S_w lower bound exceeds 0.10 (own list; reported only)",
         sprintf("%.3f", mean(K3[plant == "d" & list == "own" & quantity == "S_w_lo", value] > 0.10)))
  ts("[known] splits and bootstrap draws (written before use)")
  fwrite(data.table(split = seq_along(SPL), pseudo_aml = vapply(SPL, paste, "", collapse = ";")), wg("k3_splits.csv"))
  fwrite(as.data.table(draws(N_BOOT, OFF[["boot"]], AMLS, HEAL)), wg("boot_draws.csv"))
  CK <- put(rbindlist(CHK), "checks")
  ts(sprintf("[known] done: %d deterministic PASS; statistical bars: %s", CK[kind == "deterministic" & state == "PASS", .N],
             paste(CK[kind == "statistical", sprintf("%s %s", id, state)], collapse = ", ")))
  quit(status = 0)
}

## -- 5. shared by the readout modes --------------------------------------------------------------------------------------
if (!gate("P")) stop("a scope-P known answer has not passed")
K3FAIL <- k3_failed(); if (length(K3FAIL)) ts("NOTE: K3 statistical bar(s) failed; outputs are 'not decisive, descriptive': ", paste(K3FAIL, collapse = "; "))
STB <- c(rep("AML", length(AMLS)), rep("healthy", length(HEAL)))   # status of a draw's columns (AML block first)
fb_share <- function(x, o) { f <- x$fallback[order == o, .(fb = any(fallback)), by = label]; setNames(as.numeric(f$fb), paste0("fb_", o, ":", f$label)) }
summ <- function(x, w, B, tag) {   # long rows: estimates with CIs over the bootstrap table B
  q <- c("C_H", "C_A", "C", "S_w"); cs <- lapply(q, function(z) ci(B[[z]]))
  rbind(data.table(section = tag, quantity = c(q, "gapT", "gapC_H", "gapC_A"), label = "", sample = "", value = c(x$C_H, x$C_A, x$C, sw_of(w), x$gapT, x$gapC_H, x$gapC_A),
                   lo = c(vapply(cs, `[[`, 1, "lo"), NA, NA, NA), hi = c(vapply(cs, `[[`, 1, "hi"), NA, NA, NA),
                   n_undefined = c(vapply(cs, `[[`, 1, "n_undefined"), NA, NA, NA), note = ""),
        data.table(section = tag, quantity = "dropped_samples", label = "", sample = paste(x$dropped, collapse = ";"), value = length(x$dropped), lo = NA, hi = NA,
                   n_undefined = NA, note = "5-gene rule")) }

## -- 6. MODE primary: 3.1, 3.2, 5.6 ---------------------------------------------------------------------------------------
if (MODE == "primary") {
  ts("[primary] estimates")
  X0 <- est_C(SAMP, STAT, FOLD); XA <- est_C(SAMP, STAT, FOLD, mode = "arm_only"); W0 <- est_Sw(SAMP, STAT, FOLD)
  ts("[primary] bootstrap, ", N_BOOT, " draws on ", NC, " cores")
  DR <- as.matrix(fread(rg("boot_draws.csv"))); stopifnot(nrow(DR) == N_BOOT, ncol(DR) == 25L)
  BT <- boot_run(DR, STB, function(ids, st, b) { x <- est_C(ids, st, FOLD); w <- est_Sw(ids, st, FOLD)
    c(C_H = x$C_H, C_A = x$C_A, C = x$C, S_w = sw_of(w), gapT = x$gapT, fb_share(x, "H"), fb_share(x, "A")) })
  put(BT, "boot")
  fbt <- function(x) x$fallback[, .(folds = sum(fallback), empty = sum(empty), fold_ids = paste(fold[fallback], collapse = ";"),   # 3.1: fallback labels per fold
                                    empty_ids = paste(fold[empty], collapse = ";")), by = .(order, label)]
  FBS <- fbt(X0); FBA <- fbt(XA)
  ids_note <- function(n, ids) ifelse(n > 0, paste0(" (", ids, ")"), "")
  FBS[, draw_share := vapply(seq_len(.N), function(r) { cl <- paste0("fb_", order[r], ":", label[r]); if (cl %in% names(BT)) mean(BT[[cl]], na.rm = TRUE) else NA_real_ }, 1)]
  SC <- X0$scores
  OUT <- rbind(summ(X0, W0, BT, "3.1-3.2"),
               data.table(section = "3.1 arm-only", quantity = c("C_H", "C_A", "C"), label = "", sample = "", value = c(XA$C_H, XA$C_A, XA$C), lo = NA, hi = NA, n_undefined = NA,
                          note = "arm-only profiles (3.1)"),
               data.table(section = "3.1 scores", quantity = rep(c("c_H", "c_A", "t"), each = nrow(SC)), label = "", sample = rep(SC$sample, 3), value = c(SC$cH, SC$cA, SC$t),
                          lo = NA, hi = NA, n_undefined = NA, note = rep(sprintf("%s; genes %d up / %d down", SC$status, SC$nU, SC$nD), 3)),
               data.table(section = "3.1 fallback", quantity = paste0("folds_with_fallback_", FBS$order), label = FBS$label, sample = FBS$fold_ids, value = FBS$folds, lo = NA, hi = NA,
                          n_undefined = NA, note = sprintf("empty pools %d%s; share of draws with a fallback %.3f", FBS$empty, ids_note(FBS$empty, FBS$empty_ids), FBS$draw_share)),
               data.table(section = "3.1 arm-only fallback", quantity = paste0("folds_with_fallback_", FBA$order), label = FBA$label, sample = FBA$fold_ids, value = FBA$folds,
                          lo = NA, hi = NA, n_undefined = NA, note = sprintf("empty pools %d%s", FBA$empty, ids_note(FBA$empty, FBA$empty_ids))),
               W0[, .(section = "3.2 labels", quantity = "S_k", label, sample = "", value = gap / gapT, lo = NA, hi = NA, n_undefined = NA,
                      note = sprintf("gap %.4f gapT %.4f%s; up S_k %.3f, down S_k %.3f; eligible %d AML / %d healthy; median genes %s; dropped %s", gap, gapT,
                                     ifelse(is.finite(gapT) & gapT < X0$gapT / 4, " (gapT_k < gapT/4)", ""), gap_up / gapT_up, gap_dn / gapT_dn, n_aml, n_healthy, n_genes, dropped))],
               data.table(section = "3.2 halves", quantity = c("S_w_up", "S_w_down"), label = "", sample = "", value = c(sw_of(W0, g = "gap_up", gt = "gapT_up"), sw_of(W0, g = "gap_dn", gt = "gapT_dn")),
                          lo = NA, hi = NA, n_undefined = NA, note = ""))
  put(OUT, "primary")
  ## labels table: pi, UMI per kept cell by arm (4.6), share of kept T in the four S_w labels and kept fractions (5.6)
  upc <- function(k, s) sum(CNT[[k]][, s]) / max(sum(NCELL[k, s]), 1)
  swsh <- colSums(NCELL[SW, , drop = FALSE]) / colSums(NCELL)
  LAB <- data.table(label = TL, pi = PI[TL], sw_label = TL %in% SW, umi_per_cell_AML = vapply(TL, upc, 1, s = AMLS), umi_per_cell_healthy = vapply(TL, upc, 1, s = HEAL),
                    cells_AML = rowSums(NCELL[TL, AMLS]), cells_healthy = rowSums(NCELL[TL, HEAL]),
                    sw_share_median_AML = median(swsh[AMLS]), sw_share_median_healthy = median(swsh[HEAL]))
  if (!SMOKE) {   # kept fractions use the real per-cell table (smoke runs have synthetic label numbers and no per-cell table)
    PC <- KCELLS[fine %in% c(TL, "T Proliferating")][, status := STAT[match(sample, SAMP)]]
    kf <- dcast(PC[fine %in% TL, .(kf = mean(kept)), by = .(label = fine, status)], label ~ status, value.var = "kf")
    LAB <- merge(LAB, kf[, .(label, kept_fraction_AML = AML, kept_fraction_healthy = healthy)], by = "label", all.x = TRUE)
    tp <- PC[, .(x = mean(fine == "T Proliferating")), by = status]
    LAB[, `:=`(T_prolif_share_AML = tp[status == "AML", x], T_prolif_share_healthy = tp[status == "healthy", x])]
  }
  put(LAB, "labels")
  ts(sprintf("[primary] written (C_H %.3f, C_A %.3f, C %.3f, S_w %.3f; the category is applied by --combine)", X0$C_H, X0$C_A, X0$C, sw_of(W0)))
  quit(status = 0)
}

## -- 7. MODE sens: 4.1-4.7 ---------------------------------------------------------------------------------------------------
if (MODE == "sens") {
  P <- rtab("primary"); pC <- P[section == "3.1-3.2" & quantity == "C", value]
  bd <- if (length(pC) != 1L || !is.finite(pC)) c(NA_real_, NA_real_) else if (pC <= 1/2) c(-Inf, 1/2) else if (pC >= 2/3) c(2/3, Inf) else c(1/2, 2/3)
  DRS <- draws(N_BOOT_SENS, OFF[["sens"]], AMLS, HEAL)
  run_sens <- function(tag, ids, st, LST, DR, stb, min_arm = 1L, ...) { extra <- list(...)   # extra arguments go to est_C (umi, pois)
    x <- do.call(est_C, c(list(ids, st, LST), extra)); w <- est_Sw(ids, st, LST, min_arm = min_arm)
    B <- boot_run(DR, stb, function(ii, ss, b) { e2 <- extra; if (!is.null(e2$pois_seed)) e2$pois_seed <- e2$pois_seed + b
      y <- do.call(est_C, c(list(ii, ss, LST), e2)); c(C_H = y$C_H, C_A = y$C_A, C = y$C, S_w = sw_of(est_Sw(ii, ss, LST, min_arm = min_arm))) })
    cc <- ci(B$C)
    rbind(summ(x, w, B, tag),
          data.table(section = tag, quantity = "genes_median", label = "", sample = "", value = median(x$scores[is.finite(t), nU + nD]), lo = NA, hi = NA, n_undefined = NA,
                     note = "matched genes (up + down) per scored sample, C (5.6)"),
          w[, .(section = tag, quantity = "S_k", label, sample = "", value = gap / gapT, lo = NA, hi = NA, n_undefined = NA,   # 4.4 and 5.6: eligible counts per label
                note = sprintf("eligible %d AML / %d healthy; in S_w %s; median genes %s; dropped %s", n_aml, n_healthy, kept, n_genes, dropped))],
          data.table(section = tag, quantity = "band_reading", label = "", sample = "", value = NA, lo = NA, hi = NA, n_undefined = NA,
                     note = if (anyNA(bd)) "primary C undefined: no band" else if (cc[["lo"]] > bd[2] || cc[["hi"]] < bd[1]) "moves C out of the primary band" else "consistent within uncertainty")) }
  RES <- list(); ckp <- function() { d <- rbindlist(RES, fill = TRUE)[, script_md5 := MD5]   # checkpoint after each sensitivity (not a declared output)
    if (length(K3FAIL)) d[, k3_status := "not decisive, descriptive"]; fwrite(d, rg("sens_partial.csv")) }
  ok_ref <- gate(c("P", "refit")); if (!ok_ref) ts("[sens] 4.1, 4.2a, 4.3 and 4.4 blocked: K0b failed")
  if (ok_ref) { ts("[sens] 4.1 refit without the ambient filter (402-analogue)"); L41 <- refit_lists(WT, D5, amb = FALSE)
    RES$s41 <- run_sens("4.1", SAMP, STAT, L41, DRS, STB); ckp() }
  if (ok_ref && gate("4.2a")) { ts("[sens] 4.2a refit with the delay score in MAIN")
    cf <- function(i, smp) delay_score(logz(WT[, SAMP[-i]]), union(L41[[i]]$U, L41[[i]]$D))[smp]   # over the 24, the fold's own 402-analogue hits excluded
    RES$s42a <- run_sens("4.2a", SAMP, STAT, refit_lists(WT, D5, cov_fun = cf), DRS, STB); ckp() }
  if (gate("4.2b")) { ts("[sens] 4.2b stored lists minus the Savage genes"); sv <- fread(SAV)$gene
    RES$s42b <- run_sens("4.2b", SAMP, STAT, lapply(FOLD, function(f) list(U = setdiff(f$U, sv), D = setdiff(f$D, sv))), DRS, STB); ckp() }
  if (ok_ref) { ts("[sens] 4.3 Control4 coded female"); Df <- copy(D5)[sample == "Control4", sex := factor("female", levels(D5$sex))]
    RES$s43 <- run_sens("4.3", SAMP, STAT, refit_lists(WT, Df), DRS, STB); ckp() }
  if (ok_ref && gate("4.4")) { ts("[sens] 4.4 without Control1, Control3, Control4")
    keep <- !SAMP %in% c("Control1", "Control3", "Control4"); s4 <- SAMP[keep]; st4 <- STAT[keep]
    L44 <- refit_lists(WT, D5[keep]); x44 <- est_C(s4, st4, L44); inC <- x44$scores[!is.na(t)]
    nfail <- length(x44$dropped); nh <- inC[status == "healthy", uniqueN(sample)]; na <- inC[status == "AML", uniqueN(sample)]
    el4 <- function(arm) rowSums(NCELL[SW, s4[st4 == arm], drop = FALSE] >= MIN_CELLS)
    RES$s44 <- if (nfail <= 22 / 3 && nh >= 5 && na >= 10) {
      stb4 <- c(rep("AML", sum(st4 == "AML")), rep("healthy", sum(st4 == "healthy")))
      run_sens("4.4", s4, st4, L44, draws(N_BOOT_SENS, OFF[["sens"]], s4[st4 == "AML"], s4[st4 == "healthy"]), stb4, min_arm = 3L)
    } else rbind(data.table(section = "4.4", quantity = "not estimable", label = "", sample = paste(x44$dropped, collapse = ";"), value = NA, lo = NA, hi = NA, n_undefined = NA,
                            note = sprintf("folds failing the 5-gene rule %d; healthy in C %d; AML in C %d", nfail, nh, na)),
                 data.table(section = "4.4", quantity = "eligible", label = SW, sample = "", value = NA, lo = NA, hi = NA, n_undefined = NA,
                            note = sprintf("eligible %d AML / %d healthy", el4("AML"), el4("healthy"))))
    ckp() }
  ts("[sens] 4.5 leave one S_w label out")
  RES$s45 <- rbindlist(lapply(SW, function(k) { lb <- setdiff(SW, k)
    B <- boot_run(DRS, STB, function(ii, ss, b) c(S_w = sw_of(est_Sw(ii, ss, FOLD, labs = lb)))); cw <- ci(B$S_w)
    data.table(section = "4.5", quantity = "S_w", label = paste("without", k), sample = "", value = sw_of(est_Sw(SAMP, STAT, FOLD, labs = lb)),
               lo = cw[["lo"]], hi = cw[["hi"]], n_undefined = cw[["n_undefined"]], note = "") })); ckp()
  ts("[sens] 4.6 UMI-weighted C"); RES$s46 <- run_sens("4.6", SAMP, STAT, FOLD, DRS, STB, umi = TRUE); ckp()
  ts("[sens] 4.7 count-noise-matched C"); RES$s47 <- run_sens("4.7", SAMP, STAT, FOLD, DRS, STB, pois_seed = SEED + OFF[["pois"]])
  put(rbindlist(RES, fill = TRUE), "sens"); ts("[sens] done"); quit(status = 0)
}

## -- 8. MODE share: 5.1 (and K7's durable table) ------------------------------------------------------------------------------
if (MODE == "share") {
  if (!gate("5.1")) stop("K7 failed: 5.1 blocked")
  ho <- function(ids, M) {   # per donor: scored with its own fold's T genes, in the compartment and in whole-T (same donors and genes); NA below 5 genes per half
    cn <- inst(ids); XC <- logz(`colnames<-`(M[, ids, drop = FALSE], cn)); XT <- logz(`colnames<-`(WT[, ids, drop = FALSE], cn))
    rbindlist(lapply(unique(ids), function(i) { J <- which(ids == i); oth <- which(ids != i); ml <- mlist(FOLD[[i]], list(rownames(XC), rownames(XT)))
      if (is.null(ml) || !length(oth)) return(data.table(j = J, sample = i, s = NA_real_, t = NA_real_, cov = NA_real_))
      data.table(j = J, sample = i, s = sum(halves(XC, J[1], oth, ml$U, ml$D)), t = sum(halves(XT, J[1], oth, ml$U, ml$D)),
                 cov = (length(ml$U) + length(ml$D)) / (length(FOLD[[i]]$U) + length(FOLD[[i]]$D))) })) }
  heldout_gap <- function(S, st) { a <- st[S$j] == "AML"   # coverage = mean over the scored donors
    gs <- mean(S$s[a], na.rm = TRUE) - mean(S$s[!a], na.rm = TRUE); gt <- mean(S$t[a], na.rm = TRUE) - mean(S$t[!a], na.rm = TRUE)
    c(gap_cmp = gs, gap_T = gt, ratio = if (is.finite(gt) && gt > 0) gs / gt else NA_real_, coverage = mean(S$cov, na.rm = TRUE)) }
  SH <- rbindlist(lapply(c("B", "NK", "Mono"), function(cmp) { dc <- dsg(cmp)[sample %in% SAMP]; ids <- dc$sample; st <- as.character(dc$status)
    S0 <- ho(ids, CMP[[cmp]]); x <- heldout_gap(S0, st); a <- ids[st == "AML"]; h <- ids[st == "healthy"]; sc <- is.finite(S0$s); s0 <- st[S0$j]
    B <- boot_run(draws(N_BOOT, OFF[["share"]], a, h), c(rep("AML", length(a)), rep("healthy", length(h))), function(ii, ss, b) heldout_gap(ho(ii, CMP[[cmp]]), ss))
    cc <- ci(B$ratio); r7 <- if (SMOKE) list(R = NA_real_, n = NA_integer_) else R_of(fit_de(CMP[[cmp]][, ids], dc, MAIN), P370, min_cov = 0)
    data.table(compartment = cmp, n_AML = length(a), n_healthy = length(h), n_AML_scored = uniqueN(S0$sample[sc & s0 == "AML"]),
               n_healthy_scored = uniqueN(S0$sample[sc & s0 == "healthy"]), dropped = if (any(!sc)) paste(unique(S0$sample[!sc]), collapse = ";") else "none",
               coverage = x[["coverage"]], gap_cmp = x[["gap_cmp"]], gap_T = x[["gap_T"]],
               ratio = x[["ratio"]], ratio_lo = cc[["lo"]], ratio_hi = cc[["hi"]], n_undefined = cc[["n_undefined"]],
               reading = if (cc[["lo"]] > 0) "direction shown out of sample" else "not shown (a ratio near 0 is not evidence of T-specificity)",
               caveat = paste0("the compartment's own subset composition is not controlled", if (cmp == "Mono")   # 5.1: stated with the result
                 "; monocytes in AML samples are often clonal, so a shift may be leukaemic biology rather than a shared environment or handling effect" else ""),
               insample_R_K7 = r7$R, insample_genes_K7 = r7$n) }))
  put(SH, "share"); ts("[share] done"); quit(status = 0)
}

## -- 9. MODE regress: 5.3 -------------------------------------------------------------------------------------------------------
if (MODE == "regress") {
  tot <- colSums(NCELL)
  d <- data.table(sample = SAMP, status = factor(STAT, c("healthy", "AML")), naive = colSums(NCELL[NAIVE, ]) / tot,
                  cd8em = colSums(NCELL[c("CD8 Effector Memory 1", "CD8 Effector Memory 2"), ]) / tot)
  d[, oof := OOFT$oof_score[match(d$sample, OOFT$sample)]]
  lo <- max(d[, min(naive), by = status]$V1); hi <- min(d[, max(naive), by = status]$V1)
  n_in <- d[naive >= lo & naive <= hi, .N, by = status]; ok <- lo <= hi && nrow(n_in) == 2L && all(n_in$N >= 3)
  R <- if (ok) { f <- summary(lm(oof ~ status + naive + cd8em, data = d))$coefficients
    data.table(term = rownames(f), estimate = f[, 1], se = f[, 2], p = f[, 4], overlap = sprintf("%.3f-%.3f", lo, hi)) } else
    data.table(term = "not estimable: no overlap", estimate = NA_real_, se = NA_real_, p = NA_real_, overlap = if (lo <= hi) sprintf("%.3f-%.3f", lo, hi) else "none")
  put(R, "regress"); ts("[regress] done"); quit(status = 0)
}

## -- 10. MODE axis: 5.5 --------------------------------------------------------------------------------------------------------
if (MODE == "axis") {
  PBo <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE120221", "ws_pb_GSE120221.rds"))   # healthy donors only
  pool <- function(ks) Reduce(`+`, lapply(ks, function(k) rowSums(PBo$fine[[k]])))
  Mo <- cbind(naive = pool(NAIVE), em = pool(EM)); Lo <- logz(Mo)                               # the loading's logz universe
  cp <- cpm(calcNormFactors(DGEList(Mo[rownames(Lo), ])))
  LOAD <- log2((cp[, "em"] + 1) / (cp[, "naive"] + 1)); EG <- log2(cpm(DGEList(matrix(rowSums(Mo[rownames(Lo), ]), dimnames = list(rownames(Lo), "x"))))[, 1] + 1)
  uni <- unique(unlist(FOLD))
  AX <- rbindlist(lapply(c("CD4 Naive", "CD8 Effector Memory 2"), function(k) {
    E <- SAMP[NCELL[k, ] >= MIN_CELLS]; sE <- STAT[match(E, SAMP)]; Xk <- logz(CNT[[k]][, E]); g <- intersect(rownames(Xk), names(LOAD))
    dec <- cut(rank(EG[g], ties.method = "first"), 10, labels = FALSE)   # deciles of the status-free Oetjen expression
    stat <- function(lab) { dg <- rowMeans(Xk[g, lab == "AML", drop = FALSE]) - rowMeans(Xk[g, lab == "healthy", drop = FALSE])
      mean(vapply(1:10, function(q) suppressWarnings(cor(dg[dec == q], LOAD[g][dec == q], method = "spearman")), 1), na.rm = TRUE) }
    obs <- stat(sE); set.seed(SEED + OFF[["axis"]]); nul <- vapply(seq_len(N_RELAB), function(r) stat(sample(sE)), 1)   # sample(x) permutes x (arm sizes kept)
    gu <- intersect(uni, g); dgu <- rowMeans(Xk[gu, sE == "AML", drop = FALSE]) - rowMeans(Xk[gu, sE == "healthy", drop = FALSE])
    p <- (1 + sum(nul >= obs)) / (N_RELAB + 1)
    data.table(label = k, n_AML = sum(sE == "AML"), n_healthy = sum(sE == "healthy"), n_genes = length(g), statistic = obs, p_relabel = p,
               reading = if (obs > 0 && p < 0.05) "the within-label shift partly follows the naive-to-effector axis, so composition below the label level may contribute" else "",
               union_rho_number_only = cor(dgu, LOAD[gu], method = "spearman"), n_union = length(gu)) }))
  put(AX, "axis"); ts("[axis] done"); quit(status = 0)
}

## -- 11. MODE gdmait: 5.2 (counts written by scripts/ws/40) ---------------------------------------------------------------------
if (MODE == "gdmait") {
  cf <- file.path(TABWS, paste0("ws_a1_g01_checks_cells", SFX, ".csv")); if (!file.exists(cf) || !all(fread(cf)$state == "PASS")) stop("run scripts/ws/40 first (its checks must pass)")
  GM <- readRDS(rg("gdmait_counts.rds"))                                               # $fine (11 labels), $ncell (11 x 25), $wt_noflag, $cells
  L11 <- names(GM$fine); C11 <- lapply(GM$fine, function(M) M[G0, SAMP, drop = FALSE]); N11 <- GM$ncell[L11, SAMP]
  x11 <- est_C(SAMP, STAT, FOLD, C = C11, N = N11, W = Reduce(`+`, C11), labs = L11)
  tgap <- function(X, LST) { o <- vapply(SAMP, function(i) { ml <- mlist(LST[[i]], list(rownames(X)))
      if (is.null(ml)) NA_real_ else sum(halves(X, which(SAMP == i), which(SAMP != i), ml$U, ml$D)) }, 1)
    list(gap = mean(o[AMLS], na.rm = TRUE) - mean(o[HEAL], na.rm = TRUE), dropped = SAMP[is.na(o)]) }
  dr <- function(x) if (length(x)) paste(x, collapse = ";") else "none"
  noTRDC <- lapply(FOLD, function(f) list(U = setdiff(f$U, "TRDC"), D = setdiff(f$D, "TRDC")))
  CE <- copy(GM$cells)[, status := STAT[match(sample, SAMP)]]; sane <- CE[status == "healthy", median(gd / kept_T)]
  LW <- logz(WT0); TG <- list(tgap(LW, FOLD), tgap(logz(GM$wt_noflag[G0, SAMP]), FOLD), tgap(LW, noTRDC))
  GD <- rbind(data.table(quantity = c("C_H_11", "C_A_11", "C_11"), value = c(x11$C_H, x11$C_A, x11$C),
                         note = paste0("gamma-delta and MAIT as two extra labels (flags are lower bounds); dropped (5-gene rule): ", dr(x11$dropped))),
              data.table(quantity = c("u_gap_wholeT", "u_gap_wholeT_noflag", "u_gap_wholeT_noTRDC"), value = vapply(TG, `[[`, 1, "gap"),
                         note = paste0("fold genes fixed; dropped (5-gene rule): ", vapply(TG, function(z) dr(z$dropped), ""))),
              CE[, .(quantity = c("gd_share_median", "mait_share_median"), value = c(median(gd / kept_T), median(mait / kept_T)), note = status[1]), by = status][, -"status"],
              CE[, .(quantity = paste0("gd_share:", sample), value = gd / kept_T, note = status)],      # 5.2: flagged shares per sample, by status
              CE[, .(quantity = paste0("mait_share:", sample), value = mait / kept_T, note = status)],
              data.table(quantity = "genes_median_C_11", value = median(x11$scores[is.finite(t), nU + nD]), note = "matched genes (up + down) per scored sample (5.6)"),
              data.table(quantity = paste0("cells:", L11), value = rowSums(N11), note = sprintf("AML %d / healthy %d", rowSums(N11[, AMLS]), rowSums(N11[, HEAL]))),
              data.table(quantity = "sanity: healthy median gamma-delta share in 0.5-15%", value = sane, note = if (sane >= 0.005 && sane <= 0.15) "within range" else "OUTSIDE range"))
  put(GD, "gdmait"); ts("[gdmait] done"); quit(status = 0)
}

## -- 12. MODE combine: section 7, summary, GATES.md block ------------------------------------------------------------------------
if (MODE == "combine") {
  P <- rtab("primary"); v <- function(q) P[section == "3.1-3.2" & quantity == q]
  fail <- length(K3FAIL) > 0
  ct <- category(v("C_H")$value, v("C_A")$value, v("C")$lo, v("S_w")$value, v("S_w")$lo, fail = fail)
  q4 <- c("C_H", "C_A", "C", "S_w"); nn <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
  ci_txt <- function(q) sprintf("%.3f to %.3f (undefined draws %d)", v(q)$lo, v(q)$hi, as.integer(v(q)$n_undefined))
  md5_of <- function(f) if (!file.exists(f)) "absent" else if (!"script_md5" %in% names(fread(f, nrows = 0))) "not recorded" else
    paste(unique(fread(f, select = "script_md5")$script_md5), collapse = ";")
  ver <- rbindlist(lapply(c("checks", "k3", "primary", "boot", "labels", "sens", "share", "regress", "axis", "gdmait", "checks_cells"), function(x)   # code version of every table
    data.table(table = x, md5 = md5_of(file.path(TABWS, paste0("ws_a1_g01_", x, SFX, ".csv"))))))
  m39 <- ver[table == "checks", md5]; other <- ver[!table %in% c("checks", "checks_cells") & md5 != "absent" & md5 != m39, table]
  if (!identical(MD5, m39)) other <- c(other, "summary (this --combine run)")   # the code that applies section 7 is checked too
  SUM <- data.table(item = c("category", "wording", as.vector(rbind(q4, paste0(q4, "_95CI"))), "k3_statistical_failures",
                             "script_md5 of 39 --known", "tables written by another version of 39", "tables absent", "script_md5 of 40"),
                    value = c(ct[["category"]], ct[["wording"]], as.vector(rbind(sprintf("%.3f", vapply(q4, function(q) v(q)$value, 1)), vapply(q4, ci_txt, ""))),
                              if (fail) paste(K3FAIL, collapse = "; ") else "none", m39, nn(other), nn(ver[md5 == "absent", table]), ver[table == "checks_cells", md5]))
  if (ct[["category"]] == "composition") {   # section 7: the wording under "composition" also gives S_w's sign and CI, S9's f_age and the age confounding
    sg <- if (isTRUE(v("S_w")$value > 0)) "positive" else if (isTRUE(v("S_w")$value < 0)) "negative" else "undefined"
    SUM <- rbind(SUM, data.table(item = c("S_w sign", "f_age (S9)", "age"), value = c(sprintf("%s (95%% CI %s)", sg, ci_txt("S_w")),
                                 "0.102 (95% CI -0.085 to 0.463)", "in GSE185381 the subtype mix is confounded with age"))) }
  put(SUM, "summary"); print(SUM)
  blk <- c("<!-- A1-G01 start -->", sprintf("## T imprint: composition versus not (A1-G01, prereg-a1-g01-v1) -- %s", format(Sys.Date())), "",
           "Scripts: scripts/ws/39_a1_g01.R, scripts/ws/40_a1_g01_cells.R. 15 adult AML vs 10 healthy (GSE185381), A1-G02's stored fold lists.", "",
           "| item | value |", "|---|---|", SUM[, sprintf("| %s | %s |", item, gsub("|", "\\|", value, fixed = TRUE))], "",   # a bar name such as "median |paired C_A|" would split the cell
           sprintf("**Verdict: A1-G01 = %s%s.**", ct[["category"]], if (nzchar(ct[["wording"]])) paste0(" (", ct[["wording"]], ")") else ""), "<!-- A1-G01 end -->")
  gf <- if (SMOKE) file.path(GDIR, "gates_block_smoke.md") else here::here("scripts", "ws", "GATES.md")
  if (!SMOKE && any(grepl("A1-G01 start", readLines(gf)))) stop("GATES.md already holds an A1-G01 block")
  if (SMOKE) writeLines(blk, gf) else write(c("", blk), gf, append = TRUE)
  ts("[combine] done"); quit(status = 0)
}
