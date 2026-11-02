#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds of the 25 adult S5 T samples (raw counts);
#          /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds ($cells: kept T cells and their fine labels; $raw$T, $fine);
#          results/tables/ws/ws_s5_scores.csv; results/tables/ws/ws_a1_g01_checks.csv (scope-P known answers must have passed)
# OUTPUT : /LARGE1/.../ws_a1_g01/gdmait_counts.rds: $fine = 11 label count matrices (genes x 25: the nine T labels without
#            flagged cells, "gamma-delta", "MAIT"), $ncell (11 x 25), $wt_noflag (whole-T without flagged cells), $cells
#            (per sample: kept T cells, gamma-delta and MAIT flagged cells); results/tables/ws/ws_a1_g01_checks_cells.csv
# what it does: PREREGISTRATION_A1_G01.md section 5.2, the per-cell step only (tag prereg-a1-g01-v1). In the kept T cells:
#   gamma-delta = TRDC >= 1 and TRGC1 + TRGC2 >= 1 UMI; MAIT = TRAV1-2 >= 1 and SLC4A10 >= 1 UMI; a cell with both flags counts
#   as gamma-delta. The label counts are re-summed with the flagged cells moved into the two extra labels. Nothing is printed
#   by status; 39 --gdmait reads this output.
#   usage: sbatch scripts/ws/submit_ws_a1_g01.sbatch 40 [--smoke]
#   --smoke: no object is read; synthetic label numbers, flags and counts with real sample IDs and status only (it tests
#            39 --gdmait's path); outputs end in _smoke
# This script deletes nothing and does not modify its inputs; it reads no paediatric, post-treatment or relapse object.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(here) })
SMOKE <- "--smoke" %in% commandArgs(TRUE)
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
TABWS <- file.path(TAB_DIR, "ws"); GDIR <- file.path(LARGE1_DIR, "ws_a1_g01"); SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
guard <- function(f) { if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
CHK <- list(); chk <- function(name, expected, got) { ok <- isTRUE(all.equal(expected, got))
  CHK[[length(CHK) + 1L]] <<- data.table(check = name, expected = format(expected), got = format(got), state = if (ok) "PASS" else "FAIL")
  if (!ok) message(sprintf("  %-70s expect %s got %s ** FAIL **", name, format(expected), format(got))) }

CKF <- file.path(TABWS, paste0("ws_a1_g01_checks_cells", SFX, ".csv")); OUTF <- file.path(GDIR, paste0("gdmait_counts", SFX, ".rds"))
if (!SMOKE && (file.exists(CKF) || file.exists(OUTF))) stop("outputs exist (", CKF, " or ", OUTF, "); nothing is overwritten")   # before reading ~64 GB
ck <- file.path(TABWS, paste0("ws_a1_g01_checks", SFX, ".csv"))   # 39 --known must have passed (section 6)
if (!file.exists(ck)) stop("run 39 --known first: ", ck)
CK <- fread(ck); if (!all(CK[scope == "P" & kind == "deterministic", state] == "PASS")) stop("a scope-P known answer has not passed")
S5 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]; SAMP <- S5$sample; STAT <- S5$status
PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds")); G <- PB$genes
MAP <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE); TL <- MAP[group4 == "T", bmm_fine]
L11 <- c(TL, "gamma-delta", "MAIT"); FLAGG <- c("TRDC", "TRGC1", "TRGC2", "TRAV1-2", "SLC4A10"); stopifnot(all(FLAGG %in% G))
KC <- PB$cells[kept == TRUE & group4 == "T" & sample %in% SAMP]
RAWT <- PB$raw$T[, SAMP]   # every pseudobulk object is subset to the 25 right after loading (2.1; the file also holds the paediatric samples)
FINE0 <- lapply(setNames(TL, TL), function(k) { M <- PB$fine[[k]]; if (is.null(M)) NULL else M[, intersect(colnames(M), SAMP), drop = FALSE] }); rm(PB); invisible(gc(FALSE))
FINE <- lapply(setNames(L11, L11), function(k) matrix(0, length(G), length(SAMP), dimnames = list(G, SAMP)))
NCELL <- matrix(0, length(L11), length(SAMP), dimnames = list(L11, SAMP)); WNF <- matrix(0, length(G), length(SAMP), dimnames = list(G, SAMP))
CELLS <- list()

if (SMOKE) {   # synthetic: label numbers, flags and counts invented (39's smoke recipe: naive labels halved and a P370 shift in AML); only sample IDs and status are real
  set.seed(SEED + 98L)
  prof <- (Reduce(`+`, lapply(FINE0, rowSums)) + 0.05); prof <- prof / sum(prof)                          # pooled whole-T profile, status-free
  lprof <- lapply(FINE0, function(M) { p <- (if (is.null(M)) 0 else rowSums(M)) + 0.05; p / sum(p) })      # pooled label profiles, status-free
  IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
  AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene; MF <- fread(file.path(TABWS, "ws_s5_de_T.csv"))
  g370 <- intersect(setdiff(IMP$gene, AMB), G); lf <- setNames(numeric(length(G)), G); lf[g370] <- 0.6 * setNames(MF$logFC, MF$gene)[g370]
  sh <- setNames(rep(1 / 9, 9), TL)
  for (s in SAMP) { a <- STAT[SAMP == s] == "AML"; p <- sh * (if (a) ifelse(TL %in% c("CD4 Naive", "CD8 Naive"), 0.5, 1) else 1)
    n <- as.vector(rmultinom(1, max(60, round(rlnorm(1, log(400), 0.7))), p / sum(p))); names(n) <- TL
    gdn <- rbinom(9, n, 0.03); mtn <- rbinom(9, n - gdn, 0.02); f <- if (a) 2^lf else 1
    for (k in TL) { NCELL[k, s] <- n[[k]] - gdn[k == TL] - mtn[k == TL]; FINE[[k]][, s] <- rnbinom(length(G), size = 20, mu = NCELL[k, s] * 1500 * lprof[[k]] * f) }
    NCELL["gamma-delta", s] <- sum(gdn); NCELL["MAIT", s] <- sum(mtn)
    for (k in c("gamma-delta", "MAIT")) FINE[[k]][, s] <- rnbinom(length(G), size = 20, mu = NCELL[k, s] * 1500 * prof * f)
    WNF[, s] <- Reduce(`+`, lapply(TL, function(k) FINE[[k]][, s]))
    CELLS[[s]] <- data.table(sample = s, kept_T = sum(n), gd = sum(gdn), mait = sum(mtn)) }
  chk("smoke: the 11 labels' cells sum to kept T", sum(vapply(CELLS, `[[`, 1, "kept_T")), sum(NCELL))
} else {
  suppressPackageStartupMessages(library(SeuratObject))
  AN <- file.path(LARGE1_DIR, "02_seurat_objects", "04_annotated", "GSE185381")
  for (s in SAMP) {
    o <- readRDS(file.path(AN, paste0(s, ".rds"))); X <- SeuratObject::LayerData(o[["RNA"]], layer = "counts"); rm(o)
    k <- KC[sample == s]; chk(sprintf("%s: gene axis equals the pseudobulk's", s), TRUE, identical(rownames(X), G))
    chk(sprintf("%s: every kept T cell found in the object", s), nrow(k), sum(k$cell %in% colnames(X)))
    if (!identical(rownames(X), G) || !all(k$cell %in% colnames(X))) { rm(X); next }   # the FAIL rows are in CHK; the stamped FAILED file is written below
    X <- X[, k$cell, drop = FALSE]
    gd <- X["TRDC", ] >= 1 & (X["TRGC1", ] + X["TRGC2", ]) >= 1
    mait <- X["TRAV1-2", ] >= 1 & X["SLC4A10", ] >= 1 & !gd                     # a cell with both flags counts as gamma-delta
    lab <- ifelse(gd, "gamma-delta", ifelse(mait, "MAIT", k$fine))
    for (l in unique(lab)) { j <- which(lab == l); FINE[[l]][, s] <- as.numeric(Matrix::rowSums(X[, j, drop = FALSE])); NCELL[l, s] <- length(j) }
    WNF[, s] <- as.numeric(Matrix::rowSums(X[, !(gd | mait), drop = FALSE]))
    tot <- as.numeric(Matrix::rowSums(X))
    chk(sprintf("%s: kept-T counts equal the pseudobulk's $raw$T (max |diff|)", s), 0, max(abs(tot - RAWT[, s])))
    chk(sprintf("%s: the 11 labels sum to the kept-T counts (max |diff|)", s), 0, max(abs(Reduce(`+`, lapply(L11, function(l) FINE[[l]][, s])) - tot)))
    orig <- vapply(TL, function(kk) { j <- which(k$fine == kk); if (length(j)) sum(X[, j]) else 0 }, 1)
    chk(sprintf("%s: unflagged label sums before moving equal $fine", s), 0,
        max(abs(orig - vapply(TL, function(kk) if (!is.null(FINE0[[kk]]) && s %in% colnames(FINE0[[kk]])) sum(FINE0[[kk]][, s]) else 0, 1))))
    CELLS[[s]] <- data.table(sample = s, kept_T = ncol(X), gd = sum(gd), mait = sum(mait))
    ts(sprintf("    %s done", s)); rm(X); invisible(gc(FALSE)) }
}
C <- rbindlist(CHK)[, script_md5 := MD5]
if (!all(C$state == "PASS")) {   # a failed run leaves only a stamped file, so the rerun after the fix is not blocked by the guard
  fwrite(C, file.path(TABWS, sprintf("ws_a1_g01_checks_cells_FAILED_%s%s.csv", format(Sys.time(), "%Y%m%d%H%M"), SFX)))
  ts("[40] a check FAILED; counts not written"); quit(status = 1) }
fwrite(C, guard(CKF))
saveRDS(list(fine = FINE, ncell = NCELL, wt_noflag = WNF, cells = rbindlist(CELLS), script_md5 = MD5), guard(OUTF))
ts(sprintf("[40] done: %d checks PASS; counts written", nrow(C)))
