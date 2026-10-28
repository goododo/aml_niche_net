#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds ($raw$T); results/tables/ws/ws_s5_scores.csv (the 25 T
#            samples in S5 order), ws_s8_t_imprint_ambient_prone.csv, decontx/GSE185381/ws_s1_ambient_blacklist.csv,
#            gencode_GRCh38_gene_order.txt; for --combine also ws_a1_g02_summary_replay.csv (the replayed gap_oof and LOO AUC)
# OUTPUT : --draw    results/tables/ws/ws_a1_g02_perm_labels.csv          the 1,000 relabellings (drawn once, before any is analysed)
#          --task=N  /LARGE1/.../ws_a1_g02/perm_task_NN.csv                 relabellings 20(N-1)+1 .. 20N: gap_oof and LOO AUC
#          --combine results/tables/ws/ws_a1_g02_perm_crossfit.csv, ws_a1_g02_summary_perm.csv
# what it does: PREREGISTRATION_A1_G02_reporting.md section 4, D2, nothing else. Each relabelling of the 25 status labels runs
#   S9's complete 25-fold leave-one-out cross-fit (33_s9_age_direction.R section 5, copied) and gives gap_oof_b; the
#   one-sided p is (1 + #{gap_oof_b >= observed}) / 1,001.
#   usage: Rscript scripts/ws/36_a1_g02_perm_crossfit.R --draw | --task=N | --combine   [--smoke: one relabelling, 4 folds]
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(here) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODE <- if ("--draw" %in% .a) "draw" else if ("--combine" %in% .a) "combine" else if (any(grepl("^--task=", .a))) "task" else stop("give --draw, --task=N or --combine")
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
TABWS <- file.path(TAB_DIR, "ws"); TASKDIR <- file.path(LARGE1_DIR, "ws_a1_g02"); dir.create(TASKDIR, showWarnings = FALSE)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
N_PERM <- 1000L; PER_TASK <- 20L; FDR_CUT <- 0.10; LFC_CUT <- 0.5
S5 <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]
D5 <- S5[, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
LABFILE <- file.path(TABWS, paste0("ws_a1_g02_perm_labels", if (SMOKE) "_smoke" else "", ".csv"))

if (MODE == "draw") {   # all relabellings drawn in one process, in the declared way, and written before any is analysed
  if (file.exists(LABFILE) && !SMOKE) stop("relabellings already drawn: ", LABFILE, " (never redraw)")
  set.seed(SEED)
  lab <- vapply(seq_len(N_PERM), function(i) paste(as.integer(sample(D5$status) == "AML"), collapse = ""), "")
  true <- paste(as.integer(D5$status == "AML"), collapse = "")
  P <- data.table(perm = seq_len(N_PERM), labels = lab)
  P[, true_AML_labelled_AML := vapply(strsplit(labels, ""), function(x) sum(as.integer(x) == 1L & D5$status == "AML"), 1L)]
  fwrite(P, LABFILE)
  ts(sprintf("[draw] %d relabellings written; identity drawn %d times; overlap with the true AML set: %s", N_PERM, sum(lab == true),
             paste(names(table(P$true_AML_labelled_AML)), table(P$true_AML_labelled_AML), sep = ":", collapse = " ")))
  quit(status = 0)
}

if (MODE == "combine") {
  fs <- list.files(TASKDIR, pattern = "^perm_task_[0-9]+\\.csv$", full.names = TRUE)
  R <- rbindlist(lapply(fs, fread))[order(perm)]
  stopifnot(nrow(R) == N_PERM, identical(R$perm, seq_len(N_PERM)))
  OBS <- fread(file.path(TABWS, "ws_a1_g02_summary_replay.csv"))
  gap_obs <- OBS[quantity == "R2_gap_oof", value]; auc_obs <- OBS[quantity == "R2_pooled_LOO_AUC", value]
  p_gap <- (1 + sum(R$gap_oof >= gap_obs, na.rm = TRUE)) / (N_PERM + 1)
  p_auc <- (1 + sum(R$loo_auc >= auc_obs, na.rm = TRUE)) / (N_PERM + 1)
  fwrite(R, file.path(TABWS, "ws_a1_g02_perm_crossfit.csv"))
  S <- data.table(quantity = c("D2_gap_oof_observed", "D2_p_gap_oof", "D2_perm_gap_median", "D2_perm_gap_P95", "D2_perm_gap_max", "D2_n_perm_gap_ge_observed",
                               "D2_relabellings_with_failed_folds", "D2_failed_folds_total", "D2_loo_auc_observed", "D2_p_loo_auc", "D2_n_perm_ge_1.378"),
                  value = c(gap_obs, p_gap, median(R$gap_oof, na.rm = TRUE), quantile(R$gap_oof, .95, na.rm = TRUE), max(R$gap_oof, na.rm = TRUE),
                            sum(R$gap_oof >= gap_obs, na.rm = TRUE), sum(R$n_failed > 0), sum(R$n_failed), auc_obs, p_auc, sum(R$gap_oof >= 1.378, na.rm = TRUE)))
  fwrite(S, file.path(TABWS, "ws_a1_g02_summary_perm.csv")); print(S)
  print(R[order(-gap_oof)][1:10])
  quit(status = 0)
}

## -- task: one slice of relabellings --------------------------------------------------------------------------------
TASK <- as.integer(sub("^--task=", "", grep("^--task=", .a, value = TRUE)))
P <- fread(LABFILE, colClasses = c(labels = "character"))
rows <- if (SMOKE) 1L else ((TASK - 1L) * PER_TASK + 1L):(TASK * PER_TASK)
stopifnot(all(rows <= nrow(P)))
PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
fit_de <- function(counts, d, terms, coef) {        # 27's edgeR branch, verbatim (as in 35)
  keep_g <- !(rownames(counts) %in% BLACK)
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = d$status)
  keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(nrow(X) == nrow(d), qr(X)$rank == ncol(X))
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
  r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, FDR)]
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]; r }
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")
G_RAW <- PB$raw$T[, D5$sample]
Mf <- G_RAW[!rownames(G_RAW) %in% BLACK, , drop = FALSE]; y0 <- DGEList(Mf)                       # 33's logz: the label-free log-CPM
y0 <- calcNormFactors(y0[rowMeans(cpm(y0) > 1) >= 0.5, , keep.lib.sizes = FALSE]); L <- cpm(y0, log = TRUE, prior.count = 1)
folds <- if (SMOKE) c(1L, 2L, 16L, 17L) else seq_len(nrow(D5))
auc <- function(s, st) { a <- s[st == "AML"]; h <- s[st == "healthy"]; mean(outer(a, h, ">") + 0.5 * outer(a, h, "=="), na.rm = TRUE) }
RES <- rbindlist(lapply(rows, function(b) {
  Dp <- copy(D5); Dp[, status := factor(ifelse(strsplit(P$labels[b], "")[[1]] == "1", "AML", "healthy"), c("healthy", "AML"))]
  stopifnot(Dp[status == "AML", .N] == 15L)
  oof <- vapply(folds, function(i) tryCatch({
    d_tr <- Dp[-i]; terms <- if (d_tr[, uniqueN(sex)] < 2) setdiff(MAIN, "sex") else MAIN
    r <- fit_de(G_RAW, d_tr, terms, "statusAML")
    fu <- setdiff(r[hit == TRUE & logFC > 0, gene], AMB); fd <- setdiff(r[hit == TRUE & logFC < 0, gene], AMB)
    Ltr <- L[, -i, drop = FALSE]; z <- (L[, i] - rowMeans(Ltr)) / apply(Ltr, 1, sd); z <- z[is.finite(z)]
    (if (length(intersect(fu, names(z)))) mean(z[intersect(fu, names(z))]) else 0) - (if (length(intersect(fd, names(z)))) mean(z[intersect(fd, names(z))]) else 0)
  }, error = function(e) NA_real_), 1)
  st <- Dp$status[folds]
  x <- data.table(perm = b, gap_oof = mean(oof[st == "AML"], na.rm = TRUE) - mean(oof[st == "healthy"], na.rm = TRUE),
                  loo_auc = auc(oof, st), n_failed = sum(is.na(oof)), true_AML_labelled_AML = P$true_AML_labelled_AML[b])
  ts(sprintf("    relabelling %4d: gap_oof %.3f, LOO AUC %.3f, failed folds %d", b, x$gap_oof, x$loo_auc, x$n_failed)); x }))
fn <- file.path(TASKDIR, sprintf("perm_task_%02d%s.csv", TASK, if (SMOKE) "_smoke" else ""))
fwrite(RES, fn); ts("[task ", TASK, "] written ", fn, " (", nrow(RES), " relabellings)")
