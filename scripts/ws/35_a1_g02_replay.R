#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds ($raw$T); results/tables/ws/ws_s5_scores.csv (the 25 T
#            samples in S5 order), ws_s5_imprint_genes.csv, ws_s5_known_answers.csv, ws_s8_t_imprint_ambient_prone.csv,
#            decontx/GSE185381/ws_s1_ambient_blacklist.csv, ws_s8_scores_rep2.csv; gencode_GRCh38_gene_order.txt;
#            meta_GSE185381_v2.3.csv and the GSE185381 GEO author metadata (capture dates); the S9 log (known answers only)
# OUTPUT : results/tables/ws/ws_a1_g02_{s5_perm_draws,oof_scores,fold_genes,gene_selection,lpo_pairs,descriptive,
#            summary_replay,checks_replay}.csv
# what it does: PREREGISTRATION_A1_G02_reporting.md sections 3, 4 (D1) and 5, nothing else. R1 replays the S5 T label
#   permutation by replaying 27's random-number sequence; R2 replays S9's 25-fold cross-fit and saves what S9 did not;
#   D1 is the declared leave-pair-out AUC; section 5 recomputes three descriptive numbers from durable tables.
#   It never re-runs 27 or 33 and writes only new ws_a1_g02_* files.
#   usage: Rscript scripts/ws/35_a1_g02_replay.R [--smoke]   (--smoke: 10 permutations, 4 folds, 5 pairs; outputs *_smoke)
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(here) })
SMOKE <- "--smoke" %in% commandArgs(TRUE)
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)                                   # exactly as 27_s5_imprint.R line 27: R1 depends on this being the first draw
TABWS <- file.path(TAB_DIR, "ws"); SFX <- if (SMOKE) "_smoke" else ""
out <- function(x) file.path(TABWS, paste0("ws_a1_g02_", x, SFX, ".csv"))
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list(); chk <- function(name, expected, got) { ok <- isTRUE(all.equal(expected, got))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, expected, got, if (ok) "PASS" else "** FAIL **")) }
FDR_CUT <- 0.10; LFC_CUT <- 0.5; N_PLANT <- 200L; GAP_OOF_S9 <- 1.378
N_PERM <- if (SMOKE) 10L else 200L

## -- 1. inputs, blacklist and DE engine (27_s5_imprint.R lines 33-88, verbatim) ------------------------------------
PB  <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
BL  <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO  <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
fit_de <- function(counts, d, terms, coef) {        # 27's edgeR branch: blacklist first, filterByExpr + CPM > 1 in half, TMM, QL robust
  keep_g <- !(rownames(counts) %in% BLACK)
  y <- DGEList(counts[keep_g, d$sample, drop = FALSE], group = d$status)
  keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- y[keep, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y, method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(nrow(X) == nrow(d), qr(X)$rank == ncol(X))
  cn <- grep(coef, colnames(X), value = TRUE)[1]
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = cn)
  r <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, PValue, FDR)]
  r[, hit := FDR < FDR_CUT & abs(logFC) > LFC_CUT]; r }
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")

S5  <- fread(file.path(TABWS, "ws_s5_scores.csv"))[compartment == "T"]          # file order = 27's design order (AML first, then healthy)
D5  <- S5[, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
G_RAW <- PB$raw$T[, D5$sample]
ts(sprintf("[1] %d samples (%d AML, %d healthy); %d genes on the axis; blacklist %d", nrow(D5), D5[status == "AML", .N], D5[status == "healthy", .N], nrow(G_RAW), length(BLACK)))
print(D5[c(1:3, 23:25)])

## -- 2. R1: S5 T label permutation, replayed --------------------------------------------------------------------------
ts("[2] R1: main fit, then the planted-DE draws, then the label permutations (27's order)")
r0 <- fit_de(G_RAW, D5, MAIN, "statusAML")
up <- r0[hit == TRUE & logFC > 0, gene]; dn <- r0[hit == TRUE & logFC < 0, gene]
chk("R1 main fit reproduces the 402 S5 T imprint genes (same set)", TRUE, setequal(c(up, dn), IMP$gene))
pool <- setdiff(r0$gene, c(up, dn))
for (lfc in c(1, 1.5)) invisible(sample(pool, min(N_PLANT, length(pool))))       # the two planted-DE draws of 27 line 112 (only the RNG state matters)
meta <- fread(file.path(RAW_DIR, "metadata", "GSE185381", "meta_GSE185381_v2.3.csv"), colClasses = "character")[, sid := sub("^GSE185381__", "", sample_id)]
mf <- list.files(file.path(RAW_DIR, "public", "GEO", "GSE185381_RAW"), pattern = "_metadata.csv.gz$", full.names = TRUE)
cap <- rbindlist(lapply(mf, function(f) { d <- fread(f, fill = TRUE)                  # some author files have no inflammation_group column
  if (!"inflammation_group" %in% names(d)) d[, inflammation_group := NA_character_]
  d[, .N, by = .(sample = samples, inflammation_group)][, date := regmatches(basename(f), regexpr("20[12][0-9]-[01][0-9]-[0-3][0-9]", basename(f)))] }))
cap <- cap[sample %in% D5$sample, .(dates = paste(sort(unique(date)), collapse = ";"), inflammation_group = paste(setdiff(unique(inflammation_group), NA), collapse = "/")), by = sample]
cap[, only_2019 := !grepl("202[0-9]", dates)]
SAMPLEINFO <- merge(D5[, .(sample, status, sex, subtype_share)], cap, by = "sample", all.x = TRUE)
SAMPLEINFO[, pooled := meta$is_multiplexed[match(sample, meta$sid)] == "TRUE"]
SAMPLEINFO <- SAMPLEINFO[match(D5$sample, sample)]                               # back to design order
chk("capture date found for all 25 samples", 25L, SAMPLEINFO[!is.na(dates), .N])
print(SAMPLEINFO[, .(status, n = .N, only_2019 = sum(only_2019), pooled = sum(pooled), female = sum(sex == "female")), by = status])
perm <- rbindlist(lapply(seq_len(N_PERM), function(i) {
  dp <- copy(D5); dp[, status := sample(status)]                                   # 27 line 132, same call on the same 25-row design
  n_hit <- sum(fit_de(G_RAW, dp, MAIN, "statusAML")$hit)
  a <- dp$status == "AML"
  data.table(draw = i, n_hits = n_hit, labels = paste(as.integer(a), collapse = ""), true_AML_labelled_AML = sum(a & D5$status == "AML"),
             labelled_AML_only_2019 = sum(a & SAMPLEINFO$only_2019), labelled_AML_pooled = sum(a & SAMPLEINFO$pooled),
             labelled_AML_female = sum(a & SAMPLEINFO$sex == "female"), labelled_AML_mean_cd8 = mean(SAMPLEINFO$subtype_share[a])) }))
ts(sprintf("    permutation hit counts: median %.1f, P95 %.1f, max %d", median(perm$n_hits), quantile(perm$n_hits, .95), max(perm$n_hits)))
if (!SMOKE) {
  chk("R1 replay: median of the 200 hit counts", 0, median(perm$n_hits))
  chk("R1 replay: P95 (type 7) of the 200 hit counts", 11.1, round(unname(quantile(perm$n_hits, .95)), 1))
  chk("R1 replay: max of the 200 hit counts", 1615L, as.integer(max(perm$n_hits)))
}
P_R1 <- (1 + sum(perm$n_hits >= length(c(up, dn)))) / (N_PERM + 1)
fwrite(perm, out("s5_perm_draws")); print(perm[order(-n_hits)][1:5])

## -- 3. R2: S9 cross-fit, replayed; D1: leave-pair-out ------------------------------------------------------------
logz <- function(M) {   # 33's logz (blacklist; CPM > 1 in >= half; TMM; log CPM prior 1); only L is used here
  M <- M[!rownames(M) %in% BLACK, , drop = FALSE]; y <- DGEList(M)
  keep <- rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  cpm(y, log = TRUE, prior.count = 1) }
L <- logz(G_RAW); stopifnot(identical(colnames(L), D5$sample))
heldout_scores <- function(D, h) {   # refit without the samples h, score them with the fold imprint; z from the other samples (33 lines 192-198)
  d_tr <- D[-h]; terms <- if (d_tr[, uniqueN(sex)] < 2) setdiff(MAIN, "sex") else MAIN
  r <- fit_de(G_RAW, d_tr, terms, "statusAML")
  fu <- setdiff(r[hit == TRUE & logFC > 0, gene], AMB); fd <- setdiff(r[hit == TRUE & logFC < 0, gene], AMB)
  Ltr <- L[, -h, drop = FALSE]; mu <- rowMeans(Ltr); sg <- apply(Ltr, 1, sd)
  sc <- vapply(h, function(j) { z <- (L[, j] - mu) / sg; z <- z[is.finite(z)]
    (if (length(intersect(fu, names(z)))) mean(z[intersect(fu, names(z))]) else 0) - (if (length(intersect(fd, names(z)))) mean(z[intersect(fd, names(z))]) else 0) }, 1)
  list(scores = sc, up = fu, down = fd) }
ts("[3] R2: 25-fold leave-one-out cross-fit")
folds <- if (SMOKE) c(1L, 2L, 16L, 17L) else seq_len(nrow(D5))   # smoke: two AML and two healthy folds, so the gap is defined
CF <- lapply(folds, function(i) { x <- heldout_scores(D5, i); ts(sprintf("    fold %2d (%s): %d genes, oof %.3f", i, D5$sample[i], length(x$up) + length(x$down), x$scores)); x })
OOF <- data.table(fold = folds, sample = D5$sample[folds], status = D5$status[folds], oof_score = vapply(CF, `[[`, 1, "scores"),
                  n_up = vapply(CF, function(x) length(x$up), 1L), n_down = vapply(CF, function(x) length(x$down), 1L))
logf <- here::here("00_project", "logs", "ws_s9_age_3669091.err")
lg <- if (file.exists(logf)) grep("oof score", readLines(logf), value = TRUE) else character(0)
logged <- data.table(sample = sub("^.*fold +[0-9]+ \\(([^)]*)\\).*$", "\\1", lg), logged = as.numeric(sub("^.*oof score ", "", lg)))   # the sample is in the first parentheses after "fold N"
OOF <- merge(OOF, logged, by = "sample", all.x = TRUE, sort = FALSE)[order(fold)]
chk("R2: logged out-of-fold scores found for every replayed fold", length(folds), OOF[!is.na(logged), .N])
chk("R2: max |replayed - logged out-of-fold score| < 0.001", TRUE, isTRUE(max(abs(OOF$oof_score - OOF$logged)) < 0.001))
FG <- rbindlist(lapply(seq_along(folds), function(k) data.table(fold = folds[k], left_out = D5$sample[folds[k]],
       gene = c(CF[[k]]$up, CF[[k]]$down), direction = rep(c("up", "down"), c(length(CF[[k]]$up), length(CF[[k]]$down))))))
GAP_OOF <- OOF[status == "AML", mean(oof_score)] - OOF[status == "healthy", mean(oof_score)]
auc <- function(s, st) { a <- s[st == "AML"]; h <- s[st == "healthy"]; mean(outer(a, h, ">") + 0.5 * outer(a, h, "==")) }
if (!SMOKE) {
  chk("R2: gap_oof reproduces S9 (1.378, within 0.001)", TRUE, abs(GAP_OOF - GAP_OOF_S9) < 0.001)
  AUC_LOO <- auc(OOF$oof_score, OOF$status)
  set.seed(SEED); bt <- replicate(2000, { i <- c(sample(which(OOF$status == "AML"), replace = TRUE), sample(which(OOF$status == "healthy"), replace = TRUE)); auc(OOF$oof_score[i], OOF$status[i]) })
  SEL <- FG[, .(n_folds = .N, direction = direction[1]), by = gene][, freq := n_folds / 25][order(-n_folds)]
  SEL[, in_primary_370 := gene %in% setdiff(IMP$gene, AMB)]
  fwrite(SEL, out("gene_selection"))
  ts(sprintf("    gap_oof %.4f; pooled LOO AUC %.3f (bootstrap 95%% CI %.3f-%.3f, conditional); genes in all 25 folds: %d", GAP_OOF, AUC_LOO, quantile(bt, .025), quantile(bt, .975), SEL[n_folds == 25, .N]))
}
fwrite(OOF, out("oof_scores")); fwrite(FG, out("fold_genes")); print(OOF)

ts("[4] D1: leave-pair-out AUC over the 15 x 10 pairs")
PAIRS <- CJ(a = which(D5$status == "AML"), h = which(D5$status == "healthy")); if (SMOKE) PAIRS <- PAIRS[1:5]
LPO <- rbindlist(lapply(seq_len(nrow(PAIRS)), function(k) { x <- heldout_scores(D5, c(PAIRS$a[k], PAIRS$h[k]))
  data.table(aml = D5$sample[PAIRS$a[k]], healthy = D5$sample[PAIRS$h[k]], score_aml = x$scores[1], score_healthy = x$scores[2],
             n_genes = length(x$up) + length(x$down)) }))
LPO[, correct := (score_aml > score_healthy) + 0.5 * (score_aml == score_healthy)]
AUC_LPO <- mean(LPO$correct); fwrite(LPO, out("lpo_pairs"))
ts(sprintf("    AUC_LPO %.4f over %d pairs; pairs with an empty fold imprint: %d", AUC_LPO, nrow(LPO), LPO[n_genes == 0, .N]))

## -- 4. section 5: durable recomputation of three descriptive numbers --------------------------------------------------
R8 <- fread(file.path(TABWS, "ws_s8_scores_rep2.csv"))                            # GSE227903 diagnosis samples only (role AML_Dx)
stopifnot(all(R8$role == "AML_Dx"), all(grepl("_Dg$", R8$sample)))
R8[, cd8_share := n_cd8 / n_T]
s10 <- c("1216_Dg", "1886_Dg", "3904_Dg", "40389_Dg", "4978_Dg", "5750_Dg")
DESC <- data.table(
  quantity = c("rho_S5score_CD8share_AML", "rho_S5score_CD8share_healthy", "rho_S8rep2score_CD8share_GSE227903Dx_n7", "rho_S8rep2score_CD8share_GSE227903Dx_n6",
               "healthy_captured_only_2019", "AML_captured_only_2019", "AML_Inflammation_high_captured_only_2019"),
  value = c(S5[status == "AML", cor(score, subtype_share, method = "spearman")], S5[status == "healthy", cor(score, subtype_share, method = "spearman")],
            R8[, cor(score_primary, cd8_share, method = "spearman")], R8[sample %in% s10, cor(score_primary, cd8_share, method = "spearman")],
            SAMPLEINFO[status == "healthy", sum(only_2019)], SAMPLEINFO[status == "AML", sum(only_2019)], SAMPLEINFO[inflammation_group == "Inflammation_high", sum(only_2019)]),
  of = c(15, 10, 7, 6, 10, 15, SAMPLEINFO[inflammation_group == "Inflammation_high", .N]))
fwrite(DESC, out("descriptive")); print(DESC)

## -- 5. summary ----------------------------------------------------------------------------------------------------
SUM <- data.table(quantity = c("R1_real_imprint_size", "R1_perm_median", "R1_perm_P95", "R1_perm_max", "R1_monte_carlo_p", "R2_gap_oof", "R2_pooled_LOO_AUC",
                               "R2_pooled_LOO_AUC_boot_lo", "R2_pooled_LOO_AUC_boot_hi", "R2_genes_in_all_folds", "D1_AUC_LPO", "D1_pairs", "D1_pairs_empty_imprint"),
                  value = c(length(c(up, dn)), median(perm$n_hits), unname(quantile(perm$n_hits, .95)), max(perm$n_hits), P_R1, GAP_OOF,
                            if (SMOKE) NA else AUC_LOO, if (SMOKE) NA else quantile(bt, .025), if (SMOKE) NA else quantile(bt, .975),
                            if (SMOKE) NA else SEL[n_folds == 25, .N], AUC_LPO, nrow(LPO), LPO[n_genes == 0, .N]))
fwrite(SUM, out("summary_replay")); print(SUM)
C <- rbindlist(CHK); fwrite(C, out("checks_replay"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(C), sum(!C$pass)))
