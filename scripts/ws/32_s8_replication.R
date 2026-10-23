#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_pseudobulk/<dataset>/ws_pb_<dataset>.rds (26: per-sample T pseudobulk, kept-cell table)
#          results/tables/ws/ws_s4_pseudobulk_meta_<dataset>.csv (role, paediatric, ambient shares)
#          /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds (raw counts: CD4/CD8 pseudobulks, per-cell test)
#          results/tables/ws/ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv, decontx/GSE185381/ws_s1_ambient_blacklist.csv,
#          /LARGE1/.../reference/gencode_GRCh38_gene_order.txt; scripts/ws/PREREGISTRATION_S8_replication.md (tag prereg-s8-v1)
# OUTPUT : results/tables/ws/ws_s8_scores_<set>.csv   sample x gene set: T score, CD4 and CD8 scores, cell counts, ambient shares
#          results/tables/ws/ws_s8_tests.csv          set x gene set x test: statistic, random-set P5/P95, n, pass
#          results/tables/ws/ws_s8_known_answers.csv, ws_s8_checks.csv; scripts/ws/GATES.md block "A1-S8"
# what it does: PREREGISTRATION_S8_replication.md sections 2-4. For each cohort set (disc = GSE185381 adult AML, the positive
#   control; rep = GSE289435, decisive; rep2 = GSE227903 diagnosis; desc = GSE116256) and each gene set (primary 370, full 402):
#   T2 Spearman between the up-half and down-half means across patients (pass: < 0 and below the random-set P5);
#   T3 Spearman between CD4-based and CD8-based patient scores (pass: > 0.5 and above the random-set P95; >= 6 patients);
#   T1 per-patient coherence of the imprint genes in single T cells vs matched random sets (pass: above P95 in >= 6 of 8).
#   Scoring and the expression-matched random sets are those of S5 (27/28). Known answers: disc must pass; 100 sign-scrambled
#   imprints must pass T2 at most 5 times on rep. GATE_S8 = T2 and T3 on rep, primary set (T2 and T1 if T3 is not testable).
#   usage: Rscript scripts/ws/32_s8_replication.R [--sets=disc,rep,rep2,desc] [--smoke]   (--smoke: 100 / 20 / 10 draws)
# This script deletes nothing and modifies none of its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(edgeR); library(Seurat); library(SeuratObject); library(here) })
.a <- commandArgs(TRUE); .get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
SETS_ARG <- strsplit(.get("sets", "disc,rep,rep2,desc"), ",")[[1]]; SMOKE <- "--smoke" %in% .a
N_RAND <- if (SMOKE) 100L else 1000L; N_RAND_CELL <- if (SMOKE) 20L else 200L; N_SCR <- if (SMOKE) 10L else 100L; CAP <- 500L; MINC <- 20L; DET_MIN <- 0.10; MIN_HALF <- 20L
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "ws", "config_ws.R")); set.seed(SEED)
options(error = function() { calls <- sys.calls(); n <- length(calls); message("--- TRACEBACK ---")
  for (i in seq(max(1, n - 15), n)) message(i, ": ", substr(paste(deparse(calls[[i]]), collapse = " "), 1, 220)); q(status = 1) })
TABWS <- file.path(TAB_DIR, "ws"); GATES <- here::here("scripts", "ws", "GATES.md")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list(); chk <- function(name, expected, got) { got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-78s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **")) }
SETS <- list(disc = "GSE185381", rep = "GSE289435", rep2 = "GSE227903", desc = "GSE116256")[SETS_ARG]

## -- 1. gene sets and the blacklist (as S5) -------------------------------------------------------------------------
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))
GS <- list(primary = IMP[!gene %in% AMB$gene], full = IMP)
chk("primary set is 370 genes (156 up / 214 down)", "370/156/214", GS$primary[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))])
chk("full set is 402 genes (171 up / 231 down)", "402/171/231", GS$full[, sprintf("%d/%d/%d", .N, sum(direction == "up"), sum(direction == "down"))])
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
black_for <- function(G) unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene], "XIST", GO[chr == "chrY", gene],
                                  grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))

## -- 2. scoring helpers (27/28) ---------------------------------------------------------------------------------------
logz <- function(M) {   # M: genes x samples raw sums; blacklist removed; returns z of logCPM over the set's samples and the matching universe
  M <- M[!rownames(M) %in% black_for(rownames(M)), , drop = FALSE]; y <- DGEList(M)
  keep <- filterByExpr(y) & rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  L <- cpm(y, log = TRUE, prior.count = 1); Z <- t(scale(t(L)))
  U <- data.table(gene = rownames(L), mean_lcpm = rowMeans(L), det = rowMeans(M[rownames(L), , drop = FALSE] > 0))
  U[, bin := paste(cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE), cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf)))]
  list(Z = Z, U = U) }
score <- function(Z, up, dn) { up <- intersect(up, rownames(Z)); dn <- intersect(dn, rownames(Z))
  (if (length(up)) colMeans(Z[up, , drop = FALSE]) else 0) - (if (length(dn)) colMeans(Z[dn, , drop = FALSE]) else 0) }
matched <- function(U, genes) { tb <- table(U[gene %in% genes, bin])   # one random tested gene per imprint gene, same expression/detection bin
  unlist(lapply(names(tb), function(b) { pool <- U[bin == b & !gene %in% genes, gene]; if (!length(pool)) pool <- U[!gene %in% genes, gene]; sample(pool, tb[[b]], replace = length(pool) < tb[[b]]) })) }
t2_stat <- function(Z, up, dn) { up <- intersect(up, rownames(Z)); dn <- intersect(dn, rownames(Z)); if (length(up) < 5 || length(dn) < 5 || ncol(Z) < 4) return(NA_real_)
  cor(colMeans(Z[up, , drop = FALSE]), colMeans(Z[dn, , drop = FALSE]), method = "spearman") }
split_random <- function(U, up, dn) { g <- matched(U, c(up, dn)); list(up = head(g, length(intersect(up, U$gene))), dn = tail(g, length(intersect(dn, U$gene)))) }
coherence <- function(C, up, dn) { if (length(up) < 2 || length(dn) < 1) return(NA_real_); A <- C[up, up]; mean(A[upper.tri(A)]) - mean(C[up, dn]) }

## -- 3. per set: cells, pseudobulks, per-cell matrices ----------------------------------------------------------------
TESTS <- list(); KA <- list(); RESULT <- list()
for (set in names(SETS)) { ds <- SETS[[set]]; ts(sprintf("[%s] %s", set, ds))
  fpb <- file.path(LARGE1_DIR, "ws_pseudobulk", ds, paste0("ws_pb_", ds, ".rds")); fme <- file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", ds, ".csv"))
  if (!file.exists(fpb) || !file.exists(fme)) { ts("    pseudobulk object or meta missing; set skipped"); next }
  P <- readRDS(fpb); META <- fread(fme); M <- unique(META[, .(sample, role, paediatric, hb = hb_share_nonery, gran = gran_share_lymph)])
  print(M[, .N, by = .(role, paediatric)])
  elig <- P$cells[kept == TRUE & group4 == "T", .(n_T = .N, n_cd4 = sum(group == "CD4_T"), n_cd8 = sum(group == "CD8_T")), by = sample]
  S <- merge(M[role == "AML_Dx" & paediatric == FALSE], elig[n_T >= MINC], by = "sample")
  ts(sprintf("    adult AML_Dx with >= %d kept T cells: %d (CD4 and CD8 both >= %d: %d)", MINC, nrow(S), MINC, S[n_cd4 >= MINC & n_cd8 >= MINC, .N]))
  if (nrow(S) < 4) { ts("    fewer than 4 patients; set skipped"); next }
  AN <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", ds); G <- P$genes
  PB_T <- P$raw$T[, S$sample, drop = FALSE]; PB4 <- PB8 <- matrix(0, length(G), nrow(S), dimnames = list(G, S$sample)); CELL <- list(); maxdiff <- 0
  for (sm in S$sample) { o <- readRDS(file.path(AN, paste0(sm, ".rds"))); X <- LayerData(o, assay = "RNA", layer = "counts"); rm(o)
    cc <- P$cells[kept == TRUE & group4 == "T" & sample == sm]; X <- X[, cc$cell, drop = FALSE]; stopifnot(identical(rownames(X), G))
    maxdiff <- max(maxdiff, max(abs(Matrix::rowSums(X) - PB_T[, sm])))                                   # 26's T pseudobulk must be these cells' sum
    if (any(cc$group == "CD4_T")) PB4[, sm] <- Matrix::rowSums(X[, cc$group == "CD4_T", drop = FALSE]); if (any(cc$group == "CD8_T")) PB8[, sm] <- Matrix::rowSums(X[, cc$group == "CD8_T", drop = FALSE])
    j <- if (ncol(X) > CAP) sample.int(ncol(X), CAP) else seq_len(ncol(X)); Xc <- X[, j, drop = FALSE]                  # per-cell test: capped, log-normalised
    det <- Matrix::rowMeans(Xc > 0); keepg <- names(det)[det >= DET_MIN]; Xn <- log1p(t(t(as.matrix(Xc[keepg, , drop = FALSE])) / colSums(Xc)) * 1e4)
    CELL[[sm]] <- list(X = Xn, det = det[keepg]); rm(X, Xc); gc(FALSE) }
  chk(sprintf("[%s] 26's T pseudobulk equals the kept T cells' sum (max abs diff)", set), 0, maxdiff)
  TZ <- logz(PB_T); S4 <- S[n_cd4 >= MINC & n_cd8 >= MINC, sample]
  Z4 <- if (length(S4) >= 6) logz(PB4[, S4, drop = FALSE])$Z else NULL; Z8 <- if (length(S4) >= 6) logz(PB8[, S4, drop = FALSE])$Z else NULL
  SC <- copy(S)
  for (gsn in names(GS)) { up <- GS[[gsn]][direction == "up", gene]; dn <- GS[[gsn]][direction == "down", gene]
    cov <- mean(c(up, dn) %in% TZ$U$gene); if (gsn == "primary") chk(sprintf("[%s] primary genes in the tested universe >= 0.8", set), TRUE, cov >= 0.8)
    s <- score(TZ$Z, up, dn); SC[, (paste0("score_", gsn)) := s[sample]]
    ## T2
    obs2 <- t2_stat(TZ$Z, up, dn); nul2 <- replicate(N_RAND, { r <- split_random(TZ$U, up, dn); t2_stat(TZ$Z, r$up, r$dn) })
    p5 <- unname(quantile(nul2, 0.05, na.rm = TRUE)); pass2 <- isTRUE(obs2 < 0 && obs2 < p5)
    TESTS[[length(TESTS) + 1]] <- data.table(set = set, dataset = ds, gene_set = gsn, test = "T2_opposite_movement", statistic = obs2, null_P5 = p5, null_P95 = unname(quantile(nul2, .95, na.rm = TRUE)), n = nrow(S), pass = pass2, detail = sprintf("coverage %.2f; valid draws %d", cov, sum(!is.na(nul2))))
    ## T3
    if (!is.null(Z4)) { s4 <- score(Z4, up, dn); s8 <- score(Z8, up, dn); obs3 <- cor(s4, s8, method = "spearman")
      nul3 <- replicate(N_RAND, { r <- split_random(TZ$U, up, dn); cor(score(Z4, r$up, r$dn), score(Z8, r$up, r$dn), method = "spearman") })
      p95 <- unname(quantile(nul3, 0.95, na.rm = TRUE)); pass3 <- isTRUE(obs3 > 0.5 && obs3 > p95)
      if (gsn == "primary") SC[match(S4, sample), `:=`(score_cd4 = s4, score_cd8 = s8)]
      TESTS[[length(TESTS) + 1]] <- data.table(set = set, dataset = ds, gene_set = gsn, test = "T3_subtype_concordance", statistic = obs3, null_P5 = unname(quantile(nul3, .05, na.rm = TRUE)), null_P95 = p95, n = length(S4), pass = pass3, detail = sprintf("patients with CD4 and CD8 >= %d: %d", MINC, length(S4)))
    } else { pass3 <- NA; TESTS[[length(TESTS) + 1]] <- data.table(set = set, dataset = ds, gene_set = gsn, test = "T3_subtype_concordance", statistic = NA_real_, null_P5 = NA_real_, null_P95 = NA_real_, n = length(S4), pass = NA, detail = "not testable (< 6 patients with both subtypes)") }
    ## T1
    per <- rbindlist(lapply(names(CELL), function(sm) { cl <- CELL[[sm]]; u <- intersect(up, rownames(cl$X)); d <- intersect(dn, rownames(cl$X))
      if (length(u) < MIN_HALF || length(d) < MIN_HALF) return(data.table(sample = sm, obs = NA_real_, p95 = NA_real_, n_up = length(u), n_dn = length(d)))
      Ud <- data.table(gene = names(cl$det), det = cl$det); Ud[, bin := cut(det, unique(quantile(det, 0:10 / 10)), include.lowest = TRUE)]
      C <- cor(t(cl$X[c(u, d), ])); obs <- coherence(C, u, d)
      nul <- replicate(N_RAND_CELL, { g <- matched(Ud, c(u, d)); gu <- head(g, length(u)); gd <- tail(g, length(d)); coherence(cor(t(cl$X[c(gu, gd), ])), gu, gd) })
      data.table(sample = sm, obs = obs, p95 = unname(quantile(nul, 0.95, na.rm = TRUE)), n_up = length(u), n_dn = length(d)) }))
    per[, above := obs > p95]; n_ok <- per[!is.na(above), .N]; n_above <- per[above == TRUE, .N]; pass1 <- isTRUE(n_ok >= 1 && n_above >= min(6L, n_ok) && n_above / n_ok >= 0.75)
    if (gsn == "primary") SC[, `:=`(t1_coherence = per[match(SC$sample, sample), obs], t1_random_P95 = per[match(SC$sample, sample), p95])]
    TESTS[[length(TESTS) + 1]] <- data.table(set = set, dataset = ds, gene_set = gsn, test = "T1_within_patient_coherence", statistic = n_above / max(1, n_ok), null_P5 = NA_real_, null_P95 = NA_real_, n = n_ok, pass = pass1, detail = sprintf("%d of %d patients above their random P95 (median obs %.3f, median P95 %.3f)", n_above, n_ok, median(per$obs, na.rm = TRUE), median(per$p95, na.rm = TRUE)))
    ## covariates (informational)
    TESTS[[length(TESTS) + 1]] <- data.table(set = set, dataset = ds, gene_set = gsn, test = "info_score_vs_ambient", statistic = NA_real_, null_P5 = NA_real_, null_P95 = NA_real_, n = nrow(S), pass = NA,
      detail = sprintf("Spearman score~hb %.2f, score~gran %.2f", cor(s[S$sample], S$hb, method = "spearman", use = "complete.obs"), cor(s[S$sample], S$gran, method = "spearman", use = "complete.obs")))
    ts(sprintf("    %-7s T2 %.2f (P5 %.2f) %s | T3 %s | T1 %d/%d %s", gsn, obs2, p5, pass2, if (is.na(pass3)) "n/a" else sprintf("%.2f (P95 %.2f) %s", obs3, p95, pass3), n_above, n_ok, pass1))
    RESULT[[paste(set, gsn)]] <- list(T2 = pass2, T3 = pass3, T1 = pass1)
    ## sign-scrambled imprint (known answer 2) on the decisive set, primary gene set
    if (set == "rep" && gsn == "primary") { scr <- replicate(N_SCR, { lab <- sample(c(rep("up", length(up)), rep("down", length(dn)))); g <- c(up, dn); t2_stat(TZ$Z, g[lab == "up"], g[lab == "down"]) })
      n_pass <- sum(scr < 0 & scr < p5, na.rm = TRUE); KA[["scramble"]] <- data.table(test = "sign_scrambled_imprint_T2", value = n_pass / N_SCR, detail = sprintf("%d of %d scrambles pass T2 on %s (rho < 0 and < random P5 %.2f)", n_pass, N_SCR, ds, p5))
      chk("sign-scrambled imprint passes T2 in at most 5% of scrambles", TRUE, n_pass / N_SCR <= 0.05) }
  }
  fwrite(SC, file.path(TABWS, paste0("ws_s8_scores_", set, ".csv"))); rm(CELL); gc(FALSE)
}
## -- 4. known answers, gate, records ----------------------------------------------------------------------------------
TESTS <- rbindlist(TESTS); fwrite(TESTS, file.path(TABWS, "ws_s8_tests.csv"))
if (!is.null(RESULT[["disc primary"]])) { r <- RESULT[["disc primary"]]; ok <- isTRUE(r$T2) && (isTRUE(r$T3) || (is.na(r$T3) && isTRUE(r$T1)))
  KA[["positive"]] <- data.table(test = "discovery_cohort_positive_control", value = as.numeric(ok), detail = sprintf("GSE185381 adult AML: T2 %s, T3 %s, T1 %s", r$T2, r$T3, r$T1))
  chk("positive control: the discovery cohort passes the gate tests", TRUE, ok) }
if (length(KA)) fwrite(rbindlist(KA), file.path(TABWS, "ws_s8_known_answers.csv"))
gate <- if (!is.null(RESULT[["rep primary"]])) { r <- RESULT[["rep primary"]]; if (is.na(r$T3)) isTRUE(r$T2) && isTRUE(r$T1) else isTRUE(r$T2) && isTRUE(r$T3) } else NA
verdict <- if (isTRUE(gate)) "GATE_S8 PASS: the T imprint's structure replicates in GSE289435" else if (isFALSE(gate)) "GATE_S8 FAIL: the T imprint's structure does not replicate in GSE289435" else "GATE_S8 not evaluated (rep set missing)"
ts("    ", verdict)
if (!SMOKE) { fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))
  block <- c("<!-- A1-S8 start -->", sprintf("## T imprint, structural replication (S8, prereg-s8-v1) -- %s", format(Sys.Date())), "",
    "Scripts: scripts/ws/32_s8_replication.R. Gene sets: primary = 370 T-imprint genes without the 32 ambient-prone ones; full = 402. Sets: disc = GSE185381",
    "adult AML (positive control), rep = GSE289435 (decisive), rep2 = GSE227903 diagnosis, desc = GSE116256 (Seq-Well). Random-set nulls as in S5.", "",
    "| set | gene set | test | statistic | random P5 | random P95 | n | pass | detail |", "|---|---|---|---|---|---|---|---|---|",
    TESTS[test != "info_score_vs_ambient", sprintf("| %s | %s | %s | %s | %s | %s | %d | %s | %s |", set, gene_set, test, fmt(statistic), fmt(null_P5), fmt(null_P95), n, pass, detail)], "",
    paste0("Score vs ambient shares: ", paste(TESTS[test == "info_score_vs_ambient" & gene_set == "primary", sprintf("%s %s", set, detail)], collapse = "; ")),
    if (length(KA)) paste0("Known answers: ", paste(rbindlist(KA)$detail, collapse = "; ")) else "Known answers: not run",
    sprintf("**Verdict: %s**", verdict), "<!-- A1-S8 end -->")
  old <- readLines(GATES); s0 <- grep("<!-- A1-S8 start -->", old, fixed = TRUE); e0 <- grep("<!-- A1-S8 end -->", old, fixed = TRUE)
  writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES); ts("    GATES.md A1-S8 block written") }
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, if (SMOKE) "ws_s8_checks_smoke.csv" else "ws_s8_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
