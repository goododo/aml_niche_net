#!/usr/bin/env Rscript
# INPUT  : results/tables/ws/ws_s5_scores.csv, ws_s5_imprint_genes.csv, ws_s5_de_<compartment>.csv   (27, frozen list)
#          /LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds  (raw sums per compartment; per-cell assignment)
#          results/tables/ws/decontx/GSE185381/ws_decontx_index.csv.gz    (cell -> channel, hierarchy_bin: prop7, progenitor fraction)
#          results/tables/ws/ws_s3_gse185381_channels.csv                 (the three channels with adult AML + healthy)
#          results/tables/ws/decontx/GSE185381/ws_s1_ambient_blacklist.csv (extended ambient ranking, >= 2 channels)
#          /LARGE1/.../02_seurat_objects/04_annotated/GSE185381/<sample>.rds (only the samples of the mixed channels)
# OUTPUT : results/tables/ws/ws_s5_coordination.csv   pairwise and mean coordination, nulls, bootstrap CI (prereg 5)
#          results/tables/ws/ws_s5_same_channel.csv    per channel x compartment: sign concordance vs random sets (prereg 6)
#          results/tables/ws/ws_s5_gates.csv           GATE_A1 (1)-(4) and GATE_A2 with the frozen criteria (prereg 7)
#          results/tables/ws/ws_s5_checks2.csv; scripts/ws/GATES.md block "A1/A2"
# what it does: PREREGISTRATION_S5_imprint.md sections 5, 6, 7 and 8.4. Reads 27's imprint and scores as frozen
#   inputs; computes nothing that could change them. Coordination = mean pairwise Spearman of the residualised
#   patient scores across compartments (pairwise complete), against (1) independent per-compartment
#   permutations and (2) expression-matched random gene sets run through the same scoring and residualisation;
#   same-channel control on the three channels that hold an adult AML donor with a healthy donor.
#   usage: Rscript scripts/ws/28_s5_coordination.R [--smoke]   (--smoke: 100 permutations / random sets)
# This script deletes nothing and does not modify 27's outputs.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws"); GATES <- here::here("scripts", "ws", "GATES.md")
N_NULL <- if (SMOKE) 100L else 1000L
PB   <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", DATASET, paste0("ws_pb_", DATASET, ".rds")))
SC   <- fread(file.path(TABWS, "ws_s5_scores.csv")); IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))
BL   <- fread(file.path(TABWS, "decontx", DATASET, "ws_s1_ambient_blacklist.csv"))
IX   <- fread(file.path(TABWS, "decontx", DATASET, "ws_decontx_index.csv.gz"), select = c("cell", "sample", "unit", "hierarchy_bin", "bmm_fine", "high_error"))
CH   <- fread(file.path(TABWS, "ws_s3_gse185381_channels.csv"))
COMPS <- IMP[imprint_called == TRUE, unique(compartment)]
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-64s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
ts("[0] compartments with an imprint: ", paste(COMPS, collapse = ", "))
if (length(COMPS) < 2) { ts("fewer than two compartments carry an imprint: coordination is undefined; writing gates as not testable"); }

## -- 1. composition covariates per sample (prereg 5.2) -------------------------------------------------------------
IX[, bin := fifelse(bmm_fine %in% c("CLP", "MLP", "MLP-II", "Pre-ProB"), "EarlyLymphoid", hierarchy_bin)]
P7 <- dcast(IX[!is.na(bin) & bin != "" & !(high_error %in% 1L), .N, by = .(sample, bin)], sample ~ bin, value.var = "N", fill = 0)
bins <- setdiff(names(P7), "sample"); P7m <- as.matrix(P7[, ..bins]); P7m <- P7m / rowSums(P7m); rownames(P7m) <- P7$sample
prog_frac <- rowSums(P7m[, intersect(c("HSC_MPP", "LMPP_GMP", "EarlyLymphoid"), bins), drop = FALSE])
pcs <- prcomp(P7m[, setdiff(bins, c("Stromal", "unassigned"))], center = TRUE, scale. = FALSE)$x[, 1:2]
COV <- data.table(sample = rownames(P7m), prog_frac = prog_frac, comp_pc1 = pcs[, 1], comp_pc2 = pcs[, 2])
SC <- merge(SC, COV, by = "sample", all.x = TRUE)
SC <- merge(SC, unique(fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))[, .(sample, channels)]), by = "sample", all.x = TRUE)
chk("composition covariates for every scored sample", 0L, SC[is.na(prog_frac), .N])

## -- 2. the gene-filtered logCPM per compartment, for scoring random sets (same filter as 27) ---------------------------
G <- PB$genes
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene], "XIST", GO[chr == "chrY", gene],
                  grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))
LOGCPM <- list(); UNIV <- list()
for (g4 in COMPS) {
  d <- SC[compartment == g4]; y <- DGEList(PB$raw[[g4]][!(G %in% BLACK), d$sample], group = d$status)
  keep <- filterByExpr(y, group = d$status) & rowMeans(cpm(y) > 1) >= 0.5; y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  LOGCPM[[g4]] <- cpm(y, log = TRUE, prior.count = 1)
  U <- data.table(gene = rownames(LOGCPM[[g4]]), mean_lcpm = rowMeans(LOGCPM[[g4]]), det = rowMeans(PB$raw[[g4]][rownames(LOGCPM[[g4]]), d$sample] > 0))
  U[, bin := paste(cut(mean_lcpm, quantile(mean_lcpm, 0:10 / 10), include.lowest = TRUE), cut(det, c(-Inf, 0.5, 0.8, 0.95, Inf)))]
  UNIV[[g4]] <- U
  chk(paste(g4, "imprint genes are all in the scoring universe"), 0L, IMP[compartment == g4 & !gene %in% U$gene, .N])
}
score_set <- function(g4, up, dn) { Z <- t(scale(t(LOGCPM[[g4]]))); s <- if (length(up)) colMeans(Z[up, , drop = FALSE]) else 0
  s - (if (length(dn)) colMeans(Z[dn, , drop = FALSE]) else 0) }
matched_set <- function(g4, genes) { U <- UNIV[[g4]]; tb <- table(U[gene %in% genes, bin])
  unlist(lapply(names(tb), function(b) { pool <- U[bin == b & !gene %in% genes, gene]; if (!length(pool)) pool <- U[!gene %in% genes, gene]; sample(pool, tb[[b]], replace = length(pool) < tb[[b]]) })) }
chk("27's score reproduces from the recomputed logCPM (T, max diff < 1e-6)", TRUE, if ("T" %in% COMPS) {
  s <- score_set("T", IMP[compartment == "T" & direction == "up", gene], IMP[compartment == "T" & direction == "down", gene])
  max(abs(s[SC[compartment == "T", sample]] - SC[compartment == "T", score])) < 1e-6 } else NA)

## -- 3. coordination (prereg 5) ----------------------------------------------------------------------------------------
resid_scores <- function(scores_by_comp) {   # list comp -> named numeric over AML patients; residualise on the composition covariates
  lapply(scores_by_comp, function(s) { d <- COV[match(names(s), sample)]; x <- SC[match(names(s), sample), .(hb, gran)]   # hb/gran are sample-level
    r <- resid(lm(s ~ d$prog_frac + d$comp_pc1 + d$comp_pc2 + x$hb + x$gran)); names(r) <- names(s); r }) }
coord <- function(R) { cs <- names(R); pr <- combn(cs, 2); v <- apply(pr, 2, function(p) { a <- R[[p[1]]]; b <- R[[p[2]]]; k <- intersect(names(a), names(b)); if (length(k) < 8) NA_real_ else cor(a[k], b[k], method = "spearman") })
  list(pairs = data.table(a = pr[1, ], b = pr[2, ], rho = v, n = apply(pr, 2, function(p) length(intersect(names(R[[p[1]]]), names(R[[p[2]]]))))), mean = mean(v, na.rm = TRUE)) }
AML <- SC[status == "AML"]
real_scores <- lapply(setNames(COMPS, COMPS), function(g4) { x <- AML[compartment == g4]; setNames(x$score, x$sample) })
if (length(COMPS) >= 2) {
  RS <- resid_scores(real_scores); obs <- coord(RS)
  ts(sprintf("[3] coordination: mean pairwise Spearman %.3f over %d pairs", obs$mean, nrow(obs$pairs))); print(obs$pairs)
  null1 <- vapply(seq_len(N_NULL), function(i) coord(lapply(RS, function(r) setNames(sample(r), names(r))))$mean, 1)
  null2 <- vapply(seq_len(N_NULL), function(i) { rs <- lapply(setNames(COMPS, COMPS), function(g4) { up <- IMP[compartment == g4 & direction == "up", gene]; dn <- IMP[compartment == g4 & direction == "down", gene]
      fu <- matched_set(g4, c(up, dn)); s <- score_set(g4, head(fu, length(up)), tail(fu, length(dn))); s[AML[compartment == g4, sample]] }); coord(resid_scores(rs))$mean }, 1)
  pats <- unique(AML$sample)
  boot <- vapply(seq_len(N_NULL), function(i) { bs <- sample(pats, replace = TRUE); rs <- lapply(RS, function(r) { k <- bs[bs %in% names(r)]; v <- r[k]; names(v) <- make.unique(k); v })
    coord(rs)$mean }, 1)
  # a bootstrap draw (or a null draw) can leave a pair with fewer than 8 shared patients, which coord() scores as NA
  ok1 <- is.finite(null1); ok2 <- is.finite(null2); okb <- is.finite(boot)
  ts(sprintf("    valid draws: permutation %d, random sets %d, bootstrap %d of %d", sum(ok1), sum(ok2), sum(okb), N_NULL))
  p1 <- (1 + sum(null1[ok1] >= obs$mean)) / (sum(ok1) + 1); p2_q95 <- quantile(null2, 0.95, na.rm = TRUE); ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  COORD <- data.table(statistic = "mean_pairwise_spearman", observed = obs$mean, n_pairs = nrow(obs$pairs), perm_p = p1, perm_P95 = quantile(null1, .95, na.rm = TRUE),
                      random_set_P95 = p2_q95, random_set_median = median(null2, na.rm = TRUE), boot_ci_lo = ci[1], boot_ci_hi = ci[2], n_null = N_NULL,
                      n_valid_perm = sum(ok1), n_valid_random = sum(ok2), n_valid_boot = sum(okb))
  # sensitivity: the 11 patients with T, NK and B all eligible (prereg 5.3)
  core <- Reduce(intersect, lapply(real_scores[intersect(COMPS, c("T", "NK", "B"))], names))
  if (length(core) >= 8 && length(intersect(COMPS, c("T", "NK", "B"))) >= 2) {
    rc <- coord(lapply(RS[intersect(COMPS, c("T", "NK", "B"))], function(r) r[core]))
    COORD <- rbind(COORD, data.table(statistic = sprintf("T_NK_B_complete_n%d", length(core)), observed = rc$mean, n_pairs = nrow(rc$pairs), perm_p = NA, perm_P95 = NA, random_set_P95 = NA, random_set_median = NA, boot_ci_lo = NA, boot_ci_hi = NA, n_null = NA), fill = TRUE) }
  fwrite(cbind(COORD, obs$pairs[, .(pairs = paste(sprintf("%s-%s %.2f (n%d)", a, b, rho, n), collapse = "; "))][1]), file.path(TABWS, "ws_s5_coordination.csv"))
  ts(sprintf("    perm p %.3f (P95 %.3f) | random-set P95 %.3f (median %.3f) | bootstrap 95%% CI %.3f to %.3f", p1, quantile(null1, .95, na.rm = TRUE), p2_q95, median(null2, na.rm = TRUE), ci[1], ci[2]))
}

## -- 4. same-channel control (prereg 6) -------------------------------------------------------------------------------------
ts("[4] same-channel control")
mixed <- CH[s6_adult_AML_with_healthy == TRUE, unit]; chk("three mixed channels", 3L, length(mixed))
S6 <- list(); AN <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
cells_mixed <- IX[unit %in% mixed]; kept <- PB$cells[kept == TRUE, .(cell, group4)]; cells_mixed <- merge(cells_mixed, kept, by = "cell")
roles <- unique(SC[, .(sample, status)])
objs <- unique(cells_mixed$sample); objs <- objs[objs %in% roles$sample]
SUMS <- list()
for (sm in objs) { o <- readRDS(file.path(AN, paste0(sm, ".rds"))); X <- LayerData(o, assay = "RNA", layer = "counts")
  for (u in unique(cells_mixed[sample == sm, unit])) for (g4 in COMPS) { cc <- cells_mixed[sample == sm & unit == u & group4 == g4, cell]
    if (length(cc) >= WS_S4_MIN_CELLS) SUMS[[paste(sm, u, g4)]] <- list(sample = sm, unit = u, comp = g4, n = length(cc), sum = as.numeric(Matrix::rowSums(X[, cc, drop = FALSE]))) }
  rm(o, X); gc(FALSE) }
for (u in mixed) for (g4 in COMPS) {
  S <- Filter(function(z) z$unit == u && z$comp == g4, SUMS); if (!length(S)) next
  st <- roles[match(vapply(S, `[[`, "", "sample"), sample), status]; if (!all(c("AML", "healthy") %in% st)) next
  M <- sapply(S, `[[`, "sum"); rownames(M) <- G; L <- log2(t(t(M) / colSums(M)) * 1e6 + 1); lfc <- rowMeans(L[, st == "AML", drop = FALSE]) - rowMeans(L[, st == "healthy", drop = FALSE])
  de <- fread(file.path(TABWS, paste0("ws_s5_de_", g4, ".csv"))); imp <- IMP[compartment == g4, gene]
  conc <- mean(sign(lfc[imp]) == sign(de[match(imp, gene), logFC]))
  null <- vapply(seq_len(N_NULL), function(i) { fu <- matched_set(g4, imp); mean(sign(lfc[fu]) == sign(de[match(fu, gene), logFC])) }, 1)
  S6[[paste(u, g4)]] <- data.table(unit = u, compartment = g4, n_aml = sum(st == "AML"), n_healthy = sum(st == "healthy"), concordance = conc, random_P95 = quantile(null, .95), random_median = median(null), above_P95 = conc > quantile(null, .95))
}
S6 <- rbindlist(S6); fwrite(S6, file.path(TABWS, "ws_s5_same_channel.csv")); print(S6)
# (4) healthy donors in mixed channels vs the others
HM <- SC[status == "healthy", .(mixed = any(strsplit(channels, ";")[[1]] %in% mixed)), by = .(sample, compartment, score)]
shift <- HM[, .(d_mixed = mean(score[mixed]) - mean(score[!mixed]), n_mixed = sum(mixed)), by = compartment]
gap <- SC[, .(gap = mean(score[status == "AML"]) - mean(score[status == "healthy"])), by = compartment]
shift <- merge(shift, gap, by = "compartment")[compartment %in% COMPS]; shift[, ratio := d_mixed / gap]   # only compartments with an imprint

## -- 5. gates (prereg 7) ----------------------------------------------------------------------------------------------------
ts("[5] gates")
pres <- function(g4, cov) { d <- SC[compartment == g4]; r1 <- resid(lm(score ~ status, d)); r2 <- resid(lm(get(cov) ~ status, d)); cor(r1, r2, method = "spearman") }
A1_1 <- rbindlist(lapply(COMPS, function(g4) data.table(compartment = g4, rho_hb = pres(g4, "hb"), rho_gran = pres(g4, "gran"))))
ext <- BL[n_units >= 2, gene]
A1_2 <- rbindlist(lapply(COMPS, function(g4) { U <- UNIV[[g4]]$gene; imp <- IMP[compartment == g4, gene]; k <- sum(imp %in% ext)
  data.table(compartment = g4, imprint_in_extended_ambient = k, p_hyper = phyper(k - 1, sum(U %in% ext), sum(!U %in% ext), length(imp), lower.tail = FALSE)) }))
GT <- rbind(
  data.table(gate = "GATE_A1_1", criterion = "|partial rho(score, hb)| and |rho(score, gran)| < 0.30 in every compartment", value = paste(A1_1[, sprintf("%s hb %.2f gran %.2f", compartment, rho_hb, rho_gran)], collapse = "; "), pass = all(abs(A1_1$rho_hb) < 0.3 & abs(A1_1$rho_gran) < 0.3)),
  data.table(gate = "GATE_A1_2", criterion = "imprint not enriched in the extended ambient ranking (hypergeometric p > 0.05)", value = paste(A1_2[, sprintf("%s %d genes p %.3f", compartment, imprint_in_extended_ambient, p_hyper)], collapse = "; "), pass = all(A1_2$p_hyper > 0.05)),
  data.table(gate = "GATE_A1_3", criterion = "same-channel concordance above the random-set P95 in >= 2 of 3 channels (pooled over compartments)", value = if (nrow(S6)) paste(S6[, .(above = sum(above_P95), tested = .N), by = unit][, sprintf("%s %d/%d", unit, above, tested)], collapse = "; ") else "not testable", pass = if (nrow(S6)) S6[, .(ok = mean(above_P95) >= 0.5), by = unit][, sum(ok)] >= 2 else NA),
  data.table(gate = "GATE_A1_4", criterion = "healthy donors in mixed channels shift < 1/3 of the AML-healthy gap", value = paste(shift[, sprintf("%s %.2f", compartment, ratio)], collapse = "; "), pass = all(abs(shift$ratio) < 1/3)),
  data.table(gate = "GATE_A2", criterion = "coordination > 0, bootstrap 95% CI excludes 0, above the random-set P95", value = if (exists("COORD")) sprintf("%.3f, CI %.3f to %.3f, random P95 %.3f, perm p %.3f", COORD$observed[1], COORD$boot_ci_lo[1], COORD$boot_ci_hi[1], COORD$random_set_P95[1], COORD$perm_p[1]) else "fewer than two imprints",
             pass = if (exists("COORD")) COORD$observed[1] > 0 && COORD$boot_ci_lo[1] > 0 && COORD$observed[1] > COORD$random_set_P95[1] else NA))
fwrite(GT, file.path(TABWS, "ws_s5_gates.csv")); print(GT[, .(gate, pass, value)])
if (!SMOKE) {
  block <- c("<!-- A1A2 start -->", sprintf("## GATE_A1 / GATE_A2 -- %s (prereg-s5-v1)", format(Sys.Date())), "",
    "Scripts: scripts/ws/27_s5_imprint.R, 28_s5_coordination.R | compartments with an imprint: ", paste(COMPS, collapse = ", "), "",
    "| gate | criterion | value | pass |", "|---|---|---|---|", GT[, sprintf("| %s | %s | %s | %s |", gate, criterion, value, pass)], "",
    sprintf("**Verdict: GATE_A1 %s, GATE_A2 %s** -- any FAIL stops the line for discussion (v2.0 5.4).", if (all(GT[grepl("A1", gate), pass], na.rm = TRUE)) "PASS" else "FAIL", if (isTRUE(GT[gate == "GATE_A2", pass])) "PASS" else "FAIL"), "<!-- A1A2 end -->")
  old <- readLines(GATES); s0 <- grep("<!-- A1A2 start -->", old); e0 <- grep("<!-- A1A2 end -->", old)
  writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES); ts("    GATES.md A1/A2 block written")
}
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, "ws_s5_checks2.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
