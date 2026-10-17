#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../aml_niche_net/ws_dialogue/{rA_<run>.rds, R_<run>.rds, planted_truth.rds}   (29)
#          results/tables/ws/ws_s4_pseudobulk_meta_GSE185381.csv, decontx/GSE185381/ws_decontx_index.csv.gz (covariates)
#          scripts/ws/PREREGISTRATION_S7_dialogue.md (tag prereg-s7-v1)
# OUTPUT : results/tables/ws/ws_s7_mcps.csv           MCP x compartment: n genes up/down, empirical p per pair
#          results/tables/ws/ws_s7_scores.csv         sample x compartment x MCP: mean cell score
#          results/tables/ws/ws_s7_criteria.csv       criteria 1-3 per MCP and run, and the verdict
#          results/tables/ws/ws_s7_known_answers.csv  planted MCP recovery; shuffled-sample false-positive rate
#          results/tables/ws/ws_s7_checks.csv; scripts/ws/GATES.md block "A1b-S7"
#          /LARGE1/.../ws_dialogue/{stability,shuffle}/   DIALOGUE outputs of the resampling runs
# what it does: PREREGISTRATION_S7_dialogue.md sections 3 and 4. Criterion 1 from DIALOGUE's empirical p (the
#   largest p over the pairs of compartments the MCP spans); criterion 2 from the sample-level MCP scores
#   (mean per sample of DIALOGUE's per-cell scores) residualised on prog_frac, two composition PCs and the two
#   ambient shares; criterion 3 by rerunning DIALOGUE on 20 random 70% patient subsets, matching MCPs to the
#   full run by sample-score correlation, and the Jaccard of the top-50 genes per compartment. Known answers:
#   the planted run (29) and 20 shuffled-sample runs. Nothing here changes an MCP.
#   usage: Rscript scripts/ws/30_s7_dialogue_evaluate.R [--runs=R1,R2] [--n-resample=20] [--skip-resample]
# This script deletes nothing and does not modify 29's outputs.

suppressPackageStartupMessages({ library(data.table); library(DIALOGUE) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
RUNS <- strsplit(.get("runs", "R1,R2"), ",")[[1]]; N_RES <- as.integer(.get("n-resample", "20")); SKIP_RES <- "--skip-resample" %in% .a
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
options(error = function() { calls <- sys.calls(); n <- length(calls); message("--- TRACEBACK ---")
  for (i in seq(max(1, n - 15), n)) message(i, ": ", substr(paste(deparse(calls[[i]]), collapse = " "), 1, 220)); q(status = 1) })
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws"); GATES <- here::here("scripts", "ws", "GATES.md")
OUTD <- file.path(LARGE1_DIR, "ws_dialogue"); K <- 3L; ABN <- 15L
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-64s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 0. covariates per sample (as 28) -----------------------------------------------------------------------------
META <- fread(file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))
SMP <- unique(META[, .(sample, hb = hb_share_nonery, gran = gran_share_lymph)])
IX <- fread(file.path(TABWS, "decontx", DATASET, "ws_decontx_index.csv.gz"), select = c("sample", "hierarchy_bin", "bmm_fine", "high_error"))
IX[, bin := fifelse(bmm_fine %in% c("CLP", "MLP", "MLP-II", "Pre-ProB"), "EarlyLymphoid", hierarchy_bin)]
P7 <- dcast(IX[!is.na(bin) & bin != "" & !(high_error %in% 1L), .N, by = .(sample, bin)], sample ~ bin, value.var = "N", fill = 0)
bins <- setdiff(names(P7), "sample"); M <- as.matrix(P7[, ..bins]); M <- M / rowSums(M); rownames(M) <- P7$sample
COV <- data.table(sample = rownames(M), prog_frac = rowSums(M[, intersect(c("HSC_MPP", "LMPP_GMP", "EarlyLymphoid"), bins), drop = FALSE]),
                  pc = prcomp(M[, setdiff(bins, c("Stromal", "unassigned"))])$x[, 1:2])
setnames(COV, c("sample", "prog_frac", "comp_pc1", "comp_pc2")); COV <- merge(COV, SMP, by = "sample")

## -- helpers -------------------------------------------------------------------------------------------------------
run_dialogue <- function(rA, main, covar, rd) {
  dir.create(rd, showWarnings = FALSE, recursive = TRUE)
  tryCatch(DIALOGUE.run(rA = rA, main = main, param = DLG.get.param(k = K, results.dir = paste0(rd, "/"), conf = "cellQ", covar = covar, abn.c = ABN, seed1 = SEED, plot.flag = FALSE)),
           error = function(e) { f <- file.path(rd, paste0("DLG.full.output_", main, ".rds")); if (file.exists(f)) readRDS(f) else NULL })
}
mcp_table <- function(R) {   # one row per MCP: compartments it spans, its p (largest over its pairs), gene counts
  if (is.null(R) || is.null(R$emp.p)) return(NULL)
  E <- as.matrix(R$emp.p); rbindlist(lapply(rownames(E), function(m) {
    cts <- R$MCP.cell.types[[m]]; if (is.null(cts) || length(cts) < 2) return(data.table(mcp = m, compartments = "", n_comp = 0L, p_max = NA_real_, p_min = NA_real_))
    pairs <- combn(cts, 2); pn <- apply(pairs, 2, function(p) { a <- paste(p, collapse = "_"); b <- paste(rev(p), collapse = "_"); if (a %in% colnames(E)) E[m, a] else if (b %in% colnames(E)) E[m, b] else NA_real_ })
    data.table(mcp = m, compartments = paste(cts, collapse = "+"), n_comp = length(cts), p_max = max(pn, na.rm = TRUE), p_min = min(pn, na.rm = TRUE)) })) }
sample_scores <- function(R) {   # sample x compartment x MCP: mean of the per-cell MCP scores
  rbindlist(lapply(names(R$scores), function(ct) { S <- as.data.table(R$scores[[ct]]); mc <- grep("^MCP", names(S), value = TRUE)
    melt(S[, c(lapply(.SD, mean), .(n = .N)), by = samples, .SDcols = mc], id.vars = c("samples", "n"), variable.name = "mcp", value.name = "score")[, compartment := ct] })) }
genes_of <- function(R, m, ct, n = 50) { g <- R$MCPs.full[[m]]; if (is.null(g)) return(character(0))
  u <- g[[paste0(ct, ".up")]]; d <- g[[paste0(ct, ".down")]]; c(head(u, n), head(d, n)) }
coord_resid <- function(SS, m, cts) {   # criterion 2
  W <- dcast(SS[mcp == m & compartment %in% cts & n >= ABN], samples ~ compartment, value.var = "score")
  W <- merge(W, COV, by.x = "samples", by.y = "sample"); if (nrow(W) < 8) return(list(rho = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(W)))
  res <- sapply(cts, function(ct) resid(lm(W[[ct]] ~ prog_frac + comp_pc1 + comp_pc2 + hb + gran, data = W)))
  pr <- combn(cts, 2); f <- function(idx) mean(apply(pr, 2, function(p) cor(res[idx, p[1]], res[idx, p[2]], method = "spearman")))
  obs <- f(seq_len(nrow(W))); bt <- replicate(1000, f(sample(nrow(W), replace = TRUE)))
  list(rho = obs, lo = unname(quantile(bt, .025, na.rm = TRUE)), hi = unname(quantile(bt, .975, na.rm = TRUE)), n = nrow(W)) }
jacc <- function(a, b) if (!length(a) || !length(b)) NA_real_ else length(intersect(a, b)) / length(union(a, b))

## -- 1. criteria 1 and 2 per run ----------------------------------------------------------------------------------
MCPS <- list(); SCORES <- list(); CRIT <- list(); KA <- list()
for (run in RUNS) {
  R <- readRDS(file.path(OUTD, paste0("R_", run, ".rds"))); rA <- readRDS(file.path(OUTD, paste0("rA_", run, ".rds")))
  covar <- R$covar; ts(sprintf("[%s] %d samples in all compartments | covar %s", run, length(R$samples), paste(covar, collapse = "+")))
  MT <- mcp_table(R); if (is.null(MT)) { ts("    no MCP object"); next }
  SS <- sample_scores(R); SCORES[[run]] <- cbind(run = run, SS)
  MT[, run := run]; for (ct in R$cell.types) MT[, (paste0("n_genes_", ct)) := vapply(mcp, function(m) length(genes_of(R, m, ct, 1e6)), 1L)]
  MCPS[[run]] <- MT
  for (i in seq_len(nrow(MT))) { m <- MT$mcp[i]; cts <- strsplit(MT$compartments[i], "+", fixed = TRUE)[[1]]
    c2 <- if (length(cts) >= 2) coord_resid(SS, m, cts) else list(rho = NA, lo = NA, hi = NA, n = NA)
    CRIT[[paste(run, m)]] <- data.table(run = run, mcp = m, compartments = MT$compartments[i], crit1_p = MT$p_max[i], crit1 = isTRUE(MT$p_max[i] <= 0.05),
                                        crit2_rho = c2$rho, crit2_lo = c2$lo, crit2_hi = c2$hi, crit2_n = c2$n, crit2 = isTRUE(c2$rho > 0.3 && c2$lo > 0),
                                        crit3_jaccard_median = NA_real_, crit3_comp_pass = NA_integer_, crit3 = NA) }
  print(MT[, .(mcp, compartments, p_max = round(p_max, 3), p_min = round(p_min, 3))])
  ## -- 2. criterion 3: stability under 20 x 70% patient subsets ------------------------------------------------------
  if (!SKIP_RES) {
    pats <- R$samples; JC <- list()
    for (b in seq_len(N_RES)) {
      keep <- sample(pats, ceiling(0.7 * length(pats))); rAs <- lapply(rA, function(r) { j <- r@samples %in% keep
        make.cell.type(name = r@name, tpm = r@tpm[, j], samples = r@samples[j], X = r@X[j, ], metadata = r@metadata[j, ], cellQ = r@cellQ[j]) })
      Rb <- run_dialogue(rAs, paste0(run, "_stab", b), covar, file.path(OUTD, "stability", paste0(run, "_stab", b)))
      if (is.null(Rb) || is.null(Rb$emp.p)) { ts(sprintf("    stab %02d: no MCP object", b)); next }
      SSb <- sample_scores(Rb)
      for (m in rownames(as.matrix(R$emp.p))) { cts <- R$MCP.cell.types[[m]]; if (is.null(cts)) next
        # match: the subset MCP whose sample scores correlate best with the full MCP across the shared samples and compartments
        best <- NA_character_; bc <- -Inf
        for (mb in rownames(as.matrix(Rb$emp.p))) { v <- sapply(cts, function(ct) { a <- SS[mcp == m & compartment == ct]; bb <- SSb[mcp == mb & compartment == ct]
            k <- intersect(a$samples, bb$samples); if (length(k) < 5) NA_real_ else cor(a[match(k, samples), score], bb[match(k, samples), score], method = "spearman") })
          cc <- mean(abs(v), na.rm = TRUE); if (is.finite(cc) && cc > bc) { bc <- cc; best <- mb } }
        if (is.na(best)) next
        for (ct in cts) JC[[length(JC) + 1]] <- data.table(run = run, mcp = m, draw = b, matched = best, match_cor = bc, compartment = ct,
                                                           jaccard = jacc(genes_of(R, m, ct), genes_of(Rb, best, ct))) }
      ts(sprintf("    stab %02d done", b))
    }
    JC <- rbindlist(JC)
    if (nrow(JC)) { fwrite(JC, file.path(TABWS, paste0("ws_s7_stability_", run, ".csv")))
      for (m in unique(JC$mcp)) { j <- JC[mcp == m, .(med = median(jaccard, na.rm = TRUE)), by = compartment]
        CRIT[[paste(run, m)]][, `:=`(crit3_jaccard_median = median(JC[mcp == m, jaccard], na.rm = TRUE), crit3_comp_pass = j[med >= 0.5, .N], crit3 = j[med >= 0.5, .N] >= 2)] } }
  }
}
MCPS <- rbindlist(MCPS, fill = TRUE); SCORES <- rbindlist(SCORES); CRIT <- rbindlist(CRIT)
CRIT[, found := crit1 & crit2 & (crit3 %in% TRUE)]
fwrite(MCPS, file.path(TABWS, "ws_s7_mcps.csv")); fwrite(SCORES, file.path(TABWS, "ws_s7_scores.csv"))

## -- 3. known answers: planted MCP and shuffled samples ------------------------------------------------------------
ts("[3] known answers")
pf <- file.path(OUTD, "R_planted.rds")
if (file.exists(pf)) {
  Rp <- readRDS(pf); truth <- readRDS(file.path(OUTD, "planted_truth.rds")); MTp <- mcp_table(Rp)
  rec <- rbindlist(lapply(MTp$mcp, function(m) { cts <- Rp$MCP.cell.types[[m]]; if (is.null(cts)) return(NULL)
    rbindlist(lapply(cts, function(ct) { g <- genes_of(Rp, m, ct, 200); data.table(mcp = m, compartment = ct, p = MTp[mcp == m, p_max], n_genes = length(g), planted_found = length(intersect(g, truth$genes[[ct]])), planted_frac = length(intersect(g, truth$genes[[ct]])) / length(truth$genes[[ct]])) })) }))
  if (nrow(rec)) { best <- rec[, .(p = p[1], mean_frac = mean(planted_frac), n_comp = .N), by = mcp][order(-mean_frac)][1]
    KA[["planted"]] <- data.table(test = "planted_MCP", value = best$mean_frac, detail = sprintf("best MCP %s: p %.3f, planted genes among its top 200 per compartment %.2f (over %d compartments)", best$mcp, best$p, best$mean_frac, best$n_comp))
    chk("planted MCP recovered (p <= 0.05 and >= 50% of planted genes)", TRUE, isTRUE(best$p <= 0.05 && best$mean_frac >= 0.5)) }
  print(rec)
}
if (!SKIP_RES && "R1" %in% RUNS) {
  rA <- readRDS(file.path(OUTD, "rA_R1.rds")); R <- readRDS(file.path(OUTD, "R_R1.rds")); fp <- 0L; nd <- 0L
  for (b in seq_len(N_RES)) {
    rAs <- lapply(rA, function(r) { s <- r@samples; u <- unique(s); map <- setNames(sample(u), u)     # relabel samples within the compartment: pairing destroyed
      make.cell.type(name = r@name, tpm = r@tpm, samples = unname(map[s]), X = r@X, metadata = r@metadata, cellQ = r@cellQ) })
    Rb <- run_dialogue(rAs, paste0("R1_shuf", b), R$covar, file.path(OUTD, "shuffle", paste0("R1_shuf", b)))
    if (is.null(Rb) || is.null(Rb$emp.p)) next
    nd <- nd + 1L; mt <- mcp_table(Rb); fp <- fp + as.integer(any(mt$p_max <= 0.05, na.rm = TRUE)); ts(sprintf("    shuffle %02d: min p over MCPs %.3f", b, min(mt$p_max, na.rm = TRUE)))
  }
  KA[["shuffle"]] <- data.table(test = "shuffled_samples", value = if (nd) fp / nd else NA_real_, detail = sprintf("%d of %d shuffles produced an MCP with p <= 0.05", fp, nd))
  chk("shuffled samples: false-positive rate of criterion 1 at or below 0.05", TRUE, isTRUE(nd > 0 && fp / nd <= 0.05))
}
if (length(KA)) fwrite(rbindlist(KA), file.path(TABWS, "ws_s7_known_answers.csv"))

## -- 4. verdict ------------------------------------------------------------------------------------------------------
fwrite(CRIT, file.path(TABWS, "ws_s7_criteria.csv")); print(CRIT[, .(run, mcp, compartments, crit1_p = round(crit1_p, 3), crit1, crit2_rho = round(crit2_rho, 2), crit2_n, crit2, crit3_jaccard_median = round(crit3_jaccard_median, 2), crit3, found)])
verdict <- if (any(CRIT$found)) "MCP FOUND (candidate only; needs its own confirmatory pre-registration)" else "NO MCP meets the three criteria: H-A1b CLOSED"
ts("    ", verdict)
if (!SKIP_RES) {
  block <- c("<!-- A1b-S7 start -->", sprintf("## H-A1b, unsupervised test (DIALOGUE, prereg-s7-v1) -- %s", format(Sys.Date())), "",
    "Scripts: scripts/ws/29_s7_dialogue_build_run.R, 30_s7_dialogue_evaluate.R", "",
    "| run | MCP | compartments | crit 1 p | crit 2 rho (n) | crit 3 Jaccard | found |", "|---|---|---|---|---|---|---|",
    CRIT[, sprintf("| %s | %s | %s | %.3f | %.2f (%s) | %s | %s |", run, mcp, compartments, crit1_p, crit2_rho, crit2_n, ifelse(is.na(crit3_jaccard_median), "NA", sprintf("%.2f", crit3_jaccard_median)), found)], "",
    if (length(KA)) paste0("Known answers: ", paste(rbindlist(KA)$detail, collapse = "; ")) else "Known answers: not run",
    sprintf("**Verdict: %s**", verdict), "<!-- A1b-S7 end -->")
  old <- readLines(GATES); s0 <- grep("<!-- A1b-S7 start -->", old); e0 <- grep("<!-- A1b-S7 end -->", old)
  writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES); ts("    GATES.md A1b-S7 block written")
}
C <- rbindlist(CHK); if (nrow(C)) fwrite(C, file.path(TABWS, "ws_s7_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
