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
# what it does: PREREGISTRATION_S7_dialogue.md sections 3, 4 and 8. Criterion 1 from the exact permutation p (rank of the
#   observed canonical correlation among DIALOGUE's 1 + 99 permutations, /100, recovered by lookup from its rank-sum emp.p;
#   the largest over the pairs the MCP spans; DIALOGUE's own p is reported alongside); criterion 2 from the sample-level MCP scores
#   (mean per sample of DIALOGUE's per-cell scores) residualised on prog_frac, two composition PCs and the two
#   ambient shares; criterion 3 by rerunning DIALOGUE on 20 random 70% patient subsets, matching MCPs to the
#   full run by sample-score correlation, and the Jaccard of the top-50 genes per compartment. Known answers:
#   the planted run (29) and 20 shuffled-sample runs. Nothing here changes an MCP.
#   usage: Rscript scripts/ws/30_s7_dialogue_evaluate.R [--runs=R1,R2] [--n-resample=20] [--skip-resample] [--out-dir=DIR]
#          --out-dir (smoke test): tables and the resampling DIALOGUE outputs go under DIR; GATES.md is not touched
# This script deletes nothing and does not modify 29's outputs.

suppressPackageStartupMessages({ library(data.table); library(DIALOGUE) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
RUNS <- strsplit(.get("runs", "R1,R2"), ",")[[1]]; N_RES <- as.integer(.get("n-resample", "20")); SKIP_RES <- "--skip-resample" %in% .a; OUT_DIR <- .get("out-dir", "")
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
ws_patch_dialogue()   # DIALOGUE 1.0's MCP cell-type assignment has a row-order bug (config_ws.R section 13; prereg section 8, item 5); must precede every DIALOGUE.run
set.seed(SEED)
options(error = function() { calls <- sys.calls(); n <- length(calls); message("--- TRACEBACK ---")
  for (i in seq(max(1, n - 15), n)) message(i, ": ", substr(paste(deparse(calls[[i]]), collapse = " "), 1, 220)); q(status = 1) })
DATASET <- WS_S4_DEFINE; TABWS <- file.path(TAB_DIR, "ws"); GATES <- here::here("scripts", "ws", "GATES.md")
OUTWS <- if (nzchar(OUT_DIR)) OUT_DIR else TABWS; dir.create(OUTWS, showWarnings = FALSE, recursive = TRUE); WRITE_GATES <- !SKIP_RES && !nzchar(OUT_DIR)
OUTD <- file.path(LARGE1_DIR, "ws_dialogue"); RESD <- if (nzchar(OUT_DIR)) file.path(OUT_DIR, "dialogue") else OUTD; K <- 3L; ABN <- 15L
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
P7 <- data.table::dcast(IX[!is.na(bin) & bin != "" & !(high_error %in% 1L), .N, by = .(sample, bin)], sample ~ bin, value.var = "N", fill = 0)
bins <- setdiff(names(P7), "sample"); M <- as.matrix(P7[, ..bins]); M <- M / rowSums(M); rownames(M) <- P7$sample
COV <- data.table(sample = rownames(M), prog_frac = rowSums(M[, intersect(c("HSC_MPP", "LMPP_GMP", "EarlyLymphoid"), bins), drop = FALSE]),
                  pc = prcomp(M[, setdiff(bins, c("Stromal", "unassigned"))])$x[, 1:2])
setnames(COV, c("sample", "prog_frac", "comp_pc1", "comp_pc2")); COV <- merge(COV, SMP, by = "sample")

## -- helpers -------------------------------------------------------------------------------------------------------
run_dialogue <- function(rA, main, covar, rd) {   # warnings (lme4 rank deficiency etc.) are printed once per distinct message and DIALOGUE call site, in-line with DIALOGUE's prints
  dir.create(rd, showWarnings = FALSE, recursive = TRUE); seen <- character(0); nsf <- ls(asNamespace("DIALOGUE"))
  withCallingHandlers(
    tryCatch(DIALOGUE.run(rA = rA, main = main, plot.flag = FALSE, param = DLG.get.param(k = K, results.dir = paste0(rd, "/"), conf = "cellQ", covar = covar, abn.c = ABN, seed1 = SEED, plot.flag = FALSE)),   # plot.flag must go to DIALOGUE.run itself
             error = function(e) { f <- file.path(rd, paste0("DLG.full.output_", main, ".rds")); if (file.exists(f)) readRDS(f) else NULL }),
    warning = function(w) { fn <- unique(vapply(sys.calls(), function(cl) sub(".*:::?", "", deparse(cl[[1]])[1]), "")); fn <- fn[fn %in% nsf]
      key <- paste(conditionMessage(w), paste(tail(fn, 3), collapse = "<"))
      if (!key %in% seen) { seen <<- c(seen, key); cat(sprintf("[warn %s] %s | in %s\n", main, conditionMessage(w), paste(tail(fn, 3), collapse = " < "))) }
      invokeRestart("muffleWarning") })
}
## DIALOGUE's emp.p is a one-sided rank-sum p (normal approximation) of the observed canonical correlation against its 99 permutations
## (DIALOGUE1: n1 = 100), i.e. a deterministic function of the observed value's rank; the exact permutation p = rank / 100 is recovered by lookup (prereg section 8)
## (lookup tolerance 1e-3: ties among permuted correlations shift the rank-sum p by ~1e-5, adjacent ranks differ by >= 3e-3)
LUT <- { rt <- get("ranksum.test", asNamespace("DIALOGUE")); data.table(rank = 1:100, dlg_p = vapply(1:100, function(r) as.numeric(rt(100 - r + 0.5, 1:99, alternative = "greater")), 1)) }
exact_p <- function(p) vapply(as.numeric(p), function(x) { if (is.na(x)) return(NA_real_); i <- which.min(abs(LUT$dlg_p - x)); if (abs(LUT$dlg_p[i] - x) > 1e-3) NA_real_ else LUT$rank[i] / 100 }, 1)
draw_distinct <- function(n, f) { out <- list(); while (length(out) < n) { x <- f(); if (!any(vapply(out, identical, TRUE, x))) out[[length(out) + 1]] <- x }; out }   # DIALOGUE1.PMD calls set.seed(seed1): draw every subset and permutation BEFORE the first DIALOGUE call
mcp_table <- function(R) {   # one row per MCP: compartments it spans, its p (largest over its pairs: DIALOGUE's rank-sum p and the exact permutation p), gene counts
  if (is.null(R) || is.null(R$emp.p)) return(NULL)
  E <- as.matrix(R$emp.p); rbindlist(lapply(rownames(E), function(m) {
    cts <- R$MCP.cell.types[[m]]; if (is.null(cts) || length(cts) < 2) return(data.table(mcp = m, compartments = "", n_comp = 0L, p_max = NA_real_, p_min = NA_real_, p_max_exact = NA_real_, p_min_exact = NA_real_))
    pairs <- combn(cts, 2); pn <- apply(pairs, 2, function(p) { a <- paste(p, collapse = "_"); b <- paste(rev(p), collapse = "_"); if (a %in% colnames(E)) E[m, a] else if (b %in% colnames(E)) E[m, b] else NA_real_ }); pe <- exact_p(pn)
    data.table(mcp = m, compartments = paste(cts, collapse = "+"), n_comp = length(cts), p_max = max(pn, na.rm = TRUE), p_min = min(pn, na.rm = TRUE), p_max_exact = max(pe, na.rm = TRUE), p_min_exact = min(pe, na.rm = TRUE)) })) }
sample_scores <- function(R) {   # sample x compartment x MCP: mean of the per-cell MCP scores
  rbindlist(lapply(names(R$scores), function(ct) { S <- as.data.table(R$scores[[ct]]); mc <- grep("^MCP", names(S), value = TRUE)
    data.table::melt(S[, c(lapply(.SD, mean), .(n = .N)), by = samples, .SDcols = mc], id.vars = c("samples", "n"), variable.name = "mcp", value.name = "score")[, compartment := ct] })) }
genes_of <- function(R, m, ct, n = 50) { g <- R$MCPs.full[[m]]; if (is.null(g)) return(character(0))
  u <- g[[paste0(ct, ".up")]]; d <- g[[paste0(ct, ".down")]]; c(head(u, n), head(d, n)) }
coord_resid <- function(SS, m, cts) {   # criterion 2
  W <- data.table::dcast(SS[mcp == m & compartment %in% cts & n >= ABN], samples ~ compartment, value.var = "score")
  W <- merge(W, COV, by.x = "samples", by.y = "sample")
  W <- W[stats::complete.cases(W[, c(cts, "prog_frac", "comp_pc1", "comp_pc2", "hb", "gran"), with = FALSE])]   # a sample missing one compartment drops out
  if (nrow(W) < 8) return(list(rho = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(W)))
  res <- sapply(cts, function(ct) resid(lm(W[[ct]] ~ prog_frac + comp_pc1 + comp_pc2 + hb + gran, data = W))); res <- matrix(res, ncol = length(cts), dimnames = list(NULL, cts))
  pr <- combn(cts, 2); f <- function(idx) mean(apply(pr, 2, function(p) cor(res[idx, p[1]], res[idx, p[2]], method = "spearman")))
  obs <- f(seq_len(nrow(W))); bt <- replicate(1000, f(sample(nrow(W), replace = TRUE)))
  list(rho = obs, lo = unname(quantile(bt, .025, na.rm = TRUE)), hi = unname(quantile(bt, .975, na.rm = TRUE)), n = nrow(W)) }
jacc <- function(a, b) if (!length(a) || !length(b)) NA_real_ else length(intersect(a, b)) / length(union(a, b))

## -- 1. criteria 1 and 2 per run ----------------------------------------------------------------------------------
MCPS <- list(); SCORES <- list(); CRIT <- list(); KA <- list()
for (run in RUNS) {
  R <- readRDS(file.path(OUTD, paste0("R_", run, ".rds"))); rA <- readRDS(file.path(OUTD, paste0("rA_", run, ".rds")))
  # DIALOGUE.run's normal return lacks `samples` and `covar` (they are in DLG.full.output, which 29 reads only when the plot fails)
  if (is.null(R$covar)) R$covar <- R$param$covar
  if (is.null(R$samples)) { ff <- file.path(OUTD, run, paste0("DLG.full.output_", run, ".rds"))   # DIALOGUE's own sample set (present in every compartment after abn.c), from its full output
    R$samples <- if (file.exists(ff)) readRDS(ff)$samples else Reduce(intersect, lapply(R$scores, function(S) unique(as.character(S$samples)))) }
  covar <- R$covar; ts(sprintf("[%s] %d samples in all compartments | covar %s", run, length(R$samples), paste(covar, collapse = "+")))
  MT <- mcp_table(R); if (is.null(MT)) { ts("    no MCP object"); next }
  chk(sprintf("[%s] every emp.p entry maps to a permutation rank (exact p defined)", run), TRUE, all(!is.na(exact_p(R$emp.p))))
  SS <- sample_scores(R); SCORES[[run]] <- cbind(run = run, SS)
  MT[, run := run]; for (ct in R$cell.types) MT[, (paste0("n_genes_", ct)) := vapply(mcp, function(m) length(genes_of(R, m, ct, 1e6)), 1L)]
  MCPS[[run]] <- MT
  for (i in seq_len(nrow(MT))) { m <- MT$mcp[i]; cts <- strsplit(MT$compartments[i], "+", fixed = TRUE)[[1]]
    c2 <- if (length(cts) >= 2) coord_resid(SS, m, cts) else list(rho = NA, lo = NA, hi = NA, n = NA)
    CRIT[[paste(run, m)]] <- data.table(run = run, mcp = m, compartments = MT$compartments[i], crit1_p_dlg = MT$p_max[i], crit1_p_exact = MT$p_max_exact[i],
                                        crit1 = isTRUE(MT$p_max_exact[i] <= 0.05), crit1_prereg = isTRUE(MT$p_max[i] <= 0.05),   # exact reading decided 2026-10-19 (section 8); original reading kept alongside
                                        crit2_rho = c2$rho, crit2_lo = c2$lo, crit2_hi = c2$hi, crit2_n = c2$n, crit2 = isTRUE(c2$rho > 0.3 && c2$lo > 0),
                                        crit3_jaccard_median = NA_real_, crit3_comp_pass = NA_integer_, crit3 = NA) }
  print(MT[, .(mcp, compartments, p_max_dlg = round(p_max, 3), p_max_exact, p_min_dlg = round(p_min, 3), p_min_exact)])
  ## -- 2. criterion 3: stability under 20 x 70% patient subsets ------------------------------------------------------
  if (!SKIP_RES) {
    pats <- R$samples; JC <- list()
    KEEPS <- draw_distinct(N_RES, function() sort(sample(pats, ceiling(0.7 * length(pats)))))
    chk(sprintf("[%s] the %d stability subsets are distinct and drawn before any DIALOGUE call", run, N_RES), N_RES, length(unique(KEEPS)))
    for (b in seq_len(N_RES)) {
      keep <- KEEPS[[b]]; rAs <- lapply(rA, function(r) { j <- r@samples %in% keep
        make.cell.type(name = r@name, tpm = r@tpm[, j], samples = r@samples[j], X = r@X[j, ], metadata = r@metadata[j, ], cellQ = r@cellQ[j]) })
      Rb <- run_dialogue(rAs, paste0(run, "_stab", b), covar, file.path(RESD, "stability", paste0(run, "_stab", b)))
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
    if (nrow(JC)) { fwrite(JC, file.path(OUTWS, paste0("ws_s7_stability_", run, ".csv")))
      for (m in unique(JC$mcp)) { j <- JC[mcp == m, .(med = median(jaccard, na.rm = TRUE)), by = compartment]
        CRIT[[paste(run, m)]][, `:=`(crit3_jaccard_median = median(JC[mcp == m, jaccard], na.rm = TRUE), crit3_comp_pass = j[med >= 0.5, .N], crit3 = j[med >= 0.5, .N] >= 2)] } }
  }
}
MCPS <- rbindlist(MCPS, fill = TRUE); SCORES <- rbindlist(SCORES); CRIT <- rbindlist(CRIT)
CRIT[, found := crit1 & crit2 & (crit3 %in% TRUE)]; CRIT[, found_prereg := crit1_prereg & crit2 & (crit3 %in% TRUE)]
gcols <- grep("^n_genes_", names(MCPS), value = TRUE); MCPS[, genes := if (length(gcols)) do.call(paste, c(.SD, sep = "/")) else NA_character_, .SDcols = gcols]
CRIT <- merge(CRIT, MCPS[, .(run, mcp, genes)], by = c("run", "mcp"), all.x = TRUE, sort = FALSE)
fwrite(MCPS, file.path(OUTWS, "ws_s7_mcps.csv")); fwrite(SCORES, file.path(OUTWS, "ws_s7_scores.csv"))

## -- 3. known answers: planted MCP and shuffled samples ------------------------------------------------------------
ts("[3] known answers")
pf <- file.path(OUTD, "R_planted.rds")
if (file.exists(pf)) {
  Rp <- readRDS(pf); truth <- readRDS(file.path(OUTD, "planted_truth.rds")); MTp <- mcp_table(Rp)
  print(MTp[, .(mcp, compartments, p_max_dlg = round(p_max, 3), p_max_exact)])
  rec <- rbindlist(lapply(MTp$mcp, function(m) { cts <- Rp$MCP.cell.types[[m]]; if (is.null(cts)) return(NULL)
    rbindlist(lapply(cts, function(ct) { g <- genes_of(Rp, m, ct, 200); data.table(mcp = m, compartment = ct, p = MTp[mcp == m, p_max_exact], p_dlg = MTp[mcp == m, p_max], n_genes = length(g), planted_found = length(intersect(g, truth$genes[[ct]])), planted_frac = length(intersect(g, truth$genes[[ct]])) / length(truth$genes[[ct]])) })) }))
  if (nrow(rec)) { best <- rec[, .(p = p[1], p_dlg = p_dlg[1], mean_frac = mean(planted_frac), n_comp = .N), by = mcp][order(-mean_frac)][1]
    KA[["planted"]] <- data.table(test = "planted_MCP", value = best$mean_frac, detail = sprintf("best MCP %s: exact p %.2f (DIALOGUE p %.3f), planted genes among its top 200 per compartment %.2f (over %d compartments)", best$mcp, best$p, best$p_dlg, best$mean_frac, best$n_comp))
    chk("planted MCP recovered (exact p <= 0.05 and >= 50% of planted genes)", TRUE, isTRUE(best$p <= 0.05 && best$mean_frac >= 0.5)) }
  print(rec)
}
if (!SKIP_RES && "R1" %in% RUNS) {
  rA <- readRDS(file.path(OUTD, "rA_R1.rds")); R <- readRDS(file.path(OUTD, "R_R1.rds")); fp <- 0L; fp_dlg <- 0L; nd <- 0L
  if (is.null(R$covar)) R$covar <- R$param$covar
  MAPS <- draw_distinct(N_RES, function() lapply(rA, function(r) { u <- unique(r@samples); setNames(sample(u), u) }))   # relabel samples within each compartment: pairing destroyed
  chk(sprintf("the %d shuffles are distinct permutations drawn before any DIALOGUE call", N_RES), N_RES, length(unique(MAPS)))
  for (b in seq_len(N_RES)) {
    rAs <- lapply(names(rA), function(ct) { r <- rA[[ct]]; map <- MAPS[[b]][[ct]]
      make.cell.type(name = r@name, tpm = r@tpm, samples = unname(map[r@samples]), X = r@X, metadata = r@metadata, cellQ = r@cellQ) }); names(rAs) <- names(rA)
    Rb <- run_dialogue(rAs, paste0("R1_shuf", b), R$covar, file.path(RESD, "shuffle", paste0("R1_shuf", b)))
    if (is.null(Rb) || is.null(Rb$emp.p)) next
    nd <- nd + 1L; mt <- mcp_table(Rb); fp <- fp + as.integer(any(mt$p_max_exact <= 0.05, na.rm = TRUE)); fp_dlg <- fp_dlg + as.integer(any(mt$p_max <= 0.05, na.rm = TRUE))
    ts(sprintf("    shuffle %02d: min exact p over MCPs %.2f (DIALOGUE p %.3f)", b, min(mt$p_max_exact, na.rm = TRUE), min(mt$p_max, na.rm = TRUE)))
  }
  KA[["shuffle"]] <- data.table(test = "shuffled_samples", value = if (nd) fp / nd else NA_real_, detail = sprintf("%d of %d shuffles produced an MCP with exact p <= 0.05 (%d with DIALOGUE p <= 0.05)", fp, nd, fp_dlg))
  chk("shuffled samples: false-positive rate of criterion 1 (exact p) at or below 0.05", TRUE, isTRUE(nd > 0 && fp / nd <= 0.05))
}
if (length(KA)) fwrite(rbindlist(KA), file.path(OUTWS, "ws_s7_known_answers.csv"))

## -- 4. verdict ------------------------------------------------------------------------------------------------------
fwrite(CRIT, file.path(OUTWS, "ws_s7_criteria.csv")); print(CRIT[, .(run, mcp, compartments, p_dlg = round(crit1_p_dlg, 3), p_exact = crit1_p_exact, crit1, crit1_prereg, rho = round(crit2_rho, 2), crit2_n, crit2, genes, jacc = round(crit3_jaccard_median, 2), crit3, found, found_prereg)])
verdict <- if (any(CRIT$found)) "MCP FOUND (candidate only; needs its own confirmatory pre-registration)" else "NO MCP meets the three criteria: H-A1b CLOSED"
verdict <- sprintf("%s [criterion 1 = exact permutation p, decision 2026-10-19; under the original reading of DIALOGUE's p: %s]", verdict, if (any(CRIT$found_prereg)) "MCP FOUND" else "none found")
ts("    ", verdict)
if (WRITE_GATES) {
  block <- c("<!-- A1b-S7 start -->", sprintf("## H-A1b, unsupervised test (DIALOGUE, prereg-s7-v1 + section 8) -- %s", format(Sys.Date())), "",
    "Scripts: scripts/ws/29_s7_dialogue_build_run.R, 30_s7_dialogue_evaluate.R. Criterion 1 p is shown as DIALOGUE's rank-sum p / exact permutation p",
    "(rank of the observed canonical correlation among 1 + 99 permutations, /100); the verdict uses the exact p (decision 2026-10-19, PREREGISTRATION_S7_dialogue.md section 8).", "",
    "| run | MCP | compartments | crit 1 p DIALOGUE / exact | crit 2 rho (n) | genes T/NK/B/Mono | crit 3 Jaccard | found exact / prereg |", "|---|---|---|---|---|---|---|---|",
    CRIT[, sprintf("| %s | %s | %s | %.3f / %.2f | %.2f (%s) | %s | %s | %s / %s |", run, mcp, compartments, crit1_p_dlg, crit1_p_exact, crit2_rho, crit2_n, genes, ifelse(is.na(crit3_jaccard_median), "NA", sprintf("%.2f", crit3_jaccard_median)), found, found_prereg)], "",
    if (length(KA)) paste0("Known answers: ", paste(rbindlist(KA)$detail, collapse = "; ")) else "Known answers: not run",
    sprintf("**Verdict: %s**", verdict), "<!-- A1b-S7 end -->")
  old <- readLines(GATES); s0 <- grep("<!-- A1b-S7 start -->", old); e0 <- grep("<!-- A1b-S7 end -->", old)
  writeLines(if (length(s0) && length(e0)) c(old[seq_len(s0 - 1)], block, old[-seq_len(e0)]) else c(old, "", block), GATES); ts("    GATES.md A1b-S7 block written")
}
C <- rbindlist(CHK); if (nrow(C)) fwrite(C, file.path(OUTWS, "ws_s7_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
