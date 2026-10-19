#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../ws_dialogue/{R_<run><tag>.rds, planted_truth<tag>.rds, <run><tag>/DLG.full.output_<run>.rds}  (29 with --tag)
# OUTPUT : results/tables/ws/ws_s7_planted_toolcheck.csv   one row per check: the MCP tracking the planted factor, its compartments,
#          largest exact p over its pairs, planted-gene recovery per compartment, verdict under the prereg section 4.1 rule
# what it does: evaluates the planted known answer for DIALOGUE configurations run as tool calibration (PREREGISTRATION_S7_dialogue.md
#   section 8, item 7). No real-data criterion is computed here. The planted MCP is the one whose sample-mean step-I scores correlate
#   best (mean |cor| over compartments) with the planted factor; exact p = rank / 100 recovered from DIALOGUE's rank-sum emp.p (as 30).
#   usage: Rscript scripts/ws/31_s7_planted_toolcheck.R --checks=planted:_pc10,plantedR2:_n23,plantedR2:_n23_pc10   (run:tag pairs)
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table); library(DIALOGUE); library(here) })
.a <- commandArgs(TRUE); .get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
CHECKS <- strsplit(strsplit(.get("checks", "planted:_pc10,plantedR2:_n23,plantedR2:_n23_pc10"), ",")[[1]], ":")
source(here::here("scripts", "config", "config_paths.R")); OUTD <- file.path(LARGE1_DIR, "ws_dialogue")
LUT <- { rt <- get("ranksum.test", asNamespace("DIALOGUE")); data.table(rank = 1:100, dlg_p = vapply(1:100, function(r) as.numeric(rt(100 - r + 0.5, 1:99, alternative = "greater")), 1)) }
exact_p <- function(p) vapply(as.numeric(p), function(x) { if (is.na(x)) return(NA_real_); i <- which.min(abs(LUT$dlg_p - x)); if (abs(LUT$dlg_p[i] - x) > 1e-6) NA_real_ else LUT$rank[i] / 100 }, 1)
top200 <- function(R, m, ct) { g <- R$MCPs.full[[m]]; if (is.null(g)) return(character(0)); c(head(g[[paste0(ct, ".up")]], 200), head(g[[paste0(ct, ".down")]], 200)) }
OUT <- list()
for (ck in CHECKS) { run <- ck[1]; tag <- if (length(ck) > 1) ck[2] else ""
  fR <- file.path(OUTD, paste0("R_", run, tag, ".rds")); fT <- file.path(OUTD, paste0("planted_truth", tag, ".rds")); fF <- file.path(OUTD, paste0(run, tag), paste0("DLG.full.output_", run, ".rds"))
  if (!all(file.exists(c(fR, fT, fF)))) { message("[", run, tag, "] missing: ", paste(c(fR, fT, fF)[!file.exists(c(fR, fT, fF))], collapse = ", ")); next }
  R <- readRDS(fR); TR <- readRDS(fT); F <- readRDS(fF); f <- TR$factor; cts_all <- R$cell.types; k <- R$k["DIALOGUE"]
  # which MCP tracks the planted factor at step I (sample means of the per-cell CCA scores)
  trk <- sapply(cts_all, function(ct) { S <- F$cca.scores[[ct]][, seq_len(k), drop = FALSE]; A <- rowsum(S, F$samples.cells[[ct]]) / as.vector(table(F$samples.cells[[ct]])[rownames(rowsum(S, F$samples.cells[[ct]]))]); cor(A, f[rownames(A)])[, 1] })
  trk <- matrix(trk, nrow = k, dimnames = list(paste0("MCP", seq_len(k)), cts_all)); m <- names(which.max(rowMeans(abs(trk))))
  cts <- R$MCP.cell.types[[m]]; E <- as.matrix(R$emp.p)
  pp <- if (length(cts) >= 2) { pr <- combn(cts, 2); apply(pr, 2, function(p) { a <- paste(p, collapse = "_"); b <- paste(rev(p), collapse = "_"); if (a %in% colnames(E)) E[m, a] else E[m, b] }) } else NA_real_
  rec <- sapply(cts_all, function(ct) { g <- top200(R, m, ct); length(intersect(g, TR$genes[[ct]])) / length(TR$genes[[ct]]) })
  ngen <- sapply(cts_all, function(ct) length(top200(R, m, ct)))
  pass <- isTRUE(length(cts) >= 2 && max(exact_p(pp)) <= 0.05 && all(rec[cts] >= 0.5))
  OUT[[length(OUT) + 1]] <- data.table(check = paste0(run, tag), npcs = ncol(F$sample.PCs[[1]]), n_samples = length(F$samples), planted_mcp = m,
    track_cor = paste(sprintf("%s %+.2f", cts_all, trk[m, ]), collapse = " "), compartments = paste(cts, collapse = "+"),
    p_max_dlg = if (all(is.na(pp))) NA_real_ else max(pp), p_max_exact = if (all(is.na(pp))) NA_real_ else max(exact_p(pp)),
    genes = paste(ngen, collapse = "/"), planted_frac = paste(sprintf("%s %.2f", cts_all, rec), collapse = " "), pass_section_4_1 = pass)
  message(sprintf("[%s%s] planted MCP %s | track cor %s | compartments %s | max exact p %s | planted frac %s | PASS %s", run, tag, m, paste(round(trk[m, ], 2), collapse = "/"), paste(cts, collapse = "+"), if (all(is.na(pp))) "NA" else max(exact_p(pp)), paste(round(rec, 2), collapse = "/"), pass))
}
O <- rbindlist(OUT); if (nrow(O)) { fwrite(O, here::here("results", "tables", "ws", "ws_s7_planted_toolcheck.csv")); print(O[, .(check, npcs, n_samples, planted_mcp, compartments, p_max_exact, genes, pass_section_4_1)]) }
