#!/usr/bin/env Rscript
# input : review1/L1/04_sim_k5_stacked_raw.csv, 04_sim_k5_universe.csv (synthetic results of 04_sim_k5_stacked.R)
# output: review1/L1/04_sim_k5_stacked.log (rewritten; header lines of the failed summary kept from 04's stdout)
# what it does: applies the 3.3(b) reading to the planted sets and summarises against the K5 bars. No fitting.
suppressPackageStartupMessages(library(data.table)); options(width = 160)
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L1"
R <- fread(file.path(OUT, "04_sim_k5_stacked_raw.csv")); U <- fread(file.path(OUT, "04_sim_k5_universe.csv"))
U[, tag := fifelse(world == 99, "all20", "realistic")]
R[, rho := sX / sT]
R[, call := fifelse(FDR < 0.05 & Direction == "Up", "larger_T", fifelse(FDR < 0.05 & Direction == "Down", "larger_X",
            fifelse(rho >= 0.5, "shared", "not_resolved")))]
hdr <- readLines(file.path(OUT, "04_sim_k5_stacked.stdout"))[1:3]
sink(file.path(OUT, "04_sim_k5_stacked.log"), split = TRUE)
cat(hdr, sep = "\n")
cat("\nuniverse (genes tested; planted genes lost from the common universe after the refit):\n")
print(U[, .(nT0 = median(nT0), nT = median(nT), nX0 = median(nX0), nX = median(nX), common0 = median(common0), common = median(common), planted_lost = sum(planted_lost)), by = .(tag, X, type)])
cat("\n== 3.3(b) reading on planted sets\n")
K <- R[n_planted >= 45][, .(n = .N, larger_T = sum(call == "larger_T"), shared = sum(call == "shared"), not_res = sum(call == "not_resolved"),
        larger_X = sum(call == "larger_X"), T_flag = sum(T_FDR < 0.05 & T_Direction == "Up"), rho_med = round(median(rho), 2),
        rho_q10 = round(quantile(rho, .1), 2), rho_q90 = round(quantile(rho, .9), 2), I_p_med = signif(median(PValue), 2),
        I_p_max = signif(max(PValue), 2)), by = .(tag, X, type)]
print(K)
cat("\n== null module sets (no planted gene)\n")
print(R[n_planted == 0, .(n = .N, I_FDR05 = round(mean(FDR < 0.05), 3), I_p05 = round(mean(PValue < 0.05), 3), larger_T = sum(call == "larger_T")), by = .(tag, X, type)])
cat("\n== module sets with 1-10 planted genes (partly contaminated): larger_T calls\n")
print(R[n_planted %between% c(1, 10), .(n = .N, larger_T = sum(call == "larger_T"), I_FDR05 = round(mean(FDR < 0.05), 3)), by = .(tag, X, type)])
cat("\n== K5 bars (realistic tag): (i) larger_T >= 14/20; (ii) larger_T <= 2/20 and shared >= 14/20; (iii) larger_T > 10/20 triggers the wording\n")
print(K[tag == "realistic", .(X, type, n, larger_T, shared, not_res,
  bar = fifelse(type == "T_only", fifelse(larger_T >= 14, "(i) PASS", "(i) FAIL"),
        fifelse(type == "shared", fifelse(larger_T <= 2 & shared >= 14, "(ii) PASS", "(ii) FAIL"),
        fifelse(larger_T > 10, "(iii) wording fires", "(iii) wording not fired"))))])
sink()
