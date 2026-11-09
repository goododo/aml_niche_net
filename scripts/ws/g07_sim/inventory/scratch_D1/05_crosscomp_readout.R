#!/usr/bin/env Rscript
# INPUT : SCR/03_crosscomp_raw.csv (synthetic). OUTPUT: SCR/05_crosscomp_readout.log
# what it does: evaluates the proposed declared cross-compartment readout on the sets flagged in T
#   (T camera_est FDR < 0.05 and, as a second gate, T camera_0.01 FDR < 0.05): interaction call, ratio band, and the
#   naive "flagged in T, not in B" rule. No fitting; seconds of CPU.
suppressPackageStartupMessages(library(data.table))
SCR <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1"
R <- fread(file.path(SCR, "03_crosscomp_raw.csv"))
for (pre in c("T", "B", "I")) for (mm in c("cam01", "camest", "fry")) { pc <- paste0(pre, "_", mm, "_p"); R[, (paste0(pre, "_", mm, "_fdr")) := p.adjust(.SD[[1]], "BH"), by = rep, .SDcols = pc] }
sink(file.path(SCR, "05_crosscomp_readout.log"))
for (gate in c("camest", "cam01")) {
  F <- R[get(paste0("T_", gate, "_fdr")) < 0.05]
  F[, tdir := get(paste0("T_", gate, "_d"))]
  F[, I_larger := I_camest_fdr < 0.05 & I_camest_d == tdir]
  F[, call := fifelse(I_larger, "larger in T than B", fifelse(ratio_BT >= 0.5, "shared (ratio >= 0.5, no interaction)", "not resolved"))]
  F[, naive := B_camest_p >= 0.05]
  cat(sprintf("\n== sets flagged in T by camera %s (FDR < 0.05): n by type\n", gate)); print(F[, .N, by = type])
  cat("declared call by true type (row shares):\n"); print(dcast(F[, .N, by = .(type, call)][, share := round(N / sum(N), 3), by = type], type ~ call, value.var = "share", fill = 0))
  cat("naive rule 'not in B' (B camest p >= 0.05) by true type:\n"); print(F[, .(naive_Tspecific = round(mean(naive), 3)), by = type])
  cat("ratio B/T among flagged sets, quantiles 10/50/90 by true type:\n"); print(F[, .(q10 = round(quantile(ratio_BT, .1), 2), q50 = round(median(ratio_BT), 2), q90 = round(quantile(ratio_BT, .9), 2)), by = type])
}
sink(); cat(readLines(file.path(SCR, "05_crosscomp_readout.log")), sep = "\n")
