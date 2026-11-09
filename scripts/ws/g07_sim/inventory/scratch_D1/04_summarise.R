#!/usr/bin/env Rscript
# INPUT : SCR/02_sim_{A_noGW_rho05,B_GW_rho05,C_GW_rho20}_raw.csv (synthetic results only)
# OUTPUT: SCR/04_oc_table.csv, SCR/04_oc_by_size.csv, SCR/04_summarise.log
# what it does: builds the operating-characteristic table of D1 section 2 and adds one derived method,
#   "camera_rstar": camera with a single preset correlation r* = the replicate's median per-set estimated
#   correlation (residual-based, status-orthogonal). Its p-values are recomputed exactly from camera_0.01's output:
#   camera's statistic is delta / sqrt(varStatPooled * (vif/m + 1/(G-m))), and delta and varStatPooled do not depend
#   on the preset correlation (limma camera.default source), so t_r = t_0.01 * sqrt((vif_0.01/m + 1/m2) / (vif_r/m + 1/m2)).
suppressPackageStartupMessages(library(data.table))
SCR <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1"
sc <- c("A_noGW_rho05", "B_GW_rho05", "C_GW_rho20")
R <- rbindlist(lapply(sc, function(s) { x <- fread(file.path(SCR, sprintf("02_sim_%s_raw.csv", s))); x[, scenario := s]; x }), fill = TRUE)
## derived camera with r* -------------------------------------------------------------------------------------------
rs <- R[method == "camera_est", .(rstar = median(cor_est, na.rm = TRUE)), by = .(scenario, rep)]
c01 <- R[method == "camera_0.01"][rs, on = .(scenario, rep)]
c01[, `:=`(G = n_genes_tested, m2 = n_genes_tested - m)]
c01[, t01 := qt(p / 2, df = G - 2, lower.tail = FALSE) * dirn]
c01[, vif01 := pmax(1, 1 + (m - 1) * 0.01)][, vifr := pmax(1, 1 + (m - 1) * rstar)]
c01[, tr := t01 * sqrt((vif01 / m + 1 / m2) / (vifr / m + 1 / m2))]
c01[, `:=`(p = 2 * pt(-abs(tr), df = G - 2), dirn = sign(tr), method = "camera_rstar")]
## sanity: recomputing camera_0.01's own p from t01 must return it
stopifnot(max(abs(2 * pt(-abs(c01$t01), df = c01$G - 2) - R[method == "camera_0.01"][rs, on = .(scenario, rep)]$p)) < 1e-8)
R2 <- rbind(R, c01[, names(R), with = FALSE])
R2[, fdr := p.adjust(p, "BH"), by = .(scenario, rep, method)]
R2[, right := is.na(true_dir) | (!is.na(dirn) & dirn == true_dir)]
OC <- R2[, .(n = .N, p05 = mean(p < 0.05), fdr05 = mean(fdr < 0.05), fdr05_right = mean(fdr < 0.05 & right)), by = .(scenario, method, type)]
fwrite(OC, file.path(SCR, "04_oc_table.csv"))
W <- dcast(OC, method ~ scenario + type, value.var = c("p05", "fdr05_right"))
BYS <- R2[type == "null", .(n = .N, p05 = mean(p < 0.05)), by = .(scenario, method, size = fifelse(m < 50, "m<50", "m>=50"))]
fwrite(BYS, file.path(SCR, "04_oc_by_size.csv"))
sink(file.path(SCR, "04_summarise.log"))
cat("r* (median per-set estimated correlation) by scenario:\n"); print(rs[, .(median_rstar = round(median(rstar), 4), min = round(min(rstar), 4), max = round(max(rstar), 4)), by = scenario])
cat("\nnull sets, rate p < 0.05 / FDR < 0.05:\n")
print(dcast(OC[type == "null"], method ~ scenario, value.var = "p05"), digits = 3)
print(dcast(OC[type == "null"], method ~ scenario, value.var = "fdr05"), digits = 3)
cat("\neffect sets, rate FDR < 0.05 with the right direction:\n")
print(dcast(OC[type != "null"], method + type ~ scenario, value.var = "fdr05_right"), digits = 3)
cat("\neffect sets, rate p < 0.05 (any direction):\n")
print(dcast(OC[type != "null"], method + type ~ scenario, value.var = "p05"), digits = 3)
cat("\nnull-set p < 0.05 rate by set size:\n"); print(dcast(BYS, method ~ scenario + size, value.var = "p05"), digits = 3)
sink(); cat(readLines(file.path(SCR, "04_summarise.log")), sep = "\n")
