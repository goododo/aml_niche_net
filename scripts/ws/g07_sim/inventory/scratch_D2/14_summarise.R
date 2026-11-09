#!/usr/bin/env Rscript
# input : sim_patients_<tag>.csv.gz, sim_meta_<tag>.csv (from 13_sim.R)
# output: sum_auc.csv (uniform vs sub-population AUC per statistic), sum_null.csv (null and healthy-LODO behaviour),
#         sum_rules.csv (declared-rule error rates), printed tables for D2 sections 4-5
# what it does: summarises synthetic runs only.
suppressPackageStartupMessages(library(data.table))
fs <- list.files(".", "^sim_patients_.*csv.gz$"); R <- rbindlist(lapply(fs, function(f) fread(f)[, tag := sub("sim_patients_(.*).csv.gz", "\\1", f)]))
R[, tail_excess := mean_z - D50]; R[, neg_pihat := -pihat]
auc <- function(x, y) { x <- x[!is.na(x)]; y <- y[!is.na(y)]; if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y)); (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y)) }   # P(x > y)
ST <- c("mean_z", "F95", "ET", "logVR", "QS", "BC", "tail_excess", "LR", "neg_pihat")
A <- rbindlist(lapply(split(R[set %in% c("up", "comb") & world != "healthy_lodo"], by = c("tag", "pipe", "set", "n")), function(d) {
  rbindlist(lapply(c("sub10", "sub20", "sub40"), function(w) data.table(tag = d$tag[1], pipe = d$pipe[1], set = d$set[1], n = d$n[1], contrast = paste0(w, "_vs_uni"),
    stat = ST, auc = sapply(ST, function(s) auc(d[world == w][[s]], d[world == "uni"][[s]])))))[
  , rbind(.SD, data.table(tag = d$tag[1], pipe = d$pipe[1], set = d$set[1], n = d$n[1], contrast = "uni_vs_null", stat = "mean_z", auc = auc(d[world == "uni"]$mean_z, d[world == "null"]$mean_z)))] }))
fwrite(A, "sum_auc.csv")
options(width = 220)
cat("\n== AUC sub20 vs uniform, set up (rows: tag/pipe/n; cols: statistic) ==\n")
invisible(NULL)
cat("\n== AUC uniform vs null (mean_z), set up ==\n")
print(dcast(A[contrast == "uni_vs_null" & set == "up"], tag + pipe ~ n, value.var = "auc", fun.aggregate = function(x) round(mean(x), 2)))
# null behaviour: AML null vs healthy LODO
N <- R[set == "up" & world %in% c("null", "healthy_lodo"), .(mean_z_med = median(mean_z), mean_z_q90 = quantile(mean_z, .9), F95_med = median(F95), F95_q90 = quantile(F95, .9),
  logVR_med = median(logVR), LR_med = median(LR, na.rm = TRUE), LR_q95 = quantile(LR, .95, na.rm = TRUE), pihat_le05 = mean(pihat <= 0.5, na.rm = TRUE)), by = .(tag, pipe, world)]
fwrite(N, "sum_null.csv"); cat("\n== null AML vs healthy LODO, set up ==\n"); print(N[order(tag, pipe, world)])
# healthy LODO max per tag/pipe (the real-run reference for 'shift present')
HM <- R[world == "healthy_lodo", .(h_max_mean_z = max(mean_z), h_max_F95 = max(F95), h_max_logVR = max(logVR), h_max_LR = max(LR, na.rm = TRUE)), by = .(tag, pipe, set)]
fwrite(HM, "sum_healthy_lodo_max.csv")
# central world summaries
W <- R[set == "up" & world != "healthy_lodo", lapply(.SD, function(x) round(median(x, na.rm = TRUE), 3)), by = .(tag, pipe, world, n), .SDcols = c("mean_z", "F95", "logVR", "BC", "LR", "pihat")]
fwrite(W, "sum_world_medians.csv")
cat("\n== medians by world (E1/E06/E03, P1 and P3, up) ==\n"); print(W[pipe %in% c("P1", "P3") & n %in% c(50, 336)][order(tag, pipe, n, world)])
# ---------- declared-rule candidates (set up), error rates by world and n ----------
U <- merge(R[set == "up" & world != "healthy_lodo"], HM[set == "up"], by = c("tag", "pipe", "set"))
U[, shift := mean_z > h_max_mean_z]
rule <- function(d, which) {
  if (which == "R1_mixture") ifelse(!d$shift, "no_shift", ifelse(d$pihat <= 0.5 & d$LR >= 10, "sub", ifelse(d$pihat >= 0.9 | d$LR < 2, "uniform", "not_separable")))
  else if (which == "R2_variance") ifelse(!d$shift, "no_shift", ifelse(d$logVR > d$h_max_logVR + 0.5, "sub", ifelse(d$logVR <= d$h_max_logVR, "uniform", "not_separable")))
  else ifelse(!d$shift, "no_shift", ifelse(d$pihat <= 0.5 & d$LR >= 10 & d$logVR > d$h_max_logVR, "sub", ifelse((d$pihat >= 0.9 | d$LR < 2) & d$logVR <= d$h_max_logVR + 0.5, "uniform", "not_separable"))) }
RU <- rbindlist(lapply(c("R1_mixture", "R2_variance", "R3_both"), function(w) U[, .(tag, pipe, world, n, call = rule(.SD, w), rule = w)]))
RT <- RU[, .N, by = .(rule, tag, pipe, world, n, call)][, share := round(N / sum(N), 3), by = .(rule, tag, pipe, world, n)]
RW <- dcast(RT, rule + tag + pipe + n + world ~ call, value.var = "share", fill = 0)
fwrite(RW, "sum_rules.csv")
cat("\n== declared-rule call shares (rule R3_both and R1, P3 and P2), set up ==\n")
print(RW[rule %in% c("R1_mixture", "R3_both") & pipe %in% c("P2", "P3")][order(rule, tag, pipe, n, world)])
