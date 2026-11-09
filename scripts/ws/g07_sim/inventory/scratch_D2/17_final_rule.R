# input : sim_patients_*.csv.gz, sum_healthy_lodo_max.csv, sim_meta_*.csv, sum_auc.csv; output: final_rule_rates.csv, final_pb_match.csv, final_auc.csv
# what it does: the proposed rule (P3 reading pipeline; shift gate with sampling margin; mixture shape rule) and its
#   simulated call rates; pseudobulk matching table; AUCs for the candidate statistics.
suppressPackageStartupMessages(library(data.table)); options(width = 250)
fs <- list.files(".", "^sim_patients_.*csv.gz$"); R <- rbindlist(lapply(fs, function(f) fread(f)[, tag := sub("sim_patients_(.*).csv.gz", "\\1", f)]))
HM <- fread("sum_healthy_lodo_max.csv")
U <- merge(R[set == "up" & world != "healthy_lodo"], HM[set == "up"], by = c("tag", "pipe", "set"))
U[, se := sqrt(exp(logVR) / n)]
U[, shift := mean_z - 1.96 * se > h_max_mean_z]
U[, call := ifelse(n < 50, "too_few_cells", ifelse(!shift, "no_shift", ifelse(pihat <= 0.5 & LR >= 10, "sub", ifelse(pihat >= 0.9 | LR < 2, "uniform", "not_separable"))))]
U[, call_n20 := ifelse(!shift, "no_shift", ifelse(pihat <= 0.5 & LR >= 10, "sub", ifelse(pihat >= 0.9 | LR < 2, "uniform", "not_separable")))]
T1 <- U[, .N, by = .(tag, pipe, world, n, call = call_n20)][, share := round(N / sum(N), 3), by = .(tag, pipe, world, n)]
W1 <- dcast(T1, tag + pipe + n + world ~ call, value.var = "share", fill = 0); fwrite(W1, "final_rule_rates.csv")
cat("== final rule (gate with margin; mixture), call shares, P1/P2/P3 ==\n")
fwrite(W1[pipe %in% c("P1", "P2", "P3")][order(tag, pipe, n, world)], sep = "\t", quote = FALSE)
M <- rbindlist(lapply(list.files(".", "^sim_meta_.*csv$"), fread)); fwrite(M, "final_pb_match.csv")
cat("\n== pseudobulk log2FC of the up set (first 20 reps per world x n) ==\n"); fwrite(M[, .(tag, world, n, pb = round(pb_log2fc_up_mean, 3), pb_sd = round(pb_log2fc_up_sd, 3), aff = round(aff, 3), det_up_median, det_dn_median)][order(tag, n, world)], sep = "\t", quote = FALSE)
A <- fread("sum_auc.csv")
cat("\n== AUC P3 and P1, set up ==\n")
fwrite(dcast(A[set == "up" & pipe %in% c("P1", "P3") & stat %in% c("mean_z", "F95", "ET", "logVR", "BC", "LR", "neg_pihat") & contrast != "uni_vs_null"], tag + pipe + contrast + n ~ stat, value.var = "auc", fun.aggregate = function(x) round(mean(x), 2))[order(tag, pipe, contrast, n)], sep = "\t", quote = FALSE)
