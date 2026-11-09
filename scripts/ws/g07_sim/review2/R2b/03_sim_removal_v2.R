#!/usr/bin/env Rscript
# input : D1's synthetic engine (inventory/scratch_D1/02_sim_geneset.R up to its scenario block); nothing from the project
# output: review2/R2b/03_sim_removal_v2_{raw.csv,log}
# what it does: replays the DRAFT v2 3.2 removal rule (direction-matched universes, own p against p* = largest flagged
#   primary p, direction-signed shift s signed by the PRIMARY direction, half-shift, 10-gene floor) on synthetic worlds.
#   H is split into H_up / H_down halves plus a small shared IEG part; H covers 40 % of 12 adjacent sets, 7 % of the rest.
#   Scenarios: solo_adj  = ONE uniform 0.35-log2 set (an adjacent set; H members move like the rest), nothing else;
#              solo_non  = the same in a non-adjacent set;
#              mixed     = one adjacent set whose H_dir members move 1.2 and other members 0.15 (same sign) + 4 uniform
#                          0.35 sets elsewhere (so p* is not the mixed set's own p).
#   usage: Rscript 03_sim_removal_v2.R <n_rep> <n_cores>
src <- readLines("/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1/02_sim_geneset.R")
SCR <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1"
eval(parse(text = src[1:(grep("^## ---- scenarios", src) - 1L)]))
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review2/R2b"
a <- commandArgs(TRUE); NREP <- as.integer(a[1]); NCORE <- as.integer(a[2])
P <- list(sd_l = 0, r_TB = 0.65, delta_D = 0, tau = 0, rho = 0.05, n_eff = 0L,
          eff_types = data.table(type = character(), n = integer(), beta = numeric(), pi = numeric(), where = character()), lib_B = 0.25)
plant <- function(cnt, g, aml, b, dir) { s <- cnt[g, aml, drop = FALSE]     # K2b form: up = Poisson add, down = thinning
  cnt[g, aml] <- if (dir > 0) s + matrix(rpois(length(s), s * (2^b - 1)), nrow(s)) else matrix(rbinom(length(s), s, 2^-b), nrow(s)); cnt }
one_rep <- function(rep, SC) {
  set.seed(71000 + 1000 * match(SC, c("solo_adj", "solo_non", "mixed")) + rep); w <- make_world(P)
  adj <- sample(45L, 12L)
  H <- unlist(lapply(seq_along(w$modules), function(k) { g <- w$modules[[k]]; sample(g, round(length(g) * if (k %in% adj) 0.40 else 0.07)) }))
  bg <- setdiff(seq_len(G), unlist(w$modules)); H <- c(H, sample(bg, round(0.10 * length(bg))))
  ieg <- sample(H, round(0.04 * length(H))); rest <- setdiff(H, ieg); hu <- sample(rest, length(rest) %/% 2)
  Hup <- c(hu, ieg); Hdn <- c(setdiff(rest, hu), ieg)                      # H_up = IEG + delay-up; H_down = IEG + delay-down
  cnt <- gen_counts(w, P, "T"); aml <- w$status == 1
  eff <- switch(SC,
    solo_adj = data.table(k = sample(adj, 1), type = "uniform", dir = sample(c(-1, 1), 1)),
    solo_non = data.table(k = sample(setdiff(1:45, adj), 1), type = "uniform", dir = sample(c(-1, 1), 1)),
    mixed    = { m <- sample(adj, 1); data.table(k = c(m, sample(setdiff(1:45, m), 4)), type = c("mixed", rep("uniform", 4)),
                                                 dir = sample(c(-1, 1), 5, TRUE)) })
  for (i in seq_len(nrow(eff))) { e <- eff[i]; g <- w$modules[[e$k]]; Hd <- if (e$dir > 0) Hup else Hdn
    if (e$type == "uniform") cnt <- plant(cnt, g, aml, 0.35, e$dir)
    if (e$type == "mixed") { cnt <- plant(cnt, intersect(g, Hd), aml, 1.2, e$dir); cnt <- plant(cnt, setdiff(g, Hd), aml, 0.15, e$dir) } }
  X <- s5_design(w$status); fit <- s5_fit(cnt, X, w$status)
  y <- fit$y; lf <- fit$tab$logFC; names(lf) <- rownames(fit$tab)
  sets <- setNames(lapply(w$modules, function(i) paste0("g", i)), paste0("S", seq_along(w$modules)))
  run <- function(univ, floor_n) {
    yy <- y[univ, , keep.lib.sizes = TRUE]; idx <- lapply(sets, intersect, univ); m <- lengths(idx); ok <- names(idx)[m >= floor_n]
    r <- camera(yy, idx[ok], X, contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
    d <- data.table(set = ok, m = m[ok], dirn = ifelse(r$Direction == "Up", 1, -1), p = r$PValue, fdr = p.adjust(r$PValue, "BH"))
    d[, shift := vapply(set, function(s) mean(lf[idx[[s]]]) - mean(lf[setdiff(univ, idx[[s]])]), 1)]; d }
  A <- run(rownames(y), 15); RU <- run(setdiff(rownames(y), paste0("g", Hup)), 10); RD <- run(setdiff(rownames(y), paste0("g", Hdn)), 10)
  F <- A[fdr < 0.05]; if (!nrow(F)) return(data.table(rep = rep, SC = SC, n_flag = 0L))
  pstar <- max(F$p)
  R <- rbind(merge(F[dirn > 0], RU, by = "set", suffixes = c("_p", "_r")), merge(F[dirn < 0], RD, by = "set", suffixes = c("_p", "_r")))
  R[, k := as.integer(sub("^S", "", set))]; R <- merge(R, eff[, .(k, type)], by = "k", all.x = TRUE); R[is.na(type), type := "null"]
  R[, `:=`(rep = rep, SC = SC, pstar = pstar, n_flag = nrow(F), weakest = p_p == pstar, s_p = dirn_p * shift_p, s_r = dirn_p * shift_r)]
  R[, reading := fifelse(m_r < 10, "consists_mostly", fifelse(s_r < 0.5 * s_p, "carried",
                  fifelse(dirn_r == dirn_p & p_r <= pstar, "not_carried", fifelse(dirn_r == dirn_p, "not_resolved", "UNDEFINED_dir_flip"))))]
  R[, carried_but_still_sig := reading == "carried" & dirn_r == dirn_p & p_r <= pstar]
  R
}
t0 <- proc.time()
RL <- unlist(lapply(c("solo_adj", "solo_non", "mixed"), function(SC) mclapply(seq_len(NREP), function(i) tryCatch(one_rep(i, SC), error = function(e) e),
             mc.cores = NCORE, mc.set.seed = FALSE)), recursive = FALSE)
bad <- vapply(RL, inherits, TRUE, "error"); BADMSG <- unique(vapply(RL[bad], conditionMessage, "")); RL <- RL[!bad]
R <- rbindlist(RL, fill = TRUE); tt <- proc.time() - t0
fwrite(R, file.path(OUT, "03_sim_removal_v2_raw.csv"))
sink(file.path(OUT, "03_sim_removal_v2.log"), split = TRUE)
cat(sprintf("reps per scenario %d | failed %d %s | elapsed %.0f s | cpu children %.0f s\n", NREP, sum(bad), paste(BADMSG, collapse = ";"),
            tt[["elapsed"]], tt[["user.child"]] + tt[["sys.child"]]))
cat("\nflagged sets per replicate (primary, FDR < 0.05):\n"); print(R[, .(n_flag = n_flag[1]), by = .(SC, rep)][, .N, by = .(SC, n_flag)][order(SC, n_flag)])
cat("\n== v2 reading of flagged EFFECT sets, by scenario and true type\n")
print(dcast(R[!is.na(reading) & type != "null", .N, by = .(SC, type, reading)], SC + type ~ reading, value.var = "N", fill = 0))
cat("\n== same, only the weakest flagged set (its own p is p*)\n")
print(dcast(R[!is.na(reading) & type != "null" & weakest == TRUE, .N, by = .(SC, type, reading)], SC + type ~ reading, value.var = "N", fill = 0))
cat("\n== same, only replicates with exactly one flagged set\n")
print(dcast(R[!is.na(reading) & type != "null" & n_flag == 1, .N, by = .(SC, type, reading)], SC + type ~ reading, value.var = "N", fill = 0))
cat("\n== 'carried' calls that stay significant at p* in the same direction after removal (7.3 says 'without them it is no longer carried')\n")
print(R[reading == "carried", .(n_carried = .N, still_sig_at_pstar = sum(carried_but_still_sig), median_p_r = median(p_r), median_pstar = median(pstar)), by = .(SC, type)])
cat("\n== ratio s_r / s_p for the effect sets (median, P10, P90)\n")
print(R[!is.na(reading) & type != "null", .(med = median(s_r / s_p), p10 = quantile(s_r / s_p, .1), p90 = quantile(s_r / s_p, .9)), by = .(SC, type)])
cat("\nnull sets flagged:", R[type == "null", .N], "| their readings:\n"); print(R[type == "null", .N, by = reading])
sink()
