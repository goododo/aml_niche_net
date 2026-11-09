#!/usr/bin/env Rscript
# input : D1's synthetic engine (scratch_D1/02_sim_geneset.R up to its scenario block); nothing else from the project
# output: review1/L1/05_sim_removal_{raw.csv,log}
# what it does: tests the 3.2 removal reading on synthetic data. A handling list H covers 10 % of background genes,
#   40 % of 12 "adjacent" sets and 7 % of the other 33 (sizes = the 45 HALLMARK sizes). Effect sets per replicate:
#   5 "uniform" (every gene 0.30 log2; H members no different) and 5 "H-carried" (only the H members of an adjacent set,
#   0.6 log2). Primary = camera est. on all genes; removal = the same fit with H rows dropped (draft 2.2/3.2: no refit).
#   Draft rule: not carried = FDR < 0.05 (BH over sets with >= 15 genes after removal) + same Direction + set shift
#   >= half of primary's; not separable = < 15 genes; carried = otherwise.
#   Per-set rule (proposed fix): s = direction-signed set shift; not carried = p <= primary BH cut-off p* and s_rem >=
#   0.5 s_prim; carried = s_rem < 0.5 s_prim; "not resolved after removal" = s kept but p > p*.
#   usage: Rscript 05_sim_removal.R <n_rep> <n_cores> <scenario A|B>
src <- readLines("/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1/02_sim_geneset.R")
eval(parse(text = src[1:(grep("^## ---- scenarios", src) - 1L)]))
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L1"
a <- commandArgs(TRUE); NREP <- as.integer(a[1]); NCORE <- as.integer(a[2]); SC <- a[3]
P <- list(seed = if (SC == "A") 50500L else 50600L, sd_l = if (SC == "A") 0 else 0.10, r_TB = 0.65, delta_D = if (SC == "A") 0 else 2,
          tau = if (SC == "A") 0 else 1, rho = 0.05, n_eff = 0L,
          eff_types = data.table(type = character(), n = integer(), beta = numeric(), pi = numeric(), where = character()), lib_B = 0.25)

one_rep <- function(rep) {
  set.seed(P$seed + rep); w <- make_world(P)
  ## handling list: 40 % of 12 adjacent sets, 7 % of the other sets, 10 % of background genes
  adj <- sample(45L, 12L)
  H <- unlist(lapply(seq_along(w$modules), function(k) { g <- w$modules[[k]]; sample(g, round(length(g) * if (k %in% adj) 0.40 else 0.07)) }))
  bg <- setdiff(seq_len(G), unlist(w$modules)); H <- c(H, sample(bg, round(0.10 * length(bg))))
  ## effects: 5 uniform sets (any), 5 H-carried sets (adjacent only); random sign
  uni <- sample(setdiff(seq_len(45), adj), 3); uni <- c(uni, sample(adj, 2)); hc <- sample(setdiff(adj, uni), 5)
  eff <- data.table(k = c(uni, hc), type = rep(c("uniform", "H_carried"), each = 5), dir = sample(c(-1, 1), 10, TRUE))
  cnt <- gen_counts(w, P, "T")
  aml <- w$status == 1
  for (i in seq_len(nrow(eff))) { e <- eff[i]; g <- w$modules[[e$k]]
    if (e$type == "H_carried") g <- intersect(g, H)
    b <- if (e$type == "uniform") 0.30 else 0.60
    s <- cnt[g, aml, drop = FALSE]                     # count-level plant (K2b form): up = Poisson add, down = thinning
    cnt[g, aml] <- if (e$dir > 0) s + matrix(rpois(length(s), s * (2^b - 1)), nrow(s)) else matrix(rbinom(length(s), s, 2^-b), nrow(s)) }
  X <- s5_design(w$status); fit <- s5_fit(cnt, X, w$status)
  y <- fit$y; lf <- fit$tab$logFC; names(lf) <- rownames(fit$tab)
  sets <- setNames(lapply(w$modules, function(i) paste0("g", i)), paste0("S", seq_along(w$modules)))
  Hn <- paste0("g", H)
  run <- function(univ) {
    yy <- y[univ, , keep.lib.sizes = TRUE]; idx <- lapply(sets, intersect, univ); m <- lengths(idx)
    ok <- names(idx)[m >= 15]
    if (!length(ok)) stop(sprintf("no set >= 15 genes; genes kept by the filter %d", nrow(y)))
    r <- camera(yy, idx[ok], X, contrast = 2, inter.gene.cor = NA, allow.neg.cor = FALSE, sort = FALSE)
    d <- data.table(set = ok, m = m[ok], dirn = ifelse(r$Direction == "Up", 1, -1), p = r$PValue, fdr = p.adjust(r$PValue, "BH"))
    d[, shift := vapply(set, function(s) mean(lf[idx[[s]]]) - mean(lf[setdiff(univ, idx[[s]])]), 1)]
    rbind(d, data.table(set = setdiff(names(idx), ok), m = m[setdiff(names(idx), ok)]), fill = TRUE) }
  A <- run(rownames(y)); Rm <- run(setdiff(rownames(y), Hn))
  r <- merge(A, Rm, by = "set", suffixes = c("_p", "_r"))
  pstar <- if (any(A$fdr < 0.05, na.rm = TRUE)) max(A[fdr < 0.05, p]) else NA_real_
  r[, k := as.integer(sub("^S", "", set))]
  r <- merge(r, eff, by = "k", all.x = TRUE); r[is.na(type), type := "null"]
  r[, `:=`(adjacent = k %in% adj, fracH = vapply(k, function(kk) mean(w$modules[[kk]] %in% H), 1), rep = rep, pstar = pstar,
           n_primary_flagged = sum(A$fdr < 0.05, na.rm = TRUE), n_removal_flagged = sum(Rm$fdr < 0.05, na.rm = TRUE))]
  r
}
t0 <- proc.time()
RL <- mclapply(seq_len(NREP), function(i) tryCatch(one_rep(i), error = function(e) e), mc.cores = NCORE, mc.set.seed = FALSE)
bad <- vapply(RL, inherits, TRUE, "error"); BADMSG <- vapply(RL[bad], conditionMessage, ""); RL <- RL[!bad]   # failed reps are skipped and counted
R <- rbindlist(RL); tt <- proc.time() - t0
fwrite(R, file.path(OUT, sprintf("05_sim_removal_%s_raw.csv", SC)))
F <- R[fdr_p < 0.05]                                                 # flagged by the primary
F[, s_p := dirn_p * shift_p]; F[, s_r := dirn_p * shift_r]           # direction-signed shifts (F1 fix)
F[, draft := fifelse(m_r < 15, "not_separable", fifelse(!is.na(fdr_r) & fdr_r < 0.05 & dirn_r == dirn_p & shift_r >= 0.5 * shift_p,
                                                        "not_carried", "carried"))]
F[, perset := fifelse(m_r < 15, "not_separable", fifelse(s_r < 0.5 * s_p, "carried",
                      fifelse(p_r <= pstar & dirn_r == dirn_p, "not_carried", "not_resolved")))]
F[, kept_shift := s_r >= 0.5 * s_p]
sink(file.path(OUT, sprintf("05_sim_removal_%s.log", SC)), split = TRUE)
cat(sprintf("scenario %s | reps %d | elapsed %.0f s | cpu self %.0f + children %.0f s\n", SC, NREP, tt[["elapsed"]],
            tt[["user.self"]] + tt[["sys.self"]], tt[["user.child"]] + tt[["sys.child"]]))
cat("reps skipped:", sum(bad), if (any(bad)) paste("(", paste(unique(BADMSG), collapse = "; "), ")") else "", "\n")
cat("flagged by the primary, by true type:\n"); print(F[, .N, by = type])
cat("\n== draft rule (3.2) by true type (counts)\n"); print(dcast(F[, .N, by = .(type, draft)], type ~ draft, value.var = "N", fill = 0))
cat("\n== per-set rule (fix) by true type (counts)\n"); print(dcast(F[, .N, by = .(type, perset)], type ~ perset, value.var = "N", fill = 0))
cat("\n== draft 'carried' calls whose direction-signed shift was KEPT (>= half): power loss, not effect loss\n")
print(F[draft == "carried", .(n_carried = .N, shift_kept = sum(kept_shift), p_r_below_pstar = sum(p_r <= pstar & dirn_r == dirn_p)), by = type])
cat("\n== BH re-ranking: draft 'carried' calls with p_r <= p* (would pass at the primary's own cut-off) and shift kept\n")
print(F[draft == "carried" & kept_shift & p_r <= pstar & dirn_r == dirn_p, .N, by = type])
cat("\nnull sets flagged by the primary:", F[type == "null", .N], "| median primary flagged per rep", median(R[, n_primary_flagged[1], by = rep]$V1),
    "| after removal", median(R[, n_removal_flagged[1], by = rep]$V1), "\n")
sink()
