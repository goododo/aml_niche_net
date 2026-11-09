# input : real_background_profile_nonF402.rds, real_T_percell_umi.csv.gz, real_sample_summary.csv (from 01/02)
# output: functions (sourced by 11_*, 12_*); no file written here
# what it does: synthetic T cells = Poisson-gamma counts around a status-free pooled profile of NON-F402 genes; synthetic
#   "up"/"down" sets are random tested non-F402 genes (no imprint-gene value is used anywhere). Fast exact UCell via
#   per-cell background count histograms (validated against UCell::ScoreSignatures_UCell in 11_validate.R).
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })
RB <- readRDS("real_background_profile_nonF402.rds"); PC <- fread("real_T_percell_umi.csv.gz"); SS <- fread("real_sample_summary.csv")
LAB <- c("CD4N","CD4CM","CD4EM","Treg","CD8N","CD8CM","CD8EM1","CD8EM2","CD8TRM")
FULL <- c("CD4 Naive","CD4 Central Memory","CD4 Effector Memory","CD4 Regulatory","CD8 Naive","CD8 Central Memory","CD8 Effector Memory 1","CD8 Effector Memory 2","CD8 Tissue Resident Memory")
PC[, lab := LAB[match(fine, FULL)]]
CM <- 400L
make_sets <- function(seed, n_up = 156L, n_dn = 214L, mode = "random") {
  set.seed(seed); pool <- RB$tested_nonF402; ab <- RB$bg[pool]
  if (mode == "low") pool <- pool[ab <= median(ab)]
  if (mode == "high") pool <- pool[ab > median(ab)]
  g <- sample(pool, n_up + n_dn); list(up = g[seq_len(n_up)], dn = g[n_up + seq_len(n_dn)])
}
# background bank: per cell, G_gt[c+1] = number of background genes with count > c (c = 0..CM); last col = count > CM
make_bank <- function(L, p_bg, phi, seed, keep_full = 0L) {
  set.seed(seed); nb <- length(L); G <- length(p_bg)
  GT <- matrix(0L, nb, CM + 1L); FULLC <- vector("list", keep_full)
  for (i in seq_len(nb)) {
    lam <- L[i] * p_bg; if (phi > 0) lam <- lam * rgamma(G, shape = 1 / phi, scale = phi)
    y <- rpois(G, lam); y <- y[y > 0]
    h <- tabulate(pmin(y, CM + 1L), nbins = CM + 1L)          # h[c] = genes with count c (c <= CM), h[CM+1] = > CM
    GT[i, ] <- as.integer(rev(cumsum(rev(c(h, 0L))))[seq_len(CM + 1L)])   # GT[i, c+1] = #genes with count > c, c = 0..CM
    if (i <= keep_full) FULLC[[i]] <- rpois(0, 1)
  }
  list(L = L, GT = GT)
}
# exact UCell (UCell 2.10.1 u_stat; ties.method "average"; zeros always >= maxRank since the axis has ~36k genes)
ucell_fast <- function(Y, bankGT, bank_id, sets_idx, maxRank = 1500) {
  # Y: cells x sig-genes counts; full rank of gene with count c>0 = sig_rank_avg + bg_gt(c) + bg_eq(c)/2
  R <- matrixStats::rowRanks(-Y, ties.method = "average")
  cc <- pmin(Y, CM)                                              # counts above CM: approx (bg_gt = #>CM)
  idx_gt <- cbind(rep(bank_id, ncol(Y)), as.vector(cc) + 1L)     # col c+1 -> #bg genes with count > c
  idx_ge <- cbind(rep(bank_id, ncol(Y)), as.vector(cc))          # col c -> #bg genes with count > c-1 = >= c
  bg_gt <- bankGT[idx_gt]; bg_ge <- ifelse(as.vector(cc) >= 1L, bankGT[cbind(rep(bank_id, ncol(Y)), pmax(as.vector(cc), 1L))], 0L)
  bg_eq <- bg_ge - bg_gt
  full <- as.vector(R) + bg_gt + bg_eq / 2
  full[as.vector(Y) == 0] <- maxRank
  full <- pmin(full, maxRank); dim(full) <- dim(Y)
  sapply(sets_idx, function(j) { n <- length(j); rmin <- n * (n + 1) / 2; 1 - (rowSums(full[, j, drop = FALSE]) - rmin) / (n * maxRank - rmin) })
}
