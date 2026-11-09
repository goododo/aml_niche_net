# input : 10_engine.R; output: validate_ucell.csv
# what it does: builds 200 explicit synthetic cells (full 36k-gene count vectors), scores two synthetic sets with
#   UCell::ScoreSignatures_UCell (maxRank 1500 and 3000) and with ucell_fast; reports max |difference|.
source("10_engine.R"); suppressPackageStartupMessages({ library(Matrix); library(UCell) })
PHI <- 0.5; sets <- make_sets(11); sig <- c(sets$up, sets$dn)
p_all <- RB$bg / sum(RB$bg); bgn <- setdiff(names(p_all), sig)
set.seed(5); L <- sample(PC$umi_axis, 200)
L[1:20] <- seq(5000, 9000, length.out = 20)                       # include deep cells where the cap binds
Ybg <- sapply(L, function(l) { lam <- l * p_all[bgn] * rgamma(length(bgn), 1 / PHI, scale = PHI); rpois(length(lam), lam) })
mult <- c(rep(3, length(sets$up)), rep(1, length(sets$dn)))     # some up-shift so scores are spread
Ysig <- t(sapply(L, function(l) { lam <- l * p_all[sig] * mult * rgamma(length(sig), 1 / PHI, scale = PHI); rpois(length(lam), lam) }))
GT <- t(apply(Ybg, 2, function(y) { y <- y[y > 0]; h <- tabulate(pmin(y, CM + 1L), nbins = CM + 1L); as.integer(rev(cumsum(rev(c(h, 0L))))[seq_len(CM + 1L)]) }))
M <- rbind(Ybg, t(Ysig)); rownames(M) <- c(bgn, sig); colnames(M) <- paste0("c", 1:200)
out <- list()
for (mr in c(1500, 3000)) {
  u <- ScoreSignatures_UCell(Matrix(M, sparse = TRUE), features = list(up = sets$up, dn = sets$dn), maxRank = mr, ncores = 1)
  f <- ucell_fast(Ysig, GT, 1:200, list(up = seq_along(sets$up), dn = length(sets$up) + seq_along(sets$dn)), maxRank = mr)
  out[[as.character(mr)]] <- data.table(maxRank = mr, set = c("up", "dn"), max_abs_diff = c(max(abs(u[, 1] - f[, 1])), max(abs(u[, 2] - f[, 2]))),
                                        mean_score = c(mean(u[, 1]), mean(u[, 2])), nfeat_median = median(colSums(M > 0)), n_cells = 200)
}
V <- rbindlist(out); print(V); fwrite(V, "validate_ucell.csv")
cat("signature genes detected per cell (no shift on dn set): up median", median(rowSums(Ysig[, seq_along(sets$up)] > 0)), " dn median", median(rowSums(Ysig[, -seq_along(sets$up)] > 0)), "\n")
