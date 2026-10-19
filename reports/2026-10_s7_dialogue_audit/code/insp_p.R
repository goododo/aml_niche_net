suppressPackageStartupMessages(library(DIALOGUE))
ns <- asNamespace("DIALOGUE")
cat("######## ranksum.test.mat ########\n"); print(get("ranksum.test.mat", ns))
cat("\n######## numeric check: DIALOGUE emp.p as a function of the observed value's rank among 1 observed + 99 permutations ########\n")
# the function above is applied to one observed value vs 99 null values; reproduce it directly
set.seed(1); obs_rank <- 100:85
res <- t(sapply(obs_rank, function(rk) { null <- seq_len(99) + 0; x <- c(rk - 0.5, null); m <- matrix(x, ncol = 1)   # observed sits between ranks rk-1 and rk of the pooled 100 values
  p <- tryCatch(get("ranksum.test.mat", ns)(t(m), c(TRUE, rep(FALSE, 99)))[, "more"], error = function(e) NA_real_)
  c(rank_from_top = 101 - rk, n_null_ge_obs = sum(null >= rk - 0.5), exact_perm_p = (sum(null >= rk - 0.5) + 1) / 100, dialogue_p = as.numeric(p)) }))
print(round(res, 4))
cat("\nobserved emp.p values in R1/R2:", "0.0448 0.0482 0.0517 0.0555 0.0595 0.0637 0.0682 0.0883 0.2128", "\n")
