suppressPackageStartupMessages(library(DIALOGUE))
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
rA <- readRDS(file.path(D, "rA_planted_v1_tpm_only.rds")); TR <- readRDS(file.path(D, "planted_truth_v1_tpm_only.rds")); f <- TR$factor
vis <- function(r, f) { A <- rowsum(r@X, r@samples) / as.vector(table(r@samples)[sort(unique(r@samples))]); max(abs(cor(A, f[rownames(A)]))) }
set.seed(1)
for (ct in names(rA)) { r <- rA[[ct]]; obs <- vis(r, f); A <- rowsum(r@X, r@samples) / as.vector(table(r@samples)[sort(unique(r@samples))])
  null <- replicate(500, max(abs(cor(A, sample(f[rownames(A)])))))
  cat(sprintf("%-5s v1 (X unplanted) with its own f: %.2f | chance reference (500 random f): median %.2f  P95 %.2f  max %.2f | n samples %d\n", ct, obs, median(null), quantile(null, .95), max(null), nrow(A))) }
