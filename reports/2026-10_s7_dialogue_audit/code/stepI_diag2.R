suppressPackageStartupMessages({ library(DIALOGUE); library(data.table) })
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
F <- readRDS(file.path(D, "planted", "DLG.full.output_planted.rds")); TR <- readRDS(file.path(D, "planted_truth.rds")); f <- TR$factor
X <- F$sample.PCs; su <- F$samples; rA <- readRDS(file.path(D, "rA_planted.rds"))
cat("=== planted direction in DIALOGUE's processed space (sample medians over abundant samples, per-PC standardised, 13 samples) ===\n")
S <- list()
for (ct in names(X)) { r <- rA[[ct]]; ab <- names(which(table(r@samples) >= 15)); j <- r@samples %in% ab
  X1 <- DIALOGUE:::average.mat.rows(r@X[j, ], r@samples[j], f = matrixStats::colMedians)      # as DIALOGUE1
  sdk <- apply(X1, 2, sd); kept <- colnames(X[[ct]]); w <- TR$dir[[ct]][kept] * sdk[kept]      # direction in standardised coordinates
  s <- as.vector(X[[ct]] %*% w); names(s) <- rownames(X[[ct]]); S[[ct]] <- s
  cat(sprintf("  %-5s kept PCs %d | cor(score, f) = %.3f | share of planted-direction weight on kept PCs: %.2f\n", ct, length(kept), cor(s, f[names(s)]), sum(abs(TR$dir[[ct]][kept])) / sum(abs(TR$dir[[ct]])))) }
M <- do.call(cbind, S); cat("pairwise cross-compartment correlation of these scores:\n"); print(round(cor(M), 3))
cat("\nmean pairwise correlation reachable by the planted direction:", round(mean(cor(M)[upper.tri(cor(M))]), 3), "| DIALOGUE MCP1/2/3 mean pairwise canonical correlation:", paste(round(rowMeans(F$cca.cor$R), 3), collapse = " / "), "\n")
cat("\n=== permutation null of DIALOGUE's canonical correlations (its own routine, 30 permuted data sets) ===\n")
set.seed(3); NUL <- t(sapply(1:30, function(i) { v <- DIALOGUE:::DIALOGUE1.PMD.empirical(X, k = 3, seed = sample(1:10000, 1), n1 = 1); c(MCP1_mean = mean(v[1, ]), MCP1_max = max(v[1, ]), MCP2_mean = mean(v[2, ]), MCP3_mean = mean(v[3, ])) }))
print(round(apply(NUL, 2, quantile, c(.05, .5, .95)), 3))
cat("\nobserved MCP1/2/3 mean pairwise canonical correlation in the planted run:", paste(round(rowMeans(F$cca.cor$R), 3), collapse = " / "), "\n")
