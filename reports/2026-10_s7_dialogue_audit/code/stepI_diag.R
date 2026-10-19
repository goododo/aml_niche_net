suppressPackageStartupMessages({ library(DIALOGUE); library(data.table) })
ns <- asNamespace("DIALOGUE"); cat("######## DIALOGUE1.PMD ########\n"); cat(head(deparse(get("DIALOGUE1.PMD", ns)), 60), sep = "\n")
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
F <- readRDS(file.path(D, "planted", "DLG.full.output_planted.rds")); TR <- readRDS(file.path(D, "planted_truth.rds")); f <- TR$factor
X <- F$sample.PCs; su <- F$samples; cat("\nsamplesU:", length(su), "| features kept per compartment:", paste(sapply(X, ncol), collapse = "/"), "\n")
ws <- F$cca$ws
cat("\n=== DIALOGUE's canonical variates (X %*% ws) vs planted f, per compartment and MCP ===\n")
for (ct in names(X)) { Y <- X[[ct]] %*% ws[[ct]]; cat(sprintf("  %-5s %s | nonzero weights per MCP: %s\n", ct, paste(sprintf("MCP%d %+.2f", 1:ncol(Y), cor(Y, f[rownames(X[[ct]])])), collapse = "  "), paste(colSums(ws[[ct]] != 0), collapse = "/"))) }
cat("\n=== observed canonical correlations per MCP and pair (cca.cor$R) ===\n"); print(round(F$cca.cor$R, 3))
cat("\n=== the planted direction itself in the same processed space: sample score = X1 %*% (dir restricted to kept PCs, rescaled by each PC's sample SD) ===\n")
rA <- readRDS(file.path(D, "rA_planted.rds")); S <- list()
for (ct in names(X)) { r <- rA[[ct]]; A <- rowsum(r@X, r@samples) / as.vector(table(r@samples)[sort(unique(r@samples))]); s <- as.vector(A %*% TR$dir[[ct]]); names(s) <- rownames(A); S[[ct]] <- s[su]
  cat(sprintf("  %-5s cor(planted-direction sample score, f) over samplesU = %.3f\n", ct, cor(S[[ct]], f[su]))) }
M <- do.call(cbind, S); cat("pairwise cor of the planted-direction scores across compartments (the correlation a CCA could reach):\n"); print(round(cor(M), 2))
cat("\n=== null canonical correlations: DIALOGUE1.PMD.empirical on these X (100 perms, as DIALOGUE does) ===\n")
set.seed(1); P <- DIALOGUE:::DIALOGUE1.PMD.empirical(X, k = 3, n1 = 100); cat("emp.p recomputed (MCP x pair):\n"); print(round(P, 3))
# one permuted run with full output to see how high permuted canonical correlations get
v <- DIALOGUE:::DIALOGUE1.PMD.empirical(X, k = 3, seed = 7, n1 = 1); cat("one permuted-data run: canonical correlations per MCP x pair:\n"); print(round(v, 2))
v2 <- DIALOGUE:::DIALOGUE1.PMD.empirical(X, k = 3, seed = 99, n1 = 1); print(round(v2, 2))
