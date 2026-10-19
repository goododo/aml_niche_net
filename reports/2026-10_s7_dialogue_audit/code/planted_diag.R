suppressPackageStartupMessages({ library(DIALOGUE); library(data.table) })
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
P <- readRDS(file.path(D, "R_planted.rds")); F <- readRDS(file.path(D, "planted", "DLG.full.output_planted.rds")); TR <- readRDS(file.path(D, "planted_truth.rds"))
cat("=== MCPs.full (all) ===\n"); for (m in names(P$MCPs.full)) { g <- P$MCPs.full[[m]]; cat(m, ":", if (is.null(g)) "NULL" else paste(sprintf("%s=%d", names(g), lengths(g)), collapse = " "), "\n") }
cat("\n=== sig1 per compartment (full output) ===\n"); for (ct in names(F$sig1)) { x <- F$sig1[[ct]]; cat(sprintf("  %-5s %s\n", ct, if (!length(x)) "-" else paste(sprintf("%s=%d", names(x), lengths(x)), collapse = " "))) }
cat("\n=== step I: planted genes among the top-200 genes correlated with each MCP's cell scores (cca.sig) ===\n")
for (ct in names(F$cca.sig)) { s <- F$cca.sig[[ct]]; for (k in names(s)) cat(sprintf("  %-5s %-10s n=%3d planted=%3d\n", ct, k, length(s[[k]]), sum(s[[k]] %in% TR$genes[[ct]]))) }
cat("\n=== step I: |cor| of planted genes vs others with the MCP2 cell score (cca.gene.cor1) ===\n")
for (ct in names(F$cca.gene.cor1)) { C <- F$cca.gene.cor1[[ct]]; if (!"MCP2" %in% colnames(C)) next; v <- abs(C[, "MCP2"]); pl <- names(v) %in% TR$genes[[ct]]
  cat(sprintf("  %-5s planted median |cor| %.3f (P90 %.3f) | others median %.3f (P99 %.3f) | planted genes in the top 200 by |cor|: %d\n", ct, median(v[pl]), quantile(v[pl], .9), median(v[!pl]), quantile(v[!pl], .99), sum(pl[order(-v)][1:200]))) }
cat("\n=== step II: gene.pval, program MCP2 -- planted genes' p-values vs others ===\n")
for (ct in names(P$gene.pval)) { G <- as.data.table(P$gene.pval[[ct]]); if (!nrow(G)) { cat("  ", ct, ": no table\n"); next }
  G <- G[grepl("^MCP2", programF)]; if (!nrow(G)) { cat("  ", ct, ": no MCP2 rows\n"); next }
  G[, pmin := pmin(p.up, p.down)]; G[, planted := genes %in% TR$genes[[ct]]]
  cat(sprintf("  %-5s rows %4d | planted among candidates %3d | min p: planted median %.3g, others median %.3g | planted with p<0.05: %d, <0.001: %d | cols: %s\n", ct, nrow(G), sum(G$planted), median(G[planted == TRUE, pmin]), median(G[planted == FALSE, pmin]), sum(G[planted == TRUE, pmin] < 0.05), sum(G[planted == TRUE, pmin] < 0.001), paste(head(names(G), 12), collapse = ","))) }
cat("\n=== MCP2 sample-level score vs f per compartment (cca.scores, step I) ===\n")
for (ct in names(F$cca.scores)) { sc <- F$cca.scores[[ct]][, "MCP2"]; sm <- F$samples.cells[[ct]]; A <- tapply(sc, sm, mean); cat(sprintf("  %-5s cor(sample-mean MCP2 score, f) = %.2f over %d samples\n", ct, cor(A, TR$factor[names(A)]), length(A))) }
