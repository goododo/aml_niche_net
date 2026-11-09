#!/usr/bin/env Rscript
# input : Savage 2021 Table S4 (mmc5.xlsx) and the frozen T list (reference files only; no project data)
# output: 02_savage_rule.out (stdout)
# what it does: tries readings of the PROVENANCE T rule to see which one reproduces the frozen 792 (295/497), then
#   applies the matching reading to B, NK and Mono (draft 2.3 / K9). Gene names only.
suppressPackageStartupMessages({ library(readxl); library(data.table) })
X <- "/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/mmc5.xlsx"
rd <- function(s, e) { d <- as.data.table(suppressMessages(read_excel(X, sheet = s, skip = 1)))
  setnames(d, 1:7, c("ct", "smp", "tp", "gene", "lfc", "p", "padj")); d[, exp := e]
  if (e == 2) d[, lfc := -lfc]; d[, .(exp, ct, smp, tp, gene, lfc, padj)] }
D <- rbind(rd("Experiment 1", 1), rd("Experiment 2", 2))
D <- D[!is.na(padj) & !is.na(lfc)]
FZ <- fread("/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv")
cat("frozen:", nrow(FZ), "up", FZ[direction == "up", .N], "down", FZ[direction == "down", .N], "\n")
build <- function(types, mode = c("any", "same"), opp = c("18h_healthy", "all")) {
  mode <- match.arg(mode); opp <- match.arg(opp)
  h <- D[ct %in% types & grepl("Healthy", smp)]
  s18 <- h[tp == 18 & padj < 0.05]
  if (mode == "any") {
    a <- s18[, .(dir = unique(sign(lfc))), by = .(exp, gene)]
    both <- a[, .(ne = uniqueN(exp), nd = uniqueN(dir), dir = dir[1]), by = gene][ne == 2 & nd == 1]
  } else {
    a <- s18[, .(dir = unique(sign(lfc))), by = .(exp, ct, gene)]
    w <- dcast(a[, .(dir = if (uniqueN(dir) == 1) dir[1] else NA_real_), by = .(exp, ct, gene)], ct + gene ~ exp, value.var = "dir")
    setnames(w, c("1", "2"), c("e1", "e2"))
    both <- w[!is.na(e1) & !is.na(e2) & e1 == e2, .(dir = unique(e1)), by = gene][, .(nd = .N, dir = dir[1]), by = gene][nd == 1]
  }
  pool <- if (opp == "18h_healthy") s18 else D[ct %in% types & padj < 0.05]
  oppg <- merge(pool[, .(gene, sg = sign(lfc))], both[, .(gene, dir)], by = "gene")[sg != dir, unique(gene)]
  both[!gene %in% oppg, .(gene, direction = ifelse(dir > 0, "up", "down"))]
}
cmp <- function(L, tag) { j <- merge(L, FZ, by = "gene", all = TRUE, suffixes = c("_new", "_frozen"))
  cat(sprintf("%-40s n %4d (up %3d / down %3d) | identical to frozen: %s | only new %d, only frozen %d, dir mismatch %d\n", tag, nrow(L),
      L[direction == "up", .N], L[direction == "down", .N], identical(sort(paste(L$gene, L$direction)), sort(paste(FZ$gene, FZ$direction))),
      j[is.na(direction_frozen), .N], j[is.na(direction_new), .N], j[!is.na(direction_new) & !is.na(direction_frozen) & direction_new != direction_frozen, .N])) }
TT <- c("CD4_Naive", "CD4_Memory", "CD8_Naive", "CD8_Memory")
for (m in c("any", "same")) for (o in c("18h_healthy", "all")) cmp(build(TT, m, o), sprintf("T rule: %s-type, opposite in %s", m, o))
cat("\nper-compartment lists under each reading:\n")
CMP <- list(B = c("B_Naive", "B_Memory"), NK = c("NK_dim", "NK_bright"), Mono = c("CD14+_Monocytes", "CD16+_Monocytes"), T = TT)
for (m in c("any", "same")) for (o in c("18h_healthy", "all")) for (k in names(CMP)) { L <- build(CMP[[k]], m, o)
  cat(sprintf("%-5s %-5s %-12s n %4d up %3d down %3d | shared with frozen T list %4d (same direction %4d)\n", k, m, o, nrow(L), L[direction == "up", .N],
      L[direction == "down", .N], L[gene %in% FZ$gene, .N], merge(L, FZ, by = "gene")[direction.x == direction.y, .N])) }
