#!/usr/bin/env Rscript
# input : ws_pb_GSE185381.rds ($cells, $raw$T), 04_annotated objects of the 25 S5 T samples (meta.data + counts colSums only),
#         ws_s5_scores.csv (sample, status columns only), ws_s5_imprint_genes.csv + ws_s5_de_T.csv (gene/direction columns only)
# output: scratch_D2/real_*.csv/rds : per-sample x fine counts; per-cell UMI and nFeature of kept T cells; status-free pooled
#         background profile of non-F402 genes (F402 values never read)
# what it does: D2 section 2 (counts only). No expression score, no imprint-gene value is computed.
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D2"
TAB <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
L1 <- "/LARGE1/gr10634/gaozy/aml_niche_net"
S5 <- fread(file.path(TAB, "ws_s5_scores.csv"), select = c("sample", "status", "compartment"))[compartment == "T"]
SAMP <- S5$sample; STAT <- setNames(S5$status, S5$sample)
cat("S5 T samples:", length(SAMP), " AML", sum(STAT == "AML"), " healthy", sum(STAT != "AML"), "\n")
IMP <- fread(file.path(TAB, "ws_s5_imprint_genes.csv"), select = c("compartment", "gene", "direction", "imprint_called"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TAB, "ws_s8_t_imprint_ambient_prone.csv"), select = "gene")$gene
TESTED <- fread(file.path(TAB, "ws_s5_de_T.csv"), select = "gene")$gene
F402 <- IMP$gene; P370 <- setdiff(F402, AMB)
cat("F402 up/down:", table(IMP$direction), " P370 up/down:", table(IMP[gene %in% P370]$direction), " tested:", length(TESTED),
    " F402 in tested:", sum(F402 %in% TESTED), " ambient-prone in F402:", sum(AMB %in% F402), "\n")
PB <- readRDS(file.path(L1, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds"))
cat("PB elements:", paste(names(PB), collapse = ","), "\n"); G <- PB$genes; cat("gene axis:", length(G), "\n")
cat("$cells columns:", paste(names(PB$cells), collapse = ","), " rows:", nrow(PB$cells), "\n")
KC <- PB$cells[kept == TRUE & group4 == "T" & sample %in% SAMP]
cat("kept T cells in the 25:", nrow(KC), "\n")
tab <- dcast(KC[, .N, by = .(sample, fine)], sample ~ fine, value.var = "N", fill = 0)
tab[, status := STAT[sample]]; tab[, n_T := rowSums(.SD), .SDcols = setdiff(names(tab), c("sample", "status"))]
fwrite(tab[order(status, -n_T)], file.path(OUT, "real_T_cells_by_sample_fine.csv"))
# status-free pooled background: non-F402 genes only (F402 rows are dropped BEFORE summing; their values are never read)
RT <- PB$raw$T; stopifnot(all(SAMP %in% colnames(RT)))
bg_genes <- setdiff(rownames(RT), F402)
bg <- Matrix::rowSums(RT[bg_genes, SAMP]); rm(RT, PB); invisible(gc(FALSE))
saveRDS(list(bg = bg, tested_nonF402 = setdiff(TESTED, F402), F402_n = length(F402), axis_n = length(G),
             set_sizes = list(P370 = table(IMP[gene %in% P370]$direction), F402 = table(IMP$direction))),
        file.path(OUT, "real_background_profile_nonF402.rds"))
# per-cell UMI / nFeature of kept T cells (QC numbers; depth only)
AN <- file.path(L1, "02_seurat_objects", "04_annotated", "GSE185381")
suppressPackageStartupMessages(library(SeuratObject))
res <- list()
for (s in SAMP) {
  o <- readRDS(file.path(AN, paste0(s, ".rds"))); md <- o@meta.data
  X <- SeuratObject::LayerData(o[["RNA"]], layer = "counts"); rm(o)
  k <- KC[sample == s]$cell; stopifnot(all(k %in% colnames(X)), identical(rownames(X), G))
  X <- X[, k, drop = FALSE]
  res[[s]] <- data.table(sample = s, status = STAT[[s]], cell = k, fine = KC[sample == s]$fine[match(k, KC[sample == s]$cell)],
                         umi_axis = as.numeric(Matrix::colSums(X)), nfeat_axis = as.numeric(Matrix::colSums(X > 0)),
                         nCount_RNA = md[k, "nCount_RNA"], nFeature_RNA = md[k, "nFeature_RNA"],
                         pct_mt = if ("percent.mt" %in% names(md)) md[k, "percent.mt"] else NA_real_)
  rm(X); invisible(gc(FALSE))
}
PC <- rbindlist(res); fwrite(PC, file.path(OUT, "real_T_percell_umi.csv.gz"))
cat("meta columns of last object:", paste(names(md), collapse = ","), "\n")
cat("nCount_RNA == colSums on axis:", PC[, mean(abs(nCount_RNA - umi_axis) < 0.5)], "\n")
q <- function(x) sprintf("median %.0f IQR %.0f-%.0f (P5 %.0f, P95 %.0f)", median(x), quantile(x, .25), quantile(x, .75), quantile(x, .05), quantile(x, .95))
for (st in c("AML", "healthy")) { x <- PC[if (st == "AML") status == "AML" else status != "AML"]
  cat(st, ": cells", nrow(x), " UMI", q(x$umi_axis), " | nFeature", q(x$nfeat_axis), "\n") }
cat("all: UMI", q(PC$umi_axis), " nFeature", q(PC$nfeat_axis), "\n")
