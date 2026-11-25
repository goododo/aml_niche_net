# 05d_h2_truncal_marker_check.R
# input : results/tables/loy/04_h2/stage0_cells_<case>.csv.gz (cell groups), the Lambo 2023 Seurat object (Zeng pack)
# output: results/tables/loy/04_h2/truncal_marker_check.csv   per case x group: cells, median depth, detection fraction, mean UMI/10k
# what it does: post-hoc check (not pre-registered) that the -Y branch of the two t(8;21) cases sits downstream of the
#   truncal translocation. RUNX1T1 / POU4F1 are t(8;21)-target genes; if -Y and Y+ leukaemic cells both express them and
#   same-sample T cells (negative control) do not, -Y is a subclone of the RUNX1::RUNX1T1 clone, not an independent clone.
#   CD3E = positive control for T identity.

suppressPackageStartupMessages({ library(data.table); library(here); library(Matrix); library(Seurat); library(SeuratObject) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "loy", "config_loy.R")); source(here::here("scripts", "loy", "00_loy_common.R"))
GENES <- c("RUNX1T1", "POU4F1", "CD3E"); CASES <- c("AML15", "AML13")

stem  <- LOY_EXT_STUDIES[dataset == "ZENG_Lambo_2023", file_stem]
s_full <- readRDS(file.path(LOY_ZENG_RDS_DIR, paste0(stem, "_projected_cleaned.rds")))
res <- rbindlist(lapply(CASES, function(cs) {
  cl <- unique(fread(file.path(LOY_H2_DIR, sprintf("stage0_cells_%s.csv.gz", cs)))[, .(cell, arm, group, call, comp, nCount_RNA)])
  cl[, grp := fcase(arm == "leukaemic" & comp == "LEUK" & group == "A", "LEUK_minusY",
                    arm == "leukaemic" & comp == "LEUK" & group == "B", "LEUK_Yplus",
                    arm == "control"   & comp == "T"    & group == "B", "T_control", default = NA_character_)]
  cl <- unique(cl[!is.na(grp)], by = "cell")
  cnt <- counts_of(subset(s_full, cells = intersect(colnames(s_full), cl$cell)))
  stopifnot(all(GENES %in% rownames(cnt)))                                  # sanity: markers present on the gene axis
  lib <- Matrix::colSums(cnt); g <- cnt[GENES, , drop = FALSE]
  cl <- cl[cell %in% colnames(cnt)]
  rbindlist(lapply(split(cl$cell, cl$grp), function(cc) data.table(
    n = length(cc), median_depth = median(lib[cc]), gene = GENES,
    det_frac = round(Matrix::rowMeans(g[, cc, drop = FALSE] > 0), 3),
    umi_per10k = round(Matrix::rowMeans(sweep(g[, cc, drop = FALSE], 2, lib[cc] / 1e4, "/")), 3))), idcol = "group")[, case := cs]
}))
fwrite(res, file.path(LOY_H2_DIR, "truncal_marker_check.csv"))
print(dcast(res, case + gene ~ group, value.var = "det_frac")); print(unique(res[, .(case, group, n, median_depth)]))
