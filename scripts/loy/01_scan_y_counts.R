#!/usr/bin/env Rscript
# 01_scan_y_counts.R ----
# INPUT  : 02_seurat_objects/04_annotated/<dataset>/<sample>.rds (all objects of one dataset)
# OUTPUT : LARGE1/loy/percell/<dataset>__y_percell.csv.gz   one row per cell: metadata, nCount, UMIs per Y gene,
#            y_umi (sum over Y genes), xist_umi
#          results/tables/loy/01_scan/<dataset>__y_sample_summary.csv   one row per sample x hierarchy_bin x
#            depth bin: n cells, n with y_umi == 0 / == 1 / >= 2, median nCount, pooled Y rate, XIST-positive frac
#          results/tables/loy/01_scan/<dataset>__y_scan_checks.csv     known answers (see below)
# what it does: LOY step 0 (counts only, no calls). For every cell it records how many UMIs fall on the
#   chrY-linked genes of config_loy.R and on XIST. The three-state call needs the zero-Y background per
#   chemistry, compartment and depth, which can only be estimated once all datasets are scanned, so this
#   script deliberately stops at counting. Doublets were already removed upstream (union rule, 03_per_sample_qc),
#   which helps here: a Y+/Y- doublet would read as Y retained.
# Usage  : Rscript scripts/loy/01_scan_y_counts.R --dataset=GSE185381
#          Rscript scripts/loy/01_scan_y_counts.R --task=3     (3rd line of results/tables/loy/01_scan/datasets.txt)
# Checks : K1 every object in the directory produced rows (n files == n samples in output);
#          K2 the number of Y genes found on the gene axis is printed and >= 8 for 10x datasets;
#          K3 row count of the per-cell table equals the sum of ncol over objects.
# WHY counts and not calls: a cell with 1,500 UMIs is Y-zero by chance about 10% of the time even with an
#   intact Y (probe 2026-11-13, Control0005). Calling LOY per cell without the background would manufacture
#   LOY in low-depth compartments and in Seq-Well data. The background is estimated in 02 on cells whose
#   Y is known to be intact.

suppressPackageStartupMessages({ library(here); library(Seurat); library(Matrix); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
set.seed(SEED)

args <- commandArgs(TRUE)
getarg <- function(key) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v[1] else NA }
ds <- getarg("dataset"); task <- getarg("task")
if (is.na(ds) && !is.na(task)) {
  dsl <- readLines(file.path(LOY_SCAN_DIR, "datasets.txt")); ds <- dsl[as.integer(task)]
}
stopifnot(!is.na(ds), nzchar(ds))
dir.create(LOY_PERCELL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOY_SCAN_DIR,    showWarnings = FALSE, recursive = TRUE)

in_dir <- file.path(LOY_ANNOT_DIR, ds)
files  <- sort(list.files(in_dir, pattern = "\\.rds$", full.names = TRUE))
message(sprintf("[1] %s: %d annotated objects in %s", ds, length(files), in_dir))
stopifnot(length(files) > 0)

out_cells <- file.path(LOY_PERCELL_DIR, paste0(ds, "__y_percell.csv.gz"))
if (file.exists(out_cells) && is.na(getarg("force"))) { message("[skip] exists: ", out_cells); quit(save = "no") }

rows <- vector("list", length(files)); ncols <- integer(length(files)); ygenes_found <- NULL
for (i in seq_along(files)) {
  f <- files[i]; sm <- sub("\\.rds$", "", basename(f))
  s <- readRDS(f)
  m <- tryCatch(LayerData(s, assay = "RNA", layer = "counts"), error = function(e) GetAssayData(s, assay = "RNA", slot = "counts"))
  ncols[i] <- ncol(m)
  yg <- intersect(LOY_YGENES, rownames(m)); if (is.null(ygenes_found)) ygenes_found <- yg
  ym <- m[yg, , drop = FALSE]                                     # genes x cells, sparse
  xist <- if (LOY_XIST %in% rownames(m)) as.numeric(m[LOY_XIST, ]) else rep(NA_real_, ncol(m))
  md <- s@meta.data
  meta <- lapply(LOY_META_COLS, function(cn) if (cn %in% colnames(md)) md[[cn]] else rep(NA, ncol(m)))
  names(meta) <- LOY_META_COLS
  dt <- as.data.table(meta)
  dt[, dataset := ds]; dt[, object := sm]; dt[, cell := colnames(m)]
  if (all(is.na(dt$nCount_RNA))) dt[, nCount_RNA := Matrix::colSums(m)]
  ycounts <- as.data.table(as.matrix(Matrix::t(ym)))               # cells x Y genes (dense: <= 13 columns)
  setnames(ycounts, paste0("y_", colnames(ycounts)))
  dt <- cbind(dt, ycounts)
  dt[, y_umi := rowSums(as.matrix(ycounts))]
  dt[, xist_umi := xist]
  rows[[i]] <- dt
  message(sprintf("    %-32s cells %6d  Ygenes %2d  median nCount %6.0f  frac y_umi==0 %.3f  XIST>0 %.3f",
                  sm, ncol(m), length(yg), median(dt$nCount_RNA), mean(dt$y_umi == 0), mean(xist > 0, na.rm = TRUE)))
  rm(s, m, ym); invisible(gc())
}
cells <- rbindlist(rows, fill = TRUE)
setcolorder(cells, c("dataset", "object", "cell"))
fwrite(cells, out_cells)
message("[2] wrote ", out_cells, " rows ", nrow(cells))

## -- per sample x bin x depth summary (reporting only) ---------------------------------------------
cells[, depth_bin := cut(nCount_RNA, LOY_DEPTH_BREAKS, right = FALSE)]
cells[, bin := fifelse(is.na(hierarchy_bin), "NA", as.character(hierarchy_bin))]
summ <- cells[, .(n = .N, n_y0 = sum(y_umi == 0), n_y1 = sum(y_umi == 1), n_y2p = sum(y_umi >= 2),
                  med_ncount = median(nCount_RNA), sum_y = sum(y_umi), sum_ncount = sum(nCount_RNA),
                  xist_pos_frac = mean(xist_umi > 0, na.rm = TRUE), n_high_error = sum(high_error %in% c(1, TRUE))),
              by = .(dataset, object, Sample, Library, Chemistry, Platform, Timepoint, Disease_state, bin, depth_bin)]
summ[, y_rate := sum_y / sum_ncount]
fwrite(summ, file.path(LOY_SCAN_DIR, paste0(ds, "__y_sample_summary.csv")))

## -- known answers --------------------------------------------------------------------------------
chk <- data.table(check = c("K1_objects_with_rows", "K2_ygenes_found", "K3_rows_equal_ncol_sum"),
                  expected = c(length(files), ">=8 (10x) / report (Seq-Well)", sum(ncols)),
                  observed = c(uniqueN(cells$object), length(ygenes_found), nrow(cells)))
chk[, pass := c(observed[1] == expected[1], TRUE, observed[3] == expected[3])]
chk[2, pass := length(ygenes_found) >= 8 | grepl("SeqWell|Seq-Well", paste(unique(cells$Platform), collapse = ","))]
fwrite(chk, file.path(LOY_SCAN_DIR, paste0(ds, "__y_scan_checks.csv")))
print(chk)
message("[done] ", ds, " | Y genes on axis: ", paste(ygenes_found, collapse = ","))
