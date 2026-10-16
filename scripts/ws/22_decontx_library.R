#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds   (RNA counts; metadata GSM, Library,
#          hierarchy_bin, bmm_fine, high_error)
#          scripts/ws/config_ws.R (WS_EXCLUDE_FINE; SEED via config_paths.R)
# OUTPUT : results/tables/ws/decontx/<dataset>/ws_decontx_index.csv.gz     one row per cell: sample, unit (GSM), group
#          results/tables/ws/decontx/<dataset>/ws_decontx_units.csv         one row per unit: cells, samples (array list)
#          results/tables/ws/decontx/<dataset>/<unit>__percell.csv.gz      cell, sample, group, contamination, UMIs before/after
#          results/tables/ws/decontx/<dataset>/<unit>__ambient_top100.csv  the 100 genes with the largest removed share
#          results/tables/ws/decontx/<dataset>/<unit>__checks.csv
#          /LARGE1/.../aml_niche_net/ws_decontx/<dataset>/<unit>.rds        decontXcounts (sparse), contamination, removed share
# what it does: v2.0 S1. Ambient RNA is shared by every droplet of one 10x channel, so the unit is the channel
#   (the per-cell GSM in the annotated objects; the manifest's pool id is coarser, and most GSE185381 donors
#   span two channels). For one unit it gathers the cells of EVERY sample that carries that GSM (AML and healthy
#   donors together), runs celda::decontX with coarse cell groups (hierarchy_bin; BoneMarrowMap Early Lymphoid
#   states as their own group, see WS_EXCLUDE_FINE), and writes the per-cell contamination fraction, the
#   decontaminated counts and the unit's ambient profile (removed UMIs per gene). The S0 decision (2026-10-16):
#   the main analyses keep the raw counts and use the contamination fraction as a covariate; the decontaminated
#   counts serve sensitivity analyses only.
#   Cell groups (--zmode): "auto" (default) lets decontX cluster the channel itself; "bin" uses hierarchy_bin with
#   the Early Lymphoid states apart; "author" uses the depositors' labels where present. Decided 2026-10-16 after
#   a three-way comparison on GSM5613781 (AML009 + Control0005): with "bin", decontX called AML009's "B_Plasma"
#   cells 95% contaminated because they are author-labelled MPP/LymP/GMP, i.e. label error read as ambient RNA;
#   "auto" gives 4% there and agrees with the author grouping elsewhere (per-cell Spearman 0.75), and it is the
#   only grouping available in every dataset (GSE289435 has no depositor labels). Sample medians barely move
#   (0.03-0.04 under all three); the per-cell values and the ambient profile do.
#   usage: --dataset=GSE185381 --index          build the cell index and unit list (reads every object once)
#          --dataset=GSE185381 --unit=GSM...    run one unit;  --task=N takes the N-th unit of the unit list
#          --zmode=auto|bin|author              cell grouping handed to decontX (default auto)
#          --smoke                              subsample to 1,500 cells and write *_smoke outputs
#          --force                              recompute an existing unit
# This script deletes nothing and does not modify the 04_annotated objects.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
DATASET <- .get("dataset", "GSE185381"); UNIT <- .get("unit", ""); TASK <- as.integer(.get("task", "0"))
INDEX <- "--index" %in% .a; SMOKE <- "--smoke" %in% .a; FORCE <- "--force" %in% .a
ZMODE <- match.arg(.get("zmode", "auto"), c("auto", "bin", "author"))

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
set.seed(SEED)
AN    <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
OBJD  <- file.path(LARGE1_DIR, "ws_decontx", DATASET); dir.create(OBJD, showWarnings = FALSE, recursive = TRUE)
TABD  <- file.path(TAB_DIR, "ws", "decontx", DATASET); dir.create(TABD, showWarnings = FALSE, recursive = TRUE)
IDX   <- file.path(TABD, "ws_decontx_index.csv.gz"); UNITS <- file.path(TABD, "ws_decontx_units.csv")
# BoneMarrowMap fine states under the broad label "Early Lymphoid" (= WS_EXCLUDE_FINE): grouped apart because in
# AML samples they are largely blasts, so they must not shape the T_NK or B_Plasma native profile
EARLY_LYMPHOID_FINE <- c("CLP", "MLP", "MLP-II", "Pre-ProB")
MIN_GROUP <- 10L            # a group below this is folded into "unassigned" (decontX needs a profile per group)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
stopifnot(identical(WS_EXCLUDE_FINE, "Early Lymphoid"))     # the grouping below assumes this one label

## -- 0. index: cell -> sample, unit, group  (one pass over every object) --------------------------
if (INDEX || !file.exists(IDX)) {
  ts("[0] building the cell index for ", DATASET)
  fs <- list.files(AN, pattern = "\\.rds$", full.names = TRUE); stopifnot(length(fs) > 0)
  IX <- rbindlist(lapply(fs, function(f) {
    o <- readRDS(f); md <- o@meta.data; sm <- sub("\\.rds$", "", basename(f))
    unit <- if ("GSM" %in% names(md)) as.character(md$GSM) else sm        # no channel id -> the sample is the unit
    ak <- intersect(c("anno_Broad_cell_identity", "anno_celltype", "anno_CellType"), names(md))[1]   # depositor label, if any
    data.table(sample = sm, cell = colnames(o), unit = unit,
               library = if ("Library" %in% names(md)) as.character(md$Library) else NA_character_,
               author = if (is.na(ak)) NA_character_ else as.character(md[[ak]]),
               hierarchy_bin = as.character(md$hierarchy_bin), bmm_fine = as.character(md$bmm_fine),
               high_error = as.integer(md$high_error), nCount = as.numeric(md$nCount_RNA))
  }))
  IX[, group := fifelse(is.na(hierarchy_bin) | hierarchy_bin == "", "unassigned", hierarchy_bin)]
  IX[bmm_fine %in% EARLY_LYMPHOID_FINE, group := "EarlyLymphoid"]
  chk("every cell has a unit", 0L, IX[is.na(unit) | unit == "", .N])
  chk("cell ids unique within dataset", nrow(IX), uniqueN(IX$cell))
  fwrite(IX, IDX, compress = "gzip")
  U <- IX[, .(n_cells = .N, n_samples = uniqueN(sample), samples = paste(sort(unique(sample)), collapse = ";"),
              n_early_lymphoid = sum(group == "EarlyLymphoid")), by = unit][order(unit)]
  fwrite(U, UNITS)
  ts(sprintf("    %d cells, %d samples, %d units; EarlyLymphoid cells %d (%.1f%%)", nrow(IX), uniqueN(IX$sample), nrow(U),
             sum(IX$group == "EarlyLymphoid"), 100 * mean(IX$group == "EarlyLymphoid")))
  print(head(U[, .(unit, n_cells, n_samples)], 8))
  if (INDEX && UNIT == "" && TASK == 0L) { fwrite(rbindlist(CHK), file.path(TABD, "ws_decontx_index_checks.csv")); quit(save = "no") }
}
IX <- fread(IDX); U <- fread(UNITS)
if (TASK > 0L) UNIT <- U$unit[TASK]
if (UNIT == "") stop("give --unit=<GSM> or --task=<n> (1..", nrow(U), "), or --index")
stopifnot(UNIT %in% U$unit)
SFX <- paste0(if (ZMODE != "auto") paste0("_", ZMODE) else "", if (SMOKE) "_smoke" else "")
OUT_RDS <- file.path(OBJD, paste0(UNIT, SFX, ".rds"))
if (file.exists(OUT_RDS) && !FORCE) { ts("[skip] exists: ", OUT_RDS); quit(save = "no") }

## -- 1. gather the unit's cells from every sample that carries it ----------------------------------
ts(sprintf("[1] unit %s: %d cells from %d sample(s)", UNIT, U[unit == UNIT, n_cells], U[unit == UNIT, n_samples]))
cells <- IX[unit == UNIT]
if (SMOKE && nrow(cells) > 1500L) { set.seed(SEED); cells <- cells[sample.int(.N, 1500L)] }
mats <- lapply(unique(cells$sample), function(sm) {
  o <- readRDS(file.path(AN, paste0(sm, ".rds")))
  m <- LayerData(o, assay = "RNA", layer = "counts")[, cells[sample == sm, cell], drop = FALSE]
  rm(o); m
})
chk("one gene axis across the unit's samples", TRUE, all(vapply(mats, function(m) identical(rownames(m), rownames(mats[[1]])), TRUE)))
X <- do.call(cbind, mats); rm(mats); gc(FALSE)
X <- X[, cells$cell, drop = FALSE]                          # same order as the index rows
chk("cells gathered", nrow(cells), ncol(X))
chk("counts are integers", TRUE, all(X@x == round(X@x)))
grp <- cells$group                                           # hierarchy_bin groups, kept in the output for reading
zin <- switch(ZMODE, auto = NULL, bin = grp,
              author = { a <- cells$author; a[is.na(a) | a == ""] <- "unassigned"; if (all(a == "unassigned")) stop("no depositor labels in ", DATASET); a })
if (!is.null(zin)) { small <- names(which(table(zin) < MIN_GROUP)); zin[zin %in% small] <- "unassigned"
  if (length(small)) ts("    groups folded into unassigned (< ", MIN_GROUP, " cells): ", paste(small, collapse = ", "))
  ts("    groups handed to decontX (", ZMODE, "): ", paste(sprintf("%s=%d", names(table(zin)), table(zin)), collapse = " ")) }
ts("    hierarchy_bin groups: ", paste(sprintf("%s=%d", names(table(grp)), table(grp)), collapse = " "))

## -- 2. decontX ------------------------------------------------------------------------------------
ts("[2] decontX (zmode = ", ZMODE, ")")
t0 <- Sys.time()
res <- celda::decontX(x = X, z = zin, seed = SEED, verbose = FALSE)
ts(sprintf("    done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
contam <- as.numeric(res$contamination); DX <- res$decontXcounts; zout <- as.character(res$z)
if (ZMODE == "auto") ts("    decontX clusters: ", length(unique(zout)))
chk("contamination within [0, 1]", TRUE, all(contam >= 0 & contam <= 1))
chk("decontXcounts has the input shape", TRUE, identical(dim(DX), dim(X)))
before <- Matrix::colSums(X); after <- Matrix::colSums(DX)
chk("no cell gains UMIs", TRUE, all(after <= before + 1e-6))
removed <- Matrix::rowSums(X) - Matrix::rowSums(DX)          # UMIs attributed to ambient RNA, per gene
rem_share <- removed / sum(removed)
top <- data.table(gene = names(sort(rem_share, decreasing = TRUE))[1:100],
                  removed_share = round(sort(rem_share, decreasing = TRUE)[1:100], 5))
top[, removed_umi := round(removed[gene], 1)]
nonery <- grp != "Erythroid"; hb <- intersect(c("HBB", "HBA1", "HBA2"), rownames(X))
if (length(hb) && any(nonery)) {
  # informational, not a check: ambient haemoglobin should be removed at least as readily as the average UMI.
  # (An earlier version compared SHARES, which can rise when other genes lose more UMIs; it failed on 25 channels
  # with HB at 1e-4 and was withdrawn on 2026-10-16. The per-unit outputs are unaffected.)
  hb_rm <- 1 - sum(DX[hb, nonery]) / max(sum(X[hb, nonery]), 1); all_rm <- 1 - sum(after[nonery]) / sum(before[nonery])
  ts(sprintf("    non-erythroid cells: fraction of haemoglobin UMIs removed %.3f (HB UMIs %d) vs all UMIs removed %.3f",
             hb_rm, as.integer(sum(X[hb, nonery])), all_rm))
}
ts("    top ambient genes: ", paste(head(top$gene, 12), collapse = " "))

## -- 3. write -------------------------------------------------------------------------------------
PC <- data.table(cell = cells$cell, sample = cells$sample, unit = UNIT, group = grp, decontx_z = zout, zmode = ZMODE,
                 contamination = round(contam, 5), umi_before = before, umi_after = round(after, 2))
fwrite(PC, file.path(TABD, paste0(UNIT, SFX, "__percell.csv.gz")), compress = "gzip")
fwrite(top, file.path(TABD, paste0(UNIT, SFX, "__ambient_top100.csv")))
saveRDS(list(dataset = DATASET, unit = UNIT, samples = unique(cells$sample), cells = cells$cell, group = grp, zmode = ZMODE,
             decontx_z = zout, contamination = contam, decontXcounts = DX, removed_share = rem_share, run = res$runParams,
             date = format(Sys.Date()), smoke = SMOKE), OUT_RDS)
print(PC[, .(cells = .N, median_contamination = round(median(contamination), 3), q75 = round(quantile(contamination, .75), 3)), by = .(sample, group)][order(sample, -cells)][1:min(.N, 12)])
C <- rbindlist(CHK); fwrite(C, file.path(TABD, paste0(UNIT, SFX, "__checks.csv")))
ts(sprintf("[4] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", OUT_RDS)
