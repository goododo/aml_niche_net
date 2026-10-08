#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_integrated/ws_rpca_integrated.rds        (42 GB, read once)
#          results/tables/ws/ws_clusters_percell.csv.gz            (res.1 assignments)
# OUTPUT : results/tables/ws/ws_marker_panel_res1.csv       57 clusters x canonical panel
#          results/tables/ws/ws_markers_target_res1.csv     unbiased DE for the 4 target clusters
#          results/tables/ws/ws_cluster_qc_res1.csv         per-cluster depth / complexity / mito
#          results/tables/ws/ws_marker_checks.csv
# what it does: answers the §0C.8 question. At resolution 1.0 (the working resolution), 13 of 57
#   clusters are >=90% one dataset. Nine are >=96% AML cells -- expected, because blasts carry
#   patient- and study-specific programs and no integration should merge them. The remaining four
#   (32, 37, 42, 45; 24,584 cells, 97% Control) are all GSE253355, and healthy marrow is the thing
#   that SHOULD be conserved across studies, so those four are the anomaly.
#   Two readouts, because they fail differently:
#     (1) a canonical panel scored on ALL cells across ALL 57 clusters -- cheap (about 90 genes),
#         and it is the decisive one. The panel deliberately includes the NON-haematopoietic niche
#         compartments (MSC/stromal, endothelial, osteoblast, adipocyte, pericyte): if the four
#         clusters light those up, GSE253355 captured a compartment the other studies did not
#         sample, which is biology and not a correction failure.
#     (2) unbiased Wilcoxon markers for each target cluster against a cluster-stratified background.
#   Plus per-cluster QC, because GSE253355 has ~2x the depth of the other Control studies, so a
#   technical island would show up as a depth/complexity outlier rather than as a cell identity.
# This script deletes nothing and does not modify the integrated object.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)
RES_COL <- "res.1"
TARGETS <- c(32L, 37L, 42L, 45L)
BG_N    <- 120000L          # cluster-stratified background for the unbiased test
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-46s expect %-12s got %-12s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- the panel. Non-haematopoietic compartments are the point, not an afterthought. --------------
PANEL <- list(
  HSPC        = c("CD34","AVP","CRHBP","HLF","SPINK2","KIT","GATA2","MPO","ELANE","AZU1","PRTN3"),
  Mono_DC     = c("LYZ","CD14","FCN1","FCGR3A","VCAN","CLEC9A","IRF8","CD1C","LILRA4"),
  T_NK        = c("CD3E","CD3D","CD8A","IL7R","CCL5","GNLY","NKG7","KLRD1"),
  B_plasma    = c("MS4A1","CD79A","CD79B","JCHAIN","MZB1","DERL3"),
  Erythroid   = c("HBB","HBA1","AHSP","ALAS2","CA1","GATA1"),
  Megakaryo   = c("PF4","PPBP","ITGA2B","GP9"),
  Mast        = c("TPSAB1","TPSB2","CPA3","MS4A2"),
  MSC_stromal = c("LEPR","CXCL12","NES","PDGFRA","PDGFRB","COL1A1","COL1A2","COL3A1","LUM","DCN","THY1","NT5E","ENG"),
  Endothelial = c("PECAM1","CDH5","VWF","EMCN","KDR","CLDN5","EGFL7"),
  Osteo       = c("BGLAP","RUNX2","ALPL","IBSP","SPP1"),
  Adipo       = c("ADIPOQ","FABP4","PLIN1","LPL"),
  Pericyte    = c("ACTA2","MYH11","RGS5","TAGLN"),
  Cycling     = c("MKI67","TOP2A","PCNA"),
  Stress      = c("FOS","JUN","EGR1","FOSB","HSPA1A","ATF3"))
GENES <- unique(unlist(PANEL))

## -- 1. load -------------------------------------------------------------------------------------
ts("[1] reading the res.1 assignments")
CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell", RES_COL))
setnames(CL, RES_COL, "cl"); setkey(CL, cell)
chk("cells in the assignment table", 895228L, nrow(CL))
chk("clusters at res 1.0", 57L, uniqueN(CL$cl))
chk("target clusters present", TRUE, all(TARGETS %in% CL$cl))

PANEL_SIDE <- file.path(BASE, "ws_panel_cache.rds")   # full metadata + the panel-gene slice
load_object <- function() {
  ts("    reading the integrated object (42 GB, about 17 minutes) ...")
  x <- readRDS(file.path(BASE, "ws_rpca_integrated.rds"))
  DefaultAssay(x) <- "RNA"
  l <- SeuratObject::Layers(x[["RNA"]])
  ts("    RNA layers: ", length(l), " (", paste(head(l, 3), collapse = ", "), ", ...)")
  if (length(grep("^counts", l)) > 1L) { ts("    joining split layers"); x <- JoinLayers(x) }
  if (!"data" %in% SeuratObject::Layers(x[["RNA"]])) { ts("    no data layer -> NormalizeData"); x <- NormalizeData(x, verbose = FALSE) }
  x
}
if (file.exists(PANEL_SIDE)) {
  ts("[2] reusing the panel cache; the 42 GB object is only loaded if step 5 needs it")
  CA <- readRDS(PANEL_SIDE); o <- NULL
} else {
  ts("[2] no panel cache yet")
  o <- load_object()
  ts("    caching metadata + the panel slice so a later failure does not repay the 17 minutes")
  gkeep <- intersect(GENES, rownames(o))
  CA <- list(meta = o@meta.data,
             panel = SeuratObject::LayerData(o, assay = "RNA", layer = "data")[gkeep, , drop = FALSE],
             cells = colnames(o))
  saveRDS(CA, PANEL_SIDE)
}
chk("cached cells match the assignment table", TRUE, setequal(CA$cells, CL$cell))
cl_vec <- CL[J(CA$cells), cl]
chk("no cell left without a cluster", 0L, sum(is.na(cl_vec)))
CA$meta$cl <- factor(cl_vec)
if (!is.null(o)) o$cl <- factor(cl_vec)

## -- 3. per-cluster QC: a technical island looks like a depth outlier, not a cell type -----------
ts("[3] per-cluster QC")
QC <- as.data.table(CA$meta)[, .(cells = .N,
      med_nCount = as.numeric(median(nCount_RNA)), med_nFeature = as.numeric(median(nFeature_RNA)),
      med_pmt = round(median(percent.mt, na.rm = TRUE), 2),
      pct_control = round(100 * mean(SampleType == "Control"), 1),
      top_dataset = names(which.max(table(Dataset_ID))),
      top_dataset_share = round(max(table(Dataset_ID)) / .N, 4)), by = cl]
QC[, is_target := cl %in% TARGETS]
fwrite(QC[order(-is_target, cl)], file.path(TAB, "ws_cluster_qc_res1.csv"))
ts("    target clusters vs the rest:")
print(QC[is_target == TRUE][order(cl)])
print(QC[, .(clusters = .N, med_nCount = median(med_nCount), med_nFeature = median(med_nFeature),
             med_pmt = median(med_pmt)), by = is_target])

## -- 4. the canonical panel, on ALL cells, all 57 clusters ---------------------------------------
ts("[4] canonical panel: ", length(GENES), " genes x 57 clusters, all 895,228 cells")
have <- rownames(CA$panel)
chk("panel genes found", length(intersect(GENES, have)), length(have))
if (length(setdiff(GENES, have))) ts("    absent from the gene axis: ", paste(setdiff(GENES, have), collapse = ", "))
dat <- CA$panel
grp <- split(seq_along(CA$cells), CA$meta$cl)
P <- rbindlist(lapply(names(grp), function(k) {
  j <- grp[[k]]; sub <- dat[, j, drop = FALSE]
  data.table(cl = as.integer(k), gene = have,
             mean_expr = round(Matrix::rowMeans(sub), 4),
             pct_expr  = round(100 * Matrix::rowSums(sub > 0) / length(j), 1))
}))
ann <- rbindlist(lapply(names(PANEL), function(k) data.table(gene = PANEL[[k]], compartment = k)))
P <- merge(P, ann, by = "gene", sort = FALSE)
fwrite(P[order(cl, compartment, -mean_expr)], file.path(TAB, "ws_marker_panel_res1.csv"))

# compartment score per cluster = mean of its genes' mean_expr, then z-scored across clusters
SC <- P[, .(score = mean(mean_expr)), by = .(cl, compartment)]
SC[, z := round((score - mean(score)) / ifelse(sd(score) > 0, sd(score), 1), 2), by = compartment]
ts("    compartment z-scores, target clusters (z > 2 means this cluster owns that compartment):")
print(dcast(SC[cl %in% TARGETS], cl ~ compartment, value.var = "z"))
ts("    the single highest-z compartment per target cluster:")
print(SC[cl %in% TARGETS][order(cl, -z)][, .SD[1:3], by = cl])

## -- 5. unbiased markers, target vs a cluster-stratified background ------------------------------
ts("[5] unbiased Wilcoxon markers vs a ", BG_N, "-cell stratified background")
if (is.null(o)) { o <- load_object(); o$cl <- factor(CL[J(colnames(o)), cl]) }
stopifnot(identical(colnames(o), CA$cells))
bg_pool <- setdiff(seq_along(CA$cells), unlist(grp[as.character(TARGETS)]))
cl_bg <- CA$meta$cl[bg_pool]
per <- max(1L, floor(BG_N / uniqueN(cl_bg)))
bg <- unlist(lapply(split(bg_pool, cl_bg), function(v) if (length(v) > per) sample(v, per) else v))
ts("    background: ", length(bg), " cells from ", uniqueN(CA$meta$cl[bg]), " clusters")
MK <- rbindlist(lapply(TARGETS, function(t) {
  j <- grp[[as.character(t)]]
  sm <- subset(o, cells = CA$cells[c(j, bg)])
  sm$grp <- ifelse(colnames(sm) %in% CA$cells[j], "target", "background")
  Idents(sm) <- "grp"
  m <- tryCatch(FindMarkers(sm, ident.1 = "target", ident.2 = "background", assay = "RNA",
                            only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE),
                error = function(e) { ts("    FindMarkers failed for ", t, ": ", conditionMessage(e)); NULL })
  rm(sm); gc(FALSE)
  if (is.null(m)) return(NULL)
  d <- as.data.table(m, keep.rownames = "gene")[p_val_adj < 0.05][order(-avg_log2FC)]
  ts(sprintf("    cluster %d: %d significant markers | top: %s", t, nrow(d),
             paste(head(d$gene, 12), collapse = ", ")))
  cbind(cl = t, head(d, 40))
}), fill = TRUE)
fwrite(MK, file.path(TAB, "ws_markers_target_res1.csv"))

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_marker_checks.csv"))
ts(sprintf("[6] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", TAB)
