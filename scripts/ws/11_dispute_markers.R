#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_integrated/ws_rpca_integrated.rds   (44 GB, read once)
#          results/tables/ws/ws_anno_consolidated_res1.csv    (provenance per cluster)
#          results/tables/ws/ws_clusters_percell.csv.gz       (res.1)
# OUTPUT : results/tables/ws/ws_dispute_markers.csv        top unbiased markers per target cluster
#          results/tables/ws/ws_dispute_marker_checks.csv
# what it does: unbiased Wilcoxon markers for the clusters whose label is not settled -- the 19
#   disputed ones plus the 3 "adopted" clusters whose automatic call the independent panel
#   contradicts (27, 30, 49). The dossier so far only contains scores over the 242 marker-table
#   genes, which cannot say anything about a cluster whose identity lies outside that table; three
#   of the disputed clusters are exactly that case (e.g. cluster 13, where neither method's call is
#   supported and the panel's best compartment is Cycling at z = 2.24).
#   Each target is tested against a cluster-stratified background so no single large cluster
#   dominates the comparison.
# This script deletes nothing and does not modify the integrated object.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)
BG_N <- 120000L
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-50s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

ts("[1] targets")
CONS <- fread(file.path(TAB, "ws_anno_consolidated_res1.csv"))
TARGETS <- sort(CONS[provenance == "disputed" |
                     (provenance == "adopted" & independent_check == "contradicted"), cl])
chk("target clusters", 22L, length(TARGETS))
ts("    ", paste(TARGETS, collapse = ", "))

CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell", "res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
chk("cells", 895228L, nrow(CL))

ts("[2] reading the integrated object (44 GB, about 17 min)")
o <- readRDS(file.path(BASE, "ws_rpca_integrated.rds"))
DefaultAssay(o) <- "RNA"
if (length(grep("^counts", SeuratObject::Layers(o[["RNA"]]))) > 1L) { ts("    joining layers"); o <- JoinLayers(o) }
if (!"data" %in% SeuratObject::Layers(o[["RNA"]])) { ts("    NormalizeData"); o <- NormalizeData(o, verbose = FALSE) }
cl_vec <- CL[J(colnames(o)), cl]
chk("every cell has a cluster", 0L, sum(is.na(cl_vec)))
o$cl <- factor(cl_vec)
grp <- split(seq_len(ncol(o)), o$cl)

ts("[3] stratified background")
bg_pool <- setdiff(seq_len(ncol(o)), unlist(grp[as.character(TARGETS)]))
per <- max(1L, floor(BG_N / uniqueN(o$cl[bg_pool])))
bg  <- unlist(lapply(split(bg_pool, o$cl[bg_pool]), function(v) if (length(v) > per) sample(v, per) else v))
ts("    ", length(bg), " cells from ", uniqueN(o$cl[bg]), " clusters")

ts("[4] unbiased markers, one target at a time (written incrementally)")
OUT <- list()
for (t in TARGETS) {
  j  <- grp[[as.character(t)]]
  sm <- subset(o, cells = colnames(o)[c(j, bg)])
  sm$grp <- ifelse(colnames(sm) %in% colnames(o)[j], "target", "background")
  Idents(sm) <- "grp"
  m <- tryCatch(FindMarkers(sm, ident.1 = "target", ident.2 = "background", assay = "RNA",
                            only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE),
                error = function(e) { ts("    FindMarkers failed for ", t, ": ", conditionMessage(e)); NULL })
  rm(sm); gc(FALSE)
  if (!is.null(m)) {
    d <- as.data.table(m, keep.rownames = "gene")[p_val_adj < 0.05][order(-avg_log2FC)]
    OUT[[length(OUT)+1]] <- cbind(cl = t, n_sig = nrow(d), head(d, 40))
    ts(sprintf("    cluster %-3d %6d cells | %5d sig | top: %s", t, length(j), nrow(d),
               paste(head(d$gene, 12), collapse = ", ")))
    fwrite(rbindlist(OUT), file.path(TAB, "ws_dispute_markers.csv"))   # incremental
  }
}
M <- rbindlist(OUT)
chk("targets with markers", length(TARGETS), uniqueN(M$cl))
C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_dispute_marker_checks.csv"))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(TAB, "ws_dispute_markers.csv"))
