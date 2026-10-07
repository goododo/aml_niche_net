#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../06_reintegration/02_integrated/ws_pca_meta.rds   (250 MB: 30-dim PCA + metadata)
# OUTPUT : results/tables/ws/ws_reso_scan.csv            per-resolution metrics + 03.1's combined_score
#          results/tables/ws/ws_clusters_percell.csv.gz  cell x resolution assignments
#          results/tables/ws/ws_cluster_composition.csv  cluster x dataset counts / share / enrichment
#          results/tables/ws/ws_control_cluster_overlap.csv  the GSE253355 question, per resolution
#          results/tables/ws/ws_cluster_checks.csv
#          /LARGE1/.../02_integrated/ws_snn_graph.rds    the SNN graph (built once, reused)
# what it does: re-runs 03.1_Diff_reso.r's resolution scan on the new integration. Grid, algorithm,
#   dims, the small-cluster threshold and the combined-score formula are all copied from that script.
#   It also answers the open question from §0C.8 -- whether GSE253355's Control cells share clusters
#   with BoneMarrowMap's Control cells, or sit in clusters of their own.
#
# Two deliberate departures from 03.1, each with its reason:
#   1. It reads the small PCA sidecar instead of the 42 GB object. FindNeighbors and FindClusters
#      accept a matrix and a Graph directly, so the counts are never needed here. Saves a 12-minute
#      read and several hundred GB of memory, and changes no number.
#   2. Stability ARI is computed on a FIXED 150,000-cell subsample, not on 0.85 x 895,228. 03.1
#      rebuilds the neighbour graph inside every stability iteration, so the original form is 36
#      full FindNeighbors+FindClusters runs on 761k cells. ARI is only ever compared ACROSS
#      resolutions here, and the subsample is identical for every resolution, so the ranking the
#      combined score depends on is preserved.
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
set.seed(42)                                    # 03.1's random_seed
# future's 500 MiB globals guard refuses FindNeighbors at this cell count (job 3641585). The
# sbatch sets R_FUTURE_GLOBALS_MAXSIZE too; this line keeps the script correct when run directly.
options(future.globals.maxSize = 64 * 1024^3)   # 64 GiB
RESO     <- seq(0.2, 2.4, by = 0.2)             # 03.1: seq(0.2, 2.4, by = 0.2)
DIMS     <- 1:30                                # 03.1: dims = 1:30
ALGO     <- 1L                                  # 03.1: algorithm = 1 (Louvain)
K_PARAM  <- 20L                                 # Seurat default, as 03.1 left it
SMALL    <- 100L                                # 03.1: small_cluster_threshold = 100
N_ITER   <- 3L                                  # 03.1: n_iterations = 3
STAB_N   <- 150000L                             # see departure 2
W_SMALL  <- 0.3; W_STAB <- 0.7                  # 03.1: calculate_combined_score weights

BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
SIDE <- file.path(BASE, "ws_pca_meta.rds")
GRPH <- file.path(BASE, "ws_snn_graph.rds")
stopifnot(requireNamespace("mclust", quietly = TRUE))   # adjustedRandIndex, as 03.1 uses
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected),
                                        got = as.character(got), pass = ok)
  message(sprintf("  %-46s expect %-12s got %-12s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
flush_checks <- function() fwrite(rbindlist(CHK), file.path(TAB, "ws_cluster_checks.csv"))

## ---- 0. input ----------------------------------------------------------------------------------
ts("[0] reading the PCA sidecar")
M <- readRDS(SIDE); emb <- M$emb[, DIMS, drop = FALSE]; md <- as.data.table(M$meta); rm(M); gc(FALSE)
chk("cells", 895228L, nrow(emb))
chk("PCA dims used", 30L, ncol(emb))
chk("metadata aligned to the embedding", TRUE, identical(rownames(emb), md$cell))
chk("datasets", 9L, uniqueN(md$Dataset_ID))

## ---- 1. neighbour graph, built once ------------------------------------------------------------
if (file.exists(GRPH)) {
  ts("[1] reusing the SNN graph on disk (nothing is overwritten)")
  snn <- readRDS(GRPH)
} else {
  ts("[1] FindNeighbors(k.param = ", K_PARAM, ", dims = 1:30) on 895,228 cells")
  snn <- FindNeighbors(emb, k.param = K_PARAM, verbose = TRUE)$snn
  ts("    writing the graph so a re-run or a finer grid skips this step: ", GRPH)
  saveRDS(snn, GRPH)
}
chk("graph is square over all cells", TRUE, all(dim(snn) == nrow(emb)))

## ---- 2. clustering across the grid -------------------------------------------------------------
ts("[2] FindClusters(algorithm = ", ALGO, ") at ", length(RESO), " resolutions: ",
   paste(RESO, collapse = ", "))
CL <- FindClusters(snn, resolution = RESO, algorithm = ALGO, random.seed = 42, verbose = TRUE)
setDT(CL); CL[, cell := rownames(emb)]
chk("one column per resolution", length(RESO), ncol(CL) - 1L)
fwrite(CL, file.path(TAB, "ws_clusters_percell.csv.gz"), compress = "gzip")
rescols <- setdiff(names(CL), "cell")

## ---- 3-4. per-resolution metrics and stability -------------------------------------------------
ts("[3] metrics + stability ARI (", N_ITER, " iterations on a fixed ", STAB_N, "-cell subsample)")
stab_idx <- sort(sample(nrow(emb), min(STAB_N, nrow(emb))))
emb_s    <- emb[stab_idx, , drop = FALSE]
ts("    building the subsample graph once (shared by every resolution)")
snn_s    <- FindNeighbors(emb_s, k.param = K_PARAM, verbose = FALSE)$snn

ROWS <- list()
for (i in seq_along(RESO)) {
  res <- RESO[i]; cc <- CL[[rescols[i]]]
  sz  <- table(cc)
  ari <- vapply(seq_len(N_ITER), function(it) {
    set.seed(42 + it)
    # the sub-clustering is itself stochastic; a different seed per iteration is what 03.1 varies
    sub <- FindClusters(snn_s, resolution = res, algorithm = ALGO, random.seed = 42 + it,
                        verbose = FALSE)[[1]]
    mclust::adjustedRandIndex(cc[stab_idx], sub)
  }, numeric(1))
  ROWS[[i]] <- data.table(resolution = res, n_clusters = length(sz),
    min_size = min(sz), median_size = as.numeric(median(sz)), max_size = max(sz),
    n_small = sum(sz < SMALL), pct_small = round(100 * sum(sz < SMALL) / length(sz), 2),
    mean_ARI = round(mean(ari), 4), sd_ARI = round(sd(ari), 4))
  ts(sprintf("    res %.1f -> %3d clusters | %2d small (%.1f%%) | size %d-%d | ARI %.3f +- %.3f",
             res, length(sz), sum(sz < SMALL), 100*sum(sz < SMALL)/length(sz),
             min(sz), max(sz), mean(ari), sd(ari)))
  fwrite(rbindlist(ROWS), file.path(TAB, "ws_reso_scan.csv"))   # incremental: a timeout keeps these
}
S <- rbindlist(ROWS)

# 03.1's combined score, verbatim: min-max normalise, invert pct_small, weight 0.3 / 0.7
nrm <- function(x) if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) rep(0.5, length(x)) else
                   (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
S[, combined_score := round(W_SMALL/(W_SMALL+W_STAB) * (1 - nrm(pct_small)) +
                            W_STAB /(W_SMALL+W_STAB) * nrm(mean_ARI), 4)]
OPT <- S$resolution[which.max(S$combined_score)]
fwrite(S, file.path(TAB, "ws_reso_scan.csv"))
ts("[4] optimal resolution by 03.1's combined score: ", OPT)
print(S[order(-combined_score)])

## ---- 5. cluster x dataset composition, every resolution ----------------------------------------
ts("[5] cluster x dataset composition")
glob <- md[, .(global = .N / nrow(md)), by = Dataset_ID]
COMP <- rbindlist(lapply(seq_along(RESO), function(i) {
  d <- data.table(cluster = CL[[rescols[i]]], Dataset_ID = md$Dataset_ID)[, .N, by = .(cluster, Dataset_ID)]
  d[, size := sum(N), by = cluster][, share := N / size]
  merge(d, glob, by = "Dataset_ID", sort = FALSE)[, .(resolution = RESO[i], cluster, Dataset_ID,
    n = N, cluster_size = size, share = round(share, 4), enrichment = round(share / global, 3))]
}))
fwrite(COMP[order(resolution, cluster, -share)], file.path(TAB, "ws_cluster_composition.csv"))
# a cluster is "dataset-private" when one dataset supplies >= 90% of it
priv <- COMP[share >= 0.90, .(resolution, cluster, Dataset_ID, cluster_size, share)]
pv <- priv[, .(private_clusters = .N, private_cells = sum(cluster_size)), by = resolution]
ts("    clusters that are >= 90% one dataset:")
print(merge(S[, .(resolution, n_clusters)], pv, by = "resolution", all.x = TRUE))

## ---- 6. the §0C.8 question: do GSE253355 and BoneMarrowMap share Control clusters? -------------
ts("[6] Control-only cluster overlap (the open question from 0C.8)")
ci <- which(md$SampleType == "Control")
OV <- rbindlist(lapply(seq_along(RESO), function(i) {
  d <- data.table(cluster = CL[[rescols[i]]][ci], Dataset_ID = md$Dataset_ID[ci])
  w <- dcast(d[, .N, by = .(cluster, Dataset_ID)], cluster ~ Dataset_ID, value.var = "N", fill = 0L)
  ds <- setdiff(names(w), "cluster"); w[, tot := rowSums(.SD), .SDcols = ds]
  # for each Control dataset: what share of ITS Control cells sits in clusters where
  # BoneMarrowMap supplies < 5% of the Control cells (i.e. clusters it effectively does not enter)
  bmm_share <- if ("BoneMarrowMap" %in% ds) w$BoneMarrowMap / w$tot else rep(NA_real_, nrow(w))
  rbindlist(lapply(ds, function(x) data.table(resolution = RESO[i], Dataset_ID = x,
    control_cells = sum(w[[x]]),
    clusters_entered = sum(w[[x]] > 0),
    cells_in_bmm_free_clusters = sum(w[[x]][bmm_share < 0.05], na.rm = TRUE),
    pct_in_bmm_free = round(100 * sum(w[[x]][bmm_share < 0.05], na.rm = TRUE) / max(1, sum(w[[x]])), 1))))
}))
fwrite(OV[order(resolution, Dataset_ID)], file.path(TAB, "ws_control_cluster_overlap.csv"))
ts("    at the optimal resolution ", OPT, ":")
print(OV[resolution == OPT])
ts("    verdict basis -- a Control dataset whose cells sit mostly in BoneMarrowMap-free clusters is",
   " either a distinct compartment (biology) or uncorrected (batch); markers decide which.")

chk("composition rows cover every resolution", length(RESO), uniqueN(COMP$resolution))
chk("overlap rows cover every resolution", length(RESO), uniqueN(OV$resolution))
C <- rbindlist(CHK); flush_checks()
ts(sprintf("[7] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] optimal resolution ", OPT, " | tables in ", TAB)
