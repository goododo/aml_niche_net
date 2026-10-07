#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../06_reintegration/02_integrated/ws_rpca_integrated.rds   (41 GB)
#          -- or ws_pca_meta.rds, the small sidecar this script writes on its first run
# OUTPUT : results/tables/ws/ws_rpca_mixing_bytype.csv    mixing computed WITHIN each SampleType
#          results/tables/ws/ws_rpca_neighbor_mix.csv     from -> to neighbour shares + enrichment
#          /LARGE1/.../02_integrated/ws_pca_meta.rds      PCA embedding + metadata only (~250 MB)
# what it does: the overall mixing score cannot distinguish two very different situations, and
#   after integration the two lowest-scoring datasets were exactly the two healthy-only ones
#   (BoneMarrowMap 0.409, GSE253355 0.334):
#     (A) healthy cells sit apart from blasts            -> real biology, expected, fine
#     (B) healthy cells of one study sit apart from the
#         healthy cells of another study                 -> residual batch effect, a problem
#   So mixing is recomputed with the k-NN search RESTRICTED to one SampleType at a time. Inside
#   the Control cells, the only datasets a Control cell can neighbour are the other Control-bearing
#   studies, so a low score there is (B) and cannot be explained by (A).
#   The from -> to table answers the sharper question directly: whose neighbours are they? Shares
#   are reported next to the share expected under random mixing, so the number to read is the
#   enrichment, not the raw share (datasets differ in size by 15x).
# This script deletes nothing and does not modify the integrated object.

suppressPackageStartupMessages({ library(data.table); library(Seurat); library(SeuratObject) })
set.seed(1234)
K      <- 30L
SUB_N  <- 30000L        # cells per SampleType stratum, matching 05's subsample size
BASE   <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB    <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
OBJ    <- file.path(BASE, "ws_rpca_integrated.rds")
SIDE   <- file.path(BASE, "ws_pca_meta.rds")
stopifnot(requireNamespace("RANN", quietly = TRUE))
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)

CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected),
                                        got = as.character(got), pass = ok)
  message(sprintf("  %-46s expect %-12s got %-12s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## ---- load: the sidecar if it exists, else the 41 GB object once ---------------------------------
if (file.exists(SIDE)) {
  ts("reading the small sidecar (the 41 GB object is not touched)")
  M <- readRDS(SIDE)
} else {
  ts("reading the integrated object (41 GB, several minutes) ...")
  o <- readRDS(OBJ)
  M <- list(emb  = Embeddings(o, reduction = "pca")[, 1:30],
            meta = data.table(cell       = colnames(o),
                              Dataset_ID = as.character(o$Dataset_ID),
                              SampleType = as.character(o$SampleType),
                              Sample_Name= as.character(o$Sample_Name)))
  rm(o); gc(FALSE)
  ts("writing the sidecar so later diagnostics skip the 41 GB read: ", SIDE)
  saveRDS(M, SIDE)
}
emb <- M$emb; md <- M$meta
chk("cells", 895228L, nrow(md))
chk("embedding rows match metadata", TRUE, identical(rownames(emb), md$cell))
chk("datasets", 9L, uniqueN(md$Dataset_ID))
chk("SampleType levels", 2L, uniqueN(md$SampleType))

## ---- the mixing statistic -----------------------------------------------------------------------
# fraction of a cell's K neighbours drawn from a different dataset, divided by the fraction
# expected if neighbours were random draws from THIS stratum (1 - the dataset's share of it).
mix_within <- function(idx, label) {
  ds <- md$Dataset_ID[idx]
  if (length(idx) <= K + 1L || uniqueN(ds) < 2L) return(NULL)
  e  <- emb[idx, , drop = FALSE]
  nn <- RANN::nn2(e, k = K + 1L)$nn.idx[, -1, drop = FALSE]
  nbr_ds <- matrix(ds[nn], nrow = nrow(nn))
  frac   <- 1 - rowMeans(nbr_ds == ds)
  share  <- table(ds) / length(ds)
  own    <- as.numeric(share[ds])
  norm   <- ifelse(own < 1, frac / (1 - own), NA_real_)

  per_ds <- data.table(stratum = label, Dataset_ID = ds, mixing = norm
    )[, .(cells = .N, median_mixing = round(median(mixing, na.rm = TRUE), 4),
          mean_mixing = round(mean(mixing, na.rm = TRUE), 4)), by = .(stratum, Dataset_ID)]

  # from -> to: whose neighbours are they, against the random-mixing expectation
  long <- data.table(from = rep(ds, K), to = as.vector(nbr_ds))
  tt   <- long[, .N, by = .(from, to)]
  tt[, share := N / sum(N), by = from]
  exp_share <- data.table(to = names(share), expected = as.numeric(share))
  tt <- merge(tt, exp_share, by = "to", sort = FALSE)
  tt[, `:=`(stratum = label, share = round(share, 4), expected = round(expected, 4),
            enrichment = round(share / expected, 3))]
  list(per_ds = per_ds, pairs = tt[, .(stratum, from, to, n = N, share, expected, enrichment)])
}

strata <- list(all = seq_len(nrow(md)))
for (s in sort(unique(md$SampleType))) strata[[s]] <- which(md$SampleType == s)

PER <- list(); PAIR <- list()
for (nm in names(strata)) {
  idx <- strata[[nm]]
  if (length(idx) > SUB_N) idx <- sort(sample(idx, SUB_N))
  ts(sprintf("[%s] %d cells | %d datasets", nm, length(idx), uniqueN(md$Dataset_ID[idx])))
  r <- mix_within(idx, nm)
  if (is.null(r)) { message("    skipped (too few cells or a single dataset)"); next }
  PER[[nm]] <- r$per_ds; PAIR[[nm]] <- r$pairs
  print(r$per_ds[order(median_mixing)])
}

P <- rbindlist(PER); Q <- rbindlist(PAIR)
fwrite(P[order(stratum, median_mixing)], file.path(TAB, "ws_rpca_mixing_bytype.csv"))
fwrite(Q[order(stratum, from, -enrichment)], file.path(TAB, "ws_rpca_neighbor_mix.csv"))

## ---- the question this script exists to answer --------------------------------------------------
ts("verdict for the two healthy-only datasets")
for (d in c("BoneMarrowMap", "GSE253355")) {
  a <- P[stratum == "all" & Dataset_ID == d, median_mixing]
  c_ <- P[stratum == "Control" & Dataset_ID == d, median_mixing]
  self <- Q[stratum == "Control" & from == d & to == d, enrichment]
  message(sprintf("  %-14s all-cell %.3f | within-Control %.3f | self-neighbour enrichment %.2fx -> %s",
    d, if (length(a)) a else NA, if (length(c_)) c_ else NA, if (length(self)) self else NA,
    if (length(c_) && !is.na(c_) && c_ >= 0.5) "mixes with the other healthy studies (biology, not batch)"
    else "still isolated among healthy cells -> residual batch effect"))
}
# a Control cell's single most enriched neighbour source, per dataset
ts("most enriched neighbour source within Control")
print(Q[stratum == "Control"][order(from, -enrichment)][, .SD[1:2], by = from])

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_mixing_diag_checks.csv"))
ts(sprintf("%d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(TAB, "ws_rpca_mixing_bytype.csv"))
