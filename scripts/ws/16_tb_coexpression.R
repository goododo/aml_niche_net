#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/<ds>/<sample>.rds   (the WP1 cohort)
#          results/tables/01_preprocess/00_curated_manifest.csv           (is_multiplexed, n_donors)
# OUTPUT : results/tables/ws/ws_tb_coexpression_check.csv   per sample, the T/B co-expression counts
#          results/tables/ws/ws_tb_coexpression_summary.csv per dataset x pooled status
#          results/tables/ws/ws_tb_checks.csv
# what it does: T3b. Asks, per dataset, whether the authors had already removed doublets before
#   the data was deposited -- a T cell and a B cell captured in one droplet co-express markers that
#   no single cell expresses, so a dataset still carrying its doublets shows an excess of
#   T+B+ cells over what independence predicts.
#   Run on the SAME cohort and the SAME per-sample objects as T3d (04_annotated), so the two are
#   comparable; T3d's five studies are GSE185381, GSE239721, GSE116256, GSE227903, GSE289435.
#
#   Two specification points from the handoff, both load-bearing:
#   * judgement uses NON-POOLED samples only, and pooled libraries are reported separately. This
#     matters almost entirely for GSE185381: 40 of its 52 samples are hashing_HTO multiplexed with
#     2-5 donors, so the verdict for the dataset that serves as T3c's positive control rests on
#     12 samples. In a pooled library a T+B+ droplet may be two cells from DIFFERENT donors, which
#     demultiplexing is supposed to catch, so the same number means something different there.
#   * the RATIO is never read alone. Enrichment obs/(p1*p2) explodes when both markers are rare --
#     two niche clusters once scored 53x off a single double-positive cell (§0C.12d) -- so every
#     row carries the absolute double-positive count and the ratio is suppressed below 10.
#
#   Expected magnitude, so the result can be read: a T/B doublet needs one T and one B in the same
#   droplet, so its frequency is about 2 * P(T) * P(B) * doublet_rate. At P(T) = 0.15, P(B) = 0.05
#   and a 5% doublet rate that is 0.075% of cells. The test is therefore weak by construction, and
#   a null result is NOT evidence that doublets were removed -- it is mostly evidence that this
#   particular contrast has little power. That caveat is reported with the numbers.
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
set.seed(42)
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
AN   <- "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/04_annotated"
TAB  <- file.path(ROOT, "results/tables/ws")
STUDIES <- c("GSE185381", "GSE239721", "GSE116256", "GSE227903", "GSE289435")   # == T3d's cohort
T_MK <- c("CD3E","CD3D","TRAC","CD2")
B_MK <- c("MS4A1","CD79A","CD79B")
MIN_BOTH <- 10L            # below this the ratio is noise, not a measurement
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-52s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

ts("[1] cohort and multiplexing status")
cm <- fread(file.path(ROOT, "results/tables/01_preprocess/00_curated_manifest.csv"),
            select = c("dataset","sample","is_multiplexed","n_donors_in_library","demux_method"))
chk("multiplexing recorded for every manifest row", 0L, sum(is.na(cm$is_multiplexed)))
cm[, pooled := is_multiplexed %chin% c(TRUE, "TRUE", "yes", "Y", "1")]
ts("    pooled samples per study:")
print(cm[dataset %chin% STUDIES, .(samples = .N, pooled = sum(pooled)), by = dataset][order(dataset)])

ts("[2] per sample")
OUT <- list()
for (ds in STUDIES) {
  fs <- list.files(file.path(AN, ds), "\\.rds$", full.names = TRUE)
  if (!length(fs)) { ts("    ", ds, ": no objects"); next }
  for (p in fs) {
    smp <- sub("\\.rds$", "", basename(p))
    o <- readRDS(p)
    X <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
    tk <- intersect(T_MK, rownames(X)); bk <- intersect(B_MK, rownames(X))
    he <- if ("high_error" %in% colnames(o@meta.data)) o$high_error else rep(0L, ncol(X))
    if (!length(tk) || !length(bk)) { rm(o, X); next }
    tpos <- Matrix::colSums(X[tk, , drop = FALSE] > 0) > 0
    bpos <- Matrix::colSums(X[bk, , drop = FALSE] > 0) > 0
    keep <- is.na(he) | he != 1            # T3d drops high_error cells; report both
    f <- function(t, b) {
      n <- length(t); pt <- mean(t); pb <- mean(b); both <- sum(t & b)
      list(n = n, pct_T = round(100*pt, 2), pct_B = round(100*pb, 2),
           both_n = both, pct_both = round(100*both/n, 4),
           enrich = if (both >= MIN_BOTH && pt*pb > 0) round((both/n)/(pt*pb), 2) else NA_real_)
    }
    a <- f(tpos, bpos); b2 <- f(tpos[keep], bpos[keep])
    OUT[[length(OUT)+1]] <- data.table(dataset = ds, sample = smp,
      n_all = a$n, pct_T = a$pct_T, pct_B = a$pct_B, both_n = a$both_n, pct_both = a$pct_both, enrich = a$enrich,
      n_he = b2$n, both_n_he = b2$both_n, pct_both_he = b2$pct_both, enrich_he = b2$enrich,
      t_markers = length(tk), b_markers = length(bk))
    rm(o, X); gc(FALSE)
  }
  ts(sprintf("    %-12s %d samples", ds, sum(vapply(OUT, function(x) x$dataset == ds, TRUE))))
}
S <- rbindlist(OUT)
S <- merge(S, cm[, .(dataset, sample, pooled, n_donors_in_library, demux_method)],
           by = c("dataset","sample"), all.x = TRUE)
S[is.na(pooled), pooled := FALSE]
fwrite(S[order(dataset, -n_all)], file.path(TAB, "ws_tb_coexpression_check.csv"))
chk("samples measured", TRUE, nrow(S) > 100L)

ts("[3] per dataset, judged on NON-POOLED samples; pooled reported separately")
SUM <- S[, .(samples = .N, cells = sum(n_all),
             median_pct_T = round(median(pct_T), 2), median_pct_B = round(median(pct_B), 2),
             total_both = sum(both_n), pct_both_pooled_cells = round(100*sum(both_n)/sum(n_all), 4),
             enrich_pooled_over_cells = round((sum(both_n)/sum(n_all)) /
                 (weighted.mean(pct_T/100, n_all) * weighted.mean(pct_B/100, n_all)), 2)),
         by = .(dataset, pooled)]
setorder(SUM, dataset, pooled)
fwrite(SUM, file.path(TAB, "ws_tb_coexpression_summary.csv"))
print(SUM)
ts("    verdict basis -- non-pooled rows only:")
V <- SUM[pooled == FALSE]
V[, verdict := fcase(
  total_both < MIN_BOTH,                        "too few double positives to judge",
  enrich_pooled_over_cells > 2,                 "excess over independence -> doublets likely NOT removed",
  enrich_pooled_over_cells < 1,                 "at or below independence -> no detectable T/B doublet excess",
  default =                                     "between 1 and 2 -> inconclusive")]
print(V[, .(dataset, samples, cells, total_both, pct_both_pooled_cells, enrich_pooled_over_cells, verdict)])
ts("    NOTE: a null here is weak evidence. Expected T/B doublet frequency is about")
ts("    2 * P(T) * P(B) * doublet_rate ~ 0.07% of cells, so this contrast has little power either way.")

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_tb_checks.csv"))
ts(sprintf("[4] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(TAB, "ws_tb_coexpression_check.csv"))
