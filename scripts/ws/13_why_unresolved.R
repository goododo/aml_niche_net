#!/usr/bin/env Rscript
# INPUT  : results/tables/ws/ws_doublet_scores_percell.csv.gz   (per-library scDblFinder, 4 of 9 datasets)
#          /LARGE1/.../02_integrated/ws_panel_cache.rds          (full metadata + the 90-gene panel slice)
#          results/tables/ws/ws_clusters_percell.csv.gz, ws_adjudication_final.csv
# OUTPUT : results/tables/ws/ws_cluster_why_unresolved.csv
#          results/tables/ws/ws_why_unresolved_checks.csv
# what it does: tests, rather than assumes, the proposed explanation for the clusters that four
#   evidence streams could not settle -- that they are either doublet/admixture or poor quality.
#   Three mechanisms, each measured directly:
#     (1) doublet load      -- median scDblFinder score and flagged fraction, from the per-library
#                              run (T3a). Covers GSE185381/GSE227903/GSE239721/GSE289435.
#     (2) library quality   -- nCount / nFeature / percent.mt against the other 56 clusters.
#     (3) mixed identity    -- co-expression ENRICHMENT obs/(p1*p2) of mutually exclusive lineage
#                              pairs. This is the mechanism-specific test: a doublet or an ambient-
#                              contaminated cluster co-expresses markers no single cell should,
#                              while a genuinely novel cell state does not. The ratio is used, not
#                              the raw co-expression rate, because the raw rate rises with depth.
#   Cluster 28 (11,444 cells, still unresolved) and cluster 46 (3,234, suspected single-library
#   artifact) are the targets; all 57 are measured so the targets have a reference distribution.
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table); library(Matrix) })
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- file.path(ROOT, "results/tables/ws")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-50s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

ts("[1] cells, clusters, metadata")
CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell","res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
CA <- readRDS(file.path(BASE, "ws_panel_cache.rds"))
md <- as.data.table(CA$meta)[, .(cell = CA$cells, Dataset_ID, SampleType, Sample_Name,
                                 nCount_RNA, nFeature_RNA, percent.mt)]
setkey(md, cell); md <- md[CL]
md[, bc := sub(".*_", "", cell)]
chk("cells", 895228L, nrow(md))
chk("clusters", 57L, uniqueN(md$cl))

## -- (1) doublet load ----------------------------------------------------------------------------
ts("[2] doublet load")
DB <- fread(file.path(TAB, "ws_doublet_scores_percell.csv.gz"))
DB[, bc := sub(".*_", "", cell)]
DBK <- unique(DB[, .(dataset, ps = sample, bc, sc_score, sc_class)], by = c("dataset","ps","bc"))
# reconcile sample names by prefix, exactly as the CNV and projection joins did, and merge
# explicitly -- a keyed dt[.(a,b)] lookup resolves bare names against the target's own columns
pairs <- rbindlist(lapply(intersect(unique(DBK$dataset), unique(md$Dataset_ID)), function(d) {
  a <- unique(DBK[dataset == d, ps]); b <- unique(md[Dataset_ID == d, Sample_Name])
  rbindlist(lapply(b, function(x) {
    hit <- a[startsWith(a, x) | startsWith(x, a)]
    if (length(hit) == 1L) data.table(Dataset_ID = d, Sample_Name = x, ps = hit) else NULL
  }))
}), fill = TRUE)
ts(sprintf("    %d one-to-one sample pairs over %d datasets", nrow(pairs), uniqueN(pairs$Dataset_ID)))
md <- merge(md, pairs, by = c("Dataset_ID","Sample_Name"), all.x = TRUE)
md <- merge(md, DBK[, .(Dataset_ID = dataset, ps, bc, sc_score, sc_class)],
            by = c("Dataset_ID","ps","bc"), all.x = TRUE)
ts(sprintf("    doublet score available for %d / %d cells (%.1f%%)",
           sum(!is.na(md$sc_score)), nrow(md), 100*mean(!is.na(md$sc_score))))

## -- (3) mixed identity: co-expression enrichment over independence ------------------------------
ts("[3] mixed identity from the panel slice")
X <- CA$panel
stopifnot(identical(colnames(X), CA$cells))
X <- X[, md$cell, drop = FALSE]                 # align by NAME, never by position
PAIRS <- list(c("CD3E","LYZ"), c("CD3E","MS4A1"), c("CD3E","MPO"), c("MS4A1","LYZ"),
              c("HBB","LYZ"), c("CD3E","HBB"))
have <- vapply(PAIRS, function(p) all(p %in% rownames(X)), TRUE)
PAIRS <- PAIRS[have]
chk("exclusive marker pairs usable", length(PAIRS), length(PAIRS))
MIX <- rbindlist(lapply(PAIRS, function(p) {
  a <- X[p[1], ] > 0; b <- X[p[2], ] > 0
  data.table(cl = md$cl, pair = paste(p, collapse = "&"), a = a, b = b)
}))
# The RATIO is unusable at this cluster size and must be reported with its numerator. Clusters 42
# and 45 scored enrichment 53.7 and 52.5 on CD3E&MS4A1 -- each driven by exactly ONE double-positive
# cell out of 4,294 and 3,467, because both markers are near-absent in a stromal cluster so the
# p1*p2 denominator collapses. Across all 57 clusters the largest double-positive count is 18 of
# 23,681 (0.076%), so cross-lineage co-expression is absent everywhere and this test cannot
# discriminate at this scale. Keep `both_n` so nobody reads the ratio on its own again.
MX <- MIX[, {
  obs <- mean(a & b); e <- mean(a) * mean(b)
  .(enrich = if (e > 0) obs / e else NA_real_, both_n = sum(a & b), both_pct = round(100 * obs, 4))
}, by = .(cl, pair)]
MXW <- dcast(MX, cl ~ pair, value.var = "enrich")
MXW[, max_enrich := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = setdiff(names(MXW), "cl")]
NW <- MX[, .(max_both_n = max(both_n), tot_both_n = sum(both_n)), by = cl]
MXW <- merge(MXW, NW, by = "cl")
MXW[max_both_n < 10L, max_enrich := NA_real_]   # a ratio on <10 cells is noise, not a measurement
rm(CA, X, MIX); gc(FALSE)

## -- assemble ------------------------------------------------------------------------------------
ts("[4] per cluster")
Q <- md[, .(cluster_size = .N,
            med_nCount = as.numeric(median(nCount_RNA)), med_nFeature = as.numeric(median(nFeature_RNA)),
            med_pmt = round(median(percent.mt, na.rm = TRUE), 2),
            dbl_cov = round(mean(!is.na(sc_score)), 3),
            med_dbl_score = round(median(sc_score, na.rm = TRUE), 4),
            frac_dbl_flag = round(mean(sc_class[!is.na(sc_class)] %in% c(TRUE,"TRUE","doublet")), 4),
            top_dataset = names(which.max(table(Dataset_ID)))), by = cl]
R <- merge(Q, MXW[, .(cl, max_enrich, max_both_n, tot_both_n,
                      CD3E_LYZ = `CD3E&LYZ`, CD3E_MS4A1 = `CD3E&MS4A1`)], by = "cl")
FIN <- fread(file.path(TAB, "ws_adjudication_final.csv"))
R <- merge(R, FIN[, .(cl, proposed, basis)], by = "cl", all.x = TRUE)
# rank each cluster against the other 56 so "poor" means poor FOR THIS COHORT
for (k in c("med_nCount","med_nFeature","med_pmt","med_dbl_score","frac_dbl_flag","max_enrich"))
  R[[paste0("pct_rank_", k)]] <- round(rank(R[[k]], na.last = "keep") / sum(!is.na(R[[k]])), 3)
fwrite(R[order(-cluster_size)], file.path(TAB, "ws_cluster_why_unresolved.csv"))

ts("[5] the two targets against the cohort")
tgt <- c(28L, 46L)
cat("\n--- 57-cluster reference distribution ---\n")
print(R[, lapply(.SD, function(x) round(quantile(x, c(.1,.5,.9), na.rm = TRUE), 3)),
        .SDcols = c("med_nCount","med_nFeature","med_pmt","med_dbl_score","frac_dbl_flag","max_enrich")])
cat("\n--- targets ---\n")
print(R[cl %in% tgt, .(cl, cluster_size, top_dataset, med_nCount, med_nFeature, med_pmt,
                       dbl_cov, med_dbl_score, frac_dbl_flag, max_enrich)])
cat("\n--- their percentile rank among the 57 (1.0 = worst/highest) ---\n")
print(R[cl %in% tgt, .(cl, pct_rank_med_nCount, pct_rank_med_nFeature, pct_rank_med_pmt,
                       pct_rank_med_dbl_score, pct_rank_frac_dbl_flag, pct_rank_max_enrich)])
cat("\n--- cross-lineage double positives: the absolute counts that the ratio hides ---\n")
print(R[order(-tot_both_n)][1:8, .(cl, cluster_size, tot_both_n, max_both_n, max_enrich,
                                   proposed = substr(proposed, 1, 26))])
cat(sprintf("    largest double-positive count in any cluster: %d of %d cells -> the mixed-identity\n    test is uninformative at this scale and max_enrich is set NA below 10 cells\n",
            max(R$max_both_n), R[which.max(R$max_both_n), cluster_size]))
cat("\n--- the 6 with the highest doublet flag rate ---\n")
print(R[!is.na(frac_dbl_flag)][order(-frac_dbl_flag)][1:6,
      .(cl, cluster_size, dbl_cov, med_dbl_score, frac_dbl_flag, proposed = substr(proposed, 1, 26))])

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_why_unresolved_checks.csv"))
ts(sprintf("[6] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(TAB, "ws_cluster_why_unresolved.csv"))
