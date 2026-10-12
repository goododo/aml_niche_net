#!/usr/bin/env Rscript
# INPUT  : results/tables/ws/ws_final_labels.csv            (cluster -> single fused label)
#          results/tables/ws/ws_cluster_cnv.csv, ws_cluster_projection.csv, ws_cluster_clonality.csv
#          results/tables/ws/ws_clusters_percell.csv.gz, /LARGE1/.../ws_pca_meta.rds
# OUTPUT : results/tables/ws/ws_celltype_malignancy.csv     per cluster x sample-origin, the two axes
#          results/tables/ws/ws_percell_labels.csv.gz       per cell: cell_type + malignancy + basis
#          results/tables/ws/ws_split_checks.csv
# what it does: separates the two axes that the source vocabulary fused into one name.
#   `Blast-like: HSC/MPP` asserts a cell TYPE (progenitor) and a MALIGNANCY STATUS (blast) in one
#   string, and the second assertion failed once the cohort changed: five clusters carrying that
#   prefix are 52.6-87.0% cells from HEALTHY donors (C.12g). Splitting them makes the type
#   survive a cohort change and lets the malignancy claim be made, or withheld, on evidence.
#
#   `cell_type`  -- transcriptional identity only, no malignancy claim.
#   `malignancy` -- malignant / normal / undetermined, assigned at the CLUSTER x SAMPLE-ORIGIN
#                   level, never per cell, because no per-cell instrument here has the specificity:
#                   the inferCNV consensus runs a 12-65% healthy-donor false-positive rate
#                   (results/tables/02_malignancy/malignancy_fpr_*.csv) and the projection residual,
#                   the better of the two at 3.1%, is only meaningful as a rate over a population.
#   `malignancy_basis` -- how each call was reached, so none of them is a bare assertion.
#
#   The decision rule, and why each threshold:
#     Control donor                                   -> normal        (definitional; a healthy
#                                                                       donor's cell is not a blast)
#     AML donor, cluster CNV > 0.654                  -> malignant     (0.654 is the WORST per-bin
#                                                                       healthy false-positive rate,
#                                                                       Stromal; above it the signal
#                                                                       cannot be false positives)
#     AML donor, cluster CNV < 0.118 AND residual
#                < 0.031                              -> normal        (0.118 is the lowest per-bin
#                                                                       baseline; 0.031 is this
#                                                                       cohort's own measured Control
#                                                                       high_error rate)
#     AML donor, anything else                        -> undetermined  (the evidence sits inside the
#                                                                       false-positive band and
#                                                                       cannot separate the two)
#   A large `undetermined` share is the honest outcome, not a failure: it quantifies how strong a
#   malignancy claim this data can support. See the non-identifiability argument in section 2 --
#   a blast that acquires a monocyte program is LABELLED a monocyte, and no marker fixes that.
#   CNV-negative is NOT evidence of benignity on its own (40-50% of AML is cytogenetically normal),
#   which is why the `normal` branch requires the projection residual to agree.
# This script deletes nothing; the old fused `final` column is carried through untouched.

suppressPackageStartupMessages({ library(data.table) })
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- file.path(ROOT, "results/tables/ws")
CNV_MALIG <- 0.654; CNV_CLEAN <- 0.118; RESID_CLEAN <- 0.031
# A cluster that contains healthy-donor cells carries its own matched control, and that is a
# better reference than a cohort-wide average: the projection residual varies by CELL TYPE, while
# d15 stratified the threshold only by PLATFORM, so 0.031 (the cohort's mean Control rate) is not
# a valid per-cluster bar. Clusters 6/8/14 make the error visible -- their OWN Control cells sit at
# 0.076/0.102/0.046, all above 0.031, so the cohort-wide rule marked their AML-derived cells
# undetermined even though the within-cluster test found the two origins indistinguishable
# (d_cnv -0.007/+0.028/-0.006, d_residual +0.052/+0.019/-0.006).
# DISCLOSURE: this reference was changed AFTER seeing the first result. It is corrected because the
# cohort-wide bar was wrong by construction, not because the answer was unwelcome; both versions
# are reported side by side below so the effect of the change is visible.
MIN_CTL_FOR_INTERNAL <- 200L   # a cluster needs this many Control cells to serve as its own control
INTERNAL_RATIO <- 1.5; INTERNAL_ABS <- 0.05
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-54s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. cell_type: strip the malignancy assertion out of the name --------------------------------
ts("[1] cell_type")
L <- fread(file.path(TAB, "ws_final_labels.csv"))
CT <- c(
  "Blast-like: HSC/MPP"                                = "HSC/MPP",
  "Blast-like: GMP"                                    = "GMP",
  "Blast-like: cDC"                                    = "DC: cDC",
  "Blast-like (BAALC+, normal karyotype)"              = "HSPC (BAALC+)",
  "Blast-like primitive progenitor (BAALC+DNTT+PROM1)" = "HSPC: primitive (DNTT+PROM1+)",
  "Blast-like (basophil/mast features, PRAME+)"        = "Basophil/mast progenitor",
  "Blast-like"                                         = "HSPC (unspecified)")
L[, cell_type := ifelse(final %chin% names(CT), unname(CT[final]), final)]
chk("no cell_type still asserts malignancy", 0L, sum(grepl("Blast", L$cell_type)))
chk("clusters", 57L, nrow(L))

LIN <- function(x) fcase(
  grepl("Unresolved|Excluded|Interferon-response", x),            "Unresolved/excluded",
  grepl("^HSC/MPP|^HSPC|^GMP|CLP|LMPP|Myeloid: Generic", x),      "HSPC/progenitor",
  grepl("^T/NK|^T cell|^T:|^CD4 T|^CD8 T|^NK$", x),               "T/NK",
  grepl("^B:|Plasma", x),                                         "B/plasma",
  grepl("Erythroid|Megakaryo|Ery/Meg", x),                        "Ery/Mega",
  grepl("^Stroma|^Vascular|^Bone|^Cartilage", x),                 "Niche (stroma/vessel/bone)",
  grepl("^Mono|^DC|Granulocyte|Mast|[Bb]asophil|macrophage", x),  "Myeloid",
  default = "Unresolved/excluded")
L[, lineage := LIN(cell_type)]
ts("    lineage groups after the split:")
print(L[, .(clusters = .N, cells = sum(cluster_size)), by = lineage][order(-cells)])
chk("HSPC/progenitor is no longer empty", TRUE, "HSPC/progenitor" %in% L$lineage)

## -- 2. malignancy, at cluster x sample-origin ----------------------------------------------------
ts("[2] malignancy")
cv <- fread(file.path(TAB, "ws_cluster_cnv.csv"))[, .(cl, cnv_mal = mal_frac, cnv_called = called_frac)]
pj <- fread(file.path(TAB, "ws_cluster_projection.csv"))[, .(cl, high_err = frac_high_error, proj_cov = proj_frac)]
E  <- merge(merge(L[, .(cl, cluster_size, final, cell_type, lineage)], cv, by = "cl", all.x = TRUE),
            pj, by = "cl", all.x = TRUE)
G <- CJ(cl = E$cl, origin = c("AML","Control"), unique = TRUE)
G <- merge(G, E, by = "cl")
# the matched internal control, where one exists
W <- tryCatch(fread(file.path(TAB, "ws_within_cluster_origin.csv")), error = function(e) NULL)
INT <- if (is.null(W)) data.table(cl = integer(), internal_ok = logical()) else {
  w <- dcast(W, cl ~ SampleType, value.var = c("cells","cnv_mal","high_err"))
  w[, internal_ok := !is.na(cells_Control) & cells_Control >= MIN_CTL_FOR_INTERNAL &
       !is.na(high_err_AML) & !is.na(high_err_Control) &
       (high_err_AML <= high_err_Control * INTERNAL_RATIO | high_err_AML - high_err_Control <= INTERNAL_ABS) &
       !is.na(cnv_mal_AML) & !is.na(cnv_mal_Control) & cnv_mal_AML <= cnv_mal_Control + INTERNAL_ABS]
  w[, .(cl, internal_ok, ctl_high = high_err_Control, aml_high = high_err_AML)] }
G <- merge(G, INT, by = "cl", all.x = TRUE)
G[is.na(internal_ok), internal_ok := FALSE]
G[, malignancy_cohortbar := fcase(            # the first version, kept for comparison
  origin == "Control",                                                        "normal",
  !is.na(cnv_mal) & cnv_mal > CNV_MALIG,                                      "malignant",
  !is.na(cnv_mal) & !is.na(high_err) & cnv_mal < CNV_CLEAN & high_err < RESID_CLEAN, "normal",
  default =                                                                   "undetermined")]
G[, malignancy := fcase(
  origin == "Control",                                                        "normal",
  !is.na(cnv_mal) & cnv_mal > CNV_MALIG,                                      "malignant",
  internal_ok == TRUE,                                                        "normal",
  !is.na(cnv_mal) & !is.na(high_err) & cnv_mal < CNV_CLEAN & high_err < RESID_CLEAN, "normal",
  default =                                                                   "undetermined")]
G[, malignancy_basis := fcase(
  origin == "Control",                                 "healthy donor (definitional)",
  malignancy == "malignant",                           sprintf("cluster CNV %.3f > worst healthy FPR baseline %.3f", cnv_mal, CNV_MALIG),
  malignancy == "normal" & internal_ok,                sprintf("AML donor; indistinguishable from the SAME cluster's %.0f healthy-donor cells (residual %.3f vs %.3f)", NA_real_, aml_high, ctl_high),
  malignancy == "normal",                              sprintf("AML donor; cluster CNV %.3f < %.3f and projection residual %.3f < this cohort's Control rate %.3f", cnv_mal, CNV_CLEAN, high_err, RESID_CLEAN),
  is.na(cnv_mal) | is.na(high_err),                    "AML donor; no CNV or projection coverage for this cluster",
  default =                                            sprintf("AML donor; CNV %.3f and residual %.3f both inside the false-positive band", cnv_mal, high_err))]
fwrite(G[order(cl, origin)], file.path(TAB, "ws_celltype_malignancy.csv"))

## -- 3. per cell ----------------------------------------------------------------------------------
ts("[3] per cell")
CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell","res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
M  <- readRDS(file.path(BASE, "ws_pca_meta.rds"))
md <- as.data.table(M$meta)[, .(cell, Dataset_ID, SampleType, Sample_Name)]
setkey(md, cell); md <- md[CL]
md[, origin := SampleType]
P <- merge(md, G[, .(cl, origin, cell_type, lineage, malignancy, malignancy_cohortbar,
                     malignancy_basis, final)], by = c("cl","origin"), all.x = TRUE)
chk("cells labelled", 895228L, nrow(P))
chk("no cell without a cell_type", 0L, sum(is.na(P$cell_type)))
chk("no cell without a malignancy value", 0L, sum(is.na(P$malignancy)))
fwrite(P[order(cl)], file.path(TAB, "ws_percell_labels.csv.gz"), compress = "gzip")

ts("[4] the two axes")
cat("\n--- malignancy over all 895,228 cells ---\n")
t1 <- P[, .(cells = .N), by = malignancy][order(-cells)][, pct := round(100*cells/sum(cells),1)][]
print(t1)
cat("\n--- effect of the reference change (AML-donor cells) ---\n")
print(merge(P[origin == "AML", .(cohort_wide_bar = .N), by = .(malignancy = malignancy_cohortbar)],
            P[origin == "AML", .(matched_internal_control = .N), by = malignancy],
            by = "malignancy", all = TRUE))
cat("\n--- malignancy within AML-donor cells only ---\n")
t2 <- P[origin == "AML", .(cells = .N), by = malignancy][order(-cells)][, pct := round(100*cells/sum(cells),1)][]
print(t2)
cat("\n--- cell_type lineage x malignancy ---\n")
print(dcast(P[, .N, by = .(lineage, malignancy)], lineage ~ malignancy, value.var = "N", fill = 0L))
cat("\n--- the five clusters whose Blast-like label failed (C.12g) ---\n")
print(P[cl %in% c(6,8,14,36,56), .(cells = .N), by = .(cl, cell_type, origin, malignancy)][order(cl, origin)])

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_split_checks.csv"))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", TAB)
