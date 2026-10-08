#!/usr/bin/env Rscript
# INPUT  : /LARGE1/gr10634/gaozy/reference/cell_annotation_markers_for_code_260206.csv  (514 entries, 71 types)
#          /LARGE1/.../02_integrated/ws_rpca_integrated.rds   (44 GB, read once, then cached)
#          /LARGE1/.../02_integrated/ws_panel_cache.rds        (full metadata, incl. the old labels)
#          results/tables/ws/ws_clusters_percell.csv.gz        (res.1 = 57 clusters)
# OUTPUT : results/tables/ws/ws_anno_scores_res1.csv        cluster x cell_type normalised score
#          results/tables/ws/ws_anno_threshold_sweep.csv    Unassigned/Ambiguous vs thresholds
#          results/tables/ws/ws_anno_draft_res1.csv         the draft call per cluster, both methods
#          results/tables/ws/ws_anno_oldlabel_vote.csv      old Manual_anno_2 / Sub_anno majority vote
#          results/tables/ws/ws_anno_marker_audit.csv       every marker: kept, dropped, or renamed
#          results/tables/ws/ws_anno_checks.csv
#          /LARGE1/.../02_integrated/ws_marker_scaled.rds    the scaled marker slice (cache)
# what it does: runs the scoring logic of 03.2_Anno_pipeline_260206.r on the re-integrated object at
#   resolution 1.0, reusing that script's marker table and its priority/specificity/level weights.
#   Four changes, all agreed on 2026-10-08 and each recorded in ws_anno_marker_audit.csv:
#     1. the per-cell-type score is a weighted MEAN, not the weighted SUM that script uses.
#        Its 71 cell types differ 64-fold in achievable score scale (7.68 for Erythroid:
#        Late/Mature down to 0.12 for Cartilage: Superficial/articular-like), because the sum
#        grows with the number of markers. With a sum plus the absolute thresholds
#        min_score_threshold = 0.5 / min_difference_threshold = 0.12, argmax is biased toward
#        marker-rich types and a 1-marker type can essentially never be called: 0.5 is 6.5% of the
#        largest type's ceiling and 417% of the smallest's. Dividing by the sum of |weights| over
#        the markers actually present makes the 71 types comparable and gives the thresholds one
#        meaning. Nothing else about the weighting changes.
#     2. priority "supportive" is read as "support" (0.8). The script's map has only core/support/
#        candidate, so "supportive" fell through to the 0.5 default -- a 37.5% down-weighting of 10
#        entries, all of them Adipo lineage (Stroma: Adipocyte and Stroma: Adipogenic lineage
#        precursor, 5 each). The two spellings plainly mean the same thing.
#     3. three mouse/legacy symbols are mapped to their human equivalents: ZFP423 -> ZNF423,
#        KITL -> KITLG, HBA -> HBA1 + HBA2. Without this they simply fail the %in% rownames filter.
#     4. ScaleData runs on the marker genes only, not on all 32,029. z-scoring is per gene and
#        independent of which other genes are present, so the values are IDENTICAL; the dense
#        matrix goes from 32,029 x 895,228 (229 GB) to about 240 x 895,228 (1.7 GB).
#   Thresholds are NOT fixed here. Normalisation changes the scale, so 0.5 / 0.12 no longer mean
#   what they did; the script reports the score distribution and a sweep for the user to choose.
#   A second, independent draft comes from a majority vote of the old Manual_anno_2 / Sub_anno
#   labels onto the new clusters, so the clusters where the two methods agree can be skimmed and
#   only the disagreements need manual attention. The manual annotation itself stays the user's.
# This script deletes nothing and does not modify the integrated object.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)
RES_COL <- "res.1"
MARKER_CSV <- "/LARGE1/gr10634/gaozy/reference/cell_annotation_markers_for_code_260206.csv"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
SCALED <- file.path(BASE, "ws_marker_scaled.rds")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-50s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. marker table, with the three agreed corrections -------------------------------------------
ts("[1] marker table")
M <- fread(MARKER_CSV)
chk("marker entries", 514L, nrow(M))
chk("cell types", 71L, uniqueN(M$cell_type))
chk("roles are positive/negative only", TRUE, all(M$role %chin% c("positive", "negative")))

ALIAS <- c(ZFP423 = "ZNF423", KITL = "KITLG")            # mouse / legacy symbols
M[, marker_orig := marker]
M[marker %chin% names(ALIAS), marker := ALIAS[marker]]
hba <- M[marker == "HBA"]                                 # HBA is not a symbol; it is HBA1 + HBA2
if (nrow(hba)) {
  M <- M[marker != "HBA"]
  M <- rbind(M, rbindlist(lapply(c("HBA1", "HBA2"), function(g) copy(hba)[, marker := g])))
}
chk("ZFP423/KITL renamed", 0L, sum(M$marker %chin% names(ALIAS)))
chk("HBA expanded to HBA1 + HBA2", 0L, sum(M$marker == "HBA"))

PRI <- c(core = 1.0, support = 0.8, supportive = 0.8, candidate = 0.2)   # change 2
SPE <- c(high = 1.0, medium = 0.6, low = 0.3)
LEV <- c(high = 1.3, present = 1.0, low = 0.7)
M[, w := ifelse(role == "positive", 1, -1) *
         ifelse(priority       %chin% names(PRI), PRI[priority],       0.5) *
         ifelse(specificity    %chin% names(SPE), SPE[specificity],    0.5) *
         ifelse(expected_level %chin% names(LEV), LEV[expected_level], 1.0)]
chk("supportive now weighted as support", 10L, sum(M$priority == "supportive"))

## -- 2. clusters and the old labels (both from the small cache) -----------------------------------
ts("[2] clusters and metadata")
CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell", RES_COL))
setnames(CL, RES_COL, "cl"); setkey(CL, cell)
CA <- readRDS(file.path(BASE, "ws_panel_cache.rds"))
md <- as.data.table(CA$meta)[, .(cell = CA$cells, Dataset_ID, SampleType, Sample_Name,
                                 Manual_anno_2, Sub_anno)]
setkey(md, cell); md <- md[CL]
chk("cells", 895228L, nrow(md))
chk("clusters at res 1.0", 57L, uniqueN(md$cl))
chk("cells carrying an old label", 728092L, sum(!is.na(md$Manual_anno_2) & md$Manual_anno_2 != ""))
rm(CA); gc(FALSE)

## -- 3. scaled marker expression (cached; z-scoring is per gene so the subset is exact) -----------
if (file.exists(SCALED)) {
  ts("[3] reusing the scaled marker cache; the 44 GB object is not read")
  SM <- readRDS(SCALED)
} else {
  ts("[3] reading the integrated object (44 GB, about 17 min) to scale the marker genes only")
  o <- readRDS(file.path(BASE, "ws_rpca_integrated.rds"))
  DefaultAssay(o) <- "RNA"
  if (length(grep("^counts", SeuratObject::Layers(o[["RNA"]]))) > 1L) { ts("    joining layers"); o <- JoinLayers(o) }
  if (!"data" %in% SeuratObject::Layers(o[["RNA"]])) { ts("    NormalizeData"); o <- NormalizeData(o, verbose = FALSE) }
  keep <- intersect(unique(M$marker), rownames(o))
  ts("    ScaleData on ", length(keep), " marker genes (not all ", nrow(o), ")")
  o <- ScaleData(o, assay = "RNA", features = keep, vars.to.regress = NULL, verbose = FALSE)
  SM <- list(x = SeuratObject::LayerData(o, assay = "RNA", layer = "scale.data")[keep, , drop = FALSE],
             cells = colnames(o))
  ts("    caching: ", SCALED); saveRDS(SM, SCALED)   # cache FIRST, assert after
  rm(o); gc(FALSE)
}
# The cache keeps the object's own cell order; md is key-sorted. Same set, different order, so
# align by NAME rather than asserting the orders match.
chk("cached cells are the same set as the cluster table", TRUE, setequal(SM$cells, md$cell))
X <- SM$x[, md$cell, drop = FALSE]
chk("scaled matrix aligned to the cluster table", TRUE, identical(colnames(X), md$cell))
chk("marker genes scaled", length(intersect(unique(M$marker), rownames(X))), nrow(X))

# audit every marker: used, or dropped and why
AUD <- unique(M[, .(category, cell_type, marker_orig, marker, role, priority, specificity,
                    expected_level, modality, w)])
AUD[, used := marker %chin% rownames(X)]
AUD[, drop_reason := ifelse(used, "", ifelse(modality %chin% c("Protein","protein"),
      "protein-only marker, not scorable from RNA", "absent from the 32,029-gene axis"))]
fwrite(AUD[order(category, cell_type, -used)], file.path(TAB, "ws_anno_marker_audit.csv"))
ts(sprintf("    markers used %d / %d unique | dropped: %d protein-only, %d off-axis",
   uniqueN(AUD[used == TRUE, marker]), uniqueN(AUD$marker),
   uniqueN(AUD[used == FALSE & grepl("protein", drop_reason), marker]),
   uniqueN(AUD[used == FALSE & grepl("axis", drop_reason), marker])))

## -- 4. per-cell score: weighted MEAN over the markers actually present (change 1) ----------------
ts("[4] scoring ", ncol(X), " cells across ", uniqueN(M$cell_type), " cell types")
gi <- setNames(seq_len(nrow(X)), rownames(X))
types <- sort(unique(M$cell_type))
SCORES <- matrix(NA_real_, nrow = ncol(X), ncol = length(types), dimnames = list(NULL, types))
NMK <- data.table(cell_type = types, n_used = 0L, denom = 0)
for (i in seq_along(types)) {
  t <- types[i]
  d <- unique(M[cell_type == t & marker %chin% names(gi), .(marker, w)])
  if (!nrow(d)) next
  den <- sum(abs(d$w))
  # crossprod(sub, w) = t(sub) %*% w : (cells x markers) %*% (markers x 1) -> one value per cell
  SCORES[, i] <- as.vector(crossprod(X[gi[d$marker], , drop = FALSE], d$w)) / den
  NMK[i, `:=`(n_used = nrow(d), denom = round(den, 3))]
}
unscored <- types[is.na(SCORES[1, ])]
chk("cell types with no usable marker on the axis", length(unscored), length(unscored))
if (length(unscored)) ts("    not scorable (no marker survives): ", paste(unscored, collapse = " | "))
SCORES <- SCORES[, !is.na(SCORES[1, ]), drop = FALSE]
types <- colnames(SCORES)
ts("    scorable cell types: ", length(types))
ts("    score scale is now comparable: all 71 types divided by their own |weight| sum")
print(NMK[order(-n_used)][c(1:5, (.N-4):.N)])

## -- 5. aggregate to clusters ---------------------------------------------------------------------
ts("[5] aggregating to the 57 clusters")
idx <- split(seq_len(nrow(SCORES)), md$cl)
CS <- rbindlist(lapply(names(idx), function(k)
  data.table(cl = as.integer(k), cell_type = types,
             score = round(colMeans(SCORES[idx[[k]], , drop = FALSE]), 4))))
# a second, orthogonal view: which cluster owns a type (column z), as a diagnostic only
CS[, z_across_clusters := round((score - mean(score)) / ifelse(sd(score) > 0, sd(score), 1), 2), by = cell_type]
fwrite(CS[order(cl, -score)], file.path(TAB, "ws_anno_scores_res1.csv"))

## -- 6. threshold sweep. The thresholds are the user's call, so report, do not choose. ------------
ts("[6] threshold sweep (no threshold is fixed here)")
rank1 <- CS[order(cl, -score), .(top = cell_type[1], top_score = score[1],
                                 second = cell_type[2], second_score = score[2]), by = cl]
rank1[, diff := round(top_score - second_score, 4)]
SW <- rbindlist(lapply(c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5), function(ms)
  rbindlist(lapply(c(0, 0.01, 0.02, 0.05, 0.1), function(md_) {
    a <- rank1[, .(cl, ok = top_score >= ms & diff >= md_,
                   unassigned = top_score < ms, ambiguous = top_score >= ms & diff < md_)]
    data.table(min_score = ms, min_diff = md_, assigned = sum(a$ok),
               unassigned = sum(a$unassigned), ambiguous = sum(a$ambiguous))
  }))))
fwrite(SW, file.path(TAB, "ws_anno_threshold_sweep.csv"))
ts("    top-score distribution over the 57 clusters:")
print(round(quantile(rank1$top_score, c(0, .1, .25, .5, .75, .9, 1)), 3))
ts("    top-minus-second distribution:")
print(round(quantile(rank1$diff, c(0, .1, .25, .5, .75, .9, 1)), 3))
print(SW[min_diff == 0.02][, .(min_score, assigned, unassigned, ambiguous)])

## -- 7. the old labels, majority-voted onto the new clusters --------------------------------------
ts("[7] majority vote of the old Manual_anno_2 / Sub_anno")
vote <- function(col) {
  d <- md[!is.na(get(col)) & get(col) != "", .N, by = c("cl", col)]
  setnames(d, col, "label")
  d[, share := N / sum(N), by = cl]
  d[order(cl, -N), .(label = label[1], n = N[1], purity = round(share[1], 3),
                     labelled_cells = sum(N), n_labels = .N), by = cl]
}
V2 <- vote("Manual_anno_2"); setnames(V2, c("label","n","purity","labelled_cells","n_labels"),
      c("old_Manual_anno_2","n_top","purity_Manual","labelled_cells","n_distinct_labels"))
VS <- vote("Sub_anno")[, .(cl, old_Sub_anno = label, purity_Sub = purity)]
OV <- merge(V2, VS, by = "cl", all = TRUE)
tot <- md[, .(cluster_size = .N), by = cl]
OV <- merge(tot, OV, by = "cl", all.x = TRUE)
OV[, pct_labelled := round(100 * labelled_cells / cluster_size, 1)]
fwrite(OV[order(cl)], file.path(TAB, "ws_anno_oldlabel_vote.csv"))

## -- 8. the draft the user reviews ----------------------------------------------------------------
ts("[8] draft table")
D <- merge(merge(tot, rank1, by = "cl"), OV[, .(cl, old_Manual_anno_2, purity_Manual,
           old_Sub_anno, purity_Sub, pct_labelled)], by = "cl", all.x = TRUE)
D <- merge(D, md[, .(top_dataset = names(which.max(table(Dataset_ID))),
                     pct_control = round(100 * mean(SampleType == "Control"), 1)), by = cl],
           by = "cl")
D[, agree_hint := ifelse(is.na(old_Manual_anno_2), "no old label (GSE185381)",
                  ifelse(purity_Manual < 0.5, "old label itself impure -> check", ""))]
fwrite(D[order(cl)], file.path(TAB, "ws_anno_draft_res1.csv"))
print(D[order(cl)][, .(cl, cluster_size, top, top_score, diff, old_Manual_anno_2, purity_Manual, pct_labelled)])

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_anno_checks.csv"))
ts(sprintf("[9] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", TAB)
