#!/usr/bin/env Rscript
# INPUT  : results/tables/ws/ws_anno_draft_res1_with_mk3.csv   (3 independent calls per cluster)
#          results/tables/ws/ws_anno_scores_res1.csv           (57 x 71 score matrix)
#          results/tables/ws/ws_marker_panel_res1.csv          (the independent 14-compartment panel)
#          results/tables/ws/ws_anno_marker_audit.csv          (markers kept per cell type)
#          results/tables/ws/ws_clusters_percell.csv.gz, ws_pca_meta.rds, ws_panel_cache.rds
# OUTPUT : results/tables/ws/ws_anno_consolidated_res1.csv   one row per cluster, with provenance
#          results/tables/ws/ws_anno_dispute_dossier.csv     long-format evidence for the disputed
#          results/tables/ws/ws_anno_dispute_report.txt      the same, readable per cluster
#          results/tables/ws/ws_consolidate_checks.csv
# what it does: turns the three calls into one annotation with explicit provenance, and assembles
#   the evidence for the clusters that are still disputed.
#   Provenance, per the 2026-10-08 decision:
#     agreed      31 clusters -- automatic (>=3 marker types) and the old label agree at the
#                 coarse class; the old label's full string is kept, since it is the user's own
#                 manual work and is the more specific of the two.
#     adopted      7 clusters -- the old label was Ambiguous/Unassigned, so the automatic call is
#                 taken, AS INSTRUCTED. But each is cross-checked against the independent
#                 14-compartment panel from 08_cluster_markers.R, and 3 of the 7 FAIL that check
#                 (clusters 27, 30, 49: the compartment implied by the automatic call scores
#                 z = 0.29 / 0.36 / 0.18, where a real compartment call scores 2.6-7.2). They are
#                 adopted as instructed but flagged `contradicted` so one edit reverses them.
#     disputed    19 clusters -- the two disagree and the old label was not a give-up. Left blank.
#   The same panel cross-check is reported for EVERY cluster, because it already caught a false
#   positive outside the 7: cluster 11 (26,252 cells) is called Stroma: Adipogenic progenitor by
#   the automatic score while its highest panel z over all 14 compartments is 0.08 -- it is not
#   stromal at all, and the old label Blast-like: HSC/MPP (purity 0.933) is almost certainly right.
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table) })
TAB <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-50s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
co <- function(x) sub(":.*", "", x)

## -- 1. the three calls --------------------------------------------------------------------------
ts("[1] the three calls per cluster")
R <- fread(file.path(TAB, "ws_anno_draft_res1_with_mk3.csv"))
R[, `:=`(c_mk3 = co(mk3_top), c_old = co(old_Manual_anno_2))]
R[, `:=`(gaveup = c_old %chin% c("Ambiguous", "Unassigned"), disagree = co(mk3_top) != co(old_Manual_anno_2))]
chk("clusters", 57L, nrow(R))
chk("cells", 895228L, sum(R$cluster_size))
chk("agreed", 31L, sum(!R$disagree))
chk("old label gave up", 7L, sum(R$gaveup))
chk("all give-ups sit inside the disagreements", 7L, sum(R$gaveup & R$disagree))
chk("still disputed", 19L, sum(R$disagree & !R$gaveup))

## -- 2. the independent panel cross-check, for every cluster -------------------------------------
ts("[2] independent 14-compartment panel cross-check")
P <- fread(file.path(TAB, "ws_marker_panel_res1.csv"))
SC <- P[, .(score = mean(mean_expr)), by = .(cl, compartment)]
SC[, z := round((score - mean(score)) / ifelse(sd(score) > 0, sd(score), 1), 2), by = compartment]
# The panel is coarse (14 compartments) while the marker table has 71 types in 23 categories, so
# a 1-to-1 map mis-sends whole categories and manufactures false contradictions. First pass did
# exactly that: 13 of 23 "contradicted" flags were Blast-like (mapped to HSPC, but a blast that has
# acquired a myeloid program scores low on HSPC) or Granulocyte (the panel has NO granulocyte
# compartment at all, so 3/3 were false). Worse, cluster 42 -- Vascular: Mural (Pericyte/SMC), which
# BOTH methods agree on -- was flagged at z = -0.52 because Vascular was sent to Endothelial, while
# the panel's own Pericyte compartment scores z = 7.15 for it. Same failure for Ery/Meg prog (mapped
# to mature platelet genes though it is a progenitor) and DC: pDC (pDCs express none of LYZ/CD14/
# VCAN). So: map each category to the SET of compartments it could legitimately occupy and take the
# best z, and where the panel has no equivalent say so instead of reporting a contradiction.
MAPSET <- list(
  "CD8 T" = "T_NK", "T" = "T_NK", "CD4 T helper" = "T_NK", "T/NK" = "T_NK",
  "Innate-like T" = "T_NK", "NK" = "T_NK",
  "HSPC" = c("HSPC", "Cycling"),
  "Blast-like" = c("HSPC", "Mono_DC", "Erythroid", "Megakaryo", "Cycling"),
  "Ery/Meg prog" = c("Erythroid", "Megakaryo", "HSPC", "Cycling"),
  "Stroma" = c("MSC_stromal", "Adipo", "Pericyte"),
  "Vascular" = c("Endothelial", "Pericyte"),
  "Bone" = c("Osteo", "MSC_stromal"), "Cartilage" = c("Osteo", "MSC_stromal"),
  "Mono" = "Mono_DC", "Myeloid" = "Mono_DC",
  "DC" = c("Mono_DC", "B_plasma"),
  "B" = "B_plasma", "Plasma" = "B_plasma",
  "Erythroid" = "Erythroid", "Megakaryocyte" = "Megakaryo",
  "Granulocyte" = c("Mast", "Mono_DC"))        # Mast only; the panel has no neutrophil set
NO_EQUIV <- c("Granulocyte")                    # flagged separately below
best <- SC[order(cl, -z), .(panel_top = compartment[1], panel_top_z = z[1]), by = cl]
R <- merge(R, best, by = "cl", all.x = TRUE)
best_of_set <- function(k, cls) {
  comps <- MAPSET[[cls]]
  if (is.null(comps)) return(list(comp = NA_character_, z = NA_real_))
  s <- SC[cl == k & compartment %chin% comps][order(-z)]
  if (!nrow(s)) return(list(comp = NA_character_, z = NA_real_))
  list(comp = s$compartment[1], z = s$z[1])
}
R[, c("mapped_compartment", "z_of_mapped") := {
  r <- Map(best_of_set, cl, c_mk3)
  list(vapply(r, function(x) x$comp, character(1)), vapply(r, function(x) x$z, numeric(1)))
}]
# a neutrophil call cannot be checked against a panel with no neutrophil set -- say that, do not
# report a contradiction the evidence cannot support
R[, independent_check := ifelse(c_mk3 %chin% NO_EQUIV & !grepl("Mast", mk3_top), "no panel equivalent",
                         ifelse(is.na(mapped_compartment), "no panel equivalent",
                         ifelse(z_of_mapped >= 2, "supported",
                         ifelse(z_of_mapped >= 1, "weak", "contradicted"))))]
ts("    cross-check over all 57 clusters:")
print(R[, .N, by = independent_check][order(-N)])
chk("cluster 11 still flagged contradicted", "contradicted", R[cl == 11, independent_check])
chk("cluster 42 no longer falsely contradicted (Pericyte z 7.15)", "supported", R[cl == 42, independent_check])

## -- 3. the consolidated annotation ---------------------------------------------------------------
ts("[3] consolidating")
R[, `:=`(provenance = ifelse(!disagree, "agreed", ifelse(gaveup, "adopted", "disputed")))]
R[, label := ifelse(provenance == "agreed", old_Manual_anno_2,
             ifelse(provenance == "adopted", mk3_top, NA_character_))]
CONS <- R[order(cl), .(cl, cluster_size, label, provenance, independent_check,
  z_of_mapped, panel_top, panel_top_z,
  auto_all = all_top, auto_all_n_marker = all_n1, auto_all_diff = all_diff,
  auto_mk3 = mk3_top, auto_mk3_n_marker = mk3_n1, auto_mk3_diff = mk3_diff,
  old_label = old_Manual_anno_2, old_purity = purity_Manual, pct_labelled,
  top_dataset, pct_control)]
fwrite(CONS, file.path(TAB, "ws_anno_consolidated_res1.csv"))
chk("labelled clusters", 38L, sum(!is.na(CONS$label)))
chk("blank (disputed) clusters", 19L, sum(is.na(CONS$label)))
ts("    provenance x independent check:")
print(dcast(CONS[, .N, by = .(provenance, independent_check)], provenance ~ independent_check,
            value.var = "N", fill = 0L))
ts("    adopted clusters that FAIL the independent check (adopted as instructed, flagged):")
print(CONS[provenance == "adopted" & independent_check == "contradicted",
           .(cl, cluster_size, label, z_of_mapped, panel_top, panel_top_z)])

## -- 4. evidence dossier for the disputed clusters ------------------------------------------------
ts("[4] dossier for the ", sum(is.na(CONS$label)), " disputed clusters")
DISP <- CONS[is.na(label), cl]
S <- fread(file.path(TAB, "ws_anno_scores_res1.csv"))
A <- fread(file.path(TAB, "ws_anno_marker_audit.csv"))
nmk <- A[used == TRUE, .(n_marker = uniqueN(marker)), by = cell_type]
S <- merge(S, nmk, by = "cell_type", all.x = TRUE)

CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell", "res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
CA <- readRDS(file.path(BASE, "ws_panel_cache.rds"))
md <- as.data.table(CA$meta)[, .(cell = CA$cells, Dataset_ID, SampleType, Sample_Name,
                                 Manual_anno_2, Sub_anno, nCount_RNA, nFeature_RNA, percent.mt)]
setkey(md, cell); md <- md[CL]; rm(CA); gc(FALSE)

DOS <- list(); RPT <- character()
for (k in DISP) {
  r <- CONS[cl == k]
  sc <- S[cl == k & n_marker >= 3][order(-score)][1:8]
  oldd <- md[cl == k & !is.na(Manual_anno_2) & Manual_anno_2 != "",
             .N, by = Manual_anno_2][order(-N)][1:3]
  oldd[, share := round(N / md[cl == k & !is.na(Manual_anno_2) & Manual_anno_2 != "", .N], 3)]
  subd <- md[cl == k & !is.na(Sub_anno) & Sub_anno != "", .N, by = Sub_anno][order(-N)][1:3]
  dsd <- md[cl == k, .N, by = .(Dataset_ID, SampleType)][order(-N)][1:4]
  q  <- md[cl == k, .(med_nCount = as.numeric(median(nCount_RNA)),
                      med_nFeature = as.numeric(median(nFeature_RNA)),
                      med_pmt = round(median(percent.mt, na.rm = TRUE), 2))]
  pz <- SC[cl == k][order(-z)][1:5]
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "score_top8",
    value = paste(sprintf("%s=%.3f(%dmk)", sc$cell_type, sc$score, sc$n_marker), collapse = " | "))
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "old_Manual_anno_2_top3",
    value = paste(sprintf("%s=%.1f%%", oldd$Manual_anno_2, 100*oldd$share), collapse = " | "))
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "old_Sub_anno_top3",
    value = paste(subd$Sub_anno, collapse = " | "))
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "panel_z_top5",
    value = paste(sprintf("%s=%.2f", pz$compartment, pz$z), collapse = " | "))
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "composition",
    value = paste(sprintf("%s/%s=%d", dsd$Dataset_ID, dsd$SampleType, dsd$N), collapse = " | "))
  DOS[[length(DOS)+1]] <- data.table(cl = k, field = "qc",
    value = sprintf("nCount=%.0f nFeature=%.0f pmt=%.2f", q$med_nCount, q$med_nFeature, q$med_pmt))

  RPT <- c(RPT,
    sprintf("==== cluster %d | %d cells | %s | %.0f%% Control ====", k, r$cluster_size, r$top_dataset, r$pct_control),
    sprintf("  automatic (>=3mk) : %s   (%d markers, margin %.3f)", r$auto_mk3, r$auto_mk3_n_marker, r$auto_mk3_diff),
    sprintf("  automatic (all)   : %s   (%d markers, margin %.3f)", r$auto_all, r$auto_all_n_marker, r$auto_all_diff),
    sprintf("  old label         : %s   (purity %.3f, %.0f%% of cells labelled)", r$old_label, r$old_purity, r$pct_labelled),
    sprintf("  independent panel : %s (z=%.2f for the automatic call's compartment; panel top %s z=%.2f)",
            r$independent_check, ifelse(is.na(r$z_of_mapped), NA_real_, r$z_of_mapped), r$panel_top, r$panel_top_z),
    sprintf("  top-8 scores      : %s", DOS[[length(DOS)-5]]$value),
    sprintf("  old label spread  : %s", DOS[[length(DOS)-4]]$value),
    sprintf("  old Sub_anno      : %s", DOS[[length(DOS)-3]]$value),
    sprintf("  panel z (top 5)   : %s", DOS[[length(DOS)-2]]$value),
    sprintf("  composition       : %s", DOS[[length(DOS)-1]]$value),
    sprintf("  QC                : %s", DOS[[length(DOS)]]$value),
    sprintf("  blast ambiguity   : %s", ifelse(grepl("Blast", r$auto_mk3) | grepl("Blast", r$old_label),
            "YES -- one side is Blast-like; marker scoring cannot separate a blast from the normal cell whose program it has acquired", "no")),
    "")
}
fwrite(rbindlist(DOS), file.path(TAB, "ws_anno_dispute_dossier.csv"))
writeLines(RPT, file.path(TAB, "ws_anno_dispute_report.txt"))
chk("dossier covers every disputed cluster", length(DISP), uniqueN(rbindlist(DOS)$cl))

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_consolidate_checks.csv"))
ts(sprintf("[5] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", TAB)
