#!/usr/bin/env Rscript
# 15_pi_panels.R ----
# The six figure/table panels that answer the PI's questions in the PI's own order. Nothing new is
# computed here: this re-presents 09/10/11/12/14 against the original brief so each question has one
# artefact to point at. Where the honest answer is "the data cannot address this", the panel SHOWS
# that rather than omitting the question.
#
#   P1  which samples were used, and every CNV / genotype number behind the grouping
#   P2  myelofibrosis-related gene expression
#   P3  cytokine / growth-factor / ECM-component gene expression
#   P4  pathway activity: TGF-beta, ECM organisation & receptor interaction, EMT, hypoxia, inflammation
#   P5  CCC: leukaemia cells, HSC, macrophage, niche, MSC, adipocyte
#   P6  Nectin-4 (PVRL4)
#
# THE EFFECT MEASURE throughout is "in how many of the 9 matched pairs is TP53-mut higher", because
# that is the quantity the paired signed-rank test actually operates on. The cell-pooled odds ratio
# is NOT used for any claim: it is cell-weighted, so one large sample can set it, and it disagrees
# in sign with the paired measure for 51 of 145 panel genes.
#
# INPUT  : LCC_TAB_DIR/{09,10,11,12,14}_*.csv, 02_sample_cnv_proxy_all_datasets.csv,
#          04_detection_by_sample.csv, 04_pathway_sample.csv
# OUTPUT : LCC_FIG_DIR/P1..P6 (png + pdf) and LCC_TAB_DIR/P1..P6 (csv)
# Usage  : Rscript LCC_proj/scripts/15_pi_panels.R

suppressPackageStartupMessages({
  library(data.table); library(here); library(ggplot2); library(ggrepel); library(scales)
  library(patchwork)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
set.seed(SEED)


message("[load]")
des  <- fread(file.path(LCC_TAB_DIR, "11_matched_design.csv"))
grp  <- fread(file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))
cnv  <- fread(file.path(LCC_TAB_DIR, "02_sample_cnv_proxy_all_datasets.csv"))
gr   <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))
rb   <- fread(file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))[design == "A_sample_strict"]
det  <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))
pws  <- fread(file.path(LCC_TAB_DIR, "04_pathway_sample.csv"))
panel<- fread_commented(LCC_GENE_PANEL_TSV)
arms <- rbind(des[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              des[, .(pair_id, dataset, sample = wt_sample,  arm = "TP53-WT")])

## ===== P1. cohort and the evidence behind the grouping =====================================
message("[P1] cohort + CNV evidence")
p1 <- merge(arms, grp[, .(dataset, sample, timepoint, uid_patient, tp53_tier, geno_variants,
                          n_cells_tp53_mut, numbat_17p_loh, arm17p_frac_loh, numbat_llr)],
            by = c("dataset", "sample"), all.x = TRUE)
p1 <- merge(p1, cnv[, .(dataset, sample, icnv_max_17p_frac = max_17p_frac,
                        icnv_17p_pos = p17_pos, icnv_ck_arms = ck_max_arms)],
            by = c("dataset", "sample"), all.x = TRUE)
p1 <- merge(p1, unique(det[stratum == "all", .(dataset, sample, n_cells)]),
            by = c("dataset", "sample"), all.x = TRUE)
setorder(p1, pair_id, -arm)
fwrite_safe(p1, file.path(LCC_TAB_DIR, "P1_cohort_evidence.csv"))

# P1 SHOWS ONLY THE EVIDENCE THAT ACTUALLY DEFINED THE GROUPS -- genotype and Numbat 17p LOH.
# An earlier version also carried the two inferCNV columns. That was a presentation error: inferCNV
# was TESTED and REJECTED, so putting it beside the evidence that was used invites the reader to
# think it contributed. The inferCNV result is real and worth reporting -- the PI's first question
# was literally "use inferCNV to split the groups" -- so it now has its own figure, P1B, where it
# reads as what it is: a method check with a negative outcome.
# The tier goes in the row label as WORDS, not a letter. 809653's mutation comes from the study's
# targeted sequencing rather than from the single-cell reads, so its genotype cell count is legitimately
# blank -- with a bare "[tier A]" that blank reads as missing evidence, which is the opposite of true.
TIERLAB <- c(A_genotype   = "single-cell genotype",
             B_reported   = "reported mutation",
             B_allele_loh = "Numbat 17p LOH")
p1[, tier := fifelse(arm == "TP53-WT", "", paste0("   ", TIERLAB[tp53_tier]))]
p1[is.na(tier), tier := ""]
p1l <- melt(p1[, .(lab = paste0(sample, "  (pair ", pair_id, ")", tier), arm,
                   `TP53 mutation\nread directly\n(cells)` = as.numeric(n_cells_tp53_mut),
                   `Numbat\n17p LOH\n(fraction)` = as.numeric(arm17p_frac_loh),
                   `cells in\nsample` = as.numeric(n_cells))],
            id.vars = c("lab", "arm"), variable.name = "evidence", value.name = "v")
p1l[is.na(v), v := 0]
p1l[, vs := v / max(v, na.rm = TRUE), by = evidence]           # each column on its own scale
# per-column formatting: counts as plain separated integers (scientific notation is unreadable for a
# clinical reader), the LOH fraction to three decimals
p1l[, txt := fifelse(v == 0, "-",
              fifelse(grepl("LOH", evidence), formatC(v, format = "f", digits = 3),
                      formatC(v, format = "d", big.mark = ",")))]
ord <- p1[order(pair_id, arm != "TP53-mut")][, paste0(sample, "  (pair ", pair_id, ")", tier)]
p1l[, lab := factor(lab, levels = rev(ord))]
fP1 <- ggplot(p1l, aes(evidence, lab, fill = vs)) +
  geom_tile(colour = SURF, linewidth = 1.1) +
  geom_text(aes(label = txt), size = 4.2, colour = INK[["primary"]]) +
  facet_wrap(~ arm, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "#eef3fa", high = PAL[["blue"]], guide = "none") +
  scale_x_discrete(position = "top") +
  labs(title = "P1  Samples used, and the evidence that defined the two groups",
       subtitle = "tier A = TP53 mutation read directly; tier B = Numbat 17p LOH; WT = neither",
       x = NULL, y = NULL) +
  theme_lcc() + theme(panel.grid.major = element_blank())
save_fig(fP1, "P1_cohort_evidence", 13, 8)

## ===== P1B. the inferCNV check, reported separately because it FAILED ======================
# The PI's question 1 was "use inferCNV to split TP53-mut from TP53-WT". This is the answer, and the
# answer is that it cannot be done on this data. Reported rather than quietly dropped.
message("[P1B] inferCNV vs genotype truth")
pv <- fread(file.path(LCC_TAB_DIR, "07_proxy_vs_genotype.csv"))
# 07 writes the three evaluable samples as a single semicolon-joined `detail` string; parse it back
# rather than retyping the calls, so the figure cannot drift from the table.
ev <- rbindlist(lapply(strsplit(pv$detail, "; *")[[1]], function(s) {
  smp <- sub(":.*$", "", s)
  truth <- ifelse(grepl("mutant", s), "TP53 mutant", "TP53 wild-type")
  call  <- ifelse(grepl("17p=TRUE", s), "17p POSITIVE", "17p negative")
  data.table(sample = smp, truth = truth, call = call)
}))
ev[, verdict := fifelse(truth == "TP53 mutant" & call == "17p negative", "missed the mutation",
                 fifelse(truth == "TP53 wild-type" & call == "17p POSITIVE", "false alarm", "correct"))]
gvar <- fread(file.path(LCC_TAB_DIR, "07_tp53_genotype_sample.csv"))[, .(sample, tp53_variants)]
ev[gvar, variant := i.tp53_variants, on = "sample"]
ev[, row := paste0(sample, ifelse(nzchar(variant) & !is.na(variant), paste0("\n", variant), ""))]
fwrite_safe(ev, file.path(LCC_TAB_DIR, "P1B_infercnv_check.csv"))

evl <- melt(ev[, .(row, `genotype truth` = truth, `inferCNV call` = call, verdict)],
            id.vars = "row", variable.name = "col", value.name = "txt")
evl[, wrong := col == "verdict" & txt != "correct"]
fP1B <- ggplot(evl, aes(col, factor(row, levels = rev(sort(unique(row)))))) +
  geom_tile(aes(fill = wrong), colour = SURF, linewidth = 1.1) +
  geom_text(aes(label = txt), size = 4.4, colour = INK[["primary"]]) +
  scale_fill_manual(values = c(`TRUE` = "#fbe0da", `FALSE` = "#eef3fa"), guide = "none") +
  scale_x_discrete(position = "top") +
  labs(title = "P1B  inferCNV was tested against genotype truth, and it failed",
       subtitle = sprintf("%d samples have both a genotype and an inferCNV call; %d correct",
                          nrow(ev), sum(ev$verdict == "correct")),
       x = NULL, y = NULL) +
  theme_lcc() + theme(panel.grid.major = element_blank())
save_fig(fP1B, "P1B_infercnv_check", 13, 4.6)
cat("\n-- P1B: inferCNV vs genotype --\n"); print(ev[, .(sample, truth, call, verdict)])
cat(sprintf("   sensitivity %d/%d   specificity %d/%d\n",
            pv$true_pos, pv$true_pos + pv$false_neg, pv$true_neg, pv$true_neg + pv$false_pos))

## ===== shared machinery for P2 / P3 =========================================================
gene_panel_fig <- function(cats, title, subtitle, caption, file, h) {
  d <- merge(gr[stratum == "all" & category %in% cats,
                .(gene, category, n_pairs, frac_pairs_mut_higher, median_delta_log2cpm,
                  p_sample_higher, pct_mut, pct_WT)],
             rb[, .(gene = feature, frac_p05)], by = "gene", all.x = TRUE)
  d <- d[!is.na(frac_pairs_mut_higher) & n_pairs >= 6]
  d[, n_up := frac_pairs_mut_higher * n_pairs]
  d[, sig := fifelse(p_sample_higher <= 0.05, "p < 0.05 (uncorrected)", "not significant")]
  setorder(d, category, n_up)
  d[, gene := factor(gene, levels = unique(gene))]
  fwrite_safe(d, file.path(LCC_TAB_DIR, paste0(file, ".csv")))
  p <- ggplot(d, aes(n_up, gene)) +
    # 50% of pairs is the no-difference line; with 9 pairs the test needs 8-9 to clear p = 0.05
    geom_vline(xintercept = 4.5, linewidth = 0.4, linetype = "22", colour = INK[["muted"]]) +
    geom_segment(aes(x = 4.5, xend = n_up, yend = gene), linewidth = 1.1, colour = GRID) +
    geom_point(aes(colour = sig, size = frac_p05)) +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    scale_colour_manual(values = c(`p < 0.05 (uncorrected)` = PAL[["orange"]],
                                   `not significant` = INK[["muted"]]), name = NULL) +
    scale_size_continuous(range = c(1.2, 4.5), limits = c(0, 1),
                          name = "robust to control-set choice", labels = percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = 0:9, limits = c(0, 9)) +
    labs(title = title, subtitle = subtitle,
         x = "number of the 9 matched pairs in which TP53-mut is HIGHER", y = NULL,
         caption = caption) +
    theme_lcc()
  save_fig(p, file, 12, h)
  d
}
CAP_COMMON <- paste("Dashed line = 4.5 of 9, i.e. no difference. Point size is the share of 500",
                    "alternative valid control sets in\nwhich the gene reached p < 0.05 -- small",
                    "points depend on which controls were drawn. NOTHING here survives\nFDR",
                    "correction once patients rather than samples are the unit (6 patients; one",
                    "patient, AML328, supplies 4 of\nthe 9 mutant samples). Treat as preliminary.")

## ===== P2. myelofibrosis-related genes ======================================================
message("[P2] myelofibrosis genes")
gene_panel_fig(c("collagen_ecm_structural", "ecm_remodeling", "megakaryocyte_fibrosis_axis"),
  "P2  Myelofibrosis-related gene expression, TP53-mut vs TP53-WT",
  paste("Collagen / reticulin structural genes, the ECM-remodelling enzymes that build and cross-link",
        "them, and the\nmegakaryocyte axis that drives marrow fibrosis."),
  CAP_COMMON, "P2_myelofibrosis_genes", 9)

## ===== P3. cytokines, growth factors, ECM components ========================================
message("[P3] cytokine / GF / ECM genes")
gene_panel_fig(c("profibrotic_cytokine_gf", "inflammatory_cytokine", "msc_niche", "hypoxia", "emt"),
  "P3  Cytokine, growth-factor and ECM-component gene expression, TP53-mut vs TP53-WT",
  paste("The profibrotic cytokine / growth-factor axis (IL11, OSM, TGFB, PDGF, SPP1), inflammatory",
        "cytokines, niche\nfactors, and the hypoxia and EMT gene sets."),
  CAP_COMMON, "P3_cytokine_ecm_genes", 11)

## ===== P4. pathway activity =================================================================
message("[P4] pathway activity")
FAMILY <- c(HALLMARK_TGF_BETA_SIGNALING_UCell = "TGF-beta signalling",
            REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS_UCell = "TGF-beta signalling",
            KEGG_TGF_BETA_SIGNALING_PATHWAY_UCell = "TGF-beta signalling",
            REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_UCell = "ECM organisation / receptor",
            KEGG_ECM_RECEPTOR_INTERACTION_UCell = "ECM organisation / receptor",
            REACTOME_COLLAGEN_FORMATION_UCell = "ECM organisation / receptor",
            HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_UCell = "EMT",
            HALLMARK_HYPOXIA_UCell = "Hypoxia",
            HALLMARK_INFLAMMATORY_RESPONSE_UCell = "Inflammatory response",
            HALLMARK_TNFA_SIGNALING_VIA_NFKB_UCell = "Inflammatory response",
            HALLMARK_IL6_JAK_STAT3_SIGNALING_UCell = "Inflammatory response")
sc <- names(FAMILY)
pw <- melt(pws[unit == "sample" & group == "all"], id.vars = c("dataset", "sample"),
           measure.vars = sc, variable.name = "pathway", value.name = "score")
pw[, pathway := as.character(pathway)]
pw <- merge(pw, arms, by = c("dataset", "sample"))
pw[, family := FAMILY[pathway]]
# Same treatment as F3: the collection prefix is dropped and underscores become spaces, then the
# name is wrapped, because full MSigDB names clip the facet strips at any readable font size.
pw[, short := gsub("_", " ", sub("^(HALLMARK|REACTOME|KEGG)_", "", sub("_UCell$", "", pathway)))]
pw[, short := vapply(short, function(x) paste(strwrap(x, width = 22), collapse = "\n"), "")]
pw[, arm_s := fifelse(arm == "TP53-mut", "mut", "WT")]
fwrite_safe(pw, file.path(LCC_TAB_DIR, "P4_pathway_scores.csv"))
# One neutral colour, not one hue per pair: there are 9 pairs and the categorical palette holds 8,
# and pair identity carries no meaning here beyond joining the two points of a pair. A single series
# needs no legend -- the subtitle names it.
fP4 <- ggplot(pw[!is.na(score)], aes(arm_s, score, group = pair_id)) +
  geom_line(linewidth = 0.7, alpha = 0.5, colour = PAL[["blue"]]) +
  geom_point(size = 1.9, stroke = 0, alpha = 0.75, colour = PAL[["blue"]]) +
  facet_wrap(~ short, scales = "free_y", ncol = 4) +
  labs(title = "P4  Pathway activity, TP53-mut vs TP53-WT",
       subtitle = paste("Per-cell UCell scores averaged per sample. One line per matched pair.",
                        "Crossing lines are the result."),
       x = NULL, y = "UCell score (sample mean)",
       caption = paste("NONE of these pathways differs. Across 500 alternative valid control sets",
                       "the best of the 14 reached p < 0.05\nin 5.8% of them -- the chance rate --",
                       "and 8 of 14 never reached it in any draw. The reason is gene-set",
                       "membership,\nnot biology: IL11, C1QB, C1QC, LYVE1, MAF and MRC1 belong to",
                       "NONE of these sets, and all three TGF-beta\nsets contain none of the genes",
                       "that did move in P2/P3. A 200-gene rank score cannot see a 5-gene shift.")) +
  theme_lcc()
save_fig(fP4, "P4_pathway_activity", 13, 8)

## ===== P5. CCC, including the nodes that do not exist =======================================
message("[P5] CCC and the node inventory")
inv <- fread(file.path(LCC_TAB_DIR, "10_stromal_denovo_sample.csv"))
inv <- merge(arms, inv, by = c("dataset", "sample"), all.x = TRUE)
for (j in c("n_stromal_denovo", "n_MSC_fibroblast", "n_adipocyte", "n_endothelial", "n_pericyte",
            "n_osteolineage")) if (j %in% names(inv)) inv[is.na(get(j)), (j) := 0L]
myl <- fread(file.path(LCC_TAB_DIR, "04_myeloid_gate.csv"))
inv <- merge(inv, myl[, .(dataset, sample, n_cells_total, n_mono_dc, n_macrophage_like)],
             by = c("dataset", "sample"), all.x = TRUE)
bins <- det[stratum %in% c("HSC_MPP", "malignant") & gene == "COL1A1",
            .(dataset, sample, stratum, n = n_cells)]
bw <- dcast(bins, dataset + sample ~ stratum, value.var = "n")
inv <- merge(inv, bw, by = c("dataset", "sample"), all.x = TRUE)
NODE <- c(`Leukaemia (malignant)` = "malignant", `HSC / MPP` = "HSC_MPP",
          `Mono_DC (macrophage compartment)` = "n_mono_dc", `Macrophage-like` = "n_macrophage_like",
          `Niche / stromal (any)` = "n_stromal_denovo", `MSC / fibroblast` = "n_MSC_fibroblast",
          `Adipocyte` = "n_adipocyte", `Endothelial` = "n_endothelial")
NODE <- NODE[NODE %in% names(inv)]
p5 <- rbindlist(lapply(names(NODE), function(k)
  inv[, .(node = k, arm, sample, n = as.numeric(get(NODE[[k]])))]))
p5[is.na(n), n := 0]
p5s <- p5[, .(total = sum(n), n_samples_with = sum(n >= 10)), by = .(node, arm)]
p5s[, node := factor(node, levels = rev(names(NODE)))]
fwrite_safe(p5s, file.path(LCC_TAB_DIR, "P5_ccc_node_inventory.csv"))
a <- ggplot(p5s, aes(pmax(total, 0.5), node, fill = arm)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = comma(total)), position = position_dodge(width = 0.72),
            hjust = -0.12, size = 4.2, colour = INK[["secondary"]]) +
  scale_fill_manual(values = unname(PAL[c("orange", "blue")]), name = NULL) +
  scale_x_log10(labels = comma, expand = expansion(mult = c(0, 0.25))) +
  labs(title = "P5a  Which CCC nodes actually exist in these 18 samples",
       x = "cells (log scale)", y = NULL,
       subtitle = "The PI asked for leukaemia / HSC / macrophage / niche / MSC / adipocyte. Three of those are not present.") +
  theme_lcc() + theme(legend.position = "right")
et <- fread(file.path(LCC_TAB_DIR, "14_ccc_edge_tests.csv"))[norm == "relative" & n_pairs >= 4]
et[, edge := paste(source, "->", target)]
setorder(et, p_higher); et[, edge := factor(edge, levels = rev(edge))]
b <- ggplot(et, aes(frac_pairs_up * n_pairs / n_pairs, edge)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.4, colour = INK[["muted"]]) +
  geom_point(aes(size = n_pairs, colour = p_higher <= 0.05)) +
  scale_colour_manual(values = c(`FALSE` = INK[["muted"]], `TRUE` = PAL[["orange"]]),
                      labels = c("not significant", "p < 0.05"), name = NULL) +
  scale_size_continuous(range = c(1.5, 4), name = "usable pairs", breaks = c(4, 5, 6, 7)) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "P5b  Signalling between the compartments that do exist",
       x = "share of usable pairs in which TP53-mut is higher", y = NULL,
       subtitle = "CellChat, significant communication probability, relative normalisation.") +
  theme_lcc() + theme(legend.position = "right")
fP5 <- (a / b) + plot_layout(heights = c(1, 1.5)) 
save_fig(fP5, "P5_ccc_nodes_and_edges", 12, 11)

## ===== P6. Nectin-4 and the whole nectin / PVR receptor axis ================================
message("[P6] NECTIN4 / PVRL4 and the nectin-PVR axis")
NEC <- c("NECTIN4","PVRL4","NECTIN1","PVRL1","NECTIN2","PVRL2","NECTIN3","PVRL3",
         "PVR","TIGIT","CD226")
nd <- det[stratum == "all" & gene %in% NEC & gene_present == TRUE,
          .(dataset, sample, gene, pct = pct_nonzero, n_cells, n_nonzero, mean_lognorm)]
ALI <- c(PVRL4 = "NECTIN4", PVRL1 = "NECTIN1", PVRL2 = "NECTIN2", PVRL3 = "NECTIN3")
nd[gene %in% names(ALI), gene := ALI[gene]]
nd <- merge(nd, arms, by = c("dataset", "sample"))
nd <- nd[, .(pct = max(pct), expr = max(mean_lognorm), n_cells = n_cells[1]),
         by = .(dataset, sample, arm, pair_id, gene)]
GORD <- c("NECTIN4", "NECTIN1", "NECTIN3", "NECTIN2", "PVR", "CD226", "TIGIT")
nd <- nd[gene %in% GORD]
nd[, role := fifelse(gene %in% c("NECTIN4","NECTIN1","NECTIN2","NECTIN3","PVR"),
                     "ligand side (nectins / PVR)", "receptor side (TIGIT / CD226)")]
nd[, gene := factor(gene, GORD)]
fwrite_safe(nd, file.path(LCC_TAB_DIR, "P6_nectin_per_sample.csv"))

# per-gene paired test, so the receptor side is reported with numbers even where it is not significant
p6t <- dcast(nd, gene + role + pair_id ~ arm, value.var = "pct")
setnames(p6t, c("TP53-mut", "TP53-WT"), c("m", "w"))
p6stat <- p6t[, { ok <- !is.na(m) & !is.na(w); d <- m[ok] - w[ok]; d2 <- d[d != 0]
                  k <- length(d2)
                  p <- NA_real_
                  if (k >= 4L) { r <- rank(abs(d2)); wp <- sum(r[d2 > 0]); mx <- sum(r)
                    ways <- numeric(mx + 1L); ways[1L] <- 1
                    for (i in seq_len(k)) ways <- ways + c(numeric(r[i]), ways[seq_len(mx + 1L - r[i])])
                    p <- sum(ways[seq(floor(wp) + 1L, mx + 1L)]) / sum(ways) }
                  .(n_pairs = sum(ok), n_up = sum(m[ok] > w[ok]), p_higher = p,
                    med_mut = median(m[ok]), med_wt = median(w[ok])) }, by = .(gene, role)]
fwrite_safe(p6stat, file.path(LCC_TAB_DIR, "P6_nectin_axis_tests.csv"))
print(p6stat[order(gene)])

a6 <- ggplot(nd, aes(gene, pct + 0.01, colour = arm)) +
  geom_hline(yintercept = 1.01, linetype = "22", linewidth = 0.5, colour = INK[["muted"]]) +
  geom_point(position = position_jitterdodge(jitter.width = 0.16, dodge.width = 0.6),
             size = 2.6, alpha = 0.9, stroke = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.45,
               position = position_dodge(width = 0.6), show.legend = FALSE) +
  facet_grid(~ role, scales = "free_x", space = "free_x") +
  scale_colour_manual(values = unname(PAL[c("orange", "blue")]), name = NULL) +
  scale_y_log10(breaks = c(0, 0.1, 1, 10) + 0.01, labels = c("0", "0.1%", "1%", "10%")) +
  labs(title = "P6a  Nectin / PVR axis expression, per sample",
       subtitle = "18 samples. Dashed line = 1% of cells.",
       x = NULL, y = "% of cells with >= 1 transcript") +
  theme_lcc()

nd[, arm_s := fifelse(arm == "TP53-mut", "mut", "WT")]
b6 <- ggplot(nd, aes(arm_s, pct + 0.01, group = pair_id)) +
  geom_line(linewidth = 0.7, alpha = 0.55, colour = PAL[["blue"]]) +
  geom_point(size = 2.2, stroke = 0, alpha = 0.8, colour = PAL[["blue"]]) +
  geom_text(data = p6stat, aes(x = 1.5, y = 30,
                               label = sprintf("%d/%d up\np = %.2f", n_up, n_pairs, p_higher)),
            size = 4.6, colour = INK[["primary"]], inherit.aes = FALSE, lineheight = 0.95) +
  facet_wrap(~ gene, nrow = 2) +
  scale_y_log10(breaks = c(0, 0.1, 1, 10) + 0.01, labels = c("0", "0.1%", "1%", "10%"),
                limits = c(0.01, 60)) +
  labs(title = "P6b  The same, pair by pair",
       subtitle = "One line per matched pair; label gives pairs higher in TP53-mut and the one-sided paired p.",
       x = NULL, y = "% of cells expressing") +
  theme_lcc()

fP6 <- a6 / b6 + plot_layout(heights = c(1, 1))
save_fig(fP6, "P6_nectin4", 15, 12)

message("[done] 15_pi_panels")
