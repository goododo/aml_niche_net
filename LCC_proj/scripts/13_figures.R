#!/usr/bin/env Rscript
# 13_figures.R ----
# Figures for the TP53-mut vs TP53-WT patch, INCLUDING the negative results. The nulls are plotted
# on the same axes as the hits on purpose: a reader who sees only the surviving genes cannot judge
# them, and "14 pathways, none of which moved" is a finding that needs to be seen, not asserted.
#
# WHAT EACH FIGURE IS FOR
#   F1 effect vs robustness   every panel gene at once -- x = effect, y = how often it survives a
#                             re-draw of the control set. Null genes pile at the origin, so the
#                             handful that leave it are visible in context.
#   F2 robustness decay       the same genes across all four designs. This is the honest figure:
#                             it shows the evidence weakening when AML328's 4 timepoints stop
#                             counting as 4 patients, and again when tier C dilutes the arm.
#   F3 pathway null           all 14 pathway scores, all four designs, against the 5% chance line.
#   F4 macrophage shift       the direction split that makes this a phenotype change, not a count
#                             change -- up and down markers on one diverging axis.
#   F5 stromal reality        why the leukaemia-stroma question could not be answered.
#
# COLOUR: categorical hues assigned in fixed order and never cycled; scatter panels cap at 3 hues
# (the all-pairs limit); F4 is diverging (two poles + neutral). Text never wears a series colour.
#
# INPUT  : LCC_TAB_DIR/11_gene_results.csv, 12_robustness_{genes,pathways}.csv,
#          10_stromal_denovo_sample.csv, 10_niche_composition.csv
# OUTPUT : LCC_FIG_DIR/F1_effect_robustness.png ... F5_stromal_reality.png  (+ .pdf)
# Usage  : Rscript LCC_proj/scripts/13_figures.R

suppressPackageStartupMessages({
  library(data.table); library(here); library(ggplot2); library(ggrepel); library(scales)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
dir.create(LCC_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

## -- validated categorical palette (light mode), assigned in fixed slot order --
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))

message("[1] loading results")
gr  <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))[stratum == "all"]
rb  <- fread(file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))
rp  <- fread(file.path(LCC_TAB_DIR, "12_robustness_pathways.csv"))
DES <- c(A_sample_strict = "A  9 samples\nA+B tiers",
         B_patient_strict= "B  6 patients\nA+B tiers",
         C_sample_broad  = "C  14 samples\n+tier C",
         D_patient_broad = "D  11 patients\n+tier C")

## -- F1. effect size vs robustness, every panel gene ----
message("[2] F1 effect vs robustness")
d1 <- merge(gr[, .(gene, category, odds_ratio, pct_mut, pct_WT,
                   eff = median_delta_log2cpm, pairs_up = frac_pairs_mut_higher)],
            rb[design == "A_sample_strict", .(gene = feature, frac_p05, median_p)], by = "gene")
drop1 <- d1[is.na(eff) | is.na(frac_p05)]
if (nrow(drop1)) message("    [note] ", nrow(drop1), " genes not plottable (",
                         paste(drop1$gene, collapse = ", "),
                         ") -- deprecated aliases whose canonical symbol carries the data")
d1 <- d1[!is.na(eff) & !is.na(frac_p05)]
# Three colour slots only: this is a scatter, so the all-pairs limit applies. Categories fold into
# the two the hypothesis is about plus "other" -- folding, never a generated 4th hue.
d1[, grp := fifelse(category == "macrophage_marker", "Macrophage program",
             fifelse(category %in% c("profibrotic_cytokine_gf", "collagen_ecm_structural",
                                     "ecm_remodeling", "megakaryocyte_fibrosis_axis"),
                     "Fibrosis / ECM axis", "Other panel genes"))]
d1[, grp := factor(grp, levels = c("Macrophage program", "Fibrosis / ECM axis", "Other panel genes"))]
lab1 <- d1[frac_p05 >= 0.35]
f1 <- ggplot(d1, aes(eff, frac_p05)) +
  geom_hline(yintercept = 0.05, linewidth = 0.4, linetype = "22", colour = INK[["muted"]]) +
  geom_vline(xintercept = 0,    linewidth = 0.4, colour = INK[["muted"]]) +
  geom_point(aes(colour = grp), size = 2.2, alpha = 0.85, stroke = 0) +
  geom_text_repel(data = lab1, aes(label = gene, colour = grp), size = 4.2, seed = SEED,
                  max.overlaps = 30, min.segment.length = 0.2, show.legend = FALSE,
                  segment.colour = INK[["muted"]], segment.size = 0.3) +
  scale_colour_manual(values = unname(PAL[c("blue", "orange", "aqua")]), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Effect size vs robustness to the control-set choice",
       subtitle = "Every gene in the 160-gene panel. Genes at the origin are the null result.",
       x = "median paired difference in log2 CPM  (TP53-mut - TP53-WT, across the 9 pairs)",
       y = "share of 500 valid matched control sets reaching p < 0.05",
       caption = paste("Dashed line = 5%, the chance level. Design A: 9 matched sample pairs,",
                       "tiers A+B; one-sided exact paired\nsigned-rank, matched within dataset and",
                       "timepoint. The x axis is deliberately the SAMPLE-level effect, not the",
                       "\ncell-level odds ratio: the latter is cell-weighted, so a single large",
                       "sample can set it, and the two disagree\nin sign for 51 of 145 panel genes.",
                       "Plotting a cell-level x against a sample-level y put genes such as VSIG4",
                       "\non the wrong side of zero.")) +
  theme_lcc()
save_fig(f1, "F1_effect_robustness", 9, 6.5)

## -- F2. the honest figure: robustness decays across designs ----
message("[3] F2 robustness decay")
HITS <- c("C1QB", "SLC2A1", "SPP1", "LYVE1", "C1QC", "GP9", "OSM", "PGK1")   # 8 = the palette cap
d2 <- rb[feature %in% HITS]
# The four designs are a 2x2 (unit of analysis x mut definition), NOT a single ordered axis. An
# earlier version strung them A-B-C-D on one axis, which produced a V shape that was pure ordering
# artefact. Unit on x, definition as facets: each line now shows one real comparison.
d2[, unit := factor(fifelse(grepl("sample", design), "per SAMPLE\n(AML328 x4)",
                                                    "per PATIENT"),
                    levels = c("per SAMPLE\n(AML328 x4)", "per PATIENT"))]
d2[, defn := factor(fifelse(grepl("strict", design), "tiers A+B only  (6 patients)",
                                                     "+ tier C  (11 patients)"),
                    levels = c("tiers A+B only  (6 patients)", "+ tier C  (11 patients)"))]
d2[, feature := factor(feature, levels = HITS)]
f2 <- ggplot(d2, aes(unit, frac_p05, group = feature, colour = feature)) +
  geom_hline(yintercept = 0.05, linewidth = 0.4, linetype = "22", colour = INK[["muted"]]) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_point(size = 2.6, stroke = 0) +
  facet_wrap(~ defn) +
  scale_colour_manual(values = unname(PAL), name = NULL, guide = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "Every hit weakens once one patient stops counting as four",
       subtitle = paste("Each line is a gene. Left point: samples are the unit, so AML328's four",
                        "timepoints count four times.\nRight point: patients are the unit."),
       x = NULL, y = "share of 500 control sets reaching p < 0.05",
       caption = paste("Both panels fall from left to right, so the sample-level result was partly",
                       "AML328 counted repeatedly.\nAdding tier C (right panel) does not repair it:",
                       "tier C is expression-CNV only, and expression-CNV scored 0/3\nagainst genotype",
                       "truth, so most of what it adds is expected to be mislabelled. At 6 pairs the",
                       "smallest\nattainable p is 1/64, so part of the fall is power rather than",
                       "error -- the two cannot be fully separated here.")) +
  theme_lcc() + theme(panel.spacing.x = unit(1.4, "lines"))
save_fig(f2, "F2_robustness_decay", 10, 6.5)

## -- F3. the pathway null ----
message("[4] F3 pathway null")
d3 <- copy(rp)
# Full MSigDB names are 30-45 characters and were eating half the panel width. The collection
# prefix is dropped (it is in FIGURE_LEGENDS.md) and underscores become spaces.
d3[, pathway := sub("_UCell$", "", feature)]
d3[, pathway := gsub("_", " ", sub("^(HALLMARK|REACTOME|KEGG)_", "", pathway))]
SHORT <- c(A_sample_strict = "A  9 samples", B_patient_strict = "B  6 patients",
           C_sample_broad = "C  14 samples", D_patient_broad = "D  11 patients")
d3[, design := factor(design, levels = names(SHORT), labels = SHORT)]
ord <- d3[design == levels(d3$design)[1]][order(frac_p05)]$pathway
d3[, pathway := factor(pathway, levels = ord)]
f3 <- ggplot(d3, aes(frac_p05, pathway)) +
  geom_vline(xintercept = 0.05, linewidth = 0.4, linetype = "22", colour = INK[["muted"]]) +
  geom_segment(aes(x = 0, xend = frac_p05, yend = pathway), linewidth = 2,
               colour = PAL[["blue"]], alpha = 0.85, lineend = "round") +
  facet_wrap(~ design, nrow = 1) +
  # NO x limit. An earlier version capped this axis at 20% and silently clipped the three largest
  # bars in design D -- on a figure whose entire job is to show a negative result, that is the one
  # unacceptable error. The axis now runs to the data.
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = c(0, 0.15, 0.3),
                     expand = expansion(mult = c(0, 0.06))) +
  labs(title = "No pathway score holds up across designs",
       subtitle = paste("Dashed line = 5%, the chance level. A pathway that were real would sit",
                        "right of it in EVERY panel;\nnone does -- the ones that rise in the",
                        "patient-collapsed designs fall back in the sample-level ones."),
       x = "share of 500 valid matched control sets reaching p < 0.05", y = NULL,
       caption = paste("TNFA_SIGNALING runs 0.2% -> 7.8% -> 3.0% -> 35.4% across A/B/C/D; hypoxia",
                       "0.2% -> 12.6% -> 0.6% -> 27.0%.\nA finding that reverses with the analysis",
                       "unit is not a finding. Gene-set membership explains the floor:\nIL11, C1QB,",
                       "C1QC, LYVE1, MAF and MRC1 belong to NONE of these 14 sets, and all three",
                       "TGF-beta sets\ncontain none of the significant genes -- the scores cannot",
                       "see the signal the panel genes carry.")) +
  theme_lcc() + theme(panel.spacing.x = unit(1, "lines"))
save_fig(f3, "F3_pathway_null", 13, 6)

## -- F4. macrophage phenotype shift, not a count shift ----
message("[5] F4 macrophage shift")
# Two effect measures side by side, because they DISAGREE and the disagreement is the lesson.
d4 <- gr[category == "macrophage_marker" & !is.na(frac_pairs_mut_higher)]
d4[, n_up := round(frac_pairs_mut_higher * n_pairs)]
setorder(d4, n_up); d4[, gene := factor(gene, levels = gene)]
d4l <- melt(d4[, .(gene, `paired (9 pairs, equal weight)` = frac_pairs_mut_higher - 0.5,
                   `pooled cells (cell-weighted)` = pmax(pmin(log2(odds_ratio) / 8, 0.5), -0.5))],
            id.vars = "gene", variable.name = "measure", value.name = "v")
f4 <- ggplot(d4l, aes(v, gene, fill = measure)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = INK[["muted"]]) +
  geom_col(width = 0.62, position = position_dodge(width = 0.68)) +
  scale_fill_manual(values = unname(PAL[c("blue", "orange")]), name = NULL) +
  scale_x_continuous(breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = c("strongly\nlower", "lower", "no\ndifference", "higher",
                                "strongly\nhigher"),
                     expand = expansion(mult = 0.05)) +
  labs(title = "The two effect measures disagree, and only one matches the test",
       subtitle = paste("Blue = share of the 9 matched pairs in which TP53-mut is higher (centred",
                        "at 0.5).\nOrange = the cell-pooled odds ratio, rescaled onto the same",
                        "axis."),
       x = NULL, y = NULL,
       caption = paste("VSIG4 and CD163 read as strongly LOWER when cells are pooled (OR 0.23,",
                       "0.28) yet are higher in 6 of 9\npairs -- the pooled fraction is",
                       "cell-weighted, so a single 15,000-cell control sample sets it. On the",
                       "measure\nthat matches the test, macrophage markers trend UP broadly",
                       "(MAF 8/9, C1QB 7/9, MRC1 7/9) with STAB1\nthe one exception (2/9). This",
                       "is NOT a polarisation switch. Macrophage frequency itself does not",
                       "differ\nsignificantly (0.36% vs 0.46%, p = 0.32).")) +
  theme_lcc()
save_fig(f4, "F4_macrophage_shift", 9.5, 6)

## -- F5. why the leukaemia-stroma question could not be answered ----
message("[6] F5 stromal reality")
sd <- fread(file.path(LCC_TAB_DIR, "10_stromal_denovo_sample.csv"))
d5 <- sd[, .(n_cells = sum(n_cells), n_str = sum(n_stromal_denovo), n_samp = .N), by = library_type]
d5[, pct := 100 * n_str / n_cells]
LAB <- c(stroma_enriched = "Stroma-enriched\n(Chen2023 niche, GSE253355)",
         whole_MNC = "Whole marrow aspirate\n(the AML samples)",
         sorted_CD34 = "CD34-sorted\n(false-positive floor)")
d5[, lt := factor(LAB[library_type], levels = LAB[c("stroma_enriched", "whole_MNC", "sorted_CD34")])]
fp <- d5[library_type == "sorted_CD34"]$pct
# Lollipops, not bars. geom_col on a log axis draws from y = 0, which is -Inf in log space, so the
# bars get clipped to the panel floor and read as if they START at some non-zero value -- the 0.028%
# bar looked like it spanned 0.03% to 1%. A segment with an explicit baseline cannot lie that way.
Y0 <- 0.005
f5 <- ggplot(d5, aes(lt, pct, colour = lt)) +
  geom_hline(yintercept = fp, linetype = "22", linewidth = 0.4, colour = INK[["muted"]]) +
  geom_segment(aes(xend = lt, y = Y0, yend = pct), linewidth = 2.4, lineend = "round") +
  geom_point(size = 5, stroke = 0) +
  geom_text(aes(label = sprintf("%s%%\n%s of %s cells\n%d samples",
                                formatC(pct, format = "g", digits = 2),
                                comma(n_str), comma(n_cells), n_samp)),
            vjust = -0.6, size = 4.3, colour = INK[["secondary"]], lineheight = 1.05) +
  scale_y_log10(labels = function(x) paste0(x, "%"),
                breaks = c(0.01, 0.1, 1, 10, 100), limits = c(Y0, 300)) +
  scale_colour_manual(values = unname(PAL[c("aqua", "orange", "blue")]), guide = "none") +
  labs(title = "Marrow aspirates do not contain stroma",
       subtitle = paste("De-novo marker co-expression per cell, bypassing the projection label.",
                        "Log scale."),
       x = NULL, y = "stromal cells (% of all cells, log scale)",
       caption = paste("Threshold calibrated so the CD34-sorted libraries -- which contain no stroma",
                       "by construction -- give\n0.85 false positives per 10,000 cells (dashed line).",
                       "Whole-marrow AML samples sit at 3x that floor:\n56 stromal cells in 200,096",
                       "cells across 36 samples. No TP53-mut sample has a usable stromal\npopulation,",
                       "so the leukaemia-stroma interaction question has n = 0 on the mutant side.")) +
  theme_lcc()
save_fig(f5, "F5_stromal_reality", 8, 6)

message("[done] 13_figures -> ", LCC_FIG_DIR)
