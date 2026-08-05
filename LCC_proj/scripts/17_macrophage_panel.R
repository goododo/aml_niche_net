#!/usr/bin/env Rscript
# 17_macrophage_panel.R ----
# P7: the macrophage result, shown rather than tabulated.
#
# WHY: the C1Q macrophage finding is the one result that cross-validates against the Kao lab's DSP
# data, and until now it existed only as a table of "7/9" counts in the write-up. A count is an
# assertion; the reader cannot see whether it rests on one extreme pair or on nine consistent ones.
# This figure plots every pair, so a single driving pair would be visible.
#
# FORM: paired slope plot. The data's job is POLARITY (did this pair go up or down), so the colour
# is a diverging pair from the project palette -- orange = higher in the mutant sample, blue = lower
# -- with grey reserved for ties. Free y per facet on purpose: these genes differ ~100-fold in
# absolute expression, and the question is direction within a pair, not magnitude between genes.
#
# ALL TEN macrophage-panel genes are shown, including the two that go the other way (STAB1 2/9,
# SELENOP 4/9). Showing only the six that moved would misrepresent the panel.
#
# INPUT  : LCC_TAB_DIR/11_matched_design.csv     the 9 pairs
#          LCC_TAB_DIR/04_detection_by_sample.csv  per-sample per-gene mean lognorm
#          LCC_TAB_DIR/11_gene_results.csv       pairs-higher and the paired p
#          LCC_TAB_DIR/12_robustness_genes.csv   share of 500 alternative matchings reaching p<0.05
# OUTPUT : LCC_FIG_DIR/P7_macrophage_pairs.{png,pdf}
#          LCC_TAB_DIR/P7_macrophage_pairs.csv
suppressPackageStartupMessages({
  library(data.table); library(here); library(ggplot2); library(patchwork)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
set.seed(SEED)

GENES <- c("C1QA", "C1QB", "C1QC", "MAF", "MRC1", "LYVE1", "VSIG4", "CD163", "SELENOP", "STAB1")

## -- data ----
des <- fread(file.path(LCC_TAB_DIR, "11_matched_design.csv"))
arms <- rbind(des[, .(sample = mut_sample, dataset, pair_id, arm = "TP53-mut")],
              des[, .(sample = wt_sample,  dataset, pair_id, arm = "TP53-WT")])

det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"),
             select = c("gene", "stratum", "sample", "dataset", "mean_lognorm", "gene_present"))
d <- det[stratum == "all" & gene %in% GENES & gene_present == TRUE]
d <- merge(d, arms, by = c("dataset", "sample"))
# one row per (gene, pair, arm); 03 can emit a gene more than once per sample across strata slices
d <- d[, .(mean_lognorm = max(mean_lognorm, na.rm = TRUE)), by = .(gene, pair_id, arm, sample)]

st <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))[stratum == "all" & gene %in% GENES,
        .(gene, n_pairs, frac_pairs_mut_higher, p_sample_higher)]
st[, pairs_higher := round(frac_pairs_mut_higher * n_pairs)]
st[, lab := sprintf("%d/%d   p = %.3f", pairs_higher, n_pairs, p_sample_higher)]
rb <- fread(file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))[
        design == "A_sample_strict" & kind == "gene" & feature %in% GENES, .(gene = feature, frac_p05)]
st <- merge(st, rb, by = "gene", all.x = TRUE)

# facet order: strongest result first, so the eye reads down a gradient rather than an arbitrary list
setorder(st, -pairs_higher, p_sample_higher)
d[,  gene := factor(gene, levels = st$gene)]
st[, gene := factor(gene, levels = st$gene)]

# per-pair direction drives the colour
w <- dcast(d, gene + pair_id ~ arm, value.var = "mean_lognorm")
setnames(w, c("TP53-mut", "TP53-WT"), c("v_mut", "v_wt"))
w[, dir := fifelse(v_mut > v_wt, "higher in TP53-mut",
            fifelse(v_mut < v_wt, "lower in TP53-mut", "no change"))]
d <- merge(d, w[, .(gene, pair_id, dir)], by = c("gene", "pair_id"))
d[, arm := factor(arm, levels = c("TP53-WT", "TP53-mut"))]
fwrite(merge(w, st[, .(gene, pairs_higher, n_pairs, p_sample_higher, frac_p05)], by = "gene"),
       file.path(LCC_TAB_DIR, "P7_macrophage_pairs.csv"))

DIRCOL <- c("higher in TP53-mut" = PAL[["orange"]],
            "lower in TP53-mut"  = PAL[["blue"]],
            "no change"          = GRID)

## -- (a) every pair, every gene ----
# label position: top of each facet's own range, since y is free
lab_pos <- d[, .(y = max(mean_lognorm) * 1.18 + 1e-4), by = gene]
lab_pos <- merge(lab_pos, st[, .(gene, lab)], by = "gene")

pa <- ggplot(d, aes(arm, mean_lognorm, group = pair_id, colour = dir)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.6) +
  geom_text(data = lab_pos, aes(x = 1.5, y = y, label = lab),
            inherit.aes = FALSE, size = 4.3, colour = INK[["primary"]], vjust = 1) +
  scale_colour_manual(values = DIRCOL, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  facet_wrap(~ gene, nrow = 2, scales = "free_y") +
  labs(title = "P7a  Macrophage genes, all 9 matched pairs",
       subtitle = "each line is one matched pair",
       x = NULL, y = "mean expression per sample (log-normalised)") +
  theme_lcc() +
  theme(panel.grid.major.x = element_blank())

## -- (b) how often the result survives a different control set ----
# The 500 re-matchings answer "would another analyst, matching equally correctly, have seen this?"
# It is NOT a p-value. 5% is the chance level and is drawn so the reader cannot forget it.
pb <- ggplot(st, aes(reorder(gene, frac_p05), frac_p05)) +
  geom_segment(aes(xend = gene, y = 0, yend = frac_p05), colour = GRID, linewidth = 1.4) +
  geom_point(size = 4.2, colour = PAL[["orange"]]) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = INK[["muted"]]) +
  annotate("text", x = 1.2, y = 0.09, label = "5% = chance", hjust = 0,
           size = 4.3, colour = INK[["muted"]]) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "P7b  Share of 500 alternative control sets reaching p < 0.05",
       subtitle = "robustness to how the wild-type arm was matched",
       x = NULL, y = NULL) +
  theme_lcc()

save_fig(pa / pb + patchwork::plot_layout(heights = c(1.5, 1)),
         "P7_macrophage_pairs", 13, 12)

cat("\n-- P7 source numbers --\n")
print(st[, .(gene, pairs_higher, n_pairs, p = round(p_sample_higher, 3), robustness = frac_p05)])
message("[done]")
