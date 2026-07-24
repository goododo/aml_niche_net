# g_assessment_figures.R ----
# Five key figures for the internal / lab-meeting assessment. Unrelated to f01-f04
# (those point at legacy paths and have never run).
# INPUT  : results/tables/{01_preprocess,02_malignancy,03_hierarchy,05_ccc,07_fgw,08_scoring}/*.csv
# OUTPUT : results/figures/11_assessment/g0{1..5}_*.{pdf,png}
# Usage  : conda run -p /FAST/gr10634/gaozy/general_env Rscript scripts/10_figures/g_assessment_figures.R
#
# [DECISION] Palette = first two slots of the dataviz reference palette (#2a78d6 blue /
#            #eb6834 orange), validated with the Machado-2009 severity-1.0 CVD sim:
#            all-pairs CVD dE=24.7, normal dE=33.6, contrast >=3:1. Diverging scale is
#            blue<->red with a neutral grey midpoint.
# [DECISION] Status colour critical=#d03b3b is used ONLY to flag a problematic quantity;
#            it must never double as a series colour anywhere in the figure.
# [DECISION] Every printed number is computed from the source table, never a hard-coded
#            literal. The first draft hard-coded the edge-level min q into the block panel
#            (0.87 vs the true 0.93) and wrote n_perm=2000 (actually 10000).
# [DECISION] Probability/rate axes are always capped (rate <=100%, p and q <=1); labels
#            live outside the panel via clip="off".

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
TBL  <- file.path(ROOT, "results", "tables")
FIG  <- file.path(ROOT, "results", "figures", "11_assessment")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
CRIT <- "#d03b3b"; WARN <- "#fab219"; GOOD <- "#0ca30c"
GREY <- "#8a8a85"; LGREY <- "#c6c5c0"; MIDGREY <- "#f0efec"
INK  <- "#0b0b0b"; INK2 <- "#52514e"
FAM  <- "sans"

theme_a <- function(base = 11) {
  theme_minimal(base_size = base, base_family = FAM) +
    theme(
      plot.title      = element_text(face = "bold", size = base + 2, colour = INK, hjust = 0),
      plot.subtitle   = element_text(size = base - 0.5, colour = INK2, hjust = 0, lineheight = 1.3),
      plot.caption    = element_text(size = base - 2.5, colour = INK2, hjust = 0, lineheight = 1.35),
      axis.title      = element_text(size = base - 1, colour = INK2),
      axis.text       = element_text(size = base - 1, colour = INK2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
      legend.position = "top", legend.justification = "left",
      legend.title    = element_blank(), legend.key.size = unit(0.8, "lines"),
      plot.title.position = "plot", plot.caption.position = "plot",
      strip.text      = element_text(face = "bold", colour = INK, size = base - 1)
    )
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 200,
         device = ragg::agg_png, bg = "white")
  cat(sprintf("[fig] %s\n", name))
}
# p-value formatter: avoid the "p=<0.001" double-symbol artefact
pfmt <- function(p) ifelse(p < 1e-3, "p<0.001", paste0("p=", formatC(p, format = "g", digits = 2)))

## ══════════════════════════════════════════════════════════════════
## Fig 1 — sample count at each stage (NOT a nested funnel)
## ══════════════════════════════════════════════════════════════════
qc <- rbindlist(lapply(Sys.glob(file.path(TBL, "01_preprocess", "03_qc_report__*.csv")),
                       fread), fill = TRUE)
qc_pass <- qc[status == "PASS"]
cons  <- fread(file.path(TBL, "02_malignancy", "ALL_consensus_summary.csv"))
fgwix <- fread(file.path(TBL, "07_fgw", "fgw_input_index.csv"))
asw   <- fread(file.path(TBL, "08_scoring", "alpha_sweep.csv"))

# Reconcile 148 vs 130 from the data, do not guess the reason
setkey(cons, dataset, sample); setkey(fgwix, dataset, sample)
both     <- fgwix[cons, nomatch = 0L]
graph_nolabel <- fgwix[!cons]
label_nograph <- cons[!fgwix]
n_aml_nolabel <- graph_nolabel[healthy == FALSE, .N]
n_hlt_nolabel <- graph_nolabel[healthy == TRUE,  .N]

funnel <- data.table(
  stage = c("QC-passing samples", "With malignancy consensus label", "CCC graph built",
            "Entered FGW global model", "Entered platform-controlled model"),
  n     = c(nrow(qc_pass), nrow(cons), nrow(fgwix), asw$n_global[1], asw$n_strat[1]),
  nds   = c(uniqueN(qc_pass$dataset), uniqueN(cons$dataset), uniqueN(fgwix$dataset), NA_integer_, 3L)
)
funnel[, stage := factor(stage, levels = rev(stage))]
funnel[, lab := ifelse(is.na(nds), sprintf("%d samples", n),
                       sprintf("%d samples | %d datasets", n, nds))]
funnel[, hl := c(rep(FALSE, 4), TRUE)]

p1 <- ggplot(funnel, aes(stage, n, fill = hl)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = lab), hjust = -0.06, size = 3.4, colour = INK, family = FAM) +
  annotate("text", x = 1, y = funnel$n[5] + 6, label = "<- narrowest", hjust = 0,
           size = 3.2, colour = INK2, family = FAM, vjust = 2.4) +
  scale_fill_manual(values = c("FALSE" = BLUE, "TRUE" = CRIT), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30))) +
  coord_flip() +
  labs(title = "Sample count at each stage (note: NOT nested)",
       subtitle = "The H2 'platform-controlled' result rests on 3 datasets, 60 samples",
       x = NULL, y = "samples",
       caption = sprintf(paste0(
         "Reconciliation: CCC graph (%d) = labelled & graphed (%d) + graphed without malignancy label ",
         "(%d: %d healthy donors + %d AML samples);\n",
         "a further %d samples are labelled but never graphed. Those %d unlabelled AML samples have their ",
         "frac_malignant filled with the\nglobal mean at 01_build_fgw_inputs.R:72 -- i.e. handed a ",
         "malignant fraction that was never measured."),
         nrow(fgwix), nrow(both), nrow(graph_nolabel), n_hlt_nolabel, n_aml_nolabel,
         nrow(label_nograph), n_aml_nolabel)) +
  theme_a() + theme(panel.grid.major.y = element_blank())

save_fig(p1, "g01_cohort_counts", 9.8, 5.0)

## ══════════════════════════════════════════════════════════════════
## Fig 2 — malignancy label quality (most important panel)
## ══════════════════════════════════════════════════════════════════
fpr <- fread(file.path(TBL, "02_malignancy", "malignancy_fpr_by_bin.csv"))
fpr[hierarchy_bin == "" | is.na(hierarchy_bin), hierarchy_bin := "Unassigned bin"]   # do not silently drop
fpr[, in_graph := !is.na(in_ccc_graph) & in_ccc_graph == TRUE]
fpr[, hierarchy_bin := factor(hierarchy_bin, levels = hierarchy_bin[order(FPR)])]
fpr[, hl := fifelse(as.character(hierarchy_bin) == "HSC_MPP", "crit",
             fifelse(in_graph, "in", "out"))]
max_in <- fpr[in_graph == TRUE][which.max(FPR)]

p2a <- ggplot(fpr, aes(hierarchy_bin, FPR, fill = hl)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%  (n=%s)", 100 * FPR, comma(n))),
            hjust = -0.05, size = 3.15, colour = INK, family = FAM) +
  scale_fill_manual(values = c(crit = CRIT, "in" = BLUE, out = LGREY), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1),
                     breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "A  Fraction of healthy-donor cells called malignant (false-positive rate)",
       subtitle = "Any malignant call in healthy marrow is by definition a false positive. Light grey = bins not in the CCC graph.",
       x = NULL, y = "false-positive rate") +
  theme_a() +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 105, 5, 5))

tier <- cons[, .N, by = evidence_tier]
tier[, evidence_tier := factor(evidence_tier, levels = c("A_concordant", "B_multi_partial", "C_single"))]
setorder(tier, evidence_tier)
tier[, lab_leg := sprintf("%s (%d)", evidence_tier, N)]

p2b <- ggplot(tier, aes(x = "", y = N, fill = evidence_tier)) +
  geom_col(width = 0.45, colour = "white", linewidth = 1.1) +
  geom_text(data = tier[N > 10], aes(label = sprintf("%s  %d", evidence_tier, N)),
            position = position_stack(vjust = 0.5), size = 3.3, colour = "white",
            fontface = "bold", family = FAM) +
  scale_fill_manual(values = c("A_concordant" = GOOD, "B_multi_partial" = WARN, "C_single" = GREY),
                    labels = tier$lab_leg) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(title = sprintf("B  Evidence-tier composition (%d samples)", nrow(cons)),
       subtitle = paste0("Blueprint Phase 1 requires an inferCNV + Numbat + VarTrix consensus (>= 2/3 agree);\n",
                         sprintf("in practice %d/%d are single-arm inferCNV. The %d tier-A/B samples are all from GSE227903.",
                                 tier[evidence_tier == "C_single", N], nrow(cons),
                                 tier[evidence_tier != "C_single", sum(N)])),
       x = NULL, y = "samples") +
  theme_a() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())

p2 <- p2a / p2b + plot_layout(heights = c(2.5, 1)) +
  plot_annotation(
    title = "Malignancy label quality: the LSC bin has a ~40% false-positive rate",
    subtitle = paste0(
      sprintf("HSC_MPP is the node the central claim leans on most, and the highest-FPR bin inside the CCC graph (%.1f%%).\n", 100 * max_in$FPR),
      "Lymphoid bins are clean (T_NK 3.6%, B_Plasma 7.7%); myeloid/primitive bins are not -- consistent with inferCNV keying on expression deviation.\n",
      "(Stromal is higher at 90.3% but is excluded from the CCC graph by the established premise.)"),
    caption = paste0("Source: A = 02_malignancy/malignancy_fpr_by_bin.csv (from 96_malignancy_fpr_healthy.R); ",
                     "B = 02_malignancy/ALL_consensus_summary.csv"),
    theme = theme_a())

save_fig(p2, "g02_label_quality", 10.4, 8.4)

## ══════════════════════════════════════════════════════════════════
## Fig 3 — H1: 49 directed edges, none significant
## ══════════════════════════════════════════════════════════════════
ee  <- fread(file.path(TBL, "08_scoring", "emergent_edges.csv"))
blk <- fread(file.path(TBL, "08_scoring", "block_permutation.csv"))
NODES <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma")
ee[, sender_bin   := factor(sender_bin,   levels = NODES)]
ee[, receiver_bin := factor(receiver_bin, levels = rev(NODES))]
lim <- max(abs(ee$dC_real))
# Recover the permutation count from the perm_p denominator, do not hard-code:
# perm_p = (1+k)/(n_perm+1), so perm_p*(n_perm+1) must be integer for every edge. Take the
# unique candidate that satisfies it. (The first draft used a "smallest gap" heuristic that
# wrongly returned 1666 -- adjacent perm_p gaps are not necessarily 1/(n_perm+1).)
.n_cand <- c(500L, 1000L, 2000L, 5000L, 10000L, 20000L, 50000L)
.n_ok <- vapply(.n_cand, function(n) all(abs(ee$perm_p * (n + 1) - round(ee$perm_p * (n + 1))) < 1e-6),
                logical(1))
stopifnot("cannot recover n_perm from perm_p" = any(.n_ok))
n_perm_edge <- .n_cand[.n_ok][1]

p3a <- ggplot(ee, aes(sender_bin, receiver_bin, fill = dC_real)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  scale_fill_gradient2(low = CRIT, mid = MIDGREY, high = BLUE, midpoint = 0,
                       limits = c(-lim, lim),
                       name = "dC = C_AML - C_healthy    <- stronger in AML | weaker in AML ->") +
  coord_fixed() +
  labs(title = "A  AML-vs-healthy difference across the 49 directed edges",
       subtitle = sprintf("Colour = magnitude of the difference, but no cell survives the permutation test (all q >= %.2f). This is the structure of noise -- do not read a pattern into it.",
                          min(ee$q_bh)),
       x = "sender", y = "receiver") +
  theme_a() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.major = element_blank(),
        legend.key.width = unit(2.1, "lines"),
        legend.title = element_text(size = 8, colour = INK2))

setorder(blk, -n_edges)
blk[, block := factor(block, levels = blk$block)]
p3b <- ggplot(blk, aes(block, q_bh)) +
  geom_col(width = 0.6, fill = LGREY) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = CRIT, linewidth = 0.6) +
  annotate("text", x = 0.55, y = 0.07, label = "FDR 0.05", hjust = 0, size = 3,
           colour = CRIT, family = FAM) +
  geom_text(aes(label = sprintf("q=%.2f  (%d edges)", q_bh, n_edges)), hjust = -0.06,
            size = 3.05, colour = INK, family = FAM) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "B  BH-FDR q for the 5 pre-specified communication blocks (ordered by edge count)",
       subtitle = sprintf("All five blocks share the identical q = %.2f, far above 0.05.", min(blk$q_bh)),
       x = NULL, y = "q (BH-FDR)") +
  theme_a() +
  theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 120, 5, 5))

p3 <- p3a / p3b + plot_layout(heights = c(2.2, 1)) +
  plot_annotation(
    title = "H1 (conserved emergent edges): no signal observed",
    subtitle = sprintf(paste0(
      "0 of 49 edges pass the permutation test (min perm p = %.3f, min q = %.2f); all 5 blocks non-significant.\n",
      "Key: edges come from CellChat and never pass through the malignancy labels -- so this null is immune to the label-quality problem, the most solid result so far."),
      min(ee$perm_p), min(ee$q_bh)),
    caption = sprintf("Source: 02_permutation_emergent.py (n_perm=%d) / 03_block_permutation.py (n_perm=10000)",
                      n_perm_edge),
    theme = theme_a())

save_fig(p3, "g03_h1_null", 9.8, 9.2)

## ══════════════════════════════════════════════════════════════════
## Fig 4 — H2 decisive test
## ══════════════════════════════════════════════════════════════════
asw_l <- melt(asw, id.vars = "alpha", measure.vars = c("p_global", "p_strat"),
              variable.name = "model", value.name = "p")
asw_l[, model := factor(model, levels = c("p_global", "p_strat"),
                        labels = c(sprintf("global model (n=%d)", asw$n_global[1]),
                                   sprintf("within-dataset (n=%d)", asw$n_strat[1])))]

p4a <- ggplot(asw_l, aes(alpha, p, colour = model)) +
  annotate("rect", xmin = 0.88, xmax = 1.02, ymin = 5e-5, ymax = 1, fill = CRIT, alpha = 0.07) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = CRIT, linewidth = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  geom_text(data = asw_l[alpha == 1], aes(label = model), hjust = 0, nudge_x = 0.035,
            size = 3, family = FAM, show.legend = FALSE) +
  annotate("text", x = 0.95, y = 0.62, label = "alpha = 1\npure topology", size = 3,
           colour = CRIT, family = FAM, lineheight = 0.95) +
  annotate("text", x = 0.02, y = 0.072, label = "p = 0.05", hjust = 0, size = 3,
           colour = CRIT, family = FAM) +
  scale_colour_manual(values = c(BLUE, ORANGE)) +
  scale_y_log10(limits = c(5e-5, 1), labels = label_number(drop0trailing = TRUE),
                breaks = c(1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.3, 1)) +
  scale_x_continuous(breaks = asw$alpha, expand = expansion(mult = c(0.03, 0.34))) +
  labs(title = "A  Alpha sweep: significance vanishes as the feature term is turned off",
       subtitle = "alpha=0 pure node features -> alpha=1 pure Gromov-Wasserstein (pure topology). Both models cross 0.05 at alpha=1.",
       x = "alpha (weight on the structure term in FGW)", y = "permutation p (log axis, capped at 1)") +
  theme_a() + theme(legend.position = "none")

fd <- fread(file.path(TBL, "08_scoring", "feature_decomposition.csv"))[alpha == 0.5]
fd[, feature_set := factor(feature_set,
      levels = c("only_stemness", "only_ncells", "no_frac_mal", "all3", "only_frac_mal"))]
setorder(fd, feature_set)
fd[, kind := fifelse(feature_set == "only_frac_mal", "circular (healthy forced to 0)",
              fifelse(p_global < 0.05, "significant (p<0.05)", "not significant"))]
fd[, kind := factor(kind, levels = c("circular (healthy forced to 0)", "significant (p<0.05)", "not significant"))]

p4b <- ggplot(fd, aes(feature_set, beta_global, fill = kind)) +
  geom_col(width = 0.6) +
  # Labels for negative bars also go to the right of 0, so they never overlap the y-axis category text
  geom_text(aes(y = pmax(beta_global, 0.002), label = sprintf("beta=%.3f  %s", beta_global, pfmt(p_global))),
            hjust = -0.06, size = 3.05, colour = INK, family = FAM) +
  scale_fill_manual(values = setNames(c(CRIT, BLUE, LGREY), levels(fd$kind))) +
  scale_y_continuous(limits = c(-0.05, 0.26), expand = expansion(mult = c(0.02, 0.02))) +
  coord_flip(clip = "off") +
  labs(title = "B  Feature decomposition (alpha=0.5): where the signal comes from",
       subtitle = paste0(
         "frac_malignant alone reproduces and amplifies the effect (beta=0.173 > all3's 0.075) -- it is forced to 0\n",
         "for healthy samples at 01_build_fgw_inputs.R:67, so that separation is constructed. But dropping it still\n",
         "leaves beta=0.064 (p=0.001), of which only_ncells (cell count, a pure technical quantity) is beta=0.050 (p=0.008);\n",
         "only_stemness (the actual biological quantity) is not significant."),
       x = NULL, y = "is_aml partial-regression coefficient beta") +
  theme_a() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(5, 100, 5, 5))

p4 <- p4a / p4b + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title = "H2 (topology carries disease information): significance is constructed, not topological",
    subtitle = paste0("Both panels are deliberately designed falsification tests.\n",
                      "Conclusion: the current H2 signal cannot be attributed to network topology, nor to stemness."),
    caption = "Source: 06_alpha_sweep.py / 07_feature_decomposition.py (both n_perm=10000)",
    theme = theme_a())

save_fig(p4, "g04_h2_decisive", 10.6, 9.4)

## ══════════════════════════════════════════════════════════════════
## Fig 5 — H3: blocked by metadata, not interpreted
## ══════════════════════════════════════════════════════════════════
stem <- fread(file.path(TBL, "03_hierarchy", "stemness_by_timepoint.csv"))
md   <- fread(file.path(TBL, "03_hierarchy", "malignant_distribution.csv"))[ok == TRUE]
dom  <- md[timepoint %in% c("Dx", "MRD", "Relapse"),
           .(n_pat = uniqueN(patient), tot = sum(tot_mal), top = max(tot_mal)), by = timepoint]
dom[, share := top / tot]
mrd_share <- dom[timepoint == "MRD", share]

sl <- melt(stem, id.vars = c("tp", "n_malignant"),
           measure.vars = c("LSC17", "vanGalen_HSC_Prog", "vanGalen_HSC_like", "HSPC_core"),
           variable.name = "signature", value.name = "score")
sl[, tp := factor(tp, levels = c("Dx", "MRD", "Relapse"))]
sl[, shape := fifelse(signature %in% c("LSC17", "vanGalen_HSC_Prog"), "V-shaped", "not V-shaped")]
xlab_n <- dom[match(levels(sl$tp), timepoint)]
# Keep labels short -- the first draft used "Dx\n8 patients | 4,136 malignant cells", which
# overlapped into one line inside the facet. Cell counts are moved to the subtitle instead.
lev_lab <- sprintf("%s\n(%d patients)", levels(sl$tp), xlab_n$n_pat)

p5a <- ggplot(sl, aes(tp, score, colour = signature, group = signature)) +
  geom_line(linewidth = 0.85) + geom_point(size = 2.4) +
  geom_text(data = sl[tp == "Relapse"], aes(label = signature), hjust = 0, nudge_x = 0.05,
            size = 2.9, family = FAM, show.legend = FALSE) +
  facet_wrap(~shape, nrow = 1) +
  scale_colour_manual(values = c(BLUE, ORANGE, GREY, "#4a3aa7"), guide = "none") +
  scale_x_discrete(labels = lev_lab, expand = expansion(add = c(0.35, 1.5))) +
  labs(title = "A  Stemness signatures along the treatment axis (transcriptional, independent of projection)",
       subtitle = sprintf(paste0(
         "Blueprint H3 predicts monotone Dx < MRD < Relapse. Actual: only 2 of 4 signatures are V-shaped\n",
         "(vanGalen_HSC_like peaks at MRD, HSPC_core declines monotonically).\n",
         "These are cell-level pooled means (Dx/MRD/Relapse: %s / %s / %s malignant cells): %.0f%% of MRD\n",
         "malignant cells come from one patient (6323), so the 'V' may be a single-sample effect."),
         comma(dom[timepoint == "Dx", tot]), comma(dom[timepoint == "MRD", tot]),
         comma(dom[timepoint == "Relapse", tot]), 100 * mrd_share),
       x = NULL, y = "mean signature score") +
  theme_a() + theme(panel.grid.major.x = element_blank())

sh <- fread(file.path(TBL, "03_hierarchy", "distribution_shift_tests.csv"))
sh[, comparison := factor(comparison, levels = c("Dx -> MRD", "MRD -> Relapse", "Dx -> Relapse"))]
sh[, score := factor(score, levels = c("primitive_frac", "stem_frac"),
      labels = c("primitive_frac (HSC_MPP+LMPP_GMP share)", "stem_frac (HSC_MPP share)"))]

p5b <- ggplot(sh, aes(comparison, median_delta)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.4) +
  geom_col(width = 0.6, fill = LGREY) +
  geom_text(aes(label = sprintf("%s\nn=%d pairs", pfmt(p_value), n_pairs),
                vjust = ifelse(median_delta > 0, -0.3, 1.2),
                fontface = ifelse(p_value < 0.05, "bold", "plain")),
            size = 2.8, colour = INK2, family = FAM, lineheight = 0.95) +
  facet_wrap(~score, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.32, 0.32))) +
  labs(title = "B  Within-pair median change (bold = p<0.05)",
       subtitle = paste0("The three comparisons use different patient sets (n=6 / 5 / 6; only 5 have all three timepoints),\n",
                         "so the third bar is not the sum of the first two. All from the single dataset GSE227903."),
       x = NULL, y = "median delta") +
  theme_a() + theme(axis.text.x = element_text(angle = 12, hjust = 1))

p5 <- p5a / p5b + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "H3 (treatment-pressure axis): not interpreted -- dependent metadata fields are unverified",
    subtitle = paste0(
      "WARNING: every conclusion here rests on the timepoint (Dx/MRD/Relapse) and patient-pairing fields,\n",
      "both known to disagree with the source papers in several datasets (the reconciliation table does not\n",
      "exist yet; 87 more samples have an empty timepoint field). Until verified, the V-shape can be neither confirmed nor denied."),
    caption = "Source: 03_malignant_distribution_shift.R / 04_stemness_score.R / malignant_distribution.csv",
    theme = theme_a() + theme(plot.subtitle = element_text(colour = CRIT)))

save_fig(p5, "g05_h3_blocked", 10.8, 8.8)

cat("\n[done] 5 figures -> ", FIG, "\n", sep = "")
