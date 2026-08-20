# d_deck_figures.R ----
# Deck-optimized figures for the lab meeting: LARGE fonts, minimal text (the slide
# carries the title + conclusion; the figure shows only the visualization), and one
# NEW cell-cell communication content figure (who-talks-to-whom + top ligand-receptor).
# INPUT  : results/tables/{01_preprocess,02_malignancy,03_hierarchy,05_ccc,07_fgw,08_scoring}/*.csv
#          (CCC summary tables produced by agg_ccc_lr_summary.py)
# OUTPUT : results/figures/12_deck/d0{1..6}_*.{png,pdf}
# Usage  : conda run -p /FAST/gr10634/gaozy/general_env Rscript scripts/10_figures/d_deck_figures.R
#
# [DECISION] Same validated palette as g_assessment_figures.R (BLUE #2a78d6 / ORANGE
#            #eb6834 / status #d03b3b) so deck figures and the g-series read as one system.
# [DECISION] No long subtitles/captions on the figure -- those become slide text. Panels
#            keep only a short A/B label, axis titles, legend, and direct data labels.
# [DECISION] Sequential magnitude = single-hue blue light->dark; AML-vs-healthy difference
#            = diverging blue<->red with a grey midpoint (dataviz rules).

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
TBL  <- file.path(ROOT, "results", "tables")

# The treatment-axis vocabulary must come from the config, not be re-listed here. These figure
# scripts were standalone, so TP_AXIS_LEVELS would otherwise be undefined at runtime -- and the
# literal they used to carry, c("Dx","MRD","Relapse"), matches only "Relapse" against the rebuilt
# tables and would empty the panel without erroring.
suppressPackageStartupMessages({
  source(file.path(ROOT, "scripts", "config", "config_paths.R"))
  source(file.path(ROOT, "scripts", "config", "config_qc.R"))
  source(file.path(ROOT, "scripts", "config", "config_hierarchy.R"))
})
stopifnot(exists("TP_AXIS_LEVELS"))
FIG  <- file.path(ROOT, "results", "figures", "12_deck")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; CRIT <- "#d03b3b"; WARN <- "#fab219"
GOOD <- "#0ca30c"; GREY <- "#8a8a85"; LGREY <- "#c6c5c0"; MIDGREY <- "#f0efec"
INK  <- "#0b0b0b"; INK2 <- "#3f3e3c"; FAM <- "sans"
NODES <- c("HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma")

BASE <- 18
theme_deck <- function(base = BASE) {
  theme_minimal(base_size = base, base_family = FAM) +
    theme(
      plot.title    = element_text(face = "bold", size = base + 1, colour = INK, hjust = 0,
                                   margin = margin(b = 6)),
      axis.title    = element_text(size = base - 2, colour = INK2),
      axis.text     = element_text(size = base - 3, colour = INK2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e1", linewidth = 0.3),
      legend.position = "top", legend.justification = "left",
      legend.title  = element_text(size = base - 4, colour = INK2),
      legend.text   = element_text(size = base - 4, colour = INK2),
      legend.key.size = unit(1.0, "lines"),
      plot.title.position = "plot",
      plot.margin = margin(10, 16, 10, 10),
      strip.text = element_text(face = "bold", colour = INK, size = base - 2)
    )
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 200,
         device = ragg::agg_png, bg = "white")
  cat(sprintf("[deck-fig] %s\n", name))
}
pf <- function(p) ifelse(p < 1e-3, "p<0.001", paste0("p=", formatC(p, format = "g", digits = 2)))

## ══ d01 — cohort counts (single panel) ═════════════════════════════════════════
qc   <- rbindlist(lapply(Sys.glob(file.path(TBL,"01_preprocess","03_qc_report__*.csv")), fread), fill=TRUE)
cons <- fread(file.path(TBL,"02_malignancy","ALL_consensus_summary.csv"))
fgwix<- fread(file.path(TBL,"07_fgw","fgw_input_index.csv"))
asw  <- fread(file.path(TBL,"08_scoring","alpha_sweep.csv"))
fn <- data.table(
  stage = c("QC-passing","Malignancy label","CCC graph built","In FGW model","Platform-controlled test"),
  n     = c(nrow(qc[status=="PASS"]), nrow(cons), nrow(fgwix), asw$n_global[1], asw$n_strat[1]))
fn[, stage := factor(stage, levels = rev(stage))]
fn[, hl := c(rep(FALSE,4), TRUE)]
d1 <- ggplot(fn, aes(stage, n, fill = hl)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = paste0(n, " samples")), hjust = -0.08, size = 6, colour = INK, family = FAM) +
  scale_fill_manual(values = c("FALSE"=BLUE,"TRUE"=CRIT), guide="none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.34))) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "number of samples") +
  theme_deck() + theme(panel.grid.major.y = element_blank())
save_fig(d1, "d01_cohort", 10.5, 5.2)

## ══ d02 — label quality (A: FPR by bin, B: evidence tier) ══════════════════════
fpr <- fread(file.path(TBL,"02_malignancy","malignancy_fpr_by_bin.csv"))
fpr[hierarchy_bin=="" | is.na(hierarchy_bin), hierarchy_bin := "Unassigned"]
fpr[, ingraph := !is.na(in_ccc_graph) & in_ccc_graph == TRUE]
fpr[, hierarchy_bin := factor(hierarchy_bin, levels = hierarchy_bin[order(FPR)])]
fpr[, hl := fifelse(as.character(hierarchy_bin)=="HSC_MPP","crit", fifelse(ingraph,"in","out"))]
d2a <- ggplot(fpr, aes(hierarchy_bin, FPR, fill = hl)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", 100*FPR)), hjust = -0.12, size = 6, colour = INK, family = FAM) +
  scale_fill_manual(values = c(crit=CRIT,"in"=BLUE,out=LGREY), guide="none") +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0,1),
                     breaks = seq(0,1,.25), expand = expansion(mult=c(0,0.04))) +
  coord_flip(clip="off") +
  labs(title = "A   % of healthy cells called malignant, by cell type", x=NULL,
       y = "% of healthy cells wrongly called malignant") +
  theme_deck() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(8,60,8,8))
tier <- cons[, .N, by = evidence_tier]
tier[, evidence_tier := factor(evidence_tier, levels=c("A_concordant","B_multi_partial","C_single"),
                               labels=c("2-3 arms agree","2 arms (partial)","1 arm only (inferCNV)"))]
setorder(tier, evidence_tier)
d2b <- ggplot(tier, aes(x = "", y = N, fill = evidence_tier)) +
  geom_col(width = 0.5, colour = "white", linewidth = 1.4) +
  geom_text(data = tier[N>10], aes(label = paste0(evidence_tier, ": ", N)),
            position = position_stack(vjust = 0.5), size = 5.4, colour = "white", fontface="bold", family=FAM) +
  scale_fill_manual(values = c("2-3 arms agree"=GOOD,"2 arms (partial)"=WARN,"1 arm only (inferCNV)"=GREY), guide="none") +
  scale_y_continuous(expand = expansion(mult=c(0,0.02))) +
  coord_flip() +
  labs(title = sprintf("B   How each of %d samples was labelled", nrow(cons)), x=NULL, y="samples") +
  theme_deck() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_blank())
save_fig(d2a / d2b + plot_layout(heights = c(2.4, 1)), "d02_label_quality", 10.5, 8.6)

## ══ d03 — H1 (A: 7x7 dC heatmap, B: block q) ═══════════════════════════════════
ee  <- fread(file.path(TBL,"08_scoring","emergent_edges.csv"))
blk <- fread(file.path(TBL,"08_scoring","block_permutation.csv"))
ee[, sender_bin := factor(sender_bin, levels = NODES)]
ee[, receiver_bin := factor(receiver_bin, levels = rev(NODES))]
lim <- max(abs(ee$dC_real))
d3a <- ggplot(ee, aes(sender_bin, receiver_bin, fill = dC_real)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  scale_fill_gradient2(low = CRIT, mid = MIDGREY, high = BLUE, midpoint = 0, limits = c(-lim, lim),
                       breaks = c(-0.1, 0, 0.1), name = "AML − healthy") +
  coord_fixed() +
  labs(title = "A   AML minus healthy, per communication link", x="sender", y="receiver") +
  theme_deck() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid.major = element_blank(),
        legend.key.width = unit(2.4,"lines"))
setorder(blk, -n_edges)
blk[, block := factor(block, levels = blk$block)]
blk[, lab := gsub("_"," ", block)]
blk[, lab := factor(lab, levels = lab)]
d3b <- ggplot(blk, aes(lab, q_bh)) +
  geom_col(width = 0.62, fill = LGREY) +
  geom_hline(yintercept = 0.05, linetype = "22", colour = CRIT, linewidth = 0.8) +
  annotate("text", x = 0.6, y = 0.11, label = "significance (0.05)", hjust = 0, size = 5, colour = CRIT, family = FAM) +
  geom_text(aes(label = sprintf("q=%.2f", q_bh)), hjust = -0.15, size = 5.6, colour = INK, family = FAM) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.25), expand = expansion(mult=c(0,0.04))) +
  coord_flip(clip="off") +
  labs(title = "B   Communication families, by false-discovery rate", x=NULL, y="q (false-discovery rate)") +
  theme_deck() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(8,70,8,8))
save_fig(d3a / d3b + plot_layout(heights = c(2.2, 1.15)), "d03_h1_null", 10.2, 9.6)

## ══ d04 — H2 (A: alpha sweep, B: feature decomposition) ════════════════════════
asw_l <- melt(asw, id.vars="alpha", measure.vars=c("p_global","p_strat"), variable.name="model", value.name="p")
asw_l[, model := factor(model, levels=c("p_global","p_strat"), labels=c("all datasets","within dataset"))]
d4a <- ggplot(asw_l, aes(alpha, p, colour = model)) +
  annotate("rect", xmin=0.9, xmax=1.02, ymin=5e-5, ymax=1, fill=CRIT, alpha=0.08) +
  geom_hline(yintercept = 0.05, linetype="22", colour=CRIT, linewidth=0.8) +
  geom_line(linewidth = 1.2) + geom_point(size = 3.4) +
  annotate("text", x=0.96, y=0.42, label="only network\nshape", size=5, colour=CRIT, family=FAM, lineheight=0.95) +
  annotate("text", x=0.0, y=0.075, label="significance (0.05)", hjust=0, size=5, colour=CRIT, family=FAM) +
  scale_colour_manual(values = c(BLUE, ORANGE), name = NULL) +
  scale_y_log10(limits=c(5e-5,1), breaks=c(1e-4,1e-3,1e-2,0.05,0.1,0.3,1), labels=label_number(drop0trailing=TRUE)) +
  scale_x_continuous(breaks = asw$alpha) +
  labs(title = "A   Significance across the shape-vs-features weight",
       x = "weight on network shape  (0 = features only → 1 = shape only)", y = "p-value (log scale)") +
  theme_deck()
fd <- fread(file.path(TBL,"08_scoring","feature_decomposition.csv"))[alpha==0.5]
fd[, feature_set := factor(feature_set, levels=c("only_stemness","only_ncells","no_frac_mal","all3","only_frac_mal"),
     labels=c("stemness only","cell count only","without blast fraction","all features","blast fraction only"))]
setorder(fd, feature_set)
fd[, kind := fifelse(as.character(feature_set)=="blast fraction only","circular",
              fifelse(p_global<0.05,"significant","not significant"))]
fd[, kind := factor(kind, levels=c("circular","significant","not significant"))]
d4b <- ggplot(fd, aes(feature_set, beta_global, fill = kind)) +
  geom_col(width = 0.62) +
  geom_text(aes(y = pmax(beta_global, 0.004), label = sprintf("%.3f", beta_global)),
            hjust = -0.15, size = 5.4, colour = INK, family = FAM) +
  scale_fill_manual(values = setNames(c(CRIT,BLUE,LGREY), levels(fd$kind)), name = NULL) +
  scale_y_continuous(limits = c(-0.05, 0.24), expand = expansion(mult=c(0.02,0.06))) +
  coord_flip(clip="off") +
  labs(title = "B   Effect size, by feature set", x=NULL, y="effect on AML-vs-healthy") +
  theme_deck() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(8,55,8,8))
save_fig(d4a / d4b + plot_layout(heights = c(1.15, 1)), "d04_h2", 10.6, 9.4)

## ══ d05 — H3 (A: stemness by timepoint, B: within-pair change) ═════════════════
stem <- fread(file.path(TBL,"03_hierarchy","stemness_by_timepoint.csv"))
sl <- melt(stem, id.vars=c("tp","n_malignant"),
           measure.vars=c("LSC17","vanGalen_HSC_Prog","vanGalen_HSC_like","HSPC_core"),
           variable.name="signature", value.name="score")
# The axis comes from config (TP_AXIS_LEVELS). c("Dx","MRD","Relapse") was a third
# vocabulary and "MRD" was retired on 2026-08-04; against the rebuilt tables it matches
# only "Relapse", which would empty this figure without erroring.
sl[, tp := factor(tp, levels = TP_AXIS_LEVELS)]
sl[, vshape := signature %in% c("LSC17","vanGalen_HSC_Prog")]
d5a <- ggplot(sl, aes(tp, score, colour = signature, group = signature, linetype = vshape)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3.2) +
  geom_text(data = sl[tp=="Relapse"], aes(label = signature), hjust=0, nudge_x=0.06, size=5, family=FAM, show.legend=FALSE) +
  scale_colour_manual(values = c(LSC17=BLUE, vanGalen_HSC_Prog=ORANGE, vanGalen_HSC_like=GREY, HSPC_core="#4a3aa7"), guide="none") +
  scale_linetype_manual(values = c("TRUE"="solid","FALSE"="21"), guide="none") +
  scale_x_discrete(expand = expansion(add = c(0.35, 1.7))) +
  labs(title = "A   Stemness score by treatment timepoint", x=NULL, y="mean stemness score") +
  theme_deck() + theme(panel.grid.major.x = element_blank())
sh <- fread(file.path(TBL,"03_hierarchy","distribution_shift_tests.csv"))
sh[, comparison := factor(comparison, levels=c("Dx -> MRD","MRD -> Relapse","Dx -> Relapse"))]
sh[, score := factor(score, levels=c("primitive_frac","stem_frac"),
      labels=c("primitive fraction","stem fraction"))]
d5b <- ggplot(sh, aes(comparison, median_delta)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.5) +
  geom_col(width = 0.62, fill = LGREY) +
  geom_text(aes(label = pf(p_value), vjust = ifelse(median_delta>0,-0.4,1.3),
                fontface = ifelse(p_value<0.05,"bold","plain")), size = 4.8, colour = INK2, family = FAM) +
  facet_wrap(~score, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult=c(0.28,0.28))) +
  labs(title = "B   Median change within paired patients", x=NULL, y="median change") +
  theme_deck() + theme(axis.text.x = element_text(angle = 10, hjust = 1))
save_fig(d5a / d5b + plot_layout(heights = c(1.05, 1)), "d05_h3", 11.0, 9.0)

## ══ d06 — NEW: CCC content (A: who-talks-to-whom, B: top ligand-receptor) ══════
bs <- fread(file.path(TBL,"05_ccc","ccc_bin_strength_cohort.csv"))
bs[, sender_bin := factor(sender_bin, levels = NODES)]
bs[, receiver_bin := factor(receiver_bin, levels = rev(NODES))]
d6a <- ggplot(bs, aes(sender_bin, receiver_bin, fill = mean_weight)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  scale_fill_gradient(low = "#f2f6fc", high = BLUE, name = "communication\nstrength") +
  coord_fixed() +
  labs(title = "A   Communication strength: sender → receiver (AML)", x = "sender", y = "receiver") +
  theme_deck() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), panel.grid.major = element_blank(),
        legend.key.width = unit(2.4,"lines"))
lr <- fread(file.path(TBL,"05_ccc","ccc_lr_overall.csv"))[order(-n_samples)][1:10]
lr[, pair := paste0(ligand, " → ", gsub("_","/",receptor))]
lr[, pair := factor(pair, levels = rev(pair))]
lr[, fam := fifelse(grepl("^MIF",ligand),"MIF",
             fifelse(grepl("^LGALS",ligand),"Galectin-9", "other"))]
lr[, fam := factor(fam, levels = c("MIF","Galectin-9","other"))]
d6b <- ggplot(lr, aes(pair, frac_samples, fill = fam)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.0f%%", 100*frac_samples)), hjust = -0.14, size = 5.2, colour = INK, family = FAM) +
  scale_fill_manual(values = c("MIF"=BLUE,"Galectin-9"=ORANGE,"other"=LGREY), name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy=1), limits = c(0,1), expand = expansion(mult=c(0,0.06))) +
  coord_flip(clip="off") +
  labs(title = "B   Top ligand-receptor pairs (% of AML samples)", x=NULL, y="% of samples with this interaction") +
  theme_deck() + theme(panel.grid.major.y = element_blank(), plot.margin = margin(8,55,8,8))
save_fig(d6a / d6b + plot_layout(heights = c(1.25, 1)), "d06_ccc_content", 10.6, 10.2)

cat("\n[done] deck figures -> ", FIG, "\n", sep="")
