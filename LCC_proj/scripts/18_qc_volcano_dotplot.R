#!/usr/bin/env Rscript
# 18_qc_volcano_dotplot.R ----
# Three views the brief asked for that the panel set was missing entirely.
#
#   Q1  sample quality / QC for the 18 samples that entered the comparison
#   V1  volcano plot over the whole gene panel (145 genes with a testable paired result)
#   V2  dot plot: gene x (cell type x arm), dot size = % of cells expressing, colour = mean expression
#
# WHY THESE ARE NEW: P2/P3/P7 all show the same quantity -- "in how many of the 9 pairs is the
# mutant higher". That is the statistic the test operates on, but it is a statistical view, not a
# biological one, and repeating it across every panel left the reader with no way to see the
# expression itself. A volcano and a dot plot are the two forms a reader expects first for a
# two-group gene comparison, and neither existed. Q1 simply did not exist at all.
#
# INPUT  : results/tables/00_ingest/00_MASTER_qc_summary.csv   (main-line QC, all samples)
#          LCC_TAB_DIR/11_matched_design.csv                   the 9 pairs
#          LCC_TAB_DIR/11_gene_results.csv                     paired effect + p, per stratum
#          LCC_TAB_DIR/12_robustness_genes.csv
# OUTPUT : LCC_FIG_DIR/{Q1_sample_qc, V1_volcano, V2_dotplot}.{png,pdf}
#          LCC_TAB_DIR/Q1_sample_qc.csv
suppressPackageStartupMessages({
  library(data.table); library(here); library(Matrix); library(ggplot2); library(ggrepel); library(patchwork)
  library(scales)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
set.seed(SEED)

des  <- fread(file.path(LCC_TAB_DIR, "11_matched_design.csv"))
arms <- rbind(des[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              des[, .(pair_id, dataset, sample = wt_sample,  arm = "TP53-WT")])
arms[, arm := factor(arm, levels = c("TP53-WT", "TP53-mut"))]

## ===== Q1. sample quality =====================================================================
# The point of this panel is a NEGATIVE: the two arms must be comparable in library quality, or any
# expression difference could be a quality difference. Both arms are drawn so the reader can check
# rather than take it on trust.
#
# METRICS ARE COMPUTED ON THE CELLS THAT ENTERED THE COMPARISON, not at ingest. A first version read
# 00_ingest/00_MASTER_qc_summary.csv, which is the PRE-FILTER summary: it reported 21,731 cells for
# 809653 where 12,340 actually entered, and 84,603 cells across the 18 samples where 49,519 entered.
# A QC figure that describes a different cell set from every other figure in the project is worse
# than no QC figure. The ingest counts are kept, but as the "before" half of a retention panel.
message("[Q1] sample QC")
ing <- fread(here::here("results", "tables", "00_ingest", "00_MASTER_qc_summary.csv"))
setnames(ing, c("Dataset", "Sample", "n_cells"), c("dataset", "sample", "n_cells_ingest"))

MTPAT   <- "^MT-"
RIBOPAT <- "^RP[SL]"
HBGENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
per_sample_qc <- function(ds, sm) {
  f <- file.path(QC_RDS_DIR, ds, paste0(sm, ".rds"))
  if (!file.exists(f)) { message("   [warn] missing ", sm); return(NULL) }
  cn <- get_counts(readRDS(f))
  tot <- Matrix::colSums(cn)
  keep <- tot > 0; cn <- cn[, keep, drop = FALSE]; tot <- tot[keep]
  g <- rownames(cn)
  frac <- function(idx) if (!length(idx)) rep(0, length(tot)) else
    100 * Matrix::colSums(cn[idx, , drop = FALSE]) / tot
  nf <- Matrix::colSums(cn > 0)
  data.table(dataset = ds, sample = sm, n_cells = length(tot),
             median_nFeature = median(nf), min_nFeature = min(nf), max_nFeature = max(nf),
             median_nCount = median(tot),
             median_pct_mt   = median(frac(grep(MTPAT, g))),
             median_pct_ribo = median(frac(grep(RIBOPAT, g))),
             median_pct_hb   = median(frac(which(g %in% HBGENES))))
}
qm <- rbindlist(lapply(seq_len(nrow(arms)), function(i) {
  message("   ", i, "/", nrow(arms), "  ", arms$sample[i])
  per_sample_qc(arms$dataset[i], arms$sample[i])
}))
q <- merge(merge(arms, qm, by = c("dataset", "sample")),
           ing[, .(dataset, sample, n_cells_ingest)], by = c("dataset", "sample"), all.x = TRUE)
q[, pct_kept := 100 * n_cells / n_cells_ingest]
fwrite_safe(q, file.path(LCC_TAB_DIR, "Q1_sample_qc.csv"))
message("    ", sum(q$n_cells), " cells retained of ", sum(q$n_cells_ingest), " ingested (",
        round(100 * sum(q$n_cells) / sum(q$n_cells_ingest)), "%)")

# every metric gets a paired test, so the "arms are comparable" claim is measured, not eyeballed
for (v in c("n_cells", "median_nFeature", "median_nCount",
            "median_pct_mt", "median_pct_ribo", "median_pct_hb")) {
  w <- dcast(q, pair_id ~ arm, value.var = v); d <- w[["TP53-mut"]] - w[["TP53-WT"]]
  message(sprintf("    %-16s mut higher in %d/%d pairs, signed-rank p = %.3f",
                  v, sum(d > 0, na.rm = TRUE), sum(!is.na(d)),
                  suppressWarnings(wilcox.test(d)$p.value)))
}

METRICS <- c(n_cells          = "cells per sample",
             median_nFeature  = "genes per cell (median)",
             median_nCount    = "UMIs per cell (median)",
             median_pct_mt    = "mitochondrial %",
             median_pct_ribo  = "ribosomal %",
             median_pct_hb    = "haemoglobin %")
ql <- melt(q[, c("pair_id", "arm", "sample", names(METRICS)), with = FALSE],
           id.vars = c("pair_id", "arm", "sample"), variable.name = "metric", value.name = "v")
ql[, metric := factor(METRICS[as.character(metric)], levels = METRICS)]
# counts span 500-21,000 across samples; on a linear axis the two 15k-cell samples flatten the
# other sixteen into a single band. Percentages stay linear -- they are already on a bounded scale.
ql[, is_count := grepl("cells|genes|UMIs", metric)]

# paired: one line per pair, so a systematic arm difference would show as parallel lines
qa <- ggplot(ql, aes(arm, v, group = pair_id)) +
  geom_line(colour = GRID, linewidth = 0.9) +
  geom_point(aes(colour = arm), size = 3) +
  scale_colour_manual(values = c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]]),
                      name = NULL) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2) +
  # sqrt, not log: the haemoglobin-% facet contains exact zeros, which log cannot take. sqrt still
  # pulls in the two 15,000-cell samples that otherwise flatten the other sixteen into one band.
  scale_y_continuous(trans = "sqrt", n.breaks = 5) +
  labs(title = "Q1a  Library quality of the 18 samples that entered, pair by pair",
       subtitle = "each line is one matched pair; y on a square-root scale", x = NULL, y = NULL) +
  theme_lcc() + theme(panel.grid.major.x = element_blank())

# and the per-sample detail, because "median genes per cell" hides a sample with a bad tail
qb <- ggplot(q, aes(reorder(sample, median_nFeature), median_nFeature, fill = arm)) +
  geom_col(width = 0.68) +
  geom_errorbar(aes(ymin = min_nFeature, ymax = max_nFeature), width = 0, colour = INK[["muted"]],
                linewidth = 0.5) +
  scale_fill_manual(values = c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]]), name = NULL) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Q1b  Genes per cell, per sample",
       subtitle = "bar = median, whisker = min-max across cells", x = NULL, y = "genes per cell") +
  theme_lcc()

# Q1c -- the filtering funnel. Without it a reader cannot tell whether an arm lost more cells to QC
# than the other, which is itself a quality difference.
# ONE PANEL, NOT FACETTED BY ARM. A facetted first version put all 18 sample names on both panels'
# axes -- free_y does not drop a factor's unused levels -- so half the rows were empty and the
# "% kept" labels sat next to the wrong sample. An overlay bar (grey = ingested behind,
# arm-coloured = retained in front) needs no facet and reads in one pass.
ord <- q[order(arm, n_cells_ingest), sample]
q[, samp := factor(sample, levels = ord)]
qc <- ggplot(q, aes(samp)) +
  geom_col(aes(y = n_cells_ingest), fill = GRID, width = 0.8) +
  geom_col(aes(y = n_cells, fill = arm), width = 0.46) +
  geom_text(aes(y = n_cells_ingest, label = sprintf("%.0f%% kept", pct_kept)),
            hjust = -0.13, size = 3.9, colour = INK[["muted"]]) +
  scale_fill_manual(values = c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]]), name = NULL) +
  coord_flip() +
  scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.20))) +
  labs(title = "Q1c  Cells before and after QC filtering",
       subtitle = "grey bar = cells at ingest; coloured bar = cells that entered the comparison",
       x = NULL, y = "cells") +
  theme_lcc()

save_fig(qa / qb / qc + patchwork::plot_layout(heights = c(1, 1.05, 1.05)), "Q1_sample_qc", 13, 17)

## ===== V1. volcano ============================================================================
message("[V1] volcano")
gr <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))[stratum == "all" & !is.na(p_sample_higher)]
rb <- fread(file.path(LCC_TAB_DIR, "12_robustness_genes.csv"))[design == "A_sample_strict" & kind == "gene"]
gr[rb, robust := i.frac_p05, on = c(gene = "feature")]
gr[, pairs_up := round(frac_pairs_mut_higher * n_pairs)]

# A REAL FOLD CHANGE FOR THE X AXIS. `median_delta_log2cpm` is the median paired DIFFERENCE of mean
# log-normalised expression; for the many lowly-detected genes in this panel that difference is
# 0.001-0.05, so plotting it puts almost every gene on top of the y axis and the volcano says nothing.
# Recompute the quantity a volcano is supposed to carry: the per-pair log2 RATIO of mean expression,
# then the median across the 9 pairs. Same paired design, same 9 pairs, interpretable magnitude.
det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"),
             select = c("gene", "stratum", "sample", "dataset", "mean_lognorm", "gene_present"))
det <- det[stratum == "all" & gene_present == TRUE]
det <- merge(det, arms, by = c("dataset", "sample"))
det <- det[, .(v = max(mean_lognorm, na.rm = TRUE)), by = .(gene, pair_id, arm)]
fcw <- dcast(det, gene + pair_id ~ arm, value.var = "v")
setnames(fcw, c("TP53-mut", "TP53-WT"), c("v_mut", "v_wt"))
EPS <- 1e-3                                   # floor: both arms can be at the detection floor
fcw[, lfc := log2((v_mut + EPS) / (v_wt + EPS))]
fc <- fcw[!is.na(lfc), .(log2FC = median(lfc), n_pairs_fc = .N), by = gene]
gr[fc, `:=`(log2FC = i.log2FC, n_pairs_fc = i.n_pairs_fc), on = "gene"]
gr <- gr[!is.na(log2FC)]
# The test is ONE-SIDED ("higher in mutant"), so a gene that is lower in the mutant arm gets a p near
# 1, not a small p. For a volcano we need a two-sided reading, so the down side is scored by the
# mirror probability. Stated here because a reader who assumes a two-sided test would misread the plot.
gr[, p_two := 2 * pmin(p_sample_higher, 1 - p_sample_higher)]
gr[, p_two := pmax(p_two, 1 / 512)]                      # 9 pairs cannot resolve below 1/512
gr[, sig := fifelse(p_sample_higher < 0.05, "higher in TP53-mut (p < 0.05)",
             fifelse(p_sample_higher > 0.95, "lower in TP53-mut (p < 0.05)", "not significant"))]

V1COL <- c(`higher in TP53-mut (p < 0.05)` = PAL[["orange"]],
           `lower in TP53-mut (p < 0.05)`  = PAL[["blue"]],
           `not significant`               = "#c8c8c4")
lab <- gr[sig != "not significant"]
fV1 <- ggplot(gr, aes(log2FC, -log10(p_two))) +
  # THE LINE HAS TO MATCH THE COLOURING. Points are coloured by the project's one-sided p (the test
  # every other figure and table reports), but y is the two-sided mirror so the plot can show both
  # directions. A line drawn at two-sided 0.05 therefore sat ABOVE several coloured points -- the
  # figure contradicted its own legend. One-sided 0.05 is two-sided 0.10, so the line goes there and
  # says so.
  geom_hline(yintercept = -log10(0.10), linetype = "dashed", colour = INK[["muted"]]) +
  geom_vline(xintercept = 0, colour = GRID) +
  geom_point(aes(colour = sig, size = robust), alpha = 0.85) +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 4.3, colour = INK[["primary"]],
                           max.overlaps = 40, min.segment.length = 0.2, seed = SEED) +
  scale_colour_manual(values = V1COL, name = NULL) +
  scale_size_continuous(range = c(1.6, 6), name = "robustness", labels = percent_format(accuracy = 1)) +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  labs(title = sprintf("V1  All %d panel genes, TP53-mut vs TP53-WT", nrow(gr)),
       subtitle = "9 pairs; dashed line = one-sided p = 0.05; size = robustness",
       x = "median paired log2 fold change  (TP53-mut / TP53-WT)",
       y = "-log10(p), two-sided") +
  theme_lcc() +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 4)),
         size = guide_legend(order = 2))
save_fig(fV1, "V1_volcano", 13, 9)

## ===== V2. dot plot ===========================================================================
# The classic form: rows = genes, columns = cell compartment x arm, dot size = fraction of cells
# expressing, colour = mean expression. This is the view that shows WHERE in the marrow a gene is
# expressed, which none of the paired-count panels can show.
message("[V2] dot plot")
STRATA <- c(all = "all cells", malignant = "leukaemic blasts",
            HSC_MPP = "HSC / MPP", Mono_DC = "monocyte / DC")
SHOW <- c("C1QA", "C1QB", "C1QC", "MRC1", "LYVE1", "MAF", "VSIG4",          # macrophage
          "SPP1", "OSM", "IL11", "PDGFB", "PDGFRB", "TGFB1",                 # profibrotic
          "PLOD2", "LOX", "MMP14", "TIMP1",                                  # cross-linking
          "COL1A1", "COL1A2", "COL3A1", "DCN",                               # collagen (negative)
          "SLC2A1", "PGK1", "LDHA",                                          # hypoxia
          "GP9", "ITGA2B", "PF4",                                            # megakaryocyte
          "VIM", "SNAI2", "TWIST1")                                          # EMT (negative)
dp <- fread(file.path(LCC_TAB_DIR, "11_gene_results.csv"))[
        gene %in% SHOW & stratum %in% names(STRATA)]
d2 <- rbind(dp[, .(gene, stratum, arm = "TP53-mut", pct = pct_mut, expr = mean_lognorm_mut)],
            dp[, .(gene, stratum, arm = "TP53-WT",  pct = pct_WT,  expr = mean_lognorm_WT)])
d2 <- d2[!is.na(pct) & !is.na(expr)]
d2[, stratum := factor(STRATA[stratum], levels = STRATA)]
d2[, arm := factor(arm, levels = c("TP53-WT", "TP53-mut"))]
d2[, gene := factor(gene, levels = rev(SHOW))]
# colour is scaled WITHIN gene: absolute expression spans ~3 orders of magnitude across this panel,
# so a shared scale would render every low-expressed gene identically blank and hide the contrast the
# figure exists to show.
d2[, expr_rel := if (max(expr) > 0) expr / max(expr) else 0, by = gene]

fV2 <- ggplot(d2, aes(arm, gene)) +
  geom_point(aes(size = pct, colour = expr_rel)) +
  facet_wrap(~ stratum, nrow = 1) +
  scale_colour_gradient(low = "#e8f2ea", high = PAL[["aqua"]],
                        name = "mean expression\n(relative, within gene)", n.breaks = 3) +
  # sqrt: most panel genes are detected in 1-8% of cells, so a linear size scale renders
  # almost the whole plot as invisible specks and only TGFB1/TIMP1/PGK1/LDHA/VIM read at all.
  scale_size_continuous(range = c(0.5, 8), trans = "sqrt", name = "% of cells expressing",
                        breaks = c(1, 5, 20, 50)) +
  labs(title = "V2  Where these genes are expressed, and in which arm",
       subtitle = "dot size = % of cells detected; colour = mean expression scaled within each gene",
       x = NULL, y = NULL) +
  theme_lcc() +
  theme(panel.grid.major.x = element_blank(),
        legend.box = "horizontal", legend.key.width = unit(2.4, "lines"),
        legend.spacing.x = unit(1.6, "lines"))
save_fig(fV2, "V2_dotplot", 13, 13)

message("[done]")
