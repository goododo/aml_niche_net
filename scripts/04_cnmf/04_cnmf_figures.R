#!/usr/bin/env Rscript
# 04_cnmf_figures.R ----
# Turn the cNMF meta-program results into presentation figures + tables. Reads 02's outputs + the
# finalized labels (mp_labels.tsv) + 03's per-cell MP usage. Writes PDF+PNG to CNMF_FIG_DIR and a
# merged overview table to CNMF_TAB_DIR. Each figure is wrapped so one failure doesn't abort the rest.
#
#   Rscript scripts/04_cnmf/04_cnmf_figures.R

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

TAB <- CNMF_TAB_DIR; FIG <- CNMF_FIG_DIR; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
sig  <- fread(file.path(TAB, "malignant_metaprograms.tsv"))
met  <- fread(file.path(TAB, "malignant_mp_metrics.csv"))
ann  <- fread(file.path(TAB, "malignant_mp_annotation.csv"))
lab  <- fread(file.path(CNMF_SCRIPT_DIR, "mp_labels.tsv"))
CONF_COL <- c(robust_tumor = "#c0392b", normal = "#7f8c8d", likely_contaminant = "#e67e22", drop = "#bdc3c7")
CONF_LVL <- c("robust_tumor", "likely_contaminant", "normal", "drop")
lab[, confidence := factor(confidence, levels = CONF_LVL)]
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  message("[fig] ", name)
}
# CURATED timepoint (third and last copy of the private tp_from_name regex; see
# TP_AXIS_LEVELS in config_hierarchy.R and sample_timepoints() in utils.R). Restricting this
# figure to GSE227903 was never a scientific choice -- it is the only dataset whose sample NAMES
# encode a timepoint, which is all the regex could read. On the curated column the same axis
# covers 28 longitudinal patients across five datasets.
TPMAP <- sample_timepoints()

## ---- merged overview TABLE (CSV always; the master reference) ----
ov <- Reduce(function(a, b) merge(a, b, all.x = TRUE),
             list(lab, met[, .(MP, sampleCoverage, silhouette, numberGenes)],
                  ann[, .(MP = malignant_MP, jaccard_class = class, max_jaccard_vs_normal = max_jaccard)]))
ov <- ov[order(match(MP, paste0("MP", 1:100)))]
ov[, top_genes := sapply(MP, function(m) paste(head(sig[MP == m][order(rank)]$gene, 12), collapse = ", "))]
fwrite(ov, file.path(TAB, "mp_overview.csv"))
cat("=== MP overview ===\n"); print(ov[, .(MP, biological_label, confidence, jaccard_class,
     sampleCoverage = round(sampleCoverage, 2), silhouette = round(silhouette, 2), numberGenes)])

## ---- Fig 1: overview table as a figure ----
tryCatch({
  d <- ov[, .(MP, label = biological_label, conf = as.character(confidence),
              cover = round(sampleCoverage, 2), sil = round(silhouette, 2), nG = numberGenes)]
  d[, y := .N:1]
  p1 <- ggplot(d) +
    geom_tile(aes(x = 0, y = y, fill = conf), width = 0.35, alpha = 0.85) +
    geom_text(aes(x = 0, y = y, label = MP), size = 3, fontface = "bold") +
    geom_text(aes(x = 0.35, y = y, label = label), hjust = 0, size = 3) +
    geom_text(aes(x = 3.1, y = y, label = sprintf("cov=%.2f  sil=%.2f  n=%d", cover, sil, nG)), hjust = 0, size = 2.7, colour = "grey30") +
    scale_fill_manual(values = CONF_COL, name = "confidence", drop = FALSE) +
    scale_x_continuous(limits = c(-0.2, 4.6)) +
    labs(title = "AML malignant meta-programs (per-sample cNMF -> cross-sample)",
         subtitle = "confidence: robust_tumor (>90% Mono_DC bin) | likely_contaminant (diffuse) | normal-like | drop") +
    theme_void(base_size = 11) + theme(plot.title = element_text(face = "bold"),
         legend.position = "bottom", plot.margin = margin(10, 10, 10, 10))
  save_fig(p1, "fig1_mp_overview_table", 11, 4.2)
}, error = function(e) message("[fig1 skip] ", conditionMessage(e)))

## ---- Fig 2: gene-signature heatmap (top genes x MP, weight), MPs ordered/annotated by confidence ----
tryCatch({
  real <- lab[confidence != "drop"]$MP
  topg <- sig[MP %in% real][order(MP, rank)][, head(.SD, 10), by = MP]
  gord <- unique(topg$gene)
  mat  <- sig[gene %in% gord & MP %in% real]
  mat[, w_scaled := weight / max(weight), by = MP]                       # per-MP 0-1 for visibility
  mp_ord <- lab[confidence != "drop"][order(confidence, MP)]$MP
  mat[, MP := factor(MP, levels = mp_ord)]
  mat[, gene := factor(gene, levels = rev(gord))]
  labmap <- setNames(sprintf("%s\n%s", lab$MP, lab$confidence), lab$MP)
  p2 <- ggplot(mat, aes(MP, gene, fill = w_scaled)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "#2c3e50", name = "gene weight\n(scaled)") +
    scale_x_discrete(labels = function(x) labmap[x]) +
    labs(title = "Meta-program gene signatures (top 10 genes per program)",
         x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 6))
  save_fig(p2, "fig2_mp_gene_signatures", 10, 12)
}, error = function(e) message("[fig2 skip] ", conditionMessage(e)))

## ---- Fig 3: malignant-vs-normal Jaccard heatmap (tumor-specificity evidence) ----
tryCatch({
  sim <- fread(file.path(TAB, "malignant_vs_normal_similarity.csv"))
  sim <- merge(sim, lab[, .(malignant_MP = MP, conf = confidence)], by = "malignant_MP")
  malord <- lab[order(confidence, MP)]$MP
  sim[, malignant_MP := factor(malignant_MP, levels = malord)]
  p3 <- ggplot(sim, aes(normal_MP, malignant_MP, fill = jaccard)) + geom_tile() +
    geom_text(aes(label = ifelse(jaccard >= CNMF_MP_SIM_MAX, sprintf("%.2f", jaccard), "")), size = 2.6) +
    scale_fill_gradient2(low = "white", mid = "#f6c667", high = "#c0392b", midpoint = CNMF_MP_SIM_MAX, name = "Jaccard") +
    labs(title = "Malignant vs normal meta-program overlap",
         subtitle = sprintf("Jaccard > %.2f (labeled) = shares genes with a normal program; used to flag normal-like", CNMF_MP_SIM_MAX),
         x = "normal meta-program", y = "malignant meta-program") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p3, "fig3_malignant_vs_normal_similarity", 8, 7)
}, error = function(e) message("[fig3 skip] ", conditionMessage(e)))

## ---- load per-cell MP usage (114 samples) for Fig 4/5 ----
uf <- list.files(file.path(CNMF_RES_DIR, "malignant", "mp_usage"), pattern = "__mp_usage\\.csv$", recursive = TRUE, full.names = TRUE)
u <- if (length(uf)) rbindlist(lapply(uf, fread), fill = TRUE) else data.table()
message(sprintf("[usage] %d cells from %d samples", nrow(u), uniqueN(u$sample)))

## ---- Fig 4: dominant-MP composition per dataset (colored by biological label) ----
tryCatch({
  comp <- u[, .N, by = .(dataset, top_MP)][, frac := N / sum(N), by = dataset]
  comp <- merge(comp, lab[, .(top_MP = MP, biological_label, confidence)], by = "top_MP")
  lord <- lab[order(confidence, MP)]$biological_label                     # all labels, confidence-ordered
  comp[, biological_label := factor(biological_label, levels = lord)]
  pal <- setNames(scales::hue_pal()(length(lord)), lord)                  # named over ALL levels (robust to unused)
  p4 <- ggplot(comp, aes(dataset, frac, fill = biological_label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = pal, name = "dominant meta-program", drop = FALSE) +
    labs(title = "Malignant-cell dominant meta-program composition, per dataset",
         subtitle = sprintf("%d cells, %d samples; each cell assigned its top-scoring program", nrow(u), uniqueN(u$sample)),
         y = "fraction of malignant cells", x = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.text = element_text(size = 8))
  save_fig(p4, "fig4_mp_usage_composition", 11, 6)
}, error = function(e) message("[fig4 skip] ", conditionMessage(e)))

## ---- Fig 5: robust tumor-specific MP usage along the treatment axis (all datasets) ----
tryCatch({
  robust <- lab[confidence == "robust_tumor"]$MP
  u2 <- merge(u, TPMAP[, .(dataset, sample, tp = tp_axis)], by = c("dataset", "sample"), all.x = TRUE)
  u2 <- u2[tp %in% TP_AXIS_LEVELS]
  meas <- intersect(robust, names(u2))
  long <- melt(u2, id.vars = c("cell", "sample", "tp"), measure.vars = meas, variable.name = "MP", value.name = "score")
  long <- merge(long, lab[, .(MP, biological_label)], by = "MP")
  long[, tp := factor(tp, levels = TP_AXIS_LEVELS)]
  agg <- long[, .(mean = mean(score, na.rm = TRUE),
                  se = sd(score, na.rm = TRUE) / sqrt(.N)), by = .(tp, biological_label)]
  p5 <- ggplot(agg, aes(tp, mean, colour = biological_label, group = biological_label)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12) +
    labs(title = "Robust tumor-specific program usage along the treatment axis (GSE227903)",
         subtitle = "mean per-cell module score +/- SE; the 3 monocytic-AML programs (Mono_DC-dominant)",
         y = "mean MP module score", x = NULL, colour = "meta-program") +
    theme(legend.position = "bottom")
  save_fig(p5, "fig5_mp_usage_along_treatment", 8, 6)
}, error = function(e) message("[fig5 skip] ", conditionMessage(e)))

cat("\n[done] figures -> ", FIG, "\n[done] overview table -> ", file.path(TAB, "mp_overview.csv"), "\n", sep = "")
