#!/usr/bin/env Rscript
# 19_percell_views.R ----
# The two per-cell views the brief asked for and the panel set did not have:
#
#   V3  violin plots -- the actual expression distribution behind every "7/9 pairs" claim
#   U1  UMAP -- where the cells are, which arm they came from, and where key genes sit
#
# WHY: every earlier gene panel reports a paired COUNT ("in 7 of 9 pairs the mutant is higher").
# That is the statistic the test operates on, but a reader cannot see from it whether the underlying
# expression is a broad shift or three outlier cells. These two figures show the cells themselves.
#
# BOTH read raw counts, so they share one pass over the 18 RDS objects.
#
# Deliberately does NOT load Seurat: the user library carries Seurat 5.5.1, which needs
# promises >= 1.5.0 while the env has 1.3.3. Normalisation, HVG, PCA, Harmony and UMAP are called
# directly, which is all Seurat would have done anyway.
#
# INPUT  : LCC_TAB_DIR/11_matched_design.csv
#          QC_RDS_DIR/<ds>/<sample>.rds
#          LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz   (cell-type bin, malignant call)
# OUTPUT : LCC_FIG_DIR/{V3_violin, U1_umap}.{png,pdf}
#          LCC_TAB_DIR/U1_umap_coords.csv.gz
suppressPackageStartupMessages({
  library(data.table); library(here); library(Matrix); library(ggplot2); library(patchwork)
  library(irlba); library(uwot); library(harmony); library(ggrastr); library(scales)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
set.seed(SEED)

# genes for the violin: two from each story line, plus the collagen negative control, chosen to span
# the detection range so the reader sees both a well-detected gene and a floor-level one.
VGENES <- c("C1QB", "C1QC", "MAF", "LYVE1",        # macrophage programme
            "SPP1", "OSM", "IL11", "PDGFB",        # profibrotic effectors
            "PLOD2", "LOX",                        # collagen cross-linking
            "SLC2A1", "PGK1",                      # hypoxia
            "GP9",                                 # megakaryocyte
            "COL1A1", "COL1A2", "VIM")             # collagen / EMT -- the negatives
UGENES <- c("C1QB", "SPP1", "SLC2A1", "COL1A1")    # genes painted onto the UMAP

des  <- fread(file.path(LCC_TAB_DIR, "11_matched_design.csv"))
arms <- rbind(des[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              des[, .(pair_id, dataset, sample = wt_sample,  arm = "TP53-WT")])

## -- one pass over the 18 objects ----
message("[1] reading ", nrow(arms), " samples")
NORM_TARGET <- 1e4
grab <- function(i) {
  ds <- arms$dataset[i]; sm <- arms$sample[i]
  f <- file.path(QC_RDS_DIR, ds, paste0(sm, ".rds"))
  if (!file.exists(f)) { message("   [warn] missing ", sm); return(NULL) }
  cnt <- get_counts(readRDS(f))
  cs  <- Matrix::colSums(cnt)
  keep <- cs > 0
  cnt <- cnt[, keep, drop = FALSE]; cs <- cs[keep]
  pc <- fread(file.path(LCC_PERCELL_DIR, ds, paste0(sm, "__lcc_percell.csv.gz")),
              select = c("cell", "hierarchy_bin", "malignant", "myeloid_class"), showProgress = FALSE)
  pc <- pc[match(colnames(cnt), cell)]
  meta <- data.table(cell = colnames(cnt), dataset = ds, sample = sm,
                     pair_id = arms$pair_id[i], arm = arms$arm[i],
                     bin = pc$hierarchy_bin, malignant = pc$malignant,
                     myeloid = pc$myeloid_class, total = cs)
  list(cnt = cnt, meta = meta)
}
parts <- lapply(seq_len(nrow(arms)), function(i) { message("   ", i, "/", nrow(arms)); grab(i) })
parts <- Filter(Negate(is.null), parts)

# common gene space across the 18 objects (the studies used different references)
common <- Reduce(intersect, lapply(parts, function(p) rownames(p$cnt)))
message("    ", length(common), " genes shared by all ", length(parts), " samples")
cnt  <- do.call(cbind, lapply(parts, function(p) p$cnt[common, , drop = FALSE]))
meta <- rbindlist(lapply(parts, `[[`, "meta"))
rm(parts); invisible(gc())
message("    ", ncol(cnt), " cells x ", nrow(cnt), " genes")

# CP10K + log1p, the standard scRNA normalisation
lognorm <- function(m, cs) {
  sm <- m %*% Matrix::Diagonal(x = NORM_TARGET / cs)
  sm@x <- log1p(sm@x)
  sm
}
LN <- lognorm(cnt, meta$total)
rownames(LN) <- common

meta[, arm := factor(arm, levels = c("TP53-WT", "TP53-mut"))]
# a readable cell-type label: malignant blasts are pulled out of their haematopoietic bin, because
# "which compartment" and "leukaemic or not" are different questions and the UMAP should show both
meta[, celltype := fifelse(myeloid == "macrophage_like", "Macrophage",
                    fifelse(malignant %in% c(1L, TRUE, "1", "TRUE"), "AML blast",
                            fifelse(is.na(bin) | bin == "", "unassigned", bin)))]

## ===== V3. violin =============================================================================
message("[2] V3 violin")
gi <- intersect(VGENES, rownames(LN))
vd <- rbindlist(lapply(gi, function(g)
  data.table(gene = g, expr = as.numeric(LN[g, ]), arm = meta$arm)))
vd[, gene := factor(gene, levels = VGENES)]
# VIOLINS ARE DRAWN OVER EXPRESSING CELLS ONLY, and this is a deliberate, disclosed choice.
# A first version plotted all cells. It was honest and unreadable: 14 of these 16 genes are detected
# in 0.1-7% of cells, so every panel collapsed to a spike at zero and only PGK1 and VIM had a visible
# shape. Restricting to cells with a detected transcript lets the distribution be seen; the detection
# rate is then printed above each violin so the reader still has the base rate the violin hides.
# Both numbers are needed: a gene can move because more cells switch it on, or because the cells that
# already express it express more, and these two panels separate exactly those.
pct <- vd[, .(pct = 100 * mean(expr > 0)), by = .(gene, arm)]
vd  <- vd[expr > 0]
pct <- merge(pct, vd[, .(y = max(expr) * 1.05), by = gene], by = "gene")
pct[, lab := sprintf("%.1f%%", pct)]

fV3 <- ggplot(vd, aes(arm, expr, fill = arm)) +
  ggrastr::rasterise(ggplot2::geom_violin(scale = "width", trim = TRUE, linewidth = 0.3,
                                          colour = "#00000022"), dpi = 200) +
  stat_summary(fun = median, geom = "point", size = 2.2, colour = INK[["primary"]]) +
  geom_text(data = pct, aes(x = arm, y = y, label = lab), inherit.aes = FALSE,
            size = 4, colour = INK[["muted"]], vjust = 0) +
  scale_fill_manual(values = c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]]), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  facet_wrap(~ gene, scales = "free_y", nrow = 4) +
  labs(title = "V3  Per-cell expression behind the paired tests",
       subtitle = "expressing cells only; dot = median; % = share of ALL cells with a detected transcript",
       x = NULL, y = "log-normalised expression") +
  theme_lcc() + theme(panel.grid.major.x = element_blank())
save_fig(fV3, "V3_violin", 14, 13)

## ===== U1. UMAP ===============================================================================
if (nzchar(Sys.getenv("LCC_SKIP_UMAP"))) { message("[3] UMAP skipped (LCC_SKIP_UMAP set)"); quit(save = "no") }
message("[3] U1 UMAP: HVG -> PCA -> Harmony -> UMAP")
# HVG by variance of the log-normalised values, computed sparsely
mu  <- Matrix::rowMeans(LN)
ex2 <- Matrix::rowMeans(LN * LN)
vr  <- ex2 - mu^2
hvg <- names(sort(vr, decreasing = TRUE))[seq_len(min(2000L, length(vr)))]
X <- t(LN[hvg, , drop = FALSE])
X <- scale(as.matrix(X), center = TRUE, scale = FALSE)      # centre only; scaling amplifies noise genes
message("    PCA on ", ncol(X), " HVGs")
pc <- irlba::prcomp_irlba(X, n = 30, center = FALSE, scale. = FALSE)
emb <- pc$x
rm(X); invisible(gc())

# HARMONY IS NOT OPTIONAL HERE. The 18 samples come from 5 studies with different chemistries; an
# uncorrected UMAP separates by study, not by biology, and would say nothing about TP53 status.
# Correcting on dataset (not on arm) leaves the mut-vs-WT contrast intact -- that is the thing being
# looked at, and regressing it out would be circular.
message("    Harmony on dataset")
emb <- harmony::RunHarmony(emb, meta = as.data.frame(meta), vars_use = "dataset",
                           verbose = FALSE, ncores = 1)
message("    UMAP")
um <- uwot::umap(emb, n_neighbors = 30, min_dist = 0.3, metric = "cosine", n_threads = 4)
meta[, `:=`(UMAP1 = um[, 1], UMAP2 = um[, 2])]
fwrite(meta[, .(dataset, sample, arm, celltype, UMAP1, UMAP2)],
       file.path(LCC_TAB_DIR, "U1_umap_coords.csv.gz"))

CT <- meta[, .N, by = celltype][order(-N), celltype]
CT <- c(setdiff(CT, "unassigned"), if ("unassigned" %in% CT) "unassigned")
meta[, celltype := factor(celltype, levels = CT)]
# 8-slot categorical cap: anything past the eighth-largest cell type folds into "other" rather than
# inventing a ninth hue.
if (length(CT) > 8L) {
  keep8 <- CT[1:8]
  meta[, celltype := factor(fifelse(celltype %in% keep8, as.character(celltype), "other"),
                            levels = c(keep8, "other"))]
}
CTCOL <- setNames(c(unname(PAL)[seq_len(min(8L, nlevels(meta$celltype)))], "#c8c8c4")[
                    seq_len(nlevels(meta$celltype))], levels(meta$celltype))

base_umap <- function() list(
  ggplot2::coord_fixed(),
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2"),
  theme_lcc(),
  ggplot2::theme(axis.text = ggplot2::element_blank(),
                 panel.grid.major = ggplot2::element_blank()))

ua <- ggplot(meta[sample(.N)], aes(UMAP1, UMAP2, colour = celltype)) +
  ggrastr::rasterise(geom_point(size = 0.25, alpha = 0.6), dpi = 200) +
  scale_colour_manual(values = CTCOL, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 3)) +
  labs(title = "U1a  Cell type") + base_umap()

ub <- ggplot(meta[sample(.N)], aes(UMAP1, UMAP2, colour = arm)) +
  ggrastr::rasterise(geom_point(size = 0.25, alpha = 0.6), dpi = 200) +
  scale_colour_manual(values = c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]]), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = "U1b  TP53 group") + base_umap()

gi2 <- intersect(UGENES, rownames(LN))
gexp <- rbindlist(lapply(gi2, function(g)
  data.table(gene = g, UMAP1 = meta$UMAP1, UMAP2 = meta$UMAP2, expr = as.numeric(LN[g, ]))))
gexp[, gene := factor(gene, levels = gi2)]
setorder(gexp, expr)                                  # expressing cells drawn last, never buried
uc <- ggplot(gexp, aes(UMAP1, UMAP2, colour = expr)) +
  ggrastr::rasterise(geom_point(size = 0.25), dpi = 200) +
  scale_colour_gradient(low = "#e9e9e6", high = PAL[["aqua"]], name = "log-norm", n.breaks = 3) +
  facet_wrap(~ gene, nrow = 1) +
  labs(title = "U1c  Expression") + base_umap() +
  theme(legend.key.width = unit(2.2, "lines"))

save_fig((ua | ub) / uc + patchwork::plot_layout(heights = c(1.35, 1)), "U1_umap", 14, 12)
message("[done]")
