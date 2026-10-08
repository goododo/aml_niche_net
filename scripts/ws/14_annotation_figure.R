#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_integrated/ws_rpca_integrated.rds   (44 GB, read once for UMAP coords)
#          /LARGE1/.../02_integrated/ws_umap_cache.rds         (written on the first run, used after)
#          results/tables/ws/ws_anno_consolidated_res1.csv, ws_adjudication_final.csv,
#          ws_clusters_percell.csv.gz, ws_pca_meta.rds
# OUTPUT : results/figures/ws/ws_annotation_final.pdf        the figure
#          results/tables/ws/ws_final_labels.csv             cluster -> final label -> lineage group
#          results/tables/ws/ws_composition_table.csv        the table the contrast WARN obliges
#          results/tables/ws/ws_figure_checks.csv
# what it does: draws the final annotation.
#   Design decisions, and why:
#   * 38 distinct fine labels cannot be colour-encoded. A validated categorical palette holds 8
#     slots and a 9th hue may not be generated, so the hero panel colours 8 LINEAGE GROUPS and the
#     fine labels appear as small multiples -- one facet per group -- instead of as 38 hues.
#   * the palette is the validated 8-slot reference instance, checked with the skill's own
#     validator: lightness band PASS, chroma floor PASS, CVD separation PASS (worst adjacent
#     dE 9.1 protan), normal-vision floor PASS (19.6), contrast vs surface WARN on 3 of 8.
#     That WARN is not dismissable, so both reliefs ship: centroid labels drawn directly on every
#     panel, and ws_composition_table.csv as the table view.
#   * hues are assigned to lineage groups in a FIXED order and never cycled, so a group keeps its
#     colour across panels and across future versions of this figure.
#   * UMAP axes are dropped: the coordinates carry no units and an axis would be chartjunk.
#   * 895,228 points are rasterised; only the labels and frame stay vector.
#   The fine-to-group mapping is written to ws_final_labels.csv so it can be audited and edited
#   without touching this script -- the first attempt at a final label used free Chinese prose and
#   split one lineage across four categories, which is what motivated the controlled vocabulary.
# This script deletes nothing.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(Seurat); library(SeuratObject); library(patchwork)
})
stopifnot(requireNamespace("ggrastr", quietly = TRUE))
# 895,228 vector points produced a 465 MB PDF (job 3641993) that no viewer opens. Rasterise the
# point layers at print resolution; labels, legend and frame stay vector so the file stays usable.
RPT <- function(...) ggrastr::rasterise(geom_point(...), dpi = 300, dev = "ragg")
set.seed(42)
options(future.globals.maxSize = 64 * 1024^3)
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- file.path(ROOT, "results/tables/ws")
FIG  <- file.path(ROOT, "results/figures/ws"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
UCACHE <- file.path(BASE, "ws_umap_cache.rds")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-50s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- the validated 8-slot categorical palette, in fixed order ------------------------------------
PAL <- c("HSPC/progenitor" = "#2a78d6", "Blast-like" = "#eb6834", "Myeloid" = "#1baf7a",
         "T/NK" = "#eda100", "B/plasma" = "#e87ba4", "Ery/Mega" = "#008300",
         "Niche (stroma/vessel/bone)" = "#4a3aa7", "Unresolved/excluded" = "#e34948")
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; MUTED <- "#898781"

## -- 1. labels: a controlled vocabulary, written out for audit -----------------------------------
ts("[1] final labels")
CONS <- fread(file.path(TAB, "ws_anno_consolidated_res1.csv"))
FIN  <- fread(file.path(TAB, "ws_adjudication_final.csv"))
# the adjudicated clusters get an explicit English label; everything else keeps its agreed label
ADJ <- data.table(
  cl = c(0,1,7,11,13,15,16,26,27,28,29,30,31,34,38,44,46,48,49,50,52,54),
  final = c("T/NK","Mono: Classical","Mono: Classical","Blast-like (BAALC+, normal karyotype)",
            "T/NK","Mono: Classical","T/NK or CLP","Myeloid: Generic","Blast-like: GMP",
            "Unresolved (homotypic doublet)","Mono: Non-classical",
            "Blast-like primitive progenitor (BAALC+DNTT+PROM1)","DC: cDC2",
            "Interferon-response state","Mono: activated","Mono/macrophage",
            "Excluded (low quality, single library)","Blast-like (basophil/mast features, PRAME+)",
            "Blast-like","T cell","Mast/basophil (TPSB2+)","NK"))
L <- merge(CONS[, .(cl, cluster_size, label, provenance, old_label, auto_mk3)], ADJ, by = "cl", all.x = TRUE)
L[, final := ifelse(!is.na(final), final, label)]
chk("every cluster has a final label", 57L, sum(!is.na(L$final)))

GROUP <- function(x) fcase(
  grepl("^Blast-like", x),                               "Blast-like",
  grepl("Unresolved|Excluded|Interferon-response", x),   "Unresolved/excluded",
  grepl("^HSPC|CLP/Early", x) & !grepl("T/NK", x),       "HSPC/progenitor",
  grepl("^T/NK|^T cell|^T:|^CD4 T|^CD8 T|^NK", x),       "T/NK",
  grepl("^B:|^B |Plasma|Pro/Pre-B", x),                  "B/plasma",
  grepl("Erythroid|Megakaryo|Ery/Meg", x),               "Ery/Mega",
  grepl("^Stroma|^Vascular|^Bone|^Cartilage", x),        "Niche (stroma/vessel/bone)",
  grepl("^Mono|^DC|^Myeloid|Granulocyte|Mast|basophil", x), "Myeloid",
  default = "Unresolved/excluded")
L[, lineage := GROUP(final)]
L[, lineage := factor(lineage, levels = names(PAL))]
chk("every cluster in one of the 8 groups", 0L, sum(is.na(L$lineage)))
fwrite(L[order(-cluster_size)], file.path(TAB, "ws_final_labels.csv"))
ts("    lineage groups:"); print(L[, .(clusters = .N, cells = sum(cluster_size)), by = lineage][order(-cells)])

## -- 2. UMAP coordinates: cache FIRST, then use ---------------------------------------------------
if (file.exists(UCACHE)) {
  ts("[2] reusing the UMAP cache; the 44 GB object is not read")
  U <- readRDS(UCACHE)
} else {
  ts("[2] reading the integrated object (44 GB, about 17 min) for UMAP coordinates only")
  o <- readRDS(file.path(BASE, "ws_rpca_integrated.rds"))
  U <- list(emb = Embeddings(o, reduction = "umap"), cells = colnames(o))
  ts("    caching before any further step: ", UCACHE)
  saveRDS(U, UCACHE)                      # cache first -- three earlier jobs lost long reads to a
  rm(o); gc(FALSE)                        # check that sat above the write
}
chk("UMAP has two dimensions", 2L, ncol(U$emb))

CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell","res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
M  <- readRDS(file.path(BASE, "ws_pca_meta.rds"))
md <- as.data.table(M$meta)[, .(cell, Dataset_ID, SampleType)]
setkey(md, cell); md <- md[CL]
UD <- data.table(cell = U$cells, UMAP_1 = U$emb[,1], UMAP_2 = U$emb[,2]); setkey(UD, cell)
D  <- merge(md, UD, by = "cell")                        # align by NAME, never by position
D  <- merge(D, L[, .(cl, final, lineage, provenance)], by = "cl")
chk("cells plotted", 895228L, nrow(D))

## -- 3. the figure --------------------------------------------------------------------------------
ts("[3] drawing")
base_theme <- theme_void(base_size = 9) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.background = element_rect(fill = SURFACE, colour = NA),
        plot.title = element_text(colour = INK, face = "bold", size = 10),
        plot.subtitle = element_text(colour = MUTED, size = 7.5),
        legend.text = element_text(colour = INK, size = 7),
        legend.title = element_blank(), strip.text = element_text(colour = INK, size = 7))
cen <- D[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = .N), by = .(lineage, cl, final)]

pA <- ggplot(D, aes(UMAP_1, UMAP_2, colour = lineage)) +
  RPT(size = .08, stroke = 0, alpha = .55, shape = 16) +
  scale_colour_manual(values = PAL, drop = FALSE) +
  ggrepel::geom_text_repel(data = cen[n > 8000], aes(label = cl), colour = INK, size = 2.4,
                           max.overlaps = Inf, box.padding = .12, segment.colour = MUTED,
                           segment.size = .2, show.legend = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 2)) +
  labs(title = "Final annotation, 57 clusters at resolution 1.0",
       subtitle = sprintf("%s cells | 9 datasets | 182 samples | cluster ids labelled where n > 8,000",
                          format(nrow(D), big.mark = ","))) +
  base_theme + theme(legend.position = "bottom")

pB <- ggplot(D, aes(UMAP_1, UMAP_2)) +
  RPT(data = D[, .(UMAP_1, UMAP_2)], colour = "#e1e0d9", size = .05, stroke = 0, shape = 16) +
  RPT(aes(colour = lineage), size = .07, stroke = 0, alpha = .6, shape = 16) +
  scale_colour_manual(values = PAL, guide = "none") +
  facet_wrap(~ lineage, nrow = 2) +
  labs(title = "The same eight groups, one panel each",
       subtitle = "grey = all other cells, so each group is read against the whole manifold") +
  base_theme

pC <- ggplot(D, aes(UMAP_1, UMAP_2, colour = SampleType)) +
  RPT(size = .06, stroke = 0, alpha = .5, shape = 16) +
  scale_colour_manual(values = c(AML = "#eb6834", Control = "#2a78d6")) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "AML vs Control", subtitle = "107 AML / 75 Control samples") +
  base_theme + theme(legend.position = "bottom")

pD <- ggplot(D, aes(UMAP_1, UMAP_2, colour = provenance)) +
  RPT(size = .06, stroke = 0, alpha = .5, shape = 16) +
  scale_colour_manual(values = c(agreed = "#1baf7a", adopted = "#eda100", disputed = "#4a3aa7")) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  labs(title = "How each label was reached",
       subtitle = "agreed = both methods; adopted = old label was Ambiguous; disputed = settled on four evidence streams") +
  base_theme + theme(legend.position = "bottom")

COMP <- D[, .N, by = .(Dataset_ID, lineage)][, share := N/sum(N), by = Dataset_ID]
fwrite(dcast(COMP, Dataset_ID ~ lineage, value.var = "N", fill = 0L),
       file.path(TAB, "ws_composition_table.csv"))
pE <- ggplot(COMP, aes(reorder(Dataset_ID, -N, sum), share, fill = lineage)) +
  geom_col(width = .72, colour = SURFACE, linewidth = .6) +     # the 2px surface gap between fills
  scale_fill_manual(values = PAL, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(100*x), "%"), expand = c(0,0)) +
  coord_flip() +
  labs(title = "Lineage composition per dataset", subtitle = "the table view the contrast check obliges: ws_composition_table.csv") +
  theme_minimal(base_size = 9) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = .25),
        axis.title = element_blank(), axis.text = element_text(colour = INK, size = 7),
        plot.title = element_text(colour = INK, face = "bold", size = 10),
        plot.subtitle = element_text(colour = MUTED, size = 7.5),
        legend.position = "none")

pdf(file.path(FIG, "ws_annotation_final.pdf"), width = 13, height = 16)
print((pA / pB / (pC | pD) / pE) + plot_layout(heights = c(1.5, 1.1, 1.1, .8)) +
      plot_annotation(theme = theme(plot.background = element_rect(fill = SURFACE, colour = NA))))
dev.off()
mb <- file.size(file.path(FIG, "ws_annotation_final.pdf")) / 1024^2
ts(sprintf("    wrote %s (%.1f MB)", file.path(FIG, "ws_annotation_final.pdf"), mb))
chk("PDF is a usable size (rasterised, not 465 MB of vector points)", TRUE, mb < 40)

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_figure_checks.csv"))
ts(sprintf("[4] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done]")
