#!/usr/bin/env Rscript
# 14_ccc_compare.R ----
# Cell-cell communication, TP53-mut vs TP53-WT, on the matched sample design from 11.
#
# THE PI'S QUESTION was leukaemia-cell <-> stromal-cell interaction. That question is dead in this
# cohort and 10 measured why: whole-marrow aspirates carry 0.028% stromal cells, three times the
# false-positive floor of a library that contains none by construction, and NO TP53-mut sample has a
# usable stromal population. So the stromal node does not exist here -- consistent with the main
# line, which already sets Stromal in_ccc_graph = FALSE. What CAN be asked is the haematopoietic
# half of the question: HSC_MPP (the PI's "HSC") and Mono_DC (which carries the macrophages).
#
# PAIRING: inherited from 11's matched design, with two documented repairs.
#   pair 8  3853_R has no CellChat object -> re-matched to another GSE227903 Relapse WT that does.
#   pair 1  AML328-D0 has no CellChat object and it is the MUT side, so no substitution is possible;
#           the pair is dropped. AML328 is still represented by its other three timepoints.
#   Both missing objects were excluded as "below_min_cells" against a 500-cell floor that
#   config_ccc.R itself records as since lowered to 300 -- a stale-manifest artefact, not biology.
#   They are NOT silently re-admitted here: re-running CellChat on them would use different
#   inclusion logic from the other 16 objects, and a mixed provenance is worse than 8 clean pairs.
#
# TWO NORMALISATIONS, both reported, because neither is safe alone:
#   absolute  sum of significant communication probability. Sensitive to cell number and depth.
#   relative  the same divided by the sample's total. Removes that, but forces edges to compete --
#             a genuine global rise then reads as other edges falling. Agreement between them is
#             the only thing that should be believed.
#
# INPUT  : LCC_TAB_DIR/11_matched_design.csv, 09_tp53_groups.csv
#          CCC_GRAPH_DIR/<ds>/<sample>__cellchat.rds
# OUTPUT : LCC_TAB_DIR/14_ccc_pairs.csv        the pairing actually used, with repairs flagged
#          LCC_TAB_DIR/14_ccc_edges.csv        per sample x source x target strength
#          LCC_TAB_DIR/14_ccc_edge_tests.csv   per edge, paired test, both normalisations
#          LCC_TAB_DIR/14_ccc_pathway_tests.csv per pathway x focus node
#          LCC_FIG_DIR/G1..G4                  figures
# Usage  : Rscript LCC_proj/scripts/14_ccc_compare.R

suppressPackageStartupMessages({
  library(data.table); library(here); library(CellChat); library(ggplot2); library(ggrepel)
  library(scales); library(patchwork)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
set.seed(SEED)

CCC_GRAPH_DIR <- file.path(LARGE1_DIR, "05_ccc_graphs")
FOCUS <- c("HSC_MPP", "Mono_DC")      # the PI's HSC and macrophage compartments
PVAL  <- 0.05
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
cc_path <- function(ds, s) file.path(CCC_GRAPH_DIR, ds, paste0(s, "__cellchat.rds"))

## -- Step 1. repair the pairing ----
message("[1] pairing")
des <- fread(file.path(LCC_TAB_DIR, "11_matched_design.csv"))
grp <- fread(file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))
wt_pool <- grp[timepoint != "Healthy" & tp53_tier %in% c("WT_genotyped", "WT_presumed")]
wt_pool[, has_cc := file.exists(cc_path(dataset, sample))]

des[, `:=`(mut_cc = file.exists(cc_path(dataset, mut_sample)),
           wt_cc  = file.exists(cc_path(dataset, wt_sample)), repair = "")]
used <- des$wt_sample
for (i in which(!des$wt_cc & des$mut_cc)) {
  cand <- wt_pool[dataset == des$dataset[i] & timepoint == des$timepoint[i] &
                    has_cc == TRUE & !(sample %in% used)]
  if (!nrow(cand)) next
  setorder(cand, sample)
  des[i, `:=`(wt_sample = cand$sample[1], wt_cc = TRUE,
              repair = paste0("WT re-matched (", des$wt_sample[i], " has no CellChat)"))]
  used <- c(used, cand$sample[1])
}
des[mut_cc == FALSE, repair := paste0("DROPPED: mut side ", mut_sample, " has no CellChat")]
fwrite_safe(des, file.path(LCC_TAB_DIR, "14_ccc_pairs.csv"))
print(des[, .(pair_id, dataset, timepoint, mut_sample, wt_sample, mut_cc, wt_cc, repair)])
use <- des[mut_cc == TRUE & wt_cc == TRUE]
message("    usable pairs: ", nrow(use), " / ", nrow(des))
arms <- rbind(use[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              use[, .(pair_id, dataset, sample = wt_sample, arm = "TP53-WT")])

## -- Step 2. extract edge and pathway strength ----
message("[2] reading ", nrow(arms), " CellChat objects")
extract <- function(ds, s) {
  cc <- readRDS(cc_path(ds, s))
  pr <- cc@net$prob; pv <- cc@net$pval
  sig <- pr * (pv < PVAL)                       # significant communication probability only
  lv <- levels(cc@idents); nc <- as.integer(table(cc@idents)[lv])
  lrp <- cc@LR$LRsig$pathway_name
  edge <- as.data.table(as.table(apply(sig, c(1, 2), sum)))
  setnames(edge, c("source", "target", "strength"))
  nlr <- as.data.table(as.table(apply(sig > 0, c(1, 2), sum)))
  setnames(nlr, c("source", "target", "n_lr"))
  edge <- merge(edge, nlr, by = c("source", "target"))
  # pathway x (source,target), kept only for the focus nodes to stay small
  ps <- rbindlist(lapply(unique(lrp), function(p) {
    k <- which(lrp == p); m <- apply(sig[, , k, drop = FALSE], c(1, 2), sum)
    dimnames(m) <- list(lv, lv)
    d <- as.data.table(as.table(m)); setnames(d, c("source", "target", "strength"))
    d[, pathway := p][]
  }), fill = TRUE)
  list(edge = edge[, `:=`(dataset = ds, sample = s)],
       path = ps[, `:=`(dataset = ds, sample = s)],
       nodes = data.table(dataset = ds, sample = s, node = lv, n_cells = nc))
}
got <- lapply(seq_len(nrow(arms)), function(i) {
  if (i %% 4 == 0) message("    ", i, "/", nrow(arms))
  tryCatch(extract(arms$dataset[i], arms$sample[i]),
           error = function(e) { message("    [warn] ", arms$sample[i], ": ", conditionMessage(e)); NULL })
})
got <- got[!vapply(got, is.null, TRUE)]
edges <- rbindlist(lapply(got, `[[`, "edge"), fill = TRUE)
paths <- rbindlist(lapply(got, `[[`, "path"), fill = TRUE)
nodes <- rbindlist(lapply(got, `[[`, "nodes"), fill = TRUE)
for (d in list(edges, paths, nodes)) d[arms, `:=`(pair_id = i.pair_id, arm = i.arm), on = c("dataset", "sample")]

# relative normalisation: share of the sample's total significant communication
edges[, total := sum(strength), by = .(dataset, sample)]
edges[, rel := fifelse(total > 0, strength / total, NA_real_)]
paths[, total := sum(strength), by = .(dataset, sample)]
paths[, rel := fifelse(total > 0, strength / total, NA_real_)]
# A node absent from a sample yields structural zeros; those pairs are dropped per edge, not
# imputed, so an edge is tested only where both members of the pair could express it.
nodes_present <- nodes[n_cells >= 10L, .(dataset, sample, node)]
ok_edge <- merge(edges, nodes_present[, .(dataset, sample, source = node)],
                 by = c("dataset", "sample", "source"))
ok_edge <- merge(ok_edge, nodes_present[, .(dataset, sample, target = node)],
                 by = c("dataset", "sample", "target"))
fwrite_safe(ok_edge, file.path(LCC_TAB_DIR, "14_ccc_edges.csv"))

## -- Step 3. the paired test (same machinery as 11) ----
message("[3] paired tests")
paired <- function(dm, dw) {
  ok <- !is.na(dm) & !is.na(dw); d <- dm[ok] - dw[ok]; n <- sum(ok); d <- d[d != 0]
  k <- length(d); if (n < 4L) return(list(NA_real_, NA_real_, n, NA_real_))
  if (!k) return(list(0, 1, n, 0.5))
  r <- rank(abs(d)); w <- sum(r[d > 0]); mx <- sum(r)
  ways <- numeric(mx + 1L); ways[1L] <- 1
  for (i in seq_len(k)) ways <- ways + c(numeric(r[i]), ways[seq_len(mx + 1L - r[i])])
  list(median(d), sum(ways[seq(floor(w) + 1L, mx + 1L)]) / sum(ways), n, mean(dm[ok] > dw[ok]))
}
test_block <- function(dt, value, keys) {
  w <- dcast(dt, as.formula(paste(paste(c(keys, "pair_id"), collapse = "+"), "~ arm")),
             value.var = value)
  setnames(w, c("TP53-mut", "TP53-WT"), c("m", "v"))
  w[, { r <- paired(m, v)
        .(n_pairs = r[[3]], median_delta = r[[1]], p_higher = r[[2]], frac_pairs_up = r[[4]],
          mean_mut = mean(m, na.rm = TRUE), mean_WT = mean(v, na.rm = TRUE)) }, by = keys]
}
et <- rbind(cbind(norm = "absolute", test_block(ok_edge, "strength", c("source", "target"))),
            cbind(norm = "relative", test_block(ok_edge, "rel",      c("source", "target"))))
et[, fdr := p.adjust(p_higher, "BH"), by = norm]
setorder(et, norm, p_higher)
fwrite_safe(et, file.path(LCC_TAB_DIR, "14_ccc_edge_tests.csv"))

# pathway level, restricted to edges touching a focus node
pf <- paths[source %in% FOCUS | target %in% FOCUS]
pf[, role := fifelse(source %in% FOCUS & target %in% FOCUS, "within focus",
              fifelse(source %in% FOCUS, "focus -> other", "other -> focus"))]
pf2 <- pf[, .(strength = sum(strength), rel = sum(rel)), by = .(pair_id, arm, pathway, role)]
pt <- rbind(cbind(norm = "absolute", test_block(pf2, "strength", c("pathway", "role"))),
            cbind(norm = "relative", test_block(pf2, "rel",      c("pathway", "role"))))
pt <- pt[n_pairs >= 5]
pt[, fdr := p.adjust(p_higher, "BH"), by = .(norm, role)]
setorder(pt, norm, p_higher)
fwrite_safe(pt, file.path(LCC_TAB_DIR, "14_ccc_pathway_tests.csv"))

## -- Step 4. figures ----
message("[4] figures")
BINS <- c("HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma")

# G1 -- the 7x7 signalling map, both arms and their difference
# MAGNITUDE and DIFFERENCE need different colour jobs and therefore different scales. A single
# diverging ramp spanning both, as a first version had, washed the two magnitude panels out to near
# white because they are all small positives -- unreadable, and diverging is the wrong job for a
# quantity that cannot go negative. Sequential single hue for magnitude, diverging for the signed
# difference, combined with patchwork so each keeps its own legend.
# MEDIAN, not mean: every test in this project is rank-based, and one sample's B_Plasma edge is
# large enough to set a mean on its own.
g1d <- ok_edge[, .(v = median(rel, na.rm = TRUE)), by = .(arm, source, target)]
g1w <- dcast(g1d, source + target ~ arm, value.var = "v")
setnames(g1w, c("TP53-mut", "TP53-WT"), c("m", "w"))
g1w[is.na(m), m := 0][is.na(w), w := 0][, diff := m - w]
g1d[, `:=`(source = factor(source, BINS), target = factor(target, BINS),
           arm = factor(arm, c("TP53-mut", "TP53-WT")))]
g1w[, `:=`(source = factor(source, BINS), target = factor(target, BINS))]
base_tile <- function(d, fillvar) ggplot(d, aes(target, source, fill = .data[[fillvar]])) +
  geom_tile(colour = SURF, linewidth = 1.2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_lcc() + theme(panel.grid.major = element_blank(),
                      legend.key.width = unit(2.2, "lines"))
p_mag <- base_tile(g1d, "v") + facet_wrap(~ arm, nrow = 1) +
  # Sequential ramp is deliberately NOT blue: blue carries "mut lower" in the diverging panel
  # beside it, and one hue must not mean "high" in one panel and "negative" in the next.
  scale_fill_gradient(low = "#eaf6f1", high = PAL[["aqua"]], na.value = GRID, n.breaks = 3,
                      name = "share of sample total", labels = percent_format(accuracy = 1)) +
  labs(title = "Signalling between the seven haematopoietic compartments",
       subtitle = paste0(nrow(use), " matched pairs"),
       x = NULL, y = "source (sender)")
p_dif <- base_tile(g1w, "diff") +
  scale_fill_gradient2(low = PAL[["blue"]], mid = "#f2f1ec", high = PAL[["orange"]], midpoint = 0,
                       na.value = GRID, name = "mut - WT", n.breaks = 3,
                       labels = percent_format(accuracy = 1)) +
  labs(title = "difference (mut - WT)", subtitle = NULL, x = NULL, y = NULL)
g1 <- (p_mag | p_dif) + patchwork::plot_layout(widths = c(2, 1)) 
save_fig(g1, "G1_ccc_map", 15, 5.4)

# G2 -- per-pair trajectories for the focus nodes: the raw evidence behind any claim
g2 <- ok_edge[source %in% FOCUS | target %in% FOCUS]
g2[, role := fifelse(source %in% FOCUS & target %in% FOCUS, "within HSC / Mono_DC",
             fifelse(source %in% FOCUS, "HSC / Mono_DC  ->  other", "other  ->  HSC / Mono_DC"))]
g2s <- g2[, .(rel = sum(rel, na.rm = TRUE)), by = .(pair_id, arm, role)]
g2p <- ggplot(g2s, aes(arm, rel, group = pair_id, colour = factor(pair_id))) +
  geom_line(linewidth = 0.8, alpha = 0.85) + geom_point(size = 2.4, stroke = 0) +
  facet_wrap(~ role) +
  scale_colour_manual(values = unname(PAL), name = "pair") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Every pair, plotted: HSC and Mono_DC signalling share",
       subtitle = "One line per matched pair. Crossing lines are what a non-significant result looks like.",
       x = NULL, y = "share of the sample's total communication",
       caption = paste("Plotted rather than summarised on purpose: with 8 pairs a mean and a p-value",
                       "hide whether one pair\ncarries the result. Paired one-sided exact",
                       "signed-rank; smallest attainable p at 8 pairs is 1/256.")) +
  theme_lcc()
save_fig(g2p, "G2_focus_pairs", 11, 5.5)

# G3 -- pathway effects, with the fibrosis-relevant ones named
FIB <- c("TGFb", "PDGF", "BMP", "ACTIVIN", "IL6", "OSM", "SPP1", "ANGPT", "ANGPTL", "VEGF",
         "CXCL", "CCL", "TNF", "IL1", "COLLAGEN", "FN1", "LAMININ", "THBS", "CD226", "NECTIN", "PVR")
g3 <- pt[norm == "relative" & role != "within focus"]
g3[, fib := pathway %in% FIB]
g3[, lab := fifelse(p_higher <= 0.15 | (fib & p_higher <= 0.4), pathway, NA_character_)]
g3p <- ggplot(g3, aes(median_delta, -log10(p_higher))) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.4, colour = INK[["muted"]]) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = INK[["muted"]]) +
  geom_point(aes(colour = fib), size = 2.2, stroke = 0, alpha = 0.85) +
  geom_text_repel(aes(label = lab, colour = fib), size = 4.2, seed = SEED, max.overlaps = 25,
                  segment.colour = INK[["muted"]], segment.size = 0.3, show.legend = FALSE) +
  facet_wrap(~ role) +
  # repelled labels at the extremes were being clipped by the panel edge
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_colour_manual(values = c(`FALSE` = INK[["muted"]], `TRUE` = PAL[["orange"]]),
                      labels = c("other pathway", "fibrosis / inflammation relevant"), name = NULL) +
  labs(title = "Pathway-level signalling into and out of HSC and Mono_DC",
       subtitle = "One point per CellChat pathway. The dashed line is p = 0.05, uncorrected.",
       x = "median paired difference in signalling share (TP53-mut - TP53-WT)",
       y = expression(-log[10]~"p  (one-sided, paired)"),
       caption = paste("Relative normalisation. Nothing here survives FDR correction across",
                       "pathways -- the labels mark where to\nlook next, not findings. With 8 pairs",
                       "the floor is p = 1/256, so a pathway must move in every pair to clear",
                       "0.05\nafter correction.")) +
  theme_lcc()
save_fig(g3p, "G3_ccc_pathways", 11.5, 6)

# G4 -- the specific edges the PI asked about
g4 <- ok_edge[(source == "Mono_DC" & target == "HSC_MPP") | (source == "HSC_MPP" & target == "Mono_DC") |
              (source == "HSC_MPP" & target == "HSC_MPP") | (source == "Mono_DC" & target == "Mono_DC")]
g4[, edge := paste(source, "->", target)]
g4p <- ggplot(g4, aes(arm, rel, group = pair_id, colour = factor(pair_id))) +
  geom_line(linewidth = 0.8, alpha = 0.85) + geom_point(size = 2.4, stroke = 0) +
  facet_wrap(~ edge, nrow = 1) +
  scale_colour_manual(values = unname(PAL), name = "pair") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "The four HSC / Mono_DC edges, pair by pair",
       subtitle = "The leukaemia-stem-cell to macrophage-compartment axis, which is the part of the PI's question the data supports.",
       x = NULL, y = "share of the sample's total communication",
       caption = "Mono_DC is the compartment carrying macrophages; a dedicated macrophage node was not built because only 2,269 macrophage-like cells exist across all AML samples.") +
  theme_lcc()
save_fig(g4p, "G4_hsc_mono_edges", 13, 4.8)

message("[5] summary")
cat("\n-- edges, ranked by the paired test (relative normalisation) --\n")
print(et[norm == "relative"][order(p_higher)][1:12,
      .(source, target, n_pairs, up = paste0(round(frac_pairs_up * n_pairs), "/", n_pairs),
        d = signif(median_delta, 3), p = signif(p_higher, 3), fdr = signif(fdr, 3))])
cat("\n-- do the two normalisations agree on the top edges? --\n")
print(dcast(et[, .(source, target, norm, p_higher)], source + target ~ norm,
            value.var = "p_higher")[order(relative)][1:10])
cat("\n-- pathways touching HSC/Mono_DC, best 12 (relative) --\n")
print(pt[norm == "relative"][order(p_higher)][1:12,
      .(pathway, role, n_pairs, up = paste0(round(frac_pairs_up * n_pairs), "/", n_pairs),
        p = signif(p_higher, 3), fdr = signif(fdr, 3))])
cat("\n-- anything surviving FDR 0.10 anywhere --\n")
print(rbind(et[fdr <= 0.10, .(what = paste(source, "->", target), norm, p_higher, fdr)],
            pt[fdr <= 0.10, .(what = paste(pathway, role), norm, p_higher, fdr)]))
message("[done] 14_ccc_compare")
