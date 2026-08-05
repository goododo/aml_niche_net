#!/usr/bin/env Rscript
# 20_ccc_macrophage.R ----
# The macrophage half of the cell-cell-communication brief, and the ligand-receptor layer that 14
# stopped short of.
#
# 14 tested communication at the EDGE level (source compartment -> target compartment) and at the
# PATHWAY level. It answered "does the HSC / Mono_DC axis move", not "which signal moves", and it
# never addressed the macrophage node the brief named, because the CellChat graphs carry seven
# haematopoietic compartments and macrophage is not one of them -- it sits inside Mono_DC.
#
# This script closes both gaps and does not hide the one that cannot be closed:
#
#   C1  MACROPHAGE NODE AVAILABILITY. A strict macrophage node is not testable in this cohort and
#       the figure shows exactly why, pair by pair, against CellChat's 10-cell floor. This is a
#       negative result and it is drawn rather than written in a footnote, because "we skipped it"
#       and "we measured that it cannot be done" are different statements and only one is true.
#
#   C2  LIGAND-RECEPTOR LEVEL. Every LR pair on the macrophage-compartment axis and the HSC axis,
#       paired mut vs WT, ranked by effect. This is the "which signal" answer.
#
#   C3  AXIS TOTALS. The same evidence as a bar with every pair overlaid, so the reader sees the
#       spread the bar height is hiding.
#
# NORMALISATION: relative (share of the sample's total significant communication) is the primary
# readout, as in 14 -- absolute strength tracks cell number and sequencing depth, which differ by
# study. Absolute is computed and written to the table so the two can be checked against each other.
#
# INPUT  : LCC_TAB_DIR/14_ccc_pairs.csv        the 8 usable pairs, repairs already applied
#          LCC_TAB_DIR/04_myeloid_gate.csv     per-sample macrophage / monocyte / Mono_DC counts
#          LARGE1_DIR/05_ccc_graphs/<ds>/<sample>__cellchat.rds
# OUTPUT : LCC_TAB_DIR/20_macrophage_node_availability.csv
#          LCC_TAB_DIR/20_ccc_lr_tests.csv
#          LCC_TAB_DIR/20_ccc_axis_totals.csv
#          LCC_FIG_DIR/{C1_macrophage_node, C2_ccc_ligand_receptor, C3_ccc_axis_totals}.{png,pdf}
# Usage  : Rscript LCC_proj/scripts/20_ccc_macrophage.R

suppressPackageStartupMessages({
  library(data.table); library(here); library(CellChat); library(ggplot2); library(ggrepel)
  library(scales); library(patchwork)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
source(here::here("LCC_proj", "scripts", "theme_lcc.R"))
set.seed(SEED)

CCC_GRAPH_DIR <- file.path(LARGE1_DIR, "05_ccc_graphs")
PVAL     <- 0.05     # CellChat's own permutation p, the same cut 14 used
MIN_CELL <- 10L      # CellChat's min.cells floor: below this a node is not modelled
ARMCOL   <- c(`TP53-WT` = PAL[["blue"]], `TP53-mut` = PAL[["orange"]])
cc_path  <- function(ds, s) file.path(CCC_GRAPH_DIR, ds, paste0(s, "__cellchat.rds"))

des  <- fread(file.path(LCC_TAB_DIR, "14_ccc_pairs.csv"))
use  <- des[mut_cc == TRUE & wt_cc == TRUE]
arms <- rbind(use[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              use[, .(pair_id, dataset, sample = wt_sample,  arm = "TP53-WT")])
arms[, arm := factor(arm, levels = c("TP53-WT", "TP53-mut"))]
message("[0] ", nrow(use), " pairs, ", nrow(arms), " samples")

## ===== C1. can macrophage be its own node? ====================================================
message("[1] C1 macrophage node availability")
my <- fread(file.path(LCC_TAB_DIR, "04_myeloid_gate.csv"))
av <- merge(arms, my[, .(dataset, sample, n_cells_total, n_mono_dc,
                         n_macrophage_like, n_monocyte_like)],
            by = c("dataset", "sample"), all.x = TRUE)
setorder(av, pair_id, arm)
fwrite_safe(av, file.path(LCC_TAB_DIR, "20_macrophage_node_availability.csv"))

# THE VERDICT, computed rather than asserted: a node is testable in a PAIR only when BOTH sides
# clear the floor, because the test is paired -- one arm alone contributes nothing.
pw <- dcast(av, pair_id ~ arm, value.var = c("n_macrophage_like", "n_mono_dc"))
n_mac_pairs  <- sum(pw[["n_macrophage_like_TP53-mut"]] >= MIN_CELL &
                    pw[["n_macrophage_like_TP53-WT"]]  >= MIN_CELL, na.rm = TRUE)
n_mono_pairs <- sum(pw[["n_mono_dc_TP53-mut"]] >= MIN_CELL &
                    pw[["n_mono_dc_TP53-WT"]]  >= MIN_CELL, na.rm = TRUE)
message("    strict macrophage node usable in ", n_mac_pairs, "/", nrow(use), " pairs")
message("    Mono_DC compartment usable in    ", n_mono_pairs, "/", nrow(use), " pairs")

c1d <- melt(av, id.vars = c("pair_id", "arm", "sample"),
            measure.vars = c("n_macrophage_like", "n_mono_dc"),
            variable.name = "node", value.name = "n")
c1d[, node := factor(node, levels = c("n_macrophage_like", "n_mono_dc"),
                     labels = c("Macrophage-like (strict call)",
                                "Mono_DC compartment (what CellChat models)"))]
# zero cannot be drawn on a log axis and dropping those rows would delete the whole point of the
# figure, so zeros are pinned to a half-unit floor and labelled "0" on the axis.
FLOOR <- 0.5
c1d[, nplot := pmax(n, FLOOR)]
c1seg <- dcast(c1d, pair_id + node ~ arm, value.var = "nplot")
setnames(c1seg, c("TP53-WT", "TP53-mut"), c("wt", "mut"))

fC1 <- ggplot(c1d, aes(nplot, factor(pair_id))) +
  geom_vline(xintercept = MIN_CELL, linetype = "22", linewidth = 0.6, colour = PAL[["red"]]) +
  geom_segment(data = c1seg, aes(x = wt, xend = mut, y = factor(pair_id), yend = factor(pair_id)),
               inherit.aes = FALSE, colour = GRID, linewidth = 1.4) +
  geom_point(aes(colour = arm), size = 3.6, stroke = 0) +
  facet_wrap(~ node, nrow = 1) +
  scale_x_log10(breaks = c(FLOOR, 1, 10, 100, 1000, 10000),
                labels = c("0", "1", "10", "100", "1,000", "10,000")) +
  annotation_logticks(sides = "b", outside = FALSE, size = 0.3, colour = INK[["muted"]]) +
  scale_colour_manual(values = ARMCOL, name = NULL) +
  labs(title = "C1  Why macrophage cannot be its own communication node here",
       subtitle = paste0("Dashed line = the 10-cell floor. Both dots must clear it: ",
                         n_mac_pairs, "/", nrow(use), " pairs do."),
       x = "cells in the sample (log scale)", y = "matched pair") +
  theme_lcc() + theme(panel.grid.major.y = element_blank())
save_fig(fC1, "C1_macrophage_node", 13, 5.6)

## ===== read the graphs, ligand-receptor resolution =============================================
message("[2] reading ", nrow(arms), " CellChat objects at LR resolution")
AXES <- c(Macrophage = "Mono_DC", HSC = "HSC_MPP")
extract_lr <- function(ds, s) {
  cc  <- readRDS(cc_path(ds, s))
  pr  <- cc@net$prob; pv <- cc@net$pval
  sig <- pr * (pv < PVAL)
  lv  <- levels(cc@idents)
  dimnames(sig) <- list(lv, lv, dimnames(pr)[[3]])
  tot <- sum(sig)
  ncell <- table(cc@idents)                    # the counts the graph was actually built on
  lrmeta <- as.data.table(cc@LR$LRsig)[, .(interaction_name, pathway_name,
                                           ligand, receptor,
                                           label = interaction_name_2)]
  out <- rbindlist(lapply(names(AXES), function(ax) {
    nd <- AXES[[ax]]
    if (!nd %in% lv) return(NULL)
    # sender and receiver are kept apart on purpose: "the macrophage compartment signals more" and
    # "the macrophage compartment is signalled to more" are different biology and averaging them
    # would cancel a genuine directional shift.
    snd <- apply(sig[nd, , , drop = FALSE], 3, sum) - sig[nd, nd, ]
    rcv <- apply(sig[, nd, , drop = FALSE], 3, sum) - sig[nd, nd, ]
    slf <- sig[nd, nd, ]
    rbind(data.table(axis = ax, direction = "sends",    interaction_name = names(snd), strength = as.numeric(snd)),
          data.table(axis = ax, direction = "receives", interaction_name = names(rcv), strength = as.numeric(rcv)),
          data.table(axis = ax, direction = "self",     interaction_name = names(slf), strength = as.numeric(slf)))[
          , n_node := as.integer(ncell[[nd]])][]
  }), fill = TRUE)
  if (!nrow(out)) return(NULL)
  out <- merge(out, lrmeta, by = "interaction_name", all.x = TRUE)
  out[, `:=`(dataset = ds, sample = s, sample_total = tot)][]
}
got <- lapply(seq_len(nrow(arms)), function(i) {
  message("    ", i, "/", nrow(arms), "  ", arms$sample[i])
  tryCatch(extract_lr(arms$dataset[i], arms$sample[i]),
           error = function(e) { message("    [warn] ", conditionMessage(e)); NULL })
})
lr <- rbindlist(Filter(Negate(is.null), got), fill = TRUE)
lr[arms, `:=`(pair_id = i.pair_id, arm = i.arm), on = c("dataset", "sample")]
lr[, rel := fifelse(sample_total > 0, strength / sample_total, NA_real_)]
message("    ", uniqueN(lr$interaction_name), " LR pairs x ", uniqueN(lr$sample), " samples")
# A sample can drop out here without erroring: if neither focus compartment is an identity in its
# CellChat object, extract_lr returns NULL. Report it rather than let the pair count quietly shrink.
miss <- setdiff(arms$sample, unique(lr$sample))
if (length(miss)) {
  message("    [note] no focus node in: ", paste(miss, collapse = ", "))
  for (s in miss) {
    p <- arms[sample == s, pair_id]
    message("           -> pair ", p, " loses its ", as.character(arms[sample == s, arm]), " side")
  }
}

## -- STRUCTURAL ZEROS ARE NOT LOW VALUES -------------------------------------------------------
# A sample whose focus compartment holds fewer than MIN_CELL cells gets ~zero communication
# probability for that node BY CONSTRUCTION -- there is no node, not a quiet one. Those rows have to
# come out before any test, exactly as 14 does at the edge level.
#
# This is not a hypothetical. The first run of this script left them in and reported the macrophage
# axis at 7/7 pairs up, p = 0.008, which would have been the headline result of the whole CCC
# section. It was an artefact: two of those seven "increases" were WT samples carrying 1 and 9
# Mono_DC cells, so their zero was structural and the mutant side won by default.
before <- lr[, uniqueN(paste(axis, pair_id))]
lr <- lr[n_node >= MIN_CELL]
# a PAIR survives only if both arms survive -- an unpaired arm contributes nothing to a paired test
okp <- lr[, .(k = uniqueN(arm)), by = .(axis, pair_id)][k == 2L, .(axis, pair_id)]
lr  <- merge(lr, okp, by = c("axis", "pair_id"))
message("    node >= ", MIN_CELL, " cells on both sides: ", nrow(okp), " of ", before,
        " axis-pair combinations kept")
NPAIR <- okp[, .(n_pairs_axis = .N), by = axis]
print(NPAIR)
if (!nrow(lr)) stop("no axis has a testable pair after the cell-count filter")

## ===== the paired test (identical machinery to 11 / 14) ========================================
# Exact one-sided signed-rank over the pairs. With 8 pairs the smallest attainable p is 1/256, so
# an LR pair has to move in the same direction in every pair to reach 0.004 -- and it still will not
# survive FDR across hundreds of LR pairs. That is a property of the cohort size, not a failure of
# the analysis, and the figures say so rather than quietly reporting nominal p.
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
message("[3] paired tests, LR level")
KEYS <- c("axis", "direction", "interaction_name", "label", "pathway_name", "ligand", "receptor")
lrt <- rbind(cbind(norm = "relative", test_block(lr, "rel",      KEYS)),
             cbind(norm = "absolute", test_block(lr, "strength", KEYS)))
# An LR pair present in only a few pairs cannot be tested; those rows are kept in the table with
# their n_pairs so the reader can see the coverage, but excluded from the figure.
lrt[, fdr := p.adjust(p_higher, "BH"), by = .(norm, axis, direction)]
setorder(lrt, norm, axis, direction, p_higher)
fwrite_safe(lrt, file.path(LCC_TAB_DIR, "20_ccc_lr_tests.csv"))

## ===== C2. the ligand-receptor answer ==========================================================
message("[4] C2 ligand-receptor")
# An LR pair has to be measurable in essentially every surviving pair of its axis. The floor of 4 is
# the paired test's own minimum; there is no point plotting an effect the test declined to compute.
c2 <- merge(lrt[norm == "relative" & direction != "self" & !is.na(median_delta)],
            NPAIR, by = "axis")
c2 <- c2[n_pairs >= pmax(4L, n_pairs_axis - 1L)]
c2 <- c2[median_delta != 0]
# top by absolute effect within each axis x direction, so a panel with small effects is not crowded
# out by one with large ones
c2[, ord := -abs(median_delta)]
setorder(c2, axis, direction, ord)
c2 <- c2[, head(.SD, 12L), by = .(axis, direction)]
c2[, panel := sprintf("%s compartment %s  (%d pairs)", axis, direction, n_pairs_axis)]
setorder(c2, -axis, direction)
c2[, panel := factor(panel, levels = unique(panel))]
c2[, up := median_delta > 0]
# CellChat's interaction_name_2 is "LIGAND - RECEPTOR", but ligand names contain hyphens of their
# own (HLA-DRA). Only the spaced separator may be replaced, or "HLA-DRA - CD4" becomes "HLA -> DRA - CD4".
c2[, disp := sub(" - ", " → ", label, fixed = TRUE)]
c2[, key := paste(panel, disp)]
setorder(c2, panel, median_delta)
c2[, key := factor(key, levels = key)]
c2[, pl := sprintf("%.0f/%d", frac_pairs_up * n_pairs, n_pairs)]
# the y key has to carry the panel name to keep four free scales from sharing a level set, but the
# panel name must not be printed on every tick -- so relabel back to the bare LR pair.
C2LAB <- setNames(c2$disp, as.character(c2$key))

fC2 <- ggplot(c2, aes(median_delta, key)) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = INK[["muted"]]) +
  geom_segment(aes(x = 0, xend = median_delta, yend = key, colour = up), linewidth = 1.1) +
  geom_point(aes(colour = up), size = 3.4, stroke = 0) +
  geom_text(aes(label = pl, hjust = ifelse(up, -0.45, 1.45)), size = 3.6, colour = INK[["muted"]]) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_y_discrete(labels = function(x) unname(C2LAB[x])) +
  scale_x_continuous(labels = percent_format(accuracy = 0.01), expand = expansion(mult = 0.22)) +
  scale_colour_manual(values = c(`TRUE` = PAL[["orange"]], `FALSE` = PAL[["blue"]]),
                      labels = c(`TRUE` = "higher in TP53-mut", `FALSE` = "higher in TP53-WT"),
                      name = NULL) +
  labs(title = "C2  Which ligand-receptor pairs move, macrophage and HSC axes",
       subtitle = "Median paired difference in signalling share. Text = pairs in which TP53-mut was higher.",
       x = "median paired difference in share of total communication (mut - WT)", y = NULL) +
  theme_lcc() + theme(panel.grid.major.y = element_blank())
save_fig(fC2, "C2_ccc_ligand_receptor", 15.5, 10)

## ===== C3. axis totals, bar with every pair drawn on it ========================================
message("[5] C3 axis totals")
tot <- lr[direction != "self", .(rel = sum(rel, na.rm = TRUE)),
          by = .(pair_id, arm, axis, direction)]
tt  <- cbind(norm = "relative", test_block(tot, "rel", c("axis", "direction")))
tot <- merge(tot, tt[, .(axis, direction, n_pairs, p_higher, frac_pairs_up)],
             by = c("axis", "direction"))
tot[, panel := sprintf("%s compartment %s\n%d pairs  |  %d higher in mut  |  p = %.3f",
                       axis, direction, n_pairs, round(frac_pairs_up * n_pairs), p_higher)]
setorder(tot, -axis, direction)
tot[, panel := factor(panel, levels = unique(panel))]
bar <- tot[, .(rel = median(rel)), by = .(arm, panel)]

# POWER RUNS THE OTHER WAY, and this has to be recorded because it decides which of the two
# directions may be believed. CellChat's permutation p depends on group size, so a smaller node
# yields fewer significant edges and therefore less measured signal. If the mutant nodes are the
# smaller ones, then "mut receives more" is observed AGAINST that bias, while "mut sends less" is
# exactly what the bias would manufacture on its own.
nsz <- unique(lr[, .(axis, pair_id, arm, n_node)])
nsw <- dcast(nsz, axis + pair_id ~ arm, value.var = "n_node")
nsw[, mut_smaller := get("TP53-mut") < get("TP53-WT")]
message("    node size: mutant side smaller in ",
        nsw[, sum(mut_smaller, na.rm = TRUE)], " of ", nrow(nsw), " axis-pairs")
print(nsw[])
fwrite_safe(merge(tot, nsw[, .(axis, pair_id, n_node_mut = get("TP53-mut"),
                               n_node_wt = get("TP53-WT"))], by = c("axis", "pair_id")),
            file.path(LCC_TAB_DIR, "20_ccc_axis_totals.csv"))
print(tt)

fC3 <- ggplot(bar, aes(arm, rel, fill = arm)) +
  geom_col(width = 0.62) +
  # the bar is a median; the lines are the pairs it is a median OF. With 8 pairs the bar alone is
  # not evidence, and crossing lines are what a null result looks like.
  geom_line(data = tot, aes(arm, rel, group = pair_id), inherit.aes = FALSE,
            colour = INK[["muted"]], linewidth = 0.5, alpha = 0.75) +
  geom_point(data = tot, aes(arm, rel), inherit.aes = FALSE,
             size = 2.1, stroke = 0, colour = INK[["primary"]], alpha = 0.85) +
  # ONE y SCALE ACROSS ALL FOUR PANELS. free_y let each panel rescale to its own range, which made a
  # 3-point difference and a 10-point difference look identical -- the opposite of the comparison
  # this figure exists to support.
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = ARMCOL, name = NULL) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C3  Total signalling on each axis",
       subtitle = "Bar = median of the pairs with the node on both sides; each line is one pair.",
       x = NULL, y = "share of the sample's total communication") +
  theme_lcc() + theme(panel.grid.major.x = element_blank(),
                      strip.text = element_text(size = 12.5, lineheight = 1.15))
save_fig(fC3, "C3_ccc_axis_totals", 15, 6)

message("[done] 20_ccc_macrophage")
