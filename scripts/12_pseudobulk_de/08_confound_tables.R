#!/usr/bin/env Rscript
# 07_confound_tables.R ----
# The three standing tables the pseudobulk pre-registration promised to report alongside the DE
# result and never produced. None of them is a gate: no hit is added or removed by this script. They
# exist because PREREG 7.2, 7.3 and 6.4 each commit to reporting a confounder, and a confounder that
# is named in the design but never quantified is a claim, not a control.
#
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv               (the 67-sample roster)
#          ANNO_RECONCILED_DIR/<ds>/<sample>__anno_percell.csv (per-cell bin + QC flags)
#          PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds            (full gene set, incl. chrY/XIST)
#          results/tables/05_ccc/ccc_node_features.csv          (mp_*/pt_* per sample x bin)
#          INFERCNV_GENE_ORDER                                  (chrY identification)
# OUTPUT : DIR_PSEUDOBULK/confound_mapping_error.csv   PREREG 7.2, per bin x arm
#          DIR_PSEUDOBULK/confound_composition.csv     PREREG 7.3, per bin, mp_*/pt_* arm shift
#          DIR_PSEUDOBULK/confound_inferred_sex.csv    PREREG 6.4, per sample + arm cross-tab
# Usage  : Rscript scripts/12_pseudobulk_de/07_confound_tables.R
#
#   PREREG 7.2  "reported as a standing table -- per bin, per arm: frac_high_error, pre-filter and
#               post-filter cell counts". The repo only carries frac_high_error per SAMPLE
#               (03_hierarchy/bmm_projection_summary.csv), so the per-BIN version is computed here
#               from the reconciled per-cell files. Direction matters and is pre-stated in 7.2: the
#               !high_error filter removes proportionally more AML cells, so each bin's surviving AML
#               cells are the ones most resembling the healthy reference, which biases logFC TOWARD
#               ZERO. A null is therefore conservative and a hit is not inflated by this mechanism.
#
#   PREREG 7.3  "the per-bin arm difference in those is reported as a DIAGNOSTIC alongside the DE
#               result, and any bin where it is large has its hits labelled composition-confounded.
#               It is not a gate -- the two cannot be separated with this design." Pseudobulk
#               averages over sub-states, so a composition shift inside a bin reads as differential
#               expression. mp_* (this project's own cNMF programs) and pt_* (BMM pseudotime) are the
#               only per-(sample,bin) composition proxies on disk.
#               NOTE pt_* is read but flagged, not used for the label: 30.7% of the BMM reference
#               sits at exactly 0 (a placeholder, not a pseudotime), so pt_ means are not comparable
#               across bins. Recorded so the column is not mistaken for usable.
#
#   PREREG 6.4  "Inferred sex per sample, read off those same genes in the pseudobulk, is reported as
#               a DIAGNOSTIC cross-tab by arm -- never as a covariate, because it is derived from the
#               data being tested." No sex column exists anywhere in this repo. chrY + XIST are
#               already excluded from every tested universe, so this table cannot change a hit; it
#               exists to say whether the arms happen to be sex-imbalanced, which a reviewer will ask.
#
# A DEGREE OF FREEDOM THIS SCRIPT INTRODUCES, DECLARED. PREREG 7.3 says "any bin where it is large"
# without defining large, and 6.4 does not define the sex call. Both are fixed here, in code, before
# the tables are read: composition-confounded := any mp_* program whose |arm difference| exceeds
# 1.0 SD of that program's pooled across-sample spread in that bin; sex call := chrY sum > 1 CPM AND
# XIST < 10x chrY sum -> male, the mirrored condition -> female, otherwise ambiguous. Neither
# threshold can move a hit, because neither table gates anything.

suppressPackageStartupMessages({
  library(data.table); library(here); library(Matrix)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))

MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))[in_analysis == TRUE]
GO  <- fread(INFERCNV_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
CHRY <- GO[chr %in% c("chrY", "Y"), unique(gene)]

## ------------------------------------------- PREREG 7.2 : mapping error per bin x arm ----
## Recomputed from the per-cell files because the repo only stores this per sample. Cells are counted
## BEFORE the !high_error filter (pre) and after it (post), within each bin, so the loss is per bin.
me <- rbindlist(lapply(seq_len(nrow(MAN)), function(i) {
  f <- file.path(ANNO_RECONCILED_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__anno_percell.csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = c("cell", "hierarchy_bin", "in_ccc_graph", "high_error", "bmm_broad"))
  d[!nzchar(trimws(hierarchy_bin)), hierarchy_bin := NA_character_]
  d[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]
  # in_ccc_graph and the bin vocabulary are upstream of the mapping-error question; hold them fixed
  # so the table isolates high_error alone. CCC_EXCLUDE_FINE is applied because PREREG 4.2 makes it
  # part of the cell set whose composition is being described.
  d <- d[in_ccc_graph == TRUE & hierarchy_bin %in% CCC_NODES & !(bmm_broad %in% CCC_EXCLUDE_FINE)]
  if (!nrow(d)) return(NULL)
  d[, .(dataset = MAN$dataset[i], sample = MAN$sample[i], split = MAN$split[i], arm = MAN$arm[i],
        n_pre = .N, n_high_error = sum(high_error), n_post = sum(!high_error)),
    by = .(hierarchy_bin = as.character(hierarchy_bin))]
}))
me[, frac_high_error := n_high_error / n_pre]
ME <- me[, .(n_samples = .N, cells_pre = sum(n_pre), cells_post = sum(n_post),
             frac_high_error_median = round(median(frac_high_error), 4),
             frac_high_error_cellwise = round(sum(n_high_error) / sum(n_pre), 4)),
         by = .(hierarchy_bin, split, arm)]
W <- dcast(ME[split == "Discovery"], hierarchy_bin ~ arm,
           value.var = c("frac_high_error_cellwise", "cells_pre", "cells_post"))
W[, aml_minus_healthy := round(frac_high_error_cellwise_AML - frac_high_error_cellwise_healthy, 4)]
W[, direction := fifelse(aml_minus_healthy > 0, "AML loses more (biases logFC toward 0)",
                  fifelse(aml_minus_healthy < 0, "healthy loses more (ANTI-conservative)", "equal"))]
fwrite(ME, file.path(DIR_PSEUDOBULK, "confound_mapping_error.csv"))
message("\n=== PREREG 7.2 mapping error, Discovery, cellwise frac_high_error ===")
print(W[, .(hierarchy_bin, AML = frac_high_error_cellwise_AML, healthy = frac_high_error_cellwise_healthy,
            aml_minus_healthy, direction)])

## -------------------------------------- PREREG 7.3 : within-bin composition shift ----
NF <- fread(file.path(TAB_DIR, "05_ccc", "ccc_node_features.csv"))
prog <- grep("^mp_", names(NF), value = TRUE)
ptc  <- grep("^pt_", names(NF), value = TRUE)
prog <- prog[!grepl("_normal$|_malignant$", prog)]      # all-cell stratum only; the split is unreliable
KEY <- MAN[split == "Discovery", .(dataset, sample, arm)]
X <- NF[KEY, on = c("dataset", "sample"), nomatch = 0][hierarchy_bin %in% CCC_NODES]
comp <- rbindlist(lapply(prog, function(p) {
  X[, {
    v <- get(p); ok <- is.finite(v)
    if (sum(ok) < 4L || uniqueN(arm[ok]) < 2L) .(NULL) else {
      sdv <- sd(v[ok])
      d <- mean(v[ok & arm == "AML"]) - mean(v[ok & arm == "healthy"])
      .(program = p, n = sum(ok), arm_diff = d, pooled_sd = sdv,
        abs_diff_in_sd = if (is.finite(sdv) && sdv > 0) abs(d) / sdv else NA_real_)
    }
  }, by = hierarchy_bin]
}), fill = TRUE)
comp[, exceeds_1sd := abs_diff_in_sd > 1.0]
COMP <- comp[, .(n_programs = .N, max_abs_diff_in_sd = round(max(abs_diff_in_sd, na.rm = TRUE), 3),
                 n_programs_over_1sd = sum(exceeds_1sd, na.rm = TRUE),
                 worst_program = program[which.max(abs_diff_in_sd)]), by = hierarchy_bin]
COMP[, composition_confounded := n_programs_over_1sd > 0]
COMP[, pt_note := sprintf("pt_* read (%d cols) but EXCLUDED from the label: 30.7%% of the BMM reference is a placeholder 0", length(ptc))]
fwrite(comp, file.path(DIR_PSEUDOBULK, "confound_composition.csv"))
message("\n=== PREREG 7.3 within-bin composition shift (mp_* cNMF programs, Discovery) ===")
print(COMP[, .(hierarchy_bin, n_programs, max_abs_diff_in_sd, n_programs_over_1sd,
               worst_program, composition_confounded)])

## ------------------------------------------- PREREG 6.4 : inferred sex diagnostic ----
sx <- rbindlist(lapply(seq_len(nrow(MAN)), function(i) {
  o <- readRDS(file.path(PB_RDS_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__pseudobulk.rds")))
  tot <- sum(o$counts); if (!tot) return(NULL)
  cy <- intersect(CHRY, o$genes); xi <- intersect("XIST", o$genes)
  ycpm <- if (length(cy)) sum(o$counts[cy, , drop = FALSE]) / tot * 1e6 else 0
  xcpm <- if (length(xi)) sum(o$counts[xi, , drop = FALSE]) / tot * 1e6 else 0
  data.table(dataset = MAN$dataset[i], sample = MAN$sample[i], split = MAN$split[i], arm = MAN$arm[i],
             chrY_cpm = round(ycpm, 2), XIST_cpm = round(xcpm, 2), n_chrY_genes_present = length(cy))
}))
# Declared above: male if chrY is expressed and XIST is not dominant; female if the mirror; else ambiguous.
sx[, inferred_sex := fifelse(chrY_cpm > 1 & XIST_cpm < 10 * chrY_cpm, "male",
                      fifelse(XIST_cpm > 1 & chrY_cpm <= 1, "female", "ambiguous"))]
fwrite(sx, file.path(DIR_PSEUDOBULK, "confound_inferred_sex.csv"))
message("\n=== PREREG 6.4 inferred sex cross-tab (diagnostic only; chrY+XIST already excluded) ===")
print(dcast(sx[, .N, by = .(split, arm, inferred_sex)], split + arm ~ inferred_sex, value.var = "N", fill = 0))
message("\nFisher test, Discovery male-vs-female by arm (ambiguous dropped):")
tb <- table(sx[split == "Discovery" & inferred_sex != "ambiguous", .(arm, inferred_sex)])
print(tb); if (all(dim(tb) == c(2, 2))) print(fisher.test(tb))
