#!/usr/bin/env Rscript
# 06_reconcile_annotation.R ----
# Join the two independent cell-typing results per cell and produce ONE corrected label, without
# discarding either input.
#
#   01_bmm_project  -> hierarchy_bin from the BoneMarrowMap projection   (bin_bmm)
#   05_marker_annotate -> hierarchy_bin from curated marker scoring      (bin_marker)
#   this script     -> hierarchy_bin, plus anno_source saying WHICH produced it
#
# WHY BOTH ARE KEPT: the projection is the only method that gives a fine BMM label and a mapping
# error; the markers are the only method that survives a dataset sitting off the reference
# manifold, and the only one that resolves stroma into subtypes (GSE253355: 12 niche subtypes vs a
# single "Stromal" label). Neither is redundant, and a downstream result that changes depending on
# which was used is something we want to be able to SEE, not something to hide behind one column.
#
# THE RULE, AND WHY IT CHANGED [2026-08-14].
#
# This script USED to let a confident marker call overwrite the projection's bin. That rule was
# measured against van Galen 2019's independent cell typing (08_validate_annotation.R) and it was
# making the annotation WORSE:
#
#   agreement with the referee   BMM 0.828 | marker 0.664 | old corrected output 0.772
#   where the two disagreed      BMM right 58.6%  marker right 13.8%
#   on the rows the old rule OVERRODE the projection: BMM right 62.6%, marker right 11.3%
#
# i.e. every override was wrong about 5.5x more often than right. The failure is mechanical, not
# noise: the panel's "DC: Migratory (LAMP3+CCR7+)" shares CCR7 with naive/central-memory T cells,
# and the method types CLUSTERS, so one bad cluster call mislabels hundreds of cells at once. The
# low-confidence flag separates good calls from bad (61.8% of the wrong calls are flagged vs 23.8%
# of the right ones) but does not catch enough of them to rescue the rule.
#
# SO: THE BIN ALWAYS COMES FROM THE PROJECTION. The markers keep the job they are actually good at
# -- naming a SUBTYPE inside that bin, which the projection cannot do (it collapses all stroma to a
# single "Stromal" label, where the markers resolve 17 subtypes including MSC/perivascular).
#
#   bin_bmm present, marker agrees      -> concordant       bin = bin_bmm, subtype KEPT
#   bin_bmm present, marker disagrees   -> bmm_over_marker  bin = bin_bmm, subtype DROPPED
#   bin_bmm present, no marker call     -> bmm_only         bin = bin_bmm
#   bin_bmm absent, marker present      -> marker_fill      bin = bin_marker  (see below)
#   neither                             -> unassigned
#
# A SUBTYPE IS ONLY KEPT WHEN ITS OWN BIN MATCHES THE COMPARTMENT. Without that gate the subtype
# column re-imports the exact error the bin rule was just removed for: 46,297 cells carry a
# stromal/vascular marker label while the projection puts them in a blood bin (16,064 in LMPP_GMP
# alone), and 499 cells van Galen types as T/NK are labelled "Stroma: MSC/perivascular". Gating on
# the bin drops those and leaves the 26,342 stromal cells whose subtype is consistent with where
# the projection put them.
#
# marker_fill (3.7% of cells) is the one place a marker call still sets the bin, and only because
# the projection returned NO bin there -- against the referee it is right 37.1% of the time, which
# is poor but strictly better than the alternative of "unassigned". It carries its own anno_source
# value so downstream can drop it in one filter.
#
# Cross-method agreement is also the ONLY projection-QC signal available for the 6 datasets that
# ship no healthy donor and therefore have no measurable FPR -- so the per-dataset agreement table
# this writes is a deliverable in its own right, not just a diagnostic.
#
# OUTPUT: <HIER_PROJ_DIR>/../03_annotation_reconciled/<dataset>/<sample>__anno_percell.csv
#         HIER_TAB_DIR/annotation_agreement_by_dataset.csv
#         HIER_TAB_DIR/annotation_agreement_by_bin_dataset.csv
#
#   Rscript scripts/03_hierarchy/06_reconcile_annotation.R [--dataset X]

suppressPackageStartupMessages({ library(optparse); library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = "")
)))

R <- qc_rds_roster(on_extra = "ignore")
if (nzchar(opt$dataset)) R <- R[dataset == opt$dataset]
stopifnot(nrow(R) > 0)
outroot <- ANNO_RECONCILED_DIR
binmap  <- unique(fread(BIN_MAP_TSV)[, .(hierarchy_bin, in_ccc_graph)])

one <- function(ds, sid) {
  pf <- file.path(HIER_PROJ_DIR,  ds, paste0(sid, "__bmm_percell.csv"))
  mf <- file.path(MARKER_ANNO_DIR, ds, paste0(sid, "__marker_percell.csv"))
  if (!file.exists(pf)) return(NULL)
  P <- fread(pf)
  keep_p <- intersect(c("cell", "bmm_fine", "bmm_prob", "bmm_broad", "mapping_error",
                        "high_error", "hierarchy_bin"), names(P))
  P <- P[, ..keep_p]
  setnames(P, "hierarchy_bin", "bin_bmm")

  if (file.exists(mf)) {
    M <- fread(mf, select = c("cell", "marker_cell_type", "marker_category",
                              "hierarchy_bin", "score", "margin_bin", "low_confidence", "cluster"))
    setnames(M, c("hierarchy_bin", "cluster", "low_confidence"),
                c("bin_marker", "marker_cluster", "marker_low_conf"))
    x <- merge(P, M, by = "cell", all.x = TRUE)
  } else {
    x <- copy(P)
    x[, `:=`(marker_cell_type = NA_character_, marker_category = NA_character_,
             bin_marker = NA_character_, score = NA_real_, margin_bin = NA_real_,
             marker_low_conf = NA, marker_cluster = NA_character_)]
  }

  # The projection leaves some cells with no bin at all (empty string, ~3.7% of cells). An empty
  # string is not a bin -- normalise both sides to NA first so "missing" and "disagrees" are never
  # confused, then apply the rule.
  blank <- function(v) is.na(v) | !nzchar(trimws(v))
  x[, `:=`(bin_bmm = fifelse(blank(bin_bmm), NA_character_, bin_bmm),
           bin_marker = fifelse(blank(bin_marker), NA_character_, bin_marker))]
  x[, agree := !is.na(bin_marker) & !is.na(bin_bmm) & bin_marker == bin_bmm]
  x[, anno_source := fifelse(is.na(bin_bmm) & is.na(bin_marker), "unassigned",
                      fifelse(is.na(bin_bmm), "marker_fill",
                       fifelse(is.na(bin_marker), "bmm_only",
                        fifelse(agree, "concordant", "bmm_over_marker"))))]
  # the bin is the projection's, except where the projection has none
  x[, hierarchy_bin := fifelse(anno_source == "marker_fill", bin_marker, bin_bmm)]
  x[anno_source == "unassigned", hierarchy_bin := "unassigned"]

  # SUBTYPE: kept only where the marker's own bin matches the compartment the cell was assigned to.
  # The raw call stays in marker_cell_type so nothing is destroyed and the drop rate is auditable;
  # cell_subtype is the column downstream should read.
  x[, cell_subtype := fifelse(!is.na(bin_marker) & bin_marker == hierarchy_bin &
                              !blank(marker_cell_type), marker_cell_type, NA_character_)]
  x[, subtype_dropped := !is.na(marker_cell_type) & !blank(marker_cell_type) & is.na(cell_subtype)]

  x <- merge(x, binmap, by = "hierarchy_bin", all.x = TRUE)
  # an unassigned cell is not a CCC node; NA there would propagate into every downstream filter
  x[is.na(in_ccc_graph), in_ccc_graph := FALSE]
  x[, `:=`(dataset = ds, sample = sid)]
  setcolorder(x, c("cell", "dataset", "sample", "hierarchy_bin", "in_ccc_graph",
                   "anno_source", "cell_subtype", "agree", "bin_bmm", "bin_marker",
                   "marker_cell_type", "subtype_dropped"))
  dir.create(file.path(outroot, ds), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, file.path(outroot, ds, paste0(sid, "__anno_percell.csv")))
  x[, .(dataset, sample, bin_bmm, bin_marker, hierarchy_bin, anno_source, agree,
        high_error, marker_low_conf, cell_subtype, subtype_dropped)]
}

all <- rbindlist(lapply(seq_len(nrow(R)), function(i) {
  r <- tryCatch(one(R$dataset[i], R$sample[i]),
                error = function(e) { message("[FAIL] ", R$sample[i], ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) message(sprintf("[ok] %s :: %s (%d cells)", R$dataset[i], R$sample[i], nrow(r)))
  r
}), fill = TRUE)
stopifnot(nrow(all) > 0)

cat("\n================ annotation source mix ================\n")
print(all[, .(cells = .N, pct = round(100 * .N / nrow(all), 1)), by = anno_source][order(-cells)])

cat("\n================ subtype coverage (what the markers still contribute) ================\n")
cat(sprintf("  cells with a retained cell_subtype : %d (%.1f%%)\n",
            all[!is.na(cell_subtype), .N], 100 * all[!is.na(cell_subtype), .N] / nrow(all)))
cat(sprintf("  marker calls dropped as inconsistent with the bin : %d\n", all[subtype_dropped == TRUE, .N]))
cat("\n  stromal subtypes recovered inside the projection's single 'Stromal' label:\n")
print(all[hierarchy_bin == "Stromal" & !is.na(cell_subtype), .N, by = cell_subtype][order(-N)][1:10])

cat("\n================ cross-method agreement by dataset ================\n")
byds <- all[!is.na(bin_marker), .(cells = .N, agree = round(mean(agree), 3),
                                  corrected = round(mean(anno_source == "marker_corrected"), 3),
                                  med_highErr = round(mean(high_error), 3)), by = dataset]
setorder(byds, agree)
print(byds)
fwrite(byds, file.path(HIER_TAB_DIR, "annotation_agreement_by_dataset.csv"))

cat("\n================ agreement by bin x dataset (BMM bin) ================\n")
bybin <- all[!is.na(bin_marker), .(cells = .N, agree = round(mean(agree), 3)),
             by = .(bin_bmm, dataset)]
fwrite(bybin, file.path(HIER_TAB_DIR, "annotation_agreement_by_bin_dataset.csv"))
print(dcast(bybin[cells >= 100], bin_bmm ~ dataset, value.var = "agree"))

cat("\n[reading it] Low agreement is NOT proof the projection is wrong -- neither method is ground\n")
cat("  truth. It marks where the two disagree, and for the datasets with healthy donors that is\n")
cat("  exactly where the measured FPR is worst (GSE253355 Mono_DC: FPR 0.80, agreement 0.03).\n")
cat("  For a dataset with no healthy donor this is the only projection-QC number there is.\n")
