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
# CORRECTION RULE (markers correct the projection, per the analysis plan):
#   agree                        -> concordant       (the bin both methods give)
#   disagree, marker confident   -> marker_corrected (marker wins)
#   disagree, marker low-conf    -> bmm_kept         (do not overwrite a good projection with a
#                                                     marker call the margin says is a coin flip)
#   marker missing               -> bmm_only
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

  # The projection leaves some cells with no bin at all (empty string, ~4.7% of cells). Those are
  # not a "good projection worth protecting" -- they are no answer, so bmm_kept must not apply to
  # them: a low-confidence marker call still beats a blank. Treat missing on either side first,
  # then the agree/disagree rule.
  blank <- function(v) is.na(v) | !nzchar(trimws(v))
  x[, `:=`(bin_bmm = fifelse(blank(bin_bmm), NA_character_, bin_bmm),
           bin_marker = fifelse(blank(bin_marker), NA_character_, bin_marker))]
  x[, agree := !is.na(bin_marker) & !is.na(bin_bmm) & bin_marker == bin_bmm]
  x[, anno_source := fifelse(is.na(bin_bmm) & is.na(bin_marker), "unassigned",
                      fifelse(is.na(bin_marker), "bmm_only",
                       fifelse(is.na(bin_bmm), "marker_only",
                        fifelse(agree, "concordant",
                         fifelse(marker_low_conf %in% TRUE, "bmm_kept", "marker_corrected")))))]
  x[, hierarchy_bin := fifelse(anno_source %in% c("marker_corrected", "marker_only"),
                               bin_marker, bin_bmm)]
  x[anno_source == "unassigned", hierarchy_bin := "unassigned"]
  x <- merge(x, binmap, by = "hierarchy_bin", all.x = TRUE)
  # an unassigned cell is not a CCC node; NA there would propagate into every downstream filter
  x[is.na(in_ccc_graph), in_ccc_graph := FALSE]
  x[, `:=`(dataset = ds, sample = sid)]
  setcolorder(x, c("cell", "dataset", "sample", "hierarchy_bin", "in_ccc_graph",
                   "anno_source", "agree", "bin_bmm", "bin_marker", "marker_cell_type"))
  dir.create(file.path(outroot, ds), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, file.path(outroot, ds, paste0(sid, "__anno_percell.csv")))
  x[, .(dataset, sample, bin_bmm, bin_marker, hierarchy_bin, anno_source, agree,
        high_error, marker_low_conf)]
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
