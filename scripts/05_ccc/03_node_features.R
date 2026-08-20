# 03_node_features.R ----
# Per-(sample, hierarchy_bin) NODE FEATURES for the FGW feature term [locked: malignancy is a node
# feature, not a split]. One row per (eligible sample x 7-node vocab); absent bins -> n_cells=0, NA feats.
# Features: n_cells; frac_malignant (CNV consensus) and frac_malignant_vg (van Galen axis);
#           mean_cnv_burden; and stemness / CytoTRACE2 potency STRATIFIED by malignancy.
#
# WHY STRATIFY THE FEATURES AND NOT THE NODES. hierarchy_bin is a contaminated lineage label for
# malignant cells: measured against van Galen's own typing, 64.3% of the cells in the Erythroid bin
# are malignant MYELOID blasts, and T_NK and B_Plasma hold 299 and 112 more. So a pooled
# mean_stemness for "Erythroid" averages two different populations.
#
# The obvious fix -- split each node into normal/malignant, 7 nodes becoming 14 -- is WRONG here,
# for two measured reasons:
#   (1) FGW compares graph STRUCTURE and FGW_MASS_MODE is "ncells", so an absent node carries zero
#       mass. Healthy graphs would present a different node set from AML graphs, and the alignment
#       cost would be dominated by which nodes exist rather than by topology. HDS would separate
#       healthy from AML by construction -- the same defect FGW_ZERO_HEALTHY_MAL already has.
#   (2) The split would be made on a label that is INVERTED on its own negative controls: healthy
#       donors carry a median consensus malignant fraction of 0.113 against 0.037 for AML samples.
#       Every healthy sample would get a populated "malignant" node built from false positives.
# Stratifying the FEATURES keeps the 7-node vocabulary (so graphs stay comparable) while stopping
# the two populations from being averaged together. Nothing is lost: frac_malignant still carries
# the mixture, as a feature, which is what it was always meant to be.
#
# THE MALIGNANCY FEATURES ARE DELIBERATELY REDUNDANT. frac_malignant (CNV) is inverted on healthy
# donors; frac_malignant_vg (the van Galen axis, healthy median 0.049 vs AML 0.399, sample-level
# AUC 0.822) is not. Both are emitted, plus mean_cnv_burden as a CONTINUOUS covariate. Which of
# them earns a place in FGW_FEATURES is decided by 08_scoring/07_feature_decomposition.py, not here.
# RAW aggregates only -- feature scaling/normalization for FGW is deferred to 07.
#
# Decoupled from 02 (reads cell-level tables, not the LR tensor) and TOLERANT of an unfinished stemness
# run: mean_stemness is filled where the stemness_percell file exists, NA otherwise. Re-run with --force
# to backfill once the full-cohort stemness job completes. The run self-reports stemness coverage.
#
# CAVEATS carried from upstream:
#   - frac_malignant UNRELIABLE for T_NK/B_Plasma (inferCNV false positives) -> do not interpret there.
#   - healthy samples have no consensus file -> malignant=NA -> frac_malignant=NA (n_evaluable=0). The
#     healthy barycenter should set frac_malignant=0 explicitly downstream; 03 stays honest (NA).
#
# INPUT  : DIR_CCC/ccc_sample_manifest.csv ; CCC_BMM_DIR/<ds>/<s>__bmm_percell.csv ;
#          DIR_MALIGNANCY/<ds>/<s>__consensus_percell.csv ; CCC_BMM_DIR/<ds>/<s>__stemness_percell.csv (optional)
# OUTPUT : DIR_CCC/ccc_node_features.csv  (long: dataset,sample,timepoint,hierarchy_bin,has_graph,
#          n_cells,n_evaluable,n_malignant,frac_malignant,mean_stemness,stemness_sig,stemness_available)
# Usage  : Rscript scripts/05_ccc/03_node_features.R [--force]
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--force", action = "store_true", default = FALSE)
)))

out_feat <- file.path(DIR_CCC, "ccc_node_features.csv")
# FRESHNESS, not existence. This guard printed "[skip]" and exited 0 on a superseded cohort:
# results/tables/07_fgw/patient_scores.csv holds 148 rows of which 55 name samples that have
# left the cohort, and 47 current samples have never entered CCC at all. Re-running the chain
# hit five of these guards in a row and reported success.
# has_graph is read off the tensor glob (Step 3), so the tensors ARE an input. Leaving them out of
# .ins let this script run mid-array on 2026-08-20 and publish has_graph = TRUE for 78 of 138
# samples; re-running afterwards printed "[skip] ... is current" over the wrong answer.
.ins <- c(file.path(DIR_CCC, "ccc_sample_manifest.csv"),
          list.files(CCC_TENSOR_DIR, pattern = "__ccc_cellchat\\.csv$", recursive = TRUE, full.names = TRUE),
          list.files(CCC_BMM_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE),
          list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE),
          list.files(CCC_BMM_DIR, pattern = "__stemness_percell\\.csv$", recursive = TRUE, full.names = TRUE))
if (!is_stale(out_feat, .ins, force = opt$force)) {
  message("[skip] ", out_feat, " is current"); quit(status = 0)
}
if (file.exists(out_feat)) message("[recompute] ", stale_reason(out_feat, .ins, force = opt$force))

## -- Step 1. eligible samples (align node features to the graph cohort) ----
man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))[ccc_eligible == TRUE]
setorder(man, dataset, sample)
message("[1] eligible samples: ", nrow(man))

STEM_SIG <- CCC_STEMNESS_SIG
n_stem_ok <- 0L

## -- Step 2. per-sample per-bin aggregation ----
build_one <- function(ds, smp, tp) {
  # BINS FROM THE RECONCILED ANNOTATION, not the raw projection. This read CCC_BMM_DIR's
  # __bmm_percell.csv, which predates 06_reconcile_annotation.R: on GSE116256 alone 6.8% of cells
  # carry a different hierarchy_bin after reconciliation (the marker_fill cells the projection
  # abstained on). Node features built on the pre-reconciliation bin disagree with every other
  # stage that reports a bin.
  ann_f  <- file.path(ANNO_RECONCILED_DIR, ds, paste0(smp, "__anno_percell.csv"))
  con_f  <- file.path(DIR_MALIGNANCY,      ds, paste0(smp, "__consensus_percell.csv"))
  stem_f <- file.path(CCC_BMM_DIR,         ds, paste0(smp, "__stemness_percell.csv"))
  vg_f   <- file.path(VG_SCORE_DIR,        ds, paste0(smp, "__vgmalig_percell.csv"))
  ct_f   <- file.path(CYTOTRACE_DIR,       ds, paste0(smp, "__cytotrace2_percell.csv"))
  icnv_f <- file.path(INFERCNV_ROOT,       ds, smp, paste0(smp, "__infercnv_percell.csv"))
  if (!file.exists(ann_f)) { warning("no reconciled annotation for ", ds, "/", smp); return(NULL) }

  d <- fread(ann_f, select = c("cell", "hierarchy_bin", "in_ccc_graph", "high_error"))
  # fwrite writes NA as "" and fread returns "" -- normalise or every is.na() below is FALSE
  d[!nzchar(trimws(hierarchy_bin)), hierarchy_bin := NA_character_]
  d[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]
  d <- d[in_ccc_graph == TRUE & high_error == FALSE & hierarchy_bin %in% CCC_NODES]
  if (!nrow(d)) return(NULL)

  if (file.exists(con_f)) d <- fread(con_f, select = c("cell", "malignant"))[d, on = "cell"]
  else d[, malignant := NA_integer_]

  .join <- function(D, f, cols, newnames) {
    if (!file.exists(f)) { for (n in newnames) D[, (n) := NA_real_]; return(list(D, FALSE)) }
    hdr <- names(fread(f, nrows = 0L))
    if (!all(cols %in% hdr)) { for (n in newnames) D[, (n) := NA_real_]; return(list(D, FALSE)) }
    x <- fread(f, select = c("cell", cols)); setnames(x, cols, newnames)
    list(x[D, on = "cell"], TRUE)
  }
  r <- .join(d, stem_f, CCC_STEMNESS_SIG,   "stem_val");   d <- r[[1]]; stem_ok <- r[[2]]
  r <- .join(d, vg_f,   "vg_malignant",     "vg_val");     d <- r[[1]]; vg_ok   <- r[[2]]
  r <- .join(d, ct_f,   "CytoTRACE2_Score", "ct_val");     d <- r[[1]]; ct_ok   <- r[[2]]
  r <- .join(d, icnv_f, "score",            "burden_val"); d <- r[[1]]; bur_ok  <- r[[2]]

  # the van Galen axis is only defined where its signature was validated (HSC_MPP); scoring it
  # outside that bin would be a call in a compartment where it was never shown to work
  d[!(hierarchy_bin %in% VG_CALL_BINS), vg_val := NA_real_]

  .m <- function(v, keep) if (!any(keep & !is.na(v))) NA_real_ else mean(v[keep], na.rm = TRUE)
  agg <- d[, {
    isN <- !is.na(malignant) & malignant == 0L
    isM <- !is.na(malignant) & malignant == 1L
    .(n_cells       = .N,
      n_evaluable   = sum(!is.na(malignant)),
      n_malignant   = sum(malignant == 1, na.rm = TRUE),
      mean_stemness           = .m(stem_val, rep(TRUE, .N)),
      mean_stemness_normal    = .m(stem_val, isN),
      mean_stemness_malignant = .m(stem_val, isM),
      mean_cytotrace           = .m(ct_val, rep(TRUE, .N)),
      mean_cytotrace_normal    = .m(ct_val, isN),
      mean_cytotrace_malignant = .m(ct_val, isM),
      mean_cnv_burden = .m(burden_val, rep(TRUE, .N)),
      n_vg_scored     = sum(!is.na(vg_val)),
      n_vg_pos        = sum(vg_val > 0, na.rm = TRUE))
  }, by = .(hierarchy_bin = as.character(hierarchy_bin))]
  agg[, frac_malignant := fifelse(n_evaluable > 0L, n_malignant / n_evaluable, NA_real_)]

  # frac_malignant_vg needs the CALL, not the score; recompute it from the per-cell call column so
  # the threshold stays the dataset-matched one 81 chose rather than being re-derived here.
  if (vg_ok && file.exists(vg_f) && "vg_call" %in% names(fread(vg_f, nrows = 0L))) {
    vc <- fread(vg_f, select = c("cell", "vg_call"))
    dv <- merge(d[, .(cell, hierarchy_bin)], vc, by = "cell")
    av <- dv[!is.na(vg_call), .(n_vg_eval = .N, n_vg_mal = sum(vg_call == 1L)),
             by = .(hierarchy_bin = as.character(hierarchy_bin))]
    agg <- merge(agg, av, by = "hierarchy_bin", all.x = TRUE)
  } else agg[, `:=`(n_vg_eval = NA_integer_, n_vg_mal = NA_integer_)]
  agg[, frac_malignant_vg := fifelse(!is.na(n_vg_eval) & n_vg_eval > 0L, n_vg_mal / n_vg_eval, NA_real_)]

  full <- agg[data.table(hierarchy_bin = CCC_NODES), on = "hierarchy_bin"]
  for (cc in c("n_cells", "n_evaluable", "n_malignant", "n_vg_scored", "n_vg_pos"))
    full[is.na(get(cc)), (cc) := 0L]
  full[, `:=`(dataset = ds, sample = smp, timepoint = tp,
              stemness_sig = CCC_STEMNESS_SIG, stemness_available = stem_ok,
              vg_available = vg_ok, cytotrace_available = ct_ok, burden_available = bur_ok)]
  if (stem_ok) n_stem_ok <<- n_stem_ok + 1L
  full[]
}

feats <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  r <- man[i]
  tryCatch(build_one(r$dataset, r$sample, as.character(r$Timepoint)),
           error = function(e) { message("  [ERROR] ", r$dataset, "/", r$sample, " : ", conditionMessage(e)); NULL })
}), fill = TRUE)

## -- Step 3. flag which samples actually have a CCC graph (tensor exists; match on dataset+sample) ----
gfiles <- list.files(CCC_TENSOR_DIR, pattern = "__ccc_cellchat\\.csv$", recursive = TRUE)
graph_keys <- if (length(gfiles)) {
  data.table(dataset = basename(dirname(gfiles)),
             sample  = sub("__ccc_cellchat\\.csv$", "", basename(gfiles)))[, key := paste(dataset, sample)]$key
} else character(0)
feats[, has_graph := paste(dataset, sample) %in% graph_keys]

setcolorder(feats, c("dataset","sample","timepoint","hierarchy_bin","has_graph",
                     "n_cells","n_evaluable","n_malignant","frac_malignant",
                     "n_vg_eval","n_vg_mal","frac_malignant_vg","mean_cnv_burden",
                     "mean_stemness","mean_stemness_normal","mean_stemness_malignant",
                     "mean_cytotrace","mean_cytotrace_normal","mean_cytotrace_malignant",
                     "stemness_sig","stemness_available","vg_available",
                     "cytotrace_available","burden_available"))
setorder(feats, dataset, sample, hierarchy_bin)
fwrite_safe(feats, out_feat)
message("[3] wrote ", out_feat, "  (", nrow(feats), " rows = ", uniqueN(feats[, .(dataset,sample)]), " samples x ", length(CCC_NODES), " nodes)")

## -- Step 4. self-report: stemness coverage + treatment-axis malignant gradient ----
message("[4] stemness coverage (sig=", STEM_SIG, "): ", n_stem_ok, " / ", nrow(man), " samples have stemness")
if (n_stem_ok < nrow(man))
  message("    -> mean_stemness is partial; re-run with --force after the full-cohort stemness job finishes")

message("[4] frac_malignant by timepoint x bin (mean over samples with graph; primitive+blast bins):")
grad <- feats[has_graph == TRUE & hierarchy_bin %in% c("HSC_MPP","LMPP_GMP","Mono_DC") & n_evaluable > 0,
              .(mean_frac_mal = round(mean(frac_malignant, na.rm = TRUE), 3), n_samp = .N),
              by = .(timepoint, hierarchy_bin)]
print(dcast(grad, timepoint ~ hierarchy_bin, value.var = "mean_frac_mal"))

message("[4] mean_stemness (", STEM_SIG, ") by timepoint x bin (samples with stemness) -- ",
        "compare AML/Healthy separation vs frac_malignant above:")
sgrad <- feats[has_graph == TRUE & stemness_available == TRUE & !is.na(mean_stemness) &
               hierarchy_bin %in% c("HSC_MPP","LMPP_GMP","Mono_DC"),
               .(mean_stem = round(mean(mean_stemness, na.rm = TRUE), 3), n_samp = .N),
               by = .(timepoint, hierarchy_bin)]
print(dcast(sgrad, timepoint ~ hierarchy_bin, value.var = "mean_stem"))
message("[done]")