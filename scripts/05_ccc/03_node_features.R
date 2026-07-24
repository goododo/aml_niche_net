# 03_node_features.R ----
# Per-(sample, hierarchy_bin) NODE FEATURES for the FGW feature term [locked: malignancy is a node
# feature, not a split]. One row per (eligible sample x 7-node vocab); absent bins -> n_cells=0, NA feats.
# Features: frac_malignant (consensus), mean_stemness (CCC_STEMNESS_SIG from stemness_percell), n_cells.
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
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--force", action = "store_true", default = FALSE)
)))

out_feat <- file.path(DIR_CCC, "ccc_node_features.csv")
if (file.exists(out_feat) && !opt$force) { message("[skip] ", out_feat, " exists; --force to rebuild"); quit(status = 0) }

## -- Step 1. eligible samples (align node features to the graph cohort) ----
man <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))[ccc_eligible == TRUE]
setorder(man, dataset, sample)
message("[1] eligible samples: ", nrow(man))

STEM_SIG <- CCC_STEMNESS_SIG
n_stem_ok <- 0L

## -- Step 2. per-sample per-bin aggregation ----
build_one <- function(ds, smp, tp) {
  bmm_f  <- file.path(CCC_BMM_DIR,    ds, paste0(smp, "__bmm_percell.csv"))
  con_f  <- file.path(DIR_MALIGNANCY, ds, paste0(smp, "__consensus_percell.csv"))
  stem_f <- file.path(CCC_BMM_DIR,    ds, paste0(smp, "__stemness_percell.csv"))
  if (!file.exists(bmm_f)) { warning("no bmm for ", ds, "/", smp); return(NULL) }

  d <- fread(bmm_f, select = c("cell","hierarchy_bin","in_ccc_graph","high_error"))
  d[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]
  d <- d[in_ccc_graph == TRUE & high_error == FALSE & hierarchy_bin %in% CCC_NODES]

  # malignant (NA where no consensus, e.g. healthy)
  if (file.exists(con_f)) {
    cons <- fread(con_f, select = c("cell","malignant"))
    d <- cons[d, on = "cell"]
  } else d[, malignant := NA_integer_]

  # stemness signature (optional; NA where stemness run not yet done for this sample)
  stem_ok <- FALSE
  if (file.exists(stem_f)) {
    sh <- fread(stem_f, nrows = 0L)
    if (STEM_SIG %in% names(sh)) {
      stem <- fread(stem_f, select = c("cell", STEM_SIG))
      setnames(stem, STEM_SIG, "stem_val")
      d <- stem[d, on = "cell"]
      stem_ok <- TRUE
    }
  }
  if (!stem_ok) d[, stem_val := NA_real_]

  # aggregate per bin over QC-pass cells
  agg <- d[, .(n_cells       = .N,
               n_evaluable   = sum(!is.na(malignant)),
               n_malignant   = sum(malignant == 1, na.rm = TRUE),
               mean_stemness = if (all(is.na(stem_val))) NA_real_ else mean(stem_val, na.rm = TRUE)),
           by = .(hierarchy_bin = as.character(hierarchy_bin))]
  agg[, frac_malignant := fifelse(n_evaluable > 0L, n_malignant / n_evaluable, NA_real_)]

  # expand to the FIXED 7-node vocabulary (absent bins -> zero cells, NA feats)
  full <- data.table(hierarchy_bin = CCC_NODES)
  out  <- agg[full, on = "hierarchy_bin"]
  out[is.na(n_cells),     n_cells := 0L]
  out[is.na(n_evaluable), n_evaluable := 0L]
  out[is.na(n_malignant), n_malignant := 0L]
  out[, `:=`(dataset = ds, sample = smp, timepoint = tp,
             stemness_sig = STEM_SIG, stemness_available = stem_ok)]
  if (stem_ok) n_stem_ok <<- n_stem_ok + 1L
  out[]
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
                     "mean_stemness","stemness_sig","stemness_available"))
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