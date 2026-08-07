#!/usr/bin/env Rscript
# 01_bmm_project.R ----
# Phase 2.1 -- project every QC'd sample onto the BoneMarrowMap hematopoietic scaffold (Symphony +
# uwot), predict cell types (fine 55 + broad 24), roll broad up to the 8 hierarchy bins, and flag
# low-confidence projections. This is the SCAFFOLD for the CCC graph nodes; malignancy (from 02_
# malignancy consensus) is layered on per bin in a later step (R2: hierarchy stays orthogonal to
# malignant identity).
#
# Pipeline per sample:  map_Query -> calculate_MappingError -> predict_CellTypes(include_broad)
#
# INPUT  : QC objects under QC_RDS_DIR/<Dataset>/<sample>.rds ; BMM scaffold (BMM_SCAFFOLD_RDS + BMM_UWOT)
# OUTPUT : HIER_PROJ_DIR/<Dataset>/<sample>__bmm_percell.csv
#            cell, bmm_fine, bmm_prob, bmm_broad, hierarchy_bin, in_ccc_graph,
#            mapping_error, mapping_error_QC, high_error, author_anno, sample, dataset
#          HIER_PROJ_SUMMARY_CSV (one row per sample: bin composition + projection QC)
#
# QC flag (PRIMARY): high_error <- bmm_prob < MIN_PROJ_PROB (scale-free KNN agreement).
# mapping_error / mapping_error_QC from calculate_MappingError are kept as reference columns.
#
# Run (from PROJECT ROOT). Heavy-ish (loads the 362MB ref once); prefer a job for the full cohort:
#   Rscript scripts/03_hierarchy/01_bmm_project.R [--limit 1] [--datasets GSE227903,Chen2023]

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(BoneMarrowMap); library(data.table); library(optparse)
})
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--limit",    type = "integer",   default = 0L,  help = "process only the first N samples (0 = all; for testing)"),
  make_option("--datasets", type = "character", default = "",  help = "comma-separated dataset filter (default: all)")
)))

dir.create(dirname(HIER_PROJ_SUMMARY_CSV), recursive = TRUE, showWarnings = FALSE)

.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

# robust column resolver (BMM prob/broad column names vary slightly by version)
pick_col <- function(df, exact, pattern = NULL) {
  for (e in exact) if (e %in% names(df)) return(e)
  if (!is.null(pattern)) { hit <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE); if (length(hit)) return(hit[1]) }
  NA_character_
}

## ---- load scaffold + bin map ONCE ----
stopifnot(file.exists(BMM_SCAFFOLD_RDS), file.exists(BMM_UWOT))
message("[ref] loading BoneMarrowMap scaffold: ", basename(BMM_SCAFFOLD_RDS))
ref <- readRDS(BMM_SCAFFOLD_RDS)
ref$save_uwot_path <- BMM_UWOT           # BMM needs the on-disk uwot model for do_umap=TRUE
bin_map <- load_bin_map()                # CellType_Broad -> hierarchy_bin + in_ccc_graph
setkey(bin_map, CellType_Broad)

## ---- sample manifest (project ALL QC samples, healthy included -- they anchor the baseline) ----
# Roster from the QC report, not from ls -- see qc_rds_roster() in utils.R.
man <- qc_rds_roster(on_extra = "error")[, .(sample, dataset, rds)]
if (nzchar(opt$datasets)) man <- man[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
if (opt$limit > 0)        man <- man[seq_len(min(opt$limit, .N))]
stopifnot(nrow(man) > 0)
message(sprintf("[0] %d sample(s) across %d dataset(s) to project", nrow(man), uniqueN(man$dataset)))

project_one <- function(seu, sample_id, dataset, first) {
  cnt <- .get_counts(seu)
  md  <- tryCatch(seu@meta.data, error = function(e) as.data.frame(seu[[]]))

  q <- map_Query(exp_query = cnt, metadata_query = md, ref_obj = ref,
                 do_normalize = TRUE, do_umap = TRUE, verbose = FALSE)
  q <- calculate_MappingError(q, reference = ref, MAD_threshold = MAD_THRESHOLD_BMM)
  q <- predict_CellTypes(query_obj = q, ref_obj = ref, ref_label = BMM_REF_LABEL,
                         k = BMM_KNN_K, include_broad = TRUE)
  qm <- q@meta.data

  if (first) message("    [verify] projected metadata columns: ", paste(names(qm), collapse = ", "))

  c_fine  <- pick_col(qm, c("predicted_CellType", "final_CellType"), "^predicted_CellType$")
  c_broad <- pick_col(qm, c("predicted_CellType_Broad", "predicted_Broad"), "[Bb]road")
  c_prob  <- pick_col(qm, c("predicted_CellType_prob", "predicted_CellType_Broad_prob"), "predicted.*prob|_prob$")
  c_qc    <- pick_col(qm, c("mapping_error_QC", "mapping_QC"), "mapping.*QC")
  c_err   <- pick_col(qm, c("mapping_error", "mapping_error_score"), "mapping_error")
  c_auth  <- pick_col(md, c("anno_celltype", "anno_cellType", "CellType"), "anno.*cell")
  if (is.na(c_fine) || is.na(c_broad))
    stop("could not resolve BMM label columns; got: ", paste(names(qm), collapse = ", "))

  dt <- data.table(
    cell         = rownames(qm),
    bmm_fine     = as.character(qm[[c_fine]]),
    bmm_prob     = if (!is.na(c_prob)) as.numeric(qm[[c_prob]]) else NA_real_,
    bmm_broad    = as.character(qm[[c_broad]]),
    mapping_error    = if (!is.na(c_err)) as.numeric(qm[[c_err]]) else NA_real_,
    mapping_error_QC = if (!is.na(c_qc)) as.character(qm[[c_qc]]) else NA_character_,
    author_anno  = if (!is.na(c_auth)) as.character(md[[c_auth]][match(rownames(qm), rownames(md))]) else NA_character_,
    sample = sample_id, dataset = dataset)

  # roll broad -> 8-bin + in_ccc_graph; PRIMARY QC flag = prob-based
  dt[bin_map, `:=`(hierarchy_bin = i.hierarchy_bin, in_ccc_graph = i.in_ccc_graph), on = c(bmm_broad = "CellType_Broad")]
  dt[, high_error := as.integer(!is.na(bmm_prob) & bmm_prob < MIN_PROJ_PROB)]
  dt
}

all_summ <- vector("list", nrow(man))
for (i in seq_len(nrow(man))) {
  sid <- man$sample[i]; ds <- man$dataset[i]
  out <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  # RESUME GUARD WITH A FRESHNESS TEST. Existence alone made this unrunnable after
  # the v2 QC rebuild: all 143 projections for current samples date from 2026-07-16
  # and predate the QC objects they describe, so a re-run would have reported
  # "already projected" for every one and produced nothing. A projection is
  # reusable only if it POSTDATES the QC object it was computed from.
  # BMM_FORCE=1 recomputes regardless.
  if (file.exists(out) && !nzchar(Sys.getenv("BMM_FORCE")) &&
      file.mtime(out) >= file.mtime(man$rds[i])) {
    message(sprintf("[%d/%d] %s::%s already projected and current", i, nrow(man), ds, sid))
    all_summ[[i]] <- data.table(sample = sid, dataset = ds, status = "already_done"); next
  }
  message(sprintf("[%d/%d] %s::%s", i, nrow(man), ds, sid))
  res <- tryCatch({
    seu <- readRDS(man$rds[i])
    dt  <- project_one(seu, sid, ds, first = (i == 1))
    fwrite_safe(dt, out)
    comp <- dt[, .N, by = hierarchy_bin][order(-N)]
    data.table(sample = sid, dataset = ds, status = "ok", n_cells = nrow(dt),
               frac_high_error = round(mean(dt$high_error), 4),
               top_bin = comp$hierarchy_bin[1], top_bin_frac = round(comp$N[1] / nrow(dt), 3),
               n_stromal = dt[hierarchy_bin == "Stromal", .N])
  }, error = function(e) data.table(sample = sid, dataset = ds, status = paste0("error: ", conditionMessage(e))))
  all_summ[[i]] <- res
  rm(list = intersect(ls(), c("seu", "dt"))); gc(verbose = FALSE)
}

summary_dt <- rbindlist(all_summ, fill = TRUE)
fwrite_safe(summary_dt, HIER_PROJ_SUMMARY_CSV)
message("[done] projection summary -> ", HIER_PROJ_SUMMARY_CSV)
print(summary_dt[, .SD, .SDcols = intersect(c("dataset","sample","status","n_cells","frac_high_error","top_bin","top_bin_frac","n_stromal"), names(summary_dt))])

# ---------------------------------------------------------------------
# WHAT TO CHECK:
#   * frac_high_error: fraction of cells with bmm_prob < MIN_PROJ_PROB. High (>~0.3) -> that sample
#     projects poorly (batch/quality); inspect before trusting its bins.
#   * top_bin / composition: AML aspirates should be myeloid-skewed (HSC_MPP/LMPP_GMP/Mono_DC top);
#     healthy BM -> broad spread. Cross-check bmm bins vs author_anno for sanity.
#   * n_stromal: expected ~0 for aspirates (Stromal is B-layer, not a CCC node). Non-zero only in
#     niche-enriched samples (Chen2023 Niche_Immune).
# NEXT: per-bin malignant fraction = join this per-cell table with the 02_malignancy consensus
#   per-cell 'malignant' on `cell` (same QC barcodes -> exact join), aggregate by hierarchy_bin.
# ---------------------------------------------------------------------