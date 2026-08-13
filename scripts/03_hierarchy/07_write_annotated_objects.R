#!/usr/bin/env Rscript
# 07_write_annotated_objects.R ----
# Write per-sample Seurat objects carrying BOTH cell-typing results in meta.data, to a NEW
# directory. The QC objects under 01_per_sample_qc are read-only here and stay untouched.
#
# WHY NOT WRITE INTO THE QC OBJECTS: 44_infercnv_run_one.R decides whether a sample needs
# re-running by comparing its burden CSV's mtime against the QC .rds mtime. Re-saving 214 QC
# objects would mark every one of the 179 finished inferCNV runs stale, and the next submission
# would silently recompute the entire cohort -- tens to hundreds of machine-hours for a metadata
# edit. A separate directory costs ~4 GB and no recompute.
#
# Columns added (all per cell, joined by cell id, never by position):
#   hierarchy_bin  in_ccc_graph  anno_source  agree      <- the corrected call and its provenance
#   bin_bmm  bmm_fine  bmm_prob  mapping_error  high_error <- projection side, unchanged
#   bin_marker  marker_cell_type  marker_category  margin_bin  marker_low_conf <- marker side
#
# Malignancy is deliberately NOT written here: it is still moving (the dataset-matched reference
# evaluation is in flight), and an object that mixes a settled annotation with an unsettled label
# invites someone to use the unsettled one.
#
#   Rscript scripts/03_hierarchy/07_write_annotated_objects.R [--dataset X] [--force]

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = ""),
  make_option("--force",   action = "store_true", default = FALSE)
)))

KEEP <- c("hierarchy_bin", "in_ccc_graph", "anno_source", "agree",
          "bin_bmm", "bmm_fine", "bmm_prob", "mapping_error", "high_error",
          "bin_marker", "marker_cell_type", "marker_category", "margin_bin", "marker_low_conf")

R <- qc_rds_roster(on_extra = "ignore")
if (nzchar(opt$dataset)) R <- R[dataset == opt$dataset]
stopifnot(nrow(R) > 0)
dir.create(ANNO_OBJ_DIR, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[0] %d sample(s) -> %s", nrow(R), ANNO_OBJ_DIR))

one <- function(ds, sid, rds) {
  out <- file.path(ANNO_OBJ_DIR, ds, paste0(sid, ".rds"))
  if (file.exists(out) && !opt$force) { message("[skip] ", sid); return(NULL) }
  af <- file.path(ANNO_RECONCILED_DIR, ds, paste0(sid, "__anno_percell.csv"))
  if (!file.exists(af)) { message("[skip] no reconciled table: ", sid); return(NULL) }
  A <- fread(af)
  seu <- readRDS(rds)

  # Join by cell id. Position-based assignment is the classic way to get a table that looks fine
  # and labels every cell as its neighbour, so require an exact set match before touching anything.
  if (!setequal(A$cell, colnames(seu)))
    stop(sprintf("cell sets differ (%s): obj=%d anno=%d shared=%d",
                 sid, ncol(seu), nrow(A), length(intersect(A$cell, colnames(seu)))))
  A <- A[match(colnames(seu), cell)]
  have <- intersect(KEEP, names(A))
  md <- as.data.frame(A[, ..have]); rownames(md) <- A$cell
  seu <- AddMetaData(seu, md)

  stopifnot(identical(rownames(seu@meta.data), colnames(seu)),
            !anyNA(seu$hierarchy_bin), all(seu$hierarchy_bin %in% c(HIERARCHY_BINS, "unassigned")))
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu, out)
  message(sprintf("[ok] %s :: %s  %d cells, %d annotation columns", ds, sid, ncol(seu), length(have)))
  data.table(dataset = ds, sample = sid, cells = ncol(seu), cols = length(have))
}

res <- rbindlist(lapply(seq_len(nrow(R)), function(i)
  tryCatch(one(R$dataset[i], R$sample[i], R$rds[i]),
           error = function(e) { message("[FAIL] ", R$sample[i], ": ", conditionMessage(e)); NULL })), fill = TRUE)

cat(sprintf("\n[done] wrote %d object(s), %d cells\n", nrow(res), if (nrow(res)) sum(res$cells) else 0L))
cat("  QC objects under 01_per_sample_qc were NOT modified.\n")
