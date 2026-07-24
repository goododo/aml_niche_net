#!/usr/bin/env Rscript
# recon_node_feasibility.R ----  RECON ONLY, not part of 05_ccc pipeline
# Count cells per (dataset, sample, hierarchy_bin, malignant status) from the two
# existing per-cell tables, to ground the malignant/normal node-split decision.
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))

MAL_DIR <- DIR_MALIGNANCY                                           # FAST .../02_malignancy
BMM_DIR <- file.path(LARGE1_DIR, "02_seurat_objects/03_bmm_projected")

bmm_files <- list.files(BMM_DIR, pattern = "__bmm_percell\\.csv$",
                        recursive = TRUE, full.names = TRUE)
message("[1] found ", length(bmm_files), " bmm per-cell files")

res <- rbindlist(lapply(bmm_files, function(f) {
  bmm <- fread(f, select = c("cell","sample","dataset","hierarchy_bin",
                             "in_ccc_graph","high_error"))
  ds <- bmm$dataset[1]; smp <- bmm$sample[1]
  cf <- file.path(MAL_DIR, ds, paste0(smp, "__consensus_percell.csv"))
  if (file.exists(cf)) {
    con <- fread(cf, select = c("cell","malignant","confidence"))
    bmm <- merge(bmm, con, by = "cell", all.x = TRUE)
  } else {
    bmm[, `:=`(malignant = NA_integer_, confidence = NA_character_)]
  }
  bmm[, has_consensus := file.exists(cf)]
  bmm
}), fill = TRUE)

res[, in_ccc_graph := as.logical(in_ccc_graph)]                    # defensive coercion
res[, high_error   := as.logical(high_error)]
res[, mal_status := fifelse(is.na(malignant), "no_consensus",
                     fifelse(malignant == 1, "malignant", "normal"))]

counts <- res[in_ccc_graph == TRUE,
              .(n_cells = .N, n_high_error = sum(high_error, na.rm = TRUE)),
              by = .(dataset, sample, hierarchy_bin, mal_status)]

fwrite(counts, "recon_node_counts.csv")
message("[done] wrote recon_node_counts.csv (", nrow(counts), " rows)")