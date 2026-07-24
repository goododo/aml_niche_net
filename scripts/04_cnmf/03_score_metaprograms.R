#!/usr/bin/env Rscript
# 03_score_metaprograms.R ----
# Phase L1.3 -- score every MALIGNANT cell for the malignant meta-programs (from 02), per sample.
# Design decision (locked): score ALL meta-programs (no pre-filtering) with AddModuleScore; carry the
# per-MP annotation (tumor_specific | normal_like) + robustness metrics alongside, so filtering is a
# DOWNSTREAM choice (nothing is discarded here). The per-cell MP usage is the functional-state layer
# that later augments the CCC graph ("which malignant program is a sender/receiver engaging").
#
# Scoring is on malignant cells within the 7 CCC bins (same compartment as 01's malignant arm).
#
# INPUT : CNMF_TAB_DIR/malignant_metaprograms.tsv       (MP -> genes, from 02)
#         CNMF_TAB_DIR/malignant_mp_annotation.csv       (per MP: class + max_jaccard + coverage/sil)
#         QC objects (QC_RDS_DIR); consensus per-cell (DIR_MALIGNANCY); projection per-cell (HIER_PROJ_DIR)
# OUTPUT: CNMF_RES_DIR/malignant/mp_usage/<ds>/<sample>__mp_usage.csv
#           cell, sample, dataset, hierarchy_bin, <MP1..MPn scores>, top_MP, top_MP_score, top_MP_class
#         CNMF_TAB_DIR/mp_usage_summary.csv   (mean MP score per dataset, + dominant-MP composition)
#
#   Rscript scripts/04_cnmf/03_score_metaprograms.R [--datasets ...] [--limit N]

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", default = "", help = "comma-separated dataset filter"),
  make_option("--limit",    type = "integer",   default = 0L, help = "first N samples (0=all)")
)))

## ---- meta-program signatures + annotation (score ALL MPs; annotation travels with the output) ----
mp_tsv <- file.path(CNMF_TAB_DIR, "malignant_metaprograms.tsv")
mp_ann <- file.path(CNMF_TAB_DIR, "malignant_mp_annotation.csv")
stopifnot(file.exists(mp_tsv))
sig <- fread(mp_tsv)
mp_genes <- split(sig$gene, sig$MP)
mp_order <- unique(sig$MP)                                   # keep MP1..MPn order
ann <- if (file.exists(mp_ann)) fread(mp_ann) else data.table(malignant_MP = mp_order, class = NA_character_)
message(sprintf("[sig] %d malignant meta-programs: %s", length(mp_genes), paste(mp_order, collapse = ", ")))
message("[ann] classes: ", paste(sprintf("%s=%s", ann$malignant_MP, ann$class), collapse = ", "))

.data_mat <- function(seu) tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "data"),
                                    error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "data"))

## ---- manifest: samples with QC + consensus + projection ----
qc <- data.table(rds = list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE))
qc <- qc[!grepl("/_", rds)][, sample := sub("\\.rds$", "", basename(rds))][, dataset := basename(dirname(rds))]
cons <- data.table(cons = list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE))
cons[, sample := sub("__consensus_percell\\.csv$", "", basename(cons))][, dataset := basename(dirname(cons))]
proj <- data.table(proj = list.files(HIER_PROJ_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
proj[, sample := sub("__bmm_percell\\.csv$", "", basename(proj))][, dataset := basename(dirname(proj))]
man <- Reduce(function(a, b) merge(a, b, by = c("dataset", "sample")), list(qc, cons, proj))
if (nzchar(opt$datasets)) man <- man[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
if (opt$limit > 0) man <- man[seq_len(min(opt$limit, .N))]
message(sprintf("[0] %d samples", nrow(man)))

usage_dir <- file.path(CNMF_RES_DIR, "malignant", "mp_usage")
class_lut <- setNames(ann$class, ann$malignant_MP)

for (i in seq_len(nrow(man))) {
  ds <- man$dataset[i]; sid <- man$sample[i]
  dst <- file.path(usage_dir, ds, paste0(sid, "__mp_usage.csv"))
  if (file.exists(dst)) { message(sprintf("[%d/%d] %s::%s done", i, nrow(man), ds, sid)); next }
  seu <- readRDS(man$rds[i])
  md  <- data.table(cell = colnames(seu))
  md[fread(man$cons[i], select = c("cell", "malignant")), malignant := i.malignant, on = "cell"]
  md[fread(man$proj[i], select = c("cell", "hierarchy_bin")), hierarchy_bin := i.hierarchy_bin, on = "cell"]
  keep <- md$cell[md$malignant == 1 & md$hierarchy_bin %in% CNMF_CCC_BINS]
  keep <- keep[!is.na(keep)]
  if (length(keep) < 20) { message(sprintf("[%d/%d] %s::%s <20 malignant cells; skip", i, nrow(man), ds, sid)); rm(seu); next }
  message(sprintf("[%d/%d] %s::%s  (%d malignant cells)", i, nrow(man), ds, sid, length(keep)))

  sub <- subset(seu, cells = keep); sub <- NormalizeData(sub, verbose = FALSE)
  present <- rownames(.data_mat(sub))
  feats <- lapply(mp_genes, function(g) intersect(g, present))          # AddModuleScore needs present genes
  ok_mp <- names(feats)[sapply(feats, length) >= 3]
  sub <- AddModuleScore(sub, features = feats[ok_mp], name = "MPscore_", seed = SEED, ctrl = 50)
  sc_cols <- paste0("MPscore_", seq_along(ok_mp))

  dt <- data.table(cell = colnames(sub),
                   hierarchy_bin = md$hierarchy_bin[match(colnames(sub), md$cell)])
  for (j in seq_along(ok_mp)) dt[, (ok_mp[j]) := sub[[sc_cols[j]]][, 1]]
  # dominant MP per cell (argmax over scored MPs) + its class
  scm <- as.matrix(dt[, ..ok_mp])
  dt[, top_MP := ok_mp[max.col(scm, ties.method = "first")]]
  dt[, top_MP_score := round(scm[cbind(seq_len(.N), max.col(scm, ties.method = "first"))], 4)]
  dt[, top_MP_class := class_lut[top_MP]]
  dt[, `:=`(sample = sid, dataset = ds)]
  for (c in ok_mp) dt[, (c) := round(get(c), 4)]
  fwrite_safe(dt, dst)
  rm(seu, sub, dt); gc(verbose = FALSE)
}

## ---- summary: mean MP score + dominant-MP composition per dataset ----
uf <- list.files(usage_dir, pattern = "__mp_usage\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(uf)) {
  all <- rbindlist(lapply(uf, fread), fill = TRUE)
  mpcols <- intersect(mp_order, names(all))                 # MP score columns are named by MP (MP1..MPn)
  cat("\n================ mean MP score by dataset ================\n")
  print(all[, c(list(n_cells = .N), lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3))),
            by = dataset, .SDcols = mpcols])
  cat("\n================ dominant-MP composition (fraction of malignant cells) ================\n")
  comp <- all[, .N, by = .(dataset, top_MP, top_MP_class)][, frac := round(N / sum(N), 3), by = dataset]
  print(dcast(comp, dataset ~ top_MP, value.var = "frac", fill = 0))
  cat("\n(top_MP_class tags each program tumor_specific/normal_like; filter downstream as needed)\n")
  fwrite_safe(all[, c(list(n_cells = .N), lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3))),
                  by = dataset, .SDcols = mpcols], file.path(CNMF_TAB_DIR, "mp_usage_summary.csv"))
}
cat("\n[done] per-cell MP usage -> ", usage_dir, "/<ds>/<sample>__mp_usage.csv\n", sep = "")
