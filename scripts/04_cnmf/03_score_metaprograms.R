#!/usr/bin/env Rscript
# 03_score_metaprograms.R ----
# Phase L1.3 -- score every MALIGNANT cell for the malignant meta-programs (from 02), per sample.
# Design decision (locked): score ALL meta-programs (no pre-filtering) with AddModuleScore; carry the
# per-MP annotation (tumor_specific | normal_like) + robustness metrics alongside, so filtering is a
# DOWNSTREAM choice (nothing is discarded here). The per-cell MP usage is the functional-state layer
# that later augments the CCC graph ("which malignant program is a sender/receiver engaging").
#
# CELL SET: --cells malignant (default, the original behaviour) scores malignant cells within the 7
# CCC bins. --cells all_ccc_bins scores EVERY cell in those bins and records the malignant flag as a
# column instead of using it as a filter.
#
# WHY THE SECOND MODE EXISTS. The malignant-only output cannot serve as a CCC node feature, for four
# reasons that all dissolve at once when the filter is dropped:
#   (a) it is keyed to a malignant call that has since been rewritten -- of the cells the on-disk
#       mp_usage scored on 2026-07-21, only 56% are still called malignant by the 2026-08-18
#       consensus, and 25 of 60 checkable samples retain under half;
#   (b) 05_ccc/03 builds <feature>_normal and <feature>_malignant strata, and _normal would be
#       structurally all-NA -- which 08_scoring/07's constant-column guard does NOT catch, so it
#       would become an all-zero feature reporting a clean null;
#   (c) every healthy donor still receives MP activity, because all 23 healthy CCC-eligible samples
#       carry >= 20 inferCNV-called malignant cells (median 173, versus 171 for AML) -- so the
#       healthy values would be an average over known false positives;
#   (d) inferCNV under-calls ~9x, so "malignant" is a biased subsample of the compartment, not the
#       compartment.
# Scoring all cells and carrying the flag makes the stratification a downstream choice, which is the
# same principle this script already applies to MP filtering.
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
  make_option("--limit",    type = "integer",   default = 0L, help = "first N samples (0=all)"),
  make_option("--cells",    type = "character",  default = "malignant",
              help = "malignant (default) | all_ccc_bins -- see the header for why the second exists"),
  make_option("--force",    action = "store_true", default = FALSE,
              help = "rescore samples that already have an output file")
)))
stopifnot(opt$cells %in% c("malignant", "all_ccc_bins"))

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
# Roster from the QC report, not from ls -- see qc_rds_roster() in utils.R.
qc <- qc_rds_roster(on_extra = "error")[, .(dataset, sample, rds)]
cons <- data.table(cons = list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE))
cons[, sample := sub("__consensus_percell\\.csv$", "", basename(cons))][, dataset := basename(dirname(cons))]
proj <- data.table(proj = list.files(HIER_PROJ_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
proj[, sample := sub("__bmm_percell\\.csv$", "", basename(proj))][, dataset := basename(dirname(proj))]
man <- Reduce(function(a, b) merge(a, b, by = c("dataset", "sample")), list(qc, cons, proj))
if (nzchar(opt$datasets)) man <- man[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
if (opt$limit > 0) man <- man[seq_len(min(opt$limit, .N))]
message(sprintf("[0] %d samples", nrow(man)))

# The two cell sets write to DIFFERENT directories. Sharing one would make the file-exists skip below
# silently serve a malignant-only file to a caller that asked for all cells, and the two are not
# interchangeable -- that is the whole point of the switch.
usage_dir <- file.path(CNMF_RES_DIR, "malignant",
                       if (opt$cells == "malignant") "mp_usage" else "mp_usage_all_bins")
message("[0] cell set: ", opt$cells, "  ->  ", usage_dir)
class_lut <- setNames(ann$class, ann$malignant_MP)

for (i in seq_len(nrow(man))) {
  ds <- man$dataset[i]; sid <- man$sample[i]
  dst <- file.path(usage_dir, ds, paste0(sid, "__mp_usage.csv"))
  if (file.exists(dst) && !opt$force) { message(sprintf("[%d/%d] %s::%s done", i, nrow(man), ds, sid)); next }
  seu <- readRDS(man$rds[i])
  md  <- data.table(cell = colnames(seu))
  md[fread(man$cons[i], select = c("cell", "malignant")), malignant := i.malignant, on = "cell"]
  md[fread(man$proj[i], select = c("cell", "hierarchy_bin")), hierarchy_bin := i.hierarchy_bin, on = "cell"]
  in_bins <- md$hierarchy_bin %in% CNMF_CCC_BINS
  keep <- if (opt$cells == "malignant") md$cell[in_bins & md$malignant == 1] else md$cell[in_bins]
  keep <- keep[!is.na(keep)]
  if (length(keep) < 20) {
    message(sprintf("[%d/%d] %s::%s <20 cells in the %s set; skip", i, nrow(man), ds, sid, opt$cells))
    rm(seu); next
  }
  message(sprintf("[%d/%d] %s::%s  (%d cells, %s; %d malignant)", i, nrow(man), ds, sid,
                  length(keep), opt$cells, sum(md$malignant[match(keep, md$cell)] == 1, na.rm = TRUE)))

  sub <- subset(seu, cells = keep); sub <- NormalizeData(sub, verbose = FALSE)
  present <- rownames(.data_mat(sub))
  feats <- lapply(mp_genes, function(g) intersect(g, present))          # AddModuleScore needs present genes
  ok_mp <- names(feats)[sapply(feats, length) >= 3]
  sub <- AddModuleScore(sub, features = feats[ok_mp], name = "MPscore_", seed = SEED, ctrl = 50)
  sc_cols <- paste0("MPscore_", seq_along(ok_mp))

  # The malignant flag is CARRIED, not applied. Downstream (05_ccc/03) needs it to build the _normal
  # and _malignant strata; keeping it as a column is what lets the stratification be a downstream
  # choice rather than something baked in here against a call that keeps changing.
  dt <- data.table(cell = colnames(sub),
                   hierarchy_bin = md$hierarchy_bin[match(colnames(sub), md$cell)],
                   malignant     = md$malignant[match(colnames(sub), md$cell)])
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
