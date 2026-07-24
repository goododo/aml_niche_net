#!/usr/bin/env Rscript
# 01_nmf_per_sample.R ----
# Phase L1.1 -- per-sample NMF program discovery (GeneNMF, Gavish-style). Builds a list of per-sample
# Seurat objects subset to ONE compartment (malignant or normal cells, restricted to the 7 CCC
# hierarchy bins), then runs multiNMF across a range of k. The output (per-sample programs) is
# aggregated into cross-sample meta-programs by 02. Run ONCE per compartment.
#
# Compartment cells (join on `cell`, exact QC barcodes):
#   malignant = consensus malignant==1  AND  hierarchy_bin in CNMF_CCC_BINS
#   normal    = consensus malignant==0  AND  hierarchy_bin in CNMF_CCC_BINS
#   (Stromal/Unassigned excluded: inferCNV mis-calls non-hematopoietic cells as malignant.)
#
# INPUT : QC objects (QC_RDS_DIR); consensus per-cell (DIR_MALIGNANCY); projection per-cell (HIER_PROJ_DIR)
# OUTPUT: CNMF_RES_DIR/<compartment>/nmf_res.rds        (GeneNMF multiNMF result -> feeds 02)
#         CNMF_RES_DIR/<compartment>/nmf_samples.csv    (samples used + n_cells)
#
# Heavy: holds ~100 subset objects + runs |k| NMFs each. Submit as a job (m>=128G). Resume: skips if
# nmf_res.rds already exists (use --force to redo).
#   Rscript scripts/04_cnmf/01_nmf_per_sample.R --compartment malignant  [--limit N] [--force]

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(GeneNMF); library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_cnmf.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--compartment", type = "character", default = "malignant", help = "malignant | normal"),
  make_option("--limit",       type = "integer",   default = 0L, help = "first N qualifying samples (0=all; for testing)"),
  make_option("--force",       action = "store_true", default = FALSE, help = "recompute even if nmf_res.rds exists")
)))
stopifnot(opt$compartment %in% c("malignant", "normal"))
outdir <- file.path(CNMF_RES_DIR, opt$compartment); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
res_rds <- file.path(outdir, "nmf_res.rds")
if (file.exists(res_rds) && !isTRUE(opt$force)) { message("[skip] exists: ", res_rds, " (use --force)"); quit(save = "no") }

## ---- manifest: samples with QC + consensus + projection ----
qc <- data.table(rds = list.files(QC_RDS_DIR, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE))
qc <- qc[!grepl("/_", rds)][, sample := sub("\\.rds$", "", basename(rds))][, dataset := basename(dirname(rds))]
cons <- data.table(cons = list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE))
cons[, sample := sub("__consensus_percell\\.csv$", "", basename(cons))][, dataset := basename(dirname(cons))]
proj <- data.table(proj = list.files(HIER_PROJ_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
proj[, sample := sub("__bmm_percell\\.csv$", "", basename(proj))][, dataset := basename(dirname(proj))]
man <- Reduce(function(a, b) merge(a, b, by = c("dataset", "sample")), list(qc, cons, proj))
if (opt$limit > 0) man <- man[seq_len(min(opt$limit, .N))]
message(sprintf("[0] %d candidate samples; compartment=%s", nrow(man), opt$compartment))

## ---- build per-sample compartment objects (>= CNMF_MIN_CELLS cells) ----
obj.list <- list(); used <- list()
for (i in seq_len(nrow(man))) {
  ds <- man$dataset[i]; sid <- man$sample[i]
  seu <- readRDS(man$rds[i])
  md  <- data.table(cell = colnames(seu))
  md[fread(man$cons[i], select = c("cell", "malignant")), malignant := i.malignant, on = "cell"]
  md[fread(man$proj[i], select = c("cell", "hierarchy_bin")), hierarchy_bin := i.hierarchy_bin, on = "cell"]
  in_ccc <- md$hierarchy_bin %in% CNMF_CCC_BINS
  sel <- if (opt$compartment == "malignant") (md$malignant == 1 & in_ccc) else (md$malignant == 0 & in_ccc)
  sel[is.na(sel)] <- FALSE
  cells <- md$cell[sel]
  if (length(cells) < CNMF_MIN_CELLS) { rm(seu); next }
  sub <- subset(seu, cells = cells)
  sub <- NormalizeData(sub, verbose = FALSE)
  sub$sample_id <- sid; sub$dataset <- ds
  obj.list[[paste(ds, sid, sep = "..")]] <- sub
  used[[length(used) + 1]] <- data.table(dataset = ds, sample = sid, n_cells = length(cells))
  rm(seu, sub)
}
used_dt <- rbindlist(used, fill = TRUE)
stopifnot(nrow(used_dt) >= 2)
message(sprintf("[1] %d samples pass >= %d %s cells (median %d cells)",
                nrow(used_dt), CNMF_MIN_CELLS, opt$compartment, as.integer(median(used_dt$n_cells))))
fwrite_safe(used_dt, file.path(outdir, "nmf_samples.csv"))

## ---- HVG blocklist (keep programs off mito/ribo/hb axes) ----
allg <- unique(unlist(lapply(obj.list, rownames)))
blocklist <- grep(paste(CNMF_BLOCK_PATTERNS, collapse = "|"), allg, value = TRUE)
message(sprintf("[2] HVG blocklist: %d genes (mito/ribo/hb)", length(blocklist)))

## ---- per-sample NMF across k (this is the compute) ----
message(sprintf("[3] multiNMF: k=%s, nfeatures=%d, %d samples", paste(range(CNMF_K_RANGE), collapse = ":"),
                CNMF_NFEATURES, length(obj.list)))
nmf.res <- multiNMF(obj.list, assay = "RNA", slot = "data", k = CNMF_K_RANGE, nfeatures = CNMF_NFEATURES,
                    L1 = CNMF_L1, hvg.blocklist = blocklist, min.cells.per.sample = CNMF_MIN_CELLS, seed = SEED)
saveRDS(nmf.res, res_rds)
message(sprintf("[done] nmf_res (%d program sets) -> %s", length(nmf.res), res_rds))
cat("\nNEXT: aggregate into meta-programs -> 02_metaprogram_aggregate.R --compartment ", opt$compartment, "\n", sep = "")
