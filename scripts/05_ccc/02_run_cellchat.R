# 02_run_cellchat.R ----
# Per-sample CellChat CCC inference on the FIXED 7-bin node vocabulary. Produces one standardized
# long tensor per sample (sender_bin x receiver_bin x LR pair, with prob/pval/node-sizes) = the
# 06_distance input. Grouping is hierarchy_bin ONLY (Main graph; malignancy is a node feature, not a
# split -> computed separately in 03_node_features.R). Engine-tagged so a liana arm (02b) can add a
# schema-compatible tensor for the R3 two-method gate.
#
# API pinned to installed CellChat 2.2.0.9001 (recon_ccc_probe [4]):
#   - QC objects are Seurat v5 with a counts-only layer -> NormalizeData first (creates 'data').
#   - computeCommunProb defaults distance.use/contact.dependent=TRUE are SPATIAL -> force distance.use=FALSE.
#   - projectData ABSENT -> raw.use=TRUE (no PPI smoothing).
#   - subsetCommunication(thresh=1) keeps ALL edges + pval; significance filtered downstream (R3).
#
# INPUT  : DIR_CCC/ccc_sample_manifest.csv (from 01) ; QC rds ; bmm + consensus per-cell CSVs
# OUTPUT : CCC_TENSOR_DIR/<ds>/<sample>__ccc_cellchat.csv   (long tensor, CCC_TENSOR_COLS)
#          CCC_GRAPH_DIR/<ds>/<sample>__cellchat.rds        (full CellChat object, persistent)
# Usage  :
#   Rscript scripts/05_ccc/02_run_cellchat.R --list                       # print eligible rows + count
#   Rscript scripts/05_ccc/02_run_cellchat.R --dataset=GSE227903 --sample=3904_R --nboot=20   # smoke test one
#   Rscript scripts/05_ccc/02_run_cellchat.R --row=$SLURM_ARRAY_TASK_ID   # array element (production)
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here); library(Seurat); library(future)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample_manifest", type = "character", default = file.path(DIR_CCC, "ccc_sample_manifest.csv")),
  make_option("--dataset", type = "character", default = NA_character_),
  make_option("--sample",  type = "character", default = NA_character_),
  make_option("--row",     type = "integer",   default = NA_integer_),
  make_option("--nboot",   type = "integer",   default = CCC_NBOOT),
  make_option("--workers", type = "integer",   default = CCC_WORKERS),
  make_option("--list",    action = "store_true", default = FALSE),
  make_option("--force",   action = "store_true", default = FALSE)
)))

## -- Step 0. eligible sample list (deterministic order = stable --row indexing) ----
stopifnot(file.exists(opt$sample_manifest))
man <- fread(opt$sample_manifest)[ccc_eligible == TRUE]
setorder(man, dataset, sample)
man[, row := .I]
if (opt$list) {
  print(man[, .(row, dataset, sample, Timepoint, total_cells, n_nodes_present)])
  message("[list] eligible samples: ", nrow(man))
  quit(status = 0)
}

# select which rows to run (braces so the if/else chain parses at top level)
sel <- if (!is.na(opt$row)) {
         man[row == opt$row]
       } else if (!is.na(opt$sample)) {
         man[sample == opt$sample & (is.na(opt$dataset) | dataset == opt$dataset)]
       } else if (!is.na(opt$dataset)) {
         man[dataset == opt$dataset]
       } else {
         man
       }
if (nrow(sel) == 0L) stop("no eligible sample matched the selection")
message("[0] will process ", nrow(sel), " sample(s); nboot=", opt$nboot, " workers=", opt$workers)

suppressPackageStartupMessages(library(CellChat))
if (opt$workers > 1L) {
  plan("multisession", workers = opt$workers)
  options(future.globals.maxSize = 4 * 1024^3)
}

## -- per-sample worker ----
run_one <- function(ds, smp, tp) {
  out_tensor <- file.path(CCC_TENSOR_DIR, ds, paste0(smp, "__ccc_cellchat.csv"))
  out_obj    <- file.path(CCC_GRAPH_DIR,  ds, paste0(smp, "__cellchat.rds"))
  # freshness, not existence: the tensor must postdate the QC object and the projection it is built from
  .ins <- c(file.path(CCC_QC_OBJ_DIR, ds, paste0(smp, ".rds")),
            file.path(CCC_BMM_DIR,    ds, paste0(smp, "__bmm_percell.csv")))
  if (!is_stale(out_tensor, .ins, force = opt$force)) {
    message("  [skip] ", ds, "/", smp, " current"); return(invisible())
  }
  if (file.exists(out_tensor)) message("  [recompute] ", ds, "/", smp, " -- ", stale_reason(out_tensor, .ins, force = opt$force))
  set.seed(SEED)

  ## load QC object + attach labels by exact barcode join
  obj <- readRDS(file.path(CCC_QC_OBJ_DIR, ds, paste0(smp, ".rds")))
  bmm <- fread(file.path(CCC_BMM_DIR, ds, paste0(smp, "__bmm_percell.csv")),
               select = c("cell","hierarchy_bin","in_ccc_graph","high_error"))
  bc  <- colnames(obj); j <- match(bc, bmm$cell)
  if (mean(!is.na(j)) < 0.999) warning("  BMM join < 0.999 for ", smp, " (", round(mean(!is.na(j)), 4), ")")
  obj$hierarchy_bin <- bmm$hierarchy_bin[j]
  obj$in_ccc_graph  <- as.logical(bmm$in_ccc_graph[j])
  obj$high_error    <- as.logical(bmm$high_error[j])

  ## keep only real CCC nodes, QC-pass cells (drops Stromal/Unassigned/high_error)
  keep <- obj$in_ccc_graph %in% TRUE & obj$high_error %in% FALSE & obj$hierarchy_bin %in% CCC_NODES
  obj  <- obj[, keep]
  gs   <- as.data.table(obj@meta.data)[, .(n = .N), by = .(bin = as.character(hierarchy_bin))]
  usable <- gs[n >= CCC_MIN_CELLS_PER_NODE]
  if (nrow(usable) < 2L) { message("  [skip] ", ds, "/", smp, " : < 2 usable nodes"); return(invisible()) }
  message("  [", ds, "/", smp, "] cells=", ncol(obj), " nodes>=",
          CCC_MIN_CELLS_PER_NODE, ": ", paste(usable$bin, collapse = ","))

  ## normalize (creates 'data' layer; LogNormalize is per-cell so subset-first is fine)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  data.input <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"),
                         error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))
  meta <- obj@meta.data
  meta$hierarchy_bin <- droplevels(factor(as.character(meta$hierarchy_bin), levels = CCC_NODES))
  meta$samples <- factor(smp)   # explicit single-sample label; suppresses CellChat v2 auto-add warning

  ## build CellChat from matrix + meta (sidesteps Seurat-v5 layered-object extraction quirks)
  cc <- createCellChat(object = data.input, meta = meta, group.by = "hierarchy_bin")

  ## DB = protein L-R only (manual interaction subset; robust vs subsetDB signature drift)
  db <- CellChatDB.human
  db$interaction <- db$interaction[db$interaction$annotation %in% CCC_DB_CATEGORIES, , drop = FALSE]
  cc@DB <- db

  ## preprocess -> communication probability (non-spatial mass-action)
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc,
          type            = CCC_COMMUN_TYPE,
          trim            = CCC_TRIM,
          population.size = CCC_POPULATION_SIZE,
          raw.use         = CCC_RAW_USE,
          distance.use    = CCC_DISTANCE_USE,   # FALSE: non-spatial. If this build errors on contact.*,
                                                # add contact.dependent = FALSE here (watch the test log).
          nboot           = opt$nboot,
          seed.use        = SEED)
  cc <- filterCommunication(cc, min.cells = CCC_MIN_CELLS_PER_NODE)

  ## extract standardized long tensor (ALL edges + pval)
  df <- as.data.table(subsetCommunication(cc, thresh = CCC_SUBSET_THRESH))
  if (nrow(df) == 0L) { message("  [warn] ", ds, "/", smp, " : no LR edges after inference"); }
  setnames(df, c("source","target"), c("sender_bin","receiver_bin"))
  df[, `:=`(sample = smp, dataset = ds, timepoint = tp, engine = "cellchat",
            n_sender   = gs$n[match(as.character(sender_bin),   gs$bin)],
            n_receiver = gs$n[match(as.character(receiver_bin), gs$bin)])]
  miss <- setdiff(CCC_TENSOR_COLS, names(df)); if (length(miss)) for (m in miss) df[, (m) := NA]
  tensor <- df[, ..CCC_TENSOR_COLS]

  fwrite_safe(tensor, out_tensor)
  dir.create(dirname(out_obj), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cc, out_obj)
  message("  [write] ", out_tensor, "  (", nrow(tensor), " edges)")
  invisible()
}

## -- drive ----
for (i in seq_len(nrow(sel))) {
  r <- sel[i]
  tryCatch(run_one(r$dataset, r$sample, as.character(r$Timepoint)),
           error = function(e) message("  [ERROR] ", r$dataset, "/", r$sample, " : ", conditionMessage(e)))
}
if (opt$workers > 1L) plan("sequential")
message("[done] processed ", nrow(sel), " sample(s)")