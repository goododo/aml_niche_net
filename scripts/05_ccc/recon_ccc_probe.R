#!/usr/bin/env Rscript
# recon_ccc_probe.R ----  RECON ONLY (not part of the 05_ccc pipeline)
# Two goals before writing 02_run_cellchat.R, per project discipline (read real API first):
#   (1) On ONE L2 sample: attach hierarchy_bin + malignant to the QC object, verify the barcode
#       join rate, and report post-QC per-bin node sizes.
#   (2) Introspect the INSTALLED CellChat (2.2.x) API + CellChatDB, so the real pipeline is written
#       against ground truth, not remembered signatures.
# Also echoes the resolved path constants + SEED (so we confirm config_ccc.R loads and SEED=491638L).
# INPUT  : QC rds + bmm/consensus per-cell CSVs (all existing)
# OUTPUT : recon_ccc_probe.txt   (upload this file back)
# Usage  : Rscript scripts/05_ccc/recon_ccc_probe.R [--dataset=GSE227903] [--sample=<auto>]
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(here); library(Seurat)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))   # defines CCC_QC_OBJ_DIR / CCC_BMM_DIR / etc.

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = "GSE227903"),
  make_option("--sample",  type = "character", default = NA_character_)
)))

SINK <- "recon_ccc_probe.txt"; con <- file(SINK, open = "wt")
report <- function(...) { cat(..., "\n", file = con, append = TRUE); cat(..., "\n") }

## -- Section 0. echo resolved paths + SEED (ground truth) ----
report("[0] SEED =", SEED, "(expect 491638)")
report("[0] CCC_QC_OBJ_DIR =", CCC_QC_OBJ_DIR)
report("[0] CCC_BMM_DIR    =", CCC_BMM_DIR)
report("[0] DIR_MALIGNANCY =", DIR_MALIGNANCY)
report("[0] DIR_CCC        =", DIR_CCC)

## -- Section 1. pick one L2 sample & inspect the QC object ----
qc_dir  <- file.path(CCC_QC_OBJ_DIR, opt$dataset)
rds_all <- list.files(qc_dir, pattern = "\\.rds$", full.names = TRUE)
report("\n[1] QC objects in", opt$dataset, ":", length(rds_all), "at", qc_dir)
if (length(rds_all) == 0L) { report("[1] FATAL: no rds found -> check CCC_QC_OBJ_DIR"); close(con); quit(status = 1) }

# auto-pick largest object (most cells -> cleanest CellChat demo) unless --sample given
pick <- if (is.na(opt$sample)) rds_all[which.max(file.info(rds_all)$size)] else
        file.path(qc_dir, paste0(opt$sample, ".rds"))
smp  <- sub("\\.rds$", "", basename(pick)); report("[1] picked sample:", smp)

obj <- readRDS(pick)
report("[1] dim(genes x cells):", paste(dim(obj), collapse = " x "))
report("[1] assays:", paste(Assays(obj), collapse = ","), "| default:", DefaultAssay(obj))
report("[1] layers:", paste(tryCatch(Layers(obj[[DefaultAssay(obj)]]),
                                      error = function(e) "v4-slot"), collapse = ","))
report("[1] meta cols:", paste(colnames(obj@meta.data), collapse = ","))
report("[1] barcode head:", paste(head(colnames(obj), 3), collapse = " | "))

## -- Section 2. attach labels by exact barcode join ----
bmm_f <- file.path(CCC_BMM_DIR,      opt$dataset, paste0(smp, "__bmm_percell.csv"))
con_f <- file.path(DIR_MALIGNANCY,   opt$dataset, paste0(smp, "__consensus_percell.csv"))
report("\n[2] bmm exists:", file.exists(bmm_f), "| consensus exists:", file.exists(con_f))

cells <- colnames(obj)
bmm <- fread(bmm_f, select = c("cell","hierarchy_bin","in_ccc_graph","high_error"))
report("[2] BMM join rate (obj cells found in bmm):", round(mean(cells %in% bmm$cell), 4))
md <- bmm[data.table(cell = cells), on = "cell"]                 # keep all obj cells
if (file.exists(con_f)) {
  cons <- fread(con_f, select = c("cell","malignant","confidence"))
  report("[2] consensus join rate:", round(mean(cells %in% cons$cell), 4))
  md <- cons[md, on = "cell"]
} else report("[2] no consensus for this sample (healthy/L1) -> malignant=NA")
md[, `:=`(in_ccc_graph = as.logical(in_ccc_graph), high_error = as.logical(high_error))]

report("\n[3] per-bin node sizes (in_ccc_graph & !high_error) -> what actually seeds CCC:")
report(paste(capture.output(print(
  md[in_ccc_graph == TRUE & high_error == FALSE, .N, by = hierarchy_bin][order(-N)]
)), collapse = "\n"))
if ("malignant" %in% names(md)) {
  report("[3] malignant x bin (post-QC):")
  report(paste(capture.output(print(
    md[in_ccc_graph == TRUE & high_error == FALSE,
       .N, by = .(hierarchy_bin, malignant)][order(hierarchy_bin, malignant)]
  )), collapse = "\n"))
}

## -- Section 4. introspect installed CellChat API + DB ----
suppressPackageStartupMessages(library(CellChat))
report("\n[4] CellChat version:", as.character(packageVersion("CellChat")))
db <- CellChatDB.human
report("[4] CellChatDB.human names:", paste(names(db), collapse = ","))
report("[4] interaction cols:", paste(colnames(db$interaction), collapse = ","))
if ("annotation" %in% colnames(db$interaction))
  report("[4] categories (annotation):\n", paste(capture.output(print(
    as.data.table(db$interaction)[, .N, by = annotation][order(-N)])), collapse = "\n"))
report("[4] n interactions total:", nrow(db$interaction))
for (fn in c("createCellChat","computeCommunProb","filterCommunication",
             "subsetCommunication","identifyOverExpressedGenes",
             "identifyOverExpressedInteractions","projectData","aggregateNet",
             "netAnalysis_computeCentrality")) {
  a <- tryCatch(paste(deparse(args(get(fn, envir = asNamespace("CellChat")))),
                      collapse = " "), error = function(e) "NOT FOUND")
  report("[4] args", fn, "::", a)
}

## -- Section 5. install-status check (safe: failed installs report FALSE, no crash) ----
report("\n[5] package load check (liana / nichenetr / igraph install status)")
for (p in c("liana","nichenetr","igraph","CellChat")) {
  ok <- requireNamespace(p, quietly = TRUE)
  v  <- if (ok) as.character(packageVersion(p)) else "-"
  report(sprintf("[5] %-10s loadable=%s version=%s", p, ok, v))
}

close(con); message("[done] wrote ", SINK, " -- please upload it")