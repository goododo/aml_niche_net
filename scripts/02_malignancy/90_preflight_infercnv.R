#!/usr/bin/env Rscript
# 90_preflight_infercnv.R ----
# PRE-FLIGHT for the inferCNV arm (40 / 44). Does NOT run inferCNV -- it validates the INPUTS so a
# multi-hour array doesn't fail (or silently produce garbage) on avoidable problems:
#   1. gene_order file present + sane (the 01 output)
#   2. gene-space overlap between each spot-checked sample and gene_order (the #1 silent failure:
#      if gene names don't match, inferCNV keeps almost no genes)
#   3. autologous ref_norm_cells.txt present + barcodes intersect the QC object; external samples
#      overlap the external-ref cache gene space
#   4. STALE-OUTPUT audit: how many samples already have a burden CSV (the 40/44 resume-guard will
#      SKIP these). Under the full seed refresh those are old-QC-cell / old-route results.
#
# Run (from PROJECT ROOT):
#   Rscript scripts/02_malignancy/90_preflight_infercnv.R [--n_spotcheck 3]

suppressPackageStartupMessages({ library(optparse); library(Seurat); library(SeuratObject); library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--n_spotcheck", type = "integer", default = 3L,
              help = "samples to load per route (autologous/external) for gene-overlap check [3]")
)))

.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

problems <- character(0); note <- function(x) { message("  ", x) }
flag     <- function(x) { problems <<- c(problems, x); message("  [!] ", x) }

## ---- 1. gene_order ----
message("[1] gene_order: ", INFERCNV_GENE_ORDER)
if (!file.exists(INFERCNV_GENE_ORDER)) {
  flag("gene_order MISSING -> run 01_infercnv_make_gene_order.R first")
  go_genes <- character(0)
} else {
  go <- fread(INFERCNV_GENE_ORDER, header = FALSE)
  go_genes <- go[[1]]
  note(sprintf("%d genes; chroms: %s", length(go_genes),
               paste(sort(unique(go[[2]])), collapse = ",")))
  if (anyDuplicated(go_genes)) flag(sprintf("gene_order has %d duplicate gene names", sum(duplicated(go_genes))))
  if (length(go_genes) < 15000) flag(sprintf("gene_order only %d genes (expected ~20-36k)", length(go_genes)))
}

## ---- 2. external ref cache ----
message("[2] external-ref cache: ", INFERCNV_EXT_REF_CACHE)
ext_genes <- character(0)
if (!file.exists(INFERCNV_EXT_REF_CACHE)) {
  flag("external-ref cache MISSING -> run 44 --prebuild_ext_ref (only needed if any sample routes external)")
} else {
  ext <- readRDS(INFERCNV_EXT_REF_CACHE)
  ext_genes <- rownames(ext$counts)
  note(sprintf("%d genes x %d cells", length(ext_genes), ncol(ext$counts)))
}

## ---- 3. routing summary ----
message("[3] ref_norm_summary: ", REFNORM_SUMMARY_CSV)
if (!file.exists(REFNORM_SUMMARY_CSV)) { flag("ref_norm_summary MISSING -> run 20_refnorm_identify.R first"); quit(status = 1) }
summ <- fread(REFNORM_SUMMARY_CSV)
summ[, route := fifelse(decision == "autologous_ok", "autologous",
              fifelse(method == "sorted_guard", "skip", "external"))]
note(sprintf("routes: autologous=%d, external=%d, skip=%d",
             sum(summ$route == "autologous"), sum(summ$route == "external"), sum(summ$route == "skip")))

## ---- 4. STALE-OUTPUT audit (resume-guard will skip these) ----
message("[4] stale-output / resume-skip audit")
summ[, burden_csv := file.path(INFERCNV_BURDEN_ROOT, dataset, paste0(sample_id, "_infercnv_burden.csv"))]
summ[, has_burden := file.exists(burden_csv)]
# Apply the SAME freshness invariant 44 uses, or this audit lies in the direction
# that matters: it reported 68 samples as "will be skipped" when every one of them
# predates its QC object and will in fact be recomputed. An existing output is only
# a skip if it POSTDATES the input it was computed from.
summ[, qc_rds := file.path(QC_RDS_DIR, dataset, paste0(sample_id, ".rds"))]
summ[, burden_current := has_burden & file.exists(qc_rds) &
       file.mtime(burden_csv) >= file.mtime(qc_rds)]
n_stale <- summ[route != "skip" & has_burden & !burden_current, .N]
n_skip  <- summ[route != "skip" & burden_current, .N]
if (n_stale > 0)
  note(sprintf("%d burden CSV(s) predate their QC object -> 44 recomputes them (not a skip)", n_stale))
if (n_skip > 0) {
  flag(sprintf("%d runnable sample(s) have a burden CSV NEWER than its QC object -> 44 will SKIP them.", n_skip))
  note("    That is correct resume behaviour if the QC object is final. To force a rebuild anyway,")
  note("    run the array with INFERCNV_FORCE=1.")
  print(summ[route != "skip" & burden_current, .(dataset, sample_id, route)][1:min(20, n_skip)])
} else note(sprintf("no CURRENT burden CSVs -> clean run over all %d runnable sample(s)",
                    summ[route != "skip", .N]))

## ---- 5. missing QC rds ----
message("[5] QC rds presence")
summ[, rds_in := file.path(QC_RDS_DIR, dataset, paste0(sample_id, ".rds"))]
miss_rds <- summ[route != "skip" & !file.exists(rds_in)]
if (nrow(miss_rds)) { flag(sprintf("%d runnable sample(s) missing QC rds:", nrow(miss_rds))); print(miss_rds[, .(dataset, sample_id)]) } else note("all runnable samples have a QC rds")

## ---- 6. spot-check gene / barcode overlap ----
message(sprintf("[6] spot-check gene/barcode overlap (%d per route)", opt$n_spotcheck))
spot <- rbind(
  head(summ[route == "autologous" & file.exists(rds_in)], opt$n_spotcheck),
  head(summ[route == "external"   & file.exists(rds_in)], opt$n_spotcheck))
for (i in seq_len(nrow(spot))) {
  sid <- spot$sample_id[i]; ds <- spot$dataset[i]; route <- spot$route[i]
  seu <- readRDS(spot$rds_in[i]); g <- rownames(.get_counts(seu)); bc <- colnames(seu)
  ov <- length(intersect(g, go_genes))
  msg <- sprintf("%s::%s [%s] genes=%d, gene_order-overlap=%d (%.0f%%)",
                 ds, sid, route, length(g), ov, 100 * ov / max(length(g), 1))
  if (ov < 8000) flag(paste0(msg, "  <- LOW overlap (gene-name mismatch?)")) else note(msg)

  if (route == "autologous") {
    ref_txt <- file.path(REFNORM_REF_CELL_DIR, ds, paste0(sid, "_ref_norm_cells.txt"))
    if (!file.exists(ref_txt)) { flag(sprintf("  %s: ref_norm_cells.txt MISSING (autologous route needs it)", sid))
    } else { rc <- readLines(ref_txt); ovb <- length(intersect(rc, bc))
      if (ovb < 20) flag(sprintf("  %s: only %d/%d ref cells intersect QC barcodes", sid, ovb, length(rc)))
      else note(sprintf("  %s: ref cells %d, intersect QC barcodes=%d", sid, length(rc), ovb)) }
  } else {
    ovx <- length(intersect(g, ext_genes))
    if (ovx < 1000) flag(sprintf("  %s: external-ref gene overlap only %d (<1000 -> 40/44 stopifnot fails)", sid, ovx))
    else note(sprintf("  %s: external-ref gene overlap=%d", sid, ovx))
  }
  rm(seu); gc(verbose = FALSE)
}

## ---- verdict ----
message("\n==================== VERDICT ====================")
if (!length(problems)) {
  message("GO: inferCNV inputs look clean. Safe to build the task list and sbatch 45.")
} else {
  message(sprintf("HOLD: %d issue(s) to resolve before the full run:", length(problems)))
  for (p in problems) message("  - ", p)
}
