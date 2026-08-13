#!/usr/bin/env Rscript
# 05_marker_annotate.R ----
# Marker-based cell typing, as an ALTERNATIVE to the BoneMarrowMap projection (01) for datasets
# where the projection cannot be trusted.
#
# WHY THIS EXISTS: 01_bmm_project maps every cell onto the BMM reference, and 96/the bin x dataset
# FPR table showed that GSE253355's cells sit far off that reference -- its median mapping_error is
# ~50% above every other dataset in EVERY bin (Mono_DC 14.3 vs 7.9-9.3, T_NK 11.9 vs 6.8-8.7), so
# its hierarchy_bin labels are not reliable. That dataset is also 80% of the cohort's stromal
# cells, i.e. 80% of the niche this project is about, so dropping it is far more costly than
# typing it a different way. Marker scoring does not depend on the projection at all.
#
# METHOD: cluster first, then type the CLUSTER. Per-cell marker scores in scRNA are dominated by
# dropout; a cluster mean is not. For each cluster we z-score every marker's mean expression
# ACROSS clusters, so a ubiquitously high gene cannot win on absolute level alone, then score each
# cell type as (mean z of its positive markers) - (mean z of its negative markers).
#
# INPUT : QC objects (raw counts) + the curated marker table
# OUTPUT: <MARKER_ANNO_DIR>/<dataset>/<sample>__marker_percell.csv
#           cell, cluster, marker_cell_type, marker_category, hierarchy_bin, in_ccc_graph,
#           score, margin, low_confidence
#         Column names match 01's per-cell table where they mean the same thing, so downstream
#         code can read either source.
#
#   Rscript scripts/03_hierarchy/05_marker_annotate.R --dataset GSE253355
#   Rscript scripts/03_hierarchy/05_marker_annotate.R --dataset GSE253355 --sample GSM8019213_H21

suppressPackageStartupMessages({ library(optparse); library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset",    type = "character"),
  make_option("--sample",     type = "character", default = ""),
  make_option("--resolution", type = "double",    default = MARKER_ANNO_RESOLUTION),
  make_option("--force",      action = "store_true", default = FALSE)
)))
stopifnot(nzchar(opt$dataset))

## ---- marker table -----------------------------------------------------------------------------
mk <- fread(MARKER_TABLE_CSV, colClasses = "character")
mk[, `:=`(marker = trimws(marker), role = trimws(role), priority = trimws(priority),
          modality = trimws(modality), cell_type = trimws(cell_type), category = trimws(category))]
# Pure-protein rows (CD105, CD90, CD73 ...) carry antibody names, not gene symbols; scoring them
# against an RNA matrix would silently drop them or, worse, hit an unrelated gene of the same name.
mk <- mk[grepl("RNA", modality, ignore.case = TRUE)]
mk <- mk[role %in% c("positive", "negative") & nzchar(marker)]
mk <- mk[!category %in% MARKER_ANNO_DROP_CATEGORIES]
mk <- unique(mk, by = c("cell_type", "marker", "role"))
message(sprintf("[0] marker table: %d RNA rows, %d cell types, %d categories",
                nrow(mk), uniqueN(mk$cell_type), uniqueN(mk$category)))

ct2bin <- MARKER_CATEGORY_TO_BIN            # category -> hierarchy_bin (config_hierarchy.R)
unmapped <- setdiff(unique(mk$category), names(ct2bin))
if (length(unmapped))
  stop("category with no hierarchy_bin mapping: ", paste(unmapped, collapse = ", "))

# Cell-type override first, category second. Both live in config_hierarchy.R; see the comment on
# MARKER_CELLTYPE_TO_BIN for why the category alone is not enough.
bin_of_marker_type <- function(cell_type, category) {
  b <- unname(MARKER_CELLTYPE_TO_BIN[cell_type])
  ifelse(is.na(b), unname(ct2bin[category]), b)
}
# A bin the scoring can never emit is a mapping typo, not a biological statement -- fail loudly.
stopifnot(all(MARKER_CELLTYPE_TO_BIN %in% HIERARCHY_BINS),
          all(MARKER_CATEGORY_TO_BIN %in% HIERARCHY_BINS))
# `mk` is already filtered to RNA rows, so a name missing here is EITHER a typo OR a type the
# marker table only defines by protein -- inert either way, but say which so nobody hunts a bug
# that is really just an antibody-only cell type.
bad_override <- setdiff(names(MARKER_CELLTYPE_TO_BIN), unique(mk$cell_type))
if (length(bad_override)) {
  raw_types <- unique(trimws(fread(MARKER_TABLE_CSV, select = "cell_type",
                                   colClasses = "character")$cell_type))
  typo <- setdiff(bad_override, raw_types)
  noRNA <- intersect(bad_override, raw_types)
  if (length(noRNA))
    message("[note] override targets with no scorable RNA marker (inert): ",
            paste(noRNA, collapse = ", "))
  if (length(typo))
    warning("MARKER_CELLTYPE_TO_BIN names not in the marker table at all (typo): ",
            paste(typo, collapse = ", "))
}
binmap <- fread(BIN_MAP_TSV)
bin_ccc <- unique(binmap[, .(hierarchy_bin, in_ccc_graph)])

## ---- samples ----------------------------------------------------------------------------------
roster <- qc_rds_roster(datasets = opt$dataset, on_extra = "ignore")
if (nzchar(opt$sample)) roster <- roster[sample == opt$sample]
stopifnot(nrow(roster) > 0)
outdir <- file.path(MARKER_ANNO_DIR, opt$dataset)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

annotate_one <- function(sid, rds) {
  out <- file.path(outdir, paste0(sid, "__marker_percell.csv"))
  if (file.exists(out) && !opt$force) { message("[skip] exists: ", sid); return(NULL) }
  message("\n[run] ", opt$dataset, " :: ", sid)
  seu <- readRDS(rds)
  if (ncol(seu) < MARKER_ANNO_MIN_CELLS) {
    message(sprintf("  [skip] only %d cells (< %d)", ncol(seu), MARKER_ANNO_MIN_CELLS)); return(NULL)
  }
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  npc <- min(30L, ncol(seu) - 1L, length(VariableFeatures(seu)) - 1L)
  seu <- RunPCA(seu, npcs = npc, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = seq_len(npc), verbose = FALSE)
  seu <- FindClusters(seu, resolution = opt$resolution, verbose = FALSE)
  cl  <- as.character(Idents(seu))
  message(sprintf("  [cluster] %d clusters over %d cells", uniqueN(cl), ncol(seu)))

  expr  <- SeuratObject::LayerData(seu, assay = "RNA", layer = "data")
  genes <- rownames(expr)
  present <- intersect(unique(mk$marker), genes)
  message(sprintf("  [markers] %d of %d marker genes present (%.0f%%)",
                  length(present), uniqueN(mk$marker), 100 * length(present) / uniqueN(mk$marker)))
  m <- mk[marker %in% present]

  # cluster mean expression, then z-score each gene ACROSS clusters
  cm <- vapply(split(seq_len(ncol(expr)), cl), function(ii)
    Matrix::rowMeans(expr[present, ii, drop = FALSE]), numeric(length(present)))
  cm <- matrix(cm, nrow = length(present), dimnames = list(present, names(split(seq_len(ncol(expr)), cl))))
  z  <- t(scale(t(cm)))
  z[!is.finite(z)] <- 0      # a gene with zero variance across clusters carries no information

  types <- unique(m$cell_type)
  sc <- vapply(types, function(tp) {
    pos <- m[cell_type == tp & role == "positive"]$marker
    neg <- m[cell_type == tp & role == "negative"]$marker
    if (length(pos) < MARKER_ANNO_MIN_MARKERS) return(rep(NA_real_, ncol(z)))
    p <- colMeans(z[pos, , drop = FALSE])
    n <- if (length(neg)) colMeans(z[neg, , drop = FALSE]) else 0
    p - MARKER_ANNO_NEG_WEIGHT * n
  }, numeric(ncol(z)))
  sc <- matrix(sc, nrow = ncol(z), dimnames = list(colnames(z), types))
  keep <- colSums(is.na(sc)) == 0
  message(sprintf("  [score] %d of %d cell types scorable (>= %d positive markers present)",
                  sum(keep), length(types), MARKER_ANNO_MIN_MARKERS))
  sc <- sc[, keep, drop = FALSE]

  best  <- colnames(sc)[max.col(sc, ties.method = "first")]
  ord   <- t(apply(sc, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  ann <- data.table(cluster = rownames(sc), marker_cell_type = best, score = round(ord[, 1], 4),
                    margin_celltype = round(ord[, 1] - ord[, 2], 4))
  ann[, marker_category := m$category[match(marker_cell_type, m$cell_type)]]
  ann[, hierarchy_bin   := bin_of_marker_type(marker_cell_type, marker_category)]

  # Confidence is judged at the BIN, not at the cell type. Eleven stroma subtypes score alike by
  # construction, so a cluster can be unambiguously Stromal while its top two SUBTYPE scores sit a
  # hair apart -- scoring confidence per cell type flagged 52% of cells as uncertain for a
  # distinction the hierarchy does not even make. What downstream reads is the bin.
  # MUST go through bin_of_marker_type, exactly as the label above does. Using the category map
  # here while the label uses the override map groups the score columns by a DIFFERENT partition
  # than the one being labelled -- e.g. every "Blast-like: *" type pooled under HSC_MPP -- so the
  # top-vs-second gap is measured between bins that are not the bins in the output. The result is
  # a plausible number and a wrong low_confidence flag. (Caught by 23_verify_marker_scoring 4.4.)
  col_bin  <- bin_of_marker_type(colnames(sc), m$category[match(colnames(sc), m$cell_type)])
  bin_best <- vapply(split(seq_along(col_bin), col_bin), function(jj)
    apply(sc[, jj, drop = FALSE], 1, max), numeric(nrow(sc)))
  bin_best <- matrix(bin_best, nrow = nrow(sc),
                     dimnames = list(rownames(sc), names(split(seq_along(col_bin), col_bin))))
  bord <- t(apply(bin_best, 1, function(r) sort(r, decreasing = TRUE)[1:2]))
  ann[, margin_bin := round(bord[match(cluster, rownames(bin_best)), 1] -
                            bord[match(cluster, rownames(bin_best)), 2], 4)]
  ann[, low_confidence := margin_bin < MARKER_ANNO_MIN_MARGIN]
  ann <- merge(ann, bin_ccc, by = "hierarchy_bin", all.x = TRUE)

  res <- data.table(cell = colnames(seu), cluster = cl)
  res <- merge(res, ann, by = "cluster", all.x = TRUE, sort = FALSE)
  res[, `:=`(sample = sid, dataset = opt$dataset, method = "marker")]
  setcolorder(res, c("cell", "cluster", "marker_cell_type", "marker_category",
                     "hierarchy_bin", "in_ccc_graph", "score", "margin_bin",
                     "margin_celltype", "low_confidence"))
  fwrite_safe(res, out)
  message(sprintf("  [done] %d cells -> %s", nrow(res), out))
  message("  [bins] ", paste(sprintf("%s=%d", names(table(res$hierarchy_bin)),
                                     as.integer(table(res$hierarchy_bin))), collapse = " "))
  if (any(res$low_confidence))
    message(sprintf("  [warn] %d cells (%.1f%%) in low-confidence clusters (margin < %.2f)",
                    sum(res$low_confidence), 100 * mean(res$low_confidence), MARKER_ANNO_MIN_MARGIN))
  invisible(res)
}

for (i in seq_len(nrow(roster))) {
  tryCatch(annotate_one(roster$sample[i], roster$rds[i]),
           error = function(e) message("[FAIL] ", roster$sample[i], ": ", conditionMessage(e)))
}
message("\n[all done] ", opt$dataset)
