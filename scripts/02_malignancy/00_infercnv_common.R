# =============================================================================
# 00_infercnv_common.R
# Shared inferCNV helpers for 40_infercnv_run.R (serial driver) and
# 44_infercnv_run_one.R (SLURM-array per-sample runner).
#
# WHY THIS FILE EXISTS: the two runners each carried their own copy of
# get_external_ref() and the reference-construction branch. They had already
# drifted (44 lacked 40's cache logging), and a change made in one silently
# left the other on the old behaviour -- which for the lineage-matched
# reference fix would mean the array runner, i.e. the one that actually
# produces the cohort, kept the broken single-group reference. Sourced by both;
# defines functions only, runs nothing.
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
# BIN_MAP_TSV + load_bin_map(): CellType_Broad -> hierarchy_bin, needed to split the external
# reference into one inferCNV group per lineage (INFERCNV_REF_PER_LINEAGE).
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(infercnv); library(data.table)
})

.get_counts <- function(seu)
  tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "counts"),
           error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "counts"))

# Build (or load) the external healthy reference: a stratified BMM subsample (raw counts),
# carrying its lineage label and hierarchy bin so it can be split into per-lineage groups.
get_external_ref <- function() {
  if (file.exists(INFERCNV_EXT_REF_CACHE)) {
    ext <- readRDS(INFERCNV_EXT_REF_CACHE)
    # Schema guard: caches written before lineage-matched references carry only counts/cells.
    # Rebuild rather than silently degrading to a single reference group -- that silent
    # degradation is exactly the failure this change exists to remove.
    if (all(c("label", "bin") %in% names(ext))) {
      message(sprintf("[ext-ref] loading cached external reference: %s", INFERCNV_EXT_REF_CACHE))
      return(ext)
    }
    message("[ext-ref] cached reference predates lineage grouping (no label/bin) -> rebuilding")
  }
  message("[ext-ref] building external reference from BoneMarrowMap subsample (one-time) ...")
  stopifnot(file.exists(INFERCNV_EXT_REF_RDS))
  bbm <- readRDS(INFERCNV_EXT_REF_RDS)
  lab <- as.character(bbm@meta.data[[INFERCNV_EXT_REF_LABELCOL]])
  set.seed(INFERCNV_EXT_REF_SEED)
  idx <- unlist(lapply(split(seq_along(lab), lab), function(ii)
    if (length(ii) <= INFERCNV_EXT_REF_PER_TYPE) ii else sample(ii, INFERCNV_EXT_REF_PER_TYPE)))
  cnt <- .get_counts(bbm)[, idx, drop = FALSE]
  colnames(cnt) <- paste0("EXTREF_", colnames(cnt))   # namespace to avoid id collisions
  lab_kept <- lab[idx]                                 # aligned to cnt columns via idx
  bmap <- load_bin_map()
  bins <- bmap$hierarchy_bin[match(lab_kept, bmap$CellType_Broad)]
  if (sum(is.na(bins)) > 0)
    message(sprintf("[ext-ref] %d reference cells have no bin-map entry -> excluded from grouping",
                    sum(is.na(bins))))
  ext <- list(counts = cnt, cells = colnames(cnt), label = lab_kept, bin = bins)
  dir.create(dirname(INFERCNV_EXT_REF_CACHE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ext, INFERCNV_EXT_REF_CACHE)
  rm(bbm); gc(verbose = FALSE)
  message(sprintf("[ext-ref] %d external reference cells cached across %d lineages.",
                  ncol(cnt), length(unique(stats::na.omit(bins)))))
  ext
}

# Pick external reference cells and label them one inferCNV group per lineage.
# See INFERCNV_REF_PER_LINEAGE in config_malignancy.R: with a single reference group infercnv's
# [min,max] per-group interval collapses to a point, so the lineage component of expression is
# never cancelled and every non-reference lineage reads as aberrant.
# `drop_bins` lets the autologous route skip the lineages its own normals already cover.
.ext_ref_lineage_block <- function(ext_ref, drop_bins = character(0)) {
  b    <- ext_ref$bin
  keep <- !is.na(b) & !(b %in% drop_bins)
  tab  <- table(b[keep])
  ok   <- names(tab)[tab >= INFERCNV_REF_MIN_PER_GROUP]
  keep <- keep & b %in% ok
  list(idx = which(keep), groups = paste0("reference_external__", b[keep]))
}

# Assemble (counts, anno, ref_groups) for one sample. `ext_ref_fn` is a thunk so callers can
# cache the (large) external reference across samples.
build_infercnv_input <- function(s_cnt, route, ds, sid, ext_ref_fn = get_external_ref) {
  if (route == "autologous") {
    ref_txt <- file.path(REFNORM_REF_CELL_DIR, ds, paste0(sid, "_ref_norm_cells.txt"))
    if (!file.exists(ref_txt)) stop("ref_norm cell list missing for autologous sample: ", ref_txt)
    ref_cells <- intersect(readLines(ref_txt), colnames(s_cnt))
    grp <- ifelse(colnames(s_cnt) %in% ref_cells, "reference_normal", "observation")

    # The autologous reference is lymphoid-only BY NECESSITY (REF_KEEP_LABELS): autologous
    # myeloid/progenitor cells may themselves be malignant. But a lymphoid-only reference is
    # exactly what makes every non-lymphoid lineage read as aberrant. So augment it with HEALTHY
    # EXTERNAL cells for the lineages the autologous normals cannot cover, each as its own
    # inferCNV group. The autologous lymphoid group is kept because it also controls for patient
    # and batch; the external groups only widen the per-gene reference interval.
    if (isTRUE(INFERCNV_REF_PER_LINEAGE)) {
      ext <- ext_ref_fn()
      blk <- .ext_ref_lineage_block(ext, drop_bins = c("T_NK", "B_Plasma"))
      if (length(blk$idx)) {
        common <- intersect(rownames(s_cnt), rownames(ext$counts))
        stopifnot(length(common) > 1000)
        counts <- cbind(s_cnt[common, , drop = FALSE], ext$counts[common, blk$idx, drop = FALSE])
        anno   <- data.table(cell = colnames(counts), group = c(grp, blk$groups))
        rg     <- c("reference_normal", sort(unique(blk$groups)))
        message(sprintf("  [ref] autologous lymphoid + %d external cells across %d lineages",
                        length(blk$idx), length(rg) - 1L))
        return(list(counts = counts, anno = anno, ref_groups = rg))
      }
      message("  [ref] no external lineage group cleared INFERCNV_REF_MIN_PER_GROUP",
              " -> autologous reference only (lineage artefact NOT controlled)")
    }
    return(list(counts = s_cnt,
                anno = data.table(cell = colnames(s_cnt), group = grp),
                ref_groups = "reference_normal"))
  }

  # external route
  ext    <- ext_ref_fn()
  common <- intersect(rownames(s_cnt), rownames(ext$counts))
  stopifnot(length(common) > 1000)
  if (isTRUE(INFERCNV_REF_PER_LINEAGE)) {
    blk <- .ext_ref_lineage_block(ext)
    stopifnot(length(blk$idx) > 0)
    counts <- cbind(s_cnt[common, , drop = FALSE], ext$counts[common, blk$idx, drop = FALSE])
    anno   <- data.table(cell  = colnames(counts),
                         group = c(rep("observation", ncol(s_cnt)), blk$groups))
    rg     <- sort(unique(blk$groups))
    message(sprintf("  [ref] external reference split into %d lineage groups (%d cells)",
                    length(rg), length(blk$idx)))
    return(list(counts = counts, anno = anno, ref_groups = rg))
  }
  counts <- cbind(s_cnt[common, , drop = FALSE], ext$counts[common, , drop = FALSE])
  list(counts = counts,
       anno = data.table(cell = colnames(counts),
                         group = c(rep("observation", ncol(s_cnt)),
                                   rep("reference_external", ncol(ext$counts)))),
       ref_groups = "reference_external")
}

# Run inferCNV for one sample given a counts matrix + annotation + ref group name(s).
# ref_subtract_use_mean_bounds is left at the package default (TRUE); it is what makes multiple
# reference groups behave as a [min,max] interval rather than a single averaged reference.
run_one <- function(counts, anno_dt, ref_groups, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  anno_file <- file.path(out_dir, "annotations.tsv")
  fwrite_safe(anno_dt, anno_file, sep = "\t", col.names = FALSE)
  obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file  = anno_file,
    gene_order_file   = INFERCNV_GENE_ORDER,
    ref_group_names   = ref_groups
  )
  infercnv::run(
    obj, cutoff = INFERCNV_CUTOFF, out_dir = out_dir,
    cluster_by_groups = TRUE, analysis_mode = INFERCNV_ANALYSIS,
    denoise = INFERCNV_DENOISE, HMM = INFERCNV_HMM,
    num_threads = INFERCNV_THREADS, no_plot = FALSE, save_rds = TRUE
  )
}

# Per-cell CNV burden from the final inferCNV residual matrix (centered ~1).
burden_from_obj <- function(obj) colMeans((obj@expr.data - 1)^2)

# Route derivation, shared so 40 and 44 cannot disagree about which sample takes which path.
infercnv_routes <- function() {
  stopifnot(file.exists(REFNORM_SUMMARY_CSV))
  summ <- fread(REFNORM_SUMMARY_CSV)
  summ[, route := fifelse(decision == "autologous_ok", "autologous",
                 fifelse(method == "sorted_guard", "skip", "external"))]
  summ
}

message("[init] 00_infercnv_common loaded.")
