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
# is_healthy_sample(): the dataset-matched reference is built from this dataset's own healthy
# donors, so the healthy roster rule has to be the SAME one 96 calibrates against -- if the two
# drifted, the reference and the number judging it would describe different cohorts.
source(here::here("scripts", "config", "config_qc.R"))
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
  # Label conventions differ between the two BMM artefacts: the ANNOTATED dataset writes
  # "Plasma_Cell" while bmm_bin_map.tsv (and the Symphony projection output) write "Plasma Cell".
  # Matching raw dropped 150 plasma cells from the reference silently -- and plasma cells are
  # transcriptionally extreme (Ig-dominated), so leaving them out of the reference is exactly how
  # plasma cells in an observation sample come to read as aberrant. Normalise both sides.
  .norm_label <- function(x) trimws(gsub("[_[:space:]]+", " ", as.character(x)))
  bins <- bmap$hierarchy_bin[match(.norm_label(lab_kept), .norm_label(bmap$CellType_Broad))]
  if (sum(is.na(bins)) > 0)
    warning(sprintf("[ext-ref] %d reference cells have no bin-map entry -> EXCLUDED from grouping. Unmapped labels: %s",
                    sum(is.na(bins)),
                    paste(sort(unique(lab_kept[is.na(bins)])), collapse = ", ")))
  ext <- list(counts = cnt, cells = colnames(cnt), label = lab_kept, bin = bins)
  dir.create(dirname(INFERCNV_EXT_REF_CACHE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ext, INFERCNV_EXT_REF_CACHE)
  rm(bbm); gc(verbose = FALSE)
  message(sprintf("[ext-ref] %d external reference cells cached across %d lineages.",
                  ncol(cnt), length(unique(stats::na.omit(bins)))))
  ext
}

# Build (or load) this dataset's OWN healthy donors as a reference pool: same chemistry, depth and
# handling as the AML samples it will score, which the BMM atlas is not. Cached per dataset.
#
# Cells are pooled per hierarchy bin and the SOURCE DONOR is carried alongside, because the caller
# has to be able to drop one donor (leave-one-out) without rebuilding. We cache several times the
# per-bin target so that removing a donor still leaves enough to sample from.
get_matched_ref <- function(ds) {
  cache <- file.path(INFERCNV_MATCHED_REF_CACHE_DIR, paste0(ds, "__matched_ref.rds"))
  if (file.exists(cache)) {
    message(sprintf("[matched-ref] loading cached pool for %s", ds))
    return(readRDS(cache))
  }
  roster <- qc_rds_roster(datasets = ds, on_extra = "ignore")
  if (!nrow(roster)) return(NULL)
  healthy <- vapply(seq_len(nrow(roster)), function(i) {
    m <- tryCatch(readRDS(roster$rds[i])@meta.data, error = function(e) NULL)
    if (is.null(m)) return(FALSE)
    tp <- if ("Timepoint" %in% names(m)) as.character(m$Timepoint[1]) else NA_character_
    isTRUE(tp == "Healthy") || isTRUE(is_healthy_sample(roster$sample[i]))
  }, logical(1))
  roster <- roster[healthy]
  if (!nrow(roster)) { message(sprintf("[matched-ref] %s has no healthy donor -> BMM only", ds)); return(NULL) }
  message(sprintf("[matched-ref] building %s pool from %d healthy donor(s)", ds, nrow(roster)))

  cap <- 4L * max(INFERCNV_MATCHED_REF_PER_BIN)   # keep headroom for leave-one-out
  parts <- lapply(seq_len(nrow(roster)), function(i) {
    sid <- roster$sample[i]
    pf  <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
    if (!file.exists(pf)) return(NULL)
    prj <- fread(pf, select = c("cell", "hierarchy_bin", "high_error"))
    # only cells the projection is confident about: a mis-binned reference cell puts the wrong
    # lineage's expression into a lineage's reference interval, which is the exact artefact
    # INFERCNV_REF_PER_LINEAGE exists to prevent.
    prj <- prj[high_error == FALSE & !is.na(hierarchy_bin) & nzchar(hierarchy_bin) &
                 hierarchy_bin %in% names(INFERCNV_MATCHED_REF_PER_BIN)]
    if (!nrow(prj)) return(NULL)
    seu <- tryCatch(readRDS(roster$rds[i]), error = function(e) NULL); if (is.null(seu)) return(NULL)
    cnt <- .get_counts(seu)
    keep <- intersect(prj$cell, colnames(cnt)); if (!length(keep)) return(NULL)
    prj <- prj[match(keep, cell)]
    cnt <- cnt[, keep, drop = FALSE]
    colnames(cnt) <- paste0("MATCHREF_", sid, "_", colnames(cnt))
    rm(seu); gc(verbose = FALSE)
    list(counts = cnt, bin = prj$hierarchy_bin, donor = rep(sid, ncol(cnt)))
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)

  genes <- Reduce(intersect, lapply(parts, function(p) rownames(p$counts)))
  stopifnot(length(genes) > 1000)
  counts <- do.call(cbind, lapply(parts, function(p) p$counts[genes, , drop = FALSE]))
  bin    <- unlist(lapply(parts, `[[`, "bin"))
  donor  <- unlist(lapply(parts, `[[`, "donor"))

  # trim each bin to the cap so the cache stays small, keeping donors balanced
  set.seed(INFERCNV_MATCHED_REF_SEED)
  keep <- unlist(lapply(split(seq_along(bin), bin), function(ii)
    if (length(ii) <= cap) ii else sample(ii, cap)))
  keep <- sort(keep)
  out <- list(counts = counts[, keep, drop = FALSE], bin = bin[keep], donor = donor[keep],
              donors = unique(donor))
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache)
  message(sprintf("[matched-ref] %s: %d cells across %d bins from %d donors -> cached",
                  ds, ncol(out$counts), length(unique(out$bin)), length(out$donors)))
  out
}

# Choose, per bin, between this dataset's healthy donors and the BMM atlas. `exclude_donor` is the
# leave-one-out sample: when the sample being scored IS one of the healthy donors, its own cells
# must not sit in its own reference.
.matched_ref_lineage_block <- function(mref, drop_bins = character(0), exclude_donor = NA_character_) {
  if (is.null(mref)) return(NULL)
  ok <- !(mref$bin %in% drop_bins)
  if (!is.na(exclude_donor)) ok <- ok & mref$donor != exclude_donor
  idx_by_bin <- split(which(ok), mref$bin[ok])
  set.seed(INFERCNV_MATCHED_REF_SEED)
  sel <- unlist(lapply(names(idx_by_bin), function(b) {
    ii  <- idx_by_bin[[b]]
    tgt <- INFERCNV_MATCHED_REF_PER_BIN[[b]]
    # The floor is the BMM block size for this bin, not a flat minimum: swapping a 300-cell BMM
    # block for a 196-cell matched one would SHRINK the reference, and 98 measured that a short
    # block gives a biased-low P95 and a worse gate. Better a cross-dataset reference of the right
    # size than a matched one that is too small.
    if (length(ii) < max(tgt, INFERCNV_MATCHED_REF_MIN_PER_BIN)) return(integer(0))
    if (length(ii) <= tgt) ii else sample(ii, tgt)
  }))
  if (!length(sel)) return(NULL)
  list(idx = sel, groups = paste0("reference_matched__", mref$bin[sel]),
       bins = sort(unique(mref$bin[sel])))
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

# Assemble the foreign half of the reference: this dataset's own healthy donors where they can
# cover a bin, the BMM atlas everywhere else. Returns NULL if neither source yields a group.
# Both routes go through here so they cannot drift apart -- the drift this file exists to stop.
.assemble_foreign_ref <- function(s_cnt, ds, sid, drop_bins, ext_ref_fn, matched_ref_fn) {
  mblk <- NULL; mref <- NULL
  if (isTRUE(INFERCNV_MATCHED_REF)) {
    mref <- tryCatch(matched_ref_fn(ds), error = function(e) { message("  [matched-ref] ", conditionMessage(e)); NULL })
    mblk <- .matched_ref_lineage_block(mref, drop_bins = drop_bins, exclude_donor = sid)
  }
  covered <- if (is.null(mblk)) character(0) else mblk$bins
  ext  <- ext_ref_fn()
  eblk <- .ext_ref_lineage_block(ext, drop_bins = union(drop_bins, covered))

  cols <- list(); groups <- character(0); genes <- rownames(s_cnt)
  if (!is.null(mblk)) genes <- intersect(genes, rownames(mref$counts))
  if (length(eblk$idx)) genes <- intersect(genes, rownames(ext$counts))
  stopifnot(length(genes) > 1000)
  if (!is.null(mblk)) {
    cols[[length(cols) + 1L]] <- mref$counts[genes, mblk$idx, drop = FALSE]
    groups <- c(groups, mblk$groups)
  }
  if (length(eblk$idx)) {
    cols[[length(cols) + 1L]] <- ext$counts[genes, eblk$idx, drop = FALSE]
    groups <- c(groups, eblk$groups)
  }
  if (!length(cols)) return(NULL)
  message(sprintf("  [ref] foreign block: %d matched cells (%s) + %d BMM cells (%s)",
                  if (is.null(mblk)) 0L else length(mblk$idx),
                  if (is.null(mblk)) "-" else paste(mblk$bins, collapse = ","),
                  length(eblk$idx),
                  if (!length(eblk$idx)) "-" else paste(sort(unique(sub("^reference_external__", "", eblk$groups))), collapse = ",")))
  list(counts = do.call(cbind, cols), groups = groups, genes = genes)
}

# Assemble (counts, anno, ref_groups) for one sample. `ext_ref_fn` is a thunk so callers can
# cache the (large) external reference across samples.
build_infercnv_input <- function(s_cnt, route, ds, sid, ext_ref_fn = get_external_ref,
                                 matched_ref_fn = get_matched_ref) {
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
      fb <- .assemble_foreign_ref(s_cnt, ds, sid, drop_bins = c("T_NK", "B_Plasma"),
                                  ext_ref_fn = ext_ref_fn, matched_ref_fn = matched_ref_fn)
      if (!is.null(fb)) {
        counts <- cbind(s_cnt[fb$genes, , drop = FALSE], fb$counts)
        anno   <- data.table(cell = colnames(counts), group = c(grp, fb$groups))
        rg     <- c("reference_normal", sort(unique(fb$groups)))
        message(sprintf("  [ref] autologous lymphoid + %d foreign cells across %d lineage groups",
                        length(fb$groups), length(rg) - 1L))
        return(list(counts = counts, anno = anno, ref_groups = rg))
      }
      message("  [ref] no foreign lineage group cleared the per-group minimum",
              " -> autologous reference only (lineage artefact NOT controlled)")
    }
    return(list(counts = s_cnt,
                anno = data.table(cell = colnames(s_cnt), group = grp),
                ref_groups = "reference_normal"))
  }

  # external route -- no autologous normals at all, so the whole reference is foreign and the
  # dataset-matched donors matter MOST here (this is the frac_ref == 0 regime: 44 of 142 AML
  # samples, and the one configuration whose healthy-donor FPR was measured at 0.53-0.59).
  ext    <- ext_ref_fn()
  common <- intersect(rownames(s_cnt), rownames(ext$counts))
  stopifnot(length(common) > 1000)
  if (isTRUE(INFERCNV_REF_PER_LINEAGE)) {
    fb <- .assemble_foreign_ref(s_cnt, ds, sid, drop_bins = character(0),
                                ext_ref_fn = ext_ref_fn, matched_ref_fn = matched_ref_fn)
    stopifnot(!is.null(fb))
    counts <- cbind(s_cnt[fb$genes, , drop = FALSE], fb$counts)
    anno   <- data.table(cell  = colnames(counts),
                         group = c(rep("observation", ncol(s_cnt)), fb$groups))
    rg     <- sort(unique(fb$groups))
    message(sprintf("  [ref] foreign reference split into %d lineage groups (%d cells)",
                    length(rg), length(fb$groups)))
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
  # ROUTE ON THE DECISION, NOT ON HOW IT WAS REACHED [changed 2026-08-14].
  #
  # This used to read `method == "sorted_guard" -> "skip"`, which contradicted both the guard's own
  # documented meaning and the decision 20_refnorm_identify.R had already recorded. The guard in
  # config_malignancy.R says a sorted library "has no normal-cell population to use as an
  # AUTOLOGOUS inferCNV reference" -- that is the exact condition the external BMM reference exists
  # to cover, and 20 duly writes decision = "fallback_external" for those samples. Routing on
  # `method` then overrode that to "skip", so 33 samples (GSE185991 29, GSE147989 4) were dropped
  # from the cohort silently: 44_infercnv_run_one.R printed "[skip] sorted_guard" and exited 0, so
  # an array over them COMPLETED with no output and no error.
  #
  # 44 other samples carrying the SAME decision (fallback_external, reached via singleR finding too
  # few reference cells) ran through the external route without incident, so there was never a
  # methodological reason to treat these differently -- only the discarded `method` string.
  #
  # "skip" is kept as a legal route value: a sample with no usable reference of either kind still
  # has to be representable. Nothing produces it today.
  summ[, route := fifelse(decision == "autologous_ok", "autologous",
                 fifelse(decision == "fallback_external", "external", "skip"))]
  summ
}

message("[init] 00_infercnv_common loaded.")
