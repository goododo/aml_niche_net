#!/usr/bin/env Rscript
# 24_verify_reconciliation.R ----
# BATCH 5 of the pipeline audit: the join that produces the corrected label (06).
#
# This is the table every downstream step is meant to read, so an error here is an error in
# everything. The specific hazards:
#   * in_ccc_graph following the PROJECTION bin instead of the corrected one -- a cell moved into
#     Stromal by the marker call would then still be flagged as a CCC-graph node, and the CCC graph
#     would silently gain cells the hierarchy says are not in it.
#   * anno_source and hierarchy_bin disagreeing: the column says "marker_corrected" while the bin
#     kept is the projection's. Both are legal values, so nothing errors.
#   * cells lost or duplicated in the merge -- a per-cell table with the wrong N still aggregates.
# Everything below is recomputed from the two SOURCE tables rather than trusted from the output.
#
#   Rscript scripts/99_admin/24_verify_reconciliation.R [n_samples]

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))

NS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(NS)) NS <- 20L
FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

binmap <- unique(fread(BIN_MAP_TSV)[, .(hierarchy_bin, in_ccc_graph)])
fs <- list.files(ANNO_RECONCILED_DIR, pattern = "__anno_percell.csv$", recursive = TRUE, full.names = TRUE)
if (!length(fs)) { cat("no reconciled tables yet -- run 06 first\n"); quit(save = "no", status = 1) }
set.seed(3); fsx <- sample(fs, min(NS, length(fs)))
cat(sprintf("\n=========== auditing %d of %d reconciled tables ===========\n", length(fsx), length(fs)))

n_lost <- 0L; n_dup <- 0L; n_ccc <- 0L; n_rule <- 0L; n_src <- 0L; n_agree <- 0L; n_illegal <- 0L
n_stale <- 0L; n_blank <- 0L
for (f in fsx) {
  X  <- fread(f)
  ds <- X$dataset[1]; sid <- X$sample[1]
  P  <- fread(file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv")),
              select = c("cell", "hierarchy_bin"))
  mf <- file.path(MARKER_ANNO_DIR, ds, paste0(sid, "__marker_percell.csv"))
  M  <- if (file.exists(mf)) fread(mf, select = c("cell", "hierarchy_bin", "low_confidence")) else NULL
  # a reconciled table older than its own inputs describes a previous generation
  if (!is.null(M) && file.mtime(f) < max(file.mtime(mf))) n_stale <- n_stale + 1L

  # (a) the join must neither lose nor duplicate a projected cell
  if (nrow(X) != nrow(P)) n_lost <- n_lost + 1L
  if (anyDuplicated(X$cell)) n_dup <- n_dup + 1L

  # (b) in_ccc_graph must follow the FINAL bin
  y <- merge(X[, .(cell, hierarchy_bin, in_ccc_graph)], binmap, by = "hierarchy_bin", all.x = TRUE)
  if (any(y$in_ccc_graph.x != y$in_ccc_graph.y, na.rm = TRUE)) n_ccc <- n_ccc + 1L
  if (any(!X$hierarchy_bin %in% c(HIERARCHY_BINS, "unassigned"))) n_illegal <- n_illegal + 1L
  if (any(!nzchar(trimws(X$hierarchy_bin)))) n_blank <- n_blank + 1L

  # (c) recompute the whole rule from the two sources
  if (!is.null(M)) {
    setnames(P, "hierarchy_bin", "bb"); setnames(M, "hierarchy_bin", "bm")
    R <- merge(P, M, by = "cell", all.x = TRUE)
    bl <- function(v) is.na(v) | !nzchar(trimws(v))
    R[, `:=`(bb = fifelse(bl(bb), NA_character_, bb), bm = fifelse(bl(bm), NA_character_, bm))]
    R[, ag := !is.na(bm) & !is.na(bb) & bm == bb]
    R[, src := fifelse(is.na(bb) & is.na(bm), "unassigned",
                fifelse(is.na(bm), "bmm_only",
                 fifelse(is.na(bb), "marker_only",
                  fifelse(ag, "concordant",
                   fifelse(low_confidence %in% TRUE, "bmm_kept", "marker_corrected")))))]
    R[, fin := fifelse(src %in% c("marker_corrected","marker_only"), bm, bb)]
    R[src == "unassigned", fin := "unassigned"]
    Z <- merge(X[, .(cell, anno_source, agree, hierarchy_bin, bin_bmm, bin_marker)], R, by = "cell")
    if (any(Z$anno_source != Z$src))          n_src   <- n_src + 1L
    if (any(Z$hierarchy_bin != Z$fin))        n_rule  <- n_rule + 1L
    if (any(Z$agree != Z$ag))                 n_agree <- n_agree + 1L
    setnames(P, "bb", "hierarchy_bin"); setnames(M, "bm", "hierarchy_bin")
  }
}

cat("\n=========== 5.1 join integrity ===========\n")
chk(n_lost == 0, "reconciled rows == projected rows", sprintf("%d samples differ", n_lost))
chk(n_dup  == 0, "no duplicated cell", sprintf("%d samples", n_dup))
chk(n_stale == 0, "no reconciled table predates its marker input", sprintf("%d stale", n_stale))

cat("\n=========== 5.2 in_ccc_graph follows the CORRECTED bin ===========\n")
chk(n_ccc == 0, "in_ccc_graph matches the bin map for the final bin", sprintf("%d samples", n_ccc))
chk(n_illegal == 0, "every final bin is a legal hierarchy bin or 'unassigned'", sprintf("%d samples", n_illegal))
chk(n_blank == 0, "no cell is left with a BLANK bin", sprintf("%d samples", n_blank))

cat("\n=========== 5.3 the correction rule is what the output says ===========\n")
chk(n_src   == 0, "anno_source reproduces the documented rule", sprintf("%d samples differ", n_src))
chk(n_rule  == 0, "hierarchy_bin == marker iff marker_corrected, else projection", sprintf("%d differ", n_rule))
chk(n_agree == 0, "agree flag == (bin_bmm == bin_marker)", sprintf("%d differ", n_agree))

cat("\n=========== 5.4 rule outcomes are exhaustive and exclusive ===========\n")
A <- rbindlist(lapply(fsx, function(f) fread(f, select = c("anno_source", "agree", "hierarchy_bin",
                                                           "bin_bmm", "bin_marker"))), fill = TRUE)
chk(all(A$anno_source %in% c("concordant", "marker_corrected", "bmm_kept", "bmm_only", "marker_only", "unassigned")),
    "anno_source takes only documented values",
    paste(setdiff(unique(A$anno_source), c("concordant", "marker_corrected", "bmm_kept", "bmm_only", "marker_only", "unassigned")), collapse = ","))
chk(A[anno_source == "concordant" & bin_bmm != bin_marker, .N] == 0, "concordant implies the bins match")
chk(A[anno_source %in% c("marker_corrected","marker_only") & hierarchy_bin != bin_marker, .N] == 0,
    "marker_corrected/marker_only imply the marker bin was taken")
chk(A[anno_source %in% c("bmm_kept", "bmm_only") & hierarchy_bin != bin_bmm, .N] == 0,
    "bmm_kept/bmm_only imply the projection bin was kept")
cat(sprintf("  source mix: %s\n", paste(sprintf("%s=%.1f%%", names(table(A$anno_source)),
            100 * as.numeric(table(A$anno_source)) / nrow(A)), collapse = "  ")))

cat(sprintf("\n=========== BATCH 5: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 5 PASS\n")
