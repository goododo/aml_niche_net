#!/usr/bin/env Rscript
# 24_verify_reconciliation.R ----
# BATCH 5 of the pipeline audit: the join that produces the corrected label (06).
#
# This is the table every downstream step is meant to read, so an error here is an error in
# everything. The specific hazards:
#   * in_ccc_graph following a different bin than the one the cell ended up with -- the CCC graph
#     would silently gain or lose cells the hierarchy says are not in it.
#   * anno_source and hierarchy_bin disagreeing. Both are legal values, so nothing errors.
#   * cells lost or duplicated in the merge -- a per-cell table with the wrong N still aggregates.
#   * a marker SUBTYPE surviving on a cell whose bin the marker disagreed with. That is how the
#     rejected override rule would sneak back in through the subtype column instead of the bin.
# Everything below is recomputed from the two SOURCE tables rather than trusted from the output.
#
# UPDATED 2026-08-14 for the rule change in 06 (bin always from the projection; markers supply a
# subtype inside that bin). 5.4 carries an explicit regression guard against the old override rule.
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

n_lost <- 0L; n_dup <- 0L; n_ccc <- 0L; n_rule <- 0L; n_src <- 0L; n_agree <- 0L; n_illegal <- 0L; n_sub <- 0L
n_stale <- 0L; n_blank <- 0L
for (f in fsx) {
  X  <- fread(f)
  ds <- X$dataset[1]; sid <- X$sample[1]
  P  <- fread(file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv")),
              select = c("cell", "hierarchy_bin"))
  mf <- file.path(MARKER_ANNO_DIR, ds, paste0(sid, "__marker_percell.csv"))
  M  <- if (file.exists(mf)) fread(mf, select = c("cell", "hierarchy_bin", "low_confidence",
                                                  "marker_cell_type")) else NULL
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
                fifelse(is.na(bb), "marker_fill",
                 fifelse(is.na(bm), "bmm_only",
                  fifelse(ag, "concordant", "bmm_over_marker"))))]
    R[, fin := fifelse(src == "marker_fill", bm, bb)]
    R[src == "unassigned", fin := "unassigned"]
    # the subtype gate, recomputed independently: a marker subtype survives ONLY where the marker's
    # own bin equals the bin the cell ended up in
    R[, sub := fifelse(!is.na(bm) & bm == fin & !bl(marker_cell_type), marker_cell_type, NA_character_)]
    # fwrite writes NA as an empty field and fread reads it back as "", so a subtype-less cell
    # arrives here as "" rather than NA. Compare on the normalised value, which is what any
    # consumer must do too -- 07 does the same normalisation when it loads this table.
    Z <- merge(X[, .(cell, anno_source, agree, hierarchy_bin, bin_bmm, bin_marker,
                     cell_subtype = fifelse(bl(cell_subtype), NA_character_, cell_subtype))],
               R, by = "cell")
    if (any(Z$anno_source != Z$src))          n_src   <- n_src + 1L
    if (any(Z$hierarchy_bin != Z$fin))        n_rule  <- n_rule + 1L
    if (any(Z$agree != Z$ag))                 n_agree <- n_agree + 1L
    if (!identical(as.character(Z$cell_subtype), as.character(Z$sub))) n_sub <- n_sub + 1L
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

cat("\n=========== 5.3 the rule is what the output says ===========\n")
chk(n_src   == 0, "anno_source reproduces the documented rule", sprintf("%d samples differ", n_src))
chk(n_rule  == 0, "hierarchy_bin == projection, except marker_fill", sprintf("%d differ", n_rule))
chk(n_agree == 0, "agree flag == (bin_bmm == bin_marker)", sprintf("%d differ", n_agree))
chk(n_sub   == 0, "cell_subtype reproduces the bin-consistency gate", sprintf("%d differ", n_sub))

cat("\n=========== 5.4 rule outcomes are exhaustive and exclusive ===========\n")
LEGAL <- c("concordant", "bmm_over_marker", "bmm_only", "marker_fill", "unassigned")
A <- rbindlist(lapply(fsx, function(f) fread(f, select = c("anno_source", "agree", "hierarchy_bin",
                                                           "bin_bmm", "bin_marker", "cell_subtype",
                                                           "subtype_dropped", "marker_cell_type"))), fill = TRUE)
# same NA/"" round-trip as above: normalise before any is.na() test means what it says
for (cc in c("cell_subtype", "bin_bmm", "bin_marker", "marker_cell_type"))
  set(A, which(!nzchar(trimws(A[[cc]]))), cc, NA_character_)
chk(all(A$anno_source %in% LEGAL), "anno_source takes only documented values",
    paste(setdiff(unique(A$anno_source), LEGAL), collapse = ","))
chk(A[anno_source == "concordant" & bin_bmm != bin_marker, .N] == 0, "concordant implies the bins match")
chk(A[anno_source == "marker_fill" & hierarchy_bin != bin_marker, .N] == 0,
    "marker_fill implies the marker bin was taken")
chk(A[anno_source %in% c("bmm_over_marker", "bmm_only", "concordant") & hierarchy_bin != bin_bmm, .N] == 0,
    "every other outcome implies the PROJECTION bin was taken")
# THE REGRESSION GUARD. The old rule let a confident marker call overwrite the projection; it was
# measured against van Galen 2019 and lost (overrides wrong 5.5x more often than right). If a
# hierarchy_bin ever again differs from bin_bmm while bin_bmm exists, that rule is back.
chk(A[!is.na(bin_bmm) & nzchar(trimws(bin_bmm)) & hierarchy_bin != bin_bmm, .N] == 0,
    "NO cell's bin overrides an existing projection bin (the rejected rule has not returned)",
    sprintf("%d cells overridden", A[!is.na(bin_bmm) & nzchar(trimws(bin_bmm)) & hierarchy_bin != bin_bmm, .N]))

cat("\n=========== 5.5 the subtype gate does its job ===========\n")
chk(A[!is.na(cell_subtype) & bin_marker != hierarchy_bin, .N] == 0,
    "no retained subtype contradicts the bin it sits in")
chk(A[subtype_dropped == TRUE & !is.na(cell_subtype), .N] == 0, "dropped and retained are exclusive")
chk(A[!is.na(cell_subtype), .N] > 0 && A[subtype_dropped == TRUE, .N] > 0,
    "the gate is NOT vacuous (it both keeps and drops)",
    sprintf("kept %d dropped %d", A[!is.na(cell_subtype), .N], A[subtype_dropped == TRUE, .N]))
chk(A[hierarchy_bin == "Stromal" & !is.na(cell_subtype), uniqueN(cell_subtype)] > 1,
    "stroma is actually resolved into subtypes (the reason the marker route is kept)",
    sprintf("%d subtypes", A[hierarchy_bin == "Stromal" & !is.na(cell_subtype), uniqueN(cell_subtype)]))
cat(sprintf("  source mix: %s\n", paste(sprintf("%s=%.1f%%", names(table(A$anno_source)),
            100 * as.numeric(table(A$anno_source)) / nrow(A)), collapse = "  ")))

cat(sprintf("\n=========== BATCH 5: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 5 PASS\n")
