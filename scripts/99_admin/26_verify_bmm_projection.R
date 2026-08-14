#!/usr/bin/env Rscript
# 26_verify_bmm_projection.R ----
# BATCH 0b of the pipeline audit: the BoneMarrowMap projection (01_bmm_project.R).
#
# THIS BATCH EXISTS BECAUSE THE PROJECTION IS NOW LOAD-BEARING. Until 2026-08-14 a confident marker
# call could overwrite the projected bin, so a projection error was partly absorbed downstream.
# That rule was measured against van Galen 2019 and removed (08_validate_annotation.R): the bin now
# comes from the projection for every cell in the cohort except the 3.7% it declines to place. Six
# audit batches covered what happens AFTER the projection and none covered the projection itself.
#
# WHAT CAN GO WRONG HERE WITHOUT ANYTHING FAILING:
#
#  (1) CELL SET DRIFT. The projection is keyed to the QC object it was run from. If a QC object was
#      re-made after its projection, the per-cell table describes cells that no longer exist -- and
#      06 joins on cell id, so the mismatch shows up as missing annotation, not as an error.
#
#  (2) THE ROLL-UP FROM broad -> bin SILENTLY LOSING A TYPE. hierarchy_bin comes from bmm_bin_map.tsv
#      keyed on CellType_Broad. A broad label absent from the map yields NA, the cell ends up
#      unassigned, and the only visible symptom is a slightly smaller graph.
#
#  (3) in_ccc_graph NOT FOLLOWING THE BIN. Stromal is annotated but is NOT a CCC node. If the flag
#      is stale or independently computed, the CCC graph gains non-hematopoietic cells.
#
#  (4) high_error KEYED TO THE WRONG METRIC. config demoted absolute mapping_error in favour of
#      bmm_prob < MIN_PROJ_PROB precisely because the absolute metric is not comparable across
#      platform/depth. A flag still computed from mapping_error would look identical in the column
#      and mean something different -- and it is the flag that decides which cells are trusted.
#
#  (5) A DEGENERATE PROJECTION. Every cell mapped to one or two labels is what a broken reference
#      load looks like. It is not an error; it is a very confident wrong answer.
#
# Everything is recomputed from bmm_bin_map.tsv and the QC objects, never trusted from the output.
#
#   Rscript scripts/99_admin/26_verify_bmm_projection.R [n_samples]

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

NS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(NS)) NS <- 15L
FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

R      <- qc_rds_roster(on_extra = "ignore")
binmap <- fread(BIN_MAP_TSV)
pf_of  <- function(ds, sid) file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
R[, pf := pf_of(dataset, sample)]

cat(sprintf("\n=========== 0b.1 coverage and freshness (%d roster samples) ===========\n", nrow(R)))
chk(all(file.exists(R$pf)), "every roster sample has a projection table",
    sprintf("%d missing", sum(!file.exists(R$pf))))
S <- R[file.exists(pf)]
# (1) a projection older than the QC object it was computed from describes a previous generation
stale <- S[file.mtime(pf) < file.mtime(rds)]
chk(nrow(stale) == 0, "no projection predates its QC object", sprintf("%d stale", nrow(stale)))
if (nrow(stale)) print(head(stale[, .(dataset, sample)], 5))

cat("\n=========== 0b.2 the broad -> bin roll-up is total ===========\n")
# (2) every broad label the projection can emit must have a row in the map. Checked over the whole
# cohort's observed labels, not a sample, because a rare label is exactly the one that slips.
BR <- unique(rbindlist(lapply(S$pf, function(p) fread(p, select = c("bmm_broad", "hierarchy_bin"))), fill = TRUE))
BR <- BR[!is.na(bmm_broad) & nzchar(trimws(bmm_broad))]
unmapped <- setdiff(unique(BR$bmm_broad), binmap$CellType_Broad)
chk(length(unmapped) == 0, "every observed broad label appears in bmm_bin_map.tsv",
    sprintf("unmapped: %s", paste(head(unmapped, 6), collapse = ", ")))
# and the mapping in the output must BE the map, not something similar
J <- merge(BR, unique(binmap[, .(bmm_broad = CellType_Broad, want = hierarchy_bin)]), by = "bmm_broad")
wrong <- J[hierarchy_bin != want]
chk(nrow(wrong) == 0, "hierarchy_bin equals what bmm_bin_map.tsv says for that broad label",
    sprintf("%d label(s) rolled up differently", nrow(wrong)))
if (nrow(wrong)) print(head(wrong, 5))
chk(all(BR$hierarchy_bin %in% c(HIERARCHY_BINS, "")), "no bin outside the declared 8-bin vocabulary",
    paste(setdiff(unique(BR$hierarchy_bin), c(HIERARCHY_BINS, "")), collapse = ", "))

cat(sprintf("\n=========== 0b.3 per-cell invariants on %d of %d samples ===========\n",
            min(NS, nrow(S)), nrow(S)))
set.seed(7); pick <- S[sample(.N, min(NS, .N))]
n_cells <- 0L; n_dup <- 0L; n_ccc <- 0L; n_flag <- 0L; n_prob <- 0L; n_lbl <- 0L; n_key <- 0L
n_unscored <- 0L
ccc_map <- unique(binmap[, .(hierarchy_bin, want_ccc = in_ccc_graph)])
for (i in seq_len(nrow(pick))) {
  P <- fread(pick$pf[i])
  qc_cells <- readRDS_cellnames <- NULL
  # (1) cell set must equal the QC object's, exactly
  seu_cells <- colnames(readRDS(pick$rds[i]))
  if (!setequal(P$cell, seu_cells)) n_cells <- n_cells + 1L
  if (anyDuplicated(P$cell)) n_dup <- n_dup + 1L
  if (!all(P$sample == pick$sample[i] & P$dataset == pick$dataset[i])) n_key <- n_key + 1L

  # (3) in_ccc_graph must be the bin map's answer for THIS cell's bin
  y <- merge(P[, .(cell, hierarchy_bin, in_ccc_graph)], ccc_map, by = "hierarchy_bin", all.x = TRUE)
  if (any(y$in_ccc_graph != y$want_ccc, na.rm = TRUE)) n_ccc <- n_ccc + 1L

  # (4) high_error must be the PROB rule, not the demoted absolute-error rule.
  # NA-safe on purpose. 3.70% of cohort cells carry bmm_prob = NA -- the projection declined to
  # place them -- and 01 writes high_error = 0 for those. Read literally that says "this projection
  # is fine", which would be the wrong answer for any code selecting trustworthy cells by the flag
  # alone. It is harmless ONLY because those cells also carry no bin, so every bin-keyed filter
  # drops them anyway (asserted separately in 0b.5, which is what makes this safe). No ACTIVE
  # script consumes high_error as a trust filter today; the surviving consumers are all under
  # scripts/以前06_hierarchy/. Re-running the projection to change the flag would cost 214 samples
  # against a hazard nothing currently trips, so the invariant is asserted instead of the flag
  # being rewritten. If a consumer ever starts filtering on high_error, fix 01 first.
  scored <- !is.na(P$bmm_prob)
  if (!identical(as.logical(P$high_error[scored]), as.logical(P$bmm_prob[scored] < MIN_PROJ_PROB)))
    n_flag <- n_flag + 1L
  if (any(as.logical(P$high_error[!scored]), na.rm = TRUE)) n_flag <- n_flag + 1L

  # (5) the invariant that makes the above safe: no cell may carry a bin without a confidence score
  bl <- function(v) is.na(v) | !nzchar(trimws(as.character(v)))
  if (any(!scored & !bl(P$hierarchy_bin)) || any(scored & bl(P$hierarchy_bin))) n_unscored <- n_unscored + 1L

  if (any(P$bmm_prob < 0 | P$bmm_prob > 1, na.rm = TRUE)) n_prob <- n_prob + 1L
  # (5) a projection that emits one or two labels for a whole sample is degenerate, not confident
  if (uniqueN(P$bmm_fine) < 3L) n_lbl <- n_lbl + 1L
  rm(P); gc(verbose = FALSE)
}
chk(n_cells == 0, "projection cell set == QC object cell set (exact)", sprintf("%d samples differ", n_cells))
chk(n_dup  == 0, "no duplicated cell in a projection table", sprintf("%d samples", n_dup))
chk(n_key  == 0, "sample/dataset columns match the file's location", sprintf("%d samples", n_key))
chk(n_ccc  == 0, "in_ccc_graph follows the bin map for the cell's own bin", sprintf("%d samples", n_ccc))
chk(n_flag == 0, sprintf("high_error == (bmm_prob < MIN_PROJ_PROB = %.2f) on every SCORED cell", MIN_PROJ_PROB),
    sprintf("%d samples use a different rule", n_flag))
chk(n_prob == 0, "bmm_prob lies in [0,1]", sprintf("%d samples", n_prob))
chk(n_lbl  == 0, "no sample collapsed onto <3 fine labels (degenerate projection)", sprintf("%d samples", n_lbl))

cat("\n=========== 0b.5 unscored cells carry NO bin (what makes high_error=0 harmless) ===========\n")
chk(n_unscored == 0, "bmm_prob is NA if and only if the bin is blank -- no bin without a score",
    sprintf("%d samples break the invariant; high_error=0 would then mark a real bin as trustworthy", n_unscored))

cat("\n=========== 0b.6 the checks are not vacuous ===========\n")
# Anchoring these to ONE sample makes them fail whenever that sample happens to lack the feature
# being probed (BM4 has no stroma). Ask across every sample audited instead.
PK <- rbindlist(lapply(pick$pf, function(p)
  fread(p, select = c("hierarchy_bin", "in_ccc_graph", "bmm_prob", "bmm_fine", "high_error", "mapping_error"))),
  fill = TRUE)
chk(uniqueN(PK$hierarchy_bin) > 1, "the audited samples span several bins (so a roll-up error WOULD show)",
    sprintf("%d bins", uniqueN(PK$hierarchy_bin)))
chk(binmap[in_ccc_graph == FALSE, .N] > 0 && any(PK$in_ccc_graph == FALSE),
    "at least one non-CCC-graph cell is present (so check 0b.3 has something to catch)",
    sprintf("%d non-CCC cells across the audited samples", sum(PK$in_ccc_graph == FALSE, na.rm = TRUE)))
chk(sum(PK$bmm_prob < MIN_PROJ_PROB, na.rm = TRUE) > 0 && sum(PK$bmm_prob >= MIN_PROJ_PROB, na.rm = TRUE) > 0,
    "the audited samples contain BOTH flagged and unflagged cells",
    sprintf("%d flagged of %d", sum(PK$bmm_prob < MIN_PROJ_PROB, na.rm = TRUE), nrow(PK)))
chk(sum(is.na(PK$bmm_prob)) > 0, "unscored cells are present (so 0b.5 is not vacuous)",
    sprintf("%d unscored", sum(is.na(PK$bmm_prob))))
P <- fread(pick$pf[1])
# negative test: the flag comparison must FIRE against the demoted metric
alt <- P$mapping_error > as.numeric(quantile(P$mapping_error, MAPPING_QC_QUANTILE, na.rm = TRUE))
chk(!identical(as.logical(P$high_error), as.logical(alt)),
    "the flag test DISTINGUISHES the prob rule from the mapping_error rule")
cat(sprintf("  (ran on %s :: %s, %d cells, %d fine labels)\n",
            pick$dataset[1], pick$sample[1], nrow(P), uniqueN(P$bmm_fine)))

cat(sprintf("\n=========== BATCH 0b: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 0b PASS\n")
