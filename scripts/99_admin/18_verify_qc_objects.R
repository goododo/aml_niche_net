#!/usr/bin/env Rscript
# 18_verify_qc_objects.R ----
# BATCH 0 of the pipeline audit: the per-sample QC objects, which every later stage reads.
#
# This runs FIRST because everything downstream inherits its errors silently. inferCNV reads the
# counts matrix, the projection reads the same cells, the CCC graph counts them. If a filter did
# not actually apply, or percent.mt was computed against a gene set that matched nothing, nothing
# crashes -- the pipeline just runs to completion on cells that should not be there.
#
# THE FAILURE MODES THIS IS BUILT AROUND, all of which are silent:
#
#  (1) A THRESHOLD THAT DID NOT APPLY. config declares ABS_MIN_NFEAT/ABS_MIN_NCOUNT/ABS_MAX_MT as
#      hard bounds on top of the per-sample MAD. Whether they were actually enforced is only
#      knowable by recomputing the metric from the counts matrix and testing every retained cell.
#      A stored column that agrees with a filter proves nothing if the column itself is wrong.
#
#  (2) percent.mt COMPUTED FROM AN EMPTY GENE SET. add_qc_metrics() matches PAT_MITO ("^MT[-.]")
#      and falls back to MITO_SHORT for datasets whose MT- prefix was stripped upstream. If BOTH
#      miss, percent.mt is 0 for every cell and the ABS_MAX_MT cap passes VACUOUSLY -- no cell is
#      ever filtered, and the sample looks pristine. So the mito gene set is checked for being
#      non-empty per sample, separately from the percentage agreeing.
#
#  (3) A NORMALISED MATRIX IN THE counts LAYER. inferCNV and the CNV burden assume raw integers.
#      Normalised values there would produce a plausible-looking burden distribution that means
#      nothing. Integrality is cheap to test and impossible to notice otherwise.
#
#  (4) REPORT AND DISK DISAGREEING. The QC report is what every summary table quotes; the objects
#      are what the analysis reads. A sample that is PASS in the report with no object on disk (or
#      the reverse) makes the two permanently inconsistent.
#
# Object-level checks recompute from the counts matrix on a SAMPLE of objects (reading all 214 is
# ~25 min); report-level checks run over every sample, since they are cheap. Every threshold comes
# from config_qc.R, never re-typed here -- a check with its own copy of a constant passes forever
# after the constant changes.
#
#   Rscript scripts/99_admin/18_verify_qc_objects.R [n_objects]

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix) })

NS <- as.integer(commandArgs(trailingOnly = TRUE)[1]); if (is.na(NS)) NS <- 12L
FAIL <- 0L; N <- 0L
ok  <- function(m) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", m)) }
bad <- function(m) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", m)) }
chk <- function(cond, m, d = "") if (isTRUE(cond)) ok(m) else bad(paste0(m, if (nzchar(d)) paste0(" -- ", d)))

TAB <- file.path(FAST_DIR, "results", "tables", "01_preprocess")

## ---------------------------------------------------------------------------------------------
## PART A -- report level, every sample
## ---------------------------------------------------------------------------------------------
# 03_qc_report__ALL.csv is a COMBINER over the per-dataset reports, so globbing both double-counts
# every sample. The per-dataset files are the source of truth here; ALL is checked AGAINST them
# below rather than merged with them -- a combined table that has drifted from its parts is how two
# figures in the same paper end up quoting different cell counts.
ALL_REP   <- file.path(TAB, "03_qc_report__ALL.csv")
rep_files <- setdiff(list.files(TAB, pattern = "^03_qc_report__.*\\.csv$", full.names = TRUE), ALL_REP)
exc_files <- list.files(TAB, pattern = "^03_P2_exclude__.*\\.csv$", full.names = TRUE)
Q <- rbindlist(lapply(rep_files, fread), fill = TRUE)
E <- rbindlist(lapply(exc_files, fread), fill = TRUE)
R <- qc_rds_roster(on_extra = "ignore")

cat(sprintf("\n=========== 0.1 report and disk agree (%d report rows, %d objects) ===========\n",
            nrow(Q), nrow(R)))
pass <- Q[status == "PASS"]
chk(nrow(Q) > 0 && nrow(R) > 0, "there is something to check at all")
key <- function(D) paste(D$dataset, D$Sample %||% D$sample)
`%||%` <- function(a, b) if (is.null(a)) b else a
k_pass <- paste(pass$dataset, pass$Sample)
k_disk <- paste(R$dataset, R$sample)
chk(setequal(k_pass, k_disk), "PASS in the report <=> object on disk",
    sprintf("report-only %d, disk-only %d", length(setdiff(k_pass, k_disk)), length(setdiff(k_disk, k_pass))))
chk(!anyDuplicated(k_pass), "no sample appears twice in the reports")
if (nrow(E)) {
  k_exc <- paste(E$dataset, E$Sample)
  chk(length(intersect(k_pass, k_exc)) == 0, "no sample is both PASS and P2-excluded",
      sprintf("%d in both", length(intersect(k_pass, k_exc))))
}
if (file.exists(ALL_REP)) {
  A <- fread(ALL_REP)
  cmp <- c("n_raw", "n_after_mad", "n_doublet", "n_final", "status")
  cmp <- intersect(cmp, intersect(names(A), names(Q)))
  mm <- merge(Q[, c("dataset", "Sample", ..cmp)], A[, c("dataset", "Sample", ..cmp)],
              by = c("dataset", "Sample"), suffixes = c(".part", ".all"))
  diff <- Reduce(`|`, lapply(cmp, function(c) as.character(mm[[paste0(c, ".part")]]) !=
                                              as.character(mm[[paste0(c, ".all")]])))
  chk(nrow(mm) == nrow(Q), "the combined report covers every per-dataset row",
      sprintf("%d of %d matched", nrow(mm), nrow(Q)))
  chk(sum(diff, na.rm = TRUE) == 0, "the combined report agrees with the per-dataset reports",
      sprintf("%d rows drifted", sum(diff, na.rm = TRUE)))
}

cat("\n=========== 0.2 the report's own arithmetic ===========\n")
chk(pass[n_final > n_raw, .N] == 0, "n_final <= n_raw", sprintf("%d samples", pass[n_final > n_raw, .N]))
chk(pass[n_final != n_after_mad - n_doublet, .N] == 0, "n_final == n_after_mad - n_doublet",
    sprintf("%d samples", pass[n_final != n_after_mad - n_doublet, .N]))
chk(pass[n_final < MIN_CELLS_SAMPLE, .N] == 0,
    sprintf("every PASS sample clears MIN_CELLS_SAMPLE (%d)", MIN_CELLS_SAMPLE),
    sprintf("%d below", pass[n_final < MIN_CELLS_SAMPLE, .N]))

cat("\n=========== 0.3 the gates that were switched ON actually fired ===========\n")
if (isTRUE(FLAG_LOWCOMPLEXITY_DROP)) {
  n_lc <- pass[med_nfeat_final < FLAG_MIN_MED_NFEAT | med_ncount_final < FLAG_MIN_MED_NCOUNT, .N]
  chk(n_lc == 0, sprintf("no PASS sample is below the low-complexity floor (%d genes / %d counts)",
                         FLAG_MIN_MED_NFEAT, FLAG_MIN_MED_NCOUNT),
      sprintf("%d samples survived a filter that is switched ON", n_lc))
  # and the gate must have REMOVED something, or it is on in name only
  chk(nrow(E) > 0 && E[grepl("complex", reason, ignore.case = TRUE), .N] > 0,
      "the low-complexity gate excluded at least one sample (not on in name only)",
      sprintf("%d excluded for complexity", if (nrow(E)) E[grepl("complex", reason, ignore.case = TRUE), .N] else 0L))
}
if (identical(DOUBLET_CONSENSUS, "union") && all(c("n_dbl_sc", "n_dbl_df", "n_dbl_both") %in% names(pass))) {
  U <- pass[!is.na(n_dbl_sc) & !is.na(n_dbl_df) & !is.na(n_dbl_both)]
  chk(U[n_doublet != n_dbl_sc + n_dbl_df - n_dbl_both, .N] == 0,
      "n_doublet == |scDblFinder UNION DoubletFinder| (config says union, not intersection)",
      sprintf("%d samples match the INTERSECTION instead", U[n_doublet == n_dbl_both & n_dbl_both != n_dbl_sc + n_dbl_df - n_dbl_both, .N]))
}

## ---------------------------------------------------------------------------------------------
## PART B -- object level, recomputed from the counts matrix
## ---------------------------------------------------------------------------------------------
set.seed(11); pick <- R[sample(.N, min(NS, .N))]
cat(sprintf("\n=========== 0.4 recomputing from counts on %d of %d objects ===========\n",
            nrow(pick), nrow(R)))

n_int <- 0L; n_feat <- 0L; n_cnt <- 0L; n_mt <- 0L; n_mtgenes <- 0L
n_abs <- 0L; n_dup <- 0L; n_align <- 0L; n_sample <- 0L; n_gate <- 0L
worst <- data.table()
for (i in seq_len(nrow(pick))) {
  ds <- pick$dataset[i]; sid <- pick$sample[i]
  seu <- readRDS(pick$rds[i])
  C <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
  md <- seu@meta.data

  # (3) raw integers, not a normalised matrix parked in the counts layer
  sv <- C@x[seq_len(min(2e5, length(C@x)))]
  if (any(abs(sv - round(sv)) > 1e-8)) n_int <- n_int + 1L

  # (1) the stored metrics must be what the matrix says
  f <- Matrix::colSums(C > 0); n <- Matrix::colSums(C)
  if (!isTRUE(all.equal(as.numeric(f), as.numeric(md[[COL_NFEAT]]), tolerance = 1e-8))) n_feat <- n_feat + 1L
  if (!isTRUE(all.equal(as.numeric(n), as.numeric(md[[COL_NCOUNT]]), tolerance = 1e-8))) n_cnt  <- n_cnt  + 1L

  # (2) percent.mt against the SAME gene set add_qc_metrics uses -- and that set must be non-empty
  g  <- rownames(C)
  mg <- unique(c(grep(PAT_MITO, g, value = TRUE), intersect(g, MITO_SHORT)))
  if (!length(mg)) n_mtgenes <- n_mtgenes + 1L else {
    pm <- 100 * Matrix::colSums(C[mg, , drop = FALSE]) / pmax(n, 1)
    if (!isTRUE(all.equal(as.numeric(pm), as.numeric(md[[COL_MT]]), tolerance = 1e-6))) n_mt <- n_mt + 1L
  }

  # (1) every retained cell clears the ABSOLUTE bounds, recomputed -- not read from a column
  v <- sum(f < ABS_MIN_NFEAT) + sum(n < ABS_MIN_NCOUNT)
  if (length(mg)) v <- v + sum(100 * Matrix::colSums(C[mg, , drop = FALSE]) / pmax(n, 1) > ABS_MAX_MT)
  if (v > 0) { n_abs <- n_abs + 1L
    worst <- rbind(worst, data.table(dataset = ds, sample = sid, violating_cells = v,
                                     min_nfeat = min(f), min_ncount = min(n)), fill = TRUE) }

  if (ncol(seu) < MIN_CELLS_SAMPLE) n_gate <- n_gate + 1L
  if (anyDuplicated(colnames(seu))) n_dup <- n_dup + 1L
  if (!identical(rownames(md), colnames(seu))) n_align <- n_align + 1L
  if (!COL_SAMPLE %in% names(md) || !all(as.character(md[[COL_SAMPLE]]) == sid)) n_sample <- n_sample + 1L
  rm(seu, C, md); gc(verbose = FALSE)
}

chk(n_int == 0, "counts layer holds RAW INTEGERS (not a normalised matrix)", sprintf("%d objects", n_int))
chk(n_feat == 0, "nFeature_RNA recomputes from the matrix", sprintf("%d objects differ", n_feat))
chk(n_cnt  == 0, "nCount_RNA recomputes from the matrix", sprintf("%d objects differ", n_cnt))
chk(n_mtgenes == 0, "the mitochondrial gene set is NON-EMPTY (else ABS_MAX_MT passes vacuously)",
    sprintf("%d objects matched no MT gene via PAT_MITO or MITO_SHORT", n_mtgenes))
chk(n_mt == 0, "percent.mt recomputes from the matrix", sprintf("%d objects differ", n_mt))
chk(n_abs == 0, sprintf("every retained cell clears ABS bounds (nFeat>=%d, nCount>=%d, mt<=%g)",
                        ABS_MIN_NFEAT, ABS_MIN_NCOUNT, ABS_MAX_MT),
    sprintf("%d objects contain cells that should have been filtered", n_abs))
if (nrow(worst)) print(head(worst[order(-violating_cells)], 5))
chk(n_gate == 0, "every object still clears MIN_CELLS_SAMPLE", sprintf("%d objects", n_gate))
chk(n_dup == 0, "no duplicated cell barcode within a sample", sprintf("%d objects", n_dup))
chk(n_align == 0, "meta.data rownames align with colnames", sprintf("%d objects", n_align))
chk(n_sample == 0, "the Sample column matches the filename", sprintf("%d objects", n_sample))

## ---------------------------------------------------------------------------------------------
## PART C -- the checks are not vacuous
## ---------------------------------------------------------------------------------------------
cat("\n=========== 0.5 NEGATIVE TESTS: would these checks catch a violation? ===========\n")
# A green suite proves nothing unless a broken input turns it red. Replay the same predicates on
# deliberately corrupted copies and require them to FIRE.
seu <- readRDS(pick$rds[1])
C <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
g <- rownames(C); n <- Matrix::colSums(C); f <- Matrix::colSums(C > 0)

Cn <- C; Cn@x <- Cn@x / 3                                   # a normalised matrix
sv <- Cn@x[seq_len(min(2e5, length(Cn@x)))]
chk(any(abs(sv - round(sv)) > 1e-8), "the integrality test FIRES on a normalised matrix")

fb <- f; fb[1] <- ABS_MIN_NFEAT - 1L                        # one cell below the floor
chk(sum(fb < ABS_MIN_NFEAT) > 0, "the ABS_MIN_NFEAT test FIRES on a single sub-threshold cell")

chk(length(unique(c(grep(PAT_MITO, "FOO", value = TRUE), intersect("FOO", MITO_SHORT)))) == 0,
    "the empty-mito-set test FIRES on a feature space with no MT gene")

mg <- unique(c(grep(PAT_MITO, g, value = TRUE), intersect(g, MITO_SHORT)))
chk(length(mg) > 0 && !isTRUE(all.equal(as.numeric(100 * Matrix::colSums(C[mg, , drop = FALSE]) / pmax(n, 1)),
                                        as.numeric(seu@meta.data[[COL_MT]]) + 1, tolerance = 1e-6)),
    "the percent.mt test FIRES when the stored value is shifted by 1")
cat(sprintf("  (negative tests ran on %s :: %s, %d cells, %d MT genes)\n",
            pick$dataset[1], pick$sample[1], ncol(seu), length(mg)))

cat(sprintf("\n=========== BATCH 0: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 0 PASS\n")
