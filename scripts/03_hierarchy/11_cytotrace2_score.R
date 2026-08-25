#!/usr/bin/env Rscript
# 11_cytotrace2_score.R ----
# Per-cell developmental potency from CytoTRACE2 -- a SECOND estimate of the differentiation axis,
# independent of the BoneMarrowMap projection.
#
# WHY A SECOND ESTIMATOR. Everything the pipeline currently says about differentiation traces back
# to one Symphony projection onto one reference. That projection is already known to be wrong for a
# specific, measurable population: it places malignant HSC-like and Prog-like blasts into the T_NK
# (299 cells), B_Plasma (112) and Erythroid (1,483) bins. An estimator cannot be checked against
# itself, so the pipeline has no way to notice the next error of that kind.
#
# CytoTRACE2 predicts ABSOLUTE potency from a model trained on annotated developmental data. It
# shares no machinery with Symphony -- no reference embedding, no kNN over reference cells, no
# marker panel. Agreement between the two is therefore evidence about the biology rather than about
# a shared assumption, and disagreement localises which cells the projection should not be trusted for.
#
# WHICH COLUMN IS COMPARABLE ACROSS SAMPLES, AND WHICH IS NOT. CytoTRACE2 returns both. The
# `_Relative` column is rank-normalised WITHIN each run, so between-sample differences in it are an
# artefact of what else was in the sample -- the same trap as AddModuleScore's per-sample control
# bins. `CytoTRACE2_Score` is the absolute potency estimate. Every cross-sample statement below uses
# the absolute score; the BMM concordance is computed WITHIN sample, which avoids the question.
#
# THE DIRECTION IS A REAL TEST. BMM pseudotime increases with differentiation and CytoTRACE2 potency
# decreases with it, so the two must correlate NEGATIVELY. A positive correlation would mean one is
# inverted -- an error that is invisible while only one estimator exists.
#
# OUTPUT: CYTOTRACE_DIR/<dataset>/<sample>__cytotrace2_percell.csv
#           cell, dataset, sample, hierarchy_bin, CytoTRACE2_Score, CytoTRACE2_Potency,
#           CytoTRACE2_Relative, preKNN_CytoTRACE2_Score
#
#   Rscript scripts/03_hierarchy/11_cytotrace2_score.R [--dataset X] [--offset N --stride M] [--force]
#
# --offset/--stride slice the roster for concurrent processes: `--offset 0 --stride 4` takes samples
# 1,5,9,... Four such processes cover the cohort with no coordination and no shared state.

suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(CytoTRACE2) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", default = ""),
  make_option("--offset",  type = "integer",   default = 0L),
  make_option("--stride",  type = "integer",   default = 1L),
  make_option("--force",   action = "store_true", default = FALSE)
)))

R <- qc_rds_roster(on_extra = "ignore")
if (nzchar(opt$dataset)) R <- R[dataset %in% strsplit(opt$dataset, ",")[[1]]]
setorder(R, dataset, sample)                       # deterministic, so slices are disjoint and stable
if (opt$stride > 1L) R <- R[seq.int(opt$offset + 1L, nrow(R), by = opt$stride)]
stopifnot(nrow(R) > 0)
message(sprintf("[0] %d samples in this slice (offset %d, stride %d)", nrow(R), opt$offset, opt$stride))

dst_of  <- function(ds, sid) file.path(CYTOTRACE_DIR, ds, paste0(sid, "__cytotrace2_percell.csv"))
anno_of <- function(ds, sid) file.path(ANNO_RECONCILED_DIR, ds, paste0(sid, "__anno_percell.csv"))

KEEP <- c("CytoTRACE2_Score", "CytoTRACE2_Potency", "CytoTRACE2_Relative", "preKNN_CytoTRACE2_Score")

score_one <- function(rds, ds, sid) {
  seu <- readRDS(rds)
  # cytotrace2 chunks cells internally and throws "wrong sign in 'by' argument" for certain
  # cell-count / batch-size combinations (observed on Chen2023::AML320, 5,626 cells). It is a
  # chunking arithmetic bug, not a property of the data, and it kills the sample outright. Retry
  # with smaller chunks rather than silently dropping the sample from the cohort -- a missing
  # sample here would be invisible except as a smaller denominator.
  run <- function(bs, sbs) cytotrace2(seu, species = "human", is_seurat = TRUE, slot_type = "counts",
                                      batch_size = bs, smooth_batch_size = sbs,
                                      ncores = CYTOTRACE_NCORES, seed = CYTOTRACE_SEED)
  res <- NULL
  for (cfg in CYTOTRACE_CHUNKS) {
    res <- tryCatch(run(cfg[1], cfg[2]), error = function(e) { message("    [retry] ", conditionMessage(e)); NULL })
    if (!is.null(res)) break
  }
  if (is.null(res)) stop("cytotrace2 failed at every chunk size")
  m <- if (inherits(res, "Seurat")) res@meta.data else as.data.frame(res)
  got <- intersect(KEEP, names(m))
  if (!length(got)) stop("cytotrace2 returned none of the expected columns")
  A <- as.data.table(m[, got, drop = FALSE], keep.rownames = "cell")
  # the potency call is an ordered factor; keep it as text so the CSV round-trip cannot reorder it
  if ("CytoTRACE2_Potency" %in% names(A)) A[, CytoTRACE2_Potency := as.character(CytoTRACE2_Potency)]
  af <- anno_of(ds, sid)
  if (file.exists(af)) {
    ann <- fread(af, select = c("cell", "hierarchy_bin"))
    # fwrite writes NA as an empty string and fread returns "", not NA -- normalise on the way in,
    # or every is.na(hierarchy_bin) below is silently FALSE
    ann[!nzchar(trimws(hierarchy_bin)), hierarchy_bin := NA_character_]
    A <- merge(A, ann, by = "cell", all.x = TRUE)
  } else A[, hierarchy_bin := NA_character_]
  A[, `:=`(dataset = ds, sample = sid)]
  rm(seu, res, m); gc(verbose = FALSE)
  A[, c("cell", "dataset", "sample", "hierarchy_bin", got), with = FALSE]
}

t0 <- Sys.time()
for (i in seq_len(nrow(R))) {
  ds <- R$dataset[i]; sid <- R$sample[i]; dst <- dst_of(ds, sid)
  .ins <- c(R$rds[i], anno_of(ds, sid))
  if (!is_stale(dst, .ins, force = opt$force)) {
    message(sprintf("[%d/%d] %s::%s current", i, nrow(R), ds, sid)); next
  }
  if (file.exists(dst))
    message(sprintf("[%d/%d] %s::%s RECOMPUTE -- %s", i, nrow(R), ds, sid, stale_reason(dst, .ins, force = opt$force)))
  message(sprintf("[%d/%d] %s::%s", i, nrow(R), ds, sid))
  x <- tryCatch(score_one(R$rds[i], ds, sid),
                error = function(e) { message("  [FAIL] ", conditionMessage(e)); NULL })
  if (is.null(x)) next
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  fwrite_safe(x, dst)
  if (i == 1) message(sprintf("  [timing] first sample took %.1f min",
                              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---------------------------------------------------------------------------------------------
## validation -- only when the whole cohort is present, so a slice cannot report a cohort verdict
## ---------------------------------------------------------------------------------------------
FULL <- qc_rds_roster(on_extra = "ignore")
have <- FULL[file.exists(dst_of(dataset, sample))]
cat(sprintf("\n================ %d of %d samples scored ================\n", nrow(have), nrow(FULL)))
if (nrow(have) < nrow(FULL)) {
  cat("[hold] cohort incomplete -- validation runs when every sample is present.\n")
  cat("       A verdict computed on a partial cohort is the narrowed-check failure mode.\n")
  quit(save = "no", status = 0)
}

P <- rbindlist(lapply(seq_len(nrow(have)), function(i)
  fread(dst_of(have$dataset[i], have$sample[i]))), fill = TRUE)
cat(sprintf("%d cells across %d samples\n", nrow(P), uniqueN(P$sample)))
SC <- CYTOTRACE_SCORE_COL
stopifnot(SC %in% names(P))
FAIL <- 0L

## -- hierarchy_bin is a CONTAMINATED lineage label for malignant cells -------------------------
# Using van Galen's own labels, 64.3% of the cells in the Erythroid bin are malignant MYELOID
# blasts, and the T_NK and B_Plasma bins hold 299 and 112 more. config_malignancy.R already refuses
# to pair malignant vs normal within hierarchy_bin for exactly this reason. So "HSC_MPP > Erythroid
# in potency" can pass or fail on how much myeloid blast sits in the Erythroid bin rather than on
# whether CytoTRACE2 works -- a positive control that measures the wrong thing.
# The contrasts are therefore evaluated on NORMAL cells, where the bin means what it says. The
# malignant stratum is printed beside it because the difference is itself informative.
cons_of <- function(ds, sid) file.path(DIR_MALIGNANCY, ds, paste0(sid, "__consensus_percell.csv"))
CN <- rbindlist(lapply(seq_len(nrow(have)), function(i) {
  f <- cons_of(have$dataset[i], have$sample[i]); if (!file.exists(f)) return(NULL)
  fread(f, select = c("cell", "malignant"))[, sample := have$sample[i]][]
}), fill = TRUE)
if (nrow(CN)) P <- merge(P, CN, by = c("cell", "sample"), all.x = TRUE) else P[, malignant := NA_integer_]
cat(sprintf("\nmalignancy label joined for %d of %d cells (%.1f%%)\n",
            P[!is.na(malignant), .N], nrow(P), 100 * P[!is.na(malignant), .N] / nrow(P)))

bin_tab <- function(D) D[!is.na(hierarchy_bin) & nzchar(hierarchy_bin),
                         .(n = .N, mean_score = round(mean(get(SC), na.rm = TRUE), 4),
                           median_score = round(median(get(SC), na.rm = TRUE), 4)),
                         by = hierarchy_bin][order(-mean_score)]

cat("\n================ potency by bin, NORMAL cells (higher = less differentiated) ================\n")
B <- bin_tab(P[!is.na(malignant) & malignant == 0L])
print(B)
cat("\n================ potency by bin, MALIGNANT cells (bin label unreliable here) ================\n")
print(bin_tab(P[!is.na(malignant) & malignant == 1L]))
cat("\n================ potency by bin, ALL cells (contaminated; reference only) ================\n")
print(bin_tab(P))

cat("\n================ ORDERED CONTRASTS against settled haematopoiesis (NORMAL cells) ================\n")
if (!nrow(B)) stop("no normal-cell bin table -- the consensus join produced nothing, so the positive ",
                   "control cannot be evaluated. Refusing to print a verdict on an empty set.")
CT <- as.data.table(CYTOTRACE_CONTRASTS)
mu <- setNames(B$mean_score, B$hierarchy_bin)
CT[, hi_score := mu[hi]][, lo_score := mu[lo]][, delta := hi_score - lo_score]
CT[, pass := !is.na(delta) & delta > 0]
print(CT[, .(hi, lo, hi_score, lo_score, delta = round(delta, 4), pass, why)])
nbad <- CT[pass == FALSE | is.na(pass), .N]
if (nbad > 0) FAIL <- FAIL + 1L
cat(sprintf("  [%s] %d of %d contrasts in the expected direction\n",
            if (nbad == 0) "PASS" else "FAIL", nrow(CT) - nbad, nrow(CT)))
cat("  Ordered contrasts, not an argmax: a single unexpected peak fails an argmax for reasons that\n")
cat("  have nothing to do with whether the estimator works.\n")
cat(sprintf("  Evaluated on %d NORMAL cells; the malignant stratum is excluded because its bin label\n",
            P[!is.na(malignant) & malignant == 0L, .N]))
cat("  is a projection artefact rather than a lineage.\n")

cat("\n================ concordance with BMM pseudotime (WITHIN sample) ================\n")
pt_of <- function(ds, sid) file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
ptf <- have[file.exists(pt_of(dataset, sample))]
has_pt <- FALSE
if (nrow(ptf)) {
  hdr <- names(fread(pt_of(ptf$dataset[1], ptf$sample[1]), nrows = 0L))
  ptcol <- intersect(c("predicted_Pseudotime", "Pseudotime", "pseudotime"), hdr)
  has_pt <- length(ptcol) > 0
}
if (!has_pt) {
  cat("  [SKIP] the projection carries no pseudotime column yet (01_bmm_project.R does not compute\n")
  cat("         one). Until it does, CytoTRACE2 is an UNCHECKED second estimate, not a validated one.\n")
} else {
  ptcol <- ptcol[1]
  # OFF-TRAJECTORY CELLS ARE EXCLUDED, and this is not a convenience. BMM pseudotime exists only
  # along the HSPC trajectory; the six classes in BMM_PSEUDOTIME_OFFTRAJ carry a placeholder zero,
  # identical to HSC/MPP, despite being terminally differentiated (see config_hierarchy.R for the
  # reference measurement). Comparing a potency estimate against a placeholder is not a test of
  # either estimator. Both columns are reported so the contamination stays visible.
  # Guard on COMPLETE PAIRS, not on row count. A sample can hold plenty of on-trajectory rows and
  # still have no usable pair: predicted_Pseudotime is NA wherever mapping_error_QC == "Fail".
  .rho <- function(a, b, min_pairs = 50L) {
    keep <- is.finite(a) & is.finite(b)
    if (sum(keep) < min_pairs) return(NA_real_)
    suppressWarnings(cor(a[keep], b[keep], method = "spearman"))
  }
  CO <- rbindlist(lapply(seq_len(nrow(ptf)), function(i) {
    d <- fread(pt_of(ptf$dataset[i], ptf$sample[i]), select = c("cell", ptcol, "bmm_broad"))
    j <- merge(P[sample == ptf$sample[i], .(cell, sc = get(SC))], d, by = "cell")
    on <- j[!(bmm_broad %in% BMM_PSEUDOTIME_OFFTRAJ)]
    if (nrow(j) < 50) return(NULL)
    data.table(sample = ptf$sample[i], n = nrow(j), n_on = nrow(on),
               rho_all = .rho(j$sc,  j[[ptcol]]),
               rho     = .rho(on$sc, on[[ptcol]]))
  }), fill = TRUE)
  n_dropped <- CO[is.na(rho), .N]
  if (n_dropped > 0)
    cat(sprintf("  [note] %d sample(s) have <50 usable on-trajectory pairs and are not scored\n", n_dropped))
  CO <- CO[!is.na(rho)]
  # NON-VACUITY: if the filter removed nothing, the label spelling has drifted from the reference
  # and the gate below would be the old, broken comparison wearing a new name.
  n_removed <- sum(CO$n) - sum(CO$n_on)
  if (n_removed <= 0)
    stop("BMM_PSEUDOTIME_OFFTRAJ matched 0 cells -- bmm_broad labels no longer match the reference ",
         "class names. Fix the list; do not let the gate pass on an unfiltered comparison.")
  cat(sprintf("  samples with both estimates: %d | cells %s, on-trajectory %s (%.0f%% kept)\n",
              nrow(CO), format(sum(CO$n), big.mark = ","), format(sum(CO$n_on), big.mark = ","),
              100 * sum(CO$n_on) / sum(CO$n)))
  cat(sprintf("  median within-sample Spearman: on-trajectory %+.3f | all cells %+.3f\n",
              median(CO$rho), median(CO$rho_all, na.rm = TRUE)))
  cat(sprintf("  samples with the expected NEGATIVE correlation: %d / %d\n", CO[rho < 0, .N], nrow(CO)))
  okc <- nrow(CO) > 0 && median(CO$rho) <= CYTOTRACE_MAX_BMM_COR
  if (!okc) FAIL <- FAIL + 1L
  cat(sprintf("  [%s] median within-sample correlation <= %.2f\n",
              if (okc) "PASS" else "FAIL", CYTOTRACE_MAX_BMM_COR))
  cat("  Computed WITHIN sample on purpose: CytoTRACE2_Relative is rank-normalised per run, and\n")
  cat("  pooling cells across samples would mix that artefact into the correlation.\n")
  fwrite(CO, file.path(HIER_TAB_DIR, "cytotrace2_bmm_concordance.csv"))
}

dir.create(HIER_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
fwrite(B, file.path(HIER_TAB_DIR, "cytotrace2_by_bin.csv"))
cat(sprintf("\n[done] %s\n", file.path(HIER_TAB_DIR, "cytotrace2_by_bin.csv")))
if (FAIL > 0) { cat(sprintf("\n[!] %d gate(s) FAILED -- do not use CytoTRACE2 until resolved.\n", FAIL))
                quit(save = "no", status = 1) }
cat("all gates PASS\n")
