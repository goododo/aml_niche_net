#!/usr/bin/env Rscript
# 03_malignant_distribution_shift.R ----
# Phase 2.3 -- test the LSC-selection arm of the central hypothesis: do MALIGNANT cells redistribute
# toward primitive/stem hierarchy bins along the treatment-pressure axis (Dx -> MRD -> Relapse)?
#
# Metric (per sample): the NORMALIZED malignant distribution across the 7 CCC bins --
#   mal_dist[bin] = n_malignant[bin] / sum(n_malignant over CCC bins)  (sums to 1).
# This deliberately factors OUT blast burden and post-chemo normal-cell regeneration (both of which
# confound the raw per-bin malignant FRACTION); it asks only WHERE the malignant cells sit.
# Summary scores:  stem_frac = HSC_MPP ;  primitive_frac = HSC_MPP + LMPP_GMP.
#
# Design: within-patient PAIRED comparison across timepoints (controls for patient variation). Small
# n (a handful of complete triples) -> per-patient trajectories are PRIMARY; Wilcoxon signed-rank
# tests are secondary and explicitly underpowered.
#
# INPUT : HIER_TAB_DIR/per_bin_malignant.csv  (from 02)
# OUTPUT: HIER_TAB_DIR/malignant_distribution.csv        (per sample: full 7-bin distribution + scores)
#         HIER_TAB_DIR/distribution_shift_tests.csv      (paired test results)

suppressPackageStartupMessages({ library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--min_malignant", type = "integer", default = 50L,
              help = "samples with fewer malignant cells (over CCC bins) are too noisy -> excluded [50]")
)))

CCC_BINS   <- setdiff(HIERARCHY_BINS, "Stromal")         # 7 hematopoietic CCC nodes
STEM_BINS  <- "HSC_MPP"
PRIM_BINS  <- c("HSC_MPP", "LMPP_GMP")
# THE AXIS COMES FROM THE CONFIG, NOT FROM A LITERAL HERE. c("Dx","MRD","Relapse") was a third
# vocabulary, matching neither CANONICAL_TIMEPOINTS nor the Python stages, and "MRD" was retired on
# 2026-08-04. Combined with 02's name-derived timepoint it restricted this whole analysis to the 27
# GSE227903 samples whose names happened to carry a suffix -- and the published
# "stem_frac Dx->MRD, n_pairs=6, p=0.0313" was then described as "all longitudinal samples".
TP_LEVELS  <- TP_AXIS_LEVELS

pb <- fread(file.path(HIER_TAB_DIR, "per_bin_malignant.csv"))
ccc <- pb[in_ccc_graph == TRUE & hierarchy_bin %in% CCC_BINS]

## ---- per-sample normalized malignant distribution ----
ccc[, tot_mal := sum(n_malignant), by = .(dataset, sample)]
ccc[, mal_dist := fifelse(tot_mal > 0, n_malignant / tot_mal, NA_real_)]
# CURATED patient id and timepoint. Deriving the patient by stripping a suffix has the same
# failure as deriving the timepoint that way: it only works for samples whose names encode it, and
# it silently merges or splits patients everywhere else. uid_patient is in the manifest.
ccc <- add_timepoint(ccc)            # -> timepoint, tp_axis, patient, uid_patient (stops on poor coverage)
ccc[, patient := fifelse(!is.na(uid_patient) & nzchar(uid_patient), uid_patient, patient)]
ccc[, timepoint := tp_axis]          # this script works on the 3-level axis

dist <- dcast(ccc, dataset + patient + sample + timepoint + tot_mal ~ hierarchy_bin, value.var = "mal_dist")
for (b in setdiff(CCC_BINS, names(dist))) dist[, (b) := NA_real_]   # ensure all 7 columns exist
dist[, stem_frac      := rowSums(.SD, na.rm = TRUE), .SDcols = STEM_BINS]
dist[, primitive_frac := rowSums(.SD, na.rm = TRUE), .SDcols = PRIM_BINS]
dist[, ok := tot_mal >= opt$min_malignant]

fwrite_safe(dist, file.path(HIER_TAB_DIR, "malignant_distribution.csv"))
n_excl <- dist[timepoint %in% TP_LEVELS & !ok, .N]
message(sprintf("[dist] %d samples; %d longitudinal-timepoint samples excluded (< %d malignant cells)",
                nrow(dist), n_excl, opt$min_malignant))

## ---- drop samples whose cells cannot be attributed to ONE patient ----
# GSE185991 pools 2 patients into 3 libraries with no demux (PT01_PT10D30, PT12_PT13D14,
# PT09_PT13D30). Everything below is per-PATIENT: the composite would enter as a phantom third
# patient carrying a 50/50 mixture of two people's malignant distributions.
.man <- fread(file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv"))
if ("patient_resolved" %in% names(.man)) {
  .un <- .man[patient_resolved == FALSE, .(dataset, sample = Sample)]
  if (nrow(.un)) {
    .before <- nrow(dist)
    .hit <- dist[.un, on = .(dataset, sample), nomatch = 0L]$sample   # those actually IN the set
    dist <- dist[!.un, on = .(dataset, sample)]
    # list what was DROPPED, not every unresolved sample in the manifest. Printing 30 candidate
    # names next to "dropped 0" reads as a silent 30-sample loss.
    message(sprintf("[dist] dropped %d of %d sample(s) with patient_resolved=FALSE%s",
                    .before - nrow(dist), nrow(.un),
                    if (length(.hit)) paste0(": ", paste(head(.hit, 8), collapse = ", ")) else
                      " (none were in the analysis set)"))
  }
} else warning("manifest has no patient_resolved column -- re-run 01_dataset_roles.R", call. = FALSE)

## ---- longitudinal cohort (patients with >=2 timepoints on the treatment axis) ----
# HOW THE PAIR COUNT IS ARRIVED AT, so n_pairs is interpretable rather than merely small:
#   23  Diagnosis -> post-treatment patient pairs exist in the curated manifest
#   16  survive --min_malignant (default 50 malignant cells in BOTH samples)
#   14  survive the restriction to the 7 CCC graph bins (Stromal and Unassigned carry malignant
#       calls but are not nodes, so their cells do not count toward tot_mal)
# Until 2026-08-20 this read 6, from ONE dataset, because the timepoint came from a regex over the
# sample name that only matched GSE227903 suffixes -- and the output was described as "all
# longitudinal samples".
lon <- dist[ok & timepoint %in% TP_LEVELS]
lon[, tp := factor(timepoint, levels = TP_LEVELS)]
setorder(lon, dataset, patient, tp)
pt_tp <- lon[, .(n_tp = uniqueN(tp)), by = .(dataset, patient)]
multi <- pt_tp[n_tp >= 2]
lon <- lon[multi, on = c("dataset", "patient")]
if (!nrow(lon)) stop("no longitudinal patients with >=2 usable timepoints")

cat("\n================ per-patient stem-ward trajectory (fraction of MALIGNANT cells) ================\n")
cat("stem_frac = HSC_MPP share of malignant cells ;  primitive_frac = HSC_MPP + LMPP_GMP share\n\n")
traj_stem <- dcast(lon, dataset + patient ~ tp, value.var = "stem_frac", fun.aggregate = function(x) round(mean(x, na.rm = TRUE), 3))
traj_prim <- dcast(lon, dataset + patient ~ tp, value.var = "primitive_frac", fun.aggregate = function(x) round(mean(x, na.rm = TRUE), 3))
cat("-- stem_frac (HSC_MPP) --\n");            print(traj_stem)
cat("\n-- primitive_frac (HSC_MPP+LMPP_GMP) --\n"); print(traj_prim)

## ---- paired tests (secondary; small n, underpowered) ----
paired_test <- function(score, t1, t2) {
  w <- dcast(lon, patient ~ tp, value.var = score, fun.aggregate = function(x) mean(x, na.rm = TRUE))
  if (!all(c(t1, t2) %in% names(w))) return(NULL)
  pr <- w[is.finite(get(t1)) & is.finite(get(t2))]
  if (nrow(pr) < 2) return(data.table(score = score, comparison = paste(t1, "->", t2), n_pairs = nrow(pr),
                                       median_delta = NA_real_, mean_delta = NA_real_, p_value = NA_real_))
  d  <- pr[[t2]] - pr[[t1]]
  p  <- tryCatch(wilcox.test(pr[[t2]], pr[[t1]], paired = TRUE)$p.value, error = function(e) NA_real_)
  data.table(score = score, comparison = paste(t1, "->", t2), n_pairs = nrow(pr),
             median_delta = round(median(d), 4), mean_delta = round(mean(d), 4),
             n_up = sum(d > 0), n_down = sum(d < 0), p_value = round(p, 4))
}
# COMPARISONS DERIVED FROM THE AXIS, NOT RE-LISTED. These call sites still read
# ("Dx","MRD"),("MRD","Relapse"),("Dx","Relapse") after TP_LEVELS moved to the canonical vocabulary,
# so paired_test()'s `!all(c(t1,t2) %in% names(w))` guard returned NULL for all three and the run
# wrote an EMPTY distribution_shift_tests.csv while every other section printed normally. The guard
# did its job; the literal is what was wrong. Consecutive steps plus first-to-last.
CMP <- c(lapply(seq_len(length(TP_LEVELS) - 1L), function(i) TP_LEVELS[c(i, i + 1L)]),
         if (length(TP_LEVELS) > 2L) list(TP_LEVELS[c(1L, length(TP_LEVELS))]) else NULL)
tests <- rbindlist(lapply(c("stem_frac", "primitive_frac"), function(s)
  rbindlist(lapply(CMP, function(cc) paired_test(s, cc[1], cc[2])), fill = TRUE)), fill = TRUE)

# AN EMPTY RESULT IS A FAILURE, NOT AN OUTPUT. fwrite on a 0-column table warns and leaves the
# previous file in place, so a broken run looks like a successful one that simply changed nothing.
if (!nrow(tests) || !ncol(tests))
  stop(sprintf(paste0("paired tests produced nothing. Axis levels are %s; the trajectory table has ",
                      "columns %s. A comparison naming a level that is not on the axis returns NULL ",
                      "silently."), paste(TP_LEVELS, collapse = "/"),
               paste(setdiff(names(traj_stem), c("dataset", "patient")), collapse = "/")))
fwrite_safe(tests, file.path(HIER_TAB_DIR, "distribution_shift_tests.csv"))

cat("\n================ paired shift tests (Wilcoxon signed-rank; n small -> read as descriptive) ================\n")
cat("median_delta > 0 and n_up > n_down  =>  malignant cells shift TOWARD stem/primitive bins along the axis\n\n")
print(tests)

cat("\n================ pooled malignant distribution by timepoint (context, all longitudinal samples) ================\n")
pooled <- lon[, lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3)),
              by = tp, .SDcols = c(CCC_BINS, "stem_frac", "primitive_frac")]
print(pooled[order(tp)])

cat("\n[note] Interpretation: this tests LSC-state REDISTRIBUTION, not the CCC-network 'topological\n")
cat("       deviation' (that is the later CCC-graph phase). n is small; the per-patient trajectories\n")
cat("       above are the honest primary readout -- look for a consistent stem-ward direction.\n")
