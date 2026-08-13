#!/usr/bin/env Rscript
# 93_matched_ref_verdict.R ----
# The decision on the dataset-matched healthy reference, on both sides at once.
#
# A reference change can lower the healthy false-positive rate two ways: by removing a batch
# artefact (good), or by pulling ALL burdens down so malignant cells stop clearing the threshold
# too (worthless). The healthy cohort alone cannot tell those apart -- it has no positives. So the
# verdict needs both arms:
#
#   SPECIFICITY  healthy-donor FPR, 18 donors across 3 datasets, leave-one-out enforced
#   SENSITIVITY  van Galen single-cell genotyping, 18 GSE116256 AML samples -- the only DNA-level
#                truth in the project, and the only place a sensitivity loss would be visible
#
# The criteria below were written BEFORE the runs finished, so the outcome is read against a fixed
# rule rather than a rule chosen to fit it:
#   KEEP     sensitivity does not fall by more than SENS_TOL and mean healthy FPR falls
#   ROLLBACK sensitivity falls by more than SENS_TOL   (FPR bought with real signal)
#   DISCUSS  anything else (e.g. both flat -- the change would then be doing nothing for AML calls)
#
#   Rscript scripts/02_malignancy/93_matched_ref_verdict.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))

BASE <- file.path(LARGE1_DIR, "05_cnv_snv", "_bmmref_baseline__infercnv_burden")
RAW  <- file.path(LARGE1_DIR, "00_raw", "public", "GEO", "GSE116256_RAW")
SENS_TOL <- 0.05      # absolute sensitivity a change may cost before it is a rollback

burden_path <- function(root, ds, sid) file.path(root, ds, paste0(sid, "_infercnv_burden.csv"))

## ---------------- specificity: healthy-donor FPR ----------------
E <- fread("infercnv_tasks_matchedref_eval.tsv", header = FALSE, col.names = c("dataset", "sample"))
fpr_of <- function(p) {
  if (!file.exists(p)) return(NULL)
  d <- fread(p); ir <- grepl("reference", d$group, ignore.case = TRUE); ob <- d[group == "observation"]
  if (!nrow(ob) || !any(ir)) return(NULL)
  thr <- quantile(d$infercnv_burden[ir], INFERCNV_SCORE_Q, na.rm = TRUE)
  list(fpr = mean(ob$infercnv_burden > thr),
       matched = sum(grepl("^reference_matched__", d$group)), n_obs = nrow(ob))
}
S <- rbindlist(lapply(seq_len(nrow(E)), function(i) {
  nw <- fpr_of(burden_path(INFERCNV_BURDEN_ROOT, E$dataset[i], E$sample[i]))
  ol <- fpr_of(burden_path(BASE, E$dataset[i], E$sample[i]))
  if (is.null(nw) || is.null(ol) || nw$matched == 0) return(NULL)
  data.table(dataset = E$dataset[i], sample = E$sample[i], n_obs = nw$n_obs,
             fpr_bmm = round(ol$fpr, 4), fpr_matched = round(nw$fpr, 4),
             delta = round(nw$fpr - ol$fpr, 4))
}), fill = TRUE)

cat("\n================ SPECIFICITY: healthy-donor FPR ================\n")
if (!nrow(S)) { cat("  no re-scored healthy donor yet\n") } else {
  setorder(S, delta); print(S)
  cat(sprintf("\n  donors re-scored: %d of %d\n", nrow(S), nrow(E)))
  cat(sprintf("  mean FPR  %.4f -> %.4f   (%+.4f)\n", mean(S$fpr_bmm), mean(S$fpr_matched), mean(S$delta)))
  cat(sprintf("  improved %d | unchanged %d | worse %d\n",
              sum(S$delta < -1e-4), sum(abs(S$delta) <= 1e-4), sum(S$delta > 1e-4)))
}

## ---------------- sensitivity: single-cell genotyping ----------------
G <- fread("infercnv_tasks_sens_eval.tsv", header = FALSE, col.names = c("dataset", "sample"))
anno <- rbindlist(lapply(list.files(RAW, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE), function(f) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(f)))
  if (!s %in% G$sample) return(NULL)
  d <- fread(cmd = paste("zcat", shQuote(f)), showProgress = FALSE, colClasses = "character")
  if (!"MutTranscripts" %in% names(d)) return(NULL)
  data.table(sample = s, cell = d$Cell, ct = d$CellType,
             mut = nzchar(trimws(d$MutTranscripts)) & trimws(d$MutTranscripts) != "NA",
             nGene = suppressWarnings(as.numeric(d$NumberOfGenes)))
}), fill = TRUE)
LY <- c("T", "B", "NK", "CTL", "Plasma", "ProB")
.core <- function(x) sub("^.*_", "", x)

# How many of the sensitivity samples have actually been re-scored with the matched reference.
# Without this the two arms silently read the SAME baseline files while the AML jobs are still
# running, the deltas come out as exactly 0.0000, and the script prints a confident KEEP based on
# a comparison it never made. A verdict that cannot tell "no change" from "not measured yet" is
# worse than no verdict.
n_rescored <- sum(vapply(unique(anno$sample), function(s) {
  p <- burden_path(INFERCNV_BURDEN_ROOT, "GSE116256", s)
  file.exists(p) && any(grepl("^reference_matched__", fread(p, select = "group")$group))
}, logical(1)))

arm <- function(root, label) {
  B <- rbindlist(lapply(unique(anno$sample), function(s) {
    p <- burden_path(root, "GSE116256", s); if (!file.exists(p)) return(NULL)
    d <- fread(p); d[group == "observation", .(sample = s, cell, burden = infercnv_burden)]
  }), fill = TRUE)
  if (!nrow(B)) return(NULL)
  A <- copy(anno); A[, k := paste(sample, .core(cell))]; B[, k := paste(sample, .core(cell))]
  m <- merge(A, B[, .(k, burden)], by = "k")[!is.na(nGene)]
  pos <- m[mut == TRUE]; neg <- m[ct %in% LY & mut == FALSE]
  auc <- function(p, n) { r <- rank(c(p, n)); (sum(r[seq_along(p)]) - length(p) * (length(p) + 1) / 2) / (length(p) * length(n)) }
  # threshold per sample, exactly as 41 derives it
  thr <- rbindlist(lapply(unique(m$sample), function(s) {
    d <- fread(burden_path(root, "GSE116256", s)); ir <- grepl("reference", d$group, ignore.case = TRUE)
    data.table(sample = s, thr = as.numeric(quantile(d$infercnv_burden[ir], INFERCNV_SCORE_Q, na.rm = TRUE)))
  }), fill = TRUE)
  m <- merge(m, thr, by = "sample")
  data.table(arm = label, samples = uniqueN(m$sample),
             n_pos = nrow(pos), n_neg = nrow(neg),
             sensitivity = round(m[mut == TRUE, mean(burden > thr)], 4),
             flag_normal_lymphoid = round(m[ct %in% LY & mut == FALSE, mean(burden > thr)], 4),
             auc = round(auc(pos$burden, neg$burden), 4))
}
V <- rbindlist(list(arm(BASE, "BMM reference (baseline)"), arm(INFERCNV_BURDEN_ROOT, "dataset-matched")), fill = TRUE)

cat("\n================ SENSITIVITY: single-cell genotyping (GSE116256) ================\n")
cat(sprintf("  AML samples re-scored with the matched reference: %d of %d\n",
            n_rescored, uniqueN(anno$sample)))
if (n_rescored < uniqueN(anno$sample)) {
  cat("  NOT READY -- the live burden files are still (partly) the baseline, so the two arms\n")
  cat("  would be compared against themselves. No verdict until every sample is re-scored.\n")
  if (nrow(V) == 2) { print(V); cat("  (the identical rows above are the artefact just described)\n") }
  quit(save = "no", status = 0)
}
if (nrow(V) < 2) { cat("  sensitivity arm not ready (AML samples still re-scoring)\n") } else {
  print(V)
  d_sens <- V[arm == "dataset-matched"]$sensitivity - V[arm == "BMM reference (baseline)"]$sensitivity
  d_flag <- V[arm == "dataset-matched"]$flag_normal_lymphoid - V[arm == "BMM reference (baseline)"]$flag_normal_lymphoid
  d_auc  <- V[arm == "dataset-matched"]$auc - V[arm == "BMM reference (baseline)"]$auc
  cat(sprintf("\n  sensitivity %+.4f | normal-lymphoid flag rate %+.4f | AUC %+.4f\n", d_sens, d_flag, d_auc))

  cat("\n================ VERDICT (criteria fixed before the runs) ================\n")
  fpr_falls <- nrow(S) > 0 && mean(S$delta) < 0
  if (d_sens < -SENS_TOL) {
    cat(sprintf("  ROLLBACK: sensitivity fell %.4f (> tolerance %.2f). The FPR gain was bought\n", -d_sens, SENS_TOL))
    cat("  with real signal, not with a batch artefact.\n")
  } else if (fpr_falls) {
    cat(sprintf("  KEEP: healthy FPR falls %.4f with sensitivity within tolerance (%+.4f).\n",
                -mean(S$delta), d_sens))
    cat("  Next: re-score the remaining cohort and rebuild 41 -> 50 -> 60 -> 96.\n")
  } else {
    cat("  DISCUSS: sensitivity held but the healthy FPR did not improve -- the change is not\n")
    cat("  doing what it was adopted for; decide whether it earns its complexity.\n")
  }
  cat(sprintf("\n  [note] AUC is reported for context only. Both arms are compared on the SAME\n"))
  cat("  genotyped cells, so a change in AUC reflects the reference, not the truth set.\n")
}
