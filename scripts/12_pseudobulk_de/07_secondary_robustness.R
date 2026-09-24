#!/usr/bin/env Rscript
# 06_secondary_robustness.R ----
# Closes the two gaps that affect whether the SECONDARY tier's hits are reportable at all:
# PREREG 8.5's leave-one-out label (never run) and PREREG 8.3's permutation p-value (run at a
# resolution the design does not support). Touches nothing in the primary tier.
#
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv, PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds,
#          DIR_PSEUDOBULK/depth/*.csv, DIR_PSEUDOBULK/hitlist_frozen.csv, de_summary.csv
# OUTPUT : DIR_PSEUDOBULK/secondary_leave_one_out.csv   per frozen hit: folds survived, robust|contingent
#          DIR_PSEUDOBULK/secondary_exhaustive_perm.csv exact permutation p, all arrangements enumerated
# Usage  : Rscript scripts/12_pseudobulk_de/06_secondary_robustness.R
#
# GAP 1 -- PREREG 8.5 WAS NEVER RUN. It says: "within the secondary tier, leave-one-out over each
#   healthy sample in turn. A secondary hit is reported as ROBUST only if it remains a hit (§6.5,
#   both conditions) in EVERY leave-one-out fit; otherwise CONTINGENT, with the number of folds in
#   which it survived." LMPP_GMP carries 304 reportable hits from an n = 5 vs 4 Chen2023-only fit and
#   not one of them has that label. Until it does, none of them is reportable under the document this
#   analysis is bound by. Each fold is a FULL refit -- filterByExpr, TMM, voom, lmFit, eBayes -- on
#   the reduced sample set, because dropping a sample changes the gene universe and the
#   mean-variance trend, and reusing the full-data versions would hide exactly the fragility the
#   check is for.
#
# GAP 2 -- PREREG 8.3's PERMUTATION p FOR THIS TIER IS FALSE PRECISION. 03_permute.R drew N_PERM =
#   1000 arm labelings WITH REPLACEMENT from a space that holds only C(9,4) = 126 distinct
#   arrangements for LMPP_GMP and C(8,4) = 70 for HSC_MPP. It reported p = 0.00699 for LMPP_GMP,
#   BELOW the design's own minimum attainable p of 1/127 = 0.00787 -- an impossible value, and the
#   tell that the sampling was the problem. Here the space is ENUMERATED EXACTLY, so p is
#   (#arrangements with hits >= observed) / (#arrangements), which cannot undershoot its resolution.
#
#   A second consequence, worth stating because 03_permute.R's output implies otherwise: for this
#   tier the two pre-registered nulls are THE SAME NULL. Chen2023's library_id is 1:1 with sample
#   (verified: 9 samples, 9 libraries), so "permute samples within dataset" and "permute arm-pure
#   libraries within dataset" enumerate an identical set. 03_permute.R's 0.00699 vs 0.00899 for
#   LMPP_GMP was two noisy estimates of one quantity, not two different tests. Reported once here.
#
# BOTH SCHEMES ARE RUN AND BOTH REPORTED. Scheme "frozen_voom" reuses the observed fit's voom object
# the way 03_permute.R did, so the exhaustive p is comparable to the primary tier's. Scheme
# "full_refit" recomputes filterByExpr/TMM/voom per arrangement, which is the more correct null.
# If they agree, the freezing shortcut declared in 03_permute.R's header is vindicated for this tier;
# if they disagree, that is a finding about the shortcut and is reported as one.

suppressPackageStartupMessages({
  library(data.table); library(here); library(limma); library(edgeR); library(statmod)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))

Q_CUT <- 0.05; LFC_CUT <- 1.0                      # PREREG 6.5, unchanged
MIN_CELLS <- CCC_MIN_CELLS_PER_OCCUPIED_BIN

MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))[in_analysis == TRUE]
DEP <- rbindlist(lapply(list.files(file.path(DIR_PSEUDOBULK, "depth"), full.names = TRUE), fread))
PB  <- lapply(seq_len(nrow(MAN)), function(i)
  readRDS(file.path(PB_RDS_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__pseudobulk.rds"))))
names(PB) <- paste(MAN$dataset, MAN$sample, sep = "|")
DISC <- MAN[split == "Discovery"]
gbase <- Reduce(intersect, lapply(unique(DISC$dataset),
  function(d) PB[[paste(DISC[dataset == d][1]$dataset, DISC[dataset == d][1]$sample, sep = "|")]]$genes))
GO   <- fread(INFERCNV_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
UNIV <- setdiff(gbase, union(GO[chr %in% c("chrY", "Y"), unique(gene)], "XIST"))   # PREREG 6.4

FROZEN <- fread(file.path(DIR_PSEUDOBULK, "hitlist_frozen.csv"))
SEC <- fread(file.path(DIR_PSEUDOBULK, "de_summary.csv"))[gate == "PASS" & tier == "secondary"]
CHEN <- DISC[dataset == "Chen2023"]                                               # PREREG 5.4

bin_matrix <- function(bin) {
  keep <- DEP[hierarchy_bin == bin & n_cells >= MIN_CELLS,
              .(dataset, sample, n_cells, median_umi)][CHEN, on = c("dataset", "sample"), nomatch = 0]
  m <- vapply(seq_len(nrow(keep)), function(i)
    PB[[paste(keep$dataset[i], keep$sample[i], sep = "|")]]$counts[UNIV, bin], numeric(length(UNIV)))
  dimnames(m) <- list(UNIV, keep$sample)
  keep[, arm := relevel(factor(arm), ref = "healthy")]                            # healthy = reference
  list(counts = m, meta = keep)
}

## Full refit: everything recomputed from the counts, as Stage B does for this tier (~ arm, no block).
fit_full <- function(counts, arm) {
  md <- data.table(arm = relevel(factor(as.character(arm), levels = c("healthy", "AML")), ref = "healthy"))
  design <- model.matrix(~ arm, data = md)
  if (!"armAML" %in% colnames(design) || min(table(md$arm)) < 2L) return(NULL)
  y <- DGEList(counts); keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  fit <- eBayes(lmFit(voom(y, design), design))
  as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")
}

## Frozen-voom refit: reuse a precomputed voom object, only the design changes (03_permute.R's scheme).
fit_frozen <- function(v, arm) {
  md <- data.table(arm = relevel(factor(as.character(arm), levels = c("healthy", "AML")), ref = "healthy"))
  design <- model.matrix(~ arm, data = md)
  if (!"armAML" %in% colnames(design) || min(table(md$arm)) < 2L) return(NULL)
  fit <- eBayes(lmFit(v, design))
  as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")
}
n_hits <- function(tt) if (is.null(tt)) NA_integer_ else tt[adj.P.Val < Q_CUT & abs(logFC) >= LFC_CUT, .N]

loo_rows <- list(); ex_rows <- list()
for (bin in SEC$hierarchy_bin) {
  B <- bin_matrix(bin); md <- B$meta
  hits_frozen <- FROZEN[hierarchy_bin == bin]$gene
  obs_tt <- fit_full(B$counts, md$arm); obs_n <- n_hits(obs_tt)
  message(sprintf("\n[%s] n = %d AML vs %d healthy | observed hits %d | frozen list %d",
                  bin, sum(md$arm == "AML"), sum(md$arm == "healthy"), obs_n, length(hits_frozen)))

  ## ---- GAP 1: PREREG 8.5 leave-one-out over each HEALTHY sample ----
  hty <- md[arm == "healthy"]$sample
  if (!length(hits_frozen)) {
    message("  [8.5] not evaluable -- no frozen hits for this bin")
    loo_rows[[bin]] <- data.table(hierarchy_bin = bin, gene = NA_character_, n_folds = length(hty),
                                  folds_survived = NA_integer_, verdict = "not evaluable - no frozen hits")
  } else {
    surv <- setNames(integer(length(hits_frozen)), hits_frozen)
    for (h in hty) {
      i <- md$sample != h
      tt <- fit_full(B$counts[, i, drop = FALSE], md$arm[i])
      if (is.null(tt)) { message("  [8.5] fold dropping ", h, ": degenerate, skipped"); next }
      still <- tt[gene %in% hits_frozen & adj.P.Val < Q_CUT & abs(logFC) >= LFC_CUT]$gene
      surv[still] <- surv[still] + 1L
      message(sprintf("  [8.5] drop %-6s -> n=%d, genes tested %d, frozen hits still hits: %d / %d",
                      h, sum(i), nrow(tt), length(still), length(hits_frozen)))
    }
    loo_rows[[bin]] <- data.table(hierarchy_bin = bin, gene = hits_frozen, n_folds = length(hty),
      folds_survived = as.integer(surv[hits_frozen]),
      verdict = fifelse(surv[hits_frozen] == length(hty), "robust", "contingent"))
    message(sprintf("  [8.5] ROBUST %d / CONTINGENT %d of %d frozen hits",
                    sum(surv == length(hty)), sum(surv < length(hty)), length(hits_frozen)))
  }

  ## ---- GAP 2: PREREG 8.3 exhaustive permutation ----
  n <- nrow(md); k <- sum(md$arm == "healthy")
  combos <- combn(n, k)                                    # every way to choose which samples are healthy
  # Frozen voom built once on the OBSERVED design, matching 03_permute.R's scheme exactly.
  des0 <- model.matrix(~ arm, data = md)
  y0 <- DGEList(B$counts); keep0 <- filterByExpr(y0, des0)
  v0 <- voom(calcNormFactors(y0[keep0, , keep.lib.sizes = FALSE], method = "TMM"), des0)
  for (scheme in c("frozen_voom", "full_refit")) {
    hv <- vapply(seq_len(ncol(combos)), function(j) {
      a <- rep("AML", n); a[combos[, j]] <- "healthy"
      n_hits(if (scheme == "frozen_voom") fit_frozen(v0, a) else fit_full(B$counts, a))
    }, numeric(1))
    p_exact <- sum(hv >= obs_n, na.rm = TRUE) / sum(!is.na(hv))
    ex_rows[[paste(bin, scheme)]] <- data.table(
      hierarchy_bin = bin, scheme = scheme, n_arrangements = ncol(combos),
      n_valid = sum(!is.na(hv)), observed_hits = obs_n,
      null_median = median(hv, na.rm = TRUE), null_mean = round(mean(hv, na.rm = TRUE), 2),
      null_p95 = as.numeric(quantile(hv, 0.95, na.rm = TRUE)), null_max = max(hv, na.rm = TRUE),
      n_ge_observed = sum(hv >= obs_n, na.rm = TRUE),
      p_exact = p_exact, resolution = 1 / sum(!is.na(hv)),
      note = "both PREREG 8.3 nulls coincide here: Chen2023 library_id is 1:1 with sample")
    message(sprintf("  [8.3 %s] exact p = %d/%d = %.5f | null median %g mean %.2f max %d",
                    scheme, sum(hv >= obs_n, na.rm = TRUE), sum(!is.na(hv)), p_exact,
                    median(hv, na.rm = TRUE), mean(hv, na.rm = TRUE), max(hv, na.rm = TRUE)))
  }
}

LOO <- rbindlist(loo_rows, fill = TRUE); fwrite(LOO, file.path(DIR_PSEUDOBULK, "secondary_leave_one_out.csv"))
EX  <- rbindlist(ex_rows,  fill = TRUE); fwrite(EX,  file.path(DIR_PSEUDOBULK, "secondary_exhaustive_perm.csv"))
message("\n================ PREREG 8.5 + 8.3, secondary tier ================")
print(LOO[!is.na(gene), .(frozen_hits = .N, robust = sum(verdict == "robust"),
                          contingent = sum(verdict == "contingent")), by = hierarchy_bin])
print(EX[, .(hierarchy_bin, scheme, observed_hits, n_arrangements, n_ge_observed, p_exact, resolution)])
