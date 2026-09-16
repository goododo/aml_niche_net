#!/usr/bin/env Rscript
# 02_de_limma.R ----
# Stage B of the pseudobulk AML-vs-healthy DE. Runs the pre-registered limma-voom model per
# hierarchy bin on the Discovery arm, applies the mandatory depth-sensitivity arm, and freezes a
# hit list. Does NOT touch Validation -- that is step 7 of PREREG 9 and a separate invocation.
#
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv              (roster written by 01)
#          PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds           (genes x bins raw counts, from 01)
#          DIR_PSEUDOBULK/depth/<ds>__<sample>__pb_depth.csv  (per-bin depth, from 01)
#          INFERCNV_GENE_ORDER                                 (gencode GRCh38; chrY identification)
# OUTPUT : DIR_PSEUDOBULK/marker_gate.csv        PREREG 8.2, a GATE -- the run aborts if it fails
#          DIR_PSEUDOBULK/de/<tier>__<bin>.csv   per-bin gene table, baseline + depth arm side by side
#          DIR_PSEUDOBULK/de_summary.csv         per-bin n, universe size, MinSampleSize, hit counts
#          DIR_PSEUDOBULK/hitlist_frozen.csv     the frozen list Validation is allowed to see
# Usage  : Rscript scripts/12_pseudobulk_de/02_de_limma.R
#
# EVERYTHING HERE IS PRE-REGISTERED. PREREGISTRATION_pseudobulk_de.md was committed (a66cbd0)
# before any pseudobulk matrix existed, and amended twice after Stage A but still before any
# p-value existed. This script implements it; it does not decide anything.
#
#   PREREG 5.4  THREE TIERS, NOT SEVEN BINS. Mono_DC / T_NK / B_Plasma are PRIMARY: pooled over
#               Chen2023 + GSE185381 with dataset as a blocking fixed effect. LMPP_GMP / HSC_MPP /
#               Erythroid are SECONDARY: Chen2023 ONLY. The reason is in the cell mass -- in those
#               three bins the Discovery healthy arm is 95.0-97.5% Chen2023 cells while the AML arm
#               is 4.7-13.3%, so pooled, the `arm` coefficient would largely be "Chen2023 healthy vs
#               GSE185381 AML", a between-dataset contrast the `dataset` term cannot absorb because
#               there is almost no within-GSE185381 healthy mass to anchor it. Megakaryocyte is
#               dropped entirely. The two tiers are NEVER pooled in any denominator.
#
#   PREREG 6.1  Primary: two-pass voom + duplicateCorrelation(block = library_id). GSE185381 is
#               HTO-multiplexed and Control1/2/3/4 share ONE lane (GSE185381__lib__P02), so its
#               healthy arm is pseudoreplicated at the lane level. Without the block the healthy n
#               is a fiction. PREREG 6.2 -- secondary tier and Validation are single-dataset,
#               singleton-library designs, so they get `~ arm` with no block and no dataset term
#               (`~ dataset + arm` is rank-deficient on one dataset).
#
#   PREREG 6.4  One common base universe for the whole Discovery side: intersect(Chen2023,
#               GSE185381) minus chrY and XIST = 21843 genes. The deposits genuinely differ
#               (33694 / 36601 / 27899), so intersecting is forced; fixing it in advance is what
#               keeps it from being a denominator chosen after the fact.
#
#   PREREG 6.5  hit := q_within_bin < 0.05 AND |logFC| >= 1. BOTH conditions, everywhere, including
#               the depth gate. BH is WITHIN BIN and within tier.
#
#   PREREG 7.1  A hit that fails EITHER condition once per-bin depth enters the model is
#               depth-dependent and NOT REPORTABLE. Decided before the numbers existed.
#
#   PREREG 8.2  Lineage markers, PER ARM, top-50 by mean CPM within that arm. Pooled across arms
#               the check cannot detect arm-asymmetric bin misassignment, which is the failure that
#               would actually manufacture DE. It is a GATE: a failing primary bin is dropped.
#
# ONE IMPLEMENTATION CHOICE NOT FIXED BY THE PREREG, DECLARED HERE. PREREG 7.1 says "per-bin
# median UMI is added as a covariate" without naming a scale. This script uses log10(median_umi).
# Depth acts multiplicatively on counts, so log is the standard scale and the untransformed
# covariate would under-correct the deep tail. Declared in code, before the run, rather than chosen
# after seeing which hits survive.

suppressPackageStartupMessages({
  library(data.table); library(here); library(limma); library(edgeR); library(statmod)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))   # INFERCNV_GENE_ORDER

PRIMARY_BINS   <- c("Mono_DC", "T_NK", "B_Plasma")                  # PREREG 5.4
SECONDARY_BINS <- c("LMPP_GMP", "HSC_MPP", "Erythroid")
MIN_CELLS      <- CCC_MIN_CELLS_PER_OCCUPIED_BIN                    # 30, PREREG 5.1
Q_CUT          <- 0.05
LFC_CUT        <- 1.0                                               # PREREG 6.5
MARKERS <- list(Mono_DC = c("LYZ", "S100A8"), T_NK = c("CD3E", "CD3D"), B_Plasma = c("MS4A1", "CD79A"),
                Erythroid = c("HBB", "HBA1"), HSC_MPP = "CD34", LMPP_GMP = "CD34")

dir.create(file.path(DIR_PSEUDOBULK, "de"), recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------- roster, depth, universe ----
MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))[in_analysis == TRUE]
DEP <- rbindlist(lapply(list.files(file.path(DIR_PSEUDOBULK, "depth"), full.names = TRUE), fread))
PB  <- lapply(seq_len(nrow(MAN)), function(i)
  readRDS(file.path(PB_RDS_DIR, MAN$dataset[i], paste0(MAN$sample[i], "__pseudobulk.rds"))))
names(PB) <- paste(MAN$dataset, MAN$sample, sep = "|")

DISC <- MAN[split == "Discovery"]
base_universe <- Reduce(intersect, lapply(unique(DISC$dataset),
                                          function(d) PB[[paste(DISC[dataset == d][1]$dataset, DISC[dataset == d][1]$sample, sep = "|")]]$genes))
GO   <- fread(INFERCNV_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
sexg <- union(GO[chr %in% c("chrY", "Y"), unique(gene)], "XIST")    # PREREG 6.4 amendment
UNIV <- setdiff(base_universe, sexg)
message(sprintf("[universe] Discovery intersect = %d genes; %d sex genes removed; base = %d",
                length(base_universe), sum(base_universe %in% sexg), length(UNIV)))

## Build the (genes x samples) count matrix for one bin, over a given roster, with its metadata.
## Only (sample, bin) units clearing MIN_CELLS enter -- PREREG 5.1.
bin_matrix <- function(bin, roster) {
  keep <- DEP[hierarchy_bin == bin & n_cells >= MIN_CELLS,
              .(dataset, sample, n_cells, median_umi)][roster, on = c("dataset", "sample"), nomatch = 0]
  if (!nrow(keep)) return(NULL)
  m <- vapply(seq_len(nrow(keep)), function(i) {
    o <- PB[[paste(keep$dataset[i], keep$sample[i], sep = "|")]]
    o$counts[UNIV, bin]                      # UNIV is a subset of every Discovery sample's genes
  }, numeric(length(UNIV)))
  dimnames(m) <- list(UNIV, paste(keep$dataset, keep$sample, sep = "|"))
  # healthy is the REFERENCE level (PREREG 6.1). R would otherwise sort "AML" first and silently
  # flip the sign of every logFC in this study.
  keep[, arm := relevel(factor(arm), ref = "healthy")]
  keep[, `:=`(dataset = factor(dataset), library_id = factor(library_id),
              log_depth = log10(median_umi))]
  list(counts = m, meta = keep)
}

## ------------------------------------------------ PREREG 8.2 : marker gate ----
## Mean CPM WITHIN each arm, then rank. Pooled ranking cannot fail for the right reason.
marker_gate <- function(bin, B) {
  cpm_all <- edgeR::cpm(B$counts)
  res <- rbindlist(lapply(levels(B$meta$arm), function(a) {
    idx <- which(B$meta$arm == a)
    mu  <- rowMeans(cpm_all[, idx, drop = FALSE])
    rk  <- rank(-mu, ties.method = "min")
    mk  <- MARKERS[[bin]]
    if (identical(mk, "CD34")) {   # progenitor bins: CD34 above the bin's own median, both arms
      data.table(hierarchy_bin = bin, arm = a, marker = "CD34", n_samples = length(idx),
                 rank = rk[["CD34"]], mean_cpm = mu[["CD34"]],
                 criterion = "above bin median", pass = mu[["CD34"]] > median(mu))
    } else {
      data.table(hierarchy_bin = bin, arm = a, marker = mk, n_samples = length(idx),
                 rank = rk[mk], mean_cpm = mu[mk], criterion = "top 50 by mean CPM", pass = rk[mk] <= 50)
    }
  }))
  # "At least one marker of each named pair must clear it, in each arm" (PREREG 8.2).
  res[, arm_pass := any(pass), by = .(hierarchy_bin, arm)]
  res[]
}

## ------------------------------------------------------------ the DE fit ----
## tier "primary"   -> ~ dataset + arm, blocked on library_id, two voom passes (PREREG 6.1)
## tier "secondary" -> ~ arm, no block, no dataset term            (PREREG 6.2)
## with_depth = TRUE appends log10(median_umi)                     (PREREG 7.1)
fit_bin <- function(B, tier, with_depth = FALSE) {
  md <- B$meta
  rhs <- if (tier == "primary") "~ dataset + arm" else "~ arm"
  if (with_depth) rhs <- sub("arm$", "log_depth + arm", rhs)
  design <- model.matrix(as.formula(rhs), data = md)
  if (!"armAML" %in% colnames(design)) stop("arm coefficient missing -- relevel failed for this bin")

  y <- DGEList(B$counts)
  keep <- filterByExpr(y, design)                        # PREREG 6.4: on top of the common base
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  if (tier == "primary") {
    v  <- voom(y, design)
    c1 <- duplicateCorrelation(v, design, block = md$library_id)
    v  <- voom(y, design, block = md$library_id, correlation = c1$consensus)
    c2 <- duplicateCorrelation(v, design, block = md$library_id)
    fit <- lmFit(v, design, block = md$library_id, correlation = c2$consensus)
    consensus <- c2$consensus
  } else {
    v <- voom(y, design); fit <- lmFit(v, design); consensus <- NA_real_
  }
  fit <- eBayes(fit)
  tt  <- as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")
  list(tt = tt, n_genes = sum(keep), min_sample_size = 1 / max(hat(design)),
       consensus = consensus, design = rhs)
}

## ------------------------------------------------------------------ run ----
gate_rows <- list(); de_rows <- list(); summ <- list()
for (tier in c("primary", "secondary")) {
  bins   <- if (tier == "primary") PRIMARY_BINS else SECONDARY_BINS
  roster <- if (tier == "primary") DISC else DISC[dataset == "Chen2023"]   # PREREG 5.4
  for (bin in bins) {
    B <- bin_matrix(bin, roster)
    if (is.null(B)) { message("[skip] ", tier, "/", bin, ": no units clear ", MIN_CELLS, " cells"); next }
    g <- marker_gate(bin, B); gate_rows[[paste(tier, bin)]] <- cbind(tier = tier, g)
    n_arm <- B$meta[, .N, by = arm]
    if (!all(g$arm_pass)) {
      message(sprintf("[GATE FAIL] %s/%s -- markers not recovered in %s arm; bin DROPPED (PREREG 8.2)",
                      tier, bin, paste(unique(g[arm_pass == FALSE]$arm), collapse = "/")))
      summ[[paste(tier, bin)]] <- data.table(tier, hierarchy_bin = bin, gate = "FAIL",
                                             n_AML = n_arm[arm == "AML"]$N, n_healthy = n_arm[arm == "healthy"]$N)
      next
    }
    base <- fit_bin(B, tier, with_depth = FALSE)
    dep  <- fit_bin(B, tier, with_depth = TRUE)
    D <- merge(base$tt[, .(gene, logFC, P.Value, adj.P.Val)],
               dep$tt[,  .(gene, logFC_depth = logFC, P_depth = P.Value, q_depth = adj.P.Val)],
               by = "gene", all.x = TRUE)
    D[, hit_baseline := adj.P.Val < Q_CUT & abs(logFC) >= LFC_CUT]                 # PREREG 6.5
    D[, hit_depth    := q_depth   < Q_CUT & abs(logFC_depth) >= LFC_CUT]           # PREREG 7.1
    D[, reportable   := hit_baseline & hit_depth %in% TRUE]
    D[, `:=`(tier = tier, hierarchy_bin = bin)]
    D <- D[order(adj.P.Val, -abs(logFC))]   # setorder() takes column names only, not expressions
    fwrite(D, file.path(DIR_PSEUDOBULK, "de", paste0(tier, "__", bin, ".csv")))
    de_rows[[paste(tier, bin)]] <- D

    summ[[paste(tier, bin)]] <- data.table(
      tier, hierarchy_bin = bin, gate = "PASS",
      n_AML = n_arm[arm == "AML"]$N, n_healthy = n_arm[arm == "healthy"]$N,
      n_libraries_healthy = B$meta[arm == "healthy", uniqueN(library_id)],
      design = base$design, n_genes_tested = base$n_genes,
      min_sample_size = round(base$min_sample_size, 3),
      dupcor_consensus = round(base$consensus, 4),
      hits_baseline = sum(D$hit_baseline), hits_after_depth = sum(D$reportable),
      struck_by_depth = sum(D$hit_baseline & !D$reportable))
    message(sprintf("[%s/%s] n=%d AML vs %d healthy (%d hty libs) | genes %d | MinSS %.2f | dupcor %s | hits %d -> %d after depth",
                    tier, bin, n_arm[arm == "AML"]$N, n_arm[arm == "healthy"]$N,
                    B$meta[arm == "healthy", uniqueN(library_id)], base$n_genes, base$min_sample_size,
                    ifelse(is.na(base$consensus), "n/a", sprintf("%.4f", base$consensus)),
                    sum(D$hit_baseline), sum(D$reportable)))
  }
}

GATE <- rbindlist(gate_rows, fill = TRUE); fwrite(GATE, file.path(DIR_PSEUDOBULK, "marker_gate.csv"))
SUMM <- rbindlist(summ, fill = TRUE);      fwrite(SUMM, file.path(DIR_PSEUDOBULK, "de_summary.csv"))
ALL  <- rbindlist(de_rows, fill = TRUE)

## PREREG 6.5 sensitivity arm: BH across the three PRIMARY bins pooled, over the intersected
## universe of those bins, reported alongside and never instead of the within-bin primary.
P <- ALL[tier == "primary"]
if (nrow(P)) {
  inter <- Reduce(intersect, split(P$gene, P$hierarchy_bin))
  S <- P[gene %in% inter][, q_pooled := p.adjust(P.Value, method = "BH")]
  S[, hit_pooled := q_pooled < Q_CUT & abs(logFC) >= LFC_CUT]
  fwrite(S[, .(tier, hierarchy_bin, gene, logFC, P.Value, q_within_bin = adj.P.Val, q_pooled,
               hit_baseline, hit_pooled)],
         file.path(DIR_PSEUDOBULK, "sensitivity_pooled_bh.csv"))
  message(sprintf("[sensitivity] pooled universe %d genes x %d bins | within-bin hits %d | pooled hits %d | disagree %d",
                  length(inter), uniqueN(S$hierarchy_bin), sum(S$hit_baseline), sum(S$hit_pooled),
                  sum(S$hit_baseline != S$hit_pooled)))
}

## The frozen hit list. Validation (step 7 of PREREG 9) may see THIS and nothing else.
FROZEN <- ALL[reportable == TRUE, .(tier, hierarchy_bin, gene, logFC_discovery = logFC,
                                    q_discovery = adj.P.Val, sign_discovery = sign(logFC))]
fwrite(FROZEN, file.path(DIR_PSEUDOBULK, "hitlist_frozen.csv"))

message("\n================ PREREG 11 summary ================")
print(SUMM)
message(sprintf("\nfrozen reportable hits: %d (primary %d, secondary %d)",
                nrow(FROZEN), nrow(FROZEN[tier == "primary"]), nrow(FROZEN[tier == "secondary"])))
if (!nrow(FROZEN))
  message("Zero reportable hits. PREREG 2 named this the expected outcome. Whether it is a null or\n",
          "an underpowered test is decided by 03_permute.R's planted-effect control, NOT here.")
