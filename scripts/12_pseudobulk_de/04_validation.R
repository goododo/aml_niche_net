#!/usr/bin/env Rscript
# 04_validation.R ----
# Step 7 of PREREG 9, the last one. Directional-concordance check of the FROZEN primary-tier hit
# lists on GSE116256. Runs only after Stage B and the permutation/power step have written their
# outputs; it is a separate invocation on purpose, so that nothing here can feed back into them.
#
# INPUT  : DIR_PSEUDOBULK/pb_sample_manifest.csv              (roster; Validation rows, in_analysis)
#          DIR_PSEUDOBULK/hitlist_frozen.csv                  (from 02 -- the ONLY thing this may see)
#          DIR_PSEUDOBULK/depth/<ds>__<sample>__pb_depth.csv  (per-bin cell counts, from 01)
#          PB_RDS_DIR/<ds>/<sample>__pseudobulk.rds           (genes x bins raw counts, from 01)
#          INFERCNV_GENE_ORDER                                 (gencode GRCh38; chrY identification)
# OUTPUT : DIR_PSEUDOBULK/validation_concordance.csv          one row per primary bin
#          DIR_PSEUDOBULK/validation/<bin>.csv                per-gene Discovery vs Validation signs
# Usage  : Rscript scripts/12_pseudobulk_de/04_validation.R
#
# EVERYTHING HERE IS PRE-REGISTERED. PREREGISTRATION_pseudobulk_de.md sec 3.4(c) fixes the model,
# the universe rule, the statistic, the alpha, the multiplicity decision and the empty-list
# behaviour. This script implements that; it decides nothing.
#
#   PREREG 3.4  VALIDATION IS 9 AML vs 2 HEALTHY, AND IS NOT A REPLICATION. `BM5-34p` is dropped
#               (CD34+ sorted, fresh, against an otherwise viability-only / ficoll / cryopreserved
#               arm) and the four post-treatment marrows are dropped (a regenerating marrow is a
#               different biological state, not a second observation of the same one). Stage A
#               already wrote both exclusions into the manifest as in_analysis = FALSE; SELF-CHECK 1
#               asserts the surviving roster is exactly 9 vs 2 rather than trusting that.
#
#   PREREG 3.4  MODEL IS `~ arm`. No `dataset` term -- Validation is one dataset and
#               `~ dataset + arm` is rank-deficient there. No duplicateCorrelation -- every
#               Validation library block has size 1, so there is nothing to pool.
#
#   PREREG 9.7  PRIMARY TIER ONLY. The secondary tier (LMPP_GMP / HSC_MPP / Erythroid) has NO
#               Validation arm, because PREREG 5.3 shows GSE116256 gives it 1-2 controls. That was
#               fixed before Discovery ran, so the 304 frozen LMPP_GMP hits are deliberately not
#               tested here. Not testing them is the pre-registered behaviour, not an omission.
#
#   PREREG 3.4  SIGN AGREEMENT IS THE CLAIM. No Validation q-value is produced. Decided before the
#               arm was opened so the result cannot be upgraded to "replication" if the signs agree
#               nor downgraded to "exploratory" if they do not. alpha = 0.05, one-sided binomial
#               against p = 0.5, and NO multiplicity correction across the three bins -- all three
#               p-values are reported together.
suppressPackageStartupMessages({
  library(data.table); library(here); library(edgeR); library(limma)
})
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))   # INFERCNV_GENE_ORDER

PRIMARY_BINS <- c("Mono_DC", "T_NK", "B_Plasma")     # PREREG 5.4; secondary has no Validation arm
MIN_CELLS    <- CCC_MIN_CELLS_PER_OCCUPIED_BIN       # 30, PREREG 5.1 -- same gate as Discovery
ALPHA        <- 0.05                                 # PREREG 3.4(c)

dir.create(file.path(DIR_PSEUDOBULK, "validation"), recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------------------- roster ----
MAN <- fread(file.path(DIR_PSEUDOBULK, "pb_sample_manifest.csv"))
VAL <- MAN[in_analysis == TRUE & split == "Validation"]

## -- SELF-CHECK 1: the roster PREREG 3.4 fixed, asserted rather than assumed ----
## If Stage A's exclusions ever drift, every number below would still compute and would silently
## be a different study. This is the project's recurring failure mode, so it is checked here.
n_aml <- VAL[arm == "AML", .N]; n_hty <- VAL[arm == "healthy", .N]
message(sprintf("[SELF-CHECK 1] Validation roster: %d AML vs %d healthy, %d dataset(s)",
                n_aml, n_hty, uniqueN(VAL$dataset)))
stopifnot(n_aml == 9L, n_hty == 2L, identical(unique(VAL$dataset), "GSE116256"))
stopifnot(all(VAL[arm == "AML", timepoint] == "Diagnosis"))
stopifnot(!("BM5-34p" %in% VAL$sample))
message("               9 vs 2, GSE116256 only, all AML at Diagnosis, BM5-34p absent -- as registered")

DEP <- rbindlist(lapply(list.files(file.path(DIR_PSEUDOBULK, "depth"), full.names = TRUE), fread))
PB  <- lapply(seq_len(nrow(VAL)), function(i)
  readRDS(file.path(PB_RDS_DIR, VAL$dataset[i], paste0(VAL$sample[i], "__pseudobulk.rds"))))
names(PB) <- paste(VAL$dataset, VAL$sample, sep = "|")

## Validation's own universe. It CANNOT be Discovery's: GSE116256 is a different deposit (Seq-Well)
## with a different gene set, and PREREG 6.4 removes the same two sex-linked classes from every
## universe. Genes the frozen list loses to this intersection are counted below, not dropped quietly.
base_universe <- Reduce(intersect, lapply(PB, `[[`, "genes"))
GO   <- fread(INFERCNV_GENE_ORDER, header = FALSE, col.names = c("gene", "chr", "start", "end"))
sexg <- union(GO[chr %in% c("chrY", "Y"), unique(gene)], "XIST")
VUNIV <- setdiff(base_universe, sexg)
message(sprintf("[universe] GSE116256 intersect = %d genes; %d sex genes removed; base = %d",
                length(base_universe), sum(base_universe %in% sexg), length(VUNIV)))

## Build the (genes x samples) count matrix for one bin over the Validation roster. Same shape as
## 02's bin_matrix, but on VUNIV and with no dataset/library columns -- there is one of each.
vbin_matrix <- function(bin) {
  keep <- DEP[hierarchy_bin == bin & n_cells >= MIN_CELLS,
              .(dataset, sample, n_cells, median_umi)][VAL, on = c("dataset", "sample"), nomatch = 0]
  if (!nrow(keep)) return(NULL)
  m <- vapply(seq_len(nrow(keep)), function(i) {
    PB[[paste(keep$dataset[i], keep$sample[i], sep = "|")]]$counts[VUNIV, bin]
  }, numeric(length(VUNIV)))
  dimnames(m) <- list(VUNIV, paste(keep$dataset, keep$sample, sep = "|"))
  # healthy is the REFERENCE level (PREREG 6.1). Without relevel, R sorts "AML" first and every
  # logFC sign in this script would flip -- which is exactly the quantity being tested.
  keep[, arm := relevel(factor(arm), ref = "healthy")]
  list(counts = m, meta = keep)
}

HIT <- fread(file.path(DIR_PSEUDOBULK, "hitlist_frozen.csv"))[tier == "primary"]
message(sprintf("[frozen] primary-tier hit list: %d genes over %d bins (secondary tier excluded by PREREG 9.7)",
                nrow(HIT), uniqueN(HIT$hierarchy_bin)))

## ------------------------------------------------------- the check ----
out <- rbindlist(lapply(PRIMARY_BINS, function(bin) {
  hits <- HIT[hierarchy_bin == bin]
  B <- vbin_matrix(bin)

  ## PREREG 3.4(c): an empty frozen list is recorded as not evaluable. It is NOT pooled with other
  ## bins to manufacture an n.
  if (!nrow(hits) || is.null(B) || uniqueN(B$meta$arm) < 2L) {
    reason <- if (!nrow(hits)) "no Discovery hits" else "no Validation samples clearing the cell gate"
    message(sprintf("[%-9s] not evaluable -- %s", bin, reason))
    return(data.table(hierarchy_bin = bin, evaluable = FALSE, reason = reason,
                      n_hits_frozen = nrow(hits), n_tested = NA_integer_, n_lost = NA_integer_,
                      n_same_sign = NA_integer_, frac_same_sign = NA_real_,
                      binom_p = NA_real_, n_AML = NA_integer_, n_healthy = NA_integer_,
                      min_sample_size = NA_real_, n_genes_tested = NA_integer_))
  }

  md     <- B$meta
  design <- model.matrix(~ arm, data = md)              # PREREG 3.4: no dataset, no block
  stopifnot("armAML" %in% colnames(design))             # the coefficient the sign comes from
  y <- DGEList(counts = B$counts)
  keep_g <- filterByExpr(y, design)                     # PREREG 6.4, Validation's own design
  y <- y[keep_g, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design)
  fit <- eBayes(lmFit(v, design))
  tt <- as.data.table(topTable(fit, coef = "armAML", number = Inf), keep.rownames = "gene")

  ## PREREG 3.4(c): the tested set is the frozen list INTERSECTED with what clears the filter here.
  m <- merge(hits[, .(gene, logFC_discovery, sign_discovery)], tt[, .(gene, logFC_val = logFC)],
             by = "gene", all.x = TRUE)
  tested <- m[!is.na(logFC_val)]
  n_lost <- nrow(m) - nrow(tested)
  same   <- tested[sign(logFC_val) == sign_discovery, .N]
  bt     <- if (nrow(tested)) stats::binom.test(same, nrow(tested), p = 0.5, alternative = "greater") else NULL

  ## A bin whose frozen list survives the universe intersection with ZERO genes has no test, and
  ## must be recorded as not evaluable for the same reason PREREG 3.4(c) gives for an empty frozen
  ## list: there is nothing to pool it into. Reporting it as "evaluable" with a blank p would let a
  ## reader count it as a bin that was checked.
  evaluable <- nrow(tested) > 0L

  fwrite(m[order(-abs(logFC_discovery))], file.path(DIR_PSEUDOBULK, "validation", paste0(bin, ".csv")))
  message(sprintf("[%-9s] %d AML vs %d healthy | frozen %d -> tested %d (lost %d) | same sign %d (%.1f%%) | one-sided p = %s",
                  bin, md[arm == "AML", .N], md[arm == "healthy", .N], nrow(m), nrow(tested), n_lost,
                  same, 100 * same / max(nrow(tested), 1),
                  if (is.null(bt)) "n/a" else sprintf("%.4f", bt$p.value)))

  data.table(hierarchy_bin = bin,
             evaluable = evaluable,
             reason = if (evaluable) NA_character_ else "no frozen hit survived the Validation universe",
             n_hits_frozen = nrow(m), n_tested = nrow(tested), n_lost = n_lost,
             n_same_sign = same, frac_same_sign = same / max(nrow(tested), 1),
             binom_p = if (is.null(bt)) NA_real_ else bt$p.value,
             n_AML = md[arm == "AML", .N], n_healthy = md[arm == "healthy", .N],
             min_sample_size = 1 / max(hat(design)), n_genes_tested = sum(keep_g))
}))

fwrite(out, file.path(DIR_PSEUDOBULK, "validation_concordance.csv"))

## -- SELF-CHECK 2: no multiplicity correction is applied, and that is registered, not forgotten ----
message("\n[SELF-CHECK 2] PREREG 3.4(c): alpha = 0.05, one-sided, NO correction across the three")
message("               bins. All three p-values are reported together so the reader applies their")
message("               own if they wish. No Validation q-value is produced, by design.")
ev <- out[evaluable == TRUE]
message(sprintf("\n[result] %d of %d primary bins evaluable; %d with same-sign fraction > 0.5 at p < %.2f",
                nrow(ev), nrow(out), ev[frac_same_sign > 0.5 & binom_p < ALPHA, .N], ALPHA))
message("[result] REGISTERED READING: a bin above 0.5 at one-sided p < 0.05 is reported as")
message("         'direction replicates'; anything else as 'does not replicate'. Neither may be")
message("         described as a powered per-gene validation at 9 vs 2.")
message(sprintf("\n[done] wrote validation_concordance.csv (%d rows) and %d per-gene tables",
                nrow(out), nrow(ev)))
