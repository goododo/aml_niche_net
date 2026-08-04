#!/usr/bin/env Rscript
# 11_compare_groups.R ----
# The comparison the PI asked for: TP53-mut vs TP53-WT AML. Tasks 2, 3 and 4 of the brief.
#
# DESIGN, as instructed: no pairing by patient, no timepoint stratification -- every AML sample
# labelled TP53-mut goes in one arm, TP53-WT in the other, and CELLS ARE POOLED across samples.
#
# THE WT ARM IS DOWN-SELECTED TO MATCH, and the selection rule is the important part of this script.
# 9 mut samples against 96 WT is not just unbalanced, it is confounded: the mut samples come from
# 4 studies and 3 timepoints, the WT pool from 10 studies. Taking 9 WT samples AT RANDOM would keep
# that confound while throwing away data. Instead each mut sample is matched 1:1 to a WT sample from
# the SAME DATASET and the SAME TIMEPOINT, nearest in cell count, forced to a distinct patient.
# The pool supports this exactly (GSE116256 Dg 3>=2, Post-tx 4>=3; Petti2019 Dg 4>=1;
# GSE227903 Relapse 7>=2; GSE289435 Dg 10>=1), so study and timepoint are eliminated by design
# rather than modelled. Seeded, therefore reproducible.
#
# TWO UNITS OF ANALYSIS, BOTH REPORTED, because they answer different questions:
#   CELL level   -- cells pooled, as instructed. n is ~10^5, so p-values are driven by cell count,
#                   not patient count: this is PSEUDOREPLICATION and the p-values must be read as
#                   descriptive. What IS trustworthy here is the effect size (odds ratio, delta pct).
#   SAMPLE level -- the 9 matched pairs, exact paired signed-rank. n = 9 pairs, so the smallest
#                   attainable one-sided p is 2^-9 = 0.002. THIS is the inferential test.
#   A result that is large at cell level but absent at sample level is a one-or-two-sample effect.
#
# INPUT  : LCC_TAB_DIR/09_tp53_groups.csv            grouping
#          LCC_TAB_DIR/04_detection_by_sample.csv    per gene x stratum cell counts (exact, all cells)
#          LCC_PERCELL_DIR/<ds>/<sample>__lcc_percell.csv.gz   per-cell UCell scores
#          LCC_TAB_DIR/04_myeloid_gate.csv           macrophage calls
#          LCC_TAB_DIR/10_stromal_denovo_sample.csv  de-novo niche composition
#          LCC_BMM_DIR/<ds>/<sample>__stemness_percell.csv     LSC17 / HSPC_core stemness
# OUTPUT : LCC_TAB_DIR/11_matched_design.csv         the 9 pairs and why each WT was chosen
#          LCC_TAB_DIR/11_gene_results.csv           per gene x stratum, cell- AND sample-level
#          LCC_TAB_DIR/11_pathway_results.csv        per pathway x stratum, both levels
#          LCC_TAB_DIR/11_composition.csv            macrophage / stromal / stemness by arm
#          LCC_TAB_DIR/11_nectin_summary.csv         task 4, as a detectability statement
# Usage  : Rscript LCC_proj/scripts/11_compare_groups.R

suppressPackageStartupMessages({
  library(data.table); library(here); library(optparse)
})
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
set.seed(SEED)

MIN_CELLS <- 30L
STRATA    <- c("all", "malignant", "HSC_MPP", "Mono_DC")

## -- Step 1. arms, and the matched WT selection ----
message("[1] building the matched design")
grp <- fread(file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))
det <- fread(file.path(LCC_TAB_DIR, "04_detection_by_sample.csv"))
ncell <- unique(det[stratum == "all", .(dataset, sample, n_cells_sample = n_cells)])

aml <- grp[timepoint != "Healthy" & !is.na(uid_patient) & uid_patient != ""]
aml <- merge(aml, ncell, by = c("dataset", "sample"), all.x = TRUE)
# tier C (expression-CNV only) is in NEITHER arm. Measured at 1/3 sensitivity against genotype
# truth, so it would inject mislabelled samples into a 9-sample group.
mut  <- aml[tp53_tier %in% c("A_genotype", "B_allele_loh")]
wtp  <- aml[tp53_tier %in% c("WT_genotyped", "WT_presumed")]
message("    TP53-mut: ", nrow(mut), " samples / ", uniqueN(mut$uid_patient), " patients")
message("    TP53-WT pool: ", nrow(wtp), " samples / ", uniqueN(wtp$uid_patient), " patients")

# Matching order: scarcest stratum first, so an abundant stratum cannot consume a patient that a
# scarce one needs (matching in sample order left GSE116256 Post-treatment unmatchable).
mut[, pool_n := vapply(seq_len(.N), function(i)
  nrow(wtp[dataset == mut$dataset[i] & timepoint == mut$timepoint[i]]), 0L)]
setorder(mut, pool_n, dataset, timepoint, sample)
# Re-use rule, deliberately SYMMETRIC with the mut arm: the mut arm is 9 samples from 6 patients
# because AML328 contributes 4 timepoints, so demanding 9 distinct WT patients would impose a
# constraint the mut arm itself does not satisfy. A WT patient may therefore appear at more than one
# timepoint, but never twice within the same timepoint, and never the same sample twice.
used_ps <- character(0); used_s <- character(0); pairs <- list()
for (i in seq_len(nrow(mut))) {
  m <- mut[i]
  free <- function(d) d[!(sample %in% used_s) & !(paste(uid_patient, timepoint) %in% used_ps)]
  cand <- free(wtp[dataset == m$dataset & timepoint == m$timepoint])
  strat <- "dataset+timepoint"
  if (!nrow(cand)) {                                # documented fallback, never a silent one
    cand <- free(wtp[dataset == m$dataset]); strat <- "dataset only" }
  if (!nrow(cand)) { strat <- "UNMATCHED"; pairs[[length(pairs) + 1L]] <-
      cbind(m[, .(mut_sample = sample, dataset, timepoint, mut_patient = uid_patient,
                  mut_n_cells = n_cells_sample)],
            data.table(wt_sample = NA_character_, wt_patient = NA_character_,
                       wt_n_cells = NA_integer_, match_strata = strat)); next }
  cand[, d := abs(log2(n_cells_sample + 1) - log2(m$n_cells_sample + 1))]
  setorder(cand, d, sample)
  p <- cand[1]
  used_s <- c(used_s, p$sample); used_ps <- c(used_ps, paste(p$uid_patient, p$timepoint))
  pairs[[length(pairs) + 1L]] <- cbind(
    m[, .(mut_sample = sample, dataset, timepoint, mut_patient = uid_patient, mut_n_cells = n_cells_sample)],
    p[, .(wt_sample = sample, wt_patient = uid_patient, wt_n_cells = n_cells_sample)],
    data.table(match_strata = strat))
}
design <- rbindlist(pairs, fill = TRUE)
design[, pair_id := .I]
fwrite_safe(design, file.path(LCC_TAB_DIR, "11_matched_design.csv"))
print(design[, .(pair_id, dataset, timepoint, mut_sample, mut_n_cells, wt_sample, wt_n_cells, match_strata)])
if (any(design$match_strata == "UNMATCHED"))
  message("    [!] ", sum(design$match_strata == "UNMATCHED"), " mut samples had no WT partner")

arms <- rbind(design[, .(pair_id, dataset, sample = mut_sample, arm = "TP53-mut")],
              design[!is.na(wt_sample), .(pair_id, dataset, sample = wt_sample, arm = "TP53-WT")])

## -- Step 2. gene panel: cell-level pooled + sample-level paired ----
message("[2] gene panel")
panel <- fread_commented(LCC_GENE_PANEL_TSV)
# NECTIN4/PVRL4 symbol drift: 190 samples carry NECTIN4, 39 carry PVRL4, never both. Collapse the
# aliases so a sample is never scored zero merely for using the other symbol.
alias <- c(PVRL4 = "NECTIN4", PVRL1 = "NECTIN1", PVRL2 = "NECTIN2", PVRL3 = "NECTIN3")
g <- merge(det[stratum %in% STRATA & n_cells >= MIN_CELLS], arms, by = c("dataset", "sample"))
g[gene %in% names(alias), gene := alias[gene]]
# 03 emits a row for EVERY panel gene in EVERY sample, with gene_present = FALSE and NA counts when
# the gene is absent from that sample's feature space (4,590 such rows). A gene that a study never
# measured is MISSING DATA, not a zero -- scoring it zero would manufacture a group difference out
# of a platform difference. Such sample x gene cells are dropped, so n_cells varies by gene and is
# reported per gene rather than assumed constant.
g <- g[, .(n_cells = n_cells[1],
           n_nonzero    = if (all(is.na(n_nonzero)))    NA_integer_ else max(n_nonzero, na.rm = TRUE),
           mean_lognorm = if (all(is.na(mean_lognorm))) NA_real_    else max(mean_lognorm, na.rm = TRUE),
           gene_present = any(gene_present)),
       by = .(dataset, sample, pair_id, arm, stratum, gene)]
g <- g[gene_present == TRUE & !is.na(n_nonzero)]

# CELL LEVEL: exact pooled counts over every cell, no subsampling. The 2x2 is
# (expressing / not) x (mut / WT); the statistic is a Woolf log-odds-ratio with a Haldane 0.5
# correction so zero cells do not produce an infinite estimate.
cell <- g[, .(n_samples = .N, n_cells = sum(n_cells), n_pos = sum(n_nonzero),
              mean_lognorm = sum(mean_lognorm * n_cells) / sum(n_cells)), by = .(arm, stratum, gene)]
cw <- dcast(cell, stratum + gene ~ arm, value.var = c("n_samples", "n_cells", "n_pos", "mean_lognorm"))
setnames(cw, gsub("TP53-", "", names(cw)))
cw <- cw[!is.na(n_cells_mut) & !is.na(n_cells_WT)]
cw[, `:=`(pct_mut = 100 * n_pos_mut / n_cells_mut, pct_WT = 100 * n_pos_WT / n_cells_WT)]
cw[, `:=`(a = n_pos_mut + .5, b = n_cells_mut - n_pos_mut + .5,
          c = n_pos_WT  + .5, d = n_cells_WT  - n_pos_WT  + .5)]
cw[, `:=`(log_or = log((a / b) / (c / d)), se = sqrt(1/a + 1/b + 1/c + 1/d))]
cw[, `:=`(z_cell = log_or / se, odds_ratio = exp(log_or))]
cw[, p_cell_two_sided := 2 * pnorm(-abs(z_cell))]
cw[, `:=`(a = NULL, b = NULL, c = NULL, d = NULL)]

# SAMPLE LEVEL: exact paired signed-rank over the 9 pairs. This is the inferential test.
paired_test <- function(dm, dw) {
  ok <- !is.na(dm) & !is.na(dw); dm <- dm[ok]; dw <- dw[ok]
  n <- length(dm); if (n < 4L) return(list(NA_real_, NA_real_, n, NA_real_))
  d <- dm - dw; d <- d[d != 0]; k <- length(d)
  if (!k) return(list(0, 1, n, 0))
  r <- rank(abs(d)); wplus <- sum(r[d > 0])
  # exact null of the signed-rank statistic by DP over all 2^k sign assignments
  maxs <- k * (k + 1L) / 2L
  ways <- numeric(maxs + 1L); ways[1L] <- 1
  for (i in seq_len(k)) { sh <- c(numeric(r[i]), ways[seq_len(maxs + 1L - r[i])]); ways <- ways + sh }
  tot <- sum(ways)
  list(median(d), sum(ways[seq(floor(wplus) + 1L, maxs + 1L)]) / tot, n, wplus)
}
sw <- dcast(g, stratum + gene + pair_id ~ arm, value.var = "mean_lognorm")
setnames(sw, c("TP53-mut", "TP53-WT"), c("v_mut", "v_wt"))
# frac_pairs_mut_higher is the effect size that MATCHES the test. The cell-level odds ratio does
# not: it is cell-weighted, so one 15,000-cell sample can set the pooled fraction on its own, and
# the two measures disagree in SIGN for 43 of 145 panel genes (VSIG4: OR 0.23 pooled, yet mut is
# higher in most pairs). Any figure must plot the sample-level effect against the sample-level test.
samp <- sw[, { r <- paired_test(v_mut, v_wt)
               ok <- !is.na(v_mut) & !is.na(v_wt)
               .(n_pairs = r[[3]], median_delta_log2cpm = r[[1]], p_sample_higher = r[[2]],
                 frac_pairs_mut_higher = if (sum(ok)) mean(v_mut[ok] > v_wt[ok]) else NA_real_) },
           by = .(stratum, gene)]

gene_res <- merge(cw, samp, by = c("stratum", "gene"), all = TRUE)
gene_res[panel, category := i.category, on = "gene"]
gene_res[is.na(category), category := "other"]
# BH within (stratum x category): the panel is pre-specified and organised by hypothesis, so each
# category is its own family. Correcting across all ~170 genes at once is the wrong family.
gene_res[, fdr_sample := p.adjust(p_sample_higher, "BH"), by = .(stratum, category)]
setorder(gene_res, stratum, p_sample_higher)
fwrite_safe(gene_res, file.path(LCC_TAB_DIR, "11_gene_results.csv"))

## -- Step 3. pathway scores, same two units ----
message("[3] pathway scores")
pws <- rbindlist(lapply(seq_len(nrow(arms)), function(i) {
  f <- file.path(LCC_PERCELL_DIR, arms$dataset[i], paste0(arms$sample[i], "__lcc_percell.csv.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  sc <- grep("_UCell$", names(d), value = TRUE); if (!length(sc)) return(NULL)
  keep <- intersect(c("hierarchy_bin", "malignant", sc), names(d))
  d <- d[, ..keep]
  d[, `:=`(dataset = arms$dataset[i], sample = arms$sample[i],
           pair_id = arms$pair_id[i], arm = arms$arm[i])]
  d[]
}), fill = TRUE)
sc_cols <- grep("_UCell$", names(pws), value = TRUE)
pws[, stratum_all := "all"]
long <- rbindlist(list(
  melt(pws, id.vars = c("pair_id", "arm", "sample", "stratum_all"), measure.vars = sc_cols,
       variable.name = "pathway", value.name = "score")[, .(pair_id, arm, sample,
                                                            stratum = stratum_all, pathway, score)],
  melt(pws[hierarchy_bin %in% STRATA], id.vars = c("pair_id", "arm", "sample", "hierarchy_bin"),
       measure.vars = sc_cols, variable.name = "pathway", value.name = "score")[
         , .(pair_id, arm, sample, stratum = hierarchy_bin, pathway, score)],
  melt(pws[malignant == 1L], id.vars = c("pair_id", "arm", "sample"), measure.vars = sc_cols,
       variable.name = "pathway", value.name = "score")[, .(pair_id, arm, sample,
                                                            stratum = "malignant", pathway, score)]
), fill = TRUE)[!is.na(score)]

# CELL level: rank-sum over pooled cells. Reported with its n so the pseudoreplication is visible.
cellp <- long[, { m <- score[arm == "TP53-mut"]; w <- score[arm == "TP53-WT"]
                  if (length(m) < 30 || length(w) < 30) .(n_cells_mut = length(m), n_cells_WT = length(w),
                       auc = NA_real_, p_cell = NA_real_, mean_mut = NA_real_, mean_WT = NA_real_)
                  else { r <- rank(c(m, w)); U <- sum(r[seq_along(m)]) - length(m) * (length(m) + 1) / 2
                    a <- U / (length(m) * length(w))
                    s <- sqrt(length(m) * length(w) * (length(m) + length(w) + 1) / 12)
                    .(n_cells_mut = length(m), n_cells_WT = length(w), auc = a,
                      p_cell = 2 * pnorm(-abs((U - length(m)*length(w)/2) / s)),
                      mean_mut = mean(m), mean_WT = mean(w)) } }, by = .(stratum, pathway)]
# SAMPLE level: per-sample mean, then the same exact paired signed-rank over the 9 pairs.
sp <- long[, .(v = mean(score)), by = .(stratum, pathway, pair_id, arm)]
spw <- dcast(sp, stratum + pathway + pair_id ~ arm, value.var = "v")
setnames(spw, c("TP53-mut", "TP53-WT"), c("v_mut", "v_wt"))
sampp <- spw[, { r <- paired_test(v_mut, v_wt)
                 .(n_pairs = r[[3]], median_delta = r[[1]], p_sample_higher = r[[2]]) },
             by = .(stratum, pathway)]
pw_res <- merge(cellp, sampp, by = c("stratum", "pathway"), all = TRUE)
pw_res[, fdr_sample := p.adjust(p_sample_higher, "BH"), by = stratum]
setorder(pw_res, stratum, p_sample_higher)
fwrite_safe(pw_res, file.path(LCC_TAB_DIR, "11_pathway_results.csv"))

## -- Step 4. cell-type composition by arm: macrophage / stroma / stemness ----
message("[4] composition: macrophage, stromal niche, stemness")
comp <- copy(arms)
my <- fread(file.path(LCC_TAB_DIR, "04_myeloid_gate.csv"))
comp[my, `:=`(n_cells_total = i.n_cells_total, n_mono_dc = i.n_mono_dc,
              n_macrophage_like = i.n_macrophage_like), on = c("dataset", "sample")]
f10 <- file.path(LCC_TAB_DIR, "10_stromal_denovo_sample.csv")
if (file.exists(f10)) {
  st <- fread(f10)
  cols <- intersect(c("n_stromal_denovo", "n_MSC_fibroblast", "n_adipocyte", "n_MSC_adipo_primed",
                      "n_endothelial", "n_pericyte", "n_osteolineage"), names(st))
  comp[st, (cols) := mget(paste0("i.", cols)), on = c("dataset", "sample")]
  for (j in cols) comp[is.na(get(j)), (j) := 0L]   # screened-out samples genuinely have none
}
# stemness: per-cell LSC17 / HSPC_core from the main line, averaged over the HSC_MPP bin
stem <- rbindlist(lapply(seq_len(nrow(comp)), function(i) {
  f <- file.path(LCC_BMM_DIR, comp$dataset[i], paste0(comp$sample[i], "__stemness_percell.csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  sc <- intersect(c("LSC17", "HSPC_core", "vanGalen_HSC_like", "vanGalen_HSC_Prog"), names(d))
  if (!length(sc)) return(NULL)
  h <- if ("hierarchy_bin" %in% names(d)) d[hierarchy_bin == "HSC_MPP"] else d
  if (!nrow(h)) h <- d
  cbind(data.table(dataset = comp$dataset[i], sample = comp$sample[i], n_hsc = nrow(h)),
        h[, lapply(.SD, mean, na.rm = TRUE), .SDcols = sc])
}), fill = TRUE)
if (nrow(stem)) comp <- merge(comp, stem, by = c("dataset", "sample"), all.x = TRUE)
comp[, `:=`(pct_macrophage = 100 * n_macrophage_like / n_cells_total,
            pct_stromal    = if ("n_stromal_denovo" %in% names(comp)) 100 * n_stromal_denovo / n_cells_total else NA_real_)]
fwrite_safe(comp, file.path(LCC_TAB_DIR, "11_composition.csv"))
CVARS <- intersect(c("pct_macrophage", "pct_stromal", "LSC17", "HSPC_core",
                     "vanGalen_HSC_like", "vanGalen_HSC_Prog"), names(comp))
cres <- rbindlist(lapply(CVARS, function(v) {
  w <- dcast(comp[, .(pair_id, arm, v = get(v))], pair_id ~ arm, value.var = "v")
  if (!all(c("TP53-mut", "TP53-WT") %in% names(w))) return(NULL)
  r <- paired_test(w$`TP53-mut`, w$`TP53-WT`)
  data.table(variable = v, n_pairs = r[[3]], median_delta = r[[1]], p_sample_higher = r[[2]],
             mean_mut = mean(comp[arm == "TP53-mut"][[v]], na.rm = TRUE),
             mean_WT  = mean(comp[arm == "TP53-WT"][[v]],  na.rm = TRUE))
}), fill = TRUE)
cat("\n-- composition, matched pairs --\n"); print(cres)

## -- Step 5. task 4: NECTIN4 as a detectability statement ----
message("[5] NECTIN4 / PVRL4")
# The PI's Nectin-4 finding is IHC/protein. The question here is NOT "is Nectin-4 higher in TP53
# AML" but "can 10x 3' marrow scRNA-seq see NECTIN4 at all". At the ambient floor a group difference
# is uninterpretable and must not be reported as one.
nec <- gene_res[stratum == "all" & gene %in% c("NECTIN4", "NECTIN1", "NECTIN2", "NECTIN3", "PVR", "TIGIT", "CD226"),
                .(gene, n_cells_mut, n_pos_mut, pct_mut, n_cells_WT, n_pos_WT, pct_WT,
                  odds_ratio, p_cell_two_sided, p_sample_higher)]
nec[, interpretable := pct_mut >= 1 | pct_WT >= 1]
fwrite_safe(nec, file.path(LCC_TAB_DIR, "11_nectin_summary.csv"))
print(nec[order(-pct_WT)])

message("[6] summary")
for (s in c("all", "malignant")) {
  cat("\n== stratum:", s, "-- genes ranked by the SAMPLE-level paired test ==\n")
  print(gene_res[stratum == s][order(p_sample_higher)][1:12,
        .(gene, category, pct_mut = round(pct_mut, 2), pct_WT = round(pct_WT, 2),
          OR = round(odds_ratio, 2), d_log2cpm = round(median_delta_log2cpm, 3),
          n_pairs, p_samp = signif(p_sample_higher, 3), fdr = signif(fdr_sample, 3))])
}
cat("\n== pathways, stratum=all ==\n")
print(pw_res[stratum == "all"][order(p_sample_higher),
      .(pathway = sub("_UCell$", "", pathway), auc_cell = round(auc, 3),
        n_cells_mut, n_cells_WT, d = round(median_delta, 4), n_pairs,
        p_samp = signif(p_sample_higher, 3), fdr = signif(fdr_sample, 3))])
cat("\n== anything passing FDR 0.10 at the SAMPLE level ==\n")
print(gene_res[fdr_sample <= 0.10, .(stratum, gene, category, OR = round(odds_ratio, 2),
                                     p = signif(p_sample_higher, 3), fdr = signif(fdr_sample, 3))][order(fdr)])
message("[done] 11_compare_groups")
