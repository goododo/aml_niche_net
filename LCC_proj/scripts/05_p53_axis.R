# 05_p53_axis.R ----
# Validate, then define, the continuous p53-activity axis that replaces the failed CNV grouping.
#
# THE FALSIFIABLE PREDICTION THIS SCRIPT EXISTS TO TEST:
#   samples independently known to have lost TP53 function must sit at the LOW end of a p53
#   transcriptional-output score. The DNA-level and RNA-level readouts share no input, so the test
#   is real. If the anchors do not rank low, the axis is not measuring p53 and 06/07 must not use it.
#
# TWO ANCHOR SETS, deliberately kept separate:
#   A4 = the 4 samples called 17p-loss by 02 (CNV evidence, cohort-internal)
#   A5 = A4 + Petti2019 809653, genotype-confirmed TP53 E286G. 809653 is NOT called by the 17p track
#        (its karyotype is del(17q), no 17p loss) and was NOT used when choosing the gene set, so it
#        is a genuinely held-out anchor rather than a fifth of the same kind.
#
# WHAT HAPPENED, IN ORDER (for the methods section -- do not present this as one clean test):
#   1. HALLMARK_P53_PATHWAY (200 genes) scored per cell with UCell -> A4 percentiles 6/38/87/66,
#      permutation p = 0.47. FAIL.
#   2. Switched to 20 canonical direct p53 targets -> A4 percentiles 96/1/12/6, p = 0.067.
#      This switch is analytic flexibility and the p-value must be reported as such.
#   3. Added the held-out genotype anchor 809653 -> percentile 8, A5 p = 0.0185. This one is clean.
#
# INPUT  : LCC_PB_DIR/<ds>/<sample>__pseudobulk.rds   (03 output; raw summed counts)
#          LCC_TAB_DIR/04_pathway_sample.csv          (the HALLMARK comparator)
#          LCC_TAB_DIR/02_sample_cnv_proxy.csv        (17p anchors + primary-set flag)
# OUTPUT : LCC_TAB_DIR/05_p53_axis.csv                per-sample axis (the downstream grouping var)
#          LCC_TAB_DIR/05_anchor_validation.csv       both anchor sets, both scores, verdicts
# Usage  : Rscript LCC_proj/scripts/05_p53_axis.R
#
# WHY within-study standardisation: the 13 studies differ in protocol (3'/5', v2/v3) and blast
# content, both of which shift a whole study's score distribution. The axis is z-scored WITHIN
# study, and study also stays in the downstream model. [mirrors blueprint D2 platform-invariance]

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))
set.seed(SEED)
N_PERM <- 20000L

## -- Step 1. score every sample from pseudobulk ----
message("[1] scoring canonical p53 targets from pseudobulk")
man <- load_sample_meta()   # timepoint-CORRECTED; see config_lcc.R::load_sample_meta
cnv <- fread(file.path(LCC_TAB_DIR, "02_sample_cnv_proxy.csv"))

score_one <- function(counts, genes) {
  if (is.null(counts) || !sum(counts)) return(list(NA_real_, 0L))
  g <- intersect(LCC_P53_TARGETS, genes)
  if (length(g) < LCC_P53_MIN_TARGETS) return(list(NA_real_, length(g)))
  cpm <- log2(counts / sum(counts) * 1e6 + 1)
  list(mean(cpm[match(g, genes)]), length(g))
}
ax <- rbindlist(lapply(seq_len(nrow(man)), function(i) {
  f <- file.path(LCC_PB_DIR, man$dataset[i], paste0(man$sample[i], "__pseudobulk.rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  a <- score_one(x$all, x$genes); m <- score_one(x$malignant, x$genes)
  data.table(dataset = man$dataset[i], sample = man$sample[i],
             p53_target_all = a[[1]], p53_target_malignant = m[[1]], n_targets_found = a[[2]],
             n_cells_all = if (is.null(x$all)) 0L else as.integer(sum(x$n_cells["all"])),
             n_cells_malignant = as.integer(x$n_cells["malignant"]))
}))
ax[man, `:=`(timepoint = i.timepoint, uid_patient = i.uid_patient, study = i.study,
             in_selection = i.in_selection), on = c("dataset", "sample")]
ax <- ax[in_selection == TRUE]   # the dataset selection is locked in config_lcc.R; 06 is the evidence
ax[cnv, `:=`(cnv_proxy_class = i.cnv_proxy_class, p17_pos = i.p17_pos,
             ck_max_arms = i.ck_max_arms, in_primary_set = i.in_primary_set),
   on = c("dataset", "sample")]
ax[is.na(p17_pos), p17_pos := FALSE][is.na(in_primary_set), in_primary_set := FALSE]

# the HALLMARK score kept alongside as the documented failed comparator
hall <- fread(file.path(LCC_TAB_DIR, "04_pathway_sample.csv"))[unit == "sample" & group == "all",
              .(dataset, sample, p53_hallmark = HALLMARK_P53_PATHWAY_UCell)]
ax[hall, p53_hallmark := i.p53_hallmark, on = c("dataset", "sample")]
# Every score is validated on BOTH its raw and its within-study-standardised form. The z form is
# the declared analysis variable (see header); the raw form is reported so the reader can see how
# much of the result rests on the standardisation rather than on the biology.
ax[!is.na(p53_hallmark), p53_hallmark_z := as.numeric(scale(p53_hallmark)), by = study]
ax[!is.na(p53_target_malignant), p53_target_malignant_z := as.numeric(scale(p53_target_malignant)), by = study]

# the analysis variable: within-study z of the canonical-target score, sign flipped so that HIGH =
# more p53-deficient. Downstream reads better with the hypothesised direction pointing up.
ax[!is.na(p53_target_all), p53_z := as.numeric(scale(p53_target_all)), by = study]
ax[, p53_deficiency_axis := -p53_z]
ax[, anchor_type := fifelse(sample == LCC_GENOTYPE_ANCHOR, "genotype_TP53_E286G",
                     fifelse(p17_pos, "cnv_17p_loss", NA_character_))]
fwrite_safe(ax, file.path(LCC_TAB_DIR, "05_p53_axis.csv"))
message("    scored ", sum(!is.na(ax$p53_target_all)), "/", nrow(ax), " samples")

## -- Step 2. anchor validation ----
message("[2] anchor validation (exact permutation on the anchor rank sum)")
# With 4-5 anchors a rank test is the only honest option; a t-test would assume a shape we cannot
# check. One-sided: the prediction is specifically LOW p53 output, decided before the test.
anchor_test <- function(d, value_col, anchor_col, label) {
  d <- d[!is.na(get(value_col))]
  is_a <- d[[anchor_col]]; k <- sum(is_a); n <- nrow(d)
  if (k < 2L || k >= n) return(NULL)
  r <- rank(d[[value_col]])
  obs <- sum(r[is_a])
  perm <- replicate(N_PERM, sum(sample(r, k)))
  data.table(anchor_set = label, score = value_col, n_total = n, n_anchor = k,
             anchor_percentiles = paste(sort(round(100 * r[is_a] / n)), collapse = ","),
             rank_sum = obs, rank_sum_expected = k * (n + 1) / 2,
             p_perm_one_sided_low = (sum(perm <= obs) + 1) / (N_PERM + 1))
}
prim <- ax[in_primary_set == TRUE]
prim[, `:=`(anchor4 = p17_pos, anchor5 = p17_pos | sample == LCC_GENOTYPE_ANCHOR)]
# THE STANDARDISATION REFERENCE IS NOT INNOCENT, AND THE RESULT IS NOT ROBUST TO IT.
# p53_z centres each study on ALL its samples (healthy, relapse, post-treatment included);
# p53_z_primary centres each study on the primary set only. Both are defensible. Petti2019 is
# 4/9 healthy donors, so the choice moves the genotype anchor 809653 from the 56th percentile
# to the 8th, and the A5 p-value from ~0.10 to ~0.02. Both rows are reported; neither is "the"
# answer, and a conclusion that depends on which one is chosen is not a robust conclusion.
prim[!is.na(p53_target_all), p53_z_primary := as.numeric(scale(p53_target_all)), by = study]
SCORES <- c("p53_hallmark", "p53_hallmark_z",           # pre-specified set, raw and standardised
            "p53_target_all", "p53_z", "p53_z_primary", # canonical set; p53_z IS the declared variable
            "p53_target_malignant", "p53_target_malignant_z")
ver <- rbindlist(lapply(SCORES, function(s) rbind(
  anchor_test(prim, s, "anchor4", "A4 (17p-loss only)"),
  anchor_test(prim, s, "anchor5", "A5 (+ held-out genotype anchor)"))), fill = TRUE)
ver[, verdict := fifelse(p_perm_one_sided_low <= 0.05, "PASS",
                  fifelse(p_perm_one_sided_low <= 0.20, "WEAK (direction right, underpowered)", "FAIL"))]
ver[, is_declared_analysis_var := score == "p53_z"]
ver[, note := fifelse(grepl("hallmark", score), "pre-specified set; failed",
               fifelse(grepl("A4", anchor_set), "gene set chosen after Hallmark failed -- flexible",
                       "held-out anchor; the clean test"))]
ver[grepl("^p53_z", score), note := paste0(note, " | z reference: ",
       fifelse(score == "p53_z", "whole study", "primary set only"))]
setorder(ver, -is_declared_analysis_var, anchor_set, score)
fwrite_safe(ver, file.path(LCC_TAB_DIR, "05_anchor_validation.csv"))
print(ver[, .(anchor_set, score, anchor_percentiles, p = round(p_perm_one_sided_low, 4), verdict)])

## -- Step 3. axis usability for the downstream models ----
message("[3] axis availability by timepoint")
print(ax[!is.na(p53_deficiency_axis), .(n = .N, n_patients = uniqueN(uid_patient)), by = timepoint][order(-n)])
message("    primary set (Diagnosis, one per patient) with a usable axis: ",
        prim[!is.na(p53_deficiency_axis), .N])
message("[done] 05_p53_axis")
