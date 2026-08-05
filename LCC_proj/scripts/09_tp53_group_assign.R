# 09_tp53_group_assign.R ----
# Assemble the final TP53 group assignment from every available evidence arm, with an explicit
# confidence tier per sample. This is the table every downstream analysis groups on.
#
# EVIDENCE ARMS, ranked by how directly they measure TP53 inactivation:
#   T1  genotype        07: per-cell / reported mutation calls. The only arm that measures the
#                       MUTATION. Names the variant.
#   T2  allele LOH      08: Numbat 17p LOH from phased allele counts. Measures loss of the second
#                       allele -- the canonical partner lesion. ORTHOGONAL to expression.
#   T3  expression CNV  02: inferCNV 17p gene-level loss (+ CopyKAT when installed). Measured
#                       performance against T1 truth: 0/2 sensitivity (it missed AML916 and 809653) plus 1 false positive
#                       (21 genotype-confirmed mutant cells) and 809653. LOW trust, used only to
#                       flag candidates, never alone to assert TP53 status.
#
# TIERS (the point of this script -- a tier is a claim about how much the label can carry):
#   A_genotype       T1 present                      -> may be called "TP53-mutant"
#   B_allele_loh     T2 present, no T1               -> called "TP53-aberrant (17p LOH)"
#   C_expr_only      T3 only                         -> called "17p-loss candidate"; NOT TP53-mutant
#   WT_genotyped     T1 says wildtype                -> confident wild-type
#   WT_presumed      no arm positive, arms available -> presumed wild-type (absence of evidence)
#   unevaluable      no arm covers this sample
# [CODING_STANDARDS §9] Only tiers A (and B with the "aberrant/LOH" wording) may appear in a TP53
# claim. Tier C must be reported as a CNV finding.
#
# INPUT  : LCC_TAB_DIR/07_tp53_genotype_sample.csv, 08_numbat_17p.csv,
#          02_sample_cnv_proxy_all_datasets.csv, 04_detection_by_sample.csv (stromal counts)
# OUTPUT : LCC_TAB_DIR/09_tp53_groups.csv          per-sample group + tier + evidence detail
#          LCC_TAB_DIR/09_tp53_groups_patient.csv  per-patient rollup (a mutation at ANY timepoint
#                                                  makes the patient TP53-aberrant)
#          LCC_TAB_DIR/09_analysis_sets.csv        which comparison each sample can serve in
# Usage  : Rscript LCC_proj/scripts/09_tp53_group_assign.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("LCC_proj", "scripts", "config_lcc.R"))

message("[1] loading evidence arms")
man <- load_sample_meta(include_stroma_ref = TRUE)
g <- data.table(); nb <- data.table(); cn <- data.table()
f <- file.path(LCC_TAB_DIR, "07_tp53_genotype_sample.csv")
if (file.exists(f)) g <- fread(f)[, .(dataset, sample, geno_status = tp53_status,
                                      geno_variants = tp53_variants,
                                      n_cells_tp53_mut, n_cells_tp53_wt, geno_evidence = evidence)]
f <- file.path(LCC_TAB_DIR, "08_numbat_17p.csv")
if (file.exists(f)) nb <- fread(f)[, .(dataset, sample, numbat_17p_loh, arm17p_frac_loh,
                                       numbat_tp53_seg = tp53_seg_state, numbat_llr = best_llr)]
f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy_all_datasets.csv")
if (!file.exists(f)) f <- file.path(LCC_TAB_DIR, "02_sample_cnv_proxy.csv")
if (file.exists(f)) cn <- fread(f)[, .(dataset, sample, expr_17p_loss = p17_pos, expr_ck_arms = ck_max_arms)]

out <- man[, .(dataset, sample, timepoint, uid_patient, study, in_selection, dataset_role)]
for (x in list(g, nb, cn)) if (nrow(x)) out <- merge(out, x, by = c("dataset", "sample"), all.x = TRUE)
for (j in c("numbat_17p_loh", "expr_17p_loss")) if (j %in% names(out)) out[is.na(get(j)), (j) := FALSE]

## -- stromal capture, for the niche/MSC layer ----
f <- file.path(LCC_TAB_DIR, "04_detection_by_sample.csv")
if (file.exists(f)) {
  st <- fread(f)[gene == "COL1A1" & stratum == "Stromal", .(dataset, sample, n_stromal = n_cells)]
  out <- merge(out, st, by = c("dataset", "sample"), all.x = TRUE)
}
out[is.na(n_stromal), n_stromal := 0L]

message("[2] assigning tiers")
out[, arms_available := (!is.na(geno_status)) + (sample %in% nb$sample) + (sample %in% cn$sample)]
out[, tp53_tier := fifelse(!is.na(geno_status) & geno_status == "mutant", "A_genotype",
                    fifelse(numbat_17p_loh == TRUE, "B_allele_loh",
                     fifelse(!is.na(geno_status) & geno_status == "wildtype", "WT_genotyped",
                      fifelse(expr_17p_loss == TRUE, "C_expr_only",
                       fifelse(arms_available > 0, "WT_presumed", "unevaluable")))))]
# the label that may actually be written in text, per tier
out[, tp53_group := fifelse(tp53_tier == "A_genotype",   "TP53-mutant",
                     fifelse(tp53_tier == "B_allele_loh", "TP53-aberrant (17p LOH)",
                      fifelse(tp53_tier == "C_expr_only", "17p-loss candidate",
                       fifelse(tp53_tier %in% c("WT_genotyped", "WT_presumed"), "TP53-WT", NA_character_))))]
setorder(out, tp53_tier, dataset, sample)
fwrite_safe(out, file.path(LCC_TAB_DIR, "09_tp53_groups.csv"))

message("[3] patient rollup")
# A TP53 lesion found at ANY timepoint defines the patient: TP53 clones are frequently subclonal at
# diagnosis and expand under therapy (this cohort shows exactly that -- patient 1886 is 17p-neutral
# at diagnosis and carries 98% 17p LOH at relapse).
pat <- out[!is.na(uid_patient) & timepoint != "Healthy",
           .(dataset = dataset[1], n_samples = .N,
             timepoints = paste(sort(unique(timepoint)), collapse = "/"),
             best_tier = min(tp53_tier),
             variants = paste(sort(unique(na.omit(geno_variants[nzchar(geno_variants)]))), collapse = "+"),
             any_numbat_loh = any(numbat_17p_loh), any_expr_17p = any(expr_17p_loss),
             max_stromal = max(n_stromal)), by = uid_patient]
pat[, tp53_group_patient := fifelse(best_tier == "A_genotype", "TP53-mutant",
                             fifelse(best_tier == "B_allele_loh", "TP53-aberrant (17p LOH)",
                              fifelse(best_tier == "C_expr_only", "17p-loss candidate", "TP53-WT")))]
setorder(pat, best_tier, -max_stromal)
fwrite_safe(pat, file.path(LCC_TAB_DIR, "09_tp53_groups_patient.csv"))

message("[4] summary")
cat("\n-- samples by tier --\n");  print(out[timepoint != "Healthy", .N, by = .(tp53_tier, tp53_group)][order(tp53_tier)])
cat("\n-- AML patients by group --\n"); print(pat[, .N, by = .(best_tier, tp53_group_patient)][order(best_tier)])
cat("\n-- TP53-aberrant (tier A or B) patients, with stromal capture --\n")
print(pat[best_tier %in% c("A_genotype", "B_allele_loh"),
          .(uid_patient, dataset, timepoints, variants, any_numbat_loh, max_stromal)])
cat("\n-- samples usable for the STROMAL/niche layer (>=50 stromal cells, AML) --\n")
print(out[timepoint != "Healthy" & n_stromal >= 50,
          .(dataset, sample, timepoint, tp53_group, tp53_tier, n_stromal)][order(-n_stromal)])
cat("\n-- healthy negative control: any arm positive? --\n")
print(out[timepoint == "Healthy", .N, by = .(tp53_tier)][order(tp53_tier)])
message("[done] 09_tp53_group_assign")
