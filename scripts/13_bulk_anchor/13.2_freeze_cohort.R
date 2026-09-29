# 13.2_freeze_cohort.R ----
# INPUT  : /FAST/gr10634/gaozy/bulk_cohorts/beataml2/ (verified by 13.1, 20/20 checks)
# OUTPUT : results/tables/13_bulk_anchor/
#            beataml_clinical_tidy.csv   one row per specimen: ids, clinical, genotype, flags
#            beataml_drug_auc_long.csv   single-agent probit AUC/IC50 per (specimen, inhibitor)
#            cohort_freeze.csv           specimen membership of the frozen analysis sets
# WHAT IT DOES : freezes the BeatAML2 analysis cohorts BEFORE any association is computed, so
#          13.3's pre-registration can reference fixed sample lists:
#            SET_RNA        analysis-grade RNA specimens                    (671, known answer)
#            SET_RNA_DRUG   + analysis-grade ex vivo drug                   (520)
#            SET_TRIPLE     + analysis-grade WES                            (494 / 458 patients)
#            SET_TRIPLE_1PP one specimen per patient inside SET_TRIPLE      (458)
#          One-per-patient rule, fixed here: prefer initial-diagnosis stage, then earliest
#          collection time relative to inclusion, then bone marrow over blood, then the
#          lexicographically smallest dbgap_rnaseq_sample (fully deterministic).
#          Analysis-grade flags are OHSU's own QC (13.1 reconciliation note: the recon report's
#          542/514/471 were the loose rnaSeq-flag variant; both reproduce exactly).

suppressPackageStartupMessages({ library(data.table); library(readxl) })
DIR <- "/FAST/gr10634/gaozy/bulk_cohorts/beataml2"
OUT <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"

cl <- as.data.table(read_excel(file.path(DIR, "beataml_wv1to4_clinical.xlsx"), sheet = "summary"))
f <- function(x) !is.na(x) & (x == "y" | x == TRUE)

# ---- 1. tidy clinical ---------------------------------------------------------------------------
td <- cl[, .(dbgap_subject_id, dbgap_rnaseq_sample, dbgap_dnaseq_sample, cohort,
             specimenType, diseaseStageAtSpecimenCollection,
             time_rel_inclusion = suppressWarnings(as.numeric(timeOfSampleCollectionRelativeToInclusion)),
             isRelapse, isDenovo, isTransformed,
             ageAtDiagnosis, ELN2017, vitalStatus,
             overallSurvival = suppressWarnings(as.numeric(overallSurvival)),
             blasts_bm = suppressWarnings(as.numeric(`%.Blasts.in.BM`)),
             blasts_pb = suppressWarnings(as.numeric(`%.Blasts.in.PB`)),
             flt3_itd = `FLT3-ITD`, flt3_itd_ar = suppressWarnings(as.numeric(allelic_ratio)),
             npm1 = NPM1, runx1 = RUNX1, asxl1 = ASXL1, tp53 = TP53,
             cebpa_biallelic = CEBPA_Biallelic, fusions = consensusAMLFusions,
             has_karyotype = !is.na(karyotype) & karyotype != "",
             rnaSeq = f(rnaSeq), analysisRnaSeq = f(analysisRnaSeq),
             exomeSeq = f(exomeSeq), analysisExomeSeq = f(analysisExomeSeq),
             analysisDrug = f(analysisDrug))]
# binary genotype flags -- codings DIFFER by column (checked 2026-09-29): NPM1/FLT3-ITD are
# "positive"/"negative"; TP53/RUNX1/ASXL1 are free-text variant descriptions (empty = negative);
# CEBPA_Biallelic is "bi"/"mono"/"N/A"/empty.
td[, `:=`(npm1_pos = !is.na(npm1) & npm1 == "positive",
          flt3_itd_pos = !is.na(flt3_itd) & flt3_itd == "positive",
          tp53_pos = !is.na(tp53) & tp53 != "",
          runx1_pos = !is.na(runx1) & runx1 != "",
          asxl1_pos = !is.na(asxl1) & asxl1 != "",
          cebpa_bi = !is.na(cebpa_biallelic) & cebpa_biallelic == "bi")]
cat("[1] tidy clinical:", nrow(td), "specimens\n")
cat("    diseaseStageAtSpecimenCollection values:\n")
print(sort(table(td$diseaseStageAtSpecimenCollection), decreasing = TRUE))
cat("    NPM1 values: "); print(table(td$npm1, useNA = "ifany"))
cat("    FLT3-ITD values: "); print(table(td$flt3_itd, useNA = "ifany"))

# ---- 2. frozen sets -----------------------------------------------------------------------------
td[, set_rna := analysisRnaSeq]
td[, set_rna_drug := analysisRnaSeq & analysisDrug]
td[, set_triple := analysisRnaSeq & analysisExomeSeq & analysisDrug]
stopifnot(sum(td$set_rna) == 671L, sum(td$set_rna_drug) == 520L, sum(td$set_triple) == 494L,
          uniqueN(td$dbgap_subject_id[td$set_triple]) == 458L)     # known answers from 13.1

# one specimen per patient inside SET_TRIPLE, deterministic
tri <- td[set_triple == TRUE]
tri[, pref_stage := grepl("Initial Diagnosis", diseaseStageAtSpecimenCollection, fixed = TRUE)]
tri[, pref_bm := grepl("Marrow|BM", specimenType, ignore.case = TRUE)]
setorder(tri, dbgap_subject_id, -pref_stage, time_rel_inclusion, -pref_bm,
         dbgap_rnaseq_sample, na.last = TRUE)
tri[, one_per_patient := seq_len(.N) == 1L, by = dbgap_subject_id]
td[, set_triple_1pp := dbgap_rnaseq_sample %in% tri[one_per_patient == TRUE, dbgap_rnaseq_sample]]
stopifnot(sum(td$set_triple_1pp) == 458L)
cat(sprintf("\n[2] frozen: RNA %d | RNA+drug %d | triple %d (%d pts) | one-per-patient %d\n",
            sum(td$set_rna), sum(td$set_rna_drug), sum(td$set_triple),
            uniqueN(td$dbgap_subject_id[td$set_triple]), sum(td$set_triple_1pp)))
cat("    1pp composition: initial-diagnosis stage",
    sum(tri[one_per_patient == TRUE, pref_stage]), "| denovo",
    sum(td[set_triple_1pp == TRUE, isDenovo] == TRUE, na.rm = TRUE), "| relapse specimens",
    sum(td[set_triple_1pp == TRUE, isRelapse] == TRUE, na.rm = TRUE), "\n")
g1 <- td[set_triple_1pp == TRUE]
cat(sprintf("    genotype in 1pp: NPM1+ %d | FLT3-ITD+ %d | TP53+ %d | RUNX1+ %d | ASXL1+ %d | CEBPA-bi %d\n",
            sum(g1$npm1_pos), sum(g1$flt3_itd_pos), sum(g1$tp53_pos),
            sum(g1$runx1_pos), sum(g1$asxl1_pos), sum(g1$cebpa_bi)))
cat("    blasts_bm coverage in 1pp:", sum(!is.na(td[set_triple_1pp == TRUE, blasts_bm])), "of 458\n")

fwrite(td, file.path(OUT, "beataml_clinical_tidy.csv"))
fwrite(td[, .(dbgap_subject_id, dbgap_rnaseq_sample, set_rna, set_rna_drug, set_triple, set_triple_1pp)],
       file.path(OUT, "cohort_freeze.csv"))

# ---- 3. drug AUC long ---------------------------------------------------------------------------
pf <- fread(file.path(DIR, "beataml_probit_curve_fits_v4_dbgap.txt"))
au <- pf[type == "single-agent",
         .(dbgap_subject_id, dbgap_rnaseq_sample, dbgap_dnaseq_sample, inhibitor,
           auc, ic50, converged, status, paper_inclusion)]
fwrite(au, file.path(OUT, "beataml_drug_auc_long.csv"))
cat(sprintf("\n[3] drug AUC long: %d fits | %d inhibitors | converged %s\n",
            nrow(au), uniqueN(au$inhibitor), paste(names(table(au$converged)), collapse = "/")))

# ---- SELF-CHECKS --------------------------------------------------------------------------------
cnt_cols <- setdiff(names(fread(file.path(DIR, "beataml_waves1to4_counts_dbgap.txt"), nrows = 0)),
                    c("stable_id", "display_label", "description", "biotype"))
cat("\n===== SELF-CHECKS =====\n")
in_mat <- td[set_rna == TRUE, dbgap_rnaseq_sample] %in% cnt_cols
cat(sprintf("  SET_RNA specimens present as matrix columns: %d of %d\n", sum(in_mat), length(in_mat)))
stopifnot(all(in_mat))
extra <- setdiff(cnt_cols, td$dbgap_rnaseq_sample)
cat(sprintf("  matrix columns not in clinical: %d (controls/replicates, kept out of every set)\n",
            length(extra)))
has_drug <- td[set_triple_1pp == TRUE, dbgap_rnaseq_sample] %in%
            au[converged == TRUE | converged == "TRUE", dbgap_rnaseq_sample]
cat(sprintf("  1pp specimens with >=1 converged single-agent fit: %d of 458\n", sum(has_drug)))
stopifnot(mean(has_drug) > 0.95)
cat("[done] cohort frozen; 13.3 pre-registration may reference cohort_freeze.csv\n")
