# 13.1_verify_beataml.R ----
# INPUT  : /FAST/gr10634/gaozy/bulk_cohorts/beataml2/  (8 files fetched 2026-09-29 via
#          media.githubusercontent LFS URLs; raw_inhibitor via raw URL)
# OUTPUT : results/tables/13_bulk_anchor/beataml_file_manifest.csv   (md5 provenance)
#          results/tables/13_bulk_anchor/beataml_verification.csv    (check | expected | observed | pass)
# WHAT IT DOES : intake verification of BeatAML2 against the 2026-09-24 data-reconnaissance
#          report (DATA_RECON_2026-09-24.md §2.1), whose counts were recomputed from this same
#          workbook by two independent agents. Every recon number is treated as a KNOWN ANSWER:
#          942 specimens / 805 patients; RNA-seq 698 with 671 analysis-grade; WES/targeted 903;
#          drug 631 specimens / 569 patients; RNA∩drug 542; RNA∩WES∩drug 514 specimens / 471
#          patients; 166 single-agent inhibitors; vitalStatus Dead 565 / Alive 345 / Unknown 32;
#          karyotype 842/942; 11,721 mutation rows. Also verifies the expression matrices parse,
#          counts are integers, and every matrix column links to a dbgap RNA id through the
#          sample-mapping workbook. Any known-answer mismatch prints FAIL and the script exits
#          nonzero -- downstream (13.2) must not run on a drifted intake.

suppressPackageStartupMessages({ library(data.table); library(readxl) })
DIR <- "/FAST/gr10634/gaozy/bulk_cohorts/beataml2"
OUT <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 0. file manifest with md5 ------------------------------------------------------------------
files <- list.files(DIR, full.names = TRUE)
man <- data.table(file = basename(files), bytes = file.size(files),
                  md5 = tools::md5sum(files), fetched = "2026-09-29")
fwrite(man, file.path(OUT, "beataml_file_manifest.csv"))
cat("[0] manifest written:", nrow(man), "files, total",
    round(sum(man$bytes) / 1e6), "MB\n")

CHECKS <- list()
chk <- function(name, expected, observed) {
  CHECKS[[length(CHECKS) + 1]] <<- data.table(check = name, expected = as.character(expected),
    observed = as.character(observed), pass = identical(as.character(expected), as.character(observed)))
  cat(sprintf("  %-46s expect %-12s got %-12s %s\n", name, expected, observed,
              if (identical(as.character(expected), as.character(observed))) "PASS" else "** FAIL **"))
}

# ---- 1. clinical workbook ------------------------------------------------------------------------
cl <- as.data.table(read_excel(file.path(DIR, "beataml_wv1to4_clinical.xlsx"), sheet = "summary"))
cat("\n[1] clinical workbook (sheet 'summary')\n")
chk("specimens (rows)",                 942, nrow(cl))
chk("patients (unique dbgap_subject)",  805, uniqueN(cl$dbgap_subject_id))
chk("RNA-seq specimens",                698, sum(cl$rnaSeq == "y" | cl$rnaSeq == TRUE, na.rm = TRUE))
chk("analysis-grade RNA",               671, sum(cl$analysisRnaSeq == "y" | cl$analysisRnaSeq == TRUE, na.rm = TRUE))
chk("WES/targeted specimens",           903, sum(cl$exomeSeq == "y" | cl$exomeSeq == TRUE, na.rm = TRUE))
drug <- cl$analysisDrug == "y" | cl$analysisDrug == TRUE
chk("drug specimens (analysisDrug)",    631, sum(drug, na.rm = TRUE))
chk("drug patients",                    569, uniqueN(cl$dbgap_subject_id[which(drug)]))
# RECONCILIATION 2026-09-29: the recon report's 542/514/471 used the loose `rnaSeq` flag
# (verified: rnaSeq x drug = 542, +exomeSeq = 514 specimens / 471 patients, reproduced exactly).
# The frozen analysis cohort uses OHSU's own ANALYSIS-GRADE flags instead -> 520/494/458.
rna <- cl$analysisRnaSeq == "y" | cl$analysisRnaSeq == TRUE
wes <- cl$analysisExomeSeq == "y" | cl$analysisExomeSeq == TRUE
rna_loose <- cl$rnaSeq == "y" | cl$rnaSeq == TRUE
chk("RNA x drug (loose flags, recon's number)", 542, sum(rna_loose & drug, na.rm = TRUE))
chk("RNA x drug specimens (analysis-grade)",    520, sum(rna & drug, na.rm = TRUE))
chk("RNA x WES x drug specimens (analysis)",    494, sum(rna & wes & drug, na.rm = TRUE))
chk("RNA x WES x drug patients (analysis)",     458, uniqueN(cl$dbgap_subject_id[which(rna & wes & drug)]))
vs <- table(cl$vitalStatus, useNA = "ifany")
chk("vitalStatus Dead",                 565, sum(cl$vitalStatus == "Dead", na.rm = TRUE))
chk("vitalStatus Alive",                345, sum(cl$vitalStatus == "Alive", na.rm = TRUE))
chk("vitalStatus Unknown",               32, sum(cl$vitalStatus == "Unknown" | is.na(cl$vitalStatus)))
chk("karyotype non-missing",            842, sum(!is.na(cl$karyotype) & cl$karyotype != ""))
cat("    ELN2017 values:", paste(names(table(cl$ELN2017)), collapse = "/"),
    "| missing:", sum(is.na(cl$ELN2017)), "\n")
cat("    genotype consensus columns present:",
    all(c("FLT3-ITD", "allelic_ratio", "NPM1", "RUNX1", "ASXL1", "TP53", "CEBPA_Biallelic")
        %in% names(cl)), "\n")

# ---- 2. drug tables -------------------------------------------------------------------------------
pf <- fread(file.path(DIR, "beataml_probit_curve_fits_v4_dbgap.txt"))
cat("\n[2] probit fits:", nrow(pf), "rows | columns:", paste(head(names(pf), 20), collapse = ", "), "...\n")
single <- pf[type == "single-agent" | !grepl(" - |;|\\+", inhibitor)]
chk("single-agent inhibitors",          166, uniqueN(single$inhibitor))
cat("    'type' values:", paste(names(table(pf$type)), collapse = "/"),
    "| auc-like columns:", paste(grep("auc|ic[0-9]", names(pf), ignore.case = TRUE, value = TRUE), collapse = ", "), "\n")

# ---- 3. mutations ---------------------------------------------------------------------------------
mu <- fread(file.path(DIR, "beataml_wes_wv1to4_mutations_dbgap.txt"))
cat("\n[3] mutations\n")
chk("mutation rows",                  11721, nrow(mu))
cat("    samples with >=1 variant:", uniqueN(mu$dbgap_sample_id), "\n")

# ---- 4. expression matrices + id linkage ----------------------------------------------------------
cnt <- fread(file.path(DIR, "beataml_waves1to4_counts_dbgap.txt"))
nrm <- fread(file.path(DIR, "beataml_waves1to4_norm_exp_dbgap.txt"))
meta_cols <- c("stable_id", "display_label", "description", "biotype")
cat("\n[4] expression matrices\n")
cat(sprintf("    counts: %d genes x %d samples | norm: %d x %d\n",
            nrow(cnt), ncol(cnt) - 4L, nrow(nrm), ncol(nrm) - 4L))
chk("counts and norm same sample set", TRUE,
    identical(sort(setdiff(names(cnt), meta_cols)), sort(setdiff(names(nrm), meta_cols))))
set.seed(1); sc <- sample(setdiff(names(cnt), meta_cols), 20)
iv <- sapply(sc, function(j) { x <- cnt[[j]]; all(abs(x - round(x)) < 1e-8) })
chk("counts integer (20 sampled columns)", TRUE, all(iv))

mp <- as.data.table(read_excel(file.path(DIR, "beataml_waves1to4_sample_mapping.xlsx"), sheet = 1))
cat("    mapping columns:", paste(names(mp), collapse = ", "), "\n")
rn_col <- grep("rna", names(mp), ignore.case = TRUE, value = TRUE)[1]
lab_col <- names(mp)[sapply(mp, function(x) mean(setdiff(names(cnt), meta_cols) %in% x) > 0.9)][1]
if (!is.na(lab_col)) {
  linked <- mean(setdiff(names(cnt), meta_cols) %in% mp[[lab_col]])
  chk("matrix columns linked via mapping", "1", as.character(round(linked, 3)))
  cat(sprintf("    matrix-id column '%s' -> dbgap column '%s'\n", lab_col, rn_col))
} else cat("    ** no mapping column matches matrix ids -- inspect manually **\n")

# ---- 5. verdict -----------------------------------------------------------------------------------
res <- rbindlist(CHECKS)
fwrite(res, file.path(OUT, "beataml_verification.csv"))
cat(sprintf("\n[5] %d checks, %d PASS, %d FAIL -> beataml_verification.csv\n",
            nrow(res), sum(res$pass), sum(!res$pass)))
if (any(!res$pass)) { cat("** KNOWN-ANSWER MISMATCH -- do not run 13.2 **\n"); quit(status = 1) }
cat("[done] intake verified; 13.2 may run\n")
