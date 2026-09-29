# 13.5_tcga_intake.R ----
# INPUT  : /FAST/gr10634/gaozy/bulk_cohorts/tcga_laml/  (cBioPortal datahub laml_tcga_pub, fetched
#            2026-09-29 via media.githubusercontent LFS URLs: data_mrna_seq_v2_rsem.txt 20,531 x 173,
#            data_mutations.txt, data_clinical_patient.txt; plus tcga_laml_case_table.tsv from the
#            2026-09-29 access scout)
# OUTPUT : results/tables/13_bulk_anchor/tcga_case_table.csv    (per patient: genotype, survival, score)
#          results/tables/13_bulk_anchor/tcga_intake_checks.csv (known-answer checks)
# WHAT IT DOES : intake + genotype derivation for TCGA-LAML as the second stratum of
#          PREREGISTRATION_itdneg_conditioned.md, and it settles that document's §2.1 hard
#          precondition. FLT3-ITD is NOT a curated column anywhere open: it must be DERIVED from
#          the WashU curated MAF, where exon-14 juxtamembrane in-frame insertions (chr13
#          ~28,608,2xx GRCh37, many annotated literally as ...dup) are the ITDs. The rule, fixed
#          here before any survival is touched:
#            ITD  = FLT3 AND mutation type contains "In_Frame"
#            TKD  = FLT3 missense with protein position 820-860 (activation loop, e.g. D835)
#          Known answers from the scout's independent derivation: 200 patients, NPM1 54, ITD 30,
#          TKD 17, ITD-and-TKD 0, OS complete 200/200 with 133 deceased. The RT-PCR panel in the
#          GDC biotab (FLT3 58+/135-, NPMc 46+/150-) is a SEPARATE channel and is used only as a
#          concordance readout -- it cannot separate ITD from TKD, so it never overrides the MAF.

suppressPackageStartupMessages(library(data.table))
DIR <- "/FAST/gr10634/gaozy/bulk_cohorts/tcga_laml"
OUT <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"
CHK <- list()
chk <- function(name, expected, observed) {
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected),
    observed = as.character(observed), pass = identical(as.character(expected), as.character(observed)))
  cat(sprintf("  %-40s expect %-8s got %-8s %s\n", name, expected, observed,
              if (identical(as.character(expected), as.character(observed))) "PASS" else "** FAIL **"))
}

# ---- 1. clinical + survival ---------------------------------------------------------------------
cl <- fread(file.path(DIR, "data_clinical_patient.txt"), skip = 4L)
cat("[1] clinical:", nrow(cl), "patients |", ncol(cl), "columns\n")
cl <- cl[, .(patientId = PATIENT_ID, age = suppressWarnings(as.numeric(AGE)), sex = SEX,
             os_months = suppressWarnings(as.numeric(OS_MONTHS)), os_status = OS_STATUS,
             risk_molecular = RISK_MOLECULAR, risk_cyto = RISK_CYTO, fab = FAB,
             blast_bm = suppressWarnings(as.numeric(BM_BLAST_PERCENTAGE)))]
cl[, dead := grepl("DECEASED", os_status)]
chk("patients", 200, nrow(cl))
chk("OS complete", 200, sum(!is.na(cl$os_months)))
chk("deceased", 133, sum(cl$dead))

# ---- 2. genotype from the WashU curated MAF ------------------------------------------------------
mu <- fread(file.path(DIR, "data_mutations.txt"))
mu[, patientId := substr(Tumor_Sample_Barcode, 1, 12)]
ppos <- suppressWarnings(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", mu$HGVSp_Short)))
fl <- mu[Hugo_Symbol == "FLT3"]
fl_pos <- suppressWarnings(as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", fl$HGVSp_Short)))
itd <- unique(fl[grepl("In_Frame", Variant_Classification), patientId])
tkd <- unique(fl[Variant_Classification == "Missense_Mutation" &
                 !is.na(fl_pos) & fl_pos >= 820 & fl_pos <= 860, patientId])
npm1 <- unique(mu[Hugo_Symbol == "NPM1", patientId])
cat("\n[2] genotype derived from MAF\n")
cat("    FLT3 rows:", nrow(fl), "| classifications:", paste(names(table(fl$Variant_Classification)), collapse = "/"), "\n")
cat("    example ITD protein changes:", paste(head(unique(fl[grepl("In_Frame", Variant_Classification), HGVSp_Short]), 5), collapse = ", "), "\n")
chk("NPM1-mutant patients", 54, length(npm1))
chk("FLT3-ITD patients", 30, length(itd))
chk("FLT3-TKD patients", 17, length(tkd))
chk("ITD and TKD overlap", 0, length(intersect(itd, tkd)))

cl[, `:=`(npm1_pos = patientId %in% npm1, itd_pos = patientId %in% itd, tkd_pos = patientId %in% tkd)]
cl[, sequenced := patientId %in% mu$patientId]
cat(sprintf("    patients with any MAF row (sequenced): %d\n", sum(cl$sequenced)))

# ---- 3. expression + the Set-14 score -----------------------------------------------------------
SET14 <- c("HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10",
           "HOXB2", "HOXB3", "HOXB4", "HOXB5", "HOXB6", "MEIS1", "PBX3")  # frozen in preregistration 2 section 3
ex <- fread(file.path(DIR, "data_mrna_seq_v2_rsem.txt"))
scols <- setdiff(names(ex), c("Hugo_Symbol", "Entrez_Gene_Id"))
X <- as.matrix(ex[, ..scols]); sym <- ex$Hugo_Symbol
ok <- rowSums(!is.finite(X)) == 0
X <- X[ok, ]; sym <- sym[ok]
rm_dup <- ave(rowMeans(X), sym, FUN = function(m) seq_along(m) == which.max(m)) == 1
X <- X[rm_dup, ]; sym <- sym[rm_dup]
idx <- which(sym %in% SET14)
cat(sprintf("\n[3] expression: %d genes x %d samples | Set-14 matched %d of 14\n",
            nrow(X), ncol(X), length(idx)))
chk("Set-14 genes matched", 14, length(idx))
P <- apply(X, 2, frank) / nrow(X)
sc <- data.table(patientId = substr(scols, 1, 12), score = colMeans(P[idx, ]))
chk("expression samples", 173, nrow(sc))
stopifnot(!anyDuplicated(sc$patientId))

ct <- merge(cl, sc, by = "patientId", all.x = TRUE)
ct[, has_rna := !is.na(score)]

# ---- 4. the registered instrument gate on this cohort -------------------------------------------
sq <- ct[has_rna == TRUE & sequenced == TRUE]
auc <- as.numeric(wilcox.test(sq[npm1_pos == TRUE, score], sq[npm1_pos == FALSE, score])$statistic) /
       (sum(sq$npm1_pos) * sum(!sq$npm1_pos))
cat(sprintf("\n[4] instrument (Set-14, NPM1-mut %d vs WT %d): AUC = %.3f | registered gate >= 0.75 -> %s\n",
            sum(sq$npm1_pos), sum(!sq$npm1_pos), auc, if (auc >= 0.75) "PASS" else "FAIL"))

# the ITD-negative stratum this cohort can contribute (group sizes only, no outcome yet)
neg <- sq[itd_pos == FALSE]
thr <- quantile(neg[npm1_pos == TRUE, score], 0.25)
cat(sprintf("    ITD-negative: %d (NPM1-mut %d, WT %d)\n", nrow(neg), sum(neg$npm1_pos), sum(!neg$npm1_pos)))
cat(sprintf("    threshold (25th pct of ITD-neg mutants) = %.3f -> NPM1-like %d | WT-low %d\n",
            thr, sum(!neg$npm1_pos & neg$score >= thr), sum(!neg$npm1_pos & neg$score < thr)))
cat(sprintf("    registered group-size gate (>=10 NPM1-like) -> %s\n",
            if (sum(!neg$npm1_pos & neg$score >= thr) >= 10L) "PASS" else "FAIL"))

fwrite(ct, file.path(OUT, "tcga_case_table.csv"))
CHK[[length(CHK) + 1]] <- data.table(check = "instrument_auc_set14", expected = ">=0.75",
  observed = sprintf("%.3f", auc), pass = auc >= 0.75)
res <- rbindlist(CHK)
fwrite(res, file.path(OUT, "tcga_intake_checks.csv"))
cat(sprintf("\n[5] %d checks, %d PASS, %d FAIL\n", nrow(res), sum(res$pass), sum(!res$pass)))
if (any(!res$pass)) print(res[pass == FALSE])

## SELF-CHECK: derived ITD/NPM1 must agree with the independent RT-PCR channel where available ----
sct <- file.path(DIR, "tcga_laml_case_table.tsv")
if (file.exists(sct)) {
  s <- fread(sct)
  jc <- intersect(grep("patient|barcode", names(s), ignore.case = TRUE, value = TRUE)[1], names(s))
  rt <- grep("rtpcr|npmc|flt3_rt|clinical_flt3", names(s), ignore.case = TRUE, value = TRUE)
  cat("\n[SELF-CHECK] scout table columns available for concordance:", paste(rt, collapse = ", "), "\n")
  if (length(jc) && length(rt)) {
    m <- merge(ct[, .(patientId, npm1_pos, itd_pos)], s, by.x = "patientId", by.y = jc)
    for (r in rt) { tb <- table(m[[r]], m$npm1_pos); cat("  ", r, "vs derived NPM1:\n"); print(tb) }
  }
}
cat("[done]\n")
