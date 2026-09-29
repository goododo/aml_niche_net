# 13.8_reducibility_interrogation.R ----
# INPUT  : results/tables/13_bulk_anchor/{tcga_case_table,beataml_clinical_tidy,hox_bulk_scores}.csv
#          /FAST/gr10634/gaozy/bulk_cohorts/tcga_laml/data_clinical_patient.txt (CYTOGENETICS text)
# OUTPUT : results/tables/13_bulk_anchor/itdneg_reducibility.csv
# WHAT IT DOES : the mandatory reducibility interrogation of the ONE positive result this project
#          has produced -- 13.7's registered co-primary, HR 1.253 [1.091-1.439] per SD of HOX
#          score for overall survival in FLT3-ITD-negative NPM1-WT AML (I2 = 0%, two cohorts).
#          The rule this enforces was promoted to project policy on 2026-09-28 after SELL
#          deflated to CD34: a hit must survive conditioning on the known biology that could
#          produce it. Here the decisive objections, each turned into an arm:
#            A KMT2A/MLL rearrangement -- textbook HOX-high AND textbook adverse prognosis.
#              Tested by EXCLUSION (drop those patients, refit), not only by adjustment.
#            B other HOX-high fusions (NUP98, DEK-NUP214, MECOM) -- same treatment.
#            C cytogenetic/molecular risk class (TCGA RISK_CYTO/RISK_MOLECULAR, BeatAML ELN2017)
#              -- if the HR vanishes, the score is a risk-class proxy.
#            D blast percentage -- the composition analogue that killed earlier arms.
#            E specimen timing -- OS runs from diagnosis, but 1pp specimens are not all initial
#              diagnosis; a score measured later than the survival origin is a time-alignment
#              confound. Tested by restricting to initial-diagnosis specimens.
#            F monotonicity -- tertiles, to see whether the effect is a gradient or one tail.
#          Every arm reports the HR per SD with CI so a shrunken-but-present effect is
#          distinguishable from an abolished one. kind_of_evidence = reducibility_test.

suppressPackageStartupMessages({ library(data.table); library(survival) })
O <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"
HOX_FUSIONS <- "KMT2A|MLLT3|MLL|NUP98|NUP214|MECOM"

# ---- assemble, keeping the covariates the arms need --------------------------------------------
tc <- fread(file.path(O, "tcga_case_table.csv"))[has_rna == TRUE & sequenced == TRUE]
tcl <- fread("/FAST/gr10634/gaozy/bulk_cohorts/tcga_laml/data_clinical_patient.txt", skip = 4L)
tc <- merge(tc, tcl[, .(patientId = PATIENT_ID, cyto_text = CYTOGENETICS,
                        risk = RISK_MOLECULAR, risk2 = RISK_CYTO)], by = "patientId")
TC <- tc[, .(cohort = "TCGA-LAML", score, npm1 = npm1_pos, itd = itd_pos, age,
             os_time = os_months, dead, blast = blast_bm,
             hox_fusion = grepl(HOX_FUSIONS, cyto_text, ignore.case = TRUE) |
                          grepl("11q23|9;11", cyto_text),
             risk = risk, initial_dx = TRUE)]   # TCGA-LAML is a de-novo-at-diagnosis cohort

bb <- fread(file.path(O, "beataml_clinical_tidy.csv"))[set_triple_1pp == TRUE]
bb <- merge(bb, fread(file.path(O, "hox_bulk_scores.csv"))[, .(dbgap_rnaseq_sample, score)],
            by = "dbgap_rnaseq_sample")
BE <- bb[vitalStatus %in% c("Dead", "Alive") & !is.na(overallSurvival),
         .(cohort = "BeatAML", score, npm1 = npm1_pos, itd = flt3_itd_pos, age = ageAtDiagnosis,
           os_time = overallSurvival, dead = vitalStatus == "Dead", blast = blasts_bm,
           hox_fusion = grepl(HOX_FUSIONS, fusions, ignore.case = TRUE),
           risk = ELN2017,
           initial_dx = grepl("Initial Diagnosis", diseaseStageAtSpecimenCollection, fixed = TRUE))]
D <- rbind(TC, BE)
D <- D[itd == FALSE & npm1 == FALSE & !is.na(os_time) & !is.na(age)]     # the registered stratum
cat(sprintf("[0] registered stratum: %s\n", paste(sprintf("%s n=%d (%d HOX-fusion, %d initial-dx)",
    c("TCGA", "BeatAML"), D[, .N, cohort]$N, D[, sum(hox_fusion), cohort]$V1,
    D[, sum(initial_dx), cohort]$V1), collapse = " | ")))

# ---- does the score simply mark the HOX fusions? ------------------------------------------------
cat("\n[1] is the score a fusion detector?\n")
for (co in unique(D$cohort)) {
  d <- D[cohort == co]
  if (sum(d$hox_fusion) >= 3L) {
    cat(sprintf("    %-10s fusion+ n=%d median score %.3f | fusion- n=%d median %.3f | wilcox p %.4g\n",
        co, sum(d$hox_fusion), median(d[hox_fusion == TRUE, score]), sum(!d$hox_fusion),
        median(d[hox_fusion == FALSE, score]), wilcox.test(score ~ hox_fusion, d)$p.value))
    q <- quantile(d$score, 0.9)
    cat(sprintf("               of the top HOX decile (n=%d), %d are fusion+ (%.0f%%)\n",
                sum(d$score >= q), sum(d$score >= q & d$hox_fusion), 100 * mean(d[score >= q, hox_fusion])))
  } else cat(sprintf("    %-10s only %d fusion+ patients -- exclusion arm still run\n", co, sum(d$hox_fusion)))
}

# ---- the arms ------------------------------------------------------------------------------------
fit <- function(d, form) {
  d <- copy(d)[, z := (score - mean(score)) / sd(score)]
  s <- summary(coxph(as.formula(form), d))
  list(n = nrow(d), events = sum(d$dead), logHR = s$coefficients["z", "coef"],
       se = s$coefficients["z", "se(coef)"], HR = s$conf.int["z", 1],
       lo = s$conf.int["z", 3], hi = s$conf.int["z", 4], p = s$coefficients["z", 5])
}
meta <- function(rows) {
  w <- 1 / rows$se^2; m <- sum(w * rows$logHR) / sum(w); se <- sqrt(1 / sum(w))
  Q <- sum(w * (rows$logHR - m)^2)
  data.table(arm = rows$arm[1], cohort = "META", n = sum(rows$n), events = sum(rows$events),
             logHR = m, se = se, HR = exp(m), lo = exp(m - 1.96 * se), hi = exp(m + 1.96 * se),
             p = 2 * pnorm(-abs(m / se)), I2 = max(0, (Q - (nrow(rows) - 1)) / Q) * 100)
}
ARMS <- list(
  registered      = list(sub = quote(TRUE),                     form = "Surv(os_time, dead) ~ z + age"),
  exclude_fusions = list(sub = quote(hox_fusion == FALSE),       form = "Surv(os_time, dead) ~ z + age"),
  adjust_risk     = list(sub = quote(risk %in% c("Good", "Intermediate", "Poor", "Favorable", "Adverse")),
                                                                form = "Surv(os_time, dead) ~ z + age + risk"),
  adjust_blast    = list(sub = quote(!is.na(blast)),             form = "Surv(os_time, dead) ~ z + age + blast"),
  initial_dx_only = list(sub = quote(initial_dx == TRUE),        form = "Surv(os_time, dead) ~ z + age"),
  excl_fus_adj_risk = list(sub = quote(hox_fusion == FALSE &
                             risk %in% c("Good", "Intermediate", "Poor", "Favorable", "Adverse")),
                                                                form = "Surv(os_time, dead) ~ z + age + risk")
)
res <- list()
for (a in names(ARMS)) {
  rows <- list()
  for (co in unique(D$cohort)) {
    d <- D[cohort == co][eval(ARMS[[a]]$sub)]
    if (nrow(d) < 40L || sum(d$dead) < 15L || uniqueN(d$score) < 10L) next
    f <- ARMS[[a]]$form
    if (grepl("risk", f) && uniqueN(d$risk) < 2L) f <- sub(" \\+ risk", "", f)
    rows[[co]] <- data.table(arm = a, cohort = co, as.data.table(fit(d, f)), I2 = NA_real_)
  }
  rows <- rbindlist(rows)
  res[[a]] <- if (nrow(rows) >= 2L) rbind(rows, meta(rows)) else rows
}
R <- rbindlist(res)[, kind_of_evidence := "reducibility_test"]
fwrite(R, file.path(O, "itdneg_reducibility.csv"))
cat("\n[2] REDUCIBILITY ARMS (HR per 1 SD of HOX score; registered arm first)\n")
print(R[, .(arm, cohort, n, events, HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3),
            p = signif(p, 3), I2 = round(I2))], nrows = 40)

# ---- monotonicity --------------------------------------------------------------------------------
cat("\n[3] monotonicity: tertiles of score within cohort (HR vs lowest tertile)\n")
for (co in unique(D$cohort)) {
  d <- copy(D[cohort == co])[, ter := cut(frank(score), 3, labels = c("T1", "T2", "T3"))]
  s <- summary(coxph(Surv(os_time, dead) ~ ter + age, d))
  cat(sprintf("    %-10s T2 HR %.2f [%.2f-%.2f] | T3 HR %.2f [%.2f-%.2f]\n", co,
      s$conf.int["terT2", 1], s$conf.int["terT2", 3], s$conf.int["terT2", 4],
      s$conf.int["terT3", 1], s$conf.int["terT3", 3], s$conf.int["terT3", 4]))
  cat(sprintf("               median OS by tertile: %s\n", paste(sprintf("%s %.1f",
      c("T1", "T2", "T3"), sapply(c("T1", "T2", "T3"), function(t)
        median(d[ter == t, os_time]))), collapse = " | ")))
}

## SELF-CHECK: a scrambled score must lose the effect in both cohorts ----------------------------
set.seed(1)
cat("\n[SELF-CHECK] score permuted within cohort (5 draws), registered model:\n")
for (i in 1:5) {
  rows <- rbindlist(lapply(unique(D$cohort), function(co) {
    d <- copy(D[cohort == co])[, score := sample(score)]
    data.table(arm = "perm", cohort = co, as.data.table(fit(d, "Surv(os_time, dead) ~ z + age")), I2 = NA_real_)
  }))
  m <- meta(rows)
  cat(sprintf("  draw %d: meta HR %.3f p %.3f\n", i, m$HR, m$p))
}
cat("[done]\n")
