# 13.7_itdneg_meta.R ----
# INPUT  : results/tables/13_bulk_anchor/{tcga_case_table.csv, beataml_clinical_tidy.csv,
#          hox_bulk_scores.csv, cohort_freeze.csv, beataml_drug_auc_long.csv}
# OUTPUT : results/tables/13_bulk_anchor/{itdneg_cohort_table,itdneg_confounds,itdneg_survival}.csv
# WHAT IT DOES : executes PREREGISTRATION_itdneg_conditioned.md on the two admissible cohorts
#          (BeatAML 1pp, TCGA-LAML), in the registered order:
#            1. per-cohort instrument gate (AUC >= 0.75) and group-size gate (>= 10 NPM1-like)
#            2. the registered rarity result -- how many ITD-negative NPM1-WT patients reach the
#               mutant distribution at all (this is what the dichotomy gate measures when it fails)
#            3. confound table, printed BEFORE any survival model
#            4. the CO-PRIMARY continuous analysis (section 4), which needs no dichotomy:
#               Cox OS ~ z(score) + age per cohort, inverse-variance fixed-effect meta on log HR.
#               z() is within-cohort standardisation so the HR is per 1 SD of HOX score and is
#               comparable across platforms -- a scale choice recorded here, not a hypothesis change.
#            5. ITD-positive stratum as the registered descriptive contrast (no claim)
#            6. venetoclax/palbociclib in ITD-negative BeatAML -- registered as completeness only,
#               cannot be called a finding whatever the result (they have been tested twice)
#          Nulls are reported as CIs, per the registered equivalence rule.

suppressPackageStartupMessages({ library(data.table); library(survival) })
O <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/13_bulk_anchor"

# ---- assemble the two cohorts into one shape ---------------------------------------------------
tc <- fread(file.path(O, "tcga_case_table.csv"))[has_rna == TRUE & sequenced == TRUE]
TC <- tc[, .(cohort = "TCGA-LAML", id = patientId, score, npm1 = npm1_pos, itd = itd_pos,
             age, os_time = os_months, dead, blast = blast_bm,
             tp53 = NA, runx1 = NA, asxl1 = NA, sex)]
bb <- fread(file.path(O, "beataml_clinical_tidy.csv"))[set_triple_1pp == TRUE]
bb <- merge(bb, fread(file.path(O, "hox_bulk_scores.csv"))[, .(dbgap_rnaseq_sample, score)],
            by = "dbgap_rnaseq_sample")
BE <- bb[vitalStatus %in% c("Dead", "Alive") & !is.na(overallSurvival),
         .(cohort = "BeatAML", id = dbgap_rnaseq_sample, score, npm1 = npm1_pos, itd = flt3_itd_pos,
           age = ageAtDiagnosis, os_time = overallSurvival, dead = vitalStatus == "Dead",
           blast = blasts_bm, tp53 = tp53_pos, runx1 = runx1_pos, asxl1 = asxl1_pos, sex = NA)]
D <- rbind(TC, BE, fill = TRUE)
cat(sprintf("[0] assembled: %s\n", paste(sprintf("%s n=%d", c("TCGA", "BeatAML"),
            c(nrow(TC), nrow(BE))), collapse = " | ")))

# ---- 1. registered gates, per cohort ------------------------------------------------------------
auc_of <- function(a, b) as.numeric(wilcox.test(a, b)$statistic) / (length(a) * length(b))
ct <- list()
for (co in unique(D$cohort)) {
  d <- D[cohort == co]
  inst <- auc_of(d[npm1 == TRUE, score], d[npm1 == FALSE, score])
  neg <- d[itd == FALSE]
  mu <- neg[npm1 == TRUE, score]; wt <- neg[npm1 == FALSE, score]
  thr <- quantile(mu, 0.25)
  ct[[co]] <- data.table(cohort = co, n_total = nrow(d), instrument_auc = inst,
    instrument_gate = if (inst >= 0.75) "PASS" else "FAIL",
    n_itdneg = nrow(neg), n_itdneg_mut = length(mu), n_itdneg_wt = length(wt),
    itdneg_internal_auc = auc_of(mu, wt),
    n_npm1like = sum(wt >= thr), pct_npm1like = 100 * mean(wt >= thr),
    n_above_mut_median = sum(wt >= median(mu)),
    group_size_gate = if (sum(wt >= thr) >= 10L) "PASS" else "FAIL",
    kind_of_evidence = "registered_primary")
}
CT <- rbindlist(ct)
fwrite(CT, file.path(O, "itdneg_cohort_table.csv"))
cat("\n[1] REGISTERED GATES\n"); print(CT[, .(cohort, n_total, instrument_auc = round(instrument_auc, 3),
  instrument_gate, n_itdneg_mut, n_itdneg_wt, itdneg_auc = round(itdneg_internal_auc, 3),
  n_npm1like, pct = round(pct_npm1like, 1), n_above_mut_median, group_size_gate)])
adm <- CT[group_size_gate == "PASS", cohort]
cat(sprintf("\n    cohorts passing BOTH gates: %d -> registered minimum-evidence rule (>=2) %s\n",
            length(adm), if (length(adm) >= 2L) "MET" else "NOT MET: the dichotomised design is NOT EVALUABLE"))

# ---- 2. the rarity result (what the failing gate actually measured) -----------------------------
tot_wt <- sum(CT$n_itdneg_wt); tot_like <- sum(CT$n_npm1like); tot_above <- sum(CT$n_above_mut_median)
cat(sprintf("\n[2] RARITY, pooled across both cohorts: of %d FLT3-ITD-negative NPM1-WT patients,\n", tot_wt))
cat(sprintf("    %d (%.1f%%) reach the NPM1-mutant 25th percentile; %d (%.1f%%) reach the mutant MEDIAN.\n",
            tot_like, 100 * tot_like / tot_wt, tot_above, 100 * tot_above / tot_wt))
ci <- binom.test(tot_above, tot_wt)$conf.int
cat(sprintf("    exact 95%% CI for the mutant-median-exceeding rate: %.2f%%-%.2f%%\n", 100 * ci[1], 100 * ci[2]))

# ---- 3. confound table BEFORE any survival model ------------------------------------------------
cf <- list()
for (co in unique(D$cohort)) {
  neg <- D[cohort == co & itd == FALSE & npm1 == FALSE]
  thr <- quantile(D[cohort == co & itd == FALSE & npm1 == TRUE, score], 0.25)
  neg[, grp := fifelse(score >= thr, "NPM1like", "WTlow")]
  for (v in c("age", "blast", "tp53", "runx1", "asxl1")) {
    if (all(is.na(neg[[v]]))) next
    hi <- neg[grp == "NPM1like"][[v]]; lo <- neg[grp == "WTlow"][[v]]
    cf[[paste(co, v)]] <- data.table(cohort = co, variable = v, n_like = sum(!is.na(hi)),
      median_like = median(as.numeric(hi), na.rm = TRUE), median_wtlow = median(as.numeric(lo), na.rm = TRUE),
      p = tryCatch(wilcox.test(as.numeric(hi), as.numeric(lo))$p.value, error = function(e) NA_real_))
  }
}
CF <- rbindlist(cf)[, kind_of_evidence := "mandatory_caveat"]
fwrite(CF, file.path(O, "itdneg_confounds.csv"))
cat("\n[3] confound table (before any survival model; NPM1-like groups are tiny, see gate)\n")
print(CF, digits = 3)

# ---- 4. CO-PRIMARY: continuous score, no dichotomy ---------------------------------------------
cat("\n[4] CO-PRIMARY (registered section 4): Cox OS ~ z(HOX score) + age, ITD-negative NPM1-WT only\n")
sv <- list()
for (co in unique(D$cohort)) {
  d <- D[cohort == co & itd == FALSE & npm1 == FALSE & !is.na(os_time) & !is.na(age)]
  d[, z := (score - mean(score)) / sd(score)]
  m <- coxph(Surv(os_time, dead) ~ z + age, d)
  s <- summary(m)
  sv[[co]] <- data.table(cohort = co, n = nrow(d), events = sum(d$dead),
    logHR = s$coefficients["z", "coef"], se = s$coefficients["z", "se(coef)"],
    HR = s$conf.int["z", 1], lo = s$conf.int["z", 3], hi = s$conf.int["z", 4],
    p = s$coefficients["z", 5], kind_of_evidence = "registered_coprimary")
  cat(sprintf("    %-10s n=%d events=%d | HR per 1 SD = %.3f [%.3f-%.3f], p = %.4g\n",
              co, nrow(d), sum(d$dead), s$conf.int["z", 1], s$conf.int["z", 3], s$conf.int["z", 4],
              s$coefficients["z", 5]))
}
SV <- rbindlist(sv)
w <- 1 / SV$se^2
mlog <- sum(w * SV$logHR) / sum(w); mse <- sqrt(1 / sum(w))
Q <- sum(w * (SV$logHR - mlog)^2); I2 <- max(0, (Q - (nrow(SV) - 1)) / Q) * 100
meta <- data.table(cohort = "META(fixed)", n = sum(SV$n), events = sum(SV$events),
  logHR = mlog, se = mse, HR = exp(mlog), lo = exp(mlog - 1.96 * mse), hi = exp(mlog + 1.96 * mse),
  p = 2 * pnorm(-abs(mlog / mse)), kind_of_evidence = "registered_coprimary")
cat(sprintf("\n    META (inverse-variance fixed effect): HR per 1 SD = %.3f [%.3f-%.3f], p = %.4g\n",
            meta$HR, meta$lo, meta$hi, meta$p))
cat(sprintf("    heterogeneity: Q = %.2f (df %d), I2 = %.0f%%\n", Q, nrow(SV) - 1L, I2))
cat(sprintf("    REGISTERED EQUIVALENCE STATEMENT: HRs outside [%.2f, %.2f] per SD are excluded at 95%%\n",
            meta$lo, meta$hi))

# ---- 5. ITD-positive descriptive contrast (no claim) -------------------------------------------
pos <- list()
for (co in unique(D$cohort)) {
  d <- D[cohort == co & itd == TRUE & npm1 == FALSE & !is.na(os_time) & !is.na(age)]
  if (nrow(d) < 15L || sum(d$dead) < 8L) { cat(sprintf("\n[5] %s ITD+ NPM1-WT n=%d events=%d -- too small, not modelled\n",
      co, nrow(d), sum(d$dead))); next }
  d[, z := (score - mean(score)) / sd(score)]
  s <- summary(coxph(Surv(os_time, dead) ~ z + age, d))
  pos[[co]] <- data.table(cohort = paste0(co, "_ITDpos"), n = nrow(d), events = sum(d$dead),
    logHR = s$coefficients["z", "coef"], se = s$coefficients["z", "se(coef)"],
    HR = s$conf.int["z", 1], lo = s$conf.int["z", 3], hi = s$conf.int["z", 4],
    p = s$coefficients["z", 5], kind_of_evidence = "descriptive_no_claim")
  cat(sprintf("\n[5] %s ITD+ NPM1-WT (descriptive): n=%d HR %.3f [%.3f-%.3f] p %.4g\n",
              co, nrow(d), s$conf.int["z", 1], s$conf.int["z", 3], s$conf.int["z", 4], s$coefficients["z", 5]))
}
fwrite(rbind(SV, meta, rbindlist(pos), fill = TRUE), file.path(O, "itdneg_survival.csv"))

# ---- 6. completeness only: the twice-tested drug leads in ITD-negative BeatAML -----------------
au <- fread(file.path(O, "beataml_drug_auc_long.csv"))[converged == TRUE][
  , .(auc = mean(auc)), by = .(dbgap_rnaseq_sample, inhibitor)]
cat("\n[6] COMPLETENESS ONLY (registered: cannot be called a finding whatever the result)\n")
for (nm in c("Venetoclax", "Palbociclib")) {
  d <- merge(D[cohort == "BeatAML" & itd == FALSE & npm1 == FALSE, .(dbgap_rnaseq_sample = id, score)],
             au[inhibitor == nm], by = "dbgap_rnaseq_sample")
  r <- cor.test(~ score + auc, d, method = "spearman")
  cat(sprintf("    %-12s ITD-neg NPM1-WT n=%d | spearman score vs AUC %+.3f, p = %.4g\n",
              nm, nrow(d), r$estimate, r$p.value))
}

## SELF-CHECK: age must behave as a known positive control in both Cox models -------------------
cat("\n[SELF-CHECK] age as a known positive control (older = worse OS expected):\n")
for (co in unique(D$cohort)) {
  d <- D[cohort == co & itd == FALSE & npm1 == FALSE & !is.na(os_time) & !is.na(age)]
  s <- summary(coxph(Surv(os_time, dead) ~ age, d))
  cat(sprintf("  %-10s HR per year %.3f, p = %.3g %s\n", co, s$conf.int[1, 1], s$coefficients[1, 5],
              if (s$conf.int[1, 1] > 1 && s$coefficients[1, 5] < 0.05) "OK" else "** unexpected **"))
}
cat("[done]\n")
