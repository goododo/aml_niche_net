# 13.3_hox_transfer.R ----
# INPUT  : bulk_cohorts/beataml2/beataml_waves1to4_norm_exp_dbgap.txt
#          results/tables/13_bulk_anchor/{cohort_freeze,beataml_clinical_tidy,beataml_drug_auc_long}.csv
#          results/tables/12_pseudobulk_de/hox_axis_genes.csv   (Set L, LMPP_GMP, frozen by P2)
# OUTPUT : results/tables/13_bulk_anchor/{hox_bulk_scores,transfer_q0_gate,transfer_q1_q2,transfer_q3_screen}.csv
# WHAT IT DOES : implements PREREGISTRATION_bulk_transfer.md (committed 28befd0 BEFORE this
#          script existed), pilot arm approved 2026-09-29: venetoclax + NPM1-like + HOX.
#          Gate order is the registered order and later stages DO NOT RUN if an earlier gate
#          stops: score -> Q0 instrument gate (AUC>=0.80 pass / 0.70-0.80 exploratory / <0.70
#          STOP) -> G1 random-gene-set specificity (>=97.5%) -> NPM1-like grouping -> G2
#          confound table (printed before any outcome is joined) -> Q1 venetoclax -> Q2 OS ->
#          Q3 exploratory screen over the remaining single agents.

suppressPackageStartupMessages({ library(data.table); library(survival) })
TAB <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
D13 <- file.path(TAB, "13_bulk_anchor")
set.seed(260929)

fz <- fread(file.path(D13, "cohort_freeze.csv"))[set_triple_1pp == TRUE]
cl <- fread(file.path(D13, "beataml_clinical_tidy.csv"))[dbgap_rnaseq_sample %in% fz$dbgap_rnaseq_sample]
stopifnot(nrow(cl) == 458L)

# ---- 1. Set L, the frozen 19 ----------------------------------------------------------------
gl <- fread(file.path(TAB, "12_pseudobulk_de", "hox_axis_genes.csv"))
SETL <- sort(unique(gl[set == "L" & hierarchy_bin == "LMPP_GMP", gene]))
cat(sprintf("[1] Set L (LMPP_GMP): %d genes: %s\n", length(SETL), paste(SETL, collapse = " ")))
stopifnot(length(SETL) >= 15L, length(SETL) <= 25L)

# ---- 2. percentile-rank scores on the 458 ----------------------------------------------------
nx <- fread("/FAST/gr10634/gaozy/bulk_cohorts/beataml2/beataml_waves1to4_norm_exp_dbgap.txt")
meta_cols <- c("stable_id", "display_label", "description", "biotype")
X <- as.matrix(nx[, setdiff(names(nx), meta_cols), with = FALSE])[, cl$dbgap_rnaseq_sample]
sym <- nx$display_label
# duplicates: keep the highest-mean row per symbol (rule stated here, applied before matching)
rm_dup <- ave(rowMeans(X), sym, FUN = function(m) seq_along(m) == which.max(m)) == 1
X <- X[rm_dup, ]; sym <- sym[rm_dup]
idx <- which(sym %in% SETL)
cat(sprintf("[2] matrix %d unique symbols x %d samples | Set L matched %d of %d\n",
            nrow(X), ncol(X), length(idx), length(SETL)))
if (length(idx) < 15L) stop("REGISTERED STOP: <15 of Set L matched -- fix symbol mapping")
P <- apply(X, 2, frank) / nrow(X)                 # per-sample percentile of every gene
score <- colMeans(P[idx, ])
sc <- data.table(dbgap_rnaseq_sample = colnames(X), score = score)
sc <- merge(sc, cl, by = "dbgap_rnaseq_sample")

# ---- 3. Q0 instrument gate -------------------------------------------------------------------
auc_of <- function(x, g) { w <- wilcox.test(x[g], x[!g])$statistic; as.numeric(w) / (sum(g) * sum(!g)) }
q0 <- auc_of(sc$score, sc$npm1_pos)
cat(sprintf("\n[3] Q0: HOX score separates NPM1-mut (n=%d) from WT (n=%d): AUC = %.3f\n",
            sum(sc$npm1_pos), sum(!sc$npm1_pos), q0))
G0 <- if (q0 >= 0.80) "PASS" else if (q0 >= 0.70) "EXPLORATORY" else "STOP"

# G1 specificity: 1000 size- and expression-decile-matched random sets
dec <- cut(frank(rowMeans(P)), 10, labels = FALSE)
pool <- split(seq_len(nrow(P)), dec)
need <- table(dec[idx])
null_auc <- replicate(1000, {
  ridx <- unlist(mapply(function(d, k) sample(pool[[d]], k), names(need), need))
  auc_of(colMeans(P[ridx, , drop = FALSE]), sc$npm1_pos)
})
g1_pct <- mean(q0 > null_auc)
cat(sprintf("    G1: Set L beats %.1f%% of 1000 matched random sets (need >=97.5%%)\n", 100 * g1_pct))
G1 <- if (g1_pct >= 0.975) "PASS" else "STOP"
fwrite(data.table(check = c("Q0_auc", "G1_pct_beaten"), value = c(q0, g1_pct),
                  verdict = c(G0, G1), kind_of_evidence = "registered_primary"),
       file.path(D13, "transfer_q0_gate.csv"))
if (G0 == "STOP" || G1 == "STOP") { cat("** REGISTERED STOP -- no Q1/Q2/Q3 **\n"); quit(save = "no") }
tag <- if (G0 == "PASS") "registered_primary" else "exploratory_screen"

# ---- 4. NPM1-like grouping + G2 confound table (before outcomes) ------------------------------
thr <- quantile(sc[npm1_pos == TRUE, score], 0.25)
sc[, grp := fifelse(npm1_pos, "NPM1mut", fifelse(score >= thr, "NPM1like", "WTlow"))]
cat(sprintf("\n[4] threshold (25th pct of mut) = %.3f -> NPM1-like %d | WT-low %d | mut %d\n",
            thr, sum(sc$grp == "NPM1like"), sum(sc$grp == "WTlow"), sum(sc$grp == "NPM1mut")))
wt <- sc[grp != "NPM1mut"]
g2 <- data.table(
  variable = c("blasts_bm", "ageAtDiagnosis", "FLT3-ITD rate", "TP53 rate"),
  npm1like = c(median(wt[grp == "NPM1like", blasts_bm], na.rm = TRUE),
               median(wt[grp == "NPM1like", ageAtDiagnosis], na.rm = TRUE),
               mean(wt[grp == "NPM1like", flt3_itd_pos]), mean(wt[grp == "NPM1like", tp53_pos])),
  wtlow    = c(median(wt[grp == "WTlow", blasts_bm], na.rm = TRUE),
               median(wt[grp == "WTlow", ageAtDiagnosis], na.rm = TRUE),
               mean(wt[grp == "WTlow", flt3_itd_pos]), mean(wt[grp == "WTlow", tp53_pos])),
  p = c(wilcox.test(blasts_bm ~ grp, wt)$p.value,
        wilcox.test(ageAtDiagnosis ~ grp, wt)$p.value,
        fisher.test(table(wt$grp, wt$flt3_itd_pos))$p.value,
        fisher.test(table(wt$grp, wt$tp53_pos))$p.value))
cat("    G2 confound table (mandatory caveat, not a stop):\n"); print(g2, digits = 3)
fwrite(sc[, .(dbgap_rnaseq_sample, dbgap_subject_id, score, npm1_pos, grp)],
       file.path(D13, "hox_bulk_scores.csv"))

# ---- 5. Q1 venetoclax + Q2 survival ------------------------------------------------------------
au <- fread(file.path(D13, "beataml_drug_auc_long.csv"))[converged == TRUE]
au <- au[, .(auc = mean(auc)), by = .(dbgap_rnaseq_sample, inhibitor)]   # combine replicate fits
ve <- merge(sc, au[inhibitor == "Venetoclax"], by = "dbgap_rnaseq_sample")
cat(sprintf("\n[5] venetoclax fits: %d of 458 (like %d, WTlow %d, mut %d)\n", nrow(ve),
            sum(ve$grp == "NPM1like"), sum(ve$grp == "WTlow"), sum(ve$grp == "NPM1mut")))
cat(sprintf("    context sanity: median AUC mut %.1f vs WT %.1f (literature: mut MORE sensitive = lower)\n",
            median(ve[grp == "NPM1mut", auc]), median(ve[grp != "NPM1mut", auc])))
vew <- ve[grp != "NPM1mut"]
q1 <- wilcox.test(auc ~ grp, vew)
d_med <- median(vew[grp == "NPM1like", auc]) - median(vew[grp == "WTlow", auc])
lm1 <- summary(lm(auc ~ (grp == "NPM1like") + ageAtDiagnosis, vew))$coefficients[2, ]
lm2 <- summary(lm(auc ~ (grp == "NPM1like") + ageAtDiagnosis + blasts_bm,
                  vew[!is.na(blasts_bm)]))$coefficients[2, ]
rho_cont <- cor.test(~ score + auc, vew, method = "spearman")
cat(sprintf("    Q1 PRIMARY: NPM1-like vs WT-low venetoclax AUC: median diff %+.1f, wilcox p = %.4g\n",
            d_med, q1$p.value))
cat(sprintf("    adjusted (age): beta %+.1f p %.4g | (age+blasts, n=%d): beta %+.1f p %.4g\n",
            lm1[1], lm1[4], sum(!is.na(vew$blasts_bm)), lm2[1], lm2[4]))
cat(sprintf("    sensitivity continuous (WT only): spearman rho %+.3f p %.4g\n",
            rho_cont$estimate, rho_cont$p.value))

sv <- sc[vitalStatus %in% c("Dead", "Alive") & !is.na(overallSurvival) & grp != "NPM1mut"]
cx <- summary(coxph(Surv(overallSurvival, vitalStatus == "Dead") ~ (grp == "NPM1like") + ageAtDiagnosis, sv))
cat(sprintf("    Q2: OS Cox (n=%d, events=%d): HR %.2f [%.2f-%.2f], p = %.4g\n",
            nrow(sv), sum(sv$vitalStatus == "Dead"), cx$conf.int[1, 1], cx$conf.int[1, 3],
            cx$conf.int[1, 4], cx$coefficients[1, 5]))
fwrite(data.table(
  test = c("Q1_wilcox", "Q1_lm_age", "Q1_lm_age_blasts", "Q1_continuous_rho", "Q2_cox_HR"),
  estimate = c(d_med, lm1[1], lm2[1], rho_cont$estimate, cx$conf.int[1, 1]),
  p = c(q1$p.value, lm1[4], lm2[4], rho_cont$p.value, cx$coefficients[1, 5]),
  kind_of_evidence = c(tag, "sensitivity", "sensitivity", "sensitivity", tag)),
  file.path(D13, "transfer_q1_q2.csv"))

# ---- 6. Q3 exploratory screen ------------------------------------------------------------------
scr <- list()
for (dr in setdiff(unique(au$inhibitor), "Venetoclax")) {
  dd <- merge(vew[, .(dbgap_rnaseq_sample, grp)], au[inhibitor == dr], by = "dbgap_rnaseq_sample")
  if (sum(dd$grp == "NPM1like") < 10L || sum(dd$grp == "WTlow") < 30L) next
  scr[[dr]] <- data.table(inhibitor = dr,
    n_like = sum(dd$grp == "NPM1like"), n_wtlow = sum(dd$grp == "WTlow"),
    median_diff = median(dd[grp == "NPM1like", auc]) - median(dd[grp == "WTlow", auc]),
    p = wilcox.test(auc ~ grp, dd)$p.value)
}
q3 <- rbindlist(scr)[, q_bh := p.adjust(p, "BH")][, kind_of_evidence := "exploratory_screen"]
setorder(q3, p)
fwrite(q3, file.path(D13, "transfer_q3_screen.csv"))
cat(sprintf("\n[6] Q3 screen: %d drugs tested | q_BH<0.10: %d | top 8:\n", nrow(q3), sum(q3$q_bh < 0.10)))
print(q3[1:min(8, .N)], digits = 3)

## SELF-CHECK: scoring machinery -- a score built from the matrix's top-expressed genes must be
## ~1 for every sample (percentile ranks near the top), i.e. the rank direction is not inverted.
top_idx <- order(rowMeans(P), decreasing = TRUE)[1:19]
stopifnot(min(colMeans(P[top_idx, ])) > 0.95)
cat("\n[SELF-CHECK] rank direction OK (top-expressed control score > 0.95 in all samples)\n[done]\n")
