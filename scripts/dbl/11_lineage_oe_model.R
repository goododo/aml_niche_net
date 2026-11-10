#!/usr/bin/env Rscript
# 11_lineage_oe_model.R ----
# Doublet-pairing SMOKE, stage 11 (user decision 2026-11-10: "design the lineage-level O/E test"): power of the v1.1 Q1 /
# v1.3 secondary question -- are T-LEUK doublets more frequent than random co-encapsulation predicts? -- on SIMULATED cohorts,
# with stage 10's confirmed typing M1 (kNN labels; post-hoc exception) and its measured calibration. No real pairing count.
# Expected flagged counts per patient (stage 08/09's random co-encapsulation model on vectors; known-answer check in [1]):
#   S  = Nd * sum_k 2 f_T f_k r_k     flagged true T-LEUK doublets at R = 1 (r_k: M1 recall of T + LEUK state k)
#   Fp = Nd * sum_u pi_u r_u          other doublet types flagged (stage 08/09's pair-type list)
#   Fs = Nsing * sum_c f_c r_c        re-entered singlets flagged (pooled per compartment)
# Analysis: O_i ~ NegBin(mean = R * S_i + Fp_i + Fs_i, size phi), one shared R; H0 R = 1 vs R > 1, one-sided LRT. Reference
#   distribution by parametric draws under the fitted null with ONE shared doublet-rate factor per draw (log-normal, sd
#   DBL_OE_SD_DBR) on the doublet-derived parts S and Fp. (Refinement of the agreed plan, made before any number: each
#   patient's calibration error is independent across patients and already enters the NB dispersion; redrawing it per draw
#   would count it twice. What the dispersion cannot absorb is a doublet-rate error shared by the cohort.) 95 % interval for R
#   from draws under the fitted alternative (stage 05's studentised form in one dimension); outcome A-D by stage 05's decide()
#   with the proposal's band 2/3-1.5.
# Simulated truth: patients cycle the qualifying smoke donors (templates); yield multiplier y_i ~ LN(0, 0.5); the library's
#   synthetic-like rates are drawn from the calibration posterior (Beta, k + 1/2, n - k + 1/2) and the analyst's calibration is
#   a fresh binomial sample of the same size from them; real-doublet rates = synthetic-like x world factors; true R_i =
#   R * exp(b_i - 0.3^2/2), b_i ~ N(0, 0.3^2) (so E[R_i] = R); one doublet-rate factor per cohort c = world x LN(0, 0.2);
#   O_i ~ Poisson(y_i * (c * (R_i * S_i + Fp_i) + Fs_i)).
# Worlds (DBL_OE_COHORTS each): R 1 / 1.25 / 1.5 / 2 with 16 patients; R 1 and 1.5 with 10 and 30; R 1 with true doublet rate
#   1.3 x assumed; R 1 with real doublets typed worse than synthetic (recall x 0.85, false positives x 1.5); empty (8 patients).
# READING RULE (fixed 2026-11-10 before any stage-11 number was computed):
#   machinery (else "not read"): null_16 rejection <= 0.10; R1.5_16 mean R-hat within 1.5 +- 0.15; R1.5_16 interval covers 1.5
#     in >= 85 % of cohorts; empty world all D; no cohort error; all fits converge; the known-answer checks pass.
#   verdict at 16 patients: power(R 1.5) >= 0.8 -> "a lineage-level O/E test is adequately powered at the proposal's scale";
#     else power(R 2) >= 0.8 -> "only a doubling is detectable with 16 patients"; else "underpowered even for a doubling; close".
#   robustness (separate lines): rejection > 0.10 in null_dbr1.3_16 -> "the doublet rate needs a dbr-free normalisation (e.g.
#     hashtag cross-donor doublets) before any real run"; in null_optimistic_16 -> "the calibration needs real doublets".
# INPUT  : DBL_TAB_DIR/dbl_nbs_rates_{tune,confirm}_smoke.csv (stage 10 --rates-only); dbl_nbs_tune_smoke.csv (known answer);
#          DBL_LANE_DIR/<lane>__calib_devel_<tag>_smoke.rds (stage 04)
# OUTPUT : DBL_TAB_DIR/dbl_oe_templates_smoke.csv, dbl_oe_cohorts_smoke.csv, dbl_oe_summary_smoke.csv, dbl_oe_verdict_smoke.csv,
#          dbl_checks_11_smoke.csv
# Usage  : Rscript scripts/dbl/11_lineage_oe_model.R [--cores=8]
# Masking: simulation only; its inputs are synthetic-doublet counts and singlet counts pooled per compartment.

suppressPackageStartupMessages({ library(data.table); library(here); library(parallel) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
argv  <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-84s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
STATES <- DBL_LEUK_STATES; LIN <- c(STATES, "LEUK_other")
FPT <- c("T|Tamb", "T|NK", "T|B", "T|Ery", "T|MonoDC", "T|T", "Tamb|LEUK*", "NK|LEUK*", "B|LEUK*", "Ery|LEUK*", "MonoDC|LEUK*", "LEUK*|LEUK*")
SDD <- DBL_OE_SD_DBR; B <- DBL_SIM_BOOT

## -- [1] calibration counts, templates, known answer ----
RT <- rbindlist(lapply(c("tune", "confirm"), function(m) fread(file.path(DBL_TAB_DIR, sprintf("dbl_nbs_rates_%s_smoke.csv", m)))))
RT <- RT[variant == DBL_OE_VARIANT & labels == DBL_OE_LABELS]
blank <- function(z) is.na(z) | z == ""                               # fread reads a written NA string back as ""
RT[, key := fifelse(kind == "singlet", paste0("singlet:", comp), paste0(type, ":", fifelse(blank(truth), "", truth)))]
chk("calibration counts present for every tag x lane", TRUE, all(CJ(DBL_OE_CAL_TAGS, DBL_SMOKE_LANES)[, paste(V1, V2)] %in% RT[, paste(tag, lane)]))
counts_of <- function(tags, L) RT[tag %in% tags & lane == L, .(n = sum(n), k = sum(k)), by = key]
make_tpl <- function(cb, d, pc, L) {        # composition weights of one donor, index of every rate in its lane's key vector
  cd <- cb$comp[donor == d]; f <- setNames(cd$N / sum(cd$N), cd$compartment)
  g <- function(cp) if (cp %in% names(f)) f[[cp]] else 0
  pd <- cb$comp[, sum(N[donor == d]) / sum(N)]; Nd <- cb$dbr * cb$N * pd^2
  gT <- g("T"); fL <- sum(vapply(LIN, g, 0))
  wF <- c(2 * gT * vapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), g, 0), gT^2, 2 * fL * vapply(c("Tamb", "NK", "B", "Ery", "MonoDC"), g, 0), fL^2)
  list(lane = L, donor = d, Nd = Nd, Ns = cb$N * pd - Nd, wS = vapply(LIN, function(kk) 2 * gT * g(kk), 0), wF = setNames(wF, FPT), wSg = f,
       iS = match(ifelse(LIN %in% STATES, paste0("T|LEUK:", LIN), "T|LEUK_other:LEUK_other"), pc$key),
       iF = match(paste0(FPT, ":"), pc$key), iSg = match(paste0("singlet:", names(f)), pc$key),
       nT = cd[compartment == "T", sum(N)], nL = cd[compartment %in% STATES, sum(N)]) }
SF <- function(tp, r) {                     # expected flagged counts (S at R = 1, Fp, Fs) for rates r aligned to the lane's keys
  rk <- r[tp$iS]; ok <- !is.na(rk); ru <- r[tp$iF]; ru[is.na(ru)] <- 0; rs <- r[tp$iSg]; rs[is.na(rs)] <- 0
  c(S = tp$Nd * sum(tp$wS[ok] * rk[ok]), Fp = tp$Nd * sum(tp$wF * ru), Fs = tp$Ns * sum(tp$wSg * rs)) }
qual_donors <- function(cb) cb$comp[, .(nT = sum(N[compartment == "T"]), nL = sum(N[compartment %in% STATES])), by = donor][nT >= 30 & nL >= 50, donor]
## known answer: tuning counts + tuning calibration reproduce stage 10's stored expected counts for M1
ref10 <- fread(file.path(DBL_TAB_DIR, "dbl_nbs_tune_smoke.csv"))[variant == DBL_OE_VARIANT & labels == DBL_OE_LABELS]
TPL <- list(); PC <- list()
for (L in DBL_SMOKE_LANES) {
  cbs <- lapply(DBL_OE_CAL_TAGS, function(tg) readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_devel_", tg, "_smoke.rds"))))
  chk(sprintf("%s: doublet rate, N and composition identical across calibration tags", L), TRUE,
      all(vapply(cbs[-1], function(z) identical(z$dbr, cbs[[1]]$dbr) && identical(z$N, cbs[[1]]$N) && isTRUE(all.equal(z$comp, cbs[[1]]$comp)), TRUE)))
  cb <- cbs[[1]]
  chk(sprintf("%s: one composition row per donor x compartment", L), 0L, anyDuplicated(cb$comp[, .(donor, compartment)]))
  pc0 <- counts_of(DBL_OE_CAL_TAGS[1], L)
  for (d in qual_donors(cb)) { t0 <- SF(make_tpl(cb, d, pc0, L), pc0$k / pc0$n); r <- ref10[tag == DBL_OE_CAL_TAGS[1] & lane == L & donor == d]
    chk(sprintf("%s %s: vector model reproduces stage 10's expected signal and false positives (1e-9)", sub("_.*", "", L), d), TRUE,
        nrow(r) == 1 && abs(t0[["S"]] - r$exp_signal) < 1e-9 && abs(t0[["Fp"]] + t0[["Fs"]] - r$exp_fp) < 1e-9) }
  PC[[L]] <- counts_of(DBL_OE_CAL_TAGS, L)                                # pooled calibration: tuning folds 1-2 + fold 3
  PC[[L]][, rec := startsWith(key, "T|LEUK:") | startsWith(key, "T|LEUK_other:")]
  for (d in qual_donors(cb)) TPL[[length(TPL) + 1]] <- make_tpl(cb, d, PC[[L]], L)
}
TT <- rbindlist(lapply(TPL, function(tp) { s <- SF(tp, PC[[tp$lane]]$k / PC[[tp$lane]]$n); rk <- (PC[[tp$lane]]$k / PC[[tp$lane]]$n)[tp$iS]; ok <- !is.na(rk)
  data.table(lane = tp$lane, donor = tp$donor, n_T = tp$nT, n_LEUK = tp$nL, Nd = round(tp$Nd, 1), S = round(s[["S"]], 2), Fp = round(s[["Fp"]], 2),
             Fs = round(s[["Fs"]], 2), recall = round(sum(tp$wS[ok] * rk[ok]) / sum(tp$wS[ok]), 3), precision = round(s[["S"]] / sum(s), 3)) }))
print(TT)
chk("at least two template donors", TRUE, length(TPL) >= 2)
nT <- length(TPL)

## -- [2] model functions ----
nb_nll <- function(logR, logphi, O, S, F) -sum(dnbinom(O, size = exp(logphi), mu = pmax(exp(logR) * S + F, 1e-9), log = TRUE))
fit_nb <- function(O, S, F, alt) {
  if (alt) { o <- optim(c(0, log(20)), function(p) nb_nll(p[1], p[2], O, S, F), method = "L-BFGS-B", lower = c(-3, log(0.1)), upper = c(3, log(1e6)))
    list(ll = -o$value, logR = o$par[1], phi = exp(o$par[2]), conv = o$convergence) }
  else { o <- optim(log(20), function(p) nb_nll(0, p, O, S, F), method = "Brent", lower = log(0.1), upper = log(1e6))
    list(ll = -o$value, logR = 0, phi = exp(o$par), conv = 0L) } }
lrt1 <- function(O, S, F) { f0 <- fit_nb(O, S, F, FALSE); f1 <- fit_nb(O, S, F, TRUE)   # one-sided: 0 unless R-hat > 1
  list(T = if (f1$logR > 0) max(0, 2 * (f1$ll - f0$ll)) else 0, f0 = f0, f1 = f1) }
decide <- function(gates_ok, p, lo, hi, band = DBL_R_BAND, alpha = 0.05) {   # stage 05 / proposal A.5, verbatim
  if (!gates_ok) return("D")
  inside  <- all(lo >= band[1] & hi <= band[2])
  outside <- any(lo > band[2] | hi < band[1])
  if (p < alpha) { if (outside) "A" else if (inside) "B-small" else "B-undetermined" } else if (inside) "C" else "D" }
## known answers on the functions
set.seed(DBL_SEED_OE); Sx <- rep(50, 400); Fx <- rep(10, 400); Ox <- rpois(400, 2 * Sx + Fx)
chk("NB fit recovers R = 2 on 400 Poisson patients (|R-hat - 2| < 0.05)", TRUE, abs(exp(fit_nb(Ox, Sx, Fx, TRUE)$logR) - 2) < 0.05)
chk("one-sided LRT is 0 when the data sit below expectation", 0, lrt1(rpois(50, 30), rep(50, 50), rep(10, 50))$T)
chk("decide(): reject + interval above the band -> A", "A", decide(TRUE, 0.01, 1.6, 2.5))
chk("decide(): no reject + interval inside the band -> C", "C", decide(TRUE, 0.4, 0.8, 1.2))
chk("decide(): no reject + interval wider than the band -> D", "D", decide(TRUE, 0.4, 0.5, 2))

## -- [3] worlds ----
w <- function(R, n, dbr = 1, rec = 1, fp = 1, coh = DBL_OE_COHORTS) list(R = R, n = n, dbr = dbr, rec = rec, fp = fp, coh = coh)
WORLDS <- list(null_16 = w(1, 16), R1.25_16 = w(1.25, 16), R1.5_16 = w(1.5, 16), R2_16 = w(2, 16),
               null_10 = w(1, 10), R1.5_10 = w(1.5, 10), null_30 = w(1, 30), R1.5_30 = w(1.5, 30),
               null_dbr1.3_16 = w(1, 16, dbr = 1.3), null_optimistic_16 = w(1, 16, rec = 0.85, fp = 1.5), empty = w(1, 8, coh = 10L))
one_cohort <- function(wn, ci) {
  W <- WORLDS[[wn]]; t0 <- proc.time()[["elapsed"]]
  set.seed(DBL_SEED_OE + match(wn, names(WORLDS)) * 1000L + ci)
  base <- data.table(world = wn, cohort = ci, n_pat = W$n)
  if (W$n < DBL_OE_MIN_PAT) return(cbind(base, outcome = "D", reason = "fewer than 10 qualifying patients"))
  tmpl <- rep(seq_len(nT), length.out = W$n)
  cdbr <- W$dbr * exp(rnorm(1, 0, SDD))                                    # the cohort's true doublet-rate factor
  P <- lapply(seq_len(W$n), function(i) { tp <- TPL[[tmpl[i]]]; pc <- PC[[tp$lane]]
    y <- exp(rnorm(1, 0, DBL_OE_SD_YIELD))
    r_syn <- rbeta(nrow(pc), pc$k + 0.5, pc$n - pc$k + 0.5)                # this library's synthetic-like rates
    r_hat <- rbinom(nrow(pc), pc$n, r_syn) / pc$n                          # the analyst's calibration of it
    r_real <- r_syn * ifelse(pc$rec, W$rec, W$fp); r_real <- pmin(r_real, 1)
    st <- SF(tp, r_real) * y; sh <- SF(tp, r_hat) * y
    Ri <- W$R * exp(rnorm(1, 0, DBL_OE_SD_HET) - DBL_OE_SD_HET^2 / 2)
    list(O = rpois(1, cdbr * (Ri * st[["S"]] + st[["Fp"]]) + st[["Fs"]]), S = sh[["S"]], Fp = sh[["Fp"]], Fs = sh[["Fs"]]) })
  O <- vapply(P, `[[`, 0, "O"); S <- vapply(P, `[[`, 0, "S"); Fp <- vapply(P, `[[`, 0, "Fp"); Fs <- vapply(P, `[[`, 0, "Fs")
  ob <- lrt1(O, S, Fp + Fs)
  Tb <- vapply(seq_len(B), function(b) { cb <- exp(rnorm(1, 0, SDD))     # reference: fitted null + one shared dbr factor
    lrt1(rnbinom(length(O), size = ob$f0$phi, mu = cb * (S + Fp) + Fs), S, Fp + Fs)$T }, 0)
  p <- (1 + sum(Tb >= ob$T)) / (B + 1)
  lRb <- vapply(seq_len(B), function(b) { cb <- exp(rnorm(1, 0, SDD))    # interval: fitted alternative + shared dbr factor
    fit_nb(rnbinom(length(O), size = ob$f1$phi, mu = cb * (exp(ob$f1$logR) * S + Fp) + Fs), S, Fp + Fs, TRUE)$logR }, 0)
  hw <- quantile(abs(lRb - ob$f1$logR), 0.95, names = FALSE); lo <- exp(ob$f1$logR - hw); hi <- exp(ob$f1$logR + hw)
  cbind(base, O_total = sum(O), E_total = sum(S + Fp + Fs), precision_hat = sum(S) / sum(S + Fp + Fs), dbr_factor = cdbr,
        T = ob$T, p = p, R_hat = exp(ob$f1$logR), lo = lo, hi = hi, phi0 = ob$f0$phi, phi1 = ob$f1$phi,
        covered = lo <= W$R & W$R <= hi, conv = ob$f1$conv, outcome = decide(TRUE, p, lo, hi), reason = "",
        sec = proc.time()[["elapsed"]] - t0) }
JOBS <- rbindlist(lapply(names(WORLDS), function(wn) data.table(wn = wn, ci = seq_len(WORLDS[[wn]]$coh))))
ts(sprintf("[3] %d simulated cohorts (%d worlds) x %d reference + %d interval draws, on %d cores", nrow(JOBS), length(WORLDS), B, B, CORES))
RES <- mclapply(seq_len(nrow(JOBS)), function(j) tryCatch(one_cohort(JOBS$wn[j], JOBS$ci[j]),
                 error = function(e) data.table(world = JOBS$wn[j], cohort = JOBS$ci[j], outcome = "ERROR", reason = conditionMessage(e))),
                mc.cores = CORES, mc.preschedule = FALSE)
COH <- rbindlist(RES, fill = TRUE)
chk("every simulated cohort finished without an error", 0L, COH[outcome == "ERROR", .N])
if (COH[outcome == "ERROR", .N]) print(head(COH[outcome == "ERROR", .(world, cohort, reason)]))
chk("every alternative fit converged", 0L, COH[!is.na(conv) & conv != 0, .N])
SUM <- COH[outcome != "ERROR", .(cohorts = .N, n_pat = n_pat[1], O_median = as.numeric(median(O_total)), precision_hat = round(mean(precision_hat), 3),
           reject_rate = mean(p < 0.05), mean_R_hat = mean(R_hat), median_R_hat = as.numeric(median(R_hat)), coverage = mean(covered),
           median_width = as.numeric(median(hi / lo)), outcomes = paste(names(table(outcome)), table(outcome), sep = ":", collapse = " ")), by = world]
print(SUM[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

## -- [4] reading rule ----
g1 <- function(wn, col) SUM[world == wn, get(col)]
chk("machinery: null_16 rejection rate <= 0.10", TRUE, isTRUE(g1("null_16", "reject_rate") <= 0.10))
chk("machinery: R1.5_16 mean R-hat within 1.5 +- 0.15", TRUE, isTRUE(abs(g1("R1.5_16", "mean_R_hat") - 1.5) <= 0.15))
chk("machinery: R1.5_16 interval covers 1.5 in >= 85 % of cohorts", TRUE, isTRUE(g1("R1.5_16", "coverage") >= 0.85))
chk("machinery: empty world -> every outcome D", TRUE, COH[world == "empty", all(outcome == "D")])
MACH <- all(rbindlist(CHK)$pass)
pw15 <- g1("R1.5_16", "reject_rate"); pw2 <- g1("R2_16", "reject_rate")
verdict <- if (!MACH) "machinery check failed: not read" else if (pw15 >= 0.8) "a lineage-level O/E test is adequately powered at the proposal's scale" else
  if (pw2 >= 0.8) "only a doubling is detectable with 16 patients" else "underpowered even for a doubling; close"
rob_dbr <- g1("null_dbr1.3_16", "reject_rate"); rob_opt <- g1("null_optimistic_16", "reject_rate")
VER <- data.table(verdict = verdict, machinery_ok = MACH, power_R1.25_16 = g1("R1.25_16", "reject_rate"), power_R1.5_16 = pw15, power_R2_16 = pw2,
                  power_R1.5_10 = g1("R1.5_10", "reject_rate"), power_R1.5_30 = g1("R1.5_30", "reject_rate"),
                  size_null_16 = g1("null_16", "reject_rate"), size_dbr1.3 = rob_dbr, size_optimistic = rob_opt,
                  robust_dbr = if (rob_dbr > 0.10) "the doublet rate needs a dbr-free normalisation before any real run" else "ok",
                  robust_calibration = if (rob_opt > 0.10) "the calibration needs real doublets before any real run" else "ok")
print(t(VER))
ts("[4] VERDICT: ", verdict)
fwrite(TT, file.path(DBL_TAB_DIR, "dbl_oe_templates_smoke.csv")); fwrite(COH, file.path(DBL_TAB_DIR, "dbl_oe_cohorts_smoke.csv"))
fwrite(SUM, file.path(DBL_TAB_DIR, "dbl_oe_summary_smoke.csv")); fwrite(VER, file.path(DBL_TAB_DIR, "dbl_oe_verdict_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_11_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
