#!/usr/bin/env Rscript
# 05_pairing_model.R ----
# Doublet-pairing SMOKE, stage 5: the proposal's primary test (section 4.2; appendix A) run end to end on SIMULATED cohorts only.
# Model: per patient i, observed T-LEUK label counts y_i (K states) given their total M_i are Dirichlet-multinomial with mean
#   q_i(beta) = (1 - rho_i) * normalise(C_i[f_i * exp(beta)]) + rho_i * qFP_i         (appendix A.2)
# and a shared precision s. Null beta = 0; alternative = shared beta (sum to zero). LRT, reference distribution by parametric
# draws that also redraw C_i from its held-out counts (the "surrogate" for regenerating doublets through the pipeline);
# R_k at the reference composition fbar (A.4); simultaneous region by studentised max over states from parametric draws
# under the fitted alternative; outcome A / B-small / B-undetermined / C / D (A.5).
# Templates: the real donors of the two smoke lanes that pass the singlet floors -- their class-A singlet compositions and their
# lane's synthetic calibration (stage 04, main mode). No real pairing count exists anywhere in this script: nothing to mask.
# Worlds (typing 'ideal' = C with 0.5 on the diagonal and 0.05 to neighbouring states; 'smoke' = the C measured in stage 04):
#   null_ideal / planted_ideal      known answers on the machinery (size, R recovery, coverage)
#   null_smoke / planted_smoke      the same with today's typing
#   planted_smoke_lowyield          M_i drawn from the smoke lanes' own expected yields
#   empty                           every M_i < 5: the decision must be D, never a crash
#   (planted = R 1.5 on DBL_PLANT_STATE at fbar.) Edge cases built in: HSC-like absent (f = 0) in every fourth patient; M_i = 0.
# INPUT  : DBL_LANE_DIR/<lane>__calib_<calib>_smoke.rds (stage 04)
# OUTPUT : DBL_TAB_DIR/dbl_model_cohorts_smoke.csv    world x cohort: qualifying patients, LRT, p, R_k, region, outcome
#          DBL_TAB_DIR/dbl_model_summary_smoke.csv    world: rejection rate, mean R_k, coverage, outcomes
#          DBL_TAB_DIR/dbl_model_templates_smoke.csv  template donors: f, rho, expected yield, singlet floors
#          DBL_TAB_DIR/dbl_model_timing_smoke.csv     seconds per fit and per cohort; extrapolation to 2,000 draws
#          DBL_TAB_DIR/dbl_checks_05_smoke.csv
# Usage  : Rscript scripts/dbl/05_pairing_model.R [--calib=main|devel] [--cores=8]

suppressPackageStartupMessages({ library(data.table); library(here); library(parallel) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ARGS  <- commandArgs(TRUE)
argv  <- function(key, default) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), ARGS, value = TRUE)); if (length(v)) v else default }
CALIB <- argv("calib", "main"); CORES <- as.integer(argv("cores", Sys.getenv("SLURM_CPUS_PER_TASK", "8")))
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-78s expect %-12s got %-12s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
S <- DBL_LEUK_STATES; K <- length(S); kp <- match(DBL_PLANT_STATE, S)
DBL_BETA_BOUND <- 8   # [SMOKE] |free beta| <= 8 (R up to e^16); first BFGS run drifted to |beta| ~ 1e3 on unidentified states

## -- [1] model functions ----
dm_ll <- function(y, q, s) {   # Dirichlet-multinomial log-likelihood of counts y with mean q and precision s
  a <- s * q; M <- sum(y)
  lgamma(s) - lgamma(M + s) + sum(lgamma(y + a) - lgamma(a)) + lgamma(M + 1) - sum(lgamma(y + 1)) }
r_dm <- function(M, q, s) { g <- rgamma(length(q), shape = s * q); if (sum(g) == 0) g <- q; as.vector(rmultinom(1, M, g / sum(g))) }
q_of <- function(beta, f, Cm, rho, qfp) {   # appendix A.2: observed-label distribution of one patient
  sig <- as.vector(Cm %*% (f * exp(beta - max(beta))))   # Cm = K x K (reported x true); shift by max(beta): same after
                                                          # normalising, and exp() cannot overflow (job 3687245 did)
  sig <- if (is.finite(sum(sig)) && sum(sig) > 0) sig / sum(sig) else rep(1 / K, K)
  q <- (1 - rho) * sig + rho * qfp
  q <- pmax(q, 1e-9); q / sum(q) }
R_of <- function(beta, fbar) exp(beta) / sum(fbar * exp(beta))   # A.4: conditional O/E at the reference composition
beta_full <- function(b) c(b, -sum(b))                            # sum-to-zero
fit <- function(Y, P, alt) {   # P: list per patient (f, C, rho, qfp); returns log-lik, beta, s
  nll <- function(par) { beta <- if (alt) beta_full(par[seq_len(K - 1)]) else rep(0, K); s <- exp(par[length(par)])
    -sum(vapply(seq_along(P), function(i) dm_ll(Y[i, ], q_of(beta, P[[i]]$f, P[[i]]$C, P[[i]]$rho, P[[i]]$qfp), s), 0)) }
  start <- if (alt) c(rep(0, K - 1), log(20)) else log(20)
  o <- if (alt) optim(start, nll, method = "L-BFGS-B", lower = c(rep(-DBL_BETA_BOUND, K - 1), log(0.1)),   # bounded: a state the
                      upper = c(rep(DBL_BETA_BOUND, K - 1), log(1e5)))   # data cannot identify sits at a bound instead of drifting
       else optim(start, nll, method = "Brent", lower = log(0.1), upper = log(1e5))
  list(ll = -o$value, beta = if (alt) beta_full(o$par[seq_len(K - 1)]) else rep(0, K), s = exp(o$par[length(o$par)]), conv = o$convergence) }
draw_C <- function(Cfull, n) {   # redraw a confusion matrix from its held-out counts (Dirichlet per true state, +0.5)
  Cs <- Cfull
  for (k in seq_len(ncol(Cfull))) if (n[k] > 0) { g <- rgamma(nrow(Cfull), shape = Cfull[, k] * n[k] + 0.5); Cs[, k] <- g / sum(g) }
  Cs[seq_len(K), , drop = FALSE] }   # drop the "none" row: K x K
decide <- function(gates_ok, p, lo, hi, band = DBL_R_BAND, alpha = 0.05) {   # A.5
  if (!gates_ok) return("D")
  inside  <- all(lo >= band[1] & hi <= band[2])
  outside <- any(lo > band[2] | hi < band[1])
  if (p < alpha) { if (outside) "A" else if (inside) "B-small" else "B-undetermined" } else if (inside) "C" else "D" }

## -- [2] known answers on the functions themselves ----
chk("A.4 example: f = 0.5, weight x 1.5 -> share 0.6, R = 1.2", 1.2, R_of(c(log(1.5), 0) - mean(c(log(1.5), 0)), c(0.5, 0.5))[1])
y0 <- c(5, 3, 2, 0, 1); q0 <- c(0.4, 0.3, 0.15, 0.05, 0.1)
chk("DM log-lik -> multinomial as s -> infinity", TRUE, abs(dm_ll(y0, q0, 1e9) - dmultinom(y0, prob = q0, log = TRUE)) < 1e-4)
chk("decide(): empty / failing gates -> D", "D", decide(FALSE, 0.001, rep(2, K), rep(3, K)))
chk("decide(): reject + one interval above the band -> A", "A", decide(TRUE, 0.01, c(1.6, rep(0.8, K - 1)), c(2.5, rep(1.2, K - 1))))
chk("decide(): no reject + region inside the band -> C", "C", decide(TRUE, 0.4, rep(0.8, K), rep(1.2, K)))
chk("decide(): no reject + region wider than the band -> D", "D", decide(TRUE, 0.4, rep(0.5, K), rep(2, K)))

## -- [3] templates: six donors x their lane's calibration ----
ts("[3] templates from stage-04 calibration (", CALIB, ")")
TPL <- list()
for (L in DBL_SMOKE_LANES) {
  cb <- readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_", CALIB, "_smoke.rds")))
  Cfull <- cb$C; Cfull[is.na(Cfull)] <- 0; Cfull[nrow(Cfull), colSums(Cfull) == 0] <- 1   # unevaluable state: never reported
  nk <- setNames(rep(0, K), S); nk[cb$n_heldout$true_state] <- cb$n_heldout$N
  for (d in unique(cb$expected$donor)) {
    e  <- cb$expected[donor == d]
    cp <- cb$comp[donor == d]
    f  <- setNames(vapply(S, function(s) cp[compartment == s, sum(N)], 0), S); nL <- sum(f); f <- if (nL > 0) f / nL else f
    sig <- sum(e$exp_signal); fp <- sum(e$exp_fp)
    TPL[[length(TPL) + 1]] <- list(lane = L, donor = d, f = f, Cfull = Cfull, n = nk, rho = if (sig + fp > 0) fp / (sig + fp) else 0,
      qfp = if (fp > 0) e$exp_fp / fp else rep(1 / K, K), exp_total = sig + fp,
      n_T = cp[compartment == "T", sum(N)], n_LEUK = nL) } }
TT <- rbindlist(lapply(TPL, function(t) data.table(lane = t$lane, donor = t$donor, n_T = t$n_T, n_LEUK = t$n_LEUK,
        singlet_floors_ok = t$n_T >= 30 & t$n_LEUK >= 50, rho = round(t$rho, 3), exp_total_pairings = round(t$exp_total, 2),
        f = paste(sprintf("%s:%.3f", sub("LEUK_", "", S), t$f), collapse = " "))))
print(TT)
## simulated patients are drawn from donors that pass the proposal's singlet floors (A.5 operational minimum: >= 30 strict-T and
## >= 50 LEUK singlets); first smoke run (job 3687244) cycled all six, so every cohort had 8 qualifying patients and ended in D
TPL <- TPL[vapply(TPL, function(t) t$n_T >= 30 && t$n_LEUK >= 50, TRUE)]
ts(sprintf("    %d of %d template donors pass the singlet floors and seed the simulated cohorts: %s", length(TPL), nrow(TT),
           paste(vapply(TPL, `[[`, "", "donor"), collapse = ", ")))
chk("at least two template donors pass the singlet floors", TRUE, length(TPL) >= 2)
nT <- length(TPL)

## -- [4] worlds ----
fbar_t <- Reduce(`+`, lapply(TPL, `[[`, "f")) / nT   # expected reference composition (templates cycled evenly)
Rt <- rep((1 - DBL_PLANT_R * fbar_t[kp]) / (1 - fbar_t[kp]), K); Rt[kp] <- DBL_PLANT_R   # R_plant on one state, others equal
stopifnot(all(Rt > 0))
beta_plant <- log(Rt) - mean(log(Rt))
chk("planted beta reproduces the target R at the templates' fbar", TRUE, isTRUE(all.equal(unname(R_of(beta_plant, fbar_t)), unname(Rt))))
M_design <- c(8, 12, 20, 30, 50, 80)                   # [SMOKE] totals per patient, about the handoff appendix-A range
## Two kinds of typing (added after job 3687245): the confusion measured in stage 04 leaves some states almost never reported,
## so their beta is unidentified and, through the shared denominator of R_k, every R_k becomes unstable. The 'ideal' C (0.5 on
## the diagonal, 0.05 to each neighbouring state, the rest not retained) checks the statistical machinery itself; the 'smoke'
## worlds show what the proposal gets with today's typing.
C_ideal <- matrix(0, K + 1, K, dimnames = list(c(S, "none"), S))
for (k in seq_len(K)) { C_ideal[k, k] <- 0.5; if (k > 1) C_ideal[k - 1, k] <- 0.05; if (k < K) C_ideal[k + 1, k] <- 0.05 }
C_ideal[K + 1, ] <- 1 - colSums(C_ideal[seq_len(K), ])
WORLDS <- list(null_ideal       = list(beta = rep(0, K),  M = "design",   typing = "ideal", n = DBL_SIM_OUTER),
               planted_ideal    = list(beta = beta_plant, M = "design",   typing = "ideal", n = DBL_SIM_OUTER),
               null_smoke       = list(beta = rep(0, K),  M = "design",   typing = "smoke", n = DBL_SIM_OUTER),
               planted_smoke    = list(beta = beta_plant, M = "design",   typing = "smoke", n = DBL_SIM_OUTER),
               planted_smoke_lowyield = list(beta = beta_plant, M = "expected", typing = "smoke", n = DBL_SIM_OUTER),
               empty            = list(beta = rep(0, K),  M = "empty",    typing = "smoke", n = 10L))
one_cohort <- function(wn, ci) {
  W <- WORLDS[[wn]]
  set.seed(DBL_SEED_SIM + match(wn, names(WORLDS)) * 1000L + ci)
  tmpl <- rep(seq_len(nT), length.out = DBL_SIM_PATIENTS)
  Mi <- switch(W$M, design = rep(M_design, length.out = DBL_SIM_PATIENTS), empty = rep(c(0, 3, 4), length.out = DBL_SIM_PATIENTS),
               expected = rpois(DBL_SIM_PATIENTS, vapply(TPL[tmpl], `[[`, 0, "exp_total")))
  P <- lapply(seq_len(DBL_SIM_PATIENTS), function(i) { t <- TPL[[tmpl[i]]]
    g <- rgamma(K, shape = 30 * t$f + 0.01); f <- g / sum(g)                   # jittered composition around the template
    if (i %% 4 == 0) { f[1] <- 0; f <- f / sum(f) }                            # edge case: the HSC-like state absent
    Cf <- if (W$typing == "ideal") C_ideal else t$Cfull
    nk <- if (W$typing == "ideal") setNames(rep(120, K), S) else t$n           # held-out doublets per state behind C
    list(f = f, Ctrue = draw_C(Cf, nk), C = Cf[seq_len(K), , drop = FALSE], Cfull = Cf, n = nk,
         rho = t$rho, qfp = t$qfp, M = Mi[i], floors = t$n_T >= 30 & t$n_LEUK >= 50) })
  Y <- t(vapply(P, function(p) r_dm(p$M, q_of(W$beta, p$f, p$Ctrue, p$rho, p$qfp), 20), numeric(K)))   # true s = 20
  qual <- which(vapply(P, function(p) p$M >= 5 && p$floors, TRUE))
  base <- data.table(world = wn, cohort = ci, n_qual = length(qual), M_total = sum(Mi[qual]))
  if (length(qual) < 10) return(cbind(base, outcome = "D", reason = "fewer than 10 qualifying patients"))   # A.5 operational minimum
  Pq <- P[qual]; Yq <- Y[qual, , drop = FALSE]
  fbar <- Reduce(`+`, lapply(Pq, `[[`, "f")) / length(Pq)
  t0 <- proc.time()[["elapsed"]]
  f0 <- fit(Yq, Pq, FALSE); f1 <- fit(Yq, Pq, TRUE); T_obs <- 2 * (f1$ll - f0$ll)
  t_fit <- (proc.time()[["elapsed"]] - t0) / 2
  Tb <- vapply(seq_len(DBL_SIM_BOOT), function(b) {    # reference distribution under the fitted null, C redrawn per draw
    Yb <- t(vapply(Pq, function(p) r_dm(p$M, q_of(rep(0, K), p$f, draw_C(p$Cfull, p$n), p$rho, p$qfp), f0$s), numeric(K)))
    2 * (fit(Yb, Pq, TRUE)$ll - fit(Yb, Pq, FALSE)$ll) }, 0)
  p <- (1 + sum(Tb >= T_obs)) / (DBL_SIM_BOOT + 1)
  Rhat <- R_of(f1$beta, fbar)
  Rb <- t(vapply(seq_len(DBL_SIM_BOOT), function(b) {  # parametric draws under the fitted alternative -> simultaneous region
    Yb <- t(vapply(Pq, function(p) r_dm(p$M, q_of(f1$beta, p$f, draw_C(p$Cfull, p$n), p$rho, p$qfp), f1$s), numeric(K)))
    R_of(fit(Yb, Pq, TRUE)$beta, fbar) }, numeric(K)))
  lR <- log(Rb); sdk <- apply(lR, 2, sd); zmax <- apply(abs(sweep(lR, 2, log(Rhat))) / rep(sdk, each = nrow(lR)), 1, max)
  cq <- quantile(zmax, 0.95); lo <- exp(log(Rhat) - cq * sdk); hi <- exp(log(Rhat) + cq * sdk)
  Rtrue <- R_of(W$beta, fbar)
  t_cohort <- proc.time()[["elapsed"]] - t0
  cbind(base, T = T_obs, p = p, s_null = f0$s, s_alt = f1$s, conv = f0$conv + f1$conv,
        as.data.table(setNames(as.list(Rhat), paste0("R_", sub("LEUK_", "", S)))),
        as.data.table(setNames(as.list(lo), paste0("lo_", sub("LEUK_", "", S)))),
        as.data.table(setNames(as.list(hi), paste0("hi_", sub("LEUK_", "", S)))),
        covered = all(Rtrue >= lo & Rtrue <= hi), Rtrue_plant = Rtrue[kp],
        outcome = decide(TRUE, p, lo, hi), reason = "gates on patient counts only (typing gate reported separately)",
        sec_per_fit = t_fit, sec_cohort = t_cohort) }
JOBS <- rbindlist(lapply(names(WORLDS), function(wn) data.table(wn = wn, ci = seq_len(WORLDS[[wn]]$n))))
ts(sprintf("[4] %d simulated cohorts (%s) x %d reference draws, %d patients each, on %d cores", nrow(JOBS),
           paste(names(WORLDS), collapse = "/"), DBL_SIM_BOOT, DBL_SIM_PATIENTS, CORES))
RES <- mclapply(seq_len(nrow(JOBS)), function(j) tryCatch(one_cohort(JOBS$wn[j], JOBS$ci[j]),
                 error = function(e) data.table(world = JOBS$wn[j], cohort = JOBS$ci[j], outcome = "ERROR", reason = conditionMessage(e))),
                mc.cores = CORES, mc.preschedule = FALSE)
COH <- rbindlist(RES, fill = TRUE)
for (cn in c("p", "covered", paste0("R_", sub("LEUK_", "", DBL_PLANT_STATE)), "conv", "sec_per_fit", "sec_cohort", "n_qual", "M_total"))
  if (!cn %in% names(COH)) COH[, (cn) := NA_real_]   # empty-result branch: a world where no cohort reached the test
chk("every simulated cohort finished without an error", 0L, COH[outcome == "ERROR", .N])
if (COH[outcome == "ERROR", .N]) print(head(COH[outcome == "ERROR", .(world, cohort, reason)]))
SUM <- COH[outcome != "ERROR", .(cohorts = .N, n_qual_median = as.numeric(median(n_qual)), M_total_median = as.numeric(median(M_total)),
           reject_rate = mean(p < 0.05, na.rm = TRUE), mean_R_plant = mean(get(paste0("R_", sub("LEUK_", "", DBL_PLANT_STATE))), na.rm = TRUE),
           coverage = mean(covered, na.rm = TRUE), outcomes = paste(names(table(outcome)), table(outcome), sep = ":", collapse = " ")), by = world]
print(SUM)
nr <- SUM[world == "null_ideal", reject_rate]   # known answers on the machinery: ideal-typing worlds
chk("null_ideal: rejection rate at alpha 0.05 in [0, 0.15] (60 cohorts)", TRUE, isTRUE(nr <= 0.15))
chk("planted_ideal: mean R-hat of the planted state in [1.3, 1.7]", TRUE, isTRUE(abs(SUM[world == "planted_ideal", mean_R_plant] - DBL_PLANT_R) <= 0.2))
chk("planted_ideal: simultaneous region covers the true R in >= 85 % of cohorts", TRUE, isTRUE(SUM[world == "planted_ideal", coverage] >= 0.85))
chk("null_ideal: region covers R = 1 in >= 85 % of cohorts", TRUE, isTRUE(SUM[world == "null_ideal", coverage] >= 0.85))
chk("null_smoke: rejection rate at alpha 0.05 in [0, 0.15] (size with today's typing)", TRUE, isTRUE(SUM[world == "null_smoke", reject_rate] <= 0.15))
chk("empty world: every outcome is D", TRUE, COH[world == "empty", all(outcome == "D")])
chk("every optimiser converged", 0L, COH[!is.na(conv) & conv != 0, .N])
## -- [5] timing and extrapolation ----
RUNS <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_runs_", CALIB, "_smoke.csv")))
sec_dbl <- RUNS[kind %in% c("heldout", "real"), median(secs / (n_real + n_injected) * 1000)]   # seconds per 1,000 barcodes
TIM <- data.table(sec_per_fit_median = COH[, median(sec_per_fit, na.rm = TRUE)], sec_per_cohort_median = COH[, median(sec_cohort, na.rm = TRUE)],
                  draws_smoke = 2 * DBL_SIM_BOOT, scdbl_sec_per_1000_barcodes = sec_dbl)
TIM[, `:=`(surrogate_cpu_h_per_2000_draws = sec_per_cohort_median / draws_smoke * 2 * 2000 / 3600,
           full_pipeline_cpu_h_per_2000_draws = sec_dbl * 1.14e3 * 2000 / 3600)]   # one cohort pass ~ 1.14 M barcodes (impl map)
print(TIM)
fwrite(COH, file.path(DBL_TAB_DIR, "dbl_model_cohorts_smoke.csv")); fwrite(SUM, file.path(DBL_TAB_DIR, "dbl_model_summary_smoke.csv"))
fwrite(TT, file.path(DBL_TAB_DIR, "dbl_model_templates_smoke.csv")); fwrite(TIM, file.path(DBL_TAB_DIR, "dbl_model_timing_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_05_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
