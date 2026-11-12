#!/usr/bin/env Rscript
# 13_share_model.R ----
# Doublet-pairing redesign (share design), model stage. Only mode so far: --planning = stage 0b's planning power (plan:
# reports/2026-11_dbl_pairing_smoke/redesign_2026-11-11/STAGE0B_plan.md; user decision 2026-11-12).
# INPUT  : DBL_TAB_DIR/dbl_share_rates_smoke.csv + dbl_checks_12_smoke.csv (stage 12 --reuse-smoke; synthetic doublets only),
#          DBL_LANE_DIR/<lane>__calib_devel_<tag>_smoke.rds (class-A composition, dbr, N), DBL_TAB_DIR/dbl_oe_templates_smoke.csv
# OUTPUT : DBL_TAB_DIR/dbl_share_planning_templates_smoke.csv  per template: Nd, E[D], E[Y], pi0, slope, SE parts, binomial check
#          DBL_TAB_DIR/dbl_share_planning_power_smoke.csv      power rows (base grid, sensitivity, regression)
#          DBL_TAB_DIR/dbl_share_planning_verdict_smoke.csv    the close-rule verdict; DBL_TAB_DIR/dbl_checks_13_smoke.csv
# What it does: per template patient, the share model's expected counts from the measured synthetic rates
#   (pi(psi) = (psi S^Y + A) / (psi S^D + B); every doublet term is Nd x composition weight x rate), the Fisher SE of
#   log psi-hat plus calibration variance, then redesign/power_share.R's run() (copied verbatim) with these variances.
# FROZEN RULE (tag prereg-dbl-0b-v1, before stage 12 ran): BN-1 iff power_binding(n = DBL_0B_N, psi = 2) < DBL_0B_BAR under the
#   base assumptions: a_s = DBL_0B_A_S, b_s = a_s + DBL_0B_B_S_EXTRA, triplets DBL_0B_TRIPLET per doublet, flagged at the
#   donor's composition-weighted per-doublet rate (S + O); calibration variance by the delta method from DBL_0B_NDRAW joint
#   posterior draws (r_D ~ Beta(k_D + 1/2, n - k_D + 1/2), r_Y = r_D x Beta(k_Y + 1/2, k_D - k_Y + 1/2); Y is a subset of D)
#   on the basis DBL_0B_CAL_BASIS ("design": k/n held at the smoke's rates, n = DBL_0B_CAL_N per donor-lane, halved for P01;
#   "smoke": the smoke's own counts); tau / bias priors from config 11b. Missing pair types (capped at 1): Other|LEUK = the
#   donor-composition-weighted mean of the measured X|LEUK rates; T|Other = mean of the measured T|X rates (X in Tamb, NK,
#   B, Ery, MonoDC); non-T non-LEUK pairs: D = that T|X mean, Y = 0. Tables are written and printed only if every known
#   answer passes; otherwise only the checks and "NOT READ". Rerun: once, only for a demonstrated code bug, same seeds.
# FREEZE 1 (tag prereg-dbl-share-f1): mode --audit (Stage 3) is not written yet. Its rules are fixed in
#   scripts/dbl/prereg_share/PREREGISTRATION_SHARE_F1.md section 4 and REDESIGN_share_plan_2026-11-11.md sections 5-7: the binomial
#   GLMM with Gauss-Hermite, the one-sided LRT with parametric reference draws, the binding variants, the world list, the bars
#   and the power verdict at freeze 2; seeds DBL_SEED_SHARE_AUDIT.
# Usage  : Rscript scripts/dbl/13_share_model.R --planning [--force] | --planning --mock (stage-10 M1 counts, D = 1.5 x M1;
#          code test only, writes nothing)

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
MD5_FROZEN <- c(config_dbl = "91d7f54db724b7d87e9a82a4f3838e5e", config_paths = "b2c7ca9b7d65f2c39774e6ee2cca2523")   # configs as tagged (prereg-dbl-0b-v1)
MD5 <- unname(tools::md5sum(c(here::here("scripts", "dbl", "config_dbl.R"), here::here("scripts", "config", "config_paths.R"))))
if (!identical(MD5, unname(MD5_FROZEN))) stop("config md5 differs from the frozen literals: fix and resubmit (not a rerun)")
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))   # the dated copy that runs
ARGS <- commandArgs(TRUE); stopifnot("--planning" %in% ARGS); FORCE <- "--force" %in% ARGS; MOCK <- "--mock" %in% ARGS
OUT <- file.path(DBL_TAB_DIR, paste0("dbl_", c("share_planning_templates", "share_planning_power", "share_planning_verdict", "checks_13"), "_smoke.csv"))
names(OUT) <- c("templates", "power", "verdict", "checks")
stopifnot(all(startsWith(OUT, DBL_SMOKE_DIR)))
if (!MOCK && any(file.exists(OUT)) && !FORCE) stop("outputs exist (use --force)")
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got, tol = 1.5e-8) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE, tolerance = tol))
  fm <- function(z, d) paste(if (is.numeric(z)) signif(z, d) else z, collapse = ";")   # numbers rounded for display only
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = fm(expected, 6), got = fm(got, 6), pass = ok)
  message(sprintf("  %-84s expect %-10s got %-10s %s", name, fm(expected, 4), fm(got, 4), if (ok) "PASS" else "** FAIL **"))
  invisible(ok) }
chk(sprintf("md5 of the running copy %s (recorded)", basename(SELF)), unname(tools::md5sum(SELF)), unname(tools::md5sum(SELF)))
STATES <- DBL_LEUK_STATES; LIN <- c(STATES, "LEUK_other"); PX <- setdiff(DBL_SHARE_P, "T")   # X partners other than T
XM <- c("Tamb", "NK", "B", "Ery", "MonoDC")                        # partners injected in the smoke (Other was not)

## -- [1] measured rates, pooled over the calibration tags ----
if (MOCK) {
  r10 <- rbind(fread(file.path(DBL_TAB_DIR, "dbl_nbs_rates_tune_smoke.csv")), fread(file.path(DBL_TAB_DIR, "dbl_nbs_rates_confirm_smoke.csv")))[
    variant == "M1" & kind == "pair" & tag %in% DBL_OE_CAL_TAGS, .(tag, lane, type, truth = fifelse(truth == "", NA_character_, truth), n, k)]
  RT <- rbind(r10[, flag := "M1"], copy(r10)[, `:=`(flag = "D", k = pmin(n, round(1.5 * k) + 1L))])
} else {
  c12 <- fread(file.path(DBL_TAB_DIR, "dbl_checks_12_smoke.csv"))
  chk("stage 12's checks all PASS", TRUE, nrow(c12) > 0 && all(c12$pass))
  RT <- fread(file.path(DBL_TAB_DIR, "dbl_share_rates_smoke.csv"))[flag %in% c("M1", "D")] }
RT[, truth := fifelse(is.na(truth) | truth == "", NA_character_, truth)]          # fread may return "" or NA for an empty field
RP <- RT[, .(n = sum(n), k = sum(k)), by = .(lane, type, truth, flag)]

## -- [2] templates: stage 10/11's qualifying donors and random co-encapsulation weights ----
TPL <- rbindlist(lapply(DBL_SMOKE_LANES, function(L) { cb <- readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_devel_", DBL_OE_CAL_TAGS[1], "_smoke.rds")))
  comp <- cb$comp; qual <- comp[, .(nT = sum(N[compartment == "T"]), nL = sum(N[compartment %in% STATES])), by = donor][nT >= 30 & nL >= 50, donor]
  rbindlist(lapply(qual, function(d) { pd <- comp[, sum(N[donor == d]) / sum(N)]; Nd <- cb$dbr * cb$N * pd^2
    data.table(lane = L, donor = d, Nd = Nd, Ns = cb$N * pd - Nd, comp = list(comp[donor == d, setNames(N / sum(N), compartment)])) })) }))
o11 <- fread(file.path(DBL_TAB_DIR, "dbl_oe_templates_smoke.csv"))
chk("templates: the qualifying donors equal stage 11's", sort(o11$donor), sort(TPL$donor))
chk("templates: Nd equals stage 11's (stored to 0.1)", o11[match(TPL$donor, donor), Nd], round(TPL$Nd, 1))
chk("templates: every class-A compartment is LEUK, T or in P", TRUE, all(unlist(lapply(TPL$comp, names)) %in% c(LIN, "T", PX)))
print(TPL[, .(lane, donor, Nd = round(Nd, 1), Ns = round(Ns))])

## -- [3] expected counts: pi(psi) = (psi S^Y + A) / (psi S^D + B) ----
rate_fun <- function(L, R) {   # R: rate per row of RP (point k/n or a posterior draw); returns lookup(type, truth, flag)
  rp <- RP[lane == L]; rr <- R[RP$lane == L]
  function(type, truth = NA_character_, flag) { i <- which(rp$type == type & rp$flag == flag & ((is.na(truth) & is.na(rp$truth)) | rp$truth %in% truth))
    if (length(i) != 1) stop(sprintf("rate %s / %s / %s: %d rows", type, truth, flag, length(i))); rr[i] } }
terms <- function(tp, R, a_s = DBL_0B_A_S, b_s = DBL_0B_A_S + DBL_0B_B_S_EXTRA, trip = DBL_0B_TRIPLET, miss = 1, nd_x = 1) {
  g <- function(cp) { v <- unname(tp$comp[[1]][cp]); ifelse(is.na(v), 0, v) }   # unname keeps the c(S =, O =) names clean
  fL <- sum(g(LIN)); r <- rate_fun(tp$lane, R)
  side <- function(fl) {
    rk <- function(kk) if (kk %in% STATES) r("T|LEUK", kk, fl) else r("T|LEUK_other", "LEUK_other", fl)
    S <- sum(vapply(LIN, function(kk) 2 * g("T") * g(kk) * rk(kk), 0))
    wx <- g(XM); xl <- vapply(XM, function(x) r(paste0(x, "|LEUK*"), NA, fl), 0); tx <- vapply(XM, function(x) r(paste0("T|", x), NA, fl), 0)
    oth_l <- min(1, miss * (if (sum(wx) > 0) sum(wx * xl) / sum(wx) else mean(xl)))   # Other|LEUK (missing type; donor weights)
    oth_t <- min(1, miss * mean(tx))                                           # T|Other (missing type)
    nn    <- if (fl == "M1") 0 else min(1, miss * mean(tx))                    # non-T non-LEUK pairs (missing types)
    O <- sum(2 * g(XM) * fL * xl) + 2 * g("Other") * fL * oth_l + fL^2 * r("LEUK*|LEUK*", NA, fl) +
         sum(2 * g("T") * g(XM) * tx) + 2 * g("T") * g("Other") * oth_t + g("T")^2 * r("T|T", NA, fl) + sum(g(PX))^2 * nn
    c(S = unname(S), O = unname(O + trip * (S + O))) }   # a triplet is flagged at the donor's per-doublet rate S + O (weights sum to 1)
  y <- side("M1"); dd <- side("D"); Nd <- tp$Nd * nd_x; Ns <- tp$Ns * nd_x
  c(SY = Nd * y[["S"]], A = Nd * y[["O"]] + Ns * a_s, SD = Nd * dd[["S"]], B = Nd * dd[["O"]] + Ns * b_s) }
pi_of <- function(tm, psi) (psi * tm[["SY"]] + tm[["A"]]) / (psi * tm[["SD"]] + tm[["B"]])
info <- function(tm) { ED <- tm[["SD"]] + tm[["B"]]; p0 <- pi_of(tm, 1)   # slope = d logit pi / d log psi at psi 1 (d log pi / d log psi over 1 - pi)
  sl <- (tm[["SY"]] / (tm[["SY"]] + tm[["A"]]) - tm[["SD"]] / (tm[["SD"]] + tm[["B"]])) / (1 - p0)
  c(ED = ED, EY = ED * p0, pi0 = p0, slope = sl, var_s = 1 / (ED * p0 * (1 - p0) * sl^2)) }
R0 <- RP$k / RP$n                                                  # point rates
W <- dcast(RP, lane + type + truth ~ flag, value.var = c("n", "k"))      # Y (= M1) is a subset of D on the same synthetic doublets
chk("calibration: every pair type has M1 and D on the same n, with k(M1) <= k(D)", TRUE, W[, all(n_M1 == n_D & k_M1 <= k_D)])
iY <- match(RP[, paste(lane, type, truth, flag)], W[, paste(lane, type, truth, "M1")]); iD <- match(RP[, paste(lane, type, truth, flag)], W[, paste(lane, type, truth, "D")])
set.seed(DBL_SEED_0B)
DRAWS <- sapply(seq_len(DBL_0B_NDRAW), function(i) {   # joint posterior: r_D ~ Beta, r_Y = r_D x Beta(k_Y + 1/2, k_D - k_Y + 1/2)
  rD <- rbeta(nrow(W), W$k_D + 0.5, W$n_D - W$k_D + 0.5); rY <- rD * rbeta(nrow(W), W$k_M1 + 0.5, W$k_D - W$k_M1 + 0.5)
  ifelse(is.na(iY), rD[iD], rY[iY]) })
cal_n <- function(type) ifelse(grepl("^T\\|LEUK", type), DBL_0B_CAL_N[["T|LEUK"]], ifelse(type == "LEUK*|LEUK*", DBL_0B_CAL_N[["LEUK|LEUK"]],
                        ifelse(grepl("\\|LEUK\\*$", type), DBL_0B_CAL_N[["X|LEUK"]], DBL_0B_CAL_N[["T|X"]])))
design_draws <- function(h) { nE <- cal_n(W$type) * h; pD <- W$k_D / W$n_D; pY <- W$k_M1 / W$n_D   # same draw form at the design's n
  sapply(seq_len(DBL_0B_NDRAW), function(i) { rD <- rbeta(nrow(W), pD * nE + 0.5, (1 - pD) * nE + 0.5)
    rY <- rD * rbeta(nrow(W), pY * nE + 0.5, (pD - pY) * nE + 0.5); ifelse(is.na(iY), rD[iD], rY[iY]) }) }
DR_full <- design_draws(1); DR_half <- design_draws(0.5)          # P01 donors span 3 lanes: half the allocation per lane
DRS <- list(smoke = function(L) DRAWS, design = function(L) if (L == DBL_SMOKE_LANES[["pooled"]]) DR_half else DR_full)
stopifnot(DBL_0B_CAL_BASIS %in% names(DRS)); DRB <- DRS[[DBL_0B_CAL_BASIS]]; DRO <- DRS[[setdiff(names(DRS), DBL_0B_CAL_BASIS)]]
lgt <- function(p) log(p / (1 - p))
template_var <- function(R0, DRf, ...) rbindlist(lapply(seq_len(nrow(TPL)), function(j) { tp <- TPL[j]; tm <- terms(tp, R0, ...); f <- info(tm)
  DR <- DRf(tp$lane)
  dl <- apply(DR, 2, function(Rd) lgt(pi_of(terms(tp, Rd, ...), 1))) - lgt(f[["pi0"]])   # analyst's calibration error on logit pi0
  data.table(donor = tp$donor, t(f), var_c = var(dl / f[["slope"]])) }))[, var_t := var_s + var_c]

## -- [4] redesign/power_share.R's run(), copied verbatim (v_t is read from the calling environment) ----
run <- function(n, psi, sc, sdd, tau = 0.3, nsim = 20000, shift_comp = 0.05, info_win = 0.95, true_bias = NULL, eq_w = TRUE) {
  ds <- rep(1:3, times = round(n * c(0.45, 0.30, 0.25))); ds <- c(ds, rep(1, n - length(ds)))[1:n]
  v  <- rep(v_t, length.out = n)
  res <- replicate(nsim, {
    bc <- rnorm(1, 0, sc); bd <- rnorm(3, 0, sdd)
    if (!is.null(true_bias)) { bc <- true_bias; bd <- c(0, 0, 0) }
    li <- log(psi) + rnorm(n, 0, tau) + bc + bd[ds] + rnorm(n, 0, sqrt(v))
    w  <- 1 / (v + tau^2)
    dm <- tapply(w * li, ds, sum) / tapply(w, ds, sum); dv <- 1 / tapply(w, ds, sum)
    if (eq_w) { est <- mean(dm); se_pat <- sqrt(sum(dv) / 9); sh <- sqrt(sc^2 + sum(sdd^2) / 9) }
    else { est <- sum(w * li) / sum(w); se_pat <- sqrt(1 / sum(w)); wd <- tapply(w, ds, sum) / sum(w); sh <- sqrt(sc^2 + sum(wd^2 * sdd^2)) }
    se <- sqrt(se_pat^2 + sh^2)
    rej <- est / se > qnorm(0.95)
    ## binding window variant: correlated estimate with slightly less information (approximation: same est + small noise)
    est_w <- est + rnorm(1, 0, se_pat * sqrt(1 / info_win - 1) * 2); rej_w <- est_w / (se * sqrt(1 / info_win)) > qnorm(0.95)
    ## conservative composition: shift down by a patient-average amount
    rej_c <- (est - shift_comp) / se > qnorm(0.95)
    cons <- sum(exp(dm) > 1) >= 2
    lo <- est - qnorm(0.975) * se; hi <- est + qnorm(0.975) * se
    c(rej = rej, all = rej & rej_w & rej_c & cons, A = rej & rej_w & rej_c & cons & lo > log(1.5),
      C = !rej & lo > log(2/3) & hi < log(1.5), se = se) })
  c(n = n, psi = psi, power_primary = mean(res["rej", ]), power_binding = mean(res["all", ]), P_A = mean(res["A", ]), P_C = mean(res["C", ]),
    se_total = mean(res["se", ])) }
POW <- list(); k_row <- 0L
prow <- function(set, vt, n, psi, sc = DBL_0B_SC, sdd = DBL_0B_SDD, tau = DBL_0B_TAU) { k_row <<- k_row + 1L
  v_t <<- vt; set.seed(DBL_SEED_0B + 10L + k_row)
  POW[[length(POW) + 1]] <<- data.table(set = set, row = k_row, t(run(n, psi, sc, sdd, tau = tau, nsim = DBL_0B_NSIM))); POW[[length(POW)]] }

## -- [5] known answers (the verdict is read only if all pass) ----
ts("[5] known answers")
tm0 <- terms(TPL[1], R0, a_s = 0, b_s = 0)
chk("doublet-rate cancellation: pi(psi) identical under Nd x0.5 and x2 when singlet rates are 0 (1e-9)", TRUE,
    all(vapply(c(0.5, 1, 2), function(p) max(abs(pi_of(terms(TPL[1], R0, a_s = 0, b_s = 0, nd_x = 0.5), p) - pi_of(tm0, p)),
                                              abs(pi_of(terms(TPL[1], R0, a_s = 0, b_s = 0, nd_x = 2), p) - pi_of(tm0, p))) < 1e-9, TRUE)))
psi_hat <- function(tm, Y, D) { p <- Y / D; num <- p * tm[["B"]] - tm[["A"]]; den <- tm[["SY"]] - p * tm[["SD"]]   # solves pi(psi) = Y / D
  ps <- ifelse(den > 0, num / den, ifelse(num > 0, Inf, 0))            # share at or above S^Y / S^D: psi-hat is +Inf
  ps[D == 0] <- NA; pmin(pmax(ps, exp(-3)), exp(3)) }                  # clipped to [e^-3, e^3]; no doublets -> NA
binsim <- function(tm, nrep = 4000) { f <- info(tm); D <- rpois(nrep, f[["ED"]]); Y <- rbinom(nrep, D, f[["pi0"]]); sd(log(psi_hat(tm, Y, D)), na.rm = TRUE) }
set.seed(DBL_SEED_0B + 1L)
tmx <- terms(TPL[1], R0) * 100
chk("Fisher SE: binomial-simulation sd of log psi-hat within 5 % of it at E[D] x100", TRUE, abs(binsim(tmx) / sqrt(info(tmx)[["var_s"]]) - 1) < 0.05)
reg <- rbind(prow("regression: design 2's v_t", c(0.194, 0.348, 0.535)^2, 12, 1.5), prow("regression: design 2's v_t", c(0.194, 0.348, 0.535)^2, 12, 2))
chk("regression: design 2's v_t reproduce the redesign's binding power 0.67 / 0.98 at n 12 (+- 0.02)", TRUE,
    abs(reg$power_binding[1] - 0.67) <= 0.02 && abs(reg$power_binding[2] - 0.98) <= 0.02)
chk("regression: ... and primary power 0.78 at psi 1.5 (+- 0.02)", TRUE, abs(reg$power_primary[1] - 0.78) <= 0.02)

## -- [6] base templates, power grid, sensitivity ----
ts("[6] templates and power")
TV <- template_var(R0, DRB)
set.seed(DBL_SEED_0B + 1L); TV[, sd_binom_over_fisher := vapply(seq_len(nrow(TPL)), function(j) binsim(terms(TPL[j], R0)) / sqrt(var_s[j]), 0)]
TV <- cbind(TPL[, .(lane, Nd, Ns)], TV)[, `:=`(se_total = sqrt(var_t), cal_share_of_var = var_c / var_t, cal_basis = DBL_0B_CAL_BASIS)]
vb <- TV$var_t
for (n in c(10L, 12L, 16L)) for (p in c(1, 1.25, 1.5, 2)) prow("base", vb, n, p)
base <- rbindlist(POW)[set == "base"]
chk("size: primary rejection at psi 1, n 12, is 0.05 +- 0.01", TRUE, abs(base[n == 12 & psi == 1, power_primary] - 0.05) <= 0.01)
SENS <- list(
  list("singlet a_s 0.001", list(a_s = 0.001, b_s = 0.001 + DBL_0B_B_S_EXTRA)), list("singlet a_s 0.004", list(a_s = 0.004, b_s = 0.004 + DBL_0B_B_S_EXTRA)),
  list("singlet a_s 0.008", list(a_s = 0.008, b_s = 0.008 + DBL_0B_B_S_EXTRA)), list("missing types x0.5", list(miss = 0.5)),
  list("missing types x2", list(miss = 2)), list("triplets 0.04", list(trip = 0.04)), list("triplets 0.14", list(trip = 0.14)))
for (s in SENS) { vt <- do.call(template_var, c(list(R0, DRB), s[[2]]))$var_t; for (p in c(1.5, 2)) prow(s[[1]], vt, DBL_0B_N, p) }
v3 <- vb; i3 <- TPL$lane == DBL_SMOKE_LANES[["pooled"]]
v3[i3] <- vapply(which(i3), function(j) { tm <- terms(TPL[j], R0, nd_x = 3); f <- info(tm)
  f[["var_s"]] + var((apply(DRB(TPL$lane[j]), 2, function(Rd) lgt(pi_of(terms(TPL[j], Rd, nd_x = 3), 1))) - lgt(f[["pi0"]])) / f[["slope"]]) }, 0)
for (p in c(1.5, 2)) prow("P01 templates with three-lane Nd", v3, DBL_0B_N, p)
for (p in c(1.5, 2)) prow("tau 0.5", vb, DBL_0B_N, p, tau = 0.5)
for (p in c(1.5, 2)) prow("bias prior 0.10 + 0.10/0.15/0.15", vb, DBL_0B_N, p, sc = 0.10, sdd = c(0.10, 0.15, 0.15))
for (p in c(1.5, 2)) prow("sampling variance only (no calibration term)", TV$var_s, DBL_0B_N, p)
for (p in c(1.5, 2)) prow("binomial-simulation variance instead of Fisher", TV$var_s * TV$sd_binom_over_fisher^2 + TV$var_c, DBL_0B_N, p)
TVo <- template_var(R0, DRO)
for (p in c(1.5, 2)) prow(paste("calibration variance on the other basis:", setdiff(names(DRS), DBL_0B_CAL_BASIS)), TVo$var_t, DBL_0B_N, p)
BF <- rbindlist(lapply(seq_len(nrow(TPL)), function(j) { tp <- TPL[j]; p0 <- pi_of(terms(tp, R0), 1); DR <- DRB(tp$lane)   # re-solve psi-hat per draw
  lp <- vapply(seq_len(min(400L, ncol(DR))), function(i) log(psi_hat(terms(tp, DR[, i]), p0, 1)), 0)
  data.table(var_c_brute = var(lp[abs(lp) < 2.999]), n_not_identified = sum(abs(lp) >= 2.999)) }))
TV <- cbind(TV, BF)
for (p in c(1.5, 2)) prow("calibration variance by re-solving psi-hat (brute force)", TV$var_s + TV$var_c_brute, DBL_0B_N, p)
POW <- rbindlist(POW)

## -- [7] verdict (frozen rule) ----
CHK <- rbindlist(CHK); VALID <- all(CHK$pass)
p2 <- POW[set == "base" & n == DBL_0B_N & psi == 2, power_binding]; p15 <- POW[set == "base" & n == DBL_0B_N & psi == 1.5, power_binding]
verdict <- if (!VALID) "NOT READ: a known answer failed" else if (p2 < DBL_0B_BAR)
  sprintf("BN-1: lineage-level T-LEUK pairing has no validated null in these data (planning power %.3f at psi 2, n = %d)", p2, DBL_0B_N) else
  sprintf("pass (planning power %.3f at psi 2, n = %d); planning label: %s", p2, DBL_0B_N,
          if (p15 >= DBL_0B_BAR) "adequately powered at psi 1.5" else sprintf("only a doubling is detectable (%.3f at psi 1.5)", p15))
if (VALID) { print(TV[, lapply(.SD, function(z) if (is.numeric(z)) signif(z, 3) else z)])   # shown only when every known answer passed
  print(POW[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)]) }
VER <- data.table(power_psi2_n12 = p2, power_psi15_n12 = p15, bar = DBL_0B_BAR, valid = VALID, verdict = verdict)
if (MOCK) { ts("[mock] ", verdict, " -- nothing written"); quit(save = "no") }
if (VALID) { fwrite(TV, OUT[["templates"]]); fwrite(POW, OUT[["power"]]) }   # a failed known answer: checks and verdict only
fwrite(VER, OUT[["verdict"]]); fwrite(CHK, OUT[["checks"]])
ts(sprintf("[done] %d checks, %d FAIL. VERDICT: %s", nrow(CHK), sum(!CHK$pass), verdict))
