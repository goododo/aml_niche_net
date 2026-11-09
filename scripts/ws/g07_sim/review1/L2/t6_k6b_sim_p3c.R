#!/usr/bin/env Rscript
# input : D2 engine inventory/scratch_D2/10_engine.R (synthetic gamma-Poisson counts around a status-free NON-F402
#         profile; exact fast UCell) and real_T_percell_umi.csv.gz (QC only: healthy cells' sample, fine label, UMI)
# output: review1/L2/t6_k6b_sim_p3c_rows.csv.gz, t6_k6b_sim_p3c_out.txt (printed summary)
# variant of t2_k6b_sim.R: pipeline P3c = drop cells below 1,000 UMI, thin every remaining cell to 1,000 UMI (bank at 1,000).
# what it does: synthetic replay of the draft's K6b/K6c (PREREGISTRATION_A1_G07.md:585-586) on synthetic cells that
#   carry the 10 real healthy donors' (label, UMI) structure. Plants as written: (a) uniform additive Poisson
#   mu_hat*(2^b-1) on up genes + down genes thinned at 2^-b; (b) Bernoulli-0.2 cells get mu_hat*(2^b-1)/0.2;
#   (bx) variant: exactly round(n/5) cells with a per-gene multiplier that matches (a)'s expected increment.
#   mu_hat = cell UMI x pooled label profile of the other nine donors (synthetic counts). P3 = binomial thinning of
#   cells above 1,500 UMI. Reference = other nine; h_max two ways: nested (each e != d vs the other eight) and
#   non-nested (main-run LODO values of the other nine). Rule = draft 3.5(e) = D2 13_sim.R pstats + 17_final_rule.R.
#   No real expression value, no F402 value, no AML cell is used.
# usage: Rscript t2_k6b_sim.R <reps> <worlds comma-separated: E,Z03,DON15,DON30,LOW>
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })
D2 <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D2"
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L2"
owd <- setwd(D2); source("10_engine.R"); setwd(OUT)          # engine reads its inputs relative to D2; writes nothing
a <- commandArgs(TRUE); REPS <- as.integer(a[1]); WORLDS <- strsplit(a[2], ",")[[1]]
SEED <- 20261111L; PHI <- 0.5; TAU <- 0.5; THIN <- 1000
# ---- verbatim from D2 13_sim.R:44-67 (mk_ref, zf, f0, grids, pstats) ----
mk_ref <- function(s, lab, L, adj) {   # per label: z transform and z95; pooled z density
  ks <- sort(unique(lab)); r <- list(adj = adj, k = ks, a = numeric(9), b = numeric(9), sd = rep(1, 9), z95 = rep(Inf, 9))
  for (k in ks) { j <- lab == k
    if (adj && sum(j) > 10) { f <- lm(s[j] ~ log10(L[j])); r$a[k] <- coef(f)[1]; r$b[k] <- coef(f)[2]; e <- resid(f) } else { r$a[k] <- mean(s[j]); r$b[k] <- 0; e <- s[j] - r$a[k] }
    r$sd[k] <- max(sd(e), 1e-6); r$z95[k] <- quantile(e / r$sd[k], 0.95, type = 7) }
  z <- zf(r, s, lab, L); r$zr <- z; r$z95c <- r$z95[lab]
  d <- density(z, n = 2048, from = min(z) - 3, to = max(z) + 8); r$dx <- d$x; r$dy <- pmax(d$y, 1e-8)
  r$m50 <- median(z); r$q9010 <- diff(quantile(z, c(0.1, 0.9))); r$v <- var(z); r }
zf <- function(r, s, lab, L) (s - r$a[lab] - r$b[lab] * log10(L)) / r$sd[lab]
f0 <- function(r, x) approx(r$dx, r$dy, x, rule = 2)$y
PIG <- seq(0.1, 1, 0.1); DG <- seq(0, 20, 0.5); SG <- seq(-3, 20, 1)
pstats <- function(r, s, lab, L, mix = TRUE) {
  z <- zf(r, s, lab, L); n <- length(z); mz <- mean(z)
  F95 <- mean(z > r$z95[lab]); F95e <- mean(r$zr + mz > r$z95c)
  m3 <- mean((z - mz)^3) / sd(z)^3; m4 <- mean((z - mz)^4) / sd(z)^4 - 3
  BC <- if (n > 3) (m3^2 + 1) / (m4 + 3 * (n - 1)^2 / ((n - 2) * (n - 3))) else NA
  llU <- NA; llM <- NA; pih <- NA
  if (mix) { lu <- function(d) sum(log(f0(r, z - d))); g0 <- SG[which.max(sapply(SG, lu))]
    llU <- optimize(lu, c(g0 - 1, g0 + 1), maximum = TRUE)$objective
    F0 <- f0(r, z); FD <- sapply(DG, function(d) f0(r, z - d)); llM <- -Inf; pih <- 1
    for (p in PIG) { ll <- colSums(log((1 - p) * F0 + p * FD)); if (max(ll) > llM + 1e-9) { llM <- max(ll); pih <- p } }
    if (llM - llU < 1e-6) pih <- 1 }
  c(mean_z = mz, F95 = F95, ET = F95 - F95e, logVR = log(var(z) / r$v), D50 = median(z) - r$m50,
    QS = unname(diff(quantile(z, c(0.1, 0.9))) - r$q9010), BC = BC, LR = if (mix) max(0, llM - llU) else NA, pihat = pih) }
# ---- light reference (mean/SD per label only; same arithmetic as mk_ref's non-adjusted branch) ----
mk_light <- function(s, lab) { r <- list(a = numeric(9), b = numeric(9), sd = rep(1, 9))
  for (k in sort(unique(lab))) { j <- lab == k; r$a[k] <- mean(s[j]); r$sd[k] <- max(sd(s[j] - r$a[k]), 1e-6) }; r }
call_rule <- function(st, n, h) { se <- sqrt(exp(st[["logVR"]]) / n)
  if (!(st[["mean_z"]] - 1.96 * se > h)) return("no_shift")
  if (st[["pihat"]] <= 0.5 && st[["LR"]] >= 10) "sub" else if (st[["pihat"]] >= 0.9 || st[["LR"]] < 2) "uniform" else "not_separable" }
# ---- synthetic healthy cells with the real healthy (label, UMI) structure ----
H <- PC[status != "AML", .(sample, lab = match(lab, LAB), L = umi_axis)]; H <- H[L >= THIN]; cat("P3c healthy cells kept:", nrow(H), "\n")
donors <- unique(H$sample); NDON <- length(donors)
setsR <- make_sets(11, mode = "random"); setsL <- make_sets(11, mode = "low")
allsig <- unique(c(unlist(setsR), unlist(setsL)))
p_all <- RB$bg / sum(RB$bg); p_bg <- p_all[setdiff(names(p_all), allsig)]
set.seed(SEED + 3); Lb <- rep(THIN, 50)
bank <- make_bank(Lb, p_bg, PHI, SEED + 3); bankT <- make_bank(rep(THIN, 400), p_bg, PHI, SEED + 4)
GT <- rbind(bank$GT, bankT$GT); NB <- length(Lb)
pick <- function(L) { deep <- L > THIN; id <- integer(length(L))
  if (any(!deep)) id[!deep] <- pmin(pmax(findInterval(L[!deep], Lb) + sample(-5:5, sum(!deep), TRUE), 1L), NB)
  if (any(deep)) id[deep] <- NB + sample.int(400L, sum(deep), TRUE); id }
rows <- list()
for (world in WORLDS) {
  sets <- if (world == "LOW") setsL else setsR
  sig <- c(sets$up, sets$dn); nU <- length(sets$up); nD <- length(sets$dn); G <- nU + nD; iU <- seq_len(nU); iD <- nU + seq_len(nD)
  p_sig <- p_all[sig]
  SIGZ <- if (world == "Z03") 0.3 else 0
  SIGD <- switch(world, DON15 = 0.15, DON30 = 0.30, 0)
  score_P3 <- function(Y, L) { deep <- rep(TRUE, length(L))
    if (any(deep)) Y[deep, ] <- matrix(rbinom(sum(deep) * G, as.vector(Y[deep, , drop = FALSE]), rep(pmin(1, THIN / L[deep]), G)), sum(deep), G)
    Le <- ifelse(deep, THIN, L); ucell_fast(Y, GT, pick(Le), list(up = iU, dn = iD), maxRank = 1500)[, 1] }
  for (rp in seq_len(REPS)) {
    set.seed(SEED + 1000L * match(world, c("E", "Z03", "DON15", "DON30", "LOW")) + rp)
    Dk <- matrix(2^rnorm(G * 9, 0, TAU), G, 9)
    beta <- 0.6 * pmin(2.5, 0.5 + rexp(G, 1 / 0.4))                    # draft K6b: 0.6 x min(2.5, 0.5 + Exp(mean 0.4))
    ud <- setNames(rnorm(NDON, 0, SIGD), donors)
    MU <- outer(H$L, p_sig) * t(Dk[, H$lab, drop = FALSE])
    fz <- 2^(ud[H$sample] + rnorm(nrow(H), 0, SIGZ)); MU[, iU] <- MU[, iU] * fz; MU[, iD] <- MU[, iD] / fz
    Y0 <- matrix(rpois(length(MU), MU * rgamma(length(MU), 1 / PHI, scale = PHI)), nrow(H), G)
    S0 <- score_P3(Y0, H$L)
    jdl <- lapply(donors, function(d) which(H$sample == d)); names(jdl) <- donors
    # main-run LODO values (each donor vs the other nine): used for the non-nested h_max
    lodo10 <- sapply(donors, function(e) { je <- jdl[[e]]; r <- mk_light(S0[-je], H$lab[-je]); mean(zf(r, S0[je], H$lab[je], 1)) })
    for (d in donors) {
      jd <- jdl[[d]]; n <- length(jd); oth <- setdiff(seq_len(nrow(H)), jd)
      nest <- sapply(setdiff(donors, d), function(e) { je <- jdl[[e]]; keep <- setdiff(oth, je)
        r <- mk_light(S0[keep], H$lab[keep]); mean(zf(r, S0[je], H$lab[je], 1)) })
      h_nest <- max(nest); h_non <- max(lodo10[setdiff(donors, d)])
      R9 <- mk_ref(S0[oth], H$lab[oth], H$L[oth], FALSE)
      prof <- sapply(1:9, function(k) { jj <- oth[H$lab[oth] == k]; if (length(jj)) colSums(Y0[jj, iU, drop = FALSE]) / sum(H$L[jj]) else rep(0, nU) })
      muhat <- H$L[jd] * t(prof[, H$lab[jd], drop = FALSE]); fu <- 2^beta[iU] - 1
      Y0d <- Y0[jd, , drop = FALSE]; L0 <- H$L[jd]
      # (a) uniform
      Ya <- Y0d; addA <- matrix(rpois(n * nU, sweep(muhat, 2, fu, `*`)), n, nU); Ya[, iU] <- Ya[, iU] + addA
      Ya[, iD] <- matrix(rbinom(n * nD, as.vector(Ya[, iD]), rep(2^-beta[iD], each = n)), n, nD)
      La <- L0 + rowSums(addA) - (rowSums(Y0d[, iD, drop = FALSE]) - rowSums(Ya[, iD, drop = FALSE]))
      La_noDn <- L0 + rowSums(addA)
      # (b) Bernoulli 0.2, multiplier (2^b-1)/0.2 as written
      sb <- runif(n) < 0.2; Yb <- Y0d; addB <- matrix(0, n, nU)
      if (any(sb)) addB[sb, ] <- matrix(rpois(sum(sb) * nU, sweep(muhat[sb, , drop = FALSE], 2, fu / 0.2, `*`)), sum(sb), nU)
      Yb[, iU] <- Yb[, iU] + addB; Lbb <- L0 + rowSums(addB)
      # (bx) exactly round(n/5) cells, per-gene multiplier matching (a)'s expected increment in this donor
      sx <- seq_len(n) %in% sample.int(n, round(n / 5)); mx <- fu * colSums(muhat) / colSums(muhat[sx, , drop = FALSE])
      Yx <- Y0d; addX <- matrix(0, n, nU); addX[sx, ] <- matrix(rpois(sum(sx) * nU, sweep(muhat[sx, , drop = FALSE], 2, mx, `*`)), sum(sx), nU)
      Yx[, iU] <- Yx[, iU] + addX; Lx <- L0 + rowSums(addX)
      pb0 <- sum(Y0d[, iU]) / sum(L0)
      pbA <- log2(sum(Ya[, iU]) / sum(La) / pb0); pbA_noDn <- log2(sum(Ya[, iU]) / sum(La_noDn) / pb0)
      pbB <- log2(sum(Yb[, iU]) / sum(Lbb) / pb0); pbX <- log2(sum(Yx[, iU]) / sum(Lx) / pb0)
      arms <- list(none = list(S0[jd], L0), a = list(score_P3(Ya, La), La), b = list(score_P3(Yb, Lbb), Lbb), bx = list(score_P3(Yx, Lx), Lx))
      for (am in names(arms)) {
        st <- pstats(R9, arms[[am]][[1]], H$lab[jd], arms[[am]][[2]], TRUE)
        rows[[length(rows) + 1]] <- data.table(world, rep = rp, donor = d, n, arm = am, mean_z = st[["mean_z"]], se = sqrt(exp(st[["logVR"]]) / n),
          LR = st[["LR"]], pihat = st[["pihat"]], h_nest, h_non, h_max10 = max(lodo10), lodo_own = lodo10[[d]],
          call_nest = call_rule(st, n, h_nest), call_non = call_rule(st, n, h_non),
          pbA, pbA_noDn, pbB, pbX, frac_b = mean(sb)) }
    }
    cat(world, "rep", rp, "done", format(Sys.time()), "\n")
  }
}
R <- rbindlist(rows); fwrite(R, file.path(OUT, "t6_k6b_sim_p3c_rows.csv.gz"))
# ---- K6b / K6c as written, per world x rep (h_max nested; non-nested in brackets) ----
k6 <- function(R, cc) R[, {
  A <- get(cc)[arm == "a"]; B <- get(cc)[arm == "b"]; N0 <- get(cc)[arm == "none"]; BX <- get(cc)[arm == "bx"]
  .(a_U_or_NS = sum(A %in% c("uniform", "not_separable")), a_U = sum(A == "uniform"), a_S = sum(A == "sub"),
    b_S = sum(B == "sub"), bx_S = sum(BX == "sub"), none_noshift = sum(N0 == "no_shift"),
    k6c_med = median(abs(pbA[arm == "a"] - pbB[arm == "a"])), k6c_med_noDn = median(abs(pbA_noDn[arm == "a"] - pbB[arm == "a"])),
    k6c_med_bx = median(abs(pbA_noDn[arm == "a"] - pbX[arm == "a"])), k6c_max = max(abs(pbA[arm == "a"] - pbB[arm == "a"]))) }, by = .(world, rep)]
for (cc in c("call_nest", "call_non")) {
  K <- k6(R, cc)
  K[, pass_a := a_U_or_NS >= 9 & a_U >= 7 & a_S == 0][, pass_b := b_S >= 8][, pass_null := none_noshift >= 9][, pass_k6c := k6c_med <= 0.05]
  K[, pass_all := pass_a & pass_b & pass_null & pass_k6c]
  cat("\n== K6b/K6c per rep, h_max =", cc, "==\n"); print(K)
  cat("\n== share of reps passing, by world (", cc, ") ==\n")
  print(K[, .(reps = .N, pass_a = mean(pass_a), pass_b = mean(pass_b), pass_null = mean(pass_null), pass_k6c = mean(pass_k6c), pass_all = mean(pass_all),
             min_none_noshift = min(none_noshift), k6c_med_median = round(median(k6c_med), 3), k6c_noDn = round(median(k6c_med_noDn), 3), k6c_bx = round(median(k6c_med_bx), 3)), by = world])
}
cat("\n== per-donor calls pooled over reps (nested h_max) ==\n")
print(dcast(R[, .N, by = .(world, arm, call_nest)], world + arm ~ call_nest, value.var = "N", fill = 0))
cat("\n== h_max values (median over reps) and per-donor K6c |a-b| by donor size ==\n")
print(R[arm == "a", .(h_nest = round(median(h_nest), 3), h_max10 = round(median(h_max10), 3), absdiff_med = round(median(abs(pbA - pbB)), 3), frac_b_sd = round(sd(frac_b), 3)), by = .(world, donor, n)][order(world, n)])
