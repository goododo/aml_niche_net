#!/usr/bin/env Rscript
# input : D2 engine inventory/scratch_D2/10_engine.R (synthetic counts around a status-free non-F402 profile; exact fast
#         UCell) and real_T_percell_umi.csv.gz (QC only: healthy cells' sample, fine label, UMI). Adapted from
#         review1/L2/t6_k6b_sim_p3c.R (S5b-1000: drop < 1,000 UMI, thin the rest to exactly 1,000).
# output: review2/R2b/04_k6b_v2_{rows.csv.gz,out.txt}
# what it does: scores the DRAFT v2 GRADED K6b bars (sub-population readable if (b) "sub" >= 8/10 and (a) "sub" <= 1/10;
#   uniform readable if (a) "uniform" >= 7/10; both fail -> not readable) with (b) = the exact-match plant ("bx"), the
#   nested h_max, and the K6e detection floor (median detected up-set genes per healthy cell >= 7), in worlds:
#   E (random tested sets), MID (sets drawn below an abundance cut chosen so detection is ~7 per cell, near the floor),
#   LOW (D2 low-abundance sets), DON30 (between-donor SD 0.30 log2 on the set). No real expression, no AML cell.
# usage: Rscript 04_k6b_v2.R <reps> <worlds comma-separated: E,MID,LOW,DON15,DON30>
D2 <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D2"
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review2/R2b"
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
pa0 <- RB$bg / sum(RB$bg); edet <- function(g) sum(1 - (1 + PHI * THIN * pa0[g])^(-1 / PHI))   # expected detected genes per cell
make_mid <- function(seed, target = 7.0) { pool <- RB$tested_nonF402; ab <- RB$bg[pool]
  for (q in seq(1, 0.5, by = -0.005)) { set.seed(seed); pl <- pool[ab <= quantile(ab, q)]; g <- sample(pl, 370L)
    if (edet(g[1:156]) <= target) return(list(up = g[1:156], dn = g[157:370], q = q)) }
  stop("no cut reaches the target") }
setsM <- make_mid(11); cat(sprintf("MID cut: abundance quantile %.3f; expected detected up genes per cell at 1,000 UMI: E %.2f MID %.2f LOW %.2f\n",
  setsM$q, edet(setsR$up), edet(setsM$up), edet(setsL$up)))
allsig <- unique(c(unlist(setsR), unlist(setsL), setsM$up, setsM$dn))
p_all <- RB$bg / sum(RB$bg); p_bg <- p_all[setdiff(names(p_all), allsig)]
set.seed(SEED + 3); Lb <- rep(THIN, 50)
bank <- make_bank(Lb, p_bg, PHI, SEED + 3); bankT <- make_bank(rep(THIN, 400), p_bg, PHI, SEED + 4)
GT <- rbind(bank$GT, bankT$GT); NB <- length(Lb)
pick <- function(L) { deep <- L > THIN; id <- integer(length(L))
  if (any(!deep)) id[!deep] <- pmin(pmax(findInterval(L[!deep], Lb) + sample(-5:5, sum(!deep), TRUE), 1L), NB)
  if (any(deep)) id[deep] <- NB + sample.int(400L, sum(deep), TRUE); id }
rows <- list()
for (world in WORLDS) {
  sets <- if (world == "LOW") setsL else if (world == "MID") setsM[c("up", "dn")] else setsR
  sig <- c(sets$up, sets$dn); nU <- length(sets$up); nD <- length(sets$dn); G <- nU + nD; iU <- seq_len(nU); iD <- nU + seq_len(nD)
  p_sig <- p_all[sig]
  SIGZ <- if (world == "Z03") 0.3 else 0
  SIGD <- switch(world, DON15 = 0.15, DON30 = 0.30, 0)
  score_P3 <- function(Y, L) { deep <- rep(TRUE, length(L))
    if (any(deep)) Y[deep, ] <- matrix(rbinom(sum(deep) * G, as.vector(Y[deep, , drop = FALSE]), rep(pmin(1, THIN / L[deep]), G)), sum(deep), G)
    Le <- ifelse(deep, THIN, L); LASTDET <<- rowSums(Y[, iU, drop = FALSE] > 0); ucell_fast(Y, GT, pick(Le), list(up = iU, dn = iD), maxRank = 1500)[, 1] }
  for (rp in seq_len(REPS)) {
    set.seed(SEED + 1000L * match(world, c("E", "Z03", "DON15", "DON30", "LOW", "MID")) + rp)
    Dk <- matrix(2^rnorm(G * 9, 0, TAU), G, 9)
    beta <- 0.6 * pmin(2.5, 0.5 + rexp(G, 1 / 0.4))                    # draft K6b: 0.6 x min(2.5, 0.5 + Exp(mean 0.4))
    ud <- setNames(rnorm(NDON, 0, SIGD), donors)
    MU <- outer(H$L, p_sig) * t(Dk[, H$lab, drop = FALSE])
    fz <- 2^(ud[H$sample] + rnorm(nrow(H), 0, SIGZ)); MU[, iU] <- MU[, iU] * fz; MU[, iD] <- MU[, iD] / fz
    Y0 <- matrix(rpois(length(MU), MU * rgamma(length(MU), 1 / PHI, scale = PHI)), nrow(H), G)
    S0 <- score_P3(Y0, H$L); DET0 <- median(LASTDET)
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
          pbA, pbA_noDn, pbB, pbX, frac_b = mean(sb), det_med = DET0) }
    }
    cat(world, "rep", rp, "done", format(Sys.time()), "\n")
  }
}
R <- rbindlist(rows); fwrite(R, file.path(OUT, "04_k6b_v2_rows.csv.gz"))
K <- R[, { A <- call_nest[arm == "a"]; BX <- call_nest[arm == "bx"]; N0 <- call_nest[arm == "none"]
  .(det = det_med[1], a_U = sum(A == "uniform"), a_NS = sum(A == "not_separable"), a_S = sum(A == "sub"), a_none = sum(A == "no_shift"),
    bx_S = sum(BX == "sub"), bx_U = sum(BX == "uniform"), bx_none = sum(BX == "no_shift"), null_pass = sum(N0 != "no_shift"),
    mz_a = median(mean_z[arm == "a"]), h = median(h_nest[arm == "a"])) }, by = .(world, rep)]
K[, `:=`(floor_ok = det >= 7, sub_readable = bx_S >= 8 & a_S <= 1, uni_readable = a_U >= 7)]
K[, state := fifelse(!floor_ok, "floor fails: shape not readable", fifelse(sub_readable & uni_readable, "both readable",
           fifelse(sub_readable, "sub only (uniform -> 'uniform or not separable')", fifelse(uni_readable, "uniform only (sub-pop: no rule)", "both fail: not readable"))))]
cat("\n== per rep (nested h_max) ==\n"); print(K)
cat("\n== graded K6b states by world ==\n"); print(dcast(K[, .N, by = .(world, state)], world ~ state, value.var = "N", fill = 0))
cat("\n== world summaries ==\n")
print(K[, .(reps = .N, det_med = round(median(det), 2), floor_ok = mean(floor_ok), sub_readable = mean(sub_readable), uni_readable = mean(uni_readable),
            a_S_mean = round(mean(a_S), 2), a_S_ge1 = mean(a_S >= 1), a_U_med = median(a_U), bx_S_med = median(bx_S), null_pass_max = max(null_pass),
            h_med = round(median(h), 2), mz_a_med = round(median(mz_a), 2)), by = world])
cat("\n== pooled per-donor calls (nested) by world and arm ==\n"); print(dcast(R[, .N, by = .(world, arm, call_nest)], world + arm ~ call_nest, value.var = "N", fill = 0))
cat("\n== false 'sub' calls in arm (a) by donor (pooled reps) ==\n"); print(dcast(R[arm == "a", .(fs = sum(call_nest == "sub")), by = .(world, donor, n)], donor + n ~ world, value.var = "fs"))
