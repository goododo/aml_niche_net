#!/usr/bin/env Rscript
# input : 10_engine.R (synthetic engine; real QC numbers and a status-free non-F402 profile only)
# output: sim_patients_<tag>.csv.gz (one row per simulated patient x pipeline x set), sim_meta_<tag>.csv
# what it does: D2 sections 3-4. Synthetic AML patients (real AML label/UMI templates) in worlds null / uniform / sub10 /
#   sub20 / sub40, matched on expected pseudobulk effect of the up set; healthy reference = the 10 real healthy donors'
#   label/UMI structure, simulated with no effect. Pipelines: P1 registered (UCell maxRank 1500, per-label healthy z/P95),
#   P2 P1 + per-label depth adjustment (linear in log10 UMI, fitted on healthy), P3 thinned to 1,500 UMI then P1,
#   P4 maxRank 3000. Per-patient statistics against the healthy reference (LODO for healthy donors).
# usage: Rscript 13_sim.R <tag> <effect_scale> <sigma_z> <set_mode> <n_list> <reps_list>
source("10_engine.R")
a <- commandArgs(TRUE); TAG <- a[1]; ESC <- as.numeric(a[2]); SIGZ <- as.numeric(a[3]); SMODE <- a[4]
NS <- as.integer(strsplit(a[5], ",")[[1]]); REPS <- as.integer(strsplit(a[6], ",")[[1]])
PHI <- 0.5; TAU <- 0.5; THIN <- 1500; SEED <- 20261109L
sets <- make_sets(11, mode = SMODE); sig <- c(sets$up, sets$dn); nU <- length(sets$up); nD <- length(sets$dn); G <- nU + nD
iU <- seq_len(nU); iD <- nU + seq_len(nD)
p_all <- RB$bg / sum(RB$bg); p_bg <- p_all[setdiff(names(p_all), sig)]; p_sig <- p_all[sig]
set.seed(SEED + 1); D <- matrix(2^rnorm(G * 9, 0, TAU), G, 9)                       # subtype baselines of set genes
beta <- pmin(c(0.5 + rexp(nU, 1 / 0.4), 0.5 + rexp(nD, 1 / 0.4)), 2.5) * ESC        # synthetic |log2FC|, not S5's
F_uni <- c(2^beta[iU], 2^-beta[iD])
set.seed(SEED + 2); Lb <- sort(sample(PC$umi_axis, 3000, replace = TRUE))
bank <- make_bank(Lb, p_bg, PHI, SEED + 3); bankT <- make_bank(rep(THIN, 600), p_bg, PHI, SEED + 4)
GT <- rbind(bank$GT, bankT$GT); LB <- c(Lb, rep(THIN, 600))
pick_bank <- function(L, thin) {
  id <- pmin(pmax(findInterval(L, Lb) + sample(-5:5, length(L), TRUE), 1L), 3000L)
  if (thin) { deep <- L > THIN; id[deep] <- 3000L + sample.int(600L, sum(deep), TRUE) }
  id }
sim_scores <- function(lab, L, world, pi_aff, thin, maxRanks) {
  n <- length(L); id <- pick_bank(L, thin); Le <- LB[id]
  MU <- outer(Le, p_sig) * t(D[, lab, drop = FALSE])
  if (SIGZ > 0) { z <- 2^rnorm(n, 0, SIGZ); MU[, iU] <- MU[, iU] * z; MU[, iD] <- MU[, iD] / z }
  aff <- rep(FALSE, n)
  if (world == "uni") MU <- sweep(MU, 2, F_uni, `*`)
  if (world == "sub") { aff <- runif(n) < pi_aff
    if (any(aff)) { m <- c(1 + (F_uni[iU] - 1) / pi_aff, pmax(0, 1 - (1 - F_uni[iD]) / pi_aff)); MU[aff, ] <- sweep(MU[aff, , drop = FALSE], 2, m, `*`) } }
  Y <- matrix(rpois(n * G, MU * rgamma(n * G, 1 / PHI, scale = PHI)), n, G)
  out <- lapply(maxRanks, function(mr) ucell_fast(Y, GT, id, list(up = iU, dn = iD), maxRank = mr))
  list(sc = out, L = Le, aff = aff, det_up = rowSums(Y[, iU] > 0), det_dn = rowSums(Y[, iD] > 0), sumU = colSums(Y[, iU]), Ltot = sum(Le)) }
labi <- function(x) match(x, LAB)
# ---------- healthy reference: 10 real donors' (label, UMI) structure, no effect ----------
H <- PC[status != "AML", .(sample, lab = labi(lab), L = umi_axis)]
set.seed(SEED + 5)
hr <- sim_scores(H$lab, H$L, "null", 0, FALSE, c(1500, 3000)); ht <- sim_scores(H$lab, H$L, "null", 0, TRUE, 1500)
HS <- list(P1 = hr$sc[[1]], P2 = hr$sc[[1]], P3 = ht$sc[[1]], P4 = hr$sc[[2]]); HL <- list(P1 = hr$L, P2 = hr$L, P3 = ht$L, P4 = hr$L)
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
PIPE <- c("P1", "P2", "P3", "P4")
REF <- lapply(PIPE, function(p) lapply(c(up = 1, dn = 2, comb = 3), function(si) {
  s <- if (si == 3) HS[[p]][, 1] - HS[[p]][, 2] else HS[[p]][, si]; mk_ref(s, H$lab, HL[[p]], p == "P2") })); names(REF) <- PIPE
rows <- list()
# healthy donors, leave-one-donor-out reference
for (d in unique(H$sample)) { j <- H$sample == d
  for (p in PIPE) for (si in 1:3) {
    s <- if (si == 3) HS[[p]][, 1] - HS[[p]][, 2] else HS[[p]][, si]
    r <- mk_ref(s[!j], H$lab[!j], HL[[p]][!j], p == "P2")
    rows[[length(rows) + 1]] <- data.table(t(pstats(r, s[j], H$lab[j], HL[[p]][j], si == 1)), world = "healthy_lodo", n = sum(j), rep = 0L, donor = d, pipe = p, set = c("up", "dn", "comb")[si]) } }
# AML patients from real AML templates
A <- PC[status == "AML", .(sample, lab = labi(lab), L = umi_axis)]; TPL <- unique(A$sample)
WORLDS <- list(null = c("null", 0), uni = c("uni", 0), sub10 = c("sub", 0.1), sub20 = c("sub", 0.2), sub40 = c("sub", 0.4))
pb <- list(); set.seed(SEED + 6)
for (ni in seq_along(NS)) for (w in names(WORLDS)) for (rp in seq_len(REPS[ni])) {
  n <- NS[ni]; tp <- sample(TPL, 1); x <- A[sample == tp][sample.int(sum(A$sample == tp), n, TRUE)]
  ww <- WORLDS[[w]]; sr <- sim_scores(x$lab, x$L, ww[1], as.numeric(ww[2]), FALSE, c(1500, 3000)); st <- sim_scores(x$lab, x$L, ww[1], as.numeric(ww[2]), TRUE, 1500)
  SC <- list(P1 = sr$sc[[1]], P2 = sr$sc[[1]], P3 = st$sc[[1]], P4 = sr$sc[[2]]); LL <- list(P1 = sr$L, P2 = sr$L, P3 = st$L, P4 = sr$L)
  for (p in PIPE) for (si in 1:3) { s <- if (si == 3) SC[[p]][, 1] - SC[[p]][, 2] else SC[[p]][, si]
    rows[[length(rows) + 1]] <- data.table(t(pstats(REF[[p]][[si]], s, x$lab, LL[[p]], si == 1)), world = w, n = n, rep = rp, donor = tp, pipe = p, set = c("up", "dn", "comb")[si]) }
  if (rp <= 20) pb[[length(pb) + 1]] <- data.table(world = w, n = n, rep = rp, sumU = sum(sr$sumU), Ltot = sr$Ltot, aff = mean(sr$aff), det_up = median(sr$det_up), det_dn = median(sr$det_dn))
}
R <- rbindlist(rows); fwrite(R, sprintf("sim_patients_%s.csv.gz", TAG))
PBt <- rbindlist(pb); hU <- sum(hr$sumU); hL <- hr$Ltot
PBt[, pb_log2fc_up := log2((sumU / Ltot) / (hU / hL))]
meta <- PBt[, .(pb_log2fc_up_mean = mean(pb_log2fc_up), pb_log2fc_up_sd = sd(pb_log2fc_up), aff = mean(aff), det_up_median = median(det_up), det_dn_median = median(det_dn)), by = .(world, n)]
meta <- rbind(meta, data.table(world = "healthy_ref", n = nrow(H), pb_log2fc_up_mean = 0, pb_log2fc_up_sd = NA, aff = 0, det_up_median = median(hr$det_up), det_dn_median = median(hr$det_dn)))
meta[, `:=`(tag = TAG, effect_scale = ESC, sigma_z = SIGZ, set_mode = SMODE, mean_beta_up = mean(beta[iU]), mean_beta_dn = mean(beta[iD]))]
fwrite(meta, sprintf("sim_meta_%s.csv", TAG)); print(meta)
cat("healthy synthetic cells: detected up-set genes per cell median", median(hr$det_up), " IQR", quantile(hr$det_up, c(.25, .75)), "; down-set", median(hr$det_dn), " IQR", quantile(hr$det_dn, c(.25, .75)), "\n")
