#!/usr/bin/env Rscript
# input : nothing (synthetic)
# output: review2/R2b/07_sim_rhoX.out
# what it does: donor-level toy model of 3.3(b). Per donor a set score in T and in X (shared donor noise r), AML shift
#   mu_T in T and mu_X in X. s = AML mean - healthy mean on the paired donors. rho_X = s_X / s_T, 95 % percentile
#   interval from 1,000 unstratified donor bootstraps (as written: "donor bootstrap"). Reading as 3.3(b): interaction
#   (T - X) z > 1.96 -> "larger in T"; < -1.96 -> "larger in X"; else "shared" if rho >= 0.5; else "not resolved".
#   Only sets FLAGGED in T are read (3.3(b) "read for flagged sets"): flag = z_T >= zflag, z_T on all 25 donors, of which
#   the paired donors are a subset (winner's curse carries over). No edgeR; this isolates the ratio and the selection.
set.seed(491638 + 9707)
sim <- function(nA, nH, nAx, nHx, zt, kx, r, sdx, zflag, nsim = 1500, B = 1000) {
  se_full <- sqrt(1 / nA + 1 / nH); mu_T <- zt * se_full; mu_X <- kx * mu_T
  out <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    u <- rnorm(nA + nH); aT <- mu_T * rep(1:0, c(nA, nH)) + sqrt(r) * u + sqrt(1 - r) * rnorm(nA + nH)
    zT <- (mean(aT[1:nA]) - mean(aT[nA + 1:nH])) / se_full
    if (zT < zflag) { out[[i]] <- NULL; next }
    idx <- c(seq_len(nAx), nA + seq_len(nHx))                                  # paired donors = a subset of the T donors
    aX <- mu_X * rep(1:0, c(nA, nH)) + sdx * (sqrt(r) * u + sqrt(1 - r) * rnorm(nA + nH))
    st <- aT[idx]; sx <- aX[idx]; g <- rep(1:0, c(nAx, nHx))
    sT <- mean(st[g == 1]) - mean(st[g == 0]); sX <- mean(sx[g == 1]) - mean(sx[g == 0]); rho <- sX / sT
    d <- st - sx; sed <- sqrt(var(d[g == 1]) / nAx + var(d[g == 0]) / nHx); zi <- (mean(d[g == 1]) - mean(d[g == 0])) / sed
    bs <- replicate(B, { j <- sample.int(length(g), replace = TRUE); gj <- g[j]
      if (length(unique(gj)) < 2) return(c(NA, NA)); tt <- mean(st[j][gj == 1]) - mean(st[j][gj == 0]); c(tt, (mean(sx[j][gj == 1]) - mean(sx[j][gj == 0])) / tt) })
    ci <- quantile(bs[2, ], c(0.025, 0.975), na.rm = TRUE)
    rd <- if (zi > 1.96) "larger_T" else if (zi < -1.96) "larger_X" else if (rho >= 0.5) "shared" else "not_resolved"
    out[[i]] <- data.frame(zT, sT, rho, lo = ci[1], hi = ci[2], negden = mean(bs[1, ] <= 0, na.rm = TRUE), rd) }
  do.call(rbind, out) }
cfg <- expand.grid(X = c("B", "NK"), zt = c(3, 4.5), kx = c(1, 0.5, 0), r = 0.5, stringsAsFactors = FALSE)
for (q in seq_len(nrow(cfg))) { cc <- cfg[q, ]
  pr <- if (cc$X == "B") list(nAx = 13, nHx = 10, sdx = 1.2) else list(nAx = 12, nHx = 9, sdx = 1.4)
  R <- sim(15, 10, pr$nAx, pr$nHx, cc$zt, cc$kx, cc$r, pr$sdx, zflag = 3)
  tb <- table(factor(R$rd, c("larger_T", "shared", "not_resolved", "larger_X"))) / nrow(R)
  cat(sprintf("X=%-2s true zT %.1f, X/T %.1f | flagged %4d | median rho %.2f | CI width median %.2f, P90 %.2f | CI hi > 3: %.2f | CI lo < 0: %.2f | boot s_T<=0 median %.3f | readings: %s\n",
      cc$X, cc$zt, cc$kx, nrow(R), median(R$rho), median(R$hi - R$lo), quantile(R$hi - R$lo, 0.9), mean(R$hi > 3), mean(R$lo < 0), median(R$negden),
      paste(sprintf("%s %.2f", names(tb), tb), collapse = ", "))) }
