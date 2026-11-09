# input : 10_engine.R; output: validate_ucell.csv, calib_phi.csv
# what it does: (1) calibrates the gamma overdispersion phi so that synthetic nFeature matches real nFeature at matched UMI
#   (QC numbers only); (2) checks the fast UCell against UCell::ScoreSignatures_UCell on explicit synthetic matrices.
source("10_engine.R"); suppressPackageStartupMessages({ library(Matrix); library(UCell) })
p_all <- RB$bg / sum(RB$bg)
# (1) phi calibration on nFeature (all genes, synthetic) vs real nFeature in UMI bins
bins <- c(1000, 1500, 2000, 3000, 4500)
real <- sapply(bins, function(b) median(PC[abs(umi_axis - b) / b < 0.05]$nfeat_axis))
res <- rbindlist(lapply(c(0, 0.5, 1, 2, 4), function(phi) { set.seed(1)
  data.table(phi = phi, umi = bins, real_nfeat = real, sim_nfeat = sapply(bins, function(L) median(replicate(40, {
    lam <- L * p_all; if (phi > 0) lam <- lam * rgamma(length(lam), 1 / phi, scale = phi); sum(rpois(length(lam), lam) > 0) })))) }))
res[, ratio := round(sim_nfeat / real_nfeat, 3)]; print(res); fwrite(res, "calib_phi.csv")
