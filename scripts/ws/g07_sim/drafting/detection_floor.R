# input : D2's status-free background profile (scratch_D2/real_background_profile_nonF402.rds; no F402 gene, no status)
# output: printed expected number of detected set genes per cell at 1,000 / 1,500 / 2,000 UMI (synthetic)
# what it does: NB detection probability per gene, P(y > 0) = 1 - (1 + mu*phi)^(-1/phi), phi 0.5 as calibrated in D2,
#   summed over random 156-gene sets drawn from all tested non-F402 genes (E) or from the lower half by abundance (LOW)
suppressPackageStartupMessages(library(data.table))
B <- readRDS("/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D2/real_background_profile_nonF402.rds")
str(B, max.level = 1)
p <- B$bg / sum(B$bg)                      # status-free relative abundance on the non-F402 axis (F402 rows absent)
tg <- intersect(B$tested_nonF402, names(p)); pt <- p[tg]
low <- names(pt)[pt <= median(pt)]
phi <- 0.5
det <- function(g, N) sum(1 - (1 + p[g] * N * phi)^(-1 / phi))
set.seed(20261109)
res <- rbindlist(lapply(1:500, function(r) {
  e <- sample(tg, 156); l <- sample(low, 156)
  data.table(rep = r, N = c(1000, 1500, 2000), E = sapply(c(1000, 1500, 2000), function(N) det(e, N)),
             LOW = sapply(c(1000, 1500, 2000), function(N) det(l, N))) }))
print(res[, .(E_median = median(E), E_p10 = quantile(E, 0.1), E_p90 = quantile(E, 0.9),
              LOW_median = median(LOW), LOW_p90 = quantile(LOW, 0.9)), by = N])
