# input : sum_rules.csv, sum_healthy_lodo_max.csv; output: sum_rules_compact.csv (printed)
suppressPackageStartupMessages(library(data.table)); options(width = 250)
RW <- fread("sum_rules.csv")
g <- function(w, col) RW[world == w, c("rule", "tag", "pipe", "n", col), with = FALSE]
C <- Reduce(function(a, b) merge(a, b, by = c("rule", "tag", "pipe", "n"), all = TRUE), list(
  setnames(g("null", "no_shift"), "no_shift", "null_noshift"),
  setnames(g("uni", c("uniform", "sub", "not_separable", "no_shift")), c("uniform", "sub", "not_separable", "no_shift"), c("uni_U", "uni_S", "uni_NS", "uni_noshift")),
  setnames(g("sub20", c("sub", "uniform", "not_separable", "no_shift")), c("sub", "uniform", "not_separable", "no_shift"), c("s20_S", "s20_U", "s20_NS", "s20_noshift")),
  setnames(g("sub10", c("sub", "uniform")), c("sub", "uniform"), c("s10_S", "s10_U")),
  setnames(g("sub40", c("sub", "uniform")), c("sub", "uniform"), c("s40_S", "s40_U"))))
fwrite(C, "sum_rules_compact.csv")
print(C[pipe %in% c("P1", "P2", "P3")][order(rule, pipe, tag, n)])
H <- fread("sum_healthy_lodo_max.csv"); print(H[set == "up"])
