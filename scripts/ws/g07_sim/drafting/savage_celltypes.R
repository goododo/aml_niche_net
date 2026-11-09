# input : Savage 2021 Table S4 (mmc5.xlsx); output: cell types x sample x time point row counts (reference file only)
suppressPackageStartupMessages({ library(readxl); library(data.table) })
X <- "/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/mmc5.xlsx"
for (s in c("Experiment 1", "Experiment 2")) { d <- as.data.table(read_excel(X, sheet = s, skip = 1)); setnames(d, 1:4, c("ct", "smp", "tp", "gene"))
  cat("\n==", s, "\n"); print(dcast(d[, .N, by = .(ct, smp, tp)], ct ~ smp + tp, value.var = "N", fill = 0)) }
