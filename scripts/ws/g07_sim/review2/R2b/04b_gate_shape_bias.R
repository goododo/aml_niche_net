#!/usr/bin/env Rscript
# input : 04_k6b_v2_rows_<world>.csv.gz (synthetic K6b replays); output: 04b_gate_shape_bias.out
# what it does: is the 3.5(e) gate shape-neutral? Compares the mean z and gate pass rate of the uniform plant (a) and the
#   exact-match sub-population plant (bx), which add the same expected counts (K6c), per world.
suppressPackageStartupMessages(library(data.table))
R <- rbindlist(lapply(c("E", "MID", "LOW", "DON15", "DON30"), function(w) fread(sprintf("04_k6b_v2_rows_%s.csv.gz", w))))
print(R[arm %in% c("a", "bx"), .(median_mean_z = round(median(mean_z), 2), gate_pass = round(mean(call_nest != "no_shift"), 3),
         k6c_absdiff_med = round(median(abs(pbA_noDn - pbX)), 3)), by = .(world, arm)][order(world, arm)])
cat("\nper donor (DON30): gate pass of a vs bx\n")
print(dcast(R[world == "DON30" & arm %in% c("a", "bx"), .(pass = mean(call_nest != "no_shift")), by = .(donor, n, arm)], donor + n ~ arm, value.var = "pass"))
