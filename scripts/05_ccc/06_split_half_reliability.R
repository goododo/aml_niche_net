# 06_split_half_reliability.R ----
# INPUT  : results/tables/08_scoring/paired_roster.csv  (37 paired samples, written by
#          08_scoring/13_gw_blindness.py -- the pairing rule lives there, not here)
#          CCC_BMM_DIR/<ds>/<sample>__bmm_percell.csv   (per-cell hierarchy_bin + QC flags)
# OUTPUT : DIR_CCC/split_half/<ds>__<sample>__<A|B>.csv   (one 'cell' column per half)
#          DIR_CCC/split_half/tasks.csv                   (array manifest: row, dataset, sample, half)
# WHAT IT DOES : splits each paired sample's graph-eligible cells into two halves, stratified by
#          hierarchy_bin, so 02_run_cellchat.R --cell_subset can be run twice per sample and the
#          two C matrices compared. It does NOT call CellChat.
#
# WHY. 08_scoring D4 measured, inside the paired samples, that within-patient |dC| is as large as
# the between-patient SD (median snr 1.11 for Dx->Treatment) while the shared timepoint effect is
# only 3-5% of the variance. Large but directionless change is what real idiosyncratic biology
# looks like AND what a noisy estimate looks like. Nothing in this pipeline separates the two.
# Splitting a sample's cells in half and re-running the identical CellChat path gives the noise
# floor directly: two halves of ONE sample differ by measurement error alone.
#
# THE SPLIT IS ASSIGNMENT-BASED, NOT SELECTION-BASED. Each cell is assigned to a half up front,
# from the per-cell table. 02_run_cellchat.R then intersects with whatever cells the QC object
# actually holds. Because the assignment does not depend on which cells survive that intersection,
# each bin stays near 50/50 whatever the QC object drops, and the run prints the counts it kept.
#
# READ THE HALVES AS A LOWER BOUND. Each half has half the cells, so it detects fewer significant
# LR pairs than the full sample does. Agreement measured here therefore UNDERSTATES the reliability
# of the production C. 08_scoring/13 restricts the comparison to edges whose sender and receiver
# both clear CCC_MIN_CELLS_PER_NODE in BOTH halves, and reports the dropout separately, so the
# threshold effect is not read as noise.
#
# Usage : Rscript scripts/05_ccc/06_split_half_reliability.R
suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))

SPLIT_SEED <- 20260901L          # fixed, so the halves are reproducible run to run
OUT_DIR    <- file.path(DIR_CCC, "split_half")
ROSTER     <- here::here("results", "tables", "08_scoring", "paired_roster.csv")

stopifnot(file.exists(ROSTER))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ros <- fread(ROSTER)
message("[0] roster: ", nrow(ros), " paired samples from ", ROSTER)

tasks <- list()
skipped <- character(0)
for (i in seq_len(nrow(ros))) {
  ds <- ros$dataset[i]; smp <- ros$sample[i]
  f <- file.path(CCC_BMM_DIR, ds, paste0(smp, "__bmm_percell.csv"))
  if (!file.exists(f)) { skipped <- c(skipped, paste0(ds, "/", smp, " (no bmm csv)")); next }
  b <- fread(f, select = c("cell", "hierarchy_bin", "in_ccc_graph", "high_error"))

  # the SAME admission filter 02_run_cellchat.R applies, so a half is a half of that population
  b <- b[as.logical(in_ccc_graph) %in% TRUE &
         as.logical(high_error)   %in% FALSE &
         hierarchy_bin %in% CCC_NODES]
  if (nrow(b) == 0L) { skipped <- c(skipped, paste0(ds, "/", smp, " (0 eligible cells)")); next }

  # stratify by bin: shuffle within bin, then alternate A/B. Alternating (rather than sampling
  # n/2) keeps the two halves within one cell of each other in EVERY bin, including bins with
  # an odd count, which matters because the rare bins are the ones that fall off the threshold.
  set.seed(SPLIT_SEED + i)
  b[, ord := sample(.N), by = hierarchy_bin]
  b[, half := fifelse(ord %% 2L == 1L, "A", "B")]

  for (h in c("A", "B")) {
    out <- file.path(OUT_DIR, sprintf("%s__%s__%s.csv", ds, smp, h))
    fwrite(b[half == h, .(cell)], out)
    tasks[[length(tasks) + 1L]] <- data.table(dataset = ds, sample = smp, half = h,
                                              n_cells = b[half == h, .N], subset_file = out)
  }
  chk <- dcast(b[, .N, by = .(hierarchy_bin, half)], hierarchy_bin ~ half, value.var = "N", fill = 0L)
  worst <- max(abs(chk$A - chk$B))
  if (worst > 1L) stop("stratified split is off by ", worst, " cells in a bin for ", ds, "/", smp,
                       " -- the alternation is broken, do not launch the array")
  message("  [split] ", ds, "/", smp, "  ", b[, .N], " cells -> A ", b[half == "A", .N],
          " / B ", b[half == "B", .N], "  (max per-bin imbalance ", worst, ")")
}

if (length(skipped)) {
  message("\n[!] skipped ", length(skipped), " sample(s):")
  for (s in skipped) message("      ", s)
}
T <- rbindlist(tasks)
T[, row := .I]
setcolorder(T, c("row", "dataset", "sample", "half", "n_cells", "subset_file"))
fwrite(T, file.path(OUT_DIR, "tasks.csv"))

message("\n[done] ", nrow(T), " tasks (", uniqueN(T$sample), " samples x 2 halves) -> ",
        file.path(OUT_DIR, "tasks.csv"))
message("       cells per half: min ", min(T$n_cells), " median ", median(T$n_cells),
        " max ", max(T$n_cells))
message("       halves under ", CCC_MIN_CELLS_PER_OCCUPIED_BIN * 2L, " cells (likely to lose nodes): ",
        sum(T$n_cells < CCC_MIN_CELLS_PER_OCCUPIED_BIN * 2L))
