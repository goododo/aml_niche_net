#!/usr/bin/env Rscript
# 91_check_infercnv.R ----
# Completeness + integrity audit of the inferCNV arm, built for an interruptible session: tells you
# exactly which runnable samples FINISHED, which are CORRUPT (burden CSV half-written by a killed
# job -> resume-skip would wrongly skip them), which are PARTIAL (inferCNV dir started but no burden
# -> must be cleared before re-running), and which are MISSING (not started). Emits a ready-to-run
# resume task list + the cleanup commands for the bad ones.
#
# Universe = runnable samples from ref_norm_summary (route != skip / method != sorted_guard).
# A sample is DONE iff its burden CSV exists, reads, has the right columns, and >0 rows.
#
# Run (from PROJECT ROOT):
#   Rscript scripts/02_malignancy/91_check_infercnv.R [--min_rows 50] [--out infercnv_tasks_remaining.tsv]

suppressPackageStartupMessages({ library(optparse); library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--min_rows", type = "integer", default = 50L,
              help = "a burden CSV with fewer rows than this is treated as CORRUPT [50]"),
  make_option("--out", type = "character", default = "infercnv_tasks_remaining.tsv",
              help = "resume task list (dataset<TAB>sample) for the not-DONE samples")
)))

stopifnot(file.exists(REFNORM_SUMMARY_CSV))
summ <- fread(REFNORM_SUMMARY_CSV)
summ[, route := fifelse(decision == "autologous_ok", "autologous",
              fifelse(method == "sorted_guard", "skip", "external"))]
run <- summ[route != "skip", .(dataset, sample_id, route)]
message(sprintf("[universe] %d runnable samples (skip=%d not audited)",
                nrow(run), summ[route == "skip", .N]))

classify <- function(ds, s) {
  burden <- file.path(INFERCNV_BURDEN_ROOT, ds, paste0(s, "_infercnv_burden.csv"))
  odir   <- file.path(INFERCNV_ROOT, ds, s)
  if (file.exists(file.path(odir, "lowobs.skip")))
    return(list(status = "SKIPPED_LOWOBS", nrow = 0L, note = "too few observation cells (inferCNV N/A)"))
  if (!file.exists(burden)) {
    return(list(status = if (dir.exists(odir)) "PARTIAL" else "MISSING", nrow = 0L, note = ""))
  }
  d <- tryCatch(fread(burden), error = function(e) NULL)
  if (is.null(d))                                   return(list(status="CORRUPT", nrow=0L, note="unreadable"))
  if (!all(c("cell","infercnv_burden","group") %in% names(d)))
                                                    return(list(status="CORRUPT", nrow=nrow(d), note="bad columns"))
  if (nrow(d) < opt$min_rows)                       return(list(status="CORRUPT", nrow=nrow(d), note="too few rows"))
  if (anyNA(d$infercnv_burden) && all(is.na(d$infercnv_burden)))
                                                    return(list(status="CORRUPT", nrow=nrow(d), note="all-NA burden"))
  list(status = "DONE", nrow = nrow(d), note = "")
}

res <- rbindlist(lapply(seq_len(nrow(run)), function(i) {
  c1 <- classify(run$dataset[i], run$sample_id[i])
  data.table(dataset = run$dataset[i], sample_id = run$sample_id[i], route = run$route[i],
             status = c1$status, n_rows = c1$nrow, note = c1$note)
}))

cat("\n==================== inferCNV STATUS ====================\n")
print(res[, .N, by = status][order(-N)])
DONE_SET <- c("DONE", "SKIPPED_LOWOBS")
bad <- res[!status %in% DONE_SET]
if (nrow(bad)) { cat("\n-- samples needing attention --\n"); print(bad[order(status, dataset, sample_id)]) }

## resume task list (everything not DONE / not intentionally skipped)
todo <- res[!status %in% DONE_SET, .(dataset, sample_id)]
if (nrow(todo)) {
  fwrite(todo, opt$out, sep = "\t", col.names = FALSE)
  cat(sprintf("\n[resume] wrote %d remaining task(s) -> %s\n", nrow(todo), opt$out))
}

## cleanup commands for CORRUPT/PARTIAL (must be removed before re-run, or inferCNV resumes dirty)
dirty <- res[status %in% c("CORRUPT", "PARTIAL")]
if (nrow(dirty)) {
  cat("\n[cleanup] run these BEFORE resubmitting (stale/partial outputs block a clean re-run):\n")
  for (i in seq_len(nrow(dirty))) {
    ds <- dirty$dataset[i]; s <- dirty$sample_id[i]
    cat(sprintf("  rm -f %s/%s/%s_infercnv_burden.csv; rm -rf %s/%s/%s\n",
                INFERCNV_BURDEN_ROOT, ds, s, INFERCNV_ROOT, ds, s))
  }
}

cat("\n==================== VERDICT ====================\n")
n_done <- res[status %in% DONE_SET, .N]
n_skip_lowobs <- res[status == "SKIPPED_LOWOBS", .N]
if (n_done == nrow(res)) {
  cat(sprintf("ALL %d runnable samples resolved (%d DONE, %d skipped-low-obs). inferCNV complete -> run 41.\n",
              nrow(res), res[status == "DONE", .N], n_skip_lowobs))
} else {
  cat(sprintf("%d/%d resolved (incl. %d skipped-low-obs). To resume: run the cleanup lines above, then:\n",
              n_done, nrow(res), n_skip_lowobs))
  cat(sprintf("  N=$(wc -l < %s); sbatch --array=1-${N}%%6 scripts/02_malignancy/45_infercnv_submit.slurm %s\n",
              opt$out, opt$out))
}
