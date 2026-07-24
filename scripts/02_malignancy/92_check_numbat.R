#!/usr/bin/env Rscript
# 92_check_numbat.R ----
# Completeness + integrity audit of the Numbat arm, built for an interruptible session. Numbat is
# ~20h/sample, so the key value is knowing WHERE a killed sample stopped, to resume from the
# cheapest safe point instead of re-downloading + re-aligning:
#   DONE_CNV          : percell CSV present + readable + has clone/compartment (real CNV call)
#   DONE_DEGRADED     : percell CSV present, note == NUMBAT_DEGRADED_NOTE (valid "no CNV" -> abstain)
#   RESUMABLE_NUMBAT  : allele_counts + star_filtered present but no/bad percell -> RE-RUN 30 ONLY
#                       (pileup/phasing already done; no re-download/STARsolo needed)
#   PARTIAL_EARLY     : sample dir exists but no allele_counts -> phasing didn't finish -> full chain
#   CORRUPT           : percell exists but unreadable/empty
#   MISSING           : (only with --sheet) listed but no sample dir at all
#
# Scans CNV_ROOT/numbat/<dataset>/<sample>/ . Pass --sheet dataset<TAB>sample to also flag MISSING.
#
# Run (from PROJECT ROOT):
#   Rscript scripts/02_malignancy/92_check_numbat.R [--sheet samples_all.tsv]

suppressPackageStartupMessages({ library(optparse); library(data.table) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sheet", type = "character", default = "",
              help = "optional dataset<TAB>sample list to detect never-started (MISSING) samples")
)))

NB_ROOT <- file.path(CNV_ROOT, "numbat")
if (!dir.exists(NB_ROOT)) stop("numbat root not found: ", NB_ROOT)

# discover <dataset>/<sample> dirs on disk
ds_dirs <- list.dirs(NB_ROOT, recursive = FALSE)
disk <- rbindlist(lapply(ds_dirs, function(dd) {
  ds <- basename(dd)
  ss <- list.dirs(dd, recursive = FALSE)
  if (!length(ss)) return(NULL)
  data.table(dataset = ds, sample_id = basename(ss), dir = ss)
}), fill = TRUE)
if (!nrow(disk)) stop("no <dataset>/<sample> dirs under ", NB_ROOT)

classify <- function(dir) {
  pc     <- file.path(dir, "numbat", paste0(basename(dir), "__numbat_percell.csv"))
  allele <- length(list.files(dir, pattern = "allele_counts\\.tsv\\.gz$", full.names = TRUE)) > 0
  star   <- dir.exists(file.path(dir, "star_filtered"))
  if (file.exists(pc)) {
    d <- tryCatch(fread(pc), error = function(e) NULL)
    if (is.null(d) || !nrow(d)) return(list(status = "CORRUPT", note = "unreadable/empty percell"))
    if ("note" %in% names(d) && all(d$note == NUMBAT_DEGRADED_NOTE, na.rm = TRUE))
      return(list(status = "DONE_DEGRADED", note = "no_CNV_detected (abstains in 50)"))
    if (any(c("compartment_opt", "clone_opt", "malignant") %in% names(d)))
      return(list(status = "DONE_CNV", note = sprintf("%d cells", nrow(d))))
    return(list(status = "CORRUPT", note = "percell has no clone/compartment/note column"))
  }
  if (allele && star) return(list(status = "RESUMABLE_NUMBAT", note = "allele+matrix present; re-run 30 only"))
  if (allele)         return(list(status = "RESUMABLE_NUMBAT", note = "allele present; matrix?"))
  list(status = "PARTIAL_EARLY", note = "no allele_counts; phasing incomplete")
}

res <- rbindlist(lapply(seq_len(nrow(disk)), function(i) {
  c1 <- classify(disk$dir[i])
  data.table(dataset = disk$dataset[i], sample_id = disk$sample_id[i],
             status = c1$status, note = c1$note)
}))

# optional MISSING detection against a sheet
if (nzchar(opt$sheet) && file.exists(opt$sheet)) {
  sh <- fread(opt$sheet, header = FALSE, col.names = c("dataset", "sample_id"))
  sh[, key := paste(dataset, sample_id)]; res[, key := paste(dataset, sample_id)]
  miss <- sh[!key %in% res$key]
  if (nrow(miss)) res <- rbind(res, miss[, .(dataset, sample_id, status = "MISSING",
                                             note = "in sheet, no dir on disk")], fill = TRUE)
  res[, key := NULL]
}

cat("\n==================== Numbat STATUS ====================\n")
print(res[, .N, by = status][order(-N)])

done <- res[status %in% c("DONE_CNV", "DONE_DEGRADED")]
cat(sprintf("\nDONE (usable by 50): %d  [CNV=%d, degraded=%d]\n",
            nrow(done), res[status == "DONE_CNV", .N], res[status == "DONE_DEGRADED", .N]))

bad <- res[!status %in% c("DONE_CNV", "DONE_DEGRADED")]
if (nrow(bad)) {
  cat("\n-- samples needing attention (grouped by resume point) --\n")
  print(bad[order(status, dataset, sample_id)])
  cat("\n[resume guidance]\n")
  if (res[status == "RESUMABLE_NUMBAT", .N])
    cat("  RESUMABLE_NUMBAT -> re-run 30 only (allele/matrix kept). Per sample:\n",
        "    Rscript scripts/02_malignancy/30_numbat_run.R --sample <S> \\\n",
        "      --matrix ", file.path(CNV_ROOT, "numbat", "<ds>", "<S>", "star_filtered"), " \\\n",
        "      --allele ", file.path(CNV_ROOT, "numbat", "<ds>", "<S>", "<S>_allele_counts.tsv.gz"), " \\\n",
        "      --outdir ", file.path(CNV_ROOT, "numbat", "<ds>", "<S>", "numbat"),
        " --genome hg38 --ncores 16 --ref_rds ",
        file.path(REF_DIR, "numbat", "numbat_ref_hg38.rds"), "\n", sep = "")
  if (res[status %in% c("PARTIAL_EARLY", "CORRUPT", "MISSING"), .N])
    cat("  PARTIAL_EARLY / CORRUPT / MISSING -> re-run the full chain (10/11) for those samples\n",
        "    (clear the stale dir first: rm -rf ", file.path(CNV_ROOT, "numbat", "<ds>", "<S>"), ").\n", sep = "")
}

cat("\n==================== VERDICT ====================\n")
if (!nrow(bad)) {
  cat(sprintf("ALL %d Numbat samples on disk are DONE and readable.\n", nrow(res)))
} else {
  cat(sprintf("%d DONE, %d need attention (see resume guidance above).\n", nrow(done), nrow(bad)))
}
