#!/usr/bin/env Rscript
# =============================================================================
# 14_quarantine_stale_outputs.R
#
# Move per-sample artefacts that are NOT in the current PASS roster out of the
# working tree.
#
# WHY THIS EXISTS. Most stages build their sample list by globbing a directory,
# which makes "a file is on disk" the operative definition of "in the cohort".
# The v2 re-ingest dropped samples the v1 ingest had written (Chen2023 sorted
# fractions, the v1 GSE185381 roster, Petti2019 ND_/Normal_sorted_, E-MTAB -BLD)
# but their outputs stayed where they were, so the old definition kept returning
# them. J3 processed 267 samples for a 214-sample cohort and reported success.
# The same 62 phantoms sit in the malignancy outputs one stage further down,
# where they would enter ALL_consensus_summary -- the input to the calibration
# gate itself.
#
# NOTHING IS DELETED. Artefacts move to a sibling _stale_v1/ directory with a
# manifest, so v1 stays available for comparison and the move is one mv away
# from being undone. Selection is made against the ROSTER, never against mtime:
# a fresh timestamp on a sample that left the cohort is still stale.
#
# Usage: Rscript scripts/99_admin/14_quarantine_stale_outputs.R [--apply]
#        (default is a dry run)
# =============================================================================

suppressPackageStartupMessages({library(data.table); library(optparse); library(here)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--apply", action = "store_true", default = FALSE,
              help = "actually move the files (default: dry run)")
)))

source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

# root, filename pattern, and how to recover the sample id from the basename.
# Everything here is laid out as <root>/<dataset>/<sample><suffix>.
TARGETS <- list(
  list(id = "qc_object",        root = QC_RDS_DIR,            pat = "\\.rds$"),
  list(id = "infercnv_burden",  root = INFERCNV_BURDEN_ROOT,  pat = "_infercnv_burden\\.csv$"),
  list(id = "consensus_summary",root = DIR_MALIGNANCY,        pat = "__consensus_summary\\.csv$"),
  list(id = "consensus_percell",root = DIR_MALIGNANCY,        pat = "__consensus_percell\\.csv$"),
  list(id = "bmm_projection",   root = file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected"),
                                                              pat = "__bmm_percell\\.csv$")
)

R <- qc_rds_roster(on_extra = "ignore")          # the roster IS the QC report
R[, k := paste(dataset, sample)]
Q <- fread(file.path(DIR_PREPROCESS, "03_qc_report__ALL.csv")); Q[, k := paste(dataset, Sample)]
message("[1] roster (PASS) = ", nrow(R), " samples")

collect <- function(t) {
  if (!dir.exists(t$root)) return(data.table())
  f <- list.files(t$root, pattern = t$pat, recursive = TRUE, full.names = TRUE)
  f <- f[!grepl("/_", f)]                        # scratch / already-quarantined
  if (!length(f)) return(data.table())
  x <- data.table(kind = t$id, root = t$root, from = f,
                  dataset = basename(dirname(f)),
                  sample  = sub(t$pat, "", basename(f)),
                  mtime   = file.mtime(f), mb = round(file.size(f) / 1e6, 1))
  x[, k := paste(dataset, sample)]
  x[!k %in% R$k]                                 # keep only what is NOT in the roster
}
S <- rbindlist(lapply(TARGETS, collect), fill = TRUE)

if (!nrow(S)) { message("[done] every artefact matches the roster; nothing to do."); quit(save = "no") }

# "absent from the report" and "present but excluded by QC" are different facts.
S[, why := Q$status[match(k, Q$k)]]
S[is.na(why), why := "absent from the QC report (dropped by the v2 ingest)"]

message("\n[2] artefacts not in the roster:")
print(S[, .(n = .N, samples = uniqueN(k), oldest = min(mtime), newest = max(mtime),
            gb = round(sum(mb) / 1e3, 2)), by = .(kind, why)][order(kind, -n)])
message("\n    by dataset:")
print(dcast(S[, .N, by = .(dataset, kind)], dataset ~ kind, value.var = "N", fill = 0))

if (!opt$apply) {
  message("\n[dry run] re-run with --apply to move ", nrow(S), " file(s) across ",
          uniqueN(S$kind), " artefact type(s)")
  quit(save = "no")
}

S[, dest := file.path(dirname(root), paste0("_stale_v1__", basename(root)), dataset, basename(from))]
for (d in unique(dirname(S$dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
S[, moved := file.rename(from, dest)]
if (!all(S$moved)) stop("failed to move ", sum(!S$moved), " file(s); nothing else changed")

man_dir <- file.path(PROJECT_ROOT, "00_project", "metadata", "quarantine")
dir.create(man_dir, recursive = TRUE, showWarnings = FALSE)
fwrite_safe(S[, .(kind, dataset, sample, why, mtime, mb, from, to = dest)],
            file.path(man_dir, "stale_v1_manifest.csv"))
message("\n[3] moved ", nrow(S), " file(s); manifest -> ", file.path(man_dir, "stale_v1_manifest.csv"))

# The point is that the roster assertion stops complaining, not that files moved.
invisible(qc_rds_roster(on_extra = "error"))
message("[done] every artefact root now matches the PASS roster (", nrow(R), " samples).")
