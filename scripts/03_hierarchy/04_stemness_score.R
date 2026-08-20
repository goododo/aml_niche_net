#!/usr/bin/env Rscript
# 04_stemness_score.R ----
# Phase 2.4 -- orthogonal (transcriptional) test of the malignant-distribution V-shape from 03, and
# a check on the MRD artifact hypothesis. Scores every cell for stemness independently of the BMM
# bin projection, then asks:
#   (i)  do MALIGNANT cells' stemness scores track the V-shape (high Dx, low MRD, rebound Relapse)?
#   (ii) at MRD, are there HIGH-stemness cells in the HSC_MPP bin that were NOT called malignant?
#        -> if yes, the MRD "differentiation" is partly a detection artifact (quiescent LSCs that
#        inferCNV can't call because of low transcriptional output), not a true loss of LSCs.
#
# Signatures (from stemness_signatures.tsv):
#   LSC17     : published weighted LSC score (Ng 2016; weighted sum of log-norm expression). The
#               weighted score is an ADAPTATION to scRNA -- noisy per cell, but valid in aggregate.
#   HSPC_core : curated canonical HSPC markers, scored with Seurat AddModuleScore (dropout-robust).
#   (add van Galen 2019 Table S3 sets as more rows -> auto-picked up here.)
#
# INPUT : QC objects (QC_RDS_DIR) ; consensus per-cell (DIR_MALIGNANCY) ; projection per-cell (HIER_PROJ_DIR)
# OUTPUT: HIER_PROJ_DIR/<ds>/<sample>__stemness_percell.csv
#           cell, sample, dataset, timepoint, malignant, hierarchy_bin, LSC17, HSPC_core
#         HIER_TAB_DIR/stemness_by_timepoint.csv (summaries)

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(data.table); library(optparse) })
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", default = "", help = "comma-separated dataset filter (run GSE227903 first)"),
  make_option("--limit",    type = "integer",   default = 0L, help = "first N samples (0=all)"),
  # THERE WAS NO --force. Combined with the existence-only guard below, the 68 stemness files
  # dated 2026-07-17 could not be refreshed by any invocation of this script.
  make_option("--force",    action = "store_true", default = FALSE, help = "rebuild regardless of freshness")
)))

sig <- load_stemness()                       # signature, gene, weight
sigs <- split(sig, sig$signature)
message("[sig] ", paste(sprintf("%s(%d genes)", names(sigs), sapply(sigs, nrow)), collapse = ", "))

# CURATED timepoint (see 02_per_bin_malignant.R and TP_AXIS_LEVELS). The regex this replaces was
# not even identical to 02's copy -- 02 used "(_|^)(REL|R2?)$", this one "(_|^)REL|(_|^)R[0-9]*$" --
# so the two stages could disagree about the same sample. They happen not to on this cohort.
TPMAP <- sample_timepoints()

.data_mat <- function(seu) tryCatch(SeuratObject::LayerData(seu, assay = "RNA", layer = "data"),
                                    error = function(e) Seurat::GetAssayData(seu, assay = "RNA", slot = "data"))

## ---- samples with QC rds + consensus + projection ----
# Roster from the QC report, not from ls -- see qc_rds_roster() in utils.R.
qc <- qc_rds_roster(on_extra = "error")[, .(dataset, sample, rds)]
cons <- data.table(cons = list.files(DIR_MALIGNANCY, pattern = "__consensus_percell\\.csv$", recursive = TRUE, full.names = TRUE))
cons[, sample := sub("__consensus_percell\\.csv$", "", basename(cons))][, dataset := basename(dirname(cons))]
proj <- data.table(proj = list.files(HIER_PROJ_DIR, pattern = "__bmm_percell\\.csv$", recursive = TRUE, full.names = TRUE))
proj[, sample := sub("__bmm_percell\\.csv$", "", basename(proj))][, dataset := basename(dirname(proj))]
man <- Reduce(function(a, b) merge(a, b, by = c("dataset", "sample")), list(qc, cons, proj))
if (nzchar(opt$datasets)) man <- man[dataset %in% trimws(strsplit(opt$datasets, ",")[[1]])]
if (opt$limit > 0) man <- man[seq_len(min(opt$limit, .N))]
stopifnot(nrow(man) > 0)
message(sprintf("[0] %d sample(s) to score", nrow(man)))

score_one <- function(seu, sid, ds, first) {
  seu <- NormalizeData(seu, verbose = FALSE)
  dat <- .data_mat(seu)
  out <- data.table(cell = colnames(seu))
  for (sg in names(sigs)) {
    g <- sigs[[sg]]$gene; w <- sigs[[sg]]$weight
    present <- g %in% rownames(dat)
    if (first) message(sprintf("    [sig] %s: %d/%d genes present", sg, sum(present), length(g)))
    if (sum(present) < 3) { out[, (sg) := NA_real_]; next }
    if (any(abs(w[present] - 1) > 1e-9)) {                    # weighted (LSC17): weighted sum of log-norm
      wv <- w[present]; out[, (sg) := as.numeric(wv %*% as.matrix(dat[g[present], , drop = FALSE]))]
    } else {                                                  # unweighted: AddModuleScore (dropout-robust)
      seu <- AddModuleScore(seu, features = list(g[present]), name = paste0("MS_", sg), seed = SEED, ctrl = 50)
      out[, (sg) := seu[[paste0("MS_", sg, "1")]][, 1]]
    }
  }
  out
}

for (i in seq_len(nrow(man))) {
  ds <- man$dataset[i]; sid <- man$sample[i]
  dst <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__stemness_percell.csv"))
  # FRESHNESS, not existence. All three inputs postdate the 2026-07-17 outputs (QC objects
  # 08-05, projection 08-07, consensus 08-13/14/17) and 48,087 cells in them no longer exist.
  if (!is_stale(dst, c(man$rds[i], man$cons[i], man$proj[i]), force = opt$force)) {
    message(sprintf("[%d/%d] %s::%s current", i, nrow(man), ds, sid)); next
  }
  if (file.exists(dst))
    message(sprintf("[%d/%d] %s::%s RECOMPUTE -- %s", i, nrow(man), ds, sid,
                    stale_reason(dst, c(man$rds[i], man$cons[i], man$proj[i]), force = opt$force)))
  message(sprintf("[%d/%d] %s::%s", i, nrow(man), ds, sid))
  seu <- readRDS(man$rds[i])
  sc  <- score_one(seu, sid, ds, first = (i == 1))
  mal <- fread(man$cons[i], select = c("cell", "malignant"))
  bin <- fread(man$proj[i], select = c("cell", "hierarchy_bin"))
  sc  <- merge(merge(sc, mal, by = "cell", all.x = TRUE), bin, by = "cell", all.x = TRUE)
  .tp <- TPMAP[dataset == ds & sample == sid]
  sc[, `:=`(sample = sid, dataset = ds,
            timepoint = if (nrow(.tp)) .tp$timepoint[1] else NA_character_,
            tp_axis   = if (nrow(.tp)) .tp$tp_axis[1]   else NA_character_,
            patient   = if (nrow(.tp)) .tp$uid_patient[1] else NA_character_)]
  fwrite_safe(sc, dst)
  rm(seu, sc); gc(verbose = FALSE)
}

## ---- aggregate + the two key comparisons ----
allsc <- rbindlist(lapply(file.path(man$dataset, man$sample), function(x)
  fread(file.path(HIER_PROJ_DIR, dirname(x), paste0(basename(x), "__stemness_percell.csv")))), fill = TRUE)
# THE AXIS, NOT A LITERAL. This line still read c("Dx","MRD","Relapse") after the per-cell
# timepoint moved to the curated vocabulary, so only "Relapse" matched and
# stemness_by_timepoint.csv silently went from three rows to one. Caught by comparing the
# rebuilt table against its own snapshot, not by any error.
lon <- allsc[tp_axis %in% TP_AXIS_LEVELS]
lon[, tp := factor(tp_axis, levels = TP_AXIS_LEVELS)]

if (nrow(lon)) {
  sig_cols <- intersect(names(sigs), names(lon))   # all scored signatures present
  cat("\n================ (i) stemness of MALIGNANT cells by timepoint (does it track the V-shape?) ================\n")
  s1 <- lon[malignant == 1, c(list(n_malignant = .N),
             lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3))), by = tp, .SDcols = sig_cols][order(tp)]
  print(s1)
  cat("\n================ (ii) HSC_MPP-bin cells: stemness by timepoint x malignant call ================\n")
  cat("   (if MRD has HIGH-stemness cells that are malignant==0, the MRD 'loss of LSC' is partly a\n")
  cat("    detection artifact -- quiescent LSCs inferCNV can't call, not true differentiation)\n")
  s2 <- lon[hierarchy_bin == "HSC_MPP", c(list(n = .N),
             lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3))), by = .(tp, malignant),
             .SDcols = sig_cols][order(tp, malignant)]
  print(s2)
  fwrite_safe(s1[, grp := "malignant_by_tp"], file.path(HIER_TAB_DIR, "stemness_by_timepoint.csv"))
}

cat("\n[note] LSC17 here is a scRNA adaptation of the bulk weighted score (per-cell noisy, valid in\n")
cat("       aggregate). HSPC_core is a curated canonical set (AddModuleScore). Add van Galen 2019\n")
cat("       Table S3 signatures as rows in stemness_signatures.tsv to cross-validate.\n")
