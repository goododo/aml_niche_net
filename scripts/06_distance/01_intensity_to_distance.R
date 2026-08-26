# 01_intensity_to_distance.R ----
# Phase 5 (stage 06_distance). Intensity -> distance. Collapse each per-sample LR tensor into a
# DIRECTED 7x7 edge-weight matrix, then rank-percentile WITHIN sample to a distance matrix C (D2
# platform-invariance: each sample's edges are ranked only against its own, so depth offsets wash out).
#
# [locked] edge weight = sum(prob) over pval < DIST_PVAL_THRESH LR pairs (= CellChat aggregateNet weight);
#          n_lr_sig (count of significant LR pairs) kept as a co-product (R8 sensitivity arm).
# [locked] 7x7 = 49 DIRECTED edges INCLUDING self-loops (autocrine, e.g. HSC_MPP->HSC_MPP CD99).
#          C_ij = 1 - rank_pct(weight). Absent bin-pairs (no significant LR) -> weight 0 -> low rank -> high C.
# [locked] pval filtered EXPLICITLY here (02 kept ALL edges via subsetCommunication thresh=1).
# [boundary] sender/receiver split for GW symmetry, node mass p_i, unbalanced missing-node masking, and
#          node-feature assembly are ALL deferred to 07 (FGW-prep). This stage emits only C_dir + edge QC.
#
# Residual note: per-sample ranking leaves a mild depth dependence in how zero-edges map to C (a sparser
# sample pulls its zero block to lower C). Controlled downstream by the sparse-graph flag
# (DIST_MIN_EDGES_FLAG) + platform covariates + the count sensitivity arm. Accepted cost of per-sample rank.
#
# INPUT  : CCC_TENSOR_DIR/<ds>/<sample>__ccc_cellchat.csv   (05_ccc/02 output; all edges + pval)
#          DIR_CCC/ccc_sample_manifest.csv                   (ccc_eligible decides which tensors are admitted)
# OUTPUT : DIR_DISTANCE/edge_distance.csv  (long: dataset,sample,timepoint,sender_bin,receiver_bin,
#            weight_probsum,n_lr_sig,rank_pct,C,detected,weight_per_lr,C_perlr ;
#            49 rows/sample, node order = CCC_NODES)
#          DIR_DISTANCE/edge_qc.csv        (per-sample: n_edges_present, total_weight, sparse_flag)
# Usage  : Rscript scripts/06_distance/01_intensity_to_distance.R [--force]
suppressPackageStartupMessages({ library(optparse); library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_distance.R"))   # DIST_* params + DIR_DISTANCE; pulls CCC_NODES/CCC_TENSOR_DIR
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--force", action = "store_true", default = FALSE)
)))

out_dist <- file.path(DIR_DISTANCE, "edge_distance.csv")
out_qc   <- file.path(DIR_DISTANCE, "edge_qc.csv")
# FRESHNESS, not existence. This guard printed "[skip]" and exited 0 on a superseded cohort:
# results/tables/07_fgw/patient_scores.csv holds 148 rows of which 55 name samples that have
# left the cohort, and 47 current samples have never entered CCC at all. Re-running the chain
# hit five of these guards in a row and reported success.
.ins <- c(file.path(DIR_CCC, "ccc_node_features.csv"),
          file.path(DIR_CCC, "ccc_sample_manifest.csv"),
          list.files(CCC_TENSOR_DIR, pattern = "__ccc_cellchat\\.csv$", recursive = TRUE, full.names = TRUE))
if (!is_stale(c(out_dist, out_qc), .ins, force = opt$force)) {
  message("[skip] distance outputs are current"); quit(status = 0)
}
if (file.exists(out_dist)) message("[recompute] ", stale_reason(c(out_dist, out_qc), .ins, force = opt$force))

## -- fixed directed edge grid (49 = 7x7 incl self-loops) ----
grid <- CJ(sender_bin = CCC_NODES, receiver_bin = CCC_NODES)

tfiles <- list.files(CCC_TENSOR_DIR, pattern = "__ccc_cellchat\\.csv$", recursive = TRUE, full.names = TRUE)
message("[1] tensors on disk: ", length(tfiles))
stopifnot(length(tfiles) > 0L)

# The glob is not the cohort. On 2026-08-20 CCC_TENSOR_DIR still held 70 tensors written against a
# superseded 220-row manifest (55 pre-merge sublibrary names + 15 now-ineligible samples); nothing
# downstream would have rejected them. Admission is decided by the manifest, not by what is on disk.
MAN <- fread(file.path(DIR_CCC, "ccc_sample_manifest.csv"))
elig <- MAN[ccc_eligible == TRUE, paste(dataset, sample)]
tkey <- paste(basename(dirname(tfiles)), sub("__ccc_cellchat\\.csv$", "", basename(tfiles)))
drop <- tfiles[!tkey %in% elig]
if (length(drop)) {
  message("[1] NOT in the eligible manifest, dropped: ", length(drop))
  message("      ", paste(head(sub(paste0(CCC_TENSOR_DIR, "/"), "", drop), 10), collapse = "\n      "))
}
tfiles <- tfiles[tkey %in% elig]
miss <- setdiff(elig, tkey)
if (length(miss)) message("[1] eligible but NO tensor: ", length(miss), " -> ",
                          paste(head(miss, 10), collapse = ", "))
message("[1] tensors admitted: ", length(tfiles), " / ", length(elig), " eligible samples")
stopifnot(length(tfiles) > 0L)

process_one <- function(f) {
  t <- fread(f)
  if (nrow(t) == 0L) return(NULL)
  ds <- t$dataset[1]; smp <- t$sample[1]; tp <- as.character(t$timepoint[1])

  # explicit significance filter, then aggregate significant LR pairs -> directed bin-pair edge
  sig <- t[pval < DIST_PVAL_THRESH]
  agg <- sig[, .(weight_probsum = sum(prob), n_lr_sig = .N), by = .(sender_bin, receiver_bin)]

  # expand onto the fixed 49-edge grid; absent pairs -> 0
  e <- agg[grid, on = c("sender_bin", "receiver_bin")]
  e[is.na(weight_probsum), weight_probsum := 0]
  e[is.na(n_lr_sig),       n_lr_sig := 0L]

  # rank-percentile WITHIN sample over all 49 edges: strong -> rank_pct ~ 1 -> C ~ 0 (close)
  N <- nrow(e)                                   # 49
  e[, rank_pct := frank(weight_probsum, ties.method = "average") / N]
  e[, C := 1 - rank_pct]

  ## -- size-independent companion weight, and an explicit detection flag --
  #
  # weight_probsum CONFLATES two things, and measurement (2026-08-26, 6,762 sample-edges) says which
  # part carries the node-size dependency:
  #
  #   cor(weight_probsum, min node size)          +0.688   over ALL edges
  #   cor(weight_probsum, min node size)          +0.119   over DETECTED edges only
  #   cor(n_lr_sig,       min node size)          +0.136   over DETECTED edges only
  #   cor(weight_probsum / n_lr_sig, min size)    +0.015   <- essentially independent
  #
  # So CellChat's per-pair probability is fine. The dependency is almost entirely PRESENCE: a node
  # with ~10 cells produces an unstable triMean, nboot=100 cannot reach pval < 0.05 on any LR pair,
  # and the edge is recorded as weight 0. That zero then ranks low and becomes C ~ 0.68 -- a numeric
  # DISTANCE standing in for "we could not look". The two are not the same statement.
  #
  # weight_per_lr is the mean strength of the channels that WERE detected, which the measurement
  # above shows is not a function of how many cells were available to detect them. n_lr_sig keeps
  # the channel COUNT, which is the size-confounded part, so the two stay separable downstream.
  # `detected` exists so a consumer can mask rather than impute -- see C_perlr below.
  e[, detected := n_lr_sig > 0L]
  e[, weight_per_lr := fifelse(detected, weight_probsum / pmax(n_lr_sig, 1L), NA_real_)]

  # C_perlr ranks ONLY the detected edges among themselves; undetected edges stay NA instead of
  # being handed the rank an absent measurement would earn. A consumer that cannot take NA must say
  # so and choose an imputation explicitly, which is the point.
  nd <- sum(e$detected)
  e[, C_perlr := NA_real_]
  if (nd > 0L)
    e[detected == TRUE,
      C_perlr := 1 - frank(weight_per_lr, ties.method = "average") / nd]

  e[, `:=`(dataset = ds, sample = smp, timepoint = tp)]
  e[, sender_bin   := factor(sender_bin,   levels = CCC_NODES)]  # enforce node order for 07 matrix assembly
  e[, receiver_bin := factor(receiver_bin, levels = CCC_NODES)]
  setorder(e, sender_bin, receiver_bin)
  e[]
}

res <- rbindlist(lapply(tfiles, function(f)
  tryCatch(process_one(f),
           error = function(err) { message("  [ERROR] ", basename(f), " : ", conditionMessage(err)); NULL })
), fill = TRUE)

## -- per-sample QC: present (non-zero) edges + sparse flag ----
qc <- res[, .(n_edges_present = sum(weight_probsum > 0),
              total_weight    = sum(weight_probsum)),
          by = .(dataset, sample, timepoint)]
qc[, sparse_flag := n_edges_present < DIST_MIN_EDGES_FLAG]

setcolorder(res, c("dataset","sample","timepoint","sender_bin","receiver_bin",
                   "weight_probsum","n_lr_sig","rank_pct","C",
                   "detected","weight_per_lr","C_perlr"))
fwrite_safe(res, out_dist)
fwrite_safe(qc,  out_qc)
message("[2] wrote ", out_dist, "  (", nrow(res), " rows = ",
        uniqueN(res[, .(dataset, sample)]), " samples x 49 edges)")
message("[2] wrote ", out_qc)

## -- report: edge-count distribution + flagged sparse samples ----
message("[3] present-edge distribution (of 49) across samples:")
print(qc[, .(min = min(n_edges_present),
             q25 = quantile(n_edges_present, .25),
             median = as.numeric(median(n_edges_present)),
             q75 = quantile(n_edges_present, .75),
             max = max(n_edges_present))])
message("[3] sparse-flagged (< ", DIST_MIN_EDGES_FLAG, " present edges): ", qc[sparse_flag == TRUE, .N], " samples")
print(qc[sparse_flag == TRUE][order(n_edges_present), .(dataset, sample, timepoint, n_edges_present)])
message("[done]")
