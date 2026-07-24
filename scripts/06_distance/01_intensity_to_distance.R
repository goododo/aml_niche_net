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
# OUTPUT : DIR_DISTANCE/edge_distance.csv  (long: dataset,sample,timepoint,sender_bin,receiver_bin,
#            weight_probsum,n_lr_sig,rank_pct,C ; 49 rows/sample, node order = CCC_NODES)
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
if (file.exists(out_dist) && file.exists(out_qc) && !opt$force) {
  message("[skip] outputs exist; --force to rebuild"); quit(status = 0)
}

## -- fixed directed edge grid (49 = 7x7 incl self-loops) ----
grid <- CJ(sender_bin = CCC_NODES, receiver_bin = CCC_NODES)

tfiles <- list.files(CCC_TENSOR_DIR, pattern = "__ccc_cellchat\\.csv$", recursive = TRUE, full.names = TRUE)
message("[1] tensors found: ", length(tfiles))
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
                   "weight_probsum","n_lr_sig","rank_pct","C"))
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
