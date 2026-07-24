#!/usr/bin/env Rscript
# e20_aggregate_programs.R ----
# Phase 3 step 2 (R): pool all per-sample cNMF programs (every sample, every K), cluster similar
# programs into consensus meta-programs (MPs) Gavish-style, and compute per-MP reproducibility (SRRS).
#
# SIMILARITY (DECISION Q3a): primary = Jaccard over top-CNMF_TOP_N gene SETS; sensitivity = cosine
#   over loadings on the UNION of the two programs' top-N genes (programs store only top-N, so cosine
#   is top-N-restricted -- documented; full-loading cosine would need the spectra matrices, deferred).
# CLUSTERING (DECISION Q3b): Gavish-style -- programs are nodes, edge if Jaccard >= TAU, communities
#   (Leiden if available, else connected components) = consensus MPs.
# SRRS (DECISION Q3c): global = fraction of independent STUDIES (datasets here) in which the MP recurs
#   (>=1 member program). PILOT CAVEAT: only 2 datasets -> SRRS is mechanism-validation only, NOT a
#   real cross-study rate (blueprint requires >=3 studies for conditional, >=0.6 global hard gate).
#   Output flags this with srrs_pilot_insufficient_studies = TRUE.
#
# INPUT  : CNMF_PROG_DIR/<ds>/<sample>__programs.csv  (dataset,sample,K,program,rank,gene,loading)
# OUTPUT : results/tables/05_cnmf/
#            meta_programs.csv          program_uid -> mp_id membership (+ input_mode, dataset)
#            mp_summary.csv             per-MP: n_programs, n_samples, n_datasets, srrs_global, flags
#            mp_top_genes.csv           per-MP consensus top genes (by recurrence across members)
# Usage:
#   Rscript e20_aggregate_programs.R --datasets=GSE227903,GSE239721 [--tau=0.2]

suppressPackageStartupMessages({ library(optparse); library(data.table); library(igraph) })
source("d00_config.R")

CNMF_OUT_DIR <- file.path(FAST_DIR, "results/tables/05_cnmf")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--datasets", type = "character", help = "comma-separated dataset ids"),
  make_option("--tau", type = "double", default = 0.20, help = "Jaccard edge threshold (default 0.20)"),
  make_option("--out", type = "character", default = CNMF_OUT_DIR, help = "output dir")
)))
if (is.null(opt$datasets)) stop("--datasets is required")
datasets <- strsplit(opt$datasets, ",")[[1]]
TAU <- opt$tau

## -- load all programs, attach input_mode from e05 meta ----
read_meta_mode <- function(ds, samp) {
  mj <- file.path(CNMF_INPUT_DIR, ds, samp, "meta.json")
  if (!file.exists(mj)) return(NA_character_)
  txt <- paste(readLines(mj), collapse = " ")
  m <- regmatches(txt, regexpr("\"input_mode\"\\s*:\\s*\"[^\"]+\"", txt))
  if (length(m)) sub(".*:\\s*\"([^\"]+)\".*", "\\1", m) else NA_character_
}

prog_list <- list()
for (ds in datasets) {
  fs <- list.files(file.path(CNMF_PROG_DIR, ds), pattern = "__programs\\.csv$", full.names = TRUE)
  for (f in fs) {
    d <- fread(f)
    if (!nrow(d)) next
    d[, program_uid := paste(dataset, sample, K, program, sep = ":")]
    d[, input_mode := read_meta_mode(dataset[1], sample[1])]
    prog_list[[f]] <- d
  }
}
prog <- rbindlist(prog_list, fill = TRUE)
if (!nrow(prog)) stop("[e20] no programs found (run e10 first)")
message("[1] ", uniqueN(prog$program_uid), " programs from ",
        uniqueN(prog$sample), " samples / ", uniqueN(prog$dataset), " datasets")

## -- top-N gene sets + loading vectors per program ----
setorder(prog, program_uid, rank)
top_sets <- prog[rank <= CNMF_TOP_N, .(genes = list(gene), load = list(setNames(loading, gene))),
                 by = .(program_uid, dataset, sample, K, input_mode)]
uids  <- top_sets$program_uid
gsets <- top_sets$genes
loads <- top_sets$load
np <- length(uids)
message("[2] computing ", np, "x", np, " similarity (Jaccard primary, cosine sensitivity)")

## -- pairwise Jaccard (primary) + cosine on union-of-topN (sensitivity) ----
jac <- matrix(0, np, np, dimnames = list(uids, uids))
cos <- matrix(0, np, np, dimnames = list(uids, uids))
for (i in seq_len(np - 1)) {
  gi <- gsets[[i]]; li <- loads[[i]]
  for (j in (i + 1):np) {
    gj <- gsets[[j]]; lj <- loads[[j]]
    inter <- length(intersect(gi, gj)); uni <- length(union(gi, gj))
    jac[i, j] <- jac[j, i] <- if (uni > 0) inter / uni else 0
    u <- union(gi, gj)
    vi <- ifelse(is.na(li[u]), 0, li[u]); vj <- ifelse(is.na(lj[u]), 0, lj[u])
    den <- sqrt(sum(vi^2)) * sqrt(sum(vj^2))
    cos[i, j] <- cos[j, i] <- if (den > 0) sum(vi * vj) / den else 0
  }
}

## -- Gavish-style graph clustering on Jaccard >= TAU ----
A <- (jac >= TAU) * 1; diag(A) <- 0
g <- graph_from_adjacency_matrix(A, mode = "undirected")
V(g)$name <- uids
comm <- tryCatch(membership(cluster_leiden(g, objective_function = "modularity")),
                 error = function(e) components(g)$membership)
mp_id <- paste0("MP", sprintf("%03d", as.integer(comm)))
names(mp_id) <- uids

memb <- merge(top_sets[, .(program_uid, dataset, sample, K, input_mode)],
              data.table(program_uid = names(mp_id), mp_id = mp_id), by = "program_uid")
fwrite_safe(memb, file.path(opt$out, "meta_programs.csv"))

## -- per-MP summary + SRRS ----
n_studies_total <- uniqueN(prog$dataset)
mp_sum <- memb[, {
  ndat <- uniqueN(dataset)
  .(n_programs = .N, n_samples = uniqueN(sample), n_datasets = ndat,
    srrs_global = ndat / n_studies_total,
    n_malignant_input = sum(input_mode == "malignant"),
    n_fallback_input  = sum(input_mode == "bin_fallback"),
    frac_fallback = mean(input_mode == "bin_fallback"))
}, by = mp_id]
mp_sum[, srrs_pilot_insufficient_studies := (n_studies_total < 3)]
# a MP dominated by bin_fallback inputs is likely a NORMAL-hematopoiesis program, not malignant
mp_sum[, caution := fifelse(frac_fallback > 0.5, "likely_normal_program_fallback_dominated", "")]
setorder(mp_sum, -n_programs)
fwrite_safe(mp_sum, file.path(opt$out, "mp_summary.csv"))

## -- per-MP consensus top genes (by recurrence across member programs) ----
mp_prog <- merge(prog[rank <= CNMF_TOP_N], memb[, .(program_uid, mp_id)], by = "program_uid")
mp_top <- mp_prog[, .(recurrence = uniqueN(program_uid), mean_rank = mean(rank),
                      mean_loading = mean(loading)),
                  by = .(mp_id, gene)]
setorder(mp_top, mp_id, -recurrence, mean_rank)
mp_top <- mp_top[, head(.SD, CNMF_TOP_N), by = mp_id]
fwrite_safe(mp_top, file.path(opt$out, "mp_top_genes.csv"))

## -- console summary ----
message("\n[done] e20 -> ", opt$out)
message("  TAU=", TAU, "  MPs=", uniqueN(memb$mp_id),
        "  (singletons=", sum(table(memb$mp_id) == 1), ")")
message("  multi-dataset MPs (recur in both): ", mp_sum[n_datasets >= 2, .N],
        " / ", nrow(mp_sum))
message("  fallback-dominated (likely normal) MPs: ", mp_sum[frac_fallback > 0.5, .N])
message("  PILOT CAVEAT: ", n_studies_total, " datasets -> SRRS is mechanism-validation only, ",
        "not a real cross-study rate (need >=3 studies).")
print(head(mp_sum[, .(mp_id, n_programs, n_samples, n_datasets, srrs_global, frac_fallback, caution)], 12))
