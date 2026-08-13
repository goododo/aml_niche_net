#!/usr/bin/env Rscript
# 98_perlineage_threshold_eval.R ----
# VERDICT: per-lineage thresholding was tested here and REJECTED. Kept in the tree so the
# rejection is reproducible rather than remembered -- the pooled per-bin FPR from 96 makes a
# per-lineage threshold look like the obvious fix, and it is not.
#
#   healthy FPR, all cells : 0.3305 (global)  ->  0.3590 (per-lineage)   WORSE
#   healthy FPR, CCC bins  : 0.2334 (global)  ->  0.2438 (per-lineage)   WORSE
#   worst regressions: Erythroid +0.284, Megakaryocyte +0.087, HSC_MPP +0.041
#
# Why it loses: a lineage's own reference block is small (Erythroid 450 cells, Megakaryocyte
# 150), so its P95 is both noisy and biased low -- a short tail underestimates the quantile --
# and more observation cells clear it. The pooled reference (~2400 cells) gives a higher, more
# stable threshold. Do not re-adopt this without first enlarging the per-lineage reference.
#
# Read-only: writes nothing into the pipeline.
#
#   Rscript scripts/02_malignancy/98_perlineage_threshold_eval.R
#
# ---------------------------------------------------------------------------------------------
# Evaluate: does a PER-LINEAGE threshold beat the single global threshold?
#
# Global (current): thr = P95 over ALL reference cells pooled; applied to every observation cell.
# Per-lineage     : thr_b = P95 over the reference cells OF THAT LINEAGE; an observation cell
#                   projected to bin b is judged against thr_b.
#
# Truth set = healthy donors (every flagged cell is a false positive), so lower FPR = better.
# Read-only: writes nothing into the pipeline.

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
HIER_PROJ_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")

R <- qc_rds_roster(on_extra = "error")
info <- rbindlist(lapply(seq_len(nrow(R)), function(i) {
  m <- tryCatch(readRDS(R$rds[i])@meta.data, error = function(e) NULL); if (is.null(m)) return(NULL)
  data.table(dataset = R$dataset[i], sample = R$sample[i],
             timepoint = if ("Timepoint" %in% names(m)) as.character(m$Timepoint[1]) else NA_character_)
}), fill = TRUE)
info[, healthy := timepoint == "Healthy" | mapply(function(s) isTRUE(is_healthy_sample(s)), sample)]
H <- info[healthy == TRUE]
message(sprintf("[0] %d healthy donors", nrow(H)))

out <- rbindlist(lapply(seq_len(nrow(H)), function(i) {
  ds <- H$dataset[i]; sid <- H$sample[i]
  bf <- file.path(INFERCNV_BURDEN_ROOT, ds, paste0(sid, "_infercnv_burden.csv"))
  pf <- file.path(HIER_PROJ_DIR, ds, paste0(sid, "__bmm_percell.csv"))
  if (!file.exists(bf) || !file.exists(pf)) return(NULL)
  d <- fread(bf)
  proj <- fread(pf, select = c("cell", "hierarchy_bin", "in_ccc_graph"))

  is_ref <- grepl("reference", d$group, ignore.case = TRUE)
  # lineage of a REFERENCE cell comes from its group label; autologous normals are T/B by
  # construction (REF_KEEP_LABELS), so they back the lymphoid bins.
  ref_bin <- ifelse(grepl("^reference_external__", d$group),
                    sub("^reference_external__", "", d$group),
                    ifelse(d$group == "reference_normal", "LYMPHOID_AUTO", NA_character_))

  thr_global <- as.numeric(quantile(d$infercnv_burden[is_ref], INFERCNV_SCORE_Q, na.rm = TRUE))

  obs <- d[group == "observation", .(cell, infercnv_burden)]
  obs <- merge(obs, proj, by = "cell", all.x = TRUE)
  obs <- obs[!is.na(hierarchy_bin) & nzchar(hierarchy_bin)]
  if (!nrow(obs)) return(NULL)

  # per-lineage threshold table; lymphoid bins draw on the autologous normals too
  thr_of <- function(b) {
    src <- if (b %in% c("T_NK", "B_Plasma")) c(b, "LYMPHOID_AUTO") else b
    v <- d$infercnv_burden[is_ref & ref_bin %in% src]
    if (length(v) < INFERCNV_REF_MIN_PER_GROUP) return(NA_real_)
    as.numeric(quantile(v, INFERCNV_SCORE_Q, na.rm = TRUE))
  }
  tt <- data.table(hierarchy_bin = unique(obs$hierarchy_bin))
  tt[, thr_lin := sapply(hierarchy_bin, thr_of)]
  obs <- merge(obs, tt, by = "hierarchy_bin", all.x = TRUE)
  # a bin with too few reference cells falls back to the global threshold rather than abstaining
  obs[, thr_used := fifelse(is.na(thr_lin), thr_global, thr_lin)]

  obs[, `:=`(fp_global = as.integer(infercnv_burden > thr_global),
             fp_lin    = as.integer(infercnv_burden > thr_used))]
  obs[, .(dataset = ds, sample = sid, hierarchy_bin, in_ccc_graph,
          fp_global, fp_lin, fell_back = is.na(thr_lin))]
}), fill = TRUE)

message(sprintf("[1] %d healthy observation cells with a projected bin", nrow(out)))

cat("\n=============== healthy FPR: global vs per-lineage threshold ===============\n")
cat(sprintf("  cells                    : %d\n", nrow(out)))
cat(sprintf("  FPR global threshold     : %.4f\n", mean(out$fp_global)))
cat(sprintf("  FPR per-lineage threshold: %.4f\n", mean(out$fp_lin)))
cat(sprintf("  bins falling back to global (too few ref cells): %.1f%% of cells\n",
            100 * mean(out$fell_back)))

cat("\n-- CCC bins only --\n")
cc <- out[in_ccc_graph == TRUE]
cat(sprintf("  cells %d | global %.4f -> per-lineage %.4f\n",
            nrow(cc), mean(cc$fp_global), mean(cc$fp_lin)))

cat("\n-- by bin (CCC only) --\n")
print(cc[, .(n = .N, FPR_global = round(mean(fp_global), 4),
             FPR_perlineage = round(mean(fp_lin), 4),
             delta = round(mean(fp_lin) - mean(fp_global), 4)), by = hierarchy_bin][order(-FPR_global)])

cat("\n-- by dataset (CCC only) --\n")
print(cc[, .(n = .N, FPR_global = round(mean(fp_global), 4),
             FPR_perlineage = round(mean(fp_lin), 4),
             delta = round(mean(fp_lin) - mean(fp_global), 4)), by = dataset][order(-FPR_global)])

cat("\n-- bin x dataset delta (negative = per-lineage helps) --\n")
print(dcast(cc[, .(d = round(mean(fp_lin) - mean(fp_global), 3)), by = .(hierarchy_bin, dataset)],
            hierarchy_bin ~ dataset, value.var = "d"))
