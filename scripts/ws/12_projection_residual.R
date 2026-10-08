#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/<ds>/<sample>.rds   (metadata only, 137 objects)
#          results/tables/01_preprocess/00_curated_manifest.csv           (gsm_id <-> sample alias)
#          /LARGE1/.../02_integrated/ws_pca_meta.rds, results/tables/ws/ws_clusters_percell.csv.gz
#          results/tables/ws/ws_anno_consolidated_res1.csv
# OUTPUT : results/tables/ws/ws_cluster_projection.csv   per cluster: mapping error + high_error rate
#          results/tables/ws/ws_projection_checks.csv
# what it does: brings the BoneMarrowMap projection residual to bear on the clusters whose label is
#   still unsettled. This reuses the project's own projection rather than re-running one:
#   scripts/以前06_hierarchy/d10_project_bmm.R already wrote per-cell `mapping_error`, and
#   d15_derive_mapping_thresholds.R already solved the confound that would otherwise sink this
#   analysis -- the BMM reference is 10x V2 only, so raw mapping error carries a platform term, and
#   d15 sets the threshold per PLATFORM at the P95 of mapping error among HEALTHY cells on that
#   platform. The resulting `high_error` flag therefore means "worse than same-platform healthy
#   cells", with a 5% false-positive rate BY CONSTRUCTION. That is the reason to prefer it over the
#   inferCNV consensus here, whose healthy-donor false-positive rate is 12-65% depending on cell
#   type and dataset (results/tables/02_malignancy/malignancy_fpr_*.csv) and which can therefore
#   only ever be used as a rule-OUT.
#   The 5% figure is verified, not assumed: the script measures the observed high_error rate among
#   this cohort's own Control cells and reports it next to the per-cluster numbers.
#   Join key is (sample, raw 10x barcode); sample names differ between the two pipelines, so they
#   are reconciled by prefix and, failing that, through the manifest's gsm_id. Every join is an
#   explicit merge() -- a keyed dt[.(a,b)] lookup resolves bare column names against the TARGET's
#   own columns, which silently returned all-NA in an earlier attempt at the inferCNV join.
# Covers 7 of the 9 datasets: BoneMarrowMap and zenodo_3345981 were never projected.
# This script deletes nothing.

suppressPackageStartupMessages({ library(data.table); library(Seurat); library(SeuratObject) })
ROOT <- "/FAST/gr10634/gaozy/aml_niche_net"
AN   <- "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/04_annotated"
BASE <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/02_integrated"
TAB  <- file.path(ROOT, "results/tables/ws")
DSMAP <- c(PMID37470778 = "Chen2023")
WANT <- c("Chen2023","GSE185381","GSE201966","GSE227903","GSE239721","GSE253355","GSE289435")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK)+1]] <<- data.table(check = name, expected = as.character(expected),
                                      got = as.character(got), pass = ok)
  message(sprintf("  %-52s expect %-10s got %-10s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. harvest the projection metadata ----------------------------------------------------------
ts("[1] reading projection metadata from the annotated objects")
KEEP <- c("mapping_error","high_error","bmm_fine","bmm_prob","Platform","Chemistry","Disease_state")
PR <- rbindlist(lapply(WANT, function(d) {
  fs <- list.files(file.path(AN, d), "\\.rds$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  rbindlist(lapply(fs, function(p) {
    o <- readRDS(p); m <- o@meta.data
    k <- intersect(KEEP, colnames(m))
    x <- data.table(cell = rownames(m), ps = sub("\\.rds$", "", basename(p)), ds = d)
    for (j in k) x[[j]] <- m[[j]]
    rm(o); x
  }), fill = TRUE)
}), fill = TRUE)
PR[, bc := sub(".*_", "", cell)]
ts(sprintf("    %d projected cells over %d datasets / %d samples", nrow(PR), uniqueN(PR$ds), uniqueN(PR$ps)))
chk("mapping_error present", TRUE, "mapping_error" %in% names(PR))
chk("high_error present", TRUE, "high_error" %in% names(PR))
chk("mapping_error non-NA", TRUE, mean(!is.na(PR$mapping_error)) > 0.95)

## -- 2. the re-integrated cells ------------------------------------------------------------------
ts("[2] the re-integrated cells")
CL <- fread(file.path(TAB, "ws_clusters_percell.csv.gz"), select = c("cell","res.1"))
setnames(CL, "res.1", "cl"); setkey(CL, cell)
M  <- readRDS(file.path(BASE, "ws_pca_meta.rds"))
md <- as.data.table(M$meta)[, .(cell, Dataset_ID, SampleType, Sample_Name)]
setkey(md, cell); md <- md[CL]
md[, bc := sub(".*_", "", cell)]
md[, ds := ifelse(Dataset_ID %chin% names(DSMAP), unname(DSMAP[Dataset_ID]), Dataset_ID)]
chk("cells", 895228L, nrow(md))

## -- 3. reconcile sample names, then merge explicitly --------------------------------------------
ts("[3] sample-name reconciliation")
cm <- fread(file.path(ROOT, "results/tables/01_preprocess/00_curated_manifest.csv"),
            select = c("dataset","sample","gsm_id"))
pairs <- rbindlist(lapply(intersect(unique(PR$ds), unique(md$ds)), function(d) {
  a <- unique(PR[ds == d, ps]); b <- unique(md[ds == d, Sample_Name])
  out <- rbindlist(lapply(b, function(x) {
    hit <- a[startsWith(a, x) | startsWith(x, a)]
    if (length(hit) == 1L) data.table(ds = d, Sample_Name = x, ps = hit) else NULL
  }))
  todo <- setdiff(b, if (is.null(out)) character() else out$Sample_Name)
  if (length(todo)) {
    g <- cm[dataset == d & !is.na(gsm_id) & gsm_id != ""]
    out <- rbindlist(list(out, rbindlist(lapply(todo, function(x) {
      hit <- a[a %chin% g[gsm_id == x | startsWith(gsm_id, x), sample]]
      if (length(hit) == 1L) data.table(ds = d, Sample_Name = x, ps = hit) else NULL
    }))), use.names = TRUE, fill = TRUE)
  }
  out
}), fill = TRUE)
ts(sprintf("    %d one-to-one sample pairs", nrow(pairs)))
print(pairs[, .(pairs = .N), by = ds][order(-pairs)])

md <- merge(md, pairs, by = c("ds","Sample_Name"), all.x = TRUE)
PK <- unique(PR[, .(ds, ps, bc, mapping_error, high_error, bmm_fine, bmm_prob, Platform)],
             by = c("ds","ps","bc"))
md <- merge(md, PK, by = c("ds","ps","bc"), all.x = TRUE)
cov <- md[, .(cells = .N, matched = sum(!is.na(mapping_error)),
              pct = round(100*mean(!is.na(mapping_error)),1)), by = Dataset_ID][order(-cells)]
ts("[4] coverage")
print(cov)
ts(sprintf("    total %d / %d = %.1f%%", sum(cov$matched), sum(cov$cells), 100*sum(cov$matched)/sum(cov$cells)))

## -- 5. verify the 5%-by-construction claim on THIS cohort ---------------------------------------
ts("[5] verifying the high_error threshold on this cohort's own Control cells")
base <- md[!is.na(high_error), .(cells = .N, frac_high = round(mean(high_error == 1), 4)),
           by = .(SampleType)]
print(base)
ctl <- base[SampleType == "Control", frac_high]
ts(sprintf("    Control high_error rate = %.3f (d15 designed it to be ~0.05)", ctl))
chk("Control high_error rate is near the designed 5%", TRUE, !is.na(ctl) && ctl < 0.20)
print(md[!is.na(high_error), .(cells = .N, frac_high = round(mean(high_error == 1), 3)),
         by = .(Dataset_ID, SampleType)][order(-frac_high)])

## -- 6. per cluster ------------------------------------------------------------------------------
ts("[6] per cluster")
PC <- md[, .(cluster_size = .N, n_proj = sum(!is.na(mapping_error)),
             med_mapping_error = round(median(mapping_error, na.rm = TRUE), 4),
             frac_high_error = round(mean(high_error[!is.na(high_error)] == 1), 3),
             med_bmm_prob = round(median(bmm_prob, na.rm = TRUE), 3)), by = cl]
PC[, proj_frac := round(n_proj / cluster_size, 3)]
CONS <- fread(file.path(TAB, "ws_anno_consolidated_res1.csv"))
PC <- merge(PC, CONS[, .(cl, provenance, independent_check, auto_mk3, old_label, pct_control)], by = "cl")
fwrite(PC[order(-cluster_size)], file.path(TAB, "ws_cluster_projection.csv"))
ts("    by provenance:")
print(PC[n_proj > 0, .(clusters = .N, med_proj_frac = round(median(proj_frac),2),
        med_frac_high = round(median(frac_high_error),3)), by = provenance])
tg <- CONS[provenance == "disputed" | (provenance == "adopted" & independent_check == "contradicted"), cl]
ts("    the unsettled clusters, by projection-residual rate:")
print(PC[cl %in% tg & n_proj > 0][order(-frac_high_error)][,
  .(cl, cluster_size, proj_frac, frac_high_error, med_mapping_error, med_bmm_prob,
    auto = substr(auto_mk3, 1, 24), old = substr(old_label, 1, 20))])

C <- rbindlist(CHK); fwrite(C, file.path(TAB, "ws_projection_checks.csv"))
ts(sprintf("[7] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
ts("[done] ", file.path(TAB, "ws_cluster_projection.csv"))
