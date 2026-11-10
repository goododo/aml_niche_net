#!/usr/bin/env Rscript
# 02b_bmm_label_bc.R ----
# Doublet-pairing SMOKE, recall lever 2 (user decision 2026-11-10): label class B/C barcodes by BoneMarrowMap projection, the
# same method that made the class-A labels, instead of stage 02's kNN transfer (state accuracy 0.85 among assigned, 35 %
# unassigned in the pool). Every analysis barcode of a lane is projected; class-A cells keep their stored labels, and their
# re-projection is the known answer (it must agree with the stored bmm_fine).
# INPUT  : DBL_LANE_DIR/<lane>.rds (01), <lane>__labels.csv.gz (02); BMM scaffold + uwot model (config_hierarchy.R)
# OUTPUT : DBL_LANE_DIR/<lane>__labels_bmm.csv.gz  same columns as stage 02's labels file; B/C compartments from the projection
#          DBL_TAB_DIR/dbl_checks_02b_smoke.csv
# B/C rule (stage 02's class-A rule minus the parts that need author/marker data): bmm_fine -> dbl_label_map.tsv; a LEUK state
#   or T label with bmm_prob < 0.5 -> LEUK_other / Tamb; T also needs T-identity UMI >= 1 (no marker concordance exists for
#   B/C); projection failed (bmm_fine NA) -> "Unassigned".
# Usage  : Rscript scripts/dbl/02b_bmm_label_bc.R [--force]
# Masking: projection labels of single barcodes, no doublet typing.

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(BoneMarrowMap); library(data.table); library(Matrix); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "dbl", "config_dbl.R"))
FORCE <- "--force" %in% commandArgs(TRUE)
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-80s expect %-10s got %-10s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
MAP <- fread(here::here("scripts", "dbl", DBL_LABEL_MAP_FILE))
ref <- readRDS(BMM_SCAFFOLD_RDS); ref$save_uwot_path <- BMM_UWOT
for (L in DBL_SMOKE_LANES) {
  out <- file.path(DBL_LANE_DIR, paste0(L, "__labels_bmm.csv.gz"))
  if (file.exists(out) && !FORCE) { ts("[skip] ", out); next }
  ts("[1] lane ", L)
  x <- readRDS(file.path(DBL_LANE_DIR, paste0(L, ".rds"))); M <- x$counts
  lab <- fread(file.path(DBL_LANE_DIR, paste0(L, "__labels.csv.gz")))[match(colnames(M), cell)]
  md <- data.frame(row.names = colnames(M), class = lab$class)
  q <- map_Query(exp_query = M, metadata_query = md, ref_obj = ref, do_normalize = TRUE, do_umap = FALSE, verbose = FALSE)
  q <- calculate_MappingError(q, reference = ref, MAD_threshold = MAD_THRESHOLD_BMM)
  q <- predict_CellTypes(query_obj = q, ref_obj = ref, ref_label = BMM_REF_LABEL, k = BMM_KNN_K, include_broad = TRUE)
  qm <- as.data.table(q@meta.data, keep.rownames = "cell")
  pr <- qm[, .(cell, fine = as.character(predicted_CellType), prob = as.numeric(predicted_CellType_prob))]
  ## known answer: re-projected class-A cells agree with their stored labels
  A <- merge(lab[class == "A", .(cell)], x$ann[, .(cell, bmm_fine)], by = "cell")[pr, on = "cell", nomatch = NULL]
  agree <- A[!is.na(bmm_fine), mean(fine == bmm_fine)]
  ts(sprintf("    class-A re-projection agrees with the stored bmm_fine for %.3f of cells", agree))
  chk(sprintf("%s: class-A re-projection agrees with stored labels (>= 0.95)", L), TRUE, agree >= 0.95)
  ## B/C labels
  tid <- Matrix::colSums(M[intersect(DBL_T_ID_GENES, rownames(M)), , drop = FALSE]); names(tid) <- colnames(M)
  bc <- merge(lab[class != "A", .(cell)], pr, by = "cell")
  bc <- merge(bc, MAP[, .(fine = bmm_fine, comp0 = compartment)], by = "fine", all.x = TRUE)
  bc[, compartment := fifelse(is.na(fine), "Unassigned", comp0)]
  bc[compartment %in% DBL_LEUK_STATES & prob < DBL_MIN_BMM_PROB, compartment := "LEUK_other"]
  bc[compartment == "T_cand", compartment := fifelse(prob >= DBL_MIN_BMM_PROB & tid[cell] >= DBL_T_ID_MIN_UMI, "T", "Tamb")]
  chk(sprintf("%s: every B/C barcode got a compartment", L), 0L, bc[is.na(compartment), .N])
  ts(sprintf("    B/C: %d barcodes, %.3f without a projection label", nrow(bc), bc[, mean(compartment == "Unassigned")]))
  new <- copy(lab)
  new[bc, on = "cell", `:=`(compartment = i.compartment, bmm_fine = i.fine, source = "bmm_projection", vote = i.prob)]
  chk(sprintf("%s: class-A rows unchanged", L), TRUE, identical(new[class == "A", compartment], lab[class == "A", compartment]))
  fwrite(new, out); ts("    [write] ", out)
  rm(x, M, q); gc(FALSE)
}
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_02b_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
