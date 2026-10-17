#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds      (raw counts; bmm_fine, anno_source, high_error)
#          /LARGE1/.../aml_niche_net/ws_decontx/<dataset>/<unit>.rds              (22's decontaminated counts; sensitivity)
#          results/tables/ws/decontx/<dataset>/ws_decontx_index.csv.gz           (cell -> channel)
#          results/tables/ws/ws_doublet_scores_percell.csv.gz                     (scDblFinder score; NA where T3a never ran)
#          results/tables/01_preprocess/00_curated_manifest.csv; scripts/ws/s4_cell_groups.tsv; config_ws.R section 12
# OUTPUT : /LARGE1/.../aml_niche_net/ws_pseudobulk/<dataset>/ws_pb_<dataset>.rds
#            $raw[[group4]]   genes x samples raw-count sums        $dx[[group4]]  the same from decontaminated counts
#            $fine[[bmm_fine]] genes x samples raw sums per subtype  $cells  per-cell assignment (cell, sample, group4, fine, kept)
#          results/tables/ws/ws_s4_pseudobulk_meta_<dataset>.csv   sample x group4: cells, median UMI, sample covariates
#          results/tables/ws/ws_s4_subtype_composition_<dataset>.csv  sample x group4 x subtype: cells and share
#          results/tables/ws/ws_s4_checks_<dataset>.csv
# what it does: v2.0 S4 under the rules frozen in config_ws.R section 12. For every sample it keeps the cells whose
#   BoneMarrowMap fine label maps to one of the four groups (T, NK, B, Mono), whose marker typing agrees with the
#   projection, and that are not high_error; sums their raw counts per group and per fine subtype; sums 22's
#   decontaminated counts per group for the sensitivity arm; and records the sample-level covariates S5 will use:
#   the two raw ambient shares (section 16.1 item 13), doublet score mean, paediatric flag, age group, sex (metadata
#   or XIST/RPS4Y1), channels and pool. The per-cell assignment is kept so that S5b can score single cells and
#   summarise them per patient. Nothing is tested here.
#   usage: --dataset=GSE185381   (one dataset per run)
# This script deletes nothing and does not modify its inputs.

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(Seurat); library(SeuratObject) })
.a <- commandArgs(TRUE)
.get <- function(k, d) { h <- grep(paste0("^--", k, "="), .a, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
DATASET <- .get("dataset", "GSE185381")

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
AN   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated", DATASET)
DXD  <- file.path(LARGE1_DIR, "ws_decontx", DATASET)
OUTD <- file.path(LARGE1_DIR, "ws_pseudobulk", DATASET); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
TABWS <- file.path(TAB_DIR, "ws")
MAP <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE)
HB <- c("HBB", "HBA1", "HBA2"); GRAN <- c("MPO", "ELANE", "AZU1", "PRTN3"); EL_FINE <- c("CLP", "MLP", "MLP-II", "Pre-ProB")
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}
stopifnot(setequal(setdiff(unique(MAP$group4), "EXCLUDED"), WS_S4_GROUPS))

## -- 1. sample roles, channels, doublets --------------------------------------------------------------------
ts("[1] inputs, ", DATASET)
man <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))[dataset == DATASET]
man[, role := fcase(disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" & sorting_r %in% c("viability_only", "none"), "AML_Dx",
                    disease == "AML" & timepoint == "Diagnosis", "AML_Dx_other_prep",
                    disease == "AML" & grepl("Post_treatment|Remission|MRD", timepoint), "AML_posttreatment",
                    disease == "AML" & timepoint == "Relapse", "AML_relapse", disease != "AML", "control", default = "other")]
man[, paediatric := DATASET == "GSE185381" & disease == "AML" & !grepl("^AML[0-9]{4}$", sample)]
man[, batch := { m <- regmatches(paste(open_questions, curator_note), regexpr("nuig_batch[0-9]", paste(open_questions, curator_note))); if (length(m)) m else NA_character_ }, by = sample]
IXF <- file.path(TABWS, "decontx", DATASET, "ws_decontx_index.csv.gz"); HAVE_DX <- file.exists(IXF)
IX <- if (HAVE_DX) fread(IXF, select = c("cell", "sample", "unit")) else NULL
DBL <- fread(file.path(TABWS, WS_DBL_SCORES_TABLE))[dataset == DATASET, .(cell, sc_score)]; setkey(DBL, cell)
fs <- list.files(AN, pattern = "\\.rds$", full.names = TRUE); stopifnot(length(fs) > 0)

## -- 2. per sample ------------------------------------------------------------------------------------------
ts("[2] ", length(fs), " samples")
RAW <- list(); DX <- list(); FINE <- list(); CELLS <- list(); META <- list(); COMP <- list(); G <- NULL
dx_cache <- new.env()
get_dx <- function(unit) { if (!exists(unit, dx_cache)) assign(unit, readRDS(file.path(DXD, paste0(unit, ".rds"))), dx_cache); get(unit, dx_cache) }
for (f in fs) {
  sm <- sub("\\.rds$", "", basename(f)); o <- readRDS(f); md <- o@meta.data
  X <- LayerData(o, assay = "RNA", layer = "counts"); if (is.null(G)) G <- rownames(X)
  if (!identical(rownames(X), G)) stop("gene axis differs in ", sm)
  tot <- Matrix::colSums(X); bin <- as.character(md$hierarchy_bin); fine <- as.character(md$bmm_fine)
  m <- MAP[match(fine, bmm_fine)]
  kept <- !is.na(m$group4) & m$group4 != "EXCLUDED" & !(WS_S4_DROP_HIGH_ERROR & md$high_error %in% 1L) &
          (!WS_S4_REQUIRE_CONCORDANT | md$anno_source %in% "concordant")
  CELLS[[sm]] <- data.table(cell = colnames(X), sample = sm, fine = fine, group = m$group, group4 = m$group4, kept = kept,
                            high_error = as.integer(md$high_error), anno_source = as.character(md$anno_source))
  # sample-level covariates on ALL cells, with the S1 definitions (23_decontx_summary.R)
  grp_el <- fine %in% EL_FINE; nonery <- !is.na(bin) & bin != "Erythroid"; lymph <- bin %in% c("T_NK", "B_Plasma") & !grp_el
  hb <- intersect(HB, G); gr <- intersect(GRAN, G)
  sexg <- WS_S4_SEX_GENES; xist <- if (sexg[["female"]] %in% G) sum(X[sexg[["female"]], ]) else NA; rps <- if (sexg[["male"]] %in% G) sum(X[sexg[["male"]], ]) else NA
  sex_meta <- if ("anno_sex" %in% names(md)) names(which.max(table(md$anno_sex))) else NA_character_
  age <- if ("anno_age_group" %in% names(md)) names(which.max(table(md$anno_age_group))) else NA_character_
  d <- DBL[J(colnames(X))]$sc_score
  units <- if (HAVE_DX) IX[sample == sm, unit] else character(0)
  base <- data.table(dataset = DATASET, sample = sm, n_cells = ncol(X), n_kept = sum(kept), kept_umi = sum(tot[kept]),
                     hb_share_nonery = if (any(nonery)) sum(X[hb, nonery]) / sum(tot[nonery]) else NA_real_,
                     gran_share_lymph = if (any(lymph)) sum(X[gr, lymph]) / sum(tot[lymph]) else NA_real_,
                     n_lymph_for_share = sum(lymph), dbl_mean = if (any(!is.na(d))) mean(d, na.rm = TRUE) else NA_real_,
                     sex_meta = sex_meta, xist_umi = xist, rps4y1_umi = rps,
                     sex_inferred = fifelse(is.na(xist) | is.na(rps), NA_character_, fifelse(xist > 3 * pmax(rps, 1), "female", fifelse(rps > 3 * pmax(xist, 1), "male", "unclear"))),
                     age_group = age, channels = paste(sort(unique(units)), collapse = ";"), n_channels = length(unique(units)))
  # pseudobulk per group4 and per fine subtype on the kept cells
  dxm <- NULL
  if (HAVE_DX && length(units)) {
    parts <- lapply(unique(units), function(u) { R <- get_dx(u); j <- match(intersect(R$cells, colnames(X)), R$cells); M <- R$decontXcounts[, j, drop = FALSE]; colnames(M) <- R$cells[j]; M })
    dxm <- do.call(cbind, parts); if (!identical(rownames(dxm), G)) dxm <- dxm[match(G, rownames(dxm)), , drop = FALSE]
  }
  for (g4 in WS_S4_GROUPS) {
    j <- which(kept & m$group4 == g4)
    RAW[[g4]][[sm]] <- if (length(j)) as.numeric(Matrix::rowSums(X[, j, drop = FALSE])) else numeric(length(G))
    DX[[g4]][[sm]]  <- if (length(j) && !is.null(dxm)) { jj <- match(colnames(X)[j], colnames(dxm)); as.numeric(Matrix::rowSums(dxm[, jj[!is.na(jj)], drop = FALSE])) } else rep(NA_real_, length(G))
    META[[paste(sm, g4)]] <- cbind(base, data.table(group4 = g4, n = length(j), median_umi = if (length(j)) as.numeric(median(tot[j])) else NA_real_,
                                                    eligible = length(j) >= WS_S4_MIN_CELLS))
    if (length(j)) COMP[[paste(sm, g4)]] <- data.table(dataset = DATASET, sample = sm, group4 = g4, fine = fine[j])[, .(n = .N), by = .(dataset, sample, group4, fine)][, share := n / sum(n)]
  }
  for (fl in unique(m$group4[kept])) for (fn in unique(fine[kept & m$group4 == fl])) {
    j <- which(kept & fine == fn); FINE[[fn]][[sm]] <- as.numeric(Matrix::rowSums(X[, j, drop = FALSE])) }
  ts(sprintf("    %-14s cells %6d kept %5d | %s", sm, ncol(X), sum(kept), paste(sprintf("%s=%d", WS_S4_GROUPS, vapply(WS_S4_GROUPS, function(g) sum(kept & m$group4 == g), 1L)), collapse = " ")))
  rm(o, X, dxm); gc(FALSE)
}
tomat <- function(L) { M <- do.call(cbind, L); rownames(M) <- G; M }
PB <- list(dataset = DATASET, genes = G, raw = lapply(RAW, tomat), dx = lapply(DX, tomat), fine = lapply(FINE, tomat),
           cells = rbindlist(CELLS), rules = list(date = WS_S4_DATE, min_cells = WS_S4_MIN_CELLS, concordant = WS_S4_REQUIRE_CONCORDANT, map = MAP))
META <- rbindlist(META); META <- merge(META, man[, .(sample, role, paediatric, batch, library_id, blast_pct_clinical)], by = "sample", all.x = TRUE)
META[is.na(role), role := "not_in_manifest"]
COMP <- rbindlist(COMP)

## -- 3. checks ----------------------------------------------------------------------------------------------
ts("[3] checks")
chk("samples in the pseudobulk", length(fs), ncol(PB$raw$T))
s3 <- fread(file.path(TABWS, "ws_s3_cells_per_group.csv"))[dataset == DATASET]
s3g <- merge(s3, unique(MAP[group != "EXCLUDED", .(group, group4)]), by = "group")[group4 != "EXCLUDED", .(n_s3 = sum(n_conc)), by = .(sample, group4)]
cmp <- merge(META[, .(sample, group4, n)], s3g, by = c("sample", "group4"))
chk("kept cells per group equal S3's concordant counts", 0L, cmp[n != n_s3, .N])
tot4 <- Reduce(`+`, lapply(PB$raw, colSums)); ku <- unique(META[, .(sample, kept_umi)]); ku <- ku[match(names(tot4), sample)]
chk("the four group sums add up to the kept cells' UMIs", TRUE, max(abs(tot4 - ku$kept_umi)) < 1e-6)   # independent of the matrices' own colSums
cpm <- function(M, g) { v <- M[g, ] / pmax(colSums(M), 1) * 1e6; v[colSums(M) >= 1000] }
chk("CD3E higher in T than in Mono pseudobulk", TRUE, median(cpm(PB$raw$T, "CD3E")) > 10 * median(cpm(PB$raw$Mono, "CD3E")))
chk("CD79A higher in B than in T pseudobulk", TRUE, median(cpm(PB$raw$B, "CD79A")) > 10 * median(cpm(PB$raw$T, "CD79A")))
chk("LYZ higher in Mono than in T pseudobulk", TRUE, median(cpm(PB$raw$Mono, "LYZ")) > 10 * median(cpm(PB$raw$T, "LYZ")))
chk("NKG7 higher in NK than in B pseudobulk", TRUE, median(cpm(PB$raw$NK, "NKG7")) > 10 * median(cpm(PB$raw$B, "NKG7")))
s1 <- file.path(TABWS, "decontx", DATASET, "ws_s1_contamination_per_sample.csv")
if (file.exists(s1)) { S1 <- fread(s1); x <- merge(unique(META[, .(sample, hb_share_nonery)]), S1[, .(sample, hb_share_nonery_raw)], by = "sample")
  chk("HB share equals the S1 table", TRUE, max(abs(x$hb_share_nonery - x$hb_share_nonery_raw)) < 1e-6) }
if (HAVE_DX) chk("decontaminated sums never exceed raw sums (T)", TRUE, all(colSums(PB$dx$T, na.rm = TRUE) <= colSums(PB$raw$T) + 1e-6))
sx <- unique(META[, .(sample, sex_meta, sex_inferred)])[!is.na(sex_meta) & sex_inferred %in% c("female", "male")]
if (nrow(sx)) chk("inferred sex agrees with metadata where both exist", TRUE, mean(sx$sex_meta == sx$sex_inferred) >= 0.9)
el <- META[role == "AML_Dx" & paediatric == FALSE, .(n_samples = uniqueN(sample), eligible = sum(eligible)), by = group4]
el[, retained := eligible / n_samples >= WS_S4_RETAIN_FRAC]
ts("    adult AML Dx eligibility (>= ", WS_S4_MIN_CELLS, " concordant cells): ", paste(el[, sprintf("%s %d/%d%s", group4, eligible, n_samples, ifelse(retained, "", " (below 70%)"))], collapse = "; "))
ts("    controls: ", paste(META[role == "control", .(eligible = sum(eligible), n = .N), by = group4][, sprintf("%s %d/%d", group4, eligible, n)], collapse = "; "))

## -- 4. write -------------------------------------------------------------------------------------------------
saveRDS(PB, file.path(OUTD, paste0("ws_pb_", DATASET, ".rds")))
for (k in c("hb_share_nonery", "gran_share_lymph", "dbl_mean")) META[, (k) := signif(get(k), 5)]
fwrite(META[order(role, sample, group4)], file.path(TABWS, paste0("ws_s4_pseudobulk_meta_", DATASET, ".csv")))
COMP[, share := round(share, 4)]; fwrite(COMP[order(sample, group4, -n)], file.path(TABWS, paste0("ws_s4_subtype_composition_", DATASET, ".csv")))
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, paste0("ws_s4_checks_", DATASET, ".csv")))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL | %s", nrow(C), sum(C$pass), sum(!C$pass), file.path(OUTD, paste0("ws_pb_", DATASET, ".rds"))))
if (any(!C$pass)) print(C[pass == FALSE])
