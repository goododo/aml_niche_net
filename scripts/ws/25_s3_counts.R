#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../02_seurat_objects/04_annotated/<dataset>/<sample>.rds   (metadata only: bmm_fine, high_error,
#          anno_source, GSM, and for GSE185381 anno_age_group / anno_sex)
#          results/tables/01_preprocess/00_curated_manifest.csv              (disease, timepoint, tissue, sorting, batch notes)
#          scripts/ws/s4_cell_groups.tsv                                      (PROPOSED bmm_fine -> S4 group map; frozen at S4)
#          results/tables/ws/decontx/GSE185381/ws_decontx_index.csv.gz       (cell -> channel, for the same-channel count)
# OUTPUT : results/tables/ws/ws_s3_cells_per_group.csv        sample x group: cells (high_error dropped) and concordant cells
#          results/tables/ws/ws_s3_samples.csv                one row per sample: role, batch, age group, per-group counts (wide)
#          results/tables/ws/ws_s3_study_summary.csv          per study and group: samples with >= 30 cells; groups retained at 70%
#          results/tables/ws/ws_s3_gse185381_channels.csv     per 10x channel: adult AML / paediatric / healthy donors (S6)
#          results/tables/ws/ws_s3_gse227903_pairs.csv        the Dx -> MRD pairs: batch and per-group cells
#          results/tables/ws/ws_s3_unmapped_fine_labels.csv   fine labels inside T_NK / B_Plasma / Mono_DC not in the map
#          results/tables/ws/ws_s3_checks.csv
# what it does: v2.0 S3, the counting step before S4. Three questions, no statistics: (1) how many patients per study
#   carry >= 30 cells in each proposed S4 cell group, with and without requiring marker-projection concordance;
#   (2) how many GSE185381 channels (and pools) hold an adult AML donor together with a healthy donor (S6), and
#   what the healthy donors' age groups are; (3) for the GSE227903 Dx -> MRD pairs, the batch of each sample and
#   its per-group cells. Counting rules: high_error cells dropped (config_ws.R section 0); Early Lymphoid and the
#   other EXCLUDED rows of the map never count (section 11); "concordant" = anno_source == "concordant".
#   usage: Rscript scripts/ws/25_s3_counts.R
# This script reads metadata only and modifies nothing.

suppressPackageStartupMessages({ library(data.table); library(Seurat); library(SeuratObject) })
suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "ws", "config_ws.R"))
AN <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated")
TABWS <- file.path(TAB_DIR, "ws")
DATASETS <- c("GSE185381", "GSE227903", "GSE289435", "GSE116256", "GSE201966")
MIN_CELLS <- 30L; RETAIN_FRAC <- 0.7                      # v2.0 5.3 S4 item 4
MAP <- fread(here::here("scripts", "ws", "s4_cell_groups.tsv"), fill = TRUE); GROUPS <- setdiff(unique(MAP$group), "EXCLUDED")   # note column is often empty
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  got <- tryCatch(got, error = function(e) paste("ERROR:", conditionMessage(e)))
  ok  <- tryCatch(isTRUE(all.equal(expected, got)), error = function(e) FALSE)
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = as.character(expected), got = as.character(got), pass = ok)
  message(sprintf("  %-60s expect %-8s got %-8s %s", name, expected, got, if (ok) "PASS" else "** FAIL **"))
}

## -- 1. metadata from every object -------------------------------------------------------------------------
ts("[1] reading metadata")
MD <- rbindlist(lapply(DATASETS, function(d) {
  fs <- list.files(file.path(AN, d), pattern = "\\.rds$", full.names = TRUE); if (!length(fs)) { ts("    ", d, ": no objects"); return(NULL) }
  rbindlist(lapply(fs, function(f) { o <- readRDS(f); md <- o@meta.data
    data.table(dataset = d, sample = sub("\\.rds$", "", basename(f)), cell = colnames(o),
               bin = as.character(md$hierarchy_bin), bmm_fine = as.character(md$bmm_fine), high_error = as.integer(md$high_error),
               anno_source = as.character(md$anno_source),
               age_group = if ("anno_age_group" %in% names(md)) as.character(md$anno_age_group) else NA_character_,
               sex = if ("anno_sex" %in% names(md)) as.character(md$anno_sex) else NA_character_) }))
}), fill = TRUE)
ts(sprintf("    %d cells, %d samples, %d datasets", nrow(MD), uniqueN(MD$sample), uniqueN(MD$dataset)))
for (d in c("GSE185381", "GSE227903", "GSE289435", "GSE116256")) {
  ix <- file.path(TABWS, "decontx", d, "ws_decontx_index.csv.gz")
  if (file.exists(ix)) chk(paste("cells equal the S1 index,", d), nrow(fread(ix, select = "cell")), MD[dataset == d, .N])
}
MD <- merge(MD, MAP[, .(bmm_fine, group)], by = "bmm_fine", all.x = TRUE)
UNM <- MD[bin %in% c("T_NK", "B_Plasma", "Mono_DC") & is.na(group) & !is.na(bmm_fine) & bmm_fine != "", .N, by = .(bin, bmm_fine)][order(-N)]
fwrite(UNM, file.path(TABWS, "ws_s3_unmapped_fine_labels.csv"))
chk("every fine label in the three bins is in the map", 0L, nrow(UNM))
CELLS <- MD[!is.na(group) & group != "EXCLUDED" & !(high_error %in% 1L)]
chk("Early Lymphoid never counts", 0L, CELLS[bmm_fine %in% c("CLP", "MLP", "MLP-II", "Pre-ProB"), .N])

## -- 2. roles from the manifest --------------------------------------------------------------------------------
man <- fread(file.path(TAB_DIR, "01_preprocess", "00_curated_manifest.csv"))
R <- man[dataset %in% DATASETS, .(dataset, sample, disease, timepoint, tissue_r, sorting_r, blast_pct_clinical, library_id,
          batch = { m <- regmatches(paste(open_questions, curator_note), regexpr("nuig_batch[0-9]", paste(open_questions, curator_note))); if (length(m)) m else NA_character_ }), by = .I][, -"I"]
R[, role := fcase(disease == "AML" & timepoint == "Diagnosis" & tissue_r == "BM" & sorting_r %in% c("viability_only", "none"), "AML_Dx",
                  disease == "AML" & timepoint == "Diagnosis", "AML_Dx_other_prep",
                  disease == "AML" & grepl("Post_treatment|Remission|MRD", timepoint), "AML_posttreatment",
                  disease == "AML" & timepoint == "Relapse", "AML_relapse",
                  disease != "AML", "control", default = "other")]
R[, paediatric := dataset == "GSE185381" & disease == "AML" & !grepl("^AML[0-9]{4}$", sample)]     # config_ws.R section 10 rule
chk("GSE185381 paediatric samples", 22L, R[dataset == "GSE185381" & paediatric == TRUE, .N])

## -- 3. cells per sample and group ---------------------------------------------------------------------------
ts("[3] cells per group")
CG <- CELLS[, .(n_all = .N, n_conc = sum(anno_source %in% "concordant")), by = .(dataset, sample, group)]
grid <- CJ(sample = unique(MD$sample), group = GROUPS); grid <- merge(grid, unique(MD[, .(sample, dataset)]), by = "sample")
CG <- merge(grid, CG, by = c("dataset", "sample", "group"), all.x = TRUE); CG[is.na(n_all), `:=`(n_all = 0L, n_conc = 0L)]
fwrite(CG[order(dataset, sample, group)], file.path(TABWS, "ws_s3_cells_per_group.csv"))
W <- dcast(CG, dataset + sample ~ group, value.var = c("n_all", "n_conc"))
W[, groups_ge30_all := rowSums(.SD >= MIN_CELLS), .SDcols = paste0("n_all_", GROUPS)]
W[, groups_ge30_conc := rowSums(.SD >= MIN_CELLS), .SDcols = paste0("n_conc_", GROUPS)]
AG <- MD[!is.na(age_group), .(age_group = names(which.max(table(age_group))), sex = names(which.max(table(sex)))), by = sample]
S <- merge(R[, .(dataset, sample, role, paediatric, batch, blast_pct_clinical, library_id)], W, by = c("dataset", "sample"), all.y = TRUE)
S <- merge(S, AG, by = "sample", all.x = TRUE)
S[is.na(role), role := "not_in_manifest"]
fwrite(S[order(dataset, role, sample)], file.path(TABWS, "ws_s3_samples.csv"))

## -- 4. study summary: AML diagnosis, adults ---------------------------------------------------------------------
ts("[4] study summary")
DX <- CG[sample %in% S[role == "AML_Dx" & paediatric == FALSE, sample]]
SUM <- DX[, .(n_samples = uniqueN(sample), ge30_all = sum(n_all >= MIN_CELLS), ge30_conc = sum(n_conc >= MIN_CELLS),
              median_all = as.numeric(median(n_all)), median_conc = as.numeric(median(n_conc))), by = .(dataset, group)]
SUM[, `:=`(frac_all = round(ge30_all / n_samples, 2), frac_conc = round(ge30_conc / n_samples, 2))]
SUM[, `:=`(retained_all = frac_all >= RETAIN_FRAC, retained_conc = frac_conc >= RETAIN_FRAC)]
fwrite(SUM[order(dataset, group)], file.path(TABWS, "ws_s3_study_summary.csv"))
for (d in unique(SUM$dataset)) {
  x <- SUM[dataset == d]; ts(sprintf("    %-10s adult AML Dx n=%2d | groups >= 30 in >= 70%% of samples: all-cells {%s} | concordant {%s}", d, x$n_samples[1],
     paste(x[retained_all == TRUE, group], collapse = ","), paste(x[retained_conc == TRUE, group], collapse = ",")))
}
CT <- CG[sample %in% S[role == "control", sample]][, .(n = uniqueN(sample), ge30_all = sum(n_all >= MIN_CELLS), ge30_conc = sum(n_conc >= MIN_CELLS)), by = .(dataset, group)]
ts("    controls, samples with >= 30 cells (all/concordant): ", paste(CT[, sprintf("%s %s %d/%d", dataset, group, ge30_all, ge30_conc)], collapse = "; "))

## -- 5. GSE185381 channels and ages (S6) --------------------------------------------------------------------------
ts("[5] GSE185381 channels")
IX <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_decontx_index.csv.gz"), select = c("sample", "unit"))
IX <- merge(IX[, .N, by = .(sample, unit)], S[dataset == "GSE185381", .(sample, role, paediatric, library_id)], by = "sample")
IX[, kind := fcase(role == "control", "healthy", paediatric == TRUE, "AML_paediatric", role == "AML_Dx", "AML_adult", default = role)]
CH <- IX[, .(n_donors = .N, adult_AML = sum(kind == "AML_adult"), paediatric_AML = sum(kind == "AML_paediatric"), healthy = sum(kind == "healthy"),
             donors = paste(sort(unique(sample)), collapse = ";"), pools = paste(sort(unique(library_id)), collapse = ";")), by = unit][order(-n_donors)]
CH[, s6_adult_AML_with_healthy := adult_AML > 0 & healthy > 0]
fwrite(CH, file.path(TABWS, "ws_s3_gse185381_channels.csv"))
PL <- IX[, .(adult_AML = uniqueN(sample[kind == "AML_adult"]), healthy = uniqueN(sample[kind == "healthy"])), by = library_id]
ts(sprintf("    channels: %d | with adult AML + healthy together: %d (adult AML donors in them %d, healthy %d) | pools with both: %d",
           nrow(CH), sum(CH$s6_adult_AML_with_healthy), CH[s6_adult_AML_with_healthy == TRUE, sum(adult_AML)], CH[s6_adult_AML_with_healthy == TRUE, sum(healthy)], PL[adult_AML > 0 & healthy > 0, .N]))
ts("    age groups, healthy: ", paste(S[dataset == "GSE185381" & role == "control", sprintf("%s=%s", sample, age_group)], collapse = " "))
ts("    age groups, adult AML: ", paste(names(table(S[dataset == "GSE185381" & role == "AML_Dx" & paediatric == FALSE, age_group])), table(S[dataset == "GSE185381" & role == "AML_Dx" & paediatric == FALSE, age_group]), collapse = " "))

## -- 6. GSE227903 Dx -> MRD pairs --------------------------------------------------------------------------------
ts("[6] GSE227903 pairs")
G <- S[dataset == "GSE227903"]; G[, `:=`(patient = sub("_.*", "", sample), tp = sub(".*_", "", sample))]
pairs <- G[tp == "MRD", patient]; PR <- G[patient %in% pairs & tp %in% c("Dg", "MRD")][order(patient, tp)]
PR[, cross_batch := uniqueN(batch) > 1, by = patient]
fwrite(PR[, c("patient", "tp", "sample", "batch", "cross_batch", "groups_ge30_all", "groups_ge30_conc", paste0("n_all_", GROUPS), paste0("n_conc_", GROUPS)), with = FALSE],
       file.path(TABWS, "ws_s3_gse227903_pairs.csv"))
chk("GSE227903 Dx -> MRD pairs", 7L, length(pairs))
ts("    cross-batch pairs: ", paste(unique(PR[cross_batch == TRUE, patient]), collapse = " "), " | groups >= 30 per sample (all): ",
   paste(PR[, sprintf("%s %d", sample, groups_ge30_all)], collapse = ", "))
C <- rbindlist(CHK); fwrite(C, file.path(TABWS, "ws_s3_checks.csv"))
ts(sprintf("[done] %d checks, %d PASS, %d FAIL", nrow(C), sum(C$pass), sum(!C$pass)))
if (any(!C$pass)) print(C[pass == FALSE])
