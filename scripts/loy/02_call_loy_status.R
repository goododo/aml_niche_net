#!/usr/bin/env Rscript
# 02_call_loy_status.R ----
# INPUT  : LARGE1/loy/percell/<dataset>__y_percell.csv.gz (from 01), results/tables/01_preprocess/00_curated_manifest.csv,
#          results/tables/ws/ws_s3_samples.csv (GSE185381 sex / age group), scripts/loy/config_loy.R section 5 (frozen rules)
# OUTPUT : results/tables/loy/02_calls/loy_reference_rates.csv       r_ref per chemistry x compartment (Y-intact healthy males)
#          results/tables/loy/02_calls/loy_reference_zero_rates.csv  empirical zero-Y rate per chemistry x compartment x depth bin
#          results/tables/loy/02_calls/loy_sample_sex.csv            expression sex per sample, with karyotype sex where known
#          results/tables/loy/02_calls/loy_compartment_calls.csv     per sample x compartment: n, suspected, retained, uncertain, Z, E0, f_hat
#          results/tables/loy/02_calls/loy_patient_table.csv         one row per male sample: karyotype, drivers, LEUK class, pairing flag
#          results/tables/loy/02_calls/loy_call_checks.csv           known answers K1-K7
#          results/tables/loy/02_calls/loy_fine_label_map.csv        which bmm_fine labels became T / NK (review this)
# what it does: LOY step 1. Applies the frozen three-state rule of config_loy.R section 5 to every scanned cell,
#   estimates the per-sample LOY fraction per compartment against a background built ONLY from expression-sexed
#   healthy male donors of the same chemistry, and writes the paired-case table that decides whether the H2 design
#   (same patient, same state, -Y vs Y+ leukaemic cells) has cases.
# Usage  : Rscript scripts/loy/02_call_loy_status.R            (requires all datasets of 01_scan/datasets.txt)
#          Rscript scripts/loy/02_call_loy_status.R --partial  (whatever percell files exist)
# WHY the reference is external: estimating the zero background from the sample being called would absorb its own
#   LOY into the background and hide it; the same lesson as the inferCNV reference choice in 02_malignancy.
# Revision 2026-11-14 (external cohorts, config section 7): per-cell files named in external_datasets.txt are loaded
#   when present; their chemistry is chem_lit from LOY_EXT_STUDIES; they never enter the reference pool; lambda is
#   scaled by a per-study factor from the study's own expression-male T cells (7b) and E0 = max(exp(-lambda), e0_emp).
#   On-disk datasets are called exactly as before, so K1-K7 reproduce. New: K8 (vanGalen AML707B via the pack),
#   loy_external_calibration.csv.
# Revision 2026-11-14 (after run 1): (a) the objects' Chemistry column reads "unknown" for 13 of 16 datasets, so
#   chemistry now comes from manifest platform_fine, then a per-dataset default; (b) sex calls use XIST >= 2 UMIs on
#   cells with >= 2,000 UMIs because one XIST UMI is what a female pool-mate's ambient RNA delivers (GSE185381 pools);
#   (c) f_hat is NA when a compartment has < 50 cells or median lambda < 2 (AML1371's 23 LEUK cells are not readable).

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "loy", "config_loy.R"))
set.seed(SEED)
args <- commandArgs(TRUE); partial <- "--partial" %in% args
dir.create(LOY_CALL_DIR, showWarnings = FALSE, recursive = TRUE)
fmt <- function(x) if (length(x) && !is.na(x[1])) as.character(round(x[1], 3)) else "NA"

## [1] load ---------------------------------------------------------------------------------------
dsl   <- readLines(file.path(LOY_SCAN_DIR, "datasets.txt"))
files <- file.path(LOY_PERCELL_DIR, paste0(dsl, "__y_percell.csv.gz"))
have  <- file.exists(files)
if (!partial) stopifnot(all(have)) else files <- files[have]
message(sprintf("[1] per-cell files: %d of %d datasets%s", sum(have), length(dsl), if (partial) " (partial)" else ""))
ext_dsl   <- if (file.exists(LOY_EXT_LIST)) readLines(LOY_EXT_LIST) else character()      # config section 7
ext_files <- file.path(LOY_PERCELL_DIR, paste0(ext_dsl, "__y_percell.csv.gz")); ext_files <- ext_files[file.exists(ext_files)]
message(sprintf("    external per-cell files present: %d of %d", length(ext_files), length(ext_dsl)))
cells <- rbindlist(lapply(c(files, ext_files), fread), fill = TRUE)
cells[, is_ext := dataset %in% ext_dsl]
ycols <- setdiff(grep("^y_[A-Z0-9]+$", names(cells), value = TRUE), "y_umi")
cells[, n_ygenes := rowSums(as.matrix(.SD) > 0, na.rm = TRUE), .SDcols = ycols]
cells[, `:=`(nCount_RNA = as.numeric(nCount_RNA), y_umi = as.numeric(y_umi), xist_umi = as.numeric(xist_umi))]
cells[, Sample := fifelse(is.na(Sample) | Sample == "", object, as.character(Sample))]
message(sprintf("    cells %d | samples %d | datasets %d", nrow(cells), uniqueN(cells[, .(dataset, Sample)]), uniqueN(cells$dataset)))

## [2] chemistry ------------------------------------------------------------------------------------
man <- fread(file.path(DIR_PREPROCESS, "00_curated_manifest.csv"))
chem_of <- function(pf) fifelse(grepl("5p", pf), "5prime", fifelse(grepl("3p", pf), "3prime", fifelse(grepl("Seq", pf, ignore.case = TRUE), "SeqWell", NA_character_)))
man[, chem_man := chem_of(platform_fine)]
cells <- merge(cells, man[, .(dataset, Sample = sample, chem_man, platform_fine, karyotype, driver_mutations, fusion,
                               disease_man = disease, timepoint_man = timepoint, sorting_r, sample_state, blast_pct_clinical)],
               by = c("dataset", "Sample"), all.x = TRUE)
known_chem <- c("3prime", "5prime", "SeqWell")
cells[, chem := fifelse(as.character(Chemistry) %in% known_chem, as.character(Chemistry), chem_man)]
ds_default <- c(Chen2023 = "3prime", "E-MTAB-11536" = "5prime", GSE116256 = "SeqWell", GSE147989 = "3prime",
                GSE185381 = "3prime", GSE185991 = "3prime", GSE201966 = "3prime", GSE207356 = "3prime",
                GSE227903 = "3prime", GSE239721 = "5prime", GSE253355 = "3prime", GSE289435 = "5prime",
                Petti2019 = "5prime", GSE120221 = "3prime", GSE227690 = "3prime", GSE276360 = "3prime")
cells[is.na(chem) & dataset %in% names(ds_default), chem := ds_default[dataset]]
cells[is_ext == TRUE, chem := LOY_EXT_ALL$chem_lit[match(dataset, LOY_EXT_ALL$dataset)]]   # literature guess, calibrated in [6b]
cells[is.na(chem), chem := "unknown"]
message("[2] chemistry: ", paste(cells[, .N, by = .(dataset, chem)][order(dataset)][, paste0(dataset, "=", chem)], collapse = " "))

## [3] compartments ----------------------------------------------------------------------------------
cells[, bin := fifelse(is.na(hierarchy_bin) | hierarchy_bin == "", "NA", as.character(hierarchy_bin))]
cells[, fine := fifelse(is.na(bmm_fine), "", as.character(bmm_fine))]
cells[, comp := bin]
cells[bin %in% LOY_LEUK_BINS, comp := "LEUK"]
cells[bin == "T_NK", comp := fifelse(grepl(LOY_T_FINE_REGEX, fine) & !grepl(LOY_T_FINE_EXCLUDE, fine), "T",
                              fifelse(grepl(LOY_NK_FINE_REGEX, fine), "NK", "T_NK_other"))]
fwrite(cells[bin == "T_NK", .N, by = .(fine, comp)][order(comp, -N)], file.path(LOY_CALL_DIR, "loy_fine_label_map.csv"))
cells[, depth_bin := cut(nCount_RNA, LOY_DEPTH_BREAKS, right = FALSE)]

## [4] expression sex per sample ----------------------------------------------------------------------
sex <- cells[, .(n = .N, n_deep = sum(nCount_RNA >= 2000),
                 xist_pos_frac = mean(xist_umi[nCount_RNA >= 2000] >= 2, na.rm = TRUE),   # >= 2 UMIs: above pool-mate ambient
                 xist_any_frac = mean(xist_umi > 0, na.rm = TRUE),
                 y_rate_all = sum(y_umi) / sum(nCount_RNA), med_ncount = as.numeric(median(nCount_RNA)),
                 karyotype = karyotype[1], disease_man = disease_man[1], Disease_state = as.character(Disease_state[1]),
                 chem = chem[1], object = object[1]), by = .(dataset, Sample)]
# Two readouts, both needed (revision 2 of 2026-11-14): XIST alone mis-called two shallow 5' females and one
# Oetjen female as male (XIST >= 2 UMIs in < 5% of cells but NO Y at all); Y alone cannot tell a female from a
# male whose every cell lost Y. r_c = crude Y rate per chemistry over all samples (reference-free, used only here).
r_chem_all <- cells[, .(r_c = sum(y_umi) / sum(nCount_RNA)), by = chem]
sex <- merge(sex, r_chem_all, by = "chem", all.x = TRUE)
# Y presence is judged on the whole sample AND on the T cells alone: in a clonal -Y AML nearly every cell is a
# -Y clone member (GSM8791436: 99.8%), so the whole-sample Y rate reads "absent" although the patient is male.
t_rate <- cells[comp == "T", .(n_T_cells = .N, y_rate_T = sum(y_umi) / sum(nCount_RNA)), by = .(dataset, Sample)]
sex <- merge(sex, t_rate, by = c("dataset", "Sample"), all.x = TRUE)
sex[, y_present := y_rate_all >= 0.25 * r_c | (!is.na(n_T_cells) & n_T_cells >= 30 & y_rate_T >= 0.25 * r_c)]
sex[, y_absent  := y_rate_all < 0.10 * r_c & !(!is.na(n_T_cells) & n_T_cells >= 30 & y_rate_T >= 0.25 * r_c)]
# A female sample with hundreds of deep cells always shows XIST >= 2 UMIs in some of them (lowest seen: 0.07);
# XIST at exactly ~0 with no Y therefore means a male whose cells have all lost Y, not a female.
sex[, sex_expr := fifelse(n_deep < 30, "too_shallow",
                  fifelse(is.nan(xist_pos_frac), fifelse(y_present, "M_noXISTgene", fifelse(y_absent, "F_noXISTgene", "mixed")),
                  fifelse(xist_pos_frac >= LOY_XIST_FEMALE & !y_present, "F",
                  fifelse(xist_pos_frac <= LOY_XIST_MALE & y_present, "M",
                  fifelse(y_absent & xist_pos_frac <= 0.005 & n_deep >= 500, "M_LOY_all",       # male, Y lost in nearly every cell
                  fifelse(y_absent & xist_pos_frac <= LOY_XIST_MALE, "noY_noXIST",             # too few deep cells to decide
                  fifelse(y_absent, "F_lowdepth",                                               # XIST under-detected, Y absent
                  fifelse(y_rate_all >= 0.50 * r_c, "M_xist_ambient", "mixed"))))))))]        # Y at the male rate, XIST from pool-mates
sex[, sex_karyo := fifelse(grepl("XX", gsub(" ", "", karyotype)), "F", fifelse(grepl("XY|,-Y|X,-Y", gsub(" ", "", karyotype)), "M", NA_character_))]
# declared exclusions (config section 6): the cells are kept in the per-cell files, the sample leaves every call table
sex <- merge(sex, LOY_SAMPLE_EXCLUDE, by = c("dataset", "Sample"), all.x = TRUE)
sex[!is.na(reason), sex_expr := "excluded"]
sex[, sex_group := fifelse(sex_expr %in% c("M", "M_xist_ambient", "M_noXISTgene", "M_LOY_all"), "M",
                   fifelse(sex_expr %in% c("F", "F_lowdepth", "F_noXISTgene"), "F", "other"))]
s3 <- tryCatch(fread(file.path(TAB_DIR, "ws", "ws_s3_samples.csv"))[, .(Sample = sample, sex_s3 = sex, age_group, paediatric)], error = function(e) NULL)
if (!is.null(s3)) sex <- merge(sex, s3, by = "Sample", all.x = TRUE)
fwrite(sex, file.path(LOY_CALL_DIR, "loy_sample_sex.csv"))
review <- sex[(sex_group == "other") | (!is.na(sex_karyo) & sex_group != "other" & sex_group != sex_karyo) | !is.na(reason),
              .(dataset, Sample, sex_expr, sex_karyo, n_deep, xist_pos_frac = round(xist_pos_frac, 3), y_rate_all = signif(y_rate_all, 3),
                n_T_cells, y_rate_T = signif(y_rate_T, 3), reason)]
fwrite(review[order(dataset, Sample)], file.path(LOY_CALL_DIR, "loy_sex_review.csv"))
message("[4] samples for manual sex review: ", nrow(review))
cells <- merge(cells, sex[, .(dataset, Sample, sex_expr, sex_group)], by = c("dataset", "Sample"))
message("[4] sex calls: ", paste(names(table(sex$sex_expr)), table(sex$sex_expr), collapse = ", "))

## [5] reference background from clean healthy MALE donors ----------------------------------------------
# reference donors: the project's healthy-control datasets (E-MTAB-11536 and GSE253355 are coded non_heme_condition in the
# manifest) plus the healthy donors inside AML studies; clean "M" calls only (no ambient-XIST or no-XIST-gene calls)
is_ref_sample <- sex[sex_expr == "M" & !(dataset %in% ext_dsl) &     # external healthy donors never enter the reference (7b)
                     (Disease_state %in% LOY_REF_DISEASE | disease_man %in% c("healthy_donor", "Healthy", "healthy", "non_heme_condition")), .(dataset, Sample)]
ref <- merge(cells, is_ref_sample, by = c("dataset", "Sample"))
message(sprintf("[5] reference: %d male healthy samples, %d cells; chemistry x cells: %s", uniqueN(ref[, .(dataset, Sample)]), nrow(ref),
                paste(ref[, .N, by = chem][, paste0(chem, "=", N)], collapse = ", ")))
r_ref  <- ref[, .(r = sum(y_umi) / sum(nCount_RNA), n_ref = .N), by = .(chem, comp)]
r_chem <- ref[, .(r_chem = sum(y_umi) / sum(nCount_RNA)), by = chem]
e0_emp <- ref[, .(e0_emp = mean(y_umi == 0), n_bin = .N), by = .(chem, comp, depth_bin)]
fwrite(merge(r_ref, r_chem, by = "chem")[order(chem, -n_ref)], file.path(LOY_CALL_DIR, "loy_reference_rates.csv"))
fwrite(e0_emp[order(chem, comp, depth_bin)], file.path(LOY_CALL_DIR, "loy_reference_zero_rates.csv"))

## [6] per-cell three-state call -------------------------------------------------------------------
cells <- merge(cells, r_ref[, .(chem, comp, r)], by = c("chem", "comp"), all.x = TRUE)
cells <- merge(cells, r_chem, by = "chem", all.x = TRUE)
cells[is.na(r), r := r_chem]
cells[is.na(r), r := median(r_ref$r)]                      # chemistry with no healthy male at all: pooled fallback, flagged below
## [6b] per-study calibration for external datasets (config 7b): the study's own expression-male T cells set the scale
r_T_assigned <- cells[comp == "T", .(r_T_ref = r[1]), by = dataset]
calib <- cells[is_ext == TRUE & comp == "T" & sex_group == "M", .(n_T_male = .N, y_rate_T = sum(y_umi) / sum(nCount_RNA),
                                                                  n_male_samples = uniqueN(Sample)), by = dataset]
calib <- merge(calib, r_T_assigned, by = "dataset", all.x = TRUE)
calib[, factor := fifelse(n_T_male >= LOY_EXT_CALIB_MIN_T, y_rate_T / r_T_ref, 1)]
calib[, calib_status := fifelse(n_T_male < LOY_EXT_CALIB_MIN_T, "uncalibrated",
                        fifelse(factor < LOY_EXT_CALIB_RANGE[1] | factor > LOY_EXT_CALIB_RANGE[2], "chem_mismatch", "calibrated"))]
calib <- merge(data.table(dataset = unique(cells[is_ext == TRUE, dataset])), calib, by = "dataset", all.x = TRUE)
calib[is.na(factor), `:=`(factor = 1, calib_status = "uncalibrated")]
calib <- merge(calib, LOY_EXT_ALL[, .(dataset, chem_lit, on_disk)], by = "dataset", all.x = TRUE)
if (nrow(calib)) { fwrite(calib[order(dataset)], file.path(LOY_CALL_DIR, "loy_external_calibration.csv")); message("[6b] external calibration:"); print(calib[, .(dataset, chem_lit, n_T_male, y_rate_T = signif(y_rate_T, 3), r_T_ref = signif(r_T_ref, 3), factor = round(factor, 2), calib_status)]) }
cells <- merge(cells, calib[, .(dataset, factor, calib_status)], by = "dataset", all.x = TRUE)
cells[is.na(factor), `:=`(factor = 1, calib_status = "on_disk")]
cells[, lambda := nCount_RNA * r * factor]
cells[, call := fifelse(y_umi == 0 & lambda >= LOY_LAMBDA_MIN, "suspected_LOY",
               fifelse((y_umi >= 2 & n_ygenes >= 2) | y_umi >= 3, "Y_retained", "uncertain"))]
cells <- merge(cells, e0_emp[, .(chem, comp, depth_bin, e0_emp, n_bin)], by = c("chem", "comp", "depth_bin"), all.x = TRUE)
cells[, e0 := fifelse(!is.na(e0_emp) & n_bin >= 30, e0_emp, exp(-lambda))]
cells[is_ext == TRUE, e0 := pmax(e0, exp(-lambda))]        # 7b: conservative whichever way the chemistry guess is wrong
no_ref_chem <- setdiff(unique(cells$chem), r_chem$chem)
if (length(no_ref_chem)) message("[6] WARNING no healthy-male reference for chemistry: ", paste(no_ref_chem, collapse = ","), " -> pooled fallback rate")

## [7] per sample x compartment ---------------------------------------------------------------------
comp_tab <- cells[sex_group == "M", .(n = .N, n_susp = sum(call == "suspected_LOY"), n_ret = sum(call == "Y_retained"),
                                       n_unc = sum(call == "uncertain"), Z = mean(y_umi == 0), E0 = mean(e0),
                                       med_ncount = as.numeric(median(nCount_RNA)), med_lambda = as.numeric(median(lambda))),
                  by = .(dataset, Sample, comp)]
comp_tab[, f_hat := pmax(0, pmin(1, (Z - E0) / (1 - E0)))]
comp_tab[, readable := n >= LOY_MIN_CELLS & med_lambda >= 2]
comp_tab[readable == FALSE, f_hat := NA_real_]
comp_tab <- merge(comp_tab, unique(cells[, .(dataset, calib_status)]), by = "dataset", all.x = TRUE)
comp_tab[calib_status == "chem_mismatch", f_hat := NA_real_]     # 7b: report, do not call
# Pool-contamination guard (added 2026-11-14 after AML073 / AML051): in a MALE sample, Y-zero cells that express XIST are
# female cells carrying another donor's label (GSE185381 pool P08: 930 of AML073's deep Mono_DC cells, 89 % XIST+; P19:
# AML051's Mono_DC). Such a compartment's LOY fraction is not the patient's and is set to NA with the flag kept.
xchk <- cells[sex_group == "M" & nCount_RNA >= 3000, .(n_y0_deep = sum(y_umi == 0), xist_in_y0 = mean(xist_umi[y_umi == 0] >= 1),
                                                       n_ypos_deep = sum(y_umi > 0), xist_in_ypos = mean(xist_umi[y_umi > 0] >= 1)),
              by = .(dataset, Sample, comp)]
comp_tab <- merge(comp_tab, xchk, by = c("dataset", "Sample", "comp"), all.x = TRUE)
comp_tab[, female_cells_in_label := (n_y0_deep >= 20 & xist_in_y0 >= 0.30) %in% TRUE]
comp_tab[female_cells_in_label == TRUE, f_hat := NA_real_]
fwrite(comp_tab[female_cells_in_label == TRUE | (n_y0_deep >= 20 & xist_in_y0 >= 0.15) %in% TRUE,
                .(dataset, Sample, comp, n, n_y0_deep, xist_in_y0 = round(xist_in_y0, 3), n_ypos_deep, xist_in_ypos = round(xist_in_ypos, 3), female_cells_in_label)],
       file.path(LOY_CALL_DIR, "loy_pool_contamination.csv"))
message("[7] compartments with female cells under a male label: ", sum(comp_tab$female_cells_in_label))
comp_tab <- merge(comp_tab, sex[, .(dataset, Sample, sex_expr, chem, Disease_state, disease_man, karyotype)], by = c("dataset", "Sample"))
fwrite(comp_tab[order(dataset, Sample, comp)], file.path(LOY_CALL_DIR, "loy_compartment_calls.csv"))

## [8] patient table -------------------------------------------------------------------------------
wide <- dcast(comp_tab[comp %in% c("LEUK", "Mono_DC", "T", "NK", "B_Plasma", "Erythroid")],
              dataset + Sample ~ comp, value.var = c("n", "n_susp", "n_ret", "f_hat"), fill = NA)
pt <- merge(sex[sex_group == "M"], wide, by = c("dataset", "Sample"), all.x = TRUE)
pt <- merge(pt, cells[, .(timepoint_man = timepoint_man[1], driver_mutations = driver_mutations[1], fusion = fusion[1], platform_fine = platform_fine[1],
                          sorting_r = sorting_r[1], sample_state = sample_state[1], blast_pct_clinical = blast_pct_clinical[1],
                          Timepoint = as.character(Timepoint[1]), AgeGroup = if ("AgeGroup" %in% names(cells)) as.character(AgeGroup[1]) else NA_character_,
                          is_ext = is_ext[1], calib_status = calib_status[1], factor = factor[1]), by = .(dataset, Sample)],
            by = c("dataset", "Sample"), all.x = TRUE)
pt[is.na(timepoint_man) | timepoint_man == "", timepoint_man := Timepoint]
pt[, is_aml := disease_man == "AML" | Disease_state %in% LOY_AML_DISEASE]   # config 7d
pt[, karyo_minusY := grepl("-\\s*Y|−Y", karyotype)]
pt[, leuk_class := fifelse(!is_aml %in% TRUE, NA_character_,
                   fifelse(is.na(f_hat_LEUK), "not_readable",
                   fifelse(f_hat_LEUK >= 0.80, "clonal_LOY",
                   fifelse(f_hat_LEUK >= 0.10 & n_susp_LEUK >= LOY_MIN_CELLS & n_ret_LEUK >= LOY_MIN_CELLS, "subclonal_LOY",
                   fifelse(f_hat_LEUK >= 0.03 & n_susp_LEUK >= 20, "minor_LOY", "none")))))]
pt[, monoDC_class := fifelse(is.na(f_hat_Mono_DC), "not_readable", fifelse(f_hat_Mono_DC >= 0.5, "mostly_LOY", fifelse(f_hat_Mono_DC >= 0.1, "mixed", "mostly_retained")))]
# Monocytic AML lives in the Mono_DC bin (truth audit 2026-10-15: Mono/DC in AML samples are mostly clonal), so a -Y
# clone can sit there while the LEUK bin holds Y-intact cells. A sample is an H2 pairing candidate when ANY of the two
# leukaemia-eligible compartments holds >= LOY_MIN_CELLS suspected-LOY AND >= LOY_MIN_CELLS Y-retained cells, or when the
# two compartments disagree (one >= 0.5, the other <= 0.2) -- the latter needs a shared-driver anchor before it counts.
pt[, pair_same_comp := (n_susp_LEUK >= LOY_MIN_CELLS & n_ret_LEUK >= LOY_MIN_CELLS & f_hat_LEUK >= 0.10 & f_hat_LEUK <= 0.90) %in% TRUE |
                       (n_susp_Mono_DC >= LOY_MIN_CELLS & n_ret_Mono_DC >= LOY_MIN_CELLS & f_hat_Mono_DC >= 0.10 & f_hat_Mono_DC <= 0.90) %in% TRUE]
pt[, pair_cross_comp := (!is.na(f_hat_LEUK) & !is.na(f_hat_Mono_DC) & abs(f_hat_LEUK - f_hat_Mono_DC) >= 0.5 &
                         pmin(n_LEUK, n_Mono_DC) >= LOY_MIN_CELLS) %in% TRUE]
pt[, h2_pairable := is_aml %in% TRUE & (pair_same_comp | pair_cross_comp)]
pt[, loy_any_leuk_comp := fifelse(leuk_class %in% "clonal_LOY" | (f_hat_Mono_DC >= 0.80) %in% TRUE, "clonal",
                           fifelse((f_hat_LEUK >= 0.10) %in% TRUE | (f_hat_Mono_DC >= 0.10) %in% TRUE, "partial",
                           fifelse((f_hat_LEUK >= 0.03) %in% TRUE | (f_hat_Mono_DC >= 0.03) %in% TRUE, "minor", "none")))]
setorder(pt, -is_aml, dataset, Sample)
fwrite(pt, file.path(LOY_CALL_DIR, "loy_patient_table.csv"))

## [9] known answers -------------------------------------------------------------------------------
fem_zero <- cells[sex_group == "F" & nCount_RNA >= 2000, mean(y_umi == 0)]
g36   <- comp_tab[Sample == "GSM8791436_MLL_17844" & comp == "LEUK", f_hat]
a1371 <- comp_tab[Sample == "AML1371" & comp == "LEUK", .(f_hat, n)]
a707  <- comp_tab[Sample == "AML707B-D0" & comp == "LEUK", .(f_hat, n)]
hm_T  <- merge(comp_tab[comp == "T" & n >= 100], is_ref_sample, by = c("dataset", "Sample"))$f_hat
z707  <- comp_tab[dataset == "ZENG_vanGalen_2019" & Sample == "AML_pvg2019_AML707B" & comp == "LEUK", .(f_hat, n)]   # same cells via the pack
z36   <- comp_tab[dataset == "ZENG_Zeng_2025" & Sample == "scAML_AEL_MLL_17844" & comp == "LEUK", .(f_hat, n)]          # = GSM8791436 (K2) via the pack
z025  <- comp_tab[dataset == "ZENG_Lasry_2022" & Sample == "AML_Lasry2022_pediatric_AML025" & comp == "LEUK", .(f_hat, n)]   # = AML025 clonal case via the pack
c4    <- sex[Sample == "Control4", sex_expr]
mism  <- sex[!is.na(sex_karyo) & sex_group != "other" & sex_group != sex_karyo & Sample != "Control4" & is.na(reason),
             .(dataset, Sample, sex_expr, sex_karyo, xist_pos_frac = round(xist_pos_frac, 3), y_rate_all = signif(y_rate_all, 3))]
chk <- data.table(
  check    = c("K1_female_zeroY_frac_deep_cells", "K2_GSM8791436_LEUK_fhat", "K3_AML1371_LEUK_fhat(3of20)",
               "K4_healthy_male_T_fhat_median", "K5_Control4_expression_sex", "K6_sex_vs_karyotype_mismatches(after declared exclusions)", "K7_AML707B-D0_LEUK_fhat(19of20 metaphases lack Y)",
               "K8_ZENG_vanGalen_AML707B_LEUK_fhat(pack copy of K7)", "K9_ZENG_Zeng2025_MLL_17844_LEUK_fhat(pack copy of K2)",
               "K10_ZENG_Lasry_AML025_LEUK_fhat(pack copy of the clonal case)"),
  expected = c(">= 0.95", ">= 0.90", "0.05-0.40 if readable", "<= 0.03", "F", "0", ">= 0.85 if readable", ">= 0.85 if scanned", ">= 0.90 if scanned", ">= 0.80 if scanned"),
  observed = c(fmt(fem_zero), fmt(g36), paste0(fmt(a1371$f_hat), " n=", if (nrow(a1371)) a1371$n else 0),
               fmt(if (length(hm_T)) median(hm_T, na.rm = TRUE) else NA), if (length(c4)) c4 else "NA", as.character(nrow(mism)),
               paste0(fmt(a707$f_hat), " n=", if (nrow(a707)) a707$n else 0),
               paste0(fmt(z707$f_hat), " n=", if (nrow(z707)) z707$n else 0),
               paste0(fmt(z36$f_hat), " n=", if (nrow(z36)) z36$n else 0), paste0(fmt(z025$f_hat), " n=", if (nrow(z025)) z025$n else 0)),
  pass     = c(isTRUE(fem_zero >= 0.95), isTRUE(length(g36) > 0 && g36 >= 0.90),
               if (nrow(a1371) && !is.na(a1371$f_hat)) isTRUE(a1371$f_hat >= 0.05 && a1371$f_hat <= 0.40) else NA,
               isTRUE(length(hm_T) > 0 && median(hm_T, na.rm = TRUE) <= 0.03), isTRUE(length(c4) > 0 && c4 %in% c("F", "F_lowdepth")),
               nrow(mism) == 0, if (nrow(a707) && !is.na(a707$f_hat)) isTRUE(a707$f_hat >= 0.85) else NA,
               if (nrow(z707) && !is.na(z707$f_hat)) isTRUE(z707$f_hat >= 0.85) else NA,
               if (nrow(z36) && !is.na(z36$f_hat)) isTRUE(z36$f_hat >= 0.90) else NA, if (nrow(z025) && !is.na(z025$f_hat)) isTRUE(z025$f_hat >= 0.80) else NA))
fwrite(chk, file.path(LOY_CALL_DIR, "loy_call_checks.csv")); print(chk)
if (nrow(mism)) { message("[K6] sex mismatches (metadata vs expression):"); print(mism) }
message("[done] male samples: ", nrow(pt), " | male AML: ", sum(pt$is_aml %in% TRUE),
        " | LEUK classes: ", paste(names(table(pt$leuk_class)), table(pt$leuk_class), collapse = ", "),
        " | any leukaemic compartment: ", paste(names(table(pt[is_aml %in% TRUE, loy_any_leuk_comp])), table(pt[is_aml %in% TRUE, loy_any_leuk_comp]), collapse = ", "),
        " | H2 pairing candidates: ", sum(pt$h2_pairable))
if (any(pt$is_ext %in% TRUE)) message("[done] external: datasets ", uniqueN(pt[is_ext %in% TRUE, dataset]), " | male samples ", sum(pt$is_ext %in% TRUE),
        " | male AML ", sum(pt$is_ext %in% TRUE & pt$is_aml %in% TRUE),
        " | LEUK classes: ", paste(names(table(pt[is_ext %in% TRUE, leuk_class])), table(pt[is_ext %in% TRUE, leuk_class]), collapse = ", "),
        " | H2 pairing candidates: ", sum(pt[is_ext %in% TRUE, h2_pairable]))
print(pt[h2_pairable == TRUE, .(dataset, Sample, timepoint_man, n_LEUK, f_hat_LEUK = round(f_hat_LEUK, 3), n_Mono_DC, f_hat_Mono_DC = round(f_hat_Mono_DC, 3), n_T, f_hat_T = round(f_hat_T, 3), pair_same_comp, pair_cross_comp)])
