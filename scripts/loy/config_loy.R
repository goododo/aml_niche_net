# config_loy.R ----
# INPUT  : none (sources scripts/config/config_paths.R via here)
# OUTPUT : constants only -- everything the LOY (loss of chromosome Y) line reads
# what it does: step 0 of the LOY line (user approval 2026-11-14). Y-chromosome loss as a per-cell
#   genomic label: a male cell that lost Y has zero UMIs on every chrY-linked gene at once. These
#   constants fix the gene set and the output locations BEFORE any call is made. The three-state call
#   (suspected LOY / Y retained / uncertain) and its thresholds live in 02 and are frozen there, after
#   the scan has shown the depth-dependent zero background; nothing in this file decides a call.
# Business scripts start with:
#   source(here::here("scripts", "config", "config_paths.R"))
#   source(here::here("scripts", "loy",    "config_loy.R"))

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))   # SEED, FAST_DIR, LARGE1_DIR, TAB_DIR, ...

LOY_CONFIG_DATE <- "2026-11-14"

## -- 1. gene sets ------------------------------------------------------------------------------
# chrY non-PAR genes with broad expression in blood cells (X homologues in brackets, for the record):
# RPS4Y1 (RPS4X), DDX3Y (DDX3X), EIF1AY (EIF1AX), KDM5D (KDM5C), UTY (KDM6A), USP9Y (USP9X), ZFY (ZFX),
# TXLNGY, NLGN4Y (NLGN4X), PRKY (PRKX), TMSB4Y (TMSB4X), LINC00278, TTTY15. PAR genes are excluded on
# purpose: they are on both X and Y and do not drop to zero with Y loss. Any member missing from a
# dataset's gene axis is skipped and the number found is recorded per dataset.
LOY_YGENES <- c("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "ZFY",
                "TXLNGY", "NLGN4Y", "PRKY", "TMSB4Y", "LINC00278", "TTTY15")
LOY_XIST   <- "XIST"          # female cells express it; male cells, with or without Y, do not

## -- 2. input / output -------------------------------------------------------------------------
LOY_ANNOT_DIR   <- file.path(LARGE1_DIR, "02_seurat_objects/04_annotated")   # per-sample objects, 16 datasets
LOY_PERCELL_DIR <- file.path(LARGE1_DIR, "loy/percell")                        # <ds>__y_percell.csv.gz (large)
LOY_TAB_DIR     <- file.path(TAB_DIR, "loy")
LOY_SCAN_DIR    <- file.path(LOY_TAB_DIR, "01_scan")                           # <ds>__y_sample_summary.csv
LOY_ZENG_RDS_DIR <- file.path(ZENODO_DIR, "zenodo.18528442")                  # Zeng 2025 21-cohort pack

## -- 3. metadata columns carried from the annotated objects (missing ones become NA) ------------
LOY_META_COLS <- c("Sample", "Patient_ID", "Library", "Chemistry", "Platform", "Timepoint",
                   "Disease_state", "hierarchy_bin", "bmm_fine", "high_error", "anno_malignant",
                   "nCount_RNA", "nFeature_RNA", "percent.mt", "anno_sex")

## -- 4. depth bins for the zero-Y background (reporting only; calls are made in 02) ---------------
LOY_DEPTH_BREAKS <- c(0, 1000, 2000, 3000, 5000, 8000, Inf)

## -- 5. three-state call, FROZEN 2026-11-14 before any AML sample was read ------------------------
# Written after the smoke run on GSE276360 (healthy donors only): the male donor's zero-Y fraction fell
# from 0.33 (< 1,500 UMI) to 0.15 (1,500-3,000), 0.023 (3,000-5,000) and 0.003 (>= 5,000), i.e. Poisson
# in depth with ~1.2 Y UMIs per 1,000 UMIs, slightly over-dispersed. Rules:
#   lambda_i      = nCount_i * r_ref(chemistry, compartment); r_ref is pooled from Y-INTACT reference cells
#                   (healthy male donors; expression-sexed) -- never from the sample being called.
#   suspected_LOY : y_umi == 0 AND lambda_i >= LOY_LAMBDA_MIN  (Poisson P(0 | intact) <= 0.01)
#   Y_retained    : y_umi >= 2 AND >= 2 distinct Y genes detected, OR y_umi >= 3   (one stray UMI can be ambient)
#   uncertain     : everything else (zero at low depth, or a single Y UMI)
# Sample-level LOY fraction (all cells of a compartment, not only the called ones):
#   f_hat = (Z - E0) / (1 - E0), Z = observed zero fraction, E0 = mean over cells of the EMPIRICAL zero rate of
#   reference cells in the same chemistry x compartment x depth bin (falls back to exp(-lambda) when the bin is
#   empty). Ambient Y UMIs on true -Y cells lower Z and so bias f_hat DOWN: conservative.
# Sex of a sample from expression, not metadata: female if XIST-positive fraction >= LOY_XIST_FEMALE, male if
#   <= LOY_XIST_MALE, otherwise "mixed" (pools, transplant chimerism, mis-assembled objects) and excluded.
# A sample is classified on its LEUK compartment (HSC_MPP + LMPP_GMP in AML samples):
#   clonal_LOY     f_hat >= 0.80 ; subclonal_LOY 0.10 <= f_hat < 0.80 with >= LOY_MIN_CELLS suspected AND
#   >= LOY_MIN_CELLS retained ; minor_LOY 0.03 <= f_hat < 0.10 with >= 20 suspected ; none otherwise.
LOY_LAMBDA_MIN      <- 4.605   # -log(0.01)
LOY_MIN_CELLS       <- 50L
LOY_XIST_FEMALE     <- 0.30
LOY_XIST_MALE       <- 0.05
LOY_LEUK_BINS       <- c("HSC_MPP", "LMPP_GMP")
LOY_REF_DISEASE     <- c("Healthy", "Control", "healthy", "control", "Normal")   # Disease_state values of reference donors
LOY_T_FINE_REGEX    <- "^(CD4|CD8|T |T-|Treg|gd|gamma|MAIT|NKT)"      # strict T from bmm_fine; never the T_NK bin alone
LOY_T_FINE_EXCLUDE  <- "Lymphoid|CLP|Pre-?pro|Pro-B|Pre-B|Cycling|Progenitor"
LOY_NK_FINE_REGEX   <- "^NK"
LOY_CALL_DIR        <- file.path(LOY_TAB_DIR, "02_calls")

## -- 6. declared sample-level exclusions (2026-11-14, after run 3 of 02) ----------------------------
# Expression sex is read from the cells, never forced from metadata. Where expression and karyotype disagree the
# sample is listed here with the reason and left out of every LOY table until the identity is resolved.
#   PAUMTZ / PAXMIJ (GSE185381 pool P06, demuxed with Control0058): the cells labelled PAUMTZ have XIST in 17-58% of
#   deep cells per compartment and no Y at all; the cells labelled PAXMIJ carry Y and no XIST. PAUMTZ's karyotype is
#   46,XY and PAXMIJ's is 46,XX,t(8;21)/45,idem,-X. Each expression profile matches the OTHER donor's karyotype:
#   a label swap in the depositors' demultiplexing. Not corrected here (drivers and karyotype would have to swap too).
LOY_SAMPLE_EXCLUDE <- data.table::data.table(
  dataset = c("GSE185381", "GSE185381"),
  Sample  = c("PAUMTZ", "PAXMIJ"),
  reason  = c("suspected demux label swap in pool P06 (expression female, karyotype male)",
              "suspected demux label swap in pool P06 (expression male, karyotype female)"))

## -- 7. external cohorts for the case search (user approval 2026-11-14 "批准，继续": Zeng pack + GSE269669) -------------
# The on-disk 16 datasets gave 3 clonal -Y AML and 0 subclonal cases (report 2026-11-14). The search is extended to
# the Zeng et al. 2025 pack (Zenodo 18528442: 21 studies, 22 RDS, BoneMarrowMap labels in the SAME vocabulary as
# bmm_bin_map.tsv) and to GSE269669 (11 TP53-mutant AML, aza/ven/magrolimab, author Seurat objects).
# Chemistry below is from the source papers as remembered, NOT verified against the pack (it has no chemistry
# column); it only picks which reference zero-background is used, and the per-study calibration in 02 (section 7b)
# corrects the Y-capture rate from each study's own male T cells, so a wrong guess shows up as a factor far from 1.
# "on_disk" names the project dataset holding the same cells: those studies are reproduction checks, not new cases.
LOY_EXT_STUDIES <- data.table::fread(text = "
study,file_stem,chem_lit,on_disk,note
Abbas_2021,scAML_Abbas_2021,5prime,,10x 5' + TCR; AML on PD-1 blockade
Bailur_2020,scAML_Bailur_2020,5prime,,paediatric marrow T cells; 10x 5' (guess)
BeneytoCalabuig_2023,scAML_BeneytoCalabuig_2023,3prime,,CloneTracer; 10x 3' v3
Cobos_2023,scAML_Cobos_2023,unknown,,platform not recalled
Dufva_2020,scAML_Dufva_2020,3prime,,immunogenomic landscape; 10x 3' v2
Ennis_2023,scAML_Ennis_2023,3prime,GSE227903,on disk as GSE227903 (1216/1886/3904/40389/5750/6323 ids); 10x 3'
Gawad_ScPCA,scAML_Gawad_ScPCA,3prime,,ScPCA paediatric AML; 10x 3' v3
Granja_2019,scAML_Granja_2019,3prime,,MPAL; 10x 3' v2
Kuusanmaki_2023,scAML_Kuusanmaki_2023,3prime,,venetoclax ex vivo; 10x 3' (guess)
Lambo_2023,scAML_Lambo_2023,3prime,,paediatric AML; 10x 3' v3
Lasry_2022,scAML_Lasry_2022,3prime,GSE185381,on disk; pooled/hashed libraries
Li_2022,scAML_Li_2022,3prime,,relapsed/refractory AML; 10x 3' (guess)
Montefiori_2021,scAML_Montefiori_2021,3prime,,MPAL/T-ALL only (SJMPAL011913 + SJTALL005006); 10x (guess)
Mumme_2023_AML,scAML_Mumme_2023_AML,3prime,,paediatric AML (Emory); 10x 3' v3
Mumme_2023_MPAL,scAML_Mumme_2023_MPAL,3prime,,paediatric MPAL (Emory); 10x 3' v3
Naldini_2023,scAML_Naldini_2023,3prime,,longitudinal chemo response; 10x 3'
Petti_2019,scAML_Petti_2019,5prime,Petti2019,on disk
RodriguezMeira_2023,scAML_RodriguezMeira_2023,TARGETseq,,plate-based full-length (TARGET-seq); no 10x reference applies
Wang_2024,scAML_Wang_2024,unknown,,platform not recalled
Zeng_2025,scAML_Zeng_2025,5prime,GSE289435,PMH in-house = GSE289435 on disk (MLL_xxxxx ids; GSM8791436 = MLL_17844); 5' per manifest; 3' v3.1 for sorted fractions
vanGalen_2019,scAML_vanGalen_2019,SeqWell,GSE116256,on disk; AML707B-D0 reproduction check
")
LOY_EXT_STUDIES[, dataset := paste0("ZENG_", study)]
LOY_BMM_BIN_MAP   <- here::here("scripts", "03_hierarchy", "bmm_bin_map.tsv")   # CellType_Broad -> hierarchy_bin (same vocabulary)
LOY_EXT_LIST      <- file.path(LOY_SCAN_DIR, "external_datasets.txt")           # one dataset id per line, for --task
# 7b. calibration for external datasets (02, declared 2026-11-14 BEFORE any external sample is called):
#   factor_study = (Y rate of the study's expression-male T cells) / r_ref(chem_lit, T); lambda_i uses r_ref x factor.
#   E0_i = max(exp(-lambda_i), e0_emp(chem_lit, comp, depth bin)) -- conservative whichever way the guess is wrong.
#   Studies with < LOY_EXT_CALIB_MIN_T male T cells keep factor 1 and are flagged "uncalibrated".
#   On-disk datasets keep the frozen section-5 rule unchanged (so K1-K7 still reproduce).
LOY_EXT_CALIB_MIN_T  <- 300L
LOY_EXT_CALIB_RANGE  <- c(0.5, 2.0)   # a factor outside this range means the chemistry guess is wrong: report, do not call
# 7c. GSE269669 (Daver 2025, PMID 40198272): 11 TP53-mutant AML/AEL on azacitidine + venetoclax + magrolimab, 27 sample
#   libraries (PTnnA/B/C = patient x timepoint), author Seurat v4 objects: avm_aml (42,910 leukaemic cells, 20 libraries)
#   and avm_tme (137,169 microenvironment cells, 27 libraries), raw counts, Azimuth-style labels (predicted.celltype.l2).
#   The deposited .rds.gz is gzip of a gzip-compressed RDS (readRDS(gzcon(gzfile()))). Chemistry guessed 5' (VDJ columns).
LOY_EXT_GEO <- data.table::data.table(
  dataset = "GSE269669", chem_lit = "5prime", on_disk = NA_character_,
  note = "TP53-mutant AML/AEL aza/ven/magrolimab; author objects avm_aml + avm_tme; 10x 5' (guess, VDJ columns present)")
LOY_GSE269669_DIR <- file.path(LARGE1_DIR, "00_raw/public/GEO/GSE269669_RAW")
LOY_EXT_ALL <- rbind(LOY_EXT_STUDIES[, .(dataset, chem_lit, on_disk)], LOY_EXT_GEO[, .(dataset, chem_lit, on_disk)])
# Azimuth-style l2 labels -> hierarchy_bin, aligned to bmm_bin_map.tsv (CLP -> T_NK like Early Lymphoid; BaEoMa ->
# LMPP_GMP like EoBasoMast Precursor; EMP -> Erythroid like MEP). Unlisted labels become NA (= not in any compartment).
LOY_AZIMUTH_L2_TO_BIN <- c(
  "HSC" = "HSC_MPP", "LMPP" = "LMPP_GMP", "GMP" = "LMPP_GMP", "BaEoMa" = "LMPP_GMP",
  "EMP" = "Erythroid", "Early Eryth" = "Erythroid", "Late Eryth" = "Erythroid",
  "Prog Mk" = "Megakaryocyte", "Platelet" = "Megakaryocyte",
  "CD14 Mono" = "Mono_DC", "CD16 Mono" = "Mono_DC", "Macrophage" = "Mono_DC", "cDC1" = "Mono_DC", "cDC2" = "Mono_DC",
  "pDC" = "Mono_DC", "pre-mDC" = "Mono_DC", "pre-pDC" = "Mono_DC", "ASDC" = "Mono_DC",
  "CLP" = "T_NK", "CD4 Naive" = "T_NK", "CD4 Memory" = "T_NK", "CD4 Effector" = "T_NK", "CD8 Naive" = "T_NK",
  "CD8 Memory" = "T_NK", "CD8 Effector_1" = "T_NK", "CD8 Effector_2" = "T_NK", "CD8 Effector_3" = "T_NK",
  "T Proliferating" = "T_NK", "MAIT" = "T_NK", "ILC" = "T_NK", "NK" = "T_NK", "NK CD56+" = "T_NK", "NK Proliferating" = "T_NK",
  "pro B" = "B_Plasma", "pre B" = "B_Plasma", "transitional B" = "B_Plasma", "Naive B" = "B_Plasma", "Memory B" = "B_Plasma",
  "Plasma" = "B_Plasma", "Stromal" = "Stromal")
# 7d. which Disease_state values count as AML for the patient table (pack vocabulary: AML, sAML, AEL, AMKL, MPAL, MPN,
#   sAML_TP53, MPN_TP53_preAML, MPN_TP53_CP, AML_PDX, AML_Multiome, Healthy). MPAL and MPN stages are not AML.
LOY_AML_DISEASE <- c("AML", "sAML", "sAML_TP53", "AEL", "AMKL", "AML_PDX", "AML_Multiome")
