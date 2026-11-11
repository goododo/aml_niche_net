# config_dbl.R ----
# INPUT  : none (sources scripts/config/config_paths.R via here)
# OUTPUT : constants only -- every rule and number the doublet-pairing SMOKE test uses
# what it does: one place for the smoke test of *_project_context/PROPOSAL_doublet_pairing_external_v1.3.md, set up on
#   2026-11-09. This is NOT a pre-registration and nothing here is frozen. Each constant cites where it comes from:
#     [P x.y]  the proposal section;  [H x.y]  HANDOFF_doublet_contact_v1.1.md;
#     [D1]-[D5] the user's five decisions of 2026-11-09 (purpose, data + masking, folder, scDblFinder version, state map);
#     [SMOKE]  a default chosen for the smoke only -- the real value is set at the proposal's first freeze.
# Business scripts start with:
#   source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))   # SEED, LARGE1_DIR, TAB_DIR, GEO_RAW_DIR, RDS_INGEST_DIR, ...

DBL_CONFIG_DATE <- "2026-11-09"

## -- 0. scope, folders, masking  [D1, D2, D3] ----
DBL_DATASET    <- "GSE185381"
DBL_SMOKE_LANES <- c(single = "GSM5613785_2021-03-09-AML3948",     # one donor (AML3948); no hashtag rows
                     pooled = "GSM5613744_2019-07-01-count-1")     # P01 lane 1: 5 adult AML donors; 7 HTO rows
DBL_SMOKE_DIR  <- file.path(LARGE1_DIR, "dbl", "smoke")            # outside the repo: the 23:30 cron never sees it
DBL_LANE_DIR   <- file.path(DBL_SMOKE_DIR, "lanes")                # per-lane objects written by 01/02
DBL_TAB_DIR    <- file.path(DBL_SMOKE_DIR, "tables")               # smoke tables (technical numbers only)
DBL_CODE_DIR   <- file.path(DBL_SMOKE_DIR, "code")                 # dated read-only copies run by submit_dbl.sbatch
DBL_REPORT_DIR <- "/FAST/gr10634/gaozy/reports/2026-11_dbl_pairing_smoke"   # outside the repo
DBL_MASK_REAL_TYPING <- TRUE   # [D2] typed calls of REAL barcodes (by state, T-LEUK totals) are never written or printed

## -- 1. barcode classes  [P 4.1 "Eligibility from pre-removal objects"] ----
# A        author-assigned (Demuxed) and kept by project QC (present in 04_annotated)
# B_union  author-assigned, passes its donor's full QC thresholds, absent from 04_annotated -> removed by the doublet union
# B_mad    author-assigned, fails one of its donor's QC thresholds (the MAD/absolute filters of 03_per_sample_qc.R)
# C        deposited with counts but no author donor assignment (Sample == "UNDEMUX__<lane>")
DBL_CLASSES <- c("A", "B_union", "B_mad", "C")

## -- 2. line QC for classes B and C  [H 4.4.3] ----
# Lower tails only for nCount/nFeature and the upper tail for percent.mt: a doublet carries ~2x UMI, so the upper count tail
# that 01_preprocess cut is exactly where doublets sit. Class B uses its own donor's effective thresholds from
# 03_qc_report__GSE185381.csv (ncount_lo, nfeat_lo, mt_hi). Class C has no donor:            [SMOKE]
DBL_QC_C_RULE <- "median"      # class C threshold = median over the lane's donors of the same three thresholds

## -- 3. labels  [P 3.4, D5] ----
DBL_LABEL_MAP_FILE <- "dbl_label_map.tsv"   # bmm_fine -> compartment (scripts/dbl/)
# class-A cells whose projection failed QC carry bmm_fine = NA and a marker-filled hierarchy_bin: fall back on the bin, never
# on a state (the state is unknown)                                                                    [SMOKE; found 2026-11-09]
DBL_BIN_FALLBACK <- c(HSC_MPP = "LEUK_other", LMPP_GMP = "LEUK_other", Megakaryocyte = "LEUK_other", Mono_DC = "MonoDC",
                      Erythroid = "Ery", T_NK = "Tamb", B_Plasma = "B", Stromal = "Other")
DBL_LEUK_STATES <- c("LEUK_HSC", "LEUK_MPP", "LEUK_LMPP", "LEUK_GMP", "LEUK_ProMono")   # [P 3.4, D5] differentiation order
DBL_MIN_BMM_PROB <- 0.5        # [D5] a LEUK state or strict T needs bmm_prob >= 0.5 (= high_error == 0); else LEUK_other / Tamb
DBL_T_REQUIRE_CONCORDANT <- TRUE                     # [P 3.4] "annotated as T and confirmed": anno_source == "concordant"
DBL_T_ID_GENES   <- c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2")   # [P 3.4] "confirmed by T-identity transcripts"
DBL_T_ID_MIN_UMI <- 1L                               # [SMOKE] summed UMI over DBL_T_ID_GENES
DBL_ATLAS_EXCLUDE_T_CL <- 13L  # [P 3.4] the integrated cluster annotated T/NK that is cycling myeloid progenitors (ws atlas)
DBL_ATLAS_HOMOTYPIC_CL <- 28L  # [H 4.7.4] "Unresolved (homotypic doublet)": positive control for the caller
# kNN label transfer from the lane's class-A cells to its B and C barcodes  [H 4.4.4; SMOKE values]
DBL_TRANSFER_NHVG     <- 2000L
DBL_TRANSFER_NPC      <- 30L
DBL_TRANSFER_K        <- 15L
DBL_TRANSFER_MIN_VOTE <- 0.6   # vote share below this -> "Unassigned"

## -- 4. calling and typing  [P 4.1 "Typing with measured accuracy", D4] ----
DBL_DBR_PER_BARCODE <- 8e-6    # [H 4.4.5] 0.8 % per 1,000 DEPOSITED barcodes of the lane (pre-QC count)
DBL_DBR_CAP         <- 0.4     # scripts/ws/02_doublet_scores.R rate cap
DBL_ADJUST_SIZE     <- 0       # [D4] scDblFinder 1.20.2 workaround for issue #145 (origin labels); devel keeps its default
DBL_DEVEL_LIB       <- "/FAST/gr10634/gaozy/rlib_scdbl_devel"   # [D4] separate library: scDblFinder 1.27.12 @ 1a0831bd

## -- 5. synthetic doublets and known answers  [P 4.1, appendix C; H 4.7] ----
DBL_N_FOLDS        <- 5L       # [SMOKE] held-out cells within a lane: each run holds out one fifth of class A (removing half,
DBL_FOLDS_RUN      <- 1:2      #   as a 2-fold split would, distorts the lane); the smoke runs folds 1-2. The proposal's
                               #   "held-out patients" is infeasible in a one-donor lane (no scDblFinder model carries across
                               #   lanes); the pooled lane adds one leave-one-donor-out run.
DBL_SYN_PER_STATE  <- 30L      # [SMOKE] synthetic T-LEUK(k) doublets per state, fold and batch
DBL_SYN_PER_OTHER  <- 10L      # [SMOKE] synthetic doublets per other (false-positive source) pair type, fold and batch
DBL_SYN_SINGLETS   <- 100L     # [SMOKE] held-out singlets re-entered as known barcodes per fold and batch (equal per compartment)
DBL_SYN_BATCHES    <- 2L       # [SMOKE] injection batches per fold: about 3 % injected barcodes per run
DBL_INJECT_RATES   <- c(0.005, 0.01, 0.02)   # [H 4.7.2] planted contact doublets, fraction of the lane's barcodes
DBL_PLANT_STATE    <- "LEUK_GMP"             # [SMOKE] planted contacts all carry this state (a state-specific known answer)
DBL_AMBIENT_WINDOW <- c(0.25, 0.80)          # [H 4.5.6] LEUK-marker UMI share of a T-LEUK call; outside -> not retained
DBL_TYPING_MIN     <- 0.7                    # [P A.5] per-state recall and precision gate (reported, never decided here)

## -- 6. hashtags  [P 4.1 "Background calibration"] ----
DBL_HTO_PREFIX      <- "^HTO"   # Antibody Capture rows of the raw GEO matrix that are hashtags
DBL_HTO_BG_QUANTILE <- 0.995    # [SMOKE] a hashtag is positive above this quantile of its counts in author-assigned barcodes
DBL_HTO_MIN_POS     <- 5L       #   of OTHER donors (the background), and never below this many UMI
DBL_HTO_MAP_MIN     <- 0.8      # [SMOKE] a donor's hashtag must be the top hashtag of >= 80 % of its author-assigned barcodes

## -- 7. pairing model on a simulated cohort  [P 4.2, appendix A] ----
DBL_SIM_PATIENTS <- 16L         # [SMOKE] pseudo-patients built from the smoke lanes' singlet compositions and confusion
DBL_SIM_OUTER    <- 60L         # [SMOKE] simulated cohorts per world (null world and planted world)
DBL_SIM_BOOT     <- 199L        # [SMOKE] parametric reference draws per cohort (the proposal asks >= 2,000)
DBL_PLANT_R      <- 1.5         # [P A.4] follow-up threshold, planted on one state in the planted world
DBL_R_BAND       <- c(2 / 3, 1.5)   # [P A.4] negligible band

## -- 8. seeds: SEED (config_paths.R) plus one fixed offset per random stream ----
DBL_SEED_TRANSFER <- SEED + 801L
DBL_SEED_SYN      <- SEED + 811L   # + fold * 10 + batch
DBL_SEED_PLANT    <- SEED + 851L   # + rate index
DBL_SEED_CALL     <- SEED + 871L   # set before every scDblFinder call
DBL_SEED_SIM      <- SEED + 891L   # + world * 1000 + cohort
DBL_SEED_NBS      <- SEED + 901L   # + fold * 10 (+ 1 for the shuffled-label world); stage 10
## -- 9. Neighbor-seq-style pair classifier, stage 10 (user decision 2026-11-10) ----
# Defaults of Ghaddar & De 2022 (github.com/sjdlabgroup/Neighborseq, R/functions.R, read 2026-11-10), re-implemented without the
# package (its extra imports, multiROC/metap/ggraph, serve plots and the enrichment test, which the smoke does not use).
DBL_NBS_TOPN       <- 50L    # top positive markers per type: FindAllMarkers(logfc.threshold 1, min.pct 0.25, max 200 cells/type)
DBL_NBS_N_AM       <- 100L   # artificial barcodes per class: every pair type incl. homotypic, plus n cells of each type
DBL_NBS_TRAIN_FRAC <- 0.8    # stratified training share; the remaining 20 % is the classifier's own test (Neighbor-seq f = 0.8)
DBL_NBS_MIN_CELLS  <- 5L     # [SMOKE] a lineage type needs >= 5 training cells to become a class
DBL_NBS_NULL_MAX   <- 0.10   # [SMOKE] shuffled-label classifier: held-out synthetic T-LEUK-lineage recall must stay <= this
## -- 10. lineage-level O/E power simulation, stage 11 (user decision 2026-11-10) ----
DBL_OE_VARIANT   <- "M1"       # stage 10's confirmed typing (kNN labels, V5 AND T-identity gate; post-hoc exception)
DBL_OE_LABELS    <- "knn"
DBL_OE_CAL_TAGS  <- c("tune2_knn", "f3_knn")   # calibration counts pooled over tuning folds 1-2 and fold 3
DBL_OE_SD_YIELD  <- 0.5        # [SMOKE] log-sd of a simulated patient's yield multiplier around its template
DBL_OE_SD_HET    <- 0.3        # [SMOKE] log-sd of the true enrichment between patients (mean-centred: E[R_i] = R)
DBL_OE_SD_DBR    <- 0.2        # [SMOKE] log-sd of the doublet-rate error per patient; the analysis knows this width
DBL_OE_COHORTS   <- 200L       # [SMOKE] simulated cohorts per world
DBL_OE_MIN_PAT   <- 10L        # [H 4.8] fewer qualifying patients -> descriptive only (outcome D)
DBL_SEED_OE      <- SEED + 911L   # + world * 1000 + cohort
## -- 11a. Stage 0a, P01 hashtag desk check (reports/.../redesign_2026-11-11/STAGE0A_plan_v2.md; user decision 2026-11-11) ----
# P01 = one pool of 5 adult AML donors, 7 hashtags, three 10x lanes. count-1 IS the smoke lane "pooled" (extracted 2026-11-09)
# and is never re-extracted under a P01 key; count-2/3 come from `01 --lanes=p01_2,p01_3`, whose tables carry suffix _p01.
DBL_P01_LANES <- c(p01_1 = "GSM5613744_2019-07-01-count-1", p01_2 = "GSM5613745_2019-07-01-count-2",
                   p01_3 = "GSM5613746_2019-07-01-count-3")
stopifnot(identical(DBL_P01_LANES[["p01_1"]], DBL_SMOKE_LANES[["pooled"]]))   # never append P01 lanes to DBL_SMOKE_LANES (02-11 loop over it)
