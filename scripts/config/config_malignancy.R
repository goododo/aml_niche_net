# config_malignancy.R ----
# Stage config for 02_malignancy (Phase 1B). Sources paths, adds all CNV/SNV/consensus constants.
# Fully absorbs legacy c00_refnorm_config.R (refnorm arm) + the inferCNV-runner constants that
# used to live script-locally in c40/c46 + c50 tier vocabulary + Numbat params.
#
# Naming decisions (3a-2 config port):
#   * The BoneMarrowMap ANNOTATED dataset (counts + CellType_Broad; used by BOTH the SingleR
#     refnorm arm AND the inferCNV external reference) is BMM_ANNOTATED_RDS -- distinct from the
#     Symphony projection reference config_paths.R::BMM_REF_RDS (batch-4 hierarchy). Resolves the
#     old BMM_REF_RDS name collision.
#   * Redundant path aliases are NOT reintroduced: scripts reference the canonical config_paths.R
#     names (QC_RDS_DIR, INFERCNV_ROOT, INFERCNV_BURDEN_ROOT) rather than legacy REFNORM_QC_OBJ_DIR
#     / INFERCNV_OUT_ROOT / INFERCNV_BURDEN_DIR.
#   * ref_norm_summary / infercnv_summary now write under DIR_MALIGNANCY (02_malignancy), not the
#     legacy 03_malignancy.

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
# CANONICAL_TIMEPOINTS lives in config_qc.R and RESIDUAL_NOMINAL_SCOPE is validated against it
# below. Sourced here rather than left to the caller: 20_refnorm_identify.R sources only paths +
# malignancy, so the validation exploded with "object 'CANONICAL_TIMEPOINTS' not found" the first
# time this config was loaded on its own. A config that asserts against another config must pull
# it in itself -- the same pattern config_hierarchy.R already follows.
source(here::here("scripts", "config", "config_qc.R"))
suppressPackageStartupMessages({ library(data.table) })

## -- inferCNV: detection + burden thresholding ----
INFERCNV_CUTOFF     <- 0.1    # expression DETECTION floor (not a malignancy threshold)
REF_FLATNESS_MAX    <- 0.02   # reject a reference group if its CNV profile is too flat
INFERCNV_SCORE_Q    <- 0.95   # per-cell burden threshold = P95 of reference-cell burden

## -- Numbat ----
NUMBAT_T           <- 1e-5    # HMM transition prob
NUMBAT_MAX_ENTROPY <- 0.5     # single-cell CNV entropy filter. NOTE: raising to 0.8 was TESTED on
                              # 6323_MRD and produced catastrophic false positives (T_NK tumor_frac
                              # 44%); 0.5/0.6 give no CNV for genuinely CN-normal samples. Keep 0.5;
                              # accept inferCNV-single (tier C) when Numbat finds nothing. [D-A7]

## =====================================================================================
## Reference-normal identification (SingleR whole-profile; legacy c00_refnorm_config.R) ----
## Identifies the autologous normal reference for inferCNV only; malignancy is still called on
## ALL cells by inferCNV+Numbat+VarTrix (blueprint R2 preserved).
## =====================================================================================

## -- BoneMarrowMap ANNOTATED reference (SingleR ref AND inferCNV external ref) ----
BMM_ANNOTATED_RDS <- file.path(ZENODO_DIR, "BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds")
BMM_LABEL_COL     <- "CellType_Broad"
# Aggregated (pseudobulk) SingleR reference, built once from BMM_ANNOTATED_RDS and cached.
SINGLER_AGG_RDS   <- file.path(LARGE1_DIR, "reference", "singleR_aggref_BMM_broad.rds")
SINGLER_FINE_TUNE <- TRUE

# Mature-lymphocyte labels (from CellType_Broad) kept as the autologous normal reference.
# NK deliberately held back as optional backup (add if too many samples fall to fallback).
# This restriction is CORRECT and stays: in AML the myeloid/progenitor compartment may itself be
# malignant, so autologous myeloid cells cannot serve as normals. What is NOT correct is letting
# a lymphoid-only reference be the ONLY reference -- see INFERCNV_REF_PER_LINEAGE below.
REF_KEEP_LABELS <- c("CD4 Memory T", "CD8 Memory T", "Naive T", "B")

## -- LINEAGE-MATCHED REFERENCE GROUPS ----
# THE DOMINANT FALSE-POSITIVE MODE. On the v1 calls the per-bin healthy false-positive rate rose
# monotonically with transcriptional distance from the lymphoid reference:
#     T_NK 0.036 | B_Plasma 0.077 | LMPP_GMP 0.143 | Mono_DC 0.234 | Megakaryocyte 0.243
#     Erythroid 0.267 | HSC_MPP 0.398 | Stromal 0.903
# That ordering is lineage, not copy number. It made healthy marrow score as malignant as AML
# marrow (healthy median frac_malignant 0.135 vs diagnosis 0.115; healthy-vs-diagnosis AUC 0.546)
# and it is worst in HSC_MPP -- the bin H1's LSC-like sender axis depends on.
#
# MECHANISM: infercnv::run(ref_subtract_use_mean_bounds = TRUE, the package default) zeroes a
# gene only when it falls outside [min, max] of the PER-REFERENCE-GROUP means. With ONE reference
# group that interval collapses to a point and the step degenerates into plain mean subtraction,
# so a gene that is merely high in myeloid cells reads as amplification against a lymphoid
# reference. Splitting the reference into one group per lineage restores a real interval and the
# lineage component cancels. No other inferCNV parameter needs to change.
#
# TRADE: a wider reference interval is strictly more conservative -- it costs sensitivity to
# genuine small CNVs in exchange for removing the lineage artefact. With a healthy FPR of 0.4 in
# HSC_MPP, that is the right direction. Verified against 96_malignancy_fpr_healthy after the run.
INFERCNV_REF_PER_LINEAGE   <- TRUE
INFERCNV_REF_MIN_PER_GROUP <- 20L   # a lineage with fewer reference cells is dropped rather than
                                    # kept as a noisy group that would widen the interval by luck

## -- reference adequacy & fallback ----
REFNORM_MIN_NORMAL_CELLS <- 50L   # need at least this many autologous normals, else fall back.
                                  # (was 30 in legacy c00; raised to 50 for P95-threshold stability.
                                  # Samples in the 30-49 band now route to REF_FALLBACK_MODE -- audit
                                  # that band against external-ref platform confound at run time.)
REF_FALLBACK_MODE <- "external"   # "external" | "skip_infercnv" | "halt"

## -- lenient blast-negative gate (SingleR is the PRIMARY filter; this is a second safety) ----
# Drop a SingleR-called T/B cell only if it co-expresses MORE than this many DISTINCT blast/myeloid
# markers (>1 = real co-expression / doublet, not a single ambient molecule). Deliberately tolerant
# to avoid the ambient over-kill (the NKG7-on-blast trap) that broke the detection-based approach.
BLAST_GATE_MAX_MARKERS <- 2L
MYELOID_BLAST <- c("CD34", "KIT", "MPO", "ELANE", "PRTN3", "AZU1",
                   "CD33", "LYZ", "CD14", "FCN1")
# Lineage markers used ONLY for the self-validation printout (precision audit).
LINEAGE_VALID <- c("CD3D", "CD3E", "TRAC", "CD79A", "MS4A1", "CD19")

## -- sorted / lineage-enriched libraries (explicit fallback guard; no autologous attempt) ----
# This guard means "this library has no normal-cell population to use as an autologous inferCNV
# reference". GSE185991/GSE147989 are genuinely CD34/CD117-sorted blast libraries.
#
# REMOVED: SORTED_LIBRARY_REGEX <- "_CD34$". Chen2023's "<donor>_CD34" is not a CD34-sorted
# library -- per Methods it is pool A = HSPC + myeloid, i.e. half of a donor whose other half is
# in "<donor>_Niche_Immune". The regex therefore routed half of every Chen2023 donor to the
# fallback and, worse, computed malignant_frac from the LYMPHOID pool only. Since
# ingest_Chen2023.R now merges both pools under one donor-level Sample, no Chen2023 sample id
# ends in "_CD34" any more and the regex would be dead code as well as wrong. Left empty rather
# than deleted so the variable stays part of the guard's contract.
SORTED_LIBRARY_DATASETS <- c("GSE185991", "GSE147989")
SORTED_LIBRARY_SAMPLES  <- character(0)
SORTED_LIBRARY_REGEX    <- ""

## -- CNV-UNINFORMATIVE DATASETS: a near-zero malignant_frac here is not a low tumour burden ----
# [2026-08-14] These 33 samples were backfilled through the external BMM reference (they had been
# dropped entirely by a routing bug in infercnv_routes(); see 00_infercnv_common.R). They now have
# burden files, and their calls came out at a median malignant_frac of 0.0151 -- 4-5x below either
# other route, with 39% of samples calling under 1%. For CD34/CD117-SORTED BLAST libraries that
# looks like a caller failure, and it is not:
#
#   GSE185991 (29 of the 33) is NPM1-mutant AML. NPM1-mutant AML is the textbook CYTOGENETICALLY
#   NORMAL subtype -- roughly 85% carry a normal karyotype. inferCNV infers copy-number change.
#   A karyotypically normal leukaemia has no copy-number change to infer, so ~0 is the RIGHT
#   answer from this instrument about these cells, and it says nothing about whether they are
#   malignant. They are: they were sorted as blasts.
#
# So malignant_frac for these datasets must be read as "no CNV evidence", never as "few malignant
# cells", and they must not enter a healthy-vs-AML contrast as though their value were comparable
# to a karyotypically abnormal sample's. Their cells remain fully usable for hierarchy and CCC.
#
# There is no orthogonal check available: karyotype and mutations are empty for all 46 samples in
# these two datasets, and neither has Numbat coverage. A threshold contribution cannot be ruled
# out either -- their observation/reference burden ratio is 1.78 vs 2.95 for the 44 samples that
# took the same external route successfully -- but the NPM1 subtype explains the bulk of it and is
# the only explanation with evidence behind it.
CNV_UNINFORMATIVE_DATASETS <- c("GSE185991", "GSE147989")

## =====================================================================================
## VAN GALEN MALIGNANT-STATE AXIS (80/81) -- a malignancy call that does not use copy number
## =====================================================================================
# Built because the CNV route is inverted on negative controls (healthy median malignant_frac
# 0.108 vs autologous-AML 0.063) and structurally blind to cytogenetically normal AML. Signatures
# are learned from van Galen's PredictionRefined labels by a WITHIN-SAMPLE malignant-vs-normal
# contrast, then scored cohort-wide with UCell -- rank-based within each cell, which is what lets
# a signature learned on Seq-Well transfer to 10x at all.
VG_MIN_CELLS_PER_SIDE   <- 10L    # a (sample,bin) needs this many malignant AND normal cells.
                                  # 20 left only Mono_DC(9) and LMPP_GMP(4) testable; 10 gives
                                  # Mono_DC 10, LMPP_GMP 8, HSC_MPP 5, Erythroid 4, T_NK 4, B_Plasma 4.
                                  # Most AML samples retain almost no NORMAL myeloid cells, so this
                                  # is the binding constraint on contrast A, not a quality knob.
VG_MIN_SAMPLES_PER_BIN  <- 4L     # a bin needs this many testable samples before it gets a signature
VG_MIN_LFC              <- 0.25   # on log-normalised expression
VG_MIN_DDET             <- 0.05   # and the detection rate must move too, not just the mean
VG_MIN_FRAC_SAMPLES     <- 0.70   # CONSISTENCY, not effect size: one donor cannot install its genes
VG_MIN_LFC_B            <- 0.25   # contrast B (vs the 4 healthy donors inside GSE116256) must agree
VG_HEALTHY_REGEX        <- "^BM"  # those donors: BM3, BM4, BM5-34p, BM5-34p38n
VG_TOP_N                <- 50L    # genes per direction per bin

# The threshold is set FROM the healthy donors rather than derived from a reference quantile and
# then discovered to be miscalibrated. Pick the operating point on the negative controls, then
# report whatever sensitivity it buys on the genotyped cells -- the opposite order to the CNV gate,
# and the reason that gate ended up at FPR 0.186 against a nominal 0.05.
VG_TARGET_HEALTHY_FPR   <- 0.05

## -- refnorm I/O (QC objects come from the canonical QC_RDS_DIR) ----
REFNORM_SUMMARY_CSV  <- file.path(DIR_MALIGNANCY, "ref_norm_summary.csv")   # 02_malignancy (was 03)
REFNORM_REF_CELL_DIR <- file.path(LARGE1_DIR, "02_seurat_objects", "01b_ref_norm_cells")
REFNORM_OUT_OBJ_DIR  <- file.path(LARGE1_DIR, "02_seurat_objects", "01b_ref_norm_obj")
REFNORM_WRITE_RDS    <- FALSE      # OPTIONAL: rewrite per-sample RDS with a `ref_norm` column
REFNORM_MANIFEST_CSV <- ""         # optional manifest; "" -> auto-discover under QC_RDS_DIR
# HEALTHY-ONLY DATASETS ARE NOT SKIPPED. They used to be, on the reasoning that a
# purely-healthy cohort has no malignancy to call. That reasoning ignores what the
# healthy samples are FOR: a normal marrow contains no malignant cells, so every
# cell called malignant there is a measured false positive, and that measurement is
# the blocking calibration gate in 70_residual_stratum. Skipping E-MTAB-11536 (7)
# and GSE253355 (12) removed 19 of 37 healthy samples -- more than half the negative
# control, and specifically the two INDEPENDENT healthy cohorts, leaving the gate to
# be judged on healthy samples collected inside AML studies. The v1 gate (AUC 0.546)
# was read off 22 such samples. Running them costs 19 inferCNV jobs.
REFNORM_SKIP_DATASETS <- character(0)

## =====================================================================================
## inferCNV RUNNER (was script-local in c40/c46) ----
## Output roots come from config_paths.R: INFERCNV_ROOT, INFERCNV_BURDEN_ROOT.
## =====================================================================================
INFERCNV_GENE_ORDER  <- file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt")  # from 01
INFERCNV_ANALYSIS    <- "subclusters"   # "subclusters" (finer) or "samples" (faster)
INFERCNV_HMM         <- TRUE
INFERCNV_DENOISE     <- TRUE
INFERCNV_MIN_OBS     <- 3L               # crash-guard: with <3 observation (non-reference) cells the
                                         # residual matrix collapses to a vector (rowMeans error) and
                                         # inferCNV is meaningless anyway (near-all-normal sample) ->
                                         # skip. Does NOT invalidate completed runs (they had >=3).
INFERCNV_THREADS     <- as.integer(Sys.getenv("INFERCNV_THREADS", "4"))  # env-overridable (array sets 8)
INFERCNV_SUMMARY_CSV <- file.path(DIR_MALIGNANCY, "infercnv_summary.csv") # 02_malignancy (was 03)

# External reference for singleR-fallback / blast-packed samples: a stratified BMM subsample.
INFERCNV_EXT_REF_RDS      <- BMM_ANNOTATED_RDS   # same annotated atlas as the SingleR ref
INFERCNV_EXT_REF_LABELCOL <- BMM_LABEL_COL
INFERCNV_EXT_REF_PER_TYPE <- 150L                # cells sampled per broad type
INFERCNV_EXT_REF_CACHE    <- file.path(LARGE1_DIR, "reference", "infercnv_external_ref_BMM.rds")
INFERCNV_EXT_REF_SEED     <- 20260613L           # LEGACY cache-provenance seed, NOT the global SEED;
                                                 # changing it re-samples the external ref on rebuild.

## ---- dataset-matched healthy reference -------------------------------------------------------
# The BMM external reference is drawn from a DIFFERENT study than the sample being scored, so the
# non-lymphoid half of every reference is a cross-dataset comparison and batch effects enter the
# burden directly. Where a dataset ships its own healthy donors, those donors are the same
# chemistry, depth and handling as the AML samples beside them, so they are the better reference.
#
# Applies to 71 of 142 AML samples (GSE185381 42, GSE116256 23, Chen2023 6); the remaining 71 have
# no in-dataset healthy donor and keep the BMM block. Per BIN, matched cells are used only if at
# least INFERCNV_MATCHED_REF_MIN_PER_BIN of them exist -- otherwise that bin falls back to BMM, so
# a matched block is never SMALLER than the BMM block it replaces. That floor matters: 98 showed a
# short reference block gives a noisy, biased-low P95 and makes the gate worse, not better.
#
# HONESTY CONSTRAINT: a healthy donor must never appear in the reference used to score itself, or
# the healthy-donor FPR -- the number this change is judged by -- is measuring the reference
# against itself. build_infercnv_input() drops the target sample's own cells (leave-one-out).
# VERDICT 2026-08-14: OFF. Measured and rejected -- see 93_matched_ref_verdict.R.
# On 18 healthy donors it looked excellent: FPR 0.186 -> 0.099, 17 of 18 improved, and the two
# CD34-sorted donors (no autologous reference at all) fell from 0.59/0.53 to 0.06/0.12.
# On the 18 GSE116256 AML samples with single-cell genotyping it collapsed: sensitivity on
# genotypically malignant cells 0.492 -> 0.131. True positives fell ~3x more than false ones,
# because a same-dataset healthy reference of the SAME lineage subtracts the malignant signal
# along with the batch effect. AML blasts and healthy GMP/Mono are close in expression; using the
# latter as the reference for the former removes what we are trying to measure.
# The healthy cohort CANNOT see this: it has no positive cells. Judging a malignancy gate on
# negative controls alone systematically favours any change that lowers all burdens.
INFERCNV_MATCHED_REF          <- FALSE
INFERCNV_MATCHED_REF_MIN_PER_BIN <- 150L  # absolute floor; the EFFECTIVE floor per bin is that
                                          # bin's BMM block size (see .matched_ref_lineage_block)
INFERCNV_MATCHED_REF_CACHE_DIR   <- file.path(LARGE1_DIR, "reference", "infercnv_matched_ref")
INFERCNV_MATCHED_REF_SEED        <- 20260811L
# per-bin targets, matched to the BMM block so the two sources are the same size
INFERCNV_MATCHED_REF_PER_BIN  <- c(Mono_DC = 600L, LMPP_GMP = 600L, Erythroid = 450L,
                                   HSC_MPP = 300L, Megakaryocyte = 150L, Stromal = 150L)

## =====================================================================================
## CONSENSUS (50 is the SOLE malignancy labeler; d35 retired) ----
## Evidence TYPES: expression {inferCNV, copykat, scevan, author}, allele {Numbat}, SNV {VarTrix}.
## =====================================================================================
CONSENSUS_MIN_TOOLS <- 1L     # pilot: single-type allowed (tier C). Raise for strict >=2/3.
CONSENSUS_UNION_MODE <- TRUE   # PRODUCTION DEFAULT: with >=2 evidence types, malignant if ANY type
                              # is positive. Numbat allele evidence (del/LOH) is the ONLY way to
                              # catch CN-normal AML (inferCNV is expression-flat there); type-majority
                              # would drop those samples' malignant cells (e.g. 3853_Dg/6323_Dg went
                              # 0 -> 0.66/0.50 under union). Labels use union; the B_multi_partial
                              # tier flags cross-type disagreement for downstream QC.

# UNIFIED tier vocabulary (first letter aligns with TIER2CONF A=high/B=medium/C=low):
#   A_concordant    : >=2 evidence types AND agree
#   B_multi_partial : >=2 evidence types but only partial agreement (union legitimately mixes)
#   C_single        : a single evidence type
TIER_LEVELS <- c("A_concordant", "B_multi_partial", "C_single")
TIER2CONF   <- c(A = "high", B = "medium", C = "low")
# A_concordant vs B_multi_partial split: among cells with >=2 evidence types, the fraction whose
# types DISAGREE (not unanimous: 0 < n_types_mal < n_types). Below this -> A_concordant; at/above
# -> B_multi_partial. E.g. an inferCNV~0 / Numbat-high blind-spot sample has a high conflict
# fraction and correctly lands in B.
TIER_CONFLICT_MAX <- 0.10

# Numbat schema handling (BUG FIX carried from d35): a degraded Numbat percell file whose
# note == "no_CNV_detected" means the sample contributes NO allele evidence -> treat as NA
# (abstain), NOT malignant=0. Encoding it as 0 would falsely vote "normal" and pollute consensus.
NUMBAT_DEGRADED_NOTE <- "no_CNV_detected"

## -- SNV arm [PLANNED, NOT ACTIVE] ----
# Most datasets lack matched genotyping/VarTrix input. 50 keeps a --vartrix interface, but the
# tier logic does NOT require SNV; when absent, evidence is expr(+allele) only. 48_vartrix_run.R
# is a placeholder until a dataset with usable SNV data is added.
SNV_ARM_ACTIVE <- FALSE

## -- RESIDUAL-DISEASE STRATUM (H3 treatment-pressure axis) ----
# WHY this exists: the nominal Timepoint is the depositing author's word for the draw, and it
# conflates two different things -- GSE201966 "Complete_remission" is clinically an MRD draw,
# and GSE116256 maps every "D<n>" to the same class whether n is 14 or 171. In v1 only GSE227903
# ever landed in "MRD", purely because its authors wrote the word. So H3's stratum is derived
# from residual disease BURDEN instead of from the label. Applied ONLY to post-treatment draws:
# a low-blast diagnosis marrow is not an MRD sample.
#
# Tracks the v2 vocabulary. "MRD" and "Post_treatment" no longer exist as canonical values; the
# post-therapy phases that replaced them are listed explicitly. Diagnosis/Refractory/Relapse are
# OUT (never had a response to be residual to); On_treatment is OUT because the drug is still on
# board and the blast count is mid-transit, not a response measurement.
RESIDUAL_NOMINAL_SCOPE <- c("Post_induction", "Post_consolidation",
                            "Post_transplant", "Post_treatment_unspecified")
stopifnot(all(RESIDUAL_NOMINAL_SCOPE %in% CANONICAL_TIMEPOINTS))
RESIDUAL_CUTS   <- c(deep = 0.05, residual = 0.20)   # frac_malignant boundaries
RESIDUAL_LABELS <- c("Deep_response", "Residual", "Active_disease")
# Sensitivity sweep: the cut points are a clinical convention (CR is <5% blasts by morphology),
# not a property of this cohort, so every H3 result must be shown to survive moving them.
RESIDUAL_SWEEP_DEEP     <- c(0.02, 0.03, 0.05, 0.08, 0.10)
RESIDUAL_SWEEP_RESIDUAL <- c(0.15, 0.20, 0.25, 0.30)

## -- CALIBRATION GATE: frac_malignant must actually separate healthy from AML ----
# BLOCKING PRECONDITION, checked before any stratum is emitted. On the v1 (2026-07-15) calls the
# consensus did NOT separate disease state at all:
#     Healthy   n=22  median frac_malignant 0.135   <- HIGHER than diagnosis
#     Diagnosis n=77  median 0.115
#     Relapse   n=10  median 0.048
#     MRD       n= 6  median 0.351
# with per-sample healthy FPR up to 0.756 (GSE116256 BM5-34p38n, a CD34-sorted healthy fraction)
# and per-bin FPR of 0.398 in HSC_MPP / 0.903 in Stromal. 125/130 samples carried single-arm
# (inferCNV-only) evidence. Stratifying treatment pressure on a score with that behaviour would
# manufacture an H3 result out of caller noise, so 70_residual_stratum.R refuses to write.
# Part of the v1 signal is a labelling artefact this refactor fixes: several GSE185381 "healthy"
# entries were POOLED LIBRARIES containing AML donors (GSM5613770 = Control0082_AML052_AML022b),
# so their "false positives" were real AML cells. Donor-level Sample (ingest_GSE185381.R) removes
# that class. GSE116256's sorted healthy fractions and Chen2023 NBM are NOT explained by it.
RESIDUAL_CALIB_MAX_HEALTHY_MEDIAN <- 0.05  # healthy median frac_malignant must sit below the deep cut
RESIDUAL_CALIB_MIN_AUC            <- 0.70  # healthy vs diagnosis separation (Mann-Whitney AUC)
RESIDUAL_CALIB_MAX_HEALTHY_FPR    <- 0.20  # per-sample healthy FPR ceiling (median across healthy)
