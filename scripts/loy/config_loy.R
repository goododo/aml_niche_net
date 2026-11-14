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
