# config_ws.R ----
# INPUT  : none (sources scripts/config/config_paths.R via here)
# OUTPUT : constants only -- every threshold the within-sample line (WP1) uses
# what it does: T3c-0. Freezes, BEFORE T3c runs, every number that T3c and the WP1 main analysis
#   depend on, per §15.1: "对已冻结门槛的修订只能发生在相关结果出现之前". Each constant cites the
#   handoff section that fixed it (20261002_HANDOFF_within_sample_v1.5.md). After T3c has produced
#   any result this file is NOT edited; a change means a NEW dated constant plus an entry in GATES.md.
#   Two kinds of value live here and are kept apart on purpose:
#     FROZEN     -- sections 0-8. Decided and dated.
#     NOT FROZEN -- section 9. Suggestions that wait for numbers (GATE_W1.3) or that were proposed
#                   but never approved (the A[L,p,b] half-split). Nothing downstream may treat these
#                   as decided.
# Business scripts start with:
#   source(here::here("scripts", "config", "config_paths.R"))
#   source(here::here("scripts", "ws",     "config_ws.R"))
#   ws_assert_frozen()        # refuses to run while any GATE_W1.1 cell is still NA

suppressPackageStartupMessages(library(here))
source(here::here("scripts", "config", "config_paths.R"))   # SEED, FAST_DIR, LARGE1_DIR, TAB_DIR, ...

WS_CONFIG_DATE   <- "2026-10-12"
WS_CONFIG_FROZEN <- TRUE    # T3c-0 complete; sections 0-8 are not to be edited after T3c runs

## -- 0. cell scope  (§5.4 step 0; §0C.12j) ------------------------------------------------------
WS_STATE_UNIT          <- "hierarchy_bin"
WS_MAIN_BINS           <- c("LMPP_GMP", "Mono_DC")
WS_SUPP_BINS           <- c("HSC_MPP")                 # enough patients in GSE185381 only
WS_DROP_HIGH_ERROR     <- TRUE                          # cells with high_error == 1 are excluded
WS_MALIGNANCY_IN_MAIN  <- "none"                        # neither a filter nor a covariate (§0C.12j):
                                                        #  filter  -> 5.4% of AML cells callable
                                                        #  covariate -> constant within patient x bin
WS_MALIGNANCY_SENSITIVITY_COLS <- c("anno_malignant", "malignancy")   # sensitivity analyses only
WS_STUDIES_DISCOVERY   <- c("GSE185381", "GSE239721")
WS_STUDIES_REPLICATION <- c("GSE227903", "GSE289435")  # GSE116256 left WP1 (§5.6, W1.2b-ii)
WS_HC2_BINS            <- c("LMPP_GMP")                 # H-C2 main analysis only here, on the 22
WS_HC2_COMBOS_TABLE    <- "ws_w1_samebin_counts.csv"    #   intersection combos (§5.6); H-C1 all bins

## -- 1. GATE_W1.2b counting rule  (frozen before T3d; mirrors 01_samebin_counts.R) --------------
WS_LIG_MEASURABLE_LO    <- 0.01    # "ligand measurable" = detected in 1%-50% of the bin's cells
WS_LIG_MEASURABLE_HI    <- 0.50
WS_MIN_CELLS_BIN        <- 100L    # a bin must allow 50 + 50 at all
WS_MIN_PATIENTS         <- 10L     # per combination
WS_ARM_PRIMARY          <- 50L     # receptor-positive AND -negative cells, after removing L+ cells
WS_ARM_FALLBACK         <- 30L
WS_MIN_COMBOS_PER_STUDY <- 30L
WS_MIN_STUDIES          <- 2L

## -- 2. §15.2 item 2: response genes = genes the receiver expresses itself ----------------------
WS_RESPONSE_GENE_MIN_DETECT     <- 0.05   # detected in >= 5% of the receiver cells
WS_RESPONSE_GENE_MAX_FOLD_OTHER <- 4      # mean in any other cell type <= 4x the receiver mean
WS_RESPONSE_EXCLUDE_RECEPTOR    <- TRUE   # the receptor gene is never a response gene

## -- 3. §15.2 item 4: ligand availability A[L,p,b] --------------------------------------------
WS_LIGAND_AVAIL_PRIMARY   <- "mean_expr"    # mean expression of L over all cells of bin b, patient p
WS_LIGAND_AVAIL_SECONDARY <- "detect_frac"  # fraction of L-positive cells, always reported alongside

## -- 4. doublet  (§5.3 item 1, revised 2026-10-05) ----------------------------------------------
WS_DBL_MAIN_ARM       <- "scdbl_score_covariate"             # scDblFinder continuous score in Z
WS_DBL_SENSITIVITY    <- "union_hard_remove"                 # only on combos that get reported
WS_DBL_PROVIDER_CHECK <- c("scDblFinder", "DoubletFinder")   # H-C1 keeps BOTH callers
WS_DBL_RATE_PER_CELL  <- 8e-6                                # expected rate = n_raw * 8e-6 (pre-QC)
WS_DBL_SCORES_TABLE   <- "ws_doublet_scores_percell.csv.gz"  # per library, keyed by cell barcode
WS_DBL_MHCII_SEPARATE <- TRUE                                # MHC-II ligand combos reported apart

## -- 5. stage-1 covariates Z  (§5.4 step 4) ----------------------------------------------------
WS_Z_COVARIATES         <- c("CDR", "bmm_fine", "phase")    # + the scDblFinder score; NO malignancy
WS_Z_DROP_ZERO_VARIANCE <- TRUE   # a column constant within patient x bin duplicates the intercept
                                  # and makes Z'Z singular -- drop it for that fit and record it

## -- 6. GATE_W1.1 positive control  (§5.6, frozen 2026-10-02; rows 2 and 4 filled 2026-10-12, T3c-0) --
WS_W11_LIGANDS  <- c(IFNG = "primary", TNF = "auxiliary", TGFB1 = "auxiliary")
WS_W11_GENESETS <- c(IFNG  = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     TNF   = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     TGFB1 = "HALLMARK_TGF_BETA_SIGNALING")
WS_W11_RECEIVERS   <- c("LMPP_GMP", "Mono_DC")   # either one recovering = the pair recovered
WS_W11_MIN_PATIENTS <- 10L                        # combos below this are not run
WS_W11_COMBOS <- data.frame(
  row      = 1:6,
  ligand   = c("IFNG",    "IFNG",     "TNF",      "TNF",       "TGFB1",  "TGFB1"),
  sender   = c("T_NK",    "T_NK",     "Mono_DC",  "B_Plasma",  "T_NK",   "T_NK"),
  receiver = c("Mono_DC", "LMPP_GMP", "LMPP_GMP", "Mono_DC",   "Mono_DC","LMPP_GMP"),
  receptor = c("IFNGR2",  "IFNGR1",   "TNFRSF1A", "TNFRSF1B",  "TGFBR2", "TGFBR2"),
  patients           = c(33L, 18L, 19L, 32L, 28L, 15L),   # receptor-only: the §5.6 column (§0A口径)
  patients_intersect = c(28L, 17L, 18L, 28L, 23L, 15L),   # sender usable AND receptor usable, same patient
  stringsAsFactors = FALSE)
# Rows 2 and 4 filled 2026-10-12 by T3c-0 (job 3644488; provenance results/tables/ws/ws_t3c_preflight_fill.csv).
# Two decisions, both taken on 2026-10-12 BEFORE any T3c result exists (§15.1):
#   (1) `patients` keeps the receptor-only definition because that is what the four numbers already
#       frozen in §5.6 on 2026-10-02 are (verified: receptor-only reproduces §0A 6/6, intersection 2/6).
#       `patients_intersect` is frozen alongside as the number of patients that actually enter the
#       test, and WS_W11_MIN_PATIENTS is applied to IT (ws_w11_runnable). All six pass either way,
#       so the choice changes reported numbers, not which combos run.
#   (2) Row 4 sender = B_Plasma by the frozen rule (largest cross-patient IQR after excluding
#       Mono_DC: 0.0886). T_NK was 3% behind (0.0858) and is the biologically likelier TNF source; it
#       is recorded here as the runner-up and NOT substituted, because changing the selection rule
#       after seeing the candidates is exactly what the freeze forbids. If the user later overrides
#       to T_NK it must be a new dated constant with the reason, not an edit of this line.
WS_W11_ROW4_RUNNER_UP <- c(sender = "T_NK", IQR = 0.0858, chosen_IQR = 0.0886)
WS_W11_ELIGIBLE_SAMPLES <- 42L   # GSE185381 AML Diagnosis BM viability_only|none; settles §15.2 item 12

## -- 7. recovery criteria  (§5.6 item 4) --------------------------------------------------------
WS_W11_P_ONESIDED <- 0.05    # patient-level permutation of A[L,.], one-sided
WS_W11_FAKE_PCTL  <- 0.95    # the gene-set z must exceed this percentile of the fake-ligand z's
# fake ligands: drawn from the SAME sender compartment; not in CellChatDB; not in the tested gene set
WS_FAKE_N           <- 200L
WS_FAKE_TOL         <- 0.20  # detection rate AND inter-patient IQR within 20% of the real ligand
WS_FAKE_TOL_RELAXED <- 0.50  # used when fewer than WS_FAKE_MIN_POOL candidates pass at 0.20
WS_FAKE_MIN_POOL    <- 200L
WS_FAKE_SEED        <- SEED  # from config_paths.R

## -- 8. result tiers and the mandatory diagnostic  (§5.6 items 8-9) -----------------------------
WS_W11_TIER <- c(
  pass        = "IFNG recovered AND at least one other pair recovered, under BOTH doublet arms -> T5, T6",
  conditional = "TNF and TGFB1 recovered, IFNG not, under both arms -> STOP, report the diagnostic, user decides",
  fail        = "anything else -> check the pipeline, do not proceed",
  one_arm     = "recovered under only one doublet arm -> STOP and find out why")
WS_W11_DIAGNOSTIC <- "per patient: Spearman(mean gene-set score in receivers, A[L,p]); reported with every result, never scored"

## -- 9. NOT FROZEN -- do not treat as decided -------------------------------------------------
WS_W13_FDR_SUGGESTED        <- 0.20            # GATE_W1.3 empirical FDR; frozen only after T6 numbers (§15.2 item 1)
WS_LIGAND_AVAIL_SPLIT_PROPOSED <- "random_half" # A[L,p,b] from one random half of the bin, stage 1 on
                                                # the other, because A[L,p,b] ≡ the exclusion fraction
                                                # (Spearman 0.995). PROPOSED, never approved; must be
                                                # decided before T6.

## -- guards ----------------------------------------------------------------------------------
ws_assert_frozen <- function() {
  na <- WS_W11_COMBOS[is.na(WS_W11_COMBOS$sender) | is.na(WS_W11_COMBOS$receptor) |
                      is.na(WS_W11_COMBOS$patients) | is.na(WS_W11_COMBOS$patients_intersect), "row"]
  if (length(na)) stop("config_ws.R: GATE_W1.1 row(s) ", paste(na, collapse = ", "),
                       " still NA -- T3c-0 is not complete; T3c must not run")
  invisible(TRUE)
}
ws_w11_runnable <- function() {
  ws_assert_frozen()
  # the run/no-run rule is evaluated on the INTERSECTION -- a combo can only be tested in patients
  # where both the sender bin and the receptor are usable
  WS_W11_COMBOS[WS_W11_COMBOS$patients_intersect >= WS_W11_MIN_PATIENTS, ]
}
