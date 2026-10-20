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
#     FROZEN     -- section 10 (added 2026-10-16). GATE_B0 of HANDOFF_within_sample_v2.0.md §8.3,
#                   new dated constants written before any GATE_B0 result; sections 0-8 untouched.
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

## -- 10. GATE_B0 footprint positive control  (v2.0 §8.3; FROZEN 2026-10-16, before any GATE_B0 result) --
# Question: across AML patients, does T/NK ligand expression predict the matching CytoSig footprint in
# leukaemic progenitor cells, beyond blast proxy, ambient RNA, doublets and age group?
# Changes against the 2026-10-15 v2.0 text, all decided before any GATE_B0 number exists:
#   - receivers = progenitor bins, not Mono_DC (Petti 2026-10-15: monocytes/DC are often clone members);
#   - TNF sender = T_NK, the runner-up recorded in section 6; TNF from Mono_DC is clone-internal and is
#     reported only (WS_B0_CLONE_INTERNAL);
#   - GSE239721 cannot replicate (median 1 T_NK cell per sample); replication is reported, never scored;
#   - pass needs IFNG plus one more pair, the section 8 rule; two pairs without IFNG = conditional.
WS_B0_DATE            <- "2026-10-16"
WS_B0_GATE_DATASET    <- "GSE185381"
WS_B0_REPORT_DATASETS <- c("GSE289435", "GSE227903")     # direction only
WS_B0_RECEIVER_BINS   <- c("HSC_MPP", "LMPP_GMP")         # pooled; high_error cells dropped as in section 0
WS_B0_MIN_RECEIVER    <- 50L
WS_B0_SENDER_BIN      <- "T_NK"
WS_B0_MIN_SENDER      <- 30L
WS_B0_PAIRS <- data.frame(
  ligand    = c("IFNG", "TNF", "TGFB1"),
  cytosig   = c("IFNG", "TNFA", "TGFB1"),
  progeny   = c("JAK-STAT", "TNFa", "TGFb"),
  sender    = c("T_NK", "T_NK", "T_NK"),
  receptors = c("IFNGR1;IFNGR2", "TNFRSF1A;TNFRSF1B", "TGFBR1;TGFBR2;TGFBR3"),
  stringsAsFactors = FALSE)
WS_B0_CLONE_INTERNAL <- c(ligand = "TNF", sender = "Mono_DC")   # reported, labelled clone-internal, never scored
WS_B0_SCOPE <- "only ligands with a footprint signature (CytoSig 43 cytokines, PROGENy 14 pathways); the 2,196 detectable pairs are not the denominator"
# readout: CytoSig decides, PROGENy is reported
WS_B0_CYTOSIG_FILE   <- file.path(REF_DIR, "cytosig", "signature.centroid")
WS_B0_CYTOSIG_SHA256 <- "60a3152afe429a7dbb41fa0597ab9730141ce858af16a35f7063e64bfc239bde"
WS_B0_DECIDES        <- "CytoSig"
WS_B0_PROGENY_TOP    <- 500L
WS_B0_GENE_FILTER    <- "receiver pseudobulk CPM > 1 in at least half of the dataset's eligible samples"
WS_B0_SCORE          <- "decoupleR run_ulm score on receiver pseudobulk log2(CPM+1), each gene centred across the dataset's eligible samples; ligand and receptor genes removed from the signature"
WS_B0_IEG <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5", "ATF3",
               "DUSP1", "ZFP36", "NR4A1", "KLF2", "KLF4", "KLF6", "HSPA1A", "HSPA1B", "HSPA8", "HSPH1",
               "DNAJB1", "PPP1R15A", "SOCS3", "CEBPD")   # reported: signatures minus immediate-early/heat-shock genes
# exposure and model
WS_B0_LIGAND      <- "log1p(sender pseudobulk counts per 10k)"   # the section 3 per-cell mean, on a log scale
WS_B0_LIGAND_SENS <- "the same times the sender's fraction of binned cells"   # follows composition; reported
WS_B0_COVARIATES  <- c("prog_frac", "ambient_z", "dbl_mean", "paediatric")
#   prog_frac  = receiver-bin cells / binned cells (blast proxy; clinical blast% exists for 20 of 42)
#   ambient_z  = mean of the within-dataset z of (a) HBB+HBA1+HBA2 UMI share in non-erythroid binned
#                cells and (b) MPO+ELANE+AZU1+PRTN3 UMI share in T_NK + B_Plasma cells
#   dbl_mean   = mean scDblFinder score over the sample's cells (section 4 main arm)
#   paediatric = GSE185381 sample IDs not of the form AML + 4 digits: the 21 COG-risk samples plus AML003,
#                which has no clinical record but shares library P08 with three of them (22 vs 20, v2.0 §5.2);
#                constant elsewhere, so dropped there as zero-variance
#   receiver cell number is neither a covariate nor a weight: entry needs >= 50 cells; it moves precision, not bias
WS_B0_PERM        <- "Freedman-Lane"   # permute residuals of the covariate-only model, refit, ligand t
WS_B0_N_PERM      <- 10000L
WS_B0_P_ONESIDED  <- 0.05
WS_B0_FAKE_PCTL   <- 0.95              # the ligand t must exceed this percentile of the fake-ligand t's
WS_B0_FAKE_TOPSIG <- 200L              # fakes exclude the top 200 genes of the tested CytoSig signature
WS_B0_FAKE_RULE   <- "from the pair's own sender bin; detection rate and inter-patient IQR within WS_FAKE_TOL (WS_FAKE_TOL_RELAXED if the pool is short); not in CellChatDB; same model, same statistic"
# decision: the main doublet arm decides; the union-removal arm is reported, and a contradiction stops
WS_B0_TIER <- c(
  pass        = "IFNG passes AND at least one other pair passes (main arm) -> WP-B continues as v2.0's explanatory layer, never as the main line",
  conditional = "TNF and TGFB1 pass, IFNG does not -> STOP, report, user decides",
  fail        = "anything else -> ligand-receptor and footprint evidence are not used in this data",
  one_arm     = "the main arm passes but the union-removal arm does not -> STOP and find out why")
WS_B0_MIN_B  <- 30L
WS_B0_B_TIER <- "reported: add the B_Plasma footprint (same signature) as a covariate in samples with >= 30 B_Plasma cells; ligand still p < 0.05 -> 'leukaemia-specific', otherwise 'marrow-wide'"
WS_B0_PRIOR  <- "disclosed before running: T3c main-effect diagnostic (HALLMARK sets, n 10-18, no covariates), LMPP_GMP receivers: IFNG rho -0.08, TNF from Mono_DC 0.53, TGFB1 -0.19"

## -- 11. fine labels excluded from ws compartments  (added 2026-10-16, after the CCC audit) -------------
# BoneMarrowMap "Early Lymphoid" (CLP/MLP) is binned as T_NK by bmm_bin_map.tsv but behaves like a
# progenitor; in AML samples it is largely blasts (GSE185381, GATE_B0's 32 patients: 15% of T_NK-bin
# cells are author-labelled HSC/MPP/GMP, 8 patients >= 25%, AML009 92%). The CCC line has excluded it
# from node features since 2026-08-21 (config_ccc.R, CCC_EXCLUDE_FINE); T3c and GATE_B0 ran WITHOUT this
# exclusion (recorded in GATES.md and reports/2026-10_ccc_audit/). Every ws analysis written after this
# date drops these cells from T_NK (and from any "non-malignant" compartment) by bmm_fine / bmm_broad.
source(here::here("scripts", "config", "config_ccc.R"))   # constants only; sources config_paths.R again (idempotent)
WS_EXCLUDE_FINE      <- CCC_EXCLUDE_FINE                   # c("Early Lymphoid"), one definition for both lines
WS_EXCLUDE_FINE_DATE <- "2026-10-16"

## -- 12. S4 pseudobulk rules  (v2.0 §5.3 S4; FROZEN 2026-10-17 after the S3 counts, user decision) ----------
# Four cell groups instead of eight, a 20-cell floor instead of 30, concordance required: with eight groups and
# 30 cells only CD14 monocytes survived the 70% rule among the 20 GSE185381 adults (results/tables/ws/ws_s3_*).
# The statistical unit stays the patient (P3): every per-cell quantity is summarised within the patient before
# any test. Single-cell resolution enters as (a) subtype-level pseudobulk with fixed weights, (b) per-cell imprint
# scores against the healthy same-subtype range, summarised per patient (S5b), (c) DIALOGUE (S7).
WS_S4_DATE          <- "2026-10-17"
WS_S4_MAP_FILE      <- "s4_cell_groups.tsv"      # bmm_fine -> group (eight) -> group4 (four); EXCLUDED never counts
WS_S4_GROUPS        <- c("T", "NK", "B", "Mono")  # T = CD4 + CD8; Mono = CD14 + CD16; cDC and pDC dropped (S3)
WS_S4_MIN_CELLS     <- 20L                        # per sample and group4 (was 30 in the v2.0 text)
WS_S4_RETAIN_FRAC   <- 0.7                        # a group is kept if >= 70% of the study's samples reach the floor
WS_S4_REQUIRE_CONCORDANT <- TRUE                  # anno_source == "concordant" (marker typing agrees with the projection)
WS_S4_DROP_HIGH_ERROR    <- TRUE                  # as section 0
WS_S4_DEFINE     <- "GSE185381"                   # imprint defined here (adults; paediatric as the second cohort)
WS_S4_REPLICATE  <- "GSE289435"                   # structural replication
WS_S4_DESCRIBE   <- c("GSE227903", "GSE116256")   # too few lymphoid cells per sample (S3); description only
WS_S4_AMBIENT_COVARIATES <- c("hb_share_nonery", "gran_share_lymph")   # section 16.1 item 13: raw shares, not decontX
WS_S4_SENSITIVITY_COUNTS <- "decontX"             # decontaminated pseudobulk kept for sensitivity analyses only
WS_S4_SEX_GENES  <- c(female = "XIST", male = "RPS4Y1")   # inferred when the metadata has no sex

## -- guards ----------------------------------------------------------------------------------
ws_assert_frozen <- function() {
  na <- WS_W11_COMBOS[is.na(WS_W11_COMBOS$sender) | is.na(WS_W11_COMBOS$receptor) |
                      is.na(WS_W11_COMBOS$patients) | is.na(WS_W11_COMBOS$patients_intersect), "row"]
  if (length(na)) stop("config_ws.R: GATE_W1.1 row(s) ", paste(na, collapse = ", "),
                       " still NA -- T3c-0 is not complete; T3c must not run")
  invisible(TRUE)
}
ws_assert_b0_frozen <- function() {
  # GATE_B0 runs only against the exact signature file named in section 10
  if (!file.exists(WS_B0_CYTOSIG_FILE)) stop("config_ws.R: CytoSig file missing: ", WS_B0_CYTOSIG_FILE)
  got <- sub(" .*", "", system2("sha256sum", WS_B0_CYTOSIG_FILE, stdout = TRUE))
  if (!identical(got, WS_B0_CYTOSIG_SHA256)) stop("config_ws.R: CytoSig sha256 mismatch -- GATE_B0 must not run")
  invisible(TRUE)
}
ws_w11_runnable <- function() {
  ws_assert_frozen()
  # the run/no-run rule is evaluated on the INTERSECTION -- a combo can only be tested in patients
  # where both the sender bin and the receptor are usable
  WS_W11_COMBOS[WS_W11_COMBOS$patients_intersect >= WS_W11_MIN_PATIENTS, ]
}

## -- 13. DIALOGUE 1.0 patch (S7; PREREGISTRATION_S7_dialogue.md section 8, item 5) -------------------------------------
## DIALOGUE.identify.cell.types() builds the MCP x pair p matrix with generic.vector2mat(), whose rows are sorted
## alphabetically (B, Mono, NK, T), but its first elimination indexes R$cell.types in the original order (T, NK, B, Mono),
## so the wrong compartment is dropped. Verified 2026-10-19 on the S7 objects: R1/R2 MCP1-2 were assigned T+B instead
## of T+NK+B; the planted run's MCP2-3 were assigned B+Mono instead of T+NK. The fixed version below keeps the
## algorithm (drop the compartment with the fewest p < 0.1 links until every one links to at least half of the rest)
## and only aligns the row order. Called by 29 and 30 after library(DIALOGUE); it self-tests and then replaces the
## namespace function so that DIALOGUE1() (and hence DIALOGUE2's choice of pairs) uses it.
ws_dialogue_identify_cell_types <- function(R) {
  if (!is.null(R$param$MCP.cell.types)) return(R$param$MCP.cell.types)
  if (length(R$cell.types) == 2) { out <- lapply(as.numeric(R$emp.p), function(x) if (x < 0.1) R$cell.types else NULL); names(out) <- rownames(R$emp.p); return(out) }
  emp.p <- R$emp.p; emp.p2 <- R$emp.p
  colnames(emp.p2) <- paste(DIALOGUE:::get.strsplit(colnames(emp.p), "_", 2), DIALOGUE:::get.strsplit(colnames(emp.p), "_", 1), sep = "_")
  emp.p <- cbind(emp.p, emp.p2); out <- list()
  for (x in rownames(emp.p)) {
    p1 <- DIALOGUE:::generic.vector2mat(emp.p[x, ])[R$cell.types, R$cell.types]; diag(p1) <- 0.05   # rows and columns in R$cell.types order
    cts <- R$cell.types; lb <- min(rowSums(p1 < 0.1))
    while (lb < ceiling(length(cts) / 2)) { cts <- setdiff(cts, rownames(p1)[rowSums(p1 < 0.1) <= lb][1]); p1 <- p1[cts, cts, drop = FALSE]; lb <- min(rowSums(p1 < 0.1)) }
    out[[x]] <- cts
  }
  out
}
ws_patch_dialogue <- function() {
  stopifnot(requireNamespace("DIALOGUE", quietly = TRUE))
  toy <- list(cell.types = c("T", "NK", "B", "Mono"), param = list(),   # links T-NK and T-B only: the intended answer is T+NK+B
              emp.p = matrix(c(0.04, 0.05, 0.9, 0.2, 0.9, 0.9), nrow = 1, dimnames = list("MCP1", c("T_NK", "T_B", "T_Mono", "NK_B", "NK_Mono", "B_Mono"))))
  orig <- DIALOGUE:::DIALOGUE.identify.cell.types(toy)$MCP1; fixed <- ws_dialogue_identify_cell_types(toy)$MCP1
  stopifnot(identical(fixed, c("T", "NK", "B")))
  message(sprintf("ws_patch_dialogue: identify.cell.types on the toy case gives %s unpatched, %s patched; namespace function replaced",
                  paste(orig, collapse = "+"), paste(fixed, collapse = "+")))
  utils::assignInNamespace("DIALOGUE.identify.cell.types", ws_dialogue_identify_cell_types, ns = "DIALOGUE")
  invisible(TRUE)
}

## -- 14. S7 v2 (PREREGISTRATION_S7_dialogue.md section 9; FROZEN with tag prereg-s7-v2, 2026-10-20) ------------------------
WS_S7_COVAR <- list(R1 = c("cellQ", "hb_share_nonery", "gran_share_lymph", "prog_frac"),          # DIALOGUE gene-level covariates per run
                    R2 = c("cellQ", "hb_share_nonery", "gran_share_lymph", "prog_frac", "status"))
WS_S7_V2_NPCS <- 10L         # the only parameter changed from v1 (30 PCs); chosen on planted data only (section 8, item 7)
WS_S7_V2_TAG  <- "_v2"
WS_S7_REDUNDANT_R <- 0.8     # two MCPs whose step-I cell scores correlate above this in every compartment count as one program (reporting only)
