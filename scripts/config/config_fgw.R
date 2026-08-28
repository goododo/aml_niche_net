# config_fgw.R ----
# Stage 07_fgw config (Phase 6): FGW alignment + barycenters.
# Sources config_ccc.R for the shared CCC data model (CCC_NODES) and config_distance.R for DIR_DISTANCE
# (the C matrices it consumes). Strictly downstream -> acyclic dependency, not sideways coupling.
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))        # CANONICAL_TIMEPOINTS (FGW_BARY_GROUPS derives from it)
source(here::here("scripts", "config", "config_ccc.R"))       # CCC_NODES, CCC_NODE_FEATURES
source(here::here("scripts", "config", "config_distance.R"))  # DIR_DISTANCE (edge_distance.csv / edge_qc.csv)

## -- graph representation (LOCKED, probe-informed) ----
# [locked] DIRECTED 7x7 C (symmetric=False in POT); NO 14-node sender/receiver split. POT 0.9.7 detects
#          asymmetry natively (all GW/FGW fns have symmetric=None); matches scACCorDiON directed CCC-OT.
# [locked] FIXED 7-node vocabulary; a bin absent in a sample keeps its node with eps mass (near 0) ->
#          balanced FGW ignores it -> same effect as unbalanced, and barycenters are balanced-only in POT.
#          Unbalanced FGW (fused_unbalanced_gromov_wasserstein) reserved as a pairwise robustness arm.
FGW_NODES       <- CCC_NODES               # 7, order authoritative for all matrices
FGW_EPS_MASS    <- 1e-6                     # mass assigned to an absent node (near-0; keeps fixed vocab)

## -- node mass p (LOCKED with sensitivity switch) ----
# [locked] default proportional to n_cells (down-weights tiny noisy bins; blueprint default).
#          "uniform" (present nodes equal) = composition-free sensitivity arm (R8). Residual composition
#          confound is ALSO controlled downstream by the blast-fraction covariate (R5).
FGW_MASS_MODE   <- "ncells"                # "ncells" | "uniform"

## -- FGW parameters ----
FGW_ALPHA_MAIN  <- 0.5                      # main structure/feature trade-off
FGW_ALPHA_SWEEP <- c(0, 0.25, 0.5, 0.75, 1)# 0 = pure feature OT (Wasserstein), 1 = pure GW (topology)
FGW_LOSS        <- "square_loss"
FGW_MAX_ITER    <- 1000L
FGW_SEED        <- SEED                     # 491638 (inherited); for barycenter init / any stochastic step

## -- node feature assembly (from 03 ccc_node_features.csv) ----
# Features fed to the FGW feature term, in this fixed column order. Healthy frac_malignant forced to 0
# here (design: healthy has no malignancy; suppresses inferCNV false positives seen in 03). Cross-sample
# feature scaling applied at assembly (see 01_build_fgw_inputs.R).
FGW_FEATURES    <- CCC_NODE_FEATURES        # c("frac_malignant","mean_stemness","n_cells")
# HOW FEATURES ARE PUT ON A COMMON SCALE. This is a RECORDED CONFLICT between the design and the
# implementation, and it is left as a switch rather than settled silently.
#
#   BLUEPRINT_v1.1_PATCH.md M8, "两条不可省的实现纪律" #1:
#     "全部特征在样本内取 rank-percentile，不用归一化表达值。否则会把平台效应从边权重
#      (已由 D2 的 rank 距离压制) 重新灌进特征项，等于绕过 D2。"
#
#   01_build_fgw_inputs.R header, written later:
#     "MUST be global (per-sample scaling erases between-sample signal)"
#
# Both are right about something. The patch is right that a global z-score lets platform-level
# expression offsets into the feature term, which the rank distance C was specifically designed to
# keep out. The code is right that with only 7 nodes per sample, a within-sample rank-percentile
# hands EVERY sample the same value set {1/7, ..., 1} and can therefore represent only WHICH node
# ranks highest -- never that one sample sits higher overall. mean_stemness_normal separates AML
# from healthy by a level shift, so ranking within sample removes that finding by construction, not
# by evidence.
#
# They are different questions, so run both and report both:
#   "global_z"           z-score across all node x sample rows (level-sensitive)
#   "within_sample_rank" frank(x)/n_nodes inside each sample   (the patch, platform-invariant)
#
# ---------------------------------------------------------------------------------------------
# DECIDED 2026-08-28: global_z is primary, within_sample_rank is the sensitivity arm. Report both.
#
# RETRACTED: the 2026-08-26 table that used to stand here. It compared the two scalings by the
# p_strat OF THE HYPOTHESIS UNDER TEST and adopted whichever was smaller. That is choosing the
# preprocessing that maximises your own finding's significance, and the numbers in it were wrong
# besides -- the within_sample_rank column was computed with n_cells silently truncated to a {0,1}
# one-hot (integer sub-assignment under a grouped `:=`, fixed in 01_build_fgw_inputs.R on 08-28).
# Refitting on repaired input moved 3 of its 5 entries (alpha 0: 0.0001->0.0002, 0.50: 0.0002->0.0018,
# 0.75: 0.0200->0.1215) and reversed the alpha=0.75 verdict. alpha=0.25 and alpha=1 were unchanged.
#
# The decision is now made on two grounds, neither of which is a p-value of anything under test.
#
# 1. MECHANISM. within_sample_rank cannot see a whole class of effect, as an identity rather than an
#    estimate: adding a constant to a sample's 7 bins leaves frank() bit-identical, so its power
#    against uniform level shifts is exactly nominal (0.05) at every effect size. Measured on the
#    138 samples: the per-sample bin-mean equals 1/7 for all 123 scaled features, spread <= 8.9e-16.
#    100% of between-sample level is removed, offset and biology alike. The premise that the level
#    was mostly platform offset is false as measured: dataset explains a median 32.6% of the level
#    variance (IQR 0.23-0.44), and the level still separates AML from healthy -- median |AUC-0.5|
#    0.079 across the 123 features, with 15 of them at >= 0.20. The dataset control is already in the
#    regression as a fixed effect; the rank transform pays for it a second time, with the biology.
#
# 2. KNOWN-SIGNAL RECOVERY, on the 114-feature panel screen at n_perm=1e6 (Discovery, per-family BH).
#    This is a screen over candidate features, not the primary hypothesis, and it was run at a
#    permutation count that resolves its own threshold:
#
#      family  n    global_z                    within_sample_rank
#      pg      42   2 pass, +13.6 SE margin     0 pass, nearest miss -79.7 SE
#      cs      39   3 pass,  +4.5 SE            0 pass,             -124.6 SE
#      pt       3   2 pass, +88.7 SE            0 pass,             -349.2 SE
#      st      12   0 pass                      0 pass
#      mt      18   0 pass                      1 pass, +1.4 SE  (still a coin flip at 1e6)
#                   ---------                   ---------
#                   7 of 114                    1 of 114, and that one unresolved
#
#    At n_perm=10000 the same screen gave 7 passes run as --subsets panels and 4 run as --subsets all
#    -- same data, same seed, same scaling, beta bit-identical (max |dbeta| = 0.0), only the position
#    in the shared permutation stream differing. It was that fragile because pg and cs cleared BH by
#    just +0.39 and +0.59 MC SE. Both problems -- the shared RNG (fixed, tests T1-T3) and the
#    permutation count (raised 100x) -- had to be fixed before this table meant anything.
#    See scripts/08_scoring/PREREGISTRATION_panel_screen.md.
#
# WHAT global_z COSTS, and what must therefore be watched: the cost matrix is divided by its own
# per-sample maximum (`M = M/(M.max()+1e-9)`), and under global_z that maximum spans 49x across the
# 138 samples (2.32 to 113.71, CV 1.15, 138 distinct values) versus 1.7x under the rank (CV 0.13).
# HDS is then each sample scaled by its own worst node-pair distance, which is not obviously
# comparable across samples. Replacing it with a cohort constant is an OPEN item in
# scripts/DECISIONS_pending.md; until it is settled the rank arm is kept and reported, not dropped.
FGW_FEATURE_SCALE <- "global_z"

FGW_SCALE_FEATURES <- TRUE                  # z-score each feature across nodes-x-samples before FGW
FGW_ZERO_HEALTHY_MAL <- TRUE                # set frac_malignant = 0 for Disease_state/timepoint == Healthy

## -- CANDIDATE node features (emitted, NOT used by FGW) ----
# Carried into fgw_nodes_long.csv, z-scored the same way, so 08_scoring/07_feature_decomposition.py can
# test them. They are deliberately NOT in FGW_FEATURES: a feature earns its way in by surviving the
# decomposition first. Adding one here changes no FGW result -- only what 07 is able to ask about.
#
# WHY THESE SIX. The alpha sweep put the whole H2 signal in the FEATURE term (alpha=1 pure topology:
# within-dataset p=0.966), so "which feature" is now the question the project turns on. Each candidate
# attacks a specific alternative explanation for the surviving stemness signal:
#   frac_malignant_vg        transcriptional malignancy (van Galen axis), NOT zeroed for healthy ->
#                            the non-circular counterpart to frac_malignant. HSC_MPP only (the only
#                            compartment where the axis passed held-out AUC); elsewhere NA -> neutral.
#   mean_cnv_burden          continuous CNV signal instead of the thresholded call
#   mean_stemness_malignant  is the stemness signal just "blasts are stem-like"? (near-circular)
#   mean_stemness_normal     or is it carried by NON-malignant cells? (a microenvironment result)
#   mean_cytotrace_normal    same split under an independent potency estimate that never saw LSC17
#   mean_cytotrace_malignant
# The stemness_normal vs stemness_malignant contrast is the decisive one.
FGW_CORE_CANDIDATES <- c("frac_malignant_vg", "mean_cnv_burden",
                         "mean_stemness_normal", "mean_stemness_malignant",
                         "mean_cytotrace_normal", "mean_cytotrace_malignant")

## -- PANEL candidates: the scores 05_ccc/03 aggregates from CCC_PANELS ----
# 38 per-cell scores x {all, non-malignant, malignant} = 114 columns the pipeline already computes
# and never tested. Derived from CCC_PANELS rather than re-listed, so the two cannot drift.
#
# THESE ARE EXPLORATORY AND THE STATISTICS MUST SAY SO. Everything above this line tests a
# hypothesis fixed before the data were seen. A sweep over 114 features does not, so two rules bind
# from here on:
#   1. Correction is BH-FDR WITHIN each panel family (st / pg / cs / mt / pt), not across all 114.
#      The families are five separate questions; pooling them would spend the correction on
#      unrelated hypotheses and bury the stemness robustness arm under 14 PROGENy pathways.
#   2. Screening happens on the DISCOVERY arm only. A dataset-level 70/30 split has existed since
#      2026-08-04 in 01_preprocess/02_sample_split.csv and no analysis has respected it; that was
#      tolerable while the feature set was fixed a priori, and stops being tolerable here.
#      Only 23 healthy samples exist cohort-wide, so the healthy controls are SHARED and only the
#      AML side is genuinely held out. Say that in the Methods; do not let it pass silently.
FGW_PANEL_CANDIDATES <- unlist(lapply(names(CCC_PANELS), function(pk) {
  base <- paste0(pk, "_", gsub("[^A-Za-z0-9_]", "_", CCC_PANELS[[pk]]$cols))
  c(base, paste0(base, "_normal"), paste0(base, "_malignant"))
}), use.names = FALSE)

# Which family each candidate belongs to -- 08_scoring/07 corrects within these.
FGW_CANDIDATE_FAMILY <- c(
  setNames(rep("core", length(FGW_CORE_CANDIDATES)), FGW_CORE_CANDIDATES),
  setNames(sub("^([a-z]{2})_.*$", "\\1", FGW_PANEL_CANDIDATES), FGW_PANEL_CANDIDATES))

# Imputation ceiling. Above this a column is mostly the cohort mean wearing a feature's name.
# In-use features exceeding it are a hard stop in 01; candidates are flagged and kept, because the
# decomposition should still be able to say "not measurable here" out loud.
FGW_MAX_IMPUTED <- 0.50

FGW_CANDIDATE_FEATURES <- c(FGW_CORE_CANDIDATES, FGW_PANEL_CANDIDATES)
stopifnot(!any(FGW_CANDIDATE_FEATURES %in% FGW_FEATURES))   # a candidate is by definition not yet in use
stopifnot(!anyDuplicated(FGW_CANDIDATE_FEATURES))
stopifnot(all(FGW_CANDIDATE_FEATURES %in% names(FGW_CANDIDATE_FAMILY)))

## -- barycenter grouping ----
# Which sample groups get a consensus barycenter (built on the fixed 7-node vocab).
# DERIVED from CANONICAL_TIMEPOINTS rather than re-listed, so the v2 vocabulary change (MRD and
# Post_treatment removed; Refractory / On_treatment / Post_induction / Post_consolidation /
# Post_transplant / Post_treatment_unspecified added) cannot silently shrink B_AML. Spelling out
# c("Diagnosis","MRD","Post_treatment",...) is exactly how that would have happened: every
# post-treatment sample would have dropped out of the AML barycenter with no error.
FGW_BARY_GROUPS <- list(
  healthy = list(timepoint = "Healthy"),                                   # B_healthy
  aml     = list(timepoint = setdiff(CANONICAL_TIMEPOINTS, c("Healthy", "Unknown")))  # B_AML
)
stopifnot(length(FGW_BARY_GROUPS$aml$timepoint) == length(CANONICAL_TIMEPOINTS) - 2L)
FGW_BARY_EXCLUDE_SPARSE <- TRUE             # drop sparse_flag==TRUE samples from barycenter construction (edge_qc.csv)

## -- IO ----
# R (01) writes tidy long CSVs (NO R<->Python bridge lib needed); Python (02) reads + pivots to arrays.
DIR_FGW <- file.path(FAST_DIR, "results/tables/07_fgw")   # fgw_{nodes,edges}_long.csv + index; barycenters + patient scores