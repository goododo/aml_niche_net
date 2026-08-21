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
FGW_CANDIDATE_FEATURES <- c("frac_malignant_vg", "mean_cnv_burden",
                            "mean_stemness_normal", "mean_stemness_malignant",
                            "mean_cytotrace_normal", "mean_cytotrace_malignant")
stopifnot(!any(FGW_CANDIDATE_FEATURES %in% FGW_FEATURES))   # a candidate is by definition not yet in use

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