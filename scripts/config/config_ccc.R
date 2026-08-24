# config_ccc.R ----
# Stage 05_ccc config: per-sample CellChat CCC graphs on a FIXED 7-bin node vocab.
# Sourced AFTER config_paths.R (FAST_DIR/LARGE1_DIR/DIR_MALIGNANCY/SEED). All decisions LOCKED.
#   [D1]      niche = hematopoietic-immune + paracrine; structural stroma is B-layer / cohort-aux only.
#   [locked]  Main graph = 7 hierarchy bins; malignancy is a NODE FEATURE, not a node split.
#   [M5]      SCS = feature-modulated: LSC-primitive sender -> immune/paracrine receiver.
#   [M-keep]  Megakaryocyte stays a node; per-sample absence -> unbalanced FGW.
#   [R9]      per-patient graphs ONLY from L2-capable whole-MNC libs (sender+receiver present).
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
# NOTE: SEED inherited from config_paths.R (verified 491638L via recon_ccc_probe [0]).

## -- Upstream input locations (verified resolvable via recon probe [0]) ----
CCC_QC_OBJ_DIR    <- file.path(LARGE1_DIR, "02_seurat_objects/01_per_sample_qc")   # <ds>/<sample>.rds (v5, counts-only layer)
CCC_BMM_DIR       <- file.path(LARGE1_DIR, "02_seurat_objects/03_bmm_projected")   # <ds>/<sample>__bmm_percell.csv
# consensus per-cell: DIR_MALIGNANCY/<ds>/<sample>__consensus_percell.csv (from config_paths)
CCC_ROLE_MANIFEST <- file.path(FAST_DIR, "results/tables/01_preprocess/01_sample_role_manifest.csv")  # verified path

## -- Node vocabulary (FIXED = FGW barycenter node set) ----
CCC_NODES <- c("HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma")
CCC_NODE_FEATURES <- c("frac_malignant","mean_stemness","n_cells")
# WARNING: frac_malignant is UNRELIABLE for lymphoid receiver bins (T_NK/B_Plasma) -- inferCNV false
# positives (recon [3]: T_NK 36/194 "malignant"). Do NOT biologically interpret malignancy there.
CCC_STEMNESS_SIG <- "LSC17"   # signature column used for the mean_stemness node feature (swappable;
                              # stemness_percell has LSC17 + HSPC_core + vanGalen_*). RAW per-node means
                              # are written by 03; feature scaling/normalization for FGW is deferred to 07.
CCC_SENDER_LSC  <- c("HSC_MPP","LMPP_GMP")           # SCS senders (feature-modulated by malignant_frac)
CCC_RECV_IMMUNE <- c("T_NK","B_Plasma","Mono_DC")    # immune/paracrine receivers (incl. HSC-macrophage axis)

## -- CCC-eligibility gate [R9] (METADATA-DRIVEN; NOT composition-based) ----
# A lymphoid floor CANNOT gate: blast-dominated whole-MNC has ~0 lymphoid yet is valid (missing
# nodes biological), while a sorted CD34 library's missingness is technical. -> gate on protocol.
CCC_ELIGIBLE_REQUIRES_L2  <- TRUE
# EMPTY, was list(Chen2023 = "CD34"). That entry excluded every "<donor>_CD34" sublibrary on the
# theory that it was a CD34-enriched fraction. Per Chen 2023 Methods it is pool A = HSPC +
# myeloid, the other half of a two-library donor, and it carries 5 of these 7 nodes (HSC_MPP,
# LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte). Excluding it did not remove a biased library --
# it deleted half of every Chen2023 marrow and left a T/B-only graph. ingest_Chen2023.R now
# merges both pools into one donor-level Sample, so there is nothing left to exclude here; the
# composition bias travels as the curated covariate sorting="multi_fraction_FACS_recombined".
# Any v1 count that quoted "10 sorted_sublibrary" exclusions is void.
CCC_SORTED_SUBLIB_PATTERN <- list()   # dataset -> regex (case-insensitive) of sorted sublibs
# MUST track MIN_CELLS_SAMPLE in config_qc.R. These are two independent gates on the same
# quantity: while this one sat at 500 and MIN_CELLS_SAMPLE was lowered to 300, every sample in
# the 300-499 band was re-admitted at 01_preprocess only to be cut again here as
# "below_min_cells" -- i.e. the MIN_CELLS_SAMPLE change was a no-op end-to-end. Quality is
# enforced by FLAG_LOWCOMPLEXITY_DROP (median nFeature >= 500) plus the occupancy gate below
# (>= CCC_MIN_OCCUPIED_BINS nodes holding >= CCC_MIN_CELLS_PER_OCCUPIED_BIN cells), which is the
# constraint that actually matters for a 7-node graph.
CCC_MIN_TOTAL_CELLS       <- 300L                    # == MIN_CELLS_SAMPLE (config_qc.R)
# v1 result (stale, MIN_CELLS=500 + intersection doublets + NBM mislabelled as Diagnosis):
#   159 eligible / 220 projected (48 not_l2_capable + 10 sorted_sublibrary + 3 below_min_cells).
#   Recomputed in v2; do not quote the v1 numbers.

## -- Per-node presence threshold (WITHIN a sample's graph) ----
# Below this a node is 'missing' in that sample (NOT dropped) -> unbalanced FGW. Also = CellChat
# filterCommunication(min.cells=) floor (aligned: groups < min.cells are auto-excluded by CellChat).
## -- fine labels excluded from NODE FEATURES ----
# A projection label that does not mean what it says must not contribute to a node's features.
# `Early Lymphoid` is the live case (2026-08-21): bmm_bin_map.tsv already flags it AMBIGUOUS
# (CLP-like, upstream of BOTH T_NK and B_Plasma), and it behaves like a progenitor, not a lymphocyte --
# CytoTRACE2 potency 0.446 against 0.064-0.071 for real T/NK, LSC17 0.123 against ~0.013, the worst
# mapping_error in the node. It is 0.95% of the healthy non-malignant T_NK pool but 10.64% of the AML
# one (11x), which alone accounts for about half the T_NK stemness difference.
#
# EXCLUDED FROM FEATURES, NOT FROM THE GRAPH. CellChat has already attributed these cells' edges to
# T_NK; undoing that means changing bmm_bin_map.tsv's in_ccc_graph and re-running 05_ccc/02 for all 138
# samples, which also costs 8 samples their T_NK node (7 AML / 1 healthy) and widens the node-presence
# gap from +8.7pp to +10.4pp. Do that only if the finding survives this cheaper test. Note that node
# MASS follows n_cells, so masses do shift here even though the edges do not -- that is the intended
# direction (a node with fewer real T/NK cells should carry less mass), but it is a partial change,
# and the edge/feature inconsistency is real until the graph is rebuilt.
CCC_EXCLUDE_FINE <- c("Early Lymphoid")

CCC_MIN_CELLS_PER_NODE <- 10L

## -- Per-sample OCCUPANCY gate (the constraint that actually binds for a 7-node graph) ----
# A total-cell floor is the wrong instrument for CCC. 300 cells spread over 7 nodes can leave
# several nodes with a handful of cells each; CellChat then estimates each group's mean
# expression from almost nothing, so the graph exists but its edge weights are noise. What must
# be guaranteed is OCCUPANCY: enough nodes carrying enough cells to support a comparable
# topology. This can only be evaluated after BMM projection, which is why it lives here and not
# as a cell-count floor in 01_preprocess.
# Distinct from CCC_MIN_CELLS_PER_NODE (=10), which decides whether a node is PRESENT in a
# sample's graph (absent -> unbalanced FGW handles it); this decides whether the SAMPLE gets a
# graph at all. A sample failing it is excluded rather than contributing a 2-node graph whose
# topology distance is meaningless.
CCC_MIN_OCCUPIED_BINS          <- 3L
CCC_MIN_CELLS_PER_OCCUPIED_BIN <- 30L

## -- CellChat engine parameters (ALL verified against installed 2.2.0.9001 via probe [4]) ----
CCC_ENGINE_DEFAULT <- "cellchat"     # liana arm = 02b (schema-compatible; R3 two-method gate). liana 0.1.14 installed.
CCC_DB_SPECIES     <- "human"
# VERIFIED exact strings (probe [4] annotation table). Protein L-R only for the main run.
CCC_DB_CATEGORIES  <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")   # 1280+424+535 = 2239
CCC_DB_EXCLUDE     <- c("Non-protein Signaling")                                    # 994 metabolite/ion -> out
# [DECISION population.size=FALSE] intrinsic signaling; abundance (=blast burden) controlled DOWNSTREAM
# by composition-permutation null + rank distance + blast-fraction covariate (R3/R5). Verified default.
CCC_POPULATION_SIZE <- FALSE
CCC_COMMUN_TYPE     <- "triMean"     # conservative default
CCC_TRIM            <- 0.1           # only used if CCC_COMMUN_TYPE != "triMean"
# CRITICAL (probe [4]): computeCommunProb defaults distance.use=TRUE/contact.dependent=TRUE are for
# SPATIAL data. This is non-spatial scRNA -> force classic law-of-mass-action:
CCC_DISTANCE_USE    <- FALSE
CCC_NBOOT           <- 100L          # permutation count (CellChat default)
CCC_RAW_USE         <- TRUE          # use normalized (unsmoothed) data; projectData ABSENT in this build -> no smoothing
CCC_SUBSET_THRESH   <- 1             # subsetCommunication thresh: keep ALL edges + pval; filter significance downstream (R3)
CCC_WORKERS         <- 1L            # future multisession workers for nboot; 1 = sequential. Set = allocated cores in array.

## -- Output ----
DIR_CCC        <- file.path(FAST_DIR,   "results/tables/05_ccc")   # manifests + tensors (purgeable scratch)
CCC_TENSOR_DIR <- file.path(DIR_CCC,    "tensors")                 # <ds>/<sample>__ccc_cellchat.csv
CCC_GRAPH_DIR  <- file.path(LARGE1_DIR, "05_ccc_graphs")           # <ds>/<sample>__cellchat.rds (persistent)
# one row per directed LR edge:
CCC_TENSOR_COLS <- c("sample","dataset","timepoint","engine",
                     "sender_bin","receiver_bin","ligand","receptor",
                     "prob","pval","n_sender","n_receiver")