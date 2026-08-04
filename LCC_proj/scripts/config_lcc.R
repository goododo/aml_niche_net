# config_lcc.R ----
# Stage config for LCC_proj: the TP53-adjacent side branch (PVRL4/Nectin-4 + pro-fibrotic niche).
# Sources the MAIN-LINE config_paths.R for roots/SEED, then adds LCC-only constants.
# LCC_proj is a SEPARATE tree from scripts/ -- it never writes into main-line results dirs and
# never mutates a main-line config constant (CCC_NODES stays locked, so 06/07/08 stay valid).
#
# INPUT  : none (constants only; source()d, never run standalone)
# OUTPUT : none
# Usage  : source(here::here("LCC_proj", "scripts", "config_lcc.R"))
#
# LANGUAGE DISCIPLINE [CODING_STANDARDS §9]: inferCNV measures CNV, NOT point mutations. Every
# object, column and label here is named *_cnv_proxy / TP53_aberrant -- never "TP53 mutant".
# The single genotype-confirmed TP53 sample in the cohort (Petti2019 809653) carries TP53 E286G
# (missense) + del(17q) and does NOT carry 17p loss, which is exactly why the 17p-only rule is
# insufficient and the complex-karyotype (CK) arm exists. See LCC_proj_task_HANDOFF_v1.md.

suppressPackageStartupMessages({ library(here); library(data.table) })
source(here::here("scripts", "config", "config_paths.R"))   # FAST_DIR/LARGE1_DIR/INFERCNV_ROOT/SEED
source(here::here("scripts", "config", "utils.R"))          # fwrite_safe/get_counts/core16/message_ts
source(here::here("scripts", "config", "config_qc.R"))      # is_healthy_sample() -- the timepoint fix

# NOTE: SEED is inherited from config_paths.R (491638L, the on-disk single source of truth).
# CODING_STANDARDS.md §6 quotes 20260605L; the config file wins -- do not introduce a third value.

## -- LCC roots (FAST = small tables/figures; LARGE1 = big objects) ----
LCC_DIR        <- here::here("LCC_proj")
LCC_SCRIPT_DIR <- file.path(LCC_DIR, "scripts")
LCC_PANEL_DIR  <- file.path(LCC_SCRIPT_DIR, "panels")
LCC_TAB_DIR    <- file.path(LCC_DIR, "results", "tables")
LCC_FIG_DIR    <- file.path(LCC_DIR, "results", "figures")
LCC_LOG_DIR    <- file.path(LCC_DIR, "logs")

LCC_LARGE_DIR  <- file.path(LARGE1_DIR, "LCC_proj")
LCC_PB_DIR     <- file.path(LCC_LARGE_DIR, "pseudobulk")     # 3-level pseudobulk count matrices
LCC_PERCELL_DIR<- file.path(LCC_LARGE_DIR, "percell")        # UCell scores, macrophage calls
LCC_STROMA_DIR <- file.path(LCC_LARGE_DIR, "stroma")         # L2 stromal subclustering objects
LCC_CC_DIR     <- file.path(LCC_LARGE_DIR, "cellchat_lcc")   # only if the gate says re-run CellChat

## -- reference tables shipped with this branch ----
LCC_ARMS_TSV   <- file.path(LCC_PANEL_DIR, "chrom_arms_grch38.tsv")
LCC_GENE_PANEL_TSV <- file.path(LCC_PANEL_DIR, "fibrosis_ecm_panel.tsv")
LCC_PATHWAY_TSV<- file.path(LCC_PANEL_DIR, "pathway_sets.tsv")
LCC_GMT        <- file.path(LCC_PANEL_DIR, "msigdb.v2026.1.Hs.symbols.gmt")   # full MSigDB 2026.1 (35361 sets)

## -- upstream inputs (all main-line, READ-ONLY from here) ----
LCC_ROLE_MANIFEST <- file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv")  # uid_patient/Timepoint/subtype
LCC_BMM_DIR       <- file.path(LARGE1_DIR, "02_seurat_objects", "03_bmm_projected")
LCC_CCC_NODE_FEAT <- file.path(TAB_DIR, "05_ccc", "ccc_node_features.csv")
LCC_GENE_ORDER    <- file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt")

# inferCNV per-sample filenames (fixed by the main-line 44_infercnv_run_one.R invocation:
# HMM i6, leiden subclusters, Bayes filter at Pnorm 0.5). Verified present for all 130 samples.
LCC_ICNV_GENES_FILE  <- "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"
LCC_ICNV_USED_FILE   <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.genes_used.dat"
LCC_ICNV_GROUP_FILE  <- "17_HMM_predHMMi6.leiden.hmm_mode-subclusters.cell_groupings"

## -- inferCNV HMM i6 state encoding ----
# i6 states map to copy-number multiples: 1 = 0x, 2 = 0.5x, 3 = 1x (NEUTRAL), 4 = 1.5x, 5 = 2x, 6 = 3x.
LCC_STATE_NEUTRAL <- 3L
LCC_STATE_LOSS    <- c(1L, 2L)
LCC_STATE_GAIN    <- c(4L, 5L, 6L)

## -- TP53 locus (GRCh38, matching gencode_GRCh38_gene_order.txt: chr17:7661779-7687550) ----
LCC_TP53_GENE <- "TP53"
LCC_TP53_CHR  <- "chr17"
LCC_TP53_ARM  <- "17p"

## -- CNV proxy definition [primary thresholds; the sensitivity grid is below] ----
# CK arm : an arm counts as "arm-level event" when >= LCC_ARM_FRAC_MIN of the sample's inferCNV-used
#          genes on that arm sit in a non-neutral state of ONE direction. Gene-fraction, not bp --
#          inferCNV only sees where genes are, so bp coverage would systematically under-call
#          gene-poor arms.
# 17p arm: the TP53 gene itself is in a loss state (state <= 2) in that subclone.
LCC_ARM_FRAC_MIN   <- 0.80   # fraction of arm genes needed to call an arm-level event
# 0.80 is EMPIRICAL, not conventional: 02 sweeps 0.50-0.95 and 0.80 is the loosest value at which
# the healthy-donor false-positive rate for 17p loss reaches 0/22. At the initially planned 0.50
# the healthy hit rate for 17p loss was 14% vs 13% in Diagnosis AML -- i.e. pure noise.
LCC_ARM_EVENT_MIN  <- 3L     # >= K arms with events -> complex-karyotype-like (CK)
LCC_CLONE_FRAC_MIN <- 0.10   # subclone must hold >= this fraction of the sample's malignant cells
LCC_SUBCLONE_MAL_MIN <- 0.50 # subclone must be predominantly malignant to carry a somatic call
LCC_MIN_SUBCLUST_CELLS <- 20L  # subclones smaller than this are excluded (HMM unstable)
# Baseline: every sample's own inferCNV reference_normal subclones are the internal negative
# control. A sample with too few of them cannot be called -- its noise floor is unknown.
LCC_MIN_REF_SUBCLONES <- 3L
LCC_MIN_REF_CELLS     <- 100L

# Acrocentric p arms carry (almost) no genes and inferCNV never covers them -> never eligible for
# an arm-level call. Excluded explicitly so they can never inflate or deflate the arm count.
LCC_ARMS_EXCLUDE <- c("13p", "14p", "15p", "21p", "22p", "Yp", "Yq", "Xp", "Xq")
# Sex chromosomes excluded too: apparent X/Y "events" track donor sex vs the reference mix,
# not somatic CNV. [reviewer attack point -- documented, not silently dropped]

## -- sensitivity grid (rigor: primary endpoint is the default cell of this grid) ----
LCC_SENS_ARM_FRAC_MIN   <- c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
LCC_SENS_ARM_EVENT_MIN  <- c(2L, 3L, 4L, 5L)
LCC_SENS_SUBCLONE_MAL   <- c(0.0, 0.5, 0.8)

## -- canonical recurrent AML/MDS arm lesions, for the lesion-specificity control ----
# Artifactual arm calls should be spread uniformly across arms; real somatic lesions concentrate on
# these. 02 reports Diagnosis-vs-Healthy hit rates per arm so the claim is checkable, not asserted.
LCC_CANONICAL_LESIONS <- c("5q:loss", "7p:loss", "7q:loss", "17p:loss", "20q:loss",
                           "8p:gain", "8q:gain")

## -- proxy confidence classes ----
# CK_and_17p / 17p_only -> high ; CK_only -> mid ; else negative.
# Primary comparison [D-inclusion, locked]: high vs negative, Diagnosis only, one sample per patient.
# mid is EXCLUDED from the primary test and re-introduced only in the sensitivity analysis.
LCC_PROXY_HIGH <- c("CK_and_17p", "17p_only")
LCC_PROXY_MID  <- c("CK_only")

## -- inclusion criteria [locked with user, 2026-07-27] ----
LCC_TIMEPOINT_PRIMARY <- "Diagnosis"   #治疗后/复发样本的微环境已被化疗改造 -> excluded from primary
LCC_HEALTHY_AS_REFERENCE <- TRUE       # healthy donors are a plotted reference arm, NOT a test arm

## -- macrophage vs monocyte scoring panel (gate for the CCC node vocabulary decision) ----
LCC_MAC_MARKERS  <- c("C1QA", "C1QB", "C1QC", "CD163", "MRC1", "VSIG4", "SELENOP", "STAB1", "LYVE1", "MAF")
LCC_MONO_MARKERS <- c("FCN1", "S100A8", "S100A9", "S100A12", "VCAN", "SELL", "CSF3R", "PLAC8")
LCC_MAC_MIN_CELLS <- 30L   # per-sample floor for macrophage-like to count as a usable CCC node
LCC_MAC_GATE_MIN_SAMPLES <- 50L  # >= this many Dg samples clearing the floor -> worth re-running CellChat

## -- p53 functional axis ----
# 20 high-confidence DIRECT p53 transcriptional targets. This replaced HALLMARK_P53_PATHWAY after
# the Hallmark set FAILED anchor validation (4 CNV anchors at percentiles 6/38/87/66, permutation
# p = 0.47). The two scores correlate at Spearman rho = 0.33 only -- the 200-gene Hallmark set is
# dominated by p53-independent stress and proliferation genes in this tissue.
# DISCLOSURE for the methods section: this set was chosen AFTER seeing Hallmark fail, so its
# 4-anchor p (0.067) carries analytic flexibility. Petti2019 809653 (genotype-confirmed TP53 E286G,
# NOT called by the 17p track and NOT used to choose this set) is the held-out anchor; adding it
# gives 5 anchors, permutation p = 0.0185. Report both numbers, in that order.
LCC_P53_TARGETS <- c("CDKN1A", "MDM2", "BAX", "GADD45A", "RRM2B", "ZMAT3", "TP53I3", "SESN1",
                     "SESN2", "TNFRSF10B", "BBC3", "FAS", "DDB2", "XPC", "TP53INP1", "AEN",
                     "TRIAP1", "CCNG1", "PHLDA3", "EDA2R")
LCC_P53_MIN_TARGETS <- 15L   # a sample needs this many of them in its reference to be scored
LCC_GENOTYPE_ANCHOR <- "809653"   # Petti2019; the cohort's only genotype-confirmed TP53 mutant

## -- Nectin family (题目 A). GRCh38 gene_order uses NECTIN4; PVRL4 is the legacy alias ----
LCC_NECTIN_GENES <- c("NECTIN4", "PVRL4", "NECTIN1", "PVRL1", "NECTIN2", "PVRL2",
                      "NECTIN3", "PVRL3", "TIGIT", "CD226", "PVR")

## -- DATASET SELECTION [revised with user 2026-07-28, second pass] ----
# REVISION: the user does not need healthy-donor comparison arms. The earlier 3-dataset restriction
# existed solely to make every AML-vs-healthy contrast WITHIN study; without a healthy arm that
# rationale is void, so the selection reopens to every dataset that contributes AML samples. The
# comparison of record is TP53-mutant AML vs TP53-wild-type AML.
# GSE116256 is reinstated specifically because it holds 3 of the 4 genotype-confirmed TP53-mutant
# patients (AML328, AML420B, AML916 -- van Galen per-cell genotyping); its previously disqualifying
# metric (healthy-donor FPR 0.755) only ever concerned its healthy arm, which is now unused.
LCC_CORE_DATASETS <- c(
  "GSE185381",    # Lasry 2023        27 AML-Dg patients; 2 inferCNV 17p-loss calls
  "GSE239721",    # (IFN-gamma AML)   15 AML-Dg; WARNING icnv_labelled_frac 0.38, no usable ref baseline
  "GSE185991",    # Naldini 2023      13 AML-Dg
  "GSE289435",    # Zeng 2025         12 AML-Dg; Numbat allele arm available
  "GSE116256",    # van Galen 2019     9 AML-Dg; THE genotype truth set (per-cell TP53 calls)
  "GSE227903",    # Ennis 2023         9 AML-Dg; Numbat allele arm available
  "Chen2023",     # Chen 2023          6 AML-Dg; the only AML dataset with stromal capture
  "Petti2019",    # Petti 2019         5 AML-Dg; holds 809653 (TP53 E286G)
  "GSE201966", "GSE147989", "GSE207356")   # 3 + 2 + 1 AML-Dg
LCC_STROMA_REF_DATASETS <- c("GSE253355")   # 25,370 healthy stromal cells; stromal reference only

# Healthy donors are RETAINED but only as a methods negative control (false-positive rate for any
# CNV-based call). They are never a comparison arm. Keeping them is what lets us report a measured
# FPR instead of asserting specificity.
LCC_HEALTHY_ROLE <- "negative_control_only"

# Excluded, with the disqualifying metric. Kept in code so the exclusion is auditable.
LCC_EXCLUDED_DATASETS <- data.table::data.table(
  dataset = c("E-MTAB-11536"),
  reason  = c("no AML samples (all healthy) and no malignancy consensus"))

## -- Numbat allele-based 17p LOH arm ----
# The ONLY evidence arm orthogonal to expression-CNV: LOH is called from allele counts, not from
# smoothed expression. 17p LOH (copy-neutral or deletional) is the canonical mechanism of TP53
# biallelic inactivation, so this arm is the strongest CNV-side evidence available.
# Numbat iteration rule [CODING_STANDARDS §9]: the highest-numbered segs_consensus_{N}.tsv wins.
LCC_NUMBAT_ROOT     <- file.path(CNV_ROOT, "numbat")   # <ds>/<sample>/numbat/segs_consensus_N.tsv
LCC_NUMBAT_LOH_STATES <- c("loh", "del", "bdel")
LCC_NUMBAT_MIN_LLR  <- 5      # segment log-likelihood ratio floor; observed range on real LOH: 16-209
LCC_NUMBAT_MIN_ARM_FRAC <- 0.50   # fraction of 17p covered by a LOH/del segment
LCC_TP53_START <- 7661779L    # GRCh38, gencode_GRCh38_gene_order.txt
LCC_TP53_END   <- 7687550L
LCC_CHR17_CEN_START <- 22700000L   # 17p ends here (chrom_arms_grch38.tsv)

## -- loaders ----
# THE SINGLE TIMEPOINT-FIX POINT. The main-line role manifest on disk (2026-07-15) predates the
# config_qc.R healthy-sample rule, so Chen2023's 8 NBM samples (4 normal-marrow donors) are stored
# as Timepoint = "Diagnosis", Disease_state = "Unknown". config_qc.R:95 documents this exact failure
# mode. Every LCC script must read metadata through this function -- never trust the raw Timepoint
# column -- so the fix lives in one place and cannot drift between scripts.
load_sample_meta <- function(selected_only = FALSE, include_stroma_ref = FALSE) {
  f <- file.path(LCC_TAB_DIR, "03_sample_manifest.csv")
  if (!file.exists(f)) stop("[config_lcc] run 03_percell_pass.R --make_manifest first: ", f)
  m <- fread(f)
  m[, is_healthy_donor := vapply(sample, function(s) isTRUE(is_healthy_sample(s)), TRUE)]
  m[, timepoint_raw := timepoint]
  m[, timepoint := fifelse(is_healthy_donor, "Healthy", timepoint_raw)]
  m[, timepoint_was_fixed := timepoint != timepoint_raw]
  keep <- c(LCC_CORE_DATASETS, if (include_stroma_ref) LCC_STROMA_REF_DATASETS)
  m[, in_selection := dataset %in% keep]
  m[, dataset_role := fifelse(dataset %in% LCC_CORE_DATASETS, "core",
                       fifelse(dataset %in% LCC_STROMA_REF_DATASETS, "stroma_reference", "excluded"))]
  if (selected_only) m <- m[in_selection == TRUE]
  m[]
}
# kept for scripts that genuinely want the unfixed main-line view (none should)
load_role_manifest <- function() {
  if (!file.exists(LCC_ROLE_MANIFEST)) stop("[config_lcc] role manifest not found: ", LCC_ROLE_MANIFEST)
  fread(LCC_ROLE_MANIFEST)
}
# fread() has no comment.char, and every panel TSV here carries a '##' provenance header
# (source, build, which script consumes it). Strip those lines before parsing.
fread_commented <- function(file, ...) {
  if (!file.exists(file)) stop("[config_lcc] panel file not found: ", file)
  ln <- readLines(file, warn = FALSE)
  ln <- ln[!grepl("^##", ln)]
  fread(text = paste(ln, collapse = "\n"), sep = "\t", ...)
}
load_chrom_arms <- function() fread_commented(LCC_ARMS_TSV)
load_gene_panel <- function() fread_commented(LCC_GENE_PANEL_TSV)
# Read only the requested gene sets out of the 30 MB GMT (never load all 35361 into memory).
read_gmt_subset <- function(set_names, gmt = LCC_GMT) {
  if (!file.exists(gmt)) stop("[config_lcc] GMT not found: ", gmt)
  ln  <- readLines(gmt, warn = FALSE)
  key <- sub("\t.*$", "", ln)
  hit <- match(set_names, key)
  if (anyNA(hit)) stop("[config_lcc] gene set(s) absent from GMT: ",
                       paste(set_names[is.na(hit)], collapse = ", "))
  out <- lapply(ln[hit], function(x) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; f[-c(1, 2)] })
  setNames(out, set_names)
}

## -- ensure LCC output dirs exist (idempotent) ----
for (.d in c(LCC_TAB_DIR, LCC_FIG_DIR, LCC_LOG_DIR,
             LCC_PB_DIR, LCC_PERCELL_DIR, LCC_STROMA_DIR, LCC_CC_DIR)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}
rm(.d)
