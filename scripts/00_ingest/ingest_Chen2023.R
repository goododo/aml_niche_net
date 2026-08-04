# ============================================================================
# ingest_Chen2023.R   --   Chen 2023 Blood Cancer Discovery (paired BM niche)
# Format : nested 10x dirs
#          "<root>/filtered_feature_bc_matrix.7z_<SAMPLE>/filtered_feature_bc_matrix/"
#          SAMPLE like AML103_CD34, AML162_niche_immune, NBM1_CD34, ...
# Timepoint: NBM* -> healthy "Baseline"; AML* -> "Diagnosis".
#
# THE SAMPLE IS THE DONOR. 20 deposit dirs = 10 donors x 2 pools, and a pool is HALF A SAMPLE,
# not a sample. Per Methods each donor was FACS-split into five fractions and re-pooled into
# two libraries that were sequenced separately:
#     pool A (deposit "<donor>_CD34")          = HSPC (CD45+CD34+) + myeloid (CD45+CD34-CD117/CD33+)
#     pool B (deposit "<donor>_Niche_Immune")  = stromal + endothelial + lymphoid
# The deposit name "_CD34" describes only half of pool A's content, and v1 read it as "CD34-
# enriched, composition biased, exclude" in two independent places (SORTED_LIBRARY_REGEX in
# config_malignancy.R, CCC_SORTED_SUBLIB_PATTERN in config_ccc.R). That threw away 5 of the 7
# CCC nodes -- HSC_MPP, LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte all live in pool A -- and
# left Chen2023 as a T/B-only graph, with malignant_frac computed on the lymphoid pool alone.
# Both rules are removed. The composition bias is real but it is a COVARIATE, carried by the
# curated field sorting = "multi_fraction_FACS_recombined", not a reason to drop half the marrow.
#
# Merging: Sample = Patient_ID, so the two pools land in the same QC / CNV / CCC unit. The pools
# are independent 10x runs with colliding barcodes, so merge_samples() must disambiguate --
# obj_list is keyed by the LIBRARY dir and passes those names to add.cell.ids.
#
# Case: the niche dirs have FOUR spellings (Niche_Immune x7, niche_immune AML162,
# NICHE_IMMUNE AML320, Niche_immune NBM3). Every match here is case-insensitive.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "Chen2023"
STUDY   <- "Chen 2023 Blood Cancer Discovery"
ROOT    <- file.path(ZENODO_DIR, "10.17632_gwjh3w6ztm.2")

# ----------------------------------------------------------------------------
# [1] Discover per-sample 10x directories ----
# ----------------------------------------------------------------------------
message("[1] listing sample directories")
sample_dirs <- list.dirs(ROOT, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("filtered_feature_bc_matrix\\.7z_", basename(sample_dirs))]
message("    ", length(sample_dirs), " sample folders found")

# Derive pool identity and Patient_ID from a deposit dir name. Named for what the library
# actually CONTAINS (Methods), not for what the deposit called it.
pool_of    <- function(s) {
  if (grepl("CD34$", s, ignore.case = TRUE)) "poolA_HSPC_myeloid" else "poolB_stroma_lymphoid"
}
patient_of <- function(s) sub("_(CD34|niche_immune)$", "", s, ignore.case = TRUE)

# ----------------------------------------------------------------------------
# [2] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[2] building per-library Seurat objects")
obj_list <- list()

for (d in sample_dirs) {
  library_id  <- sub("^filtered_feature_bc_matrix\\.7z_", "", basename(d))
  inner       <- file.path(d, "filtered_feature_bc_matrix")
  patient     <- patient_of(library_id)
  pool        <- pool_of(library_id)
  tp          <- if (grepl("^NBM", patient)) "Baseline" else "Diagnosis"
  message("    - ", library_id, "  (patient=", patient, ", pool=", pool, ")")

  counts <- read_10x_dir(inner)
  s <- make_seurat(
    counts    = counts,
    sample    = patient,       # THE SAMPLE IS THE DONOR: both pools share one Sample key
    dataset   = DATASET,
    patient   = patient,
    timepoint = tp,
    study     = STUDY
  )
  s$Library     <- library_id  # provenance: which of the donor's two runs this cell came from
  s$Pool        <- pool
  s$Compartment <- pool        # kept for backward compatibility with existing readers
  obj_list[[library_id]] <- s  # names -> add.cell.ids, so colliding pool barcodes stay distinct
}

# Every donor must contribute both pools; a donor with one pool is a half sample and its node
# set would be structurally incomplete in a way FGW must not read as biology.
.n_pool <- table(vapply(obj_list, function(x) x$Patient_ID[1], character(1)))
if (any(.n_pool != 2L))
  warning("Chen2023: donor(s) not contributing exactly 2 pools: ",
          paste(names(.n_pool)[.n_pool != 2L], collapse = ", "), call. = FALSE)

# ----------------------------------------------------------------------------
# [3] Merge into one dataset object (20 libraries -> 10 donor-level samples) ----
# ----------------------------------------------------------------------------
message("[3] merging ", length(obj_list), " libraries into ", length(.n_pool), " donor samples")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [4] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[4] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] Chen2023 finished")
