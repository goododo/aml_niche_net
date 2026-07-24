# =============================================================================
# 02_study_split.R   (Task 2)
# Study-level (NOT sample-level) 70/30 split, stratified by genetic subtype.
#   - Randomisation acts ONLY on the cross-sectional Discovery-AML candidate
#     studies (ROLE_TABLE$discovery_candidate == TRUE & !is_longitudinal).
#   - Reference / Healthy studies are assigned deterministically.
#   - ALL longitudinal / relapse studies are FORCED to the Validation/treatment
#     side (and never enter the random pool) -> prevents same-patient leakage.
#
# Input : results/tables/01_qc/01_sample_role_manifest.csv  (from Task 1)
# Output: results/tables/01_qc/02_study_split.csv           (study-level)
#         results/tables/01_qc/02_sample_split.csv          (sample-level join)
# Run   : Rscript 02_study_split.R
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
suppressPackageStartupMessages({ library(data.table) })

set.seed(SEED)   # reproducible split

manifest <- fread(file.path(DIR_PREPROCESS, "01_sample_role_manifest.csv"))

# ----------------------------------------------------------------------------
# Step 1. Study-level table with the variables needed for the split.
# ----------------------------------------------------------------------------
studies <- unique(manifest[, .(study_role = first(study_role),
                               is_longitudinal = first(is_longitudinal),
                               discovery_candidate = first(discovery_candidate),
                               subtype_stratum = first(subtype_stratum),
                               n_discovery_samples = sum(in_discovery_pool)),
                           by = dataset])

# ----------------------------------------------------------------------------
# Step 2. Deterministic assignment for everything outside the random pool.
# ----------------------------------------------------------------------------
studies[, split := NA_character_]
studies[study_role == "Healthy-control",                  split := "Healthy"]
studies[grepl("Reference", study_role),                   split := "Reference"]
studies[is_longitudinal == TRUE,                          split := "Validation"]  # forced
studies[study_role == "Exploratory",                      split := "Exploratory"]

# Random pool = cross-sectional discovery-candidate studies still unassigned.
pool <- studies[discovery_candidate == TRUE &
                is_longitudinal == FALSE &
                is.na(split), dataset]

# ----------------------------------------------------------------------------
# Step 3. Stratified 70/30 split within the pool, stratum = subtype.
#   With only ~5 candidate studies, stratification is intentionally coarse;
#   we therefore (i) split within each stratum at 70% (round to >=1 discovery),
#   and (ii) enforce a global guard: at least 1 pool study must land in
#   Validation, so cross-cohort validation is never empty.
# ----------------------------------------------------------------------------
disc_frac <- 0.70
pool_dt   <- studies[dataset %in% pool]

split_one_stratum <- function(ds_vec) {
  n  <- length(ds_vec)
  if (n == 0) return(character(0))
  ds_vec <- sample(ds_vec)                     # shuffle
  n_disc <- max(1L, round(n * disc_frac))      # at least 1 discovery per stratum
  n_disc <- min(n_disc, n)                      # cannot exceed stratum size
  res <- setNames(rep("Validation", n), ds_vec)
  res[seq_len(n_disc)] <- "Discovery"
  res
}

assign_dt <- rbindlist(lapply(
  split(pool_dt$dataset, pool_dt$subtype_stratum),
  function(ds_vec) {
    res <- split_one_stratum(ds_vec)
    data.table(dataset = names(res), split = unname(res))
  }))
# Join by dataset (robust: no reliance on unlist() name mangling).
studies[assign_dt, on = "dataset", split := i.split]

# Global guard: ensure >=1 pool study in Validation (demote a random Discovery
# pool study if the stochastic draw sent them all to Discovery).
pool_split <- studies[dataset %in% pool, split]
if (!any(pool_split == "Validation", na.rm = TRUE)) {
  demote <- sample(studies[dataset %in% pool & split == "Discovery", dataset], 1)
  studies[dataset == demote, split := "Validation"]
  message_ts("guard: demoted ", demote, " to Validation to keep validation non-empty")
}

# ----------------------------------------------------------------------------
# Step 4. Persist study- and sample-level split tables.
# ----------------------------------------------------------------------------
setorder(studies, split, dataset)
fwrite_safe(studies, file.path(DIR_PREPROCESS, "02_study_split.csv"))

sample_split <- merge(manifest, studies[, .(dataset, split)], by = "dataset", all.x = TRUE)
# Sample-level split label: longitudinal AML -> "Validation/Treatment-axis";
# healthy samples in any study -> "Healthy" regardless of their study's split.
sample_split[, split_sample := split]
sample_split[sample_role == "Healthy-control", split_sample := "Healthy"]
sample_split[grepl("^Exclude", sample_role),   split_sample := "Excluded"]
fwrite_safe(sample_split, file.path(DIR_PREPROCESS, "02_sample_split.csv"))

message_ts("Task 2 done. Seed = ", SEED)
cat("\n--- Study-level split ---\n")
print(studies[, .(dataset, study_role, subtype_stratum, is_longitudinal, split)])
cat("\n--- Study counts by split ---\n")
print(studies[, .N, by = split])