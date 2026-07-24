# config_distance.R ----
# Stage 06_distance config (Phase 5): intensity -> distance (rank-percentile).
# Sources config_ccc.R for the SHARED CCC data model it consumes (CCC_NODES node vocabulary +
# CCC_TENSOR_DIR). distance is strictly DOWNSTREAM of CCC, so this is an acyclic dependency, not the
# sideways cross-stage coupling the standards forbid. (If strict isolation is preferred, promote
# CCC_NODES + CCC_TENSOR_DIR to config_paths.R as project-wide truth and drop the config_ccc source.)
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_ccc.R"))   # CCC_NODES, CCC_TENSOR_DIR (upstream data model)

## -- parameters ----
DIST_PVAL_THRESH    <- 0.05      # explicit significance filter when aggregating LR pairs into a bin-pair edge
DIST_EDGE_WEIGHT    <- "probsum" # primary edge weight = sum(prob) over pval<DIST_PVAL_THRESH (= aggregateNet weight); count = co-product (R8)
DIST_MIN_EDGES_FLAG <- 10L       # < this many present (non-zero) directed bin-pair edges (of 49) -> sparse-graph QC flag (flag only; 07 decides)

## -- output ----
DIR_DISTANCE <- file.path(FAST_DIR, "results/tables/06_distance")   # (ideally promote to config_paths alongside other DIR_*)
