#!/usr/bin/env Rscript
# 50_consensus_malignancy.R ----
# Type-level malignancy consensus, layered onto QC cells (R2). The blueprint's ">=2/3" means
# >=2 of THREE evidence TYPES (expression-CNV / allele-CNV / SNV), NOT >=2 of N tools -- because
# expression-CNV callers (inferCNV, CopyKAT, SCEVAN, author CopyKAT) are correlated and must not
# be counted as independent votes.
#
# Arms (standardized schema: cell, malignant, score, method, sample):
#   EXPRESSION-CNV : --infercnv (41) | --copykat (42) | --scevan (42) | --author (43)
#   ALLELE-CNV     : --numbat   (30 native percell; mapped here)
#   SNV            : --vartrix  (48 [planned])
#
# Voting:
#   expr_call   = strict-majority of available expression callers (tie -> NA/ambiguous)
#   allele_call = Numbat ;  snv_call = VarTrix
#   n_types     = how many of {expr, allele, snv} have a call
#   malignant   = production default (union): >=2 types -> malignant if ANY type positive
#                 (--majority forces the stricter >=2-of-present-types rule instead)
#                 (pilot, --min_tools 1) single available type is taken, confidence=low
#   evidence_tier (unified): C_single (1 type) | A_concordant (>=2 types, conflict<X) |
#                            B_multi_partial (>=2 types, conflict>=X)   X = --conflict_max
#
# Barcodes differ across tools; we join on the 16bp DNA core, per sample.
#
# 3b logic changes (vs the 3a rename):
#   * SCHEMA FIX: a degraded Numbat file (note == NUMBAT_DEGRADED_NOTE) or one with no usable
#     malignant/compartment column ABSTAINS (read_arm returns NULL) instead of voting malignant=0.
#     It stays out of `present`, so an inferCNV-only sample lands tier C, not a false tier A / a
#     false "confident normal". [d35 parity; consistent with D-A7]
#   * --union_mode (default CONSENSUS_UNION_MODE): >=2 types -> malignant if ANY type is positive
#     (blind-spot complementarity, e.g. inferCNV~0 but Numbat high), vs the default type-majority.
#   * UNIFIED tier vocabulary A_concordant / B_multi_partial / C_single, split by the sample's
#     cross-type conflict fraction (--conflict_max, default TIER_CONFLICT_MAX = 0.10).
#   * core16() now comes from utils.R (single implementation); --min_tools defaults to
#     CONSENSUS_MIN_TOOLS; the degraded-note literal is the shared NUMBAT_DEGRADED_NOTE.
#
# Usage (pilot, single sample; run from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/50_consensus_malignancy.R --sample 6323_R --dataset GSE227903 \
#     --qc_rds  .../6323_R.rds \
#     --numbat  .../6323_R__numbat_percell.csv \
#     --infercnv .../6323_R__infercnv_percell.csv \
#     --outdir  .../02_malignancy/GSE227903 [--union_mode] [--majority] [--min_tools 1] [--conflict_max 0.10]
#   (voting defaults to union for production, config CONSENSUS_UNION_MODE=TRUE; no flag needed.
#    --union_mode just re-affirms union (bare, safe); --majority forces type-majority (bare).)

suppressPackageStartupMessages({ library(optparse); library(data.table) })

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--sample",   type = "character"),
  make_option("--dataset",  type = "character", default = NA),
  make_option("--qc_rds",   type = "character", help = "per-sample QC Seurat .rds (canonical cells)"),
  make_option("--numbat",   type = "character", default = ""),   # allele
  make_option("--vartrix",  type = "character", default = ""),   # snv
  make_option("--infercnv", type = "character", default = ""),   # expression
  make_option("--copykat",  type = "character", default = ""),   # expression
  make_option("--scevan",   type = "character", default = ""),   # expression
  make_option("--author",   type = "character", default = ""),   # expression (author CopyKAT etc.)
  make_option("--outdir",   type = "character"),
  make_option("--min_tools",type = "integer", default = CONSENSUS_MIN_TOOLS,
              help = "1 = pilot (allow a single evidence type); 2 = production (>=2/3 types) [%default]"),
  make_option("--union_mode", action = "store_true", default = FALSE,
              help = "affirm union voting (>=2 types -> malignant if ANY positive). Already the production default via config CONSENSUS_UNION_MODE=TRUE, so safe/optional to pass. Bare flag, no value."),
  make_option("--majority", action = "store_true", default = FALSE,
              help = "force type-majority voting (>=2 of the present types) instead of union. Bare flag, no value. (Overrides union if both are given.)"),
  make_option("--conflict_max", type = "double", default = TIER_CONFLICT_MAX,
              help = "cross-type conflict fraction splitting A_concordant vs B_multi_partial [%default]")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

EXPR_ARMS   <- c("infercnv", "copykat", "scevan", "author")
ALLELE_ARMS <- c("numbat")
SNV_ARMS    <- c("vartrix")
# core16() is sourced from utils.R (single implementation with barcode-length validation).

## ---- 1. canonical cells from the QC object ----
message("[1] Loading QC cells: ", opt$qc_rds)
qc <- readRDS(opt$qc_rds)
qc_cells <- tryCatch(colnames(qc), error = function(e) NULL)
if (is.null(qc_cells)) qc_cells <- rownames(qc@meta.data)
base <- data.table(cell = qc_cells, key16 = core16(qc_cells))
message("    ", length(qc_cells), " QC cells")

## ---- 2. ingest each available arm -> per-arm malignant column ----
read_arm <- function(path, method) {
  if (!nzchar(path) || !file.exists(path)) { message("    [", method, "] absent"); return(NULL) }
  d <- fread(path)
  if (method == "numbat") {
    # degraded Numbat (note == NUMBAT_DEGRADED_NOTE) contributes NO allele evidence -> ABSTAIN.
    # Returning NULL keeps numbat OUT of `present`, so a sample left with only inferCNV correctly
    # lands tier C (single) instead of a false "confident normal" (the old malignant=0 bug) or a
    # false tier A. [3b schema fix; d35 parity; consistent with D-A7]
    if ("note" %in% names(d) && all(d$note == NUMBAT_DEGRADED_NOTE, na.rm = TRUE)) {
      message("    [numbat] degraded (", NUMBAT_DEGRADED_NOTE, ") -> abstain"); return(NULL)
    }
    if (!"malignant" %in% names(d)) {
      if ("compartment_opt" %in% names(d)) {
        d[, malignant := as.integer(compartment_opt == "tumor")]
      } else {
        message("    [numbat] no malignant/compartment column -> abstain"); return(NULL)
      }
    }
  }
  if (!"cell" %in% names(d)) { message("    [", method, "] no 'cell' col"); return(NULL) }
  d[, key16 := core16(cell)]; d <- d[!is.na(malignant)]
  message("    [", method, "] ", nrow(d), " cells, ", sum(d$malignant == 1), " malignant")
  unique(d[, .(key16, m = malignant)], by = "key16")
}
arm_paths <- list(numbat = opt$numbat, vartrix = opt$vartrix, infercnv = opt$infercnv,
                  copykat = opt$copykat, scevan = opt$scevan, author = opt$author)
present <- c()
for (nm in names(arm_paths)) {
  a <- read_arm(arm_paths[[nm]], nm)
  if (!is.null(a)) { base[a, paste0("m_", nm) := i.m, on = "key16"]; present <- c(present, nm) }
}
if (!length(present)) stop("No malignancy arms available for ", opt$sample)
message("[2] Arms present: ", paste(present, collapse = ", "))

## ---- 3. collapse to per-TYPE calls, then vote across types ----
union_mode <- (isTRUE(opt$union_mode) || isTRUE(CONSENSUS_UNION_MODE)) && !isTRUE(opt$majority)
message("[3] Type-level voting (min_tools=", opt$min_tools,
        ", union_mode=", union_mode, ")")
expr_cols <- intersect(paste0("m_", EXPR_ARMS), names(base))
strict_majority <- function(dt, cols) {
  if (!length(cols)) return(rep(NA_integer_, nrow(dt)))
  n   <- rowSums(!is.na(dt[, ..cols]))
  mal <- rowSums(dt[, ..cols] == 1, na.rm = TRUE)
  fifelse(n == 0, NA_integer_,
          fifelse(mal * 2 > n, 1L, fifelse((n - mal) * 2 > n, 0L, NA_integer_)))  # tie -> NA
}
base[, expr_call   := strict_majority(.SD, expr_cols)]
base[, allele_call := if ("m_numbat"  %in% names(base)) m_numbat  else NA_integer_]
base[, snv_call    := if ("m_vartrix" %in% names(base)) m_vartrix else NA_integer_]

tcols <- c("expr_call", "allele_call", "snv_call")
base[, n_types := rowSums(!is.na(.SD)), .SDcols = tcols]
base[, n_types_mal := rowSums(.SD == 1, na.rm = TRUE), .SDcols = tcols]

base[, malignant := NA_integer_]
if (union_mode) {
  # any evidence type positive -> malignant (blind-spot complementarity; former d35 union behavior)
  base[n_types >= 2, malignant := as.integer(n_types_mal >= 1)]
} else {
  # default: type-level majority (>=2 of the present types agree malignant)
  base[n_types >= 2, malignant := as.integer(n_types_mal >= 2)]
}
if (opt$min_tools <= 1) base[n_types == 1, malignant := as.integer(n_types_mal >= 1)]

base[, confidence := fifelse(n_types >= 2 & n_types_mal >= 2, "high",
                     fifelse(n_types >= 2 & n_types_mal == 1, "conflicting",
                     fifelse(n_types == 1, "low_single_type",
                     fifelse(n_types == 0, "no_evidence", "normal_consensus"))))]
base[n_types >= 2 & n_types_mal == 0, confidence := "high"]   # confident normal

if (!is.na(opt$dataset)) base[, dataset := opt$dataset]
base[, sample := opt$sample]

## ---- 4. write per-cell + per-sample summary ----
keep <- c("cell","sample","malignant","confidence","n_types","n_types_mal",
          "expr_call","allele_call","snv_call", intersect(paste0("m_",names(arm_paths)), names(base)))
out_cell <- file.path(opt$outdir, paste0(opt$sample, "__consensus_percell.csv"))
fwrite_safe(base[, ..keep], out_cell)

# sample-level evidence-type presence (expr / allele / snv) and cross-type conflict
type_present <- c(expr   = length(intersect(EXPR_ARMS, present)) > 0,
                  allele = any(ALLELE_ARMS %in% present),
                  snv    = any(SNV_ARMS   %in% present))
n_types_sample <- sum(type_present)
# among cells with >=2 types, fraction whose types DISAGREE (not unanimous: 0 < n_types_mal < n_types)
multi <- base[n_types >= 2]
conflict_frac <- if (nrow(multi)) mean(multi$n_types_mal > 0 & multi$n_types_mal < multi$n_types) else 0
evidence_tier <- if (n_types_sample <= 1) "C_single" else
                 if (conflict_frac < opt$conflict_max) "A_concordant" else "B_multi_partial"
n_called <- base[!is.na(malignant), .N]
summ <- data.table(
  sample = opt$sample, dataset = opt$dataset, n_qc_cells = nrow(base),
  arms = paste(present, collapse = "+"), n_types = n_types_sample,
  evidence_tier = evidence_tier, conflict_frac = round(conflict_frac, 4),
  n_called = n_called, n_malignant = base[malignant == 1, .N],
  malignant_frac = round(base[malignant == 1, .N] / max(n_called, 1), 4),
  n_high = base[confidence == "high", .N],
  n_conflicting = base[confidence == "conflicting", .N],
  n_low_single = base[confidence == "low_single_type", .N])
out_summ <- file.path(opt$outdir, paste0(opt$sample, "__consensus_summary.csv"))
fwrite_safe(summ, out_summ)

message("[done] ", opt$sample, ": tier=", evidence_tier, " malignant_frac=", summ$malignant_frac,
        " (", summ$n_malignant, "/", n_called, "), arms=", summ$arms,
        "\n        per-cell -> ", out_cell, "\n        summary  -> ", out_summ)
