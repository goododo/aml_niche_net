# ============================================================================
# ingest_GSE276360.R   --   Rodriguez-Sevilla 2025 Nat Commun (CCUS + healthy marrow MNC)
# INPUT  : GEO_RAW_DIR/GSE276360_RAW/GSE276360_bm_scrna_expr.rds   dgCMatrix 27,022 x 22,194, raw integer UMI,
#            author-filtered (>= 500 UMI); the GEO file is double-gzipped and this is the copy gunzipped once
#            (2026-10-23); cell ids "<prefix>_<barcode>", prefixes CCUS-50, CCUS-53, CCUS-PR, HD-2, HD-3, HD-36, HD-95
#          RAW_DIR/metadata/GSE276360/donor_ages_GSE276360_v1.0.csv (gsm, donor, fraction, age, sex; authoritative)
# OUTPUT : DIR_INGEST/GSE276360_qc_percell.csv.gz + _qc_summary.csv, RDS_INGEST_DIR/GSE276360.rds
# DOES   : one Seurat object per healthy LIBRARY (one per donor; Sample = "<GSM>_<donor>", Patient_ID = donor), merged
#          UNFILTERED. The three CCUS libraries leave through INGEST_EXCLUDE (key ^CCUS-): the dataset is in
#          TP_HEALTHY_DATASETS, whose rule labels every ingested cell Healthy.
# Role   : S9 descriptive age/direction check (PREREGISTRATION_S9_age_direction.md, tag prereg-s9-v1, sections 2.2 and
#          6 step 3).
#
# THE DONOR MAPPING. Two young donors appear under matrix prefixes that differ from their GEO titles. Resolved on
# 2026-10-23 from XIST and RPS4Y1: HD-36 is male -> YHD658 (23 y, M), HD-95 is female -> YHD636 (24 y, F); HD-2 and
# HD-3 are the elderly donors HD2 (61, F) and HD3 (70, F). The mapping is asserted again below with the S4 sex rule
# (female if XIST > 3 max(RPS4Y1, 1), male if RPS4Y1 > 3 max(XIST, 1)); a mismatch stops the ingest (prereg K3).
# Barcodes lose the matrix prefix (it is kept in the Library column), so merged cell names read "<Sample>_<barcode>".
#
#   Rscript scripts/00_ingest/ingest_GSE276360.R
# This script writes only the outputs listed above; it deletes nothing and does not modify its inputs.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE276360"
STUDY   <- "Rodriguez-Sevilla 2025 Nat Commun"
RDS_IN  <- file.path(GEO_RAW_DIR, "GSE276360_RAW", "GSE276360_bm_scrna_expr.rds")
DONORS  <- file.path(RAW_DIR, "metadata", DATASET, "donor_ages_GSE276360_v1.0.csv")
PREFIX_TO_DONOR <- c("HD-95" = "YHD636", "HD-36" = "YHD658", "HD-2" = "HD2", "HD-3" = "HD3")
N_CELL_EXPECTED <- c("HD-2" = 1921L, "HD-3" = 1798L, "HD-36" = 2995L, "HD-95" = 3499L)   # counted 2026-10-23

age_group_of <- function(a) {   # meta_v2.3 vocabulary, same breaks as build_meta_GSE276360.py
  a <- as.numeric(a)
  ifelse(a < 18, "pediatric_lt18", ifelse(a < 40, "AYA_18_39", ifelse(a < 60, "adult_40_59", "older_60plus")))
}

# ----------------------------------------------------------------------------
# [1] Donor table (healthy BM MNC rows only) ----
# ----------------------------------------------------------------------------
message("[1] donor table")
don <- data.table::fread(DONORS, colClasses = "character")[fraction == "BM MNC"]
stopifnot(setequal(don$donor, PREFIX_TO_DONOR), all(don$sex %in% c("F", "M")))
print(don[, .(gsm, donor, age, sex)])

# ----------------------------------------------------------------------------
# [2] Load the matrix ----
# ----------------------------------------------------------------------------
message("[2] loading ", basename(RDS_IN))
m <- readRDS(RDS_IN)
stopifnot(inherits(m, "dgCMatrix"), all(m@x == round(m@x)), !anyDuplicated(rownames(m)))
message("    ", nrow(m), " genes x ", ncol(m), " cells")
cell_prefix <- sub("_[ACGTN]+-[0-9]+$", "", colnames(m))              # per cell: "HD-36", "CCUS-50", ...
print(table(cell_prefix))

# ----------------------------------------------------------------------------
# [3] Exclusions and the donor mapping ----
# ----------------------------------------------------------------------------
message("[3] exclusions and donor mapping")
kept <- ingest_keep(DATASET, sort(unique(cell_prefix)), "prefix")     # drops CCUS-50, CCUS-53, CCUS-PR
stopifnot(all(grepl("^HD-", kept)), setequal(kept, names(PREFIX_TO_DONOR)))
sex_call <- function(p) {                                             # the S4 rule (26_s4_pseudobulk.R)
  i <- cell_prefix == p; x <- sum(m["XIST", i]); y <- sum(m["RPS4Y1", i])
  if (x > 3 * max(y, 1)) "F" else if (y > 3 * max(x, 1)) "M" else "unclear"
}
n_by_prefix <- table(cell_prefix)
map <- data.table::data.table(prefix = kept, donor = unname(PREFIX_TO_DONOR[kept]))
map[, gsm := don$gsm[match(donor, don$donor)]]                        # inside map[...], `donor` is map's column
map[, table_sex := don$sex[match(donor, don$donor)]]
map[, inferred_sex := vapply(prefix, sex_call, "")]
map[, n_cells := as.integer(n_by_prefix[prefix])]
print(map)
stopifnot(all(map$inferred_sex == map$table_sex),                     # prereg K3 at ingest: the mapping must hold
          all(map$n_cells == N_CELL_EXPECTED[map$prefix]))

# ----------------------------------------------------------------------------
# [4] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[4] building per-library Seurat objects")
obj_list <- list()
for (k in seq_len(nrow(map))) {
  p   <- map$prefix[k]
  d   <- don[donor == map$donor[k]]
  key <- paste0(d$gsm, "_", d$donor)                                  # "GSM8497031_YHD658"
  counts <- m[, cell_prefix == p, drop = FALSE]
  colnames(counts) <- sub("^.*_", "", colnames(counts))              # "HD-36_AAAC...-1" -> "AAAC...-1"
  stopifnot(!anyDuplicated(colnames(counts)))
  n <- ncol(counts)
  extra <- data.frame(age_years = rep(d$age, n), age_group = rep(age_group_of(d$age), n),
                      sex = rep(ifelse(d$sex == "F", "female", "male"), n),   # the words 26_s4_pseudobulk.R expects
                      row.names = colnames(counts), stringsAsFactors = FALSE)
  message(sprintf("    - %-18s prefix %-6s age %s  sex %s  cells %d", key, p, d$age, d$sex, ncol(counts)))
  s <- make_seurat(counts = counts, sample = key, dataset = DATASET, patient = d$donor,
                   timepoint = "Baseline", study = STUDY, disease = "Healthy", extra_meta = extra)
  s$GSM     <- d$gsm        # decontX channel unit (ws/22)
  s$Library <- p            # the matrix prefix in the deposit
  obj_list[[key]] <- s
}
rm(m); gc()

# ----------------------------------------------------------------------------
# [5] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[5] merging ", length(obj_list), " libraries")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [6] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[6] writing outputs")
smry <- write_qc_outputs(seu, DATASET)

# ----------------------------------------------------------------------------
# [7] Sanity checks on the output ----
# ----------------------------------------------------------------------------
message("[7] sanity checks")
print(as.data.frame(smry)[, c("Sample", "Patient_ID", "Timepoint", "n_cells", "median_nFeature", "max_pct_mt")])
stopifnot(nrow(smry) == 4L, sum(smry$n_cells) == sum(N_CELL_EXPECTED), all(smry$Timepoint == "Healthy"),
          all(smry$max_pct_mt > 0))          # an all-zero percent.mt would mean the gene names were read wrongly

message("[ok] GSE276360 finished")
