# ============================================================================
# ingest_GSE120221.R   --   Oetjen 2018 JCI Insight (healthy marrow, 20 donors 24-84 y)
# INPUT  : GEO_RAW_DIR/GSE120221_RAW/  -- 25 Cell Ranger 2.0 filtered triplets
#            "<GSM>_{barcodes,genes,matrix}_<lib>.{tsv,mtx}.gz"; genes.tsv = Ensembl id + CR 2.0 symbol
#          RAW_DIR/metadata/GSE120221/donor_ages_GSE120221_v1.0.csv  (donor, age, sex per library; authoritative)
#          REF_FEATURES_2020A  (GRCh38 2020-A features table of GSE289435 = the project's symbol set)
#          results/tables/ws/ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv  (read only: coverage assertion)
# OUTPUT : DIR_INGEST/GSE120221_qc_percell.csv.gz + _qc_summary.csv, RDS_INGEST_DIR/GSE120221.rds
#          DIR_INGEST/GSE120221_gene_remap.csv  (per gene: Ensembl id, CR 2.0 symbol, 2020-A name, final name, rule)
# DOES   : one Seurat object per LIBRARY (Sample = "<GSM>_<lib>", Patient_ID = donor letter), merged UNFILTERED.
#          Genes are renamed by Ensembl id to the 2020-A names every other 10x ingest carries.
# Role   : S9 age/direction control, the DECISIVE cohort (scripts/ws/PREREGISTRATION_S9_age_direction.md,
#          tag prereg-s9-v1, sections 2.2-2.3 and 6 step 3). Healthy throughout; no timepoints.
#
# WHY SAMPLE = LIBRARY, NOT DONOR. Donors C (C1, C2, Ck) and S (S1, S2, Sk1, Sk2) have several libraries.
# 03_per_sample_qc.R sets the doublet prior from the Sample's cell count, so a donor-level Sample would triple
# (C) or quadruple (S) it, and C is one of the five donors >= 60 the age test rests on. Libraries are pooled
# per donor later, inside scripts/ws/33_s9_age_direction.R (prereg section 2.3).
#
# WHY THE GENE REMAP. Cell Ranger 2.0 used Ensembl-84-era symbols: 10,784 of 33,694 differ from 2020-A, and
# 44 of the 370 primary imprint genes would be missing by symbol. Rule (prereg 2.3): a gene whose Ensembl id is
# in 2020-A takes 2020-A's name for that id (the make.unique'd name, as ReadMtx gives every 2020-A ingest, so
# MATR3 / MATR3.1 point to the same ids as in GSE185381); an unmapped id keeps its CR 2.0 symbol unless that
# symbol is already used by a mapped gene on this axis, in which case it becomes "<symbol>.<Ensembl id>";
# underscores become dashes (Seurat would do it silently); remaining duplicates get make.unique.
#
#   Rscript scripts/00_ingest/ingest_GSE120221.R
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

DATASET <- "GSE120221"
STUDY   <- "Oetjen 2018 JCI Insight"
RAW     <- file.path(GEO_RAW_DIR, "GSE120221_RAW")
DONORS  <- file.path(RAW_DIR, "metadata", DATASET, "donor_ages_GSE120221_v1.0.csv")
TABWS   <- file.path(TAB_DIR, "ws")
N_LIB_EXPECTED   <- 25L
N_CELL_EXPECTED  <- 90653L     # barcodes in the 25 files (counted 2026-10-23)
MIN_PRIMARY_COV  <- 0.95       # prereg 6 step 3: >= 0.95 of the 370 primary genes on the axis

age_group_of <- function(a) {   # meta_v2.3 vocabulary, same breaks as build_meta_GSE120221.py
  a <- as.numeric(a)
  ifelse(a < 18, "pediatric_lt18", ifelse(a < 40, "AYA_18_39", ifelse(a < 60, "adult_40_59", "older_60plus")))
}

# ----------------------------------------------------------------------------
# [1] Donor table (all columns as character: donor letters "T" and "F" must not become logicals) ----
# ----------------------------------------------------------------------------
message("[1] donor table")
don <- data.table::fread(DONORS, colClasses = "character")
stopifnot(nrow(don) == N_LIB_EXPECTED, !anyDuplicated(don$gsm), all(don$sex %in% c("F", "M")))
print(head(don[, .(gsm, library, donor, age, sex)], 5))

# ----------------------------------------------------------------------------
# [2] Discover library prefixes ----
# ----------------------------------------------------------------------------
message("[2] listing libraries")
bc_files <- list.files(RAW, pattern = "_barcodes_.*\\.tsv\\.gz$", full.names = FALSE)
libs <- data.table::data.table(
  gsm = sub("_barcodes_.*$", "", bc_files),                           # "GSM3396161"
  lib = sub("^.*_barcodes_(.*)\\.tsv\\.gz$", "\\1", bc_files))       # "A", "C1", "Sk2", ...
libs[, key := paste0(gsm, "_", lib)]                                  # the pipeline Sample id
libs <- libs[key %in% ingest_keep(DATASET, libs$key, "library")]     # no GSE120221 rows today; kept for the audit trail
stopifnot(nrow(libs) == N_LIB_EXPECTED,
          setequal(libs$gsm, don$gsm),
          all(don$library[match(libs$gsm, don$gsm)] == libs$lib))     # file library label == donor-table label
message("    ", nrow(libs), " libraries; donors: ", data.table::uniqueN(don$donor))

# ----------------------------------------------------------------------------
# [3] Gene remap: Ensembl id -> the project's 2020-A names (rule in the header) ----
# ----------------------------------------------------------------------------
message("[3] gene remap to the 2020-A symbol set")
gfile <- function(gsm, lib) file.path(RAW, paste0(gsm, "_genes_", lib, ".tsv.gz"))
g <- data.table::fread(gfile(libs$gsm[1], libs$lib[1]), header = FALSE, col.names = c("ens", "sym_cr2"))
for (i in seq_len(nrow(libs))[-1]) {                                   # every library must carry the same gene table
  gi <- data.table::fread(gfile(libs$gsm[i], libs$lib[i]), header = FALSE, col.names = c("ens", "sym_cr2"))
  stopifnot(identical(gi$ens, g$ens), identical(gi$sym_cr2, g$sym_cr2))
}
stopifnot(!anyDuplicated(g$ens))
ref <- data.table::fread(REF_FEATURES_2020A, header = FALSE, select = 1:2, col.names = c("ens", "sym"))
ref[, sym_unique := make.unique(sym)]                                  # what ReadMtx(unique.features = TRUE) gives a 2020-A ingest
g[, sym_2020 := ref$sym_unique[match(ens, ref$ens)]]
g[, rule := data.table::fifelse(!is.na(sym_2020), "mapped_by_ensembl", "kept_cr2")]
g[, final := data.table::fifelse(rule == "mapped_by_ensembl", sym_2020, sym_cr2)]
g[rule == "kept_cr2" & grepl("_", final), `:=`(final = gsub("_", "-", final), rule = "kept_cr2_underscore_to_dash")]
coll <- g$rule != "mapped_by_ensembl" & g$final %in% g[rule == "mapped_by_ensembl", final]
g[coll, `:=`(final = paste0(final, ".", ens), rule = "collision_suffixed")]
g[, final_unique := make.unique(final)]
g[final_unique != final, rule := paste0(rule, "+made_unique")]
print(g[, .N, by = rule])
stopifnot(!anyDuplicated(g$final_unique), !any(grepl("_", g$final_unique)))

# coverage of the primary imprint genes (prereg 2.3 / K4 input): present, and never on an ambiguous row
imp  <- data.table::fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
amb  <- data.table::fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))
prim <- setdiff(imp$gene, amb$gene)
stopifnot(length(prim) == 370L)
hit  <- g[final_unique %in% prim]
cov  <- nrow(hit) / length(prim)
message(sprintf("    primary genes on the axis: %d / %d (%.3f); rules: %s", nrow(hit), length(prim), cov,
                paste(names(table(hit$rule)), table(hit$rule), sep = " ", collapse = "; ")))
message("    primary genes absent: ", paste(setdiff(prim, g$final_unique), collapse = " "))
stopifnot(cov >= MIN_PRIMARY_COV, all(hit$rule %in% c("mapped_by_ensembl", "kept_cr2")))
fwrite_safe(g[, .(ens, sym_cr2, sym_2020, final = final_unique, rule)], file.path(DIR_INGEST, paste0(DATASET, "_gene_remap.csv")))

# ----------------------------------------------------------------------------
# [4] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[4] building per-library Seurat objects")
obj_list <- list(); sex_list <- list()
for (i in seq_len(nrow(libs))) {
  gsm <- libs$gsm[i]; lib <- libs$lib[i]; key <- libs$key[i]
  d   <- don[gsm == libs$gsm[i]]                                      # inside don[...], `gsm` is don's column
  counts <- Seurat::ReadMtx(mtx      = file.path(RAW, paste0(gsm, "_matrix_", lib, ".mtx.gz")),
                            cells    = file.path(RAW, paste0(gsm, "_barcodes_", lib, ".tsv.gz")),
                            features = gfile(gsm, lib),
                            feature.column = 1)                       # Ensembl ids (unique), renamed just below
  stopifnot(all(rownames(counts) %in% g$ens))
  rownames(counts) <- g$final_unique[match(rownames(counts), g$ens)]
  n <- ncol(counts)
  extra <- data.frame(age_years = rep(d$age, n), age_group = rep(age_group_of(d$age), n),
                      sex = rep(ifelse(d$sex == "F", "female", "male"), n),   # the words 26_s4_pseudobulk.R compares with XIST/RPS4Y1
                      row.names = colnames(counts), stringsAsFactors = FALSE)
  sex_list[[key]] <- data.table::data.table(Sample = key, table_sex = d$sex,
                                            xist = sum(counts["XIST", ]), rps4y1 = sum(counts["RPS4Y1", ]))
  message(sprintf("    - %-14s donor %s  age %s  sex %s  cells %d", key, d$donor, d$age, d$sex, n))
  s <- make_seurat(counts = counts, sample = key, dataset = DATASET, patient = d$donor,
                   timepoint = "Baseline", study = STUDY, disease = "Healthy", extra_meta = extra)
  s$GSM     <- gsm          # decontX channel unit (ws/22)
  s$Library <- lib          # the raw library label of the deposit
  obj_list[[key]] <- s
}

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
print(head(as.data.frame(smry)[, c("Sample", "Patient_ID", "Timepoint", "n_cells", "median_nFeature", "max_pct_mt")], 6))
stopifnot(nrow(smry) == N_LIB_EXPECTED, sum(smry$n_cells) == N_CELL_EXPECTED,
          all(smry$Timepoint == "Healthy"),
          all(smry$max_pct_mt > 0))          # an all-zero percent.mt would mean the gene names were read wrongly
# sex genes per library vs the donor table (printed; the binding identity check is K3 in 33_s9_age_direction.R)
sx <- data.table::rbindlist(sex_list)
sx[, inferred := data.table::fifelse(xist > 3 * pmax(rps4y1, 1), "F", data.table::fifelse(rps4y1 > 3 * pmax(xist, 1), "M", "unclear"))]
print(sx)
message("    sex genes agree with the donor table in ", sum(sx$inferred == sx$table_sex), " of ", nrow(sx), " libraries")
if (any(sx$inferred != sx$table_sex)) warning("sex-gene call differs from the donor table: ", paste(sx[inferred != table_sex, Sample], collapse = ", "))

message("[ok] GSE120221 finished")
