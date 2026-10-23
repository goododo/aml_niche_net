# ============================================================================
# ingest_GSE227690.R   --   Griffin 2023 Nature (healthy marrow BM1-BM5 of a BPDCN study)
# INPUT  : GEO_RAW_DIR/GSE227690_RAW/GSE227690_BM_counts.tsv.gz    dense genes x cells, 19,939 x 20,411, author-filtered
#          GEO_RAW_DIR/GSE227690_RAW/GSE227690_BM_metadata.tsv.gz  per cell: replicate (BM1.1 ...), tech, CellType, ...
#          RAW_DIR/metadata/GSE227690/donor_ages_GSE227690_v1.0.csv (gsm, library BM1_BioRep1, donor, age, sex; authoritative)
# OUTPUT : DIR_INGEST/GSE227690_qc_percell.csv.gz + _qc_summary.csv, RDS_INGEST_DIR/GSE227690.rds
# DOES   : splits the series-level matrix by replicate library, one Seurat object per LIBRARY
#          (Sample = "<GSM>_<donor>_BioRep<k>", Patient_ID = donor), merged UNFILTERED. Author columns kept as anno_*.
# Role   : S9 descriptive age/direction check (PREREGISTRATION_S9_age_direction.md, tag prereg-s9-v1, sections 2.2 and
#          6 step 3). BM6 (replicate BM6.1) is Seq-Well and leaves through INGEST_EXCLUDE (key ^BM6).
#
# WHAT THE DEPOSIT IS. Raw integer UMI counts after the authors' own cell QC (>= 2,000 UMI per cell) and with genes
# removed: all MT- genes, XIST, MALAT1 and NEAT1 are absent (checked 2026-10-23), so per-cell column sums are a few
# UMI below the authors' nCount_RNA and percent.mt is 0 for every cell (add_qc_metrics warns once per library; this is
# a property of the deposit, not a reading error). Cell ids are "<barcode>-<donor>.<rep>", e.g. CACTTCGTCATTCATC-BM1.1;
# BM1.1 is BioRep1 of donor BM1 = GSM7105550 (the donor table carries the GSM). Cell ids are kept as deposited so that
# they stay joinable to the author metadata. Sex can only be read from Y-linked genes here (no XIST).
#
# ONE MALFORMED ROW IN THE DEPOSIT. Line 2451 of BM_counts.tsv.gz (gene CCDC60) carries 23,679 values for 20,411 cells;
# every other line has exactly 20,412 fields. A plain fread stops there and silently returns the first 2,449 genes only
# (first run, job 3658416, 2026-10-23: median 272 UMI per cell instead of ~5,000; caught by the marker check). The reader
# below keeps only well-formed lines (awk NF == 20412), so CCDC60 is dropped (it is not an imprint gene); the checks in
# [2] stop the script before any output is written if the gene count or the per-cell totals do not match the deposit.
#
#   Rscript scripts/00_ingest/ingest_GSE227690.R
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

DATASET <- "GSE227690"
STUDY   <- "Griffin 2023 Nature"
RAW     <- file.path(GEO_RAW_DIR, "GSE227690_RAW")
COUNTS  <- file.path(RAW, "GSE227690_BM_counts.tsv.gz")
META    <- file.path(RAW, "GSE227690_BM_metadata.tsv.gz")
DONORS  <- file.path(RAW_DIR, "metadata", DATASET, "donor_ages_GSE227690_v1.0.csv")
N_CELL_ALL      <- 20411L    # cells in the deposit, incl. BM6 (2,474 Seq-Well)
N_GENE_LINES    <- 19939L    # gene lines in the deposit (counted 2026-10-23)
MALFORMED_GENES <- "CCDC60"  # the one line with the wrong field count (see the header)
N_LIB_EXPECTED  <- 13L       # 10x libraries of BM1-BM5
ANNO_COLS <- c("orig.ident2", "replicate", "tech", "Phase", "S.Score", "G2M.Score", "CellType", "bm_involvement", "is_malignant")
Y_GENES   <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y")

age_group_of <- function(a) {   # meta_v2.3 vocabulary, same breaks as build_meta_GSE227690.py
  a <- as.numeric(a)
  ifelse(a < 18, "pediatric_lt18", ifelse(a < 40, "AYA_18_39", ifelse(a < 60, "adult_40_59", "older_60plus")))
}

# ----------------------------------------------------------------------------
# [1] Donor table and author metadata ----
# ----------------------------------------------------------------------------
message("[1] donor table and author metadata")
don <- data.table::fread(DONORS, colClasses = "character")
don <- don[!grepl("Seq-Well", platform)]                              # BM6 row: not ingested (INGEST_EXCLUDE ^BM6)
stopifnot(nrow(don) == N_LIB_EXPECTED, !anyDuplicated(don$library), all(don$sex %in% c("F", "M")))
meta <- read_meta_txt(META, cell_col = 1)                             # rownames = author cell ids
stopifnot(nrow(meta) == N_CELL_ALL, all(ANNO_COLS %in% colnames(meta)))
print(table(meta$replicate, meta$tech))

# ----------------------------------------------------------------------------
# [2] Dense counts -> sparse genes x cells ----
# ----------------------------------------------------------------------------
message("[2] reading the dense count matrix (", basename(COUNTS), "), well-formed lines only")
all_genes <- data.table::fread(cmd = sprintf("zcat %s | cut -f1", shQuote(COUNTS)), header = TRUE)[[1]]
stopifnot(length(all_genes) == N_GENE_LINES)
# awk keeps the header and every line with exactly 1 + N_CELL_ALL fields; the output goes through a temp file (TMPDIR)
dt <- data.table::fread(cmd = sprintf("zcat %s | awk -F'\\t' 'NF == %d'", shQuote(COUNTS), N_CELL_ALL + 1L),
                        sep = "\t", header = TRUE, data.table = FALSE)
genes <- as.character(dt[[1]])
mat   <- as.matrix(dt[, -1, drop = FALSE]); rm(dt); gc()
rownames(mat) <- genes
M <- Matrix::Matrix(mat, sparse = TRUE); rm(mat); gc()
message("    ", nrow(M), " genes x ", ncol(M), " cells; dropped malformed: ", paste(setdiff(all_genes, rownames(M)), collapse = " "))
stopifnot(ncol(M) == N_CELL_ALL, setequal(colnames(M), rownames(meta)),
          nrow(M) == N_GENE_LINES - length(MALFORMED_GENES),
          identical(setdiff(all_genes, rownames(M)), MALFORMED_GENES),
          !anyDuplicated(rownames(M)), all(M@x == round(M@x)))        # raw integer UMI, as the deposit claims
message("    absent by deposit: MT- genes ", sum(grepl("^MT-", rownames(M))), "; XIST ", "XIST" %in% rownames(M),
        "; MALAT1 ", "MALAT1" %in% rownames(M))
d_ncount <- meta[colnames(M), "nCount_RNA"] - Matrix::colSums(M)
rel <- d_ncount / meta[colnames(M), "nCount_RNA"]
message(sprintf("    author nCount_RNA minus column sum: median %.0f (%.4f of the cell's UMI), 99th pct %.0f, range %.0f to %.0f",
                median(d_ncount), median(rel), quantile(d_ncount, 0.99), min(d_ncount), max(d_ncount)))
# a misread (as in the first run) removes most of each cell's UMI; the deposit itself differs by a few UMI per cell
stopifnot(all(d_ncount >= 0), median(rel) < 0.01,
          all(c("CD3E", "CD8A", "MS4A1", "CD14", "HBB", "RPS4Y1") %in% rownames(M)))

# ----------------------------------------------------------------------------
# [3] Replicate libraries: exclusions, and the join to the donor table ----
# ----------------------------------------------------------------------------
message("[3] replicate libraries")
rep_of <- sub("^.*-", "", colnames(M))                                # "BM1.1"
stopifnot(all(rep_of == meta[colnames(M), "replicate"]))              # cell-id suffix == the author's replicate column
reps <- ingest_keep(DATASET, sort(unique(rep_of)), "replicate")       # drops BM6.1 (Seq-Well)
lib_of_rep <- paste0(sub("\\..*$", "", reps), "_BioRep", sub("^.*\\.", "", reps))   # BM1.1 -> BM1_BioRep1
stopifnot(length(reps) == N_LIB_EXPECTED, setequal(lib_of_rep, don$library))
keep <- rep_of %in% reps
stopifnot(all(meta[colnames(M)[keep], "tech"] == "TenX"))
message("    ", length(reps), " libraries, ", sum(keep), " cells kept; dropped ", sum(!keep), " cells of ",
        paste(setdiff(unique(rep_of), reps), collapse = ", "))

# ----------------------------------------------------------------------------
# [4] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[4] building per-library Seurat objects")
obj_list <- list(); sex_list <- list()
yg <- intersect(Y_GENES, rownames(M))
for (k in seq_along(reps)) {
  r   <- reps[k]
  d   <- don[library == lib_of_rep[k]]                                # inside don[...], `library` is don's column
  key <- paste0(d$gsm, "_", d$library)                                # "GSM7105550_BM1_BioRep1"
  cells  <- colnames(M)[rep_of == r]
  counts <- M[, cells, drop = FALSE]
  sex_list[[key]] <- data.table::data.table(Sample = key, table_sex = d$sex,
                                            y_umi_per_cell = sum(counts[yg, ]) / length(cells))
  extra  <- meta[cells, ANNO_COLS, drop = FALSE]
  extra$age_years <- d$age; extra$age_group <- age_group_of(d$age)
  extra$sex <- ifelse(d$sex == "F", "female", "male")                 # the words 26_s4_pseudobulk.R expects
  message(sprintf("    - %-24s replicate %-6s donor %s  age %s  sex %s  cells %d", key, r, d$donor, d$age, d$sex, length(cells)))
  s <- make_seurat(counts = counts, sample = key, dataset = DATASET, patient = d$donor,
                   timepoint = "Baseline", study = STUDY, disease = "Healthy", extra_meta = extra)
  s$GSM     <- d$gsm        # decontX channel unit (ws/22)
  s$Library <- r            # the replicate label inside the deposit
  obj_list[[key]] <- s
}
rm(M); gc()

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
print(head(as.data.frame(smry)[, c("Sample", "Patient_ID", "Timepoint", "n_cells", "median_nFeature", "median_nCount")], 6))
stopifnot(nrow(smry) == N_LIB_EXPECTED, sum(smry$n_cells) == sum(keep), all(smry$Timepoint == "Healthy"),
          all(c("CD3E", "CD8A", "MS4A1", "CD14", "HBB", "RPS4Y1") %in% rownames(seu)),   # wrong-column read would lose these
          all(!is.na(seu$anno_CellType)))
# Y-linked UMI per cell per library vs the donor table (printed; no XIST in this deposit, see the header and prereg K3)
print(data.table::rbindlist(sex_list))
message("    Y-linked genes used: ", paste(yg, collapse = " "))

message("[ok] GSE227690 finished")
