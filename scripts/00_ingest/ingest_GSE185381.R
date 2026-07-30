# ============================================================================
# ingest_GSE185381.R   --   Lasry 2023 Nat Cancer (inflammatory AML immune atlas)
# Format : CITE-seq 10x mtx, prefixed in flat GEO_RAW folder
#          "<GSM>_<date-count-N>_{barcodes,features,matrix}.tsv/.mtx .gz"
#          plus per-library "<prefix>_metadata.csv.gz"
# RNA only: features 3rd column "Antibody Capture" rows (ADT) are dropped;
#           ADT / soupx / feature_reference / vdj_t files are NOT read.
# Metadata: the real cell id lives in the "cell" column as "<lib>:<core>"
#           (NO -1 suffix); column 1 is a stray global row index. Cells are
#           aligned to raw barcodes by 16bp core (verified overlap = 1).
# Patient (strategy C): per-cell donor column is "samples".
#           assigned cell -> Patient_ID = donor, Demuxed = TRUE
#           unassigned    -> Patient_ID = library composition string (fallback)
# Sample  : the DONOR, not the library. Libraries pool up to 5 donors and mix AML with Healthy,
#           so the library is not a biological sample; Library/GSM/Platform keep the provenance.
#           50 libraries -> 53 donors (42 AML Diagnosis + 11 Healthy); 41 donors span >1 library.
# Guard  : if a library's metadata barcodes barely map to its matrix
#          (source-file mismatch, e.g. GSM5613798), warn and use fallback.
# ============================================================================

# ----------------------------------------------------------------------------
# [0] Load config layer + shared readers (via here anchor) ----
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_qc.R"))
source(here::here("scripts", "config", "utils.R"))
source(here::here("scripts", "00_ingest", "00_common_readers.R"))

DATASET <- "GSE185381"
STUDY   <- "Lasry 2023 Nat Cancer"
RAW     <- file.path(GEO_RAW_DIR, "GSE185381_RAW")

# ----------------------------------------------------------------------------
# [1] GSM -> library donor-composition map (fallback labels) ----
# ----------------------------------------------------------------------------
# Trailing a/b/c are lanes of the same pool; kept verbatim for provenance and
# used as Patient_ID only for cells the author metadata did not demux.
SAMPLE_MAP <- c(
  GSM5613744="AML0612_AML3762_AML0160_AML0310_AML3133a", GSM5613745="AML0612_AML3762_AML0160_AML0310_AML3134b",
  GSM5613746="AML0612_AML3762_AML0160_AML0310_AML3135c", GSM5613747="Control1_Control2_Control3_Control4a",
  GSM5613748="Control1_Control2_Control3_Control4b",     GSM5613749="Control1_Control2_Control3_Control4c",
  GSM5613750="Control1_Control2_Control3_Control4d",     GSM5613751="Control5_AML2910_AML3050_AML0361a",
  GSM5613752="Control5_AML2910_AML3050_AML0361b",        GSM5613753="PAWWEEa",
  GSM5613754="PAWWEEb",                                  GSM5613755="PAWWEEc",
  GSM5613756="AML2451_Control0004",                      GSM5613757="Control0058_PAXMIJ_PAUMTZa",
  GSM5613758="Control0058_PAXMIJ_PAUMTZb",               GSM5613759="Control0058_PAXMIJ_PAUMTZc",
  GSM5613760="AML038_AML008_AML043_AML028_AML056a",      GSM5613761="AML038_AML008_AML043_AML028_AML057b",
  GSM5613762="AML073_AML006_AML025_AML003a",             GSM5613763="AML073_AML006_AML025_AML003b",
  GSM5613764="AML012_AML005_AML055_AML048a",             GSM5613765="AML012_AML005_AML055_AML048b",
  GSM5613766="AML0048-c1",                               GSM5613767="AML0048-c2",
  GSM5613768="AML3266",                                  GSM5613769="Control0082_AML052_AML022a",
  GSM5613770="Control0082_AML052_AML022b",               GSM5613771="Control0082_AML052_AML022c",
  GSM5613772="AML2123a",                                 GSM5613773="AML2123b",
  GSM5613774="Control4003a",                             GSM5613775="Control4003b",
  GSM5613776="Control4003c",                             GSM5613777="AML1371_AML4340a",
  GSM5613778="AML1371_AML4340b",                         GSM5613779="AML0024",
  GSM5613780="AML4897",                                  GSM5613781="Control0005_AML009a",
  GSM5613782="AML026_AML051",                            GSM5613783="AML001",
  GSM5613784="AML0693",                                  GSM5613785="AML3948",
  GSM5613786="AML3730",                                  GSM5613787="Control0005_AML009b",
  GSM5613788="AML0114a",                                 GSM5613789="AML0114b",
  GSM5613790="AML0114c",                                 GSM5613791="Control0004",
  GSM5613792="AML052",                                   GSM5613793="Control0082",
  GSM5613794="Control0182a",                             GSM5613795="Control0182b",
  GSM5613796="Control4003a",                             GSM5613797="Control4003b",
  GSM5613798="AML051",                                   GSM5613799="Control0005_AML009",
  GSM5613800="AML0102",                                  GSM5613801="AML0134",
  GSM5613802="AML0693",                                  GSM5613803="AML1133",
  GSM5613804="AML2910",                                  GSM5613805="AML2975",
  GSM5613806="AML3948"
)

# Per-cell donor column candidates ("samples" is the real one here; do NOT
# include orig.ident, which holds the library name).
DONOR_COL_CANDIDATES <- c(
  "samples","sample","donor","donor_id","donorID","best_hashtag","hash.ID",
  "HTO_classification","MULTI_ID","assignment","souporcell","vireo",
  "patient","patient_id","subject","individual","SNG.BEST.GUESS"
)

# ----------------------------------------------------------------------------
# [2] CITE reader: keep RNA (Gene Expression) rows only ----
# ----------------------------------------------------------------------------
read_cite_rna <- function(dir, prefix) {
  mtx <- file.path(dir, paste0(prefix, "_matrix.mtx.gz"))
  bc  <- file.path(dir, paste0(prefix, "_barcodes.tsv.gz"))
  ft  <- file.path(dir, paste0(prefix, "_features.tsv.gz"))
  stopifnot(file.exists(mtx), file.exists(bc), file.exists(ft))

  mat <- tryCatch(
    ReadMtx(mtx = mtx, cells = bc, features = ft, feature.column = 2),
    error = function(e) ReadMtx(mtx = mtx, cells = bc, features = ft, feature.column = 1)
  )

  # Drop Antibody Capture rows if the features file carries a type column.
  feats <- data.table::fread(ft, header = FALSE, data.table = FALSE)
  if (ncol(feats) >= 3 && any(feats[[3]] == "Antibody Capture")) {
    rna_idx <- which(feats[[3]] == "Gene Expression")
    mat <- mat[rna_idx, , drop = FALSE]
  }
  mat
}

# ----------------------------------------------------------------------------
# [3] Align author metadata to raw barcodes by 16bp core ----
# ----------------------------------------------------------------------------
# Returns list(aligned, n_meta, n_mapped):
#   aligned  : data.frame, rownames == full_bc, NA rows for unannotated cells
#   n_meta   : author cells in this library
#   n_mapped : how many of them land in the matrix (mismatch -> ~0)
align_lasry_meta <- function(full_bc, meta_file) {
  if (!file.exists(meta_file)) return(NULL)
  md <- data.table::fread(meta_file, sep = ",", data.table = FALSE)
  if (!"cell" %in% colnames(md)) return(NULL)

  meta_core <- sub("-\\d+$", "", sub("^.*:", "", md$cell))   # "<lib>:<core>" -> core
  drop_cols <- intersect(colnames(md), c("V1", "", "cell"))
  md_anno   <- md[, setdiff(colnames(md), drop_cols), drop = FALSE]
  rownames(md_anno) <- make.unique(meta_core)

  core_bc <- sub("-\\d+$", "", full_bc)
  idx     <- match(core_bc, rownames(md_anno))
  aligned <- md_anno[idx, , drop = FALSE]
  rownames(aligned) <- full_bc

  list(aligned = aligned, n_meta = nrow(md_anno), n_mapped = sum(!is.na(idx)))
}

# Find the first matching donor column (case-insensitive).
detect_donor_col <- function(meta_df) {
  if (is.null(meta_df) || ncol(meta_df) == 0) return(NA_character_)
  for (cand in DONOR_COL_CANDIDATES) {
    hit <- which(tolower(colnames(meta_df)) == tolower(cand))
    if (length(hit)) return(colnames(meta_df)[hit[1]])
  }
  NA_character_
}

# Rough disease state from a donor / composition label.
state_of <- function(lab) {
  has_aml  <- grepl("AML|^PA", lab)
  has_ctrl <- grepl("Control", lab, ignore.case = TRUE)
  if (has_aml && has_ctrl) "Mixed" else if (has_aml) "AML" else if (has_ctrl) "Healthy" else "Unknown"
}

# ----------------------------------------------------------------------------
# [4] Discover libraries (anchored on _matrix.mtx.gz) ----
# ----------------------------------------------------------------------------
message("[4] listing RNA libraries")
mtx_files <- list.files(RAW, pattern = "_matrix\\.mtx\\.gz$", full.names = FALSE)
prefixes  <- sub("_matrix\\.mtx\\.gz$", "", mtx_files)   # "GSM5613744_2019-07-01-count-1"
message("    ", length(prefixes), " RNA libraries found")

# ----------------------------------------------------------------------------
# [5] Build one Seurat object per library (no filtering) ----
# ----------------------------------------------------------------------------
message("[5] building per-library Seurat objects")
obj_list <- list()

for (prefix in prefixes) {
  gsm        <- stringr::str_extract(prefix, "^GSM\\d+")
  comp_label <- SAMPLE_MAP[[gsm]]
  if (is.null(comp_label) || is.na(comp_label)) comp_label <- prefix   # safety
  is_pooled  <- grepl("_", comp_label)
  message("    - ", prefix, "  (gsm=", gsm, ", composition=", comp_label, ")")

  counts  <- read_cite_rna(RAW, prefix)
  full_bc <- colnames(counts)

  meta_file <- file.path(RAW, paste0(prefix, "_metadata.csv.gz"))
  ali       <- align_lasry_meta(full_bc, meta_file)
  aligned   <- if (is.null(ali)) NULL else ali$aligned   # rownames == full_bc or NULL
  donor_col <- detect_donor_col(aligned)

  # Guard: source-file barcode mismatch (matrix vs metadata = different cells),
  # e.g. GSM5613798 (1/3132). Warning only; fallback label is used downstream.
  if (!is.null(ali) && ali$n_meta > 0 && ali$n_mapped / ali$n_meta < 0.5) {
    warning("Library ", prefix, ": only ", ali$n_mapped, "/", ali$n_meta,
            " metadata cells map to the matrix barcodes - likely a matrix/metadata ",
            "mismatch in the source files. Using fallback label '", comp_label, "'.")
  }

  s <- make_seurat(
    counts     = counts,
    sample     = prefix,        # PROVISIONAL: replaced by the donor id below
    dataset    = DATASET,
    patient    = comp_label,    # fallback default; overridden per cell below
    timepoint  = NA_character_, # cross-sectional
    study      = STUDY,
    extra_meta = aligned
  )
  # Library provenance moves into its own column now that Sample is the DONOR. A donor spans up
  # to 5 libraries in this deposit, so Library/GSM are the only remaining handle on batch.
  s$Library <- prefix
  s$GSM     <- gsm
  s$Pooled  <- is_pooled
  .plat       <- platform_of(DATASET, prefix)
  s$Platform  <- .plat$platform
  s$Chemistry <- .plat$chemistry     # 3prime vs 5prime: resolved per LIBRARY, not per dataset

  # Strategy C: per-cell demux on the donor column when available.
  if (!is.na(donor_col)) {
    donor_vals <- as.character(s@meta.data[[paste0("anno_", donor_col)]])
    has_call   <- !is.na(donor_vals) & donor_vals != "" & donor_vals != "NA"
    pid        <- ifelse(has_call, donor_vals, comp_label)
    s$Patient_ID    <- pid
    s$orig.ident    <- pid
    s$Demuxed       <- has_call
    s$Disease_state <- vapply(pid, state_of, character(1), USE.NAMES = FALSE)
    # THE SAMPLE IS THE DONOR, NOT THE LIBRARY. A library here pools up to 5 donors and mixes
    # AML with Healthy in one barcode set, so splitting on the library key (03_per_sample_qc.R
    # splits by Sample) yielded "samples" that are not biological samples: one QC unit, one CNV
    # run and one CCC graph per POOL. 50 libraries -> 53 donors (42 AML Diagnosis + 11 Healthy),
    # all 53 clearing the cell-count and complexity gates. Un-demuxed cells get a marked id so
    # they can never be mistaken for a donor; DEMUX_PREFILTER drops them before the split anyway.
    s$Sample <- ifelse(has_call, donor_vals, paste0("UNDEMUX__", prefix))
    message("        demuxed on '", donor_col, "': ",
            sum(has_call), "/", length(has_call), " cells assigned -> ",
            length(unique(donor_vals[has_call])), " donor-level samples")
  } else {
    s$Demuxed       <- FALSE
    s$Disease_state <- state_of(comp_label)
    s$Sample        <- paste0("UNDEMUX__", prefix)
    message("        no donor column -> fallback to composition label")
  }
  # make_seurat stamped uid_patient from the SCALAR fallback label, so it must be rebuilt per
  # cell after demux; otherwise every cell in a pooled library keeps the composition string.
  s$uid_patient <- paste0(DATASET, ":", s$Patient_ID)

  # C1: GSE185381 is disease-driven AND per-cell (a library pools AML+Healthy donors), so
  # canonical Timepoint must be set per cell from Disease_state, not via the sample-level
  # make_seurat path. AML->Diagnosis, Healthy->Healthy, Mixed->Unknown (DISEASE_TO_TP);
  # ignore Demuxed status (demux filtering happens later in 01_preprocess).
  s$Timepoint_detail <- s$Disease_state
  s$Timepoint        <- unname(DISEASE_TO_TP[s$Disease_state])
  s$Timepoint[is.na(s$Timepoint)] <- "Unknown"

  obj_list[[prefix]] <- s
}

# ----------------------------------------------------------------------------
# [6] Merge into one dataset object ----
# ----------------------------------------------------------------------------
message("[6] merging ", length(obj_list), " libraries")
seu <- merge_samples(obj_list)
rm(obj_list); gc()

# ----------------------------------------------------------------------------
# [7] Write QC tables + unfiltered RDS ----
# ----------------------------------------------------------------------------
message("[7] writing outputs")
write_qc_outputs(seu, DATASET)

message("[ok] GSE185381 finished")
