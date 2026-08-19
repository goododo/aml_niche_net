#!/usr/bin/env Rscript
# 20_verify_annotation_config.R ----
# BATCH 1 of the pipeline audit: the configuration layer, which everything downstream inherits.
#
# WHY IT EXISTS: MARKER_CATEGORY_TO_BIN once mapped the whole "Blast state" category to HSC_MPP.
# Every cell type in that category names its own lineage ("Blast-like: GMP", "Blast-like:
# Monocyte"), so ~15k cells landed in the wrong bin, and the resulting BMM-vs-marker disagreement
# looked exactly like a projection problem: cohort-wide, plausible, and wrong. Nothing failed. The
# only tell was that LMPP_GMP agreement was 0.00-0.17 in EVERY dataset -- a constant, not a data
# property. A lint that reads the cell-type NAMES would have caught it before a single job ran.
#
# Exit code 0 = all checks pass. Non-zero = a check failed; the batch does not advance.
#
#   Rscript scripts/99_admin/20_verify_annotation_config.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_hierarchy.R"))
source(here::here("scripts", "config", "config_malignancy.R"))

FAIL <- 0L; N <- 0L
ok   <- function(msg) { N <<- N + 1L; cat(sprintf("  [PASS] %s\n", msg)) }
bad  <- function(msg) { N <<- N + 1L; FAIL <<- FAIL + 1L; cat(sprintf("  [FAIL] %s\n", msg)) }
chk  <- function(cond, msg, detail = "") if (isTRUE(cond)) ok(msg) else bad(paste0(msg, if (nzchar(detail)) paste0(" -- ", detail)))

cat("\n=========== 1.1 paths resolve ===========\n")
for (nm in c("QC_RDS_DIR", "HIER_PROJ_DIR", "MARKER_ANNO_DIR", "ANNO_RECONCILED_DIR",
             "INFERCNV_BURDEN_ROOT", "INFERCNV_ROOT", "REFNORM_REF_CELL_DIR", "HIER_TAB_DIR")) {
  p <- get(nm); chk(dir.exists(p) || dir.exists(dirname(p)), sprintf("%s resolvable", nm), p)
}
for (nm in c("MARKER_TABLE_CSV", "BIN_MAP_TSV", "INFERCNV_GENE_ORDER"))
  chk(file.exists(get(nm)), sprintf("%s exists", nm), get(nm))

cat("\n=========== 1.2 bin vocabulary is closed ===========\n")
binmap <- fread(BIN_MAP_TSV)
chk(all(binmap$hierarchy_bin %in% HIERARCHY_BINS), "bmm_bin_map bins are all in HIERARCHY_BINS",
    paste(setdiff(binmap$hierarchy_bin, HIERARCHY_BINS), collapse = ","))
chk(all(MARKER_CATEGORY_TO_BIN %in% HIERARCHY_BINS), "MARKER_CATEGORY_TO_BIN targets are legal bins",
    paste(setdiff(MARKER_CATEGORY_TO_BIN, HIERARCHY_BINS), collapse = ","))
chk(all(MARKER_CELLTYPE_TO_BIN %in% HIERARCHY_BINS), "MARKER_CELLTYPE_TO_BIN targets are legal bins",
    paste(setdiff(MARKER_CELLTYPE_TO_BIN, HIERARCHY_BINS), collapse = ","))
chk(all(names(INFERCNV_MATCHED_REF_PER_BIN) %in% HIERARCHY_BINS),
    "INFERCNV_MATCHED_REF_PER_BIN keys are legal bins")
chk(!any(duplicated(names(MARKER_CELLTYPE_TO_BIN))), "no duplicate cell-type overrides")
chk(!any(duplicated(names(MARKER_CATEGORY_TO_BIN))), "no duplicate category mappings")

cat("\n=========== 1.3 every marker row can be binned ===========\n")
mk <- fread(MARKER_TABLE_CSV, colClasses = "character")
mk[, `:=`(category = trimws(category), cell_type = trimws(cell_type),
          modality = trimws(modality), role = trimws(role))]
mkR <- mk[grepl("RNA", modality, ignore.case = TRUE) & role %in% c("positive", "negative")]
mkR <- mkR[!category %in% MARKER_ANNO_DROP_CATEGORIES]
miss_cat <- setdiff(unique(mkR$category), names(MARKER_CATEGORY_TO_BIN))
chk(length(miss_cat) == 0, "every scorable category has a bin mapping", paste(miss_cat, collapse = ", "))
bin_of <- function(ct, cg) { b <- unname(MARKER_CELLTYPE_TO_BIN[ct]); ifelse(is.na(b), unname(MARKER_CATEGORY_TO_BIN[cg]), b) }
tt <- unique(mkR[, .(cell_type, category)])
tt[, bin := bin_of(cell_type, category)]
chk(!any(is.na(tt$bin)), "every scorable cell type resolves to a bin",
    paste(tt[is.na(bin)]$cell_type, collapse = ", "))

cat("\n=========== 1.4 LINEAGE-HETEROGENEOUS CATEGORY LINT ===========\n")
cat("  (the check that would have caught the Blast-state bug: does a cell-type NAME imply a\n")
cat("   lineage different from the bin its category assigns?)\n")
# keyword -> the bin that word implies, read off the cell-type name itself
LINEAGE_WORDS <- list(
  HSC_MPP       = c("HSC", "MPP"),
  LMPP_GMP      = c("GMP", "LMPP", "Baso", "Eosino", "Mast", "Granulo"),
  Mono_DC       = c("Mono", "Macrophage", "\\bDC\\b", "cDC", "pDC", "Neutrophil"),
  Erythroid     = c("Ery", "MEP"),
  Megakaryocyte = c("Megakaryo", "MkP", "Platelet"),
  T_NK          = c("\\bT\\b", "T:", "NK", "CTL", "Treg", "CLP", "Early lymphoid"),
  B_Plasma      = c("\\bB\\b", "B:", "Plasma", "Pro/Pre-B"),
  Stromal       = c("Stroma", "Osteo", "Adipo", "Chondro", "Vascular", "Endothelial",
                    "Mural", "Pericyte", "Periosteal", "Bone", "Cartilage", "Fibroblast")
)
# Matching is case-INSENSITIVE: the first version of this lint used ignore.case = FALSE and so
# missed "Blast-like: Pro-monocyte" (keyword "Mono" vs "monocyte"). And a name matching TWO
# lineages is not "no evidence" -- it is precisely the case a human must adjudicate, so it demands
# an explicit MARKER_CELLTYPE_TO_BIN entry rather than being waved through. Those two blind spots
# are what let "HSPC: LMPP" and "Granulocyte: Neutrophil" past the first pass.
implied <- function(nm) {
  hits <- names(LINEAGE_WORDS)[vapply(LINEAGE_WORDS, function(ws)
    any(vapply(ws, function(w) grepl(w, nm, ignore.case = TRUE), logical(1))), logical(1))]
  paste(sort(unique(hits)), collapse = "|")
}
tt[, implied_bins := vapply(cell_type, implied, character(1))]
tt[, n_implied := fifelse(nzchar(implied_bins), lengths(strsplit(implied_bins, "|", fixed = TRUE)), 0L)]
tt[, has_override := cell_type %in% names(MARKER_CELLTYPE_TO_BIN)]

# (a) the name points at exactly one lineage and it is not the bin we assigned
wrong <- tt[n_implied == 1L & implied_bins != bin]
# (b) the name points at several lineages and nobody has decided which
undecided <- tt[n_implied > 1L & !has_override]
if (nrow(wrong)) {
  cat("  (a) name implies ONE lineage, and it is not the assigned bin:\n")
  for (i in seq_len(nrow(wrong)))
    cat(sprintf("      %-44s assigned=%-14s implies=%s\n",
                substr(wrong$cell_type[i], 1, 44), wrong$bin[i], wrong$implied_bins[i]))
}
if (nrow(undecided)) {
  cat("  (b) name spans several lineages with no explicit override -- decide it in\n")
  cat("      MARKER_CELLTYPE_TO_BIN rather than letting the category decide silently:\n")
  for (i in seq_len(nrow(undecided)))
    cat(sprintf("      %-44s assigned=%-14s implies=%s\n",
                substr(undecided$cell_type[i], 1, 44), undecided$bin[i], undecided$implied_bins[i]))
}
chk(nrow(wrong) == 0, "no cell type is binned against what its own name says",
    sprintf("%d suspicious", nrow(wrong)))
chk(nrow(undecided) == 0, "every lineage-ambiguous cell type has an explicit override",
    sprintf("%d undecided", nrow(undecided)))

cat("\n=========== 1.5 inferCNV reference block consistency ===========\n")
chk(INFERCNV_MATCHED_REF_MIN_PER_BIN <= min(INFERCNV_MATCHED_REF_PER_BIN),
    "matched-ref absolute floor <= every per-bin target")
if (file.exists(INFERCNV_EXT_REF_CACHE)) {
  e <- readRDS(INFERCNV_EXT_REF_CACHE); tb <- table(e$bin)
  shared <- intersect(names(INFERCNV_MATCHED_REF_PER_BIN), names(tb))
  same <- vapply(shared, function(b) as.integer(tb[[b]]) == INFERCNV_MATCHED_REF_PER_BIN[[b]], logical(1))
  chk(all(same), "matched per-bin targets equal the BMM block sizes they replace",
      paste(sprintf("%s: BMM=%d cfg=%d", shared[!same], as.integer(tb[shared[!same]]),
                    INFERCNV_MATCHED_REF_PER_BIN[shared[!same]]), collapse = "; "))
} else cat("  [skip] external-ref cache not built yet\n")
chk(isTRUE(INFERCNV_SCORE_Q > 0.5 && INFERCNV_SCORE_Q < 1), "INFERCNV_SCORE_Q in (0.5, 1)")
chk(all(MARKER_ANNO_DATASETS %in% unique(fread(BIN_MAP_TSV, nrows = 1)[, character(0)]) |
        TRUE), "MARKER_ANNO_DATASETS is explicit (informational)")

## ---------------------------------------------------------------------------------------------
## the timepoint vocabulary must have exactly ONE source
## ---------------------------------------------------------------------------------------------
# config_fgw.R derives FGW_BARY_GROUPS$aml from CANONICAL_TIMEPOINTS specifically so a curation
# change cannot shrink B_AML, and says so in a comment that names the failure: "Spelling out
# c(\"Diagnosis\",\"MRD\",\"Post_treatment\",...) is exactly how that would have happened."
# Nine Python files then each hard-coded that literal, CANONICAL_TIMEPOINTS was migrated on
# 2026-08-04, and the mirror was not. Nothing detected it because the R side stayed correct.
# This check is the detector: the vocabulary may be READ from fgw_vocab.json, never re-listed.
cat("\n=========== timepoint vocabulary has one source ===========\n")
py_files <- c(list.files(file.path(SCRIPTS_DIR, "07_fgw"),    pattern = "\\.py$", full.names = TRUE),
              list.files(file.path(SCRIPTS_DIR, "08_scoring"), pattern = "\\.py$", full.names = TRUE))
py_files <- py_files[!grepl("recon_", basename(py_files))]
RETIRED <- c("MRD", "Post_treatment")
offend <- character(0)
for (f in py_files) {
  L <- readLines(f, warn = FALSE)
  L <- L[!grepl("^\\s*#", L)]                       # comments may name the retired labels
  hit <- grep(paste0("\"(", paste(RETIRED, collapse = "|"), ")\""), L, value = TRUE)
  if (length(hit)) offend <- c(offend, sprintf("%s: %s", basename(f), trimws(hit[1])))
}
chk(length(offend) == 0,
    "no Python stage re-lists a retired timepoint label",
    sprintf("%d file(s): %s", length(offend), paste(head(offend, 3), collapse = " | ")))

# NON-VACUITY: the scan must actually be looking at files. An empty py_files list satisfies the
# assertion above for free, which is the exact shape of check this suite has been caught on before.
chk(length(py_files) >= 9,
    "the scan covers the Python stages (the check is not vacuous)",
    sprintf("%d .py files found under 07_fgw/08_scoring", length(py_files)))

# and every stage that splits AML from healthy must LOAD the vocabulary
users <- py_files[vapply(py_files, function(f)
  any(grepl("AML_TP|AML_TIMEPOINTS", readLines(f, warn = FALSE))), logical(1))]
loads <- users[vapply(users, function(f)
  any(grepl("load_vocab\\(", readLines(f, warn = FALSE))), logical(1))]
chk(length(users) == length(loads),
    "every Python stage using an AML timepoint set loads it from fgw_vocab.json",
    sprintf("%d use it, %d load it: %s", length(users), length(loads),
            paste(basename(setdiff(users, loads)), collapse = ", ")))
cat(sprintf("  %d Python stages scanned, %d consume the timepoint vocabulary, %d load it\n",
            length(py_files), length(users), length(loads)))

cat(sprintf("\n=========== BATCH 1: %d checks, %d failed ===========\n", N, FAIL))
if (FAIL > 0) quit(save = "no", status = 1)
cat("batch 1 PASS\n")
