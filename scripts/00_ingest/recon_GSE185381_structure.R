#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# recon_GSE185381_structure.R
#
# PURPOSE  Recon probe (no pipeline changes). Answers three questions before we
#          rewrite the GSE185381 sample key:
#            Q1  Does any donor (Patient_ID) span more than one library?
#            Q2  How many cells survive per donor after the split?
#            Q3  Does any library mix Healthy and AML donors (i.e. does
#                TP_DISEASE_DRIVEN actually create mixed-timepoint libraries)?
#
# USAGE    Rscript recon_GSE185381_structure.R
#          (login node is fine; only meta.data is retained per object)
#
# OUTPUT   3 CSVs + a console summary, written to OUT_DIR
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

# --- config ---------------------------------------------------------------
QC_DIR  <- "/LARGE1/gr10634/gaozy/aml_niche_net/02_seurat_objects/01_per_sample_qc/GSE185381"
OUT_DIR <- "."          # change to FAST_DIR/tables if you prefer
MIN_CELLS_DONOR   <- 500   # candidate gate: cells per donor
MIN_MED_NFEATURE  <- 500   # candidate gate: median genes per cell

message_ts <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# --- collect meta.data only ------------------------------------------------
fs <- list.files(QC_DIR, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(fs) > 0)
message_ts("[1] found ", length(fs), " library-level RDS under ", QC_DIR)

want <- c("Sample", "Patient_ID", "uid_patient", "Disease_state",
          "Timepoint", "Timepoint_detail", "Pooled", "Demuxed",
          "nCount_RNA", "nFeature_RNA", "percent.mt")

meta <- rbindlist(lapply(seq_along(fs), function(i) {
  f <- fs[i]
  o <- readRDS(f)
  m <- o@meta.data
  keep <- intersect(want, colnames(m))
  dt <- as.data.table(m[, keep, drop = FALSE])
  dt[, lib_file := basename(f)]
  rm(o, m); gc(verbose = FALSE)
  if (i %% 5 == 0) message_ts("    ... ", i, "/", length(fs))
  dt
}), fill = TRUE)

message_ts("[2] pooled meta.data: ", nrow(meta), " cells x ", ncol(meta), " cols")

# --- Q1: does a donor span multiple libraries? -----------------------------
donor_lib <- meta[, .(n_cells = .N), by = .(Patient_ID, Sample)]
donor_span <- donor_lib[, .(n_libraries = uniqueN(Sample),
                            libraries   = paste(sort(unique(Sample)), collapse = " | "),
                            n_cells     = sum(n_cells)),
                        by = Patient_ID][order(-n_libraries, -n_cells)]

message_ts("[3] Q1 donors spanning >1 library: ",
           donor_span[n_libraries > 1, .N], " / ", nrow(donor_span))
print(head(donor_span[n_libraries > 1], 20))

# --- Q2: per-donor cell counts and complexity ------------------------------
donor_qc <- meta[, .(
  n_cells        = .N,
  n_libraries    = uniqueN(Sample),
  med_nFeature   = as.numeric(median(nFeature_RNA, na.rm = TRUE)),
  med_nCount     = as.numeric(median(nCount_RNA,  na.rm = TRUE)),
  disease_states = paste(sort(unique(Disease_state)), collapse = ","),
  timepoints     = paste(sort(unique(Timepoint)),     collapse = ",")
), by = Patient_ID][order(n_cells)]

donor_qc[, pass_cells   := n_cells      >= MIN_CELLS_DONOR]
donor_qc[, pass_complex := med_nFeature >= MIN_MED_NFEATURE]
donor_qc[, pass_both    := pass_cells & pass_complex]

message_ts("[4] Q2 donors total = ", nrow(donor_qc),
           " | pass cells>=", MIN_CELLS_DONOR, ": ", donor_qc[, sum(pass_cells)],
           " | pass medFeat>=", MIN_MED_NFEATURE, ": ", donor_qc[, sum(pass_complex)],
           " | pass both: ", donor_qc[, sum(pass_both)])
print(donor_qc)

# --- Q3: does any library mix disease states / timepoints? -----------------
lib_mix <- meta[, .(
  n_cells          = .N,
  n_donors         = uniqueN(Patient_ID),
  donors           = paste(sort(unique(Patient_ID)), collapse = ","),
  n_disease_states = uniqueN(Disease_state),
  disease_states   = paste(sort(unique(Disease_state)), collapse = ","),
  n_timepoints     = uniqueN(Timepoint),
  timepoints       = paste(sort(unique(Timepoint)), collapse = ","),
  n_uid_patient    = uniqueN(uid_patient)
), by = Sample][order(-n_disease_states, -n_donors)]

message_ts("[5] Q3 libraries mixing >1 Disease_state: ",
           lib_mix[n_disease_states > 1, .N], " / ", nrow(lib_mix),
           " | mixing >1 Timepoint: ", lib_mix[n_timepoints > 1, .N])
print(lib_mix)

# --- uid_patient sanity ----------------------------------------------------
uid_check <- meta[, .(n_donors = uniqueN(Patient_ID),
                      donors   = paste(sort(unique(Patient_ID)), collapse = ",")),
                  by = uid_patient][order(-n_donors)]
message_ts("[6] uid_patient values collapsing >1 donor: ",
           uid_check[n_donors > 1, .N], " / ", nrow(uid_check))
print(head(uid_check, 20))

# --- write -----------------------------------------------------------------
fwrite(donor_span, file.path(OUT_DIR, "recon_GSE185381_donor_span.csv"))
fwrite(donor_qc,   file.path(OUT_DIR, "recon_GSE185381_donor_qc.csv"))
fwrite(lib_mix,    file.path(OUT_DIR, "recon_GSE185381_library_mix.csv"))
message_ts("[7] wrote 3 CSVs to ", normalizePath(OUT_DIR))
