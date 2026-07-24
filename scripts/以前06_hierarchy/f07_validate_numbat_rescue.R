#!/usr/bin/env Rscript
# f07_validate_numbat_rescue.R ----
# Sanity check for a Numbat run that was rescued by loosening max_entropy (or any re-tuned run):
# does "tumor" concentrate on myeloid/primitive bins (expected, real signal) or bleed into T_NK/
# B_Plasma (red flag -- the looser threshold is calling noise as CNV)?
#
# NOT a substitute for f05 (f05 = BMM projection-mapping quality; this = CNV-caller malignancy
# quality). The two are orthogonal QC axes.
#
# Reads the Numbat percell csv directly (full schema, has compartment_opt) + the d20 binned csv,
# joins on the 16bp barcode core, reports tumor-fraction by hierarchy_bin.
#
# Usage:
#   Rscript f07_validate_numbat_rescue.R --dataset=GSE227903 --sample=6323_MRD

suppressPackageStartupMessages({ library(optparse); library(data.table) })
source("d00_config.R")

PERCELL_BINNED_DIR <- file.path(FAST_DIR, "results/tables/04_hierarchy/percell_binned")

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dataset", type = "character", help = "dataset id"),
  make_option("--sample",  type = "character", help = "sample id"),
  make_option("--numbat_dir", type = "character", default = NULL,
              help = "override dir (default: NUMBAT_PERCELL_DIR/<ds>)")
)))
if (is.null(opt$dataset) || is.null(opt$sample)) stop("--dataset and --sample are required")

nb_dir <- if (!is.null(opt$numbat_dir)) opt$numbat_dir else file.path(NUMBAT_PERCELL_DIR, opt$dataset)
nb_f   <- file.path(nb_dir, opt$sample, "numbat", paste0(opt$sample, "__numbat_percell.csv"))
bn_f   <- file.path(PERCELL_BINNED_DIR, opt$dataset, paste0(opt$sample, "_binned.csv"))
if (!file.exists(nb_f)) stop("[f07] no numbat percell at ", nb_f)
if (!file.exists(bn_f)) stop("[f07] no binned csv at ", bn_f)

core16 <- function(x) sub("-\\d+$", "", sub("^.*_", "", x))

nb <- fread(nb_f)
if (!"compartment_opt" %in% names(nb)) stop("[f07] ", opt$sample, " is degraded-schema Numbat (no compartment_opt) -- nothing to validate")
nb[, core := core16(cell)]
nb[, tumor := as.integer(compartment_opt == "tumor")]

bn <- fread(bn_f, select = c("cell", "hierarchy_bin"))
bn[, core := core16(cell)]

m <- merge(bn, nb[, .(core, tumor, clone_opt, GT_opt)], by = "core")
message("[1] joined ", nrow(m), " / ", nrow(bn), " binned cells to Numbat calls (coverage ",
        round(100 * nrow(m) / nrow(bn), 1), "%)")

by_bin <- m[, .(n_cells = .N, n_tumor = sum(tumor), tumor_frac = round(mean(tumor), 4)),
           by = hierarchy_bin]
setorder(by_bin, -tumor_frac)
message("\n== tumor fraction by hierarchy_bin ==")
print(by_bin)

tnk_frac <- by_bin[hierarchy_bin == "T_NK", tumor_frac]
mye_frac <- by_bin[hierarchy_bin %in% c("HSC_MPP", "LMPP_GMP"), sum(n_tumor) / sum(n_cells)]
message("\nT_NK tumor_frac = ", ifelse(length(tnk_frac), tnk_frac, "NA (bin absent/too small)"),
        "   myeloid(HSC_MPP+LMPP_GMP) tumor_frac = ", round(mye_frac, 4))
message(if (length(tnk_frac) && !is.na(tnk_frac) && tnk_frac > 0.10)
  "[!] WARNING: T_NK tumor_frac > 10% -- looser max_entropy may be calling noise as CNV in normal cells."
  else "[ok] T_NK tumor_frac looks low -- signal concentrates on myeloid bins, consistent with real malignancy.")

# clone composition, for a quick look at subclonal structure
message("\n== clone_opt x GT_opt (segment) composition ==")
print(m[, .N, by = .(clone_opt, GT_opt)][order(clone_opt)])
