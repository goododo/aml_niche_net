# 00_loy_common.R ----
# input : nothing on its own; sourced after scripts/config/config_paths.R and scripts/loy/config_loy.R
# output: four functions shared by the H2 scripts (04, 05, 06)
# what it does: holds the LOY stage's shared pieces so 04-06 do not each carry a copy. 02 and 03 keep their own copies
#   of the same code (they are frozen steps; editing them would change nothing but the diff). Stage-scoped, not a
#   general utility file: everything here is about chrY counts, the three-state call and the H2 case objects.

# counts_of(): raw UMI matrix from a Seurat object (v5 layer or v4 slot); refuses ambient-corrected non-integer counts
# (Gawad / Lasry in the Zeng pack) because the Poisson model behind the -Y call needs raw UMIs.
counts_of <- function(s) {
  m <- tryCatch(SeuratObject::LayerData(s, assay = "RNA", layer = "counts"),
                error = function(e) Seurat::GetAssayData(s, assay = "RNA", slot = "counts"))
  i <- seq_len(min(1e5, length(m@x)))                       # m@x = the non-zero values of the sparse matrix
  stopifnot(all(m@x[i] == round(m@x[i])))
  m
}

# read_gse269669(): the deposited .rds.gz is a gzip of a gzip-compressed RDS; gzcon() peels the outer layer, readRDS the inner.
read_gse269669 <- function(obj = c("avm_aml", "avm_tme")) {
  obj <- match.arg(obj)
  f <- file.path(LOY_GSE269669_DIR, paste0("GSE269669_", obj, ".rds.gz")); stopifnot(file.exists(f))
  con <- gzcon(gzfile(f, "rb")); on.exit(close(con)); readRDS(con)
}

# loy_comp_of(): hierarchy bin + fine label -> compartment, exactly as 02_call_loy_status.R section [3]
loy_comp_of <- function(hierarchy_bin, bmm_fine) {
  bin  <- data.table::fifelse(is.na(hierarchy_bin) | hierarchy_bin == "", "NA", as.character(hierarchy_bin))
  fine <- data.table::fifelse(is.na(bmm_fine), "", as.character(bmm_fine))
  comp <- bin
  comp[bin %in% LOY_LEUK_BINS] <- "LEUK"
  tnk  <- bin == "T_NK"
  comp[tnk] <- data.table::fifelse(grepl(LOY_T_FINE_REGEX, fine[tnk]) & !grepl(LOY_T_FINE_EXCLUDE, fine[tnk]), "T",
                data.table::fifelse(grepl(LOY_NK_FINE_REGEX, fine[tnk]), "NK", "T_NK_other"))
  comp
}

# loy_cell_call(): the frozen three-state rule (config section 5); lambda = nCount * r(chem, comp) * factor
loy_cell_call <- function(y_umi, n_ygenes, lambda) {
  data.table::fifelse(y_umi == 0 & lambda >= LOY_LAMBDA_MIN, "suspected_LOY",
  data.table::fifelse((y_umi >= 2 & n_ygenes >= 2) | y_umi >= 3, "Y_retained", "uncertain"))
}

# loy_h2_is_leuk(): which cells count as leukaemic in an H2 library -- the depositors' own call, never ours:
# Zeng pack = NormalvsMalignant (anno_malignant); GSE269669 = membership of the authors' AML object.
loy_h2_is_leuk <- function(dataset, object, anno_malignant) {
  data.table::fifelse(dataset == "GSE269669", object == "avm_aml",
                      !is.na(anno_malignant) & tolower(as.character(anno_malignant)) == "malignant")
}
