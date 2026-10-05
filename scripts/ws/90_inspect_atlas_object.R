#!/usr/bin/env Rscript
# 90_inspect_atlas_object.R ----
# INPUT  : /LARGE1/gr10634/gaozy/leukemia_ml/03.2_anno/inte_rough_anno.RDS  (42.5 GB)
# OUTPUT : stdout only -- what the object still contains
# WHAT IT DOES : the leukemia_ml integration inputs were deleted from LARGE0, so this object is the
#          only surviving copy of the QC'd cells for 161 samples including BoneMarrowMap. Whether a
#          re-integration is possible at all depends on whether RAW COUNTS survive in it, or only
#          the integrated/scaled values. Read-only; decides nothing.
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
f <- "/LARGE1/gr10634/gaozy/leukemia_ml/03.2_anno/inte_rough_anno.RDS"
message("[1] reading ", f, " ...")
s <- readRDS(f)
message(sprintf("[2] class %s | cells %d | default assay %s", class(s)[1], ncol(s), DefaultAssay(s)))
for (a in names(s@assays)) {
  A <- s[[a]]
  ly <- tryCatch(SeuratObject::Layers(A), error = function(e) "(v3 assay)")
  message(sprintf("    assay %-12s class %-10s features %6d | layers/slots: %s",
                  a, class(A)[1], nrow(A), paste(ly, collapse = ", ")))
  if (inherits(A, "Assay")) {
    for (sl in c("counts", "data", "scale.data")) {
      d <- tryCatch(dim(methods::slot(A, sl)), error = function(e) NULL)
      message(sprintf("      slot %-11s %s", sl, if (is.null(d)) "absent" else paste(d, collapse = " x ")))
    }
  }
}
message(sprintf("[3] reductions: %s", paste(names(s@reductions), collapse = ", ")))
message(sprintf("[4] graphs: %s", paste(names(s@graphs), collapse = ", ")))
md <- s@meta.data
message("[5] metadata columns:"); print(colnames(md))
message("\n[6] cells per dataset")
dsc <- grep("Dataset|dataset|orig.ident", colnames(md), value = TRUE)[1]
print(as.data.table(md)[, .N, by = c(dsc)][order(-N)])
message("\n[7] annotation columns, level counts")
for (k in grep("anno|Anno|CellType|cluster|snn", colnames(md), value = TRUE))
  message(sprintf("    %-28s %d levels", k, length(unique(md[[k]]))))
message("\n[8] is a raw counts matrix usable for re-integration? ")
ok <- FALSE
for (a in names(s@assays)) {
  A <- s[[a]]
  cnt <- tryCatch(SeuratObject::LayerData(s, assay = a, layer = "counts"), error = function(e) NULL)
  if (!is.null(cnt) && nrow(cnt) > 0 && max(cnt[, 1:min(50, ncol(cnt))]) %% 1 == 0) {
    message(sprintf("    YES -- assay %s has integer counts (%d x %d)", a, nrow(cnt), ncol(cnt))); ok <- TRUE }
}
if (!ok) message("    NO integer counts found -- a re-integration would need the data re-downloaded")
message("[done]")
