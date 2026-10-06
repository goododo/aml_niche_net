#!/usr/bin/env Rscript
# INPUT  : /LARGE1/.../06_reintegration/01_input/<Dataset_ID>.rds  (9 objects)
#          results/tables/ws/ws_reintegration_common_genes.csv     (32,042-gene intersection)
# OUTPUT : results/tables/ws/ws_preflight_integration.csv          per-dataset pre-flight numbers
# what it does: answers the three questions that can silently break or bias the RPCA step --
#   (A) does vars.to.regress have a usable percent.mt in every dataset?   -> hard blocker
#   (B) do the cell-cycle gene sets survive the gene intersection?        -> hard blocker
#   (C) is residual doublet content comparable across datasets?           -> bias, not a blocker
# (C) is measured as co-expression ENRICHMENT over independence, obs/(p1*p2), because the raw
# co-expression rate rises with sequencing depth and would confound the comparison.

suppressPackageStartupMessages({ library(data.table); library(Seurat); library(Matrix) })
IN   <- "/LARGE1/gr10634/gaozy/leukemia_ml/06_reintegration/01_input"
TAB  <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
GENES <- fread(file.path(TAB, "ws_reintegration_common_genes.csv"))$gene

# (B) cell-cycle sets, checked once against the intersection
S_G <- Seurat::cc.genes$s.genes; G2_G <- Seurat::cc.genes$g2m.genes
message(sprintf("[B] cell cycle: s.genes %d/%d in intersection | g2m.genes %d/%d",
                sum(S_G %in% GENES), length(S_G), sum(G2_G %in% GENES), length(G2_G)))
message("    missing s.genes  : ", paste(setdiff(S_G, GENES), collapse = ", "))
message("    missing g2m.genes: ", paste(setdiff(G2_G, GENES), collapse = ", "))

# (C) mutually exclusive lineage pairs: a doublet is the only ordinary way to get both
PAIRS <- list(c("CD3E", "LYZ"), c("CD3E", "MS4A1"), c("CD3E", "MPO"), c("MS4A1", "LYZ"))
stopifnot(all(unlist(PAIRS) %in% GENES))

OUT <- rbindlist(lapply(sort(list.files(IN, "\\.rds$")), function(f) {
  ds <- sub("\\.rds$", "", f); o <- readRDS(file.path(IN, f))
  md <- o@meta.data
  cnt <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
  n <- ncol(cnt)
  r <- data.table(
    Dataset_ID = ds, cells = n,
    genes_in_file = nrow(cnt),
    genes_lost_to_intersection = nrow(cnt) - length(intersect(rownames(cnt), GENES)),
    # (A)
    pmt_present = "percent.mt" %in% colnames(md),
    pmt_na      = if ("percent.mt" %in% colnames(md)) sum(is.na(md$percent.mt)) else NA_integer_,
    pmt_median  = if ("percent.mt" %in% colnames(md)) median(md$percent.mt, na.rm = TRUE) else NA_real_,
    med_ncount  = median(md$nCount_RNA), med_nfeat = median(md$nFeature_RNA))
  # (C)
  for (p in PAIRS) {
    a <- cnt[p[1], ] > 0; b <- cnt[p[2], ] > 0
    obs <- mean(a & b); exp <- mean(a) * mean(b)
    r[[paste0("coex_", p[1], "_", p[2])]]    <- round(100 * obs, 3)
    r[[paste0("enrich_", p[1], "_", p[2])]] <- if (exp > 0) round(obs / exp, 3) else NA_real_
  }
  rm(o, cnt); gc(FALSE)
  message(sprintf("  %-15s %7d cells | percent.mt NA %s | CD3E&LYZ enrich %.2f",
                  ds, n, r$pmt_na, r$enrich_CD3E_LYZ))
  r
}), fill = TRUE)

fwrite(OUT, file.path(TAB, "ws_preflight_integration.csv"))
message("\n[A] datasets with an unusable percent.mt (absent or any NA): ",
        OUT[pmt_present == FALSE | pmt_na > 0, paste(Dataset_ID, collapse = ", ")],
        if (nrow(OUT[pmt_present == FALSE | pmt_na > 0]) == 0L) "(none) -> PASS" else "  ** BLOCKER **")
print(OUT[, .(Dataset_ID, cells, pmt_na, pmt_median, med_ncount,
              enrich_CD3E_LYZ, enrich_CD3E_MS4A1, enrich_CD3E_MPO, enrich_MS4A1_LYZ)])
message("[done] ", file.path(TAB, "ws_preflight_integration.csv"))
