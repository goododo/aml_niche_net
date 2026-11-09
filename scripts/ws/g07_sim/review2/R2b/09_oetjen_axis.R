#!/usr/bin/env Rscript
# input : rownames only of the Oetjen pseudobulk (ws_pb_GSE120221.rds) and of the GSE185381 T pseudobulk; S5 gene column;
#         ambient table; the GTF lncRNA rule of draft 2.3. No count value is read.
# output: 09_oetjen_axis.out
# what it does: how many U_A genes (and U_A lncRNA genes) are absent from Oetjen's gene axis ("not assessable", 3.4)
suppressPackageStartupMessages(library(data.table))
O <- readRDS("/LARGE1/gr10634/gaozy/aml_niche_net/ws_pseudobulk/GSE120221/ws_pb_GSE120221.rds")
cat("Oetjen object names:", paste(names(O), collapse = ", "), "\n")
go <- if (!is.null(O$raw$T)) rownames(O$raw$T) else rownames(O$raw[[1]]); rm(O); invisible(gc())
U <- fread("/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/ws_s5_de_T.csv", select = "gene")$gene
ext <- fread("/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/decontx/GSE185381/ws_s1_ambient_blacklist.csv")[n_units >= 2, gene]
UA <- setdiff(U, ext)
gtf <- "/LARGE1/gr10634/gaozy/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
if (!file.exists(gtf)) gtf <- Sys.glob("/LARGE1/gr10634/gaozy/**/refdata-gex-GRCh38-2020-A/genes/genes.gtf")[1]
cat("GTF:", gtf, "\n")
L <- if (!is.na(gtf) && file.exists(gtf)) { x <- fread(cmd = sprintf("awk -F'\\t' '$3==\"gene\"' %s", gtf), header = FALSE, sep = "\t")$V9
  data.table(name = sub('.*gene_name "([^"]+)".*', "\\1", x), type = sub('.*gene_type "([^"]+)".*', "\\1", x)) } else NULL
lnc <- if (!is.null(L)) L[, .(only = all(type == "lncRNA")), by = name][only == TRUE, name] else character()
cat(sprintf("Oetjen axis genes %d | U_A %d, absent from Oetjen %d (%.1f %%) | U_A lncRNA %d, absent %d (%.1f %%) | U_A non-lncRNA absent %d\n",
    length(go), length(UA), sum(!UA %in% go), 100 * mean(!UA %in% go), sum(UA %in% lnc), sum(UA %in% lnc & !UA %in% go),
    100 * mean(!UA[UA %in% lnc] %in% go), sum(!UA %in% lnc & !UA %in% go)))
cat("examples absent:", head(setdiff(UA, go), 25), "\n")
