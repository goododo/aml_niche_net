# input:  project MP table, curated stemness table, van Galen Table S3, pinned HALLMARK GMT
# output: scripts/ws/gene_sets_clonal_env_v1.tsv (set, role, gene, weight, source), read by the CLONAL_ENV scripts
# what it does: assembles every gene list the clonal-vs-environment prereg names, before any expression is read (scripts/ws/PREREGISTRATION_CLONAL_ENV.md)
suppressPackageStartupMessages({ library(data.table); library(readxl) })
R  <- "/FAST/gr10634/gaozy/aml_niche_net"
OUT <- file.path(R, "scripts/ws/gene_sets_clonal_env_v1.tsv")
add <- function(set, role, genes, src, w = NA_real_) data.table(set, role, gene = genes, weight = w, source = src)
mp <- fread(file.path(R, "results/tables/04_cnmf/malignant_metaprograms.tsv"))
L <- lapply(c("MP5", "MP6", "MP10"), function(m) add(m, "primary", mp[MP == m, gene],
            "this project cNMF malignant_metaprograms.tsv (2026-07-17)"))
st <- fread(file.path(R, "scripts/以前06_hierarchy/stemness_signatures.tsv"), skip = "signature")
lab <- c(LSC17 = "F1_LSC17", vanGalen_HSC_Prog = "F2_vanGalen_HSC_Prog",
         vanGalen_HSCprog_like = "F3_vanGalen_HSCprog_like", HSPC_core = "F4_HSPC_core")
for (s in names(lab)) L[[length(L) + 1]] <- add(lab[[s]], "panel", st[signature == s, gene],
         "stemness_signatures.tsv (Ng 2016 / van Galen 2019 Table S3 / curated HSPC)", st[signature == s, as.numeric(weight)])
g <- strsplit(readLines("/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt"), "\t")
names(g) <- vapply(g, `[`, "", 1)                     # element 1 = set name, 3+ = genes
hm <- c(F5_IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE", F6_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        F7_TNFA_NFKB = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
for (s in names(hm)) L[[length(L) + 1]] <- add(s, "panel", g[[hm[[s]]]][-(1:2)], paste("MSigDB", hm[[s]], "v2024.1 (pinned GMT)"))
t3 <- as.data.frame(read_excel(file.path(R, "scripts/03_hierarchy/van Galen 2019 TableS3.xlsx"), col_names = FALSE, .name_repair = "minimal"))
j <- which(t3[2, ] == "Monocyte-like")                # row 2 holds the column labels
mono <- na.omit(t3[-(1:2), j]); mono <- mono[mono != ""]
L[[length(L) + 1]] <- add("F9_vanGalen_Monocyte_like", "panel", mono, "van Galen 2019 Table S3, tumor-derived per cell type, Monocyte-like")
L[[length(L) + 1]] <- add("F10_inhibitory_ligands", "panel", c("LILRB4", "CD274", "PDCD1LG2", "LGALS9", "CD47", "VSIR"), "curated (prereg draft v0)")
L[[length(L) + 1]] <- add("K1_HOX", "known_answer", c("HOXA5", "HOXA7", "HOXA9", "HOXA10", "HOXB3", "HOXB5", "MEIS1", "PBX3"), "curated NPM1/HOX (prereg draft v0)")
D <- rbindlist(L)
stopifnot(!anyDuplicated(D[, .(set, gene)]))           # sanity: no gene listed twice in one set
fwrite(D, OUT, sep = "\t")
print(D[, .(n_genes = .N, weighted = any(!is.na(weight))), by = .(set, role)])
cat("md5:", tools::md5sum(OUT), "\n")
