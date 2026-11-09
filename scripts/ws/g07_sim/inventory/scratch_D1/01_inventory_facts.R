#!/usr/bin/env Rscript
# INPUT : gene columns ONLY of results/tables/ws/ws_s5_de_{T,B,NK,Mono}.csv (no logFC/P read);
#         HALLMARK GMT; ws_s1_ambient_blacklist.csv; Savage delay list; config_ws.R WS_B0_IEG;
#         ws_pb_GSE185381.rds ($cells: cell counts only); ws_s4_pseudobulk_meta_GSE185381.csv (covariate layout only)
# OUTPUT: scratch_D1/01_inventory_facts.txt, scratch_D1/01_hallmark_coverage.csv
# what it does: outcome-independent inventory (gene-universe sizes, set coverage, cell counts). Blind-safe:
#   never reads any statistic of the S5 contrast and never touches the 370/402 gene lists.
suppressPackageStartupMessages({ library(data.table) })
REPO <- "/FAST/gr10634/gaozy/aml_niche_net"; TABWS <- file.path(REPO, "results/tables/ws")
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/inventory/scratch_D1"
sink(file.path(OUT, "01_inventory_facts.txt"), split = TRUE)
U <- list()
for (cp in c("T", "B", "NK", "Mono")) U[[cp]] <- fread(file.path(TABWS, paste0("ws_s5_de_", cp, ".csv")), select = "gene")$gene
cat("tested genes per compartment:", paste(names(U), lengths(U), collapse = "; "), "\n")
cat("T∩B", length(intersect(U$T, U$B)), " T∩NK", length(intersect(U$T, U$NK)), " T∩Mono", length(intersect(U$T, U$Mono)),
    " T∩B∩NK", length(Reduce(intersect, U[c("T","B","NK")])), "\n")
# gmt
gl <- strsplit(readLines("/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt"), "\t")
H <- setNames(lapply(gl, function(x) unique(x[-(1:2)])), vapply(gl, `[`, "", 1))
cat("HALLMARK sets:", length(H), "; md5", tools::md5sum("/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt"), "\n")
BL <- fread(file.path(TABWS, "decontx/GSE185381/ws_s1_ambient_blacklist.csv"))
EXT <- BL[n_units >= 2, gene]; B5 <- BL[n_units >= 5, gene]
cat("ambient blacklist rows", nrow(BL), "; n_units>=5:", length(B5), "; extended (n_units>=2):", length(EXT),
    "; extended in T universe:", sum(EXT %in% U$T), "\n")
SV <- fread("/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv")
cat("Savage md5", tools::md5sum("/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv"),
    "; genes", nrow(SV), "; up", SV[direction=="up", .N], "; down", SV[direction=="down", .N],
    "; in T universe up", sum(SV[direction=="up", gene] %in% U$T), "down", sum(SV[direction=="down", gene] %in% U$T), "\n")
src <- readLines(file.path(REPO, "scripts/ws/config_ws.R")); i0 <- grep("^WS_B0_IEG <-", src)
eval(parse(text = paste(src[i0:(i0+2)], collapse = "\n")))
cat("IEG list", length(WS_B0_IEG), "; in T universe", sum(WS_B0_IEG %in% U$T), "\n")
UT <- U$T; UTx <- setdiff(UT, EXT)
cat("T universe minus extended ambient:", length(UTx), "\n")
cov <- rbindlist(lapply(names(H), function(s) { g <- H[[s]]
  data.table(set = s, size = length(g), inT = sum(g %in% UT), inT_noAmb = sum(g %in% UTx),
             inT_noAmb_noSavage = sum(g %in% setdiff(UTx, SV$gene)), inT_noAmb_noSavage_noIEG = sum(g %in% setdiff(UTx, c(SV$gene, WS_B0_IEG))),
             savage_in = sum(g %in% intersect(UTx, SV$gene)), ieg_in = sum(g %in% intersect(UTx, WS_B0_IEG)),
             inTB = sum(g %in% intersect(UTx, U$B)), inTNK = sum(g %in% intersect(UTx, U$NK)), inTMono = sum(g %in% intersect(UTx, U$Mono))) }))
cov[, frac_handling := round((savage_in + ieg_in - sum(H[[set]] %in% intersect(intersect(UTx, SV$gene), WS_B0_IEG))) / inT_noAmb, 3), by = set]
fwrite(cov[order(-inT_noAmb)], file.path(OUT, "01_hallmark_coverage.csv"))
cat("covered size (T universe minus extended ambient): min", min(cov$inT_noAmb), "q25", quantile(cov$inT_noAmb, .25), "median", median(cov$inT_noAmb),
    "q75", quantile(cov$inT_noAmb, .75), "max", max(cov$inT_noAmb), "\n")
cat("sets with >= 15 covered genes:", sum(cov$inT_noAmb >= 15), "; >= 10:", sum(cov$inT_noAmb >= 10), "\n")
cat("sets with >= 15 covered after Savage+IEG removal:", sum(cov$inT_noAmb_noSavage_noIEG >= 15), "\n")
cat("distinct genes in union of HALLMARK covered:", length(intersect(unique(unlist(H)), UTx)), "of", length(UTx), "\n")
cat("sets with handling-gene fraction >= 0.25:", paste(cov[frac_handling >= 0.25, sprintf("%s(%.2f)", set, frac_handling)], collapse = ", "), "\n")
cat("sets with handling-gene fraction >= 0.10:", cov[frac_handling >= 0.10, .N], "\n")
# pseudobulk structure (no expression)
PB <- readRDS(file.path("/LARGE1/gr10634/gaozy/aml_niche_net/ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds"))
cat("PB names:", paste(names(PB), collapse = ", "), "\n")
cat("PB$raw names:", paste(names(PB$raw), collapse = ", "), "\n")
cat("PB$fine names:", paste(names(PB$fine), collapse = " | "), "\n")
META <- fread(file.path(TABWS, "ws_s4_pseudobulk_meta_GSE185381.csv"))
META[, status := fifelse(role == "control", "healthy", fifelse(role == "AML_Dx" & paediatric == FALSE, "AML", NA_character_))]
T25 <- META[group4 == "T" & eligible == TRUE & !is.na(status), .(sample, status)]
cat("T eligible:", T25[, paste(status, .N, collapse = " / "), by = status]$V1, "\n")
# are covariates sample-level or sample x compartment?
cv <- META[sample %in% T25$sample, .(n_hb = uniqueN(hb_share_nonery), n_gran = uniqueN(gran_share_lymph), n_dbl = uniqueN(dbl_mean)), by = sample]
cat("covariate values per sample across group4 rows: hb max", max(cv$n_hb), " gran max", max(cv$n_gran), "\n")
for (g4 in c("B", "NK", "Mono")) { e <- META[group4 == g4 & eligible == TRUE & !is.na(status)]
  cat(g4, "eligible: AML", e[status=="AML", .N], "healthy", e[status=="healthy", .N], "; of them in T25:", sum(e$sample %in% T25$sample), "\n") }
CE <- PB$cells[kept == TRUE & sample %in% T25$sample & group4 == "T", .N, by = .(sample, fine)]
CE <- merge(CE, T25, by = "sample")
tab <- CE[, .(healthy_ge20 = sum(N >= 20 & status == "healthy"), AML_ge20 = sum(N >= 20 & status == "AML"),
             healthy_cells = sum(N[status == "healthy"]), AML_cells = sum(N[status == "AML"])), by = fine][order(-healthy_cells)]
print(tab)
fwrite(tab, file.path(OUT, "01_fine_label_cells.csv"))
sink()
