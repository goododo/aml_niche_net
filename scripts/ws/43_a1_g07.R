#!/usr/bin/env Rscript
# INPUT  : PREREGISTRATION_A1_G07.md section 2 (tag prereg-a1-g07-v1), every file md5-recorded by --known (K0): ws_pb_GSE185381.rds
#          ($raw, $fine, $cells), ws_pb_GSE120221.rds + donor_ages_GSE120221_v1.0.csv (Oetjen); results/tables/ws/ ws_s5_scores.csv,
#          ws_s5_de_{T,B,NK,Mono}.csv, ws_s5_imprint_genes.csv, ws_s8_t_imprint_ambient_prone.csv, ws_a1_g02_{fold_genes,oof_scores}.csv,
#          decontx/GSE185381/ws_s1_ambient_blacklist.csv; gene sets and lists under /LARGE1/gr10634/gaozy/reference/
# OUTPUT : results/tables/ws/ws_a1_g07_{checks,sets,sets_removal,sets_sens,opposite,compartments,cross,tags,classes,overlap,summary}.csv;
#          work files under /LARGE1/.../ws_a1_g07/; GATES.md block A1-G07 (--combine). Smoke runs write only to /LARGE1/.../ws_a1_g07/smoke/.
# what it does: PREREGISTRATION_A1_G07.md 3.1-3.4, 3.6, 4, 6 (K0-K5, K7, K9) and 7, nothing else (3.5 S5b is scripts/ws/44).
#   --prep     pins the 22 Gavish sets; derives the B/NK/Mono delay lists with K9      --known   K0-K5, K7, K9 (re-derived)
#   --sets     3.1, 3.2 and section 4 for them      --cross  3.3 and its 4.1/4.2 views (bootstrap for the primary only)      --tags  3.4
#   --classes  3.6(a)      --overlap 3.6(b)      --combine  section 7 readings, summary table, GATES block (refuses to run before 44 --score)
#   --smoke    any mode on synthetic counts with the real IDs (unknown arguments stop the script, so a mistyped --smoke cannot run a real mode)
# This script deletes nothing, modifies no input, and never runs or sources another analysis script (functions are copied, with sources).

options(warn = 1)   # warnings print where they happen (smoke review)
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma); library(Matrix); library(here); library(parallel) })
.a <- commandArgs(TRUE); SMOKE <- "--smoke" %in% .a
MODES <- c("prep", "known", "sets", "cross", "tags", "classes", "overlap", "combine")
MODE <- MODES[paste0("--", MODES) %in% .a]
if (length(MODE) != 1L) stop("give exactly one of ", paste0("--", MODES, collapse = " "), " [--smoke]")
bad_arg <- setdiff(.a, c(paste0("--", MODES), "--smoke"))   # a mistyped --smoke (e.g. --smok) would otherwise run the mode on real data
if (length(bad_arg)) stop("unknown argument(s): ", paste(bad_arg, collapse = " "))
source(here::here("scripts", "config", "config_paths.R"))    # SEED and the roots
source(here::here("scripts", "ws", "config_ws.R"))           # constants only: WS_B0_IEG (K0), WS_S4_MAP_FILE
NC <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))); if (is.na(NC) || NC < 1L) NC <- 2L
TABWS <- file.path(TAB_DIR, "ws"); REF <- "/LARGE1/gr10634/gaozy/reference"
GDIR <- file.path(LARGE1_DIR, "ws_a1_g07"); SDIR <- file.path(GDIR, "smoke")
WDIR <- if (SMOKE) SDIR else GDIR; ODIR <- if (SMOKE) SDIR else TABWS   # smoke tables never go under results/ (the 23:30 cron)
dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)
SFX <- if (SMOKE) "_smoke" else ""
SELF <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
MD5 <- if (length(SELF) && !is.na(SELF) && file.exists(SELF)) unname(tools::md5sum(SELF)) else NA_character_
if (!SMOKE && is.na(MD5)) stop("run the script with Rscript (its md5 is needed for every output)")
if (!SMOKE && file.access(SELF, 2) == 0) stop("a real run must use the dated read-only copy: sbatch scripts/ws/submit_ws_a1_g07.sbatch 43 --", MODE)   # an edit cannot reach a running job
ts <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
## write guards (42's pattern): a real (non-smoke) file is never overwritten; every table carries the script md5
wtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g07_", x, SFX, ".csv")); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }
put  <- function(d, x) { d <- as.data.table(d); if (!ncol(d)) d <- data.table(script_md5 = character(0)) else d[, script_md5 := MD5]; fwrite(d, wtab(x)); invisible(d) }   # an empty table is written as a header
rtab <- function(x) { f <- file.path(ODIR, paste0("ws_a1_g07_", x, SFX, ".csv")); if (!file.exists(f)) stop("missing ", f, " (run its mode first)"); fread(f) }
wg   <- function(x) { f <- file.path(WDIR, sub("(\\.[a-z.]+)$", paste0(SFX, "\\1"), x)); if (file.exists(f) && !SMOKE) stop("refusing to overwrite ", f); f }   # x.csv.gz -> x_smoke.csv.gz
rg   <- function(x) file.path(WDIR, sub("(\\.[a-z.]+)$", paste0(SFX, "\\1"), x))
ts(sprintf("[%s%s] %d worker processes", MODE, if (SMOKE) " smoke" else "", NC))

## -- 0. constants (section numbers refer to the pre-registration) --------------------------------------------------------------
Q <- 0.05; MINSET <- 15L; MINREM <- 10L; HALF <- 0.5; RHO_SHARED <- 0.5                         # 3.1-3.3
N_K2A <- if (SMOKE) 60L else 500L; N_K2B <- if (SMOKE) 3L else 10L; N_K2C <- if (SMOKE) 3L else 20L; N_K5 <- if (SMOKE) 3L else 20L
N_RAND <- if (SMOKE) 50L else 1000L; N_BOOT <- if (SMOKE) 12L else 1000L                           # 3.2(1), 3.4, 3.3(b)
OFF <- c(k2a = 701L, k2b = 702L, k5 = 703L, tags_rand = 704L, fgsea = 705L, h_rand = 706L, boot = 707L, smoke = 709L)   # 8.3
T15 <- paste0("Control", 1:5)                                                                       # the T-ALL lanes (2.1, 4.3)
OPP12 <- c("ANXA1", "APOBEC3G", "CCL5", "CKLF", "GZMA", "GZMH", "GZMK", "HOPX", "KLRD1", "KLRG1", "NKG7", "TIGIT")   # 2.3
CLOCK <- c("ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "TEF", "HLF", "CIART", "BHLHE40", "BHLHE41")
GCLIST <- c("FKBP5", "TSC22D3", "ZBTB16", "KLF9", "DDIT4", "PER1", "SGK1", "DUSP1", "TXNIP", "CXCR4")
K4N <- c("CCR7", "LEF1", "TCF7", "SELL"); K4E <- c("GZMB", "GNLY", "NKG7", "PRF1", "GZMH")         # K4
CLASSES <- list(naive = c("CD4 Naive", "CD8 Naive"), central_memory = c("CD4 Central Memory", "CD8 Central Memory"),
                effector_memory = c("CD4 Effector Memory", "CD8 Effector Memory 1", "CD8 Effector Memory 2"))
CD4L <- c("CD4 Naive", "CD4 Central Memory", "CD4 Effector Memory", "CD4 Regulatory")
PAIRS_EXPECTED <- c("CD4_CELL_CYCLE|CD8_CELL_CYCLE", "CD4_CYTOTOXIC|CD8_CYTOTOXIC", "CD4_GLYCOLYSIS_MYC|CD8_GLYCOLYSIS_MYC", "CD4_INTERFERON|CD8_INTERFERON",
                    "CD4_NAIVE_1|CD8_MEMORY_NAIVE_1", "CD4_NAIVE_2|CD8_NAIVE_2", "CD4_STRESS_HSP|CD8_HEAT_SHOCK")   # 2.3
FULL_BELOW15 <- c("CD4_CELL_CYCLE", "CD4_CYTOTOXIC", "CD4_NAIVE_1", "CD4_STRESS_HSP", "CD8_CELL_CYCLE", "CD8_CYTOTOXIC", "CD8_MEMORY_NAIVE_1")   # 1.4
K9_SIZES <- list(T = c(792L, 295L, 497L), B = c(849L, 328L, 521L), NK = c(717L, 320L, 397L), Mono = c(826L, 300L, 526L))   # Mono: dated pre-run note (text read)
DELAY_TYPES <- list(T = c("CD4_Naive", "CD4_Memory", "CD8_Naive", "CD8_Memory"), B = c("B_Naive", "B_Memory"), NK = c("NK_dim", "NK_bright"),
                    Mono = c("CD14+_Monocytes", "CD16+_Monocytes"))
SAV_DIR <- file.path(REF, "savage2021_processing_delay"); SAV_T <- file.path(SAV_DIR, "savage2021_T_delay18h_consistent.csv"); SAV_X <- file.path(SAV_DIR, "mmc5.xlsx")
HM_F <- file.path(REF, "h.all.v2024.1.Hs.symbols.gmt"); GAV_SRC <- here::here("LCC_proj", "scripts", "panels", "msigdb.v2026.1.Hs.symbols.gmt")
GAV_PIN <- if (SMOKE) file.path(SDIR, "gavish2023_3ca_T.v2026.1_smoke.gmt") else file.path(REF, "gavish2023_3ca_T.v2026.1.gmt")
DELAY_F <- function(cmp) if (cmp == "T") SAV_T else if (SMOKE) file.path(SDIR, sprintf("savage2021_%s_delay18h_consistent_smoke.csv", cmp)) else
  file.path(SAV_DIR, sprintf("savage2021_%s_delay18h_consistent.csv", cmp))
LASRY_X <- file.path(REF, "lasry2023_natcancer", "43018_2022_480_MOESM2_ESM.xlsx"); RUTELLA_X <- file.path(REF, "rutella2022_jci", "JCI159579.sdt1-10.xlsx")
MAZZ_DIR <- file.path(REF, "mazziotta2024_blood"); GTF <- file.path(REF, "refdata-gex-GRCh38-2020-A", "genes", "genes.gtf")
CNMF_F <- file.path(TAB_DIR, "04_cnmf", "normal_metaprograms.tsv")
PROV_F <- if (SMOKE) file.path(SDIR, "PROVENANCE_G07_smoke.txt") else file.path(SAV_DIR, "PROVENANCE_G07_derived.txt")   # next to the T set (2.3); PROVENANCE.txt is left as it is
if (!SMOKE) {   # every fixed target of the mode, work files included, is checked before any work
  tg <- switch(MODE, prep = c(vapply(c("B", "NK", "Mono"), DELAY_F, ""), GAV_PIN, PROV_F, file.path(GDIR, "prep_checks.csv")),
               known = c(file.path(TABWS, "ws_a1_g07_checks.csv"), file.path(GDIR, c("input_md5.csv", "k2_plants.csv", "k5_plants.csv"))),
               sets = file.path(TABWS, paste0("ws_a1_g07_", c("sets", "sets_removal", "sets_sens", "opposite"), ".csv")),
               cross = file.path(TABWS, paste0("ws_a1_g07_", c("compartments", "cross"), ".csv")),
               tags = c(file.path(TABWS, "ws_a1_g07_tags.csv"), file.path(GDIR, "tags_per_gene.csv.gz")),
               classes = file.path(TABWS, "ws_a1_g07_classes.csv"), overlap = file.path(TABWS, "ws_a1_g07_overlap.csv"),
               combine = file.path(TABWS, "ws_a1_g07_summary.csv"))
  if (any(file.exists(tg))) stop("refusing to overwrite ", paste(tg[file.exists(tg)], collapse = ", "), " (a rerun needs a dated note and the old file renamed)")
  if (MODE == "combine" && any(grepl("A1-G07 start", readLines(here::here("scripts", "ws", "GATES.md"))))) stop("GATES.md already holds an A1-G07 block")
}

## known-answer bookkeeping (section 6: PASS, FAIL_deterministic, FAIL_statistical; "reported" rows carry no state)
CHK <- list()
chk <- function(id, scope, name, expected, got, tol = NULL) {
  ok <- if (is.null(tol)) isTRUE(all.equal(expected, got)) else length(got) == length(expected) && isTRUE(all(abs(expected - got) <= tol))
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "deterministic", check = name, expected = paste(format(expected), collapse = ";"),
                                         got = paste(format(got), collapse = ";"), state = if (ok) "PASS" else "FAIL_deterministic")
  message(sprintf("  %-4s %-90s expect %-16s got %-16s %s", id, name, paste(format(expected), collapse = ";"), paste(format(got), collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
chk_stat <- function(id, scope, name, bar, got, ok) {
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "statistical", check = name, expected = bar, got = paste(got, collapse = ";"),
                                         state = if (isTRUE(ok)) "PASS" else "FAIL_statistical")
  message(sprintf("  %-4s %-90s bar %-22s got %-16s %s", id, name, bar, paste(got, collapse = ";"), if (isTRUE(ok)) "PASS" else "** FAIL (statistical) **")) }
rep_ <- function(id, scope, name, value) {
  CHK[[length(CHK) + 1L]] <<- data.table(id = id, scope = scope, kind = "reported", check = name, expected = "", got = paste(value, collapse = ";"), state = "REPORTED")
  message(sprintf("  %-4s %-90s %s [REPORTED]", id, name, paste(value, collapse = ";"))) }
checks <- function() { f <- file.path(ODIR, paste0("ws_a1_g07_checks", SFX, ".csv")); if (!file.exists(f)) stop("run --known first: ", f); fread(f) }
gate <- function(sc) {   # a mode runs only when every deterministic check of its scope passed, under this script version
  C <- checks(); if (!SMOKE && !identical(unique(C$script_md5), MD5)) stop("the checks table was written by another version of this script; rerun --known")
  bad <- C[grepl(sc, scope) & kind == "deterministic" & state != "PASS"]
  if (nrow(bad)) stop("deterministic known answers failed for scope ", sc, ": ", paste(bad$id, collapse = ", ")); invisible(C) }
stat_ok <- function(C, ids) all(C[id %in% ids & kind == "statistical", state] == "PASS")
syn <- function(M) { attr(M, "synthetic") <- TRUE; M }
need_syn <- function(M, what) if (SMOKE && !isTRUE(attr(M, "synthetic"))) stop("smoke run: a real ", what, " matrix would enter a computation")
cols <- function(M, j) { X <- M[, j, drop = FALSE]; if (isTRUE(attr(M, "synthetic"))) syn(X) else X }
mcl <- function(X, FUN) { r <- mclapply(X, FUN, mc.cores = NC, mc.preschedule = FALSE)
  bad <- vapply(r, function(z) is.null(z) || inherits(z, "try-error"), TRUE); if (any(bad)) stop(sum(bad), " of ", length(r), " parallel workers failed: ", paste(unique(unlist(lapply(r[bad], as.character))), collapse = " | ")); r }
recheck_inputs <- function() { f <- rg("input_md5.csv"); if (!file.exists(f)) stop("run --known first: ", f); I <- fread(f)
  if (!SMOKE && !identical(unique(I$script_md5), MD5)) stop("input_md5.csv was written by another version of this script; rerun --known")
  now <- unname(tools::md5sum(I$path)); bad <- I$path[is.na(now) | now != I$md5]; if (length(bad)) stop("inputs changed since --known: ", paste(basename(bad), collapse = ", ")); invisible(I) }
read_gmt <- function(f, pat = NULL) { L <- strsplit(readLines(f), "\t"); S <- lapply(L, function(x) unique(x[-(1:2)])); names(S) <- vapply(L, `[`, "", 1)   # 20_t3c:85
  if (!is.null(pat)) S <- S[grepl(pat, names(S))]; S }
short <- function(x) sub("^GAVISH_3CA_METAPROGRAM_", "", sub("_T_CELLS_", "_", x))
f3 <- function(x) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3))
fq <- function(q) ifelse(is.na(q), "NA", ifelse(q < 0.001, "< 0.001", formatC(q, format = "fg", digits = 2)))   # FDR and p in sentences (never "0.000")

## -- 1. --prep: the pinned Gavish file (9.4) and the per-compartment delay lists (9.14; K9) ----------------------------------------
savage_rule <- function(types) {   # review2/R2b/02_savage_rule.R build(types, "any", "18h_healthy"), the reading of 2.3
  suppressPackageStartupMessages(library(readxl))
  rd <- function(s, e) { d <- as.data.table(suppressMessages(read_excel(SAV_X, sheet = s, skip = 1, col_types = "text")))   # all text: 163 tiny p values
    setnames(d, 1:7, c("ct", "smp", "tp", "gene", "lfc", "p", "padj")); d[, exp := e]                                       # are stored as text (dated note)
    d[, `:=`(tp = as.numeric(tp), lfc = as.numeric(lfc), padj = as.numeric(padj))]; if (e == 2) d[, lfc := -lfc]          # Experiment 2 sign flipped (PROVENANCE)
    d[, .(exp, ct, smp, tp, gene, lfc, padj)] }
  D0 <- rbind(rd("Experiment 1", 1), rd("Experiment 2", 2))[!is.na(padj) & !is.na(lfc)]
  s18 <- D0[ct %in% types & grepl("Healthy", smp) & tp == 18 & padj < 0.05]
  a <- s18[, .(dir = unique(sign(lfc))), by = .(exp, gene)]
  both <- a[, .(ne = uniqueN(exp), nd = uniqueN(dir), dir = dir[1]), by = gene][ne == 2 & nd == 1]
  oppg <- merge(s18[, .(gene, sg = sign(lfc))], both[, .(gene, dir)], by = "gene")[sg != dir, unique(gene)]
  both[!gene %in% oppg, .(gene, direction = ifelse(dir > 0, "up", "down"))][order(gene)] }
if (MODE == "prep") {   # every K9 row is evaluated first; nothing is written unless all pass; files are locked only after all are written
  ts("[prep] K9: input md5s; the 2.3 reading of Table S4 reproduces the frozen T list; B, NK and Mono sizes; the 22 Gavish T-cell sets")
  chk("K9", "P", "md5 of mmc5.xlsx, the Savage T list (prefix) and the Gavish source (as recorded)", c("60fac77c3d75689f8abdc63bbd2f2a26", "5ee444a9", "30893650e0c0735a39819f9f7691f0d6"),
      c(unname(tools::md5sum(SAV_X)), substr(unname(tools::md5sum(SAV_T)), 1, 8), unname(tools::md5sum(GAV_SRC))))
  FZ <- fread(SAV_T)
  LST <- lapply(DELAY_TYPES, savage_rule)
  chk("K9", "P", "T rule reproduces savage2021_T_delay18h_consistent.csv (gene + direction, identical)", TRUE,
      identical(sort(paste(LST$T$gene, LST$T$direction)), sort(paste(FZ$gene, FZ$direction))))
  for (k in names(K9_SIZES)) chk("K9", "P", sprintf("%s list: n, up, down", k), K9_SIZES[[k]], c(nrow(LST[[k]]), LST[[k]][direction == "up", .N], LST[[k]][direction == "down", .N]))
  gl <- readLines(GAV_SRC); gv <- gl[grepl("^GAVISH_3CA_METAPROGRAM_CD[48]_T_CELLS_", gl)]   # extracted verbatim from the repo's MSigDB v2026.1 file
  chk("K9", "P", "Gavish T-cell sets extracted", 22L, length(gv))
  C <- rbindlist(CHK)
  if (!all(C$state == "PASS")) { ff <- file.path(WDIR, sprintf("prep_FAILED_%s%s.csv", format(Sys.time(), "%Y%m%d%H%M%S"), SFX)); if (file.exists(ff)) stop("refusing to overwrite ", ff)
    fwrite(C[, script_md5 := MD5], ff); stop("K9 failed: nothing written") }
  OUTF <- c(vapply(c("B", "NK", "Mono"), DELAY_F, ""), GAV_PIN, PROV_F)
  for (k in c("B", "NK", "Mono")) fwrite(LST[[k]], DELAY_F(k))
  writeLines(gv, GAV_PIN)
  szs <- vapply(c("B", "NK", "Mono"), function(k) sprintf("%s %d (%d / %d)", k, nrow(LST[[k]]), LST[[k]][direction == "up", .N], LST[[k]][direction == "down", .N]), "")
  prov <- c(sprintf("Written %s by scripts/ws/43_a1_g07.R --prep (script md5 %s), A1-G07 pre-registration decisions 9.4 and 9.14 (user approval 2026-11-09). It extends this folder's PROVENANCE.txt, which is an input and is left unchanged.", format(Sys.time()), MD5),
            sprintf("%s: the 22 lines starting GAVISH_3CA_METAPROGRAM_CD4_T_CELLS_ / _CD8_T_CELLS_, copied verbatim from %s (md5 %s; committed in aac665f, MSigDB v2026.1, CC BY 4.0); md5 of this file %s.",
                    GAV_PIN, GAV_SRC, unname(tools::md5sum(GAV_SRC)), unname(tools::md5sum(GAV_PIN))),
            "Gavish M, et al. 2023. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. Nature 618:598-606 (3CA metaprograms).",
            sprintf("savage2021_{B,NK,Mono}_delay18h_consistent.csv: derived from mmc5.xlsx (md5 %s) with the T list's rule as read in the G07 pre-registration 2.3", unname(tools::md5sum(SAV_X))),
            "(any of the compartment's cell types significant at 18 h vs 2 h, adjusted p < 0.05, in Pilot01_Healthy and in Pilot02_Healthy with Experiment 2's sign flipped; one direction; no significant opposite-sign 18 h healthy call).",
            "Table S4 was read with readxl col_types = 'text' and converted with as.numeric: 163 tiny adjusted p values are stored as text, and a default read turns them into NA (G07 dated pre-run note, section 11).",
            "Cell types: B = B_Naive, B_Memory; NK = NK_dim, NK_bright; Mono = CD14+_Monocytes, CD16+_Monocytes. The same rule reproduces savage2021_T_delay18h_consistent.csv exactly (K9).",
            sprintf("Sizes, genes (up / down): %s. 2.3 and K9 of the pre-registration quote Mono 822 (299 / 523) from a default read; the text read gives the numbers here (dated pre-run note).", paste(szs, collapse = "; ")),
            sprintf("md5: %s", paste(sprintf("%s %s", c("B", "NK", "Mono"), unname(tools::md5sum(vapply(c("B", "NK", "Mono"), DELAY_F, "")))), collapse = "; ")))
  writeLines(prov, PROV_F)
  if (!SMOKE) for (f in OUTF) Sys.chmod(f, "0444")   # locked only once every file is written
  fwrite(rbindlist(CHK)[, script_md5 := MD5], wg("prep_checks.csv")); ts("[prep] done"); quit(status = 0)
}

## -- 2. data: designs, lists, universes (27:28-67; 37:45-80) -------------------------------------------------------------------------
BL <- fread(file.path(TABWS, "decontx", "GSE185381", "ws_s1_ambient_blacklist.csv"))
GO <- fread(file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), header = FALSE, select = 1:2, col.names = c("gene", "chr"))
PB <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE185381", "ws_pb_GSE185381.rds")); G <- PB$genes
BLACK <- unique(c("HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "MPO", "ELANE", "AZU1", "PRTN3", BL[n_units >= 5, gene],
                  "XIST", GO[chr == "chrY", gene], grep("^MT-", G, value = TRUE), grep("^RP[LS][0-9]", G, value = TRUE)))   # 27:47-49
S5ALL <- fread(file.path(TABWS, "ws_s5_scores.csv"))   # design columns only; the score column is never read (2.1)
dsg <- function(cmp) S5ALL[compartment == cmp, .(sample, status = factor(status, c("healthy", "AML")), hb, gran, sex = factor(sex), subtype_share)]
CMPS <- c(T = "T", B = "B", NK = "NK", Mono = "Mono"); D <- lapply(CMPS, dsg)
S25 <- D$T$sample; AMLS <- D$T[status == "AML", sample]; HEAL <- D$T[status == "healthy", sample]
IMP <- fread(file.path(TABWS, "ws_s5_imprint_genes.csv"))[compartment == "T" & imprint_called == TRUE]
AMB <- fread(file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"))$gene
F402 <- IMP[, .(gene, s = sign(logFC))]; P370 <- F402[!gene %in% AMB]   # the frozen lists (gene and direction only)
TESTED <- lapply(CMPS, function(k) fread(file.path(TABWS, sprintf("ws_s5_de_%s.csv", k)), select = "gene")$gene)   # gene names only
UT <- TESTED$T; EXT <- BL[n_units >= 2, gene]; UA <- setdiff(UT, EXT)                                                # 2.2
stopifnot(length(S25) == 25L, length(AMLS) == 15L, nrow(F402) == 402L, nrow(P370) == 370L, all(P370$gene %in% UA))
MAP <- fread(here::here("scripts", "ws", WS_S4_MAP_FILE), fill = TRUE); TL <- MAP[group4 == "T", bmm_fine]; stopifnot(length(TL) == 9L)
CD8L <- setdiff(TL, CD4L)
HM <- read_gmt(HM_F)
GV <- if (file.exists(GAV_PIN)) read_gmt(GAV_PIN) else NULL
if (is.null(GV) && MODE != "known") stop("run --prep first: ", GAV_PIN)
if (!is.null(GV)) { names(GV) <- short(names(GV)); UNASSIGNED <- grep("UNASSIGNED", names(GV), value = TRUE); GV19 <- GV[!names(GV) %in% UNASSIGNED] }
SV <- fread(SAV_T); IEG <- WS_B0_IEG
H_up <- union(IEG, SV[direction == "up", gene]); H_dn <- union(IEG, SV[direction == "down", gene]); H <- union(H_up, H_dn)
gt <- fread(cmd = sprintf("grep -v '^#' %s | awk -F'\\t' '$3==\"gene\"'", GTF), header = FALSE, sep = "\t", quote = "", select = 9)
gt[, `:=`(gname = sub('.*gene_name "([^"]+)".*', "\\1", V9), gtype = sub('.*gene_type "([^"]+)".*', "\\1", V9))]
LNC <- gt[, .(only_lnc = all(gtype == "lncRNA")), by = gname][only_lnc == TRUE, gname]; rm(gt)
if (!SMOKE) {
  RAW <- lapply(CMPS, function(k) PB$raw[[k]][, D[[k]]$sample])
  FINE <- lapply(setNames(TL, TL), function(k) { M <- PB$fine[[k]]; out <- matrix(0, length(G), length(S25), dimnames = list(G, S25))   # zero-filled (39:120-122)
    if (!is.null(M)) { h <- intersect(colnames(M), S25); out[, h] <- M[, h, drop = FALSE] }; out })
  NCELL <- matrix(0, length(TL), length(S25), dimnames = list(TL, S25))
  tmp <- PB$cells[kept == TRUE & sample %in% S25 & fine %in% TL, .N, by = .(fine, sample)]; NCELL[cbind(tmp$fine, tmp$sample)] <- tmp$N
}
rm(PB); invisible(gc(FALSE))

## -- 3. smoke data: synthetic counts, covariates, labels and an Oetjen stand-in on the real IDs ------------------------------------------
if (SMOKE) {
  set.seed(SEED + OFF[["smoke"]]); SMS <- sample.int(.Machine$integer.max, 30L)
  ng <- length(G); set.seed(SMS[1]); mu <- setNames(numeric(ng), G); on <- sample.int(ng, 12000L); mu[on] <- exp(rnorm(12000L, 0, 1.8))
  mu[UT] <- pmax(mu[UT], exp(rnorm(length(UT), 2, 1)))                                     # the universe is expressed
  PUP <- unique(c(unlist(HM[c("HALLMARK_MTORC1_SIGNALING", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")]), intersect(HM$HALLMARK_TNFA_SIGNALING_VIA_NFKB, H_up)))   # the last part
  PDN <- unique(unlist(HM["HALLMARK_INTERFERON_ALPHA_RESPONSE"]))                                     # plants a set carried by handling genes
  if (!is.null(GV)) PDN <- union(PDN, GV[["CD8_NAIVE_2"]])
  synth <- function(ids, depth, aml, fu, fd) { M <- vapply(ids, function(s) { m <- mu * depth / sum(mu)
      if (s %in% aml) { m[names(m) %in% PUP] <- m[names(m) %in% PUP] * fu; m[names(m) %in% PDN] <- m[names(m) %in% PDN] * fd }
      rnbinom(ng, size = 10, mu = m) }, numeric(ng)); dimnames(M) <- list(G, ids); syn(M) }
  set.seed(SMS[2]); RAW <- list(T = synth(D$T$sample, 2e6, AMLS, 2^0.9, 2^-0.9), B = synth(D$B$sample, 5e5, AMLS, 2^0.7, 2^-0.2),
                                NK = synth(D$NK$sample, 5e5, AMLS, 1, 1), Mono = synth(D$Mono$sample, 2e6, character(0), 1, 1))
  set.seed(SMS[3]); for (k in CMPS) D[[k]][, `:=`(hb = runif(.N, 0, 0.02), gran = runif(.N, 0, 0.01),
                                                    sex = factor(sample(c("female", "male"), .N, TRUE), levels = c("female", "male")), subtype_share = runif(.N, 0.2, 0.6))]
  for (k in CMPS) { i <- match(D[[k]]$sample, D$T$sample); j <- !is.na(i); D[[k]][j, `:=`(hb = D$T$hb[i[j]], gran = D$T$gran[i[j]], sex = D$T$sex[i[j]])] }   # donor-level covariates agree
  lprof <- function(boost) { p <- mu; p[boost] <- p[boost] * 8; p }
  LPROF <- lapply(setNames(TL, TL), function(k) lprof(if (k %in% CLASSES$naive) K4N else if (k %in% CLASSES$effector_memory) K4E else character(0)))
  set.seed(SMS[4]); NCELL <- matrix(rpois(length(TL) * 25L, 120), length(TL), 25L, dimnames = list(TL, S25))
  FINE <- lapply(setNames(TL, TL), function(k) syn(vapply(S25, function(s) { m <- NCELL[k, s] * 1500 * LPROF[[k]] / sum(LPROF[[k]])   # the AML plant goes into
    if (s %in% AMLS) { m[names(m) %in% PUP] <- m[names(m) %in% PUP] * 2^0.9; m[names(m) %in% PDN] <- m[names(m) %in% PDN] * 2^-0.9 }    # the labels, so whole T
    rnbinom(ng, size = 20, mu = m) }, numeric(ng))))                                                                                    # (their sum) carries it
  for (k in TL) dimnames(FINE[[k]]) <- list(G, S25)
  RAW$T[, S25] <- Reduce(`+`, FINE); RAW$T <- syn(RAW$T)                                    # whole T = sum of the labels in smoke
}

## -- 4. engines -------------------------------------------------------------------------------------------------------------------------
fml <- function(terms) as.formula(paste("~", paste(terms, collapse = " + ")))
MAIN <- c("status", "hb", "gran", "sex", "subtype_share")
fit_s5 <- function(cnt, d, terms = MAIN) {   # 27:70-87 edgeR branch, return_y = TRUE, plus the QL F (4.4 fgsea)
  need_syn(cnt, "count")
  if ("sex" %in% terms && uniqueN(d$sex) < 2) terms <- setdiff(terms, "sex")
  y <- DGEList(cnt[!rownames(cnt) %in% BLACK, d$sample, drop = FALSE], group = d$status)
  keep <- filterByExpr(y, group = d$status) & (rowMeans(cpm(y) > 1) >= 0.5)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE], method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(qr(X)$rank == ncol(X), colnames(X)[2] == "statusAML")
  if (nrow(X) - ncol(X) < 2) stop("residual df < 2: not estimable")
  y <- estimateDisp(y, X, robust = TRUE); f <- glmQLFit(y, X, robust = TRUE); t <- glmQLFTest(f, coef = "statusAML")
  tt <- as.data.table(topTags(t, n = Inf, sort.by = "none")$table, keep.rownames = "gene")[, .(gene, logFC, logCPM, F, PValue, FDR)]
  list(tt = tt, y = y, X = X, lfc = setNames(tt$logFC, tt$gene)) }
fit_rows <- function(cnt, d, terms = MAIN) {   # K2b/K2c/K5/3.3(b): fixed rows (no re-filtering), TMM recomputed (review 2)
  need_syn(cnt, "count"); if ("sex" %in% terms && uniqueN(d$sex) < 2) terms <- setdiff(terms, "sex")
  y <- calcNormFactors(DGEList(cnt[, d$sample, drop = FALSE], group = d$status), method = "TMM")
  X <- model.matrix(fml(terms), data = as.data.frame(d)); stopifnot(qr(X)$rank == ncol(X)); y <- estimateDisp(y, X, robust = TRUE); list(y = y, X = X) }
lfc_rows <- function(cnt, d, terms = MAIN) { r <- fit_rows(cnt, d, terms); f <- glmQLFit(r$y, r$X, robust = TRUE)   # logFC on fixed rows (rho_X)
  setNames(f$coefficients[, "statusAML"] / log(2), rownames(r$y)) }
cam <- function(y, X, sets, univ, contrast = "statusAML", cor = NA, min_n = MINSET) {   # camera.DGEList on the rows in univ (2.2)
  yy <- y[rownames(y) %in% univ, , keep.lib.sizes = TRUE]
  idx <- ids2indices(sets, rownames(yy), remove.empty = FALSE); n <- lengths(idx); ok <- n >= min_n
  out <- data.table(set = names(sets), NGenes = as.integer(n), tested = ok, Correlation = NA_real_, Direction = NA_character_, PValue = NA_real_)
  if (any(ok)) { r <- as.data.table(camera(yy, idx[ok], X, contrast = contrast, inter.gene.cor = cor, allow.neg.cor = FALSE, sort = FALSE), keep.rownames = "set")
    j <- match(r$set, out$set); out[j, `:=`(Direction = r$Direction, PValue = r$PValue)]; if ("Correlation" %in% names(r)) out[j, Correlation := r$Correlation] }
  out[, FDR := NA_real_][tested == TRUE, FDR := p.adjust(PValue, "BH")][] }
shift_of <- function(lfc, set, univ) { u <- intersect(names(lfc), univ); a <- intersect(set, u); if (!length(a) || length(a) == length(u)) return(NA_real_)
  mean(lfc[a]) - mean(lfc[setdiff(u, a)]) }
sgn <- function(dirn) ifelse(dirn == "Up", 1, -1)
plant <- function(M, genes, samples, lfc) {   # count-level plant (K2b, K2c, K5): Poisson add for up, binomial thinning for down
  sub <- M[genes, samples, drop = FALSE]
  sub <- if (lfc > 0) sub + matrix(rpois(length(sub), as.vector(sub) * (2^lfc - 1)), nrow(sub)) else matrix(rbinom(length(sub), as.vector(sub), 2^lfc), nrow(sub))
  M[genes, samples] <- sub; M }
bins_of <- function(e, nb = 10L) setNames(cut(rank(e, ties.method = "first"), nb, labels = FALSE), names(e))
match_draw <- function(target, pool, bins, n, seed) {   # random sets copying the target's bin composition (37:234-238, without replacement)
  target <- intersect(target, names(bins)); pool <- setdiff(intersect(pool, names(bins)), target); P <- split(pool, bins[pool]); tb <- table(bins[target]); set.seed(seed)
  lapply(seq_len(n), function(i) unlist(lapply(names(tb), function(b) { p <- P[[b]]; k <- tb[[b]]; if (is.null(p)) character(0) else if (length(p) <= k) p else p[sample.int(length(p), k)] }), use.names = FALSE)) }
pct_rank <- function(obs, ref) (1 + sum(ref < obs) + 0.5 * sum(ref == obs)) / (length(ref) + 1)
stream <- function(key, n) { set.seed(SEED + OFF[[key]]); sample.int(.Machine$integer.max, n) }   # one offset per stream (8.3); child seeds for each draw or job
logz <- function(M) {   # 37:65-68 (label-free log-CPM over the 25)
  need_syn(M, "count"); M <- M[!rownames(M) %in% BLACK, , drop = FALSE]; y <- DGEList(M); keep <- rowMeans(cpm(y) > 1) >= 0.5
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); cpm(y, log = TRUE, prior.count = 1) }
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
gav_pairs <- function(sets) { nm <- names(sets); out <- character(0)
  for (i in seq_along(nm)) for (j in seq_along(nm)) if (i < j && jac(intersect(sets[[i]], UA), intersect(sets[[j]], UA)) >= 0.3) out <- c(out, paste(sort(c(nm[i], nm[j])), collapse = "|"))
  sort(out) }

## stacked T x X model (3.3(b)); donor fixed effects absorb every donor-level main effect
paired <- function(X) intersect(D$T$sample, D[[X]]$sample)
stack_design <- function(X, dons, dT = NULL, dX = NULL) {   # dT / dX given only for the 4.2 recode
  if (is.null(dT)) dT <- D$T[match(dons, sample)]; if (is.null(dX)) dX <- D[[X]][match(dons, sample)]; n <- length(dons)
  donor <- factor(rep(dons, 2), levels = dons); cT <- rep(c(1, 0), each = n); aml <- as.numeric(rep(dT$status == "AML", 2))
  M <- cbind(model.matrix(~ 0 + donor), cT = cT, `cT:AML` = cT * aml, `cT:hb` = cT * rep(dT$hb, 2), `cT:gran` = cT * rep(dT$gran, 2),
             `cT:sexmale` = cT * as.numeric(rep(dT$sex == "male", 2)), `cT:CD8share` = cT * rep(dT$subtype_share, 2), `cT:Xshare` = cT * rep(dX$subtype_share, 2))
  if (uniqueN(dT$sex) < 2) M <- M[, colnames(M) != "cT:sexmale"]
  stopifnot(qr(M)$rank == ncol(M), nrow(M) - ncol(M) >= 2); M }
stack_counts <- function(X, MT, MX, dons, genes) { M <- cbind(MT[genes, dons, drop = FALSE], MX[genes, dons, drop = FALSE])
  colnames(M) <- c(paste0("T__", dons), paste0(X, "__", dons)); if (isTRUE(attr(MT, "synthetic")) && isTRUE(attr(MX, "synthetic"))) syn(M) else M }   # marked only if both halves are
stack_fit <- function(M, Xm) { need_syn(M, "stacked count"); y <- calcNormFactors(DGEList(M), method = "TMM"); estimateDisp(y, Xm, robust = TRUE) }
rsub <- function(M, i) { Y <- M[i, , drop = FALSE]; if (isTRUE(attr(M, "synthetic"))) syn(Y) else Y }   # row subset that keeps the smoke mark
rho_pair <- function(MT, MX, dT, dX, genes) list(T = lfc_rows(rsub(MT, genes), dT), X = lfc_rows(rsub(MX, genes), dX))

## -- 5. gene identity in healthy marrow (3.4): two references -----------------------------------------------------------------------------
tag_ref <- function(F, W, NCL, donors, genes, thr = 1, cmin = 0.8, min_don = 5L) {   # F: label -> genes x donors; W: whole-T; NCL: labels x donors
  need_syn(W, "whole-T count"); for (k in names(F)) need_syn(F[[k]], "label count")
  rbindlist(lapply(names(F), function(k) {
    el <- donors[NCL[k, donors] >= 20]; if (length(el) < min_don) return(NULL)
    Fk <- F[[k]][, el, drop = FALSE]; Rk <- W[, el, drop = FALSE] - Fk
    ck <- t(t(Fk) / colSums(Fk)) * 1e6; cr <- t(t(Rk) / colSums(Rk)) * 1e6; g <- intersect(genes, rownames(Fk))
    L <- log2((ck[g, , drop = FALSE] + 1) / (cr[g, , drop = FALSE] + 1))
    data.table(gene = g, label = k, M = matrixStats_rowMedians(L), c = rowMeans(L > 0), medCPM = matrixStats_rowMedians(ck[g, , drop = FALSE]), n_donors = length(el)) }))[
    , high := M >= thr & c >= cmin & medCPM >= 10][] }
matrixStats_rowMedians <- function(M) apply(M, 1, median)
class_of <- function(TG) { if (is.null(TG) || !nrow(TG)) return(data.table(gene = character(0)))
  out <- TG[, .(naive = any(high & label %in% CLASSES$naive), central_memory = any(high & label %in% CLASSES$central_memory),
                effector_memory = any(high & label %in% CLASSES$effector_memory)), by = gene]; out }
load_oetjen <- function() {   # 33:113-126 donor pooling; Oetjen's own gene axis (map by symbol)
  if (SMOKE) {   # synthetic stand-in: 20 donors, the same label profiles
    set.seed(SMS[5]); dn <- sprintf("O%02d", 1:20); nc <- matrix(rpois(length(TL) * 20L, 150), length(TL), 20L, dimnames = list(TL, dn))
    Fo <- lapply(setNames(TL, TL), function(k) { M <- vapply(dn, function(s) rnbinom(length(G), size = 20, mu = nc[k, s] * 1500 * LPROF[[k]] / sum(LPROF[[k]])), numeric(length(G))); dimnames(M) <- list(G, dn); syn(M) })
    return(list(F = Fo, W = syn(Reduce(`+`, Fo)), NCL = nc, donors = dn, axis = G)) }
  PO <- readRDS(file.path(LARGE1_DIR, "ws_pseudobulk", "GSE120221", "ws_pb_GSE120221.rds"))
  don <- fread(file.path(RAW_DIR, "metadata", "GSE120221", "donor_ages_GSE120221_v1.0.csv"), colClasses = "character")
  libs <- colnames(PO$raw$T); dmap <- don$donor[match(sub("_.*$", "", libs), don$gsm)]; stopifnot(!anyNA(dmap)); dn <- sort(unique(dmap))
  pool <- function(M) { out <- matrix(0, length(PO$genes), length(dn), dimnames = list(PO$genes, dn)); if (is.null(M)) return(out)
    for (j in seq_along(dn)) { l <- intersect(colnames(M), libs[dmap == dn[j]]); if (length(l)) out[, j] <- rowSums(M[, l, drop = FALSE]) }; out }
  Fo <- lapply(setNames(TL, TL), function(k) pool(PO$fine[[k]])); W <- pool(PO$raw$T)
  cc <- PO$cells[kept == TRUE & fine %in% TL]; cc[, donor := dmap[match(sample, libs)]]
  nc <- matrix(0, length(TL), length(dn), dimnames = list(TL, dn)); tmp <- cc[, .N, by = .(fine, donor)]; nc[cbind(tmp$fine, tmp$donor)] <- tmp$N
  list(F = Fo, W = W, NCL = nc, donors = dn, axis = PO$genes) }
tags_both <- function(OJ, thr = 1, cmin = 0.8, min_don = 5L, gdon = HEAL) {
  FG <- lapply(FINE, cols, gdon); WG <- cols(RAW$T, gdon)   # cols() keeps the smoke mark
  TGg <- tag_ref(FG, WG, NCELL, gdon, UT, thr, cmin, min_don); TGo <- tag_ref(OJ$F, OJ$W, OJ$NCL, OJ$donors, intersect(UT, OJ$axis), thr, 0.8, 5L)
  list(G = TGg, O = TGo, cG = class_of(TGg), cO = class_of(TGo)) }

## -- 6. the cross-fitted score with genes removed (3.6(a); A1-BATCH post-hoc a1_posthoc_d8_shared_genes.R:19-33) -------------------
oof_gap <- function(excl = character(0)) {
  FG <- fread(file.path(TABWS, "ws_a1_g02_fold_genes.csv")); L25 <- logz(cols(RAW$T, S25)); a <- S25 %in% AMLS
  o <- vapply(S25, function(s) { g <- FG[left_out == s & !gene %in% excl]; tr <- setdiff(S25, s)
    mu <- rowMeans(L25[, tr]); sg <- apply(L25[, tr], 1, sd); z <- (L25[, s] - mu) / sg; z <- z[is.finite(z)]
    u <- intersect(g[direction == "up", gene], names(z)); d <- intersect(g[direction == "down", gene], names(z))
    (if (length(u)) mean(z[u]) else 0) - (if (length(d)) mean(z[d]) else 0) }, 1)
  list(oof = o, gap = mean(o[a]) - mean(o[!a])) }

## -- 7. --known (section 6) --------------------------------------------------------------------------------------------------------------
if (MODE == "known") {
  ts("[K0] input md5s")
  INPUTS <- c(file.path(LARGE1_DIR, "ws_pseudobulk", c("GSE185381/ws_pb_GSE185381.rds", "GSE120221/ws_pb_GSE120221.rds")),
              file.path(RAW_DIR, "metadata", "GSE120221", "donor_ages_GSE120221_v1.0.csv"),
              file.path(TABWS, c("ws_s5_scores.csv", paste0("ws_s5_de_", CMPS, ".csv"), "ws_s5_imprint_genes.csv", "ws_s8_t_imprint_ambient_prone.csv",
                                 "ws_a1_g02_fold_genes.csv", "ws_a1_g02_oof_scores.csv", "decontx/GSE185381/ws_s1_ambient_blacklist.csv")),
              file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"), HM_F, GAV_SRC, SAV_T, SAV_X, LASRY_X, RUTELLA_X, GTF, CNMF_F,
              here::here("scripts", "config", "config_paths.R"), here::here("scripts", "ws", "config_ws.R"), here::here("scripts", "ws", WS_S4_MAP_FILE),
              file.path(TABWS, "ws_a1_g02_gene_selection.csv"))   # appended: G02's 25-gene core, read by --tags and --overlap
  MZ <- list.files(MAZZ_DIR, "\\.csv$", full.names = TRUE)   # Mazziotta enters 3.6(b) only as a csv on disk at --known (md5-recorded; dated note)
  m0 <- unname(tools::md5sum(INPUTS)); chk("K0", "S;X;G;C;O", sprintf("all %d declared inputs exist", length(INPUTS)), length(INPUTS), sum(!is.na(m0)))
  RECF <- c(ws_pb = INPUTS[1], scores = file.path(TABWS, "ws_s5_scores.csv"), de_T = file.path(TABWS, "ws_s5_de_T.csv"), imprint = file.path(TABWS, "ws_s5_imprint_genes.csv"),
            amb = file.path(TABWS, "ws_s8_t_imprint_ambient_prone.csv"), bl = file.path(TABWS, "decontx/GSE185381/ws_s1_ambient_blacklist.csv"), hm = HM_F, gav = GAV_SRC, sav = SAV_T,
            mmc5 = SAV_X, oetjen = INPUTS[2], oetjen_donors = INPUTS[3], gene_order = file.path(LARGE1_DIR, "reference", "gencode_GRCh38_gene_order.txt"),
            s4_map = here::here("scripts", "ws", WS_S4_MAP_FILE))
  REC <- c(ws_pb = "6c2fb685c49d4a23ef20a6842aa7aa67", scores = "58374fc60aabc53ba020c3149764fe8f", de_T = "b13de32b8eaeb400f00a0c32ca736c28",
           imprint = "977f58b06da076cc39a561ab7445b69d", amb = "acb2b043049a4286e4052d04cea15feb", bl = "acebed4dbeed2ef6ecd5a60791e28fb9",
           hm = "d3606c32", gav = "30893650e0c0735a39819f9f7691f0d6", sav = "5ee444a9", mmc5 = "60fac77c3d75689f8abdc63bbd2f2a26",   # mmc5: PREREGISTRATION_A1_batch.md:241
           oetjen = "b4ad90bd5c587d48cc66e02a537c5079", oetjen_donors = "178ced1c73e531b67172ed413b009c07",                                 # PREREGISTRATION_A1_G04.md:220-224
           gene_order = "c0950d592084a8b713de34ab30a7ddf8", s4_map = "f5a5a0283b243e5efaf17b1ddfdda5f1")
  stopifnot(identical(names(RECF), names(REC)), all(RECF %in% INPUTS))
  got_k <- substr(unname(tools::md5sum(RECF)), 1, nchar(REC))   # hm and sav are recorded as 8-character prefixes
  chk("K0", "S;X;G;C;O", sprintf("md5 of the %d inputs with recorded values (%s)", length(REC), paste(names(REC), collapse = ", ")), unname(REC), got_k)
  chk("K0", "O", "sha256 of the downloaded Lasry and Rutella files (PROVENANCE)", c("c2239a6107f43eb67ae6baaf67245312c7a7a94cfb237fdeb6bbc3c7ad91e5ef", "000eaaf485e777a525c1f207ba8b7d6d9591869585835210360bd07a83f0a0aa"),
      unname(vapply(c(LASRY_X, RUTELLA_X), function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1], "")))
  IEG26 <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "IER5", "ATF3", "DUSP1", "ZFP36", "NR4A1", "KLF2", "KLF4", "KLF6",
             "HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "DNAJB1", "PPP1R15A", "SOCS3", "CEBPD")   # frozen names: A1-BATCH post-hoc a1_posthoc_d8_shared_genes.R:21
  chk("K0", "S;X;G;C;O", "the IEG list WS_B0_IEG equals the 26 frozen names", sort(IEG26), sort(IEG))
  rep_("K0", "O", "Mazziotta 2024 csv files at --known, with md5 (only these can enter 3.6(b))", if (length(MZ)) paste(basename(MZ), unname(tools::md5sum(MZ))) else "none")
  if (!is.null(GV)) chk("K0", "S;X;G;O", "pinned Gavish file: 22 sets, each equal to the source", TRUE,
                        { src <- read_gmt(GAV_SRC, "^GAVISH_3CA_METAPROGRAM_CD[48]_T_CELLS_"); names(src) <- short(names(src)); length(GV) == 22L && identical(GV[sort(names(GV))], src[sort(names(src))]) })
  else chk("K0", "S;X;G;O", "pinned Gavish file exists (run --prep first)", TRUE, FALSE)
  MD5IN <- data.table(path = c(INPUTS, MZ, if (!is.null(GV)) c(GAV_PIN, DELAY_F("B"), DELAY_F("NK"), DELAY_F("Mono"))), md5 = NA_character_)
  MD5IN[, `:=`(md5 = unname(tools::md5sum(path)), script_md5 = MD5)]

  ts("[K9] the T list's identity and the list sizes; the derived delay lists equal a fresh derivation")
  LT <- savage_rule(DELAY_TYPES$T); FZ <- fread(SAV_T)
  chk("K9", "S;X;C;O", "T rule reproduces savage2021_T_delay18h_consistent.csv (gene + direction, identical)", TRUE,
      identical(sort(paste(LT$gene, LT$direction)), sort(paste(FZ$gene, FZ$direction))))
  chk("K9", "S;X;C;O", "T list: n, up, down", K9_SIZES$T, c(nrow(FZ), FZ[direction == "up", .N], FZ[direction == "down", .N]))
  for (k in c("B", "NK", "Mono")) { f <- DELAY_F(k); L1 <- if (file.exists(f)) fread(f) else data.table(gene = character(0), direction = character(0)); L2 <- savage_rule(DELAY_TYPES[[k]])
    chk("K9", "X", sprintf("%s delay list on disk equals the rule (n %d)", k, nrow(L2)), TRUE, identical(paste(L1$gene, L1$direction), paste(L2$gene, L2$direction)))
    chk("K9", "X", sprintf("%s list on disk: n, up, down", k), K9_SIZES[[k]], c(nrow(L1), L1[direction == "up", .N], L1[direction == "down", .N])) }

  ts("[K3] universes and set sizes (gene names only)")
  chk("K3", "S", "U_A, U_A - H_up, U_A - H_down, U_A - H", c(6957L, 6683L, 6517L, 6259L), c(length(UA), length(setdiff(UA, H_up)), length(setdiff(UA, H_dn)), length(setdiff(UA, H))))
  HM45 <- names(HM)[vapply(HM, function(s) length(intersect(s, UA)) >= MINSET, TRUE)]
  chk("K3", "S", "HALLMARK sets with >= 15 genes in U_A", 45L, length(HM45))
  if (!is.null(GV)) {
    mn <- function(sets, U) min(vapply(sets, function(s) length(intersect(s, U)), 1L))
    chk("K3", "S", "min genes left: HALLMARK(45) and Gavish(22) in U_A - H_up and U_A - H_down (>= 10)", TRUE,
        min(mn(HM[HM45], setdiff(UA, H_up)), mn(HM[HM45], setdiff(UA, H_dn)), mn(GV, setdiff(UA, H_up)), mn(GV, setdiff(UA, H_dn))) >= MINREM)
    chk("K3", "S", "Gavish sets below 15 genes after full removal (1.4)", sort(FULL_BELOW15), sort(names(GV)[vapply(GV, function(s) length(intersect(s, setdiff(UA, H))) < 15L, TRUE)]))
    chk("K3", "S", "near-duplicate pairs among the 19 interpreted sets (2.3)", sort(PAIRS_EXPECTED), gav_pairs(GV19))
    chk("K3", "S", "lncRNA genes in U_T (GTF gene_type)", 305L, length(intersect(UT, LNC)))
  }

  ts("[K1-K1b] S5 refits")
  if (!SMOKE) {
    st <- fread(file.path(TABWS, "ws_s5_de_T.csv")); r <- fit_s5(RAW$T, D$T)
    h <- r$tt[FDR < 0.10 & abs(logFC) > 0.5]
    chk("K1", "S;X;C", "T refit: genes tested; hits; up; down", c(7379L, 402L, 171L, 231L), c(nrow(r$tt), nrow(h), h[logFC > 0, .N], h[logFC < 0, .N]))
    chk("K1", "S;X;C", "T refit: hit set equals the frozen F402", TRUE, setequal(h$gene, F402$gene))
    chk("K1", "S;X;C", "T refit: max |logFC - ws_s5_de_T|", 0, max(abs(r$tt$logFC - st$logFC[match(r$tt$gene, st$gene)])), tol = 1e-9)
    chk("K1", "S;X;C", "design column 2 is statusAML", "statusAML", colnames(r$X)[2])
    for (k in c("B", "NK", "Mono")) { sk <- fread(file.path(TABWS, sprintf("ws_s5_de_%s.csv", k))); rk <- fit_s5(RAW[[k]], D[[k]])
      chk("K1b", "X", sprintf("%s refit: genes tested", k), c(B = 4292L, NK = 4151L, Mono = 9252L)[[k]], nrow(rk$tt))
      chk("K1b", "X", sprintf("%s refit: max |logFC - stored|", k), 0, max(abs(rk$tt$logFC - sk$logFC[match(rk$tt$gene, sk$gene)])), tol = 1e-9) }
  } else r <- fit_s5(RAW$T, D$T)

  ts("[K2a-K2c] camera on the T fit: random sets, plants, real-fit power")
  e <- setNames(r$y$AveLogCPM, rownames(r$y)); BIN <- bins_of(e[intersect(names(e), UA)]); HMs <- HM[HM45]
  ks <- stream("k2a", N_K2A + 1L); set.seed(ks[1]); tmpl <- sample(names(HMs), N_K2A, replace = TRUE)
  RS <- lapply(seq_len(N_K2A), function(i) match_draw(HMs[[tmpl[i]]], names(BIN), BIN, 1L, ks[i + 1L])[[1]]); names(RS) <- sprintf("rand%03d", seq_len(N_K2A))
  ca <- cam(r$y, r$X, RS, UA)
  chk_stat("K2a", "S", "camera (estimated correlation) on random matched sets: share with p < 0.05", "<= 0.10", f3(mean(ca$PValue < 0.05, na.rm = TRUE)), mean(ca$PValue < 0.05, na.rm = TRUE) <= 0.10)
  rep_("K2a", "S", "random sets: Up / Down split; median p", c(ca[Direction == "Up", .N], ca[Direction == "Down", .N], f3(median(ca$PValue, na.rm = TRUE))))
  rows <- rownames(r$y); base <- rsub(cols(RAW$T, D$T$sample), rows)   # rsub/cols keep the smoke mark; plant() keeps attributes
  plant_run <- function(i, lfc, frac, seed) { set.seed(seed); pg <- sample(intersect(UA, rows), 50L); pp <- if (frac < 1) pg[seq_len(round(frac * 50))] else pg
    M <- plant(base, pp, AMLS, lfc)
    ff <- fit_rows(M, D$T); yk <- DGEList(M, group = D$T$status); fails <- sum(!(filterByExpr(yk, group = D$T$status) & rowMeans(cpm(yk) > 1) >= 0.5)[pp])
    cc <- cam(ff$y, ff$X, c(HMs, list(PLANT = pg)), UA); p <- cc[set == "PLANT"]   # the stored p is the plant's own; its FDR (BH with the real sets) only sets the flags
    data.table(i = i, lfc = lfc, frac = frac, Direction = p$Direction, PValue = p$PValue, flagged_right = p$FDR < Q & p$Direction == if (lfc > 0) "Up" else "Down",
               flagged_wrong = p$FDR < Q & p$Direction != if (lfc > 0) "Up" else "Down", filter_fails = fails) }
  ps <- stream("k2b", 2L * N_K2B + 2L * N_K2C)
  jobs <- c(lapply(seq_len(N_K2B), function(i) list(i, 0.6, 1, ps[i])), lapply(seq_len(N_K2B), function(i) list(i, -0.6, 1, ps[N_K2B + i])),
            lapply(seq_len(N_K2C), function(i) list(i, 0.3, 1, ps[2L * N_K2B + i])), lapply(seq_len(N_K2C), function(i) list(i, 0.5, 0.3, ps[2L * N_K2B + N_K2C + i])))
  PL <- rbindlist(mcl(jobs, function(j) plant_run(j[[1]], j[[2]], j[[3]], j[[4]])))
  okb <- PL[frac == 1 & abs(lfc) == 0.6, .(right = sum(flagged_right), wrong = sum(flagged_wrong)), by = .(up = lfc > 0)]
  chk_stat("K2b", "S", "code check: 0.6 log2 plants flagged in their direction (up; down), never the wrong way", sprintf(">= %d/%d each; 0 wrong", ceiling(0.8 * N_K2B), N_K2B),
           c(okb[up == TRUE, right], okb[up == FALSE, right], sum(okb$wrong)), all(okb$right >= ceiling(0.8 * N_K2B)) && sum(okb$wrong) == 0)
  rep_("K2b", "S", "planted genes that would fail the S5 filter (total over plants)", PL[abs(lfc) == 0.6, sum(filter_fails)])
  rep_("K2c", "S", "real-fit power: whole-set 0.3 log2 plants flagged; 30 % of genes at 0.5 flagged", c(PL[lfc == 0.3, sum(flagged_right)], PL[frac < 1, sum(flagged_right)], N_K2C))
  fwrite(PL[, script_md5 := MD5], wg("k2_plants.csv"))

  ts("[K4] healthy-subtype tags: positive control in both references")
  OJ <- load_oetjen(); TB <- tags_both(OJ)
  for (rf in c("G", "O")) { cl <- TB[[paste0("c", rf)]]; nv <- intersect(K4N, cl$gene); ef <- intersect(K4E, cl$gene)
    a1 <- cl[gene %in% nv, sum(naive)]; a2 <- cl[gene %in% ef, sum(effector_memory)]; w1 <- cl[gene %in% nv, sum(effector_memory)]; w2 <- cl[gene %in% ef, sum(naive)]
    chk_stat("K4", "G", sprintf("%s: naive genes tagged naive (of %d); effector genes tagged EM (of %d); opposite tags", rf, length(nv), length(ef)), ">= 3; >= 3; 0",
             c(a1, a2, w1 + w2), length(nv) >= 3 && length(ef) >= 3 && a1 >= 3 && a2 >= 3 && w1 + w2 == 0) }

  ts("[K5] stacked T x X plants (code checks; one set per refit; rho_X point estimates)")
  K5 <- list(); k5s <- stream("k5", 2L * 5L * N_K5)
  for (X in c("B", "NK")) {
    dons <- paired(X); com <- setdiff(intersect(UT, TESTED[[X]]), EXT); if (SMOKE) com <- intersect(com, G)
    Xm <- stack_design(X, dons); dT <- D$T[match(dons, sample)]; dX <- D[[X]][match(dons, sample)]
    MT <- cols(RAW$T, dons); MX <- cols(RAW[[X]], dons); aml <- dT[status == "AML", sample]
    HMx <- HM[vapply(HM, function(s) length(intersect(s, com)) >= MINSET, TRUE)]; stopifnot(length(HMx) == c(B = 41L, NK = 43L)[[X]])   # K5 BH family (section 6)
    types <- list(i = c(0.6, 0), ii = c(0.6, 0.6), iii = c(0.6, 0.3), iv = c(0.3, 0), v = c(0.3, 0.3))
    jobs <- unlist(lapply(names(types), function(ty) lapply(seq_len(N_K5), function(i) list(ty = ty, i = i))), recursive = FALSE)
    res <- rbindlist(mcl(jobs, function(j) { set.seed(k5s[(match(X, c("B", "NK")) - 1L) * 5L * N_K5 + (match(j$ty, names(types)) - 1L) * N_K5 + j$i])
      pg <- sample(com, 50L); tl <- types[[j$ty]]; mt <- plant(MT, pg, aml, tl[1]); mx <- if (tl[2] > 0) plant(MX, pg, aml, tl[2]) else MX   # plant() keeps the smoke mark
      y <- stack_fit(stack_counts(X, mt, mx, dons, com), Xm); cc <- cam(y, Xm, c(HMx, list(PLANT = pg)), com, contrast = "cT:AML"); p <- cc[set == "PLANT"]
      rp <- rho_pair(mt, mx, dT, dX, com); sT <- shift_of(rp$T, pg, com); sX <- shift_of(rp$X, pg, com); rho <- if (!is.na(sT) && sT > 0) sX / sT else NA_real_
      rd <- if (p$FDR < Q && p$Direction == "Up") "larger in T" else if (p$FDR < Q) "larger in X" else if (!is.na(rho) && rho >= RHO_SHARED) "shared" else "not resolved"
      data.table(X = X, type = j$ty, i = j$i, int_dir = p$Direction, int_p = p$PValue, rho = rho, reading = rd) }))   # the plant's own p (its FDR depends on the real sets)
    K5[[X]] <- res
    n1 <- res[type == "i", sum(reading == "larger in T")]; n2 <- res[type == "ii", sum(reading == "larger in T")]; s2 <- res[type == "ii", sum(reading == "shared")]
    chk_stat("K5", "X", sprintf("%s: (i) T-only 0.6 'larger in T'; (ii) shared 0.6 'larger in T' and 'shared'", X),
             sprintf(">= %d; <= %d and >= %d", ceiling(0.7 * N_K5), floor(0.1 * N_K5), ceiling(0.7 * N_K5)), c(n1, n2, s2),
             n1 >= ceiling(0.7 * N_K5) && n2 <= floor(0.1 * N_K5) && s2 >= ceiling(0.7 * N_K5))
    rep_("K5", "X", sprintf("%s: readings of (iii) half-shared, (iv) T-only 0.3, (v) shared 0.3", X),
         vapply(c("iii", "iv", "v"), function(ty) paste(names(table(res[type == ty, reading])), table(res[type == ty, reading]), sep = ":", collapse = ","), "")) }
  fwrite(rbindlist(K5)[, script_md5 := MD5], wg("k5_plants.csv"))

  ts("[K7] gene-class removal machinery (3.6(a))")
  if (!SMOKE) { OOFT <- fread(file.path(TABWS, "ws_a1_g02_oof_scores.csv")); g0 <- oof_gap()
    chk("K7", "C", "rebuilt out-of-fold scores equal ws_a1_g02_oof_scores.csv (max |diff|)", 0, max(abs(g0$oof - OOFT[match(S25, sample), oof_score])), tol = 1e-9)
    chk("K7", "C", "gap: none; without Savage; without Savage + IEG", c(1.3782, 1.3184, 1.3228), round(c(g0$gap, oof_gap(SV$gene)$gap, oof_gap(c(SV$gene, IEG))$gap), 4), tol = 1e-9)
  } else rep_("K7", "C", "smoke: gap without removal (synthetic)", f3(oof_gap()$gap))

  C <- rbindlist(CHK)[, script_md5 := MD5]
  fwrite(C, wtab("checks")); fwrite(MD5IN, wg("input_md5.csv"))
  ts(sprintf("[known] %d rows; deterministic FAIL %d; statistical FAIL %d", nrow(C), C[state == "FAIL_deterministic", .N], C[state == "FAIL_statistical", .N])); quit(status = 0)
}

## -- 8. --sets: 3.1, 3.2, section 4 for them ------------------------------------------------------------------------------------------
removal <- function(fam, sets, P, y, X, lfc, univ) {   # 3.2(2)-(5) for the flagged sets of one family
  fl <- P[tested == TRUE & FDR < Q]; if (!nrow(fl)) return(NULL)
  pstar <- nrow(fl) * Q / P[tested == TRUE, .N]                                                    # k x 0.05 / m (3.2(3))
  U <- list(Up = setdiff(univ, H_up), Down = setdiff(univ, H_dn)); UI <- setdiff(univ, IEG); UF <- setdiff(univ, H)
  R <- lapply(U, function(u) cam(y, X, sets, u, min_n = MINREM)); RI <- cam(y, X, sets, UI, min_n = MINREM); RF <- cam(y, X, sets, UF, min_n = MINREM)
  read_one <- function(d0, s0, z, s2) {
    if (is.na(s0) || s0 <= 0) return("shift opposes the camera direction")
    if (is.na(z$NGenes) || z$NGenes < MINREM) return("consists mostly of the listed handling genes")
    if (is.na(s2) || s2 < HALF * s0) return("carried by the listed handling genes")
    if (!is.na(z$PValue) && z$PValue <= pstar && identical(z$Direction, d0)) return("not carried by the listed handling genes")
    "not resolved after removal" }
  rbindlist(lapply(fl$set, function(sn) {
    d0 <- fl[set == sn, Direction]; s0 <- sgn(d0) * shift_of(lfc, sets[[sn]], univ); u <- U[[d0]]
    if (is.na(s0) || s0 <= 0) return(data.table(family = fam, set = sn, Direction = d0, s_primary = s0, pstar = pstar, opposes = TRUE,   # 3.1: marked, no further reading;
                                                reading = "shift opposes the camera direction", bound = FALSE))                   # not a 7.3 unit (dated note)
    z <- R[[d0]][set == sn]; zi <- RI[set == sn]; zf <- RF[set == sn]
    rd <- read_one(d0, s0, z, sgn(d0) * shift_of(lfc, sets[[sn]], u)); rdi <- read_one(d0, s0, zi, sgn(d0) * shift_of(lfc, sets[[sn]], UI))
    rdf <- read_one(d0, s0, zf, sgn(d0) * shift_of(lfc, sets[[sn]], UF))
    sv <- intersect(intersect(sets[[sn]], univ), SV$gene); agree <- if (length(sv) >= 8) mean(sign(lfc[sv]) == ifelse(SV$direction[match(sv, SV$gene)] == "up", 1, -1)) else NA_real_
    dchk <- if (is.na(agree)) "too few delay members" else if (agree >= 0.75) "runs with the delay direction" else if (agree <= 0.25) "runs against it" else "mixed"
    carrier <- if (rd == "carried by the listed handling genes") { if (rdi == "carried by the listed handling genes") "immediate-early genes" else "delay genes" } else ""
    full_rd <- rd; if (rd == "carried by the listed handling genes") full_rd <- if (dchk == "runs against it")   # precedence (3.2(3)): the fixed text, no carrier
      "carried by the listed handling genes; its delay-list members as a whole run against the delay direction" else sprintf("carried by the listed handling genes (%s)", carrier)
    data.table(family = fam, set = sn, Direction = d0, s_primary = s0, pstar = pstar, opposes = FALSE, carrier = carrier, n_left = z$NGenes, p_after = z$PValue, dir_after = z$Direction,
               s_after = sgn(d0) * shift_of(lfc, sets[[sn]], u), reading = full_rd, bound = rd %in% c("carried by the listed handling genes", "consists mostly of the listed handling genes"),
               reading_ieg_only = rdi, reading_full_removal = rdf, n_savage = length(sv), delay_agreement = agree, direction_check = dchk) }), fill = TRUE) }   # an opposing set's row is shorter
if (MODE == "sets") {
  C <- gate("S"); recheck_inputs()
  ts("[sets] primary T fit (S5 as run)")
  r <- fit_s5(RAW$T, D$T); lfc <- r$lfc; e <- setNames(r$y$AveLogCPM, rownames(r$y))
  FAM <- list(`P-H` = HM, `P-T` = GV19); OUT <- list(); REM <- list(); HS <- stream("h_rand", length(HM) + length(GV19))
  BIN <- bins_of(e[intersect(names(e), UA)])
  for (fm in names(FAM)) {
    P <- cam(r$y, r$X, FAM[[fm]], UA)
    P[, `:=`(family = fm, size = lengths(FAM[[fm]])[set], in_UA = vapply(FAM[[fm]][set], function(s) length(intersect(s, UA)), 1L))]
    P[, s := vapply(set, function(z) sgn(Direction[set == z]) * shift_of(lfc, FAM[[fm]][[z]], UA), 1)]
    P[, `:=`(coverage = in_UA / size, opposes = !is.na(s) & s <= 0)][, mark := fifelse(opposes, "shift opposes the camera direction", "")]   # 3.1: such a set gets no further reading
    P[, fH := vapply(set, function(z) { g <- intersect(FAM[[fm]][[z]], UA); length(intersect(g, H)) / length(g) }, 1)]
    P[, c("fH_rand_median", "fH_rand_p95", "fH_pct") := rbindlist(lapply(set, function(z) { g <- intersect(FAM[[fm]][[z]], UA)
      rr <- vapply(match_draw(g, names(BIN), BIN, N_RAND, HS[(if (fm == "P-H") 0L else length(HM)) + match(z, set)]), function(x) length(intersect(x, H)) / max(1, length(x)), 1)
      data.table(median(rr), quantile(rr, 0.95), pct_rank(length(intersect(g, H)) / length(g), rr)) }))]
    for (arm in c("up", "down")) { gp <- P370[s == if (arm == "up") 1 else -1, gene]; gf <- F402[s == if (arm == "up") 1 else -1, gene]
      P[, paste0("P370_", arm) := vapply(set, function(z) length(intersect(intersect(FAM[[fm]][[z]], UA), gp)), 1L)]
      P[, paste0("P370_", arm, "_exp") := in_UA * length(gp) / length(UA)]
      P[, paste0("F402_", arm) := vapply(set, function(z) length(intersect(intersect(FAM[[fm]][[z]], UT), gf)), 1L)]
      P[, paste0("F402_", arm, "_exp") := vapply(set, function(z) length(intersect(FAM[[fm]][[z]], UT)), 1L) * length(gf) / length(UT)] }
    if (fm == "P-T") P[, absent_from_UA := vapply(set, function(z) paste(setdiff(FAM[[fm]][[z]], UA), collapse = " "), "")]
    OUT[[fm]] <- P; REM[[fm]] <- removal(fm, FAM[[fm]], P, r$y, r$X, lfc, UA) }
  UNA <- cam(r$y, r$X, GV[UNASSIGNED], UA)[, `:=`(family = "UNASSIGNED (reported, outside every BH family)", FDR = NA_real_)]
  UNA[, `:=`(size = lengths(GV[UNASSIGNED])[set], in_UA = vapply(GV[UNASSIGNED][set], function(z) length(intersect(z, UA)), 1L))][, coverage := in_UA / size]
  UNA[, s := vapply(seq_len(.N), function(i) sgn(Direction[i]) * shift_of(lfc, GV[[set[i]]], UA), 1)]
  UNA[, `:=`(fH = vapply(set, function(z) { g <- intersect(GV[[z]], UA); if (length(g)) length(intersect(g, H)) / length(g) else NA_real_ }, 1),
             absent_from_UA = vapply(set, function(z) paste(setdiff(GV[[z]], UA), collapse = " "), ""))]
  ts("[sets] 3.2(6) opposite-sign class; section 4 views and sensitivities")
  OPP <- rbind(cam(r$y, r$X, list(opposite12 = OPP12), UA, min_n = 5L), cam(r$y, r$X, list(opposite10 = setdiff(OPP12, c("KLRD1", "NKG7"))), UA, min_n = 5L))
  OPP[, reading := ifelse(Direction == "Up", "runs against the delay direction, as an effector shift would", "runs with the delay direction")]
  SENS <- list()
  for (fm in names(FAM)) {
    SENS[[paste(fm, "4.1 U_T")]] <- cam(r$y, r$X, FAM[[fm]], UT)[, view := "4.1 all 7,379 genes"]
    SENS[[paste(fm, "4.4 camera 0.01")]] <- cam(r$y, r$X, FAM[[fm]], UA, cor = 0.01)[, view := "4.4 camera 0.01"]
    yy <- r$y[rownames(r$y) %in% UA, , keep.lib.sizes = TRUE]; idx <- ids2indices(FAM[[fm]], rownames(yy)); idx <- idx[lengths(idx) >= MINSET]
    fr <- as.data.table(fry(yy, idx, r$X, contrast = "statusAML", sort = "none"), keep.rownames = "set")
    SENS[[paste(fm, "4.4 fry")]] <- fr[, .(set, NGenes, Direction, PValue, FDR, view = "4.4 fry directional")]
    st <- r$tt[gene %in% UA, setNames(sign(logFC) * sqrt(F), gene)]; set.seed(SEED + OFF[["fgsea"]])
    fg <- as.data.table(fgsea::fgseaMultilevel(lapply(FAM[[fm]], intersect, UA), st, minSize = MINSET, maxSize = 500, eps = 0))
    SENS[[paste(fm, "4.4 fgsea")]] <- fg[, .(set = pathway, NGenes = size, Direction = ifelse(NES > 0, "Up", "Down"), PValue = pval, FDR = padj, NES,
                                            leading_edge = vapply(leadingEdge, paste, "", collapse = " "), view = "4.4 fgsea preranked")] }
  J <- rbind(OUT$`P-H`[tested == TRUE, .(family, set, PValue)], OUT$`P-T`[tested == TRUE, .(family, set, PValue)])[, FDR_joint := p.adjust(PValue, "BH")]
  SENS$joint <- J[, .(set, family, PValue, FDR = FDR_joint, view = "4.4 joint BH over P-H and P-T")]
  d42 <- copy(D$T); stopifnot(identical(levels(d42$sex), c("female", "male"))); if (!SMOKE) d42[sample == "Control4", sex := factor("female", levels(sex))]
  r42 <- fit_s5(RAW$T, d42); d43 <- D$T[!sample %in% T15]; r43 <- fit_s5(cols(RAW$T, d43$sample), d43)
  view_rem <- function(fm, p, f, univ, lab) { rr <- removal(fm, FAM[[fm]], p, f$y, f$X, f$lfc, univ)   # 3.2 under a section-4 view (reported)
    if (is.null(rr)) data.table(set = character(0)) else rr[, .(set, Direction, reading, view = lab)] }
  for (fm in names(FAM)) { SENS[[paste(fm, "4.2")]] <- cam(r42$y, r42$X, FAM[[fm]], setdiff(rownames(r42$y), EXT))[, view := "4.2 Control4 coded female"]
                           SENS[[paste(fm, "4.3")]] <- cam(r43$y, r43$X, FAM[[fm]], setdiff(rownames(r43$y), EXT))[, view := "4.3 without Control1-5"]
    SENS[[paste(fm, "4.1 removal")]] <- view_rem(fm, SENS[[paste(fm, "4.1 U_T")]], r, UT, "4.1 removal reading")
    SENS[[paste(fm, "4.2 removal")]] <- view_rem(fm, SENS[[paste(fm, "4.2")]], r42, setdiff(rownames(r42$y), EXT), "4.2 removal reading")
    SENS[[paste(fm, "4.3 removal")]] <- view_rem(fm, SENS[[paste(fm, "4.3")]], r43, setdiff(rownames(r43$y), EXT), "4.3 removal reading") }
  SN <- rbindlist(lapply(names(SENS), function(k) SENS[[k]][, family_view := k]), fill = TRUE)
  LE <- rbindlist(lapply(names(FAM), function(fm) SENS[[paste(fm, "4.4 fgsea")]][, .(set, leading_edge)]))
  ALL <- rbind(rbindlist(OUT, fill = TRUE), UNA, fill = TRUE); ALL[, leading_edge := LE$leading_edge[match(set, LE$set)]]
  put(ALL, "sets"); put(if (length(REM)) rbindlist(REM, fill = TRUE) else data.table(), "sets_removal"); put(SN, "sets_sens"); put(OPP, "opposite")
  ts(sprintf("[sets] flagged: P-H %d, P-T %d", OUT$`P-H`[FDR < Q, .N], OUT$`P-T`[FDR < Q, .N])); quit(status = 0)
}

## -- 9. --cross: 3.3 ----------------------------------------------------------------------------------------------------------------------
if (MODE == "cross") {
  C <- gate("X"); recheck_inputs(); K5ok <- stat_ok(C, "K5")
  S <- rtab("sets"); FL <- S[family %in% c("P-H", "P-T") & FDR < Q & opposes == FALSE]   # an opposing set gets no further reading (3.1)
  FAM <- c(HM, GV)
  rc4 <- function(d) { d <- copy(d); if (!SMOKE && "Control4" %in% d$sample) {   # 4.2: Control4 coded female (cells are female; memory note)
    d[, sex := factor(as.character(sex), levels = c("female", "male"))]; d[sample == "Control4", sex := factor("female", levels = c("female", "male"))] }; d }
  ts("[cross] 3.3(a): each compartment alone, and the delay lists (primary; views 4.1 all tested genes and 4.2 Control4 coded female)")
  DL <- lapply(c(T = "T", B = "B", NK = "NK", Mono = "Mono"), function(k) fread(DELAY_F(k)))
  CMPT <- list()
  for (k in CMPS) { f <- fit_s5(RAW[[k]], D[[k]]); f42 <- fit_s5(RAW[[k]], rc4(D[[k]]))
    dl <- list(own_up = DL[[k]][direction == "up", gene], own_down = DL[[k]][direction == "down", gene], T_up = DL$T[direction == "up", gene], T_down = DL$T[direction == "down", gene])
    if (k == "T") dl <- dl[1:2]
    VW <- list(primary = list(f = f, U = setdiff(rownames(f$y), EXT)), `4.1 all tested genes` = list(f = f, U = rownames(f$y)),
               `4.2 Control4 coded female` = list(f = f42, U = setdiff(rownames(f42$y), EXT)))
    for (vw in names(VW)) { ff <- VW[[vw]]$f; U <- VW[[vw]]$U
      for (fm in c("P-H", "P-T")) CMPT[[paste(k, fm, vw)]] <- cam(ff$y, ff$X, if (fm == "P-H") HM else GV19, U)[, `:=`(compartment = k, family = fm, view = vw)]
      CMPT[[paste(k, "delay", vw)]] <- cam(ff$y, ff$X, dl, U, min_n = MINREM)[, `:=`(compartment = k, family = "delay lists", view = vw)]
      if (k != "Mono") CMPT[[paste(k, "UNA", vw)]] <- cam(ff$y, ff$X, GV[UNASSIGNED], U)[, `:=`(compartment = k, family = "UNASSIGNED (own p only)", FDR = NA_real_, view = vw)] } }
  CM <- rbindlist(CMPT, fill = TRUE)
  CM[, reading := ""][family == "delay lists" & view == "primary", reading := "reported (not read)"]   # Mono and the T lists inside B/NK/Mono: reported, not read (3.3(a))
  CM[family == "delay lists" & view == "primary" & compartment %in% c("T", "B", "NK") & set %in% c("own_up", "own_down"),   # 7.6 wording, a compartment and its own list
     reading := fifelse(!is.na(FDR) & FDR < Q & ((set == "own_up" & Direction == "Up") | (set == "own_down" & Direction == "Down")),
                        sprintf("%s aligned with its own delay list in the delay direction (consistent with tube-level delay)", compartment),
                        sprintf("%s: alignment with its delay list not detected (power not established)", compartment))]
  ts("[cross] 3.3(b): stacked interaction, rho_X with its donor bootstrap (primary); 4.1 and 4.2 views (interaction camera and rho_X point estimate)")
  CR <- list(); set.seed(SEED + OFF[["boot"]])   # one stream, drawn before any parallel work: B's resamples first, then NK's
  BSALL <- lapply(c(B = "B", NK = "NK"), function(X) { dT <- D$T[match(paired(X), sample)]; a <- dT[status == "AML", sample]; h <- dT[status == "healthy", sample]
    lapply(seq_len(N_BOOT), function(b) c(sample(a, replace = TRUE), sample(h, replace = TRUE))) })
  inter <- function(y, Xd, U) rbind(cam(y, Xd, HM, U, contrast = "cT:AML")[, family := "P-H"], cam(y, Xd, GV19, U, contrast = "cT:AML")[, family := "P-T"],
                                    cam(y, Xd, GV[UNASSIGNED], U, contrast = "cT:AML")[, `:=`(family = "UNASSIGNED", FDR = NA_real_)])
  for (X in c("B", "NK")) {
    dons <- paired(X); com <- setdiff(intersect(UT, TESTED[[X]]), EXT); com41 <- intersect(UT, TESTED[[X]]); if (SMOKE) { com <- intersect(com, G); com41 <- intersect(com41, G) }
    Xm <- stack_design(X, dons); dT <- D$T[match(dons, sample)]; dX <- D[[X]][match(dons, sample)]; MT <- cols(RAW$T, dons); MX <- cols(RAW[[X]], dons)
    INT <- inter(stack_fit(stack_counts(X, MT, MX, dons, com), Xm), Xm, com)
    stopifnot(INT[family == "P-H", sum(tested)] == c(B = 41L, NK = 43L)[[X]], INT[family == "P-T", sum(tested)] == c(B = 15L, NK = 18L)[[X]])   # 3.3(b) BH families
    rp <- rho_pair(MT, MX, dT, dX, com)
    dT42 <- rc4(dT); dX42 <- rc4(dX); Xm42 <- stack_design(X, dons, dT42, dX42)
    VW <- list(`4.1 all common genes` = list(INT = inter(stack_fit(stack_counts(X, MT, MX, dons, com41), Xm), Xm, com41), rp = rho_pair(MT, MX, dT, dX, com41), com = com41),
               `4.2 Control4 coded female` = list(INT = inter(stack_fit(stack_counts(X, MT, MX, dons, com), Xm42), Xm42, com), rp = rho_pair(MT, MX, dT42, dX42, com), com = com))
    rho_of <- function(lt, lx, s, d0, cc = com) { sT <- sgn(d0) * shift_of(lt, FAM[[s]], cc); sX <- sgn(d0) * shift_of(lx, FAM[[s]], cc); if (is.na(sT) || sT <= 0) NA_real_ else sX / sT }
    BS <- BSALL[[X]]
    BR <- mcl(seq_len(N_BOOT), function(b) { s <- BS[[b]]; nm <- paste0(s, "_", seq_along(s))
      mt <- rsub(cols(MT, s), com); mx <- rsub(cols(MX, s), com); colnames(mt) <- colnames(mx) <- nm   # rsub/cols keep the smoke mark
      dt <- copy(dT[match(s, sample)])[, sample := nm]; dx <- copy(dX[match(s, sample)])[, sample := nm]
      out <- try(rho_pair(mt, mx, dt, dx, com), silent = TRUE); if (inherits(out, "try-error")) list(T = NULL, X = NULL) else out })
    for (i in seq_len(nrow(FL))) { sn <- FL$set[i]; d0 <- FL$Direction[i]; z <- INT[set == sn]
      n_common <- length(intersect(FAM[[sn]], com)); rho <- rho_of(rp$T, rp$X, sn, d0)
      bt <- vapply(BR, function(o) if (is.null(o$T)) NA_real_ else rho_of(o$T, o$X, sn, d0), 1); nund <- sum(is.na(bt))
      bt[is.na(bt)] <- -Inf   # a failed refit or s_T <= 0 counts against "shared" (dated note): the bound is over all 1,000 draws
      lo <- unname(quantile(bt, 0.025)); hi <- unname(quantile(bt, 0.975))
      rd <- if (!K5ok) "not readable (K5 failed)" else if (n_common < MINSET) sprintf("not testable in %s", X) else if (!is.na(z$FDR) && z$FDR < Q && z$Direction == d0) sprintf("larger in T than in %s", X) else
        if (!is.na(z$FDR) && z$FDR < Q) sprintf("larger in %s than in T", X) else if (!is.na(rho) && rho >= RHO_SHARED && lo > 0) sprintf("shared with %s", X) else "not resolved"
      CR[[paste(X, sn, "primary")]] <- data.table(view = "primary", X = X, family = FL$family[i], set = sn, Direction = d0, n_common = n_common, int_dir = z$Direction, int_p = z$PValue,
                                                 int_FDR = z$FDR, rho = rho, rho_lo = lo, rho_hi = hi, boot_undefined = nund, reading = rd)
      for (vw in names(VW)) { v <- VW[[vw]]; zz <- v$INT[set == sn]   # views: reported, no reading word
        CR[[paste(X, sn, vw)]] <- data.table(view = vw, X = X, family = FL$family[i], set = sn, Direction = d0, n_common = length(intersect(FAM[[sn]], v$com)),
                                            int_dir = zz$Direction, int_p = zz$PValue, int_FDR = zz$FDR, rho = rho_of(v$rp$T, v$rp$X, sn, d0, v$com), reading = "") } }
    CR[[paste(X, "all")]] <- INT[, `:=`(view = "primary", X = X, set_row = "interaction, all sets")]
    for (vw in names(VW)) CR[[paste(X, "all", vw)]] <- VW[[vw]]$INT[, `:=`(view = vw, X = X, set_row = "interaction, all sets")] }
  put(CM, "compartments"); put(rbindlist(CR, fill = TRUE), "cross"); ts("[cross] done"); quit(status = 0)
}

## -- 10. --tags: 3.4 ---------------------------------------------------------------------------------------------------------------------
if (MODE == "tags") {
  C <- gate("G"); recheck_inputs(); K4ok <- stat_ok(C, "K4")
  need_syn(RAW$T, "whole-T count"); for (k in TL) need_syn(FINE[[k]], "label count")
  OJ <- load_oetjen(); main <- tags_both(OJ)
  core <- fread(file.path(TABWS, "ws_a1_g02_gene_selection.csv"))[n_folds == 25L, gene]; stopifnot(length(core) == 25L, all(core %in% P370$gene))
  ARMS <- list(P370_up = P370[s > 0, gene], P370_down = P370[s < 0, gene], core25 = core, F402_up = F402[s > 0, gene], F402_down = F402[s < 0, gene])
  assessable <- intersect(UT, OJ$axis)
  Lh <- log2(rowSums(RAW$T[, HEAL]) / sum(RAW$T[, HEAL]) * 1e6 + 1); BINT <- bins_of(Lh[assessable])   # pooled logCPM of the healthy donors only (no AML expression enters 3.4; dated note)
  ## composition-only expectation (3.4): UMI-weighted label mix at a common CD8 share, centred on the U_A median
  Pk <- vapply(TL, function(k) { v <- rowSums(FINE[[k]][, HEAL, drop = FALSE]); v / max(1, sum(v)) * 1e6 }, numeric(length(G))); rownames(Pk) <- G
  ub <- vapply(TL, function(k) sum(FINE[[k]][, HEAL]) / max(1, sum(NCELL[k, HEAL])), 1)
  shr <- function(sm) rowMeans(sweep(NCELL[, sm, drop = FALSE], 2, colSums(NCELL[, sm, drop = FALSE]), "/"))
  cbar <- mean(colSums(NCELL[CD8L, S25]) / colSums(NCELL[, S25]))
  wts <- function(s) { s[CD4L] <- s[CD4L] / sum(s[CD4L]) * (1 - cbar); s[CD8L] <- s[CD8L] / sum(s[CD8L]) * cbar; w <- s * ub; w / sum(w) }
  wA <- wts(shr(AMLS)); wH <- wts(shr(HEAL)); comp <- log2((Pk %*% wA + 1) / (Pk %*% wH + 1))[, 1]; comp <- comp - median(comp[intersect(names(comp), UA)])
  readout <- function(TB, tag_view) {
    cons <- merge(TB$cG, TB$cO, by = "gene", suffixes = c("_G", "_O"))
    for (cl in names(CLASSES)) cons[, (cl) := get(paste0(cl, "_G")) & get(paste0(cl, "_O"))]
    rbindlist(lapply(names(ARMS), function(an) {
      uni <- intersect(if (grepl("^F402", an)) UT else UA, assessable); g <- intersect(ARMS[[an]], uni)
      RS <- match_draw(g, uni, BINT, N_RAND, stream("tags_rand", length(ARMS))[match(an, names(ARMS))])
      rbindlist(lapply(names(CLASSES), function(cl) {
        tagged <- cons[get(cl) == TRUE, gene]; obs <- mean(g %in% tagged); base <- mean(uni %in% tagged); rr <- vapply(RS, function(x) mean(x %in% tagged), 1)
        tagG <- cons[get(paste0(cl, "_G")) == TRUE, gene]; tagO <- cons[get(paste0(cl, "_O")) == TRUE, gene]   # each reference alone, on the genes assessable in both
        tg <- intersect(g, tagged); cg <- comp[tg]; sgt <- F402$s[match(tg, F402$gene)]   # each gene's own direction (the core mixes up and down genes)
        data.table(view = tag_view, arm = an, class = cl, n_arm = length(ARMS[[an]]), n_assessable = length(g), not_assessable = length(setdiff(ARMS[[an]], uni)),
                   n_consensus = length(tg), share = obs, base_rate = base, random_pct = pct_rank(obs, rr),
                   n_tag_G = sum(g %in% tagG), base_rate_G = mean(uni %in% tagG), random_pct_G = pct_rank(mean(g %in% tagG), vapply(RS, function(x) mean(x %in% tagG), 1)),
                   n_tag_O = sum(g %in% tagO), base_rate_O = mean(uni %in% tagO), random_pct_O = pct_rank(mean(g %in% tagO), vapply(RS, function(x) mean(x %in% tagO), 1)),
                   comp_p10 = if (length(cg)) quantile(cg, 0.1) else NA_real_, comp_p50 = if (length(cg)) median(cg) else NA_real_, comp_p90 = if (length(cg)) quantile(cg, 0.9) else NA_real_,
                   comp_share_arm_sign = if (length(cg)) mean(sign(cg) == sgt & abs(cg) >= 0.25) else NA_real_,
                   genes = if (grepl("^F402", an)) "" else paste(tg, collapse = " ")) })) })) }   # F402 values appear only as numbers (2.1): no gene names
  TG <- readout(main, "primary (threshold 1, both references)")
  TG <- rbind(TG, readout(tags_both(OJ, thr = 0.585), "sensitivity: threshold 0.585"))
  g5 <- setdiff(HEAL, T15); TG <- rbind(TG, readout(tags_both(OJ, cmin = 0.75, min_don = 4L, gdon = g5), "sensitivity: G without Control1-5"))
  TG[, consensus_readable := K4ok]
  stopifnot(!any(nzchar(TG[grepl("^F402", arm), genes])))   # no F402 gene name reaches the public table
  fwrite(rbind(main$G[, ref := "G"], main$O[, ref := "O"])[, script_md5 := MD5], wg("tags_per_gene.csv.gz")); put(TG, "tags")   # work file first (/LARGE1)
  ts("[tags] done"); quit(status = 0)
}

## -- 11. --classes: 3.6(a) ------------------------------------------------------------------------------------------------------------
if (MODE == "classes") {
  C <- gate("C"); recheck_inputs()
  CG <- union(CLOCK, GCLIST)
  rows <- list(none = character(0), Savage = SV$gene, `Savage + IEG` = c(SV$gene, IEG), lncRNA = LNC, `clock + glucocorticoid` = CG, `lncRNA + clock/glucocorticoid + H` = c(LNC, CG, H))
  g0 <- oof_gap()$gap
  CL <- rbindlist(lapply(names(rows), function(k) { g <- oof_gap(rows[[k]])$gap
    data.table(removed = k, gap_oof = g, ratio_to_none = g / g0, sentence = if (k == "none") sprintf("The cross-fitted gap with no genes removed is %.3f.", g) else
      sprintf("The cross-fitted gap is %.3f without %s (ratio %.2f to %.3f).", g, k, g / g0, g0)) }))
  LS <- rbindlist(lapply(c("up", "down"), function(arm) { gp <- P370[s == if (arm == "up") 1 else -1, gene]; gf <- F402[s == if (arm == "up") 1 else -1, gene]
    data.table(arm = arm, P370_lnc = sum(gp %in% LNC), P370_n = length(gp), P370_lnc_expected = length(gp) * mean(UA %in% LNC),
               F402_lnc = sum(gf %in% LNC), F402_n = length(gf), F402_lnc_expected = length(gf) * mean(UT %in% LNC)) }))
  CL[, clock_gc_note := fifelse(removed == "clock + glucocorticoid", "this cannot exclude time of day or steroids, which act beyond these genes and are unrecorded", "")]
  put(rbind(CL[, table := "cross-fitted gap without a gene class"], LS[, table := "lncRNA share per arm"], fill = TRUE), "classes"); ts("[classes] done"); quit(status = 0)
}

## -- 12. --overlap: 3.6(b) ------------------------------------------------------------------------------------------------------------
if (MODE == "overlap") {
  C <- gate("O"); I <- recheck_inputs()
  suppressPackageStartupMessages(library(readxl))
  core <- fread(file.path(TABWS, "ws_a1_g02_gene_selection.csv"))[n_folds == 25L, gene]; stopifnot(length(core) == 25L)
  LISTS <- c(list(`Savage 2021 delay up` = SV[direction == "up", gene], `Savage 2021 delay down` = SV[direction == "down", gene], `project IEG list` = IEG),
             setNames(HM, paste("HALLMARK:", names(HM))), setNames(GV, paste("Gavish 3CA:", names(GV))))
  cn <- fread(CNMF_F); LISTS[["cNMF normal MP10 (derived from this project's cells)"]] <- cn[MP == "MP10", gene]; LISTS[["cNMF normal MP4 (derived from this project's cells)"]] <- cn[MP == "MP4", gene]
  la <- as.data.table(read_excel(LASRY_X, sheet = "Supp6. Adult inflammation sig")); li <- as.data.table(read_excel(LASRY_X, sheet = "Supp11. Adult iScore genes"))
  LISTS[["Lasry 2023 adult inflammation signature (malignant cells; positioning only)"]] <- unique(la$gene_id); LISTS[["Lasry 2023 adult iScore (bulk; positioning only)"]] <- unique(li$`Gene name`)
  ru <- as.data.table(suppressMessages(read_excel(RUTELLA_X, sheet = "Supplementary Table 2", skip = 1)))
  for (j in seq_len(ncol(ru))) { v <- unique(na.omit(as.character(ru[[j]]))); v <- v[grepl("^[A-Za-z0-9][A-Za-z0-9.-]*$", v)]; if (length(v)) LISTS[[paste("Rutella 2022:", names(ru)[j])]] <- v }
  mz_disk <- list.files(MAZZ_DIR, "\\.csv$", full.names = TRUE); mz <- intersect(mz_disk, I$path); mz_late <- setdiff(mz_disk, mz)   # only csvs md5-recorded at --known
  mz_note <- paste(Filter(nzchar, c(if (length(mz)) "" else if (length(list.files(MAZZ_DIR, "\\.pdf$"))) "Mazziotta 2024: PDF present but no extracted gene list (csv) at --known; omitted" else
                                      "Mazziotta 2024: no extracted gene list (csv) at --known (pending manual download); omitted",
                                    if (length(mz_late)) sprintf("Mazziotta 2024: %s added after --known; omitted (dated note)", paste(basename(mz_late), collapse = ", ")) else "")), collapse = "; ")
  for (f in mz) LISTS[[paste("Mazziotta 2024:", basename(f))]] <- unique(fread(f)[[1]])
  up <- P370[s > 0, gene]; dn <- P370[s < 0, gene]; fu <- F402[s > 0, gene]; fd <- F402[s < 0, gene]
  OV <- rbindlist(lapply(names(LISTS), function(k) { g <- LISTS[[k]]; ua <- intersect(g, UA); ut <- intersect(g, UT)
    data.table(list = k, size = length(g), in_UT = length(ut), in_UA = length(ua), P370_up = length(intersect(ua, up)), P370_up_expected = length(ua) * length(up) / length(UA),
               P370_down = length(intersect(ua, dn)), P370_down_expected = length(ua) * length(dn) / length(UA), core25 = length(intersect(g, core)),
               core25_expected = length(ua) * length(core) / length(UA),
               F402_up = length(intersect(ut, fu)), F402_up_expected = length(ut) * length(fu) / length(UT), F402_down = length(intersect(ut, fd)),
               F402_down_expected = length(ut) * length(fd) / length(UT), P370_genes = paste(sort(intersect(ua, P370$gene)), collapse = " ")) }))
  OV[, note := mz_note][]
  put(OV, "overlap"); ts("[overlap] done"); quit(status = 0)
}

## -- 13. --combine: section 7 ---------------------------------------------------------------------------------------------------------
if (MODE == "combine") {
  C <- gate("S|X|G|C|O"); recheck_inputs()
  S <- rtab("sets"); RM <- rtab("sets_removal"); SN <- rtab("sets_sens"); OP <- rtab("opposite"); CM <- rtab("compartments"); CR <- rtab("cross")
  TG <- rtab("tags"); CL <- rtab("classes"); OV <- rtab("overlap")
  md5s <- unique(c(S$script_md5, RM$script_md5, SN$script_md5, OP$script_md5, CM$script_md5, CR$script_md5, TG$script_md5, CL$script_md5, OV$script_md5, C$script_md5))
  if (!SMOKE && !identical(md5s, MD5)) stop("tables from different script versions: ", paste(md5s, collapse = " "))
  f5 <- file.path(ODIR, paste0("ws_a1_g07_s5b_patients", SFX, ".csv")); f5c <- file.path(ODIR, paste0("ws_a1_g07_s5b_checks", SFX, ".csv"))
  if (!file.exists(f5) || !file.exists(f5c)) stop("run 44 --known and 44 --score first: the GATES block is appended once and must hold 3.5")
  S5B <- fread(f5); S5C <- fread(f5c); MD544 <- unique(S5B$script_md5)
  if (length(MD544) != 1L || !identical(unique(S5C$script_md5), MD544) || any(S5C$state == "FAIL_deterministic")) stop("the S5b tables are not from one version of 44 with every deterministic check passed")
  K2ok <- stat_ok(C, c("K2a", "K2b")); FAMG <- c(HM, GV)
  FLA <- S[family %in% c("P-H", "P-T") & FDR < Q]; FL <- FLA[opposes == FALSE]; FLX <- FLA[opposes == TRUE]   # opposing sets: the mark only (3.1)
  pw <- C[id == "K2c", got][1]
  NOTCAR <- " (the lists represent 18 h ex-vivo delay of fresh blood from two donors and immediate-early stress; freezing, thawing, shipping and sorting are not represented)"
  COMPX <- " (two compartments of the same donors; monocytes and other cells are not compared this way, so this is not a statement of T-specificity)"
  ## 7.1 / 7.3 / 7.6 per flagged set
  name_of <- function(p) if (p$family == "P-T") sprintf(paste("the Gavish program \"%s\" (a program defined in tumour-infiltrating T cells; as tested here %d of %d genes, %d of them in the handling lists;",
                                                             "tested on whole-T pseudobulks, so the CD4/CD8 prefix does not localise the shift)"), p$set, p$in_UA, p$size, round(p$fH * p$in_UA)) else p$set
  lost_of <- function(p) {   # section 4: one parenthesis naming any of 4.1-4.3 under which the set loses FDR < 0.05 or flips direction; 4.3 has its fixed text
    sens <- SN[set == p$set & view %in% c("4.1 all 7,379 genes", "4.2 Control4 coded female", "4.3 without Control1-5")]; out <- character(0)
    for (v in c("4.1 all 7,379 genes", "4.2 Control4 coded female")) { z <- sens[view == v][1]; if (is.na(z$set)) next
      if (is.na(z$FDR) || z$FDR >= Q) out <- c(out, sprintf("loses FDR < 0.05 under %s", v)) else if (!identical(z$Direction, p$Direction)) out <- c(out, sprintf("direction flips under %s", v)) }
    z <- sens[view == "4.3 without Control1-5"][1]
    if (!is.na(z$set)) { if (is.na(z$FDR) || z$FDR >= Q) out <- c(out, "not significant without Control1-5, 5 healthy donors") else
      if (!identical(z$Direction, p$Direction)) out <- c(out, "direction flips without Control1-5, 5 healthy donors") }
    if (length(out)) sprintf(" (%s)", paste(out, collapse = "; ")) else "" }
  hand_of <- function(p) { r <- RM[set == p$set & family == p$family]   # 7.3: the 3.2 reading, the H fraction with its random reference, the direction check
    hf <- sprintf("H fraction %s (matched-random median %s, P95 %s; percentile %s)", f3(p$fH), f3(p$fH_rand_median), f3(p$fH_rand_p95), f3(p$fH_pct))
    if (!nrow(r)) return(hf)
    paste0(r$reading, if (r$reading == "not carried by the listed handling genes") NOTCAR else "", "; ", hf, "; direction check: ", r$direction_check) }
  comp_of <- function(p) { x <- CR[set == p$set & view == "primary" & !is.na(reading) & nzchar(reading)]; if (!nrow(x)) return("")   # fread returns "" for an empty field
    rr <- ifelse(grepl("^shared with", x$reading), sprintf("%s (\u03c1 %.2f, 95 %% CI %.2f to %.2f)", x$reading, x$rho, x$rho_lo, x$rho_hi), x$reading)
    paste0(paste(sprintf("T versus %s: %s", x$X, rr), collapse = "; "), COMPX) }
  set_sent <- function(i) { p <- FL[i]; cx <- comp_of(p)
    sprintf(paste("In the S5 T-cell contrast (15 adult AML vs 10 healthy donors, GSE185381), %s had a %s average AML-versus-healthy statistic than the average tested gene",
                  "(camera with estimated inter-gene correlation, FDR %s; same contrast that defined the axis; descriptive)%s. Handling (7.3): %s.%s"),
            name_of(p), if (p$Direction == "Up") "higher" else "lower", fq(p$FDR), lost_of(p), hand_of(p), if (nzchar(cx)) sprintf(" Compartments (7.6): %s.", cx) else "") }
  opp_sent <- function(i) { p <- FLX[i]
    sprintf("%s was flagged by camera (%s, FDR %s), but its set shift opposes the camera direction, so it gets no further reading (3.1).", name_of(p), p$Direction, fq(p$FDR)) }
  SETS <- if (!K2ok) "Set-level description not interpretable (K2a or K2b failed); tables are reported as computed." else if (!nrow(FLA))
    sprintf(paste("No HALLMARK or T-cell-program set differed from the average tested gene at FDR 0.05 in a competitive test that accounts for inter-gene correlation (n = 25).",
                  "On this fit, planted whole-set shifts of 0.3 log2 were flagged in %s of %s, so this is not evidence that the axis lacks set structure."), strsplit(pw, ";")[[1]][1], strsplit(pw, ";")[[1]][3]) else
    paste(c(vapply(seq_len(nrow(FL)), set_sent, ""), vapply(seq_len(nrow(FLX)), opp_sent, "")), collapse = " ")
  ## 7.2 secondary views
  only44 <- setdiff(SN[grepl("^4\\.4", view) & !is.na(FDR) & FDR < Q, unique(set)], FLA$set)
  SEC <- if (length(only44)) sprintf("Flagged only in a 4.4 view, so not detected by the primary: %s.", paste(only44, collapse = ", ")) else "No set was flagged only in a 4.4 view."
  ## 7.3 consequence rule: per direction over the non-opposing flagged sets, a CD4/CD8 Gavish pair as one unit
  unit_of <- function(s) { for (p in PAIRS_EXPECTED) if (s %in% strsplit(p, "\\|")[[1]]) return(p); s }
  CONS <- if (!K2ok) "not evaluated (3.1 not interpretable: K2a or K2b failed)" else if (!nrow(FL)) "not triggered (no flagged set whose shift runs in its camera direction)" else
    paste(unlist(lapply(c("Up", "Down"), function(d0) { fs <- FL[Direction == d0, set]; if (!length(fs)) return(NULL)
      r <- RM[set %in% fs & opposes == FALSE]; r[, unit := vapply(set, unit_of, "")]
      u <- r[, .(bound = all(bound), mostly = all(reading == "consists mostly of the listed handling genes")), by = unit]
      arm <- if (d0 == "Up") "up" else "down"; Hd <- if (d0 == "Up") H_up else H_dn; gs <- unique(unlist(lapply(FAMG[fs], intersect, UA)))
      side <- sprintf(paste("Reported beside it: %d of %d flagged %s-arm units bound (a CD4/CD8 Gavish pair counts once); %s of the %d U_A genes in the union of the flagged %s-arm sets",
                            "are in H_%s, against a base rate of %s in U_A."), sum(u$bound), nrow(u), arm, f3(mean(gs %in% Hd)), length(gs), arm, arm, f3(mean(UA %in% Hd)))
      if (nrow(u) && sum(u$bound) >= nrow(u) / 2)
        paste0("Main limitation (7.3): \"The set-level description of the ", arm, " arm rests on genes that also change with 18 h ex-vivo blood processing delay (Savage 2021) or are immediate-early genes; without them it is ",
               if (all(u$mostly[u$bound])) "not testable" else "no longer carried",
               ". These lists cover at least a quarter of the genes tested here (outside the ambient ranking) of 20 of the 22 T-cell programs used (9 of them NAIVE, CYTOTOXIC, T_REG or DYSFUNCTION programs), ",
               "and most naive-identity genes also fall with delay, so a handling reading and a T-cell-state reading of this arm cannot be separated with these data.\" ", side) else
        sprintf("Consequence rule (7.3), %s arm: not triggered. %s", arm, side) })), collapse = " ")
  OPPS <- sprintf("Opposite-sign class (12 genes): %s (p %s); without KLRD1 and NKG7: %s (p %s).", OP[set == "opposite12", reading], fq(OP[set == "opposite12", PValue]),
                  OP[set == "opposite10", reading], fq(OP[set == "opposite10", PValue]))
  ## 7.6 delay lists: T, B and NK against their own lists only
  DLS <- paste(c(CM[family == "delay lists" & view == "primary" & !is.na(reading) & !reading %in% c("", "reported (not read)"), sprintf("%s list: %s", sub("^own_", "delay-", set), reading)],
                 "Monocytes, and the T lists inside B, NK and Mono, are reported, not read (3.3(a))."), collapse = "; ")
  ## 7.4 gene identity, then A1-G01's sentence verbatim and the triggered sentences
  G01S <- paste("In GSE185381 (15 adult AML, 10 healthy donors), the mix of BoneMarrowMap T subtypes accounts for about 5% of the cross-fitted T-imprint gap",
                "(C 0.045, 95% CI \u22120.100 to 0.191; C_H 0.051, \u22120.107 to 0.193; C_A 0.039, \u22120.104 to 0.196), and a shift is detected inside the four main T labels",
                "(S_w 0.664, 95% CI 0.497 to 0.772; a sign check, not a share of the gap). The intervals are descriptive, with the fold lists held fixed.",
                "With gamma-delta and MAIT cells as separate labels the composition share is 0.057 (the flags are lower bounds).")
  t1 <- TG[view == "primary (threshold 1, both references)"]; READ4 <- all(t1$consensus_readable); tt <- t1[arm %in% c("P370_up", "P370_down")]
  cls <- function(cl) gsub("_", "-", cl)
  TAGS <- if (!READ4) "Gene-identity consensus not readable (K4 failed in a reference); each reference is reported alone in the tags table." else
    paste(c(tt[, ifelse(n_consensus > 0,
                         sprintf(paste("%d of the %s genes are %s-high in healthy marrow T cells in both references (GSE185381 healthy donors and Oetjen; paired within donor; base rate %s;",
                                       "matched-random percentile %s). Such genes are expected both from a change within cells and from a %s %s share in AML (composition-only expectation,",
                                       "UMI-weighted at a common CD8 share: median log2 %s, P10 to P90 %s to %s); the tags cannot tell these apart."),
                                 n_consensus, sub("P370_", "", arm), cls(class), f3(base_rate), f3(random_pct), ifelse(arm == "P370_up", "higher", "lower"), cls(class), f3(comp_p50), f3(comp_p10), f3(comp_p90)),
                         sprintf("0 of the %s genes are %s-high in healthy marrow T cells in both references (base rate %s).", sub("P370_", "", arm), cls(class), f3(base_rate)))],
            sprintf("Not assessable in both references: %d up and %d down genes.", tt[arm == "P370_up", not_assessable][1], tt[arm == "P370_down", not_assessable][1])), collapse = " ")
  TAGS <- paste(TAGS, G01S, "C is a share of the score gap, not of gene classes.")
  if (READ4 && tt[class == "naive", any(n_consensus > 0)])   # a naive-like description (7.4)
    TAGS <- paste(TAGS, "The healthy donors are 31.5 years younger in median; in healthy marrow the weak age association of the score runs through the naive-T share (S9).")
  if ((READ4 && t1[arm %in% c("P370_up", "P370_down", "core25") & class == "effector_memory", any(n_consensus > 0)]) || any(grepl("CYTOTOXIC", FL$set)) ||
      any(CR[view == "primary", reading] == "shared with NK", na.rm = TRUE))   # the effector-like trigger (7.4)
    TAGS <- paste(TAGS, "CMV serostatus is unrecorded in both arms; CMV expands TEMRA-like cells inside the CD8 effector-memory labels and adaptive NK cells, and rises with age.")
  ## 7.5 S5b (44 writes the sentence on the primary rows; fwrite doubles inner quotes and fread keeps them)
  S5BS <- paste(unique(gsub('""', '"', S5B[!is.na(cohort_sentence) & nzchar(cohort_sentence), cohort_sentence], fixed = TRUE)), collapse = " ")
  fails <- function(T0, kind) { x <- T0[state == kind, unique(id)]; if (length(x)) paste(x, collapse = ", ") else "none" }
  KA <- sprintf("43: deterministic FAIL %d, statistical FAIL %s. 44: deterministic FAIL %d, statistical FAIL %s.", C[state == "FAIL_deterministic", .N], fails(C, "FAIL_statistical"),
                S5C[state == "FAIL_deterministic", .N], fails(S5C, "FAIL_statistical"))
  ## 7.7 context, verbatim
  CTX <- paste(
    "A1-BATCH (GATES.md:293): \"Allowed wording: none of the six within-study processing factors, tested one at a time, explained or weakened the frozen direction; sample source, storage, sorting method and chemistry version are constant within arms or unrecorded, and were not tested.\"",
    "A1-BATCH 12.5: \"The imprint direction is shared with B and NK cells of the same donors (D10: R 0.69 and 0.93), so it is not T-cell-specific. It overlaps a blood-processing-delay program in the same direction (D8), and adding the delay score to MAIN lowers R to 0.82 (fhat 0.30; descriptive, no reference). Dropping the Control1-5 group lowers R to 0.81 (fhat 0.31; percentile 0.10), a share G1 could not have flagged.\"",
    "A1-G01 (GATES.md:319): \"beyond T-subset composition (detected inside the main labels)\", together with A1-G01 11.5's sentence (gene-identity row) and A1-G01's section-0 limits.",
    "Monocytes (A1-G01 11.4): \"monocytes in AML samples are often clonal, so a monocyte shift may be leukaemic biology rather than a shared environment or handling effect.\"",
    "A1-G04 (GATES.md:346): \"The direction was not shown in Chen2023 (exact p = 17/126); a common shift larger than 1.36 of GSE185381's selection-corrected gap is excluded (one-sided 95.2%, common-shift model).\" \"5 vs 4 cannot exclude an effect of the GSE185381 size.\" The direction stays \"shown in one study\".",
    "S9 (GATES.md:241 and the S9 block): \"below one third at the point estimate; upper bound 0.463 (worded as a bound)\"; \"The S5 label 'AML-or-age' is narrowed to 'AML, with an age share bounded below one third in healthy marrow (one protocol)'\"; \"a weak age association exists (f_slope 0.102 above the random P95 0.046)\".",
    "G07: \"Samples were taken at diagnosis; acute illness, infection, transfusion and medication are unrecorded, and no non-AML patient comparator was run, so 'AML-associated' includes being acutely ill. The commercial-versus-biobank source of the arms is not separated from disease.\"")
  SUM <- data.table(item = c("set-level description (7.1)", "secondary views (7.2)", "consequence rule (7.3)", "opposite-sign class (3.2(6))", "delay lists (3.3(a), 7.6)",
                             "gene identity (3.4, 7.4)", "S5b (3.5, 7.5)", "known answers (section 6)", "context (7.7)", "scripts"),
                    value = c(SETS, SEC, CONS, OPPS, DLS, TAGS, S5BS, KA, CTX, sprintf("scripts/ws/43_a1_g07.R md5 %s; scripts/ws/44_a1_g07_cells.R md5 %s", MD5, MD544)))
  put(SUM, "summary")
  BLOCK <- c("<!-- A1-G07 start -->", sprintf("## What the T axis is: a pre-registered description (A1-G07, prereg-a1-g07-v1) -- %s", format(Sys.Date())), "",
             sprintf("Scripts: scripts/ws/43_a1_g07.R (md5 %s), scripts/ws/44_a1_g07_cells.R (md5 %s). GSE185381, 15 adult AML vs 10 healthy, S5's T fit; P370. A description, not a gate.", MD5, MD544), "",
             "| item | reading |", "|---|---|", SUM[item != "scripts", sprintf("| %s | %s |", item, gsub("\\|", "/", value))], "",
             "**Reading: see the table; every set sentence re-describes the contrast that defined the axis and is not independent evidence (PREREGISTRATION_A1_G07.md section 7).**",
             "<!-- A1-G07 end -->")
  if (SMOKE) writeLines(BLOCK, rg("gates_block.md")) else cat(paste0("\n", paste(BLOCK, collapse = "\n"), "\n"), file = here::here("scripts", "ws", "GATES.md"), append = TRUE)
  ts("[combine] done"); writeLines(BLOCK)
}
