#!/usr/bin/env Rscript
# input : gene NAMES only of ws_s5_de_{T,B,NK}.csv (fread select = "gene"); the ambient table; HALLMARK v2024.1; the
#         repo MSigDB v2026.1 GMT (Gavish sets); the Savage set; WS_B0_IEG (parsed from config_ws.R);
#         design covariates only of ws_s5_scores.csv (select = sample, status, hb, gran, sex, subtype_share, compartment;
#         the score column is never read)
# output: /FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L1/03_membership_and_design.log
# what it does: (A) set sizes after removing H from U_A (draft 3.2, K3) for HALLMARK and Gavish, and in the T x X common
#   universes (3.3(b)); (B) the stacked T x X design of 3.3(b) built from the REAL sample-level covariates: rank,
#   residual df, conditioning, and the equivalent paired-difference regression. No expression value is read.
suppressPackageStartupMessages({ library(data.table) })
OUT <- "/FAST/gr10634/gaozy/reports/2026-11_a1_g07_prereg/review1/L1"
TAB <- "/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws"
sink(file.path(OUT, "03_membership_and_design.log"), split = TRUE)

## ---- (A) memberships ----------------------------------------------------------------------------------------
UT <- fread(file.path(TAB, "ws_s5_de_T.csv"), select = "gene")$gene
UB <- fread(file.path(TAB, "ws_s5_de_B.csv"), select = "gene")$gene
UN <- fread(file.path(TAB, "ws_s5_de_NK.csv"), select = "gene")$gene
BL <- fread(file.path(TAB, "decontx/GSE185381/ws_s1_ambient_blacklist.csv"))
EXT <- BL[n_units >= 2, gene]
UA <- setdiff(UT, EXT)
cfg <- readLines("/FAST/gr10634/gaozy/aml_niche_net/scripts/ws/config_ws.R")
i0 <- grep("^WS_B0_IEG <- c\\(", cfg); i1 <- i0 + 2L                        # the three-line definition (lines 160-162)
e_ieg <- new.env(); eval(parse(text = paste(cfg[i0:i1], collapse = " ")), envir = e_ieg); IEG <- e_ieg$WS_B0_IEG  # evaluates only the 3-line definition
stopifnot(length(IEG) == 26L)
SV <- fread("/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv", select = "gene")$gene
H <- intersect(union(SV, IEG), UA)
gmt <- function(p, pat = NULL) { L <- strsplit(readLines(p), "\t"); s <- setNames(lapply(L, `[`, -(1:2)), vapply(L, `[`, "", 1))
  if (!is.null(pat)) s <- s[grepl(pat, names(s))]; s }
HM <- gmt("/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt")
GV <- gmt("/FAST/gr10634/gaozy/aml_niche_net/LCC_proj/scripts/panels/msigdb.v2026.1.Hs.symbols.gmt", "^GAVISH_3CA_METAPROGRAM_CD[48]_T_CELLS_")
cat(sprintf("U_T %d | ext in U_T %d | U_A %d | H in U_A %d | U_A minus H %d | HALLMARK %d sets | Gavish %d sets\n",
            length(UT), length(intersect(UT, EXT)), length(UA), length(H), length(setdiff(UA, H)), length(HM), length(GV)))
UAH <- setdiff(UA, H)
UTB <- setdiff(intersect(UT, UB), EXT); UTN <- setdiff(intersect(UT, UN), EXT)
cat(sprintf("common universes minus ext: T&B %d, T&NK %d\n", length(UTB), length(UTN)))
sz <- function(S, fam) rbindlist(lapply(names(S), function(k) { s <- S[[k]]
  data.table(family = fam, set = sub("^(HALLMARK|GAVISH_3CA_METAPROGRAM)_", "", k), inUA = sum(s %in% UA), H_in = sum(s %in% H),
             afterH = sum(s %in% UAH), inTB = sum(s %in% UTB), inTN = sum(s %in% UTN)) }))
M <- rbind(sz(HM, "P-H"), sz(GV, "P-T"))
M[, fracH := round(H_in / inUA, 2)]
M[, adjacent := fracH >= 0.25]
cat("\n== P-T (Gavish) sets: U_A size, H members, size after H removal, T&B / T&NK sizes\n")
print(M[family == "P-T"][order(afterH)], row.names = FALSE)
cat("\n== counts of sets with >= 15 genes\n")
print(M[inUA >= 15, .(n_tested = .N, keep15_afterH = sum(afterH >= 15), lost_afterH = sum(afterH < 15),
                      lost_and_not_adjacent = sum(afterH < 15 & !adjacent),
                      TB15 = sum(inTB >= 15), TN15 = sum(inTN >= 15)), by = family])
cat("\nP-H sets with the smallest size after H removal:\n"); print(M[family == "P-H" & inUA >= 15][order(afterH)][1:6], row.names = FALSE)
fwrite(M, file.path(OUT, "03_set_sizes.csv"))

## ---- (B) stacked design from the real covariates ----------------------------------------------------------------
D <- fread(file.path(TAB, "ws_s5_scores.csv"), select = c("sample", "status", "hb", "gran", "sex", "subtype_share", "compartment"))
cat("\nrows per compartment:\n"); print(D[, .(n = .N, AML = sum(status == "AML"), healthy = sum(status == "healthy")), by = compartment])
cat("hb, gran identical across compartment rows of a sample:",
    D[, .(u = uniqueN(hb) == 1 && uniqueN(gran) == 1 && uniqueN(sex) == 1), by = sample][, all(u)], "\n")
for (X in c("B", "NK")) {
  T_ <- D[compartment == "T"]; X_ <- D[compartment == X]
  dn <- intersect(T_$sample, X_$sample)
  cat(sprintf("\n==== X = %s: T donors %d, X donors %d, paired %d; T-only donors: %s\n", X, nrow(T_), nrow(X_), length(dn),
              paste(setdiff(T_$sample, X_$sample), collapse = " ")))
  P <- merge(T_[sample %in% dn, .(sample, status, hb, gran, sex, CD8share = subtype_share)],
             X_[sample %in% dn, .(sample, Xshare = subtype_share)], by = "sample")
  cat(sprintf("paired: %d AML / %d healthy; sex male %d of %d (AML %d, healthy %d)\n", P[status == "AML", .N], P[status == "healthy", .N],
              P[sex == "male", .N], nrow(P), P[status == "AML" & sex == "male", .N], P[status == "healthy" & sex == "male", .N]))
  ## stacked rows: T then X
  S <- rbind(P[, .(sample, cT = 1, AML = as.integer(status == "AML"), hb, gran, sexmale = as.integer(sex == "male"), CD8share, Xshare)],
             P[, .(sample, cT = 0, AML = as.integer(status == "AML"), hb, gran, sexmale = as.integer(sex == "male"), CD8share, Xshare)])
  S[, donor := factor(sample)]
  Xs <- cbind(model.matrix(~ 0 + donor, S), cT = S$cT, `cT:AML` = S$cT * S$AML, `cT:hb` = S$cT * S$hb, `cT:gran` = S$cT * S$gran,
              `cT:sexmale` = S$cT * S$sexmale, `cT:CD8share` = S$cT * S$CD8share, `cT:Xshare` = S$cT * S$Xshare)
  cat(sprintf("stacked design: %d rows x %d columns; rank %d; residual df %d\n", nrow(Xs), ncol(Xs), qr(Xs)$rank, nrow(Xs) - qr(Xs)$rank))
  ## equivalent paired-difference regression (2 rows per donor): T - X on [1, AML, hb, gran, sexmale, CD8share, Xshare]
  Dm <- model.matrix(~ AML + hb + gran + sexmale + CD8share + Xshare, P[, .(AML = as.integer(status == "AML"), hb, gran,
                     sexmale = as.integer(sex == "male"), CD8share, Xshare)])
  cat(sprintf("difference regression: %d donors x %d columns, rank %d, residual df %d\n", nrow(Dm), ncol(Dm), qr(Dm)$rank, nrow(Dm) - qr(Dm)$rank))
  Zs <- scale(Dm[, -1]); cat(sprintf("condition number of the centred-scaled covariates: %.1f\n", kappa(cbind(1, Zs), exact = TRUE)))
  vif <- sapply(seq_len(ncol(Zs)), function(j) 1 / (1 - summary(lm(Zs[, j] ~ Zs[, -j]))$r.squared)); names(vif) <- colnames(Zs)
  cat("VIF:", paste(sprintf("%s %.2f", names(vif), vif), collapse = " | "), "\n")
  cat("cor(AML, covariates):", paste(sprintf("%s %+.2f", colnames(Zs)[-1], cor(Zs[, 1], Zs[, -1])), collapse = " | "), "\n")
  ## how well can AML be told apart from the other covariates (R^2 of AML on the rest): the share of AML information left
  r2 <- summary(lm(Zs[, 1] ~ Zs[, -1]))$r.squared
  cat(sprintf("R^2 of AML on the other 5 covariates: %.2f (information on cT:AML left: %.0f %%)\n", r2, 100 * (1 - r2)))
  ## same quantities for the single-compartment S5 designs (25 / X rows)
  for (cc in c("T", X)) { Q <- D[compartment == cc]
    Qm <- model.matrix(~ status + hb + gran + sex + subtype_share, transform(Q, status = factor(status, c("healthy", "AML")), sex = factor(sex)))
    q2 <- summary(lm(Qm[, 2] ~ Qm[, -(1:2)]))$r.squared
    cat(sprintf("  own S5 design %s: %d x %d, rank %d, residual df %d, R^2(AML ~ rest) %.2f\n", cc, nrow(Qm), ncol(Qm), qr(Qm)$rank, nrow(Qm) - qr(Qm)$rank, q2)) }
  fwrite(P[, .(sample, status, hb, gran, sex, CD8share, Xshare)], file.path(OUT, sprintf("03_paired_design_%s.csv", X)))
}
sink()
