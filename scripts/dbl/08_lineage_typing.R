#!/usr/bin/env Rscript
# 08_lineage_typing.R ----
# Doublet-pairing SMOKE, stage 8 (user decision 2026-11-10: fall back to the lineage level, v1.1's Q1 -- are T-LEUK doublets
# more frequent than random co-encapsulation predicts, state ignored). First question: is LINEAGE-level typing good enough?
# Re-scores the EXISTING synthetic doublets of stages 04 (scDblFinder origins) and 06 (NNLS weights) at the lineage level;
# runs no caller. LEUK (lineage) = the five states + LEUK_other (v1.3: "stem/progenitor and blast-like compartment");
# Mono/DC stays outside, as in the proposal.
#   scdbl_main / scdbl_devel   stage 04 origin {T, one LEUK-lineage compartment}, called, not ambiguous, inside the ambient window
#   nnls_devel                 stage 06 NNLS: two heaviest components T + one LEUK-lineage compartment (06's margins), among the
#                              doublets devel called (stage 06's primary scope); no ambient window, as in 06
# Precision is read at NATURAL abundance, not at the arbitrary injection mix: for each qualifying patient (singlet floors),
# expected reported T-LEUK from signal and from every false-positive source under random co-encapsulation (stage 04's
# expectation model, lineage level): same-donor doublets dbr * N * p_d^2, pair probabilities from the donor's class-A
# composition, rates from the synthetic doublets. Unsimulated partners borrow the nearest simulated type: X|LEUK_other uses
# X|LEUK* (review R-04 F5). Recall is abundance-weighted over the donor's LEUK-lineage compartments.
# MASKING [D2 + decision 2026-11-10]: the singlet rows of stage 04's 2026-11-09 tables are opened ONLY inside
#   singlet_fp_pooled(), which returns per-compartment counts of "typed as lineage T-LEUK" (pooled, the allowed summary) and
#   nothing per barcode. NNLS singlets carry only stage 06's masked five-state flag: used as a lower bound (stated).
# READING RULE (fixed 2026-11-10 before any run): a method passes if, in BOTH lanes and for EVERY qualifying patient, natural
#   recall >= 0.7 and natural precision >= 0.7. Verdict "lineage-level typing is good enough to design the O/E test" if any
#   method passes; otherwise "NK/Tamb confusion has to be solved first".
# INPUT  : DBL_TAB_DIR/dbl_typing_injected_{main,devel}_smoke.csv.gz (04, jobs 3687239/3687240), dbl_nnls_injected_smoke.csv.gz (06),
#          DBL_LANE_DIR/<lane>__calib_main_smoke.rds (04: dbr, N, class-A composition by donor; nothing typed is read from it)
# OUTPUT : DBL_TAB_DIR/dbl_lineage_{rates,natural,verdict}_smoke.csv, dbl_checks_08_smoke.csv
# Usage  : Rscript scripts/dbl/08_lineage_typing.R

suppressPackageStartupMessages({ library(data.table); library(here) })
source(here::here("scripts", "config", "config_paths.R")); source(here::here("scripts", "dbl", "config_dbl.R"))
ts  <- function(...) message(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ...)
CHK <- list()
chk <- function(name, expected, got) {
  ok <- isTRUE(all.equal(expected, got, check.attributes = FALSE))
  CHK[[length(CHK) + 1]] <<- data.table(check = name, expected = paste(expected, collapse = ";"), got = paste(got, collapse = ";"), pass = ok)
  message(sprintf("  %-84s expect %-8s got %-8s %s", name, paste(expected, collapse = ";"), paste(got, collapse = ";"), if (ok) "PASS" else "** FAIL **")) }
STATES <- DBL_LEUK_STATES
LIN    <- c(STATES, "LEUK_other")                         # LEUK at the lineage level
WIN    <- DBL_AMBIENT_WINDOW
nz <- function(z) fifelse(z == "", NA_character_, z)
lin_origin <- function(origin) {   # TRUE when an scDblFinder origin is exactly {T, one LEUK-lineage compartment}
  p <- strsplit(ifelse(is.na(origin), "", origin), "+", fixed = TRUE)
  vapply(p, function(z) length(z) == 2 && "T" %in% z && sum(z %in% LIN) == 1, TRUE) }
scdbl_lin <- function(d) d$called %in% TRUE & d$ambiguous %in% FALSE & lin_origin(d$origin) & !is.na(d$share) & d$share >= WIN[1] & d$share <= WIN[2]
singlet_fp_pooled <- function(path) {   # MASKED: per-compartment pooled counts only; nothing per barcode leaves this function
  x <- fread(path)[startsWith(type, "singlet:")]
  x[, origin := nz(origin)]
  r <- x[, .(n = .N, lin_any = sum(scdbl_lin(.SD)), called_any = sum(called %in% TRUE)), by = .(lane, comp = sub("^singlet:", "", type))]
  rm(x); r }
## -- [1] synthetic pair rows (non-singlet), lineage calls ----
S04 <- rbindlist(lapply(c("main", "devel"), function(m) {
  x <- fread(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_", m, "_smoke.csv.gz")))[!startsWith(type, "singlet:")]
  x[, `:=`(origin = nz(origin), true_state = nz(true_state))]
  x[, .(lane, run, kind, syn, type, true_state, method = paste0("scdbl_", m), lin = scdbl_lin(.SD))] }))
N06 <- fread(file.path(DBL_TAB_DIR, "dbl_nnls_injected_smoke.csv.gz"))[method == "nnls" & kind %in% c("heldout", "lodo") & !startsWith(type, "singlet:")]
N06[, true_state := nz(true_state)]
N06[, lin := called_devel %in% TRUE & ((c1 == "T" & c2 %in% LIN) | (c2 == "T" & c1 %in% LIN)) & pmin(w1, w2) >= 0.15 & (w2 - w3) >= 0.05]
ALL <- rbind(S04, N06[, .(lane, run, kind, syn, type, true_state, method = "nnls_devel", lin)])
## the true lineage partner of each synthetic type: T|LEUK -> its state; T|LEUK_other -> LEUK_other; everything else is a false source
ALL[, truth := fifelse(type == "T|LEUK", true_state, fifelse(type == "T|LEUK_other", "LEUK_other", NA_character_))]
chk("every T|LEUK row carries its true state", 0L, ALL[type == "T|LEUK" & is.na(true_state), .N])
## -- [2] per-type rates (held-out runs; lodo reported) ----
RATES <- ALL[, .(n = .N, rate = mean(lin)), by = .(lane, method, kind, type, truth)]
## singlet false-positive rates: scDblFinder from the 11-09 tables (pooled inside the masked function); NNLS lower bound
SG <- rbindlist(lapply(c("main", "devel"), function(m) singlet_fp_pooled(file.path(DBL_TAB_DIR, paste0("dbl_typing_injected_", m, "_smoke.csv.gz")))[, method := paste0("scdbl_", m)]))
SG[, rate := lin_any / n]
## NNLS singlets: stage 06 decomposed every re-entered singlet without a detection step and kept only its masked five-state flag.
## Detect-then-type rate per compartment = devel's pooled called rate x NNLS's typed rate (independence assumed; five-state flag,
## so T + LEUK_other singlets are missed: a lower bound on that part)
SGn <- fread(file.path(DBL_TAB_DIR, "dbl_nnls_injected_smoke.csv.gz"))[method == "nnls" & kind == "heldout" & startsWith(type, "singlet:")][
  , .(n = .N, typed = sum(tleuk_any %in% TRUE)), by = .(lane, comp = sub("^singlet:", "", type))]
SGn <- merge(SGn, SG[method == "scdbl_devel", .(lane, comp, called_rate = called_any / n)], by = c("lane", "comp"), all.x = TRUE)
SGn[, `:=`(rate = (typed / n) * fifelse(is.na(called_rate), 0, called_rate), lin_any = NA_integer_, called_any = NA_integer_, method = "nnls_devel")]
SG <- rbind(SG, SGn[, .(lane, comp, n, lin_any, called_any, method, rate)], use.names = TRUE)
ts("[2] singlet false positives (lineage, pooled over compartments):"); print(SG[, .(n = sum(n), mean_rate = round(mean(rate), 4)), by = .(lane = sub("_.*", "", lane), method)])
## -- [3] natural abundance per qualifying patient ----
NAT <- list()
for (L in DBL_SMOKE_LANES) {
  cb <- readRDS(file.path(DBL_LANE_DIR, paste0(L, "__calib_main_smoke.rds")))   # dbr, N and composition only
  comp <- cb$comp; dbr <- cb$dbr; N_tot <- cb$N                                  # (comp has its own column N: cell counts)
  qual <- comp[, .(nT = sum(N[compartment == "T"]), nL = sum(N[compartment %in% STATES])), by = donor][nT >= 30 & nL >= 50, donor]
  for (m in unique(ALL$method)) for (d in qual) {
    cd <- comp[donor == d][, f := N / sum(N)]; g <- function(cp) cd[compartment == cp, if (.N) f else 0]
    pd <- comp[, sum(N[donor == d]) / sum(N)]; Nd_dbl <- dbr * N_tot * pd^2; Nd_sing <- N_tot * pd - Nd_dbl
    R  <- RATES[lane == L & method == m & kind == "heldout"]
    rk <- function(k) { r <- if (k %in% STATES) R[type == "T|LEUK" & truth == k, rate] else R[type == "T|LEUK_other", rate]; if (length(r)) r else NA_real_ }
    fL <- sum(vapply(LIN, g, 0))
    sig_k <- vapply(LIN, function(k) 2 * g("T") * g(k) * rk(k), 0)        # NA when a compartment had no synthetic doublets
    ru <- function(u) { r <- R[type == u, rate]; if (length(r)) r else 0 }
    fp_pairs <- 2 * g("T") * g("Tamb") * ru("T|Tamb") + 2 * g("T") * g("NK") * ru("T|NK") + 2 * g("T") * g("B") * ru("T|B") +
                2 * g("T") * g("Ery") * ru("T|Ery") + 2 * g("T") * g("MonoDC") * ru("T|MonoDC") + g("T")^2 * ru("T|T") +
                2 * g("Tamb") * fL * ru("Tamb|LEUK*") + 2 * g("NK") * fL * ru("NK|LEUK*") + 2 * g("B") * fL * ru("B|LEUK*") +
                2 * g("Ery") * fL * ru("Ery|LEUK*") + 2 * g("MonoDC") * fL * ru("MonoDC|LEUK*") + fL^2 * ru("LEUK*|LEUK*")
    sg <- SG[lane == L & method == m]
    fp_sing <- sum(vapply(cd$compartment, function(cp) { r <- sg[comp == cp, rate]; if (length(r)) g(cp) * r else 0 }, 0))
    signal <- Nd_dbl * sum(sig_k, na.rm = TRUE); fp <- Nd_dbl * fp_pairs + Nd_sing * fp_sing
    w <- vapply(LIN, g, 0); ok <- !is.na(sig_k)
    NAT[[length(NAT) + 1]] <- data.table(lane = L, method = m, donor = d, f_T = g("T"), f_LEUK = fL,
      recall_nat = sum(w[ok] * vapply(LIN[ok], rk, 0)) / sum(w[ok]), exp_signal = signal, exp_fp_pairs = Nd_dbl * fp_pairs,
      exp_fp_singlets = Nd_sing * fp_sing, precision_nat = signal / (signal + fp)) }
}
NAT <- rbindlist(NAT)
NAT[, pass := recall_nat >= 0.7 & precision_nat >= 0.7]
## -- [4] verdict ----
byM <- NAT[, .(lanes = uniqueN(lane), all_pass = all(pass)), by = method]
winners <- byM[lanes == length(DBL_SMOKE_LANES) & all_pass, method]
verdict <- if (length(winners)) paste0("lineage-level typing is good enough to design the O/E test (", paste(winners, collapse = ", "), ")") else
  "NK/Tamb confusion has to be solved first"
ts("[3] natural recall / precision per qualifying patient:")
print(NAT[, .(lane = sub("_.*", "", lane), method, donor, recall_nat = round(recall_nat, 2), precision_nat = round(precision_nat, 2),
              exp_signal = round(exp_signal, 1), exp_fp_pairs = round(exp_fp_pairs, 1), exp_fp_singlets = round(exp_fp_singlets, 1), pass)])
ts("[3] main false-positive sources (held-out, lineage rate):")
fpv <- RATES[kind == "heldout" & is.na(truth) & type %in% c("NK|LEUK*", "Tamb|LEUK*", "T|NK", "T|Tamb", "T|Ery", "T|B", "T|MonoDC", "LEUK*|LEUK*")]
fpv[, col := paste(method, sub("_.*", "", lane), sep = ".")]
print(dcast(fpv, type ~ col, value.var = "rate", fun.aggregate = function(z) round(z, 2)))
ts("[3] lineage recall, leave-one-donor-out (pool):"); print(ALL[kind == "lodo" & !is.na(truth), .(n = .N, recall = round(mean(lin), 2)), by = method])
ts("[3] VERDICT by the pre-registered reading rule: ", verdict)
chk("every qualifying patient has a natural precision", 0L, NAT[is.na(precision_nat), .N])
chk("qualifying patients: AML3948 (single lane) + AML0612, AML3762 (pool)", 3L, NAT[, uniqueN(donor)])
fwrite(RATES, file.path(DBL_TAB_DIR, "dbl_lineage_rates_smoke.csv")); fwrite(NAT, file.path(DBL_TAB_DIR, "dbl_lineage_natural_smoke.csv"))
fwrite(data.table(verdict = verdict, rule = "a method passes if natural recall and natural precision are >= 0.7 for every qualifying patient in both lanes"),
       file.path(DBL_TAB_DIR, "dbl_lineage_verdict_smoke.csv"))
CHK <- rbindlist(CHK); fwrite(CHK, file.path(DBL_TAB_DIR, "dbl_checks_08_smoke.csv"))
ts(sprintf("[done] %d checks, %d FAIL", nrow(CHK), sum(!CHK$pass)))
