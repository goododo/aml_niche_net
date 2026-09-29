# Pre-registration: first sc→bulk transfer experiment (HOX program legibility in BeatAML2)

**Status: DRAFT — awaiting user approval. No association in §4–§6 has been computed.**
Written 2026-09-29, after 13.1 (intake verified 20/20) and 13.2 (cohort frozen) and before any
score touches an outcome. The frozen sample lists are `results/tables/13_bulk_anchor/cohort_freeze.csv`
(committed 4b3ba95); this document must be committed before the implementing script exists.

## 0. What was already seen (honesty ledger)

- The sc instrument: LMPP_GMP HOX score separates NPM1-mut at AUC 0.975 (n=36) in single-cell
  pseudobulk (P2, `12_pseudobulk_de/05_hox_axis.R`); held-out GSE116256 arm p = 0.0079; 0.2% of
  1000 matched random gene sets do as well. NPM1→HOX itself is textbook — the *instrument* is
  ours, the association is not.
- From 13.2's freeze log (seen, no outcome attached): 1pp cohort n=458; NPM1+ 128, FLT3-ITD+ 113,
  TP53+ 42; blasts_bm available for 356/458. Nothing else.
- Known prior art to cite, not compete with: Mer 2021 (ITD-like transcriptome in FLT3-WT NPM1-mut);
  the phenocopy/discordance precedents and the DeltaNp73 caution collected in
  `DATA_RECON_2026-09-24.md`; Zeng 2022 owns composition→drug response.

## 1. Questions, in test order

- **Q0 (instrument gate).** Is the HOX program legible in bulk at all: does the Set L score
  separate NPM1-mut from NPM1-WT in the 458?
- **Q1 (discordance, primary).** Among NPM1-WT patients, do "NPM1-like" cases (defined §3)
  differ from NPM1-WT/HOX-low in **venetoclax ex vivo AUC**?
- **Q2 (discordance, survival).** Do NPM1-like cases differ in overall survival?
- **Q3 (exploratory screen).** The remaining 165 single-agent inhibitors, BH-corrected,
  `kind_of_evidence = exploratory_screen`.

## 2. Cohort and data (all frozen upstream)

- Samples: `SET_TRIPLE_1PP` (458 patients, one specimen each; selection rule in 13.2 header).
- Expression: `beataml_waves1to4_norm_exp_dbgap.txt` (22,843 genes, GRCh37 symbols), matched by
  `display_label`. Counts matrix is not used here.
- Genotype: 13.2's binary flags (consensus workbook columns, codings documented there).
- Drugs: `beataml_drug_auc_long.csv`, converged single-agent fits only.
- Survival: `overallSurvival` + `vitalStatus` (Unknown excluded from Q2).

## 3. Frozen definitions

- **Score**: Set L exactly as frozen in `12_pseudobulk_de/PREREGISTRATION_hox_axis.md` (19
  locus-defined genes). Per-sample score = mean percentile rank of the matched genes within that
  sample's 22,843-gene expression ranking. Report n matched; **if <15 of 19 match, STOP** (symbol
  mapping inadequate; fix mapping, do not substitute genes).
- **NPM1-like** (primary): NPM1-WT patient whose score ≥ the **25th percentile of the NPM1-mut
  score distribution**. Sensitivity: (a) top decile of NPM1-WT scores; (b) continuous score, no
  dichotomy. The comparison group is all remaining NPM1-WT.
- **Primary endpoint**: venetoclax AUC (one named drug, no multiplicity). Q3 covers the rest.
- **Models**: Q1 primary = Wilcoxon NPM1-like vs rest-of-WT on venetoclax AUC. Adjusted
  sensitivity = linear model AUC ~ group + ageAtDiagnosis (+ blasts_bm in the 356-patient
  subset — the bulk analogue of the composition rung). Q2 = Cox OS ~ group + ageAtDiagnosis;
  ELN2017 as sensitivity only (it partially encodes genotype, adjusting for it over-controls).

## 4. Decision rules, fixed before running

| gate | rule | on failure |
|---|---|---|
| **G0 instrument** | Q0 AUC ≥ 0.80 (direction: mut high) | 0.70–0.80: report, everything downstream becomes exploratory. < 0.70: **STOP**, record "HOX program not legible in bulk", no Q1–Q3 |
| **G1 specificity** | Set L beats ≥ 97.5% of 1000 size- and expression-decile-matched random gene sets on Q0 AUC | fail → STOP (score is not the program, it is expression load) |
| **G2 confound check** | NPM1-like vs rest-of-WT compared on blasts_bm, age, FLT3-ITD and TP53 rates BEFORE outcomes are read; any imbalance is reported next to every result | not a stop — a mandatory caveat table |
| **Q1/Q2 significance** | two-sided p < 0.05 on the primary definition, direction free | report as is; no re-definition of NPM1-like after seeing results |
| **external gate** | any Q1/Q3 drug hit must show the same direction in FIMM/Malani sDSS (Zenodo 10.5281/zenodo.7274740) before it may be called a finding; Q2 direction checked in TCGA-LAML (151 RNA) | inconsistent → stays `exploratory_screen`, stated |

## 5. Wording ceiling (fixed now)

Best case licenses: *"in NPM1-WT AML, an NPM1-like HOX expression state identifies patients whose
ex vivo drug response / survival resembles NPM1-mut disease"* — a phenotype statement.
It does NOT license: any mechanism; "the program causes sensitivity"; ruling out an unmeasured
lesion producing both (the DeltaNp73 objection is named here, in advance, as the reviewer's
first question).

## 6. Outputs

`results/tables/13_bulk_anchor/`: `hox_bulk_scores.csv` (per patient: score, genotype, group),
`transfer_q0_gate.csv`, `transfer_q1_q2.csv`, `transfer_q3_screen.csv` — every row carrying
`kind_of_evidence`. Implementing script: `13.3_hox_transfer.R`, to be written only after this
document is approved and committed.
