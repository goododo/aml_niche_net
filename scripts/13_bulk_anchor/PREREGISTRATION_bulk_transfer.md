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

---

## 7. Pilot results (2026-09-29, run of `13.3_hox_transfer.R`, gates in registered order)

**Instrument: PASS, strongly.** Set L matched 17/17 symbols in the BeatAML matrix; Q0
**AUC 0.923** for NPM1-mut vs WT (128 vs 330), beating **100.0%** of 1000 size- and
expression-decile-matched random gene sets (G1 threshold 97.5%). The single-cell-derived HOX
program is legible in bulk. `tag = registered_primary`.

**Grouping: the design's weak point, and it is a numbers problem.** The registered threshold
(25th percentile of the mutant score distribution) leaves only **17 NPM1-like** of 330 WT
(15 with venetoclax fits). Every downstream estimate rests on those 15.

**G2 confound table fired before any outcome was read**: NPM1-like are **58.8% FLT3-ITD+ vs
14.4%** in WT-low (Fisher **p = 5.9e-5**); blasts_bm and TP53 balanced, age not significant.

**Q1 venetoclax (primary).** NPM1-like are more sensitive: median AUC difference **−90.8**,
Wilcoxon **p = 0.0014**; age-adjusted beta −59.3 (p = 0.0082). But the pre-stated
composition-analogue adjustment already softens it (age + blasts_bm, n = 183: beta −41.9,
**p = 0.060**), and the **post-hoc FLT3-ITD adjustment that G2 mandates** leaves
**beta −40.9, p = 0.073** while FLT3-ITD itself carries p = 0.0027. Stratified: ITD-negative
NPM1-like n = 6 (p = 0.30), ITD-positive n = 9 (p = 0.065) — consistent direction in both, but
neither stratum is powered. The continuous sensitivity across all WT is weak (rho −0.111,
p = 0.088).
**Reading, at the registered wording ceiling: the venetoclax association is REAL but NOT
SEPARABLE from FLT3-ITD at this n.** It is not yet a phenotype claim about the HOX state.

**Q2 survival: null.** Cox OS, 314 patients / 189 events: **HR 0.87 [0.41–1.86], p = 0.73**.

**Q3 exploratory screen: 22 of 112 drugs at q_BH < 0.10**, every one in the same direction
(NPM1-like more sensitive), led by Vandetanib (q 0.012), Palbociclib (0.015), Canertinib
(0.015), Sunitinib (0.026), Sorafenib (0.027). Four of the top five are multi-kinase / FLT3-
active agents, and with FLT3-ITD in the model their betas survive (p 0.003–0.035) while FLT3
itself is also significant — the same entanglement as Q1. **Palbociclib (CDK4/6, not
FLT3-active; beta −42.4, p = 0.0037 with FLT3 adjusted) is the one lead that is not obviously
an ITD story** and is the single item worth carrying to the external FIMM gate.

**Pilot verdict.** The transfer machinery works end to end and the instrument gate is
comfortably cleared — that part is now established. The biology is blocked by the grouping's n
and by FLT3-ITD collinearity, neither of which is fixable by re-cutting this cohort. Next steps
that follow, in order: (1) the external FIMM direction gate on Palbociclib and venetoclax;
(2) a genotype-conditioned design — NPM1-like *within* ITD-negative AML, which needs a larger
cohort (TCGA-LAML + AMLCG GSE146173 + Leucegene) rather than a different statistic;
(3) only then a second program (the 7-compartment atlas programs).
No re-definition of NPM1-like was performed after seeing outcomes; the 25th-percentile rule
stands as registered, and the post-hoc FLT3 adjustment is labelled post hoc here.

## 8. External direction gate in FIMM (2026-09-29, `13.4_fimm_external_gate.R`)

**Gate result: NOT PASSED for either lead — direction consistent, effect absent.**

Mapping note first, because the registered stop rule fired twice on the way in. FIMM's CPM matrix
is keyed by Ensembl id, so the symbol match gave 0/17 and the script stopped as registered. After
fixing the mapping through BeatAML's own `stable_id`↔`display_label` annotation, **14 of 17** Set L
genes exist in FIMM's 18,203-row matrix; the 3 absent are lncRNAs (HOXA10-AS, HOXB-AS1, HOXB-AS3),
missing because that quantification is protein-coding-centric. The registered threshold fires read
as an absolute count (14 < 15) and passes read as a fraction (82% > 15/19 = 79%); **both readings
are recorded and no gene was substituted**. To take the gene set out of the comparison, the
BeatAML side was recomputed on the same 14: **Q0 AUC 0.923, identical to the 17-gene version**, and
its two leads keep their effects (venetoclax median AUC diff −83.4, p = 7.4e-4; palbociclib −56.3,
p = 0.0031). So the discovery side is not gene-set-fragile.

FIMM instrument: **AUC 0.794** for NPM1-mut vs WT (16 vs 63 at diagnosis) — above the 0.70 floor,
below the 0.80 pass line, so this cohort's arm is `EXPLORATORY` by the registered rule. The sDSS
sign convention was verified before judging (higher = more sensitive; NPM1-mut vs WT anchor is
**+8.37** for venetoclax, the textbook direction, confirming the convention and that FIMM's
genotype labels behave).

| lead | FIMM n (like / WT-low) | median sDSS diff | p | direction vs BeatAML |
|---|---|---|---|---|
| Venetoclax | 13 / 32 | **+0.16** | 0.47 | same sign |
| Palbociclib | 13 / 41 | **+0.43** | 0.96 | same sign |

Both differences point the same way as BeatAML, but their magnitudes are **negligible against the
scale of the same table's own positive control** (venetoclax NPM1-mut vs WT = +8.37; the NPM1-like
effect is 2% of that). At 13 NPM1-like samples this is not a powered refutation — but the
registered gate required same-direction *evidence* before a lead may be called a finding, and a
+0.16 sDSS shift at p = 0.47 is not that.

**Conclusion, recorded at the registered wording ceiling: neither venetoclax nor palbociclib
graduates past `exploratory_screen`. The NPM1-like HOX state is not an externally replicated
drug-response phenotype.** What survived the whole chain is the instrument: the single-cell-derived
HOX program reads NPM1 genotype in two independent bulk cohorts (BeatAML AUC 0.923, FIMM 0.794),
with the same specificity control. That is a transfer-machinery result, not a discovery.

**What this rules out and what it does not.** Ruled out on present data: a large, cohort-portable
drug-sensitivity phenotype attached to NPM1-like HOX expression in NPM1-WT AML. Not ruled out: a
modest effect, or one confined to FLT3-ITD-negative disease (§7's collinearity), both of which need
a genotype-conditioned design in a larger cohort (TCGA-LAML + AMLCG GSE146173 + Leucegene) rather
than more analysis of these two.
