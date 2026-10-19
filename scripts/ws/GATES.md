# GATES.md -- one entry per gate: date, numbers, pass or fail

<!-- W1.1 start -->
## GATE_W1.1 positive control -- 2026-10-12

Script: scripts/ws/20_t3c_positive_control.R | config: scripts/ws/config_ws.R (frozen 2026-10-12)

| row | ligand | sender -> receiver | receptor | arm | n | z | p_perm | fake P95 | recovered | diag rho |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | IFNG | T_NK -> Mono_DC | IFNGR2 | main | 28 | -0.55 | 0.5534 | 3.55 | FALSE | 0.38 |
| 1 | IFNG | T_NK -> Mono_DC | IFNGR2 | sens | 28 | 0.10 | 0.4316 | 3.56 | FALSE | 0.38 |
| 2 | IFNG | T_NK -> LMPP_GMP | IFNGR1 | main | 17 | 1.14 | 0.1588 | 1.71 | FALSE | -0.08 |
| 2 | IFNG | T_NK -> LMPP_GMP | IFNGR1 | sens | 17 | 1.64 | 0.0899 | 1.61 | FALSE | -0.04 |
| 3 | TNF | Mono_DC -> LMPP_GMP | TNFRSF1A | main | 18 | 0.35 | 0.4216 | 2.45 | FALSE | 0.53 |
| 3 | TNF | Mono_DC -> LMPP_GMP | TNFRSF1A | sens | 18 | -1.09 | 0.6743 | 2.99 | FALSE | 0.52 |
| 4 | TNF | B_Plasma -> Mono_DC | TNFRSF1B | main | 28 | 0.93 | 0.3886 | 3.90 | FALSE | -0.35 |
| 4 | TNF | B_Plasma -> Mono_DC | TNFRSF1B | sens | 27 | 0.55 | 0.4565 | 4.98 | FALSE | -0.35 |
| 5 | TGFB1 | T_NK -> Mono_DC | TGFBR2 | main | 16 | -1.37 | 0.8871 | 2.52 | FALSE | 0.51 |
| 5 | TGFB1 | T_NK -> Mono_DC | TGFBR2 | sens | 15 | -1.43 | 0.9451 | 1.77 | FALSE | 0.54 |
| 6 | TGFB1 | T_NK -> LMPP_GMP | TGFBR2 | main | 10 | -0.10 | 0.4885 | 1.08 | FALSE | -0.19 |
| 6 | TGFB1 | T_NK -> LMPP_GMP | TGFBR2 | sens | 10 | 0.15 | 0.4386 | 0.90 | FALSE | -0.19 |

Pairs recovered -- IFNG: FALSE | TNF: FALSE | TGFB1: FALSE | arms agree: TRUE
**Verdict: FAIL** -- anything else -> check the pipeline, do not proceed
<!-- W1.1 end -->

<!-- B0 prereg start -->
## GATE_B0 footprint positive control -- PRE-REGISTERED 2026-10-16 (no GATE_B0 result exists)

Config: scripts/ws/config_ws.R section 10 (frozen 2026-10-16) | script: scripts/ws/21_gate_b0_footprint.R |
handoff: HANDOFF_within_sample_v2.0.md §8.3

- **Scope.** Only ligands with a footprint signature (CytoSig 43 cytokines, PROGENy 14 pathways). The 2,196
  detectable ligand-receptor pairs are not the denominator of anything this gate or WP-B reports.
- **Question.** Across GSE185381 AML diagnosis patients, does T/NK ligand expression predict the matching
  CytoSig footprint in leukaemic progenitor cells (HSC_MPP + LMPP_GMP bins, >= 50 cells; T_NK >= 30 cells),
  beyond a blast proxy, ambient RNA, doublets and age group?
- **Pairs.** IFNG, TNF, TGFB1, all sent by T_NK. TNF sent by Mono_DC is reported and labelled clone-internal.
- **Per pair.** Freedman-Lane permutation (10,000), one-sided p < 0.05, AND the ligand t above the 95th
  percentile of 200 fake ligands (same sender, matched detection and inter-patient IQR, not in CellChatDB,
  not in the top 200 genes of the tested signature; same model, same statistic).
- **Verdict** (main doublet arm): pass = IFNG plus one more pair; conditional = TNF and TGFB1 without IFNG
  (stop, user decides); fail = anything else; one_arm = the union-removal arm contradicts a pass (stop).
- **A pass means the method works, not that WP-B is new.** WP-B then continues as v2.0's explanatory layer;
  it is not promoted to the main line, and whether it deserves months is decided separately.
- **Reported, never scored.** PROGENy readout; union-removal doublet arm; ligand x sender fraction;
  signatures without immediate-early / heat-shock genes; B-cell tier ('leukaemia-specific' if the ligand
  survives adding the B_Plasma footprint, else 'marrow-wide'); GSE289435 and GSE227903 direction; TNF from
  Mono_DC (clone-internal).
- **Prior evidence, disclosed before running.** T3c main-effect diagnostic (HALLMARK sets, n 10-18, no
  covariates), LMPP_GMP receivers: IFNG rho -0.08, TNF from Mono_DC 0.53, TGFB1 -0.19.
<!-- B0 prereg end -->

<!-- B0 result start -->
## GATE_B0 footprint positive control -- RESULT 2026-10-16

Script: scripts/ws/21_gate_b0_footprint.R | config: scripts/ws/config_ws.R section 10 (frozen 2026-10-16, commit 4c54f6f)

| analysis | ligand | arm | n | t | p (Freedman-Lane) | fake P95 | pass |
|---|---|---|---|---|---|---|---|
| primary | IFNG | main | 32 | -0.20 | 0.5707 | 1.32 | FALSE |
| progeny | IFNG | main | 32 | 0.81 | 0.1889 | 2.77 | FALSE |
| ligand_x_fraction | IFNG | main | 32 | -0.33 | 0.6241 | 1.32 | FALSE |
| no_IEG | IFNG | main | 32 | -0.25 | 0.5855 | 1.42 | FALSE |
| B_tier | IFNG | main | 31 | -1.15 | 0.8710 | 2.13 | FALSE |
| primary | TNF | main | 32 | 0.35 | 0.3519 | 1.56 | FALSE |
| progeny | TNF | main | 32 | 1.32 | 0.1046 | 1.40 | FALSE |
| ligand_x_fraction | TNF | main | 32 | -0.18 | 0.5492 | 1.56 | FALSE |
| no_IEG | TNF | main | 32 | 0.24 | 0.3961 | 1.56 | FALSE |
| B_tier | TNF | main | 31 | 0.01 | 0.4897 | 1.07 | FALSE |
| primary | TGFB1 | main | 32 | 3.18 | 0.0023 | 2.76 | TRUE |
| progeny | TGFB1 | main | 32 | 1.93 | 0.0321 | 2.08 | FALSE |
| ligand_x_fraction | TGFB1 | main | 32 | 2.52 | 0.0085 | 2.76 | FALSE |
| no_IEG | TGFB1 | main | 32 | 3.27 | 0.0017 | 2.73 | TRUE |
| B_tier | TGFB1 | main | 31 | 2.01 | 0.0301 | 1.20 | TRUE |
| clone_internal_MonoDC | TNF | main | 30 | -0.21 | 0.5839 | 2.82 | FALSE |
| primary | IFNG | sens | 30 | -0.13 | 0.5469 | 1.43 | FALSE |
| primary | TNF | sens | 30 | 0.32 | 0.3678 | 1.43 | FALSE |
| primary | TGFB1 | sens | 30 | 2.81 | 0.0062 | 2.73 | TRUE |

Primary, main arm -- IFNG: FALSE | TNF: FALSE | TGFB1: TRUE | sens-arm tier: fail
B-cell tier (reported): IFNG marrow-wide; TNF marrow-wide; TGFB1 leukaemia-specific
**Verdict: FAIL** -- anything else -> ligand-receptor and footprint evidence are not used in this data
<!-- B0 result end -->

<!-- B0 audit note -->
Audit note (2026-10-16, after the FAIL; details in reports/2026-10_ccc_audit/ outside the repo): the T_NK
sender bin of T3c and GATE_B0 contained BoneMarrowMap "Early Lymphoid" cells, which in AML samples are
largely blasts (15% of T_NK-bin cells across the 32 GATE_B0 patients; 8 patients >= 25%). Re-fitting with
six corrected sender definitions leaves IFNG (t -0.05 to -0.48) and TNF (t -0.07 to 0.15) null and lowers
TGFB1 from t 3.18 to 1.77-2.00. The verdict stands. config_ws.R section 11 now carries the exclusion for
every later ws analysis.

<!-- A1A2 start -->
## GATE_A1 / GATE_A2 -- 2026-10-17 (prereg-s5-v1)

Scripts: scripts/ws/27_s5_imprint.R, 28_s5_coordination.R | compartments with an imprint: 
T, NK, Mono

| gate | criterion | value | pass |
|---|---|---|---|
| GATE_A1_1 | |partial rho(score, hb)| and |rho(score, gran)| < 0.30 in every compartment | T hb 0.00 gran -0.04; NK hb -0.26 gran 0.28; Mono hb 0.18 gran -0.04 | TRUE |
| GATE_A1_2 | imprint not enriched in the extended ambient ranking (hypergeometric p > 0.05) | T 32 genes p 0.035; NK 8 genes p 0.107; Mono 31 genes p 0.000 | FALSE |
| GATE_A1_3 | same-channel concordance above the random-set P95 in >= 2 of 3 channels (pooled over compartments) | GSM5613751 1/1; GSM5613752 1/1; GSM5613756 2/2 | TRUE |
| GATE_A1_4 | healthy donors in mixed channels shift < 1/3 of the AML-healthy gap | Mono 0.17; NK 0.07; T 0.28 | TRUE |
| GATE_A2 | coordination > 0, bootstrap 95% CI excludes 0, above the random-set P95 | 0.123, CI -0.352 to 0.702, random P95 0.452, perm p 0.263 | FALSE |

**Verdict: GATE_A1 FAIL, GATE_A2 FAIL** -- any FAIL stops the line for discussion (v2.0 5.4).
<!-- A1A2 end -->

<!-- A1b-S7 start -->
## H-A1b, unsupervised test (DIALOGUE, prereg-s7-v1 + section 8) -- 2026-10-19

Scripts: scripts/ws/29_s7_dialogue_build_run.R, 30_s7_dialogue_evaluate.R. Criterion 1 p is shown as DIALOGUE's rank-sum p / exact permutation p
(rank of the observed canonical correlation among 1 + 99 permutations, /100); the verdict uses the exact p (decision 2026-10-19, PREREGISTRATION_S7_dialogue.md section 8).

| run | MCP | compartments | crit 1 p DIALOGUE / exact | crit 2 rho (n) | genes T/NK/B/Mono | crit 3 Jaccard | found exact / prereg |
|---|---|---|---|---|---|---|---|
| R1 | MCP1 | T+NK+B | 0.223 / 0.28 | 0.39 (12) | 222/12/0/0 | 0.26 | FALSE / FALSE |
| R1 | MCP2 | T+NK+B | 0.596 / 0.57 | 0.78 (12) | 228/69/0/0 | 0.42 | FALSE / FALSE |
| R1 | MCP3 | B+Mono | 0.912 / 0.89 | NA (10) | 0/0/0/0 | NA | FALSE / FALSE |
| R2 | MCP1 | T+NK+B | 0.106 / 0.14 | 0.66 (22) | 244/159/0/0 | 0.23 | FALSE / FALSE |
| R2 | MCP2 | T+NK+B | 0.193 / 0.25 | 0.68 (22) | 279/145/0/0 | 0.11 | FALSE / FALSE |
| R2 | MCP3 | T+NK+B+Mono | 0.213 / 0.27 | 0.75 (18) | 149/152/149/207 | 0.17 | FALSE / FALSE |

Known answers: best MCP MCP1: exact p 0.75 (DIALOGUE p 0.807), planted genes among its top 200 per compartment 0.05 (over 2 compartments); 10 of 20 shuffles produced an MCP with exact p <= 0.05 (6 with DIALOGUE p <= 0.05)
**Verdict: NO MCP meets the three criteria: H-A1b CLOSED [criterion 1 = exact permutation p, decision 2026-10-19; under the original reading of DIALOGUE's p: none found]**
<!-- A1b-S7 end -->

Note added 2026-10-19 on the first run of this block (job 3650058, dated 2026-10-17; not a result): that run's verdict
was **withdrawn as computed**. Its criteria 1 and 2 stand as printed, but the two known answers and criterion 3 were not valid runs:
(i) the planted program was added to the expression matrix only, not to the PCs that DIALOGUE step I uses, so
the planted run's empirical p-values are identical to R1's and the control tested nothing; (ii) DIALOGUE1.PMD
calls set.seed(seed1), so stability draws 2-20 and all 20 shuffles were the same draw (identical sample sets
and p matrices on disk). Also recorded: DIALOGUE's emp.p is a rank-sum p of the observed value against 99
permutations, so "p <= 0.05" means top 2 of 100 (exact p <= 0.02); R1 MCP1 and R2 MCP1 sit at ranks 4 and 3
(exact p 0.04 and 0.03), and neither has a gene signature in two compartments. The corrected known answers and
criterion 3 are rerun under section 8 of PREREGISTRATION_S7_dialogue.md, which also records the user's decision of
2026-10-19 to read criterion 1 as the exact permutation p; 30_s7_dialogue_evaluate.R regenerates the block, so a
block dated 2026-10-19 or later above this note is the corrected run.
Addendum, later on 2026-10-19: DIALOGUE 1.0 also misassigns compartments to MCPs (row-order bug in
DIALOGUE.identify.cell.types; PREREGISTRATION_S7_dialogue.md section 8, item 5). With it fixed, R1 and R2 MCP1-2 span
T+NK+B rather than T+B and their NK-B pair fails criterion 1 under both readings. R1, R2 and the planted run are rebuilt
with the patched function before the corrected evaluation.

Reading of the 2026-10-19 block (added after job 3650855 completed): the verdict line is the script's fixed wording
for "no MCP passes all three criteria". Both pre-registered known answers FAILED in this run, so the block does not
close H-A1b. (i) The planted program (100 genes, 0.5 SD, all four compartments) was not recovered: step I linked it
in T+NK+B at exact p 0.04 but step II returned no genes for it (5% of planted genes in the best MCP's lists); the
audit in reports/2026-10_s7_dialogue_audit shows step I converging to noise-level optima at n = 13 (its own
permutation null reaches canonical correlations of ~0.86 while the planted direction offers 0.97). (ii) 10 of 20
sample-relabelled shuffles produced an MCP with exact p <= 0.05 (6 of 20 with DIALOGUE's p) against the
pre-registered ceiling of 1 of 20, so criterion-1 values of 0.03-0.04 on the real data carry no weight.
**Recorded outcome: S7 cannot test H-A1b in GSE185381 adults with DIALOGUE 1.0 as configured (30 PCs, k = 3,
13 patients present in all four compartments). H-A1b stays open; it is not a must-have (handoff v2.0 §16.1 item 15).**
Descriptive facts kept for the record: R1 and R2 each return two T+NK+B programs with T and NK gene signatures and
none in B; their sample scores are moderately reproducible across 70% patient subsets (match correlation 0.65-0.78)
but their top-50 gene lists are not (median Jaccard 0.07-0.43, below 0.5 in every compartment).
