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
Tool checks, 2026-10-20 (planted data only; prereg section 8, item 7): with 10 PCs per compartment DIALOGUE recovers
the planted program at n = 13 (recovery 0.88-0.94, largest exact p 0.05) and n = 23 (0.96-1.00, 0.01); with 30 PCs it
fails at n = 23 as it did at n = 13. The "untestable" outcome above therefore applies to the 30-PC configuration; a v2
protocol with 10 PCs is pending approval, with the shuffled-sample specificity gate to be passed before any real result is read.

v2 protocol (PREREGISTRATION_S7_dialogue.md section 9, tag prereg-s7-v2), 2026-10-20: the specificity gate FAILED before
any real-data run. With 10 PCs, 10 of 20 sample-relabelled shuffles of the R1 input produced an MCP with largest-pair
exact p <= 0.05 (5 of 20 with DIALOGUE's p; ceiling 1 of 20), the same rate as with 30 PCs; the R2 gate was at 6 of 11
when this note was written (final count below once complete). No A1b-S7_v2 block is generated. **Final outcome for S7:
untestable in this cohort (10 PCs fix sensitivity, not specificity). H-A1b stays open and is not a must-have.**
R2 final count (2026-10-20, job 3651199): 9 of 20 relabelled shuffles with largest-pair exact p <= 0.05 (7 of 20 with
DIALOGUE's p). Both v2 gates fail; S7 closed as recorded above.

<!-- A1-S8 start -->
## T imprint, structural replication (S8, prereg-s8-v1) -- 2026-10-23

Scripts: scripts/ws/32_s8_replication.R. Gene sets: primary = 370 T-imprint genes without the 32 ambient-prone ones; full = 402. Sets: disc = GSE185381
adult AML (positive control), rep = GSE289435 (decisive), rep2 = GSE227903 diagnosis, desc = GSE116256 (Seq-Well). Random-set nulls as in S5.

| set | gene set | test | statistic | random P5 | random P95 | n | pass | detail |
|---|---|---|---|---|---|---|---|---|
| disc | primary | T2_opposite_movement | -0.29 | 0.41 | 0.78 | 15 | TRUE | coverage 1.00; valid draws 1000 |
| disc | primary | T3_subtype_concordance | 0.78 | 0.22 | 0.84 | 14 | FALSE | patients with CD4 and CD8 >= 20: 14; share of random sets >= observed 0.12; null median 0.59 |
| disc | primary | T1_within_patient_coherence | 0.87 | NA | NA | 15 | TRUE | 13 of 15 patients above their random P95 (median obs 0.021, median P95 0.003) |
| disc | full | T2_opposite_movement | -0.26 | 0.44 | 0.79 | 15 | TRUE | coverage 1.00; valid draws 1000 |
| disc | full | T3_subtype_concordance | 0.74 | 0.23 | 0.83 | 14 | FALSE | patients with CD4 and CD8 >= 20: 14; share of random sets >= observed 0.20; null median 0.58 |
| disc | full | T1_within_patient_coherence | 0.87 | NA | NA | 15 | TRUE | 13 of 15 patients above their random P95 (median obs 0.025, median P95 0.003) |
| rep | primary | T2_opposite_movement | -0.55 | 0.38 | 0.88 | 8 | TRUE | coverage 0.99; valid draws 1000 |
| rep | primary | T3_subtype_concordance | 1.00 | 0.03 | 0.94 | 6 | TRUE | patients with CD4 and CD8 >= 20: 6; share of random sets >= observed 0.03; null median 0.71 |
| rep | primary | T1_within_patient_coherence | 1.00 | NA | NA | 8 | TRUE | 8 of 8 patients above their random P95 (median obs 0.021, median P95 0.004) |
| rep | full | T2_opposite_movement | -0.57 | 0.38 | 0.88 | 8 | TRUE | coverage 0.99; valid draws 1000 |
| rep | full | T3_subtype_concordance | 1.00 | 0.09 | 0.94 | 6 | TRUE | patients with CD4 and CD8 >= 20: 6; share of random sets >= observed 0.04; null median 0.77 |
| rep | full | T1_within_patient_coherence | 1.00 | NA | NA | 8 | TRUE | 8 of 8 patients above their random P95 (median obs 0.024, median P95 0.004) |
| rep2 | primary | T2_opposite_movement | -0.43 | 0.21 | 0.89 | 7 | TRUE | coverage 1.00; valid draws 1000 |
| rep2 | primary | T3_subtype_concordance | NA | NA | NA | 4 | NA | not testable (< 6 patients with both subtypes) |
| rep2 | primary | T1_within_patient_coherence | 1.00 | NA | NA | 7 | TRUE | 7 of 7 patients above their random P95 (median obs 0.024, median P95 0.004) |
| rep2 | full | T2_opposite_movement | -0.71 | 0.32 | 0.89 | 7 | TRUE | coverage 1.00; valid draws 1000 |
| rep2 | full | T3_subtype_concordance | NA | NA | NA | 4 | NA | not testable (< 6 patients with both subtypes) |
| rep2 | full | T1_within_patient_coherence | 1.00 | NA | NA | 7 | TRUE | 7 of 7 patients above their random P95 (median obs 0.029, median P95 0.003) |
| desc | primary | T2_opposite_movement | 0.60 | 0.54 | 0.94 | 6 | FALSE | coverage 0.84; valid draws 1000 |
| desc | primary | T3_subtype_concordance | NA | NA | NA | 2 | NA | not testable (< 6 patients with both subtypes) |
| desc | primary | T1_within_patient_coherence | 0.60 | NA | NA | 5 | FALSE | 3 of 5 patients above their random P95 (median obs 0.012, median P95 0.012) |
| desc | full | T2_opposite_movement | 0.60 | 0.54 | 0.94 | 6 | FALSE | coverage 0.83; valid draws 1000 |
| desc | full | T3_subtype_concordance | NA | NA | NA | 2 | NA | not testable (< 6 patients with both subtypes) |
| desc | full | T1_within_patient_coherence | 0.50 | NA | NA | 6 | FALSE | 3 of 6 patients above their random P95 (median obs 0.011, median P95 0.011) |

Score vs ambient shares: disc Spearman score~hb 0.07, score~gran -0.35; rep Spearman score~hb 0.14, score~gran -0.24; rep2 Spearman score~hb 0.43, score~gran 0.04; desc Spearman score~hb -0.68, score~gran 0.23
Known answers: 0 of 100 scrambles pass T2 on GSE289435 (rho < 0 and < random P5 0.38); GSE185381 adult AML: T2 TRUE, T3 FALSE, T1 TRUE
**Verdict: GATE_S8 PASS: the T imprint's structure replicates in GSE289435**
<!-- A1-S8 end -->

Reading of the A1-S8 block (user decision 2026-10-23, PREREGISTRATION_S8_replication.md section 9): T3 is uninformative
in this design (random gene sets are 0.6-0.8 concordant between CD4 and CD8 of the same patient; the positive control fails
it; the replication passes it at the ceiling). GATE_S8 is read as T2 AND T1: **PASS** in GSE289435, consistent in GSE227903,
not reproduced on Seq-Well GSE116256. Direction versus healthy donors and age remain untested in independent data.

<!-- A1-S9 start -->
## T imprint, age/direction control on older healthy donors (S9, prereg-s9-v1) -- 2026-10-26

Scripts: scripts/ws/33_s9_age_direction.R. Primary set = 370 genes. Decisive cohort GSE120221 (Oetjen, 20 of 20 donors eligible, 24-84 y); GSE227690 and GSE276360 descriptive.
Yardstick: gap_in 1.9714 (3.84 healthy-SD; registered 1.970 = 3.87 healthy-SD, computed on 10-23 without the blacklist step; K6 |diff| 0.0014), cross-fitted gap_oof 1.378, r 1.43; D_res_J (Oetjen) 1.968; gate denominator D_res_J / r 1.376. Delta_age 31.5 y.

| test | cohort | gene set | statistic | threshold | random P5 | random P95 | n | pass | detail |
|---|---|---|---|---|---|---|---|---|---|
| A1_f_slope | GSE120221 | primary | 0.102 | component of f_age | -0.032 | 0.046 | 20 | NA | beta_age x 31.5 / (D_res_J/r); sex-adjusted; unadjusted 0.164; rank among random sets 1.000 |
| A1_f_contrast | GSE120221 | primary | 0.073 | component of f_age | -0.018 | 0.065 | 15 | NA | old (n 5, mean age 67.4) minus young (n 10, mean age 39.0) |
| A1_f_age | GSE120221 | primary | 0.102 | < 1/3 (point estimate) | NA | NA | 20 | TRUE | 95% donor-bootstrap CI -0.085 to 0.463 (conditional on the yardstick); r 1.430; LOO range 0.045 to 0.185; LOO flip FALSE |
| A1_f_age_full402 | GSE120221 | full | 0.121 | sensitivity | NA | NA | 20 | NA | f_slope 0.121, f_contrast 0.097 over the 402 set's own D_res_J/r 1.373 (with the primary set's denominator: f_age 0.120) |
| A2b_alignment_b_over_lambda | GSE120221 | primary | 0.050 | below random P95 -> wording | -0.050 | 0.206 | 358 | TRUE | b 0.046, lambda 0.92 (random median lambda 0.51, median b/lambda 0.074 = genome-wide rate); Spearman 0.167; sign concordance 0.601 |
| A2a_removal_ratio | GSE120221 | primary | 1.007 | reported | 0.992 | 1.008 | 29 | NA | 29 of 370 removed at raw p < 0.05 (0 at FDR < 0.10); removed genes' age sign = AML sign in 0.62; random-removal median ratio 1.000 |
| A3_direction | GSE120221 | primary | 0.073 | A_E/(D_res/r) < 1/3 and <= 1 old donor above the midpoint | NA | NA | 15 | TRUE | S_E 0.275 (\|S_E\| < D_res/2 = 0.984: TRUE), sd ratio 0.586, interpretable TRUE, descriptive (<= 3 young) FALSE; D_dir 1.593; old above midpoint 0 of 5 |
| A3_direction | GSE227690 | primary | -0.022 | A_E/(D_res/r) < 1/3 and <= 1 old donor above the midpoint | NA | NA | 5 | NA | S_E -0.120 (\|S_E\| < D_res/2 = 0.985: TRUE), sd ratio 0.734, interpretable TRUE, descriptive (<= 3 young) TRUE; D_dir 2.120; old above midpoint 0 of 2 |
| A3_direction | GSE276360 | primary | 0.296 | A_E/(D_res/r) < 1/3 and <= 1 old donor above the midpoint | NA | NA | 4 | NA | S_E 0.838 (\|S_E\| < D_res/2 = 0.981: TRUE), sd ratio 1.526, interpretable TRUE, descriptive (<= 3 young) TRUE; D_dir 0.718; old above midpoint 0 of 2 |
| W_rho_score_age | GSE120221 | primary | 0.277 | reported | -0.195 | 0.366 | 20 | NA | two-sided p 0.237; Fisher-z 95% CI -0.189 to 0.641; slope 0.122 W-units per decade; young-donor SD 0.568; \|rho\| rank among random sets 0.847 |
| A1_f_slope_naive_share_adjusted | GSE120221 | primary | -0.047 | sensitivity | NA | NA | 20 | NA | lm(score ~ age + sex + naive share among kept T cells), J space |
| W_rho_subtype_CD4Naive | GSE120221 | primary | 0.321 | sensitivity | NA | NA | 19 | NA | CD4 Naive pseudobulk, donors with >= 20 kept cells (computed only with >= 6 such donors, an unregistered minimum); two-sided p 0.180 |
| W_rho_subtype_CD8EM | GSE120221 | primary | 0.043 | sensitivity | NA | NA | 18 | NA | CD8 EM pseudobulk, donors with >= 20 kept cells (computed only with >= 6 such donors, an unregistered minimum); two-sided p 0.864 |
| info_score_vs_ambient_partial | GSE120221 | primary | NA | informational | NA | NA | 20 | NA | partial Spearman of the J score with hb_share_nonery -0.571 and gran_share_lymph -0.020 given age |
| info_GSE185381_internal_age | GSE185381 | primary | 0.903 | informational | NA | NA | 25 | NA | sign concordance of the 370 genes between 4 AML < 60 and 5 healthy 39-53: 0.903; rho(score, age) within AML -0.074 (p 0.794), within healthy 0.261 (p 0.466); the four AML < 60 score 1.07/0.48/0.83/0.53 vs healthy range -2.27 to -0.43 |

Known answers: K7 libraries whose meta_<GSE> row equals the donor table (age, sex, donor) = 42/42 (PASS); K7 Oetjen donors whose age and sex equal the BoneMarrowMap reference metadata = 20/20 (PASS); K3 donors whose every library's sex call equals the donor table (29 donors) = 29/29 (PASS); K0 eligible Oetjen donors (>= 20 kept T cells after pooling) / of which >= 60 y = 20/5 (PASS); K4 GSE120221: primary genes in U_J = 358/370 (0.968) (PASS); K4 GSE227690: primary genes in U_J = 329/370 (0.889) (PASS); K4 GSE276360: primary genes in U_J = 355/370 (0.959) (PASS); 3.5 gap_oof > 0 and r <= 2 (yardstick cross-validates and is resolved) = gap_oof 1.378, r 1.430 (PASS); K6 S9 scoring vs S8 (15 AML, z over the 15): Spearman / max |diff| = 1.0000 / 0.00000 (PASS); K6 yardstick on the 25 reproduces the registered 1.970 (|gap_in - 1.970|) = 0.0014 (PASS); K2 technical replicates C1-C2 and S1-S2: |score difference| in J-Oetjen = 0.053 / 0.033 (PASS); K5 full plant: slope recovery (f_slope planted - unplanted, undeflated) = 0.978 (PASS); K5 full plant: contrast recovery vs (mean old age - mean young age)/31.5 = 0.885 vs 0.902 (PASS); K5 third plant: slope recovery = 0.327 (PASS); K5 third plant: contrast recovery vs expected = 0.299 vs 0.301 (PASS); K5 planted A2b: b/lambda of the imprint vs the random-set P95 = 0.976 vs 0.034 (PASS); K5 100 sign-scrambled imprints (permuted ages): f_slope P5 / P95, denominator D_res_J/r (calibration report) = -0.036 / 0.030 (reported); K1 Oetjen: Spearman of the CD8 Naive share of projected T cells with age (one-sided, rho < 0) = rho -0.767, p 0.0000 (PASS); K1 Oetjen: the same for the CD4 Naive share (reported) = rho -0.245, p 0.1493 (reported); K1 GSE227690: both old donors' CD8 Naive share below every young donor's (descriptive) = BM1(31 y) 0.234; BM2(37 y) 0.101; BM3(45 y) 0.060; BM4(75 y) 0.031; BM5(74 y) 0.039 (reported); K1 GSE276360: both old donors' CD8 Naive share below every young donor's (descriptive) = HD2(61 y) 0.048; HD3(70 y) 0.047; YHD636(24 y) 0.089; YHD658(23 y) 0.149 (reported); K1 GSE227690 author labels vs ours on projected T cells: CD8Naive agreement (of author-labelled) / CD4Naive = 0.791 / 0.910 (reported); K1b Oetjen T pseudobulk: mean age coefficient (log2 per year) over naive genes / cytotoxic genes = -0.0131 / 0.0212 (reported)

**Verdict: GATE_S9 PASS: In one healthy marrow cohort (Oetjen 2018: 20 of 20 donors, 24-84 y, 10x 3' v2, one protocol, part of the BoneMarrowMap reference), the T imprint's age shift over the GSE185381 age gap is 0.102 (95% CI -0.085 to 0.463) of the cross-fitted AML-healthy gap, below one third at the point estimate; upper bound 0.463 (worded as a bound). The age alternative to the direction is bounded, not excluded. Direction versus healthy remains shown in the defining cohort only (no independent AML + healthy cohort exists in this work); the 32 ambient-prone genes stay excluded from every statement. The S5 label 'AML-or-age' is narrowed to 'AML, with an age share bounded below one third in healthy marrow (one protocol)'.**
- f_age 0.102 (f_slope 0.102, f_contrast 0.073; 95% donor-bootstrap CI -0.085 to 0.463, conditional on the yardstick); r 1.43; LOO range 0.045 to 0.185
- A2b: imprint genes' age effects not more aligned with their AML effects than the random P95 (b/lambda 0.050, P95 0.206)
- a weak age association exists (f_slope 0.102 above the random P95 0.046)
- older healthy donors in GSE120221 do not score like AML (A_E/(D_res/r) 0.073, n young 10, n old 5)
- GSE227690 is descriptive (n young 3, n old 2): A_E/(D_res/r) -0.022, 0 of 2 older donors above the AML-healthy midpoint after study correction
- GSE276360 is descriptive (n young 2, n old 2): A_E/(D_res/r) 0.296, 0 of 2 older donors above the AML-healthy midpoint after study correction
<!-- A1-S9 end -->

<!-- A1-G02 start -->
## T imprint, out-of-sample statistics (A1-G02 reporting amendment) -- 2026-10-28

Scripts: scripts/ws/35_a1_g02_replay.R (R1, R2, D1), scripts/ws/36_a1_g02_perm_crossfit.R (D2). Declaration: PREREGISTRATION_A1_G02_reporting.md (committed before any run; section 7 before any D1/D2 result was read). GSE185381, 15 adult AML vs 10 healthy, S5 design with metadata sex.

| test | statistic | threshold | n | pass | detail |
|---|---|---|---|---|---|
| R1_monte_carlo_p_402 | 0.025 | reported | 200 relabellings | NA | 4 of 200 reach >= 402 genes (563, 707, 1,166, 1,615); median 0, P95 11.1 |
| R2_pooled_LOO_AUC | 0.927 | reported, conditional on the folds | 25 | NA | bootstrap 0.80-1.00; 25 genes selected in all 25 folds, 299 of the 370 in >= 13 |
| D1_AUC_LPO | 0.927 | the AUC a manuscript quotes | 150 pairs | NA | no empty fold imprint; misordered pairs involve AML1371, AML2451, AML2123 |
| D2_p_gap_oof | 0.008 | < 0.05, else stop for discussion | 1,000 relabellings | TRUE | observed 1.378; null median 0, P95 0.660, max 1.622; no failed folds |
| D2_p_loo_auc | 0.004 | secondary | 1,000 relabellings | NA | |

Known answers: 8/8 PASS (402 genes reproduced; 200 hit counts median 0 / P95 11.1 / max 1,615; 25 out-of-fold scores within 0.001; gap_oof 1.378).

**Verdict: A1-G02 PASS (stop rule not triggered). The T imprint exists out of sample in GSE185381: cross-fitted gap 1.378 (1.98 out-of-fold healthy-SD), full-pipeline permutation p 0.008, leave-pair-out AUC 0.93. Batch, sample source and T-subset composition are not addressed here; the next step is the batch amendment.**
<!-- A1-G02 end -->

<!-- A1-BATCH start -->
## T imprint, processing factors inside GSE185381 (A1-BATCH, prereg-a1-batch-v1) -- 2026-10-31

Scripts: scripts/ws/37_a1_batch.R, scripts/ws/38_a1_batch_cells.R. 15 adult AML vs 10 healthy, S5 design, frozen P370 direction.

| readout | value | reference | n | percentile | verdict | note |
|---|---|---|---|---|---|---|
| G1 drop Control1-5: R | 0.812 | 40 design-matched drops | 40 | 0.100 | no effect | fhat 0.31 (0.26-0.38); z-percentile over 252: 0.135 |
| G2 block-respecting p of the cross-fitted gap | 0.017 | block-respecting relabellings | 1000 | 0.017 | no effect | null median 0.000, P95 1.092; free p (A1-G02) 0.008 |
| W1 period: R | 1.051 | 1,000 residual-permutation draws | 1000 | 0.772 | no effect | fhat -0.08; covariate cannot explain (t cap below 2.78) |
| W2 depth: R | 1.039 | 1,000 residual-permutation draws | 1000 | 0.761 | no effect | fhat -0.06; covariate cannot explain (t cap below 2.78) |
| W3 mito: R | 1.052 | 1,000 residual-permutation draws | 1000 | 0.900 | no effect | fhat -0.09; covariate cannot explain (t cap below 2.78) |
| W4 t_dbl: R | 0.997 | 1,000 residual-permutation draws | 1000 | 0.426 | no effect | fhat 0.00; covariate cannot explain (t cap below 2.78) |
| D3 CF without Control1-5 / median matched CF | 0.157 | 39 matched drops | 39 | NA | descriptive |  |
| D3 drop P02: R | 0.825 | 13 matched drops of 4 healthy | 13 | 0.077 | descriptive |  |
| D7 drop P01: R | 0.979 | matched drops of 4 AML | 1222 | 0.339 | descriptive |  |
| D6 calendar year: R | 1.025 | 1,000 residual-permutation draws | 1000 | 0.652 | descriptive |  |
| D6 early stratum: R | 1.027 | 500 subsets (P01, P02 whole or absent) | 500 | NA | descriptive |  |
| D2 within-healthy: share of P370 with Control1-5 sign opposite to the imprint | 0.710 | 252 splits of the 10 healthy | 252 | 0.873 | descriptive | p_ge 0.143 |
| D5 leave-pair-out AUC with Control4 coded female | 0.913 | 150 pairs (A1-G02 D1 = 0.927) | 150 | NA | descriptive |  |

Calibration (healthy-ward plant on Control1-5, x main logFC): 0.25: R 0.705, percentile 0.050; 0.50: R 0.623, percentile 0.025; 1.00: R 0.506, percentile 0.025

**Verdict: GATE_A1B PASS.**
<!-- A1-BATCH end -->
