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
