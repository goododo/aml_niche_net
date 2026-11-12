# Amendment 1 to freeze 1: the SE(logit π0) gate of Stage 2 (2026-11-12)

Status: **declared post-hoc.** This change was decided after the Stage 2 verdict had been seen.

## What was seen first
- Gates job 3688071 ran with commit a2e0b47 and the pigeonhole SE.
- Its verdict was **BN-1, with 7 eligible patients (the bar is 10)**. The outputs carry no suffix in DBL_TAB_DIR (`dbl_share_stage2_*.csv`) and are kept unchanged.
- Without the SE gate, 12 patients would be eligible (post-hoc sensitivity). The SE gate decided the verdict.

## The defect
- Freeze 1 (section 4) asks for SE(logit π0) ≤ 0.10 "from 2,000 parent-cluster bootstrap replicates". Before any full calibration number existed, Stage 2 implemented this as a pigeonhole bootstrap: every parent cell draws a Poisson(1) weight, and a synthetic doublet carries the product of its two parents' weights (STAGE2_plan.md section 3c).
- With product weights, the barcode-level variance enters three times: Var(W_a W_b) = 3 for a pair that is seen once (Owen 2007, "The pigeonhole bootstrap").
- Post-hoc diagnostic on 4 lanes, 500 replicates:

  | Lane | Pigeonhole | One-way parent clusters | Barcode-only |
  |---|---|---|---|
  | AML3948 | 0.086 | 0.046-0.053 | 0.047 |
  | AML025 | 0.072 | 0.039-0.047 | 0.036 |
  | AML0612 | 0.176 | 0.094-0.110 | 0.094 |
  | AML055 | 0.101 | 0.058-0.059 | 0.058 |

  Parent clustering adds little, and the pigeonhole is about √3 times a cluster-robust SE.
- So the implementation measured a quantity about 1.7 times the calibration SE that the gate is meant to bound.

## The decision
- The user decided on 2026-11-12: treat this as an implementation defect and recompute 2b only, with a parent-cluster-robust SE.
- The freeze says "no third exception; a machinery rerun only for a demonstrated code bug". This amendment is recorded as that rerun, **made after the result was known.** It is labelled so wherever Stage 2 is reported.

## The amended rule (config_dbl.R section 11e; 12_share_calibrate.R --gates)
- **SE(logit π0)** is the square root of the dyadic cluster-robust variance (Aronow, Samii & Assenova 2015, *Political Analysis*).
  - Take the delta-method influence of each synthetic doublet on the patient-lane's doublet-level logit π0.
  - Sum φ_i φ_j over every pair of doublets that share at least one parent cell, including i = j.
  - Shared cells are counted in either role.
- **Unchanged:** every other gate, the gate unit (patient × lane), the bar 0.10, the calibration files, all seeds, and the pigeonhole bootstrap replicates. The replicates are still stored, and their SE is reported as `se_logit_pi0_pigeonhole`.
- **Checks:**
  - unit tests: unique parents give Σφ², and the result equals a brute-force pair sum;
  - regression: every 2b column except the SE equals the original run, and the stored pigeonhole SE equals the original SE.
- **Outputs** carry the suffix `_amend1`. 12 --replica reads the amended tables.

## Reporting
Every Stage 2 statement gives both verdicts: the registered one (pigeonhole, BN-1, n = 7) and the amended one, labelled "amendment 1, decided after the registered verdict was seen". Any later result that rests on the amended eligibility carries the same label.
