# Pre-registration: the paired-sample gate on the FGW transport score

Written and committed **before any of the runs below has produced a number**.
That is the only thing that makes this a pre-registration rather than a
description. The August screen was amended after a result existed and every
surviving feature came from the amendment; this document exists so that does
not happen twice.

## Why this exists

Every topological result in this project is negative, and the planted-effect
harness showed the alpha = 1 score has no power against the class of difference
this pipeline makes. Two explanations remain and they have not been separated:

1. the **graph representation** carries no usable signal, or
2. the **edge-to-cost transform** destroys it, because 46.9% of every cost
   matrix is a placeholder whose value is set by how sparse the sample is.

Paired samples are the cleanest material we have for (1): same donor, same
laboratory, same protocol. Four distance variants address (2). The two are
crossed, so each gate is evaluated under each variant.

## What is being tested, in order

### GATE 1 . identity. Does the representation carry the patient?

For a patient with graphs at two timepoints, ask whether the two graphs are
closer to each other than to graphs from other patients.

- distance is **FGW2(graph_i, graph_j)** computed directly between two patient
  graphs. The healthy barycentre is not involved, so nothing in this gate
  depends on the healthy arm, the malignant label, or any covariate.
- **The candidate pool is restricted to the same dataset.** Without this the
  test is won by batch: a patient's two samples are trivially from the same
  study. This restriction is the whole reason the gate means anything.
- statistic: for each paired patient, the **percentile rank** of the true
  partner among all same-dataset AML graphs. Under the null this is uniform on
  [0, 1] with median 0.5.
- **PASS**: one-sided Wilcoxon signed-rank of those percentiles against 0.5,
  p < 0.05, in the direction of the true partner ranking closer.
- reported alongside, as description and not as the criterion: top-1 and top-5
  retrieval rate against the chance rate implied by each dataset's pool size.

This gate assumes no biology at all. It asks only whether the representation is
reproducible within a donor. If it fails, nothing downstream can work and no
amount of distance engineering repairs it.

### GATE 2 . the positive control. Diagnosis to post-treatment

11 pairs. Between diagnosis and post-treatment the marrow goes from
blast-packed to regenerating, so a difference is close to guaranteed by
biology. This is the gate the phrase "if it cannot do it here, stop looking"
actually applies to.

- statistic: HDS at each timepoint, then the within-pair difference.
- **the sign is pre-specified**: treatment moves the marrow toward normal, so
  **HDS_treatment < HDS_diagnosis**. One-sided.
- **depth is not removed by pairing.** Measured within-pair depth ratios run
  0.23x to 3.54x. So the test is a regression of the within-pair HDS difference
  on the within-pair difference in log10 median UMI, and the quantity tested is
  the **intercept**.
- **PASS**: intercept negative with p < 0.05.

### GATE 3 . the question. Diagnosis to relapse

11 pairs. This is a scientific question with an unknown answer, **not a gate**.
A null here is a finding about relapse, not evidence that the pipeline failed,
and it will not be used to stop anything.

- pre-specified direction, from the original blueprint H3: relapse moves away
  from healthy, **HDS_relapse > HDS_diagnosis**. Reported two-sided.
- same depth handling as GATE 2.

## The four distance arms, all reported

Every gate is run under all four. The point is not to choose a winner but to
find out whether the alpha = 1 null is a property of the data or of the
transform.

| arm | what an undetected edge receives |
| --- | --- |
| `rank` | the current production rule, `C = 1 - rank_pct(weight)`, so absent edges tie at the bottom and take `1 - (k+1)/98` |
| `mask` | nothing. Rank only the detected edges among themselves and leave absent ones out of the cost. `C_perlr` and `detected` are already written by `06_distance/01` and are consumed by nothing today |
| `const` | a fixed `C = 1`, independent of how many are missing |
| `logs` | `C = -log(S + eps)`, normalised then min-max scaled. Specified in the design blueprint and never implemented |

**No arm is designated primary in advance and none may be selected after the
fact.** All four numbers go into the report for every gate. If they disagree,
that disagreement is the result.

## Decision rules, fixed now

1. **GATE 1 fails in all four arms** -> the graph representation does not carry
   the patient. The FGW transport score is retired as a primary readout. The
   project continues with per-edge regression and per-feature screening, both of
   which have demonstrated power in the planted-effect harness. The network
   question is not abandoned; the instrument is.
2. **GATE 1 passes but GATE 2 fails in all four arms** -> the representation is
   reproducible but the score cannot see a difference that is certainly there.
   Same consequence for the transport score, and it is a stronger result than
   (1) because it localises the failure to the scoring step.
3. **GATE 1 and GATE 2 pass in at least one arm** -> that arm becomes the
   production distance, the alpha sweep and the planted-effect harness are rerun
   under it, and the August topological null is **withdrawn and recomputed**.
4. **GATE 3 is reported whatever it says** and never used to stop or continue
   anything.
5. No gate is re-run with a different statistic, a different covariate set, or a
   different candidate pool to change a verdict. The arms above are the whole
   analysis.

## Two free parameters, fixed here before any run

Neither is set by the blueprint, so both are pinned now rather than chosen later.

- **The alpha grid.** GATE 1 runs at **alpha 0, 0.5 and 1**, and all three are
  reported. Identity is a property of the representation as a whole, and running
  the grid says *where* the identity lives: in the node features, in the
  topology, or in both. GATE 2 and GATE 3 run at **alpha 0.5**, the production
  setting, and at **alpha 1**, the term actually under suspicion.
- **The epsilon in the `logs` arm.** `eps = 1e-3`, applied to a within-sample
  max-normalised weight, so an undetected edge lands at `-log(1e-3)` and, after
  the min-max rescale, at exactly `C = 1` in every sample regardless of how many
  are missing. That constancy is the entire point of the arm.

## What the four arms do to an undetected edge, stated as formulas

Let `w` be `weight_probsum` over the 49 directed edges of one sample, and `k`
the number with `w = 0`.

| arm | detected edges | undetected edges |
| --- | --- | --- |
| `rank` | `1 - rank(w over all 49)/49` | `1 - (k+1)/98`, which is why the value moves with `k` |
| `mask` | `1 - rank(weight_per_lr over detected only)/n_detected` | `1.0`, fixed |
| `const` | `1 - rank(w over all 49)/49` | `1.0`, fixed |
| `logs` | min-max of `-log(w/max(w) + 1e-3)` | `1.0`, fixed |

`mask` and `const` differ in the detected edges, not the absent ones. Under
`rank` and `const`, a sparse sample's detected edges are compressed into the top
of the ranking; under `mask` they are spread across the full range. That
difference in geometry is what the two arms separate.

## Reproduction self-checks, both must pass or the run is void

1. The `rank` arm must reproduce the `C` column of
   `results/tables/06_distance/edge_distance.csv` bit for bit. This proves the
   re-implemented transform is faithful.
2. The `rank` arm at alpha 0.5 must reproduce `HDS` in
   `results/tables/07_fgw/patient_scores.csv`. This proves the re-implemented
   barycentre and scoring are faithful.

If either fails the script aborts and no gate result is written.

## What this cannot settle, and must not be used to claim

- **Sample size.** 11 pairs. A paired test at 80% power detects Cohen d of about
  0.81. A GATE 3 null is uninterpretable on its own and only becomes
  interpretable if GATE 2 passed.
- **Which dataset carries it.** 9 of the 11 relapse pairs and 7 of the 11
  treatment pairs are GSE227903. GSE201966 contributes 2 relapse pairs against
  only 3 same-dataset distractors, so its contribution to GATE 1 is close to
  nil. Per-dataset numbers are reported so this is visible.
- **The arm.** All paired samples sit in the Validation split. These gates test
  the instrument, not a screened feature, so this is not double-dipping on the
  screen hypothesis. It does mean the Validation arm is no longer naive for
  anything involving the transport score, and that goes in the Methods.
- **Anything about the malignant label.** GATE 1 never reads it. GATE 2 and 3
  inherit it only through `frac_malignant` inside the FGW feature term, which is
  a known circularity recorded elsewhere and is not repaired here.
