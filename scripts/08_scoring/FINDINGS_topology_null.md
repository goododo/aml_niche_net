# Findings: the communication-graph topology line

Status as of 2026-09-01. Every number is reproduced by a named script on a named sample set.
Nothing here is an estimate or a recollection. Where a result is mixed or unresolved it is written
as mixed or unresolved.

The hypothesis was that AML bone marrow differs from healthy marrow in the **topology** of its
cell-cell communication graph, and that this topology shifts with treatment and relapse.

The work splits into three parts. **Part I** is settled and negative: the graphs change within a
patient, but with no direction shared across patients, and AML does not differ from healthy in
topology. **Part II** diagnoses why the statistic could not have found anything even if there had
been something, and localises the cause to the within-sample rank transform. **Part III** is the
attempted correction, run as a 2x2 ablation. It settles three things - the damped-walk cost is
rejected on its own positive control, the planted-effect detection is an interaction of both
changes rather than either one, and the AML-vs-healthy topology question is closed - and it opens
one: replacing the cell-count node mass with LR-signal mass takes the patient-identity test from
5 of 24 cells to 5 of 6, at the cost of the positive control falling from p 0.027 to 0.068. That
last one is a pattern, not a result; section 14 says exactly how much weight it can carry.

**Where the whole line stands**: the original hypothesis is answered and the answer is no. What
survives is a methodological finding about why a common OT formulation is blind here, and one
untested lead about node mass.

---

## 0. Sample sets

Naming the sample set is part of every claim, because several earlier results in this project were
read as cohort-wide when the code had silently used a subset.

| Set | n | Composition |
|---|---|---|
| Full graph cohort | 138 | 115 AML, 23 healthy, 10 datasets |
| Healthy barycentre | 19 | healthy donors, sparse-flagged excluded |
| **Paired samples** | **37** | 22 pairs, 15 patients, 3 datasets, all AML timepoints |

Paired breakdown: GSE227903 16 pairs (9 relapse, 7 treatment), GSE116256 4 pairs (all treatment),
GSE201966 2 pairs (both relapse). Contrasts are Dx-to-Treatment n=11 and Dx-to-Relapse n=11.

**All 22 pairs sit in the Validation arm** of `01_preprocess/02_study_split.csv`: every longitudinal
dataset was assigned there by design. There is no held-out paired cohort, so the paired analyses are
exploratory and cannot be called validated.

---

# PART I - What was established

## 1. C is a reliable measurement, not noise

Each paired sample's cells were split in half, stratified by hierarchy bin (every bin balanced to
within one cell), and the identical CellChat path run on each half. 74 runs, all completed.

| Quantity | Value |
|---|---|
| **Per-sample** Spearman(w_A, w_B), summarised over the 37 samples | **median 0.858**, IQR [0.696, 0.915] |
| Edge detected/absent agreement | **93.8%** |
| Bins surviving in both halves | median 5 of 7 |
| Testable edges | median 25 of 49 |

**Which Spearman, and why it matters.** The reliability figure is the *per-sample* correlation:
one Spearman computed inside each sample across its own testable edges, then summarised over the
37 samples. Pooling all 863 testable sample-edges into a single correlation instead gives
**0.867**, and that number is **not** the reliability: samples differ in overall signal level, so a
strong sample's edges sit high on both axes and lift a pooled correlation without being evidence
that any one sample reproduces itself. The two happen to be close here, which says the inflation is
small in this cohort, not that the distinction is unimportant.

**What both summaries hide.** The per-sample values run from **-0.021 to 0.983**. At least one
sample's two halves do not correlate at all. A median of 0.858 is the right headline but it is not
a guarantee about any individual sample, and any per-sample claim needs that sample's own value.

Restricted to edges whose sender and receiver both clear `CCC_MIN_CELLS_PER_NODE = 10` in **both**
halves, so node dropout from halving the cells is reported separately, not counted as disagreement.

## 2. Within-patient change is real, but has no shared direction

**Real.** Comparing like with like, the same edge set on both sides of each ratio:

| Contrast | pairs used | edges | within-patient median abs dC | split-half noise floor | ratio |
|---|---|---|---|---|---|
| Dx-to-Treatment | 10 of 11 | 16 | 0.1582 | 0.0408 | **3.88** |
| Dx-to-Relapse | 11 of 11 | 25 | 0.1122 | 0.0459 | **2.44** |

One Dx-to-Treatment pair is dropped: fewer than 3 edges were testable in both halves of both its
timepoints. The noise floor is measured at half the cells, so it overstates the noise relevant to a
full-sample comparison and the ratios are conservative. No correction factor is given: C is a rank
transform and its noise does not scale as a clean 1/sqrt(n).

**No shared direction.** Four independent measurements agree:

1. Variance decomposition (approximate, not a formal ANOVA): patient identity accounts for
   51.7-74.1% of the variance of C, the timepoint main effect for **2.6-5.1%**.
2. Per-edge paired test over 49 directed edges, two statistics (Wilcoxon signed-rank; and the
   intercept of `dC ~ d(logdepth)`, the form GATE 2 uses): **0 of 49 survive BH at q<0.05** in
   either contrast.
3. Family-level sign-flip permutation, 10,000 draws. Decisive, because it does not depend on BH:

   | Contrast | observed edges at raw p<.05 | chance | p |
   |---|---|---|---|
   | Dx-to-Treatment | 3 / 49 | 2.11 | **0.280** |
   | Dx-to-Relapse | 2 / 49 | 2.26 | **0.495** |

   A common effect too weak to survive BH would still push the raw count above chance. It does not.
4. Cross-dataset consistency, GSE227903 (7 pairs) vs the rest (4 pairs), Dx-to-Treatment: Spearman
   of the 49 per-edge median dC = **-0.002**; 16 of 49 edges agree in sign, below the 24.5 expected
   by chance.

**Power bound, stated explicitly.** With n=11 the smallest attainable two-sided Wilcoxon p is
0.00098 and the BH threshold for the most significant of 49 edges is 0.00102, so passing BH needs
essentially all 11 patients moving the same way on one edge. The strongest observed edges reached
9 of 11 (raw p 0.0039-0.0059, q 0.19-0.24). Point 3 is what makes the conclusion stronger than this
bar alone would support.

## 3. No AML vs healthy difference in topology

Pure topology (alpha = 1), healthy barycentre, within-dataset permutation, **mass mode `ncells`**:

| Distance arm | mean healthy | mean AML | beta | p |
|---|---|---|---|---|
| rank (production) | 0.0539 | **0.0445** | -0.00150 | 0.966 |
| mask | 0.0590 | 0.0731 | +0.00736 | 0.389 |
| const | 0.0700 | **0.0645** | +0.00517 | 0.799 |
| logs | 0.0561 | 0.0607 | +0.01582 | 0.524 |

A transport-free baseline agrees: the plain mass-weighted squared difference between cost matrices
gives healthy 0.0590 vs AML 0.0611, beta +0.0149, permutation p = 0.679 (20,000 draws).

**The direction is not consistent and is partly an artefact of the scoring rule.** Two arms put AML
closer to the healthy barycentre than healthy samples are, two put it farther. Healthy samples are
scored **leave-one-out** (`06_alpha_sweep.py:243`), which is the right thing to do - it removes the
circularity of scoring a sample against a barycentre it helped build - but it also costs healthy
samples a distance that AML samples never pay. Measured at alpha=1:

| arm | reported gap, healthy(LOO) - AML | LOO overhead on healthy | gap if healthy also used the full barycentre |
|---|---|---|---|
| rank | +0.0093 | +0.0038 | +0.0056 |
| const | +0.0052 | +0.0048 | **+0.0004** |
| mask | -0.0134 | +0.0032 | -0.0166 |
| logs | -0.0047 | +0.0037 | -0.0084 |

For `const` the LOO overhead accounts for essentially the whole gap; for `rank` it accounts for
about 40% and AML is still nearer without it. None of this changes any conclusion, because every
one of these gaps is null (p 0.39-0.97) - but the sign of a null gap should not be read as a
biological direction, and no figure or sentence should present it as one.

Pure features (alpha = 0) do separate: healthy 0.1161 vs AML 0.3337, beta +0.122, p = 0.016 under
mass mode `ncells`, and healthy 0.0778 vs AML 0.1653, beta +0.094, p = 0.0009 under `uniform`. Both
are correct; they are different models, and any figure quoting one must say which.

**That separation is not a discovery.** Of the three features actually in the FGW distance,
`mean_stemness` has node-level AUC 0.516 and `n_cells` 0.520, both at chance, while `frac_malignant`
is set to a single constant for every healthy row by the deliberate `FGW_ZERO_HEALTHY_MAL` switch
(`07_fgw/01_build_fgw_inputs.R:104`) - even though the same 23 donors carry measured malignant
fractions from 0.0101 to 0.5299 elsewhere in the pipeline.

**Five weight-to-cost formulas were tried, but not all of them on this question.** The four
registered arms - `rank` (production), `mask`, `const`, `logs` - have all now been swept against
AML vs healthy in both mass modes; the `mask` sweep was missing from an earlier version of this
document and was run on 2026-09-02 (alpha=1, `uniform` mass: healthy 0.0607, AML 0.0601, beta
-0.00246, p 0.697, also null). The fifth form, raw `log1p(weight_probsum)`, was tested only at the
per-edge level, never as an FGW cost: it finds 5 of 49 edges at q<0.05 pooled across datasets and
**0 of 49 within dataset** - those five are platform effects, which is why the rank transform
exists.

## 4. The pre-registered gates: GATE 1 is where pure topology did work

`PREREGISTRATION_paired_gate.md` was written and frozen before these were run. GATE 1 asks whether a
patient's own second timepoint ranks closer to their diagnosis graph than other patients' graphs
from the same dataset (statistic = median percentile of the true partner; 0.5 = chance, lower is
better), Dx-to-Treatment, n=11:

| arm | alpha=0 | alpha=0.5 | alpha=1 |
|---|---|---|---|
| rank | 0.273 (p .062) | **0.182 (p .026)** | 0.364 (p .084) |
| mask | 0.273 (p .062) | 0.273 (p .124) | 0.636 (p .719) |
| const | 0.273 (p .062) | **0.130 (p .024)** | **0.273 (p .021)** |
| logs | 0.273 (p .062) | **0.130 (p .022)** | **0.261 (p .0015)** |

Dx-to-Relapse fails in every cell. GATE 2, the positive control (Dx-to-Treatment, direction
pre-specified as intercept < 0), passes at alpha=0.5 in **all four arms** (p .027 / .030 / .029 /
.031, intercepts -0.075 / -0.074 / -0.077 / -0.078) and fails at alpha=1 in all four. GATE 3 is null
throughout.

The pre-registered retirement condition - GATE 1 failing in all four arms - did **not** occur, so
the transport score was not retired; decision rule 3 fired instead.

**Disclosed gap in the registration**: GATE 1 spans 24 grid cells and the registration specified no
multiplicity correction. Under BH only `logs` at alpha=1 survives (q = 0.035). That cell is the
finding; the uncorrected cells are not.

## 5. Resolution is not the limiting factor

Counting node types clearing 10 cells across the 37 paired samples:

| Resolution | types | in 100% of samples | in >=90% | in >=75% |
|---|---|---|---|---|
| `hierarchy_bin` | 7 | 1 | **3** | 5 |
| `bmm_broad` | 23 | 0 | **3** | 7 |
| `bmm_fine` | 54 | 0 | **2** | 6 |

Per-sample usable nodes rise with resolution (median 6, 13, 19) but the shared backbone does not: at
any resolution only 2-3 node types are reliably present across samples, because different samples
contain different cell populations, a composition effect already traced to sample preparation.
**Caveat added 2026-09-01**: this uses a hard >=10-cell presence threshold. The pipeline already
tolerates absent nodes through `EPS_MASS`, and scACCorDiON forces a common node set the same way, so
finer resolution is a design choice rather than something the data forbids. An earlier version of
this document stated it too strongly.

---

# PART II - Why the statistic could not have found anything

## 6. The information was encoded in the one dimension GW discards

Measured as the between-sample SD of the sorted spectrum (value information) against the SD at a
fixed edge position (value + arrangement):

| arm | SD of sorted C | SD of C at edge | ratio |
|---|---|---|---|
| **rank (production)** | 0.0833 | 0.1998 | **2.40** |
| const | 0.1511 | 0.2972 | 1.97 |
| mask | 0.1271 | 0.2786 | 2.19 |
| **logs** | **0.2264** | 0.3377 | **1.49** |

Gromov-Wasserstein re-optimises its coupling, which is exactly the operation that discards
arrangement. **The pipeline encoded its signal in the one dimension its statistic cannot read.**

One independent observation lines up with this and was not used to derive it: `logs` has the lowest
ratio, i.e. preserves the most value spectrum, **and** is the only arm passing GATE 1 at alpha=1
(p=0.0015) and the only arm with a monotone alpha=1 trend in the paired planting.

**Two claims that an earlier version of this document made here have been withdrawn**, both after
being checked rather than argued:

1. *"Rank leaves every sample holding nearly the same 49 numbers."* **Overstated.** As a fraction
   of its own value range the production spectrum varies by 8.5%, and the corrected geometry's by
   8.1% - the same. The raw SDs (0.083 vs 0.395) differ only because the value ranges do (1.0 vs
   4.9). What is true is narrower and was measured separately: the production spectrum's variation
   is substantially a restatement of **how many edges were detected**, regressing the spectrum on
   the detected-edge count gives mean R^2 = **0.421** across the 49 positions (8 positions above
   0.9), against **0.171** for the corrected geometry (0 positions above 0.9). The production
   spectrum carries less *independent* information, not less variation.
2. *"The sorted spectrum separates AML from healthy in 0 of 49 positions for production but 15 of 49
   for the corrected geometry."* **Does not survive dataset stratification.** Within dataset the
   counts are **2 of 49** and **4 of 49** against a chance expectation of 2.45, i.e. neither
   separates. The unstratified 15 was a dataset effect. Only 3 of 10 datasets contain both AML and
   healthy samples, so the stratified test runs on 72 of 138 samples.

## 7. The consequence, measured

**The statistic cannot see an effect planted into it.** Planting into 6 of 49 edges, additive mode,
largest tested amplitude:

| arm | targets moved | mean abs dC | planted edges recovered per-edge | FGW omnibus p |
|---|---|---|---|---|
| logs | 97.7% | 0.639 | **6 / 6 at q<0.05** | **0.555** |
| const | 97.7% | 0.672 | **6 / 6 at q<0.05** | **0.871** |

**The coupling has collapsed toward the uninformative solution.** With 7 nodes the product coupling
p (x) q has diagonal mass 1/7 = 0.1429 and the identity coupling has 1.0. Measured:

| Block | n | diagonal mass | distance to identity | distance to product |
|---|---|---|---|---|
| healthy | 23 | 0.147 | 1.013 | 0.456 |
| paired | 37 | 0.153 | 1.011 | 0.384 |
| all AML | 115 | 0.172 | 0.998 | 0.380 |

**Confirmed by measurement, not inference.** The same planting scored two ways:

| delta | GW (alpha=1, optimal T) | Frobenius (no transport) |
|---|---|---|
| 0.25 | +7.7% | +11.9% |
| 1.0 | **+14.3%** | +28.7% |
| 4.0 | **+5.9%** | +36.9% |

GW is **non-monotone**: it peaks at delta=1 then falls back as the planted effect grows, because the
coupling re-optimises to absorb the imposed structure. The transport-free statistic on the same data
rises monotonically. (Effect sizes, not significance tests; the omnibus p-values are above.)

**Publishable form:** on a 7-node communication graph, the GW term of FGW loses sensitivity to
edge-level change through coupling degeneracy, and that sensitivity *decreases* as the perturbation
grows.

## 8. The cost matrix is not symmetrised - verified, not assumed

POT 0.9.7 with `symmetric=False`: FGW(C, Cb) = 0.050927 and FGW(C^T, Cb) = 0.051522, which differ;
the transpose-equivariance identity FGW(C, Cb) = FGW(C^T, Cb^T) holds to 6.9e-17, the correct
property of an asymmetric GW loss; explicitly symmetrising the input changes the value from 0.0509
to 0.0289. `symmetric=None` (auto-detect) gives the same asymmetric answer.

Direction is therefore used. But its share is modest: relative asymmetry
norm(C - C^T) / norm(C - mean(C)) has median **0.680**, and the two directions of a cell-type pair
correlate at rho = **+0.752**, so most of the production C's variation is a pair-level component
tracking how abundant both cell types are rather than who signals to whom. The corrected geometry of
Part III raises that asymmetry share to **1.486** and drops the redundancy to rho = -0.192.

---

# PART III - The corrected formulation: what was tried, and what it settled

## 9. Why the formulation was questioned at all

Part II established *how* the statistic fails but not *what to do instead*. The prompt came from
scACCorDiON (Nagai et al., *Bioinformatics* 2025, 41(5):btaf288), which compares CCC graphs with
optimal transport successfully. Reading its methods against ours showed four differences, and one
of them is not a tuning choice but a category error:

| | scACCorDiON | this pipeline |
|---|---|---|
| distance | Wasserstein | Gromov-Wasserstein |
| points mass moves between | edges (via a line graph) | nodes |
| **cost** | **graph structure (hitting time on the Markov chain)** | **LR signal (rank of edge weight)** |
| **mass** | **LR signal** | **cell counts** |
| node correspondence | given, `"all samples k have the same cell types"` | solved for by GW |
| question | clustering into patient subgroups | two-group hypothesis test |

**The two roles were swapped.** In an OT comparison of graphs the structure supplies the cost and
the signal supplies the mass. This pipeline used the signal it wanted to test as the geometry, and
used cell composition - the quantity already traced here to sample preparation - as the signal.

A second point follows from the same table and matters for how GW was used at all: our seven bins
are labelled and shared across every sample, so the correspondence is **known**. GW exists to
recover a correspondence that is unknown. Using it here paid the full price of its relabelling
invariance and bought nothing, which is exactly the degeneracy Part II measured.

## 10. What was changed, and the decision to change two things at once

`config/graph_geometry.py`, one function, imported by both the builder and the planting harness so
the two cannot drift:

- **cost** = `-log R`, where `R = (1 - alpha) * inv(I - alpha * P)` is the damped random-walk
  reachability of the directed graph. `P` is the row-normalised weight matrix, mixed with a small
  uniform teleport. Continuation 0.85, teleport 0.05.
- **mass** = node signal strength, sent plus received, normalised to sum 1.

**Why a random walk rather than scACCorDiON's hitting time**: both turn local edge weights into a
global geometry; reachability is bounded in (0, 1) and needs no arbitrary scaling, and it was the
smaller change to a pipeline whose downstream code expects a bounded cost.

**Why no within-sample normalisation, which is the whole point**: row-normalising `W` into a
transition matrix already removes the absolute depth scale, which is the job the rank transform was
introduced to do. The reachability matrix is a probability, so it is comparable across samples with
no further rescaling - and unlike ranking, the value spectrum survives. Part II identified the rank
transform as the root cause; this replaces it without giving up what it was for.

**Two implementation choices that had to be made and are reported, not hidden**: 25.5% of node-rows
carry no outgoing signal, so their transition row is undefined and is completed to uniform (the
standard dangling-node convention). With a pure walk, 15.9% of ordered pairs then turn out to be
unreachable and their cost lands on an arbitrary `-log` floor that dominates the scale. The teleport
term fixes both: it makes the chain irreducible so every entry is strictly positive and no floor is
needed, and it turns the dangling convention into a limiting case of one rule rather than a special
case. Measured after the fix: floor entries 0 of 6762.

**This changes two things at once (cost and mass), which was a deliberate risk.** The two belong to
one conceptual correction, and the planting harness was expected to give a clean yes/no. It did not,
which is what forced the ablation in section 14.

**None of this is pre-registered.** The registration covers four weight-to-cost transforms; this
changes which quantity *is* the cost.

## 11. Step one: validate on a planted effect before touching biology

**Why this order.** A formulation that cannot see an effect deliberately put into it cannot be
trusted to report that no effect exists. Part II had just shown the production statistic failing
exactly this test, so running biology first would have repeated the mistake.

**How.** The identical harness, the same six pinned edges, the same delta grid, additive mode. The
delta=0 row is an unplanted control that must stay null.

**Result, k=1, the first time anything has passed:**

| delta | mean abs dC | mean healthy | mean AML | **omnibus p** |
|---|---|---|---|---|
| 0 (control) | 0 | 0.3345 | 0.3403 | **0.5332** |
| 0.25 | 0.817 | 0.3345 | 0.3123 | 0.1264 |
| 0.5 | 1.033 | 0.3345 | 0.3028 | 0.0675 |
| 1 | 1.259 | 0.3345 | 0.2952 | **0.0270** |
| 2 | 1.468 | 0.3345 | 0.2862 | **0.0110** |
| 4 | 1.638 | 0.3345 | 0.2928 | **0.0060** |

Monotone, control still null. The production geometry on the same edges is flat at p 0.52-0.61 and
non-monotone.

**Limits, stated:** k=3 reaches only p=0.072 and k=6 only p=0.093 at the largest delta. This detects
a large, concentrated, consistent change, not a diffuse one. `mean abs dC` is on a `-log` scale here
and is **not** comparable to Part II's [0, 1] numbers.

## 12. Step two: establishing that the sign of an HDS shift means nothing

**Why.** Planting the same six edges one at a time flipped the direction of the HDS shift between
k=2 and k=3. Carrying an unexplained sign flip into real data would have meant not knowing how to
read a result if one appeared.

**How.** Two hypotheses were tested and both were dropped before the right one was found, which is
recorded here because the failures were informative:

1. *Walk concentration.* Refuted: k=1 and k=6 have nearly identical entropy shifts (0.0914 vs
   0.0845) and opposite HDS directions.
2. *Alignment with the toward-barycentre direction.* Only 3 of 6 correct on its own.
3. *The full geometric criterion.* For a squared distance to a fixed reference, dist^2 falls only if
   `cos(theta) > norm(displacement) / (2 * norm(toward))`. **6 of 6 correct.**

| k | cosine | norm(disp) | threshold | predicted | observed dHDS |
|---|---|---|---|---|---|
| 1 | +0.243 | 2.743 | 0.243 | fall | **-0.0475** |
| 2 | +0.098 | 3.080 | 0.273 | rise | +0.0172 |
| 3 | -0.060 | 3.674 | 0.325 | rise | +0.1831 |
| 4 | -0.004 | 4.585 | 0.406 | rise | +0.1715 |
| 5 | +0.026 | 4.201 | 0.372 | rise | +0.1801 |
| 6 | +0.107 | 3.949 | 0.350 | rise | +0.2407 |

k=1 sits exactly on its own threshold, which is why it is the only one that falls.

**Consequence.** The planted directions are near-orthogonal to the AML-to-healthy axis
(abs cosine <= 0.24), and orthogonal displacement always increases a distance. **HDS answers "how
different", never "different in which direction".** Intrinsic to any distance-to-a-reference score,
not to this implementation, and it must be carried into any write-up that reports an HDS sign.

## 13. Step three: real data, and the positive control breaks

**AML vs healthy, walk cost + signal mass:** alpha=1 gives healthy 0.3334, AML 0.3368,
within-dataset beta -0.0524, **p = 0.507**. Still null.

**The gates:** GATE 1 passes only at alpha=0 (p .026 / .038) and is *worse* than the registered arms
at alpha=1 (p 0.115 against `logs` 0.0015). **GATE 2, the pre-specified positive control, fails
(p = 0.380)** where all four registered arms pass (p .027-.031).

**A formulation that fails its positive control cannot support a negative result.** The corrected
geometry's alpha=1 null is therefore not stronger evidence than the old one, and is not claimed as
such anywhere in this document.

## 14. Step four: the 2x2 ablation

**Why.** Section 10's deliberate risk came due. The change moved the cost and the mass together, so
neither the planted-effect success nor the GATE 2 failure could be attributed. The hypothesis on
record before running was: *GATE 2's signal is 100% feature-driven (Part II: planting topology left
its intercept at -0.0750 across every delta), the feature term is transported by the mass, so
replacing the cell-count mass is what broke it.*

**How.** Two further arms completing the square, run through the identical harness and gates:
`walk_ncells` = walk cost with the production cell-count mass, and `rank_signal` = production rank
cost with the signal mass. Total compute for all four configurations across three tests: about four
minutes (gates 9-14 s per arm, planting 66-68 s, alpha sweep 13-14 s).

**Two built-in consistency checks passed before any result was read.** The alpha=0 row of the sweep
must group by *mass*, not by cost, because the cost matrix does not enter the objective at alpha=0.
It does: `walk+ncells` reproduces production exactly (0.1161 / 0.3337 / p 0.0161) and `walk+signal`
reproduces `rank+signal` exactly (0.0953 / 0.1954 / p 0.0205). And the default path still reproduces
the frozen sweep and all 40 rows of `paired_gate.csv` to 0.000e+00.

### The complete square

| cost | mass | GATE 1 cells at raw p<.05 | GATE 2 (alpha=0.5) | planted k=1 best p | alpha=1 AML vs healthy |
|---|---|---|---|---|---|
| rank | ncells (**production**) | 5 of **24** (4 arms) | **0.027-0.031** | 0.524 (flat) | 0.966 |
| rank | signal | **5 of 6** | 0.068 | 0.111 (monotone) | 0.672 |
| walk | ncells | 2 of 6 | 0.508 | 0.387 (monotone) | 0.717 |
| walk | signal | 2 of 6 | 0.380 | **0.006** | 0.507 |

### What it settles

**1. The registered hypothesis is refuted, and the walk cost is out.** Restoring the cell-count mass
did **not** bring GATE 2 back; it made it worse (0.508 against walk+signal's 0.380). The walk cost
breaks the positive control under **both** mass modes, so the cost change is responsible and the
mass change is not. The damped-walk cost is therefore rejected on its own positive control.

**2. The planted-effect detection is an interaction, not attributable to either change.** rank+signal
reaches only p=0.111 and walk+ncells only p=0.387; both are monotone, both are large improvements on
production's flat 0.524, and only the two together cross 0.05. An earlier reading of this document
attributed the detection to the mass alone; that was wrong and is corrected here.

**3. An unlooked-for result: the signal mass substantially improves the identity test.** With the
production cost unchanged, replacing the cell-count mass with signal strength takes GATE 1 from
5 of 24 cells to **5 of 6**:

| alpha | contrast | percentile of the true partner | p |
|---|---|---|---|
| 0 | Dx-to-Relapse | 0.333 | **0.0259** |
| 0 | Dx-to-Treatment | 0.348 | **0.0376** |
| 0.5 | Dx-to-Relapse | 0.333 | **0.0195** |
| 0.5 | Dx-to-Treatment | 0.217 | **0.0322** |
| 1 | Dx-to-Relapse | 0.217 | 0.1387 |
| 1 | Dx-to-Treatment | 0.261 | **0.0488** |

**Dx-to-Relapse passes here, and it never passed in any registered arm at any alpha.** The same
configuration is also the only one significant at every intermediate alpha of the AML-vs-healthy
sweep (p .0205 / .0186 / .0151 / .0170 at alpha 0 / .25 / .5 / .75, against production dropping out
at alpha=0.75 with p 0.143).

**The cost of it**: GATE 2 falls from 0.027 to 0.068, just past the line.

**How much weight this can carry.** No single GATE 1 cell survives BH over the six (smallest
q approximately 0.056). The six cells share the same 11 pairs and the alphas are nested, so they are
**not independent** and a naive family-level binomial test would be anticonservative; none is
reported. This is a pattern worth pursuing, not an established result.

**4. AML vs healthy in pure topology is closed.** alpha=1 is null in all four configurations
(p 0.507 / 0.672 / 0.717 / 0.966). Adding the four registered weight-to-cost arms in both mass
modes (rank p 0.966, mask 0.389, const 0.799, logs 0.524 at `ncells`; all null at `uniform` too) and
the transport-free baseline (p 0.679), every route to this question that has been run returns null.

### One bug found and fixed during this step

`rank_signal`'s first alpha sweep reproduced the production sweep to the last digit. That is not a
result, it is a tell: the new branch was gated on `startswith("walk")`, which `rank_signal` does not
satisfy, so it fell silently back to the production path. Fixed, rerun, and the default-path guard
re-checked at 0.000e+00. Recorded here because this class of failure - a new arm silently scoring
the old model and exiting 0 - is the one this project keeps meeting.

---

## 15. What can and cannot be claimed

**Can be claimed**
- C is a reliable measurement (split-half Spearman 0.858).
- Within-patient graph change between diagnosis and treatment is real, 3.88x the measurement noise.
- In this cohort there is no common direction to that change detectable at n=11, and the edge-level
  signal is indistinguishable from chance at the family level.
- There is **no** AML vs healthy difference in graph topology. The four registered weight-to-cost
  arms in both mass modes, two mass definitions, a transport-free baseline, and a corrected geometry
  validated against a planted effect all agree. (The fifth weight form, raw `log1p`, was only ever
  tested per-edge, not as an FGW cost; its 5 of 49 edges vanish within dataset.)
- The production FGW alpha=1 statistic is insensitive to edge-level change, by a mechanism that has
  been localised and measured, with the root cause traced to the within-sample rank transform.
- The damped-walk cost is rejected: it fails the pre-specified positive control under both mass
  definitions.
- A pure-topology distance does carry patient identity information (GATE 1), and the LR-signal node
  mass improves it markedly.

**Cannot be claimed**
- That no common change exists. n=11 per contrast; BH significance needs near-unanimity.
- That the paired findings are validated. Every paired dataset is in the Validation arm.
- That the alpha=0 feature separation is a biological discovery. One of its three features is
  constant-by-construction in healthy samples; the other two are at chance.
- **That `rank_signal` is established.** 5 of 6 GATE 1 cells at raw p<.05 is a pattern; no cell
  survives BH, the cells are not independent, and it costs the positive control (0.027 to 0.068).
- Any directional reading of an HDS shift (section 12).
- Anything about structural stroma. The graph is hematopoietic-immune only.

## 16. Limitations to carry into the write-up

1. n=11 pairs per contrast, 3 datasets, 16 of 22 pairs from one study.
2. "Post-treatment" is not one biological state: GSE116256 contributes day 20-113
   induction/consolidation samples, GSE227903 contributes MRD.
3. The split-half noise floor is measured at half depth and is therefore conservative.
4. `pt_predicted_Pseudotime` carries the off-trajectory placeholder zero in 138/138 samples for the
   T_NK bin. It is **not** among the three FGW features, so nothing above is affected, but
   `08_scoring/07_feature_decomposition.py` does consume it and any conclusion from that script
   needs rerunning with `BMM_PSEUDOTIME_OFFTRAJ` filtered.
5. The Discovery/Validation split has not been respected by any analysis in this line, because the
   longitudinal data is entirely on one side of it.
6. The walk geometry has two free parameters (damping 0.85, teleport 0.05) and no sensitivity
   analysis was run on either. Since the walk cost is now rejected on its positive control, that
   analysis is moot unless the cost is revisited.
7. 25.5% of node-rows carry no outgoing signal and are completed to uniform. The teleport term makes
   that a limiting case rather than a special rule, but a quarter of rows is a lot for a convention
   to be carrying.
8. The signal mass has not been tested for the confound it was introduced to avoid. Cell-count mass
   tracks sample preparation; whether signal-strength mass tracks sequencing depth instead has not
   been measured, and should be before `rank_signal` is taken further.

## 17. Provenance

| Claim | Script | Output |
|---|---|---|
| Paired roster, edge occupancy, variance decomposition, spectrum diagnostic | `08_scoring/13_gw_blindness.py` (D4) | `gw_blindness_D4.csv` |
| Split-half reliability, noise floor | `05_ccc/06_split_half_reliability.R` + `06_split_half.sbatch` + `13_gw_blindness.py` (D5) | `gw_blindness_D5.csv` |
| Coupling geometry, Frobenius baseline | `08_scoring/13_gw_blindness.py` (D1, D3) | `gw_blindness_D1.csv` |
| Per-edge paired test, sign-flip permutation | `08_scoring/14_paired_edge_test.py` | `paired_edge_test.csv` |
| alpha sweeps, all four registered arms x both mass modes | `06_alpha_sweep.py --distance` | `alpha_sweep__{arm}_{mass}.csv` |
| alpha sweeps, the 2x2 | `06_alpha_sweep.py --geometry {walk,walk_ncells,rank_signal}` | `alpha_sweep__{geometry}.csv` |
| Planted-effect power, registered arms | `10_planted_effect_power.py --distance` | `planted_effect_power__{arm}.csv` |
| Planted-effect power, the 2x2 | `10_planted_effect_power.py --geometry {walk,walk_ncells,rank_signal}` | `planted_effect_power__{geometry}.csv` |
| Registered gates | `11_paired_gate.py` | `paired_gate.csv` |
| Gates, the 2x2 | `11_paired_gate.py --arms {walk,walk_ncells,rank_signal}` | `paired_gate__{arm}.csv` |
| Paired planting | `12_paired_planting.py` | `paired_planting.csv` |
| Raw vs rank weight | `08_raw_vs_rank.py` | `raw_vs_rank.csv` |
| Walk geometry itself | `06_distance/02_graph_geometry.py`, `config/graph_geometry.py` | `graph_geometry.csv` |
| Registered arms and gates | `PREREGISTRATION_paired_gate.md`, `config/distance_variants.py` | - |

**Self-checks that ran and passed.** The rank transform reproduces the stored C to 4.7e-16 wherever
it is rebuilt. The Python re-aggregation of production LR tensors reproduces `edge_distance.csv` to
4.7e-16, so half-sample tensors are on the same footing. The planted delta=0 control is null in
every configuration. The sign-flip permutation scores observed and null with the same hand-rolled
statistic, so tie handling cannot bias the comparison. The alpha=0 sweep row groups by mass and not
by cost, as the objective requires. Extracting the walk kernel into `config/graph_geometry.py` left
`graph_geometry.csv` byte-identical. After every change above, the default paths were re-run and
reproduce the frozen results bit for bit: the alpha sweep to 0.000e+00 and all 40 rows of
`paired_gate.csv` to 0.000e+00.
