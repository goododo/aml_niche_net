# Pre-registration — running scACCorDiON's real construction on this cohort

**Written 2026-09-16, before any label was scored.** Steps 1–5 are enforced in that order by
`15_scaccordion_benchmark.py`, which refuses to run a step until the previous step's output file
exists. Nothing in section 6 or 7 had been computed when this was written.

Method under test: Nagai JS, Maié T, Schaub MT, Costa IG. *scACCorDiON: a clustering approach for
explainable patient level cell–cell communication graph analysis.* Bioinformatics 2025;
41(5):btaf288. Code vendored read-only at `/FAST/gr10634/gaozy/external/scACCorDiON`, upstream
commit `9d8e7d65fa335b2186473de8c07b69b930b96c59` (2026-08-21, "fix np.alltrue").

---

## 1. The question, and why it is worth the run

`FINDINGS_topology_null.md` section A closes AML-vs-healthy topology across eleven configurations
of this project's own FGW pipeline (α=1 p from 0.305 to 0.966). What that cannot settle is whether
the null belongs to **the cohort** or to **the statistic we chose**. scACCorDiON is the cleanest
available instrument for that question: a published optimal-transport method for exactly this data
type that recovered disease labels on seven cohorts.

Both outcomes are useful, which is why this is worth pre-registering rather than exploring:

- **It separates the groups** → the FGW null was our statistic, and the comparison of the two
  becomes the finding.
- **It is also null** → a much stronger negative result, provided the positive control passes.

## 2. Sample set

Fixed before any analysis. **All 138 samples that have a per-sample ligand-receptor tensor on
disk** (`05_ccc/tensors/<dataset>/<sample>__ccc_cellchat.csv`), across 10 datasets.

| | |
|---|---|
| disease | 115 AML / **23 healthy** |
| timepoint (4-level) | Diagnosis 88, Healthy 23, Treatment 14, Relapse 13 |
| datasets | GSE185381 47, GSE227903 26, GSE116256 16, GSE289435 11, GSE239721 10, Chen2023 9, E-MTAB-11536 6, GSE201966 5, Petti2019 5, GSE207356 3 |
| healthy samples come from **only 4** datasets | GSE185381 10, E-MTAB-11536 6, Chen2023 4, GSE116256 3 |

**Deliberate departure from the standing paired-samples rule.** scACCorDiON is a cohort clustering
method scored by ARI; the 37-sample paired roster cannot support it (see §5 for the numbers that
rule it out). The paired set is not abandoned — it returns in §8 as a separate, non-ARI analysis.

**E-MTAB-11536 is 6/6 healthy.** One whole dataset *is* a disease class. This is registered here
because it is the single largest reason a positive result could be batch.

## 3. The ablation ladder (the six arms)

Built by `06_distance/03_scaccordion_distance.py`. Each rung moves exactly one ingredient, so a
separation can be attributed rather than merely observed.

| # | arm | mass | ground cost | transport |
|---|---|---|---|---|
| 1 | `prop7` | 7 bin cell counts, CLR | — | — (Euclidean) |
| 2 | `tabular` | 49→39 edge slots | — | — (Euclidean) |
| 3 | `tv` | same | — | — (½·L₁) |
| 4 | `corrot` | same | correlation between slots | `ot.emd2` |
| 5 | `dwot_shipped` | same | hitting time, **their code verbatim** | `ot.emd2` |
| 6 | `dwot_fixed` | same | hitting time, **as the paper describes** | `ot.emd2` |

Rung 1 is the **mandatory negative control**: it contains no communication at all. If it matches
anything above it, the ladder has no methodological finding and the signal is blast fraction.

Rungs 5 and 6 differ because the shipped code and the published text disagree; all three
disagreements plus one ordering defect are documented in `scripts/config/scaccordion_geometry.py`.

## 4. Disclosed pre-knowledge

Measured **before** this registration, all label-free. Disclosed because it constrains what the
run can claim, and because a reader must be able to tell what was known in advance.

**(a) The published ground cost is nearly inert, and not because of our resolution.**

| arm | off-diag CV | max/min | ρ(d, total variation) |
|---|---|---|---|
| `dwot_shipped` | 0.086 | 1.571 | **0.9968** |
| `dwot_fixed` | 0.407 | 3.838 | 0.9122 |
| `corrot` | 0.294 | 17.183 | 0.9565 |
| `tabular` | — | — | 0.9395 |
| `prop7` | — | — | 0.4294 |

The authors' **own** published Peng PDAC cohort (35 samples, 10 cell types, 80 slots) gives
ρ = **0.99954**, cost max/min 1.494 — *worse* than ours. The degeneracy is a property of the
construction (a near-complete line graph has near-uniform stationary distribution), not of our
7-bin resolution. A finer-resolution re-run would not fix it.

**(b) Registered inertness rule.** An arm with `cost_max_over_min < 2.0` **or** `ρ(d, TV) > 0.98`
is declared inert and **may not be reported as an optimal-transport result**; it is reported as
total variation. On the observed values this fires for `dwot_shipped` and not for `dwot_fixed`.

**(c) Their default variance filter removes 10 of 49 slots**, including **all seven Megakaryocyte
edges**. Registered as disclosed, not corrected.

**(d) The shipped implementation misorders its own cost matrix** — 33 of 39 positions differ
between `expgraph` node order and `p.index` order, so `ot.emd2` receives a permuted cost. The
`dwot_fixed` arm asserts alignment; `dwot_shipped` reproduces the defect deliberately.

**(e) The adapter is exact.** Their group-by-sum reproduces production `weight_probsum` to
max |Δ| = 5.3e-15, so these distances are comparable to section A's eleven configurations.

**(f) This distance is largely a detection-pattern statistic.** Pre-flight: ρ(d, support Jaccard)
≈ 0.86, and binarising the weights barely changes it. Hence the mandatory companion readout in §7.

## 5. What chance looks like (simulated, 5000 draws, label structure only)

*The authoritative values are the ones step 1 writes to `08_scoring/scaccordion_null_ari.csv`.
The table below is the pre-run estimate from an independent draw of the same simulation; the two
differ only by Monte Carlo noise (full-cohort disease max-over-k p95: 0.0251 here, 0.0240 there).
Thresholds in §7 are evaluated against the file, not against this table.*

| cell | label | k=2 p95 | max over k∈2..7, p95 |
|---|---|---|---|
| full cohort n=138 | disease | 0.0184 | **0.0251** |
| full cohort n=138 | tp4 | 0.0197 | 0.0313 |
| GSE185381 n=47 | disease | 0.0594 | **0.0818** |

**Batch ceiling** — the best ARI reachable by *any* clustering that knows only `dataset`,
by exhaustive enumeration of all 115,179 coarsenings of 10 datasets into ≤7 blocks:

| label | ceiling | achieved by |
|---|---|---|
| disease | **0.3813** | Chen2023 \| E-MTAB-11536 \| the other eight |
| tp4 | 0.3763 | 5 blocks |

**The entire interpretable window for full-cohort disease ARI is 0.025 → 0.381.** At or above
0.381 the result is indistinguishable from perfect recovery of study of origin.

**Why the paired roster is excluded from ARI.** At n=37 (and n=22 for Dx→Relapse) the max-over-k
null p95 exceeds 0.12, so no ARI below roughly 0.2 would mean anything. ARI is not estimable
there and will not be computed there.

## 6. Locked analysis order

1. **Simulated null** — label structure only, no distance touched.
2. **Positive control** — see §7. Nothing downstream may be interpreted if this fails.
3. **`dataset` label** — a label known to exist. If no arm recovers it, the distances are at the
   noise floor and step 4 cannot separate "no effect" from "no measurement".
4. **`disease` label** — the primary outcome. Last, on purpose.
5. **Support-sparsity companion** — see §7.

## 7. Positive control and decision rules

### 7.1 Positive control

Plant a **multiplicative** effect on the **same six edges** `10_planted_effect_power.py` uses
(`B_Plasma→Erythroid`, `HSC_MPP→HSC_MPP`, `Megakaryocyte→HSC_MPP`, `T_NK→Erythroid`,
`B_Plasma→Mono_DC`, `Megakaryocyte→T_NK`), in a random half of the 138 samples, over
δ ∈ {0, 0.25, 0.5, 1, 2, 4}. Planting is at the **ligand-receptor level**, so the entire published
construction — group-by-sum, variance filter, shared topology graph, hitting-time cost, transport
— is rebuilt from the perturbed input. Multiplicative rather than additive so that the support
pattern is unchanged and the control tests **strength**, not detection.

**Four of the six pinned edges do not survive their default variance filter.** Registered in
advance, with this consequence: the control is run at `filter_q = 0.2` (their default) **and**
`filter_q = 0`.

**Pass rule:** `dwot_fixed` at δ = 4, `filter_q = 0`, must reach **ARI ≥ 0.30** against the planted
split, and ARI must be non-decreasing in δ over the grid up to numerical noise.

- Fails → the run is reported as **inconclusive — the instrument failed its own control**, never
  as evidence that the cohort has no signal. This is section 11's rule applied to somebody else's
  method.
- Passes at `filter_q = 0` but not at 0.2 → a reportable finding in its own right: their default
  filter destroys the planted effect.

### 7.2 Primary outcome

**One** number: ARI of k-medoids (`fasterpam`, their clusterer) on `dwot_fixed`, full cohort,
against **disease**, at k = 2 and max over k ∈ 2..7.

**Declared positive only if all four hold:**

1. the positive control passed;
2. ARI > 0.0251 (the max-over-k null p95);
3. ARI exceeds the same arm's ARI against the support median-split (§7.3);
4. it is not batch — either ARI < 0.3813 *and* the within-dataset cell GSE185381 also exceeds its
   own null (0.0818), or ARI(dataset) is itself at chance.

**Declared null** only if the positive control passed and (2) fails.

### 7.3 Mandatory companion readouts

- ARI of the same clustering against the **support median-split**. An arm whose disease ARI does
  not exceed its support ARI has not shown a communication result.
- ARI against **`dataset`**, reported beside every disease ARI. Never report one without the other.

### 7.4 Multiplicity

Primary is a single pre-specified cell: `dwot_fixed` × full cohort × disease. Everything else
(6 arms × 2 cells × 2 labels) is secondary and BH-FDR corrected within the arm family. Secondary
results may not be promoted to primary after the fact.

## 8. Registered secondary analysis — the paired set

Separate from the ARI programme and obeying the standing paired-samples rule. The 37-sample paired
roster (`08_scoring/paired_roster.csv`) is scored with **this project's registered paired test**,
not ARI: the scACCorDiON distance replaces the FGW distance in the GATE 1 within-patient
comparison, making it directly comparable as a **twelfth row** in FINDINGS section A's table.

## 9. What is *not* being done, and why

- **No `bmm_broad` re-run at finer resolution.** 7 bins is below every cohort scACCorDiON was
  tested on, so a null is arguable on those grounds. But §4(a) shows finer resolution does not
  rescue the geometry, and the re-run costs 138 CellChat jobs. Deferred to an explicit decision
  after the first numbers, not folded in silently.
- **No bug fixes pushed upstream as part of this analysis.** The defects in §4(b,d) are
  reproduced in `dwot_shipped` and corrected in `dwot_fixed`; both are reported.
- **No k chosen by silhouette.** k is either the number of true classes or swept 2..7 with the
  max reported against a max-over-k null. Choosing k by a data-driven criterion and then scoring
  ARI at that k would not have a matching null.

## 10. Provenance

| | |
|---|---|
| construction | `scripts/config/scaccordion_geometry.py` |
| distances | `scripts/06_distance/03_scaccordion_distance.py` |
| scoring | `scripts/08_scoring/15_scaccordion_benchmark.py` |
| vendored method | `/FAST/gr10634/gaozy/external/scACCorDiON` @ `9d8e7d6` |
| extra deps | `/FAST/gr10634/gaozy/external/pylibs` (kmedoids, genieclust, pydiffmap + 2 transitive); the production `general_env` is untouched |
| seed | 491638 (`config_paths.sh`) |
