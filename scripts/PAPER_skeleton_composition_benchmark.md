# Paper skeleton — the composition benchmark

Drafted 2026-09-24; **revised same day after external review of the C2 chain** (the product-coupling
narration was wrong and is replaced below; the reviewer's controls are now
`08_scoring/18_solver_controls.py` and `19_gate1_equivalence.py`). Every figure panel names the file
it is built from. **Status**: ✓ on disk · ◐ decision needed · ✗ needs a run.

---

## 1. Abstract, draft 1

> Optimal-transport (OT) distances over cell–cell communication (CCC) graphs are increasingly
> proposed for comparing patients from single-cell RNA-seq. **We ask a bounded question: on
> cell-type-labelled graphs, do patient-level OT/CCC distances carry information beyond cell-type
> composition?** Across three pre-specified readouts — patient identity, AML vs healthy, and
> within-patient treatment change — on 13 public AML bone-marrow datasets (244 curated samples,
> 138 in the CCC cohort), they do not. On the identity readout — **11 paired patients, 9 from one
> dataset** — a seven-number composition vector is the strongest arm, and bootstrap equivalence
> bounds exclude median improvements larger than ~4.5 percentile points for every communication
> distance tested. The reason is structural: Gromov–Wasserstein minimises over all couplings and is
> therefore relabelling-invariant by construction; on our graphs the optimiser always finds a
> label-ignoring matching that makes two graphs look more similar than matching by label does
> (median 44% lower cost, in 138 of 138 samples), so perturbations planted on specific labelled
> edges are absorbed by re-matching — recovered 6/6 by per-edge tests, invisible to the FGW omnibus
> (p = 0.55). The published scACCorDiON construction fares no better: its hitting-time ground
> metric is statistically indistinguishable from total variation on its own published cohort
> (ρ ≥ 0.9995). These nulls are informative, not under-powered: split-half retesting places
> between-sample distances at 3.9–5.6× measurement noise, and planted-effect power gives 87–96%
> recovery at |logFC| = 1.5 in the primary pseudobulk tiers. The same pipeline recovers real signal
> where it exists — differential-expression directions replicate across platforms (163/211
> sign-concordant, p = 3.8e-16) and an annotation-free HOX/MEIS1 axis separates NPM1-mutant AML at
> AUC 0.89–0.98 held-out. The claim is not that CCC carries no information; it is that
> **patient-level graph distances on cell-type nodes are composition in disguise**.

## 2. Title candidates

1. **Counting cells suffices: patient-level optimal-transport distances on cell-type-labelled
   communication graphs are composition in disguise in AML bone marrow**
2. Counting cells suffices: a benchmark of OT distances on AML communication graphs
3. Composition in disguise: what OT adds to patient-level single-cell comparison, measured

The AML scope stays in the title: C1 is an empirical claim about this disease; only C2/C3 are
structural. Removing the scope invites the exact review the scope answers.

"Pre-registered" enters the title ONLY after the Zenodo deposit (§8); until then the abstract says
"time-stamped, pre-specified analysis plans".

## 2.5 Framing decisions (external review 2026-09-24, second round)

**Not a strawman, and the receipts go in the Introduction.** "FGW was never appropriate for
labelled graphs" is half-right — and that it is right is the *finding*, not the flaw, because the
field is answering this question by publication: (i) this project's own pre-specified plans chose
FGW as the primary method after explicit deliberation; (ii) scACCorDiON ships OT distances on CCC
graphs in Bioinformatics (2025); (iii) current reviews list OT-based comparison as a category of
CCC methodology. This paper is the first to test the premise directly.

**The paper must not rest on FGW — and does not.** The correct reviewer counter is "on a labelled
graph the right tool is an identity-pinned per-edge distance; you never tested it." It is tested:
that is the `tabular` rung (per-edge Euclidean, no transport) and the `tv` rung, and on the
identity readout both lose to composition (observed median −4.2 points) with the same equivalence
bound as every OT arm (rule out > 4.5 points). The load-bearing sentence is therefore:

> **From the most flexible OT (FGW, re-matching allowed) to the plainest identity-pinned per-edge
> distance, every rung of the ladder fails to exceed cell-type composition. Label-invariance
> explains only why the OT rungs could not have won; the per-edge rungs lose without it.**

**Label-invariance is written as a trade-off, not a defect.** When node correspondence is unknown
(cross-species or multi-omics alignment), invariance is the point; when correspondence is fixed
and meaningful, it is a cost — and the cost is now quantified (median +44% for matching by label,
138/138). The one-hot control closes the circle: pin identity in the features and FGW's feature
end collapses to TV(composition) exactly (3.9e-16); pin the coupling on the structure end and it
collapses to the per-edge comparison — and both of those are rungs already in the ladder. The
precedent for this framing is Errica et al.: not "GNNs are useless" but "structure is not being
used on these datasets."

**One concession made before a reviewer makes it:** AML marrow is a best case for composition —
blast expansion shifts it massively — so "composition wins" is least surprising exactly here.
That is what §8's path C answers (a low-composition-shift, non-haematological cohort from the
authors' own benchmark).

---

## 3. The claims — each with its evidence

| claim | evidence | source | status |
|---|---|---|---|
| **C1** Composition matches or beats every tested communication distance; improvements beyond ~4.5 percentile points are excluded | GATE 1 strongest arm = 7 proportions; paired p 0.23–0.56; **equivalence bounds +4.3…+4.5 pts (Dx→Relapse, all 5 distance arms)**; survives dataset drop | `17_scaccordion_gate1.py`, `19_gate1_equivalence.py` | ✓ |
| **C2** GW-family distances are relabelling-invariant **by construction**, and on these graphs re-matching absorbs label-anchored signal | label-respecting coupling costs **median +44% (IQR +23..+77%)** over the optimum, strictly worse in **138/138**; planted 6/49 edges: per-edge 6/6, omnibus p 0.55/0.87; solver ruled out (exact CG, one-hot control 100% of cap, 7-init spread 1e-16) | `18_solver_controls.py`, `13_gw_blindness.py` | ✓ |
| **C3** The published hitting-time ground metric ≈ total variation, incl. the authors' own data — and it worsens with resolution, as predicted | ρ(d,TV) 0.9968 at 7 types; **0.99954** Peng PDAC (10); **0.99972 at 22 types** (max/min 1.57 → 1.32) | `03_scaccordion_distance.py` (+ `--suffix=bmm_broad`) | ✓ |
| **C4** The nulls are informative | split-half SNR 3.9–5.6×; DE LFC80 = 1.5, recovery 87–96% primary tier (HSC_MPP 46% labelled weak) | `16_…noise_floor.py`, `power.csv` | ✓ |
| **C5** The same pipeline detects real signal | 77.3% cross-platform sign concordance p 3.8e-16; HOX/NPM1 AUC 0.891/0.975, held-out floor-p; random-set null 0.2% | `04_validation.R`, `05_hox_axis.R` | ✓ |
| **C6** Worked example: a replicating compartment program reducible to another compartment | Mono_DC adj. LMPP_GMP AUC 0.490 vs reverse 0.733 (p 0.021) | `06_hox_compartment.R` | ✓ |
| (C7) No working malignancy call — motivates the annotation-free design | inferCNV ρ −0.069 (n=59); 6× under-call; 23/59 blind | `05_undercall_contamination.R` | ✓ |

**Demoted to a remark** (Methods or supplement): the product-coupling moment identity. It holds
only at p⊗q, and the optimum is measurably elsewhere (product coupling costs +56% over optimum).
It no longer carries C2.

**Disclosed as open leads, not buried**: (i) the two node-mass arms (`node_feat`,
`node_feat150`): positive observed medians (+8.7 / +4.3 pts), equivalence bounds up to 33 pts;
(ii) `corrot` at 22 types: best fine-bin point estimates (median pct 0.091, p = 0.016, paired
median +4.5 pts) but CI spans zero. Improvement over composition can NOT be excluded for either.
`gate1_equivalence{,__bmm_broad}.csv` report both as unresolved.

---

## 4. Figures

### F1 — Cohort, design, and the question (scope panel first)
(a) the three pre-specified readouts and where each is answered — the "wrong question" defence,
promoted to Figure 1 ✓ (b) cohort/split schematic ✓ (c) plan-timestamp timeline (git; Zenodo DOI
when deposited ◐) ✓

### F2 — The benchmark (C1)
(a) GATE 1 percentile heatmap, 11 patients × arms, **with per-dataset annotation (9 GSE227903 / 2
GSE201966)** ✓ (b) **equivalence bounds per arm** — observed median + bootstrap CI, the "rule out
> X pts" panel · `gate1_equivalence.csv` ✓ (c) **GSE227903-only: composition survives, every
communication arm does not — promoted from robustness to a main panel** ✓ (d) ablation ladder ✓ (e) the 22-type rerun: inertness worsens, composition still passes, bounds widen · `scaccordion_{qc,gate1}__bmm_broad.csv`, `gate1_equivalence__bmm_broad.csv` ✓

### F3 — Why: relabelling-invariance, made visible (C2)
(a) per-sample cost of the label-respecting coupling vs the optimum (the +44% panel; product
coupling shown as a reference line, +56%) · `solver_controls_C.csv` ✓
(b) planted-effect: per-edge 6/6 vs omnibus p by amplitude · `gw_blindness` ✓
(c) **how the optimal coupling re-arranges when an effect is planted on labelled edges — the
mechanism picture** ✗ (one run: planting harness + coupling readout, ~half day)
(d) ground-metric inertness incl. Peng PDAC (C3) ✓
(e) solver controls as an inset or supplement: one-hot snap (100% of cap), 7-init cost spread
1e-16, TV identity check 3.9e-16 · `solver_controls_{A,B}.csv` ✓

### F4 — The nulls are informative (C4) — unchanged
(a) split-half floor ✓ (b) power curves, tiers separated ✓ (c) permutation nulls ✓

### F5 — Positive controls and the worked example (C5+C6) — unchanged
(a) cross-platform concordance ✓ (b) HOX×NPM1 scores ✓ (c) random-set null ✓ (d) compartment
adjustment ✓

Supplementary adds: **S7 solver controls (full)**, **S8 α→0 identity**, **S9 synthetic sweeps** (inertness surface; exact-invariance dichotomy; heterogeneity erosion — with the gap-vs-h miss stated): with one-hot features and
the uniform barycentre marginal, the feature end of the production distance **equals TV(composition,
uniform) exactly** (verified 3.9e-16) — once identity is pinned, HDS's feature term literally is a
composition distance ✓.

---

## 5. Methods items that pre-empt review (new)

- **`population.size = FALSE`** (`config_ccc.R:154`, a recorded decision): CellChat edge strengths
  are NOT scaled by cell-type abundance, so "composition wins" is not built into the edge weights
  by construction. The rank arm (magnitude removed) is reported as the additional control. The
  open question this setting leaves — whether the `logs` arm's α=1 GATE-1 pass reflects
  composition's shadow or patient-specific mean expression — is stated, not resolved.
- **The scACCorDiON node-order finding is handled as an upstream report, not a selling point**: a
  GitHub issue is filed before submission; Methods carries one neutral sentence ("we noted a
  discrepancy between X and the paper's description, reported in issue #N; our reimplementation
  follows the paper") and S1 documents the reimplementation.
- The GSE185381 NPM1 peek is disclosed verbatim in Methods (registered arm = GSE116256 only).

## 6. External validation of the HOX score (unchanged, decision ◐)

BeatAML / TCGA-LAML bulk RNA-seq, NPM1 status, n in the hundreds; freeze Set L verbatim,
pre-register one page, one script. Converts C5's strongest number into a registered test on data
this project never touched.

## 7. The reviewer table — updated

| objection | answer |
|---|---|
| **"CCC methods aren't for re-identifying patients — wrong question"** | Introduction ends on this: three pre-specified readouts (identity, disease, within-patient change); OT/CCC exceeds composition on none; meanwhile the same tensors are split-half reproducible and the same pipeline detects DE and genotype signal. Claim scoped to "composition in disguise", not "no information" |
| "n=11, 9 from one dataset" | disclosed in the abstract; equivalence bounds instead of bare non-significance; GSE227903-only panel promoted |
| "Coupling collapse is a solver artifact (ε / CG init)" | no entropic call anywhere (exact CG); one-hot control reaches 100% of attainable diagonal; 7 inits agree to 1e-16; label-init walks away from the diagonal — the optimum, not the init, is label-ignoring |
| "CCC = composition by construction" | `population.size = FALSE`; rank arm as magnitude-free control |
| "7 bins too coarse" | answered empirically: the full 22-type rerun makes the ground metric MORE inert (ρ 0.99972), composition still passes GATE 1, and no communication arm demonstrably beats it (all paired CIs span 0). Disclosed limit: fine-bin equivalence bounds are wide (33 pts vs 4.5 at 7 bins), so fine-resolution improvements are undemonstrated, not excluded; `corrot`-at-22-types is listed as an unresolved lead |
| "GW label-invariance is obvious" | written as a trade-off, not a discovery (§2.5): an advantage when correspondence is unknown, a quantified cost (+44%, 138/138) when it is fixed — the Errica-style framing |
| "On labelled graphs the right method is an identity-pinned per-edge distance — you never tested it" | tested: the `tabular` and `tv` rungs; both lose to composition (obs median −4.2 pts, rule out > 4.5 pts) — same bound as the OT arms |
| "AML is composition's best case — of course it wins" | conceded up front (§2.5); answered by §8 path C: add the composition baseline to a mild-composition-shift cohort from the authors' own benchmark |
| "Just a limits report" | C4 + C5 |

## 8. Gaps before submission — revised order (external review 2026-09-24)

| # | item | size | status |
|---|---|---|---|
| 1 | C2 narrative fixed; label-cost number added; TV identity added | done today | ✓ (F3c panel still ✗) |
| 2 | GATE 1 equivalence bounds | done today | ✓ |
| 3 | bmm_broad re-run | **DONE 2026-09-24**: 138/138 clean; C3 strengthened, C1 unchanged, limits disclosed (see §7 row 5); F2 gains a 22-type panel | ✓ |
| 4 | GitHub issue on the node-order finding — file now, responses take time | 1 h | ◐ draft ready |
| 5 | Zenodo deposit of the six plans + git history → DOI; until then no "pre-registered" in title | 1 h | ◐ needs account |
| 6 | External HOX validation (§6) | 2–3 days | ◐ |
| 7 | Clean-clone replication of the three load-bearing numbers (Peng ρ; label vs optimum cost; GATE 1 paired p) by a second person | half day for them | ◐ |
| 8 | Write per §10 order; bioRxiv at submission | — | after 1–7 |

**Scale-up paths (external review #2), cost-ordered — n=11 with honest bounds reads like a pilot
unless at least one of these lands:**

| path | what | cost | status |
|---|---|---|---|
| **A. simulation sweep** | **DONE 2026-09-24** (`20_synthetic_sweep.py`, 3 design iterations logged): inertness surface monotone in K × density (anchor matches production); relabelling-invariance exact (3e-17) vs per-edge 91–99%; GW not spectrum-blind, detection erodes with heterogeneity 25→16/30. One sub-prediction missed and disclosed (gap-vs-h); the real-data planted experiment stays the blindness demonstration | done | ✓ → S9 |
| **B. AML-vs-healthy as co-primary** | "after regressing out composition, topology adds nothing" on the disease readout (n=72 across 3 dual-arm datasets) — a second, larger-n reading of C1 | needs a designed test + its own pre-specified plan before running; the existing §A grid is close but was not built as an equivalence test | ◐ design first |
| **C. the authors' own benchmark + one baseline** | rerun scACCorDiON's public cohorts with THEIR readout (disease-state clustering) plus one added arm: composition TV. Ahlmann-Eltze style — no new benchmark, one added baseline. Answers "one disease" and "one lab's pipeline" at once; pick one mild-composition-shift, non-haematological cohort to answer §2.5's concession | ~1 week incl. data fetch | ◐ **the decisive one** |

Peng PDAC at ρ ≥ 0.9995 already constitutes one instance of path C.

Dropped from the paper: N5 (distance-to-healthy probe) — not load-bearing, currently not citable.

## 9. Venues — with the fit caveats

| venue | fit | caveat |
|---|---|---|
| **Bioinformatics (OUP)** | the evaluated method's own venue; publishes evaluations | tight page limits — C2's structural analysis likely compressed into supplement; some referees will call label-invariance "obvious" |
| **Genome Biology** | benchmark home; strongest reach | needs the positive-control story front and centre |
| **PLOS Comp Bio** | rigor + negative-result tolerance | slower |
| **GigaScience** | reproducibility-first | lower visibility |
| **Briefings in Bioinformatics** | benchmark track | scope stretch |

## 10. Order of execution (= §8 rows 1–8, then write)

Methods from the six plans → Results tracking C1–C6 → Discussion = FINDINGS §4+§5 expanded.
