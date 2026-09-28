# Findings: where the project stands

Last result produced **2026-09-24** (`12_pseudobulk_de/06_hox_compartment.R`).

**What this document is.** A synthesis across two lines of work that each have their own detailed
findings document. It exists at `scripts/` root rather than inside a stage because it spans stages;
it adds no analysis of its own and replaces neither of the two documents it points into:

| line | where the reasoning lives |
|---|---|
| communication-graph topology / OT | `scripts/08_scoring/FINDINGS_topology_null.md` (1022 lines) |
| pseudobulk expression / HOX axis | `scripts/12_pseudobulk_de/PREREGISTRATION_{pseudobulk_de,hox_axis,hox_compartment}.md` |

**Standard.** Every number below names the file that produces it. **One item (N5) does not meet
that standard and is marked so** — its output was written to a scratchpad that has since been
cleared. It is listed as negative because that is what was observed, and as not citable because it
cannot currently be reproduced from disk.

---

## 1. The whole thing in one table

### Six negatives — all pre-registered or proven, none of them a failure to measure

| | claim | headline number | file |
|---|---|---|---|
| **N1** | Cell composition matches or beats every communication distance tried | paired Wilcoxon p **0.23–0.56**; composition alone is the strongest GATE 1 arm | `08_scoring/17_scaccordion_gate1.py` |
| **N2** | GW-family distances are relabelling-invariant by construction; re-matching absorbs label-anchored signal | matching BY LABEL costs median **+44%** over the optimum, strictly worse in **138/138**; FGW omnibus p **0.555** on an effect per-edge tests recover 6/6 | `08_scoring/18_solver_controls.py`, `13_gw_blindness.py` |
| **N3** | scACCorDiON's hitting-time ground cost is inert — on their own published data too | ρ(d, total variation) = **0.9968** ours, **0.99954** on their Peng PDAC cohort | `06_distance/03_scaccordion_distance.py` |
| **N4** | No working cell-level malignancy call exists in this project | inferCNV vs clinical blast %, n=59: ρ = **−0.069**, p = 0.605 | `05_ccc/05_undercall_contamination.R` |
| **N5** | A reference-free "distance to healthy" also fails against blast % | ρ ≈ 0 — **NOT CITABLE, see §2.5** | *(output cleared)* |
| **N6** | The Mono_DC expression signal is reducible to the blast compartment | Mono_DC adjusted for LMPP_GMP: AUC **0.490**, perm p = 0.549 | `12_pseudobulk_de/06_hox_compartment.R` |

### Four positives — two of them are what make the negatives mean something

| | claim | headline number | file |
|---|---|---|---|
| **P1** | The pseudobulk DE pipeline is a working measurement that replicates across platform | **163 / 211 = 77.3%** sign concordance, one-sided binomial p = **3.8e-16** | `12_pseudobulk_de/04_validation.R` |
| **P2** | A patient-genotype classifier built with **zero** malignancy annotation | LMPP_GMP HOX score, NPM1 status: raw AUC **0.975** (n=36, 9 mut) | `12_pseudobulk_de/05_hox_axis.R` |
| **P3** | The distance nulls are not detection-limit nulls | between-sample distances at **3.9–5.6×** measurement noise, within-patient at **2.9–4.2×** | `08_scoring/16_scaccordion_noise_floor.py` |
| **P4** | The **primary-tier** DE nulls are not detection-limit nulls | **LFC80 = 1.5** in all three primary bins, 87.5–96.0% recovery of a planted \|logFC\| 1.5 — but HSC_MPP at n=8 reaches only 46%, see §3.4 | `12_pseudobulk_de/03_permute.R` |

---

## 2. The six negatives

### 2.1 N1 — composition matches or beats communication, four independent times

The comparison ran four times on different formulations and returned the same answer each time
(`FINDINGS_topology_null.md` §A.2, §A.3, §B, §B.7):

- **GATE 1, the within-patient identity test** (§B.7). Passes for Dx-to-Relapse in all six
  representation arms; a 1000-draw permutation of the whole table confirms it (6 of 12 cells
  against a null mean of 0.55, p < 0.001). **But the strongest arm is seven cell-type proportions
  with no communication in them.** No transport arm improves on it — paired Wilcoxon of per-patient
  percentiles, 11 patients, **p = 0.23–0.56**, and each CCC arm beats composition in only 3–4 of 11.
- **Dropping one dataset removes every communication arm and composition survives.** Restricted to
  GSE227903 alone (9 of the 11 patients), the communication arms go and `prop7` stays.
- **The scACCorDiON ablation ladder** (§B): `prop7` (composition only) scores 0.1118 / **0.0569** /
  0.498 — at or better than every distance rung above it.

**The n, stated where it belongs:** the identity readout has **11 paired patients, 9 of them from
GSE227903**. At that n, "p = 0.23–0.56" is absence of evidence; the honest form is an equivalence
bound (`08_scoring/19_gate1_equivalence.py`, 10,000 bootstrap draws): for every communication
distance arm on Dx→Relapse, **median improvements over composition larger than +4.3…+4.5
percentile points are excluded at 95%**. The two node-mass arms are the exception and are reported
as unresolved: observed medians +8.7/+4.3 with upper bounds 26–33 points.

**The fine-bin rerun does not change the answer, and its limits are stated (2026-09-24).**
GATE 1 rerun on the 22-type distances (`17 --suffix=bmm_broad`): every arm weakens — the graphs
are sparser and noisier at fine bins — and only two cells pass on Dx→Relapse: composition
(p = 0.0254) and `corrot` (p = 0.0161); the family-level table permutation gives 2/6 relapse
cells vs null mean 0.29 (p = 0.029). **No arm demonstrably beats composition at 22 types
either**: every paired bootstrap CI spans zero (`corrot`'s observed median improvement is +4.5
percentile points, CI −17.4…+33.3). But honesty cuts both ways: at this n the fine-bin
equivalence bounds are wide (**33 points**, against 4.5 at 7 bins), so improvements at fine
resolution are *undemonstrated*, not *excluded*. `corrot`-at-fine-bins joins the node-mass arms
on the unresolved-leads list.

**Not built in by construction:** CellChat ran with `population.size = FALSE`
(`config_ccc.R:154`, a recorded decision), so edge strengths are not scaled by cell-type
abundance; the rank arm additionally removes magnitude.

**What it licenses:** what is reproducible about a patient is their cell composition, and the
communication graph adds nothing measurable on top of it — with the exclusion bound above as the
quantitative form.
**What it does not:** that cell-cell communication carries no information in AML marrow generally.
It says these distances, on these 7 bins, on these 138 samples, do not beat counting cells.

### 2.2 N2 — the statistic is invariant to the thing the signal was in, by construction

**Correction, 2026-09-24.** An earlier version of this section (and of
`FINDINGS_topology_null.md` §7) narrated this as "the coupling collapses to the product coupling,
where the moment identity applies". That narration was wrong twice: the converged optimum costs
**56% less** than the product coupling, so the optimum is not the product and the moment identity
(which holds only at that point) cannot carry the claim; and diagonal mass near 1/7 says "does not
favour identity", not "is near p ⊗ q". The moment identity is demoted to a remark. The defensible
chain, each link measured (`08_scoring/solver_controls_{A,B,C}.csv`):

1. **By construction.** GW minimises over all couplings, so it is invariant to node relabelling.
   This is definitional, not empirical.
2. **What that costs here, measured.** The label-respecting coupling (maximal diagonal mass) is
   strictly worse than the label-ignoring optimum in **138 of 138 samples**, by **median +44%
   (IQR +23–77%)**; in 43/138 it is worse than even the product coupling. The optimiser always
   finds a re-matching that makes two marrows look more alike than aligning cell types does.
3. **What gets absorbed.** An effect planted into 6 of 49 labelled edges moves 97.7% of targets and
   is recovered **6/6 by per-edge tests**, while the **FGW omnibus reads p = 0.555 / 0.871**.
4. **Not a solver artifact.** Every solve in this project is POT's exact conditional gradient
   (no entropic term, no ε anywhere). One-hot features pull the coupling to **100% of the
   attainable diagonal cap** at α = 0.1 *and* at the production α = 0.5; seven initialisations per
   sample (product, max-diagonal, 5 random vertices) agree in final cost to median relative spread
   **1.3e-16** (4/138 samples above 1e-4, max 1.3e-2); started **on** the diagonal, the solver
   walks off it.

Two corollaries worth keeping: with one-hot features and the uniform barycentre marginal, the
α→0 end of the production distance **equals TV(composition, uniform)** exactly (checked to
3.9e-16) — once identity is pinned, the feature term literally is a composition distance. And the
value-vs-arrangement split stands: sorted-spectrum SD over fixed-position SD = **2.40** for the
production rank transform, and the `logs` arm — the one preserving the most value spectrum — is
the only arm passing GATE 1 at α=1 (p = 0.0015): GW reads value spectra, not arrangements.

**Synthetic sweeps (2026-09-24, `08_scoring/20_synthetic_sweep.py`, design iterations logged in
the script header).** Three results and one disclosed miss:

1. **The inertness surface (C3).** Over K ∈ {5..30} node types × density ∈ {0.35..0.95}, the
   registered inertness rule fires monotonically earlier as either axis grows; hitting-time
   geometry exists only on small sparse graphs (K=5, d=0.35: max/min 44.7, ρ_TV 0.67) and is gone
   in the regime CCC graphs occupy (K≥22: ρ_TV ≥ 0.996 everywhere). The K=7 anchor (ρ 0.991)
   matches production (0.9968).
2. **The invariance, measured directly.** median |GW(C, πC′πᵀ) − GW(C, C′)| = **2.8–3.5e-17** at
   every K — node-relabelled group differences are *exactly* invisible to GW while per-edge tests
   at fixed labels recover 91–99% of moved positions. This is the definitional dichotomy as a
   measurement.
3. **GW is not spectrum-blind — and the real-data blindness is a regime phenomenon.** Under the
   rank transform every sample's value multiset is identical by construction, yet GW separates
   additive-planted groups (p = 0.001) when between-sample heterogeneity is low; its detection
   erodes monotonically as heterogeneity h grows (25/30 → 21/30 → 16/30 cells at h = 0/0.5/1.0).
4. **Disclosed miss:** the sub-prediction that a per-edge-detects/GW-misses *gap* would open
   monotonically with h was NOT confirmed at this grid (gap cells 2/1/1) — per-edge power erodes
   with h too at small amplitudes. The real-data blindness demonstration therefore remains the
   real-data planted experiment (per-edge 6/6, omnibus p 0.55), with the sweep supplying the
   invariance and the direction of the heterogeneity effect, not a full regime map.

**What it licenses:** on a labelled space where node identity is meaningful — which a fixed
cell-type panel is — a GW-family distance cannot use the labels, and what the labels anchored is
absorbed by re-matching at measurably lower cost.
**What it does not:** any claim that FGW is a bad method, or that this is surprising to an OT
audience; the contribution is quantifying what the invariance discards on real data.

### 2.3 N3 — the published method's ground metric is inert on its own data

A linear OT whose cost is constant off the diagonal is exactly proportional to total variation, so
the question is how much spread the cost has.

| arm | off-diag CV | max/min | ρ(d, total variation) |
|---|---|---|---|
| `dwot_shipped` | 0.086 | 1.571 | **0.9968** |
| `dwot_fixed` | 0.407 | 3.838 | 0.9122 |
| `corrot` | 0.294 | 17.183 | 0.9565 |

The registered inertness rule (`max/min < 2.0` **or** ρ > 0.98 ⇒ may not be reported as an
optimal-transport result) fires for the shipped configuration.

**This is not a property of our 7 bins.** Re-running the authors' own published Peng PDAC cohort
(35 samples, 10 cell types, 80 line-graph slots, shipped with their repo) gives ρ = **0.99954** and
cost max/min 1.494 — *worse* than ours. The cause is structural: the shared topology graph is a
near-complete line graph with near-uniform out-degree, so the stationary distribution is
near-uniform and −log of the hitting probability is near-constant. **More cell types makes this
worse, not better.**

**And that prediction is now confirmed on our own data (2026-09-24).** The full CellChat + ladder
rerun at the fine 22-type `bmm_broad` vocabulary (138/138 jobs, 0 failures; 484 line-graph slots,
387 after the variance filter) gives, for the shipped construction: off-diagonal CV **0.0234**,
max/min **1.320**, ρ(d, TV) = **0.99972** — each strictly more inert than at 7 types (0.086 /
1.571 / 0.9968). The node-order misalignment also generalises: **372 of 387** positions at 22
types (33 of 39 at 7). The "finer bins would fix it" escape is closed empirically, not just
structurally. (`06_distance/03 --suffix=bmm_broad` → `scaccordion_qc__bmm_broad.csv`.)

### 2.4 N4 — there is no working malignancy call, and the gap is not small

inferCNV `malignant_frac` against clinical blast %, 59 samples with both:

- **pooled ρ = −0.069, p = 0.605**
- per dataset, mixed signs and none significant: GSE116256 +0.393 (n=13), GSE185381 +0.125 (n=18),
  GSE227903 −0.293 (n=9), GSE289435 +0.182 (n=11), Petti2019 −0.500 (n=5)
- the van Galen signature caller on the same ruler: **ρ = −0.165**
- **median `malignant_frac` 0.105 against median clinical blast fraction 0.650** — a systematic
  ~6× under-call, not just a ranking failure
- **23 of 59 samples are blind**: `malignant_frac` < 5% while clinical blast ≥ 20%

**Consequence that propagates everywhere else in this document:** no analysis in this project can
separate malignant from non-malignant cells. `PREREGISTRATION_pseudobulk_de.md` §10 states this as
a design decision rather than working around it, and N6 exists because of it.

### 2.5 N5 — reference-free distance-to-healthy: negative, and NOT CITABLE

A leave-one-dataset-out k-NN distance from each cell to the nearest healthy cells, computed both
pooled and within hierarchy bin, scored against clinical blast %. **Observed outcome: ρ ≈ 0 in both
forms**, i.e. the same failure as N4 without needing an annotation.

**This item does not meet this document's standard.** The probe ran in a scratchpad
(`probe_far_from_healthy.py`) and wrote `far_from_healthy_per_sample.csv` there; that directory has
since been cleared, so no number here is reproducible from disk. Two further caveats were recorded
at the time and survive: the probe operated **only in the 24-signature-score space**, which is
precisely the representation the collaborator's objection called inadequate; and an apparent
inversion ("healthy further from healthy than AML", AUC 0.279) was traced to a single dataset
(GSE253355, `Reference-scaffold+AuxStroma`) and vanished on its exclusion (AUC 0.477, p = 0.71).

**Open item O1** in §5 is to restore the script into the repo and re-run it so this row becomes
citable or is withdrawn.

### 2.6 N6 — the Mono_DC signal is real, replicates, and is a shadow of the blast compartment

This is the newest result and the most consequential, because it closes the only line that had
produced a positive.

**The signal is real.** Pseudobulk AML-vs-healthy DE in the Mono_DC bin, pre-registered before any
p-value existed: 537 hits at baseline, **385 after the mandatory depth-sensitivity arm**,
permutation p = **0.002** (`sample_within_dataset`) and 0.001 (`library_arm_pure`). On the held-out
Validation arm — GSE116256, Seq-Well, a different platform from Discovery's 10x — **163 of 211
testable hits keep their sign, 77.3%, one-sided binomial p = 3.8e-16**.

**What the genes are.** Up in AML: `HOXB-AS3` +6.6, `HOXA3` +4.6, `HOXB5` +4.4, `HOXB3` +4.3,
`HOXB6` +4.3, `MEIS1` +4.2, `HOXA9` +4.0, `HOXA10` +3.4, `HOXB4` +3.3, `HOXB7` +3.0, plus `MYB`
+1.9 and `FLT3` +1.7. Down in AML: `FCER1A` −3.4, `CLEC10A` −3.3, `CLEC4D` −3.4, `CLEC4E` −3.3,
`FCGR3A` −2.8, `FPR1`, `FPR2`, `G0S2`, `RETN`, `RBP7` −4.9. The up side is the leukemic HOX/MEIS1
program; the down side is mature monocyte and cDC2 identity.

**It is not blast burden.** `PREREGISTRATION_hox_axis.md` Q1 returned **INDETERMINATE** under its
own pre-registered rule: GSE185381 ρ = **−0.434** (n=17) and GSE116256 ρ = **+0.500** (n=7) —
opposite signs. On the same samples inferCNV gives +0.142 and +0.036.

**But it is not its own compartment either.** `PREREGISTRATION_hox_compartment.md`, 36 samples /
9 NPM1-mutants, GSE185381 + GSE116256, Chen2023 excluded for having zero within-dataset contrast
(5/5 mutant):

| | AUC for NPM1 | permutation p |
|---|---|---|
| Mono_DC, raw *(reference)* | 0.893 | — |
| LMPP_GMP, raw *(reference)* | **0.975** | — |
| **Mono_DC adjusted for LMPP_GMP** | **0.490** | **0.549** |
| LMPP_GMP adjusted for Mono_DC | 0.733 | 0.021 |

Adjusted for the blast compartment, Mono_DC lands at **exactly chance**. The reverse direction
survives. Three controls agree: the Set D repeat (12 Discovery-frozen genes instead of 19
locus-defined) gives 0.523 / p 0.36 and 0.753 / **p 0.0075**; the T_NK negative control does not
add (0.429, p 0.71), so the "adding" is not generic; and cell count is not a confound (ρ −0.225
and −0.170 against a 0.40 threshold). Compartment collinearity was ρ 0.822 / 0.857, below the
pre-registered 0.90 stop.

**What it licenses**, at the wording ceiling `PREREGISTRATION_hox_compartment.md` §6 fixed in
advance: the Mono_DC compartment's HOX signal is **not** separable from the LMPP/GMP compartment's.
The 385 genes are not a microenvironment finding.
**What it does not:** it does not establish that the cells carrying the signal are malignant. That
would need the cell-level call N4 says does not exist.

---

## 3. The four positives

### 3.1 P1 — the pseudobulk pipeline is a working measurement

Discovery is 42 AML + 14 healthy (Chen2023 + GSE185381, 10x); Validation is 9 AML + 2 healthy
(GSE116256, Seq-Well), dataset-level and fixed on 2026-08-04 before any of this ran.

- Mono_DC: **163 / 211 = 77.3%** sign concordance, p = **3.8e-16**
- B_Plasma: 11 / 13 = 84.6%, p = 0.0112 — but at **3 AML vs 1 healthy**, which carries no weight
- T_NK: its single frozen hit did not survive the Validation universe; recorded not evaluable
- 174 of 385 and 15 of 28 frozen hits were lost to the Validation universe intersection, counted
  and reported as the pre-registration required

Registered in advance as sign agreement and **not** a powered per-gene replication, so the result
could be neither upgraded if the signs agreed nor downgraded if they did not.

### 3.2 P2 — a patient-genotype classifier that uses no malignancy annotation

- **Held-out and pre-registered:** GSE116256, all 4 NPM1-mutants above all NPM1-wild-types,
  one-sided Wilcoxon p = **0.0079 = 1/C(9,4)**, the minimum attainable at that n
- **Exploratory and larger:** GSE185381, 7 mutants vs 25, AUC **0.891**, p = **0.0004**; the
  mutants' HOX ranks are 1, 4, 5, 6, 7, 11, 13 of 32 against a random-expectation mean rank of 16.5
- **Best single compartment:** LMPP_GMP raw AUC **0.975** (n=36, 9 mutants) — exploratory
- **Not any 19 genes:** of 1000 size- and expression-decile-matched random gene sets, **0.2%**
  separate NPM1 as well (p = 0.002)
- **Two independent gene-set definitions agree:** locus-defined Set L (19 genes) against
  Discovery-frozen Set D (12), Spearman ρ = **0.984**
- Genotype source: `00_curated_manifest.csv` `driver_mutations`, populated for **166 of 244**
  samples — 44 of the scored AML carry a call, 16 of them NPM1-mutant

**This is a validation, not a discovery.** That NPM1-mutant AML expresses HOXA/HOXB/MEIS1 highly
has been textbook for two decades. What is shown is that this pipeline recovers it from raw counts
with no malignancy annotation, across platforms, at near-perfect separation — which is a statement
about the instrument, not about the biology.

Confounds that cannot be removed at this n and are attached to every statement of it: 2 of 4
GSE116256 mutants are also FLT3-ITD; `DNMT3A` co-occurs on both sides; the wild-type group is
heterogeneous (TP53 ×2, KRAS/NRAS, KIT/RAD21, DNMT3A/RUNX1), so the contrast is really
"NPM1-mutant vs everything else".

### 3.3 P3 — the distance nulls are informative nulls

Split-half retesting of 37 samples, five different distances including a published one:

- between-sample distances at **3.9–5.6×** measurement noise
- within-patient distances at **2.9–4.2×**
- Wilcoxon p = **7.3e-12**
- the split is stratified by hierarchy bin, and composition total variation between halves has
  median **0.0013**, max 0.0118 — so the halves are matched on the thing that N1 says carries the
  signal

**This is the load-bearing positive for N1–N3.** Without it, every one of those negatives has the
alternative explanation "the measurement is too noisy to see anything". With it, the measurement
is 3–6× above its own noise and what it reproduces is large idiosyncratic between-patient
variation that is not disease.

### 3.4 P4 — the DE nulls are informative nulls

Planted-effect recovery, 200 planted genes, `12_pseudobulk_de/power.csv`:

| bin | tier | n | recovery at \|logFC\| 1.0 | 1.5 | 2.0 | LFC80 |
|---|---|---|---|---|---|---|
| Mono_DC | primary | 54 | 0.485 | 0.875 | 0.965 | **1.5** |
| T_NK | primary | 53 | 0.505 | 0.905 | 0.995 | **1.5** |
| B_Plasma | primary | 54 | 0.475 | 0.960 | 0.995 | **1.5** |
| LMPP_GMP | secondary | 9 | 0.570 | 0.855 | 0.950 | **1.5** |
| HSC_MPP | secondary | 8 | 0.000 | **0.460** | 0.825 | **2.0** |

**The two tiers must not be quoted at one power.** T_NK's 1 hit is a null at **90.5% power for
\|logFC\| ≥ 1.5** — that is a real null. **HSC_MPP's 0 hits are not**: at n=8 its LFC80 is **2.0**
and it recovers only **46%** of a planted \|logFC\| 1.5, so its zero is weak evidence and is
recorded as such. `power.csv` carries the per-bin classification string; the secondary tier's floor
is stated in `PREREGISTRATION_pseudobulk_de.md` §8.3 and is not the 14-control figure that applies
to the three primary bins.

---

## 4. Why this is a benchmark result and not an inconclusive one

The six negatives would individually each admit the reading "we could not measure it". P3 and P4
remove that reading for all of them: the distance measurements sit 3–6× above their own split-half
noise, and the expression measurements recover 87–90% of a planted |logFC| 1.5. P1 and P2 then
supply what a benchmark needs and a limits paper lacks — **a positive control that works**: the
same pipeline that returns these nulls replicates across platforms at p = 3.8e-16 and recovers a
textbook genotype association at AUC 0.975 with no annotation.

The shape this supports is a statement with both halves measured:

> **When node correspondence is fixed and meaningful, direct distances — and in the limit, counting
> cells — are sufficient; the optimal-transport formulations tested here add nothing on top, and
> §2.2 gives the structural reason why they cannot.**

The second half of the original conjecture — that FGW earns its place once disease drives
structural drift in population state so that cell states no longer correspond one-to-one — is
**not supported on this data**. N6 is the closest this project came to testing it, and the drifted
compartment turned out to be a shadow of the undrifted one.

---

## 5. Open items

| | item | why it is open |
|---|---|---|
| **O1** | Restore `probe_far_from_healthy.py` into the repo and re-run it | N5 is currently not citable (§2.5) |
| ~~O2~~ | ~~The `bmm_broad` re-run~~ **DONE 2026-09-24** | 138/138 jobs clean; N3 confirmed more inert at 22 types (ρ 0.99972); GATE 1: composition still passing, no arm demonstrably better, fine-bin equivalence bounds wide (33 pts) — see §2.1/§2.3 |
| **O3** | LMPP_GMP has 304 frozen hits and **no Validation arm** | pre-registered: GSE116256 gives that bin 1–2 controls (`PREREGISTRATION_pseudobulk_de.md` §5.3). Its permutation null max equals its observed count (431 = 431) |
| **O4** | The P2 numbers outside GSE116256 are exploratory | the GSE185381 ranking was seen before the p-value was computed; only the GSE116256 arm is a registered test |
| **O5** | Whether the LMPP_GMP HOX→NPM1 result can be made a registered test | `driver_mutations` exists for 166/244 samples, which was not known when `PREREGISTRATION_hox_axis.md` was written; a fresh held-out split is possible but has not been designed |

---

## 6. Provenance

| result | script | output |
|---|---|---|
| N1 GATE 1, ablation ladder | `08_scoring/17_scaccordion_gate1.py`, `15_scaccordion_benchmark.py` | `results/tables/08_scoring/` |
| N2 GW blindness, planting | `08_scoring/13_gw_blindness.py` | `gw_blindness_D4.csv` |
| N2 solver controls, label-cost, TV identity | `08_scoring/18_solver_controls.py` | `solver_controls_{A,B,C}.csv` |
| N1 equivalence bounds | `08_scoring/19_gate1_equivalence.py` | `gate1_equivalence.csv` |
| N2/N3 synthetic sweeps | `08_scoring/20_synthetic_sweep.py` | `synthetic_sweep_{inertness,planting}.csv` |
| N3 ground metric, Peng PDAC re-run | `config/scaccordion_geometry.py`, `06_distance/03_scaccordion_distance.py` | `results/tables/06_distance/` |
| N4 malignancy callers | `05_ccc/05_undercall_contamination.R` | `results/tables/05_ccc/undercall_contamination.csv` |
| N5 | *(scratchpad, cleared — see O1)* | — |
| N6, P2 | `12_pseudobulk_de/05_hox_axis.R`, `06_hox_compartment.R` | `hox_axis_*.csv`, `hox_compartment_*.csv` |
| P1 | `12_pseudobulk_de/04_validation.R` | `validation_concordance.csv`, `validation/*.csv` |
| P3 noise floor | `08_scoring/16_scaccordion_noise_floor.py` | `results/tables/08_scoring/` |
| P4 power, permutation | `12_pseudobulk_de/03_permute.R` | `power.csv`, `perm_summary.csv`, `perm_null.csv` |

Pre-registrations, all committed before the code they govern existed:
`08_scoring/PREREGISTRATION_{scaccordion,paired_gate,panel_screen}.md`,
`12_pseudobulk_de/PREREGISTRATION_{pseudobulk_de,hox_axis,hox_compartment}.md`.
