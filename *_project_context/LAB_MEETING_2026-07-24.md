# AML Niche-Topology Pipeline — Phase-1 Progress & Stress-Test Report

**Lab meeting · 2026-07-24**
*Extracted from the internal assessment (`ASSESSMENT_2026-07-24.md`). One page ≈ one slide; each slide names the backing figure (`results/figures/11_assessment/g0*`).*

**Thesis in one line:** the pipeline is built end-to-end and we ran the falsification tests most groups skip. They gave us one clean, label-independent negative (H1) and — before any niche claim went out — localized a fixable foundation problem in the malignant-labelling layer. The fix is already scoped.

---

## Slide 1 — The question and the design

**Question.** Beyond *how many* blasts a patient has, does the *communication network* between bone-marrow cell types carry AML disease information? Three pre-registered hypotheses:

- **H1** — are there conserved "emergent" communication edges that appear in AML but not healthy marrow?
- **H2** — does network **topology** discriminate disease *beyond* blast composition?
- **H3** — does a treatment-pressure axis (Dx → MRD → Relapse) deepen stemness/primitiveness?

**Design.** 13 public AML + healthy bone-marrow scRNA-seq datasets → consensus malignant labels → projection onto a frozen hematopoiesis hierarchy (7-node niche: HSC_MPP, LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte, T_NK, B_Plasma) → per-sample cell–cell communication graphs → Fused Gromov-Wasserstein barycenter alignment → permutation-based hypothesis tests.

> **Backing figure — g01 (cohort counts).** Samples at each stage: QC-passing → labelled → 148 CCC graphs → FGW model. Not a nested funnel; the platform-controlled H2 test ultimately rests on 3 datasets / 60 samples.

---

## Slide 2 — Status: the machine is built, and it is deliberately self-critical

**Built and running end-to-end.** 148 cell–cell communication graphs; FGW barycenter alignment of healthy vs AML; and a testing layer that is unusually adversarial by design:

- permutation nulls on every edge and block,
- a **pre-registered α-sweep** that dials the test from pure node-features (α=0) to pure topology (α=1),
- a platform / batch control,
- a feature-decomposition diagnostic that asks *which* feature actually carries any signal.

**One phase passes clean — Phase 2 (hierarchy projection).** Per-sample projection onto the frozen reference, full 220/220 coverage, malignant identity kept strictly orthogonal to the projection. The 7-node vocabulary is correctly frozen (stroma excluded by the established premise — it is the highest false-positive bin).

**The point for today:** we did not stop at "we got a result." We built the tests that try to *break* the result — and that is exactly what surfaced everything below.

---

## Slide 3 — Headline result: the network null (H1) is clean and label-independent

**No conserved emergent edge exists.**

- 0 of 49 directed edges pass the permutation test (best single edge HSC_MPP → T_NK, perm p = 0.052; **min BH-q = 0.87**).
- All 5 pre-specified communication blocks share q = 0.93 — far above 0.05.
- The pre-registered primitive→immune axis is null (mean ΔC = −0.009, p = 0.86).

**Why this one is solid and worth stating plainly:** H1's edges come straight from CellChat and **never pass through the malignant labels**. So this null is immune to the labelling problem on the next slide — it is a genuine, publishable-as-is negative, not a borderline miss. *The colored heatmap is the structure of noise; there is no pattern to read into it.*

> **Backing figure — g03 (H1 null).** (A) AML−healthy difference across all 49 edges; (B) the 5 block q-values, all at 0.93.

---

## Slide 4 — What the stress tests localized: the label layer is the foundation crack

The blueprint calls for a **three-way consensus** malignant call (inferCNV + Numbat + VarTrix, malignant if ≥ 2/3 agree). In practice:

- **125 / 130 samples are single-arm inferCNV.** Only 5 (all GSE227903) reach two arms; **0 reach three.**
- Our own negative control (call malignancy in *healthy* marrow, where every positive is by definition false) quantifies the cost **by bin**:

| Niche bin | False-positive rate |
|---|---|
| **HSC_MPP (the LSC bin)** | **≈ 40%** |
| Erythroid / Megakaryocyte / Mono_DC | 23–27% |
| LMPP_GMP | 14% |
| B_Plasma | 7.7% |
| T_NK | 3.6% |

**Why HSC_MPP matters most:** an LSC *is* a malignant high-stemness cell in the HSC/MPP bins — so HSC_MPP is the exact node the whole niche hypothesis leans on, and it is the dirtiest node in the graph. Lymphoid bins are clean, primitive/myeloid bins are not — the signature of inferCNV keying on primitive-cell expression, i.e. a **method artifact, not rare pre-leukemic clones.**

> **Backing figure — g02 (label quality).** (A) false-positive rate per bin, HSC_MPP flagged; (B) evidence-tier composition, 125/130 single-arm.

---

## Slide 5 — How that propagates: the H2 "signal" is composition, not topology

The decisive pre-registered test is the α-sweep. **At pure topology (α = 1) the effect is null** (global p = 0.11, within-dataset p = 0.27). The significance lives entirely at the *feature* end. Decomposing the features (α = 0.5):

| Feature | β | p | Reading |
|---|---|---|---|
| `frac_malignant` alone | 0.173 | <0.001 | **circular** — forced to 0 in healthy by construction, so it re-encodes the AML/healthy label |
| drop `frac_malignant` | 0.064 | 0.001 | residual… |
| cell count alone | 0.050 | 0.008 | …is largely **technical** (raw cell count) |
| stemness alone | −0.015 | 0.56 | the **intended biology — null, and wrong-signed** |

**Conclusion:** whatever discrimination exists is blast fraction (circular) + cell count (technical) — **not graph topology, and not stemness.** This is precisely the artifact the α = 1 test was built to expose, and it traces straight back to the label-quality problem on Slide 4.

> **Backing figure — g04 (H2 decisive).** (A) α-sweep: significance vanishes as the feature term is turned off; (B) feature decomposition.

---

## Slide 6 — H3 is blocked, and we are not over-reading it

The treatment-pressure axis cannot be adjudicated yet:

- The trend is **not monotone** — it is a V-shape (stemness dips at MRD, rebounds at Relapse), and only 2 of 4 stemness signatures even show the V.
- It is **patient-confounded**: 62% of all MRD malignant cells come from a single sample (6323_MRD), whose malignant compartment was labelled almost entirely non-primitive.
- It rests almost entirely on **one dataset** (GSE227903, the only clean Dx/MRD/Relapse triplets — 23 of 130 samples carry any timepoint at all).
- The timepoint and patient-pairing fields are **not yet reconciled against the source papers.**

**We are deliberately not interpreting H3** until the metadata reconciliation is done. The V-shape can be neither confirmed nor denied today.

> **Backing figure — g05 (H3 blocked).** (A) stemness signatures along Dx/MRD/Relapse; (B) within-pair median change. Header flagged as unverified.

---

## Slide 7 — The fix is scoped: a second, cohort-wide malignancy arm

**Highest-leverage single fix = a second independent evidence arm for the malignant call.**

- Numbat can't be it (needs BAM/FASTQ; only 2–3 datasets qualify).
- **CopyKAT / SCEVAN can** — they read the per-sample QC object directly, **no BAM needed**, so they add an independent expression-CNV arm **across the whole cohort**. This is the most realistic approximation of the blueprint's ≥2/3 consensus under our data, and it directly attacks the HSC_MPP false-positive rate.

**Sequencing (agreed):**
1. Reconcile metadata (timepoint / patient-pairing in the affected datasets).
2. Add the CopyKAT/SCEVAN arm → re-run consensus labelling.
3. Re-run everything downstream **once** (folding in one pending allelic-evidence merge at the same time).
4. Only then add the blueprint pieces skipped so far: the reproducibility screen (Phase 3), the two-method CCC robustness gate (Phase 4), feature discipline (Phase 5–6), and the full robustness stage (Phase 8).

**Take-home:** we used the falsification tests most groups skip, found and *localized* a foundation problem before it could produce a false-positive niche claim, and the one label-independent result (H1) is a clean null. The fix is defined and on the calendar.

---

## Backup — verdict table (if asked for detail)

| Phase / hypothesis | Verdict | One-line reason |
|---|---|---|
| Phase 1 — data + consensus malignancy | 🔴 does-not-stand | 125/130 single-arm inferCNV; HSC_MPP FPR ≈ 40% |
| Phase 2 — hierarchy projection | 🟢 stands | clean reference-only projection, orthogonal labels |
| Phase 3 — meta-programs / SRRS | 🟡 partial | programs extracted; reproducibility screen not computed |
| Phase 4 — CCC graphs | 🟡 partial | 148 graphs built; two-method robustness gate absent |
| Phase 5–6 — distance + FGW | 🟡 partial | runs; uses un-amended features + circular healthy frac_malignant |
| Phase 7 — emergent edges + scoring | 🟡 partial | machinery correct; result strongly negative |
| Phase 8 — robustness | ⚫ not-run | no code yet |
| **H1** conserved emergent edges | 🔴 **null** (label-independent) | 0/49 edges, all blocks q = 0.93 |
| **H2** topology carries info | 🔴 **not topological** | α=1 null; signal = composition + cell count |
| **H3** treatment-pressure axis | ⚫ **blocked** | V-shape, single-patient confound, unverified metadata |

**Figure index:** g01 cohort counts · g02 label quality · g03 H1 null · g04 H2 decisive · g05 H3 blocked — all in `results/figures/11_assessment/` (`.pdf` + `.png`).
