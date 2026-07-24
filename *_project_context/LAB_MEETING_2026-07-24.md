# AML Niche-Topology Pipeline — Phase-1 Progress & Stress-Test Report

**Lab meeting · 2026-07-24**
*Extracted from the internal assessment (`ASSESSMENT_2026-07-24.md`). One page ≈ one slide; each slide names the backing figure (`results/figures/11_assessment/g0*`). Detailed answers to likely questions are in the "Anticipated Q&A" backup at the end.*

**Thesis in one line:** the pipeline is **preliminarily built** end-to-end and we ran the falsification tests most groups skip. They gave us one label-independent negative (H1) and — before any niche claim went out — localized a fixable foundation problem in the malignant-labelling layer. The fix is already scoped.

---

## Slide 1 — The question and the design

**Question.** Beyond *how many* blasts a patient has, does the *communication network* between bone-marrow cell types carry AML disease information? Three pre-registered hypotheses (with their falsification criteria):

- **H1 — Convergence.** There is a set of **LSC-like → microenvironment communication axes** that reproduce in **≥ 60% of independent studies** and are enriched across multiple genetic subtypes (NPM1, TP53, KMT2A-rearranged, monocytic, …), significantly above the permutation null.
  *Falsified if:* conserved axes do not exceed the null, **or** conservation is fully explained by platform/batch.
- **H2 — Topology > best single pair.** Network-topology deviation (HDS/ATS via FGW/GW distance), under **strict cross-cohort validation**, discriminates AML vs healthy **better than the best single L–R model and better than a composition-only model**.
  *Falsified if:* a single-L–R or composition-only model matches it (topology adds no incremental information).
- **H3 — Monotonic intensification.** **Within paired patients**, topology deviation increases **monotonically Dx < MRD/residual < Relapse**.
  *Falsified if:* no monotone trend, **or** the trend is explained by blast fraction alone.

**Design.** 13 public AML + healthy bone-marrow scRNA-seq datasets → consensus malignant labels → projection onto a frozen hematopoiesis hierarchy (7-node niche: HSC_MPP, LMPP_GMP, Mono_DC, Erythroid, Megakaryocyte, T_NK, B_Plasma) → per-sample cell–cell communication graphs → Fused Gromov-Wasserstein barycenter alignment → permutation-based hypothesis tests.

> **Backing figure — g01 (cohort counts).** Samples at each stage: QC-passing → labelled → 148 CCC graphs → FGW model. Not a nested funnel; the platform-controlled H2 test ultimately rests on 3 datasets / 60 samples.

---

## Slide 2 — Status: the machine is preliminarily built, and it is deliberately self-critical

**Running end-to-end (first pass).** 148 cell–cell communication graphs; FGW barycenter alignment of healthy vs AML; and a testing layer that is unusually adversarial by design:

- permutation nulls on every edge and block,
- a **pre-registered α-sweep** that dials the test from pure node-features (α=0) to pure topology (α=1),
- a platform / batch control,
- a feature-decomposition diagnostic that asks *which* feature actually carries any signal.

**One phase passes clean — Phase 2 (hierarchy projection).** Per-sample projection onto the frozen reference, full 220/220 coverage, malignant identity kept strictly orthogonal to the projection. The 7-node vocabulary is correctly frozen (stroma excluded by the established premise — it is the highest false-positive bin).

**One caveat to state up front:** the communication layer is currently **single-method (CellChat only)**. The blueprint's *primary* engine is LIANA+, with CellChat as the sensitivity control, plus a two-method robustness gate (and we will also want NicheNet). So **every CCC-based result below — including H1 — is provisional on method choice** until the multi-method gate is in.

**The point for today:** we did not stop at "we got a result." We built the tests that try to *break* the result — and that is what surfaced everything below.

---

## Slide 3 — Headline: the network null (H1) is clean *against the label problem* (still single-method)

**No conserved emergent edge exists — in CellChat.**

- **Edge level:** 0 of 49 directed edges pass permutation (best edge HSC_MPP → T_NK, perm p = 0.052; min BH-q = 0.87).
- **Block level:** the 5 pre-specified communication blocks have raw perm p ranging **0.34 → 0.93** and all collapse to q = 0.93. *(This identical q is correct BH-FDR behavior, not a computation error — when even the best block is far from significant, the largest p becomes the shared FDR ceiling. Verified against the code.)*
- The **primary block — the SCS axis, primitive→immune ({HSC_MPP, LMPP_GMP} → {T_NK, B_Plasma, Mono_DC}) — is exactly the hypothesized LSC→microenvironment interaction.** It comes out *nominally* stronger in AML (mean ΔC = −0.009) but **p = 0.86, q = 0.93 — not significant.** A hint in the predicted direction, drowned in noise.

**Why this is our most defensible result today:** H1's edges come straight from CellChat and **never pass through the malignant labels**, so this null is immune to the label-quality problem on the next slide. **But** it is CellChat-only — adding LIANA+/NicheNet is required before calling it a final null. *(The colored heatmap is the structure of noise; there is no pattern to read into it.)*

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

**Why HSC_MPP matters most:** an LSC *is* a malignant high-stemness cell in the HSC/MPP bins — so HSC_MPP is the exact node the whole niche hypothesis leans on, and it is the dirtiest node in the graph. Lymphoid clean, primitive/myeloid dirty — the signature of inferCNV keying on primitive-cell expression, i.e. a **method artifact, not rare pre-leukemic clones.**

**On stroma (FPR 90.3%):** this is *not* an error we ship — stroma is **excluded from the CCC graph** precisely *because* of this rate (the established premise). Its 90% is a **diagnostic** that inferCNV is unreliable on non-hematopoietic/primitive cells; it never enters the analysis. The in-graph problem to fix is HSC_MPP.

> **Backing figure — g02 (label quality).** (A) false-positive rate per bin, HSC_MPP flagged; (B) evidence-tier composition, 125/130 single-arm.

---

## Slide 5 — How that propagates: the H2 "signal" is composition, not topology

FGW mixes two costs by a weight **α**: α weights **topology** (Gromov-Wasserstein structure term), (1−α) weights **node features** (Wasserstein term). So α=1 is pure topology, α=0 is pure node features.

The decisive pre-registered test is the α-sweep. **At pure topology (α = 1) the effect is null** (global p = 0.11, within-dataset p = 0.27). The significance lives entirely at the *feature* end (α=0: p = 0.0009). Decomposing the features (α = 0.5):

| Feature | β | p | Reading |
|---|---|---|---|
| `frac_malignant` alone | 0.173 | <0.001 | **circular** — this is the sample's blast/malignant fraction, but it is *forced to 0* for healthy samples, so it just re-encodes the AML/healthy label |
| drop `frac_malignant` | 0.064 | 0.001 | residual… |
| cell count alone | 0.050 | 0.008 | …is largely **technical** (raw cell count) |
| stemness alone | −0.015 | 0.56 | the **intended biology — null, and wrong-signed** |

**Conclusion:** whatever discrimination exists is blast fraction (circular) + cell count (technical) — **not graph topology, and not stemness.** This is precisely the artifact the α = 1 test was built to expose, and it traces straight back to the label-quality problem on Slide 4.

> **Backing figure — g04 (H2 decisive).** (A) α-sweep: significance vanishes as the feature term is turned off; (B) feature decomposition.

---

## Slide 6 — H3 is blocked; we test neither claim yet

Two separate things must not be conflated:

- **Stemness along the timepoint axis** (transcriptional, measured): **not monotone** — a V-shape (dips at MRD, rebounds at Relapse), and only 2 of 4 stemness signatures even show the V. It is also **patient-confounded**: 62% of all MRD malignant cells come from a single sample (6323_MRD). So this is likely a single-sample effect, not a cohort trend.
- **Topology deviation along the timepoint axis** (the *actual* H3 variable): **never tested** — the pre-registered three-model (A/B/C) competition was not fit; there is no result table. So we cannot say topology deviation is monotone *or* non-monotone — only that it is **not yet computed**.

Both point to the same conclusion — **H3 cannot be adjudicated today** — but for different reasons: stemness is measured-but-confounded, topology is un-run. And it all rests on **one dataset** (GSE227903, the only clean Dx/MRD/Relapse triplets; 23 of 130 samples carry any timepoint), whose timepoint/pairing fields are **not yet reconciled against the source papers.**

> **Backing figure — g05 (H3 blocked).** (A) stemness signatures along Dx/MRD/Relapse; (B) within-pair median change. Header flagged as unverified.

---

## Slide 7 — The fix is scoped: harden the malignant call, then re-run once

**Highest-leverage single fix = a second, cohort-wide evidence arm for the malignant call: CopyKAT / SCEVAN.** They read the per-sample QC object directly (**no BAM needed**), so they can run across the whole cohort. Stated honestly:

- **What it buys (real, and testable):** CopyKAT (Bayesian, self-calibrating baseline) and SCEVAN (reference-free) segment the genome and define "who is normal" *differently* from inferCNV (which depends on an external reference — our weakest link). Agreement across all three reduces **method-specific** false positives. We will *measure* the gain, not assert it: run the same healthy-marrow FPR control on CopyKAT/SCEVAN — **if they also show ~40% HSC_MPP FPR, the confounder is shared/biological and consensus won't fix it; if not, it genuinely helps.**
- **What it does NOT buy (say this before you're asked):** all three read the *same* signal (expression), so they are **not** modality-orthogonal like allelic Numbat or SNV VarTrix; they can share the HSC_MPP failure mode. And **copy-number-neutral clones (e.g. NPM1-mutated normal-karyotype AML) are invisible to all expression-CNV methods.** The arm meant to catch those is **SNV/VarTrix on known-mutation datasets (Petti2019)** — currently an inactive stub (`SNV_ARM_ACTIVE=FALSE`, 0/130 samples), a separate reactivation track.

**Sequencing (agreed):**
1. Reconcile metadata (timepoint / patient-pairing).
2. Add CopyKAT/SCEVAN → re-run consensus labelling (measure the HSC_MPP FPR gain).
3. Re-run everything downstream **once** (folding in one pending allelic-evidence merge).
4. Then add the blueprint pieces skipped so far: multi-method CCC gate (LIANA+/NicheNet, Phase 4), the SRRS reproducibility screen (Phase 3), feature discipline (Phase 5–6), the full robustness stage (Phase 8), and — as its own track — the SNV/VarTrix arm.

**Take-home:** we used the falsification tests most groups skip, found and *localized* a foundation problem before it could produce a false-positive niche claim, and the one label-independent result (H1) is a clean null in CellChat. The fix is defined and on the calendar.

---

## Backup A — verdict table (if asked for detail)

| Phase / hypothesis | Verdict | One-line reason |
|---|---|---|
| Phase 1 — data + consensus malignancy | 🔴 does-not-stand | 125/130 single-arm inferCNV; HSC_MPP FPR ≈ 40% |
| Phase 2 — hierarchy projection | 🟢 stands | clean reference-only projection, orthogonal labels |
| Phase 3 — meta-programs / SRRS | 🟡 partial | programs extracted; reproducibility screen not computed |
| Phase 4 — CCC graphs | 🟡 partial | 148 graphs built; **single-method (CellChat); LIANA+/NicheNet + two-method gate absent** |
| Phase 5–6 — distance + FGW | 🟡 partial | runs; uses un-amended features + circular healthy frac_malignant |
| Phase 7 — emergent edges + scoring | 🟡 partial | machinery correct; result strongly negative |
| Phase 8 — robustness | ⚫ not-run | no code yet |
| **H1** convergence (LSC→niche axes) | 🔴 **null** (label-independent, CellChat-only) | 0/49 edges, all blocks q = 0.93, SCS axis p = 0.86 |
| **H2** topology > single pair | 🔴 **not topological** | α=1 null; signal = composition + cell count |
| **H3** monotonic intensification | ⚫ **blocked** | stemness V-shape (1-patient confound); topology test un-run; metadata unverified |

**Figure index:** g01 cohort counts · g02 label quality · g03 H1 null · g04 H2 decisive · g05 H3 blocked — all in `results/figures/11_assessment/` (`.pdf` + `.png`).

---

## Backup B — Anticipated Q&A (speaker notes)

**Q: All 5 block q-values are 0.93 — did the FDR correction break?**
No. The raw p-values differ (0.34, 0.45, 0.84, 0.86, 0.93). Benjamini-Hochberg is a step-up: `q_(i) = min over j≥i of (p_(j)·m/j)`. Because every block is non-significant, the largest p (0.93, rank 5) gives q = 0.93, and every smaller block's inflated value p·5/j exceeds 1, so the monotonicity rule pulls them all down to 0.93. The shared q just equals the least-significant p — the correct signature of a uniformly null result. n_perm = 10,000; code is a textbook hand-rolled BH, verified line-by-line.

**Q: What exactly is the "primitive→immune" / SCS axis?**
The pre-registered primary block: senders = LSC/primitive {HSC_MPP, LMPP_GMP}, receivers = immune/paracrine {T_NK, B_Plasma, Mono_DC} (6 directed edges). Locked in `config_ccc.R` before results existed. It *is* the hypothesized LSC→microenvironment interaction. Result: nominally stronger in AML (ΔC = −0.009) but p = 0.86 — not significant. The other blocks: immune→primitive (reverse), immune→immune (9 edges), primitive-autocrine (LSC internal, 4 edges), other-lineages (anything touching Erythroid/Megakaryocyte, 24 edges).

**Q: Isn't CopyKAT/SCEVAN just inferCNV again (all expression-based)?**
Algorithmically no — different segmentation (moving-average vs Bayesian-KS vs variational) and, crucially, different ways of defining the normal baseline (inferCNV needs an external reference; the other two self-calibrate). So consensus removes *method-specific* artifacts. But honestly: same modality (expression), so they can share *biological* errors like the HSC_MPP over-call, and adding them keeps the sample at one evidence "type," not the blueprint's cross-modality 2/3. We will test whether the HSC_MPP FPR actually drops by re-running the healthy-marrow FPR control on the new callers — turning the argument into a measurement.

**Q: What about malignant cells with no CNV?**
Correct limitation — NPM1-mutated normal-karyotype AML carries no CNV, so all expression-CNV methods (and Numbat's CNV core) miss it. That is the SNV/VarTrix arm's job (Petti2019 is the natural first dataset), and that arm is currently an inactive stub — a separate reactivation track, not something CopyKAT/SCEVAN solves.

**Q: What is the SRRS reproducibility screen (Phase 3)?**
SRRS = Study-level Reproducibility Score, dual-layer. **Global SRRS** = fraction of *independent studies* in which a meta-program is recovered (recovered = a per-sample program with top-gene Jaccard ≥ τ=0.2); hard gate **≥ 0.6** (recovered in ≥ 60% of studies). **Conditional layer**: within strata of ≥ 3 studies, an evidence weight `w_cond = (n_recovered/n_total)·log(1+n_studies)`, plus HSI (hierarchy-specificity) and healthy-match annotation. Note: the 60%-of-studies bar is the *same* reproducibility standard as H1 — SRRS is essentially H1 operationalized at the meta-program level. None of it is implemented yet; current "robust" programs rest on a hand-curated `mp_labels.tsv`.

**Q: Why is healthy `frac_malignant` = 0?**
A config switch (`FGW_ZERO_HEALTHY_MAL`, `01_build_fgw_inputs.R:67`) sets it deterministically, which is what makes that feature circular for the AML/healthy contrast. It is fixable (drop the switch / impute honestly); it's on the deferred bug list for the re-run.
