# AML Niche-Topology Pipeline — Phase-1 Progress & Stress-Test Report

**Lab meeting · 2026-07-24**
*Extracted from the internal assessment (`ASSESSMENT_2026-07-24.md`), read against the design blueprint **v1.0 + v1.1 patch** (v1.1 dated 2026-07-22). One page ≈ one slide; each slide names the backing figure (`results/figures/11_assessment/g0*`). Detailed answers to likely questions are in the "Anticipated Q&A" backup.*

**Thesis in one line:** the pipeline is **preliminarily built** end-to-end and we ran the falsification tests most groups skip. The current results are the **un-amended v1.0 pipeline**; the design has since (v1.1) already anticipated the interpretation pitfalls these first-pass results reveal. Two things are solid today: a **label-independent H1 null** (in CellChat), and a **localized foundation problem in the malignant-labelling layer** that sits upstream of every downstream claim. The fix is scoped.

---

## Slide 1 — The question and the design (hypotheses as of v1.1)

**Question.** Beyond *how many* blasts a patient has, does the *communication network* between bone-marrow cell types — its **topology and the direction of its deviation** from healthy — carry AML disease information? Three pre-registered, falsifiable hypotheses:

- **H1 — Convergence.** A set of **LSC-like → microenvironment communication axes** reproduces in **≥ 60% of independent studies** and is enriched under **three orthogonal stratifications**: (i) genetic subtype (NPM1/TP53/KMT2A-r/monocytic), (ii) epigenetic subgroup (eCHROMA's 16 ATAC subgroups, RNA-assigned), (iii) transcriptional-state class (our cNMF meta-programs).
  *Falsified if:* conserved axes don't exceed the null, **or** conservation vanishes under *any* stratification, **or** it's explained by platform/study.
- **H2 — Topology (+ direction) > single pair.** Topology deviation discriminates AML vs healthy better than **(a)** the best single L–R model, **(b)** a composition-only model, **and (c)** a node-feature-only model (the α=0 endpoint); *and* the **direction** of deviation (cosθ, orthogonal residual HDS⊥) adds information beyond the scalar distance.
  *Falsified if:* any of (a)/(b)/(c) matches it, **or** direction adds nothing after composition-matching.
- **H3 — Trajectory model (not assumed monotone).** A pre-registered **competition of three trajectory models** for how topology deviation moves along Dx → MRD → Relapse: **A** monotone intensification, **B** rebound/V-shape (MRD trough, HDS⊥ elevated at MRD), **C** hysteresis (relapse enters a *new* topological direction, cosθ_Rel < cosθ_Dx). Fit by a mixed-effects model (patient/study/platform random effects, frac_malignant covariate), compared by AIC/BIC + nested LRT.
  *Falsified if:* none beats the null (timepoint uninformative), **or** the best model's effect vanishes after regressing out blast fraction.

**Design.** 13 datasets → consensus malignant labels → projection onto a frozen 7-node hematopoiesis hierarchy → per-sample cell–cell communication graphs → Fused Gromov-Wasserstein alignment → permutation tests.

> **Backing figure — g01 (cohort counts).** Samples per stage; the platform-controlled test rests on 3 datasets / 60 samples.

---

## Slide 2 — Status: preliminarily built, deliberately self-critical, and the design is self-correcting

**Running end-to-end (first pass).** 148 CellChat graphs; FGW barycenter alignment; and an adversarial testing layer (per-edge/block permutation nulls, a pre-registered α-sweep, a platform control, feature decomposition).

**One phase passes clean — Phase 2 (hierarchy projection):** full 220/220 coverage, malignant identity kept orthogonal to the projection, 7-node vocabulary frozen (stroma excluded — it is the highest-FPR bin).

**Two framing facts for today:**
- The communication layer is **single-method (CellChat only)** so far. The blueprint's *primary* engine is LIANA+, with CellChat as the control, plus a two-method gate and NicheNet. So every CCC-based result — including H1 — is **provisional on method** until that gate is in. *(This is an implementation gap, not a design gap.)*
- The **v1.1 patch (two days ago) already scoped the fixes** for most of what these first-pass results show — it moved H3 to a three-model competition, formalized the α=0 pure-feature control for H2, and mandated honest richer node features (M8). What we're presenting is the pipeline *before* those amendments.

---

## Slide 3 — Headline: the network null (H1) is clean *against the label problem* (still CellChat-only)

**No conserved emergent axis is detectable — in CellChat.**

- **Edge level:** 0/49 directed edges pass permutation (best edge HSC_MPP → T_NK, p = 0.052; min BH-q = 0.87).
- **Block level:** the 5 pre-specified blocks have raw p 0.34 → 0.93, all collapsing to q = 0.93. *(Correct BH-FDR behavior, not a bug — when even the best block is far from significant, the largest p becomes the shared FDR ceiling.)*
- The **primary block — the SCS axis, primitive→immune ({HSC_MPP, LMPP_GMP} → {T_NK, B_Plasma, Mono_DC}) — is exactly the hypothesized LSC→microenvironment interaction.** It is *nominally* stronger in AML (ΔC = −0.009) but **p = 0.86 — not significant.**

**Reading it against v1.1:** external work (eCHROMA, 1,563 cases) has established the *precondition* for H1 — LSC-like cells sit at the same early hierarchical position across all 16 epigenetic subgroups. Our job is the next step: does the *communication dependency* also converge? In CellChat, **the conserved axis does not surface** (the SCS axis is null). Two honest limits: this is **CellChat-only** (LIANA+/NicheNet pending), and v1.1's **triple stratification has not been run** — the current null is already there even before stratifying.

**Why it's still our most defensible result:** H1's edges come straight from CellChat and never pass through the malignant labels, so this null is immune to the label-quality problem on the next slide.

> **Backing figure — g03 (H1 null).**

---

## Slide 4 — What the stress tests localized: the label layer is the foundation crack

The blueprint calls for a **three-way consensus** call (inferCNV + Numbat + VarTrix, ≥ 2/3 agree). In practice: **125/130 samples are single-arm inferCNV**; only 5 reach two arms; **0 reach three.** Our healthy-marrow negative control (every malignant call there is false) quantifies the cost by bin:

| Niche bin | False-positive rate |
|---|---|
| **HSC_MPP (the LSC bin)** | **≈ 40%** |
| Erythroid / Megakaryocyte / Mono_DC | 23–27% |
| LMPP_GMP | 14% |
| B_Plasma / T_NK | 7.7% / 3.6% |

**Why HSC_MPP matters most:** an LSC *is* a malignant high-stemness cell in the HSC/MPP bins — the exact node the niche hypothesis leans on, and the dirtiest in the graph. Lymphoid clean, primitive/myeloid dirty — the signature of inferCNV keying on primitive-cell expression (a method artifact, not rare clones). v1.1's R10 adds an independent mechanism for the *MRD* version of this blind spot — quiescent LSCs transcribe little **and** sit in the metaphysis that aspirates under-sample ("dual invisibility").

**On stroma (FPR 90.3%):** not an error we ship — stroma is **excluded from the graph** precisely *because* of this rate. Its 90% is a **diagnostic** that inferCNV is unreliable on non-hematopoietic/primitive cells; it never enters the analysis. The in-graph problem to fix is HSC_MPP.

> **Backing figure — g02 (label quality).**

---

## Slide 5 — H2: the signal sits exactly on the pre-registered control arm (topology adds nothing *yet*)

FGW mixes two costs by weight **α**: α weights **topology** (Gromov-Wasserstein), (1−α) weights **node features** (Wasserstein). α=1 = pure topology, α=0 = pure node features — and **v1.1 makes α=0 H2's formal control arm (c).**

The α-sweep result: **at pure topology (α=1) the effect is null** (p = 0.11); the discrimination lives entirely at the **feature end** (α=0: p = 0.0009). So by H2's own (c) criterion, **topology currently adds no increment over node features.** Feature decomposition (α=0.5):

| Feature | β | p | Reading |
|---|---|---|---|
| `frac_malignant` alone | 0.173 | <0.001 | **circular** — the blast fraction, *forced to 0* in healthy, so it re-encodes the AML/healthy label |
| drop `frac_malignant` | 0.064 | 0.001 | residual… |
| cell count alone | 0.050 | 0.008 | …largely **technical** |
| stemness alone | −0.015 | 0.56 | the intended biology — **null, wrong-signed** |

**But this is not yet the decisive v1.1 test.** The α=0 winner is itself **circular** (frac_malignant), and the pieces that would give topology a fair chance — **M8's honest 7-block node features** (within-sample rank, de-circularized) and **M6's direction components** (cosθ, HDS⊥) — **are not implemented.** So today's verdict is "not topological *with the un-amended v1.0 features*"; the real H2 test is pending those.

> **Backing figure — g04 (H2 decisive).**

---

## Slide 6 — H3: the V-shape is a *pre-registered model*, not a failure — but the test is un-run

Two things must not be conflated:

- **What we observed (stemness along the axis):** a **V-shape** (dips at MRD, rebounds at Relapse). Under v1.0 this looked like "failing the monotone prediction." **It is not** — v1.1's M6 replaced monotone with a **three-model competition**, and the V-shape is the pre-registered **Model B**, motivated by our own 03_hierarchy observation *and* external non-monotone dynamics (and mechanistically expected via the MRD dual-invisibility above). So the direction of the early signal is **consistent with a pre-specified alternative.**
- **What we have NOT done (the actual H3 test):** the three-model (A/B/C) competition was never fit, and it *cannot* be until **M6's direction components (HDS∥ / HDS⊥ / cosθ) are computed** — those are the quantities that distinguish B and C from A, and they don't exist yet. On top of that: the signal is **patient-confounded** (62% of MRD malignant cells from one sample, 6323_MRD), rests on **one dataset** (GSE227903), and its timepoint/pairing metadata are **not yet reconciled** against the source papers.

**Net:** H3 is **not testable today** — not because the trend "failed to be monotone," but because the machinery that adjudicates the three models is un-implemented and the metadata is unverified. The early stemness V is a *hint toward Model B*, nothing more.

> **Backing figure — g05 (H3 blocked).**

---

## Slide 7 — The fix is scoped: harden the label layer first, then implement v1.1

**Highest-leverage single fix = a second, cohort-wide evidence arm for the malignant call: CopyKAT / SCEVAN** (read the QC object directly, **no BAM**). Honestly:
- **Buys (real, testable):** they define "who is normal" differently from inferCNV (which needs an external reference — our weakest link), so agreement cuts **method-specific** false positives. We will *measure* the gain — re-run the healthy-marrow FPR control on them; **if HSC_MPP FPR stays ~40%, the confounder is biological and consensus won't fix it.**
- **Does NOT buy:** same modality (expression) → not orthogonal like Numbat/VarTrix, may share the HSC_MPP failure; and **copy-neutral clones (NPM1 normal-karyotype) are invisible to all expression-CNV methods** — that's the **SNV/VarTrix arm's** job (Petti2019), currently an inactive stub (`SNV_ARM_ACTIVE=FALSE`), a separate track.

**Sequencing (refined against v1.1):**
1. Reconcile metadata (timepoint / patient-pairing).
2. Add CopyKAT/SCEVAN → re-run consensus labelling; measure the HSC_MPP FPR gain.
3. Re-run **only Phase 6–8** (node features → FGW → scoring). **Phase 5 / CellChat need not re-run** — edges are label-independent (which is also why H1 is), per v1.1 §11.
4. Then implement the v1.1 backlog, in this order of leverage: **M8** honest node features → **M6** direction components (unlocks the real H2 and H3 tests) → multi-method CCC gate (LIANA+/NicheNet) → **M7** functional axes, **M9** NSD validation, **D4** ATAC subgroups, Phase 8 robustness.

**Take-home:** the design is sound and self-correcting — v1.1 already anticipates the pitfalls our first pass exposed. The results are early: one clean label-independent null (H1, in CellChat), and a foundation label problem that must be fixed *before* the richer v1.1 tests are worth running (else we compute elegant features on 40%-contaminated labels).

---

## Backup A — verdict table

| Phase / hypothesis | Verdict | One-line reason |
|---|---|---|
| Phase 1 — data + consensus malignancy | 🔴 does-not-stand | 125/130 single-arm inferCNV; HSC_MPP FPR ≈ 40% |
| Phase 2 — hierarchy projection | 🟢 stands | clean reference-only projection, orthogonal labels |
| Phase 3 — meta-programs / SRRS | 🟡 partial | programs extracted; SRRS screen (now triple-stratified in v1.1) not computed |
| Phase 4 — CCC graphs | 🟡 partial | 148 graphs; single-method (CellChat); LIANA+/NicheNet + two-method gate absent |
| Phase 5–6 — distance + FGW | 🟡 partial | runs; un-amended v1.0 3-dim features (M8 pending) + circular healthy frac_malignant |
| Phase 7 — emergent edges + scoring | 🟡 partial | machinery correct; result strongly negative; M6 direction components not implemented |
| Phase 8 — robustness | ⚫ not-run | M9/S9/S10 absent |
| **H1** convergence | 🔴 **null** (label-independent, CellChat-only) | 0/49 edges, blocks q = 0.93, SCS axis p = 0.86; triple stratification un-run |
| **H2** topology (+direction) > pair | 🔴 **not topological (v1.0 features)** | α=1 null; signal = circular composition; M8/M6 decisive test pending |
| **H3** trajectory-model competition | ⚫ **blocked / un-run** | V-shape = pre-registered Model B; A/B/C fit needs M6 direction + metadata |

**Figure index:** g01 cohort counts · g02 label quality · g03 H1 null · g04 H2 decisive · g05 H3 blocked — in `results/figures/11_assessment/`.

---

## Backup B — Anticipated Q&A (speaker notes)

**Q: What changed between v1.0 and v1.1, and does it invalidate these results?**
No — no number changes; the *interpretation* of H2/H3 sharpens. v1.1 (2026-07-22): H3 monotone → **A/B/C three-model competition** (M6, adds HDS direction components); H2 gains a formal **α=0 pure-feature control arm** + a direction-increment test (M8); node features expand 3-dim → **7 semantic blocks, within-sample rank** (M8); plus functional-axis decomposition (M7), NSD receiver validation (M9), RNA-assigned ATAC subgroups (D4), and three new reviewer defenses (R10 spatial collapse / R11 epigenetic stratification / R12 post-translational ligand processing). Our first-pass results are the pre-amendment pipeline.

**Q: The stemness dips at MRD — isn't that the opposite of what you predicted?**
Under the *old* monotone framing, yes. Under v1.1 it is **Model B (rebound/V)**, a pre-registered alternative — motivated by our own hierarchy data and by external non-monotone state dynamics, and mechanistically expected (at MRD the quiescent LSCs both transcribe little and are under-sampled by aspiration). So it's a hint toward B, not a failure. But we can't *adjudicate* A vs B vs C until the direction components (HDS⊥, cosθ) are computed.

**Q: All 5 block q-values are 0.93 — did the FDR break?**
No. Raw p differ (0.34, 0.45, 0.84, 0.86, 0.93). BH step-up: the largest p (0.93) sets q, and every smaller block's inflated p·5/j exceeds 1, so monotonicity pulls them all to 0.93 — the correct signature of a uniformly null result. n_perm = 10,000; verified line-by-line.

**Q: What is the SCS axis (the "primitive→immune" block)?**
The pre-registered primary block: LSC/primitive {HSC_MPP, LMPP_GMP} → immune/paracrine {T_NK, B_Plasma, Mono_DC}, 6 edges; locked in `config_ccc.R` before results. It *is* the hypothesized LSC→microenvironment interaction. Result: nominally stronger in AML (ΔC = −0.009) but p = 0.86.

**Q: Isn't CopyKAT/SCEVAN just inferCNV again?**
Different algorithms (Bayesian-KS / variational vs moving-average) and, crucially, different ways of defining the normal baseline (they self-calibrate; inferCNV needs an external reference). So consensus removes *method-specific* artifacts — but same modality, so it may share *biological* errors like the HSC_MPP over-call, and it stays one evidence "type" (not the cross-modality 2/3). We test the actual HSC_MPP FPR drop rather than assert it.

**Q: What about copy-neutral malignant cells (NPM1, normal karyotype)?**
Invisible to all expression-CNV methods. That's the SNV/VarTrix arm's job (Petti2019 first) — currently an inactive stub, a separate reactivation track.

**Q: Do we have to re-run CellChat after re-labelling?**
No. CCC edges depend on hierarchy bins + expression, not on malignant labels (that's why H1 is label-independent). Re-labelling invalidates node features + FGW + scoring (Phase 6–8) only — consistent with v1.1 §11.

**Q: Why is healthy `frac_malignant` = 0?**
A config switch (`FGW_ZERO_HEALTHY_MAL`, `01_build_fgw_inputs.R:67`) sets it deterministically — the source of that feature's circularity. Fixable; on the deferred bug list, and superseded by M8's feature discipline.
