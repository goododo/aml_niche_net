# Findings: pseudobulk AML-vs-healthy differential expression, per hierarchy bin

Last result produced 2026-09-24. The analysis was pre-registered in
`PREREGISTRATION_pseudobulk_de.md`, committed 2026-09-16 (`a66cbd0`) before any pseudobulk matrix
existed. Every number below is reproduced by a named script on a named sample set (§16). Nothing
here is an estimate or a recollection. Where a result is mixed or unresolved it is written as mixed
or unresolved.

**Why this line exists.** Two questions were raised about the analysis in September 2026: how much
of the single-cell data actually enters it, given that everything is reduced to 7 nodes, and whether
node vectors should be built from real expression rather than from derived scores. Both were
well founded. The node vector is 150 features, ~141 of them bin-level means of signature, pathway,
cNMF-program, pseudotime or CNV scores, and only 9 gene-expression features (3 genes x 3 strata). A
repo-wide search confirmed **no cohort-wide, full-transcriptome, model-based DE had ever been run
here** — no `edgeR`, `limma`, `DESeq2`, `muscat` or `FindMarkers` call existed in `scripts/`. This
closes that gap.

**Where the line stands.** The pre-registration predicted a null (§2 of that document) and the
prediction was wrong in exactly one place. **Mono_DC carries a real, cross-platform-replicated
transcriptional difference**: 385 reportable hits, exact permutation p = 0.002 and 0.001 under the
two pre-registered nulls, and 77.3% sign concordance on a held-out arm run on a different sequencing
platform (one-sided binomial p = 3.8e-16). Three bins return clean nulls with quantified detection
limits. One bin was dropped by its own marker gate, for the reason the design had predicted in
advance.

**And the content of the signal argues that the bin label is the variable.** The up side is the
HOXA/HOXB clusters plus `MEIS1`; the down side is 255 genes of mature monocyte / cDC2 identity. The
most economical reading is that the cells a *healthy*-reference projection assigns to "Mono_DC" in
AML marrow are substantially not mature monocytes but HOX-high immature cells — i.e. this is a
composition difference inside the bin as much as a per-cell expression difference. **That reading is
an interpretation, not a result** (§7), and the design cannot separate the two — the
pre-registration said so before the numbers existed.

**It also explains the topology null rather than contradicting it.**
`08_scoring/FINDINGS_topology_null.md` closes the AML-vs-healthy question through graph topology as
negative across 11 cost x mass configurations (p 0.305-0.966). If the signal lives in *within-node
composition*, then a 7-node graph summary discards exactly the dimension that carries it. The two
results are two sides of one story, and that is worth more than either alone (§8).

---

## 0. Sample sets

Only the **3 of 10 datasets that contain both arms** enter. For the other 66 samples dataset and
disease status are perfectly collinear, and no covariate separates perfectly collinear factors.

| | Discovery | Validation |
| --- | --- | --- |
| datasets | Chen2023, GSE185381 (10x) | GSE116256 (Seq-Well) |
| AML | 42, all `Diagnosis` | 9, all `Diagnosis` |
| healthy | 14 | 2 |
| patients | all distinct, zero repeats | 9 distinct |

72 both-arm samples exist; **67 enter**. Five were removed by pre-registered rule before any model
ran: `BM5-34p` (CD34+ magnetic-bead sorted and fresh, against an otherwise viability-only /
Ficoll-mononuclear / cryopreserved cohort) and four post-treatment marrows from patients already
contributing their own Diagnosis sample. The Validation control arm is therefore **n = 2**, and that
was fixed in advance rather than discovered at write-up.

Controls are **disjoint** between arms — a genuine improvement on
`08_scoring/PREREGISTRATION_panel_screen.md`, where only 23 healthy samples exist cohort-wide and
the controls had to be shared. It is also why Validation is this small.

**The Discovery/Validation split is respected here for the first time in this project**
(`FINDINGS_topology_null.md` §16 limitation 5 records that no prior analysis did).

### 0.1 Stage A verified against an independent path

| check | result |
| --- | --- |
| samples aggregated | 67 |
| (sample, bin) units | 449 |
| sum over bins == filtered count total | **exact, all 67** |
| sum over bins == filtered cell count | **exact, all 67** |
| barcode join fraction | **1.000, every sample** |
| `Early Lymphoid` cells removed (`CCC_EXCLUDE_FINE`) | 2585 |
| units reproducing `ccc_node_features.csv` `n_cells` | **449 / 449** |

The last row is the strongest check in this document. Stage A reaches the production cell set by its
own route — reconciled annotations, its own filter chain, its own aggregation — and lands on the same
449 numbers. It is a cross-check, not a restatement.

---

# PART I — What was established

## 1. The whole result in one table

| tier | bin | n (AML v healthy) | healthy libs | genes | reportable hits | exact permutation p | LFC80 | Validation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary | **Mono_DC** | 40 v 14 | 11 | 11029 | **385** | **0.002 / 0.001** | 1.5 | **163/211 = 77.3%**, p = 3.8e-16 |
| primary | **B_Plasma** | 40 v 14 | 11 | 7734 | **28** | **0.025 / 0.023** | 1.5 | 11/13 = 84.6%, p = 0.011 |
| primary | T_NK | 39 v 14 | 11 | 9281 | 0 (null) | 0.102 / 0.320 | 1.5 | not evaluable |
| secondary | **LMPP_GMP** | 5 v 4 | 4 | 9497 | **175 robust** (+129 contingent) | **0.016** | 1.5 | no arm by design |
| secondary | HSC_MPP | 4 v 4 | 4 | 8981 | 0 (null) | 1.0 | 2.0 | — |
| — | Erythroid | 3 v 3 | — | — | **gate failure** | — | — | — |

Primary-tier p-values are the two pre-registered nulls (`sample_within_dataset` /
`library_arm_pure`) at `N_PERM = 1000`. The secondary tier's is **exact**, enumerated over all 126
arrangements — see §4 for why the earlier sampled value was withdrawn. `Megakaryocyte` was dropped
before any model ran: zero GSE185381 controls clear 30 cells, so its healthy arm is one dataset.

## 2. Mono_DC is real, and it replicates across sequencing platform

Three independent lines, none of which is the others restated:

1. **Against its own null.** 537 baseline hits, 385 surviving the mandatory depth arm, against a
   null whose median is 0 and whose mean is 3.3 (sample null) / 4.5 (library null). p = 0.002 and
   0.001. Both nulls agree.
2. **Across platform.** Discovery is 10x (Chen2023 + GSE185381); Validation is Seq-Well
   (GSE116256), different donors, different library chemistry. 163 of 211 testable frozen hits keep
   their Discovery sign. One-sided binomial p = **3.8e-16**.
3. **Asymmetrically, and the asymmetry is informative.** Reported here for the first time — the
   pre-registration only required the pooled figure:

   | direction | concordant | one-sided binomial p |
   | --- | --- | --- |
   | **up in AML** | **63 / 71 = 88.7%** | **5.1e-12** |
   | down in AML | 100 / 140 = 71.4% | 2.1e-07 |

   The up side — the HOX/`MEIS1` side (§6) — is close to fully reproducible. The down side, mature
   myeloid identity, is weaker. Both clear 0.5 decisively.

**The honest caveat on Validation size.** 174 of the 385 hits (**45%**) never reach the test: they
are absent from GSE116256's 27899-gene deposit. And the Validation arm is **7 AML vs 2 healthy** in
this bin. A 77.3% concordance on 211 genes against 2 control samples is real evidence of direction
and is not evidence about effect size.

## 3. Three clean nulls — and "clean" here has a number attached

A null is only a result if the study could have found something. The pre-registered planted-effect
control (§5 of the pre-registration, as amended) gives each bin an **LFC80**: the smallest planted
logFC at which 80% of 200 spiked middle-tertile genes are recovered under the full hit rule.

| bin | LFC80 | the null it licenses |
| --- | --- | --- |
| T_NK | 1.5 | no effect of \|logFC\| >= 1.5 at this prevalence |
| HSC_MPP | 2.0 | no effect of \|logFC\| >= 2.0 at this prevalence |
| Mono_DC / B_Plasma / LMPP_GMP | 1.5 | (not null — see §2, §4) |

**No bin is power-limited.** T_NK's single baseline hit does not survive either null (p = 0.102,
0.320) and is reported as zero. HSC_MPP returns 0 hits against a null whose maximum over all 70
arrangements is 1.

T_NK deserves one sentence on its own, because it is the cleanest comparison in the study: it has
**the same n, the same datasets, the same model, the same library blocking and the same LFC80 as
Mono_DC**, and returns nothing. That rules out the pipeline manufacturing hits from design
structure, which no permutation test can rule out by itself.

## 4. LMPP_GMP is real but fragile, and it is smaller than first reported

Two corrections were applied to this bin after its first pass, both recorded in the
pre-registration:

**(a) The reportable count is 175, not 304.** `PREREG 8.5` requires every secondary hit to be
labelled robust or contingent by leave-one-out over each healthy donor, and it had never been run.
Run as a full refit per fold:

| fold | frozen hits surviving |
| --- | --- |
| drop NBM1 | 231 / 304 |
| drop NBM2 | 280 / 304 |
| **drop NBM3** | **203 / 304** |
| drop NBM4 | 220 / 304 |

**175 robust, 129 contingent.** 42% of this bin's hits depend on one donor being present; NBM3 alone
carries 101 of them. With 4 controls that is what the design can support.

**(b) The permutation p was false precision and is withdrawn.** `03_permute.R` drew 1000 labellings
*with replacement* from a space holding only `C(9,4) = 126` distinct arrangements and reported
p = 0.00699 — **below the design's own minimum attainable p of 1/127 = 0.00787**, an impossible
value. Enumerated exactly:

| scheme | p |
| --- | --- |
| frozen voom (comparable to the primary tier) | 1/126 = 0.0079 |
| **full refit (the correct one, and the one to quote)** | **2/126 = 0.0159** |

The two schemes differ by a factor of two, which is itself a finding about the voom-freezing
shortcut `03_permute.R` declared in its header: the more correct scheme is weaker, and one
arrangement out-hits the observed data (458 vs 431). **Do not read `perm_summary.csv` for this tier;
read `secondary_exhaustive_perm.csv`.**

**(c) The two nulls are one null here.** Chen2023's `library_id` is 1:1 with sample (9 samples, 9
libraries), so permuting samples and permuting arm-pure libraries enumerate an identical set.
`perm_summary.csv`'s 0.00699 vs 0.00899 was one quantity measured twice with noise. The distinction
stays real for the **primary** tier, where GSE185381's shared lanes make it bite.

## 5. Erythroid failed its own marker gate, for the predicted reason

`PREREG 5.4` predicted from the Chen2023 protocol — five FACS fractions, every one gating out
CD235a+/CD71+ cells — that this bin could not be read as erythroid biology. The gate then found it
in the expression data independently:

| | healthy arm | AML arm | log2 ratio |
| --- | --- | --- | --- |
| `HBB` | 330 CPM | 2 CPM | **7.1** |
| `HBA1` | 10 CPM | 0 CPM | inf |

Two independent routes — a sorting protocol recorded in `00_curated_manifest.csv`, and haemoglobin
expression — reach the same conclusion. The bin is excluded and nothing from it is reported.

---

# PART II — What the signal is made of

## 6. HOX up, mature myeloid identity down

Descriptive characterisation of the frozen 385-hit list. **Not a hypothesis test** — no p-value in
this section is claimed as evidence. (A formal MSigDB enrichment was attempted and abandoned: the
installed `msigdbr` fetches gene sets from Zenodo at run time and the compute nodes have no
outbound network. Family counting needs no network and no p-value.)

**Up in AML — 130 genes (34%)**

| family | n | median logFC | testable in Validation | sign-concordant |
| --- | --- | --- | --- | --- |
| HOXB cluster | 5 | +4.31 | 5 | **5 / 5** |
| HOXA cluster | 5 | +4.02 | 3 | **3 / 3** |
| HOX antisense | 2 | +5.33 | 1 | **1 / 1** |
| `MEIS1` | 1 | +4.20 | 1 | **1 / 1** |

**12 distinct HOX/`MEIS1` genes; 10 are testable in Validation and all 10 keep their sign.** The two
that are not testable (`HOXA7`, `HOXA10-AS`) are absent from GSE116256's deposit, not discordant —
a distinction worth making explicitly, because collapsing "not tested" into "failed" understates
this by two genes.

`HOXB-AS3` +6.57, `HOXA3` +4.62, `HOXB5` +4.35, `HOXB3` +4.34, `HOXB6` +4.31, `MEIS1` +4.20,
`HOXA10-AS` +4.10, `HOXA9` +4.02, `HOXA7` +3.69, `HOXA10` +3.40, `HOXB4` +3.28, `HOXB7` +2.98 —
**13 genes, half of the top 25 on the up side.**

**Down in AML — 255 genes (66%)**

| family | n | median logFC | Validation concordant |
| --- | --- | --- | --- |
| cDC2 identity | 6 | −3.33 | 5 / 6 |
| mature monocyte | 5 | −2.75 | 4 / 5 |

`FCER1A` −3.42, `CLEC4D` −3.37, `CLEC10A` −3.35, `CLEC4E` −3.30, `CD1C` −3.16, `CLEC4A`;
`FCGR3A` −2.75, `FPR1`, `FPR2`, `G0S2`, `RETN`. The largest-magnitude down genes overall are
`LYPD2` −5.60, `NRG1` −5.11, `RBP7` −4.90, `STEAP4` −4.82, `HES1` −4.46, `WLS` −4.33.

HOXA9/`MEIS1` co-expression is the canonical signature of NPM1-mutant and KMT2A-rearranged AML. The
down side is not a pathway being modulated; it is a differentiation state being absent.

## 7. The composition reading — and why this design cannot settle it

**The reading.** A pseudobulk profile averages over sub-states within its bin. If the AML "Mono_DC"
compartment contains HOX-high immature cells where the healthy one contains mature monocytes and
cDC2, then a composition shift presents as differential expression, and the 385 genes are partly a
description of *which cells are in the bin* rather than *how the same cells differ*.

**What supports it:** the up side is a progenitor/leukaemic-stemness programme, the down side is 255
genes of mature identity loss, and the `PREREG 7.3` diagnostic labels Mono_DC
composition-confounded.

**What limits it, stated plainly:**

- The composition label for Mono_DC rests on **one** cNMF programme at **1.016 SD**, against a
  1.0 SD line **this analysis chose, not the pre-registration**. HSC_MPP sits at 0.968 on the other
  side. Both labels flip under a 4% change of threshold. LMPP_GMP at 1.470 (4 programmes) does not.
- `PREREG 7.3` fixed in advance that this is a **diagnostic, not a gate**, precisely because
  "the two cannot be separated with this design". That remains true. Nothing here upgrades it.
- The obvious external anchor — blast burden — has been tried and is **unresolved** (§14).

**So: the composition reading is the most economical explanation of the gene content, and this
document cannot establish it.** It was the right target for the next pre-registration — and that
pre-registration has since been written and run.

### 7.1 Update 2026-09-24 — `PREREGISTRATION_hox_compartment.md` answered it

The question above is no longer open. `06_hox_compartment.R`, on 36 samples / 9 NPM1-mutants
(GSE185381 + GSE116256; Chen2023 excluded for having zero within-dataset contrast at 5/5 mutant),
asked whether the Mono_DC HOX signal is separable from the LMPP/GMP blast compartment:

| | AUC for NPM1 | permutation p |
| --- | --- | --- |
| Mono_DC, raw (reference) | 0.893 | — |
| LMPP_GMP, raw (reference) | 0.975 | — |
| **Mono_DC adjusted for LMPP_GMP** | **0.490** | **0.549** |
| LMPP_GMP adjusted for Mono_DC | 0.733 | 0.021 |

**Adjusted for the blast compartment, Mono_DC lands at exactly chance, and the reverse direction
survives.** Three pre-registered controls agree: a Set D repeat using 12 Discovery-frozen genes
instead of 19 locus-defined ones gives 0.523 / p 0.36 and 0.753 / p 0.0075; a T_NK negative control
does not add (0.429, p 0.71), so the "adding" is not generic; and cell count is not a confound
(rho −0.225 and −0.170 against a 0.40 threshold). Compartment collinearity was rho 0.822 / 0.857,
under the pre-registered 0.90 stop.

**What this settles, at the wording ceiling that document fixed in advance:** the Mono_DC
compartment's HOX signal is **not separable** from the LMPP/GMP compartment's. The 385 genes are
**not a microenvironment finding**. The composition reading of §7 is therefore supported rather than
merely economical — and §13's first bullet should be read with that constraint attached.

**What it still does not settle:** that the cells carrying the signal are malignant. That needs the
cell-level malignancy call which `scripts/FINDINGS_project_status.md` N4 records as not existing in
this project (inferCNV vs clinical blast %, n=59, rho = −0.069, p = 0.605).

See `scripts/FINDINGS_project_status.md` §2.6 for this result in its cross-line context; it is the
project-level synthesis and this document is the line-level one.

## 8. How this relates to the topology null

`FINDINGS_topology_null.md` closes AML-vs-healthy through graph topology as negative: 11 cost x mass
configurations at alpha=1 give p 0.305-0.966, the transport-free Frobenius baseline p = 0.679, and
the per-edge regression finds 0 of 49 edges within dataset. Its §6 localises the cause to
information being encoded in the one dimension the statistic discards.

This line adds a second instance of the same mechanism, one level down. If the AML-vs-healthy signal
sits in **composition within a node**, then reducing a sample to 7 nodes deletes it by construction —
whatever is done with the graph afterwards. That is consistent with the split-half result in
`FINDINGS_topology_null.md` B.7, where the strongest patient-identity arm was **seven cell-type
proportions with no communication in them**.

Both lines therefore point at composition, from opposite directions. Stated as a convergence, not as
a test: no analysis here was designed to compare the two, and none is claimed to.

---

# PART III — The controls, and what each one actually did

## 9. Depth — the mandatory arm removed 288 hits, and 8 of them prove the gate's shape mattered

AML libraries are systematically deeper cohort-wide (`med_ncount_final` 4874 vs 3006.5, AUC 0.695,
p = 0.0033; `DECISIONS_pending.md` O6). Stage A therefore computed **per-(sample, bin)** depth,
which did not previously exist anywhere in `results/tables/` — and the per-bin picture does **not**
inherit the sample-level one. Per-bin median UMI on the Discovery set the model is fit on, units
clearing 30 cells:

| bin | AML | healthy | deeper arm |
| --- | --- | --- | --- |
| Mono_DC | 4474 | 3213 | AML |
| T_NK | 2534 | 2158 | AML |
| Megakaryocyte | 2297 | 1762 | AML |
| B_Plasma | 2612 | 2798 | **healthy** |
| Erythroid | 5064 | 5976 | **healthy** |
| HSC_MPP | 6760 | 7273 | **healthy** |
| LMPP_GMP | 6526 | 7758 | **healthy** |

The cohort-level confound has AML deeper (4874 vs 3006.5). **Per bin it reverses in four of seven**,
including both bins where the secondary tier found hits. Using the sample-level value as the
covariate would have been wrong in sign for the majority of bins — which is the whole reason
`PREREG 4.3` insisted the depth covariate be computed per bin at aggregation time, when it is free
and before it is unrecoverable.

288 hits were struck. Which half of the hit rule failed:

| failed condition | n |
| --- | --- |
| `q` only | 256 |
| both | 24 |
| **\|logFC\| only** | **8** |

Those 8 Mono_DC genes keep `q < 0.05` with depth in the model but their `|logFC|` falls below 1. The
handoff specified a gate on significance alone; the pre-registration changed it to test **both**
halves of the hit definition. Without that change, 8 depth-dependent genes would have been reported
as depth-independent. **The amendment caught a real instance, not a hypothetical one.**

## 10. Mapping error — the confound runs in the conservative direction, in all seven bins

Bin labels come from projection onto a *healthy* reference, and the `!high_error` filter removes
proportionally more AML cells. Recomputed per bin from the reconciled per-cell files (the repo only
stored this per sample), Discovery, cellwise:

| bin | AML | healthy |
| --- | --- | --- |
| Mono_DC | 1.73% | 0.55% |
| B_Plasma | 0.62% | 0.53% |
| T_NK | 7.76% | 6.65% |
| LMPP_GMP | 7.44% | 4.66% |
| HSC_MPP | 21.6% | 9.4% |
| Erythroid | 4.17% | 1.69% |
| Megakaryocyte | 0.24% | 0.00% |

**AML loses more in every bin.** So the surviving AML cells are the ones most resembling the healthy
reference, and `logFC` is biased **toward zero**. A null here is conservative, and the Mono_DC hits
are not inflated by this mechanism. `PREREG 7.2` pre-stated this direction; it now has numbers.

## 11. Sex — no imbalance, and it could not have mattered anyway

No sex or age column exists anywhere in this repo. chrY (111 genes) and `XIST` were therefore
excluded from every tested universe *a priori*, and sex was inferred post hoc from those same genes
purely as a diagnostic:

| | male | female | ambiguous |
| --- | --- | --- | --- |
| Discovery AML | 25 | 9 | 8 |
| Discovery healthy | 9 | 2 | 3 |

Fisher exact p = **0.70**, odds ratio 1.60. No detectable imbalance. Since the genes were already
excluded, this table cannot have changed a hit either way; it exists because a reviewer will ask.

## 12. The multiple-testing sensitivity arm behaves asymmetrically

`PREREG 6.5` made BH **within bin** primary — reversing the handoff — because Chen2023's bins come
from two physically separate sequencing libraries and the per-bin gene universes differ. The
across-bin pooled arm was reported alongside, over the 7323-gene intersection of the three primary
bins:

| bin | within-bin hits | pooled hits | disagree |
| --- | --- | --- | --- |
| Mono_DC | 149 | 70 | 79 |
| B_Plasma | 15 | 15 | 0 |
| T_NK | 1 | **4** | 3 |

Pooling is **more** conservative for Mono_DC (149 → 70) and **less** conservative for T_NK (1 → 4),
because pooled BH lets T_NK borrow from Mono_DC's strong tail. That is the concrete reason the
within-bin correction is primary: a pooled denominator moves hits into a bin that has none of its
own. The conclusion does not change under either arm — T_NK's 4 pooled hits still fail the
permutation null.

---

# PART IV

## 13. What can be claimed

- AML and healthy bone marrow **differ transcriptionally in the Mono_DC compartment**, reproducibly
  and across sequencing platform. 385 genes, permutation p 0.002 / 0.001, 77.3% sign concordance on
  a held-out platform (p = 3.8e-16).
- The difference is **HOXA/HOXB + `MEIS1` up, mature monocyte / cDC2 identity down**.
- **B_Plasma** carries a smaller real difference: 28 genes, p 0.025 / 0.023. Its Validation arm is
  3 AML vs **1** healthy sample; that arm is directional support at best.
- **LMPP_GMP** carries 175 leave-one-out-robust genes at exact p = 0.016, from a 5-vs-4
  single-dataset contrast.
- **T_NK and HSC_MPP are null**, at \|logFC\| >= 1.5 and >= 2.0 respectively.
- **Erythroid cannot be assessed** in this cohort.
- Per-bin sequencing depth does not inherit the cohort-level depth confound and must be computed per
  bin (§9). This is a reusable methodological point for the project.

## 14. What cannot be claimed

- **Not resolved: whether the HOX axis is blast burden.** `PREREGISTRATION_hox_axis.md` Q1 returned
  **INDETERMINATE** — rho = **−0.43** in GSE185381 (n = 17) against **+0.50** in GSE116256 (n = 7),
  **opposite signs**, neither reaching the pre-set 0.6. Its §2.1 fixed the reading in advance for
  exactly this outcome: *"No claim. n is reported and the question is left open."* Q2 (NPM1) passed
  at p = 0.029, **n = 7**, with the FLT3-ITD confound named in its §6. The axis is therefore neither
  established as blast contamination nor as a stratification axis.
- **Not a per-cell expression finding.** §7 — composition and expression are not separable here.
- **Not a niche or microenvironment finding.** `Stromal` is outside `CCC_NODES` and outside this
  analysis. The repo convention that "microenvironment" excludes structural stroma applies.
- **Not monocyte biology**, if the composition reading holds.
- **Not a cell-level malignancy caller**, not risk group, not treatment response, not survival.
- **Not an effect-size estimate.** The hit rule cuts at \|logFC\| >= 1 while LFC80 is 1.5, so the
  1.0-1.5 band is detected at roughly half power and the hit lists are incomplete there. This is a
  sensitivity limitation, not a false-positive one.

## 15. Limitations to carry into the write-up

1. **The pre-registration was amended twice and corrected once, all after seeing results.** Every
   change is a dated block in that document stating it was post hoc. They are: the §8.2 marker-gate
   criterion (top-50 by mean CPM is, on pseudobulk, a test of whether a marker outranks the ribosome
   — 88% of T_NK's top 50 is `RP*`/`MT-`); the §8.3 planted control (planted **at** the \|logFC\|>=1
   decision boundary, pinning recovery near 0.50 by construction regardless of n — 0.485/0.505/0.475
   at n = 54/53/54, while the same Mono_DC bin recovers 0.965 at logFC 2.0); and the §8.3
   secondary-tier p (§4b). **Each has a defensible reason and the cost is still real.** A document
   amended three times after the fact is weaker evidence than one that was not. For the marker gate
   the mitigation is that a structurally different repair — keeping absolute rank and widening to
   top 2% — returns identical verdicts for all four bins; two repairs agreeing is weaker than
   pre-specification and is not claimed to be equivalent.
2. **Discovery's healthy arm is 14 samples in 11 libraries at best.** GSE185381 is HTO-multiplexed
   and `Control1/2/3/4` share one lane (`lib__P02`), which is why `library_id` enters the primary
   model via `duplicateCorrelation` (consensus 0.169-0.224). Without it the healthy n is a fiction.
3. **Chen2023's two-pool nesting is not in any model.** Each donor is two separately sequenced pools
   (HSPC+myeloid | stromal+lymphoid) merged, so within Chen2023 library is nested in bin. The
   `library_id` column cannot express it.
4. **45% of Mono_DC's hits never reach Validation** — absent from GSE116256's 27899-gene deposit.
5. **Validation is 2 control samples** (1 in B_Plasma). It is a sign-concordance check, fixed as
   such before the arm was opened, and is not an independent replication.
6. **Three of six analysed bins came from a single dataset** (the secondary tier), because in
   LMPP_GMP / HSC_MPP / Erythroid the Discovery healthy arm is 95.0-97.5% Chen2023 cells against an
   AML arm that is 4.7-13.3%.
7. **No ambient-RNA correction exists anywhere in `scripts/`** (no SoupX, CellBender, decontX).
8. **inferCNV under-calls ~9x**, so no malignant/non-malignant split was used; every filtered cell
   in a bin enters its pseudobulk regardless of CNV call.
9. **Donor sex, age, karyotype and treatment response are not usable** at this cohort's coverage.

## 16. Provenance

| Claim | Script | Output |
| --- | --- | --- |
| Design, all decision rules, three amendments | `PREREGISTRATION_pseudobulk_de.md` | — |
| Roster (67 of 72), per-(sample,bin) counts + depth, §8.1 conservation | `01_aggregate.R` | `pb_sample_manifest.csv`, `depth/*.csv`, `PB_RDS_DIR/**` |
| §8.2 marker gate, per-bin DE, depth arm, frozen hit list, pooled-BH arm | `02_de_limma.R` | `marker_gate.csv`, `de/*.csv`, `de_summary.csv`, `hitlist_frozen.csv`, `sensitivity_pooled_bh.csv` |
| Both §8.3 nulls, planted-effect power curve, LFC80 | `03_permute.R` | `perm_null.csv`, `perm_summary.csv`, `power.csv` |
| Validation sign concordance | `04_validation.R` | `validation_concordance.csv`, `validation/*.csv` |
| §8.5 leave-one-out, §8.3 exhaustive permutation | `06_secondary_robustness.R` | `secondary_leave_one_out.csv`, `secondary_exhaustive_perm.csv` |
| §7.2 / §7.3 / §6.4 confounder tables | `07_confound_tables.R` | `confound_mapping_error.csv`, `confound_composition.csv`, `confound_inferred_sex.csv` |
| §11 item 6 depth-struck list | inline, from `de/*.csv` | `depth_struck_hits.csv` |
| HOX axis, separate pre-registration | `PREREGISTRATION_hox_axis.md`, `05_hox_axis.R` | `hox_axis_*.csv` |

All outputs under `results/tables/12_pseudobulk_de/`. All R run through `ENV_PREFIX`
(`/FAST/gr10634/gaozy/general_env`, R 4.4.3, limma 3.62.2, edgeR 4.4.2) — the login-node `Rscript`
is R 4.5 and cannot see limma.

**Self-checks that ran and passed.** Both §8.1 conservation equalities are exact on all 67 samples.
Barcode join is 1.000 on all 67. Stage A's 449 units reproduce `ccc_node_features.csv`'s `n_cells`
449/449 by an independent path. `filterByExpr`'s `MinSampleSize` is reproduced analytically at 7.043
for the three primary designs, matching a synthetic-data check run before the real fit. The
`armAML` coefficient is asserted present in every design, since R's default factor ordering would
silently invert the sign of every logFC in the study. The secondary tier's permutation space is
enumerated exhaustively rather than sampled, and the two enumeration schemes are both reported
because they disagree.

**One failure worth recording.** The first Stage B run crashed on
`setorder(D, adj.P.Val, -abs(logFC))` — `setorder` takes column names, not expressions. It crashed
rather than producing a wrong ordering, which is the good failure mode and the opposite of the
family of bugs `FINDINGS_topology_null.md` §14 warns about.
