# Phase-1 Foundation Assessment — aml_niche_net

**Date:** 2026-07-24
**Scope:** Every pipeline phase assessed against the project blueprint (v1.0 + v1.1 patch), each conclusion tagged **stands / partial / does-not-stand / blocked / not-run**.
**Audience:** internal (detailed). A lab-meeting-length version is extracted separately.
**Provenance:** every number below traces to a source table cell (cited inline). The numbers were produced by the pipeline and independently re-checked against the tables under `results/tables/`; 62 distinct figures were fact-checked and 1 overstatement corrected (noted in §H1).

---

## 0. Executive summary

The pipeline runs end-to-end and its statistical machinery is, if anything, unusually self-critical (permutation nulls, FWL residualization, a pre-registered alpha sweep, a platform control, and feature-decomposition diagnostics). But **none of the three headline hypotheses is supported by the current results. H2 and H3 come up short along a single causal chain that starts at the label layer; H1 fails independently of it:**

> Single-arm inferCNV labelling → **~40% false-positive rate in the HSC_MPP (LSC) bin** → the `frac_malignant` node feature is both unreliable and, for healthy samples, **forced to 0 by construction** → the only statistically significant term behind H2 is that circular feature plus raw cell count → **at pure topology (α=1) H2 is null** → and H3's only nominal signal is a single-patient, label-dependent artifact on metadata that is itself unverified (**H3 blocked**). H1 is *not* part of this chain — its emergent-edge test never uses the labels, and it is null on its own (see below).

**The one result that stands on its own is the H1 null.** H1's edges come straight from CellChat and never pass through the malignancy labels, so it is immune to the label-quality problem below — a genuine negative, not a borderline miss. One honest limit: it is CellChat-only. The blueprint's main CCC engine is LIANA+ (CellChat is only the sensitivity check), so H1 still has to be re-tested under LIANA+/NicheNet before we call it final. But the label problem cannot rescue it either way — precisely because the labels never touch it.

**Frame for the lab meeting:** this is not "our analysis failed." It is "we built the falsification tests most groups skip, and they told us the foundation (malignant labelling) needs a second evidence arm before any niche claim can stand." The fix is already scoped (CopyKAT/SCEVAN, see §Next).

### Verdict table

| Phase / hypothesis | Verdict | One-line reason |
|---|---|---|
| **Phase 1** — data curation + consensus malignancy | 🔴 **does-not-stand** | 125/130 single-arm inferCNV (spec: 3-way ≥2/3); HSC_MPP FPR 39.8%; ≥200-malignant gate not implemented |
| **Phase 2** — hierarchy projection / node vocabulary | 🟢 **stands** | Reference-only Symphony projection, 220/220 coverage, orthogonal malignant layer, correct 7-bin vocabulary |
| **Phase 3** — meta-programs / cNMF + SRRS | 🟡 **partial** | MP extraction + tumor/normal flag real; the defining **dual-layer SRRS reproducibility screen is entirely absent** |
| **Phase 4** — CCC graph construction | 🟡 **partial** | 148 CellChat graphs built, but the **R3 two-method (LIANA+ ∧ CellChat) robustness gate does not exist** |
| **Phase 5–6** — distance + FGW barycenter | 🟡 **partial** | Runs and is directionally sane, but implements the **un-amended v1.0 feature vector (M8 ignored)** and a **circular** healthy `frac_malignant`=0 |
| **Phase 7** — emergent edges + patient scoring | 🟡 **partial** | Machinery correct and self-falsifying; **scientific output strongly negative** (H1 null, H2 not topological) |
| **Phase 8** — robustness & cross-validation | ⚫ **not-run** | `09_robustness/` empty, no script dir; perturbation-recovery, bootstrap CI, external validation all absent |
| **H1** (conserved emergent edges) | 🔴 **null result** | 0/49 edges emergent (min q=0.865); 5 blocks all q=0.927 |
| **H2** (topology + direction carry information) | 🔴 **not topological (v1.0 features)** | α=1 (pure topology) null (p=0.105/0.271); signal = circular `frac_malignant` + cell count. This is exactly the α=0 result on v1.1's own control arm (c); the decisive test (M8 features + M6 direction) has not been run yet |
| **H3** (treatment-pressure trajectory) | ⚫ **blocked / un-run** | The V-shape is v1.1's pre-registered Model B, not a failed monotone prediction; the A/B/C model competition can't run until the M6 direction scores are built (they aren't); 62% of MRD malignant cells from one patient; metadata unverified |

Figures backing this: `results/figures/11_assessment/g01`–`g05`.

---

## 1. The central finding — label quality is the foundation crack

The blueprint's Phase 1 specifies a **three-way consensus** (inferCNV + Numbat + VarTrix, malignant if ≥2/3 arms agree). What was realized:

- **125/130 samples are single-arm inferCNV** (`tier=C_single`); only 5 (all GSE227903) reach two arms; **0/130 reach the specified three arms**. VarTrix/SNV never contributed to any sample, despite Petti2019 being a known-mutation dataset.
  *(ALL_consensus_summary.csv: evidence_tier → C_single=125, B_multi_partial=4, A_concordant=1)*
- Even the 5 two-arm samples disagree heavily: 3 have `conflict_frac ≥ 0.40` (3853_Dg=0.76, 6323_Dg=0.50, 1886_Dg=0.40) — inferCNV and Numbat disagree on 40–76% of cells.

The negative control the team built (`96_malignancy_fpr_healthy.R`) then quantifies the cost, by bin:

| Bin (in CCC graph) | False-positive rate | n (healthy cells) |
|---|---|---|
| **HSC_MPP** | **0.398** | 2,469 |
| Erythroid | 0.267 | 864 |
| Megakaryocyte | 0.243 | 1,855 |
| Mono_DC | 0.234 | 21,331 |
| LMPP_GMP | 0.143 | 3,075 |
| B_Plasma | 0.077 | 10,867 |
| T_NK | 0.036 | 43,708 |

*(malignancy_fpr_by_bin.csv)*

**Why HSC_MPP matters most:** the operational definition of an LSC is a *malignant, high-stemness cell in the HSC/MPP/LMPP bins*. HSC_MPP is therefore the exact node the entire niche hypothesis leans on — and it is the highest-FPR node in the graph. A ~40% false-positive contamination corrupts the `frac_malignant` node feature precisely where the biology lives. The pattern (lymphoid clean, primitive/myeloid dirty) is the signature of inferCNV keying on primitive-cell expression deviation, not clonal biology — so it is a method artifact, not rare pre-leukemic clones.

**Two compounding issues:**
- The blueprint's **≥200-malignant-cell sample gate is not implemented anywhere** (only the ≥500-total-cell gate exists, and that one is satisfied — min 528). 38/130 samples (29%) have <200 malignant cells; 19 have <25; 3853_R has exactly 0. All pass into downstream stages.
- **Timestamp inversion:** per-sample labels are dated 07-14 but ingest/QC were refreshed 07-15, so on-disk labels are keyed to a superseded QC cell batch, and `ref_norm_summary.csv` (needed to reproduce the routing) is gone. Original-run reproducibility is lost (see bug list).

**Verdict: does-not-stand.** The ≥500-cell gate and the study-split machinery work, so this is a substantive *foundation failure*, not a total non-execution. But the malignancy layer the downstream causal claims rest on is not the one the blueprint specified.

---

## 2. Per-phase detail

### Phase 2 — hierarchy projection / node vocabulary — 🟢 stands

- Pure per-sample Symphony projection onto a frozen BoneMarrowMap reference with **no in-cohort re-integration** (as specified). Coverage complete: 220/220 `__bmm_percell` tables across 13 datasets.
- Malignant identity is **genuinely orthogonal** (R2 holds): projection computes only `hierarchy_bin`/`in_ccc_graph`; malignancy enters via an exact barcode join in step 02 (join_rate=1.000 for all 130 joined samples).
- The 7-bin CCC vocabulary is frozen and correctly excludes Stromal (90.3% FPR — this is *why* stroma is excluded, premise 2).
- **Deviations (minor):** Palantir pseudotime / diffusion entropy (the R7 mitigation) not implemented — only stemness scoring exists. A 9th "Unassigned" category (~3.9% of evaluable cells) is created but not in the original 8-bin map. Four AMBIGUOUS bin placements left unresolved in `bmm_bin_map.tsv` (e.g. "Cycling Progenitor"→HSC_MPP folds lineage-agnostic cycling cells into the LSC bin — worth noting given §1). Megakaryocyte bin is fed by a single broad type (thin as a node).
- The HSC_MPP FPR issue lands *on* this phase's key node but is inherited from Phase 1, not caused here.

### Phase 3 — meta-programs / cNMF + SRRS — 🟡 partial

- **What stands:** 10 malignant + 10 normal meta-programs extracted (malignant MP1 is a degenerate single-gene JUN cluster, leaving 9 substantive). A `tumor_specific`/`normal_like` flag exists (5/5 split) and the malignant-vs-normal Jaccard is a usable co-opted-vs-emergent surrogate. Per-cell MP usage scored.
- **What does not stand:** the deliverable that *defines* Phase 3 as a reproducibility screen — the **dual-layer SRRS (Global SRRS ≥ 0.6 hard gate) + M4 conditional evidence-weighting + HSI** — was **never computed**. `grep` for SRRS/HSI/w_cond/n_recovered across `04_cnmf` returns nothing. The retained "robust" programs rest on subjective hand labels (`mp_labels.tsv`) whose hierarchy-specificity percentages are not traceable to any table.
- **Note:** the per-sample cNMF cell floor is `CNMF_MIN_CELLS=100`, a 2× relaxation of the blueprint's R8 ≥200 stability bound — on the exact low-cell instability axis the blueprint flagged. Only 9/13 datasets contribute malignant cells; even if SRRS were computed, |S|=9 with 3 studies contributing ≤3 samples.
- **Notable biology:** no LSC-like/primitive/stemness meta-program survived on the malignant side — the 3 robust_tumor MPs are all monocytic (MP5/MP6/MP10). Consistent with the §1 concern that primitive-bin malignant calls are unreliable.

### Phase 4 — CCC graph construction — 🟡 partial

- **What stands:** a usable single-method CellChat cohort was genuinely built — 220 projected → 159 metadata-eligible → **148 graphs** (edge counts 2–49, mean 25.5; 28/148 flagged sparse). The R9 eligibility gate is correctly implemented (159 = 48 not_l2_capable + 10 sorted_sublibrary + 3 below_min_cells). Node presence is uneven (Megakaryocyte in only 52 samples vs T_NK 193), confirming the unbalanced-FGW design is load-bearing.
- **What does not stand — the critical deviation:** the blueprint's **R3 two-method robustness gate does not exist at all.** No LIANA+ code, tensor, or agreement filter anywhere (`find -iname '*liana*'` empty; `config_ccc.R:45` names a "liana arm 02b" that was never written). **Double inversion:** LIANA+ was the blueprint's *primary* engine and CellChat only the sensitivity control; the delivered pipeline runs the control alone and drops the primary. The composition-permutation null (second half of R3) is also absent. So edges enter the barycenter with no cross-method agreement check and no significance filter.
- M7 functional-axis annotation has no substrate in the delivered tensor (only ligand/receptor kept, no pathway map).

### Phase 5–6 — distance + FGW barycenter — 🟡 partial

- **What stands:** the pipeline runs end-to-end and is directionally sane (healthy HDS < AML HDS). Sparse-graph flag/exclude implemented. Barycenter membership: of 148 graphs, 120 non-sparse; **Healthy barycenter = 20 samples**, **AML barycenter = 93 samples** (113 total carry mass; 35 are scored but contribute none). This is a **small, asymmetric** basis — especially 20 healthy graphs.
- **What does not stand:**
  - **M8 not implemented.** The node-feature vector is still the v1.0 3-dim (`frac_malignant, mean_stemness, n_cells`) that M8 explicitly rejects. Features are **globally z-scored**, not within-sample rank-percentile — which M8 mandates precisely so platform effects are not re-injected after the D2 rank transform removed them from the edges.
  - **Unbalanced FGW not used** — balanced FGW with an eps-mass (1e-6) substitute, despite Phase 6 and R4 both calling for unbalanced.
  - **R4 safeguards absent:** no per-edge barycenter sample-support, no bootstrap CI, no permutation null, no simple-distance benchmark. Uniform per-sample weights with no study balancing → barycenters can be dominated by the largest cohort.
- **Circularity confirmed (feeds H2):** `01_build_fgw_inputs.R:67` forces `frac_malignant := 0` for all healthy samples. Post global-z-score, every healthy node's `frac_malignant` becomes a constant, so an HDS separation on this feature is *imposed by the pipeline, not observed.* Additionally, `:72` imputes NA `frac_malignant` to the global mean — giving **19 unlabelled AML samples a fabricated malignant fraction** (see bug list).

### Phase 7 — emergent edges + patient scoring — 🟡 partial (machinery), 🔴 (result)

The machinery is correctly implemented and admirably self-falsifying (FWL permutation, LOO barycenter, pre-specified alpha sweep, platform control). Its scientific output is strongly negative — see H1/H2 below.

### Phase 8 — robustness & cross-validation — ⚫ not-run

`results/tables/09_robustness/` is empty and `scripts/09_robustness/` does not exist. The three load-bearing pillars are entirely absent:
- **Topology perturbation** (30% edge dropout + 20% node mask, 100 bootstraps, ≥70% recovery) — absent.
- **Sample-level bootstrap edge CIs** — absent.
- **External Discovery-vs-Validation validation** — absent (the 70/30 split is defined at Phase 1 but never used at the barycenter; `fgw_input_index.csv` has no discovery/validation column).

A minority of robustness logic was displaced into 08_scoring (platform control, alpha sensitivity, feature decomposition, raw-vs-rank) — but these are stratification/diagnostic checks, not the blueprint's formal robustness stage. PATCH items M9 (receiver-side NSD validation) and D4 (ATAC-subgroup) scripts do not exist.

---

## 3. Hypothesis verdicts

### H1 — conserved emergent edges — 🔴 null result *(does NOT depend on labels)*

No edge is emergent. `emergent_edges.csv`: `is_emergent=False` for all 49; best single edge HSC_MPP→T_NK has perm_p=0.052 (above 0.05); **min q_bh = 0.865**. The block test is flatter still: all 5 blocks q_bh = 0.927, min block perm_p = 0.341; the pre-specified primitive→immune SCS axis is null (mean_dC=−0.009, perm_p=0.857).

**Honest caveat (and the one fact-check correction):** two isolated edges reach q≈0.034 in *alternative* coordinate systems and only within-dataset — raw log-weight T_NK→B_Plasma (`raw_vs_rank.csv`) and CLR B_Plasma→LMPP_GMP (`clr_redistribution.csv`). These are single edges in a stratified re-parameterization, not concordant emergent edges, and the **rank-transformed** weights (the pipeline's actual C) kill even those — the smallest `rank_C_q_strat` is 0.093 (B_Plasma→LMPP_GMP), still above 0.05. *(A draft of this assessment overstated that floor as 0.493; corrected here — the conclusion is unchanged.)*

**Because H1's edges come straight from CellChat and never touch the malignancy labels, this null is immune to the §1 label-quality problem — the most solid result the project currently has.** The one caveat is method, not labels: it is CellChat-only, and the blueprint's primary CCC engine is LIANA+ (CellChat is the sensitivity control). So H1 must be re-checked under LIANA+/NicheNet before it is a final null. Nothing about the labelling problem can move this result, because it never uses the labels.

### H2 — topology (and direction) carry disease information — 🔴 not topological, with the v1.0 features *(depends on labels)*

The key pre-registered test is the alpha sweep. At the **pure-topology endpoint α=1**: p_global=0.105, p_strat=0.271, beta_global=0.012 — **null.** The signal lives entirely at the feature end (α=0: p_global=0.0009). In v1.1 that α=0 endpoint is no longer just a diagnostic — it is H2's **formal control arm (c)** (the node-feature-only model that topology has to beat). So the plain reading is: with the current features, topology does not beat the feature model. Decomposing the features:

| Feature set (α=0.5) | β (global) | p | Reading |
|---|---|---|---|
| only_frac_mal | 0.173 | <0.001 | **circular** — `frac_malignant` mean_healthy = 0.0 by construction; re-encodes the AML/healthy label |
| all3 | 0.075 | <0.001 | |
| no_frac_mal | 0.064 | 0.0014 | residual after dropping the circular term |
| only_ncells | 0.050 | 0.008 | **technical** — raw cell count |
| only_stemness | −0.015 | 0.56 | the *intended* biological quantity — **null, wrong-signed** |

*(feature_decomposition.csv)*

So whatever discrimination exists is blast fraction (circular) + cell count (technical), **not graph topology, and not stemness.** Two qualifiers matter, though. First, the α=0 model that "wins" is itself the circular `frac_malignant`, so even the control arm is contaminated. Second, this is **not yet the decisive v1.1 test**: the two things that would give topology a fair chance — M8's honest 7-block node features (within-sample rank, no circular term) and M6's direction scores (cosθ, HDS⊥) — are not built. So read this as *"not topological with the un-amended v1.0 features,"* with the real test still pending. (The earlier H2-supporting scripts `01_h2_blast_regression` / `05_h2_platform_control`, dated Jul 21–22, report H2 as *supported*, but they are fixed at α=0.5 and predate the Jul 24 diagnostics — their "support" is exactly the composition artifact the α=1 test exposes.)

### H3 — treatment-pressure trajectory — ⚫ blocked / un-run *(needs the M6 machinery + verified metadata)*

The pre-registered three-model (A/B/C) competition cannot be run yet — and the reason is structural, not just a missing table. Telling A (monotone) from B (V-shape) from C (new direction) needs the **M6 direction scores — HDS∥, HDS⊥, cosθ** — and those are not built (placeholder scripts only), so there is no A/B/C fit. What we *do* have is a shape, measured on stemness: stem_frac and primitive_frac fall Dx→MRD (both p=0.031, n=6, all 6 pairs down), rebound MRD→Relapse (p=0.063, n.s.), and net Dx→Relapse is null (p=0.156 / 0.688) — a **V-shape**. Under v1.0 this looked like "failing the monotone prediction." It is not: v1.1's M6 pre-registers this exact V as **Model B (rebound)**, so the shape is *consistent with a pre-specified alternative*, not a miss. (Only 2 of 4 stemness signatures even show the V — LSC17 and vanGalen_HSC_Prog do; HSPC_core is monotone-down, vanGalen_HSC_like peaks at MRD.)

**Decisively, the axis is patient-confounded and label-dependent:** 62% of all MRD malignant cells (5324 of ~8634) come from a single sample, 6323_MRD, which is 98.7% Mono_DC with primitive_frac=0.006 — so the MRD stem-collapse is essentially one patient whose malignant compartment was labelled almost entirely non-primitive. The comparison rests almost entirely on GSE227903 (the only dataset with clean Dx/MRD/Relapse triplets — 23 of 130 samples carry any timepoint; the other 107 are cross-sectional/baseline by design, not broken metadata). Timepoints in other datasets are unverified against the source papers. **So H3 is blocked on two things, and neither of them is the shape of the trend: (1) the M6 direction machinery has to be built, and (2) the timepoint/pairing metadata has to be reconciled. The stemness V we see is a hint toward Model B — nothing that can be adjudicated yet.**

---

## 4. Deferred bug list (for the post-metadata rerun — NOT to fix this week)

Ordered by severity. None of these change the three verdicts above; they matter for the rerun.

| # | Where | Problem | Severity |
|---|---|---|---|
| 1 | `config_malignancy.R:73` → `ref_norm_summary.csv` (absent); guarded by 44:82 / 91:28 / 90:61 | Routing table missing; any inferCNV rerun-from-disk hard-aborts at `stopifnot`. Rebuildable via `20_refnorm_identify.R` but not identical to the 07-14 run | **blocks rerun** |
| 2 | `10_figures/f00_fig_config.R:18-19` (f01/f04) | Hardcode legacy `01_qc/` & `03_malignancy/` paths + a never-existing `03_qc_report__ALL.csv`; abort at `stopifnot`. (The English `g_assessment_figures.R` is unaffected) | **blocks rerun** |
| 3 | `07_fgw/01_build_fgw_inputs.R:67` | `frac_malignant := 0` for healthy → one of 3 FGW features made deterministic by group label → H2 partly circular | **distorts result** |
| 4 | `07_fgw/01_build_fgw_inputs.R:72` | NA `frac_malignant` imputed to global mean → **19 unlabelled AML samples get a fabricated malignant fraction** feeding the barycenter and 08 tests | **distorts result** |
| 5 | GSE289435 `__consensus_percell.csv` (07-14) vs 5 new `__numbat_percell.csv` (07-16..23) | Allelic evidence produced but never merged; all 12 labels remain `arms=infercnv`. Merging changes labels → must be folded into the joint rerun, not triggered alone | **distorts result** |
| 6 | `side_tp53_pvrl4/{infercnv_proxy_manifest,tp53_locus_proxy_infercnv}.csv` | Orphan outputs; tree-wide grep finds no producing script → provenance unknown | reproducibility |
| 7 | Labels (07-14) predate QC refresh (07-15) | On-disk labels keyed to a superseded QC cell batch; every downstream stage inherits it | reproducibility |
| 8 | `07_fgw/02_fgw_align.py:30-32` | Re-declares `config_fgw.R` node/feature constants inline ("keep in sync") — one-sided edit silently desyncs R inputs from the Python barycenter | reproducibility |
| 9 | `05_ccc/ccc_edge_distance.csv`, `ccc_edge_qc.csv` | Byte-identical stale duplicates of `06_distance/*` (matching md5); no script reads them → stale-read trap after any 06 rerun. Safe to delete | cosmetic |
| 10 | `config_paths.R:42` `PROJ_OBJ_DIR` | Dead constant (`04_bmm_projected` vs real `03_bmm_projected`); only archived scripts use it | cosmetic |
| 11 | `09_robustness/` | Empty dir, no script dir — Phase 8 has no code | cosmetic |

---

## 5. What this means, and the next step

**The single highest-leverage fix is a second, cohort-wide malignancy evidence arm.** Numbat can't be it (needs FASTQ/BAM, only 2–3 datasets qualify). But `42_exprcnv_run.R` (CopyKAT/SCEVAN) reads the per-sample QC `.rds` directly — **no BAM needed** — so it can add an independent expression-CNV arm across the whole cohort, lifting most samples out of `tier=C_single` and directly attacking the HSC_MPP false-positive rate. This is the most realistic approximation of the blueprint's "≥2/3 three-way consensus" under the current data. It is **decided, and scheduled for after the lab meeting** (it needs compute).

**Sequencing (agreed, refined against v1.1):** fix metadata (timepoint/pairing) → add the CopyKAT/SCEVAN arm → re-run `50_consensus` labelling. Re-labelling invalidates everything that *uses* the malignant labels — the malignant-distribution / stemness-by-timepoint tables, the malignant-fraction node features, and the FGW + scoring downstream (Phase 6–8) — **but not the CellChat graphs (Phase 5) or the hierarchy projection (Phase 2)**: the CCC edges come from hierarchy bins and expression, not from the labels, which is exactly why H1 is label-independent. So the CellChat array does **not** need to re-run (v1.1 §11 says the same). Re-run **Phase 6–8 once**, folding in the GSE289435 Numbat merge at the same time. Only then does it make sense to add the pieces the blueprint specified but the pipeline skipped, in order of leverage: M8's honest node features and M6's direction scores (these unlock the real H2 and H3 tests), the two-method CCC gate (LIANA+/NicheNet, Phase 4), the SRRS reproducibility screen (Phase 3), and the full Phase 8 robustness stage.

**For the lab meeting**, the honest and strong message is: *we used the falsification tests most groups skip and found — and localized — a foundation problem before it could produce a false-positive conclusion. The one result that is label-independent (H1) is a clean null. The fix is scoped.*
