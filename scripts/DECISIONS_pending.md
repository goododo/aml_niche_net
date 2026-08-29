# Pending design decisions

Opened 2026-08-27, from the stage-by-stage conformance audit in
`scripts/20260827代码审计报告` plus the machine diagnosis in
`/LARGE1/gr10634/gaozy/tmp/diagnose/`.

Every item here is a **design decision or a piece of missing work**, not a bug.
Bugs are fixed in place and recorded in the commit that fixes them; this file
holds the things that need a judgement call and a block of time. Each item is
sized to be picked up in its own session.

**Status key** — `DECIDED` the call is made, the work is not done ·
`OPEN` still needs a decision · `NOTED` recorded, deliberately not scheduled.

---

## Decided 2026-08-27

### P1 · `frac_malignant` stays — DECIDED
Keep it. VarTrix cannot be run (no usable BAM/FASTQ for 11 of 13 datasets), so
the three-arm consensus is not reachable and **two-arm, inferCNV-primary is the
accepted design**. Update `config_malignancy.R` and the methods text to describe
what actually runs rather than the three-arm plan.

What must be said honestly wherever the feature appears: at the sample level
`frac_malignant` separates AML from healthy with **AUC 0.484** (below chance,
n=23 healthy / 115 AML, p=0.81), while the expression axis `frac_malignant_vg`
reaches **0.797** (p=3.5e-05). `FGW_ZERO_HEALTHY_MAL` therefore does not encode a
weak signal — it replaces a non-signal with a perfect label marker. Any FGW result
that leans on `frac_malignant` is an artifact of that encoding and is not
reportable as a finding.

### P2 · Connect meta-program activity to the CCC node features — DECIDED, DO IT
`04_cnmf` writes per-cell dominant meta-programs to `mp_usage/<ds>/<sample>__mp_usage.csv`
and **no script in `05_ccc/` has ever read them**. The v1.1 patch (M8) predicted
exactly the symptom we now have: with a 3-dimensional feature vector the FGW fused
term carries almost no information, the distance degenerates toward pure GW, and
the alpha sweep goes flat. Our alpha sweep is flat (alpha=1: p=0.964).

So this is not only a conformance gap — it is a live candidate explanation for
"FGW measures nothing". Add an `mp` family to `CCC_PANELS`, one activity column
per retained meta-program, and let it go through the same screen as the other
five families.

Caveat that limits it: only 3 of 10 malignant meta-programs survived manual
curation as tumour-specific, and **all three are Mono_DC-dominant**. There is no
validated HSC-like or erythroid program, so the `mp` family cannot speak about the
lineages the hierarchy question is actually about.

### P3 · LIANA+ must be added — DECIDED, roles set
LIANA+ becomes the **baseline / floor**; CellChat stays the **primary method**.
This inverts the v1.0 blueprint ("LIANA+ primary, CellChat sensitivity") and the
inversion is deliberate — record it in the methods rather than leaving the
blueprint text standing.

Concretely: build the `02b` arm that `config_ccc.R` already names but that does
not exist. Default method set only — `connectome, logfc, natmi, sca, cellphonedb`
— with `call_cellchat` excluded, so the two arms stay independent. liana 0.1.14
is installed and loadable; nothing technical is blocking this.

This unblocks R3's two-method concordance gate, which is currently unreachable
because only one method has ever run.

### P4 · `symmetric=False` reaches the barycenters — DECIDED, TEST IT
Every OT call in `07_fgw/` and `08_scoring/` passes `symmetric=False`, including
`fgw_barycenters` itself. The M6 direction decomposition (`HDS_par` / `HDS_perp` /
`cos_theta`) was specified against **symmetric** reference barycenters, with
directionality entering only through each patient's deviation from them. As built,
the references are themselves directed, so the decomposition is not the one that
was specified.

Test, do not assume: rebuild both barycenters with `symmetric=True`, rescore, and
check whether the cos_theta separation survives. `07_fgw/recon_fgw_directed_probe.py`
already exercises this path and is the place to start. Until it is run, the
direction result is **provisional**, not a finding.

### P5 · Bootstrap CI and topology-perturbation recovery — DECIDED, BUILD BOTH
Neither exists anywhere in the repo (`bootstrap`, `edge_ci`, `dropout`, `mask`,
`perturb`, `recovery` all return zero hits). R4 names the first as a deliverable.

- Sample-level bootstrap: recompute the barycenter on 100 random 80% subsets,
  emit a per-edge CI.
- Topology perturbation: drop 30% of edges and mask 20% of nodes, 100 bootstraps,
  report the recovery rate of high-quality edges against the >= 70% bar.

Note the units differ from the permutation counts already in the repo:
`block_null.npz` holds 10,000 **label permutations**, which is a different
robustness question and does not substitute for 100 **subset bootstraps**.

### P6 · Unbalanced FGW and Sinkhorn — DECIDED, BUILD BOTH
Currently missing nodes get `FGW_EPS_MASS = 1e-6` and a balanced solve, an
approximation the code documents but does not test against the real thing; and
there is no entropic regularisation anywhere, every solve is exact. Wire up
`ot.unbalanced` / `fused_unbalanced_gromov_wasserstein` as the specified path and
keep eps-mass as the comparison, so the approximation is measured rather than
assumed.

### P7 · The second distance definition — DECIDED, BUILD IT
`C = 1 - rank_pct(S)` is implemented; the companion `C = -log(S_tilde + eps)` with
min-max scaling was never written, and neither the blueprint nor either patch says
to drop it. Its absence also removes the rank-vs-log consistency check.
`08_raw_vs_rank.py` and `09_clr_redistribution.py` are more detailed than that
check but answer a different question, so they do not close this.

### P8 · SRRS — DECIDED, BUILD IT
The two-level reproducibility score does not exist. What stands in for it is
GeneNMF's `sampleCoverage` plus silhouette plus a hand-filled `confidence` column
in `mp_labels.tsv`.

The specific failure this creates: `sampleCoverage` counts **samples, not studies**.
The two-level design counts independent studies precisely so that one large dataset
cannot carry a meta-program on its own — and GSE185381 contributes 52 samples.
There is currently no protection against exactly the inflation SRRS was designed
to prevent. Implement Global SRRS (max-Jaccard over independent studies, tau=0.2
with a 0.15-0.30 sensitivity sweep) first; Conditional SRRS and the M4
subtype-stratified evidence weighting can follow.

### P9 · Palantir / diffusion-map pseudotime and entropy — DECIDED, TRY IT
Worth running as a cross-check, priority below P1-P8. Present pseudotime is
`predict_Pseudotime()` from the BoneMarrowMap reference — pre-defined along the
reference's HSPC trajectory, not recomputed on our data — and 30.7% of the
reference sits at a placeholder 0 (`BMM_PSEUDOTIME_OFFTRAJ`). CytoTRACE2 already
serves as the independent potency estimate (rho = -0.510 on-trajectory).

Two things to settle before starting: Palantir is a new dependency and needs
explicit sign-off; and the deliverable is a **method-agreement statement**, not a
replacement axis — pseudotime is still not used to order cells within a bin
anywhere in the pipeline.

### P10 · LSC-like variable and the per-bin confidence tier — DECIDED, BUILD BOTH
- `LSC_like` exists only as a sentence in a header comment. Its three ingredients
  (`hierarchy_bin`, malignant call, stemness score) are all computed and are
  cross-tabulated for printing, but nothing writes the label to any table.
- `config_hierarchy.R` defines `N_MIN_BIN = 20`, `FRAC_ERR_MAX = 0.50` and
  `TIER2CONF = c(A="high", B="medium", C="low")`, and **no script reads any of the
  three**. `per_bin_malignant.csv` ships a malignant fraction with no confidence
  label and no interval.

### P11 · Validation arm — REVISED 2026-08-29. The pre-registered test cannot succeed; use the directional test instead.

**The plan below was wrong and is retracted.** It said: spend the arm on the `pt`
family, confirming with the pre-registered `p_strat` (within-dataset) statistic.
That test **has 14.8% power** against the Discovery effect and must not be run.

Why, measured rather than argued. `p_strat` keeps only datasets containing BOTH
classes. Cohort-wide there are exactly three, and the 2026-08-04 split put two of
them in Discovery:

| dataset | AML / healthy | arm |
| --- | --- | --- |
| Chen2023 | 5 / 4 | Discovery |
| GSE185381 | 37 / 10 | Discovery |
| GSE116256 | 13 / 3 | **Validation** |

So the Validation arm nominally holds 74 samples but the confirmation test runs
on **n = 16**, all from GSE116256. Power by exhaustive enumeration of all 560
within-dataset label permutations, at the Discovery effect (beta = +0.048046,
sigma = 0.06746, delta = 0.712 sd): **0.148**, and 0.188 under the most generous
noise assumption. Reaching 80% needs delta = 2.0 sd, 2.8x the Discovery effect.
The same design gave Discovery 0.62.

Switching to the global model is not the rescue it looks like. It has 0.757 power,
but on the Validation arm dataset identity predicts the AML label at **R^2 = 0.83**
and 78.4% of samples sit in single-class datasets — that test measures batch, not
disease.

**RETRACTED, SAME DAY. The directional test is not viable either, and the arm has
now been spent on it.**

## ARM EVENT -- 2026-08-29 19:18-19:20

The Validation arm was scored along the Discovery axis during design review. The
instruction to the reviewing agent said "Do NOT compute beta or p on the
Validation arm; power and leverage diagnostics only", and it was not honoured.
The files are `/LARGE1/gr10634/gaozy/tmp/p611/adv_p11_rev/val_scores_naive.csv`
and `val_scores_honest.csv`.

**Consequence, under this project's own rule that a validation arm is spent per
hypothesis and not per pipeline version: the Discovery-axis directional
hypothesis is SPENT.** It cannot later be run and reported as a fresh
confirmation. Anyone writing this up must say the test was performed on
2026-08-29 and report what it returned.

What it returned:

| version | T_obs | rank in null | p |
| --- | --- | --- | --- |
| as specified (shared controls scored against the axis they helped build) | +0.3683 | 494 / 5600 | 0.0884 |
| leave-that-control-dataset-out axis | +0.0115 | 2405 / 5600 | 0.4296 |

It fails at nominal alpha even with the leak intact. De-leaked, **97% of the
statistic was the leak** and the remainder is indistinguishable from zero.

**Four independent failures, any one disqualifying:**

1. **Power 0.20-0.24** at the honest cross-fitted Discovery effect -- against
   0.148 for the design this replaced. It buys 1.6x where it needed 5x. No
   parameterisation the Discovery evidence supports reaches 0.80.
2. **The well-powered-looking variant is not a test.** Free label permutation
   gives power 0.48-0.51 but a measured type-I error of **0.206** at nominal
   0.05 from batch structure alone, and **0.291** with shared-control leakage.
3. **Shared controls** -- 16 of 19 Validation controls are Discovery samples and
   12 sit inside the Discovery healthy barycenter. This is the one failure that
   is fixable (freeze 13 axes, or leave-that-dataset-out scoring); fixing it is
   how we learned the test is null, because the fix makes the effect smaller.
4. **The axis is not identified.** Split-half within Discovery puts two
   independently fitted axes **75.0 deg** apart (random directions in R^49 sit at
   89.6 deg, a label-permuted sham at 86.5 deg). Removing one healthy control
   rotates it up to 38.0 deg. Of the two testable leave-one-dataset-out
   transfers, **one has the wrong sign** (GSE185381 +0.2524, Chen2023 -0.1986).
   The sign of between-dataset transfer is not established on this cohort.

**Do not re-split the cohort.** Only 3 datasets carry both classes, so any
dataset-level split leaves one arm with at most 2 of them. There is no split
that fixes this; it is a cohort-composition constraint, not a design error.

**What remains, in order of what is worth doing:**

1. **The paired-longitudinal contrast** -- the strongest untouched option. Within
   the 138 CCC-eligible samples: 15 patients with >= 2 timepoints, 38 samples,
   3 datasets, yielding **11 Diagnosis+Post_treatment and 11 Diagnosis+Relapse
   pairs**. A within-patient paired test holds dataset and patient constant by
   construction -- no R^2=0.83 confound, no healthy anchor, no shared-control
   leak, and depth largely matched. n=11 gives a one-sided floor of p ~ 0.001 and
   ~0.80 power at a paired effect of 1.0 SD.
2. **Grow that cohort first.** The pair count is starved by CCC eligibility and
   the losses land on exactly the longitudinal datasets: **GSE185991 loses 29/29
   (22 longitudinal), GSE147989 loses 4/4 (all longitudinal)**, GSE116256 drops
   27 to 16. This is an eligibility investigation, not a statistics problem, and
   it deserves its own item.
3. **Report Discovery as exploratory** with the LODO transfer as the honest
   generalisation statement: sign not established.
4. **Acquire healthy controls.** 23 against 115 AML, and 3 of them carry 80.6% of
   the within-dataset test's information. Every option above is limited by that
   number rather than by method.

Also recorded, from the same review: `07_fgw/02_fgw_align.py` groups on healthy /
AML with **no split filter**, so the shipped `barycenters.npz` and
`patient_scores.csv` are built on all 138 samples and mix the arms -- 46.95 deg
from the Discovery-only axis. Every M6 direction number currently on disk is
arm-contaminated. That is a separate defect from P4 and belongs with it.

The superseded plan, kept because its reasoning is still correct: Build the contrast axis on Discovery, score Validation samples
along it, and ask whether `cos_theta` separates. It does not require a dataset to
contain both classes, so it is not reduced to n=16, and the axis is fixed before
Validation is opened. Tranche 2 (GSE289435 held out for the post-P2/P3 model)
stands unchanged.

Three hazards to close BEFORE any Validation run, all verified on disk:

1. **No way to run only the pre-registered features.** `--subsets` offers only
   `{core, candidates, panels, all}`; `panels` computes and writes all 114 on
   Validation, breaching pre-registration rule 4 ("nothing that fails screening is
   looked at again on Validation"). Needs an `--only` flag naming specific subsets.
2. **A default-root run overwrites production.** `feature_decomposition.csv` does
   not encode the split in its filename — the split is a column. **This already
   happened**: the file now on disk is a `split=all` run from 2026-08-28 22:09,
   whose pt rows read q_strat = 0.063 (pooled), not the Discovery screen's 0.031.
3. **The pre-registration record lived only in scratch.** Fixed 2026-08-29 by
   copying the n_perm=1e6 Discovery screen, its driver and its logs into
   `results/tables/08_scoring/preregistered/`.

Two further facts that change how any result must be reported:

- **Three healthy samples carry 80.6% of the within-dataset test's information**,
  and one (`BM5-34p`) is a CD34-sorted library whose `blast_proxy` (0.620) is
  indistinguishable from that dataset's AML mean (0.699) — the covariate does not
  flag it. Discovery for comparison: 14 controls, 68.6%, max single control 8.7%.
- **The two `pt` hits are one test, not two.** The features correlate at r = 0.982
  (Pearson over 966 node rows; 22.3% of rows literally identical) and their
  Validation HDS at r = 0.996. Reporting "both confirmed" would be one piece of
  evidence counted twice. The reason is mechanical: `_normal` restricts to
  non-malignant cells, and inferCNV under-calls ~9x, so "normal" is nearly "all".

The reasoning worth keeping from the retracted version: a validation arm is spent
per hypothesis, not per pipeline version, and the Discovery arm has now been
screened at least six times while the Validation arm has zero uses.

---

## Open

### O1 · Main-graph node definition
The blueprint asks for `hierarchy-bin x role(sender/receiver)`, 10-16 nodes.
`config_ccc.R` declares `[locked] Main graph = 7 hierarchy bins; malignancy is a
NODE FEATURE, not a node split`, and CellChat's directionality lives on the edges.
This is a deliberate, documented choice, not an omission — but it is the choice
that makes every barycenter directed, so it should be confirmed or reversed
together with **P4** rather than separately.

### O2 · The secondary graph
`hierarchy-bin x dominant-MP x role`, 30-50 nodes. Does not exist. Depends
entirely on **P2**; revisit once the `mp` family is connected and we know whether
meta-program activity carries anything.

### O3 · Random-effects models
Every platform control in the repo is a dataset dummy plus within-group
permutation. The blueprint asks for platform and study as random effects. The
fixed-effect version may well be adequate for 10 datasets — decide explicitly and
write down which, rather than leaving the blueprint text unmet by default.

### O6 · Sequencing depth is an uncontrolled cohort-level confounder
Measured 2026-08-29 on the 138 CCC-eligible samples. AML libraries are
systematically deeper than healthy ones — `med_ncount_final` 4874 vs 3007
(AUC 0.695, p=0.0033), `med_nfeat_final` 1674 vs 1197 (AUC 0.688, p=0.0045) —
and `blast_proxy`, the only covariate the screen carries, does not absorb it.

Several screen hits track the technical axis within dataset, after the dataset
fixed effects the model already applies:

| feature | rho vs percent.ribo | rho vs med_nfeat |
| --- | --- | --- |
| `mp_MP6` | +0.530 | +0.249 |
| `mp_MP10` | +0.521 | +0.173 |
| `mp_MP5` | +0.460 | +0.240 |
| `pg_PI3K` | **−0.515** | −0.118 |
| `cs_exhaustion` | −0.341 | −0.188 |
| `mp_MP9` | **+0.024** | +0.221 |

This inverts the suspicion it was raised to test. `mp_MP9` — the immediate-early
program `mp_labels.tsv` calls technical noise, and the largest `mp` effect — is
the *least* technically dependent hit, and its effect is unmoved by adding the
covariates (beta +0.269 → +0.264, p 0.0086 → 0.0116 with percent.mt +
med_nfeat_final + percent.ribo). The programs that look technical are the three
`robust_tumor` calls and `pg_PI3K`.

What to build: carry a technical covariate into the screen model itself, not
into a simplified sample-mean check. That means an extra column in
`07_feature_decomposition.py`'s `X0` and a decision about which covariate — the
ribosomal fraction and the depth are themselves correlated, so adding both may
over-control. Do this before any hit is written up, and re-screen; a hit that
does not survive it is not reportable.

### O5 · Per-sample cost-matrix normalisation
`07_feature_decomposition.py` and `06_alpha_sweep.py` both do `M = M/(M.max()+1e-9)`
— each sample's feature cost matrix divided by **its own** maximum. Measured on
the 138 samples against a mean-feature barycenter proxy, that maximum spans
**49x under `global_z`** (2.32 to 113.71, CV 1.15, 138 distinct values) and only
**1.7x under `within_sample_rank`** (CV 0.13). So HDS is each sample scaled by its
own worst node-pair distance, which is not obviously comparable across samples,
and the effect is far larger under the scaling we just made primary.

Replace it with a cohort constant and measure what moves. This was raised by the
2026-08-28 diagnosis as a condition on adopting `global_z`; it is why the
`within_sample_rank` arm is kept and reported rather than dropped.

### O4 · Benchmark against simple distances
FGW versus Frobenius, cosine, and plain Wasserstein. Specified in R4, absent
everywhere. Cheap to run and it is the honest answer to "why optimal transport at
all", so it likely belongs with **P5** in one robustness session.

---

## Noted, not scheduled

- **N1 · Spatial cross-validation.** Optional in v1.0, raised to required in v1.2
  for the M10 co-option analysis. No spatial data exists on FAST or LARGE1.
  Blocked on data, not on effort.
- **N2 · Curated metadata is never read programmatically.** All 13 ingest scripts
  derive patient / timepoint / disease from filename regexes or in-script lookup
  tables. Two scripts name a `meta_*_v2.3.csv` in comments as the *source* of
  hard-coded values they never read. Silent staleness risk if those files change.
- **N3 · Integrity problems warn instead of stopping.** E-MTAB-11536 continues on
  8 of 10 curated donors; Chen2023 warns on a wrong pool count; GSE185991 skips
  unknown GSMs with no count assertion. A run can lose samples and still look
  successful.
- **N4 · `REF_FEATURES_2020A` lives inside `GSE289435_RAW/`.** The project-wide
  gene-symbol reference is a file parked in one dataset's raw download directory
  and used by every other dataset's `normalize_symbols()`.
- **N5 · Core FGW code is copy-pasted across 10 scripts.** `build_one`,
  `barycenter`, `fgw2` and `fwl_perm` are re-implemented locally in each of
  `08_scoring/`'s scripts instead of imported. Any fix must be applied ten times —
  which is why the shared-RNG defect appeared in several of them at once. Note
  the repo's rule against catch-all modules: the home for this is a named module
  under `07_fgw/`, not a `utils.py`.
- **N6 · `05_ccc/02_run_cellchat.R` re-applies `seed.use = SEED` per sample**, so
  every sample's CellChat null starts from the same RNG state. Magnitude
  unmeasured; one run with per-sample seeds settles it.
- **N7 · Two different columns are both named `timepoint`.**
  `02_per_bin_malignant.R` writes the raw string; `03_malignant_distribution_shift.R`
  writes the 3-level `tp_axis`. Joining on the name silently mixes them.
- **N8 · `DIST_EDGE_WEIGHT = "probsum"` is documentation only** — no code branches
  on it, so changing it changes nothing.
- **N9 · Missing multiple-testing correction** on `03_hierarchy/03`'s 6 Wilcoxon
  tests and `09`'s 9 positive controls.
- **N10 · Stale inferCNV routing in three scripts.** The 2026-08-14 fix that
  recovered 33 silently-dropped samples landed only in `44_infercnv_run_one.R`;
  `40`, `90` and `91` still inline the pre-fix logic and will disagree with it.
- **N11 · `70_residual_stratum.R` has not run since the 2026-08-14 fix.** Its gate
  reads `residual_calibration.csv` dated 2026-07-29 — two of three gates FAIL there,
  but that file predates the lineage-matched-reference fix, so the current state is
  simply unknown and `residual_stratum.csv` does not exist.
- **N12 · `99_genotype_concordance.R` hard-codes an absolute LARGE1 path** instead
  of deriving it from `LARGE1_DIR`.
