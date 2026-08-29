# Pre-registration: node-feature panel screen

Written and committed **2026-08-26, before the screening run produced any
result**. That is the only thing that makes it a pre-registration rather than a
description.

## What is being screened

114 candidate node features, produced by `05_ccc/03_node_features.R` from
`CCC_PANELS`, in five families:

| family | n | question |
| --- | --- | --- |
| `st` | 12 | does the stemness result depend on which signature is used? |
| `pg` | 42 | do PROGENy pathway activities separate AML from healthy? |
| `cs` | 39 | do immune-tone / retention / senescence programs separate them? |
| `mt` | 18 | do metabolic programs or the venetoclax axis separate them? |
| `pt` | 3 | does BMM pseudotime (on-trajectory cells only) separate them? |

Each is tested alone, at alpha = 0, in the within-dataset model with
`blast_proxy` as a covariate and the label permuted within dataset.

## The arms

Screening runs on **Discovery only**: 57 AML + 20 healthy, from Chen2023,
GSE185381, GSE239721, Petti2019. Confirmation runs on **Validation**: 55 AML
from GSE116256, GSE201966, GSE227903, GSE289435.

Only 23 healthy samples exist cohort-wide, so **the healthy controls are shared
between the arms**. This validates against AML heterogeneity, not against
control heterogeneity. It is a limitation of the cohort, not a choice, and it
goes in the Methods.

## Decision rule, fixed in advance

1. A feature **passes screening** if `q_strat < 0.05` after BH-FDR **within its
   family** on the Discovery arm. Families are corrected separately because they
   are five separate questions; pooling them would spend the correction on
   unrelated hypotheses.
2. Every feature that passes is carried to Validation with its **sign fixed** to
   whatever Discovery gave. A Validation result with the opposite sign is a
   failure, not a finding.
3. A feature **confirms** if the Validation estimate has the pre-registered sign
   and raw `p_strat < 0.05`. No FDR on the confirmation step: the hypotheses were
   fixed before the arm was opened, so there is nothing left to correct for.
4. Nothing that fails screening is looked at again on Validation. If a family
   yields nothing, that is the result for that family.
5. If a feature passes screening but fails confirmation, it is reported as
   failed. It does not get re-tested with a different covariate set, a different
   stratum, or the pooled cohort.

## What this screen cannot settle, and must not be used to claim

- **Edge weight tracks node cell count.** `cor(weight_probsum, min node size)
  = +0.688`; edges between nodes in the bottom size quintile are present 1% of
  the time versus 94% in the top. `population.size = FALSE` keeps abundance out
  of the probability, but detection power puts it back: with ~10 cells the
  triMean is unstable and `nboot = 100` cannot reach p < 0.05. Any topological
  readout is therefore partly a readout of composition, and the block-flow lead
  (`imm_to_imm`, p = 0.035 after a sample-level immune-fraction control) is
  **not** cleared by that control, which is coarser than the mechanism.
- **The malignant / non-malignant split is unreliable.** inferCNV under-calls
  ~9x and is uncorrelated with clinical blast burden. The `_normal` and
  `_malignant` strata inherit that. `05_ccc/05_undercall_contamination.R` shows
  the stemness contrast survives the split being wrong; nothing shows that for
  any other feature, and each one that passes needs the same test.
- **MSC / stromal biology.** Median 0 stromal cells per sample; only 4 healthy
  samples carry >= 20. No claim about mesenchymal stroma is supportable here,
  whatever `cs_senescence` does in the haematopoietic compartment.
- **Anything about treatment response.** The field does not exist in the
  metadata.

## What has already been used

Every analysis before 2026-08-26 pooled all 138 samples, so the Validation arm
has been used for **estimation** — the alpha sweep, the per-edge regression, the
planted-effect harness. It was never used for **selection**: `FGW_FEATURES` was
fixed a priori in `config_fgw.R`. Existing results are therefore not overfitted,
but the Validation arm is no longer naive, and this pre-registration is the first
time it is being held back on purpose.

The consequence: a feature confirming here is evidence against overfitting to
the Discovery arm, not evidence from a fully naive cohort. Say it that way.

---

## Amendment 2026-08-29 — the technical-covariate gate, declared before the run

Two things changed the screen after the original registration, and one confounder
was found that invalidates part of it. The rule below is written **before** the
re-screen is run, so it cannot be chosen after seeing which hits survive.

### What changed

- A sixth family, `mp` (9 cNMF meta-programs x 3 strata = 27 features), joined the
  five originally registered. The screen is 141 features, not 114. Both versions
  are frozen under `results/tables/08_scoring/preregistered/` with input
  checksums. The v141 screen returns 18 hits: cs 3, mp 11, pg 2, pt 2, mt 0, st 0.
- Meta-program activity is now scored on all cells in the 7 CCC bins rather than
  on the inferCNV-malignant subset.

### The confounder (O6)

AML libraries are systematically deeper than healthy ones: `med_ncount_final`
4874 vs 3006 (AUC 0.695, p=0.0033). On the **Discovery arm**, after the dataset
fixed effects the model already applies, depth still tracks the label at
r = +0.284 (p = 0.012), and 68.2% of the label variance survives dataset
adjustment. Several hits track `percent.ribo` within dataset: `mp_MP6` +0.53,
`mp_MP10` +0.52, `pg_PI3K` -0.52, `cs_exhaustion` -0.34, while `mp_MP9` -- the
largest `mp` effect -- does not (+0.02).

`blast_proxy` is the only covariate the screen carries and it does not absorb
this.

### The rule, fixed in advance

1. The re-screen runs three arms, all `--split discovery --subsets panels
   --alphas 0 --n_perm 1000000`, differing only in flags:
   - **baseline** `--covar none --absent_mask off`
   - **depth**    `--covar depth --absent_mask off`
   - **both**     `--covar depth --absent_mask on`
2. **The baseline arm is a gate, not a result.** It must return 18 hits with the
   same family split and `beta` matching the frozen v141 screen to 0.0e+00. If it
   does not, the edit has changed the default path and arms 2 and 3 are
   uninterpretable. Stop and fix.
3. **A hit that loses q < 0.05 under `--covar depth` is NOT REPORTABLE.** It is
   recorded as depth-dependent. This is decided now, before the numbers exist.
4. A hit that survives `--covar depth` but changes under `--absent_mask on` is
   flagged, not dropped: the mask fixes a real defect (the cost matrix is
   normalised by a maximum taken over absent nodes carrying ~1e-10 mass, which
   moves HDS by up to 12.8% on 5 of 138 samples and hits healthy samples more
   often than AML), so the masked arm is the more correct one and any change is
   reported as such.
5. No hit is re-tested with a different covariate set to recover significance. The
   three arms above are the whole analysis.

### What this amendment cannot fix

Depth and ribosomal fraction are themselves correlated, so controlling for both
may over-control and remove real biology along with the artifact. The depth arm
is the primary; `percent.ribo` is not added to it for that reason. This is a
limitation, not a solved problem.
