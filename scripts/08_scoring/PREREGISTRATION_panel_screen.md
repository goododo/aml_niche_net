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
