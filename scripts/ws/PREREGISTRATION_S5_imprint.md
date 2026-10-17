# PRE-REGISTRATION — S5: the diagnosis-time imprint on non-malignant cells (WP-A1)

Draft 2026-10-17, written before any S5 code exists and before any AML-vs-healthy differential test has
been run on the S4 pseudobulks. Frozen inputs: `config_ws.R` sections 0, 10-12; the S1 tables under
`results/tables/ws/decontx/`; the S3 tables `ws_s3_*`; the S4 pseudobulks under
`/LARGE1/.../ws_pseudobulk/<dataset>/`. Handoff: `HANDOFF_within_sample_v2.0.md` §5.3 S5-S6, §5.4, §16.

Status: **BINDING from 2026-10-17** (user approval the same day; tag `prereg-s5-v1`). From this tag on,
changes are new dated sections appended at the end, never edits of the text above them.

---

## 1. Hypotheses (from v2.0 §5.1, unchanged)

- **H-A1a** At diagnosis, the non-malignant T, NK, B and monocyte compartments of adult AML patients differ
  in expression from the same compartments of healthy donors, within compartment.
- **H-A1b** Within a patient, the strength of that difference is coordinated across compartments after
  composition, blast burden, ambient RNA and doublets are removed.
- **H-A1c** Direction (and, where testable, coordination) replicates in another study.

What this file does NOT claim: that the imprint is caused by AML (age is confounded, §3.4), that monocytes
are non-malignant (Petti 2019: monocytes are often clone members; the Mono compartment is reported as
"mixed"), or anything at the single-cell level that is not first summarised within the patient.

## 2. Data and units

| item | value | source |
|---|---|---|
| defining study | GSE185381, adult AML at diagnosis (role AML_Dx, not paediatric) vs healthy donors | S4 meta |
| n | 20 adult AML, 10 healthy | S3 |
| compartments | T (CD4 + CD8), NK, B, Mono (CD14 + CD16) | `s4_cell_groups.tsv`, section 12 |
| per-sample entry | >= 20 cells whose marker typing agrees with the projection, not high_error | section 12 |
| eligible adults per compartment | T 15, NK 12, B 13, Mono 16 (all four kept: each >= 12, 16.1 item 14) | S4 |
| eligible healthy | T 10, NK 9, B 10, Mono 10 | S4 |
| unit of every test | the patient. Cells never enter a test unsummarised | P3 |
| paediatric GSE185381 (22) | second cohort; same pipeline, run AFTER the adult imprint is frozen; not pooled | 16.1 item 3 |
| replication | GSE289435 (10 AML, no healthy): direction and structure only (S8) | section 12 |
| description only | GSE227903, GSE116256 | section 12 |

## 3. Differential expression per compartment (defines the imprint)

3.1 **Tool.** edgeR quasi-likelihood (`glmQLFit`, robust = TRUE) on the raw pseudobulk sums; TMM
normalisation; genes kept by `filterByExpr` with the AML/healthy factor, and additionally CPM > 1 in at
least half of the compartment's samples (v2.0 S4 item 5). Decided for 16.2 item 4: edgeR QL, not limma-voom
(small n, count data); limma-voom is run once as a sensitivity arm and reported, never used to pick genes.

3.2 **Design.** `~ status + hb_share_nonery + gran_share_lymph + sex + subtype_share`, where
- `status` = AML vs healthy (the coefficient of interest);
- `hb_share_nonery`, `gran_share_lymph` = the two raw ambient shares (16.1 item 13), centred;
- `sex` = metadata where present, else inferred from XIST/RPS4Y1 (S4: 20/20 agreement);
- `subtype_share` = one within-compartment composition term: T: CD8 share; Mono: CD16 share; NK: CD56high
  share; B: immature-B share. This replaces the fixed-weight subtype composite of v2.0 S4 item 3, which is
  under-powered at n 12-16 per compartment; subtype-level DE (CD4 T, CD8 T, CD14 Mono) is a sensitivity arm.
- NOT in the design: age (§3.4), doublet mean (no T3a scores for controls; added as a sensitivity arm once
  T3a has run on the 10 control libraries), channel/pool (29 of 42 AML and 9 of 10 controls are pooled;
  `multiplexed` is tested as a sensitivity covariate), clinical blast % (missing for half).

3.3 **Imprint genes.** FDR < 0.10 (Benjamini-Hochberg within compartment) and |log2FC| > 0.5, up and down
kept separately. Removed beforehand, not after: HB genes; MPO, ELANE, AZU1, PRTN3; the S1 ambient
blacklist (genes in the top-100-by-removed-fraction of >= 5 GSE185381 channels: 79 genes,
`ws_s1_ambient_blacklist.csv`); XIST and Y-linked genes; MT-, RPL/RPS genes. Immediate-early and heat-shock
genes (`WS_B0_IEG`) stay in but are flagged; a sensitivity arm removes them.
**Minimum imprint:** a compartment with fewer than 10 genes (up + down) has no imprint, and H-A1a is
recorded as not supported for that compartment. No threshold is relaxed to reach 10.

3.4 **Age.** Healthy donors are 19-26 (6) and 39-53 (4) years old; 16 of 20 adult AML patients are >= 60.
Age cannot be a covariate in this contrast. Every AML-vs-healthy statement from S5 is therefore an
"AML-or-age" statement. Two things are pre-specified to bound it: (a) a descriptive comparison of the 4 AML
patients < 60 against the 4 healthy donors aged 39-53 (direction concordance of the imprint genes; no test);
(b) a literature-derived ageing signature of T cells / monocytes, if one is adopted, is applied as a
sensitivity covariate and declared before use. Resolving the confound needs age-matched healthy marrow
(see memory `healthy-marrow-controls-with-ages`); that is a separate decision.

## 4. Scores

4.1 **Patient-compartment score (primary).** log-CPM of the compartment pseudobulk, each gene z-scored
across all eligible samples of the study (AML + healthy); score = mean z(up) − mean z(down).

4.2 **Single-cell layer (secondary, summarised per patient; S5b).** UCell scores of the up and down sets
per cell; per patient and compartment: the mean, and the fraction of cells above the 95th percentile of
the healthy donors' cells of the same fine subtype. These describe whether the imprint is a uniform shift
or a sub-population, and never replace 4.1 in a gate.

## 5. Coordination across compartments (H-A1b)

5.1 Matrix: adult AML patients × compartments (score 4.1), compartments with >= 20 cells only.
5.2 Residualise each compartment's score on: progenitor fraction (blast proxy), the 7-bin composition
(`prop7`, first two principal components), `hb_share_nonery`, `gran_share_lymph`.
5.3 Statistic: mean of the pairwise Spearman correlations between compartments, each pair over its complete
patients (n per pair 10-14; T+NK+B complete in 11, used as a sensitivity arm).
5.4 Null 1: permute each compartment's residual scores across patients independently, 1,000 times.
5.5 Null 2: 1,000 random gene sets per compartment matched to the imprint on mean log-CPM and detection
(bins), scored and residualised the same way; the real coordination must exceed their 95th percentile.
5.6 CI: patient bootstrap, 1,000 draws.

## 6. Same-channel control (S6; GATE_A1 items 3-4)

Three channels hold an adult AML donor with a healthy donor (7 AML, 3 healthy; `ws_s3_gse185381_channels.csv`).
6.1 Within each channel, for each compartment with >= 20 cells on both sides: sign concordance of the
imprint genes' AML-minus-healthy log-fold-changes with the study-wide imprint; compared with the same for
1,000 random gene sets.
6.2 Healthy donors in mixed channels (3) vs healthy donors in healthy-only channels (7): difference in
imprint score, as a fraction of the AML-minus-healthy difference.
Three channels is exactly the v2.0 minimum; the result counts in GATE_A1 but is reported with its n.

## 7. Gates (numbers fixed here, for 16.2 item 2)

| gate | criterion | passes if |
|---|---|---|
| GATE_A1 (1) | partial Spearman of the compartment score with `hb_share_nonery` and with `gran_share_lymph`, given status | both |rho| < 0.30 in every compartment with an imprint |
| GATE_A1 (2) | imprint genes vs the extended ambient ranking (top-100-by-removed-fraction in >= 2 channels) | hypergeometric p > 0.05 |
| GATE_A1 (3) | S6 sign concordance (6.1) | above the random-set 95th percentile in >= 2 of 3 channels, pooled over compartments |
| GATE_A1 (4) | S6 healthy-in-mixed-channel shift (6.2) | < one third of the AML-healthy difference |
| GATE_A2 | coordination (5.3) | > 0, bootstrap 95% CI excludes 0, above random-set P95 (5.5) |
| GATE_A3 | S8, GSE289435: sign concordance of imprint genes (AML vs the GSE185381 healthy reference is NOT available there; direction uses the within-study rank of AML patients' scores vs their coordination) | coordination > 0 and above its random-set P95; direction test as specified in the S8 addendum |

Any gate failing stops the line for discussion (P5). A failed GATE_A1 (1) for one compartment removes that
compartment from WP-A2, it does not relax the threshold.

## 8. Known answers the S5 code must pass before the real contrast is read

1. Planted DE: 200 random genes multiplied by 2^(±1) in the AML samples of one compartment; recovery at
   FDR 0.10 reported as the compartment's power statement (as `12_pseudobulk_de/power.csv` did).
2. Sex contrast: with `sex` as the coefficient of interest, XIST (female) and RPS4Y1/DDX3Y (male) must be
   the top hits in every compartment.
3. Label permutation: 200 permutations of AML/healthy; the number of imprint genes under permutation gives
   the empirical false-discovery count; the real count must exceed its 95th percentile for the
   compartment to carry an imprint (this is in addition to the FDR threshold).
4. Random-gene-set coordination null calibrated: at r = 0 the permutation p is uniform (checked on
   simulated scores).

## 9. Outputs

`results/tables/ws/ws_s5_de_<compartment>.csv` (all genes), `ws_s5_imprint_genes.csv` (frozen list, tagged),
`ws_s5_scores.csv` (patient × compartment, 4.1 and 4.2 summaries), `ws_s5_coordination.csv`,
`ws_s5_same_channel.csv`, `ws_s5_gates.csv`, `ws_s5_checks.csv`; GATES.md block "A1/A2".

## 10. Not done, by decision

No gene set is changed after the result is seen; no sample is added or dropped after the result is seen;
no cell-level p-value is reported; the paediatric cohort and the replication study are not opened until
`ws_s5_imprint_genes.csv` is tagged.

---

## 11. Outcome and decision (appended 2026-10-17, after the run; nothing above this line was edited)

Run: `27_s5_imprint.R` (job 3649941) and `28_s5_coordination.R` (job 3650033); tables `results/tables/ws/ws_s5_*`.
- Imprints: T 402 genes, Mono 275, NK 91; B 9 (not called). Known answers passed (section 8).
- GATE_A1: (1) pass, (3) pass, (4) pass; **(2) FAIL** (T 32 imprint genes in the extended ambient ranking,
  p 0.035; Mono 31, p < 0.001; NK 8, p 0.11).
- **GATE_A2 FAIL**: coordination 0.123, permutation p 0.26, random-set P95 0.45, bootstrap CI -0.35 to 0.70.
- Sensitivity: a multiplexed-library covariate removes the NK imprint (91 -> 2) and leaves T intact (402 -> 405);
  the Mono imprint contains HOXB2, HOXB3 and NKG7 (clonal or doublet admixture).

Decision (user, 2026-10-17): the "coordinated imprint" claim is dropped. The claim is narrowed to: adult AML
T cells differ at diagnosis from (younger) healthy donors; the difference is robust to the pre-registered
controls except that 32 of its 402 genes are ambient-prone, which is acknowledged; replication in other
studies is pending. WP-A2 may ask the remission question with the T imprint only. S7 (DIALOGUE) proceeds
unsupervised; if it finds no cross-compartment program, H-A1b is closed. Any further test on the T imprint
(replication, remission) gets its own pre-registered section before it is run.

## 12. S7 DIALOGUE settings (appended 2026-10-17, before S7 is run; nothing above this line was edited)

Purpose (handoff 16.1 item 15): an unsupervised look for cross-compartment programs. If none is found, H-A1b is
closed; if one is found, it gets its own pre-registered test before any claim.

- Tool: DIALOGUE 1.0 (`DIALOGUE.run`, PMD on sample-level averaged PCs, then mixed-effects gene calls).
- Cells: the S4 kept cells (concordant, not high_error, map group4) of GSE185381 adult AML patients and healthy
  donors; a sample enters a compartment with >= 20 kept cells (section 12 rule).
- Primary run: compartments T, NK, B, Mono on the samples eligible in all four. Secondary run: T and Mono on the
  samples eligible in both (the two compartments whose imprints were stable in S5).
- Per compartment: log-normalised expression (counts per 10k, log1p), 2,000 most variable genes within the
  compartment, 15 PCs (`X`); `cellQ` = log10 UMIs per cell; sample-level covariates in `covar`: cellQ, the two
  raw ambient shares. k = 3 programs; DIALOGUE defaults otherwise (`n.genes` 200, `p.anova` 0.05).
- Reading rule, fixed here: a "cross-compartment program" exists if some MCP has (a) DIALOGUE's cross-cell-type
  correlation R positive in every compartment pair of the run with its mixed-effects p < 0.05, (b) that R above
  the 95th percentile of 10 runs in which sample labels are shuffled independently within each compartment
  (breaks cross-compartment pairing, keeps everything else), and (c) stability: median Jaccard >= 0.5 of its
  top-50 genes per compartment over 10 runs on random 80% subsets of samples. Programs failing any of the three
  are reported as not found. Association of MCP sample scores with status and age group is reported, not tested.
- Overlap of MCP genes with the S5 T imprint is reported (hypergeometric), not used in the reading rule.

**Superseded (2026-10-17, same day, before any run under this section produced a result):** S7 is governed by
`PREREGISTRATION_S7_dialogue.md` (tag `prereg-s7-v1`), written in a parallel session with a fuller design (two
runs, a planted known answer, 20 shuffles, 20 resamples, 500-cell cap). This section and `29_s7_dialogue.R`
are kept as a record and are not run.

## 13. DRAFT — the T imprint beyond GSE185381 (appended 2026-10-17; NOT binding until the user approves and it is tagged)

Acknowledgement carried into every statement below: 32 of the 402 T-imprint genes are in the extended ambient
ranking (GATE_A1 (2) failed); the imprint is "AML-or-age" (section 3.4); B carries no imprint; NK's imprint is
fragile; Mono is clone-mixed. Only the T imprint (`ws_s5_imprint_genes.csv`, compartment T, tag to be set) is used.

13.1 **Coherence replication in GSE289435** (10 AML samples with >= 20 kept T cells; no healthy donors there, so
direction against healthy cannot be tested). Statistic: within the T pseudobulks, the mean pairwise correlation
of the up genes plus that of the down genes, minus the mean up-vs-down correlation ("program coherence").
Null: 1,000 expression-matched random gene sets of the same sizes. Passes if the imprint's coherence exceeds the
random P95. Secondary, reported: the same in GSE227903's 7 diagnosis samples with >= 20 T cells.

13.2 **Remission (GSE227903, 6 Dx -> MRD pairs with >= 20 kept T cells at both times: 1216, 1886, 3904, 40389,
4978, 5750; 6323 has 12 T cells at diagnosis and is excluded before any result is seen).**
- Score: log-CPM of the T pseudobulk, genes z-scored across the 13 GSE227903 T pseudobulks (Dx + MRD), score =
  mean z(up) − mean z(down); covariates hb and gran shares regressed out before testing.
- Test 1, rank preservation: Spearman of the six patients' Dx vs MRD scores, one-sided (positive). With n = 6,
  p < 0.05 needs rho >= 0.829: this is a direction statement, as v2.0 6.4 already said for adults.
- Test 2, self-recognition: for each MRD sample, correlation of its imprint-gene z-profile with each of the six Dx
  profiles; statistic = mean rank of the own diagnosis sample (1 = best). Nulls: all 720 pairings; and 1,000
  matched random gene sets. Counts only if it beats both.
- Reported: MRD score relative to the Dx scores (no healthy reference in this study); batch (all six pairs are
  within one batch after 6323's exclusion).
13.3 **Not done:** no gene is added or removed from the T imprint; no sample is added after a result; no paediatric
run until the adult statements above are written down with their outcome.
