# REPORTING AMENDMENT A1-G02 — durable, out-of-sample statistics of the T imprint (WP-A1 gap G02)

Written 2026-10-28 and committed **before** any analysis below is run. Authorised by the user on 2026-10-28 (start the
WP-A1 gap work in the order of `reports/2026-10_wpa1_paper_gap/WPA1_PAPER_GAP_ASSESSMENT_2026-10-28.md` §5, step 1 = G02;
handoff v2.0 §16.1 item 16, §17 S15). Nothing in S5, S8 or S9 changes: no gene set, threshold, gate or registered wording.
This document only (1) replays two registered computations whose per-sample or per-draw values were never saved, with
known answers, and (2) declares, before seeing them, the out-of-sample statistics a WP-A1 manuscript will quote for "the
imprint exists". Later changes are appended as dated sections.

## 1. Why

- The in-sample gap (1.971 = 3.84 healthy-SD, AUC 1.0) is circular: the genes were selected on the same 25 samples.
- The cross-fitted gap (S9 §3.5: gap_oof 1.378 = 1.98 out-of-fold healthy-SD, `ws_s9_yardstick.csv`) is honest, but its 25
  out-of-fold scores exist only as log lines (`00_project/logs/ws_s9_age_3669091.err` lines 21-45, git-ignored, on /FAST),
  and the per-fold gene lists were never written.
- The S5 label permutation (`ws_s5_known_answers.csv`: real 402; median 0, P95 11, max 1,615 over 200 draws) saved no
  draws. From the table alone its Monte Carlo p for >= 402 genes lies between 0.010 and 0.055, and the null is bimodal
  (near-inversions of the labels recreate the contrast). A signed out-of-sample statistic avoids that bimodality.

## 2. Inputs (all on disk, all frozen)

`ws_pb_GSE185381.rds` (`$raw$T`), `results/tables/ws/ws_s5_scores.csv` (the 25 eligible T samples in their S5 order, with
status, hb, gran, sex, subtype_share = CD8 share), `ws_s5_imprint_genes.csv`, `ws_s5_de_T.csv`, `ws_s8_t_imprint_ambient_prone.csv`
(32 genes), `ws_s5_known_answers.csv`, `decontx/GSE185381/ws_s1_ambient_blacklist.csv`, `gencode_GRCh38_gene_order.txt`,
the S9 log above (known answers only), and for section 5: `ws_s8_scores_rep2.csv`, `ws_s4_subtype_composition_GSE227903.csv`,
`ws_s4_subtype_composition_GSE185381.csv`, `meta_GSE185381_v2.3.csv`, the GSE185381 GEO file names and author per-cell
metadata (capture date, `samples`). SEED = 491638 (`config_paths.R`).

## 3. Replays (must reproduce; each has a known answer)

**R1 — S5 T label permutation.** Regenerate the 200 relabellings of `27_s5_imprint.R` for T by replaying its random-number
sequence: `set.seed(SEED)`; T is the first compartment; the only draws before T's permutations are the two
`sample(pool, 200)` calls of the planted-DE arm (pool = tested genes minus the 402 imprint genes, in the order of the main
fit); then 200 × `sample(status)` on the 25-row design in S5 order, each refitted with the S5 design. The 27 code (blacklist,
`design_table`, `fit_de`) is copied verbatim. Known answers: the main fit reproduces the 402 imprint genes exactly; the 200
hit counts reproduce median 0, P95 (type 7) 11.1 (printed 11) and max 1,615. Saved per draw: hit count, the permuted label
vector, how many true AML samples are labelled AML, and for the samples labelled AML: number captured only in 2019, number
in hashed pools, number female, mean CD8 share. Reported: the Monte Carlo p = (1 + #{count >= 402}) / 201. If the replay
does not reproduce, that is recorded and no redraw with another seed is reported as the S5 null.

**R2 — S9 cross-fitted yardstick.** The 25 leave-one-sample-out folds exactly as `33_s9_age_direction.R` section 5 (same
`fit_de`, design `~ status + hb + gran + sex + subtype_share` with sex dropped if a fold has one sex, imprint = FDR < 0.10
and |log2FC| > 0.5 minus the 32 ambient-prone genes, label-free log-CPM of the 25, z with the 24 other samples' mean and SD,
empty half = 0). Known answers: every logged out-of-fold score reproduced within 0.001 and gap_oof 1.378 within 0.001.
Saved: the 25 out-of-fold scores; every fold's gene list with direction; per-gene selection frequency over the 25 folds and
the always-selected core. Reported, descriptive: the pooled leave-one-out AUC of the 25 out-of-fold scores with a 2,000-draw
stratified bootstrap CI, labelled "conditional on the fitted folds" (pooling across folds mixes scales; Parker 2007).

## 4. Declared statistics (new; computed once)

**D1 — leave-pair-out AUC** (Airola et al. 2011). For each of the 15 x 10 = 150 AML-healthy pairs: refit the S5 design on
the other 23 samples (R2's fold rules), take that fold's imprint, z the label-free log-CPM with the 23 samples' mean and SD,
and score the two held-out samples. AUC_LPO = mean over the 150 pairs of 1[S_AML > S_healthy] + 0.5 · 1[equal]. Reported
with the number of pairs whose fold imprint is empty. This is the AUC a manuscript quotes.

**D2 — full-pipeline label-permutation p of gap_oof** (the out-of-sample "imprint exists" test; Ojala & Garriga 2010).
1,000 relabellings of the 25 status labels (covariates stay with their sample), drawn once in one process after
`set.seed(SEED)` as 1,000 calls of `sample(status)` and written to a table before any of them is analysed. Each
relabelling runs R2's complete 25-fold cross-fit with the relabelled status and gives gap_oof_b = mean out-of-fold score of
the samples labelled AML minus that of the samples labelled healthy. **p_D2 = (1 + #{gap_oof_b >= 1.378}) / 1,001**,
one-sided. A fold that fails to fit gives an NA score, is counted, and is left out of that relabelling's means. Secondary,
from the same runs: the permutation p of the pooled leave-one-out AUC. Run as a SLURM array (50 tasks x 20 relabellings).
**Consequence: p_D2 < 0.05 → the manuscript may state that the imprint exists out of sample, with p_D2 and AUC_LPO.
p_D2 >= 0.05 → the existence claim is not supported out of sample; the line stops for discussion (P5) before G01.**

## 5. Durable recomputation of numbers quoted on 2026-10-28 (descriptive)

- Spearman of the S5 T score with the CD8 share (S5's subtype_share) within the 15 AML and within the 10 healthy.
- Spearman of the S8 rep2 score with the kept-T CD8 share in the 7 GSE227903 diagnosis samples, and in the 6 that formed
  the paused S10 pairs. Diagnosis rows only: post-treatment and relapse rows of GSE227903 are never read.
- Capture period of the 25 samples (from GEO file names) and Lasry's `inflammation_group` (author per-cell metadata).

Not here: per-patient cell-bootstrap score uncertainty (moves to G09); anything about composition (G01).

## 6. Rules, code, outputs

- New scripts only: `scripts/ws/35_a1_g02_replay.R` (R1, R2, D1, section 5) and `scripts/ws/36_a1_g02_perm_crossfit.R` (D2:
  `--draw` writes the relabellings, `--task=N` runs one array slice, `--combine` summarises), each with a SLURM wrapper.
  Helpers are copied from 27/33, not sourced. `27_s5_imprint.R` and `33_s9_age_direction.R` are never re-run (they would
  overwrite frozen tables and S9 would read its age relation a second time). No new dependency.
- Outputs, all new, prefix `results/tables/ws/ws_a1_g02_`: `s5_perm_draws.csv`, `oof_scores.csv`, `fold_genes.csv`,
  `gene_selection.csv`, `lpo_pairs.csv`, `perm_labels.csv`, `perm_crossfit.csv`, `descriptive.csv`, `summary.csv`,
  `checks.csv`; a GATES.md block `A1-G02`; the outcome is appended here as a dated section.

## 7. Dated note, 2026-10-28 15:49, written before any D1 or D2 result is read

Written after replay job 3677853 finished (15:38) and while the D2 array 3677854 is still running. No D1 output
(`ws_a1_g02_lpo_pairs.csv`, `ws_a1_g02_summary_replay.csv`) and no D2 output has been opened. User decision 2026-10-28.

1. **Control4's sex.** Curated and author metadata say male; its cells are female (XIST 7,797 vs RPS4Y1 56 UMI;
   `ws_s4_pseudobulk_meta_GSE185381.csv`, `sex_inferred`). S5, S9 and this amendment use the metadata sex. R1 and R2 are
   replays and keep it. **D1 and D2 stay exactly as declared in section 4, with the metadata sex, so they test the imprint
   as S5 defined it; they remain the primary out-of-sample statistics.** A cross-fitted gap and D1 recomputed with Control4
   coded female are declared as a sensitivity in the next amendment (G06 core). They are reported next to D1/D2 whatever
   their values and do not replace them. D2 is not re-run with the corrected sex.
2. **Also known before any D1/D2 result:** five of the ten healthy donors (Control1-4 in pool P02, Control5 in pool P03)
   were captured in hashed pools together with T-ALL/ETP-ALL samples that are not in the GEO deposit (feature references;
   T-ALL/ETP hashtags are the top hashtag of 8.5-22 % of called barcodes per lane; no cell kept by the authors carries
   one). Their ambient RNA was not removed. This does not change D1/D2 or how they are read; it is tested in the next
   amendment.
3. **Order after G02 (user decision 2026-10-28):** the next step is one batch amendment (capture day and pool, including
   the T-ALL pools; the library-blocking rungs of the September/S5 bridge; the decontX, per-library doublet and
   corrected-sex arms; handling-artefact scores; Lasry's label as one comparison with capture period), before the
   composition step (G01). The D2 stop rule (p_D2 >= 0.05 -> stop for discussion) is unchanged.
