# PREREGISTRATION A1-G03: S5's T imprint and the September T_NK null, one factor at a time

**Status: v1, APPROVED by the user 2026-11-05 (all section-9 items as recommended), after reviews 1 and 2 (section
10); tag `prereg-a1-g03-v1`.** Nothing declared below was computed before the tag.

**Sources.** It is written from two read-only inventories, kept outside the repo:
`/FAST/gr10634/gaozy/reports/2026-11_a1_g03_prereg/inventory/A_september_spec.md` (A) and
`B_s5_bridge_obligations.md` (B). No September per-gene value was read; 1.4 lists what was seen.

**Place in the sequence.** Order of the WP-A1 steps (user decisions 2026-10-28 to 2026-11-05): A1-G02 (PASS) →
A1-BATCH (PASS with qualifications) → A1-G01 → A1-G04 ((a) "not shown", (b) "not readable") → **this**. Handoff v2.0
§17 S15; G04 handover (`PREREGISTRATION_A1_G04.md` sections 5.7 and 7).

## 0. Question, scope, and what each outcome lets the paper say

Two analyses of largely the same samples disagree:

| | September (2026-09-16) | S5 |
|---|---|---|
| source | `scripts/12_pseudobulk_de/PREREGISTRATION_pseudobulk_de.md`, commit a66cbd0 | `PREREGISTRATION_S5_imprint.md` |
| design | pooled pseudobulk DE on the T_NK bin | per-compartment pseudobulk DE |
| samples | 39 AML vs 14 healthy | 25 GSE185381 adults |
| result | no reportable hit beyond its permutation nulls | 402 genes; frozen primary list P370 = 156 up / 214 down |

24 of S5's 25 samples are in the September fit, so a reviewer will read the pair as a forking path.

**G03 asks: which design factors separate the two results?** It answers in three parts:
- **The ladder (section 3).** Start from S5's design and change one factor per rung until the design is September's.
  R0 is S5 on its 25 samples; R1 to R6 run on the 24 shared adults.
- **Two extension rungs fixed by G04 5.7.** They add the children, then Chen2023.
- **A recomputed prediction.** The prediction on record (1.1) is recomputed with G04's estimates.

**What G03 is not.**
- **It does not test again whether the imprint exists.** A1-G02 tested that out of sample (gap 1.378, D2 p 0.008), and
  A1-BATCH G2 tested it against lane structure (p 0.017). Both use a sample-level statistic; G03 counts genes.
- **It does not reinterpret the September result.** That result stays reported as a null for |logFC| ≥ 1.5 per gene.
  No G03 result is described as showing that the September analysis was underpowered, used the wrong model, or is
  explained by its pooled samples (section 7, "never written"). This repeats the September promise
  (`PREREGISTRATION_pseudobulk_de.md`:41-44) and G04 1.1's binding wording.

**What each ladder outcome lets the paper say.** The categories and wording are fixed in section 7; the classes are
decision 9.7.
- **"held through R6".** S5's count stays beyond both of its permutation nulls (sample and library) at every rung, up
  to the September cell definition, engine with library block, and gene universe. It is reported with R6's September-rule count beside it.
  The paper shows the ladder table and the September null side by side.
- **"lost at a specification rung"** (R2 covariates, R3 T+NK, R4a engine, R5 cell definition, R6 universe). The
  paper names the S5 choice the count needs. The per-gene list keeps its status with that qualifier.
- **"lost at a validity rung"** (R1 AML0310, R4b library block). The paper calls the per-gene list exploratory. The
  sample-level results of A1-G02 and A1-BATCH G2 are a different statistic and are not changed.
- **"borderline"**, or **"S5's count is not beyond its own null at R0"**. These are reported as such (section 7).

### 0.1 What G03 owes, where it is handled, and declared departures

| source | obligation | here |
|---|---|---|
| assessment 2.1; addendum 1.3 | move S5's design towards September's one factor at a time; per rung, the hit count and its permutation p | 3 |
| addendum 1.3; verify_G03 | the verifier's rung order: covariates → T+NK → limma + duplicateCorrelation → projection T_NK bin → universe and hit rule | 3.1, which keeps it. The engine step is split (R4a, R4b), and the hit rule is a second count at every rung (departure 1) |
| verify_G03 (corrections 3) | not to be omitted: the cell set, library blocking and the arm-pure-library null, the depth arm, the universe filters, the 30- vs 20-cell gate, T vs T+NK | 3.1-3.2; the 30-cell gate bites only on AML0310 among S5's samples, so it acts at R1 |
| verify_G03 (best_option 2) | "If the adult contrast still beats the library-arm-pure null at the last rung…" | 3.2, 7: the library null decides jointly with the sample null at every rung (decision 9.3) |
| G04 5.7 | rungs "+children" and "+Chen2023": September specification, permutation procedure, PAUMTZ/PAXMIJ out; readout includes P370 hits and sign share | 3.3. R6 → E1 → E2 share September's frozen procedure (4.1). The `library_id` block in the one-dataset rungs R6 and E1 is the declared reading of "September design" (departure 2) |
| G04 5.7 and 7 | G03 states its reading rules before it runs, and opens the September per-gene tables only in the run | 7; 5.3 |
| G04 7 | relate the result to 1.1's prediction; take (b)'s A_b and φ and their CIs under any outcome | 5.4 |
| assessment 2.1 | read the 370 genes in the September table against matched random sets | 5.3: descriptive, without random sets (departure 3; decision 9.9) |
| addendum 1.3; B 3(c) | disclosures | 1.5 |
| A1-BATCH, "Moved out" | the hit-count ladder goes to G03 after G04's tag | this document |

Declared departures:
1. The hit rule is not a rung. Both rules are counted at every rung, and R6 is read on both (7). The verifier read
   "the hit count under that rung's own rule" at the last rung only.
2. R6 and E1 keep the `library_id` block. The September pre-registration's single-dataset form (§6.2) drops the
   block, arguing that libraries are singletons there. GSE185381's are not. 4.5 runs R6 and E1 without the block.
3. The 370-in-September read uses no matched random sets. P370 was selected on the shared samples, and unselected
   genes cannot calibrate that selection.

## 1. Already known before this draft (disclosed)

**1.1 The September T_NK result and the prediction on record.** G04 1.1 holds the full text.

| item | value |
|---|---|
| samples | 53: 39 AML (5 Chen2023, 20 GSE185381 children incl. PAUMTZ and PAXMIJ, 14 GSE185381 adults) vs 14 healthy (4 Chen2023, 10 GSE185381) |
| model | `~ dataset + arm`, two-pass voom with duplicateCorrelation on `library_id`, consensus 0.224 |
| genes tested | 9,281 |
| hit rule | q < 0.05 and \|logFC\| ≥ 1 |
| hits | 1 (AC007952.4), inside both Monte Carlo nulls: p 0.102 and 0.320, 1,000 draws each |
| planted recovery | 0.005 / 0.505 / 0.905 / 0.995 at \|logFC\| 0.5 / 1 / 1.5 / 2; LFC80 1.5 |

The prediction (addendum 1.3, 2026-10-28): about 156 of the 370 genes would have been called if all 39 AML carried
S5's sizes; about 4 (about 1 if shrunk by r) if only the 14 adults did. The fit called 1.

**1.2 S5 as run** (B section 1).

| item | value |
|---|---|
| samples | 15 adult AML vs 10 healthy |
| cell floor | ≥ 20 kept T cells |
| cell set | marker-concordant T (CD4 and CD8 fine labels), `high_error` dropped |
| model | edgeR QL robust; `~ status + hb + gran + sex + CD8 share`; no library term |
| gene universe | GSE185381 axis minus HB, granule, 79 ambient, XIST, chrY, MT and RP genes; then `filterByExpr(group)` and CPM > 1 in half the samples; 7,379 genes |
| hit rule | BH FDR < 0.10 and \|logFC\| > 0.5 |
| hits | 402 (171 up / 231 down) |
| permutation | 200 free relabellings, full refit; 4 of 200 reached 402, p = 5/201 = 0.025 (95 % interval of the exceedance rate 0.005-0.050) |

Other facts on record:
- The four draws that reached 402 placed 7, 9, 10 and 6 of the 15 true AML in the AML label; chance expects 9
  (`ws_a1_g02_s5_perm_draws.csv`, read in review 1). So the null's tail comes from another strong axis of the data,
  not from near-inverted labels.
- Of the 370 genes, 252 have |log2FC| ≤ 1 and 17 have |log2FC| ≥ 1.5 (verify_G03, data_claims_checked).
- S5's unblocked limma-voom sensitivity gave 678 hits, 383 of them in the 402.

**1.3 G04's results** (G04 section 11).
- (a) Chen2023: "not shown"; D 0.712 J, exact p 17/126; ρ 0.518 with one-sided 95.2 % bounds [−0.247, 1.361]. Without
  AML162, D is 0.225 (recorded in G04 11; its qualifier does not apply under "not shown").
- (b) children: "not readable"; A_b 0.372 [0.161, 0.594]; φ −0.244 [−0.827, 0.092] (95 % percentile intervals).

**1.4 Seen while writing this draft.**
- Read in full: the September roster (`pb_sample_manifest.csv`, sample-level) and `power.csv`.
- Every other file under `results/tables/12_pseudobulk_de/`: names, sizes, headers, row counts and md5 sums only.
- T_NK cell counts per sample, from the reconciled label columns (A 3.2; no expression).
- S5-side outputs (B), including by accident the S5 logFC of the 32 ambient-prone genes, and S5's stored permutation
  draws (review 1).
- No September gene name, logFC, p or q, except AC007952.4's, which are already on record (G04 1.1).

**1.5 Disclosures** (each goes into the paper's limitations in the form decided at writing).
1. The September pre-registration promised not to reframe its null as "underpowered" or "needs a different model"
   (:41-44). G03 is bound by that promise (section 0, section 7).
2. S5's hit rule (FDR 0.10, |logFC| 0.5) differs from September's (q 0.05, |logFC| ≥ 1). Handoff v2.0 (:137; §16.2
   item 5) proposed S5's rule before S5 ran, but gave no reason.
3. S5 has no library term. The September write-up says that without one "the healthy n is a fiction"
   (`FINDINGS_pseudobulk_de.md`:461-463). A1-BATCH G2 later tested lane clustering on the cross-fitted gap, not on the
   hit count.
4. Each S5 difference has a pre-S5 documented reason (verify_G03, best_option 1):
   - the 2026-10-16 Early Lymphoid / CCC audit led to marker-concordant T (`config_ws.R`:191-200);
   - no age column existed, so children could not be separated;
   - ambient RNA, first measured in S1, led to the ambient covariates (handoff §16.1 item 13).
5. S5 knew of the September line: the S5 pre-registration (:123) cites `12_pseudobulk_de/power.csv`.
6. T_NK entered the September analysis only through the post-hoc marker-gate amendment (September §8.2). The September
   write-up records its pre-registration as amended twice and corrected once, after seeing results.
7. The September fit pooled children with adults because no age column existed (September :422-424). Handoff v2.0
   forbids that pooling.
8. September p values are Monte Carlo over 1,000 draws; the September write-up calls them "exact".
9. Reading T_NK's single hit as zero is a FINDINGS reading (:141-142), not a pre-registered rule. The gene stays in
   `hitlist_frozen.csv`.
10. The September universe is a symbol intersection. AC007952.4 pools two Ensembl loci when Chen2023 columns are
    present.
11. S5 deviations recorded elsewhere: the `multiplexed` coding (assessment 6.2 item 4); GATE_A1(4) registered 3 vs 7
    but run 2 vs 8; the healthy age miscount; Control4 is female by expression and male in the metadata.
12. Healthy lanes P02 (Control1-4) and P03 (Control5) were captured together with undeposited T-ALL/ETP-ALL samples.
    They are flagged in every T contrast (4.4).
13. The blind-safe rungs were not run before G04, as the addendum had recommended. They are specified here with G04's
    result known.
14. 24 of S5's 25 samples are in the September fit, and both per-gene tables are public. "Unread" is attested, not
    enforced.

## 2. Data

**2.1 Samples.**

| set | AML | healthy |
|---|---|---|
| R0: S5's 25 | the 14 below plus AML0310 | the 10 below |
| R1-R6: the 24 shared adults | AML0048, AML0114, AML0612, AML0693, AML1371, AML2123, AML2451, AML3133, AML3266, AML3730, AML3762, AML3948, AML4340, AML4897 (14) | Control0004, Control0005, Control0058, Control0082, Control1, Control2, Control3, Control4, Control4003, Control5 (10) |

AML0310 is not in the September roster. It is CCC-ineligible (`below_min_occupied_bins`), and its T_NK bin holds 29
cells under the September filter, below the 30-cell gate (A 3.3).

Library structure of the 24 (17 libraries):
- P01 holds AML0612, AML3133 and AML3762 (AML0310 is also in P01, at R0).
- P02 holds Control1-4.
- P05 holds AML2451 and Control0004; it is the only mixed library.
- P15 holds AML1371 and AML4340.
- Every other sample is the only rung member of its library. Control0005 (P18), Control0058 (P06), Control0082 (P12)
  and Control5 (P03) share their lanes with samples outside the rung.

**2.2 Inputs** (read-only; md5 checked in K0).

| input | used in |
|---|---|
| S5 pseudobulks `/LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds` (`$raw$T`, `$raw$NK`), subset to the rung samples at load | R0-R4b |
| S5 design table (`27_s5_imprint.R`:55-67; rows in `ws_s5_scores.csv`) | R0, R1, 5.8 |
| September per-sample pseudobulks `/LARGE1/.../12_pseudobulk_de/<dataset>/<sample>__pseudobulk.rds`, T_NK column, opened by explicit name (no glob). These are the stored files, not re-aggregated: the Chen2023 QC objects were rewritten on 2026-10-05 (A 3.1). | R5, R6, E1, E2 |
| September depth tables `results/tables/12_pseudobulk_de/depth/<dataset>__<sample>__pb_depth.csv`, by explicit name | R5, R6, E1, E2 |
| `library_id` and capture period from `results/tables/ws/ws_a1_batch_covariates.csv`. It holds all 25 samples, AML0310 included, in September's `GSE185381__…` format | R0-R4b, 5.7 |
| `library_id` from `pb_sample_manifest.csv`. It agrees with the line above on the 24 (K2d) | R5, R6, E1, E2 |

**2.3 Extension samples** (G04 5.7).
- E1 adds the 18 children of G04(b): AML001, AML003, AML005, AML006, AML009, AML012, AML022, AML025, AML026, AML028,
  AML038, AML043, AML048, AML051, AML052, AML055, AML056, AML073.
- E2 adds AML103, AML162, AML24, AML304, AML320 and NBM1-4. These are September's pseudobulks: both pools, projection
  bin, symbol axis.
- PAUMTZ and PAXMIJ stay closed (G04 2.3): no column of theirs is read. ws_pb holds their columns and the September
  directory holds their files, so ws_pb is subset immediately and files are read by explicit name.

## 3. The ladder

### 3.1 Rungs (each changes one factor; order of the verifier, decision 9.1)

| rung | samples | cells | design | engine | gene universe |
|---|---|---|---|---|---|
| R0 | S5's 25 | S5: concordant T, ≥ 20 | `~ status + hb + gran + sex + CD8 share` | edgeR QL robust | S5 (blacklist; `filterByExpr(group)`; CPM > 1 in half) |
| **R1** | **24** (AML0310 out) | as R0 | as R0 | as R0 | as R0 |
| **R2** | 24 | as R0 | **`~ status`** | as R0 | as R0 |
| **R3** | 24 | **concordant T + NK** (`$raw$T + $raw$NK`) | `~ status` | as R0 | as R0 |
| **R4a** | 24 | as R3 | `~ arm` | **limma-voom** (`filterByExpr` as S5; TMM; voom; lmFit; eBayes) | as R0 |
| **R4b** | 24 | as R3 | `~ arm` + **`library_id` block** | **two-pass voom and duplicateCorrelation** (September 02:157-167) | as R0 |
| **R5** | 24 | **the September cell definition**: the reconciled `hierarchy_bin` T_NK (includes `marker_fill` cells), `in_ccc_graph`, not `high_error`, Early Lymphoid out, no concordance requirement, counts from `01_per_sample_qc`, ≥ 30 cells | as R4b | as R4b | as R0 |
| **R6** | 24 | as R5 | as R4b | as R4b | **September base**: the pinned symbol intersection Chen2023 ∩ GSE185381 minus 91 sex genes (21,843), then `filterByExpr(design)` only |

Notes:
- R1's floor, cell set and engine are S5's. The R3 floor (≥ 20 on T, as S5) is moot for the 24; all pass on either
  definition. The 30-cell floor at R5 is a no-op for the 24.
- R6 is "the September specification on the 24 S5 samples" (G04 5.7), with the September rule read as the second count
  (3.2).
- The September §5.2 bin rules fail by construction for a one-dataset fit and are not applied. The §8.2 marker gate is
  reported, not applied (5.6).

### 3.2 Statistics at every rung

**Counts.** Both come from one fit.
- **c_S5**: genes with BH q < 0.10 and |logFC| > 0.5.
- **c_Sep**: genes with q < 0.05 and |logFC| ≥ 1 (the September baseline rule and permutation statistic).
- At R6, E1 and E2 each count is also given restricted to genes outside S5's blacklist (B 1, 27:48-49). These are
  **c_S5′** and **c_Sep′**. The q values come from the rung's fit over its whole tested universe; blacklisted genes are
  then dropped from the count, not from BH. Every permutation draw is treated the same way, at no extra compute. R6
  re-admits the HB, granule, ambient, MT and RP genes, and S5's ambient gate failed in T. So the category at R6 uses
  c_S5′ (section 7), and c_S5 is reported beside it.

**Nulls.** Each count gets two Monte Carlo p values, p = (1 + #{null ≥ observed}) / (1 + valid draws):
- **p_free**: free relabelling of samples, arm sizes kept (S5's null; the September primary null, within dataset).
- **p_lib**: whole arm-pure libraries relabelled, mixed libraries held fixed (September 03:116-127). Library units are
  permuted, so arm sizes in samples may vary.
- Arrangements: 3,268,760 free and 8,008 library for the 25; 1,961,256 and 8,008 for the 24 (K5).

**Permutation procedure.**
- **Full refit at every ladder rung, R0-R6.** Each draw reruns the whole rung pipeline: filtering, TMM, dispersion or
  voom, both duplicateCorrelation passes where present, and the fit (decision 9.4).
- At R4b and R6, September's frozen scheme (03:149-160) is a sensitivity (4.1). The only comparison on record shows the
  frozen scheme can give a smaller p than a full refit (LMPP_GMP 1/126 vs 2/126; September :700-702), so neither
  scheme is assumed conservative.
- E1 and E2 use September's frozen scheme, as G04 5.7 binds (3.3).

**Draws.**
- Draws are stored as sample-keyed labellings and joined by ID at every rung (K7).
  - R0: 1,000 free (`SEED + 601`) and 1,000 library (`SEED + 602`) labellings of the 25.
  - R1-R6: 1,000 free (`SEED + 603`) and 1,000 library (`SEED + 604`) labellings of the 24, the same draws at every
    rung.
- Common draws reduce Monte Carlo noise between rungs; they do not remove it. R0 and R1 cannot share draws, because
  their arm sizes differ. The R0 → R1 step therefore compares two independent Monte Carlo estimates.
- A draw whose fit fails (for example rank deficiency) is caught, counted as NA, excluded from "valid draws", and the
  number of such draws is reported.

**The deciding count.** It is c_S5 at R0-R5 and c_S5′ at R6.

**State of each null for the deciding count.**
- **p_free** (Monte Carlo). Take the exceedance count x of n valid draws, with its Clopper-Pearson 95 % interval
  [lo, hi] for x/n:
  - **held**: hi ≤ 0.05;
  - **lost**: lo > 0.05;
  - **borderline**: otherwise.
  - At n = 1,000 this means held ⇔ x ≤ 36 and lost ⇔ x ≥ 65.
- **p_lib**. Read the same way at 1,000 Monte Carlo draws. Its space is small (8,008 labellings for the 25 and for the
  24), so a borderline p_lib is resolved by exact enumeration, not by more random draws. Fits already computed for
  labellings drawn in the first 1,000 are reused. The exact p is held if ≤ 0.05 and lost otherwise. This follows
  September's own 2026-09-24 correction for small spaces (September :693-702).

**State of a rung.** The worse of its two nulls' states, at every rung from R0 to R6 (lost > borderline > held;
decision 9.3). The deciding null therefore does not change along the ladder.

**Extension.** It is fixed now and applies to the deciding count only.
- A borderline p_free is extended to 10,000 draws: the first 1,000 plus 9,000 new draws (`SEED + 605` for the 25,
  `SEED + 606` for the 24, common across R1-R6). The state is then re-read and is final.
- A borderline p_lib is enumerated (above).
- A null is not extended if the rung's state is already fixed by its other null, i.e. that null is lost.
- Priority order:
  1. rungs up to and including the first lost rung, in rung order;
  2. then the later rungs, in order;
  3. nothing else. The R6 second axis (c_Sep′), the non-deciding counts and the extension rungs are never extended.
- Cap: 300 CPU-h in total. Before each extension starts, its cost is projected from K8's timing. It runs whole if the
  projection fits the remaining budget, and is skipped otherwise; no partial extension is used. A skipped extension
  leaves the rung at its 1,000-draw state, and the wording says so (section 7).

**Also per rung (descriptive).**
- genes tested and the duplicateCorrelation consensus;
- c_Sep after the depth arm at R5, R6, E1 and E2: a hit both in the baseline fit and with `log10(median_umi)` added
  (September 02:190-197). No depth arm is run at R0-R4b;
- the draws that exceed the observed count: their overlap with the true labels, their library and capture-period
  alignment (5.7);
- 5.1's P370 readouts.

### 3.3 Extension rungs (fixed by G04 5.7; nothing else changes)

| rung | samples | design | libraries |
|---|---|---|---|
| E1 "+children" | R6's 24 + the 18 children = 42 (32 AML v 10 healthy) | `~ arm` + `library_id` block | 22 (19 pure, 3 mixed: P05, P12, P18) |
| E2 "+Chen2023" | E1 + 9 Chen2023 = 51 (37 v 14) | `~ dataset + arm` + block (exactly September 02:147-163) | 31 (28 pure, 3 mixed) |

- Universe: the pinned 21,843, then `filterByExpr(design)`.
- Permutations: September's procedure. That is 1,000 draws per null with the frozen scheme. Seeds use September's
  `.perm_seed()` with labels `perm|sample_within_dataset|g03_<rung>|T_NK` and `perm|library|g03_<rung>|T_NK` (03:74-78,
  180).
- Readout: c_Sep and c_Sep′ with p_free and p_lib; c_Sep after the depth arm; c_S5 and c_S5′ with both p values from the
  same draws; 5.1's P370 readouts.
- The September fit (53, 1 hit) is the reference. E2 differs from it only by PAUMTZ and PAXMIJ, which also changes
  P06's library purity.
- So that "nothing else changes" between R6 and the extension rungs, the extension sentence (section 7) takes R6's p
  from 4.1, which uses September's frozen scheme on R6's draws. R6 → E1 → E2 → September then share one permutation
  procedure. R6's full-refit p stays in the rung table.

## 4. Sensitivities (each reported; none changes a category)

- **4.1 Frozen scheme at R4b and R6.** September's frozen voom weights and consensus, using the same draws as 3.2. At
  R6, it also supplies the p of the extension sentence (3.3), for c_Sep and c_Sep′ under both nulls.
- **4.2 Control4 coded female at R0 and R1.** Only those rungs have a sex covariate.
- **4.3 E1 without AML073.** AML073 has 210 of 605 kept T cells labelled malignant by the authors; the September cell
  set kept them. c_Sep and c_S5, no p.
- **4.4 R6 and E1 without Control1-5.** These are the T-ALL lanes P02 and P03. c_S5′ and c_Sep′, no p.
- **4.5 R6 and E1 without the `library_id` block.** This is the September single-dataset form, under the frozen
  scheme, with the same draws.

Path dependence is a limitation, not a sensitivity. A factor's effect can depend on the order in which factors change,
and the reverse ladder is not feasible for every factor: a T-only projection-bin pseudobulk does not exist without new
per-cell aggregation (B 5.2).

## 5. Readouts (descriptive; no category)

**5.1 P370 per rung.** Among the P370 genes in the rung's tested universe:
- the number tested;
- hits under each rule;
- the share with S5's sign;
- the median of `sign_S5 × logFC_rung / |logFC_S5|` (attenuation).

At R0-R6 this is circular, because P370 was selected on the 25, so it is flagged. E1 and E2 add samples P370 never
saw.

**5.2 AC007952.4 per rung.** Its logFC, q and universe membership. It is one locus (ENSG00000262202) at R0-E1, and two
pooled loci at E2.

**5.3 The stored September table** (`de/primary__T_NK.csv`, opened here for the first time; md5 checked first). Among
the P370 genes in its 9,281, report:
- the number of genes and of hits;
- the share with S5's sign;
- the median `sign_S5 × logFC`;
- the median signed z, `sign_S5 × sign(logFC) × qnorm(P.Value / 2, lower.tail = FALSE)`.

There is no test. It is labelled circular, because 24 of the 53 samples are shared.

**5.4 The prediction recomputed with G04's estimates.** It uses only numbers on record, no September per-gene value.
- Method, as in the addendum: `power.csv`'s full-rule recovery curve, interpolated linearly, applied to each P370
  gene's |f × log2FC_S5|.
- ρ and φ are shares of the selection-corrected gap and of the out-of-fold gap, so both are relative to about S5's size
  divided by r. f is therefore defined by two rows:
  - **unshrunk** (adults carry S5's sizes): f = (14 + (5ρ + 20φ) / r) / 39;
  - **shrunk** (adults carry S5 / r): f = (14 + 5ρ + 20φ) / (39 r).
- r = 1.43. φ is applied to all 20 children, PAUMTZ and PAXMIJ included. |f| is used when f < 0; the curve is
  symmetric in sign.
- Rows are computed at the point estimates and at the corners of ρ's one-sided bounds and φ's interval. Each row gives
  the expected number of September calls among the 370, against the observed 1.
- The addendum's caveat applies: these are rough upper estimates, because the curve comes from mid-tertile genes
  planted into every AML sample.
- A_b, φ and their CIs go into G03's GATES block as G04 handed them over.

**5.5 Planted recovery at R6, E1 and E2.** September's grid: 200 middle-tertile genes, with logFC 0.5, 1, 1.5, 2 and 3
added to the voom E of every AML sample, frozen weights (03:202-224; seed `.perm_seed("planted|g03_<rung>|T_NK")`).
Reported: recovery under the September rule, and LFC80. LFC80 is the licence of that rung's null: "no effect of
|logFC| ≥ LFC80 in every AML sample".

**5.6 Marker gate at R5, R6, E1 and E2.** September §8.2's amended gate (CD3E/CD3D: mean CPM ≥ 100 in both arms,
|log2 arm ratio| ≤ 1). Reported only.

**5.7 Exceeding draws per rung.** For every draw at or above the observed c_S5 (c_S5′ at R6):
- the overlap of its AML label with the true AML;
- its alignment with P01, P02 and P15;
- its alignment with the capture period (`ws_a1_batch_covariates.csv`).

**5.8 References for R1 and R2** (c_S5, no p).
- R1: leave-one-AML-out on R0's design, for each of the other 14 AML, so that "depends on AML0310" has a reference.
- R2: R1 with each covariate dropped singly (hb, gran, sex, CD8 share), so that a loss at R2 can be traced.

## 6. Known answers (run first; three states as in G04)

| ID | scope | check |
|---|---|---|
| K0 | all | md5 of the declared inputs: A section 5's September files, the 51 rung pseudobulk RDS files, `ws_pb_GSE185381.rds` (6c2fb685…), S5's gene, score and draw tables, `ws_a1_batch_covariates.csv` |
| K1 | L | R0's observed fit reproduces S5: 7,379 tested, 402 hits (171 up / 231 down), P370 membership |
| K1b | L | the copied edgeR pipeline, run on S5's 200 stored labellings, reproduces every draw's n_hits exactly, including 1,615 / 1,166 / 707 / 563. `ws_a1_g02_s5_perm_draws.csv` stores each draw as a 25-character string. Its positions follow the T rows of `ws_s5_scores.csv` in file order: 27's design order, AML first, then healthy (35:48-49, 76-80) |
| K2 | L, E | for every rung sample, the stored T_NK column sum equals the depth table's `total_counts`, and `n_cells` ≥ 30 equals `passes_30`. The base universe rebuilt from gene names is 21,843 |
| K2b | L | `identical()` of the gene axes of `ws_pb` and of every GSE185381 September RDS used. For two samples, R3's cell set is re-aggregated from `01_per_sample_qc` into memory and compared with `identical()` against `$raw$T + $raw$NK`, without printing. This shows the count source is not a hidden factor at R5 |
| K2c | L | per sample, R3 vs R5 cell counts and barcode overlap (metadata only; reported) |
| K2d | L | on the 24, `library_id` from `ws_a1_batch_covariates.csv` equals `library_id` from `pb_sample_manifest.csv` |
| K3 | L, E | the copied September engine and permutation code on September fits that hold no closed sample. Main target: the Chen2023 secondary bin LMPP_GMP (5 v 4, 304 baseline hits), reproducing `de/secondary__LMPP_GMP.csv` logFC (max \|Δ\| < 1e-8), the hit count, and `perm_null.csv`'s per-draw n_hits for its sample null from `.perm_seed("perm|sample_within_dataset|secondary|LMPP_GMP")`. Second target: HSC_MPP (4 v 4) |
| K3c | L | the full-refit limma loop (rerunning filterByExpr, TMM, voom, lmFit and eBayes per labelling, `~ arm`, September filter) over LMPP_GMP's 126 labellings reproduces the `full_refit` row of `secondary_exhaustive_perm.csv`: `observed_hits`, `n_ge_observed` and `null_max` (07:84-91, 144-159) |
| K4 | L, E | the blocked two-pass fit equals September 02's `fit_bin`, copied verbatim and run side by side, on synthetic counts with shared libraries |
| K5 | L, E | permutable units: 3,268,760 / 8,008 (25); 1,961,256 / 8,008 (24); A 3.4's counts for E1 and E2 |
| K6 | L | each rung's per-gene logFC vector differs from the previous rung's (September §8.6's "differs from baseline" tell) |
| K7 | L | stored draws are keyed by sample. At every rung the AML label set of each draw is the same set of IDs |
| K8 | all | a smoke run of every mode on synthetic counts, with no real covariate written next to a score. It also times one full refit per engine for the compute note (8) |

A deterministic failure is fixed, a dated note is added, and everything is rerun. A failed scope-L known answer makes
the ladder "not decisive". A failed scope-E known answer makes the extension readout "not run".

## 7. Reading rules (fixed now; first match wins)

**Ladder.**
1. **"not decisive"**: a scope-L known answer failed.
2. **"S5's count is not beyond its own null"**: R0 is lost.
3. **"lost at <rungs>"**: at least one rung in R1-R6 is lost.
   - Every rung Rk that is lost while R(k−1) is held is named as a loss of its class. A validity loss so named carries
     its consequence, whatever happens at later rungs.
   - A rung lost after a lost or borderline rung is reported as "lost, not attributable to <factor>", because the ladder
     is cumulative.
   - Borderline rungs are named. A rung held again after a loss is reported as "held again at Rj".
4. **"borderline"**: no rung is lost, and R0 or at least one rung is borderline. The borderline rungs are named.
5. **"held through R6"**: R0-R6 are all held.

**Second axis, R6 under the September rule.** The state of c_Sep′ (worse of p_free and p_lib) at R6: beyond its null,
borderline, or inside its null. Every ladder category is reported with it.

**Wording.**
- **5** → "On the 24 adults that S5 and the September fit share, S5's rule keeps c_S5′ hits beyond both permutation
  nulls (p_free, p_lib) through the September cell definition, engine with library block, and gene universe. Under
  the September rule the count is c_Sep′ (p_free, p_lib; state). The two analyses count different effects: S5 counts
  |logFC| > 0.5 at FDR 0.10, September |logFC| ≥ 1 at q 0.05. The September T_NK result stays reported as a null for
  |logFC| ≥ 1.5 per gene."
- **3**, k ∈ {R2, R3, R4a, R5, R6} (specification rungs) → "S5's per-gene count needs <factor; for R5, the S5 cell
  definition: marker concordance, the fine-label system and the `in_ccc_graph` filter>. Without it the count falls
  inside its permutation null (rung table). The per-gene list is reported with this qualifier."
- **3**, named loss at R1 → "S5's per-gene count does not survive removing AML0310 (5.8: removing any other single AML
  gives c_S5 from a to b)."
- **3**, named loss at R4b → "S5's per-gene count does not survive the library-level correlation of a blocked model."
- **3**, any named validity loss (R1 or R4b) → "The per-gene list is reported as exploratory. The sample-level
  separation (A1-G02, A1-BATCH G2) is a different statistic and stands as reported there."
- **3**, every sentence names the null or nulls under which the rung was lost.
- **4** → "S5's per-gene count is borderline at <rungs>; no rung loses it." For each rung, the wording adds either
  "(Monte Carlo interval straddles 0.05 after 10,000 draws)" or "(1,000 draws; extension not funded under the cap)".
  The per-gene list carries this qualifier.
- **2** → "S5's own count is not beyond its permutation null (p_free, p_lib, number of draws). The per-gene list is
  reported as exploratory."
- The rung classes and their consequences are decision 9.7.

**Extension (E1, E2).** No category. Fixed sentence, on September's count c_Sep, with c_Sep′ in parentheses, all p
from the frozen scheme (3.3):

"On the 24 adults the September specification finds c_Sep (c_Sep′) hits (p_free, p_lib). Adding the 18 children gives
c_Sep (c_Sep′) (p_free, p_lib); adding Chen2023 gives c_Sep (c_Sep′) (p_free, p_lib). The September fit on 53 found 1,
inside its nulls. G03 does not test which samples account for any difference."

The P370 hits at R6, E1 and E2 follow.

**Prediction (5.4).** Fixed sentence: "With G04's estimates, the addendum's method expects X (unshrunk) and Y (shrunk)
September calls among the 370, ranging from a to b over the interval corners; the fit called 1. These are rough upper
estimates from a planted-curve approximation. They are not a test, and they are not read as explaining the September
result."

**Never written, under any outcome:**
- "September was underpowered", "used the wrong model", or "is explained by, or traces to, pooling the children or
  Chen2023";
- "the imprint is confirmed by the September data";
- a reading of 5.1 at R0-R6, or of 5.3, as evidence (both circular);
- a reading of 5.4 as explaining the September null;
- a factor named from a borderline rung, or from a loss that follows a lost or borderline rung, as if it were an
  attributable loss.

**Handover.**
- `--combine` appends a GATES.md block with the category, the second axis, the rung table and the fixed sentences.
- A dated note goes into handoff S15.
- G03 has no pass/fail gate of its own. It sets the status of the per-gene list (section 0) for G07 and for writing.

## 8. Code, outputs, compute (plan; the script is written after the tag)

**Script.**
- One script, `scripts/ws/42_a1_g03.R`, with modes `--known`, `--ladder`, `--extend`, `--sens`, `--read` (5.3, 5.4)
  and `--combine`. Each mode can also run with `--smoke`.
- A wrapper, `submit_ws_a1_g03.sbatch`, runs dated read-only copies of the script, as in G04. The ladder's permutation
  draws are split over a SLURM array.
- Functions are copied, never sourced. September 01/02/03 and ws/26, 27 and 37 write frozen public tables and are not
  run (B 6.1).

**Outputs.**
- `results/tables/ws/ws_a1_g03_{checks,ladder,extend,sens,p370,read,summary}.csv`.
- Draws and null counts under `/LARGE1/.../ws_a1_g03/`.
- Every target is checked at mode start, and the script refuses to overwrite. A later mode stops unless the checks
  table carries the same script md5.

**Compute** (estimate; the K8 timings go into a dated pre-run note before any real mode runs).
- edgeR rungs R0-R3: about 8,000 full refits, about 9-12 CPU-h.
- R4a: about 2,000 voom refits, under 2 CPU-h.
- R4b, R5, R6: about 6,000 full refits with two duplicateCorrelation passes, about 50-100 CPU-h.
- Frozen-scheme rungs and sensitivities: about 1-3 CPU-h.
- Extensions are capped at 300 CPU-h (3.2).
  - An edgeR extension (R0-R3) costs about 10-14 CPU-h.
  - A blocked extension (R4b-R6) costs about 25-60 CPU-h at September's implied rate (03:44-48, about 22 s per
    blocked refit at n = 53) and up to 75-170 CPU-h at this plan's conservative rate.
  - p_lib enumeration costs about 7,000 refits per rung.
  - The priority order in 3.2 decides which extensions are funded.
- Every job sets `OMP_PROC_BIND=false`.

**Seeds.** `SEED = 491638`. Ladder offsets are +601 to +607; `scripts/ws` uses offsets up to +534 (39's 4.7 uses
+35 to +534). E1, E2 and the planted grid use September's `.perm_seed()` labels.

## 9. Decisions (the user, 2026-11-05: all nine as recommended)

Real choices are marked ★.

| # | decision | recommended | alternative |
|---|---|---|---|
| 9.1 | rung order | the verifier's order, engine split into R4a/R4b, AML0310 out first (R1) | AML0310 leaves at R5, where the September cell set first excludes it |
| 9.2 | the hit rule is a second count at every rung, read at R6 as the second axis | yes | a final rung R7 with the September rule only |
| 9.3 ★ | which null decides the state | the worse of p_free and p_lib at every rung, R0-R6. This is the verifier's library-null rule, applied throughout so that no rung changes the deciding null | p_free throughout, with the qualifier "held under the sample null only (library p = …)" when p_lib is not held |
| 9.4 ★ | permutation procedure on the ladder, and compute | full refit at every rung (about 60-115 CPU-h before extensions; extensions capped at 300 CPU-h); frozen as a sensitivity at R4b and R6, and for the extension sentence | September's frozen scheme for R4b-R6 (about 3 CPU-h; the procedure then changes together with the engine) |
| 9.5 | E1 keeps the `library_id` block | yes; 4.5 runs it without the block | the September single-dataset form, no block |
| 9.6 | R6 universe | September's pinned symbol base (21,843) | Ensembl matching (not September's specification) |
| 9.7 ★ | rung classes | validity = R1 and R4b; a validity loss makes the per-gene list "exploratory". Specification = R2, R3, R4a, R5 and R6; a loss is named as a qualifier | no classes: a loss only names the factor, with no consequence fixed now |
| 9.8 | 4.3-4.5 and 5.8 | include (descriptive, cheap) | drop |
| 9.9 | the stored September table (5.3) | descriptive only, no random sets | matched random sets, as in the assessment (circular either way) |

## 10. Review history

- **DRAFT v1**, 2026-11-05: written from inventories A and B.
- **Review 1**, 2026-11-05 (one read-only reviewer: 4 blocking, 10 should-fix, 12 minor; notes in
  `/FAST/gr10634/gaozy/reports/2026-11_a1_g03_prereg/review1/review.md`). DRAFT v2 folds in all of the blocking and
  should-fix findings and the minor ones:
  1. R0 is a rung with its own draws. A three-state rule (Clopper-Pearson) with one extension to 10,000 draws
     replaces the bare 0.05 crossing.
  2. The first-match order is rewritten, so that a later rung cannot override a loss. R6 is read on blacklist-
     restricted counts.
  3. R5 of v1 is split into R4a (engine) and R4b (library block). Full refits are used on the whole ladder; the frozen
     direction claim is withdrawn.
  4. The sentence attributing the September null to the added samples is replaced by a descriptive one and added to
     "never written".
  5. The library null decides jointly from R4b on (decision 9.3).
  6. The verifier's rung order is restored.
  7. R6 has a two-axis reading.
  8. R5 is named as the S5 cell definition, with gene-axis and count-source known answers (K2b, K2c).
  9. The edgeR permutation path is checked on S5's stored draws (K1b).
  10. K3 targets LMPP_GMP.
  11. 5.3 uses logFC and a signed z.
  12. Disclosures 4, 5, 11 and 12 are added, and the AML162 wording is corrected.
  13. Minor findings: seeds moved to +601 and up; library wording; the AC007952.4 locus; 5.4's scale and |f|;
      sample-keyed draws (K7); explicit-name loading; failed-draw handling; the compute figures; the citation; K7 of
      v1 replaced; 5.8's references.
- **Review 2**, 2026-11-05 (one read-only reviewer: 2 blocking, 7 should-fix, 2 minor; notes in `…/review2/review.md`).
  All 14 of review 1's blocking and should-fix findings were confirmed resolved, and the counts, arrangement numbers and
  seed offsets were checked. DRAFT v3 folds in all of review 2:
  1. The extension rule is fixed. Only the deciding count is extended. p_free: the first 1,000 draws plus 9,000 new.
     p_lib: exact enumeration. A null is not extended once the other null fixes the state. Priority order as in 3.2.
     The cap is tested before each extension on a projected cost, all or nothing.
  2. R4b no longer changes the deciding null: both nulls decide at every rung (9.3).
  3. Every loss that follows a held rung is named and classified. A validity loss keeps its consequence. A loss after a
     lost or borderline rung is "not attributable". R1's sentence cites 5.8.
  4. The extension sentence is written on c_Sep, with R6's p from the frozen scheme, so R6 → E1 → E2 share one
     procedure. Departure 2 and 4.5 now cover R6.
  5. 5.4 has a fixed sentence, and "never written" excludes reading it as an explanation.
  6. p_lib uses exact enumeration instead of more Monte Carlo draws.
  7. R0's `library_id` comes from `ws_a1_batch_covariates.csv`, with K2d checking agreement on the 24.
  8. The BH set of the restricted counts is fixed (whole universe, then drop).
  9. K3c checks the full-refit limma path on LMPP_GMP's 126 labellings.
  10. K1b's positional mapping is stated, and R0 and R1 are declared not to share draws.
- **v1 APPROVED**, 2026-11-05: the user decided all nine section-9 items as recommended ("按推荐"). Tag
  `prereg-a1-g03-v1`. The script is written after the tag (section 8).
- **Code review 1 and dated pre-run note**, 2026-11-05. Written before any real-data mode of 42 was run. Before this
  note, 42 had touched no real count, apart from 5.4, which has no synthetic input (item 16).
  `scripts/ws/42_a1_g03.R` and `submit_ws_a1_g03.sbatch` were written in three stages, each run end to end on
  synthetic data. Two smoke patterns were run: signal in every rung, and signal on S5's side only. The second exercises
  the "lost" and "not attributable" branches.

  Three read-only reviewers checked the code (notes in `reports/2026-11_a1_g03_prereg/code_review1/`, outside the repo):
  - conformance with this document;
  - correctness on real data;
  - blindness and output safety.

  Two of them found the same blocking bug. 4.2 recoded Control4 as "F", while the S5 table codes sex as
  `female`/`male`, so Control4 would have got a sex level of its own. It is fixed and asserted, and the smoke covariates
  now use the real level names. All should-fix findings were fixed.

  Declared now; none of these can change a category:
  1. E1 and E2 are built in manifest order (September's roster order), so `.perm_seed` gives September's procedure.
     Their library IDs come from the manifest.
  2. K3 and K3c run the production functions (`fit_rung`, `refit_frozen`, `frozen_null`), not inline copies. K3 keeps
     only the secondary-bin rows of `perm_null.csv` and `de_summary.csv`.
  3. K2c is reported, as declared. September's barcode join needs to reach only 0.999, so its cell counts are not a
     deterministic target.
  4. K7b (added, scope L): every draw carries one labelling across all rungs that share its draw set.
  5. K9 (added, scope "R"): 5.4's method reproduces the addendum's rows 155.9 / 63.3 / 16.8 / 3.8 / 1.2. The recovery
     curve is linear, anchored at (0, 0) and flat above 3. A K9 failure suppresses only 5.4's sentence.
  6. Every mode after `--known` re-checks the input md5s that `--known` recorded. These now include `config_paths.R`,
     the two K2b per-cell objects, the 24 K2c label files and the two G04 tables read by 5.4. Chunk files and the
     observed-fit file carry the script md5. Every parallel worker must return; a dead worker stops the mode.
  7. `--extend` needs scope E only. It refits R5 and R6 itself, and runs R6's frozen nulls on the stored 24-sample
     draws (the same as 4.1 at R6). The extension sentence takes R6's p from there.
  8. `--sens` needs scopes L and E.
  9. `--read` opens `de/primary__T_NK.csv` only after the ladder, extension and sensitivity tables exist under the
     same script md5.
  10. `--combine` writes nothing if a table that its gate outcome requires is missing.
  11. The rung state "not computed" counts as the worst state. "Held through R6" requires every rung to be held.
  12. One output is added: `ws_a1_g03_checks_obs.csv` (K6). `ws_a1_g03_p370.csv` also carries genes tested and the
      duplicateCorrelation consensus per rung.
  13. 4.2's R0f and R1f get full permutation nulls on R0's and R1's draws (about 4,000 more edgeR refits, about 3.5
      CPU-h). They are never extended.
  14. 5.7 is read on the round-1 draws (the first 1,000 of each null). Failed fits are counted (`n_failed`) and
      excluded from the valid draws.
  15. Synthetic-only seeds are SEED + 697 to 699, above the declared +601 to +607.
  16. 5.4 was computed in the smoke runs from the real on-record numbers (it has no synthetic input). Point estimate:
      2.4 expected September calls unshrunk and 0.7 shrunk; 0.1 to 19.0 over the interval corners; against the
      observed 1. No category depends on it.
  17. Timings on synthetic data, per refit: edgeR 3.0 s, voom 0.85 s, blocked (two passes) 31.6 s, frozen 0.77 s.
      Round 1 is 172 array tasks, about 60 CPU-h, plus about 3.5 CPU-h for item 13. The extension cap is projected from
      the real `--known` timings, including a blocked refit on R6's universe.
  18. Each array task copies the script when it starts. The script is not edited while an array is queued; `gate()`
      stops any task whose copy differs.
