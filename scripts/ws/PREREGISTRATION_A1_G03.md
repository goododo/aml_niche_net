# PREREGISTRATION A1-G03: S5's T imprint and the September T_NK null, one factor at a time

**Status: DRAFT v1, 2026-11-05. Not approved; nothing declared below has been computed.**

- Written from two read-only inventories, kept outside the repo:
  `/FAST/gr10634/gaozy/reports/2026-11_a1_g03_prereg/inventory/A_september_spec.md` (A) and
  `B_s5_bridge_obligations.md` (B).
- No September per-gene value was read. Section 1.4 lists what was seen.
- Order of the WP-A1 steps (user decisions 2026-10-28 to 2026-11-05): A1-G02 (PASS) → A1-BATCH (PASS with
  qualifications) → A1-G01 → A1-G04 ((a) "not shown", (b) "not readable") → **this**.
- Sources: handoff v2.0 §17 S15; the G04 handover (`PREREGISTRATION_A1_G04.md` sections 5.7 and 7).

## 0. Question, scope, and what each outcome lets the paper say

Two analyses of largely the same samples disagree:

- **September, 2026-09-16** (`scripts/12_pseudobulk_de/PREREGISTRATION_pseudobulk_de.md`, commit a66cbd0). A pooled
  pseudobulk DE on the T_NK bin, 39 AML vs 14 healthy. It found no reportable hit beyond its permutation nulls.
- **S5.** On 25 GSE185381 adults, it found 402 genes; the frozen primary list P370 has 156 up and 214 down.
- **Overlap.** 24 of S5's 25 samples are in the September fit. A reviewer will read the pair as a forking path.

G03 asks one question: **which design factors separate the two results?** It answers it in three ways:

- **The ladder (section 3, rungs R1-R6).** It starts from S5's design on the 24 shared adults and changes one factor
  per rung until the design is September's.
- **The extension rungs (E1, E2), fixed by G04 5.7.** They add back the 18 children, then the 9 Chen2023 donors.
- **Relation to the prior prediction.** The result is compared with the prediction on record (1.1), using G04's
  estimates for the children and Chen2023.

What G03 is not:
- **Not a new test of whether the imprint exists.** A1-G02 tested that out of sample (gap 1.378, D2 p 0.008), and
  A1-BATCH G2 tested it against lane structure (p 0.017).
- **Not a reinterpretation of the September result.** That result stays reported as a null for |logFC| ≥ 1.5 per
  gene. No G03 result is described as showing that the September analysis was underpowered or used the wrong model.
  This repeats the September pre-registration's promise (`PREREGISTRATION_pseudobulk_de.md`:41-44) and G04's binding
  wording (G04 1.1).

What each ladder outcome lets the paper say. The wording is fixed in section 7; which factor classes count as
"validity" is decision 9.7.
- **"held to the September design".** On the 24 adults S5 shares with the September fit, S5's contrast keeps a hit
  count beyond its permutation null when the September cell set, engine and gene universe replace S5's. The paper shows
  the ladder table and the September null side by side.
- **"lost at a specification rung" (R2, R3, R4, R6).** The count needs one of S5's specification choices (covariates,
  T only, marker-concordant cells, S5's gene universe). The paper names that choice. The per-gene list keeps its status
  and carries the qualifier.
- **"lost at a validity rung" (R1 or R5).** The count depends on one sample (AML0310), or does not survive the
  library-level correlation of September's engine. The paper calls the per-gene list exploratory. The sample-level
  result (A1-G02, A1-BATCH G2) is a different statistic and is not changed by G03.

### 0.1 What G03 owes, and where it is handled

| source | obligation | here |
|---|---|---|
| assessment 2.1; addendum 1.3 | move S5's design towards September's one factor at a time; per rung, the hit count under the rung's rule and its permutation p | 3.1-3.3 |
| addendum 1.3 | the blind-safe rungs use S5's own samples (24 of 25 in the September fit) | 2.1, 3.1 |
| addendum 1.3 (verifier) | not to be omitted: cell set, library blocking and the arm-pure-library null, the depth arm, the gene-universe filters, the 30- vs 20-cell gate, T vs T+NK as its own rung | 3.1 (R3-R6), 3.2, 5.2 |
| G04 5.7 | rung "+children" and rung "+Chen2023", September specification, PAUMTZ/PAXMIJ out, readout incl. P370 hits and sign share | 3.3 |
| G04 5.7, 7 | G03 states its reading rules before it runs; it opens the September per-gene tables only there | 7; 5.3 |
| G04 7 | relate the result to 1.1's prediction; take (b)'s A_b, φ and CIs under any outcome | 5.4 |
| assessment 2.1 | read the 370 genes in the September table (against matched random sets) | 5.3 (descriptive; decision 9.9) |
| addendum 1.3 | disclose the September promise and the hit-rule difference | 1.5 |
| A1-BATCH, "Moved out" | the September/S5 hit-count ladder goes to G03 after G04's tag | this document |

## 1. Already known before this draft (disclosed)

**1.1 The September T_NK result and the prediction on record.** G04 1.1 holds the full text; it is not repeated here.

The fit:

| item | value |
|---|---|
| samples | 53: 39 AML (5 Chen2023, 20 GSE185381 children incl. PAUMTZ and PAXMIJ, 14 GSE185381 adults) vs 14 healthy (4 Chen2023, 10 GSE185381) |
| model | `~ dataset + arm`, two-pass voom with duplicateCorrelation on `library_id` (consensus 0.224) |
| genes tested | 9,281 |
| hit rule | q < 0.05 and \|logFC\| ≥ 1 |
| hits | 1 (AC007952.4), inside both Monte Carlo nulls (p 0.102 / 0.320, 1,000 draws each) |
| planted recovery | 0.005 / 0.505 / 0.905 / 0.995 at \|logFC\| 0.5 / 1 / 1.5 / 2; LFC80 1.5 |

The prediction (addendum 1.3, 2026-10-28): about 156 of the 370 genes would have been called if all 39 AML carried
S5's sizes, and about 4 (about 1 shrunk by r) if only the 14 adults did. The fit called 1.

**1.2 S5 as run** (B section 1).

| item | value |
|---|---|
| samples | 15 adult AML vs 10 healthy |
| cell floor | ≥ 20 kept T cells |
| cell set | marker-concordant T cells (CD4 and CD8 fine labels), `high_error` dropped |
| model | edgeR QL, robust; `~ status + hb + gran + sex + CD8 share`; no library term |
| gene universe | GSE185381 axis minus HB, granule, 79 ambient, XIST, chrY, MT and RP genes; then `filterByExpr(group)` and CPM > 1 in half the samples (7,379 genes) |
| hit rule | BH FDR < 0.10 and \|logFC\| > 0.5 |
| hits | 402 (171 up / 231 down) |
| permutation | 200 free relabellings with full refit; Monte Carlo p 0.025 |

Two further facts:
- 252 of the 370 genes have |log2FC| ≤ 1, and 17 have |log2FC| ≥ 1.5 (addendum 1.3).
- S5's limma-voom sensitivity, without blocking, found 678 hits, 383 of them in the 402.

**1.3 G04's results** (G04 section 11, 2026-11-05).

| arm | category | numbers |
|---|---|---|
| (a) Chen2023 | "not shown" | D 0.712 J, exact p 17/126; ρ 0.518 [−0.247, 1.361]; the positive D rests on AML162 |
| (b) children | "not readable" | A_b 0.372 [0.161, 0.594]; φ −0.244 [−0.827, 0.092] |

**1.4 Seen while writing this draft.**
- The September roster file (`pb_sample_manifest.csv`, sample-level) and `power.csv` (bin-level) were read.
- For every other file under `results/tables/12_pseudobulk_de/`, only names, sizes, headers, row counts and md5 sums
  were read.
- T_NK cell counts per sample were counted from the reconciled per-cell label columns (inventory A 3.2; no
  expression).
- S5-side outputs were read (inventory B). By accident, they include the S5 logFC of the 32 ambient-prone genes.
- No September gene name, logFC, p or q was read, except AC007952.4's values, which are already on record (G04 1.1).

**1.5 Disclosures** (inventory B 3(c); each goes into the paper's limitations in the form decided at writing).
1. The September pre-registration promised not to reframe its null as "underpowered" or as "needs a different model"
   (`PREREGISTRATION_pseudobulk_de.md`:41-44). G03 is bound by it (section 0).
2. S5's hit rule (FDR 0.10, |logFC| 0.5) differs from September's (q 0.05, |logFC| ≥ 1). The handoff (v2.0 :137,
   §16.2 item 5) suggested S5's rule before S5 ran, but gave no reason.
3. S5 has no library term, although the September write-up notes that without one "the healthy n is a fiction"
   (`FINDINGS_pseudobulk_de.md`:461-463). A1-BATCH G2 later tested lane clustering on the cross-fitted gap.
4. T_NK entered the September analysis only through the post-hoc marker-gate amendment (September pre-registration
   §8.2). The September write-up says its pre-registration was amended twice and corrected once, each time after seeing
   results.
5. The September fit pooled children with adults because no age column existed then (September pre-registration
   :422-424). Handoff v2.0 forbids that pooling.
6. September p values are Monte Carlo over 1,000 draws. The September write-up calls them "exact".
7. Reading T_NK's single hit as zero is a FINDINGS reading (`FINDINGS_pseudobulk_de.md`:141-142), not a pre-registered
   rule. The gene stays in `hitlist_frozen.csv`.
8. The September universe is a symbol intersection. AC007952.4 pools two Ensembl loci.
9. The blind-safe rungs were not run before G04, as the addendum had recommended. They are specified here with G04's
   result known.
10. 24 of S5's 25 samples are in the September fit, and both per-gene tables are public. "Unread" is attested, not
    enforced.

## 2. Data

**2.1 The 24 shared adults** (B 1.1; A 3.3). These are S5's 25 without AML0310.
- AML0310 is not in the September roster: it is CCC-ineligible (`below_min_occupied_bins`). Under the September filter
  its T_NK bin also holds 29 cells, below the 30-cell gate.

| arm | samples |
|---|---|
| AML (14) | AML0048, AML0114, AML0612, AML0693, AML1371, AML2123, AML2451, AML3133, AML3266, AML3730, AML3762, AML3948, AML4340, AML4897 |
| healthy (10) | Control0004, Control0005, Control0058, Control0082, Control1, Control2, Control3, Control4, Control4003, Control5 |

Libraries: 17.
- Multi-sample libraries: P01 (AML0612, AML3133, AML3762); P02 (Control1-4, a T-ALL pool); P05 (AML2451 and
  Control0004, the only mixed one); P15 (AML1371, AML4340).
- Every other sample is its own library.
- Control4 is female by expression and male in the metadata.

**2.2 Inputs** (read only; md5 checked in K0).

| input | used in | notes |
|---|---|---|
| S5 pseudobulks, `/LARGE1/.../ws_pseudobulk/GSE185381/ws_pb_GSE185381.rds` (`$raw$T`, `$raw$NK`) | R0-R3 | |
| S5 design table (`27_s5_imprint.R`:55-67 logic; rows in `ws_s5_scores.csv`) | R0, R1 | |
| September per-sample pseudobulks, `/LARGE1/.../12_pseudobulk_de/<dataset>/<sample>__pseudobulk.rds` (T_NK bin column) | R4-R6, E1, E2 | the stored files, not re-aggregated (the Chen2023 QC objects were rewritten 2026-10-05; A 3.1) |
| September depth tables `results/tables/12_pseudobulk_de/depth/*.csv` | R4-R6, E1, E2 | n_cells, median_umi, total_counts |
| `library_id` from `pb_sample_manifest.csv` | R4-R6, E1, E2 | for R1-R3 the same IDs, joined by sample |

**2.3 Extension samples** (G04 5.7).
- E1 adds the 18 children of G04(b): AML001, AML003, AML005, AML006, AML009, AML012, AML022, AML025, AML026, AML028,
  AML038, AML043, AML048, AML051, AML052, AML055, AML056, AML073.
- E2 adds Chen2023's AML103, AML162, AML24, AML304 and AML320, and NBM1-4. The pseudobulks are September's (both pools,
  projection bin, symbol axis).
- PAUMTZ and PAXMIJ stay closed: no column of theirs is read (G04 decision 9.4). Their September files exist, so the
  script drops them by name before loading.

## 3. The ladder

### 3.1 Rungs (each changes one factor; decision 9.1)

| rung | samples | cells | design | engine (and its permutation procedure, 3.2) | gene universe |
|---|---|---|---|---|---|
| R0 | S5's 25 | S5: concordant T, ≥ 20 | `~ status + hb + gran + sex + CD8 share` | edgeR QL robust | S5 (blacklist; filterByExpr(group); CPM > 1 in half) |
| **R1** | **24 (AML0310 out)** | as R0 | as R0 | as R0 | as R0 |
| **R2** | 24 | as R0 | **`~ status`** | as R0 | as R0 |
| **R3** | 24 | **concordant T + NK** (`$raw$T + $raw$NK`), ≥ 20 | `~ status` | as R0 | as R0 |
| **R4** | 24 | **September T_NK bin** (reconciled bin, `in_ccc_graph`, not `high_error`, Early Lymphoid out), ≥ 30 | `~ status` | as R0 | as R0 |
| **R5** | 24 | as R4 | `~ arm` + `library_id` block | **limma-voom, two passes of duplicateCorrelation** (September 02:157-167) | as R0 |
| **R6** | 24 | as R4 | as R5 | as R5 | **September base: the pinned symbol intersection Chen2023 ∩ GSE185381 minus 91 sex genes (21,843), then `filterByExpr(design)` only** |

Notes:
- R0 is a known answer, not a rung (K1). Its result is on record.
- R6 is "the September specification on the 24 S5 samples" (G04 5.7), with the September hit rule read as the second
  count (3.2).
- R5 and R6 use `~ arm` because a one-level `dataset` factor fails (A 3.4). They keep the `library_id` block, because
  GSE185381 has shared lanes (decision 9.5). The September pre-registration's single-dataset form drops the block, but
  only on the grounds that libraries are singletons there.
- The September §5.2 bin rules fail by construction for a one-dataset fit and are not applied. The §8.2 marker gate is
  reported, not applied (5.6).
- The order follows the addendum verifier's ladder, except that the gene universe and the hit rule are no longer one
  rung. The universe is R6; the hit rule becomes a second count at every rung (3.2; decision 9.2). AML0310 leaves first,
  so that R1-R6 hold the samples fixed (decision 9.1).

### 3.2 What is computed at every rung

**Two counts from one fit.**
- **c_S5**: genes with BH q < 0.10 and |logFC| > 0.5 (S5's rule).
- **c_Sep**: genes with q < 0.05 and |logFC| ≥ 1 (the September baseline rule; this is also the September permutation
  statistic).

Each count gets two Monte Carlo p values, `(1 + #{null ≥ observed}) / (1 + valid draws)`, 1,000 draws each:
- **p_free**: free relabelling of the samples, arm sizes kept. This is S5's null and the September primary null (within
  dataset; for a one-dataset rung it is the same thing).
- **p_lib**: whole arm-pure libraries relabelled, mixed libraries held fixed (September 03:116-127). It permutes
  library units, so arm sizes in samples may vary.
- For the 24 there are 1,961,256 free arrangements and 8,008 = C(16, 6) library arrangements (A 3.4).

Permutation procedure (decision 9.4):
- **edgeR rungs (R1-R4): full refit.** Each draw reruns `filterByExpr`, TMM, `estimateDisp` and `glmQLFit` (S5's
  procedure).
- **limma rungs (R5, R6, E1, E2): September's frozen scheme.** The gene set, voom weights and duplicateCorrelation
  consensus are kept from the observed fit; only `lmFit` and `eBayes` are rerun (September 03:149-160). The September
  code itself notes that this scheme leans towards calling the observed count unremarkable.
- 4.1 adds full refits at R5 and R6.

Common random numbers:
- One set of 1,000 free draws (seed `SEED + 61`) and one set of 1,000 library draws (`SEED + 62`) is drawn once for the
  24, in the S5 design order, and reused at R1-R6.
- Rung-to-rung differences therefore do not include draw noise.

Also per rung, descriptive:
- genes tested;
- duplicateCorrelation consensus (R5, R6);
- c_Sep after the depth arm (R4-R6, E1, E2; September 02:190-197). "Reportable" means a hit in both the baseline fit
  and the fit with `log10(median_umi)` added.
- 5.1's P370 readouts.

### 3.3 Extension rungs (fixed by G04 5.7; nothing else changes)

| rung | samples | design | libraries |
|---|---|---|---|
| E1 "+children" | R6's 24 + the 18 children = 42 (32 AML v 10 healthy) | `~ arm` + block | 22 (19 pure, 3 mixed: P05, P12, P18) |
| E2 "+Chen2023" | E1 + 9 Chen2023 = 51 (37 v 14) | `~ dataset + arm` + block (exactly September 02:147-163) | 31 (28 pure, 3 mixed) |

- Gene universe: the pinned 21,843, then `filterByExpr(design)`.
- Permutations: September's procedure, 1,000 draws per null (E1: free `SEED + 63`, library `SEED + 64`; E2: free within
  dataset `SEED + 65`, library `SEED + 66`).
- Readout: c_Sep with p_free and p_lib, c_Sep after the depth arm, and 5.1's P370 readouts (G04 5.7). c_S5 is reported
  beside them.
- The September fit (53, 1 hit) is the reference. E2 differs from it only by PAUMTZ and PAXMIJ, and therefore P06's
  library purity.

## 4. Sensitivities (each reported; none changes a category)

- **4.1 Full refits at R5 and R6.** voom and both duplicateCorrelation passes are rerun per draw, with the same draws
  as 3.2. This separates the engine from the frozen permutation scheme.
- **4.2 Control4 coded female at R1.** The sex covariate exists only at R0 and R1.
- **4.3 The children's author-malignant T cells (E1).** AML073 has 210 of 605 kept T cells labelled malignant by the
  authors. The September cell set did not remove them, and E1 keeps them ("nothing else changes"). As a descriptive
  check, E1 is rerun without AML073: c_Sep and c_S5, with no p.

Path dependence is a limitation, not a sensitivity. A factor's effect can depend on the order in which factors change.
The reverse ladder (from September's end) is not feasible for every factor: a T-only projection-bin pseudobulk does not
exist without new per-cell aggregation (B 5.2).

## 5. Readouts (descriptive; no category)

**5.1 P370 per rung.** Among the P370 genes in the rung's tested universe:
- the number tested;
- hits under each rule;
- the share with S5's sign;
- the median of `sign_S5 × logFC_rung / |logFC_S5|` (attenuation).

At R1-R6 this is circular: P370 was selected on 25 samples that include these 24. Flag it as circular. E1 and E2 add
samples that P370 never saw.

**5.2 AC007952.4 per rung.** Its logFC and q, plus whether it is in the universe. On R0-R5's GSE185381 axis it is one
locus. At R6 and after, the symbol intersection pools two loci.

**5.3 The stored September table** (`de/primary__T_NK.csv`, opened here for the first time; md5 checked first).
- Among P370 genes in its 9,281: the number of genes, hits, the share with S5's sign, and the median `sign_S5 × t`.
- No test; labelled circular, because 24 of 53 samples are shared (decision 9.9).
- The matched-random-set comparison suggested in the assessment is not run: P370 was selected on the shared samples,
  and unselected random genes cannot calibrate that selection.

**5.4 The prediction updated with G04's estimates** (computed from on-record numbers; no September per-gene value).
- The addendum's calculation is repeated: `power.csv`'s full-rule recovery curve, interpolated linearly, applied to each
  P370 gene's S5 |log2FC| scaled by the AML arm's carrier fraction f.
- f = (14 + 5ρ + 20φ) / 39, with ρ and φ from G04 (point, and the four corners of their CIs; φ is applied to all 20
  children, PAUMTZ and PAXMIJ included).
- The rows are the expected number of September calls among the 370, unshrunk and shrunk by r = 1.43, against the
  observed 1.
- The same rough-upper-estimate caveat as the addendum applies (mid-tertile genes planted into every AML sample).
- A_b, φ and their CIs are carried into G03's GATES block as G04 handed them over.

**5.5 Planted recovery per rung (R6, E1, E2).**
- September's planted grid: 200 middle-tertile genes, logFC 0.5/1/1.5/2/3 added to the voom E of every AML sample,
  frozen weights (September 03:202-224; seeds `SEED + 67`, `+68`, `+69`).
- Reported: recovery under the September rule and LFC80 per rung, i.e. the licence that rung's null carries ("no effect
  of |logFC| ≥ LFC80 in every AML sample").

**5.6 Marker gate per rung (R4-R6, E1, E2).** The September §8.2 amended gate (CD3E/CD3D: mean CPM ≥ 100 in both arms,
|log2 arm ratio| ≤ 1), reported only.

## 6. Known answers (run first; three states as in G04)

| ID | scope | check |
|---|---|---|
| K0 | all | md5 of the declared inputs: inventory A section 5's September files, the 51 rung pseudobulk RDS files, `ws_pb_GSE185381.rds` (md5 6c2fb685…), S5's gene and score tables |
| K1 | L | R0 reproduces S5: 7,379 tested, 402 hits (171 up / 231 down), P370 membership |
| K2 | L, E | for every rung sample, the stored T_NK column sum equals the depth table's `total_counts`, and `n_cells` ≥ 30 equals `passes_30`. The pinned base universe rebuilt from gene names is 21,843 |
| K3 | L, E | engine and permutation code on a September fit that holds no closed sample: the Chen2023 secondary bin HSC_MPP (5 v 4, `~ arm`). Genes tested and baseline hits equal `de_summary.csv`; max \|ΔlogFC\| vs `de/secondary__HSC_MPP.csv` < 1e-8; the sample-null summary row of `perm_summary.csv` is reproduced from its own seed label |
| K4 | L, E | the blocked two-pass fit equals September 02's `fit_bin` (copied verbatim, run side by side) on synthetic counts with a shared-library structure |
| K5 | L, E | permutable units: 1,961,256 / 8,008 (24), and A 3.4's counts for E1 and E2 |
| K6 | L | each rung's per-gene logFC vector differs from the previous rung's (a no-op rung signals a bug; September §8.6's "differs from baseline" tell) |
| K7 | L | R4's 24 T_NK columns equal the September files' columns, i.e. nothing is re-aggregated |
| K8 | all | smoke run of every mode on synthetic counts (no real covariate written next to a score) |

A deterministic failure is fixed, a dated note is added, and everything is rerun. A failed known answer of scope L
makes the ladder "not decisive"; of scope E, it makes the extension readout "not run".

## 7. Reading rules (fixed now; first match wins)

**Ladder (R1-R6).** "Held at rung k" means c_S5's p_free ≤ 0.05 at k (decision 9.3). p_lib is reported beside it in
every sentence.
1. **"not decisive"**: a scope-L known answer failed.
2. **"held to the September design"**: held at R6.
3. **"lost at Rk"**: otherwise, where k is the first rung not held. If a later rung is held again, add "regained at Rj".

Wording:
- 2 → "On the 24 adults S5 shares with the September fit, S5's rule keeps c_S5 hits beyond its permutation null under
  the September cell set, engine and gene universe (p_free, p_lib). Under the September rule the count is c_Sep
  (p_free, p_lib). The two analyses count different effects: S5 counts |logFC| > 0.5 at FDR 0.10, September |logFC| ≥ 1
  at q 0.05. The September T_NK result stays reported as a null for |logFC| ≥ 1.5 per gene."
- 3, k ∈ {R2, R3, R4, R6} (specification rungs) → "S5's per-gene count needs <factor>: without it the count falls
  inside its permutation null (rung table). The per-gene list is reported with this qualifier."
- 3, k ∈ {R1, R5} (validity rungs) → R1: "S5's per-gene count depends on one sample (AML0310)". R5: "S5's per-gene
  count does not survive library-level correlation". In both cases: "the per-gene list is reported as exploratory; the
  sample-level separation (A1-G02, A1-BATCH G2) is a different statistic and stands as reported there."
- The rung classes and their consequences are decision 9.7.

**September rule on the 24 (R6, c_Sep).** If p_free ≤ 0.05: "the September specification finds hits beyond its null
on the 24 adults; the September null then traces to the added samples (E1, E2)". Otherwise: "the September
specification finds none beyond its null on the 24 adults, as on the 53".

**Extension (E1, E2).** No category. The sentence reports c_Sep (p_free, p_lib) and the P370 hits at R6, E1 and E2,
then 5.4's updated prediction beside the observed 1.

**Never written, under any outcome:**
- "September was underpowered" or "used the wrong model";
- "the imprint is confirmed by the September data";
- a reading of 5.1 at R1-R6, or of 5.3, as evidence (both circular).

**Handover.**
- `--combine` appends a GATES.md block with the ladder category, the rung table and the fixed sentences.
- A dated note goes in handoff S15.
- G03 has no pass/fail gate of its own. It decides the status of the per-gene list (section 0) for G07 and for writing.

## 8. Code, outputs, compute (plan; the script is written after the tag)

- One script, `scripts/ws/42_a1_g03.R`, with modes `--known`, `--ladder`, `--extend`, `--sens`, `--read` (5.3, 5.4)
  and `--combine`, each also runnable with `--smoke`. A wrapper `submit_ws_a1_g03.sbatch` runs dated read-only copies
  (as in G04).
- Functions are copied, never sourced. September 02/03 and ws/26, 27 and 37 write frozen public tables and must not be
  run (B 6.1).
- Outputs:
  - `results/tables/ws/ws_a1_g03_{checks,ladder,extend,sens,p370,read,summary}.csv`;
  - permutation draws and null counts under `/LARGE1/.../ws_a1_g03/`.
- Every target is checked at mode start, and a mode refuses to overwrite. A later mode stops unless the checks table
  carries the same script md5.
- Compute, to be confirmed by timing in the smoke run:
  - R1-R4: about 8,000 edgeR refits, about 9 CPU-h;
  - the frozen limma rungs: under 1 CPU-h;
  - 4.1: 4,000 full limma refits, 20-70 CPU-h depending on duplicateCorrelation time.
  - All mclapply jobs set `OMP_PROC_BIND=false`. A SLURM array is used if 4.1 exceeds 12 h on one node.
- Seeds: `SEED = 491638`, offsets +61 to +69, all unused in `scripts/ws`.

## 9. Decisions for the user

Real choices are marked ★; the rest have a conventional default.

| # | decision | recommended | alternative |
|---|---|---|---|
| 9.1 | rung order; AML0310 leaves first (R1) | as in 3.1 | AML0310 leaves at R4, where the September cell set first excludes it |
| 9.2 | the hit rule is a second count at every rung, not its own rung | yes | a final rung R7 with the September rule only |
| 9.3 ★ | the null that decides "held" | p_free (S5's and September's primary null); p_lib reported | p_lib (stricter; it also tests lane structure, which A1-BATCH G2 already tested on the gap) |
| 9.4 | permutation procedure | each engine's own (edgeR full refit; limma frozen as September); 4.1 full refits at R5 and R6 | full refit everywhere (20-70 CPU-h more) |
| 9.5 | one-dataset rungs keep the `library_id` block | yes | `~ arm` without block (the September single-dataset form) |
| 9.6 | R6 universe | September's pinned symbol base (21,843) | Ensembl matching (not September's specification) |
| 9.7 ★ | rung classes: validity = R1, R5; specification = R2, R3, R4, R6; a validity loss makes the per-gene list "exploratory" | as stated | no classes: every loss only names the factor, with no consequence fixed now |
| 9.8 | 4.3 (E1 rerun without AML073): descriptive | include | drop |
| 9.9 | the stored September table (5.3): descriptive only, no random-set test | yes | matched random sets as in the assessment (circular either way) |

## 10. Review history

- **DRAFT v1**, 2026-11-05: written from inventories A and B. Not yet reviewed.
