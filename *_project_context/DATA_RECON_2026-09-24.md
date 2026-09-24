# Data reconnaissance and competition check — new project direction

**Date:** 2026-09-24. **Purpose:** decide whether the "genotype legibility" direction is viable before
committing to it, and before pitching it.

**Candidate direction.** *A hierarchy-resolved map of where AML driver genotypes are transcriptionally
legible.* Three layers:

| | question | instrument |
|---|---|---|
| **L1** | For each driver lesion, which hematopoietic compartment carries an **independent** transcriptional imprint, and which compartment signals are shadows of the blast compartment — measured annotation-free at compartment pseudobulk | `12_pseudobulk_de` (P1), HOX scoring (P2), compartment adjustment (N6) |
| **L2** | Transfer compartment programs to bulk cohorts; does genotype–program **discordance** carry prognostic / drug-response information beyond genotype | new; open bulk cohorts |
| **L3** | Does the compartment program follow the lesion | mouse driver models (human isogenic does not exist — §4) |

**How this document was produced.** Eight parallel web-research agents, each output independently
re-verified by a second agent that re-searched from scratch and re-downloaded primary files, then a
completeness critic and four gap-fill agents. 21 agents, 791 fetches. Every number below was
confirmed against a primary source in that session; items that could not be confirmed are marked.
Verifiers overturned four first-pass conclusions — those corrections are flagged **[corrected]**.

---

## 1. Verdict

**The direction is viable. The lane is open but narrow, and it is closing.** Two 2025–2026 papers
now flank it on both sides, and one 2023 paper has already executed the whole L1→L2 shape in mouse.
What remains unoccupied is a specific composite, and the pitch has to be pinned to exactly that
composite — not to any one of its parts, each of which is taken.

**Still unoccupied, as of 2026-09-24:**

> A per-driver × per-compartment legibility matrix, measured **annotation-free** at cohort scale,
> with explicit **independent-signal vs blast-shadow adjudication**, plus **systematic multi-driver,
> both-direction, composition-corrected discordance** tested against survival **and** ex vivo drug
> response.

**Already taken, and each must be cited rather than rediscovered:**

| claim | owner |
|---|---|
| genotype → hierarchy **composition**, cohort scale, on our own reference | Zeng 2025 *Blood Cancer Discov* |
| composition → prognosis + drug response | Zeng 2022 *Nat Med*; Karakaslar/seAMLess 2024 |
| which driver programs are "legible", **reported per AML cell type** | **scOPE** bioRxiv 2026-08 |
| single-driver discordance → survival **and** venetoclax resistance | **Lee 2024** *Blood Cancer J* (TP53) |
| within-blast program made legible in bulk, clinically consequential | Lilljebjörn 2025 *Nat Commun* (NPM1) |
| recovering **within-cell-type expression** from bulk, with prognosis | **BLUE** *PLoS Comput Biol* 2026 |
| per-driver program in a compartment → human bulk survival | **Isobe 2023** *Cell Genomics* (mouse, preleukemic) |
| genotype-from-expression classifiers | commodity — see §5.4 |

Three consequences that change the plan, not just the wording:

1. **P2 (NPM1 AUC 0.975) is not a contribution.** Genotype-from-expression in AML is a solved
   commodity: IDH AUROC 0.994 (Jung 2026), NPM1 0.971 (scOPE 2026), a seven-model benchmark
   (Silva 2024). Only *where the signal lives* and *that it is independent of composition* are
   defensible.
2. **Zeng 2025's composition models are the mandatory null.** Any program claim a composition model
   reproduces is dead on arrival. L1 must be framed "programs **given** composition", never "which
   compartments expand".
3. **Do not headline a new prognostic score.** That genre is crowded (Zhang 2025 ACEsig,
   Severens 2024, Karakaslar 2024, Bertoline 2026 and more). Discordance must be the object, and the
   comparison must be against genotype plus composition, not against other scores.

---

## 2. L2 anchors — all open, no application needed

This is the strongest part of the picture. Everything needed for L2 downloads from a Japanese HPC
with zero data-access applications.

### 2.1 Primary: Beat AML waves 1–4 (Bottomly 2022 *Cancer Cell*)

`https://biodev.github.io/BeatAML2/` — CC-BY-4.0, flat files, no registration. Counts below were
recomputed from the official clinical workbook by two independent agents and agree exactly.

- **942 specimens / 805 patients.** RNA-seq 698 (671 analysis-grade), WES/targeted 903,
  ex vivo drug 631 specimens / 569 patients.
- **Usable overlaps: RNA∩drug 542; RNA∩WES∩drug 514 specimens / 471 patients.**
- Drug panel **166 single-agent inhibitors**, probit-fit AUC + IC10/25/50/75/90.
- Survival in all 942 rows (Dead 565 / Alive 345 / Unknown 32); ELN2017 all rows; karyotype 842/942.
- Genotype: full mutation table (11,721 variants) **plus** consensus per-specimen FLT3-ITD with
  allelic ratio, NPM1, RUNX1, ASXL1, TP53, CEBPA-biallelic and fusions inside the one open clinical
  workbook.
- Two practical traps: several files are **Git-LFS** — `/raw/` GitHub URLs silently return 133-byte
  pointers, use `media.githubusercontent.com/media/biodev/beataml2.0_data/main/<file>`. And the
  biodev matrix is **GRCh37** while GDC's open quantifications are GRCh38 — pick one.
- cBioPortal `aml_ohsu_2022` and the AWS Open Data mirror are alternative channels. Raw reads need
  dbGaP phs001657 (not needed for L2). No wave-5 release exists as of today.

### 2.2 Second drug cohort: FIMM / Helsinki (Malani 2022 *Cancer Discov*)

**Zenodo 10.5281/zenodo.7274740** — open CC-BY, one 30.9 MB zip, *not mentioned in the paper's own
data-availability statement*. Verified by downloading and opening every file:

- sDSS **515 compounds × 181 samples**; log2CPM 18,202 genes × 167 samples; VAF 340 genes × 225;
  binary calls for 57 recurrent genes; clinical table for 186 patients.
- **Survival fields are present** (`os_time`, `os_event`, `efs_*`, `rfs_*`, `eln2017_risk_class`).
- **RNA∩DSRT = 132 samples (109 patients).**
- Relapsed/refractory-weighted; **sDSS is healthy-normalized, not the Beat AML probit AUC** — needs
  per-cohort rescaling, never pooling.

### 2.3 Survival cohorts

| cohort | n | genotype openly per sample | survival | notes |
|---|---|---|---|---|
| **GSE146173** (AMLCG) | 246 | **246/246** — NPM1 82, FLT3-ITD 49, FLT3-TKD 22, TP53 19, RUNX1 48, DNMT3A 78, IDH1 19, IDH2 31, ASXL1 34, TET2 41, CEBPA 9, EZH2 4, SRSF2 26, U2AF1 11, SF3B1 7 + hotspots | **OS 246/246**, EFS 246/246, RFS 155/246 | **the richest open cohort found.** Cytogenetics 238/246, ELN2017 238/246 |
| TCGA-LAML | 200 cases | open masked somatic MAFs | verified populated | **RNA-seq only 151** (183 are array) |
| GSE37642 (AMLCG-1999) | 562 patients / 984 GSMs | RUNX1 only | OS 553/562 | microarray, two platform generations |
| TARGET-AML | 2,492 cases | mostly controlled (phs000465) | OS verified | expression + clinical fully open; **pediatric spectrum** |
| HOVON GSE6891 / GSE14468 | 537 / 526 | **NPM1, FLT3-ITD/TKD, CEBPA, IDH1/2, N/KRAS, karyotype, risk group** | **not in GEO** | **[corrected]** — first pass wrongly called genotype absent, from one unrepresentative sample |
| Leucegene GSE232130 | 691 | none in GEO (BCLQ-controlled) | none in GEO | **[corrected]** 691, not 1,038 (452 re-analyzed + 239 new). Expression-only |

**Cohort-independence trap, resolved by direct patient matching:** the three AMLCG series are not
three validations. **237 of GSE146173's 246 patients are the same individuals as in GSE106291**
(identical GEO ids, and gender/age/OS/status matching 237/237). GSE12417 sits inside the GSE37642
family. GSE37642 shows **no** detectable overlap with the RNA-seq series, so it can serve as a
second, platform-independent AMLCG cohort.

> **Recommendation:** BeatAML2 (drug) + GSE146173 (genotype + survival) + TCGA-LAML, with FIMM as the
> independent drug replication and GSE37642 as an optional microarray arm. Never count GSE106291 and
> GSE146173 separately. German AMLCG patients cannot overlap US cohorts — that independence is
> structural.

**GSE146173 caveat to carry into analysis:** resistance-enriched by construction (all 41
AMLCG-1999 add-ons are non-responders), so it is valid for genotype-conditioned association and
invalid for absolute survival calibration. Its Lexogen 3'-biased chemistry also means transfer
to/from polyA TruSeq cohorts needs rank-based scoring.

### 2.4 Cell-line drug anchors (optional)

**[corrected] — GDSC was wrongly reported unavailable.** `GDSC2_fitted_dose_response_27Oct23.csv`
downloads openly from `cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/` (the old cancerrxgene
download page now returns 410). 242,036 curves, 969 lines, 286 drugs, **26 AML lines with multi-dose
LN_IC50 and AUC** — the same functional form as Beat AML AUC, unlike PRISM's single-dose LFC
(21 AML lines). Union of both: 35 AML lines. Sanger Project Score adds nothing (only 2 AML lines).

---

## 3. L1 anchors — cell-level truth is better than we thought

Our 244-sample cohort stays the discovery set. What the recon adds is external truth for validating
annotation-free compartment calls.

- **GSE230559** (Cell Stem Cell 2025) — **[corrected]**: the first pass concluded "no GoT-style
  genotype+transcriptome cohort in de novo AML exists"; the verifier refuted it. **20 IDH-mutant de
  novo AML patients** (6 IDH1, 13 IDH2, 1 both), integrated single-cell genotyping + transcriptomics,
  clone-specific programs for NPM1 / NRAS / SRSF2 co-mutations, clone-level inhibitor response.
  Open GEO.
- **nanoranger GSE243227** — per-cell mutation calls reported **per hematopoietic compartment**
  (HSC 847, LMPP 1,198, GMP 154, MK 815, Ery 3,223). This is literally our question, at n≈5 patients.
- **TARGET-seq GSE226340 + Zenodo 8060602** — 17,517 HSPCs, 14 TP53-sAML donors; per-cell curated
  genotype labels openly downloadable. Lin⁻CD34⁺ only, so it anchors HSPC compartments only.
- **LOTR-Seq E-MTAB-15981** — multi-locus per-cell genotype covering IDH1/2, TP53, **SRSF2, U2AF1**
  (splicing-factor coverage is otherwise absent), and it already projects onto BoneMarrowMap.
- **GSE158067** — DNMT3A R882 clonal hematopoiesis, 5 donors, CD34⁺, per-cell genotype +
  transcriptome + methylome, donor as own control. CH not AML, but the cleanest internal-control
  design available.
- **CloneTracer** (already in our 13) and **MutaSeq** (4 patients, figshare) complete the set.
- **E-MTAB-16321** — fully open adult CITE-seq atlas, 26 patients × diagnosis+relapse = 52 samples,
  81 surface antigens, per-patient cytogenetics/mutations/survival. The best-access new adult cohort.
- **Lambo GSE235063** — 28 pediatric patients, 684,031 cells; **[corrected]** the per-cell metadata
  ships `Expected_Driving_Aberration`, `Subgroup`, `Known_CNVs`, `Malignant` and
  `Clinical_Blast_Percent` openly, no paywalled table needed. Pediatric spectrum, separate stratum.

**Pre-registered expectation worth noting:** CloneTracer found the disease-defining aberrant
compartment to be downstream myeloid progenitors, with near-normal active LSCs. If our legibility
matrix concentrates signal in LMPP/GMP-like compartments, that is a replication, not a discovery —
plan the claim accordingly.

---

## 4. L3 — the human isogenic route does not exist; mouse does

**Twice-checked negative: there is no public human isogenic NPM1c, FLT3-ITD or IDH knock-in with
transcriptomics.** Searched from several angles by both agents. The human isogenic literature covers
the splicing-factor / epigenetic-modifier / CH axis (Papapetrou GSE163034: ASXL1+SRSF2+NRAS;
GSE164666: SRSF2/U2AF1) and not the NPM1c–FLT3–IDH axis — which is exactly where our strongest
findings sit.

**What exists instead:**

- **Mouse driver-model scRNA is strong and fully open.** Isobe **GSE227026** covers 8 lesions
  including Flt3-ITD, Npm1c, Dnmt3a-R882H and Idh1-R132H — 38 animals, 269,048 cells, Lin⁻cKit⁺ only.
  **[corrected]** the split is 19 mutant / 19 WT per GEO, not 21/17.
- **STRACK GSE266232** (Cell Stem Cell 2025) — LARRY-barcoded HSCs tracked **before and after**
  acquiring Dnmt3a-R878H or Npm1c. The single cleanest public "does the program follow the lesion"
  design for our exact lesions. Missed on the first pass.
- **SPLINTR GSE161676** — **[corrected]**: in vivo leukemic-stage MLL-AF9 mouse scRNA with clonal
  barcodes, including a chemotherapy arm, does exist (first pass said it did not).
- Also: GSE318348 (Npm1c;Flt3-ITD LARRY + CROPseq + chemo), Izzo GSE124822 (Tet2/Dnmt3a KO,
  Idh2-R140Q), GSE272266 (Dnmt3a R878H point mutant LSK).
- **NPM1c degron systems** test maintenance, not compartment: four independent systems
  (GSE111178/111180, GSE251919, GSE306285, and the missed **GSE197387** SuperSeries with CUT&RUN /
  ChIP / Bru-seq that addresses whether HOX is a *direct* NPM1c effect). All built on the **only two
  NPM1c-mutant human AML lines in existence** — OCI-AML3 and IMS-M2 — so cross-dataset replication is
  weaker than it looks, and OCI-AML3 is also DNMT3A-R882C. **[corrected]** GSE251919 ships **no**
  processed matrix; raw reads only.
- **DepMap is a dependency anchor, not a hierarchy anchor**, and its AML panel cannot support
  per-driver contrasts: of 61 AML lines, NPM1c n=1, canonical IDH1/IDH2 **n=0**, SRSF2 n=1,
  U2AF1 n=2, FLT3-ITD 4 lines but 3 donors; only TP53 (33/61) has n. Two landmines:
  the hotspot convenience matrix **has no NPM1 column at all** and misses FLT3-ITD, and the full
  mutations CSV has **embedded newlines inside quoted fields** — line-based filtering produces
  confident false zeros (one verifier hit exactly this and caught itself).
- **Replogle Perturb-seq is unusable for our drivers**: CRISPRi knocks *down* wild-type genes while
  NPM1c/FLT3-ITD/IDH/DNMT3A-R882 are gain-of-function, neomorphic or dominant-negative; FLT3 is not
  in the library; K562 is CML blast crisis, and carries ASXL1-truncating on top of TP53-null.
  Tahoe-100M has **zero** hematopoietic lines.

> **L3 decision:** route NPM1c/FLT3-ITD/IDH through **mouse** (GSE227026 + GSE266232), and treat the
> degron systems as maintenance-level support only. Keep L3 optional — it is a strengthening layer,
> not a gate.

---

## 5. Competition, in detail

### 5.1 The two flanking papers

**Zeng 2025** (*Blood Cancer Discov* 6(4):307–324) is the closest competitor and it is built on **our
own reference**: BoneMarrowMap, 1,223,411 cells from 318 samples, deconvolution of **1,224 bulk
samples across 5 cohorts**, >45 drivers mapped to 12 aberrant-differentiation patterns, plus NK-AML
survival and BeatAML drug associations. Verified by two independent full-text reads: the
genotype–phenotype unit is **state abundance**, and there is **no within-state differential
expression by genotype**. Reuse their associations as positive controls and nulls; do not re-derive.

**scOPE** (Ashford, Lapadat & Demir, bioRxiv 2026-07-26 / v2 2026-08-12) is the nearer threat to L1's
*rhetoric*. Full methods read: bulk-trained per-driver classifiers, frozen, projected onto single
cells across 7 cancers with AML as flagship (NPM1 AUROC 0.971, SRSF2 0.960, CEBPA 0.929, IDH2 0.910);
158 driver–cancer models audited, only 11 passed. It already uses **"legibility"** language, already
residualizes against same-state healthy medians, and **Fig 6D is literally "Mean residual score by
AML cell type and driver"**. What it does **not** do: pseudobulk anywhere (0 occurrences), any test
of independent-vs-shadow signal per compartment, any survival or drug layer (0 occurrences of
survival/prognosis/venetoclax/ex vivo), and its "discordance" means score-vs-CNV, not
genotype-vs-program. Scale is one cohort — van Galen, 16 patients, 35,843 cells.

### 5.2 Discordance is no longer virgin ground

**[corrected]** — the first pass concluded "discordance-predicts-outcome has not been done in AML".
The verifier refuted it with **Lee, Baughn, Myers & Sachs 2024** (*Blood Cancer J* 14:80): a ridge
classifier of TP53 status on Beat AML, top-10% TP53-**wild-type** scorers (n=40) defined as
"TP53Mut-like", **median survival 204 vs 861 days**, broad ex vivo drug resistance including
venetoclax, validated in TCGA. Single driver, whole-bulk, one direction, no composition adjustment —
but it is a published instance of the core idea and must be cited as such.

Supporting precedent: Mer 2021 (37% of FLT3-WT NPM1-mutant samples carry an ITD-like transcriptome
with equally poor prognosis — a by-product, never the modeled variable); Mosquera Orgueira 2021
(FLT3-like WT, no outcome); Bakhtiar 2022 (pan-cancer phenocopies improve drug-response prediction in
68% of 165 combinations — concept taken, no AML, no survival); Helzer 2024; the Ph-like ALL
literature as the leukemia-wide precedent.

**Two cautions to build in rather than discover in review.** First, **DeltaNp73**
(Pereira-Martins, *Cell Rep Med* Jan 2026): TP53-WT AML with high DeltaNp73 phenocopies TP53-mutant
outcome *with an identified molecular cause* — reviewers will ask whether our discordance is just an
unmeasured lesion. Second, **Lilljebjörn 2025** treats classifier-discordant cases as a technical
artifact of low blast content (one third of their external NPM1 cases were unclassifiable for that
reason). Our framing has to do their correction first and treat only the residual as biology.

### 5.3 The competitor that did our whole shape, in mouse

**Isobe 2023** (*Cell Genomics* 3:100426) — found only by a verifier, missed by two agents. Eight
mouse driver models, 269,048 HSPC transcriptomes, **per-mutation programs within the HSPC compartment
transferred to TCGA / BeatAML / TARGET for survival** (the Stem11 signature). This is the L1→L2 shape
executed end to end. It is mouse, preleukemic, one compartment, and has no discordance analysis — but
it is simultaneously our best L3 dataset and the proof that this shape publishes. Cite it as both.

### 5.4 Two more boundary stones

**BLUE** (*PLoS Comput Biol* 2026, Zhu & Qiu) recovers **per-sample cell-type-specific expression**
from bulk with AML as its primary case study and derives survival-stratifying subtypes — with **no
genotype link**. This falsifies any sentence of the form "nobody recovers within-compartment programs
from bulk". It is also a candidate L2 instrument.

**Karakaslar / seAMLess 2024** — full methods read closes a flagged unknown **in our favour**: pure
MuSiC composition on 1,350 samples, the word "differential" occurs once in the paper and it is in a
bibliography entry, and the "joint mutation and maturation modelling" is an unimplemented proposal.
L1 is untouched; but it and Zeng 2022 firmly own composition→drug/survival, so L2 must survive
conditioning on composition.

### 5.5 Conference and preprint sweep

ASH 2025 / EHA 2026 / preprints 2025-06→2026-09: **no new occupant in either lane.** One new
single-driver entrant — a Schuringa-lab preprint (2026-06) tying a "RAS mutant-like" program in
RAS-**wild-type** patients to venetoclax resistance. Commercial phenocopy *detection* is heating up
(two Genomic Testing Cooperative ASH posters: KMT2Ar-like for menin-inhibitor selection, IDH1/2-like)
— motivation support for L2, and a signal that trial-selection applications may move fast. Sweep
limitation: `ashpublications.org` blocks scripted fetches, so the ~7,000-abstract supplement was
searched only through indexed pages.

---

## 6. eCHROMA — the epigenetic subgroup layer is half-open

Ochi et al., *Nature* 2026;656:752–762 (1,563 cases, 16 ATAC subgroups; Sweden 1,040 / Kyoto 523).

**Open, verified by download and cell-by-cell parse:**

- **Supplementary Table 3** — all **1,563 samples** with ATAC subgroup, WHO/ICC/ELN class,
  **per-sample driver-gene string**, and per-assay availability flags. Free, no login.
- **Supplementary Table 11** — the ClaNC classifier gene lists: 240 genes (15 per subgroup) plus four
  30-gene models.

**Blocked:**

- **ClaNC is not runnable from open files.** ST11 contains **gene symbols only, zero numeric cells**
  — verified in two independent copies (journal supplement and bioRxiv supplement). ClaNC classifies
  by standardized distance to class centroids, so gene names alone are unusable, and retraining is
  impossible because the labels are open but the expression matrix is not.
- The runnable model is **Zenodo 10.5281/zenodo.17585383 — restricted**. The processed-data record
  carried a July 2026 notice pausing new access requests; its 2026-08-24 version drops that wording
  and asks for name/affiliation/PI/purpose, which suggests requests are being processed again.
  Corresponding author: Yotaro Ochi, Kyoto University.
- External-cohort predictions (n=1,079 across four adult cohorts) are published as **summary figures
  only** — no per-sample table.
- Raw data: EGA EGAS00001008315 (DAC review ~2 weeks) and WGS in the **Japan-hosted G-CARD**
  (G-CARDS000002) — likely the easiest controlled tier for an Osaka-U applicant. **Not** JGA, **not**
  dbGaP.

**Open fallback: AMLmapR** (`github.com/jeppeseverens/AMLmapR`, CC BY-NC-SA) — trained SVMs ship
inside the repo (46.9 MB), runs offline, predicts the **Severens 2024 17 transcriptional clusters**
(not eCHROMA's 16 ATAC subgroups). Requires a full 60,660-gene GENCODE v36 STAR count matrix, so
microarray cohorts need re-quantification. Verified by code review, not execution.

> **Action:** file the Zenodo request and email Ochi now, since the lead time is weeks; use AMLmapR in
> the meantime; and note that ST3 already gives open subgroup labels for the eCHROMA samples
> themselves. eCHROMA stratification is a **strengthening layer, not a gate** — do not let L1 depend
> on it.

---

## 7. What to do next

1. **Pilot L1 on 2–3 drivers** using the existing pseudobulk pipeline: extend P2's design to
   FLT3-ITD and TP53 (or IDH2), with the N6 compartment-adjustment arm as the shadow test and the
   2026-08-04 discovery/validation split respected for the first time. Deliverable: a 3-driver ×
   7-compartment legibility matrix with shadow verdicts.
2. **Stand up L2 as a Snakemake pipeline** — download BeatAML2 + GSE146173 + TCGA-LAML, harmonize,
   score compartment programs, model discordance conditioned on genotype **and composition**. This is
   also the cleanest route to a working pipeline repo.
3. **Read before writing the pitch:** scOPE Supplementary Figs S7/S8 (which AML drivers localize
   where — the closest thing to a scoop of L1), Isobe 2023 in full (the shape competitor), and the
   Schuringa RAS preprint (the drug-response half of L2).
4. **Positioning, fixed in advance:** *programs given composition* against Zeng's composition;
   *compartment-resolved and composition-adjudicated* against scOPE's whole-cell axes;
   *systematic and multi-driver and both-direction* against Lee's single-driver TP53;
   *human and multi-compartment with discordance* against Isobe's mouse HSPC.

## 8. Known gaps in this recon

Honest list of what remains unchecked, so none of it gets mistaken for verified:

- The ASH 2025 supplement was keyword-searched through indexed pages only (site blocks scripted
  fetches); a low-profile abstract using different vocabulary could have been missed.
- Whether Zenodo access requests on the eCHROMA records are actually being **granted** post-August
  2026 is inferred from version history, not confirmed.
- scOPE Supplementary Tables S1/S2 and Figs S7/S8 were not opened.
- GEN-PHEN-VEN (Lachowiez 2025, *Blood Cancer Discov*) — cohort size and endpoints unread (publisher
  blocked); the strongest peer-reviewed 2025 statement that phenotype adds to genotype for
  venetoclax response. Read before finalizing the competition section.
- Severens 2024's full cluster–genotype table is paywalled; the AMLmapR cluster correspondence is
  abstract-level only.
- GEO-vs-SRA raw-access contradictions in GSE158067 and GSE268962 were reproduced but not resolved —
  matters only if raw reads are ever needed.
- DepMap 26Q1's AML count (30) is a lower bound; 16 of its models are absent from the 24Q4 model
  table and the portal blocks scripted enumeration.
