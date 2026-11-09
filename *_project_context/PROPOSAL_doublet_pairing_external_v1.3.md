# Mapping state-dependent T cell–leukemia pairing in acute myeloid leukemia using heterotypic scRNA-seq doublets

## Project summary

Acute myeloid leukemia (AML) is not one cell state but a hierarchy of leukemic states that mirrors, imperfectly, normal myeloid differentiation. Whether T cells in the marrow engage all of these states equally, or whether some states are systematically under-represented among T-cell partners, is not known. Expression-based inference of ligand–receptor signalling cannot answer this, because it reports what cells could signal, not which cells were together; and cell-type abundance reports only how many cells of each kind were present. Droplet single-cell RNA-seq contains a third, routinely discarded, class of information: barcodes that captured two cells. Some of these heterotypic doublets arise from random co-encapsulation, but some capture cells that were physically associated when the suspension was loaded. Treated carefully, the population of T cell–leukemic cell doublets can be tested for structure that random co-encapsulation would not produce.

We hypothesize that, after accounting for random co-encapsulation, ambient RNA, and the detection and typing biases of the doublet-calling pipeline, a subset of leukemic cell states in AML is reproducibly over- or under-represented among T–leukemia candidate pairings, and that this pairing preference is associated with specific expression programs on the leukemic and T-cell sides. The study proceeds in three stages: establish candidate pairings and calibrate their technical background against synthetic doublets and confirmed cross-donor doublets (Aim 1); test, with the patient as the unit of inference, whether the allocation of pairings across leukemic states departs from state availability, and identify associated programs (Aim 2); and, only if the first two aims meet pre-specified entry conditions, test candidate leukemic factors in co-culture with separate readouts for binding, T-cell response, and killing (Aim 3).

The expected contribution is a map of which stem/progenitor and blast-like leukemic states T cells preferentially pair with in patient marrow, with the programs that accompany pairing, derived from public data with bias controls designed for the question and tested for reproducibility in an independent cohort. Outcomes are pre-defined so that a positive result identifies factors for mechanistic follow-up, a bounded negative result states what the data can exclude, and an inconclusive result is reported as such.

---

## 1. Background, knowledge gap and rationale

### 1.1 Leukemic cell states and the immune question

Single-cell transcriptomics has established that AML blasts occupy a range of differentiation-like states, from stem/progenitor-like to promonocyte-like, within one patient and across patients, and that the distribution of these states carries prognostic and therapeutic information (van Galen et al., 2019; Zeng et al., 2022). The immune side has been studied in parallel: T cells in AML marrow show altered composition and exhaustion-like features, and immune-escape mechanisms have been described at relapse. What has not been addressed is the joint question: given the leukemic states available in a patient's marrow, are T cells found with all of them in proportion to their abundance, or do some states appear disproportionately often, or rarely, among T-cell partners? If pairing is state-dependent, the states that are under-paired and the programs that distinguish them are the natural starting point for understanding how leukemic cells limit T-cell engagement; the states that are over-paired, and the T subtypes that pair with them, are where functional T–leukemia interaction is most likely to be examined.

### 1.2 What different kinds of single-cell evidence can and cannot say

Three kinds of evidence are commonly extracted from a single-cell dataset, and they answer different questions.

*Cell-type abundance* says how many cells of each kind were captured. It is the strongest single descriptor of patient-level variation in most single-cell cohorts, but it says nothing about which cells were associated.

*Ligand–receptor inference* scores the co-expression of a ligand in one cell type and its receptor in another. It is a statement about potential: two populations could be signalling. It is sensitive to abundance and to ambient RNA, and it cannot say whether two particular cells were together. For the question posed here, a different class of evidence is needed.

*Heterotypic doublets* are a record of an event in the sample: two cells shared a droplet. For any single doublet this cannot be distinguished from random co-encapsulation. For the population of doublets, however, random co-encapsulation makes a specific prediction, namely that pair types occur in proportion to the product of their constituents' frequencies, modulated by how well the pipeline detects each pair type. Departures from that prediction, if they survive calibration, are information that neither abundance nor ligand–receptor inference contains.

### 1.3 Doublets as evidence: what prior work has established

Several groups have used doublets in this way. Neighbor-seq reconstructs physical interaction networks from standard droplet scRNA-seq by training a classifier on synthetic doublets and comparing observed heterotypic doublets with the random expectation (Ghaddar and De, 2022). CIcADA (Cell type-specific Interaction Analysis using Doublets) applies a related strategy to identify interacting cell-type pairs and their expression features (Schiebout et al., 2023). PIC-seq took the experimental route, sorting physically interacting cells and sequencing them, and showed that paired-cell transcriptomes report contact-dependent T-cell states in immune tissues (Giladi et al., 2020). In AML specifically, García-Guerrero et al. (2018) formed T cell–leukemic cell conjugates from patient samples in culture, sorted them by flow cytometry, and found the conjugated T cells enriched for tumor-reactive cytotoxic T lymphocytes. That study demonstrates that T–leukemia conjugates in AML are biologically informative. It does not establish that such associations survive Ficoll separation, cryopreservation, thawing, and loading in the public datasets we will use. This proposal cannot determine the fate of any particular association; what it can test is whether detectable candidate-pairing structure remains after processing.

To our knowledge, heterotypic doublets have not been analyzed as pairing evidence in AML single-cell data; a literature check will be repeated before publication. The methodological base is therefore established, and the contribution sought here is biological: whether AML shows reproducible state-dependent pairing structure, which states and T subtypes are involved, and which programs accompany it.

### 1.4 Terminology

A T–leukemia doublet that passes calling and origin typing is a **candidate pairing** ("pairing" for short). *Physical contact* and *stable conjugate* are used only for studies that measured contact directly and for the conditional experimental aim. Under-representation of a state among pairings is reported as such, not as immune evasion; programs associated with pairing are associations, not contact-induced changes; cytotoxic expression or TCR clonality in paired T cells does not by itself establish tumor specificity.

---

## 2. Central hypothesis, objectives and innovation

### 2.1 Hypothesis

After controlling for random co-encapsulation, ambient RNA, and detection and typing bias, a subset of leukemic cell states in AML is reproducibly over- or under-represented among T–leukemia candidate pairings, and this pairing preference is associated with specific expression programs on the leukemic side and on the T-cell side.

The hypothesis is direction-free with respect to which states are enriched, and it does not presuppose that pairing reflects cytotoxic engagement. Activated, cytotoxic, regulatory, and dysfunctional T-cell states are all admissible findings; the T-side readout is designed to distinguish them rather than to assume one.

### 2.2 Objectives

**Aim 1 — Establish reliable candidate pairings and their technical background.** Determine which public AML libraries retain doublets, type the doublets by origin with measured accuracy, and calibrate the random-co-encapsulation background, including false-positive assignments, against synthetic doublets and against confirmed cross-donor doublets where multiplexed libraries allow.

**Aim 2 — Identify state-dependent pairing structure and associated programs.** Test, conditional on each patient's total number of T–leukemia pairings, whether their allocation across leukemic states departs systematically from state availability (primary test). Characterize which T subtypes pair with which leukemic states. Identify expression programs by two complementary analyses: deviations of real pairings from matched synthetic pairings, and programs associated with pairing propensity at the patient-by-state level.

**Aim 3 (conditional) — Validate candidate factors.** If Aim 2 yields a statistically supported, replicated, above-threshold result and a candidate factor that replicates in held-out patients, perturb that factor in AML cell lines co-cultured with T cells and measure conjugate formation, T-cell response, and killing as separate readouts.

Secondary and exploratory analyses are kept distinct from the aims: overall enrichment of T–leukemia pairings over background (secondary; not required for the primary result), the cytotoxic/activation state of paired T cells (exploratory readout), patient-level heterogeneity in pairing, and association of a paired-T program with outcome in bulk cohorts (exploratory). In vivo experiments are contemplated only if Aim 3 supports a candidate and are not committed here.

### 2.3 Where the contribution lies

The methods used — doublet calling, origin typing, synthetic-doublet comparators, non-negative mixture decomposition — are established tools and are not claimed as novel. The contribution is (i) the biological finding, if the data support it: a reproducible, state-resolved map of T–leukemia pairing in patient marrow with its associated programs; (ii) the bias-control design built for this question, in particular a null model that propagates state availability through the pipeline's full confusion structure and a decision framework that distinguishes supported, bounded-negative, and inconclusive outcomes; and (iii) a pre-specified path from a computational association to a functional test that reads out binding, response, and killing separately.

---

## 3. Preliminary foundation and data resources

### 3.1 Assets already in hand

*Curated cohort inventory.* Thirteen public AML single-cell datasets have been curated into a single metadata schema (72 fields per sample) covering donor, timepoint, tissue, sorting, processing chain, and platform. This is what makes library-level eligibility decisions possible.

*Compartment audit against genotype.* Using datasets with per-cell somatic mutation calls (Petti et al., 2019) and targeted genotyping (van Galen et al., 2019), we determined which annotated compartments are clonal (preliminary, unpublished). In the genotype-covered T-cell reference (395 cells), none carried the queried AML mutations, which supports T cells as a reference compartment within the limits of that check. The hematopoietic stem/progenitor compartment in AML samples is predominantly clonal (92–99% in genotyped cells). Monocyte-like and dendritic-cell populations are largely clonal in adult AML, and one integrated cluster annotated as T/NK was found to consist of cycling myeloid progenitors. These findings define which cells may serve as the non-malignant partner and which may not, and they motivate defining leukemic cells by compartment rather than by copy-number-based malignancy calls, which fail in normal-karyotype AML.

*Recoverable pre-removal barcodes.* The project's ingest pipeline stores the deposited matrices before any project-side quality control, including, for multiplexed libraries, barcodes that carry counts but no author donor assignment. Doublets removed by project QC are therefore recoverable from these objects.

*Multimodal cross-checks.* Per-cell mutation calls (Petti et al., 2019) and genotype-based demultiplexing in multiplexed libraries provide partial, independent checks on doublet origin that pure expression typing lacks.

### 3.2 What remains to be audited

Whether a given library retains doublets depends on whether the depositing authors removed them before deposition. This must be verified per library from the processing chain, the deposited barcode set, and the count distribution, not assumed. Reliable-pairing yield per patient and per leukemic state, and therefore the achievable precision, will be estimated by simulation before any real pairing counts are examined. The number of patients meeting eligibility is not yet known; the sample sizes below are nominal.

### 3.3 Cohorts under consideration

| Dataset | Nominal content | Platform | Proposed role | Eligibility status |
|---|---|---|---|---|
| GSE185381 (Lasry et al., 2023) | 42 AML patients (20 adult, 22 pediatric), 10 healthy donors; multiplexed libraries demultiplexed by genotype | 10x 3′ | Proposed primary cohort | To be audited per library; undemultiplexed barcodes present and recoverable |
| GSE289435 (Zeng et al., 2025) | Adult AML, diagnostic; same study as the reference map | 10x 5′ | Candidate replication | To be audited |
| Petti et al., 2019 | Adult AML with per-cell mutation calls | 10x 5′ | Candidate replication; genotype cross-check | To be audited |
| GSE227903 (Ennis et al., 2023) | Diagnosis, residual disease, relapse | 10x 3′ | Excluded | Authors removed doublets |
| van Galen et al., 2019 | Diagnosis and post-treatment | Seq-Well | Excluded | Non-droplet doublet mechanism |
| GSE185991 (Naldini et al., 2023); Chen et al., 2023 | Sorted progenitor fractions | 10x | Excluded | No T cells |
| GSE239721 (Wang et al., 2024) | Adult AML | 10x 5′ | Excluded | Essentially no T/NK cells captured |

Final assignment of primary and replication roles follows the eligibility audit and the simulated yield; a library may be partially eligible. No post-treatment comparison is planned, because the post-treatment datasets available are either author-doublet-removed or sorted without T cells.

### 3.4 Cell compartments and states

Non-malignant partner: **strict T cells**, annotated as T and confirmed by T-identity transcripts, excluding the cycling-progenitor cluster and cells with NK-only markers. NK cells and unresolved T/NK cells are analyzed separately and do not enter the primary test. Leukemic partner (**LEUK**): the stem/progenitor and blast-like compartment of AML samples, defined by projection onto a reference bone-marrow map. B cells and erythroid cells serve as reference pair partners (B with the caveat that it is partly clonal in some genotypes); monocyte-like and dendritic cells are excluded from tests because they are predominantly clonal, and are reported separately.

Leukemic **states** are differentiation bins on the projection axis (HSC-like, MPP-like, LMPP-like, GMP-like, promonocyte-like). **Cycling** is an independent axis, assigned to every leukemic cell as a flag and used for stratified sensitivity analysis, not as a state in the differentiation series. T subtypes (naive, central/effector memory, cytotoxic effector, regulatory, γδ) are projection labels.

Clonal evidence applies at the compartment level; it does not mean that every individual cell carries a confirmed malignant genotype.

---

## 4. Research strategy

### 4.1 Aim 1 — Reliable candidate pairings and their technical background

*Question.* Do eligible libraries contain T–leukemia doublets whose origin can be typed with measured accuracy, and what is the frequency of such doublets expected from random co-encapsulation and from mis-assignment?

*Eligibility from pre-removal objects.* For each library we return to the deposited matrix and classify barcodes into three classes: A, author-assigned and retained by project QC; B, author-assigned but removed by project QC or doublet calling; C, present with counts but without an author donor assignment (in multiplexed libraries). Barcodes absent from the deposited matrix are unrecoverable. Classes A and B are available for analysis; class C is retained and judged per library, depending on whether the donor can be confirmed downstream by SNP or hashtag information. Class C is neither discarded wholesale nor used as a cross-donor control. A library is eligible when the deposited matrix demonstrably predates author doublet removal, judged from the processing chain, the barcode set relative to the Cell Ranger filtered set, and the upper tail of the count distribution.

*Typing with measured accuracy.* Doublets are called and their origin typed with scDblFinder at a version in which the documented origin-assignment error is fixed, using cluster labels that include the strict-T compartment, the other compartments, and the leukemic states, so that typing resolves leukemic state as well as lineage. Accuracy is measured on synthetic doublets built from each library's singlets, passed through the same calling and typing, and evaluated on synthetic doublets from held-out patients, so that the typing model is not scored on the cells that trained it. The result is, per library, a confusion matrix giving the probability that a true pair of each type is retained and reported under each label, and a separate false-positive background: the rate at which singlets and other pair types are reported as T–leukemia. The genes used by typing are recorded and excluded from all later functional gene sets.

*Background calibration.* In multiplexed libraries, doublets confirmed by hashtag or SNP to contain cells from two donors serve as a technical-background calibration set: the two cells did not share a marrow, although adhesion after pooling cannot be excluded. Their allocation across states and their overall frequency check that the modelled background is calibrated; deviation is reported as a bound on residual technical background, not corrected away. Where genotype information exists, the strict-T singlet reference is checked for the absence of mutant reads, and candidate pairings for mutant evidence compatible with a leukemic member. Reads in a mixed barcode cannot be assigned to one member, and absence of a detected mutation is recorded as unknown, not wild-type.

*What Aim 1 delivers.* A per-library eligibility and barcode-class table; confusion matrices and false-positive backgrounds; a background-calibration report; and, from simulation at the expected number of pairings per patient, the number of reliable pairings per state and the precision attainable. These quantities set the state resolution for Aim 2 and are frozen before any real pairing count is examined.

*Risks.* Few libraries may retain doublets; typing may be inaccurate at fine state resolution; sample processing may have dissociated associations. Each is addressed by a pre-specified response (§6).

### 4.2 Aim 2 — State-dependent pairing and associated programs

*The primary question, stated precisely.* Relative to the leukemic-state composition available in the same patient and library, and to the pipeline's state-specific detection and typing behaviour, does the allocation of T–leukemia pairings across leukemic states show a systematic departure across patients?

*Why conditioning matters.* The test is conditioned on each patient's total number of T–leukemia pairings. A patient with many pairings does not by that fact contribute evidence; what matters is how those pairings are distributed across the states that patient has. This separates the primary question from overall enrichment, which is a secondary test. A state can be relatively under-represented among pairings without any claim that its absolute pairing probability is below random.

*The null and the alternative share the same bias model.* The expected allocation under the null is the patient's state availability pushed through the library's confusion matrix and false-positive background. The alternative multiplies true-state availability by a shared vector of pairing propensities before the same confusion structure is applied. Because both models pass through the same measurement process, propensities are interpretable on the true-state scale, and effects are never read off observed-label counts directly. A state absent from a patient's singlets contributes no true pairings in that patient but can still receive mistyped counts; its observed category is retained and is never recorded as "unpaired".

*Inference.* The model is a Dirichlet–multinomial over observed-label counts per patient, with the patient as the unit and overdispersion absorbing patient-level variation, tested by likelihood ratio against the null, direction-free. Because counts are sparse and the confusion matrices are estimated, the reference distribution of the statistic and all confidence regions are obtained by parametric simulation that regenerates doublets from each patient's singlets through the full pipeline, preserving patient and library structure. Effect sizes are reported on a conditional observed/expected scale, R, defined at a frozen reference composition; R is not the exponentiated model coefficient, and the follow-up threshold of R ≥ 1.5 or R ≤ 2/3 is defined on R. Per-state estimates with simultaneous intervals are reported for all states regardless of the global result; statements naming specific states are made only after a positive global test, with multiplicity control. Appendix A gives the model, the reference composition, and the region construction.

*Outcomes.* Four pre-defined outcomes (Appendix A.5) separate statistical support from follow-up magnitude: supported above the follow-up threshold (A), supported with magnitude not established (B), not supported at adequate precision (C), or insufficient (D). Eligibility, typing accuracy and identifiability gates are applied before any of these. Replication in an independent cohort is recorded separately.

*Pairing by T subtype.* Each T–leukemia doublet is decomposed into one T-subtype profile and one leukemic-state profile using patient-specific references (non-negative least squares; an optimal-transport variant as a sensitivity analysis), yielding the most likely pair, the mixing fraction, and a fit margin. Decomposition accuracy is calibrated on held-out synthetic doublets, and ambiguous assignments are dropped. This gives the state-by-subtype pairing map.

*Two program analyses, two questions.* Two synthetic control sets are kept distinct. The random null for the primary test draws pairs according to singlet state availability and never matches the state proportions of real pairings, which would copy the effect under test into the null. The matched reference for program analysis draws synthetic pairs matched to each real pairing on patient, T subtype, leukemic state, sequencing depth, and estimated mixing ratio.

Analysis A asks what is present in real pairings beyond what mixing two singlets predicts: real pairings are compared with their matched synthetic pairings, within patient and state pair, on genes pre-selected on singlets and never used in typing. This is the only analysis that interrogates the doublet transcriptome itself. The T-side readout (activation and cytotoxic genes expressed in T singlets and silent in leukemic singlets, scored after normalization by the estimated T mixing fraction) is one instance of Analysis A and is exploratory with respect to the primary question.

Analysis B asks which singlet programs track pairing propensity: at the patient-by-state level, pairing propensity is related to pseudobulk program scores of that patient's singlets in that state, within similar differentiation backgrounds (same or adjacent states, with pseudotime as a covariate). Two estimation routes are pre-specified — shrunken patient-by-state deviations with propagated uncertainty, or entering the program score as an explanatory variable in the count model — and one will be selected and frozen before unblinding. Singlets are never labelled "paired" or "under-paired", because the pairing map does not identify which singlet participated. Results are reported as programs associated with pairing propensity.

Both analyses require replication in held-out patients and, where an eligible cohort exists, in an independent cohort. A candidate list of leukemic factors — surface, adhesion, and antigen-presentation genes with available reagents, supported by both analyses where possible — is frozen before Aim 3.

*Exploratory extensions.* Patient-level anomalies in pairing allocation are listed as exploratory and do not enter the primary decision. A paired-T program drawn from Analysis A, restricted to genes absent from leukemic singlets, is scored in BeatAML and TCGA-LAML bulk RNA-seq with Cox models adjusted for ELN risk, age and blast percentage; this is reported as association, with the explicit caveat that it cannot validate pairing and may reflect T-cell abundance or cytotoxic tone.

### 4.3 Aim 3 — Conditional validation of candidate factors

*Entry conditions.* Outcome A in the primary cohort; replication status "replicated" in at least one independent cohort with adequate yield; and a candidate factor that replicates in held-out patients.

*Design.* The laboratory maintains human AML cell lines spanning CD34-positive, promyelocytic-like, myelomonocytic and monocytic states (e.g., TF-1, HL-60, OCI-AML3, MOLM-13); two to three matching the states implicated by Aim 2 will be selected once that result is known. Imaging flow cytometry is available. T cells are not currently held in the laboratory and are a planned condition: commercially sourced healthy-donor T cells or a collaboration, under applicable ethics approval. Either source is allogeneic and its reactivity is not tumor-specific; antigen-matched T cells will be used only if secured. Candidate factors are perturbed by CRISPR knockout or lentiviral overexpression.

*Three separate readouts.* Conjugate formation is measured by flow-cytometric doublet gating and imaging flow cytometry, with LFA-1/ICAM-1 blockade as the reduction control; T-cell response by CD69, IFN-γ, and degranulation; killing by target viability. Each candidate is scored by which step it affects, because a factor may change binding without changing killing, or vice versa. Conjugated T cells are sorted and assayed for the paired-T program.

*Criteria.* Support requires an effect in the direction implied by Aim 2, exceeding blockade-normalized noise, in at least three biological replicates in at least two lines. In vivo work is contemplated only if an in vitro effect is clear and is not committed here.

---

## 5. Expected outcomes and significance

If pairing is state-dependent in the primary cohort (outcome A or B), the study delivers a map of which stem/progenitor and blast-like leukemic states T cells preferentially pair with in patient marrow, resolved by T subtype, with the programs that accompany pairing on both sides; the map is described as reproducible only when an independent cohort agrees. Outcome A with replication yields candidate factors for Aim 3. If the data exclude pairing differences of the pre-set size (C), the result is a bounded negative restricted to these datasets, this state resolution and this model. If precision is inadequate (D), the study reports what was attainable and what would be needed.

How a positive result would be read depends on its pattern; the following are conditional readings, not predictions. If under-paired states differ from their adjacent, better-paired states mainly in adhesion or antigen-presentation programs on the leukemic side, the candidates for Aim 3 are those surface factors, and the first experimental question is whether restoring them changes conjugate formation rather than killing. If the difference lies on the T-cell side — for instance, pairings with primitive states carrying regulatory or dysfunctional subtypes while pairings with committed states carry cytotoxic subtypes — the finding concerns T-cell state rather than leukemic surface biology, and functional follow-up would require antigen-matched T cells before any causal claim. If pairing tracks cycling status more than differentiation state, the result concerns proliferation-linked programs and is reported as such. Each reading is testable within the design; none is assumed.

The direct deliverables are the eligibility audit of public AML libraries for retained doublets, the calibrated null model, the pairing map, the program analyses, and, conditionally, functional evidence on candidate factors. Downstream uses such as TCR discovery or biomarker development depend on results the design cannot anticipate and are not claimed.

---

## 6. Feasibility, risks and alternatives

| Risk | Control or response within the design |
|---|---|
| Few libraries retain doublets | Audit from pre-removal objects decides; a single eligible cohort still yields an A–D outcome with replication recorded as unavailable |
| Processing dissociated associations | Yields outcome C or D depending on attainable precision; conclusion stated within scope |
| Typing inaccurate at fine state resolution | Pre-registered merging along the differentiation axis with frozen trigger, order and stopping rule; fewer than two comparable states yields D |
| Confusion matrices biased | Held-out evaluation; cross-donor calibration bounds residual background |
| Result driven by one patient | Patient-level overdispersion in the model; anomalies reported separately |
| Programs reflect differentiation rather than pairing | Matched synthetic pairs of the same states; similar-background comparison with pseudotime covariate; held-out replication; Aim 3 tests causality |
| Too few pairings per state | Simulation audit before unblinding; outcome D if precision inadequate |

Alternatives within the design are merging states, declaring insufficiency, and waiting for additional eligible cohorts.

---

## 7. Milestones and timeline

| Milestone | Target |
|---|---|
| Eligibility audit of libraries and barcode classes completed; design and thresholds registered | November 2026 |
| Typing accuracy, confusion matrices, background calibration and power audit completed; state resolution and numeric parameters frozen | December 2026 |
| Primary test of state-dependent pairing unblinded; outcome and replication status recorded | December 2026 |
| Pairing map by T subtype; program analyses; candidate factor list frozen (if outcome A or B) | January–February 2027 |
| Exploratory bulk association; Phase 2 report; experimental design | February 2027 |
| Functional validation, if entry conditions are met | March–September 2027 |
| Manuscript and thesis chapter | late 2027 – early 2028 |

Experimental entry is conditional and is not assumed in the timeline.

---

## References

- García-Guerrero E, et al. Selection of tumor-specific cytotoxic T lymphocytes in acute myeloid leukemia patients through the identification of T-cells capable to establish stable interactions with the leukemic cells ("doublet technology"). *Front Immunol.* 2018;9:1971. doi:10.3389/fimmu.2018.01971
- Ghaddar B, De S. Reconstructing physical cell interaction networks from single-cell data using Neighbor-seq. *Nucleic Acids Res.* 2022;50(14):e82.
- Giladi A, et al. Dissecting cellular crosstalk by sequencing physically interacting cells. *Nat Biotechnol.* 2020;38:629–637. doi:10.1038/s41587-020-0442-2
- Schiebout C, et al. Cell type-specific interaction analysis using doublets in scRNA-seq (CIcADA). *Bioinform Adv.* 2023;3(1):vbad120.
- Lasry A, et al. An inflammatory state remodels the immune microenvironment and improves risk stratification in acute myeloid leukemia. *Nat Cancer.* 2023. (GSE185381)
- Petti AA, et al. A general approach for detecting expressed mutations in AML cells using single cell RNA-sequencing. *Nat Commun.* 2019;10:3660.
- van Galen P, et al. Single-cell RNA-seq reveals AML hierarchies relevant to disease progression and immunity. *Cell.* 2019;176:1265–1281.
- Zeng AGX, et al. A cellular hierarchy framework for understanding heterogeneity and predicting drug response in acute myeloid leukemia. *Nat Med.* 2022;28:1212–1223.
- Zeng AGX, et al. Single-cell transcriptional atlas of human hematopoiesis reveals genetic and hierarchy-based determinants of aberrant AML differentiation. *Blood Cancer Discov.* 2025;6(4):307–324. (reference map; GSE289435)
- Ennis S, et al. Cell–cell interactome of the hematopoietic niche and its changes in acute myeloid leukemia. *iScience.* 2023. (GSE227903)
- Naldini MM, et al. Longitudinal single-cell profiling of chemotherapy response in acute myeloid leukemia. *Nat Commun.* 2023;14:1285. (GSE185991)
- Wang B, et al. Comprehensive characterization of IFNγ signaling in acute myeloid leukemia reveals prognostic and therapeutic strategies. *Nat Commun.* 2024;15:1821. doi:10.1038/s41467-024-45916-6 (GSE239721)
- Chen Y, et al. A single-cell taxonomy predicts inflammatory niche remodeling to drive tissue failure and outcome in human AML. *Blood Cancer Discov.* 2023;4(5):394–417. PMID 37470778. Data: Mendeley Data, doi:10.17632/gwjh3w6ztm.2
- Germain PL, Lun A, Garcia Meixide C, Macnair W, Robinson MD. Doublet identification in single-cell sequencing data using scDblFinder. *F1000Res.* 2021;10:979.

---

# Technical appendix

## A. Statistical design

### A.1 Notation
For patient *i*: *M_i* reliable T–LEUK candidate pairings; true leukemic states *k* = 1…K at the frozen resolution; *f_i* the vector of state shares among available LEUK singlets; *C_i* the K×K confusion matrix for that patient's library, with C_i[j,k] the probability that a true T–LEUK(k) doublet is retained by calling and filtering **and** reported as label *j* (column sums may be < 1); ρ_i and q_i^FP the false-positive fraction and label distribution among barcodes reported as T–LEUK (Aim 1); *β* a K-vector of log pairing propensities with a sum-to-zero constraint.

### A.2 Generative model
Observable allocation is a proportion mixture of a signal component and a false-positive component:

q_i(β) = (1 − ρ_i) · q_i^signal(β) + ρ_i · q_i^FP,  with q_i^signal(β) ∝ C_i [ f_i ⊙ exp(β) ],

where ρ_i is the estimated fraction of barcodes reported as T–LEUK in patient *i* that are false positives (singlets or other pair types), and q_i^FP is the normalized label distribution of those false positives, both estimated from synthetic data in Aim 1. Uncertainty in ρ_i and q_i^FP is propagated through the parametric simulation. Null: β = 0. Alternative: shared β with Dirichlet–multinomial overdispersion. A true state with f_ik = 0 contributes no true pairings in patient *i*; its observed label is retained because it can receive mistyped counts.

### A.3 Library-to-patient aggregation
Expectations are computed forward per library (C and f at library level) and aggregated to the patient with structure preserved. C and f are never averaged across libraries before the forward computation.

### A.4 Conditional O/E scale and reference composition
R_k = [exp(β_k)/Σ_j f̄_j exp(β_j)] / [1/Σ_j f̄_j], evaluated at one pre-specified primary reference composition f̄ (equal weights across qualifying patients in the primary cohort), frozen before unblinding and used unchanged for replication cohorts. Age-stratified references are a sensitivity analysis only; they do not enter the A–D decision. R_k ≠ exp(β_k): with f_k = 0.5, multiplying the state's weight by 1.5 yields a share of 0.6 and R_k = 1.2. Negligible band: R_k ∈ [2/3, 1.5] for all k. Follow-up threshold: R_k ≥ 1.5 or ≤ 2/3 for at least one k. Power simulations generate the target R_k directly by solving for β at f̄ and pushing through C_i.

### A.5 Test, regions and outcomes
Likelihood-ratio test of β = 0, direction-free, with reference distribution from ≥ 2,000 parametric simulations regenerating doublets from singlets through the full pipeline. A simultaneous confidence region for (R_1…R_K) is constructed by one pre-specified method (likelihood-based or simulation-calibrated) whose coverage is verified on the parametric simulation before use.

Decision priority: eligibility, typing accuracy (per-state recall and precision ≥ 0.7 on held-out synthetic doublets), state identifiability (≥ 2 comparable states with sufficient within-patient overlap), and the operational minimum (≥ 10 qualifying patients with ≥ 30 strict-T singlets, ≥ 50 LEUK singlets, ≥ 5 typed T-containing doublets) are checked first; failure yields D. Then:

- **A.** Global test rejects at α = 0.05; after localization with multiplicity control, at least one state's simultaneous interval lies entirely outside the negligible band.
- **B.** Global test rejects; no state's interval lies entirely outside the band. B-small if the simultaneous region lies entirely within the band; B-undetermined otherwise.
- **C.** Global test does not reject **and** the simultaneous region lies entirely within the band.
- **D.** Any gate failure; or the global test does not reject and the region is not contained in the band.

Per-state R_k and simultaneous intervals are reported in all outcomes.

### A.6 Replication and experimental entry
Replication status, recorded separately: replicated (direction of localized effects agrees in ≥ 1 independent cohort with adequate yield), not replicated, or unavailable. Outcome A with replication unavailable is "discovery evidence established; cross-cohort reproducibility not yet established". Entry to Aim 3 requires A, replicated, and a candidate factor replicating in held-out patients.

### A.7 Secondary test
Overall T–LEUK observed/expected per patient, one row per patient; summarized by patient-bootstrap interval and a permutation test against the simulated null; T–B and T–Ery as reference pair types. Not required for the primary result.

## B. State resolution and merging
Differentiation states: HSC-like, MPP-like, LMPP-like, GMP-like, promonocyte-like. Cycling is a separate flag. Merging trigger: a state fails the typing threshold or has fewer than the frozen minimum of reliable pairings in more than half of qualifying patients. Order: HSC+MPP → +LMPP → +GMP; promonocyte-like merges last. Stop when all remaining states meet both thresholds; the common level applies to all cohorts. Fewer than two states, or insufficient within-patient overlap, yields D. The minimum pairings per state and the overlap threshold are set by the simulation audit and frozen before unblinding.

## C. Synthetic control sets
Random null: pairs drawn by singlet state availability through detection, typing and filtering; never matched to real pairing proportions. Matched reference: synthetic pairs matched to each real pairing on patient, T subtype, leukemic state (same or adjacent), depth, and estimated mixing ratio. Typing genes excluded from all functional gene sets; typing evaluated on held-out patients.

## D. Program analyses
Decomposition: NNLS over patient-specific T-subtype and leukemic-state references (entropic OT as sensitivity), returning best pair, mixing fraction and margin; calibrated on held-out synthetic doublets. T-side readout: genes with mean CPM < 1 in LEUK singlets and > 5 in strict-T singlets per patient, not used in typing; score normalized by estimated T mixing fraction; minimum 20 T-identity UMIs; pseudo-count 1; 200 expression-matched random gene sets as null. Analysis B estimation routes (one to be frozen): (1) patient × state random deviations with shrinkage and propagated uncertainty, related to pseudobulk program scores within similar backgrounds; (2) program score entered as a covariate on β in the count model, q_i ∝ C_i[f_i ⊙ exp(β + γ·score_i)], testing γ.

## E. Pre-registration plan
Design, model, thresholds and gates are recorded in a version-controlled, timestamped document before eligibility results are inspected (first freeze). Numeric parameters that depend on the simulation audit (minimum pairings per state, overlap threshold, final state resolution, region construction, Analysis B route) are recorded in a second freeze after the audit and before any real pairing count is examined. Both freezes are preserved with the results.

---

# Concept figure: design note and legend

**Layout (left to right, five panels).**

1. *Question.* A schematic AML marrow with leukemic cells drawn along a differentiation gradient (five shades) and T cells of several subtypes. A question mark over the T–leukemia interface: "Do T cells pair with all leukemic states in proportion to their abundance?"
2. *Data opportunity.* Two droplet icons: a singlet and a heterotypic doublet. Beneath the doublet, two latent origins (random co-encapsulation; pre-existing association) drawn as dashed arrows converging on the same observed barcode, with the caption: individual doublets cannot be assigned to an origin; the population of candidate pairings can be tested against the random expectation.
3. *Calibration.* Two parallel arrows feeding one box labelled "same measurement model": (top) state availability → confusion matrix + false-positive background → expected allocation under the null; (bottom) state availability × pairing propensity → same confusion matrix → alternative. Side inputs: synthetic doublets from singlets; confirmed cross-donor doublets.
4. *Analysis.* A stacked-bar schematic of observed versus expected allocation across states (no values drawn), with "R" on an axis and the negligible band shaded; beside it, two small boxes: "Analysis A: real pairing vs matched synthetic pairing" and "Analysis B: pairing propensity vs singlet programs".
5. *Conditional validation.* A gated arrow (labelled "entry conditions: outcome A, replicated, candidate replicates") to co-culture icons with three readout labels: binding, response, killing.

A footer row shows the four outcomes A–D as labelled boxes, with B split into B-small and B-undetermined.

**Legend (English).** *Study design for recovering state-dependent T cell–leukemia pairing from heterotypic doublets.* (1) The question: whether T cells in AML marrow pair with leukemic cell states in proportion to their availability. (2) Droplet single-cell RNA-seq contains heterotypic doublets whose origin — random co-encapsulation or pre-existing association — cannot be assigned per barcode; all T–leukemia doublets enter one candidate pool, and inference is made at the population level. (3) Null and alternative share one measurement model: state availability is propagated through the library-specific confusion matrix and false-positive background, calibrated with synthetic doublets and confirmed cross-donor doublets. (4) The primary test asks, conditional on each patient's pairing total, whether allocation across leukemic states departs from availability; effects are reported on a conditional observed/expected scale with simultaneous intervals. Two program analyses compare real pairings with matched synthetic pairings and relate pairing propensity to singlet programs. (5) Functional validation of candidate factors proceeds only when pre-specified entry conditions are met, with separate readouts for conjugate formation, T-cell response, and killing. The figure depicts design, not results.
