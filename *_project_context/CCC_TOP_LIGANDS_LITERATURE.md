# Top CCC ligands in AML — literature grounding, druggability, and reliability plan

*Prepared 2026-07-27, in response to the lab-meeting requests: (1) references for Galectin-9 / CD99 / MIF in AML — who sends, who receives, what mechanism; (2) identify druggable molecules for clinical relevance; (3) how we make the results reliable enough to justify wet-lab validation.*

> All references below were retrieved and their citation metadata verified against PubMed/PMC/publisher pages. Honesty flags (novel-vs-established, clinical status, unverified author lists) are kept inline — please re-check flagged items against source PDFs before manuscript submission.

---

## 0. Why these three, and how they connect to our own data

These are the three most frequent ligands in our per-sample AML CCC tensors (LGALS9 > MIF > CD99 by significant-edge count). Our co-option probe (`results/tables/05_ccc/cooption/`, see `CANDIDATE_M10_cooption_lens.md`) found a clear split in *who sends them*:

| Ligand | primitive/LSC sender share, AML vs healthy | verdict from our data |
|---|---|---|
| **Galectin-9** | 0.47 vs 0.19 (pooled p=8.6e-4; within-study p=0.07, trend) | shifts to LSC sender in AML |
| **CD99** | 0.62 vs 0.18 (pooled p=2.9e-5; within-study p=0.14) | shifts to LSC sender in AML |
| **MIF** | 0.48 vs 0.47 (n.s.) | broad/unchanged — negative control |

**The literature independently corroborates exactly this split** (Sections 1–3): Galectin-9 and CD99 are established LSC-associated signals with autocrine/LSC-maintenance roles, whereas MIF is a broadly produced niche factor (blasts + stroma + macrophages), not LSC-restricted. So our unbiased computation both **recovers known biology** (a built-in positive control) and **extends it** (the explicit sender-compartment shift is, as far as the reviewers found, not yet published as a cell-cell-communication result).

---

## 1. Galectin-9 (gene *LGALS9*)

- **Sender (AML):** AML blasts and the LSC compartment (CD34+CD38−TIM-3+) secrete Galectin-9; serum Gal-9 is elevated in AML patients and in AML xenografts, and per-cell output far exceeds healthy leukocytes (~5,980 pg/10⁶ cells). [Kikushige 2015; Yasinska 2020; Zhang 2021]
- **Receptor & receiver:** principal receptor **TIM-3/HAVCR2**, present both on the **LSCs themselves (autocrine loop)** and on **T/NK cells**; **VISTA** on T cells is a second AML-verified receptor. PD-1 / CD44 / CD45(PTPRC) are documented Gal-9 glycoreceptors in general immunology but lack strong AML-specific functional evidence — cite as candidate. [Kikushige 2015; Gonçalves Silva 2017; Yasinska 2020]
- **Mechanism:** autocrine Gal-9→TIM-3 co-activates **NF-κB and β-catenin → LSC self-renewal** (present in LSCs, not normal HSCs); immune evasion via blocked NK granzyme-B delivery and T-cell suppression/apoptosis. Nuance: *exogenous/high-dose* Gal-9 is directly cytotoxic to AML incl. CD34+ AML stem cells (autophagy-inhibition, non-apoptotic) — a "double-edged" biology worth stating. [Kikushige 2010/2015; Gonçalves Silva 2017; Choukrani 2023]
- **Prognosis:** high LGALS9/Gal-9 → worse OS in AML (esp. at post-HSCT relapse); plasma Gal-9 independent adverse marker in MDS. [Zhang 2021; Asayama 2017]

---

## 2. CD99 (*MIC2*)

- **Sender/expression (AML):** CD99 is aberrantly overexpressed on AML blasts and, critically, on **LSCs (~5× vs normal HSCs)**; high surface CD99 enriches functional LSC activity and marks an LMPP-like (CD34+CD38−CD90−CD45RA+) leukemic compartment, enriched in *FLT3*-ITD. Consistent with our "primitive/LSC compartment carries more CD99 in AML." [Chung 2017; Vaikari 2020; Travaglini/Ottone 2022]
- **Interaction partner / receiver:** AML-validated axis is **homophilic CD99–CD99** (on neighboring blasts/LSCs and marrow endothelium). Heterophilic partners **PILR** and **CD99L2** are established in normal leukocyte–endothelial trafficking, not shown in AML — candidate only. [Schenkel 2002; Goswami 2017]
- **Mechanism:** CD99 sustains the differentiation-arrested LSC state and adhesion/retention. Agonistic **anti-CD99 mAb directly kills AML/MDS cells incl. LSCs via SRC-family kinase (effector/complement-independent), forces differentiation, downregulates MDM2 → ROS/apoptosis, and spares normal CD34+ cells.** [Chung 2017; Vaikari 2020; Travaglini 2022]
- **Prognosis (nuanced):** as an LSC surface marker, CD99-high = adverse LSC-enriched compartment; but *bulk transcript* high CD99 correlated with **better** OS (ECOG-1900) — cite with the isoform/compartment nuance, not as "high = worse." [Vaikari 2020]

---

## 3. MIF (Macrophage Migration Inhibitory Factor)

- **Sender (AML):** broadly produced — AML **blasts** and the **stromal/macrophage niche** both secrete MIF; single-cell/CellChat work finds it ubiquitously expressed across samples. **No paper claims MIF is LSC-restricted** — this matches our result that MIF's sender composition barely shifts (our negative control). [Abdul-Aziz 2017; Spertini 2024; Anderson-Crannage 2025]
- **Receptors & receiver:** canonical **CD74** with **CD44** co-receptor (reliably high on both AML blasts — autocrine — and BM macrophages); CXCR2/CXCR4 heterophilic but heterogeneous in AML; blast-derived MIF also signals to **BM-MSCs via CD74**. [Leng 2003; Shi 2006; Bernhagen 2007; Abdul-Aziz 2017; Spertini 2024]
- **Mechanism:** CD74/CD44 → ERK1/2, Src, PI3K/AKT → anti-apoptosis/pro-survival; blast→MSC MIF drives **PKCβ→IL-8** feedback that protects blasts; sustains M2-like macrophages (immunosuppression). [Abdul-Aziz 2017; Spertini 2024; Pantouris 2025]
- **Prognosis:** high **CD74** (the MIF receptor) → inferior outcome in CN-AML and pediatric AML. Circulating MIF as an AML prognostic marker is **undocumented**. [Attar 2009; Menssen 2024]

---

## 4. Druggability summary (professor's request #2)

| Ligand / axis | Drug target | Agent(s) | Most advanced stage | Honest status note |
|---|---|---|---|---|
| **Galectin-9 / TIM-3** | TIM-3 (receptor) | **Sabatolimab (MBG453)** anti-TIM-3 | Phase 2/3 in MDS; Phase 2 AML (STIMULUS-AML1, +aza+ven) | **Failed pivotal MDS endpoints (STIMULUS-MDS1/-MDS2); program incl. AML discontinued (2024, press-level).** Targets the receptor, not Gal-9. |
| Galectin-9 (direct) | Gal-9 (ligand) | anti-Gal-9; recombinant Gal-9 (cytotoxic) | Preclinical | No anti-Gal-9 agent in AML trials verified. |
| **CD99** | CD99 (homophilic) | anti-CD99 mAb; FLT3/CD99 bispecific NP; CD99 CAR-T; CD99 ADC (NV101/NV103) | **Preclinical in AML** (clinical in T-ALL context) | mAb kills LSCs sparing normal CD34+; ADC "planned Phase I" is press-level/unverified. |
| **MIF** | MIF (tautomerase/allosteric) | ISO-1, 4-IPP, allosteric "cmpd 1" | Preclinical in AML | Research tool compounds. |
| MIF → CD74 | CD74 (receptor) | Ibudilast (MN-166); Milatuzumab (anti-CD74); STRO-001 (anti-CD74 ADC) | Ibudilast Phase 1b/2a GBM; Milatuzumab Phase 1 B-cell; STRO-001 preclinical peds-AML | None in AML clinical use; NCT03918655 is an observational MIF-in-AML cohort. |

**Bottom line for clinical relevance:** all three are druggable, and two (Galectin-9/TIM-3, CD99) have **LSC-directed** agents — which is exactly the compartment our data implicates. The most clinically mature axis (TIM-3/sabatolimab) has **negative pivotal data in MDS**, so if we build a therapeutic argument, CD99 (LSC-selective mAb, spares normal HSC) and the MIF→CD74 niche loop are the cleaner near-term stories, with Galectin-9/TIM-3 framed as biologically validated but clinically unproven.

---

## 5. Reliability plan (professor's request #3) — what we do before betting experiments

A candidate ligand→niche interaction should clear these gates before it justifies wet-lab work. Current status in brackets.

1. **Multi-method CCC agreement.** Reproduce every headline interaction with **LIANA+ and NicheNet**, not CellChat alone. A real edge should survive ≥2 independent inference methods. *[pending — currently CellChat only]*
2. **Within-study replication, not just pooled.** Our co-option shifts are strong pooled but only trends within-study (study-confounded). Require the sender-shift to hold **within the datasets that contain both AML and healthy** (currently 3), and report a Study-level Reproducibility Score (≥60% of studies). *[partially done; underpowered — needs more both-arm datasets]*
3. **Proper null + confound control.** (a) permutation/label-shuffle null for the sender-composition shift; (b) check the residual abundance confound flagged for the H2 uniform-mass result — does the number of *present* nodes differ AML vs healthy, and does controlling for it change the result. *[pending]*
4. **Robustness safeguards (R4).** Bootstrap CIs, per-edge barycenter support, and a simple-distance benchmark for the FGW/HDS results (e.g., the α=1 uniform-mass topological signal, p=0.006 but small). *[pending]*
5. **Malignant-cell-specific sender (the strongest test).** The sender-shift claim is most convincing if the conserved ligand is shown to be sent disproportionately **by malignant cells** in AML — i.e., the bin×malignant-state×role secondary graph. This needs the label fix (Numbat/VarTrix on the BAM subset; consensus expression callers elsewhere). *[pending — depends on Phase-1 label work]*
6. **Biological-plausibility check (already positive).** The literature in Sections 1–3 independently recovers the Gal-9/CD99-are-LSC vs MIF-is-broad split — a de-facto external validation of the pipeline's output. *[done — see this document]*

**Recommended order:** 6 (done) → 3b present-node confound → 1 LIANA+/NicheNet on the top ligands → 3a null + 4 bootstrap → 2 within-study/SRRS → 5 malignant-specific sender (after labels). Only edges surviving 1–5 go to the experimentalists.

---

## References (verified 2026-07-27)

**Galectin-9 / TIM-3**
- Kikushige Y, et al. (2010). TIM-3 is a promising target to selectively kill AML stem cells. *Cell Stem Cell* 7(6):708–717. PMID 21112565 / DOI 10.1016/j.stem.2010.11.014.
- Kikushige Y, et al. (2015). A TIM-3/Gal-9 autocrine stimulatory loop drives self-renewal of human myeloid leukemia stem cells. *Cell Stem Cell* 17(3):341–352. PMID 26279267 / DOI 10.1016/j.stem.2015.07.011. **(landmark)**
- Gonçalves Silva I, et al. (2017). The Tim-3–galectin-9 secretory pathway in immune escape of human AML. *EBioMedicine* 22:44–57. PMID 28750861 / DOI 10.1016/j.ebiom.2017.07.018.
- Yasinska IM, et al. (2020). Galectin-9 and VISTA suppress human T-lymphocyte cytotoxicity. *Front Immunol* 11:580557. PMID 33329552 / DOI 10.3389/fimmu.2020.580557.
- Asayama T, et al. (2017). Tim-3 on blasts and prognostic impact of galectin-9 in MDS. *Oncotarget* 8(51):88904–88917. PMID 29179486 / DOI 10.18632/oncotarget.21492.
- Zhang Y, et al. (2021). Galectin-9 and PSMB8 overexpression predict unfavorable prognosis in AML. *J Cancer* 12(14):4257–4263. PMID 34093826 / DOI 10.7150/jca.53686.
- Choukrani G, et al. (2023). Galectin-9 has non-apoptotic cytotoxic activity toward AML. *Cell Death Discov* 9. PMID 37407572 / DOI 10.1038/s41420-023-01515-w.
- Zeidan AM, et al. (2024). Sabatolimab + HMA in higher-risk MDS (STIMULUS-MDS1), phase 2. *Lancet Haematol* 11(1):e38–e50. DOI 10.1016/S2352-3026(23)00333-2. *(negative)*
- STIMULUS-MDS2 (NCT04266301): missed primary OS endpoint (EHA 2024, congress/press-level — not a full paper).
- Reviews (metadata-verified, full text not opened): *Mol Biol Rep* (2024) 51:571, DOI 10.1007/s11033-024-09563-w; "Galectin-9: a double-edged sword in AML," *Ann Hematol* (2025), DOI 10.1007/s00277-025-06387-x. *(author lists unverified)*

**CD99**
- Chung SS, et al. (2017). CD99 is a therapeutic target on disease stem cells in myeloid malignancies. *Sci Transl Med* 9(374):eaaj2025. PMID 28123069 / DOI 10.1126/scitranslmed.aaj2025. **(landmark)**
- Vaikari VP, et al. (2020). Clinical and preclinical characterization of CD99 isoforms in AML. *Haematologica* 105(4):999–1012. PMID 31371417 / DOI 10.3324/haematol.2018.207001.
- Travaglini S, Ottone T, et al. (2022). CD99 as a novel therapeutic target on leukemic progenitor cells in FLT3-ITDmut AML. *Leukemia* 36(6):1685–1688. PMID 35422094 / DOI 10.1038/s41375-022-01566-5.
- Schenkel AR, et al. (2002). CD99 in monocyte transendothelial migration. *Nat Immunol* 3(2):143–150. PMID 11812991. *(homophilic CD99, general)*
- Goswami D, et al. (2017). Endothelial CD99 binds neutrophil PILRs. *Blood* 129(13):1811–1822. PMID 28223280 / DOI 10.1182/blood-2016-08-733394. *(partner biology, general)*
- Ali A, et al. (2024). FLT3/CD99 bispecific nanoparticles for AML. *Cancer Res Commun* 4(8):1946–1962. PMID 39007347 / DOI 10.1158/2767-9764.CRC-24-0096. *(preclinical)*
- Shi J, et al. (2021). CAR-T targeting CD99 in T-ALL. *J Hematol Oncol* 14:162. PMID 34627328 / DOI 10.1186/s13045-021-01178-z. *(T-ALL context)*
- Kotemul K, et al. (2024). CD99 as a target for antibody therapy of T-ALL (review). *Explor Target Antitumor Ther* 5(1):96–107. PMID 38468825 / DOI 10.37349/etat.2024.00207.
- Guerzoni C, et al. (2015). CD99 triggering in Ewing sarcoma (p53 reactivation, +doxorubicin). *Clin Cancer Res* 21(1):146–156. PMID 25501132 / DOI 10.1158/1078-0432.CCR-14-0492. *(shared mechanism)*

**MIF**
- Leng L, et al. (2003). MIF signal transduction via CD74. *J Exp Med* 197(11):1467–1476. PMID 12782713 / DOI 10.1084/jem.20030286. **(receptor)**
- Shi X, et al. (2006). CD44 is the signaling component of the MIF–CD74 complex. *Immunity* 25(4):595–606. PMID 17045821 / DOI 10.1016/j.immuni.2006.08.020.
- Bernhagen J, et al. (2007). MIF is a noncognate ligand of CXCR2/CXCR4. *Nat Med* 13(5):587–596. PMID 17435771 / DOI 10.1038/nm1567.
- Abdul-Aziz AM, et al. (2017). MIF-induced stromal PKCβ/IL-8 is essential in human AML. *Cancer Res* 77(2):303–311. PMID 27872094 / DOI 10.1158/0008-5472.CAN-16-1095. **(key AML)**
- Spertini C, et al. (2024). MIF blockade reprograms macrophages and disrupts prosurvival signaling in AML. *Cell Death Discov* 10:157. PMID 38548753 / DOI 10.1038/s41420-024-01924-5. *(author spellings to reconfirm)*
- Pantouris G, et al. (2025). Allosteric MIF inhibitor triggers cell-cycle arrest in AML. *ACS Omega* 10:17441–17452. DOI 10.1021/acsomega.4c10969. *(authors/pages to reconfirm)*
- Menssen AJ, et al. (2024). CD74 in pediatric AML — target for therapy (COG). *Haematologica* 109(10):3182–3193. PMID 38299667 / DOI 10.3324/haematol.2023.283757.
- Martin P, et al. (2015). Phase I anti-CD74 milatuzumab in B-cell lymphoma. *Leuk Lymphoma* 56(11):3065–3070. PMID 25754579 / DOI 10.3109/10428194.2015.1028052. *(not AML)*
- Attar EC, et al. (2009). CD74 and inferior outcome in CN-AML (CALGB). *Blood* 114(22):1616 (ASH abstract). DOI 10.1182/blood.V114.22.1616.1616. *(abstract)*
- Anderson-Crannage M, et al. (2025). MIF–CD74/CD44 axis in pediatric AML. *Blood* 146(Suppl 1):6786 (ASH abstract). *(abstract, single-cell/CellChat)*
- NCT03918655 — "MIF: a new target to eradicate the pre-leukemic clone in AML" (AP-HP, observational).
