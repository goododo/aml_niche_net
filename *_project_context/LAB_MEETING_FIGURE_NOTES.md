# Lab-meeting deck — figure conclusions (speaker notes)

*The deck figures (`results/figures/12_deck/d0*.png`) show only the visualization — big fonts, no baked-in sentences. The conclusion for each figure is here, so you can read/edit them separately and use them as talking points. Plain language throughout.*

---

**d01 — Samples remaining at each step.**
We start with 220 QC-passing samples. 130 get a cancer-cell label, 148 have a communication network, 141 enter the model, and the strictest (platform-controlled) test rests on just 60 samples from 3 datasets.
→ *Takeaway:* the headline tests lean on a modest, uneven slice of the data — worth stating up front.

**d06 — Who signals to whom, and with which molecules. (NEW)**
Across the 121 AML samples: the busiest senders and receivers are the immature stem/progenitor cells and the monocytes/dendritic cells. A handful of ligand-receptor pairs dominate — **MIF** (→ CD74/CXCR4/CD44, in ~90–98% of samples), **CD99**, and **Galectin-9** (LGALS9 → PTPRC/CD44, ~79%). On the stem→immune route specifically, MIF and Galectin-9 are again the leads.
→ *Takeaway:* this is the actual content of the networks (we hadn't shown it before). MIF and Galectin-9 are both classic immune-dampening signals in AML — a plausible immune-evasion channel, and exactly what H1 tests next.

**d03 — H1: AML vs healthy, per link.**
No communication link is reliably different between AML and healthy marrow. 0 of 49 sender→receiver links pass the significance test; grouped into 5 families, all sit at q = 0.93 (far above the 0.05 line). The stem→immune route leans slightly stronger in AML, but not significantly (p = 0.86).
→ *Takeaway:* a genuine "no signal," not a near-miss. And because it never uses the cancer labels, the labelling problem (next figure) can't distort it — so this is our most trustworthy result. Caveat: CellChat only so far; repeat with LIANA+/NicheNet before final.

**d02 — Label quality.**
125 of 130 samples rest on a single method (inferCNV); the plan was three methods agreeing. Tested on healthy marrow (where any "cancer" call is wrong by definition), the stem-cell type is mislabelled ~40% of the time — exactly the cell type our hypothesis needs most. Immune cells are clean (3–8%). Stroma is 90% wrong, which is why we exclude it.
→ *Takeaway:* the foundation — deciding which cells are cancerous — is the weak spot, and it sits right where the biology matters.

**d04 — H2: network shape vs cell features.**
When we rely on the network shape alone, the AML-vs-healthy signal disappears (p = 0.11, not significant). All the signal comes from the cell features — and mostly from the blast-fraction feature, which we hard-set to 0 in healthy samples, so it just re-states "AML vs healthy." Stemness (the feature we actually care about) does nothing.
→ *Takeaway:* today the signal is cell counts, not network shape. The fair test needs the planned upgrades (honest richer features; a "direction" score) — not yet built.

**d05 — H3: treatment trajectory.**
Stemness dips at residual disease (MRD) and rebounds at relapse — a V-shape. This is not a failure: the updated plan pre-registered this V as one of three allowed patterns, and it's the one expected at MRD. But we can't judge it yet — the scores that tell the three patterns apart aren't built, one patient contributes 62% of the residual-disease cells, and the timepoint labels still need checking against the source papers.
→ *Takeaway:* on hold, for the right reasons — not because the trend "failed."
