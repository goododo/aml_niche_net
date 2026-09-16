#!/usr/bin/env python
# 17_scaccordion_gate1.py ----
# INPUT  : <root>/06_distance/scaccordion_distance__<arm>.csv  (06_distance/03, six arms)
#          <root>/05_ccc/ccc_sample_manifest.csv               (uid_patient, Timepoint, dataset)
#          <root>/07_fgw/fgw_vocab.json                        (the locked timepoint vocabulary)
#          <root>/08_scoring/scaccordion_noise_floor.csv       (16's retest floor, for context)
# OUTPUT : <root>/08_scoring/scaccordion_gate1.csv             (one row per arm x contrast)
#          <root>/08_scoring/scaccordion_gate1_ranks.csv       (one row per patient per cell)
# WHAT IT DOES : runs this project's pre-registered GATE 1 -- identity -- with the scACCorDiON
#          distance substituted for FGW2. Registered in PREREGISTRATION_scaccordion.md section 8.
#
# WHY THIS CELL AND NOT ANOTHER. Section B.6 established that the measurement is real (between
# sample distances are 3.9-5.6x the split-half noise floor) but that the variation it captures is
# between patients and is not disease. That leaves exactly one design whose variance structure is
# favourable: a within-patient comparison, which is what GATE 1 is. Every other cell in section A
# and section B has now been closed, so this is the last registered question standing.
#
# THE RULE IS COPIED, NOT REINVENTED. GATE 1 is defined in PREREGISTRATION_paired_gate.md and
# implemented in 11_paired_gate.py lines 213-239. Reproduced here exactly:
#   - for a patient with graphs at two timepoints, is the pair closer to each other than to other
#     patients' graphs;
#   - candidate pool = SAME DATASET, OTHER PATIENTS, AML timepoints only. Without the same-dataset
#     restriction the test is won by batch and means nothing;
#   - statistic = percentile rank of the true partner in that pool, uniform on [0,1] under the null;
#   - PASS = one-sided Wilcoxon signed-rank of the percentiles against 0.5, p < 0.05;
#   - pool < 3 distractors skips the patient; fewer than 5 patients skips the cell.
# The only substitution is the distance. FGW's alpha grid has no analogue here, so the six
# representation arms take its place as the axis that is swept.
#
# Usage : python scripts/08_scoring/17_scaccordion_gate1.py
import json
import os
import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
ARMS = ["prop7", "tabular", "tv", "corrot", "dwot_shipped", "dwot_fixed"]
REL_TP = {"Relapse", "Relapse2"}                       # 11_paired_gate.py lines 59-61, copied
TRT_TP = {"On_treatment", "Post_induction", "Post_consolidation",
          "Post_treatment_unspecified", "Refractory"}
MIN_POOL = 3
MIN_PATIENTS = 5
SEED = 491638

root = DEFAULT_ROOT
D_DST = os.path.join(root, "06_distance")
D_SCO = os.path.join(root, "08_scoring")
rng = np.random.default_rng(SEED)

AML_TP = set(json.load(open(os.path.join(root, "07_fgw", "fgw_vocab.json")))["aml_timepoints"])

D0 = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__%s.csv" % ARMS[0]), index_col=0)
SAMPLES = list(D0.index)
man = pd.read_csv(os.path.join(root, "05_ccc", "ccc_sample_manifest.csv")).set_index("sample")
idx = man.reindex(SAMPLES)[["dataset", "uid_patient", "Timepoint"]]
idx.columns = ["dataset", "patient", "timepoint"]
print("[0] %d samples, %d datasets, %d patients" %
      (len(SAMPLES), idx.dataset.nunique(), idx.patient.nunique()))

# ---- pairs, built exactly as 11_paired_gate.py builds them ------------------------------------
pairs = []
for (ds, pt), g in idx.groupby(["dataset", "patient"]):
    tps = {t: s for s, t in zip(g.index, g.timepoint)}
    if "Diagnosis" not in tps:
        continue
    for kind, want in [("Dx_to_Relapse", REL_TP), ("Dx_to_Treatment", TRT_TP)]:
        hit = [t for t in g.timepoint if t in want]
        if hit:
            pairs.append(dict(kind=kind, dataset=ds, patient=pt,
                              a=tps["Diagnosis"], b=tps[sorted(hit)[0]]))
P = pd.DataFrame(pairs)
print("[1] pairs: %s" % P.kind.value_counts().to_dict())
print("[1] SAMPLE SET: %d distinct samples enter as a pair member, from %d patients"
      % (len(set(P.a) | set(P.b)), P.patient.nunique()))

floor = pd.read_csv(os.path.join(D_SCO, "scaccordion_noise_floor.csv")).set_index("arm")


# Pools are a property of the roster, not of the distance, so they are built once. Doing this
# inside the loop made the 1000-draw family null take minutes of pandas .at lookups.
POS = {s: i for i, s in enumerate(SAMPLES)}
POOLS = {}
for _, r in P.iterrows():
    POOLS[(r.kind, r.patient)] = np.array(
        [POS[s] for s in SAMPLES
         if idx.at[s, "dataset"] == r.dataset and idx.at[s, "patient"] != r.patient
         and idx.at[s, "timepoint"] in AML_TP], dtype=int)


def gate1(Dnp, grp, shuffle=False):
    """Percentile rank of the true partner among same-dataset other-patient AML graphs."""
    pcts, top1, top5, rows = [], 0, 0, []
    for _, r in grp.iterrows():
        pool = POOLS[(r.kind, r.patient)]
        if len(pool) < MIN_POOL:
            continue
        ia = POS[r.a]
        ib = int(rng.choice(pool)) if shuffle else POS[r.b]    # negative control
        d_true = float(Dnp[ia, ib])
        d_pool = Dnp[ia, pool[pool != ib]]
        if len(d_pool) < MIN_POOL:
            continue
        # Rounded, and the rounding is load-bearing. A percentile is a ratio of small integers, so
        # the signed-rank test sees exact ties: with a pool of 23, pct 4/23 and pct 19/23 sit the
        # same distance either side of 0.5 (both 15/46). In binary those two distances differ by
        # one ULP, which makes scipy treat a mathematically real tie as two distinct ranks. On this
        # data that single ULP moves the prop7 Dx-to-Relapse p between 0.0034 and 0.0049 -- the
        # verdict is unchanged, but the number is not reproducible across a CSV round-trip unless
        # the tie is made deterministic here.
        pct = round(float((d_pool < d_true).sum()) / len(d_pool), 12)
        pcts.append(pct)
        top1 += int(pct == 0.0)
        top5 += int((d_pool < d_true).sum() < 5)
        rows.append(dict(dataset=r.dataset, patient=r.patient, n_pool=len(d_pool),
                         d_true=d_true, d_pool_median=float(np.median(d_pool)), pct=pct))
    return np.array(pcts), top1, top5, rows


# ---- the registered test ----------------------------------------------------------------------
res, rankrows = [], []
print("\n[2] %-13s %-16s %3s %9s %9s %8s %8s %10s"
      % ("arm", "contrast", "n", "med pct", "p", "top1", "top5", "d_true/floor"))
for arm in ARMS:
    Dnp = pd.read_csv(os.path.join(D_DST, "scaccordion_distance__%s.csv" % arm),
                      index_col=0).reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
    for kind, grp in P.groupby("kind"):
        pcts, t1, t5, rows = gate1(Dnp, grp)
        if len(pcts) < MIN_PATIENTS:
            print("    [skip] %s %s: only %d patients" % (arm, kind, len(pcts)))
            continue
        w = stats.wilcoxon(pcts - 0.5, alternative="less")
        # chance rate for top-1 is the mean of 1/pool_size, not a single number, because pools differ
        chance1 = float(np.mean([1.0 / r["n_pool"] for r in rows]))
        chance5 = float(np.mean([min(5, r["n_pool"]) / r["n_pool"] for r in rows]))
        # how big is the within-patient distance relative to THIS arm's split-half noise floor?
        fl = float(floor.loc[arm, "retest_half_median"]) if arm in floor.index else np.nan
        ratio = float(np.median([r["d_true"] for r in rows]) / fl) if fl and fl > 0 else np.nan
        res.append(dict(arm=arm, contrast=kind, n=len(pcts), median_pct=float(np.median(pcts)),
                        p_raw=float(w.pvalue), top1=t1, top1_chance=chance1 * len(pcts),
                        top5=t5, top5_chance=chance5 * len(pcts),
                        d_true_median=float(np.median([r["d_true"] for r in rows])),
                        retest_floor=fl, d_true_over_floor=ratio))
        for r in rows:
            rankrows.append(dict(arm=arm, contrast=kind, **r))
        print("    %-13s %-16s %3d %9.3f %9.4f %4d/%-3d %4d/%-3d %10s"
              % (arm, kind, len(pcts), np.median(pcts), w.pvalue, t1, len(pcts), t5, len(pcts),
                 "%.2f" % ratio if ratio == ratio else "n/a"))

R = pd.DataFrame(res)
# BH-FDR over the 12 registered cells (6 arms x 2 contrasts). PREREGISTRATION_scaccordion.md 7.4.
o = np.argsort(R.p_raw.to_numpy())
m = len(R)
adj = np.empty(m)
adj[o] = np.minimum.accumulate((R.p_raw.to_numpy()[o] * m / (np.arange(m) + 1))[::-1])[::-1]
R["p_bh"] = np.minimum(adj, 1.0)
R.to_csv(os.path.join(D_SCO, "scaccordion_gate1.csv"), index=False)
pd.DataFrame(rankrows).to_csv(os.path.join(D_SCO, "scaccordion_gate1_ranks.csv"), index=False)

## -- FAMILY-LEVEL NULL. Two problems the per-cell p cannot handle, both handled here by the same
## -- permutation. (1) The six arms are six versions of ONE distance on ONE dataset -- their
## -- correlation against total variation runs 0.89 to 1.00 -- so BH over twelve cells treats
## -- highly dependent tests as independent and is far too conservative. (2) A single shuffled
## -- draw per cell is a weak control. Replacing the true partner with a random same-dataset AML
## -- graph, many times, gives the null distribution of the WHOLE table at once, which is the
## -- statistic actually being read ("Dx-to-Relapse passes in all six arms"). This is the same
## -- device section 2 used to settle the per-edge family count.
N_PERM = 1000
print("\n[SELF-CHECK] family-level null: true partner replaced by a random same-dataset AML")
print("             graph, %d draws. Counts how often the shuffled table looks like the real one."
      % N_PERM)
DBANK = {a: pd.read_csv(os.path.join(D_DST, "scaccordion_distance__%s.csv" % a),
                        index_col=0).reindex(index=SAMPLES, columns=SAMPLES).to_numpy(float)
         for a in ARMS}
CELLS = [(a, k) for a in ARMS for k, _ in P.groupby("kind")]
obs_pass = int((R.p_raw < 0.05).sum())
obs_rel = int(((R.p_raw < 0.05) & (R.contrast == "Dx_to_Relapse")).sum())
null_pass, null_rel = np.zeros(N_PERM, int), np.zeros(N_PERM, int)
for t in range(N_PERM):
    for arm in ARMS:
        for kind, grp in P.groupby("kind"):
            pcts, _, _, _ = gate1(DBANK[arm], grp, shuffle=True)
            if len(pcts) < MIN_PATIENTS:
                continue
            if float(stats.wilcoxon(pcts - 0.5, alternative="less").pvalue) < 0.05:
                null_pass[t] += 1
                null_rel[t] += int(kind == "Dx_to_Relapse")
p_fam = float((null_pass >= obs_pass).mean())
p_rel = float((null_rel >= obs_rel).mean())
print("    cells passing raw p<0.05 : observed %d of 12 ; null mean %.2f ; p = %.4f"
      % (obs_pass, null_pass.mean(), p_fam))
print("    of those, Dx_to_Relapse  : observed %d of 6  ; null mean %.2f ; p = %.4f"
      % (obs_rel, null_rel.mean(), p_rel))
R["family_p_all_cells"] = p_fam
R["family_p_relapse"] = p_rel
R.to_csv(os.path.join(D_SCO, "scaccordion_gate1.csv"), index=False)

print("\n[3] REGISTERED READING")
print("    PASS = one-sided Wilcoxon p < 0.05 with the true partner ranking closer.")
print("    Cells passing raw: %d of %d ; surviving BH over all %d: %d"
      % (int((R.p_raw < 0.05).sum()), len(R), len(R), int((R.p_bh < 0.05).sum())))
print("    d_true_over_floor < 1 means the within-patient distance is BELOW this arm's own")
print("    split-half measurement noise, and the cell cannot be believed whatever its p.")
print("\n[done] %d cells" % len(R))
