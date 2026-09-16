#!/usr/bin/env python
# 16_scaccordion_noise_floor.py ----
# INPUT  : <root>/05_ccc/tensors__split_{A,B}/<ds>/<sample>__ccc_cellchat.csv (37 samples, 2 halves)
#          <root>/05_ccc/tensors/<ds>/<sample>__ccc_cellchat.csv               (same 37, full depth)
#          <root>/06_distance/scaccordion_cost__<arm>.csv                      (cohort geometry)
#          <root>/06_distance/scaccordion_pmat.csv                             (the 39 kept slots)
#          <root>/05_ccc/ccc_sample_manifest.csv                               (uid_patient)
# OUTPUT : <root>/08_scoring/scaccordion_noise_floor.csv        (one row per arm)
#          <root>/08_scoring/scaccordion_noise_floor_pairs.csv  (per-sample retest distances)
# WHAT IT DOES : answers the one question section B left open -- is the ~0.56 median between-sample
#          distance MEASUREMENT NOISE, or real between-sample variation that simply is not disease?
#          Splitting one sample's cells in two and re-running the identical CellChat path gives the
#          measurement noise directly: two halves of ONE sample differ by measurement error alone.
#
# WHY THIS IS THE DECIDING NUMBER. Section B measured that a 5x planted effect on 6 of 49 edges
# moves the between-group distance by 3.4%, because sample-to-sample distances are already ~0.56 on
# a 0-1 scale. That is a ceiling, but it does not say what the 0.56 IS. Two readings, opposite
# consequences:
#   (a) it is measurement noise  -> the whole approach sits below its detection limit; nothing
#       downstream can work, and the line closes.
#   (b) it is real between-sample variation that happens not to be disease -> the measurement is
#       sound and the fix is a design that removes between-patient variance, i.e. the paired
#       within-patient comparison GATE 1 already runs.
#
# THE DEPTH CONFOUND, AND HOW IT IS CONTROLLED. Each half holds half the cells, so CellChat's
# permutation reaches significance on fewer LR pairs: measured here, the halves occupy a median
# of 21 of the cohort's 39 slots against the full sample's 25. A retest distance computed at half
# depth therefore MAY NOT be divided by a between-sample distance computed at full depth -- that
# ratio would credit the depth difference to reliability. This project has already made exactly
# that error once (the D5 edge-set mismatch, which moved the published ratios from 2.75/1.50 to
# 3.88/2.44). So every ratio below has BOTH of its terms computed at half depth, and full-depth
# numbers are carried only as a separately labelled reference column.
#
# Usage : PYTHONPATH=/FAST/gr10634/gaozy/external/pylibs \
#         python scripts/08_scoring/16_scaccordion_noise_floor.py
import glob
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "config"))
from scaccordion_geometry import tables_from_tensors, pmat                   # noqa: E402

CCC_NODES = ["HSC_MPP", "LMPP_GMP", "Mono_DC", "Erythroid", "Megakaryocyte", "T_NK", "B_Plasma"]
DEFAULT_ROOT = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
# prop7 is deliberately excluded: 05_ccc/06 splits cells STRATIFIED BY hierarchy_bin, so the two
# halves have identical composition by construction. Its retest distance would be ~0 by design
# rather than by measurement. That fact is instead used as SELF-CHECK 2, a control on the split.
ARMS = ["tv", "tabular", "corrot", "dwot_shipped", "dwot_fixed"]
COST_ARMS = {"corrot", "dwot_shipped", "dwot_fixed"}

root = DEFAULT_ROOT
D_CCC = os.path.join(root, "05_ccc")
D_DST = os.path.join(root, "06_distance")
D_SCO = os.path.join(root, "08_scoring")

# ---- 1. load the three depths onto ONE fixed slot set -----------------------------------------
SLOTS = list(pd.read_csv(os.path.join(D_DST, "scaccordion_pmat.csv"), index_col=0).index)
print("[0] cohort geometry: %d line-graph slots (from the 138-sample production run)" % len(SLOTS))


def mass(tensor_glob, keep=None):
    """LR tensors -> (39 slots x samples) mass matrix on the FIXED cohort slot set, zero-padded.

    Zero-padding onto the cohort's own slot set is scACCorDiON's stated convention for a sample
    that lacks an interaction, and it is what keeps all three depths in one geometry.
    """
    t = {k: v for k, v in tables_from_tensors(sorted(glob.glob(tensor_glob))).items() if len(v)}
    if keep is not None:
        t = {k: v for k, v in t.items() if k in keep}
    P = pmat(t, nodes=CCC_NODES).reindex(SLOTS).fillna(0.0)
    return P


PA = mass(os.path.join(D_CCC, "tensors__split_A", "*", "*__ccc_cellchat.csv"))
PB = mass(os.path.join(D_CCC, "tensors__split_B", "*", "*__ccc_cellchat.csv"))
IDS = sorted(set(PA.columns) & set(PB.columns))
PF = mass(os.path.join(D_CCC, "tensors", "*", "*__ccc_cellchat.csv"), keep=set(IDS))
IDS = [s for s in IDS if s in PF.columns]
PA, PB, PF = PA[IDS], PB[IDS], PF[IDS]
n = len(IDS)
print("[1] %d samples with both halves and a full-depth run" % n)

man = pd.read_csv(os.path.join(D_CCC, "ccc_sample_manifest.csv")).set_index("sample")
pat = man.reindex(IDS)["uid_patient"].to_numpy()
print("[1] %d patients; %d longitudinal same-patient pairs"
      % (len(set(pat)), sum(1 for i in range(n) for j in range(i + 1, n) if pat[i] == pat[j])))

## -- SELF-CHECK 1: the halves must live on the production slot set, or nothing is comparable ----
assert list(PA.index) == SLOTS == list(PF.index), "slot sets diverged between depths"
detA = (PA.to_numpy() > 0).sum(0)
detF = (PF.to_numpy() > 0).sum(0)
print("\n[SELF-CHECK 1] depth effect, the confound being controlled")
print("               edges detected, median: full %.0f | half A %.0f | half B %.0f"
      % (np.median(detF), np.median(detA), np.median((PB.to_numpy() > 0).sum(0))))
print("               every ratio below divides a half-depth number by a half-depth number.")


def normed(P):
    A = P.to_numpy(float)
    return A / np.where(A.sum(0, keepdims=True) > 0, A.sum(0, keepdims=True), 1.0)


def dist(X, Y, arm, C=None):
    """Pairwise distances between the columns of X and of Y, in the cohort geometry."""
    import ot
    if arm == "tv":
        return np.abs(X[:, :, None] - Y[:, None, :]).sum(0) / 2.0
    if arm == "tabular":
        return np.sqrt(((X[:, :, None] - Y[:, None, :]) ** 2).sum(0))
    C = np.ascontiguousarray(C)
    D = np.zeros((X.shape[1], Y.shape[1]))
    for i in range(X.shape[1]):
        for j in range(Y.shape[1]):
            D[i, j] = ot.emd2(X[:, i].copy(), Y[:, j].copy(), C)
    return D


## -- SELF-CHECK 2: a control on the SPLIT itself, not on the distance --------------------------
## 05_ccc/06 stratifies the split by hierarchy_bin, so the two halves must hold the same cell
## composition. If they do not, the retest distance is partly a composition difference and every
## ratio below is measuring the wrong thing. Cell counts per bin are read from the halves' own
## tensors (n_sender), not from ccc_node_presence.csv, which only covers full-depth samples.
def half_composition(sub):
    out = {}
    for f in sorted(glob.glob(os.path.join(D_CCC, sub, "*", "*__ccc_cellchat.csv"))):
        d = pd.read_csv(f)
        if not len(d):
            continue
        a = d.groupby("sender_bin")["n_sender"].first().reindex(CCC_NODES).fillna(0.0)
        out[str(d["sample"].iloc[0])] = a.to_numpy(float)
    return pd.DataFrame(out).T.reindex(IDS).fillna(0.0).to_numpy(float)


cA, cB = half_composition("tensors__split_A") + 0.5, half_composition("tensors__split_B") + 0.5
fa, fb = cA / cA.sum(1, keepdims=True), cB / cB.sum(1, keepdims=True)
comp_tv = np.abs(fa - fb).sum(1) / 2.0
print("\n[SELF-CHECK 2] did the bin-stratified split equalise composition?")
print("               composition TV between halves: median %.4f  max %.4f"
      % (np.median(comp_tv), comp_tv.max()))
print("               cell-count ratio A/B: median %.3f" % np.median(cA.sum(1) / cB.sum(1)))
if comp_tv.max() > 0.05:
    raise SystemExit("the halves differ in composition; the retest distance is confounded")
print("               -> the halves differ in WHICH cells landed where, not in composition, so")
print("                  the retest distance is pure ligand-receptor detection noise.")

# ---- 2. the three quantities, all at matched depth --------------------------------------------
XA, XB, XF = normed(PA), normed(PB), normed(PF)
iu = np.triu_indices(n, 1)
same_pat = np.array([[pat[i] == pat[j] for j in range(n)] for i in range(n)])
rows, pairrows = [], []
print("\n[2] %-13s %10s %12s %12s %10s %10s" % ("arm", "retest", "between", "within-pt",
                                                "SNR_betw", "SNR_within"))
for arm in ARMS:
    C = None
    if arm in COST_ARMS:
        C = pd.read_csv(os.path.join(D_DST, "scaccordion_cost__%s.csv" % arm), index_col=0)
        C = C.reindex(index=SLOTS, columns=SLOTS).to_numpy(float)

    DAB = dist(XA, XB, arm, C)                                   # every A against every B
    retest = np.diag(DAB)                                        # d(A_i, B_i), half depth
    DAA = dist(XA, XA, arm, C)                                   # d(A_i, A_j), half depth
    DFF = dist(XF, XF, arm, C)                                   # full depth, reference only

    ## -- SELF-CHECK 3: the retest is an A-vs-B comparison while the baseline is A-vs-A. If the
    ## -- two halves differed systematically, the retest would be inflated against its own
    ## -- baseline and every SNR below would be understated. Compare like with like.
    off = ~np.eye(n, dtype=bool) & ~same_pat
    ab_between, aa_between = float(np.median(DAB[off])), float(np.median(DAA[off]))
    if aa_between > 0 and abs(ab_between - aa_between) / aa_between > 0.05:
        raise SystemExit("%s: systematic A/B offset (%.4f vs %.4f); retest is not a fair baseline"
                         % (arm, ab_between, aa_between))

    betw_half = DAA[iu][~same_pat[iu]]
    within_half = DAA[iu][same_pat[iu]]
    betw_full = DFF[iu][~same_pat[iu]]
    within_full = DFF[iu][same_pat[iu]]

    med_r = float(np.median(retest))
    med_b = float(np.median(betw_half))
    med_w = float(np.median(within_half))

    # Paired, per sample: is sample i further from OTHER samples than from its own other half?
    othermed = np.array([np.median(DAA[i, [j for j in range(n) if j != i and pat[j] != pat[i]]])
                         for i in range(n)])
    w = stats.wilcoxon(othermed, retest, alternative="greater")

    rows.append(dict(arm=arm, n_samples=n, retest_half_median=med_r,
                     between_half_median=med_b, within_patient_half_median=med_w,
                     snr_between=med_b / med_r if med_r > 0 else np.inf,
                     snr_within_patient=med_w / med_r if med_r > 0 else np.inf,
                     between_full_median=float(np.median(betw_full)),
                     within_patient_full_median=float(np.median(within_full)),
                     ab_vs_aa_between=ab_between / aa_between if aa_between > 0 else np.nan,
                     composition_tv_max=float(comp_tv.max()),
                     wilcoxon_stat=float(w.statistic), wilcoxon_p=float(w.pvalue)))
    for i, s in enumerate(IDS):
        pairrows.append(dict(arm=arm, sample=s, uid_patient=pat[i], retest=float(retest[i]),
                             median_to_other_patients=float(othermed[i]),
                             n_edges_half=int(detA[i])))
    print("    %-13s %10.4f %12.4f %12.4f %10.2f %10.2f"
          % (arm, med_r, med_b, med_w, rows[-1]["snr_between"], rows[-1]["snr_within_patient"]))

R = pd.DataFrame(rows)
R.to_csv(os.path.join(D_SCO, "scaccordion_noise_floor.csv"), index=False)
pd.DataFrame(pairrows).to_csv(os.path.join(D_SCO, "scaccordion_noise_floor_pairs.csv"), index=False)

# ---- 3. the verdict ---------------------------------------------------------------------------
print("\n[3] paired Wilcoxon, per sample: median distance to other patients > own retest distance?")
for _, r in R.iterrows():
    print("    %-13s V=%8.0f  p=%.2e  %s"
          % (r.arm, r.wilcoxon_stat, r.wilcoxon_p,
             "above noise" if r.wilcoxon_p < 0.05 else "AT NOISE"))
print("\n[3] READING (registered in advance of looking):")
print("    SNR_between near 1.0  -> the 0.56 IS measurement noise; the approach is below its")
print("                             detection limit and the line closes there.")
print("    SNR_between >> 1.0    -> between-sample variation is REAL but is not disease; the")
print("                             measurement is sound and a within-patient design is the fix.")
print("    SNR_within_patient    -> whether longitudinal change clears the same floor, at the")
print("                          DISTANCE level (section 2 answered this per-edge: 3.88/2.44).")
print("\n[done] %d arms x %d samples" % (len(ARMS), n))
