#!/usr/bin/env python
# distance_variants.py ----
# The four edge-weight-to-cost transforms pre-registered in
# scripts/08_scoring/PREREGISTRATION_paired_gate.md. ONE function, four arms, nothing else.
#
# WHY THIS FILE EXISTS. Three call sites need these transforms: 08_scoring/11_paired_gate.py,
# 08_scoring/06_alpha_sweep.py and 08_scoring/10_planted_effect_power.py. Re-listing them in
# three places is exactly how the hard-coded timepoint vocabulary and the hard-coded family
# map both went wrong in this project -- a rule was changed in one place and the other two
# silently kept scoring the old model. Same reason fgw_vocab.py exists, and this file sits
# next to it for the same reason.
#
# This is NOT a utils module. It holds the four registered transforms and will not be given
# a second responsibility.
import numpy as np
from scipy import stats

ARMS = ("rank", "mask", "const", "logs")
LOG_EPS = 1e-3   # pre-registered. See the amendment in PREREGISTRATION_paired_gate.md.


def weights_to_C(w, arm, n_lr_sig=None, eps=LOG_EPS):
    """Edge weights -> costs, for one sample.

    w         (n,) weight_probsum over that sample's directed edges, 0 where undetected.
    arm       one of ARMS.
    n_lr_sig  (n,) count of significant ligand-receptor pairs. Required by 'mask' only,
              which needs the size-independent weight_per_lr = w / max(n_lr_sig, 1).

    Returns a (n,) cost vector in [0, 1]. Strongest edge -> 0 (close), weakest -> 1 (far).

    An undetected edge is defined as w == 0. Under the planted-effect harness an edge that
    was undetected can be pushed above 0; it then counts as detected, and for 'mask' its
    weight_per_lr uses max(n_lr_sig, 1) = 1. That convention is stated here rather than
    left implicit at the three call sites.
    """
    w = np.asarray(w, float)
    n = len(w)
    det = w > 0

    if arm == "rank":
        # The production rule from 06_distance/01: rank all edges together, absent ones tie
        # at the bottom and take the average rank, so their cost moves with how many there are.
        return 1.0 - stats.rankdata(w, method="average") / n

    if arm == "const":
        # Detected edges keep the production ranking; absent edges take a fixed cost.
        c = 1.0 - stats.rankdata(w, method="average") / n
        c[~det] = 1.0
        return c

    if arm == "mask":
        # Rank only the detected edges, among themselves, on the size-independent weight.
        # Absent edges take the same fixed cost as 'const'. The two arms therefore differ
        # in the DETECTED edges, not the absent ones.
        if n_lr_sig is None:
            raise ValueError("arm 'mask' needs n_lr_sig")
        c = np.ones(n)
        if det.any():
            wpl = w[det] / np.maximum(np.asarray(n_lr_sig, float)[det], 1.0)
            c[det] = 1.0 - stats.rankdata(wpl, method="average") / det.sum()
        return c

    if arm == "logs":
        # The blueprint's second formula, never implemented until now:
        # normalise within sample, take -log, then min-max to [0, 1]. An absent edge lands at
        # -log(eps), which is the maximum, so it takes cost exactly 1 in every sample
        # regardless of how many are missing. That constancy is the point of the arm.
        mx = w.max()
        s = w / mx if mx > 0 else w
        raw = -np.log(s + eps)
        rng = raw.max() - raw.min()
        return (raw - raw.min()) / rng if rng > 0 else np.zeros(n)

    raise ValueError("unknown arm %r; expected one of %s" % (arm, ARMS))
