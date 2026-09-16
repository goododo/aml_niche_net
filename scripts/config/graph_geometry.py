#!/usr/bin/env python
# graph_geometry.py ----
# The damped-random-walk geometry of a directed communication graph. ONE function, nothing else.
#
# WHY THIS FILE EXISTS. Two call sites need the identical construction: 06_distance/02 builds it
# for the cohort, and 08_scoring/10 must rebuild it from PLANTED weights to ask whether the new
# geometry can detect an effect the old one could not. If the two drifted, the planting harness
# would be testing a formula the pipeline does not use, and would do so silently. This is the same
# reason distance_variants.py and fgw_vocab.py exist, and it sits next to them for the same reason.
#
# This is NOT a utils module. It holds the walk geometry and will not be given a second
# responsibility.
import numpy as np

WALK_ALPHA = 0.85      # probability the walk continues at each step
WALK_TELEPORT = 0.05   # uniform jump mixed into every row, so the chain is irreducible
EPS_MASS = 1e-6        # matches config_fgw.R's treatment of an absent node


def walk_geometry(Wm, alpha=WALK_ALPHA, teleport=WALK_TELEPORT, eps_mass=EPS_MASS):
    """One sample's directed weight matrix -> (R, C_geo, mass).

    Wm  (n, n) directed edge weights, senders on rows. Zeros allowed.

    R       (n, n) damped reachability. R[i, j] is the probability that a walk starting at i,
            continuing with probability `alpha` and stopping with probability (1 - alpha), stops
            at j. R = (1 - alpha) * inv(I - alpha * P). Row-stochastic, so it is on the same
            scale in every sample however deeply that sample was sequenced. Row-normalising W
            into P is what removes the depth scale here, which is the job the within-sample rank
            transform used to do -- and unlike ranking it leaves the value spectrum intact.
    C_geo   -log(R). The GW structure matrix: a GLOBAL geometry (how hard is it for signal to
            reach j from i through the whole network) rather than a per-edge readout.
    mass    node mass from total signal strength (sent + received), normalised to sum 1. This is
            the half of the swap that matters as much as the cost: the signal belongs in the
            mass, and cell counts -- known here to track sample preparation -- do not.

    A node with no outgoing signal has an undefined transition row; it is completed to uniform,
    the standard dangling-node convention. The teleport term makes that a limiting case of the
    same rule rather than a special one, and guarantees every R entry is strictly positive so
    -log needs no arbitrary floor. Both are properties of this construction, not tuned choices;
    06_distance/02 reports how many rows and entries they touch.
    """
    Wm = np.asarray(Wm, float)
    n = Wm.shape[0]
    if Wm.shape != (n, n):
        raise ValueError("walk_geometry expects a square matrix, got %r" % (Wm.shape,))
    s = Wm.sum(1)
    P = np.where(s[:, None] > 0, Wm / np.where(s[:, None] > 0, s[:, None], 1.0), 1.0 / n)
    P = (1.0 - teleport) * P + teleport / n
    R = (1.0 - alpha) * np.linalg.inv(np.eye(n) - alpha * P)
    stre = Wm.sum(1) + Wm.sum(0)
    m = np.where(stre > 0, stre, eps_mass)
    return R, -np.log(R), m / m.sum()
