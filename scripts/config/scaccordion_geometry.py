#!/usr/bin/env python
# scaccordion_geometry.py ----
# The scACCorDiON construction (Nagai, Maie, Schaub, Costa; Bioinformatics 2025, 41(5), btaf288),
# in the two forms this project has to keep apart: exactly as their code ships it, and exactly as
# their paper describes it. ONE construction, nothing else.
#
# WHY THIS FILE EXISTS. Two call sites need the identical construction: 06_distance/03 builds it
# for the cohort, and 08_scoring/15 must rebuild it from PLANTED ligand-receptor tables to ask
# whether the published method can see an effect that was deliberately put there. If those two
# drifted, the positive control would be certifying a construction the benchmark does not use,
# and would do so silently. That is the same failure this project has already had twice, and the
# same reason graph_geometry.py and distance_variants.py exist. This file sits next to them.
#
# This is NOT a utils module. It holds the scACCorDiON construction and will not be given a
# second responsibility.
#
# WHY TWO ARMS -- the shipped code and the published text disagree in three places. All three
# were verified on 2026-09-16 against upstream commit 9d8e7d65fa335b2186473de8c07b69b930b96c59
# (vendored, read-only, at SCACC_REPO below):
#
#   (1) SHARED-TOPOLOGY EDGE WEIGHT, INVERTED.
#       Paper S3.2: "we defined the weight w_{u'v'} ... as the proportion of graphs containing
#       both edges u' to v'", whose stated purpose is "This makes transport of masses between
#       common edges ... more likely than transport between rare edges".
#       Code (Accordion.py): score = (n_samples + 1) - n_co_present ; weight = score / n_samples.
#       That is HIGH when two edges RARELY co-occur -- the opposite sign. Because the weight
#       becomes a transition probability, the shipped code makes transport between RARE edge
#       pairs cheap.
#
#   (2) LINE-GRAPH EDGE RULE, WIDENED.
#       Paper S3.2: an arc exists "if the target of u' is the source of v'" -- head-to-tail only.
#       Code: `if tmpi[1] == tmpj[0] or tmpi[0] == tmpj[1]` -- head-to-tail OR tail-to-head.
#
#   (3) PAGERANK REGULARISATION, ABSENT.
#       Paper S3.2 adds "a low-rank regularization term, as popularized within the context of the
#       well-known PageRank algorithm ... this guarantees the global reachability of all nodes in
#       the STG, which is required to compute the distance between nodes in a graph."
#       Code (compute_cost) adds only `tmp += 1e-10 * ones`. Measured consequence on this cohort:
#       with the paper's narrower edge rule (2) and no teleport, the hitting-time cost returns 12
#       non-finite entries. The teleport is load-bearing exactly as the paper says.
#
# A fourth defect is an ordering bug rather than a paper disagreement, and is NOT reproduced in
# the `paper` arm: Accordion.compute_cost builds the cost with nx.to_numpy_array(self.expgraph),
# whose node order is networkx insertion order, while compute_wassestein passes marginals in
# self.p.index order. On this cohort 33 of 39 positions differ, so ot.emd2 receives a cost matrix
# permuted with respect to its own marginals. stg_paper() below returns an adjacency already
# indexed by p.index, and 06_distance/03 asserts the alignment.
import os
import numpy as np
import pandas as pd

SCACC_REPO = "/FAST/gr10634/gaozy/external/scACCorDiON"   # vendored, pinned; see .VENDORED_SHA
SCACC_PYLIBS = "/FAST/gr10634/gaozy/external/pylibs"      # kmedoids/genieclust/pydiffmap only

PSEUDO = 1e-10          # Accordion's own `pseudo` default; an edge counts as present above this
VAR_FILTER_Q = 0.20     # Accordion's own filter=0.2, filter_mode='edge' default
HTD_BETA = 0.5          # the arm their compute_cost_all lists first and the paper describes
TELEPORT = 0.05         # the paper's PageRank term. Matches graph_geometry.WALK_TELEPORT so the
                        # two hand-built geometries in this project use one teleport, not two.


def scaccordion_on_path():
    """Put the vendored package and its three extra dependencies on sys.path. Idempotent.

    Kept as a function rather than executed at import so that a caller which only needs the
    paper-faithful arm (which is implemented here, in numpy) does not require the install.
    """
    import sys
    for p in (SCACC_REPO, SCACC_PYLIBS):
        if p not in sys.path:
            sys.path.insert(0, p)


# ----------------------------------------------------------------------------------------------
# 1. adapter: our per-sample ligand-receptor tensors -> their `tbls` input
# ----------------------------------------------------------------------------------------------
def tables_from_tensors(paths, pval_thresh=0.05):
    """Per-sample CellChat tensors -> {sample: DataFrame[source, target, lr_means]}.

    paths         iterable of <sample>__ccc_cellchat.csv written by 05_ccc/02.
    pval_thresh   keep ligand-receptor rows below this. This is NOT a scACCorDiON choice; it is
                  06_distance/01's locked production filter (`sig <- t[pval < DIST_PVAL_THRESH]`),
                  reused so that their groupby-sum reproduces our weight_probsum exactly. That
                  identity is asserted by 06_distance/03 SELF-CHECK 1, and it is the only reason
                  a distance built by their code is comparable to the eleven configurations
                  already in FINDINGS_topology_null.md section A.

    Samples whose significant-LR table is empty are returned as empty frames rather than dropped,
    so the caller decides what to do with them instead of losing them silently.
    """
    out = {}
    for f in paths:
        d = pd.read_csv(f)
        sample = str(d["sample"].iloc[0]) if len(d) else os.path.basename(f).split("__")[0]
        d = d.loc[d["pval"] < pval_thresh, ["sender_bin", "receiver_bin", "prob"]]
        out[sample] = d.rename(columns={"sender_bin": "source", "receiver_bin": "target",
                                        "prob": "lr_means"}).reset_index(drop=True)
    return out


def pmat(tbls, nodes=None):
    """{sample: LR table} -> (line-graph slot x sample) edge-weight matrix.

    Reproduces Accordion.__init__ + utils.graphs_to_pmat: sum lr_means over the LR pairs
    supporting each ordered (source, target), then lay the result out over the full Cartesian
    product of the node universe with absent slots zero-padded. Slot labels are "source$target",
    their convention, because the line-graph edge rule parses that separator.

    nodes  fix the node universe explicitly. Pass it when planting, so a planted edge cannot
           enlarge the universe and change the geometry underneath the comparison.
    """
    agg = {}
    seen = set()
    for k, v in tbls.items():
        if len(v):
            g = v.groupby(["source", "target"], sort=False)["lr_means"].sum()
            agg[k] = {f"{a}${b}": float(w) for (a, b), w in g.items()}
            seen.update(v["source"]); seen.update(v["target"])
        else:
            agg[k] = {}
    uni = sorted(seen) if nodes is None else list(nodes)
    slots = [f"{i}${j}" for i in uni for j in uni]
    p = pd.DataFrame(0.0, index=slots, columns=list(tbls))
    for k, d in agg.items():
        for s, w in d.items():
            if s in p.index:
                p.at[s, k] = w
    return p.sort_index()


def variance_filter(p, q=VAR_FILTER_Q):
    """Accordion's filter_mode='edge': drop slots whose across-sample variance is below the
    q-quantile, then drop all-zero slots. Reported, never silent -- 06_distance/03 prints which
    slots this removes, because on a 49-slot graph it removes ten of them.
    """
    p = p.loc[p.T.var() > np.quantile(p.T.var(), q=q), :]
    return p.loc[p.sum(axis=1) != 0, :]


# ----------------------------------------------------------------------------------------------
# 2. the shared topology graph, both ways
# ----------------------------------------------------------------------------------------------
def _composable(a, b, rule):
    """Is there a line-graph arc from slot a to slot b?"""
    ta, tb = a.split("$"), b.split("$")
    if rule == "paper":
        return ta[1] == tb[0]                          # head of a is tail of b
    return ta[1] == tb[0] or ta[0] == tb[1]            # the shipped code's widened rule


def stg(p, arm, teleport=TELEPORT):
    """Line-graph slots -> row-stochastic transition matrix of the shared topology graph.

    arm 'paper'   co-presence proportion weight, head-to-tail arcs, PageRank teleport.
    arm 'shipped' ((n+1) - co-presence) weight, widened arcs, no teleport (a flat 1e-10 is
                  added by the caller's cost step, as Accordion.compute_cost does).

    Returns (P, index) with P row-stochastic and index == list(p.index), so the cost matrix it
    produces is aligned to the marginals by construction. This is the one place where the
    `paper` arm deliberately does not reproduce the shipped code -- see the header, defect (4).
    """
    idx = list(p.index)
    n_s = p.shape[1]
    B = (p.to_numpy() > PSEUDO).astype(np.int64)
    CO = B @ B.T                                       # co-presence counts between slots
    W = np.zeros((len(idx), len(idx)))
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            if a == b or not _composable(i, j, "paper" if arm == "paper" else "code"):
                continue
            W[a, b] = CO[a, b] / n_s if arm == "paper" else ((n_s + 1) - CO[a, b]) / n_s
    rs = W.sum(1, keepdims=True)
    rs[rs == 0] = 1.0                                  # a slot with no composable partner
    P = W / rs
    if arm == "paper":
        # The paper's low-rank PageRank term. Without it the head-to-tail-only line graph is not
        # strongly connected, the chain is reducible, and -log of the hitting probability is not
        # finite -- measured, 12 entries on this cohort.
        P = (1.0 - teleport) * P + teleport / len(idx)
    else:
        P = P + 1e-10                                  # Accordion.compute_cost's `d` default
    return P, idx


def htd_cost(P, beta=HTD_BETA):
    """Hitting Time Distance cost from a transition matrix, via THEIR getCTD.

    Deliberately calls the vendored implementation rather than reimplementing Boyd et al.: the
    whole point of this comparison is that the published method is run, not approximated. The
    known hazard inside it -- `-np.log10(Aht, where=Aht != 0)` has no `out=`, so entries where
    the symmetrised matrix is exactly zero are uninitialised memory -- is not patched here. It
    is instead checked for by the caller (06_distance/03 SELF-CHECK 2 rejects any |C| > 1e6),
    so that a defect in the published code surfaces as a failed check rather than as a number.
    """
    scaccordion_on_path()
    import importlib.util as iu
    spec = iu.spec_from_file_location(
        "scacc_distances", os.path.join(SCACC_REPO, "scaccordion", "tools", "distances.py"))
    mod = iu.module_from_spec(spec)
    spec.loader.exec_module(mod)                       # bypasses the package __init__ chain
    C = np.asarray(mod.getCTD(np.asarray(P, float), beta=beta), float)
    C = (C + C.T) / 2.0                                # getCTD already symmetrises; idempotent
    np.fill_diagonal(C, 0.0)
    return C


# ----------------------------------------------------------------------------------------------
# 3. the transport itself
# ----------------------------------------------------------------------------------------------
def emd_matrix(p, C):
    """Balanced OT distance between every pair of samples, with ONE cohort-level cost matrix.

    This is the line that separates scACCorDiON from this project's own FGW pipeline, and the
    reason it was worth running at all. C is the same matrix for every pair, so slot i means the
    same interaction in both samples and the linear program only moves mass between labelled
    bins. FGW instead solves for a coupling between two per-sample cost matrices, which is a
    node-matching problem -- and Gromov-Wasserstein is by definition invariant to relabelling
    the nodes, so it cannot use the fact that the seven bins are already matched. See
    FINDINGS_topology_null.md section 13.
    """
    import ot
    A = p.to_numpy(float)
    tot = A.sum(0, keepdims=True)
    if (tot <= 0).any():
        raise ValueError("a sample has zero total edge weight; balanced OT is undefined for it")
    A = A / tot
    C = np.ascontiguousarray(np.asarray(C, float))
    n = A.shape[1]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = ot.emd2(A[:, i].copy(), A[:, j].copy(), C)
    return D


# ----------------------------------------------------------------------------------------------
# 4. the positive control
# ----------------------------------------------------------------------------------------------
def plant(tbls, edges, delta, samples):
    """Multiply the ligand-receptor signal on `edges` by (1 + delta) in `samples`.

    Planting happens at the LIGAND-RECEPTOR level, before any aggregation, so the entire
    published construction downstream -- groupby-sum, variance filter, shared topology graph,
    hitting-time cost, transport -- is rebuilt from the perturbed input exactly as it would be
    from real data. Planting the aggregated edge weight instead would bypass the variance
    filter and the co-presence geometry, which are two of the three things being tested.

    edges    iterable of (source, target) pairs.
    delta    0 gives back an unmodified copy, which is the blank control arm.
    samples  the subset that receives the effect.

    A multiplicative effect is used, not additive, because an additive effect on a slot a sample
    never expressed would also change that slot's PRESENCE, and presence is what drives most of
    this distance (measured: Spearman 0.86 against support Jaccard). A multiplicative effect
    leaves the support pattern untouched, so the control tests strength, not detection.
    """
    want = set(samples)
    pairs = {(str(a), str(b)) for a, b in edges}
    out = {}
    for k, v in tbls.items():
        v = v.copy()
        if k in want and len(v) and delta != 0:
            hit = [(s, t) in pairs for s, t in zip(v["source"], v["target"])]
            v.loc[hit, "lr_means"] = v.loc[hit, "lr_means"] * (1.0 + delta)
        out[k] = v
    return out
