#!/usr/bin/env python
# 24_mantel_ceiling.py ----
# INPUT  : leukemia_ml/05_CCC_per_sample/{cohort_manifest.csv, sample_disease_corrected.csv,
#            02_cellchat/*__mat.csv, 02_cellchat__half_{A,B}/*__mat.csv}
#          aml_niche_net/results/tables/05_ccc/{ccc_node_presence.csv, ccc_sample_manifest.csv,
#            tensors/*/*__ccc_cellchat.csv}; 01_preprocess/00_curated_manifest.csv
# OUTPUT : aml_niche_net/results/tables/08_scoring/mantel_ceiling.csv
# WHAT IT DOES : the three numbers the 2026-09-30 prior-art verdict named as the cheapest decisive
#          test, before any four-arm benchmark is attempted:
#            (1) Mantel rho between D_comp and D_ccc -- are composition and communication the same
#                measurement? D_comp is CLR/Aitchison (Halter's baseline, not raw proportions, so
#                the baseline cannot be called a straw man); D_ccc is cosine on the 7x7 CellChat
#                weight vector with population.size = FALSE.
#            (2) The RELIABILITY CEILING of each matrix, and the disattenuated Mantel
#                rho / sqrt(rel_comp * rel_ccc). This is the only number that answers the
#                "your null is a power problem" objection.
#            (3) Whether D_ccc still tracks within-AML genotype once D_comp is partialled out
#                (partial Mantel, within-dataset pairs only, so batch cannot carry it).
#            (4) THE ALTERNATIVE EXPLANATION, tested rather than assumed: is the part of D_ccc
#                that is independent of composition simply BATCH? Partial Mantel of D_ccc against
#                a dataset-mismatch matrix given D_comp, plus the comp-vs-ccc Mantel recomputed on
#                within-dataset pairs only.
#          TWO RELIABILITY CAVEATS, both load-bearing and both stated rather than hidden:
#            - composition reliability CANNOT come from the existing split halves: that split is
#              stratified by compartment, which forces both halves to the same composition by
#              construction (06.4 measured retest 5.8e-4). It is estimated here instead from
#              N_SPLIT unstratified random half-splits of each sample's cells, simulated exactly
#              by hypergeometric sampling from the per-compartment counts.
#            - CCC reliability DOES come from the stratified halves, which hold composition fixed,
#              so it is an OPTIMISTIC estimate for the molecular term (it excludes composition
#              sampling noise). The asymmetry is interpretable in both directions: if the observed
#              Mantel rho is far below the ceiling even with CCC flattered, the two are genuinely
#              different measurements; if it is at the ceiling, the true ceiling is lower still.
# Usage : general_env/bin/python scripts/08_scoring/24_mantel_ceiling.py
import os, glob, json
import numpy as np, pandas as pd
from scipy import stats

LM  = "/FAST/gr10634/gaozy/leukemia_ml/05_CCC_per_sample"
AN  = "/FAST/gr10634/gaozy/aml_niche_net/results/tables"
OUT = os.path.join(AN, "08_scoring", "mantel_ceiling.csv")
N_PERM, N_SPLIT = 999, 20
rng = np.random.default_rng(260930)
LM_CP = ["Blast_HSPC","Mono_DC","Granulocyte","Ery_Meg","T_NK","B_Plasma","Stroma_niche"]
AN_CP = ["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
fn = lambda s: "".join(c if c.isalnum() or c in "._-" else "_" for c in s)
UT = lambda D: D[np.triu_indices(len(D), 1)]


def clr_d(counts):
    """Aitchison distance on CLR-transformed compositions. +0.5 pseudocount for structural zeros."""
    p = counts + 0.5
    p = p / p.sum(1, keepdims=True)
    L = np.log(p)
    C = L - L.mean(1, keepdims=True)
    return np.sqrt(np.maximum((C**2).sum(1)[:,None] + (C**2).sum(1)[None,:] - 2*C@C.T, 0))


def cos_d(X):
    n = np.sqrt((X**2).sum(1, keepdims=True)); n[n == 0] = 1.0
    Y = X/n; D = 1.0 - Y@Y.T
    D = (D+D.T)/2; np.fill_diagonal(D, 0.0); return D


def mantel(D1, D2, nperm=N_PERM):
    a, b = UT(D1), UT(D2)
    r = stats.spearmanr(a, b).statistic
    n = len(D1); null = np.empty(nperm)
    for i in range(nperm):
        q = rng.permutation(n)
        null[i] = stats.spearmanr(a, UT(D2[np.ix_(q, q)])).statistic
    return float(r), float((np.abs(null) >= abs(r)).mean())


def partial_mantel(Dx, Dy, Dz, mask, nperm=N_PERM):
    """Spearman partial correlation of Dx,Dy given Dz on the masked pairs, permutation p."""
    iu = np.triu_indices(len(Dx), 1)
    m = mask[iu]
    x, y, z = [stats.rankdata(M[iu][m]) for M in (Dx, Dy, Dz)]
    def pr(xx):
        rx = xx - np.poly1d(np.polyfit(z, xx, 1))(z)
        ry = y - np.poly1d(np.polyfit(z, y, 1))(z)
        return stats.spearmanr(rx, ry).statistic
    r = pr(x); n = len(Dx)
    null = np.empty(nperm)
    for i in range(nperm):
        q = rng.permutation(n)
        null[i] = pr(stats.rankdata(Dx[np.ix_(q, q)][iu][m]))
    return float(r), float((np.abs(null) >= abs(r)).mean())


rows = []
# =============================== COHORT 1: leukemia_ml (67) ====================================
mf = pd.read_csv(f"{LM}/cohort_manifest.csv").query("included == True").sort_values("Sample_Name")
S = list(mf.Sample_Name)
cnt = mf[[f"n_{c}" for c in LM_CP]].to_numpy(float)
Dc = clr_d(cnt)

grid = [(a, b) for a in LM_CP for b in LM_CP]
def ccc_vec(d):
    V = np.zeros((len(S), 49))
    for i, s in enumerate(S):
        t = pd.read_csv(f"{d}/{fn(s)}__mat.csv")
        k = {(r.source, r.target): r.weight for r in t.itertuples()}
        V[i] = [k.get(g, 0.0) for g in grid]
    return V
Vp = ccc_vec(f"{LM}/02_cellchat"); Dk = cos_d(Vp)
r_obs, p_obs = mantel(Dc, Dk)
print(f"[1] leukemia_ml n={len(S)}  Mantel rho(D_comp, D_ccc) = {r_obs:+.3f}  (perm p {p_obs:.3f})")

# reliability: CCC from the stratified halves (optimistic), composition from unstratified splits
Dk_A, Dk_B = cos_d(ccc_vec(f"{LM}/02_cellchat__half_A")), cos_d(ccc_vec(f"{LM}/02_cellchat__half_B"))
rel_ccc = stats.spearmanr(UT(Dk_A), UT(Dk_B)).statistic
rel_c = []
for _ in range(N_SPLIT):
    N = cnt.sum(1); half = np.floor(N/2).astype(int)
    A = np.array([rng.multivariate_hypergeometric(cnt[i].astype(int), half[i]) for i in range(len(S))], float)
    rel_c.append(stats.spearmanr(UT(clr_d(A)), UT(clr_d(cnt - A))).statistic)
rel_comp = float(np.mean(rel_c))
ceil = np.sqrt(max(rel_comp, 0) * max(rel_ccc, 0))
print(f"    reliability: comp {rel_comp:.3f} (unstratified, {N_SPLIT} splits) | ccc {rel_ccc:.3f} (stratified halves, optimistic)")
print(f"    ceiling sqrt(rel*rel) = {ceil:.3f}  ->  disattenuated Mantel = {r_obs/ceil:+.3f}")

# genotype, within dataset
lab = pd.read_csv(f"{LM}/sample_disease_corrected.csv").set_index("Sample_Name").reindex(S)
mal = (lab.disease_class == "malignant") & lab.driver_mutations.notna() & (lab.driver_mutations != "")
npm1 = mal & lab.driver_mutations.fillna("").str.contains("NPM1") & ~lab.driver_mutations.fillna("").str.contains("NPM1_WT")
ds = lab.Dataset_ID.to_numpy()
G = (npm1.to_numpy()[:,None] != npm1.to_numpy()[None,:]).astype(float)
mask = (ds[:,None] == ds[None,:]) & mal.to_numpy()[:,None] & mal.to_numpy()[None,:]
print(f"    genotype: {int(npm1.sum())} NPM1-mut of {int(mal.sum())} with driver calls | within-dataset malignant pairs {int(np.triu(mask,1).sum())}")
if np.triu(mask, 1).sum() >= 20:
    r_g_ccc, p_g_ccc = partial_mantel(Dk, G, Dc, mask)
    r_g_c, p_g_c = partial_mantel(Dc, G, Dk, mask)
    print(f"    partial Mantel  D_ccc~genotype | D_comp : {r_g_ccc:+.3f} (p {p_g_ccc:.3f})")
    print(f"                    D_comp~genotype | D_ccc : {r_g_c:+.3f} (p {p_g_c:.3f})")
else:
    r_g_ccc = p_g_ccc = r_g_c = p_g_c = np.nan
rows.append(dict(cohort="leukemia_ml", n=len(S), mantel_rho=r_obs, mantel_p=p_obs,
                 rel_comp=rel_comp, rel_ccc=rel_ccc, ceiling=ceil, mantel_disattenuated=r_obs/ceil,
                 partial_ccc_genotype=r_g_ccc, partial_ccc_genotype_p=p_g_ccc,
                 partial_comp_genotype=r_g_c, partial_comp_genotype_p=p_g_c))

# =============================== COHORT 2: aml_niche_net (138) =================================
pr7 = pd.read_csv(f"{AN}/06_distance/scaccordion_distance__prop7.csv", index_col=0)
S2 = list(pr7.index)
pres = pd.read_csv(f"{AN}/05_ccc/ccc_node_presence.csv")
c2 = (pres[pres.hierarchy_bin.isin(AN_CP)]
      .pivot_table(index="sample", columns="hierarchy_bin", values="n_cells", aggfunc="sum")
      .reindex(index=S2, columns=AN_CP).fillna(0.0).to_numpy(float))
Dc2 = clr_d(c2)
bi = {b: i for i, b in enumerate(AN_CP)}
V2 = np.zeros((len(S2), 49))
for i, s in enumerate(S2):
    h = glob.glob(f"{AN}/05_ccc/tensors/*/{s}__ccc_cellchat.csv")
    if not h: continue
    t = pd.read_csv(h[0]); t = t[(t.pval <= 0.05) & t.sender_bin.isin(AN_CP) & t.receiver_bin.isin(AN_CP)]
    for (a, b), g in t.groupby(["sender_bin", "receiver_bin"]):
        V2[i, bi[a]*7 + bi[b]] = g.prob.sum()
Dk2 = cos_d(V2)
r2, p2 = mantel(Dc2, Dk2)
print(f"\n[2] aml_niche_net n={len(S2)}  Mantel rho(D_comp, D_ccc) = {r2:+.3f}  (perm p {p2:.3f})")
nf = pd.read_csv(f"{AN}/08_scoring/scaccordion_noise_floor.csv")
print(f"    no half tensors on disk here; published SNR_between for edge arms: "
      f"{', '.join('%s %.1f' % (r.arm, r.snr_between) for r in nf.itertuples())}")
man2 = pd.read_csv(f"{AN}/05_ccc/ccc_sample_manifest.csv").set_index("sample").reindex(S2)
cur = pd.read_csv(f"{AN}/01_preprocess/00_curated_manifest.csv")
dm = cur.set_index("sample").reindex(S2).driver_mutations.fillna("")
has = dm != ""
n1 = has & dm.str.contains("NPM1") & ~dm.str.contains("NPM1_WT")
ds2 = man2.dataset.to_numpy()
G2 = (n1.to_numpy()[:,None] != n1.to_numpy()[None,:]).astype(float)
m2 = (ds2[:,None] == ds2[None,:]) & has.to_numpy()[:,None] & has.to_numpy()[None,:]
print(f"    genotype: {int(n1.sum())} NPM1-mut of {int(has.sum())} with calls | within-dataset pairs {int(np.triu(m2,1).sum())}")
if np.triu(m2, 1).sum() >= 20:
    rg, pg = partial_mantel(Dk2, G2, Dc2, m2)
    rc, pc = partial_mantel(Dc2, G2, Dk2, m2)
    print(f"    partial Mantel  D_ccc~genotype | D_comp : {rg:+.3f} (p {pg:.3f})")
    print(f"                    D_comp~genotype | D_ccc : {rc:+.3f} (p {pc:.3f})")
else:
    rg = pg = rc = pc = np.nan
rows.append(dict(cohort="aml_niche_net", n=len(S2), mantel_rho=r2, mantel_p=p2,
                 rel_comp=np.nan, rel_ccc=np.nan, ceiling=np.nan, mantel_disattenuated=np.nan,
                 partial_ccc_genotype=rg, partial_ccc_genotype_p=pg,
                 partial_comp_genotype=rc, partial_comp_genotype_p=pc))

# =============================== (4) is the independent part batch? ============================
print("\n[4] is D_ccc's composition-independent variation just batch?")
for name, Dc_, Dk_, dsv in [("leukemia_ml", Dc, Dk, ds), ("aml_niche_net", Dc2, Dk2, ds2)]:
    B = (dsv[:, None] != dsv[None, :]).astype(float)
    Wm = (dsv[:, None] == dsv[None, :])
    rb, pb = partial_mantel(Dk_, B, Dc_, np.ones_like(Wm))
    rcb, pcb = partial_mantel(Dc_, B, Dk_, np.ones_like(Wm))
    iu = np.triu_indices(len(Dc_), 1); m = Wm[iu]
    rw = stats.spearmanr(Dc_[iu][m], Dk_[iu][m]).statistic
    print("    %-14s D_ccc~batch|D_comp %+.3f (p %.3f) | D_comp~batch|D_ccc %+.3f (p %.3f) | "
          "comp~ccc within-dataset %+.3f" % (name, rb, pb, rcb, pcb, rw))
    for r in rows:
        if r["cohort"] == name:
            r.update(partial_ccc_batch=rb, partial_ccc_batch_p=pb,
                     partial_comp_batch=rcb, mantel_within_dataset=rw,
                     mantel_within_disattenuated=(rw / r["ceiling"] if r["ceiling"] == r["ceiling"] else np.nan))

pd.DataFrame(rows).assign(kind_of_evidence="decisive_cheap_test").to_csv(OUT, index=False)

## SELF-CHECK: the CLR composition distance must track the conventional TV composition distance ---
p7 = c2 / c2.sum(1, keepdims=True)
TV = 0.5*np.abs(p7[:,None,:] - p7[None,:,:]).sum(2)
rs = stats.spearmanr(UT(Dc2), UT(TV)).statistic
print(f"\n[SELF-CHECK] CLR/Aitchison vs TV composition distance: rho {rs:.3f} (must be > 0.7)")
assert rs > 0.7, "CLR composition distance does not track TV -- construction bug"
print(f"[done] wrote {OUT}")
