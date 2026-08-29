#!/usr/bin/env python
# 07_feature_decomposition.py ----
# Phase 7 (stage 08_scoring). DIAGNOSTIC 1, following the alpha sweep.
#
# CONTEXT: 06_alpha_sweep showed a perfectly monotone dose-response -- separation scales with the FEATURE
# weight and vanishes at alpha=1 (pure topology: healthy 0.048 vs AML 0.052, p=0.271). So whatever signal
# HDS carries lives in the node-FEATURE term, not the edge structure. That raises a circularity worry:
# one of the three features is frac_malignant, which we FORCE to 0 for healthy samples by design
# (FGW_ZERO_HEALTHY_MAL). If that single feature drives everything, the significance is an artifact of our
# own encoding and cannot be reported as a finding at all.
#
# THIS SCRIPT re-runs the whole HDS pipeline (LOO healthy barycenter + scoring + the two regression models)
# separately for each FEATURE SUBSET:
#   all3          : frac_malignant, mean_stemness, n_cells      (the current definition)
#   no_frac_mal   : mean_stemness, n_cells                      <- KEY: signal here is NOT circular
#   only_frac_mal : frac_malignant                              <- KEY: isolates the encoded label
#   only_stemness : mean_stemness                               (transcriptional, CNV-independent)
#   only_ncells   : n_cells                                     (pure composition; already covered by blast_proxy)
#
# READING THE RESULT:
#   only_frac_mal strong AND no_frac_mal null  -> CIRCULAR. The result is an artifact of our encoding.
#   no_frac_mal still significant              -> a real (if non-topological) feature signal survives;
#                                                 then check WHICH: stemness (biology) vs n_cells (composition).
#   only_ncells carries it                     -> trivially compositional, and blast_proxy already covers it.
# Run at alpha=0 (pure feature, cleanest attribution) and alpha=0.5 (the default setting).
#
# INPUT  : <root>/07_fgw/fgw_{nodes,edges,input_index}_long.csv
# OUTPUT : <root>/08_scoring/feature_decomposition.csv          (--split all)
#          <root>/08_scoring/feature_decomposition__<split>.csv (otherwise; override with --out)
# Usage  : python scripts/08_scoring/07_feature_decomposition.py [--n_perm 10000]
import argparse, hashlib, os, warnings
import sys
import numpy as np
import pandas as pd
import ot

# The timepoint vocabulary is LOADED, not literal. A hard-coded set here silently deleted 34 of
# 214 samples (64% of the treated arm) after CANONICAL_TIMEPOINTS changed on 2026-08-04.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'config'))
from fgw_vocab import load_vocab, load_features, assert_index_covered

# SUPPRESS THE POT alpha=0 DIVIDE-BY-ZERO AT THE CALL SITE, NOT PROCESS-WIDE. The previous line here
# was a module-level warnings.filterwarnings("ignore", RuntimeWarning), which silenced the expected
# ot/gromov/_gw.py log-term 0/0 at alpha=0 AND every other RuntimeWarning in the process -- overflow,
# invalid value in our own arithmetic, a numpy division we did not intend. A blanket filter is how an
# honest "expected warning" note turns into a blindfold. This context manager covers only the two POT
# entry points that actually raise it; anything else still reaches the log, where the chain's
# gate_no_warnings can act on it.
from contextlib import contextmanager
@contextmanager
def pot_quiet():
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="divide by zero encountered")
        warnings.filterwarnings("ignore", category=RuntimeWarning,
                                message="invalid value encountered")
        yield

FGW_NODES=["HSC_MPP","LMPP_GMP","Mono_DC","Erythroid","Megakaryocyte","T_NK","B_Plasma"]
ALL_FEATURES=["frac_malignant","mean_stemness","n_cells"]
BLAST_BINS=["HSC_MPP","LMPP_GMP","Mono_DC"]
AML_TP = None  # set from fgw_vocab.json below -- see scripts/config/fgw_vocab.py
DEFAULT_ROOT="/FAST/gr10634/gaozy/aml_niche_net/results/tables"; SEED=491638
# A node is PRESENT if its normalised mass exceeds this. FGW_EPS_MASS is 1e-6 before the within-sample
# renormalisation, so an absent node lands near 1e-7 and a present one is orders of magnitude above.
_PRESENT_EPS=1e-5

# CORE: the three features actually in FGW_FEATURES. This is the original circularity test.
CORE_SUBSETS={
 "all3":          ["frac_malignant","mean_stemness","n_cells"],
 "no_frac_mal":   ["mean_stemness","n_cells"],
 "only_frac_mal": ["frac_malignant"],
 "only_stemness": ["mean_stemness"],
 "only_ncells":   ["n_cells"],
}

# CANDIDATES: emitted by 07_fgw/01 as FGW_CANDIDATE_FEATURES but NOT part of the FGW distance. A feature
# earns its way into FGW_FEATURES by surviving here first.
#
# WHY EACH ONE EXISTS. The alpha sweep put the entire H2 signal in the feature term (alpha=1 pure
# topology: within-dataset p=0.966), and only_stemness survived the circularity test (p=0.0063) while
# only_ncells did not (p=0.410). So the live question is WHICH stemness -- and whether a malignancy axis
# we did NOT encode by hand behaves like the one we did.
#
#   only_vg_mal          transcriptional malignancy (van Galen), NEVER zeroed for healthy. If this
#                        separates, malignancy is a real signal and not just our encoding. This is the
#                        honest counterpart to only_frac_mal.
#   only_cnv_burden      continuous CNV instead of the thresholded call
#   only_stem_normal     <- DECISIVE. Stemness of NON-malignant cells. If this carries the signal, the
#                        finding is about the microenvironment.
#   only_stem_malignant  <- DECISIVE. Stemness of malignant cells. If the signal is only here, the
#                        finding reduces to "blasts are stem-like" -- known, and near-circular with
#                        malignancy itself.
#   only_cyto_normal / only_cyto_malignant   the same split under CytoTRACE2, which never saw LSC17.
#                        Agreement across two independent potency estimates is the reproducibility arm.
#   noncircular          everything we did NOT encode by hand, together.
CANDIDATE_SUBSETS={
 "only_vg_mal":         ["frac_malignant_vg"],
 "only_cnv_burden":     ["mean_cnv_burden"],
 "only_stem_normal":    ["mean_stemness_normal"],
 "only_stem_malignant": ["mean_stemness_malignant"],
 "only_cyto_normal":    ["mean_cytotrace_normal"],
 "only_cyto_malignant": ["mean_cytotrace_malignant"],
 "noncircular":         ["frac_malignant_vg","mean_stemness_normal"],
}

# PANEL_SUBSETS is built from the vocabulary at run time, not listed here: 05_ccc/03 decides which
# scores exist and config_fgw.R derives the candidate names from the same CCC_PANELS declaration.
# Re-listing ~114 names in a third place is how the timepoint bug happened (see fgw_vocab.py).
def panel_subsets(cand):
    """One single-feature subset per panel candidate, tagged with its family for FDR."""
    out = {}
    for c in cand:
        if c in ("frac_malignant_vg","mean_cnv_burden","mean_stemness_normal",
                 "mean_stemness_malignant","mean_cytotrace_normal","mean_cytotrace_malignant"):
            continue
        out["only_" + c] = [c]
    return out

def family_of(name, fam_map):
    """Family key for BH correction, LOOKED UP in the vocabulary rather than pattern-matched.

    This used to test the name's two-letter prefix against a hard-coded ("st","pg","cs","mt","pt").
    That is a silent-failure shape: adding a sixth family to CCC_PANELS -- the `mp` family P2 calls
    for -- would not have raised anything. Every mp_* feature would have been labelled "core", the
    BH loop skips "core", and the whole family would have been reported with q = NaN, i.e.
    uncorrected, while looking exactly like the other families in the output.

    fam_map comes from fgw_vocab.json, which 07_fgw/01 writes from config_ccc.R's CCC_PANELS. An
    unrecognised candidate is a hard error: the alternative is to guess, and guessing here means
    silently skipping a multiple-testing correction.
    """
    f = name[len("only_"):] if name.startswith("only_") else name
    if f in fam_map:
        return fam_map[f]
    if name in CORE_SUBSETS or name in CANDIDATE_SUBSETS or name in ("all3","no_frac_mal","noncircular"):
        return "core"
    raise SystemExit(
        f"cannot assign a multiple-testing family to '{name}'.\n"
        f"  It is not in fgw_vocab.json's candidate_family map ({len(fam_map)} entries) and is not a\n"
        f"  core/candidate subset. Re-run scripts/07_fgw/01_build_fgw_inputs.R so the vocabulary\n"
        f"  carries the family for every candidate CCC_PANELS declares.")

ap=argparse.ArgumentParser()
ap.add_argument("--root",default=DEFAULT_ROOT)
ap.add_argument("--n_perm",type=int,default=10000)
ap.add_argument("--alphas",default="0,0.5")
ap.add_argument("--split",default="all",choices=["all","discovery","validation"],
    help="which arm of the dataset-level 70/30 split to use. 'all' reproduces every analysis before "
         "2026-08-26 and is the right choice ONLY for a hypothesis fixed in advance. A sweep over the "
         "panel families is selection and must run on 'discovery'.")
ap.add_argument("--split_csv",default=None,
    help="path to 01_preprocess/02_sample_split.csv; defaults to <root>/01_preprocess/02_sample_split.csv")
ap.add_argument("--subsets",default="all",choices=["core","candidates","panels","all"],
                help="core = the 3 features in FGW_FEATURES; candidates = FGW_CANDIDATE_FEATURES")
ap.add_argument("--only",default="",
    help="comma-separated subset names to run, e.g. 'only_pt_predicted_Pseudotime'. Applied AFTER "
         "--subsets builds the candidate set, so --subsets panels --only <two names> tests exactly "
         "those two. This exists for CONFIRMATION runs: the pre-registration says nothing that "
         "failed screening may be looked at again on Validation, and without this flag the only way "
         "to test a screened hit was --subsets panels, which computes, writes AND prints all 114 -- "
         "including the 107 that failed. An unknown name is a hard error, never a silent no-op.")
ap.add_argument("--out",default=None,
    help="output basename under <root>/08_scoring. Defaults to feature_decomposition.csv for "
         "--split all, and feature_decomposition__<split>.csv otherwise -- see the note at the "
         "write site for why the split MUST be in the path and not only in a column.")
ap.add_argument("--covar",default="none",choices=["none","depth","ribo","both"],
    help="technical covariate added to BOTH model matrices. AML libraries are systematically deeper "
         "than healthy (med_ncount 4874 vs 3006, AUC 0.695), and on the Discovery arm depth still "
         "tracks the label at r=+0.284 (p=0.012) AFTER the dataset fixed effects. blast_proxy does "
         "not absorb it. Default none reproduces the frozen screen bit-for-bit; see the 2026-08-29 "
         "amendment in PREREGISTRATION_panel_screen.md for the rule fixed before this was run.")
ap.add_argument("--qc_csv",default=None,
    help="path to 01_preprocess/03_qc_report__ALL.csv for --covar; defaults to <root>/01_preprocess/")
ap.add_argument("--absent_mask",action="store_true",
    help="normalise the feature cost matrix by its maximum over PRESENT nodes only. Absent nodes "
         "carry ~1e-6 mass, transport nothing, and hold imputed features -- yet a phantom row can "
         "set the scale for the whole sample. Changes HDS by up to 12.8%% on 5 of 138 samples and "
         "hits healthy samples (52.2%%) more often than AML (35.7%%), so the leak points toward the "
         "hypothesis. Off by default so the baseline arm stays a reproduction.")
ap.add_argument("--max_iter",type=int,default=1000)
args=ap.parse_args()
# ONE RNG PER (alpha, subset, model), KEYED ON A STABLE LABEL. A single shared stream made a subset's
# permutation p depend on how many subsets ran before it: --subsets panels and --subsets all disagreed
# on 113 of 114 panel p_strat (median |dp|=0.0030, max 0.0135) while beta was bit-identical, and that
# flipped 3 BH-FDR verdicts on 2026-08-26. Keying on the LOOP INDEX would fix the symptom and
# reintroduce the bug the moment SUBSETS is reordered or a subset is inserted ahead of another; keying
# on the label survives both. Do NOT "fix" this by reseeding one rng inside the loop -- that makes
# every subset share the same permutations, which is a worse bug (see test T3 below).
def _lab(s): return int.from_bytes(hashlib.blake2b(s.encode(),digest_size=8).digest(),"big")
def _rng(alpha,name,model): return np.random.default_rng([SEED,_lab(f"{alpha:g}|{name}|{model}")])
ALPHAS=[float(x) for x in args.alphas.split(",")]
D_FGW=os.path.join(args.root,"07_fgw"); D_OUT=os.path.join(args.root,"08_scoring"); os.makedirs(D_OUT,exist_ok=True)
_VOCAB = load_vocab(D_FGW)
# The subset definitions below are DELIBERATELY literal: "all3" names a specific model we compare
# candidates against. But if the inputs were built with a different feature set, that name lies --
# so say so loudly rather than silently reporting "all3" for a model that is not all of anything.
_IN_MODEL = load_features(D_FGW)
if sorted(_IN_MODEL) != sorted(["frac_malignant", "mean_stemness", "n_cells"]):
    print("[!] these inputs were built with features %s, not the production triple." % _IN_MODEL)
    print("    The subset named 'all3' below is therefore NOT the model that produced this root's")
    print("    alpha sweep. Read the subset rows as standalone feature tests, not as a decomposition")
    print("    of the fitted model.")
AML_TP = _VOCAB["aml_timepoints"]

SUBSETS = dict(CORE_SUBSETS) if args.subsets in ("core","all") else {}
if args.subsets in ("candidates","all"): SUBSETS.update(CANDIDATE_SUBSETS)


edges=pd.read_csv(os.path.join(D_FGW,"fgw_edges_long.csv"))
nodes=pd.read_csv(os.path.join(D_FGW,"fgw_nodes_long.csv"))

# Panel subsets are added HERE, after nodes is read, because they are defined by what the vocabulary
# lists AND what the file actually carries -- a candidate named in the vocab but absent from the
# table would otherwise become a hard failure in the column check below rather than a skipped test.
if args.subsets in ("panels","all"):
    _cand = [c for c in _VOCAB.get("candidate_features", []) if c in nodes.columns]
    _absent = [c for c in _VOCAB.get("candidate_features", []) if c not in nodes.columns]
    SUBSETS.update(panel_subsets(_cand))
    if _absent:
        print(f"[0] {len(_absent)} vocab candidate(s) absent from fgw_nodes_long.csv, skipped: "
              f"{', '.join(_absent[:6])}{' ...' if len(_absent) > 6 else ''}")

# --only is applied HERE, after every source of subsets has contributed, so a name from any of them
# resolves. Unknown names abort: a confirmation run that silently tested nothing, or tested a
# mistyped neighbour, is the worst possible outcome for an arm that can only be spent once.
if args.only:
    _want_only = [x.strip() for x in args.only.split(",") if x.strip()]
    _unknown = [x for x in _want_only if x not in SUBSETS]
    if _unknown:
        raise SystemExit(
            f"--only names subset(s) that do not exist: {', '.join(_unknown)}\n"
            f"  Available under --subsets {args.subsets}: {len(SUBSETS)} names, e.g. "
            f"{', '.join(sorted(SUBSETS)[:6])} ...\n"
            f"  (A hit from the panel screen needs --subsets panels or all.)")
    SUBSETS = {k: SUBSETS[k] for k in _want_only}
    print(f"[0] --only: restricted to {len(SUBSETS)} subset(s): {', '.join(_want_only)}")
idx=pd.read_csv(os.path.join(D_FGW,"fgw_input_index.csv"))

# Fail loudly on a label the vocabulary does not cover. Without this an unrecognised
# timepoint is not an error -- it is a sample that quietly stops being AML and stops
# being healthy, and therefore stops existing for every test below.
assert_index_covered(idx, _VOCAB)

# A subset naming a column that fgw_nodes_long.csv does not carry must abort, not silently produce a
# zero-variance feature. Candidate columns only appear once 07_fgw/01 has been re-run with
# FGW_CANDIDATE_FEATURES set, and a missing column would otherwise read as "this feature has no signal".
_want = sorted({c for f in SUBSETS.values() for c in f})
_missing = [c for c in _want if c not in nodes.columns]
if _missing:
    raise SystemExit(
        f"fgw_nodes_long.csv lacks column(s): {', '.join(_missing)}\n"
        f"  Re-run scripts/07_fgw/01_build_fgw_inputs.R (it emits FGW_CANDIDATE_FEATURES).\n"
        f"  Columns present: {', '.join(nodes.columns)}")
_const = [c for c in _want if float(nodes[c].std(skipna=True) or 0.0) == 0.0]
if _const:
    raise SystemExit(f"column(s) constant across all nodes, cannot be tested: {', '.join(_const)}")
print(f"[0] subsets={args.subsets} -> {len(SUBSETS)} sets over {len(_want)} distinct feature columns")

_cache={}
def build_one(ds,smp,feats):
    k=(ds,smp,tuple(feats))
    if k in _cache: return _cache[k]
    e=edges[(edges["dataset"]==ds)&(edges["sample"]==smp)]
    C=e.pivot(index="sender_bin",columns="receiver_bin",values="C").reindex(FGW_NODES,columns=FGW_NODES).to_numpy(float)
    C=np.nan_to_num(C,nan=1.0)
    nd=nodes[(nodes["dataset"]==ds)&(nodes["sample"]==smp)].set_index("hierarchy_bin").reindex(FGW_NODES)
    F=np.nan_to_num(nd[feats].to_numpy(float),nan=0.0)
    if F.ndim==1: F=F.reshape(-1,1)
    p=np.nan_to_num(nd["mass"].to_numpy(float),nan=1e-6); p=p/p.sum()
    _cache[k]=(C,F,p); return _cache[k]

def barycenter(keys,alpha,feats):
    Cs=[];Fs=[];ps=[]
    for k in keys:
        C,F,p=build_one(*k,feats); Cs.append(C);Fs.append(F);ps.append(p)
    m=len(Cs); n=7
    with pot_quiet():
        out=ot.gromov.fgw_barycenters(n,Fs,Cs,ps,lambdas=[1.0/m]*m,alpha=alpha,loss_fun="square_loss",
            symmetric=False,max_iter=args.max_iter,p=np.ones(n)/n,init_C=np.mean(np.stack(Cs),0),
            init_X=np.mean(np.stack(Fs),0),random_state=SEED,log=True)
    return np.asarray(out[1]),np.asarray(out[0]),np.ones(n)/n

def fgw2(C,F,p,Cb,Fb,pb,alpha):
    M=ot.dist(F,Fb)
    # THE SCALE MUST COME FROM NODES THAT TRANSPORT. An absent node carries FGW_EPS_MASS (1e-6 before
    # renormalisation) and an imputed feature vector, so its row of M is a distance to something that
    # is not there -- yet under a plain M.max() it can set the normalisation for the entire sample.
    mx = M[p > _PRESENT_EPS].max() if (args.absent_mask and (p > _PRESENT_EPS).any()) else M.max()
    M=M/(mx+1e-9) if mx>0 else M
    with pot_quiet():
        return float(ot.gromov.fused_gromov_wasserstein2(M,C,Cb,p,pb,loss_fun="square_loss",alpha=alpha,symmetric=False))

def fwl_perm(y,X0,a,groups,n_perm,rng):
    Q,_=np.linalg.qr(X0)
    def resid(v): return v-Q@(Q.T@v)
    ry=resid(y); ra=resid(a); den=float(ra@ra)
    if den<=0: return np.nan,np.nan
    beta=float(ry@ra/den)
    gidx=None if groups is None else [np.where(groups==g)[0] for g in np.unique(groups)]
    exceed=0
    for _ in range(n_perm):
        if gidx is None: ap_=rng.permutation(a)
        else:
            ap_=a.copy()
            for gi in gidx: ap_[gi]=rng.permutation(ap_[gi])
        rp=resid(ap_); dd=float(rp@rp)
        if dd<=0: continue
        exceed+=(abs(float(ry@rp/dd))>=abs(beta))
    return beta,(1.0+exceed)/(n_perm+1.0)

## -- cohort ----
d=idx.copy()
d["grp"]=np.where(d.timepoint=="Healthy","Healthy",np.where(d.timepoint.isin(AML_TP),"AML","other"))
d=d[d.grp!="other"].copy(); d["is_aml"]=(d.grp=="AML").astype(float)
nc=nodes.pivot_table(index=["dataset","sample"],columns="hierarchy_bin",values="n_cells_raw",aggfunc="first").reset_index()
for b in FGW_NODES:
    if b not in nc: nc[b]=0
nc["total"]=nc[FGW_NODES].sum(axis=1); nc["blast_proxy"]=nc[BLAST_BINS].sum(axis=1)/nc["total"].clip(lower=1)
d=d.merge(nc[["dataset","sample","blast_proxy"]],on=["dataset","sample"],how="left").dropna(subset=["blast_proxy"]).reset_index(drop=True)

# TECHNICAL COVARIATES. Joined on Sample (capital S in the QC report) and z-scored so the column
# scale cannot dominate the QR. log10 of depth because library size is multiplicative.
COVAR_COLS=[]
if args.covar!="none":
    qc_p = args.qc_csv or os.path.join(args.root,"01_preprocess","03_qc_report__ALL.csv")
    if not os.path.exists(qc_p):
        raise SystemExit(f"--covar {args.covar} needs the QC report, not found: {qc_p}\n"
                         f"  Pass --qc_csv, or run with --covar none.")
    qc = pd.read_csv(qc_p).rename(columns={"Sample":"sample"})
    qc = qc[["dataset","sample","med_ncount_final"]].drop_duplicates(subset=["dataset","sample"])
    qc["depth"]=np.log10(qc["med_ncount_final"].clip(lower=1))
    if args.covar in ("ribo","both"):
        raise SystemExit("--covar ribo/both needs a per-sample ribosomal fraction, which is not in "
                         "03_qc_report__ALL.csv. It would have to be aggregated from the per-cell "
                         "ingest QC first; do that as its own step rather than inline here.")
    d = d.merge(qc[["dataset","sample","depth"]],on=["dataset","sample"],how="left")
    miss=int(d["depth"].isna().sum())
    if miss:
        # Dropping silently would change the cohort between arms and make the comparison meaningless.
        raise SystemExit(f"{miss} of {len(d)} sample(s) have no depth in the QC report. The arms must "
                         f"run on the SAME cohort or the comparison is not a comparison.")
    d["depth"]=(d["depth"]-d["depth"].mean())/d["depth"].std(ddof=0)
    COVAR_COLS=["depth"]
    print(f"[1] covariate: {args.covar} -> columns {COVAR_COLS} (z-scored, log10 for depth)")
# THE SPLIT. A dataset-level 70/30 Discovery/Validation assignment has existed since 2026-08-04 and
# no analysis respected it, which was tolerable only while the feature set was fixed a priori. A
# sweep over the panel families IS selection, so it must run on one arm.
# Healthy donors carry their own split label ("Healthy") because there are only 23 of them
# cohort-wide; splitting those would leave both arms underpowered, so controls are SHARED and only
# the AML side is genuinely held out. That limitation belongs in the Methods, not in a footnote.
if args.split != "all":
    sp_csv = args.split_csv or os.path.join(args.root, "01_preprocess", "02_sample_split.csv")
    if not os.path.exists(sp_csv):
        raise SystemExit(f"--split {args.split} needs the split table, not found: {sp_csv}\n"
                         f"  Pass --split_csv, or run with --split all and say so when reporting.")
    sp = pd.read_csv(sp_csv)[["dataset","Sample","split_sample"]].rename(columns={"Sample":"sample"})
    sp = sp.drop_duplicates(subset=["dataset","sample"])
    before = len(d)
    d = d.merge(sp, on=["dataset","sample"], how="left")
    unlabelled = int(d["split_sample"].isna().sum())
    if unlabelled:
        raise SystemExit(f"{unlabelled} sample(s) carry no split label. Assigning them here would be "
                         f"choosing the arm after seeing the data; fix 02_sample_split.csv instead.")
    keep = {"discovery": {"Discovery","Healthy"}, "validation": {"Validation","Healthy"}}[args.split]
    d = d[d["split_sample"].isin(keep)].reset_index(drop=True)
    n_aml = int((d.is_aml==1).sum()); n_h = int((d.is_aml==0).sum())
    print(f"[1] SPLIT={args.split}: {len(d)} of {before} samples | AML={n_aml} healthy={n_h} "
          f"(controls are shared between arms)")
    if n_aml < 20 or n_h < 8:
        raise SystemExit(f"arm too small to test: AML={n_aml}, healthy={n_h}")

heal_all=[(r.dataset,r["sample"]) for _,r in d[d.is_aml==0].iterrows()]
heal_bary=[(r.dataset,r["sample"]) for _,r in d[(d.is_aml==0)&(~d["sparse_flag"].fillna(False))].iterrows()]
aml_all=[(r.dataset,r["sample"]) for _,r in d[d.is_aml==1].iterrows()]
print(f"[1] healthy={len(heal_all)} (barycenter {len(heal_bary)}) | AML={len(aml_all)} | alphas={ALPHAS}")

rows=[]
for alpha in ALPHAS:
    for name,feats in SUBSETS.items():
        Cb_full,Fb_full,pb=barycenter(heal_bary,alpha,feats)
        hds={}
        for k in heal_all:
            loo=[x for x in heal_bary if x!=k]
            Cb,Fb,_=barycenter(loo,alpha,feats)
            C,F,p=build_one(*k,feats); hds[k]=fgw2(C,F,p,Cb,Fb,pb,alpha)
        for k in aml_all:
            C,F,p=build_one(*k,feats); hds[k]=fgw2(C,F,p,Cb_full,Fb_full,pb,alpha)
        dd=d.copy(); dd["HDS"]=[hds[(r.dataset,r["sample"])] for _,r in dd.iterrows()]
        y=dd["HDS"].to_numpy(float); a=dd["is_aml"].to_numpy(float)
        X0=np.column_stack([np.ones(len(dd)),dd["blast_proxy"].to_numpy(float)]
                          +[dd[c].to_numpy(float) for c in COVAR_COLS])
        bA,pA=fwl_perm(y,X0,a,None,args.n_perm,_rng(alpha,name,"global"))
        both=[ds for ds,g in dd.groupby("dataset") if g.is_aml.nunique()==2]
        db=dd[dd.dataset.isin(both)].reset_index(drop=True)
        yb=db["HDS"].to_numpy(float); ab=db["is_aml"].to_numpy(float)
        dums=pd.get_dummies(db["dataset"],drop_first=True).to_numpy(float)
        X0b=np.column_stack([np.ones(len(db)),db["blast_proxy"].to_numpy(float),dums]
                           +[db[c].to_numpy(float) for c in COVAR_COLS])
        bB,pB=fwl_perm(yb,X0b,ab,db["dataset"].to_numpy(),args.n_perm,_rng(alpha,name,"strat"))
        rows.append(dict(alpha=alpha,feature_set=name,features="+".join(feats),
                         mean_healthy=float(y[a==0].mean()),mean_aml=float(y[a==1].mean()),
                         beta_global=bA,p_global=pA,beta_strat=bB,p_strat=pB))
        print(f"[2] alpha={alpha:<4g} {name:14s} | healthy={y[a==0].mean():.4f} AML={y[a==1].mean():.4f} "
              f"| global b={bA:+.5f} p={pA:.5f} | within-ds b={bB:+.5f} p={pB:.5f}",flush=True)

res=pd.DataFrame(rows)

# BH-FDR WITHIN each panel family, not across all of them. The families are five separate questions
# (stemness robustness / PROGENy pathways / cell state / metabolism+drug target / pseudotime);
# pooling them would spend the correction on unrelated hypotheses and bury the stemness arm under 14
# PROGENy pathways. The core subsets are a pre-registered hypothesis and are NOT corrected here --
# correcting a pre-specified test against an exploratory sweep would be the same mistake inverted.
def _bh(pv):
    pv = np.asarray(pv, float); n = len(pv); o = np.argsort(pv)
    q = np.empty(n); q[o] = np.minimum.accumulate((pv[o] * n / np.arange(1, n+1))[::-1])[::-1]
    return np.minimum(q, 1.0)

# The family map is LOADED from the vocabulary, never inferred from the name. jsonlite writes a
# length-1 R list element as a bare string, so normalise to str either way.
_FAM_MAP = {k: (v if isinstance(v, str) else v[0])
            for k, v in (_VOCAB.get("candidate_family") or {}).items()}
if not _FAM_MAP:
    raise SystemExit(
        "fgw_vocab.json carries no candidate_family map.\n"
        "  Re-run scripts/07_fgw/01_build_fgw_inputs.R --force: it emits the map from CCC_PANELS,\n"
        "  and without it the BH correction cannot know which features form a family.")
res["family"] = [family_of(n, _FAM_MAP) for n in res["feature_set"]]
res["q_strat"] = np.nan
for (al, fam), g in res.groupby(["alpha","family"]):
    if fam == "core" or len(g) < 2:
        continue
    res.loc[g.index, "q_strat"] = _bh(g["p_strat"].to_numpy())
res["split"] = args.split
# THE SPLIT GOES IN THE FILENAME, NOT ONLY IN A COLUMN. It used to be a column alone, so every run
# landed on the same feature_decomposition.csv regardless of --split, --subsets or --n_perm. That is
# not hypothetical: on 2026-08-28 a --split all run silently overwrote the file that a reader would
# take to be "the screen", leaving pt at q_strat=0.063 (pooled cohort) where the pre-registered
# Discovery screen has 0.031. A Validation confirmation run overwriting the Discovery screen that
# defines its own hypotheses would be unrecoverable.
_out = args.out or ("feature_decomposition.csv" if args.split == "all"
                    else f"feature_decomposition__{args.split}.csv")
_only_tag = f" [--only {args.only}]" if args.only else ""
res.to_csv(os.path.join(D_OUT,_out),index=False)
print(f"[done] wrote {os.path.join(D_OUT,_out)}  "
      f"({len(res)} rows, split={args.split}, n_perm={args.n_perm}){_only_tag}")

fams = [f for f in res["family"].unique() if f != "core"]
if fams:
    print("\n[3b] PANEL FAMILIES -- BH-FDR within each family, alpha=0, within-dataset model")
    a0 = res[res.alpha == min(ALPHAS)]
    for fam in sorted(fams):
        g = a0[a0.family == fam].sort_values("p_strat")
        n_sig = int((g["q_strat"] < 0.05).sum())
        print(f"  {fam}: {len(g)} features tested | q<0.05: {n_sig}")
        for _, r in g.head(3).iterrows():
            mark = " *" if r["q_strat"] < 0.05 else ""
            print(f"      {r['feature_set'][:44]:46s} beta={r['beta_strat']:+.4f} "
                  f"p={r['p_strat']:.4f} q={r['q_strat']:.4f}{mark}")
    if not (a0[a0.family != "core"]["q_strat"] < 0.05).any():
        print("  nothing survives FDR in any family. Report that, not the smallest raw p.")
print("\n[3] FEATURE DECOMPOSITION")
print(res[["alpha","feature_set","mean_healthy","mean_aml","beta_global","p_global","beta_strat","p_strat"]].round(5).to_string(index=False))

print("\n[4] VERDICT (alpha=0, within-dataset model)")
z=res[(res.alpha==0.0)].set_index("feature_set")
def g(k,c):
    return float(z.loc[k,c]) if k in z.index else float("nan")
p_no=g("no_frac_mal","p_strat"); p_only=g("only_frac_mal","p_strat")
print(f"    only_frac_mal : p={p_only:.5f}   (the feature we force to 0 for healthy)")
print(f"    no_frac_mal   : p={p_no:.5f}   (stemness + n_cells, no encoded label)")
if p_only<0.05 and p_no>=0.05:
    print("    -> CIRCULAR: the separation comes from the feature we encoded by design.")
    print("       It cannot be reported as a finding.")
elif p_no<0.05:
    print("    -> A non-circular feature signal survives. Check only_stemness vs only_ncells below to see")
    print("       whether it is biological (stemness) or merely compositional (n_cells, already covered by blast_proxy).")
    print(f"       only_stemness p={g('only_stemness','p_strat'):.5f} | only_ncells p={g('only_ncells','p_strat'):.5f}")
else:
    print("    -> neither subset is significant on its own; the signal needs the full feature combination.")

# no_frac_mal bundles a non-encoded feature (stemness) with a null one (n_cells), so it answers two
# questions at once and answers neither well. The single-feature rows are what the decision rests on.
if any(k in z.index for k in CANDIDATE_SUBSETS):
    print("\n[5] CANDIDATE FEATURES (alpha=0, within-dataset model)")
    for k in CANDIDATE_SUBSETS:
        if k in z.index:
            print(f"    {k:<20} beta={g(k,'beta_strat'):+.5f}  p={g(k,'p_strat'):.5f}"
                  f"   healthy={g(k,'mean_healthy'):.4f} AML={g(k,'mean_aml'):.4f}")
    pn, pm = g("only_stem_normal","p_strat"), g("only_stem_malignant","p_strat")
    if np.isfinite(pn) and np.isfinite(pm):
        print("\n    DECISIVE CONTRAST -- which cells carry the stemness signal?")
        if pn<0.05 and pm>=0.05:
            print("    -> NON-MALIGNANT cells. The residual hematopoiesis is shifted; this is a")
            print("       microenvironment finding and is not reducible to blast content.")
        elif pm<0.05 and pn>=0.05:
            print("    -> MALIGNANT cells only. This reduces to 'blasts are stem-like' -- already known,")
            print("       and near-circular with malignancy itself. Weak as a finding.")
        elif pn<0.05 and pm<0.05:
            print("    -> BOTH compartments. Report the normal-cell arm as the novel part; the malignant")
            print("       arm is expected and should be presented as a positive control, not a result.")
        else:
            print("    -> NEITHER survives once split. The pooled stemness signal was carried by the")
            print("       mixture itself (composition), not by either compartment.")
    # A feature imputed on most nodes is near-constant after scaling, so a null result measures the
    # ENCODING, not the biology. frac_malignant_vg is the live case: the van Galen axis passed held-out
    # AUC in HSC_MPP only, so 86% of its nodes are imputed and its post-z sd is ~0.37 against ~0.95 for
    # the others. Reporting "malignancy does not separate" off that row would be wrong -- the same axis
    # separates at SAMPLE level (AUC 0.822, scripts/02_malignancy/81). Say so instead of concluding.
    IMPUTED_MAX = 0.50
    imp = {c: float((nodes[c] == 0.0).mean()) for c in _want}
    pv = g("only_vg_mal","p_strat")
    if np.isfinite(pv):
        iv = imp.get("frac_malignant_vg", 0.0)
        print(f"\n    only_vg_mal p={pv:.5f} -- malignancy WITHOUT the hand-coded healthy zero.")
        if iv > IMPUTED_MAX:
            print(f"    UNINFORMATIVE: {iv:.1%} of its nodes are imputed (axis validated in HSC_MPP only),")
            print( "    so as a 7-node feature it is near-constant. This tests the encoding, not malignancy.")
            print( "    Draw no conclusion here; use the sample-level test in 02_malignancy/81 instead.")
        elif pv < 0.05:
            print("    Malignancy separates on its own merits; only_frac_mal was not pure artifact.")
        else:
            print("    Does not separate. The only_frac_mal result was the encoding, not malignancy.")
    heavy = {c: v for c, v in imp.items() if v > IMPUTED_MAX}
    if heavy:
        print("\n    imputation warning (>50% of nodes imputed -> null results uninformative):")
        for c, v in sorted(heavy.items(), key=lambda kv: -kv[1]): print(f"      {c:<26} {v:.1%}")

