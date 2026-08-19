"""Single source for the FGW/scoring vocabulary on the Python side.

WHY THIS MODULE EXISTS. The R config derives the AML timepoint set from CANONICAL_TIMEPOINTS
(config_fgw.R:48) specifically so a curation change cannot shrink B_AML, and its comment names the
failure mode outright: "Spelling out c(\"Diagnosis\",\"MRD\",\"Post_treatment\",...) is exactly how
that would have happened." Nine Python files then each hard-coded exactly that literal.

CANONICAL_TIMEPOINTS was migrated on 2026-08-04 -- MRD and Post_treatment retired, Refractory /
On_treatment / Post_induction / Post_consolidation / Post_transplant / Post_treatment_unspecified
added -- and the Python mirror was not. Measured against the 214-sample roster, the stale literal
sends 34 samples to grp == "other", where every downstream filter deletes them: 17 Post_induction,
8 Post_treatment_unspecified, 7 On_treatment, 1 Post_consolidation, 1 Refractory. That is 64% of
the treated arm, which is the arm H3 is about, and no stage prints an "other" count -- the row
count of patient_scores.csv is unchanged because the alignment keeps one row per sample.

So the vocabulary is READ, never written here. 01_build_fgw_inputs.R emits fgw_vocab.json next to
fgw_input_index.csv; load_vocab() fails loudly if it is absent rather than falling back to a
literal, because a silent fallback is the bug this module replaces.
"""
import json
import os


class VocabError(RuntimeError):
    pass


def load_vocab(d_fgw):
    """Load fgw_vocab.json from the 07_fgw output directory. Raises if missing or malformed."""
    p = os.path.join(d_fgw, "fgw_vocab.json")
    if not os.path.exists(p):
        raise VocabError(
            "missing %s -- run scripts/07_fgw/01_build_fgw_inputs.R first. Refusing to fall back to "
            "a hard-coded timepoint list: the stale literal silently deleted 34 of 214 samples "
            "(64%% of the treated arm) and that is the failure this file exists to prevent." % p)
    with open(p) as fh:
        v = json.load(fh)
    for k in ("aml_timepoints", "healthy_timepoints", "all_timepoints", "nodes", "features"):
        if k not in v:
            raise VocabError("%s is missing key %r" % (p, k))
    # jsonlite's auto_unbox writes a length-1 vector as a bare STRING, not a list. set("Unknown")
    # is then {'U','n','k','o','w'} -- a set of characters -- so the label never matches and every
    # Unknown-timepoint sample reads as uncovered. Normalise scalars to lists before making a set.
    def _as_set(x):
        if x is None:
            return set()
        return set(x if isinstance(x, (list, tuple)) else [x])

    v["aml_timepoints"] = _as_set(v["aml_timepoints"])
    v["healthy_timepoints"] = _as_set(v["healthy_timepoints"])
    v["excluded_timepoints"] = _as_set(v.get("excluded_timepoints"))
    return v


def assert_index_covered(idx, vocab, tp_col="timepoint"):
    """Every timepoint in the index must be a label the vocabulary knows.

    Without this, an unrecognised label is not an error -- it is a sample that quietly stops being
    AML and stops being healthy, and therefore stops existing for every test downstream.
    """
    known = vocab["aml_timepoints"] | vocab["healthy_timepoints"] | vocab.get("excluded_timepoints", set())
    seen = set(str(x) for x in idx[tp_col].dropna().unique())
    unknown = seen - known
    if unknown:
        raise VocabError(
            "timepoint(s) in the index that the vocabulary does not cover: %s. These would be "
            "dropped from every AML/healthy split without a message." % ", ".join(sorted(unknown)))
    n_aml = int(idx[tp_col].isin(vocab["aml_timepoints"]).sum())
    n_heal = int(idx[tp_col].isin(vocab["healthy_timepoints"]).sum())
    n_other = len(idx) - n_aml - n_heal
    print("[vocab] %d AML, %d healthy, %d excluded (of %d rows in the index)"
          % (n_aml, n_heal, n_other, len(idx)))
    return n_aml, n_heal, n_other
