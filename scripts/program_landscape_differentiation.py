"""Cross-dataset replication, step 1: does ANTXR2 rise with epithelial
differentiation at the mean level? (planning_summaries/program_coexpression.md,
"Cross-dataset replication").

Both tests compare ANTXR2 with 50 genes of its genomic length matched on raw
mean + detection (span window; candidates above LOW_FLOOR), so a shift every
long or every similarly expressed gene shows is not counted.

A. ladder   each epithelial population against the atlas's stem population
            (LADDER), per donor with >= MIN_GROUP_CELLS cells of both, size
            factors computed jointly. Anchors are matched in the stem
            population only (matching in both would select genes with
            ANTXR2's fold change). lfc = log2 mean(X/sf) population / stem.
B. within   inside each absorptive population x state, per donor, top vs
            bottom quartile of the fixed maturation score
            (program_config.MATURATION_MARKERS: mature markers minus OLFM4),
            anchors matched in that population x state.

excess = lfc_ANTXR2 - mean lfc_anchor, z = excess / anchor SD, pct = share of
anchors below ANTXR2, donors_above = share of donors where ANTXR2's lfc
exceeds the median anchor's.

Writes LANDSCAPE_OUT/differentiation/{ladder,within}.csv.
ANTXR2_LANDSCAPE=gca runs it on the Gut Cell Atlas.

Usage:
    python program_landscape_differentiation.py
"""
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import (                       # noqa: E402
    LANDSCAPE_ATLAS, LANDSCAPE_CELLS, LANDSCAPE_OUT, LANDSCAPE_STATES, MATURATION_MARKERS,
    MIN_GROUP_CELLS, N_ANCHORS, TARGET, gene_spans, pop_slug,
)
from program_landscape_cycling import LOW_FLOOR, lfc, normalise, summarise  # noqa: E402
from program_scores import gene_stats, load_cells  # noqa: E402

OUT = LANDSCAPE_OUT / "differentiation"
OUT.mkdir(parents=True, exist_ok=True)
STATE_OF = {v: k for k, v in LANDSCAPE_STATES.items()}

# stem population, then the populations compared with it (least to most
# differentiated within each branch); `absorptive` get test B
LADDER = {
    "ibdverse": dict(
        stem="Enterocyte stem OLFM4lo LGR5hi",
        absorptive=["Enterocyte stem MKI67", "Enterocyte progenitor crypt OLFM4hi",
                    "Enterocyte progenitor/precursor crypt OLFM4lo", "Enterocyte top villus"],
        secretory=["Ileum goblet MUC2lo BCAS1int", "Ileum goblet MUC2int",
                   "Ileum goblet MUC2hi SPINK4hi", "Ileum goblet MUC2hi BCAS1hi"]),
    "gca": dict(
        stem="Epithelial_stem",
        absorptive=["TA", "Enterocyte"],
        secretory=["Goblet_cycling", "Goblet_progenitor", "Goblet", "Paneth"]),
}[LANDSCAPE_ATLAS]


def load(pop):
    return load_cells(LANDSCAPE_CELLS / f"cells_{pop_slug(pop)}.npz")


def anchors_in(X, genes, rows):
    """N_ANCHORS span-matched genes nearest ANTXR2 on raw mean + detection in `rows`."""
    from program_landscape import span_matched   # binds LANDSCAPE_OUT; import late
    raw_mean, det = gene_stats(X, rows)
    f = pd.DataFrame({"log10_mean": np.log10(raw_mean + 1e-6), "det": det}, index=genes)
    anchors, _ = span_matched({"c": f}, [set(genes[raw_mean >= LOW_FLOOR])],
                              np.log10(gene_spans()), N_ANCHORS)
    return anchors, raw_mean[genes.get_loc(TARGET)]


def summary(L, **extra):
    r = summarise(L, **extra)
    r["donors_above"] = 1 - r.pop("donors_below")
    return r


def ladder():
    stem = LADDER["stem"]
    Xs, os_, genes = load(stem)
    rows = []
    for branch in ("absorptive", "secretory"):
        for pop in LADDER[branch]:
            Xp, op, gp = load(pop)
            assert gp.equals(genes)
            X = sp.vstack([Xs, Xp]).tocsr()
            obs = pd.concat([os_, op], ignore_index=True)
            obs["is_stem"] = np.arange(len(obs)) < len(os_)
            obs["state"] = obs.disease.map(STATE_OF)
            Y, _ = normalise(X, genes)
            for st, o in obs[obs.state.notna()].groupby("state"):
                n = o.groupby(["donor", "is_stem"]).size().unstack(fill_value=0)
                if n.shape[1] < 2:
                    continue
                donors = n.index[(n[True] >= MIN_GROUP_CELLS) & (n[False] >= MIN_GROUP_CELLS)]
                if len(donors) < 3:
                    continue
                oo = o[o.donor.isin(donors)]
                anchors, stem_mean = anchors_in(X, genes, oo.index[oo.is_stem].values)
                T = np.array([genes.get_loc(g) for g in [TARGET] + anchors])
                L = np.array([lfc(Y, oo.index[(oo.donor == d) & ~oo.is_stem].values,
                                  oo.index[(oo.donor == d) & oo.is_stem].values, T)
                              for d in donors])
                pop_mean = gene_stats(X, oo.index[~oo.is_stem].values)[0][genes.get_loc(TARGET)]
                rows.append(summary(L, branch=branch, population=pop, state=st,
                                    antxr2_raw_stem=stem_mean, antxr2_raw_pop=pop_mean))
            print(f"ladder {pop}: done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "ladder.csv", index=False)
    return out


def maturation_score(Y, scale, genes, rows):
    """Mean z of log1p mature markers minus mean z of the immature ones, over
    markers detected in >= 1% of `rows`; None if either side has none."""
    s = []
    for side in ("mature", "immature"):
        idx = [genes.get_loc(g) for g in MATURATION_MARKERS[side] if g in genes]
        L = Y[rows][:, idx].toarray() * scale
        L = np.log1p(L[:, (L > 0).mean(0) >= 0.01])
        if L.shape[1] == 0:
            return None
        s.append(((L - L.mean(0)) / (L.std(0) + 1e-9)).mean(1))
    return s[0] - s[1]


def within():
    rows = []
    for pop in [LADDER["stem"]] + LADDER["absorptive"]:
        X, obs, genes = load(pop)
        Y, scale = normalise(X, genes)
        state = obs.disease.map(STATE_OF)
        for st in LANDSCAPE_STATES:
            ok = (state == st).values
            if ok.sum() == 0:
                continue
            anchors, raw = anchors_in(X, genes, np.where(ok)[0])
            T = np.array([genes.get_loc(g) for g in [TARGET] + anchors])
            L = []
            for d, ix in obs[ok].groupby("donor").groups.items():
                ix = np.asarray(ix)
                if len(ix) < MIN_GROUP_CELLS:
                    continue
                s = maturation_score(Y, scale, genes, ix)
                if s is None:
                    continue
                r = pd.Series(s).rank(method="first").values
                L.append(lfc(Y, ix[r > 0.75 * len(ix)], ix[r <= 0.25 * len(ix)], T))
            if len(L) >= 3:
                rows.append(summary(np.array(L), population=pop, state=st, antxr2_raw_mean=raw))
        print(f"within {pop}: done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "within.csv", index=False)
    return out


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    print(ladder().round(3).to_string())
    print(within().round(3).to_string())
