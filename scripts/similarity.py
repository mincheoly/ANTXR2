"""Population x population similarity of the ANTXR2 co-expression profile.

Input: a directory of <population>.csv files written by antxr2_coexpr.py
(columns population, gene, corr_coef, corr_se, corr_pval, qval, raw_mean,
det_rate).

Similarity definitions (all computed on the genes tested in BOTH populations,
which differ between tissues because the raw-mean floor is applied per
population):

  pearson_all     Pearson r of corr_coef over all shared genes
  spearman_all    rank version, robust to the heavy tails of the estimator
  pearson_sig     Pearson r restricted to genes reaching q<0.05 in at least
                  one of the two populations -- the informative subset, since
                  most shared genes are noise around zero
  jaccard_sig     signed Jaccard of the q<0.05 gene sets: |shared same-sign| /
                  |union|, a set-level view that ignores effect size

The detection control matrix uses det_rate (fraction of cells with a nonzero
count) in place of corr_coef: two tissues whose co-expression profiles look
alike only because their detectable-gene sets look alike will score high there
too, so the co-expression matrix should be read against it.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import rankdata

SIG = 0.05


def load_profiles(indir, min_genes=2000):
    """-> dict population -> DataFrame indexed by gene."""
    out = {}
    for p in sorted(glob.glob(os.path.join(indir, "*.csv"))):
        pop = os.path.splitext(os.path.basename(p))[0]
        df = pd.read_csv(p)
        if len(df) < min_genes:
            continue
        df = df.dropna(subset=["corr_coef"]).drop_duplicates("gene").set_index("gene")
        out[pop] = df
    return out


def _pearson(a, b):
    if len(a) < 50:
        return np.nan
    sa, sb = a.std(), b.std()
    if not np.isfinite(sa) or not np.isfinite(sb) or sa == 0 or sb == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def similarity(profiles, value="corr_coef", sig_col="qval"):
    pops = list(profiles)
    n = len(pops)
    keys = ["pearson_all", "spearman_all", "pearson_sig", "jaccard_sig", "n_shared", "n_sig_union"]
    M = {k: pd.DataFrame(np.nan, index=pops, columns=pops, dtype=float) for k in keys}
    for i, pi in enumerate(pops):
        for j in range(i, n):
            pj = pops[j]
            A, B = profiles[pi], profiles[pj]
            g = A.index.intersection(B.index)
            a, b = A.loc[g, value].to_numpy(), B.loc[g, value].to_numpy()
            ok = np.isfinite(a) & np.isfinite(b)
            a, b, g = a[ok], b[ok], g[ok]
            res = {"n_shared": len(g)}
            res["pearson_all"] = _pearson(a, b)
            res["spearman_all"] = _pearson(rankdata(a), rankdata(b)) if len(a) >= 50 else np.nan
            if sig_col in A.columns and sig_col in B.columns:
                sa = A.loc[g, sig_col].to_numpy() < SIG
                sb = B.loc[g, sig_col].to_numpy() < SIG
                u = sa | sb
                res["n_sig_union"] = int(u.sum())
                res["pearson_sig"] = _pearson(a[u], b[u]) if u.sum() >= 50 else np.nan
                both = sa & sb
                same = both & (np.sign(a) == np.sign(b))
                res["jaccard_sig"] = float(same.sum() / u.sum()) if u.sum() else np.nan
            for k, v in res.items():
                M[k].iloc[i, j] = M[k].iloc[j, i] = v
    return M


def reliability(df, value="corr_coef", se="corr_se", genes=None):
    """Fraction of the observed variance in a population's co-expression
    profile that is signal rather than bootstrap noise: 1 - E[se^2]/Var(coef).

    This is the quantity that attenuates a between-population correlation, and
    it varies ~200-fold across these populations (cells, donors and ANTXR2
    abundance all differ), so the RAW similarity matrix is dominated by which
    populations were well powered rather than by biology.
    """
    d = df if genes is None else df.loc[genes]
    a = d[value].to_numpy()
    s = d[se].to_numpy()
    ok = np.isfinite(a) & np.isfinite(s)
    if ok.sum() < 50:
        return np.nan
    v = np.var(a[ok])
    return float(max(0.0, 1 - np.mean(s[ok] ** 2) / v)) if v > 0 else np.nan


def disattenuated(profiles, value="corr_coef", se="corr_se"):
    """Spearman's correction for attenuation: r_ij / sqrt(rel_i * rel_j),
    both reliabilities computed on the genes the pair shares. Returns
    (corrected matrix, pairwise sqrt(rel_i*rel_j)); the diagonal of the second
    is each population's own reliability."""
    pops = list(profiles)
    n = len(pops)
    R = pd.DataFrame(np.nan, index=pops, columns=pops, dtype=float)
    K = pd.DataFrame(np.nan, index=pops, columns=pops, dtype=float)
    for i, pi in enumerate(pops):
        for j in range(i, n):
            pj = pops[j]
            A, B = profiles[pi], profiles[pj]
            g = A.index.intersection(B.index)
            if len(g) < 50:
                continue
            ra = reliability(A, value, se, g)
            rb = reliability(B, value, se, g)
            r = _pearson(A.loc[g, value].to_numpy(), B.loc[g, value].to_numpy())
            d = np.sqrt(ra * rb) if np.isfinite(ra) and np.isfinite(rb) else np.nan
            K.iloc[i, j] = K.iloc[j, i] = d
            if np.isfinite(d) and d > 0:
                R.iloc[i, j] = R.iloc[j, i] = float(np.clip(r / d, -1, 1))
    return R, K


def cluster_order(M):
    """Average-linkage leaf order on 1 - similarity."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    V = M.to_numpy(dtype=float).copy()
    V[~np.isfinite(V)] = 0.0
    D = 1 - V
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2
    Z = linkage(squareform(D, checks=False), method="average")
    return [M.index[i] for i in leaves_list(Z)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("outprefix")
    ap.add_argument("--value", default="corr_coef")
    a = ap.parse_args()
    prof = load_profiles(a.indir)
    print(f"{len(prof)} populations", flush=True)
    M = similarity(prof, value=a.value)
    for k, m in M.items():
        m.to_csv(f"{a.outprefix}_{k}.csv")
    print("wrote", a.outprefix, flush=True)


if __name__ == "__main__":
    main()
