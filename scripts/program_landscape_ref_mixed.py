"""Prototype: mixed-model-type test of reference-scaled coupling differences.

Per cell (population x state) the coupling is a mean over cells of a product
z = S~ * T~ / (sd_T * sigma_ref) (group-centred program score times target,
or target minus its matched-gene background for `excess`). A random-intercept
model on z,  z_ij = mu + b_i + e_ij  (donor i, cell j), combines the
between-donor variance sb2 and the within-donor variance se2:
    mu_hat = sum_i w_i zbar_i / sum_i w_i,   w_i = 1 / (sb2 + se2 / n_i),
    se(mu_hat) = 1 / sqrt(sum_i w_i).
sb2 = 0 gives the pooled estimate with the cells-only standard error; large
sb2 gives equal donor weights. sb2 is the unbalanced one-way ANOVA moment
estimate, truncated at 0. sd_T is pooled over the cell and treated as fixed.

compare    A vs B (donor effects of the two cells treated as independent)
           -> LANDSCAPE_OUT/test_mixed/{A}__vs__{B}.csv
calibrate  random donor halves of one cell against each other
           -> LANDSCAPE_OUT/test_mixed/calibration.csv
"""
import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm

import program_landscape_ref_memento as M
from memento.program import _prepare

OUT = M.LANDSCAPE_OUT / "test_mixed"
OUT.mkdir(parents=True, exist_ok=True)
STATS = ("beta", "excess")


def donor_stats(cell, progs, units, sigma):
    """Per donor: n, sum z, sum z^2 for beta and excess (donors x programs)."""
    X, o, genes, anchors = M.load_cell(cell)
    X = sp.csr_matrix(X)
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    rr = np.concatenate([gpos[progs[u]["genes"]].values for u in units])
    cc = np.concatenate([np.full(len(progs[u]["genes"]), k) for k, u in enumerate(units)])
    vv = np.concatenate([progs[u][M.SCHEME] / progs[u]["s_ref"] for u in units])
    W = sp.csc_matrix((vv, (rr, cc)), shape=(len(genes), len(units)))
    W2 = W.multiply(W).tocsc()
    cols = gpos[[M.TARGET] + list(anchors)].values
    member = np.asarray(W[cols].T.toarray() != 0)            # programs x targets
    key = (o.donor.astype(str) + "^" + o.fine.astype(str)).values
    prepared = []
    var_t, N = np.zeros(len(cols)), 0
    for g in pd.unique(key):
        idx = np.where(key == g)[0]
        s, _, t, t2 = _prepare(X[idx], o.memento_size_factor.values[idx].astype(float), M.Q, W, W2, cols)
        prepared.append((o.donor.values[idx[0]], s, t))
        var_t += t2.sum(0)
        N += len(idx)
    sd_t = np.sqrt(np.where(var_t > 0, var_t / N, np.nan))
    keep = (~member[:, 1:]).astype(float) * np.isfinite(sd_t[1:])   # anchors usable per program
    count = keep.sum(1)
    donors = sorted(set(o.donor))
    di = {d: i for i, d in enumerate(donors)}
    n = np.zeros(len(donors))
    acc = {k: (np.zeros((len(donors), len(units))), np.zeros((len(donors), len(units)))) for k in STATS}
    sg = sigma.values[None, :]
    for d, s, t in prepared:
        u = t / sd_t
        bg = np.nan_to_num(u[:, 1:]) @ keep.T / count
        z = {"beta": s * u[:, :1] / sg, "excess": s * (u[:, :1] - bg) / sg}
        i = di[d]
        n[i] += len(s)
        for k in STATS:
            acc[k][0][i] += z[k].sum(0)
            acc[k][1][i] += (z[k] ** 2).sum(0)
    dead = member[:, 0]
    for k in STATS:
        acc[k][0][:, dead] = np.nan
    return np.array(donors), n, acc


HET = False   # --het: donor-specific residual variances (DerSimonian-Laird) instead of one pooled se2


def fit_het(n, s1, s2):
    N, k = n.sum(), len(n)
    zbar = s1 / n[:, None]
    pooled = s1.sum(0) / N
    v = (s2 / n[:, None] - zbar ** 2) / (n[:, None] - 1)          # variance of each donor mean
    w = 1 / v
    fixed = (w * zbar).sum(0) / w.sum(0)
    q = (w * (zbar - fixed) ** 2).sum(0)
    sb2 = np.clip((q - (k - 1)) / (w.sum(0) - (w ** 2).sum(0) / w.sum(0)), 0, None)
    w = 1 / (v + sb2[None, :])
    mu = (w * zbar).sum(0) / w.sum(0)
    se2 = (v * n[:, None] * n[:, None]).sum(0) / N               # cell-weighted per-cell variance
    return mu, 1 / np.sqrt(w.sum(0)), pooled, np.sqrt(se2 / N), sb2, se2


def fit(n, s1, s2):
    """Random-intercept fit per program. Returns mu, se, pooled mean, cells-only se, sb2, se2."""
    if HET:
        return fit_het(n, s1, s2)
    N, k = n.sum(), len(n)
    zbar = s1 / n[:, None]
    pooled = s1.sum(0) / N
    se2 = (s2.sum(0) - (s1 ** 2 / n[:, None]).sum(0)) / (N - k)
    msb = (n[:, None] * (zbar - pooled) ** 2).sum(0) / (k - 1)
    n0 = (N - (n ** 2).sum() / N) / (k - 1)
    sb2 = np.clip((msb - se2) / n0, 0, None)
    w = 1 / (sb2[None, :] + se2[None, :] / n[:, None])
    mu = (w * zbar).sum(0) / w.sum(0)
    return mu, 1 / np.sqrt(w.sum(0)), pooled, np.sqrt(se2 / N), sb2, se2


def contrast(a, b):
    out = {}
    for k in STATS:
        fa, fb = fit(a[1], *a[2][k]), fit(b[1], *b[2][k])
        d, se = fa[0] - fb[0], np.hypot(fa[1], fb[1])
        out[k] = dict(d=d, se=se, p=2 * norm.sf(np.abs(d / se)), d_pooled=fa[2] - fb[2],
                      se_cells=np.hypot(fa[3], fb[3]), mu_a=fa[0], mu_b=fb[0], pooled_a=fa[2], pooled_b=fb[2],
                      icc_a=fa[4] / (fa[4] + fa[5]), icc_b=fb[4] / (fb[4] + fb[5]))
    return out


def cmd_compare(args):
    progs, units, sigma = M.setup()
    A, B = donor_stats(args.a, progs, units, sigma), donor_stats(args.b, progs, units, sigma)
    out = pd.DataFrame({"unit": units})
    for cell, tag, S in ((args.a, "a", A), (args.b, "b", B)):
        act, full = M.active(cell, units)
        out[f"active_{tag}"] = act
        pooled = S[2]["beta"][0].sum(0) / S[1].sum()
        assert np.allclose(pooled, full.beta.values, rtol=1e-3, atol=1e-5, equal_nan=True), cell
    c = contrast(A, B)
    for k in STATS:
        for name, v in c[k].items():
            out[f"{k}_{name}"] = v
        out[f"{k}_q"] = M.bh(out[f"{k}_p"])
        out[f"{k}_q_cells"] = M.bh(2 * norm.sf(np.abs(out[f"{k}_d_pooled"] / out[f"{k}_se_cells"])))
    out.to_csv(OUT / f"{args.a}__vs__{args.b}{'__het' if HET else ''}.csv", index=False)
    e = out[out.active_a | out.active_b]
    print(f"{args.a} ({len(A[0])} donors, {int(A[1].sum())} cells) vs {args.b} ({len(B[0])} donors, {int(B[1].sum())} cells); "
          f"{len(e)} programs")
    for k in STATS:
        print(f"d_{k}: mixed q<0.05 {(e[f'{k}_q'] < .05).sum()}, p<0.05 {(e[f'{k}_p'] < .05).sum()}, median se {e[f'{k}_se'].median():.4f} | "
              f"sb2=0 (cells-only analytic) q<0.05 {(e[f'{k}_q_cells'] < .05).sum()}, median se {e[f'{k}_se_cells'].median():.4f} | "
              f"mixed vs pooled estimate r {np.corrcoef(e[f'{k}_d'], e[f'{k}_d_pooled'])[0, 1]:.3f} | "
              f"programs with sb2>0: a {(e[f'{k}_icc_a'] > 0).mean():.2f}, b {(e[f'{k}_icc_b'] > 0).mean():.2f}; "
              f"median donor share of variance where >0: a {e[f'{k}_icc_a'][e[f'{k}_icc_a'] > 0].median():.4f}, "
              f"b {e[f'{k}_icc_b'][e[f'{k}_icc_b'] > 0].median():.4f}")


def cmd_calibrate(args):
    progs, units, sigma = M.setup()
    rng = np.random.default_rng(M.SEED)
    rows = []
    for cell in args.cells:
        donors, n, acc = donor_stats(cell, progs, units, sigma)
        act, _ = M.active(cell, units)
        for i in range(args.splits):
            h = np.zeros(len(donors), bool)
            h[rng.permutation(len(donors))[: len(donors) // 2]] = True
            sub = lambda m: (donors[m], n[m], {k: (acc[k][0][m], acc[k][1][m]) for k in STATS})
            c = contrast(sub(h), sub(~h))
            rec = dict(cell=cell, split=i, donors=len(donors))
            for k in STATS:
                z = (c[k]["d"] / c[k]["se"])[act]
                zc = (c[k]["d_pooled"] / c[k]["se_cells"])[act]
                rec.update({f"{k}_p05": (np.abs(z) > 1.96).mean(), f"{k}_z_sd": np.nanstd(z),
                            f"{k}_q05": (M.bh(2 * norm.sf(np.abs(z))) < .05).mean(),
                            f"{k}_cells_p05": (np.abs(zc) > 1.96).mean(), f"{k}_cells_z_sd": np.nanstd(zc)})
            rows.append(rec)
        d = pd.DataFrame([r for r in rows if r["cell"] == cell])
        print(cell, len(donors), "donors;", d.drop(columns=["cell", "split", "donors"]).median().round(3).to_dict(),
              "| max excess p05", d.excess_p05.max().round(3), "max excess q05", d.excess_q05.max().round(3), flush=True)
    pd.DataFrame(rows).to_csv(OUT / f"calibration{'_het' if HET else ''}.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare"); c.add_argument("--a", required=True); c.add_argument("--b", required=True)
    k = sub.add_parser("calibrate"); k.add_argument("--cells", nargs="+", required=True)
    k.add_argument("--splits", type=int, default=200)
    for q in (c, k):
        q.add_argument("--het", action="store_true")
    a = ap.parse_args()
    HET = a.het
    {"compare": cmd_compare, "calibrate": cmd_calibrate}[a.cmd](a)
