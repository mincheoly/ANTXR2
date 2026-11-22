"""Reference-scaled coupling tests through memento.ht_program_correlation
(memento branch feature/program-correlation): the same estimator as
program_landscape_ref.score, with memento's bootstrap over donors (default),
over donors and cells within each donor x subtype group, or over cells only
(--mode). Supersedes program_landscape_ref_test.py, whose donor-only numbers
it reproduces.

Each population x state cell becomes one memento condition. Its size factors
are multiplied by the cell's per-10k `scale`, so one weight matrix
(v / s_ref) serves every condition; `scale=` is the reference sigma, and the
population's anchors are the `references` (statistic `excess`).

compare     A vs B                 -> LANDSCAPE_OUT/test_memento/{A}__vs__{B}.csv
calibrate   random donor halves of one cell against each other (truth = 0)
                                   -> LANDSCAPE_OUT/test_memento/calibration.csv

Usage (memento worktree first on PYTHONPATH):
    PYTHONPATH=/home/ubuntu/Github/memento-program-correlation ANTXR2_LANDSCAPE_RUN=ref \
        python program_landscape_ref_memento.py compare --a CELL --b CELL [--mode donors|both|cells] [--backend cpu]
"""
import argparse
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import LANDSCAPE_OUT, SEED, TARGET  # noqa: E402
import program_landscape as PL                     # noqa: E402
import program_landscape_ref as REF                # noqa: E402
from program_scores import Q                       # noqa: E402
import memento                                     # noqa: E402

assert hasattr(memento, "ht_program_correlation"), "put the memento program-correlation worktree on PYTHONPATH"
OUT = LANDSCAPE_OUT / "test_memento"
OUT.mkdir(parents=True, exist_ok=True)
N_BOOT = 2000
SCHEME = "pc1"
MIN_ACTIVITY_SD = 0.2


def setup():
    progs = REF.load_programs()
    units = list(progs)
    sigma = pd.Series([progs[u][f"sigma_{SCHEME}"] for u in units], index=units)
    return progs, units, sigma


def load_cell(cell):
    """Counts, obs (donor, fine, condition) and adjusted size factors of one filled cell."""
    cells = REF.filled_cells().set_index("cell")
    pop, scale = cells.population[cell], cells.scale[cell]
    X, Y, sf, obs, genes, groups = PL.load_population(pop)
    gc = groups[groups.cell_type == cell]
    rows = np.concatenate(gc.rows.values)
    o = pd.DataFrame({"donor": obs.donor.values[rows], "fine": obs.fine.values[rows], "condition": cell,
                      "memento_size_factor": sf[rows] * scale})
    anchors = json.load(open(LANDSCAPE_OUT / "score_meta.json"))[pop]["anchors"]
    return X[rows], o, genes, anchors


def weight_frame(progs, units, genes):
    """genes x programs weights v / s_ref (the per-10k scale lives in the size factors)."""
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    rr = np.concatenate([gpos[progs[u]["genes"]].values for u in units])
    cc = np.concatenate([np.full(len(progs[u]["genes"]), k) for k, u in enumerate(units)])
    vv = np.concatenate([progs[u][SCHEME] / progs[u]["s_ref"] for u in units])
    W = sp.csc_matrix((vv, (rr, cc)), shape=(len(genes), len(units)))
    used = np.unique(rr)
    return pd.DataFrame.sparse.from_spmatrix(W[used], index=genes[used], columns=units).sparse.to_dense()


def memento_adata(X, obs, genes):
    """AnnData carrying the fields ht_program_correlation reads, with our size
    factors, constant q and donor x condition x subtype groups (floors already applied)."""
    obs = obs.reset_index(drop=True)
    obs.index = [f"c{i}" for i in range(len(obs))]
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=pd.DataFrame(index=genes))
    labels = ["donor", "condition", "fine"]
    a.obs["memento_group"] = "sg^" + obs.donor.astype(str) + "^" + obs.condition + "^" + obs.fine.astype(str)
    groups = a.obs.memento_group.drop_duplicates().tolist()
    code = a.obs.memento_group.values
    a.uns["memento"] = dict(label_columns=labels, label_delimiter="^", groups=groups,
                            group_cells={g: a.X[code == g].tocsc() for g in groups},
                            group_q=dict.fromkeys(groups, Q), estimator_type="hyper_relative")
    return a


def test(a, W, sigma, refs, contrast, args, seed=SEED):
    r = memento.ht_program_correlation(
        a, W, [TARGET], "condition", None if args.mode == "cells" else "donor", contrast=contrast,
        references=refs, scale=sigma, num_boot=args.boot, resample_cells=args.mode != "donors",
        random_state=seed, backend=args.backend)
    return r


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    o = np.argsort(p[ok])
    r = p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1)
    q[np.where(ok)[0][o]] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def active(cell, units):
    full = pd.read_csv(LANDSCAPE_OUT / "coupling.csv")
    full = full[(full.scheme == SCHEME) & (full.cell == cell)].set_index("unit").reindex(units)
    return full.activity_sd.values >= MIN_ACTIVITY_SD, full


def wide(r, kind, condition):
    x = r[(r.kind == kind) & (r.condition == condition)]
    return x.pivot(index="program", columns="statistic", values=["coef", "se", "ci_lo", "ci_hi", "pval"])


def cmd_compare(args):
    progs, units, sigma = setup()
    Xa, oa, genes, anc_a = load_cell(args.a)
    Xb, ob, gb, anc_b = load_cell(args.b)
    assert genes.equals(gb)
    if args.shared:   # donors as strata: keep only donors that carry both cells
        both = set(oa.donor) & set(ob.donor)
        ka, kb = oa.donor.isin(both).values, ob.donor.isin(both).values
        Xa, oa, Xb, ob = Xa[ka], oa[ka], Xb[kb], ob[kb]
    a = memento_adata(sp.vstack([Xa, Xb]), pd.concat([oa, ob]), genes)
    t0 = time.time()
    r = test(a, weight_frame(progs, units, genes), sigma, {args.a: anc_a, args.b: anc_b}, (args.a, args.b), args)
    secs = time.time() - t0
    out = pd.DataFrame({"unit": units})
    for cell, tag in ((args.a, "a"), (args.b, "b")):
        act, full = active(cell, units)
        w = wide(r, "estimate", cell).reindex(units)
        # the observed row must reproduce the scored coupling
        assert args.shared or np.allclose(w[("coef", "coef")].values, full.beta.values, rtol=1e-4, atol=1e-6,
                                          equal_nan=True), cell
        out[f"active_{tag}"] = act
        out[f"beta_{tag}"], out[f"excess_{tag}"] = w[("coef", "coef")].values, w[("coef", "excess")].values
    w = wide(r, "contrast", f"{args.a} - {args.b}").reindex(units)
    for stat, name in (("coef", "d_beta"), ("excess", "d_excess")):
        out[name] = w[("coef", stat)].values
        for k in ("se", "ci_lo", "ci_hi", "pval"):
            out[f"{name}_{k}"] = w[(k, stat)].values
        out[f"{name}_q"] = bh(out[f"{name}_pval"])
    tag = ("__shared" if args.shared else "") + ("" if args.mode == "donors" else f"__{args.mode}")  # donors-mode files keep the bare name
    path = OUT / f"{args.a}__vs__{args.b}{tag}.csv"
    out.to_csv(path, index=False)
    either = out[out.active_a | out.active_b]
    print(f"{args.a} vs {args.b}: {a.n_obs} cells, {len(a.uns['memento']['groups'])} groups, "
          f"{r.attrs['memento']['n_replicates']} replicates, {args.boot} draws, resampled: {args.mode}, "
          f"{args.backend}, {secs:.0f}s -> {path}")
    for k in ("d_beta", "d_excess"):
        print(f"{k}: {len(either)} programs active in either; q < 0.05: {(either[f'{k}_q'] < 0.05).sum()}, "
              f"p < 0.05: {(either[f'{k}_pval'] < 0.05).sum()}; median se {either[f'{k}_se'].median():.4f}")


def cmd_estimate(args):
    """One cell: coupling and excess over its anchors, with bootstrap SE and p."""
    progs, units, sigma = setup()
    if args.units:   # a pre-specified subset; q is then over this subset
        units = [u for u in units if u in set(args.units)]
        sigma = sigma[units]
    X, o, genes, anchors = load_cell(args.cell)
    a = memento_adata(X, o, genes)
    r = test(a, weight_frame(progs, units, genes), sigma, anchors, None, args)
    act, full = active(args.cell, units)
    w = wide(r, "estimate", args.cell).reindex(units)
    assert np.allclose(w[("coef", "coef")].values, full.beta.values, rtol=1e-4, atol=1e-6, equal_nan=True)
    out = pd.DataFrame({"unit": units, "active": act, "activity_sd": full.activity_sd.values,
                        "activity_mean": full.activity_mean.values, "rho": full.rho.values, "pct": full.pct.values})
    for stat, name in (("coef", "beta"), ("excess", "excess")):
        out[name] = w[("coef", stat)].values
        for k in ("se", "pval"):
            out[f"{name}_{k}"] = w[(k, stat)].values
        out[f"{name}_q"] = np.nan
        out.loc[act, f"{name}_q"] = bh(out.loc[act, f"{name}_pval"])
    path = OUT / f"{args.cell}__estimate{'__' + args.tag if args.tag else ''}__{args.mode}.csv"
    out.to_csv(path, index=False)
    e = out[out.active]
    print(f"{args.cell}: {a.n_obs} cells, {o.donor.nunique()} donors, {len(e)} active programs; excess q<0.05 "
          f"{(e.excess_q < .05).sum()}, beta q<0.05 {(e.beta_q < .05).sum()}, median se {e.excess_se.median():.4f}, "
          f"var(excess)/mean se^2 {e.excess.var() / (e.excess_se ** 2).mean():.2f} -> {path}")


COLLAGENS = None   # filled from the matrisome collagen program


def cmd_genes(args):
    """ANTXR2 vs individual genes (default: matrisome collagens) in one cell:
    noise-corrected correlation pooled within donor x subtype, its excess over
    the cell's 200 matched anchors (each anchor's own correlation with the same
    gene), ANTXR2's percentile among them, cell-bootstrap SE and p."""
    from program_config import FILTER_MEAN_THRESH
    from program_scores import load_library
    lib = load_library()
    genes_req = args.genes or sorted(g for g in lib[lib.program_id == "matrisome:NABA_COLLAGENS"].gene
                                     if not g.startswith("COMPLEX:"))
    X, o, genes, anchors = load_cell(args.cell)
    raw = np.asarray(X.mean(0)).ravel()
    det = np.asarray((X > 0).mean(0)).ravel()
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    have = [g for g in genes_req if g in gpos.index]
    keep = [g for g in have if raw[gpos[g]] >= FILTER_MEAN_THRESH]
    W = pd.DataFrame(np.eye(len(keep)), index=keep, columns=keep)
    a = memento_adata(X, o, genes)
    rep = None if args.mode == "cells" else "donor"
    r = memento.ht_program_correlation(a, W, [TARGET], "condition", rep, references=list(anchors),
                                       num_boot=args.boot, resample_cells=args.mode != "donors",
                                       random_state=SEED, backend=args.backend)
    # the anchors' own correlations with each gene (observed values only)
    ra = memento.ht_program_correlation(a, W, list(anchors), "condition", rep, num_boot=2,
                                        resample_cells=args.mode != "donors", random_state=SEED,
                                        backend=args.backend)
    e = r[r.kind == "estimate"]
    t = e.pivot(index="program", columns="statistic", values=["coef", "se", "pval"]).reindex(keep)
    anc = ra[ra.kind == "estimate"].pivot(index="program", columns="target", values="coef").reindex(keep)
    corr = pd.DataFrame({TARGET: t[("coef", "coef")].values}, index=keep)
    out = pd.DataFrame({
        "gene": keep, "raw_mean": raw[gpos[keep]], "det": det[gpos[keep]],
        "corr": t[("coef", "coef")].values, "corr_se": t[("se", "coef")].values, "corr_p": t[("pval", "coef")].values,
        "anchor_mean": anc.mean(1).values, "anchor_sd": anc.std(1).values,
        "excess": t[("coef", "excess")].values, "excess_se": t[("se", "excess")].values,
        "excess_p": t[("pval", "excess")].values,
        "pct": ((anc.values < corr[[TARGET]].values).sum(1) / np.isfinite(anc.values).sum(1) * 100)})
    out["corr_q"], out["excess_q"] = bh(out.corr_p), bh(out.excess_p)
    path = OUT / f"{args.cell}__genes__{args.mode}.csv"
    out.to_csv(path, index=False)
    print(f"{args.cell}: {a.n_obs} cells, {o.donor.nunique()} donors, ANTXR2 raw mean {raw[gpos[TARGET]]:.3f}; "
          f"{len(keep)} of {len(have)} genes clear the floor -> {path}")
    print(out.round(4).to_string(index=False))


def cmd_calibrate(args):
    progs, units, sigma = setup()
    rng = np.random.default_rng(SEED)
    rows = []
    for cell in args.cells:
        X, o, genes, anchors = load_cell(cell)
        W = weight_frame(progs, units, genes)
        act, _ = active(cell, units)
        donors = np.array(sorted(o.donor.unique()))
        for i in range(args.splits):
            half = set(rng.permutation(donors)[: len(donors) // 2])
            oo = o.assign(condition=np.where(o.donor.isin(half), "h1", "h2"))
            r = test(memento_adata(X, oo, genes), W, sigma, anchors, ("h1", "h2"), args, seed=SEED + i)
            w = wide(r, "contrast", "h1 - h2").reindex(units)[act]
            rec = dict(cell=cell, split=i, donors=len(donors), programs=int(act.sum()),
                       mode=args.mode)
            for stat, name in (("coef", "d_beta"), ("excess", "d_excess")):
                p, est, se = w[("pval", stat)], w[("coef", stat)], w[("se", stat)]
                rec.update({f"{name}_p05": (p < 0.05).mean(), f"{name}_p01": (p < 0.01).mean(),
                            f"{name}_q05": (bh(p) < 0.05).mean(), f"{name}_z_sd": (est / se).std(),
                            f"{name}_se": se.median()})
            rows.append(rec)
        d = pd.DataFrame([x for x in rows if x["cell"] == cell])
        print(cell, d.drop(columns=["cell", "split"]).median(numeric_only=True).round(4).to_dict(),
              "max p05:", d.d_excess_p05.max().round(3), flush=True)
    out = pd.DataFrame(rows)
    path = OUT / "calibration.csv"
    if path.exists():
        old = pd.read_csv(path)
        same = old.cell.isin(args.cells) & (old["mode"] == args.mode)
        out = pd.concat([old[~same], out], ignore_index=True)
    out.to_csv(path, index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--a", required=True)
    c.add_argument("--b", required=True)
    c.add_argument("--shared", action="store_true", help="keep only donors present in both cells")
    e = sub.add_parser("estimate")
    e.add_argument("--cell", required=True)
    e.add_argument("--units", nargs="*", help="only these programs (q over them)")
    e.add_argument("--tag", default="", help="file-name tag for a --units run")
    gn = sub.add_parser("genes")
    gn.add_argument("--cell", required=True)
    gn.add_argument("--genes", nargs="*")
    k = sub.add_parser("calibrate")
    k.add_argument("--cells", nargs="+", required=True)
    k.add_argument("--splits", type=int, default=20)
    for p in (c, k, e, gn):
        p.add_argument("--boot", type=int, default=N_BOOT)
        p.add_argument("--mode", default="cells", choices=["donors", "both", "cells"],
                       help="what is resampled: cells within donor x subtype strata, all donors (default, agreed 2026-11-15), donors, or both")
        p.add_argument("--backend", default="gpu", choices=["gpu", "cpu"])
    a = ap.parse_args()
    {"compare": cmd_compare, "calibrate": cmd_calibrate, "estimate": cmd_estimate, "genes": cmd_genes}[a.cmd](a)
