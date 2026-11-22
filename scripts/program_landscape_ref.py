"""Reference-scaled ANTXR2 x program coupling: every program gets a number in
every population x state cell, near 0 where the program is not expressed
(planning_summaries/program_coexpression.md, "Reference-scaled coupling").

The landscape (program_landscape.py) scores a program only where >= 80% of a
matrix-wide core clears the floor, so lineage-restricted programs (collagens)
never appear and a shallower atlas loses most programs. Here:

    weights  a gene enters a program if it clears the floor in >=
             REF_MIN_CELLS filled cells. Its reference SD s_g is the 90th
             percentile, over those cells, of its corrected SD in per-10k
             units (the gene where it varies most). PC1 (and
             equal) weights come from the gene-gene correlation averaged over
             the cells where each pair is measured (0 if never). The
             program's reference SD is sigma_ref = sqrt(w' Rbar w).
    score    in a cell the program score is S = sum_g w_g Y_g / (scale s_g)
             over ALL its genes, so a gene silent here adds ~nothing.
                 beta  = cov(T, S) / (sd(T) * sigma_ref)
                       = rho * activity_sd,  activity_sd = sd(S) / sigma_ref
             for T = ANTXR2 and N_ANCHORS (200) span-matched anchors;
             excess = beta_ANTXR2 - mean beta_anchor, z = excess / anchor SD,
             pct = share of anchors below ANTXR2. Alongside: activity_mean
             (mean per-10k expression over the program's genes) and frac_pass
             (share of them clearing the floor here).

ANTXR2 must still clear its own floor in a cell (its SD is in beta). No
bootstrap yet: point estimates and the anchor distribution only.

ANTXR2_LANDSCAPE=gca: IBDverse weights and Rbar (never re-fitted), GCA's own
reference SDs (genes clearing the floor in >= 1 GCA cell; sigma_ref over
those), ANTXR2 floor 0.05.

Outputs in LANDSCAPE_OUT (landscape_ref/ or landscape_gca_ref/): cells.csv,
prep/, weights.csv, programs.csv, rbar.npz, gene_ref.csv, coupling.csv,
score_meta.json.

Usage:
    ANTXR2_LANDSCAPE_RUN=ref python program_landscape_ref.py {prep,weights,score,all}
        [--pops SLUG ...]
"""
import argparse
import json
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import (                       # noqa: E402
    LANDSCAPE_OUT, LANDSCAPE_RUN, LANDSCAPE_WEIGHTS_FROM, MIN_GENES, MIN_GENES_DISCOVERY,
    N_ANCHORS, REF_MIN_CELLS, TARGET, gene_spans, pop_slug,
)
assert LANDSCAPE_RUN == "ref", "run with ANTXR2_LANDSCAPE_RUN=ref"
import program_landscape as PL                     # noqa: E402
from program_scores import gene_stats, group_stats, load_library  # noqa: E402

LOW_FLOOR = 0.01              # anchor candidates: raw mean floor in every state
SCHEMES = ("pc1", "equal")
REF_QUANTILE = 0.9


def filled_cells():
    c = pd.read_csv(LANDSCAPE_OUT / "cells.csv")
    return c[c.filled]


def gene_reference(min_cells):
    """Reference SD per gene (per-10k units): REF_QUANTILE of its corrected SD
    over the filled cells where it clears the floor, i.e. the gene where it
    varies most (a median would let a gene dwarf its reference in the one
    lineage that expresses it highly); genes in >= min_cells cells."""
    sd = {}
    for c, scale in filled_cells()[["cell", "scale"]].itertuples(index=False):
        z, genes = PL.load_prep(c)
        for g, v in zip(genes[z["idx"]], np.sqrt(z["var"]) / scale):
            sd.setdefault(g, []).append(v)
    ref = pd.DataFrame({"s_ref": {g: np.quantile(v, REF_QUANTILE) for g, v in sd.items()},
                        "n_pass": {g: len(v) for g, v in sd.items()}})
    ref = ref[ref.n_pass >= min_cells].rename_axis("gene")
    ref.to_csv(LANDSCAPE_OUT / "gene_ref.csv")
    return ref


def fit_weights(lib):
    ref = gene_reference(REF_MIN_CELLS)
    members = {}
    for pid, df in lib.groupby("program_id", sort=True):
        k = df[df.gene.isin(ref.index) & (df.gene != TARGET)].drop_duplicates("gene")
        fixed = pid.startswith(("discovery:", "bracq:"))
        if len(k) >= (MIN_GENES_DISCOVERY if fixed else MIN_GENES):
            members[pid] = k
    U = sorted(set().union(*[set(k.gene) for k in members.values()]))
    upos = pd.Series(np.arange(len(U)), index=U)
    Ssum = np.zeros((len(U), len(U)), np.float32)
    Scnt = np.zeros((len(U), len(U)), np.float32)
    for c in filled_cells().cell:
        z, genes = PL.load_prep(c)
        rg = genes[z["ridx"]]
        keep = np.where(rg.isin(U))[0]
        u = upos[rg[keep]].values
        Ssum[np.ix_(u, u)] += z["R"][np.ix_(keep, keep)]
        Scnt[np.ix_(u, u)] += 1
    wrows, prows, rbar = [], [], {}
    for i, (pid, k) in enumerate(members.items()):
        u = upos[k.gene].values
        with np.errstate(invalid="ignore"):
            Rbar = np.nan_to_num(Ssum[np.ix_(u, u)] / Scnt[np.ix_(u, u)]).astype(np.float64)
        np.fill_diagonal(Rbar, 1.0)
        e = k.sign.values.astype(float)
        lam, V = np.linalg.eigh(Rbar)
        v1 = V[:, -1] if V[:, -1] @ Rbar @ e >= 0 else -V[:, -1]
        eq = e / np.linalg.norm(e)
        rbar[f"p{i}"] = Rbar.astype(np.float32)
        prows.append(dict(program_id=pid, collection=k.collection.iloc[0], rbar_key=f"p{i}",
                          n_genes=len(k), pc1_ve=lam[-1] / len(k), cos_pc1_equal=float(v1 @ eq),
                          never_measured_pairs=float((Scnt[np.ix_(u, u)] == 0).mean())))
        wrows.append(pd.DataFrame(dict(program_id=pid, gene=k.gene.values, sign=e, v_equal=eq,
                                       v_pc1=v1, n_pass=ref.n_pass[k.gene].values)))
    pd.DataFrame(prows).to_csv(LANDSCAPE_OUT / "programs.csv", index=False)
    pd.concat(wrows).to_csv(LANDSCAPE_OUT / "weights.csv", index=False)
    np.savez_compressed(LANDSCAPE_OUT / "rbar.npz", **rbar)
    print(f"weights: {len(members)} programs, {len(U)} genes", flush=True)


def weights(lib):
    if LANDSCAPE_WEIGHTS_FROM is None:
        fit_weights(lib)
    else:                                   # a replication atlas: fixed weights
        for f in ("weights.csv", "programs.csv", "rbar.npz"):
            shutil.copy(LANDSCAPE_WEIGHTS_FROM / f, LANDSCAPE_OUT / f)
        gene_reference(1)


def load_programs():
    """Per program: genes, weights per scheme restricted to genes with a
    reference SD in this atlas, and sigma_ref over those."""
    W = pd.read_csv(LANDSCAPE_OUT / "weights.csv")
    P = pd.read_csv(LANDSCAPE_OUT / "programs.csv").set_index("program_id")
    ref = pd.read_csv(LANDSCAPE_OUT / "gene_ref.csv", index_col=0).s_ref
    rbar = np.load(LANDSCAPE_OUT / "rbar.npz")
    out = {}
    for pid, w in W.groupby("program_id", sort=False):
        ok = w.gene.isin(ref.index).values
        if not ok.any():
            continue
        R = rbar[P.rbar_key[pid]][np.ix_(ok, ok)].astype(np.float64)
        d = dict(genes_all=w.gene.values, genes=w.gene.values[ok], s_ref=ref[w.gene[ok]].values)
        for s in SCHEMES:
            v = w[f"v_{s}"].values[ok]
            d[s] = v
            d[f"sigma_{s}"] = float(np.sqrt(max(v @ R @ v, 0)))
            d[f"w2_{s}"] = float((v ** 2).sum())      # share of the weight vector kept
        out[pid] = d
    return out


def score(pops):
    cells = filled_cells()
    progs = load_programs()
    keys = [(u, s) for u in progs for s in SCHEMES]
    sigma = np.array([progs[u][f"sigma_{s}"] for u, s in keys])
    span = np.log10(gene_spans())
    rows, meta = [], {}
    for pop in pops:
        t0 = time.time()
        fc = cells[cells.population == pop]
        if fc.empty:
            print(f"score {pop}: no filled cell", flush=True)
            continue
        X, Y, sf, obs, genes, groups = PL.load_population(pop)
        gpos = pd.Series(np.arange(len(genes)), index=genes)
        feats, cand, stats = {}, [], {}
        for c in fc.cell:
            gr = np.concatenate(groups[groups.cell_type == c].rows.values)
            raw_mean, det = gene_stats(X, gr)
            stats[c] = (raw_mean, np.asarray(Y[gr].mean(0)).ravel())
            feats[c] = pd.DataFrame({"log10_mean": np.log10(raw_mean + 1e-6), "det": det},
                                    index=genes)
            cand.append(set(genes[raw_mean >= LOW_FLOOR]))
        anchors, n_cand = PL.span_matched(feats, cand, span, N_ANCHORS)
        targets = [TARGET] + anchors
        T = gpos[targets].values
        f0 = feats[fc.cell.iloc[0]]
        meta[pop] = dict(anchors=anchors, n_candidates=n_cand,
                         antxr2_log10_mean=float(f0.log10_mean[TARGET]),
                         anchor_log10_mean_range=[float(f0.log10_mean[anchors].min()),
                                                  float(f0.log10_mean[anchors].max())])
        # targets that are members of a program are not scored against it
        member = np.array([[g in set(progs[u]["genes_all"]) for g in targets] for u, _ in keys])
        for c, state, scale in fc[["cell", "state", "scale"]].itertuples(index=False):
            gc = groups[groups.cell_type == c].reset_index(drop=True)
            raw_mean, ymean = stats[c]
            passing = set(genes[PL.load_prep(c)[0]["idx"]])
            rr, cc, vv = [], [], []
            for k, (u, s) in enumerate(keys):
                p = progs[u]
                rr.append(gpos[p["genes"]].values)
                cc.append(np.full(len(p["genes"]), k))
                vv.append(p[s] / (p["s_ref"] * scale))
            Wm = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))),
                               shape=(len(genes), len(keys)))
            (n, s_sum, s_sq, s_noise), (t_sum, t_sq, t_noise), ST = group_stats(X, Y, sf, gc, Wm, T)
            N = n.sum()
            vs = (s_sq - s_sum ** 2 / n[:, None] - s_noise).sum(0) / N
            vt = (t_sq - t_sum ** 2 / n[:, None] - t_noise).sum(0) / N
            cov = (ST - s_sum[:, :, None] * t_sum[:, None, :] / n[:, None, None]).sum(0) / N
            with np.errstate(invalid="ignore", divide="ignore"):
                sdt = np.sqrt(np.where(vt > 0, vt, np.nan))
                beta = np.where(member, np.nan, cov / (sdt[None, :] * sigma[:, None]))
                rho = cov[:, 0] / (np.sqrt(np.where(vs > 0, vs, np.nan)) * sdt[0])
                am, asd = np.nanmean(beta[:, 1:], 1), np.nanstd(beta[:, 1:], 1)
                pct = (np.nansum(beta[:, 1:] < beta[:, :1], 1)
                       / np.isfinite(beta[:, 1:]).sum(1) * 100)
            act_sd = np.sqrt(np.clip(vs, 0, None)) / sigma
            for k, (u, s) in enumerate(keys):
                ga = gpos[[g for g in progs[u]["genes_all"] if g in gpos.index]].values
                rows.append(dict(
                    population=pop, state=state, cell=c, unit=u, scheme=s,
                    n_genes=len(progs[u]["genes"]), weight_kept=progs[u][f"w2_{s}"],
                    frac_pass=np.mean([g in passing for g in progs[u]["genes_all"]]),
                    activity_mean=ymean[ga].mean() / scale, activity_sd=act_sd[k],
                    rho=rho[k], beta=beta[k, 0], anchor_mean=am[k], anchor_sd=asd[k],
                    excess=beta[k, 0] - am[k], excess_z=(beta[k, 0] - am[k]) / asd[k],
                    pct=pct[k]))
        print(f"score {pop}: {len(fc)} cells ({time.time() - t0:.0f}s)", flush=True)
    L = pd.DataFrame(rows)
    path = LANDSCAPE_OUT / "coupling.csv"
    if path.exists():
        old = pd.read_csv(path)
        L = pd.concat([old[~old.population.isin(pops)], L], ignore_index=True)
    L.to_csv(path, index=False)
    mpath = LANDSCAPE_OUT / "score_meta.json"
    allmeta = json.load(open(mpath)) if mpath.exists() else {}
    allmeta.update(meta)
    json.dump(allmeta, open(mpath, "w"), indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prep", "weights", "score", "all"])
    ap.add_argument("--pops", nargs="*", help="population slugs (default: all)")
    a = ap.parse_args()
    pops = [p for p in PL.POPS if not a.pops or pop_slug(p) in a.pops]
    lib = load_library()
    if a.stage in ("prep", "all"):
        PL.prep(pops, lib)
    if a.stage in ("weights", "all"):
        weights(lib)
    if a.stage in ("score", "all"):
        score(pops)


if __name__ == "__main__":
    main()
