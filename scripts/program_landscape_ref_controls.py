"""Positive-control calibration of the matched-gene (anchor) test
(user, 2026-11-18: "maybe the negative control is a bit too hard to beat").

For ANTXR2 and genes with established fibroblast ECM roles that are not
members of the matrisome collagen / core programs (MRC2, ITGA11, DDR2, FAP,
PDGFRB, LRP1, ANTXR1), each gene gets its own 200 anchors, chosen exactly as
ANTXR2's are (genomic span within ANCHOR_SPAN_WINDOW in log10, then nearest on
log10 raw mean + detection in this cell). Per program: coupling (beta), excess
over the gene's anchors with cell-bootstrap SE, and the gene's percentile
among its anchors. If the controls cannot beat their anchors on ECM programs
either, the anchor test is too strict for this purpose.

    -> <run>/test_memento/{cell}__controls.csv

Usage:
    PYTHONPATH=<memento worktree> ANTXR2_LANDSCAPE_RUN=ref [ANTXR2_LANDSCAPE=skin] \\
        python program_landscape_ref_controls.py --cell CELL [--genes ...] [--backend cpu]
"""
import argparse
import numpy as np
import pandas as pd

import program_landscape_ref_memento as M
from program_config import ANCHOR_SPAN_WINDOW, FILTER_MEAN_THRESH, gene_spans
from program_scores import match_anchors
import memento

CONTROLS = ["ANTXR2", "MRC2", "ITGA11", "DDR2", "FAP", "PDGFRB", "LRP1", "ANTXR1"]
LOW_FLOOR = 0.01
N_ANCHORS = 200


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--genes", nargs="*", default=CONTROLS)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    a = ap.parse_args()
    progs, units, sigma = M.setup()
    X, o, genes, _ = M.load_cell(a.cell)
    raw = np.asarray(X.mean(0)).ravel()
    det = np.asarray((X > 0).mean(0)).ravel()
    feats = pd.DataFrame({"log10_mean": np.log10(raw + 1e-6), "det": det}, index=genes)
    feats = feats[~feats.index.duplicated()]
    span = np.log10(gene_spans())
    span = span[~span.index.duplicated()]
    pool = [g for g in feats.index[raw[:len(feats)] >= LOW_FLOOR] if g in span.index]
    W = M.weight_frame(progs, units, genes)
    ad = M.memento_adata(X, o, genes)
    act, _ = M.active(a.cell, units)
    rows = []
    for g in a.genes:
        if g not in feats.index or 10 ** feats.log10_mean[g] < FILTER_MEAN_THRESH or g not in span.index:
            print(f"{g}: below floor or no span, skipped")
            continue
        cand = [c for c in pool if c != g and c not in a.genes and abs(span[c] - span[g]) <= ANCHOR_SPAN_WINDOW]
        anchors = match_anchors(feats, feats.loc[g], cand, N_ANCHORS)
        kw = dict(scale=sigma, resample_cells=True, random_state=M.SEED, backend=a.backend)
        r = memento.ht_program_correlation(ad, W, [g], "condition", None, references=anchors,
                                           num_boot=a.boot, **kw)
        ra = memento.ht_program_correlation(ad, W, anchors, "condition", None, num_boot=2, **kw)
        e = r[r.kind == "estimate"].pivot(index="program", columns="statistic",
                                          values=["coef", "se", "pval"]).reindex(units)
        anc = ra[ra.kind == "estimate"].pivot(index="program", columns="target", values="coef").reindex(units)
        beta = e[("coef", "coef")].values
        with np.errstate(invalid="ignore"):
            pct = (anc.values < beta[:, None]).sum(1) / np.isfinite(anc.values).sum(1) * 100
        d = pd.DataFrame({"gene": g, "unit": units, "active": act, "raw_mean": 10 ** feats.log10_mean[g],
                          "beta": beta, "beta_se": e[("se", "coef")].values, "beta_p": e[("pval", "coef")].values,
                          "anchor_mean": anc.mean(1).values, "excess": e[("coef", "excess")].values,
                          "excess_se": e[("se", "excess")].values, "excess_p": e[("pval", "excess")].values,
                          "pct": np.where(np.isfinite(beta), pct, np.nan)})
        d.loc[d.active, "excess_q"] = M.bh(d.loc[d.active, "excess_p"])
        d.loc[d.active, "beta_q"] = M.bh(d.loc[d.active, "beta_p"])
        rows.append(d)
        x = d[d.active]
        print(f"{g}: raw mean {d.raw_mean.iloc[0]:.3f}; excess q<0.05 {(x.excess_q < .05).sum()} "
              f"(+{((x.excess_q < .05) & (x.excess > 0)).sum()}), beta q<0.05 {(x.beta_q < .05).sum()}, "
              f"median excess se {x.excess_se.median():.4f}", flush=True)
    out = pd.concat(rows, ignore_index=True)
    path = M.OUT / f"{a.cell}__controls.csv"
    out.to_csv(path, index=False)
    print("->", path)


if __name__ == "__main__":
    main()
