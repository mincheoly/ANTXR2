"""ANTXR2 and interferon in stimulated airway epithelium (user, 2026-11-18: "Let's
try memento's HTECs"). Kim et al., Cell 2024 (memento), GEO GSE274751
(GSE274751_HTEC_counts.h5ad): primary human tracheal epithelial cells, 2 donors,
IFN-alpha / beta / gamma / lambda at 3, 6, 9, 24, 48 h plus unstimulated
controls (time 0).

1. mean   ANTXR2 per 10k by cell type x donor x stim x time, ratio to control
          (UMI-pooled per donor).
2. couple Within ciliated cells (the only type above the 0.07 raw-mean floor),
          noise-corrected correlation of ANTXR2 with each program (IBDverse
          weights, scale=None), per stim (strata donor x time; control: donor),
          excess over 200 anchors matched as in the landscape (span, then raw
          mean + detection in these cells), cell bootstrap; anchors' own
          correlations give ANTXR2's percentile.

    -> OUT_DIR/htec_ifn/{mean.csv, coupling_<celltype>.csv}

Usage:
    PYTHONPATH=<memento worktree> ANTXR2_LANDSCAPE_RUN=ref python program_htec_ifn.py [--boot N] [--units ...]
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import program_landscape_ref_memento as M
from program_config import ANCHOR_SPAN_WINDOW, FILTER_MEAN_THRESH, OUT_DIR, ROOT, TARGET, gene_spans
from program_scores import match_anchors
import memento

H5AD = ROOT / "raw" / "htec" / "GSE274751_HTEC_counts.h5ad"
OUT = OUT_DIR / "htec_ifn"
OUT.mkdir(parents=True, exist_ok=True)
CT = {"basal": "basal", "basal/club": "club", "ciliated": "ciliated"}
LOW_FLOOR = 0.01
N_ANCHORS = 200


def load():
    a = ad.read_h5ad(H5AD)
    o = a.obs.astype({"cell_type": str, "stim": str, "donor": str}).copy()
    o["time"] = o.time.astype(int)
    o["ct"] = o.cell_type.map(CT)
    keep = o.ct.notna().values
    X = sp.csr_matrix(a.X)[keep]
    return X, o[keep].reset_index(drop=True), np.asarray(a.var_names)


def cmd_mean(X, o, genes):
    j = list(genes).index(TARGET)
    o = o.assign(umi=np.asarray(X[:, j].todense()).ravel())
    g = o.groupby(["ct", "donor", "stim", "time"], observed=True)
    t = g.agg(cells=("umi", "size"), umis=("umi", "sum"), depth=("total_counts", "sum")).reset_index()
    t["per10k"] = t.umis / t.depth * 1e4
    ctl = t[t.stim == "control"].set_index(["ct", "donor"]).per10k
    t["ratio"] = t.per10k / ctl.reindex(pd.MultiIndex.from_frame(t[["ct", "donor"]])).values
    t.to_csv(OUT / "mean.csv", index=False)
    print(t[t.stim != "control"].pivot_table(index=["ct", "stim"], columns="time", values="ratio").round(2).to_string())


def cmd_couple(X, o, genes, ct, boot, units_req, backend):
    m = (o.ct == ct).values
    X, o = X[m], o[m].reset_index(drop=True)
    # memento_adata groups by donor x condition x fine: condition = stim, fine = time
    obs = pd.DataFrame({"donor": o.donor.values, "condition": o.stim.values, "fine": o.time.astype(str).values,
                        "memento_size_factor": o.total_counts.values / 1e4})
    out = couple(X, obs, genes, ["control", "alpha", "beta", "gamma", "lambda"], boot, units_req, backend, ct)
    tag = "" if not units_req else "_subset"
    out.to_csv(OUT / f"coupling_{ct}{tag}.csv", index=False)
    print(f"-> {OUT / f'coupling_{ct}{tag}.csv'}")
    return out


def couple(X, obs, genes, conditions, boot, units_req, backend, label, anchors=None, extra=None):
    """ANTXR2 noise-corrected correlation with each program per condition (strata
    donor x fine within it), excess over 200 anchors matched in these cells (or
    the given `anchors`), and ANTXR2's percentile among the anchors' own
    correlations. `extra`: {name: genes} gene sets added as equal-weight programs
    (each gene scaled by its reference SD). Also returns each program's mean and
    SD of per-cell activity (sum_g w_g x_g / (size factor s_ref_g); size factors
    are per 10k)."""
    raw = np.asarray(X.mean(0)).ravel()
    det = np.asarray((X > 0).mean(0)).ravel()
    feats = pd.DataFrame({"log10_mean": np.log10(raw + 1e-6), "det": det}, index=genes)
    feats = feats[~feats.index.duplicated()]
    print(f"{label}: {X.shape[0]} cells, ANTXR2 raw mean {10 ** feats.log10_mean[TARGET]:.3f} (floor {FILTER_MEAN_THRESH})")
    if anchors is None:
        span = np.log10(gene_spans())
        span = span[~span.index.duplicated()]
        pool = [g for g in feats.index[(10 ** feats.log10_mean >= LOW_FLOOR).values]
                if g in span.index and g != TARGET and abs(span[g] - span[TARGET]) <= ANCHOR_SPAN_WINDOW]
        anchors = match_anchors(feats, feats.loc[TARGET], pool, N_ANCHORS)
    else:
        anchors = [g for g in anchors if g in feats.index]

    progs, units, _ = M.setup()
    if units_req:
        units = [u for u in units if u in set(units_req)]
    if extra:   # a list of genes (equal weights) or a Series gene -> weight
        sref = pd.read_csv(M.LANDSCAPE_OUT / "gene_ref.csv", index_col=0).s_ref
        for name, gs in extra.items():
            if not isinstance(gs, pd.Series):
                gs = [g for g in dict.fromkeys(gs) if g in sref.index and g != TARGET]
                gs = pd.Series(1 / np.sqrt(len(gs)), index=gs)
            gs = gs[gs.index.isin(sref.index) & (gs.index != TARGET)]
            progs[name] = {"genes": gs.index.values, "s_ref": sref[gs.index].values, M.SCHEME: gs.values}
            units.append(name)
    # keep each program's genes present in this matrix (weights of the rest dropped)
    gset = set(genes)
    for u in units:
        ok = np.array([g in gset for g in progs[u]["genes"]])
        if not ok.all():
            progs[u] = dict(progs[u], genes=progs[u]["genes"][ok], s_ref=progs[u]["s_ref"][ok],
                            **{M.SCHEME: progs[u][M.SCHEME][ok]})
    units = [u for u in units if len(progs[u]["genes"]) >= 3]
    W = M.weight_frame(progs, units, genes)
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    gpos = gpos[~gpos.index.duplicated()]
    Y = sp.diags(1 / obs.memento_size_factor.values) @ X[:, gpos[W.index].values]
    act = pd.DataFrame(np.asarray(Y @ W.values), columns=units)
    act["condition"] = obs.condition.values
    a = M.memento_adata(X, obs, genes)
    kw = dict(scale=None, resample_cells=True, random_state=M.SEED, backend=backend)
    r = memento.ht_program_correlation(a, W, [TARGET], "condition", None, references=list(anchors),
                                       num_boot=boot, **kw)
    ra = memento.ht_program_correlation(a, W, list(anchors), "condition", None, num_boot=2, **kw)
    rows = []
    for cond in conditions:
        e = M.wide(r, "estimate", cond).reindex(units)
        anc = ra[(ra.kind == "estimate") & (ra.condition == cond)].pivot(
            index="program", columns="target", values="coef").reindex(units)
        rho = e[("coef", "coef")].values
        with np.errstate(invalid="ignore"):
            pct = (anc.values < rho[:, None]).sum(1) / np.isfinite(anc.values).sum(1) * 100
        d = pd.DataFrame({"cell_type": label, "stim": cond, "unit": units, "rho": rho,
                          "rho_se": e[("se", "coef")].values, "rho_p": e[("pval", "coef")].values,
                          "anchor_mean": anc.mean(1).values, "anchor_sd": anc.std(1).values,
                          "excess": e[("coef", "excess")].values, "excess_se": e[("se", "excess")].values,
                          "excess_p": e[("pval", "excess")].values, "pct": pct})
        d["excess_q"] = M.bh(d.excess_p)
        ac = act[act.condition == cond][units]
        d["activity_mean"], d["activity_sd"] = ac.mean().values, ac.std().values
        rows.append(d)
    out = pd.concat(rows, ignore_index=True)
    out.attrs["activity"] = act
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["mean", "couple"])
    ap.add_argument("--cell-type", default="ciliated")
    ap.add_argument("--units", nargs="*")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    a = ap.parse_args()
    X, o, genes = load()
    if a.cmd == "mean":
        cmd_mean(X, o, genes)
    else:
        cmd_couple(X, o, genes, a.cell_type, a.boot, a.units, a.backend)


if __name__ == "__main__":
    main()
