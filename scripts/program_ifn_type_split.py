"""Which interferon does ANTXR2 track? (user, 2026-11-21: "do the gene set
separation ... and fibroblast type 1 interferon next").

1. sets        Type I-biased / IFN-gamma-biased / shared gene sets from the HTEC
               stimulation (GSE274751): pseudobulk per cell type x donor x stim
               (6-48 h pooled) vs unstimulated, median log2FC over the 6
               cell type x donor pairs (pseudocount 0.1 per 10k).
               - type I-biased: alpha and beta each >= 2, mean(alpha, beta) -
                 gamma >= 1.5, >= 0.5 per 10k under alpha and beta.
               - gamma-biased: gamma >= 1, gamma - max(alpha, beta) >= 0.8,
                 >= 0.2 per 10k under gamma. (In airway epithelium gamma also
                 induces the ISGF3 genes ~2x, so strict "only" sets are 2-15
                 genes; bias is the usable split.)
               - shared: alpha, beta, gamma each >= 1.5, |mean(alpha, beta) -
                 gamma| < 0.8.
2. tissue      Per landscape cell (ANTXR2_LANDSCAPE=ibdverse|gca): ANTXR2's
               noise-corrected correlation (scale None) with the three sets and
               the 8 pre-specified IFN programs, excess over the population's
               landscape anchors, cell bootstrap within donor x subtype strata
               (program_htec_ifn.couple); per-cell set activity (mean, SD) and
               the correlation between the sets' activities.
3. fibroblast  Dermal fibroblasts, Hagai et al. 2018:
               - dsrna: 10x, unstimulated (E-MTAB-5988) vs poly I:C 6 h
                 (E-MTAB-5989); 3 donors pooled, not demultiplexed.
               - ifnb:  Smart-seq2, IFN-beta 2 h / 6 h, 5 donors
                 (E-MTAB-7051); no unstimulated cells. Read counts (Salmon,
                 rounded), not UMIs: the noise correction is near zero at this
                 depth, so rho is close to a plain correlation; the anchors
                 carry the technical baseline.
               ANTXR2 per 10k by condition; coupling as in 2, anchors matched
               in these cells.

    -> OUT_DIR/ifn_type_split/{sets.csv, htec_lfc.csv, tissue/<atlas>/<cell>.csv,
       fibroblast/{dsrna,ifnb}.csv, fibroblast/mean.csv}

Usage:
    PYTHONPATH=<memento worktree> ANTXR2_LANDSCAPE_RUN=ref python program_ifn_type_split.py sets
    ... tissue --cells CELL [CELL ...] [--boot N]
    ... fibroblast [--boot N]
"""
import argparse

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

from program_config import OUT_DIR, ROOT, TARGET
import program_htec_ifn as H
import program_landscape_ref_memento as M

OUT = OUT_DIR / "ifn_type_split"
OUT.mkdir(parents=True, exist_ok=True)
IFN8 = ["hallmark:INTERFERON_ALPHA_RESPONSE", "hallmark:INTERFERON_GAMMA_RESPONSE", "cytosig:IFN1", "cytosig:IFNG",
        "collectri:IRF1", "collectri:IRF2", "collectri:STAT1", "progeny:JAK-STAT"]
SETS = ["IFN_I_biased", "IFN_G_biased", "IFN_shared"]
RESID = ["I_resid_G", "G_resid_I"]
PC = 0.1
FIB_DSRNA = ROOT / "raw" / "fib_dsrna_hagai"
FIB_IFNB = ROOT / "raw" / "fib_ifnb_hagai"
FIB_MIN_UMI = 1000
FIB_MAX_MT = 0.2


def htec_pseudobulk():
    path = OUT / "htec_pseudobulk.parquet"
    if path.exists():
        return pd.read_parquet(path)
    X, o, genes = H.load()
    late = np.where(o.stim == "control", "0", np.where(o.time >= 6, "late", "3"))
    key = (o.ct + "|" + o.donor + "|" + o.stim + "|" + late).values
    keys = pd.unique(key)
    G = sp.csr_matrix((np.ones(len(o)), (pd.Categorical(key, categories=keys).codes, np.arange(len(o)))),
                      shape=(len(keys), len(o)))
    tot = np.asarray(G @ o.total_counts.values).ravel()
    pb = pd.DataFrame(np.asarray((G @ X).todense()) / tot[:, None] * 1e4, index=keys, columns=genes)
    pb = pb.loc[:, ~pb.columns.duplicated()]
    pb.to_parquet(path)
    return pb


def cmd_sets():
    pb = htec_pseudobulk()
    idx = pb.index.to_series().str.split("|", expand=True)
    idx.columns = ["ct", "donor", "stim", "t"]
    lfc = {}
    for s in ["alpha", "beta", "gamma", "lambda"]:
        cols = []
        for (ct, d), g in idx.groupby(["ct", "donor"]):
            c = pb.loc[g.index[g.stim == "control"][0]]
            t = pb.loc[g.index[(g.stim == s) & (g.t == "late")][0]]
            cols.append(np.log2((t + PC) / (c + PC)))
        lfc[s] = pd.concat(cols, axis=1).median(1)
    med = pd.DataFrame(lfc)
    expr = pd.DataFrame({s: pb.loc[(idx.stim == s).values & (idx.t == "late").values].median()
                         for s in ["alpha", "beta", "gamma"]})
    one = med[["alpha", "beta"]].mean(1)
    sets = {
        "IFN_I_biased": (med[["alpha", "beta"]].min(1) >= 2) & (one - med.gamma >= 1.5)
                        & (expr[["alpha", "beta"]].min(1) >= .5),
        "IFN_G_biased": (med.gamma >= 1) & (med.gamma - med[["alpha", "beta"]].max(1) >= .8) & (expr.gamma >= .2),
        "IFN_shared": (med[["alpha", "beta", "gamma"]].min(1) >= 1.5) & ((one - med.gamma).abs() < .8)
                      & (expr[["alpha", "beta", "gamma"]].min(1) >= .5),
    }
    rows = [pd.DataFrame({"set": k, "gene": med.index[v.values]}) for k, v in sets.items()]
    t = pd.concat(rows, ignore_index=True)
    t = t[t.gene != TARGET]
    t.to_csv(OUT / "sets.csv", index=False)
    med.assign(expr_alpha=expr.alpha, expr_beta=expr.beta, expr_gamma=expr.gamma).to_csv(OUT / "htec_lfc.csv")
    for k, g in t.groupby("set"):
        print(f"{k} ({len(g)}): {' '.join(g.gene)}")
    print(med.loc[[TARGET]].round(2).to_string())


def load_sets():
    t = pd.read_csv(OUT / "sets.csv")
    return {k: list(t.gene[t.set == k]) for k in SETS}


def residual_sets(X, obs, genes, sets):
    """Add I_resid_G = I - b G and G_resid_I = G - b' I (b from the within-stratum
    regression of one set's per-cell activity on the other's, Poisson noise
    removed from the regressor's variance). Both stay linear in the genes, so
    memento scores them like any program; ANTXR2's coupling to them is the part
    of each response the other does not explain."""
    sref = pd.read_csv(M.LANDSCAPE_OUT / "gene_ref.csv", index_col=0).s_ref
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    gpos = gpos[~gpos.index.duplicated()]
    sf = obs.memento_size_factor.values
    grp = (obs.donor.astype(str) + "^" + obs.condition.astype(str) + "^" + obs.fine.astype(str)).values
    v, A, N = {}, {}, {}
    for k in ("IFN_I_biased", "IFN_G_biased"):
        gs = [g for g in dict.fromkeys(sets[k]) if g in sref.index and g in gpos.index and g != TARGET]
        v[k] = pd.Series(1 / np.sqrt(len(gs)), index=gs)
        w = (v[k] / sref[gs]).values
        Xk = X[:, gpos[gs].values]
        A[k] = np.asarray(Xk @ w).ravel() / sf
        A[k] = A[k] - pd.Series(A[k]).groupby(grp).transform("mean").values
        N[k] = (np.asarray(Xk @ w ** 2).ravel() / sf ** 2).mean()
    aI, aG = A["IFN_I_biased"], A["IFN_G_biased"]
    cov = (aI * aG).mean()
    vI = max((aI ** 2).mean() - N["IFN_I_biased"], 1e-3 * (aI ** 2).mean())
    vG = max((aG ** 2).mean() - N["IFN_G_biased"], 1e-3 * (aG ** 2).mean())
    info = dict(b_I_on_G=cov / vG, b_G_on_I=cov / vI, r_IG=cov / np.sqrt(vI * vG),
                noise_share_I=N["IFN_I_biased"] / (aI ** 2).mean(), noise_share_G=N["IFN_G_biased"] / (aG ** 2).mean())
    I, G = v["IFN_I_biased"], v["IFN_G_biased"]
    extra = dict(sets)
    extra["I_resid_G"] = I.add(-info["b_I_on_G"] * G, fill_value=0)
    extra["G_resid_I"] = G.add(-info["b_G_on_I"] * I, fill_value=0)
    print("residualization: " + ", ".join(f"{k} {x:.3f}" for k, x in info.items()), flush=True)
    return extra, info


def summarize(out, label):
    act = out.attrs["activity"]
    cols = ["unit", "rho", "rho_se", "excess", "excess_se", "excess_q", "pct", "activity_mean", "activity_sd"]
    for cond, d in out.groupby("stim", sort=False):
        print(f"--- {label} / {cond}")
        print(d[cols].round(3).to_string(index=False))
        a = act[act.condition == cond][[s for s in SETS + RESID if s in act]]
        print("set activity correlation (raw, per cell):")
        print(a.corr().round(2).to_string())


def cmd_tissue(cells, boot, backend):
    sets = load_sets()
    d = OUT / "tissue" / M.LANDSCAPE_OUT.name
    d.mkdir(parents=True, exist_ok=True)
    for cell in cells:
        X, o, genes, anchors = M.load_cell(cell)
        extra, info = residual_sets(X, o, genes, sets)
        out = H.couple(X, o, genes, [cell], boot, IFN8, backend, cell, anchors=anchors, extra=extra)
        out = out.assign(**info)
        j = list(genes).index(TARGET)
        out["antxr2_per10k"] = float((X[:, j].toarray().ravel() / o.memento_size_factor.values).mean())
        out["n_cells"], out["n_donors"] = X.shape[0], o.donor.nunique()
        out.to_csv(d / f"{cell}.csv", index=False)
        summarize(out, cell)


def cmd_htec(boot, backend, ct="ciliated"):
    """Check on the source data: the residual programs, fitted per stimulus
    (b from that stimulus' cells), should carry ANTXR2's coupling under alpha /
    beta and nothing under gamma. Circular for the raw sets (they were chosen
    here), not for the question which residual carries the coupling."""
    sets = load_sets()
    X, o, genes = H.load()
    m = (o.ct == ct).values
    X, o = X[m], o[m].reset_index(drop=True)
    rows = []
    for stim in ["alpha", "beta", "gamma", "lambda"]:
        k = (o.stim == stim).values
        obs = pd.DataFrame({"donor": o.donor.values[k], "condition": stim, "fine": o.time.astype(str).values[k],
                            "memento_size_factor": o.total_counts.values[k] / 1e4})
        extra, info = residual_sets(X[k], obs, genes, sets)
        out = H.couple(X[k], obs, genes, [stim], boot, IFN8, backend, f"htec_{ct}_{stim}", extra=extra)
        summarize(out, f"htec {ct}")
        out.attrs = {}   # the activity frame in attrs breaks concat
        rows.append(out.assign(**info))
    d = OUT / "htec"
    d.mkdir(parents=True, exist_ok=True)
    pd.concat(rows).to_csv(d / f"{ct}.csv", index=False)


def sample_scores(X, sf, genes, keys, sets):
    """Per sample (keys): log2 ANTXR2 per 10k and each set's score (mean over
    genes of log2(per10k + 0.1), centred per gene over samples)."""
    keys = np.asarray(keys)
    u = pd.unique(keys)
    G = sp.csr_matrix((np.ones(len(keys)), (pd.Categorical(keys, categories=u).codes, np.arange(len(keys)))),
                      shape=(len(u), len(keys)))
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    gpos = gpos[~gpos.index.duplicated()]
    depth = np.asarray(G @ sf).ravel()
    out = pd.DataFrame({"n": np.asarray(G.sum(1)).ravel()}, index=u)
    out["antxr2"] = np.log2(np.asarray(G @ X[:, gpos[TARGET]].toarray()).ravel() / depth + 0.01)
    for k in ("IFN_I_biased", "IFN_G_biased", "IFN_shared"):
        gs = [g for g in sets[k] if g in gpos.index]
        L = np.log2(np.asarray((G @ X[:, gpos[gs].values]).todense()) / depth[:, None] + 0.1)
        out[k] = (L - L.mean(0)).mean(1)
    return out


def fit(d, covars, label):
    """OLS of ANTXR2 on the two biased set scores (+ covariates), standardized
    predictors, donor-clustered SEs when donors repeat."""
    import statsmodels.formula.api as smf
    d = d.copy()
    for k in ("IFN_I_biased", "IFN_G_biased"):
        d[k] = (d[k] - d[k].mean()) / d[k].std()
    f = "antxr2 ~ IFN_I_biased + IFN_G_biased" + "".join(f" + C({c})" for c in covars)
    cluster = "donor" not in covars and d.donor.duplicated().any()   # HTEC: 2 donors as fixed effects
    kw = dict(cov_type="cluster", cov_kwds={"groups": pd.factorize(d.donor)[0]}) if cluster else {}
    m = smf.ols(f, d).fit(**kw)
    rows = []
    for k in ("IFN_I_biased", "IFN_G_biased"):
        rows.append(dict(model=label, term=k, coef=m.params[k], se=m.bse[k], p=m.pvalues[k], n=len(d),
                         r_IG=d.IFN_I_biased.corr(d.IFN_G_biased),
                         marginal_r=d.antxr2.corr(d[k])))
    return rows


def cmd_donors(cells, min_cells=20):
    """Between samples: does ANTXR2's mean follow type I- or IFN-gamma-biased
    activity? Means discriminate the cytokine (HTEC: alpha / beta raise ANTXR2,
    gamma does not); within-cell coupling does not (cmd_htec). HTEC first as the
    positive control, then each IBDverse population over donor x state samples."""
    sets = load_sets()
    (OUT / "donors").mkdir(exist_ok=True)
    rows = []
    X, o, genes = H.load()
    for ct in ["ciliated", "basal", "club"]:
        m = (o.ct == ct).values
        s = sample_scores(X[m], o.total_counts.values[m] / 1e4, genes,
                          (o.donor + "|" + o.stim + "|" + o.time.astype(str)).values[m], sets)
        s[["donor", "stim", "time"]] = s.index.to_series().str.split("|", expand=True).values
        rows += fit(s[s.n >= min_cells], ["donor"], f"htec_{ct}")
    del X
    pops = {}
    for cell in cells:
        Xc, oc, gc, _ = M.load_cell(cell)
        pop, state = cell.split("@")
        s = sample_scores(Xc, oc.memento_size_factor.values, gc, oc.donor.astype(str).values, sets)
        pops.setdefault(pop, []).append(s.assign(state=state, donor=s.index))
    for pop, ss in pops.items():
        d = pd.concat(ss)
        d = d[d.n >= min_cells]
        rows += fit(d, ["state"], pop)
        for st, g in d.groupby("state"):
            if len(g) >= 10:
                rows += fit(g, [], f"{pop}@{st}")
        d.to_csv(OUT / "donors" / f"{pop}.csv")
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "donors" / "fits.csv", index=False)
    pd.set_option("display.width", 200)
    print(t.round(3).to_string(index=False))


AMBIENT_IG = ["IGKC", "IGHA1", "IGHA2", "JCHAIN", "IGLC2", "IGLC3", "IGHM", "IGHG1", "MZB1"]   # plasma-cell soup
IMMUNE = ["PTPRC", "CD52", "CORO1A", "LCP1", "CD3E", "CD3D", "CD2", "CD79A", "MS4A1", "LYZ", "TYROBP",
          "FCER1G", "CD14", "C1QA"]   # pan-immune, no IFN-inducible genes
MIN_STRATUM = 10


def cmd_contam(cells, boot, backend):
    """Is the tissue IFN coupling carried by immune-derived reads (user,
    2026-11-21: "we don't have to be super hardcore about it")? Per cell, the
    IFN set re-estimated (i) on all cells, (ii) without likely immune doublets
    (>= 2 UMI of IMMUNE genes), (iii) without each donor's top quarter by
    ambient Ig share (AMBIENT_IG UMI / size factor); plus ANTXR2's coupling to
    the two contamination programs themselves. Strata left under MIN_STRATUM
    cells after a filter are dropped."""
    sets = load_sets()
    extra = {**sets, "ambient_Ig": AMBIENT_IG, "immune_doublet": IMMUNE}
    d = OUT / "contam" / M.LANDSCAPE_OUT.name
    d.mkdir(parents=True, exist_ok=True)
    for cell in cells:
        X, o, genes, anchors = M.load_cell(cell)
        o = o.reset_index(drop=True)
        gpos = pd.Series(np.arange(len(genes)), index=genes)
        imm = np.asarray(X[:, gpos[[g for g in IMMUNE if g in gpos.index]].values].sum(1)).ravel()
        amb = np.asarray(X[:, gpos[[g for g in AMBIENT_IG if g in gpos.index]].values].sum(1)).ravel() \
            / o.memento_size_factor.values
        hi_amb = amb > pd.Series(amb).groupby(o.donor.values).transform(lambda v: v.quantile(.75)).values
        filters = {"all": np.ones(len(o), bool), "no_immune_doublet": imm < 2, "low_ambient": ~hi_amb}
        rows = []
        for f, m in filters.items():
            grp = (o.donor.astype(str) + "^" + o.fine.astype(str)).values
            n = pd.Series(grp[m]).value_counts()
            m = m & np.isin(grp, n.index[n >= MIN_STRATUM])
            out = H.couple(X[m], o[m], genes, [cell], boot, IFN8, backend, f"{cell}:{f}", anchors=anchors, extra=extra)
            out.attrs = {}
            rows.append(out.assign(filter=f, n_cells=int(m.sum()), kept=m.mean(),
                                   immune_umi_mean=imm[m].mean(), ambient_per10k_mean=amb[m].mean()))
        t = pd.concat(rows, ignore_index=True)
        t.to_csv(d / f"{cell}.csv", index=False)
        show = t[t.unit.isin(IFN8 + SETS + ["ambient_Ig", "immune_doublet"])]
        print(f"=== {cell}: kept " + ", ".join(f"{f} {k:.2f}" for f, k in t.groupby("filter", sort=False).kept.first().items()))
        print(show.pivot(index="unit", columns="filter", values="excess")[list(filters)].round(3).to_string())


BULK = {   # Hagai 2018 bulk (E-MTAB-5919): ANTXR2 and MX1 / IFIT1 / RSAD2 Ensembl ids per species
    "human": ("ENSG00000163297", ["ENSG00000157601", "ENSG00000185745", "ENSG00000134321"]),
    "rhesus": ("ENSMMUG00000003333", ["ENSMMUG00000015329", "ENSMMUG00000017533", "ENSMMUG00000015969"]),
    "mouse": ("ENSMUSG00000029338", ["ENSMUSG00000000386", "ENSMUSG00000034459", "ENSMUSG00000020641"]),
    "rat": ("ENSRNOG00000000081", ["ENSRNOG00000001959", "ENSRNOG00000019050", "ENSRNOG00000007539"]),
}


def cmd_bulk():
    """Untreated vs IFN-beta 8 h dermal fibroblasts, paired by donor (bulk,
    E-MTAB-5919 `<species>_IFNB_unst.txt`). The count columns are Sanger lane ids
    with no donor / treatment key in the SDRF: treatment is called from the ISG
    score (top half), and each treated sample is paired with the untreated one it
    correlates with best on the 1,000 most variable non-responsive genes (always
    the adjacent lane id). Median-of-ratios size factors (ISGs take read share)."""
    rows = []
    for sp, (ant, isg) in BULK.items():
        m = pd.read_csv(FIB_IFNB / "bulk" / f"{sp}_IFNB_unst.txt", sep=r"\s+", index_col=0)
        lg = np.log(m[(m > 0).all(1)])
        sf = np.exp(lg.sub(lg.mean(1), axis=0).median())
        L = np.log2(m / sf / m.sum().mean() * 1e6 + 1)
        score = L.reindex(isg).mean()
        treated = (score > score.median()).values
        expd = L[(L > np.log2(11)).all(1)]
        d = expd.loc[:, treated].mean(1) - expd.loc[:, ~treated].mean(1)
        st = expd[d.abs() < 0.2]
        st = st.sub(st.mean(1), axis=0)
        C = st.loc[st.var(1).nlargest(1000).index].corr()
        for t in m.columns[treated]:
            u = C.loc[t, m.columns[~treated]].idxmax()
            rows.append(dict(species=sp, treated=t, untreated=u, pair_r=C.loc[t, u], isg_lfc=score[t] - score[u],
                             antxr2_unst=2 ** L.loc[ant, u] - 1, antxr2_ifnb=2 ** L.loc[ant, t] - 1,
                             antxr2_lfc=L.loc[ant, t] - L.loc[ant, u]))
    t = pd.DataFrame(rows)
    (OUT / "fibroblast").mkdir(exist_ok=True)
    t.to_csv(OUT / "fibroblast" / "bulk_ifnb.csv", index=False)
    print(t.round(3).to_string(index=False))
    print(t.groupby("species").antxr2_lfc.agg(["size", "mean", "min", "max"]).round(3).to_string())


def ens_to_symbol():
    g = pd.read_csv(FIB_DSRNA / "pIC" / "genes.tsv", sep="\t", header=None, names=["ens", "sym"])
    return g.drop_duplicates("ens").set_index("ens").sym


def load_dsrna():
    mats, obs = [], []
    for sub, cond in (("unst", "control"), ("pIC", "dsRNA_6h")):
        X = sio.mmread(FIB_DSRNA / sub / "matrix.mtx").tocsr().T.tocsr()
        g = pd.read_csv(FIB_DSRNA / sub / "genes.tsv", sep="\t", header=None)
        genes = g[1].values.astype(str)
        tot = np.asarray(X.sum(1)).ravel()
        mt = np.asarray(X[:, np.char.startswith(genes, "MT-")].sum(1)).ravel() / np.maximum(tot, 1)
        keep = (tot >= FIB_MIN_UMI) & (mt <= FIB_MAX_MT)
        print(f"dsrna {cond}: {keep.sum()} / {len(tot)} cells kept, median UMI {np.median(tot[keep]):.0f}")
        mats.append(X[keep])
        obs.append(pd.DataFrame({"donor": "mix", "condition": cond, "fine": "fib", "total": tot[keep]}))
    X = sp.vstack(mats).tocsr()
    o = pd.concat(obs, ignore_index=True)
    dup = pd.Index(genes).duplicated()
    return X[:, ~dup], o, genes[~dup]


def load_ifnb():
    m = pd.read_csv(FIB_IFNB / "expr_count_matrix.csv", index_col=0)
    m = m[m.index.str.startswith("ENSG")]
    sym = ens_to_symbol().reindex(m.index)
    m = m[sym.notna().values]
    m.index = sym.dropna().values
    m = m.groupby(level=0).sum()
    s = pd.read_csv(FIB_IFNB / "E-MTAB-7051.sdrf.txt", sep="\t").drop_duplicates("Source Name")
    s = s.set_index("Source Name")
    cells = m.columns.str.replace("#", "_")
    o = pd.DataFrame({"donor": s["Characteristics[individual]"].reindex(cells).values,
                      "condition": "IFNb_" + s["Factor Value[time]"].reindex(cells).astype(str).values + "h"})
    X = sp.csr_matrix(np.rint(m.values.T.clip(min=0)))
    tot = np.asarray(X.sum(1)).ravel()
    keep = o.donor.notna().values & (tot >= 1e5)
    print(f"ifnb: {keep.sum()} / {len(o)} cells kept (>= 1e5 reads), median reads {np.median(tot[keep]):.0f}")
    o = o[keep].reset_index(drop=True).assign(fine="fib", total=tot[keep])
    return X[keep], o, m.index.values.astype(str)


def cmd_fibroblast(boot, backend, all_programs=False, only=None):
    """all_programs: every library program plus the sets (is IFN the top of
    ANTXR2's profile, or does it couple to everything that changes?)."""
    sets = load_sets()
    d = OUT / "fibroblast"
    d.mkdir(parents=True, exist_ok=True)
    means = []
    for name, loader in (("dsrna", load_dsrna), ("ifnb", load_ifnb)):
        if only and name not in only:
            continue
        X, o, genes = loader()
        j = list(genes).index(TARGET)
        o["antxr2"] = X[:, j].toarray().ravel()
        t = o.groupby(["condition", "donor"]).agg(cells=("antxr2", "size"), counts=("antxr2", "sum"),
                                                  depth=("total", "sum"), det=("antxr2", lambda v: (v > 0).mean()))
        t["per10k"] = t.counts / t.depth * 1e4
        t = t.reset_index().assign(dataset=name)
        print(t.round(4).to_string(index=False))
        means.append(t)
        obs = o[["donor", "condition", "fine"]].assign(memento_size_factor=o.total.values / 1e4)
        conds = list(dict.fromkeys(o.condition))
        extra, info = residual_sets(X, obs, genes, sets)   # b pooled over conditions, within strata
        out = H.couple(X, obs, genes, conds, boot, None if all_programs else IFN8, backend,
                       f"fibroblast_{name}", extra=extra).assign(**info)
        if all_programs:
            out.to_csv(d / f"{name}_all.csv", index=False)
            for cond, g in out.groupby("stim", sort=False):
                g = g.assign(rank=g.excess.rank(ascending=False).astype(int))
                print(f"--- {name} / {cond}: {len(g)} programs, q < 0.05 {(g.excess_q < .05).sum()} "
                      f"(positive {((g.excess_q < .05) & (g.excess > 0)).sum()}), median excess {g.excess.median():.3f}")
                print(g.nsmallest(25, "rank")[["rank", "unit", "excess", "excess_q", "pct"]].round(3).to_string(index=False))
                print(g[g.unit.isin(IFN8 + SETS + RESID)][["rank", "unit", "excess", "excess_q"]].round(3).to_string(index=False))
            continue
        out.to_csv(d / f"{name}.csv", index=False)
        summarize(out, name)
    if not all_programs and not only:
        pd.concat(means).to_csv(d / "mean.csv", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["sets", "tissue", "fibroblast", "htec", "donors", "bulk", "contam"])
    ap.add_argument("--cells", nargs="*")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"])
    ap.add_argument("--all-programs", action="store_true", help="fibroblast: every library program too")
    ap.add_argument("--only", nargs="*", choices=["dsrna", "ifnb"])
    a = ap.parse_args()
    if a.cmd == "sets":
        cmd_sets()
    elif a.cmd == "tissue":
        cmd_tissue(a.cells, a.boot, a.backend)
    elif a.cmd == "htec":
        cmd_htec(a.boot, a.backend)
    elif a.cmd == "donors":
        cmd_donors(a.cells)
    elif a.cmd == "bulk":
        cmd_bulk()
    elif a.cmd == "contam":
        cmd_contam(a.cells, a.boot, a.backend)
    else:
        cmd_fibroblast(a.boot, a.backend, a.all_programs, a.only)


if __name__ == "__main__":
    main()
