"""Difference between two cell populations (program_config.CONTRASTS; e.g.
fibroblast - enterocyte) in ANTXR2 x gene-program association, with one fixed
weight vector per program shared by both populations.

Estimator (planning_summaries/program_coexpression.md, "Estimator"): memento's
hypergeometric model corrects only the diagonal of the covariance, so for a
fixed weight vector a on the X/sf scale and a target gene g outside the program,
within one group

    Var(s)   = sample var(sum_h a_h X_h/sf) - (1-q) sum_h a_h^2 mean(X_h/sf^2)
    Cov(s,g) = sample cov(sum_h a_h X_h/sf, X_g/sf)

and the cell-type estimate pools these within-group moments over groups (i.e.
centres on the group, Decision 4). Groups are donor x subtype, and donor x
maturation bin for a population flagged `maturation` (enterocytes;
program_config.MATURATION_MARKERS), so neither population's association is
driven by its differentiation structure. A population may also be one disease
state of a cell type, so A - B can be a state contrast (inflamed - uninflamed).
Weights (Decision 3) are on genes standardised within each cell type:
a_c = v / sigma_c, sigma_c the pooled corrected SD, v unit-norm and shared:

    equal      published signs / norm                     (Decision 5)
    pc1        PC1 of the equal-weight average of the two cell types' pooled
               corrected correlation matrices, oriented so the score
               correlates positively with the signed equal-weight score
    published  PROGENy / CytoSig weights / norm, where they exist

The statistic is delta = rho_a - rho_b (populations in CONTRASTS order), and
its excess over 50 anchor genes matched to ANTXR2 in BOTH populations:
excess = delta_ANTXR2 - mean(delta_anchor). Everything reduces to per-group sufficient statistics, so
the donor bootstrap (one donor resample shared by both cell types) is a
reweighting of groups and carries ANTXR2 and the anchors together. When the
two populations share no donor (a state contrast between people), donors are
resampled within each population instead.

`--validate` first checks the per-group single-gene correlation from these
formulas against memento's own compute_2d_moments on the same size factors.

--contrast picks the slice (source superset npz + age range, applied here
before size factors are computed) and the two populations.

Outputs in OUT_DIR/<contrast>_<grouping>[_mat<N>][_gm][_covNN]/:
    differential.csv   one row per program x scheme: delta, anchor delta, excess
    associations.csv   one row per program x scheme x population (the two rhos,
                       score reliability = corrected / raw score variance)
    group_rho.csv      per donor group ANTXR2 rho, per program x scheme
    programs.csv       one row per program: coverage, status, coherence, PC1 VE
    loadings.csv       one row per program x gene: weights per scheme
    anchors.csv        anchor genes and their matching features
    maturation.csv     enterocyte maturation bins: score, ANTXR2 mean, cells
    run.json           parameters, group table, sanity checks

Usage:
    python program_scores.py --contrast <name> [--grouping fine|broad]
                             [--maturation-bins N] [--min-coverage F]
                             [--estimator pooled|groupmean] [--validate]
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from program_config import (
    CONTRASTS, DISCOVERY_PROGRAMS, FILTER_MEAN_THRESH, GENESET_DIR, HYPOTHESIS_PROGRAMS,
    MATURATION_MARKERS,
    MIN_COVERAGE, MIN_GENES, MIN_GENES_DISCOVERY,
    MIN_GROUP_CELLS, N_ANCHORS, N_BOOT, N_MATURATION_BINS, OUT_DIR, Q, SEED, SHRINKAGE,
    SLICES, TARGET, TRIM_PERCENT, cells_npz,
)

OBS_KEYS = ["donor", "cell_type", "fine", "study", "age", "region", "assay", "disease"]


# ------------------------------------------------------------------ inputs

def load_cells(npz_path):
    z = np.load(npz_path, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.DataFrame({k: z[k].astype(str) if k in z else np.full(z["shape"][0], "")
                        for k in OBS_KEYS})
    obs["total_counts"] = z["total_counts"]
    genes = pd.Index(z["genes"].astype(str))
    if genes.has_duplicates:                     # keep the first occurrence
        keep = ~genes.duplicated()
        X, genes = X[:, np.where(keep)[0]], genes[keep]
    return X, obs, genes


def load_library():
    """programs_long.csv plus the fixed discovery / hypothesis modules (equal
    positive weights)."""
    lib = pd.read_csv(GENESET_DIR / "programs_long.csv")
    return pd.concat([lib] + [pd.DataFrame(dict(collection=col, program=k, gene=v,
                                                sign=1.0, weight=np.nan,
                                                program_id=f"{col}:{k}"))
                              for col, fixed in [("discovery", DISCOVERY_PROGRAMS),
                                                 ("bracq", HYPOTHESIS_PROGRAMS)]
                              for k, v in fixed.items()], ignore_index=True)


def size_factors(X, genes):
    """memento's size factors on fibroblasts + enterocytes jointly, so both
    cell types share one normalisation."""
    import anndata as ad
    import memento
    A = ad.AnnData(X=X.copy(), obs=pd.DataFrame(index=[f"c{i}" for i in range(X.shape[0])]),
                   var=pd.DataFrame(index=genes))
    A.obs["q"] = Q
    memento.setup_memento(A, q_column="q", filter_mean_thresh=FILTER_MEAN_THRESH,
                          shrinkage=SHRINKAGE, trim_percent=TRIM_PERCENT)
    return A.obs["memento_size_factor"].values.astype(np.float64)


def maturation_bins(Y, genes, obs, nbins, pops):
    """Enterocyte maturation score (mean within-donor z of log1p normalised
    mature markers minus that of the immature markers) and within-donor
    quantile bins m1 (least mature) .. mN, for the cells of `pops`. Other
    cells get ""."""
    obs["mat_bin"], obs["mat_score"] = "", np.nan
    if nbins <= 0 or not pops:
        return
    ent = np.where(obs.cell_type.isin(pops).values)[0]
    sides = {}
    for side, gl in MATURATION_MARKERS.items():
        missing = [g for g in gl if g not in genes]
        assert not missing, f"maturation markers absent from data: {missing}"
        L = pd.DataFrame(np.log1p(Y[ent][:, [genes.get_loc(g) for g in gl]].toarray()),
                         columns=gl)
        d = obs.donor.values[ent]
        z = L.groupby(d).transform(lambda x: (x - x.mean()) / (x.std(ddof=0) or 1))
        sides[side] = z.mean(1).values
    score = sides["mature"] - sides["immature"]
    obs.loc[ent, "mat_score"] = score
    rank = pd.Series(score).groupby(obs.donor.values[ent]).rank(pct=True).values
    obs.loc[ent, "mat_bin"] = [f"m{b}" for b in np.minimum((rank * nbins).astype(int), nbins - 1) + 1]


def make_groups(obs, grouping, floors):
    """floors: population -> minimum cells per group."""
    key = obs.cell_type + "|" + obs.donor
    if grouping == "fine":
        key = key + "|" + obs.fine
    key = key + np.where(obs.mat_bin != "", "|" + obs.mat_bin, "")
    n = key.value_counts()
    pop_of = key.str.split("|").str[0].groupby(key).first()
    keep = {k for k, c in n.items() if c >= floors[pop_of[k]]}
    obs["group"] = np.where(key.isin(keep), key, "")
    gt = (obs[obs.group != ""].groupby("group")
          .agg(cell_type=("cell_type", "first"), donor=("donor", "first"),
               fine=("fine", lambda s: s.iloc[0] if grouping == "fine" else "*"),
               mat_bin=("mat_bin", "first"), n=("donor", "size")).reset_index())
    gt["rows"] = [np.where(obs.group.values == g)[0] for g in gt.group]
    return gt


# ------------------------------------------------------------------ moments

def gene_stats(X, rows):
    Xc = X[rows]
    raw_mean = np.asarray(Xc.sum(0)).ravel() / len(rows)
    det = np.asarray((Xc > 0).sum(0)).ravel() / len(rows)
    return raw_mean, det


def noise_sums(X, sf, rows):
    """sum over the group's cells of X/sf^2, per gene (memento's diagonal term)."""
    w = 1.0 / sf[rows] ** 2
    return np.asarray(X[rows].T @ w).ravel()


def pooled_var(Y, X, sf, groups, gidx):
    """Pooled within-group corrected variance (the diagonal only), gidx genes."""
    v = np.zeros(len(gidx))
    N = 0
    for rows in groups.rows:
        Yg = Y[rows][:, gidx].astype(np.float64)
        s = np.asarray(Yg.sum(0)).ravel()
        v += np.asarray(Yg.multiply(Yg).sum(0)).ravel() - s ** 2 / len(rows) \
            - (1 - Q) * noise_sums(X, sf, rows)[gidx]
        N += len(rows)
    return v / N


def pooled_cov(Y, X, sf, groups, gidx, block=20000):
    """Pooled within-group corrected covariance over genes gidx: one dense
    product over within-group-centred cells, in blocks of whole groups."""
    p = len(gidx)
    C = np.zeros((p, p))
    noise = np.zeros(p)
    N = 0
    chunk, size = [], 0

    def flush():
        B = np.vstack(chunk).astype(np.float32)
        C[:] += (B.T @ B).astype(np.float64)

    for rows in groups.rows:
        Yg = Y[rows][:, gidx].toarray().astype(np.float64)
        chunk.append(Yg - Yg.mean(0))
        size += len(rows)
        noise += (1 - Q) * noise_sums(X, sf, rows)[gidx]
        N += len(rows)
        if size >= block:
            flush()
            chunk, size = [], 0
    if chunk:
        flush()
    C[np.diag_indices(p)] -= noise
    return C / N


def cov_to_corr(C):
    d = np.sqrt(np.clip(np.diag(C), 1e-300, None))
    R = C / np.outer(d, d)
    np.fill_diagonal(R, 1.0)
    return R


# ------------------------------------------------------------------ programs

def build_units(lib, genes, passing, pops, min_coverage=MIN_COVERAGE):
    """Comparable programs (Decision 6): genes clearing the floor in BOTH cell
    types, >= MIN_GENES of them and >= min_coverage of the measured members.
    Programs passing in only one cell type cannot carry a difference; they are
    listed with their status and not scored."""
    present = set(genes)
    shared = passing[pops[0]] & passing[pops[1]]
    prog_rows, units = [], []
    for pid, df in lib.groupby("program_id", sort=True):
        meas = df[df.gene.isin(present)]
        n_meas = len(meas)
        in_c = {c: meas[meas.gene.isin(passing[c])] for c in pops}
        both = meas[meas.gene.isin(shared)]
        fixed = pid.startswith(("discovery:", "bracq:"))
        min_genes = MIN_GENES_DISCOVERY if fixed else MIN_GENES
        ok = lambda d: len(d) >= min_genes and n_meas > 0 and len(d) / n_meas >= min_coverage
        row = dict(program_id=pid, collection=df.collection.iloc[0], n_members=len(df),
                   n_measured=n_meas, n_pass_a=len(in_c[pops[0]]),
                   n_pass_b=len(in_c[pops[1]]), n_pass_both=len(both))
        if ok(both):
            units.append(dict(unit=pid, genes=both))
            row["status"] = "comparable"
        else:
            solo = [c for c in pops if ok(in_c[c])]
            row["status"] = ("one_cell_type:" + "+".join(solo)) if solo else "insufficient"
        prog_rows.append(row)
    return pd.DataFrame(prog_rows), units


def learn_weights(unit, R, pos, pops):
    """Unit-norm weight vectors per scheme + diagnostics. R[c] is the cell
    type's corrected correlation matrix over its program genes, pos[c] maps
    gene -> row in it."""
    g = unit["genes"]
    e = g.sign.values.astype(float)
    sub = {c: R[c][np.ix_(pos[c][g.gene.values], pos[c][g.gene.values])] for c in pops}
    Rbar = np.mean([sub[c] for c in pops], axis=0)
    lam, V = np.linalg.eigh(Rbar)
    v1 = V[:, -1]
    if v1 @ Rbar @ e < 0:
        v1 = -v1
    p = len(g)
    schemes = {"equal": e / np.linalg.norm(e), "pc1": v1}
    if g.weight.notna().all():
        w = g.weight.values.astype(float)
        schemes["published"] = w / np.linalg.norm(w)
    off = ~np.eye(p, dtype=bool)
    ee = np.outer(e, e)
    diag = dict(n_genes=p, pc1_ve=lam[-1] / p,
                cos_pc1_equal=float(v1 @ schemes["equal"]))
    for c in pops:
        diag[f"pc1_ve_{c}"] = float(v1 @ sub[c] @ v1 / p)
        diag[f"coherence_{c}"] = float((ee * sub[c])[off].mean())
    return schemes, diag


def match_anchors(feats, target_feats, candidates, k):
    """k nearest genes to ANTXR2 in standardised feature space."""
    F = feats.loc[candidates]
    sd = F.std(ddof=0).replace(0, 1)
    d = (((F - target_feats) / sd) ** 2).sum(1)
    return d.sort_values().index[:k].tolist()


# ------------------------------------------------------------------ pooling

def pool_rho(m, S, Tt, ST):
    """Pooled correlation from per-group sufficient statistics.

    m: (B, G) group multiplicities. S = (n, s_sum, s_sq, s_noise) with
    s_* shape (G, K); Tt = (t_sum, t_sq, t_noise) shape (G, T);
    ST shape (G, K, T). Returns (B, K, T)."""
    n, s_sum, s_sq, s_noise = S
    t_sum, t_sq, t_noise = Tt
    N = m @ n
    with np.errstate(invalid="ignore", divide="ignore"):
        vs = (m @ (s_sq - s_sum ** 2 / n[:, None] - s_noise)) / N[:, None]
        vt = (m @ (t_sq - t_sum ** 2 / n[:, None] - t_noise)) / N[:, None]
        cross = ST - s_sum[:, :, None] * t_sum[:, None, :] / n[:, None, None]
        cov = np.tensordot(m, cross, axes=(1, 0)) / N[:, None, None]
        vs = np.where(vs > 0, vs, np.nan)
        vt = np.where(vt > 0, vt, np.nan)
        return cov / np.sqrt(vs[:, :, None] * vt[:, None, :])


def group_mean_rho(w, rg):
    """Cell-weighted mean of per-group correlations, each clipped to [-1, 1]
    (memento's own get_2d_moments / one-sample convention). w: (B, G) weights
    (multiplicity x cells), rg: (G, K, T) per-group rho. Returns (B, K, T).

    Unlike pool_rho's ratio of sums, a group whose noise-corrected variance is
    over-subtracted cannot drag every other group's estimate past |1|."""
    ok = np.isfinite(rg)
    r0 = np.where(ok, np.clip(rg, -1, 1), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.einsum("bg,gkt->bkt", w, r0) / np.einsum("bg,gkt->bkt", w, ok.astype(float))


def group_stats(X, Y, sf, groups, W, T):
    """Per-group sufficient statistics for scores Y @ W against targets T."""
    G, K = len(groups), W.shape[1]
    n = groups.n.values.astype(float)
    s_sum, s_sq, s_noise = (np.zeros((G, K)) for _ in range(3))
    t_sum, t_sq, t_noise = (np.zeros((G, len(T))) for _ in range(3))
    ST = np.zeros((G, K, len(T)))
    W2 = W.multiply(W).tocsc()
    for i, rows in enumerate(groups.rows):
        Yg = Y[rows].astype(np.float64)
        S = np.asarray((Yg @ W).todense())
        YT = Yg[:, T].toarray()
        z = noise_sums(X, sf, rows)
        s_sum[i], s_sq[i] = S.sum(0), (S ** 2).sum(0)
        s_noise[i] = (1 - Q) * (W2.T @ z)
        t_sum[i], t_sq[i] = YT.sum(0), (YT ** 2).sum(0)
        t_noise[i] = (1 - Q) * z[T]
        ST[i] = S.T @ YT
    return (n, s_sum, s_sq, s_noise), (t_sum, t_sq, t_noise), ST


# ------------------------------------------------------------------ validation

def validate(X, genes, sf, obs, groups, passing, pops, rng, n_genes=30, n_groups=4):
    """Single-gene programs vs memento.compute_2d_moments, per group, same sf."""
    import anndata as ad
    import memento
    both = sorted((passing[pops[0]] & passing[pops[1]]) - {TARGET})
    pick = list(rng.choice(both, size=min(n_genes, len(both)), replace=False))
    gsel = groups.sort_values("n").iloc[np.linspace(0, len(groups) - 1, n_groups).astype(int)]
    rows = np.concatenate(gsel.rows.values)
    cols = [genes.get_loc(g) for g in [TARGET] + pick]
    A = ad.AnnData(X=X[rows][:, cols].copy(), obs=obs.iloc[rows][["group"]].copy(),
                   var=pd.DataFrame(index=[TARGET] + pick))
    A.obs["q"] = Q
    A.obs.index = [f"c{i}" for i in range(len(rows))]
    memento.setup_memento(A, q_column="q", filter_mean_thresh=0.0,
                          shrinkage=SHRINKAGE, trim_percent=TRIM_PERCENT)
    A.obs["memento_size_factor"] = sf[rows]       # use the full-data size factors
    memento.create_groups(A, label_columns=["group"])
    memento.compute_1d_moments(A, min_perc_group=0.0, filter_genes=False)
    memento.compute_2d_moments(A, [(TARGET, g) for g in pick])
    mom, _ = memento.get_2d_moments(A)
    out = []
    Y = sp.diags(1 / sf) @ X
    for _, gr in gsel.iterrows():
        col = f"sg^{gr.group}"
        if col not in mom.columns:
            continue
        r = gr.rows
        yt = Y[r][:, genes.get_loc(TARGET)].toarray().ravel()
        z = noise_sums(X, sf, r)
        vt = yt.var() - (1 - Q) * z[genes.get_loc(TARGET)] / len(r)
        for g in pick:
            j = genes.get_loc(g)
            yh = Y[r][:, j].toarray().ravel()
            vh = yh.var() - (1 - Q) * z[j] / len(r)
            mine = np.cov(yh, yt, bias=True)[0, 1] / np.sqrt(vh * vt) if vh > 0 and vt > 0 else np.nan
            mem = mom.loc[mom.gene_2 == g, col].values[0]
            out.append(dict(group=gr.group, gene=g, ours=mine, memento=mem))
    v = pd.DataFrame(out)
    ok = v.dropna()
    unclipped = ok[ok.ours.abs() < 1]
    err = float((unclipped.ours - unclipped.memento).abs().max()) if len(unclipped) else np.nan
    print(f"validate: {len(ok)} group x gene pairs, max |ours - memento| = {err:.2e} "
          f"(on {len(unclipped)} unclipped)")
    return v, err


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contrast", required=True, choices=sorted(CONTRASTS))
    ap.add_argument("--grouping", choices=["fine", "broad"], default="fine")
    ap.add_argument("--maturation-bins", type=int, default=N_MATURATION_BINS)
    ap.add_argument("--min-coverage", type=float, default=MIN_COVERAGE,
                    help="sensitivity only; Decision 6 fixes the primary at MIN_COVERAGE")
    ap.add_argument("--estimator", choices=["pooled", "groupmean"], default="pooled",
                    help="pooled: ratio of summed within-group moments; groupmean: "
                         "cell-weighted mean of clipped per-group rho (memento convention)")
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    con = CONTRASTS[a.contrast]
    pops = list(con["pops"])
    mat_pops = [p for p, d in con["pops"].items() if d.get("maturation")]
    nbins = a.maturation_bins if mat_pops else 0
    suffix = (f"_mat{nbins}" if mat_pops else "") + ("_gm" if a.estimator == "groupmean" else "") + \
        ("" if a.min_coverage == MIN_COVERAGE else f"_cov{round(100 * a.min_coverage)}")
    out = OUT_DIR / f"{a.contrast}_{a.grouping}{suffix}"
    out.mkdir(parents=True, exist_ok=True)
    sl = SLICES[con["slice"]]
    npz = cells_npz(sl["source"])
    meta = dict(npz=str(npz), contrast=a.contrast, pops=pops, pop_defs=con["pops"],
                slice=con["slice"], grouping=a.grouping,
                estimator=a.estimator, maturation_pops=mat_pops, maturation_bins=nbins, maturation_markers=MATURATION_MARKERS,
                q=Q, shrinkage=SHRINKAGE, trim_percent=TRIM_PERCENT,
                filter_mean_thresh=FILTER_MEAN_THRESH, min_group_cells=MIN_GROUP_CELLS,
                min_genes=MIN_GENES, min_coverage=a.min_coverage, n_anchors=N_ANCHORS,
                n_boot=N_BOOT, seed=SEED)

    X, obs, genes = load_cells(npz)
    pop = np.full(len(obs), "", dtype=object)
    for name, d in con["pops"].items():
        m = obs.cell_type.values == d["cell_type"]
        if "fine" in d:
            m &= obs.fine.isin(d["fine"]).values
        if "disease" in d:
            m &= obs.disease.isin(d["disease"]).values
        pop[m] = name
    inslice = pop != ""
    for key, col in [("ages", "age"), ("disease", "disease"), ("regions", "region")]:
        if key in sl:
            inslice &= obs[col].isin(sl[key]).values
    keep = np.where(inslice)[0]
    X, obs = X[keep], obs.iloc[keep].reset_index(drop=True)
    obs["cell_type"] = pop[keep]
    sf = size_factors(X, genes)
    Y = (sp.diags(1 / sf) @ X).tocsr().astype(np.float32)
    maturation_bins(Y, genes, obs, nbins, mat_pops)
    floors = {p: d.get("min_group_cells", MIN_GROUP_CELLS) for p, d in con["pops"].items()}
    meta["min_group_cells"] = floors
    groups = make_groups(obs, a.grouping, floors)
    print(f"{X.shape[0]} cells x {X.shape[1]} genes; {len(groups)} groups "
          f"({time.time() - t0:.0f}s)", flush=True)
    meta["groups"] = groups.drop(columns="rows").to_dict("records")
    tix = genes.get_loc(TARGET)

    if nbins > 0:
        e = obs[(obs.group != "") & obs.cell_type.isin(mat_pops)]
        yt = np.asarray(X[e.index.values][:, tix].todense()).ravel()
        mt = (e.assign(antxr2=yt, antxr2_det=yt > 0).groupby(["cell_type", "mat_bin"])
              .agg(cells=("donor", "size"), donors=("donor", "nunique"),
                   score=("mat_score", "mean"), antxr2_raw_mean=("antxr2", "mean"),
                   antxr2_det=("antxr2_det", "mean"))
              .reset_index())
        mt.to_csv(out / "maturation.csv", index=False)
        print(mt.round(3).to_string(index=False), flush=True)

    # per-cell-type gene stats, floor, pooled corrected moments. Variances for
    # every floor gene (anchors need them); the correlation matrix only over
    # floor genes that belong to some program.
    lib = load_library()
    in_lib = genes.isin(set(lib.gene))
    feats, passing, R, pos, sigma = {}, {}, {}, {}, {}
    for c in pops:
        gc = groups[groups.cell_type == c]
        rows = np.concatenate(gc.rows.values)
        raw_mean, det = gene_stats(X, rows)
        feats[c] = pd.DataFrame({"log10_mean": np.log10(raw_mean + 1e-6), "det": det,
                                 "raw_mean": raw_mean}, index=genes)
        floor = np.where(raw_mean >= FILTER_MEAN_THRESH)[0]
        var = pooled_var(Y, X, sf, gc, floor)
        idx = floor[var > 0]
        sigma[c] = pd.Series(np.sqrt(var[var > 0]), index=genes[idx])
        passing[c] = set(genes[idx])
        ridx = idx[in_lib[idx]]
        R[c] = cov_to_corr(pooled_cov(Y, X, sf, gc, ridx)).astype(np.float32)
        pos[c] = pd.Series(np.arange(len(ridx)), index=genes[ridx])
        # Decision 8 sanity check: random-pair median correlation
        iu = rng.integers(0, len(ridx), size=(20000, 2))
        iu = iu[iu[:, 0] != iu[:, 1]]
        rp = R[c][iu[:, 0], iu[:, 1]]
        meta[f"random_pair_{c}"] = dict(median=float(np.median(rp)), mean=float(rp.mean()),
                                        n=int(len(rp)))
        meta[f"target_{c}"] = dict(raw_mean=float(raw_mean[tix]), det=float(det[tix]),
                                   passes_floor=TARGET in passing[c],
                                   n_groups=int(len(gc)), n_donors=int(gc.donor.nunique()),
                                   n_cells=int(len(rows)), n_genes_floor=int(len(idx)),
                                   n_genes_corr=int(len(ridx)))
        print(f"{c}: {len(rows)} cells, {len(gc)} groups, {gc.donor.nunique()} donors, "
              f"{len(idx)} genes at floor; {TARGET} raw mean {raw_mean[tix]:.3f} "
              f"det {det[tix]:.3f}; random-pair median r = {np.median(rp):+.4f} "
              f"({time.time() - t0:.0f}s)", flush=True)

    if a.validate:
        v, err = validate(X, genes, sf, obs, groups, passing, pops, rng)
        v.to_csv(out / "validate.csv", index=False)
        meta["validate_max_abs_err"] = err

    # programs and weights (ANTXR2 already removed from the library)
    for c in pops:
        passing[c].discard(TARGET)
    progs, units = build_units(lib, genes, passing, pops, a.min_coverage)
    diag_rows, load_rows = [], []
    for u in units:
        u["schemes"], d = learn_weights(u, R, pos, pops)
        diag_rows.append(dict(program_id=u["unit"], **d))
        L = u["genes"][["gene", "sign", "weight"]].copy()
        for s, v in u["schemes"].items():
            L[f"v_{s}"] = v
        for c in pops:
            L[f"sigma_{c}"] = sigma[c].loc[L.gene].values
        L.insert(0, "unit", u["unit"])
        load_rows.append(L)
    pd.concat(load_rows).to_csv(out / "loadings.csv", index=False)
    progs = progs.merge(pd.DataFrame(diag_rows), on="program_id", how="left")
    progs.to_csv(out / "programs.csv", index=False)
    print(progs.status.str.split(":").str[0].value_counts().to_string(), flush=True)

    # anchors matched to ANTXR2 in both populations at once
    jf = pd.concat({c: feats[c][["log10_mean", "det"]] for c in pops}, axis=1)
    jf.columns = [f"{a_}_{b}" for a_, b in jf.columns]
    anchors = match_anchors(jf, jf.loc[TARGET],
                            sorted(passing[pops[0]] & passing[pops[1]]), N_ANCHORS)
    targets = [TARGET] + anchors
    pd.DataFrame([dict(set="target" if g == TARGET else "anchor", gene=g,
                       **{f"{c}_{k}": feats[c].loc[g, k] for c in pops
                          for k in ["raw_mean", "det"]}) for g in targets]
                 ).to_csv(out / "anchors.csv", index=False)
    T = np.array([genes.get_loc(g) for g in targets])

    # one score column per program x scheme; an anchor inside a program is not
    # a null for it
    keys = [(u, s) for u in units for s in u["schemes"]]
    mask = np.ones((len(keys), len(targets)), bool)
    for k, (u, s) in enumerate(keys):
        for j, g in enumerate(targets):
            if g in set(u["genes"].gene):
                mask[k, j] = False

    donors = sorted(groups.donor.unique())
    dpos = {d: i for i, d in enumerate(donors)}
    pop_donors = [set(groups.donor[groups.cell_type == c]) for c in pops]
    paired = bool(pop_donors[0] & pop_donors[1])
    meta["bootstrap"] = "paired" if paired else "stratified by population"
    strata = [np.arange(len(donors))] if paired else \
        [np.array([dpos[d] for d in sorted(D)]) for D in pop_donors]
    boot_counts = np.zeros((N_BOOT, len(donors)), int)
    for b in range(N_BOOT):
        for st in strata:
            np.add.at(boot_counts[b], st[rng.integers(0, len(st), len(st))], 1)
    point, boot, assoc, reliability, group_rows = {}, {}, [], {}, []
    for c in pops:
        gc = groups[groups.cell_type == c].reset_index(drop=True)
        rr, cc, vv = [], [], []
        for k, (u, s) in enumerate(keys):
            gi = np.array([genes.get_loc(g) for g in u["genes"].gene])
            rr.append(gi)
            cc.append(np.full(len(gi), k))
            vv.append(u["schemes"][s] / sigma[c].loc[u["genes"].gene].values)
        W = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))),
                          shape=(len(genes), len(keys)))
        S, Tt, ST = group_stats(X, Y, sf, gc, W, T)
        m = boot_counts[:, [dpos[d] for d in gc.donor]].astype(float)
        if a.estimator == "pooled":
            point[c] = np.where(mask, pool_rho(np.ones((1, len(gc))), S, Tt, ST)[0], np.nan)
            boot[c] = np.where(mask[None], pool_rho(m, S, Tt, ST), np.nan)
        else:
            w = np.vstack([np.ones(len(gc)), m]) * S[0][None]
            r = group_mean_rho(w, pool_rho(np.eye(len(gc)), S, Tt, ST))
            point[c] = np.where(mask, r[0], np.nan)
            boot[c] = np.where(mask[None], r[1:], np.nan)
        lo, hi = np.nanpercentile(boot[c][:, :, 0], [2.5, 97.5], axis=0)
        # reliability: share of the pooled within-group score variance left
        # after memento's noise correction. Near 0, the corrected variance in
        # the rho denominator is mostly subtraction and rho is unstable.
        n, s_sum, s_sq, s_noise = S
        raw = (s_sq - s_sum ** 2 / n[:, None]).sum(0)
        rel = (raw - s_noise.sum(0)) / raw
        t_raw = (Tt[1][:, 0] - Tt[0][:, 0] ** 2 / n).sum()
        reliability[c] = rel
        # per-group ANTXR2 rho (identity multiplicities), to see whether a
        # pooled value is carried by one donor group
        pg = pool_rho(np.eye(len(gc)), S, tuple(t[:, :1] for t in Tt), ST[:, :, :1])[:, :, 0]
        group_rows.append(pd.DataFrame(pg.T, columns=gc.group).assign(
            unit=[u["unit"] for u, _ in keys], scheme=[s_ for _, s_ in keys], pop=c))
        meta[f"target_{c}"]["reliability"] = float((t_raw - Tt[2][:, 0].sum()) / t_raw)
        for k, (u, s) in enumerate(keys):
            assoc.append(dict(unit=u["unit"], scheme=s, cell_type=c, n_genes=len(u["genes"]),
                              reliability=rel[k], rho=point[c][k, 0], ci_lo=lo[k], ci_hi=hi[k],
                              anchor_mean=np.nanmean(point[c][k, 1:]),
                              anchor_sd=np.nanstd(point[c][k, 1:])))
        print(f"{c}: {len(keys)} program x scheme scores ({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(assoc).to_csv(out / "associations.csv", index=False)
    pd.concat([g.melt(id_vars=["unit", "scheme", "pop"], var_name="group", value_name="rho")
               for g in group_rows]).merge(groups[["group", "n"]], on="group") \
        .to_csv(out / "group_rho.csv", index=False)

    # the difference, and ANTXR2's excess over the anchors' difference; every
    # bootstrap draw uses the same donor resample in both populations
    pa, pb = pops
    dP = point[pa] - point[pb]                                 # (K, T)
    dB = boot[pa] - boot[pb]                                   # (B, K, T)
    with np.errstate(invalid="ignore"):
        exP = dP[:, 0] - np.nanmean(dP[:, 1:], 1)
        exB = dB[:, :, 0] - np.nanmean(dB[:, :, 1:], 2)
        asd = np.nanstd(dP[:, 1:], 1)
        n_a = np.isfinite(dP[:, 1:]).sum(1)
        pct = 100 * np.nansum(dP[:, 1:] < dP[:, :1], 1) / n_a
    dlo, dhi = np.nanpercentile(dB[:, :, 0], [2.5, 97.5], axis=0)
    elo, ehi = np.nanpercentile(exB, [2.5, 97.5], axis=0)
    fpos = np.nanmean(exB > 0, 0)
    diff = pd.DataFrame(dict(
        unit=[u["unit"] for u, _ in keys], scheme=[s for _, s in keys],
        n_genes=[len(u["genes"]) for u, _ in keys],
        rho_a=point[pa][:, 0], rho_b=point[pb][:, 0],
        reliability_a=reliability[pa], reliability_b=reliability[pb],
        delta=dP[:, 0], delta_ci_lo=dlo, delta_ci_hi=dhi,
        anchor_delta_mean=np.nanmean(dP[:, 1:], 1), anchor_delta_sd=asd,
        anchor_delta_pct=pct, n_anchors=n_a,
        excess=exP, excess_ci_lo=elo, excess_ci_hi=ehi, excess_z=exP / asd,
        p_boot_excess=2 * np.minimum(fpos, 1 - fpos)))
    diff.sort_values(["scheme", "unit"]).to_csv(out / "differential.csv", index=False)
    meta["runtime_s"] = round(time.time() - t0, 1)
    json.dump(meta, open(out / "run.json", "w"), indent=1, default=str)
    print(f"done -> {out} ({meta['runtime_s']}s)")


if __name__ == "__main__":
    main()
