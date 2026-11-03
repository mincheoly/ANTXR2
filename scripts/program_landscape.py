"""ANTXR2 co-expression landscape over IBDverse ileum populations x states
(planning_summaries/program_coexpression.md, "Landscape design").

The matrix: one cell per population x state (program_config.LANDSCAPE_POPS
+ pooled fibroblasts; LANDSCAPE_STATES). In every filled cell, ANTXR2's
corrected association with each program -- one weight vector per program for
the whole matrix -- and with every floor gene, each calibrated against 50
genes matched to ANTXR2 in that population. Estimator, grouping and memento
settings are program_scores.py's; this module reuses its functions.

Three stages, each reading the previous one's output in LANDSCAPE_OUT:

    prep     per population: size factors (all states jointly), donor x
             subtype groups, and per state the floor, pooled corrected
             variances and the correlation matrix over library genes
             -> prep/<pop>@<state>.npz, cells.csv
    weights  core genes (floor in >= CORE_PASS_FRAC of filled cells) and one
             PC1 per program from the core correlation matrix averaged over
             filled cells (each entry over the cells where both genes pass);
             the shared donor bootstrap (donors resampled within state)
             -> weights.csv, programs.csv, boot_counts.npz
    score    per population: anchors within ANCHOR_SPAN_WINDOW of ANTXR2's
             genomic span, matched on raw mean + detection jointly across its
             filled states (the per-cell calibration), plus the lineage's
             anchors (LANDSCAPE_LINEAGES: matched jointly across every filled
             cell of the lineage, scored in all of them -- same-gene nulls for
             lineage-wide consistency); per cell the program scores (core genes passing there,
             >= CELL_CORE_COVERAGE of the core, weights renormalised) and the
             genome-wide gene correlations
             -> landscape.csv (point estimates, CIs), boot/<pop>.npz (excess
                draws, for the CI of any difference between cells, and every
                target's point rho for anchor-based calibration),
                genes/<pop>@<state>.csv.gz

Cells with < LANDSCAPE_MIN_DONORS donors or ANTXR2 below the raw-mean floor
are listed in cells.csv as not measurable and never scored.

Usage:
    python program_landscape.py {prep,weights,score,all} [--pops SLUG ...] [--validate]
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from program_config import (
    ANCHOR_SPAN_WINDOW, CELL_CORE_COVERAGE, CORE_PASS_FRAC, FILTER_MEAN_THRESH, LANDSCAPE_CELLS,
    LANDSCAPE_FIB, LANDSCAPE_FIB_FLOOR, LANDSCAPE_LINEAGES, LANDSCAPE_MIN_DONORS, LANDSCAPE_OUT,
    LANDSCAPE_POPS, LANDSCAPE_STATES,
    MIN_GENES, MIN_GENES_DISCOVERY, MIN_GROUP_CELLS, N_ANCHORS, N_BOOT, SEED, TARGET, gene_spans,
    pop_slug,
)
from program_scores import (
    cov_to_corr, gene_stats, group_stats, load_cells, load_library, make_groups, match_anchors,
    pool_rho, pooled_cov, pooled_var, size_factors, validate,
)

POPS = LANDSCAPE_POPS + [LANDSCAPE_FIB]
STATE_OF = {v: k for k, v in LANDSCAPE_STATES.items()}
PREP, BOOT, GENES = (LANDSCAPE_OUT / d for d in ("prep", "boot", "genes"))
for _d in (PREP, BOOT, GENES):
    _d.mkdir(parents=True, exist_ok=True)


def cell_id(pop, state):
    return f"{pop_slug(pop)}@{state}"


def load_population(pop):
    """Counts, size factors, normalised counts and groups for one population;
    obs.cell_type becomes the matrix cell (population@state)."""
    X, obs, genes = load_cells(LANDSCAPE_CELLS / f"cells_{pop_slug(pop)}.npz")
    state = obs.disease.map(STATE_OF)
    keep = np.where(state.notna().values)[0]
    X, obs = X[keep], obs.iloc[keep].reset_index(drop=True)
    obs["cell_type"] = [cell_id(pop, s) for s in state.iloc[keep]]
    obs["mat_bin"] = ""
    sf = size_factors(X, genes)
    Y = (sp.diags(1 / sf) @ X).tocsr().astype(np.float32)
    floor = LANDSCAPE_FIB_FLOOR if pop == LANDSCAPE_FIB else MIN_GROUP_CELLS
    groups = make_groups(obs, "fine", {c: floor for c in obs.cell_type.unique()})
    return X, Y, sf, obs, genes, groups


# ------------------------------------------------------------------ prep

def prep(pops, lib):
    rows = []
    for pop in pops:
        t0 = time.time()
        X, Y, sf, obs, genes, groups = load_population(pop)
        in_lib = genes.isin(set(lib.gene))
        tix = genes.get_loc(TARGET)
        for state in LANDSCAPE_STATES:
            c = cell_id(pop, state)
            gc = groups[groups.cell_type == c]
            r = dict(population=pop, state=state, cell=c, n_cells=int((obs.cell_type == c).sum()),
                     n_groups=len(gc), n_donors=int(gc.donor.nunique()), cells_in_groups=0)
            if r["n_donors"] < LANDSCAPE_MIN_DONORS:
                rows.append(dict(r, filled=False, reason=f"< {LANDSCAPE_MIN_DONORS} donors"))
                continue
            gr = np.concatenate(gc.rows.values)
            raw_mean, det = gene_stats(X, gr)
            r.update(cells_in_groups=len(gr), antxr2_raw_mean=raw_mean[tix], antxr2_det=det[tix])
            if raw_mean[tix] < FILTER_MEAN_THRESH:
                rows.append(dict(r, filled=False, reason="ANTXR2 below floor"))
                continue
            floor = np.where(raw_mean >= FILTER_MEAN_THRESH)[0]
            var = pooled_var(Y, X, sf, gc, floor)
            idx, var = floor[var > 0], var[var > 0]
            ridx = idx[in_lib[idx]]
            R = cov_to_corr(pooled_cov(Y, X, sf, gc, ridx)).astype(np.float32)
            np.savez(PREP / f"{c}.npz", raw_mean=raw_mean.astype(np.float32),
                     det=det.astype(np.float32), idx=idx, var=var, ridx=ridx, R=R,
                     genes=np.asarray(genes, dtype=str))
            rows.append(dict(r, filled=True, reason="", n_genes_floor=len(idx),
                             n_genes_corr=len(ridx)))
        print(f"prep {pop}: {len(obs)} cells, {len(groups)} groups "
              f"({time.time() - t0:.0f}s)", flush=True)
    cells = pd.DataFrame(rows)
    path = LANDSCAPE_OUT / "cells.csv"
    if path.exists():                      # --pops reruns replace their own rows
        old = pd.read_csv(path)
        cells = pd.concat([old[~old.population.isin(pops)], cells], ignore_index=True)
    cells.to_csv(path, index=False)
    return cells


def load_prep(c):
    z = np.load(PREP / f"{c}.npz", allow_pickle=False)
    genes = pd.Index(z["genes"].astype(str))
    return z, genes


# ------------------------------------------------------------------ weights

def weights(lib):
    cells = pd.read_csv(LANDSCAPE_OUT / "cells.csv")
    filled = cells[cells.filled].cell.tolist()
    passing, genes = {}, None
    for c in filled:
        z, genes = load_prep(c)
        passing[c] = set(genes[z["idx"]]) - {TARGET}
    frac = pd.Series({g: np.mean([g in passing[c] for c in filled])
                      for g in set(lib.gene) & set(genes)})
    prog_rows, core = [], {}
    for pid, df in lib.groupby("program_id", sort=True):
        df = df[df.gene.isin(genes) & (df.gene != TARGET)].drop_duplicates("gene")
        k = df[df.gene.map(frac).fillna(0) >= CORE_PASS_FRAC]
        min_genes = MIN_GENES_DISCOVERY if pid.startswith(("discovery:", "bracq:")) else MIN_GENES
        prog_rows.append(dict(program_id=pid, collection=df.collection.iloc[0] if len(df) else "",
                              n_measured=len(df), n_core=len(k),
                              status="core" if len(k) >= min_genes else "too few core genes"))
        if len(k) >= min_genes:
            core[pid] = k
    U = sorted(set().union(*[set(k.gene) for k in core.values()]))
    upos = pd.Series(np.arange(len(U)), index=U)
    Ssum = np.zeros((len(U), len(U)))
    Scnt = np.zeros((len(U), len(U)), np.float32)
    for c in filled:
        z, genes = load_prep(c)
        rg = genes[z["ridx"]]
        keep = np.where(rg.isin(U))[0]
        u = upos[rg[keep]].values
        Ssum[np.ix_(u, u)] += z["R"][np.ix_(keep, keep)]
        Scnt[np.ix_(u, u)] += 1
    wrows = []
    diag = {}
    for pid, k in core.items():
        u = upos[k.gene].values
        with np.errstate(invalid="ignore"):
            Rbar = np.nan_to_num(Ssum[np.ix_(u, u)] / Scnt[np.ix_(u, u)])
        e = k.sign.values.astype(float)
        lam, V = np.linalg.eigh(Rbar)
        v1 = V[:, -1] if V[:, -1] @ Rbar @ e >= 0 else -V[:, -1]
        eq = e / np.linalg.norm(e)
        diag[pid] = dict(pc1_ve=lam[-1] / len(u), cos_pc1_equal=float(v1 @ eq),
                         min_pair_cells=int(Scnt[np.ix_(u, u)].min()))
        wrows.append(pd.DataFrame(dict(program_id=pid, gene=k.gene.values, sign=e,
                                       v_equal=eq, v_pc1=v1, pass_frac=k.gene.map(frac).values)))
    programs = pd.DataFrame(prog_rows)
    programs = programs.join(pd.DataFrame(diag).T, on="program_id")
    programs.to_csv(LANDSCAPE_OUT / "programs.csv", index=False)
    pd.concat(wrows).to_csv(LANDSCAPE_OUT / "weights.csv", index=False)

    # one donor bootstrap for the whole matrix: donors are resampled within
    # their state (each IBDverse ileum donor is in exactly one state)
    dstate = {}
    for pop in POPS:
        z = np.load(LANDSCAPE_CELLS / f"cells_{pop_slug(pop)}.npz", allow_pickle=False)
        for d, s in set(zip(z["donor"].astype(str), z["disease"].astype(str))):
            if s in STATE_OF:
                assert dstate.setdefault(d, s) == s, f"donor {d} in two states"
    donors = np.array(sorted(dstate))
    rng = np.random.default_rng(SEED)
    counts = np.zeros((N_BOOT, len(donors)), np.int32)
    strata = [np.where(np.array([dstate[d] for d in donors]) == s)[0] for s in LANDSCAPE_STATES.values()]
    for b in range(N_BOOT):
        for st in strata:
            np.add.at(counts[b], st[rng.integers(0, len(st), len(st))], 1)
    np.savez(LANDSCAPE_OUT / "boot_counts.npz", donors=donors, counts=counts)
    print(f"weights: {len(core)} programs with a core of >= MIN_GENES genes "
          f"({len(U)} genes); bootstrap over {len(donors)} donors", flush=True)


# ------------------------------------------------------------------ score

def pooled_cross(Y, groups, T, gidx):
    """Pooled within-group covariance of targets T with genes gidx. Off the
    diagonal memento's estimator is the plain sample covariance."""
    C = np.zeros((len(T), len(gidx)))
    N = 0
    for rows in groups.rows:
        Yg = Y[rows]
        A = Yg[:, T].toarray().astype(np.float64)
        A -= A.mean(0)
        C += (Yg[:, gidx].T @ A).T
        N += len(rows)
    return C / N


def span_matched(feats, passing, span, k):
    """k genes of ANTXR2's genomic length (ANCHOR_SPAN_WINDOW), nearest to it
    on raw mean + detection across every cell in feats / passing."""
    jf = pd.concat(feats, axis=1)
    jf.columns = [f"{a}_{b}" for a, b in jf.columns]
    cand = [g for g in sorted(set.intersection(*passing) - {TARGET})
            if g in span.index and abs(span[g] - span[TARGET]) <= ANCHOR_SPAN_WINDOW]
    return match_anchors(jf, jf.loc[TARGET], cand, k), len(cand)


def lineage_anchors(cells, span):
    out = {}
    for lin, members in LANDSCAPE_LINEAGES.items():
        fc = cells[cells.population.isin(members) & cells.filled].cell
        if fc.empty:
            continue
        feats, passing = {}, []
        for c in fc:
            z, genes = load_prep(c)
            feats[c] = pd.DataFrame({"log10_mean": np.log10(z["raw_mean"] + 1e-6),
                                     "det": z["det"]}, index=genes)
            passing.append(set(genes[z["idx"]]))
        out[lin], n = span_matched(feats, passing, span, N_ANCHORS)
        print(f"lineage anchors {lin}: {len(fc)} cells, {n} candidates", flush=True)
    return out


def score(pops, do_validate):
    cells = pd.read_csv(LANDSCAPE_OUT / "cells.csv")
    W_all = pd.read_csv(LANDSCAPE_OUT / "weights.csv")
    bc = np.load(LANDSCAPE_OUT / "boot_counts.npz", allow_pickle=True)
    dpos = {d: i for i, d in enumerate(bc["donors"].astype(str))}
    counts = bc["counts"].astype(float)
    units = list(W_all.program_id.unique())
    span = np.log10(gene_spans())
    keys = [(u, s) for u in units for s in ("pc1", "equal")]
    rows, meta = [], {}
    rng = np.random.default_rng(SEED)
    lin_anchors = lineage_anchors(cells, span)
    lineage_of = {p: lin for lin, ps in LANDSCAPE_LINEAGES.items() for p in ps}
    for pop in pops:
        t0 = time.time()
        filled = cells[(cells.population == pop) & cells.filled]
        if filled.empty:
            print(f"score {pop}: no filled cell", flush=True)
            continue
        X, Y, sf, obs, genes, groups = load_population(pop)
        prep_z = {c: load_prep(c)[0] for c in filled.cell}
        passing = {c: set(genes[z["idx"]]) for c, z in prep_z.items()}
        sigma = {c: pd.Series(np.sqrt(z["var"]), index=genes[z["idx"]]) for c, z in prep_z.items()}

        # anchors: genes of ANTXR2's genomic length (the length axis, planning
        # note "Correction"), matched to it on raw mean + detection jointly
        # across the population's filled states
        anchors, n_cand = span_matched(
            {c: pd.DataFrame({"log10_mean": np.log10(z["raw_mean"] + 1e-6), "det": z["det"]},
                             index=genes) for c, z in prep_z.items()},
            list(passing.values()), span, N_ANCHORS)
        lin = lineage_of[pop]
        targets = [TARGET] + anchors + lin_anchors[lin]    # point rho for all;
        nb = 1 + len(anchors)                              # excess + bootstrap on the first nb
        meta[pop] = dict(anchors=anchors, n_candidates=n_cand, lineage=lin,
                         lineage_anchors=lin_anchors[lin],
                         anchor_log10_span_median=float(span[anchors].median()),
                         antxr2_log10_span=float(span[TARGET]))
        T = np.array([genes.get_loc(g) for g in targets])

        boot_out = {}
        for c in filled.cell:
            state = c.split("@")[1]
            gc = groups[groups.cell_type == c].reset_index(drop=True)
            if do_validate:
                v, err = validate(X, genes, sf, obs, gc, {c: passing[c]}, [c, c], rng)
                meta[pop][f"validate_{state}"] = err
            # program scores restricted to the core genes passing here
            rr, cc, vv, info = [], [], [], []
            for k, (u, s) in enumerate(keys):
                w = W_all[W_all.program_id == u]
                ok = w.gene.isin(passing[c]).values
                n_ok, cov = int(ok.sum()), ok.mean()
                min_genes = MIN_GENES_DISCOVERY if u.startswith(("discovery:", "bracq:")) else MIN_GENES
                scored = cov >= CELL_CORE_COVERAGE and n_ok >= min_genes
                info.append(dict(n_genes=n_ok, coverage=cov, scored=scored,
                                 genes=set(w.gene[ok])))
                if not scored:
                    continue
                v = w[f"v_{s}"].values[ok]
                v = v / np.linalg.norm(v)
                gi = np.array([genes.get_loc(g) for g in w.gene[ok]])
                rr.append(gi)
                cc.append(np.full(len(gi), k))
                vv.append(v / sigma[c].loc[w.gene[ok]].values)
            cells.loc[cells.cell == c, "n_programs"] = sum(i["scored"] for i in info) // 2
            if not vv:                   # no program reaches core coverage here
                print(f"  {c}: no program scored; gene level only", flush=True)
            else:
                Wm = sp.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))),
                                   shape=(len(genes), len(keys)))
                mask = np.array([[inf["scored"] and g not in inf["genes"] for g in targets]
                                 for inf in info])
                S, Tt, ST = group_stats(X, Y, sf, gc, Wm, T)
                m = counts[:, [dpos[d] for d in gc.donor]]
                point_all = np.where(mask, pool_rho(np.ones((1, len(gc))), S, Tt, ST)[0], np.nan)
                point, mask = point_all[:, :nb], mask[:, :nb]
                Tt, ST = tuple(t[:, :nb] for t in Tt), ST[:, :, :nb]
                boot = np.where(mask[None], pool_rho(m, S, Tt, ST), np.nan)
                with np.errstate(invalid="ignore"):
                    exP = point[:, 0] - np.nanmean(point[:, 1:], 1)
                    exB = boot[:, :, 0] - np.nanmean(boot[:, :, 1:], 2)
                boot_out[f"{state}_excess"] = exB.astype(np.float32)
                boot_out[f"{state}_rho"] = boot[:, :, 0].astype(np.float32)
                boot_out[f"{state}_point"] = point_all.astype(np.float32)
                n, s_sum, s_sq, s_noise = S
                raw = (s_sq - s_sum ** 2 / n[:, None]).sum(0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    rel = (raw - s_noise.sum(0)) / raw
                t_raw = (Tt[1][:, 0] - Tt[0][:, 0] ** 2 / n).sum()
                cells.loc[cells.cell == c, "antxr2_reliability"] = (t_raw - Tt[2][:, 0].sum()) / t_raw
                lo, hi = np.nanpercentile(boot[:, :, 0], [2.5, 97.5], axis=0)
                elo, ehi = np.nanpercentile(exB, [2.5, 97.5], axis=0)
                fpos = np.nanmean(exB > 0, 0)
                asd = np.nanstd(point[:, 1:], 1)
                for k, (u, s) in enumerate(keys):
                    if not info[k]["scored"]:
                        continue
                    rows.append(dict(population=pop, state=state, cell=c, unit=u, scheme=s,
                                     n_genes=info[k]["n_genes"], coverage=info[k]["coverage"],
                                     reliability=rel[k], rho=point[k, 0], rho_ci_lo=lo[k],
                                     rho_ci_hi=hi[k], anchor_mean=np.nanmean(point[k, 1:]),
                                     anchor_sd=asd[k], excess=exP[k], excess_ci_lo=elo[k],
                                     excess_ci_hi=ehi[k], excess_z=exP[k] / asd[k],
                                     p_boot_excess=2 * min(fpos[k], 1 - fpos[k])))

            # genome-wide: ANTXR2 and anchors against every floor gene
            z = prep_z[c]
            gidx = z["idx"]
            vt = pooled_var(Y, X, sf, gc, T)
            C = pooled_cross(Y, gc, T, gidx)
            with np.errstate(invalid="ignore", divide="ignore"):
                Rg = C / np.sqrt(np.where(vt > 0, vt, np.nan)[:, None] * z["var"][None, :])
            Rg[np.asarray(gidx)[None, :] == T[:, None]] = np.nan     # a gene with itself
            am, asd_g = np.nanmean(Rg[1:nb], 0), np.nanstd(Rg[1:nb], 0)
            pd.DataFrame(dict(gene=genes[gidx], raw_mean=z["raw_mean"][gidx], det=z["det"][gidx],
                              rho=Rg[0], anchor_mean=am, anchor_sd=asd_g,
                              z=(Rg[0] - am) / asd_g, is_anchor=np.isin(gidx, T[1:nb]))) \
                .dropna(subset=["rho"]).to_csv(GENES / f"{c}.csv.gz", index=False)
        np.savez_compressed(BOOT / f"{pop_slug(pop)}.npz", units=np.array([u for u, _ in keys]),
                            targets=np.array(targets), n_pop_anchors=len(anchors),
                            schemes=np.array([s for _, s in keys]), **boot_out)
        print(f"score {pop}: {len(filled)} cells ({time.time() - t0:.0f}s)", flush=True)

    L = pd.DataFrame(rows)
    path = LANDSCAPE_OUT / "landscape.csv"
    if path.exists():
        old = pd.read_csv(path)
        L = pd.concat([old[~old.population.isin(pops)], L], ignore_index=True)
    L.to_csv(path, index=False)
    cells.to_csv(LANDSCAPE_OUT / "cells.csv", index=False)       # + antxr2_reliability
    mpath = LANDSCAPE_OUT / "score_meta.json"
    allmeta = json.load(open(mpath)) if mpath.exists() else {}
    allmeta.update(meta)
    json.dump(allmeta, open(mpath, "w"), indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["prep", "weights", "score", "all"])
    ap.add_argument("--pops", nargs="*", help="population slugs (default: all)")
    ap.add_argument("--validate", action="store_true")
    a = ap.parse_args()
    pops = [p for p in POPS if not a.pops or pop_slug(p) in a.pops]
    lib = load_library()
    if a.stage in ("prep", "all"):
        prep(pops, lib)
    if a.stage in ("weights", "all"):
        weights(lib)
    if a.stage in ("score", "all"):
        score(pops, a.validate)


if __name__ == "__main__":
    main()
