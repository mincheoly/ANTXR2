"""Donor-bootstrap test for the reference-scaled coupling: is ANTXR2's coupling
to a program different between two population x state cells of one atlas?

Statistic, per program (PC1 scheme), from the same pooled within-group
estimator as program_landscape_ref.score:
    beta      ANTXR2's reference-scaled coupling in a cell
    excess    beta minus the mean beta of the cell's population anchors
    d_beta    beta_A - beta_B          d_excess   excess_A - excess_B

Bootstrap: donors are the sampling unit. Each draw resamples the union of the
two cells' donors with replacement (a donor present in both cells enters both
with the same multiplicity, so shared donors stay paired) and re-sums the
stored per-group (donor x subtype) sufficient statistics; no cell is resampled
within a donor. se = SD over draws; p is the normal approximation
2 * Phi(-|estimate| / se); ci = 2.5th / 97.5th percentile of the draws; q =
Benjamini-Hochberg over programs. `active_*` = activity_sd >= 0.2 in the cell.

compare      A vs B                  -> LANDSCAPE_OUT/test/{A}__vs__{B}.csv
calibrate    random donor halves of ONE cell tested against each other (the
             true difference is zero): share of programs with p < 0.05 per
             split                   -> LANDSCAPE_OUT/test/calibration.csv

Usage:
    ANTXR2_LANDSCAPE_RUN=ref [ANTXR2_LANDSCAPE=gca] python program_landscape_ref_test.py \
        compare --a "Enterocyte_top_villus@normal" --b "Enterocyte_stem_OLFM4lo_LGR5hi@normal"
    ... calibrate --cells "Enterocyte_top_villus@normal" [--splits 20]
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from program_config import LANDSCAPE_OUT, SEED, TARGET  # noqa: E402
import program_landscape as PL                     # noqa: E402
import program_landscape_ref as REF                # noqa: E402
from program_scores import group_stats             # noqa: E402

OUT = LANDSCAPE_OUT / "test"
OUT.mkdir(parents=True, exist_ok=True)
N_BOOT = 2000
SCHEME = "pc1"
MIN_ACTIVITY_SD = 0.2
CHUNK = 100


class Cell:
    """Per-group sufficient statistics of one population x state cell."""

    def __init__(self, cell, progs, units, sigma):
        cells = REF.filled_cells().set_index("cell")
        assert cell in cells.index, f"{cell} is not a filled cell of {LANDSCAPE_OUT.name}"
        pop, scale = cells.population[cell], cells.scale[cell]
        anchors = json.load(open(LANDSCAPE_OUT / "score_meta.json"))[pop]["anchors"]
        targets = [TARGET] + anchors
        X, Y, sf, obs, genes, groups = PL.load_population(pop)
        gc = groups[groups.cell_type == cell].reset_index(drop=True)
        gpos = pd.Series(np.arange(len(genes)), index=genes)
        rr = np.concatenate([gpos[progs[u]["genes"]].values for u in units])
        cc = np.concatenate([np.full(len(progs[u]["genes"]), k) for k, u in enumerate(units)])
        vv = np.concatenate([progs[u][SCHEME] / (progs[u]["s_ref"] * scale) for u in units])
        Wm = sp.csr_matrix((vv, (rr, cc)), shape=(len(genes), len(units)))
        (n, s_sum, _, _), (t_sum, t_sq, t_noise), ST = group_stats(X, Y, sf, gc, Wm, gpos[targets].values)
        self.cell, self.sigma, self.donor, self.n = cell, sigma, gc.donor.values, n
        # per-group contributions: sums over groups give N * cov and N * var
        self.C = (ST - s_sum[:, :, None] * t_sum[:, None, :] / n[:, None, None]).reshape(len(n), -1)
        self.V = t_sq - t_sum ** 2 / n[:, None] - t_noise
        self.member = np.array([[g in set(progs[u]["genes_all"]) for g in targets] for u in units])
        full = pd.read_csv(LANDSCAPE_OUT / "coupling.csv")
        full = full[(full.scheme == SCHEME) & (full.cell == cell)].set_index("unit").reindex(units)
        self.active = full.activity_sd.values >= MIN_ACTIVITY_SD
        self.full_beta = full.beta.values

    def estimate(self, M):
        """M: draws x groups multiplicities -> (beta, excess), each draws x programs."""
        K, T = self.member.shape
        beta, excess = (np.empty((len(M), K)) for _ in range(2))
        for s in range(0, len(M), CHUNK):
            m = M[s:s + CHUNK]
            N = (m @ self.n)[:, None]
            vt = m @ self.V / N
            cov = (m @ self.C).reshape(len(m), K, T) / N[:, :, None]
            with np.errstate(invalid="ignore", divide="ignore"):
                b = cov / (np.sqrt(np.where(vt > 0, vt, np.nan))[:, None, :] * self.sigma[None, :, None])
                b = np.where(self.member[None], np.nan, b)
                beta[s:s + CHUNK] = b[:, :, 0]
                excess[s:s + CHUNK] = b[:, :, 0] - np.nanmean(b[:, :, 1:], 2)
        return beta, excess


def setup():
    progs = REF.load_programs()
    units = list(progs)
    return progs, units, np.array([progs[u][f"sigma_{SCHEME}"] for u in units])


def multiplicities(donors_all, donor_of_group, draws):
    """draws x donors counts -> draws x groups; row 0 of `draws` is the observed sample."""
    pos = pd.Series(np.arange(len(donors_all)), index=donors_all)
    return draws[:, pos[donor_of_group].values].astype(np.float64)


def draw_donors(n_donors, rng, n_boot=N_BOOT):
    d = rng.multinomial(n_donors, np.full(n_donors, 1 / n_donors), size=n_boot)
    return np.vstack([np.ones(n_donors, dtype=int), d])


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    o = np.argsort(p[ok])
    r = p[ok][o] * ok.sum() / np.arange(1, ok.sum() + 1)
    q[np.where(ok)[0][o]] = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    return q


def infer(x, prefix):
    """x: (1 + draws) x programs, row 0 observed."""
    est, boot = x[0], x[1:]
    se = np.nanstd(boot, 0)
    lo, hi = np.nanquantile(boot, [0.025, 0.975], axis=0)
    p = 2 * norm.sf(np.abs(est) / se)
    return {prefix: est, f"{prefix}_se": se, f"{prefix}_lo": lo, f"{prefix}_hi": hi, f"{prefix}_p": p,
            f"{prefix}_q": bh(p)}


def compare(A, B, units, rng, draws=None, Ma=None, Mb=None):
    if Ma is None:
        donors = np.array(sorted(set(A.donor) | set(B.donor)))
        draws = draw_donors(len(donors), rng)
        Ma, Mb = multiplicities(donors, A.donor, draws), multiplicities(donors, B.donor, draws)
    # a draw that leaves a cell without donors is undefined for it
    ok = (Ma.sum(1) > 0) & (Mb.sum(1) > 0)
    (ba, ea), (bb, eb) = A.estimate(Ma[ok]), B.estimate(Mb[ok])
    out = dict(unit=units, active_a=A.active, active_b=B.active)
    for k, v in (("beta_a", ba), ("beta_b", bb), ("excess_a", ea), ("excess_b", eb),
                 ("d_beta", ba - bb), ("d_excess", ea - eb)):
        out.update(infer(v, k))
    return pd.DataFrame(out), int(ok.sum() - 1)


def cmd_compare(a):
    progs, units, sigma = setup()
    rng = np.random.default_rng(SEED)
    A, B = Cell(a.a, progs, units, sigma), Cell(a.b, progs, units, sigma)
    for c in (A, B):   # the observed row must reproduce the scored coupling
        assert np.allclose(c.estimate(np.ones((1, len(c.n))))[0][0], c.full_beta, equal_nan=True), c.cell
    df, n_boot = compare(A, B, units, rng)
    shared = len(set(A.donor) & set(B.donor))
    path = OUT / f"{a.a}__vs__{a.b}.csv"
    df.to_csv(path, index=False)
    print(f"{a.a} ({len(set(A.donor))} donors) vs {a.b} ({len(set(B.donor))} donors), {shared} shared; "
          f"{n_boot} draws -> {path}")
    for k in ("d_beta", "d_excess"):
        both = df[df.active_a | df.active_b]
        print(f"{k}: {len(both)} programs active in either; q < 0.05: {(both[f'{k}_q'] < 0.05).sum()}, "
              f"p < 0.05: {(both[f'{k}_p'] < 0.05).sum()}")
    cols = ["unit", "beta_a", "beta_b", "d_beta", "d_beta_se", "d_beta_p", "d_beta_q", "d_excess",
            "d_excess_se", "d_excess_q"]
    s = df.sort_values("d_beta_p")
    print(s[cols].head(25).round(4).to_string(index=False))


def cmd_calibrate(a):
    progs, units, sigma = setup()
    rng = np.random.default_rng(SEED)
    rows = []
    for cell in a.cells:
        C = Cell(cell, progs, units, sigma)
        donors = np.array(sorted(set(C.donor)))
        for i in range(a.splits):
            half = set(rng.permutation(donors)[: len(donors) // 2])
            ina = np.isin(C.donor, list(half))
            draws = draw_donors(len(donors), rng, a.boot)
            M = multiplicities(donors, C.donor, draws)
            df, _ = compare(C, C, units, rng, Ma=M * ina, Mb=M * ~ina)
            df = df[C.active]
            r = dict(cell=cell, split=i, donors=len(donors), programs=len(df))
            for k in ("d_beta", "d_excess"):
                z = df[k] / df[f"{k}_se"]
                r.update({f"{k}_p05": (df[f"{k}_p"] < 0.05).mean(), f"{k}_p01": (df[f"{k}_p"] < 0.01).mean(),
                          f"{k}_q05": (df[f"{k}_q"] < 0.05).mean(), f"{k}_z_sd": z.std(),
                          f"{k}_cover": ((df[f"{k}_lo"] <= 0) & (df[f"{k}_hi"] >= 0)).mean()})
            rows.append(r)
        d = pd.DataFrame([r for r in rows if r["cell"] == cell])
        print(cell, d.drop(columns=["cell", "split"]).median().round(3).to_dict(), flush=True)
    out = pd.DataFrame(rows)
    path = OUT / "calibration.csv"
    if path.exists():
        old = pd.read_csv(path)
        out = pd.concat([old[~old.cell.isin(a.cells)], out], ignore_index=True)
    out.to_csv(path, index=False)


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare")
    c.add_argument("--a", required=True)
    c.add_argument("--b", required=True)
    k = sub.add_parser("calibrate")
    k.add_argument("--cells", nargs="+", required=True)
    k.add_argument("--splits", type=int, default=20)
    k.add_argument("--boot", type=int, default=500)
    a = ap.parse_args()
    {"compare": cmd_compare, "calibrate": cmd_calibrate}[a.cmd](a)
