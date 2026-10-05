"""Pick expression/detection-matched anchor genes per population for ANTXR2.

Matching is within each population (compartment x tissue), never pooled, per the
per-compartment rule in the detection-matched-coexpression skill. Tolerances
are relative on the raw mean (+/- 20%) plus an absolute detection-rate
window (+/- 2 pp), following antxr2_insert_run.py: an absolute log1p
tolerance admits genes spanning orders of magnitude at ANTXR2's low expression.

Candidates are genes tested in that population's scan (so they clear the same
memento filters), excluding ribosomal/mitochondrial genes and the project's
own panel genes. If fewer than N_ANCHORS candidates match, the tolerances are
widened stepwise and the widening is recorded.

With --npz, candidates and ANTXR2's own mean/detection are computed directly
from population matrices (e.g. the filtered output of endothelial_decontam.py),
for populations that have no scan yet. Candidates must then clear
CANDIDATE_FLOOR on the population raw mean, an approximation of memento's
per-group filter; an anchor failing memento's filter simply drops out of the
scan. Populations whose ANTXR2 raw mean is below --antxr2-floor are skipped
and reported, not rescued (house reliability floor 0.07).

Usage: python coexpr_anchor_select.py <out.csv> --pops <scan.csv> ... --npz-dirs <dir> ...
       python coexpr_anchor_select.py <out.csv> --npz <pop.npz> ... [--antxr2-floor 0.07]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coexpr_pathway_gsea import population_name

N_ANCHORS = 5
CANDIDATE_FLOOR = 0.07
WIDEN = [(0.20, 0.02), (0.30, 0.03), (0.50, 0.05)]
SEED = 20261003
EXCLUDE_RE = r"^(RP[LS]\d|MRP[LS]\d|MT-)"
PANEL = {"ANTXR1", "ANTXR2", "MRC2", "CTSB", "CTSK", "MMP14", "TIMP2", "LAMP1",
         "COL6A1", "COL6A2", "COL6A3"}


def candidates(scan_csv):
    """-> population name, ANTXR2 raw mean (from the run metadata, since ANTXR2
    is the target rather than a row of the scan), candidate anchor table."""
    df = pd.read_csv(scan_csv).dropna(subset=["raw_mean", "det_rate"])
    df = df.drop_duplicates("gene").set_index("gene")
    meta = pd.read_json(scan_csv[:-len(".csv")] + ".json", typ="series")
    cand = df[~df.index.str.match(EXCLUDE_RE) & ~df.index.isin(PANEL)]
    return population_name(scan_csv), float(meta["antxr2_raw_mean"]), cand


def candidates_from_npz(npz):
    """-> population name, ANTXR2 raw mean, ANTXR2 detection, candidate table."""
    z = np.load(npz, allow_pickle=False)
    n, g = int(z["shape"][0]), int(z["shape"][1])
    genes = z["genes"].astype(str)
    raw_mean = np.bincount(z["indices"], weights=z["data"], minlength=g) / n
    det_rate = np.bincount(z["indices"], minlength=g) / n
    df = pd.DataFrame({"raw_mean": raw_mean, "det_rate": det_rate}, index=genes)
    df = df[~df.index.duplicated()]
    t = df.loc["ANTXR2"]
    cand = df[(df.raw_mean >= CANDIDATE_FLOOR) & ~df.index.str.match(EXCLUDE_RE) & ~df.index.isin(PANEL)]
    return os.path.basename(npz)[:-4], float(t.raw_mean), float(t.det_rate), cand


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--pops", nargs="+")
    ap.add_argument("--npz-dirs", nargs="+",
                    help="dirs holding the population .npz files (for ANTXR2 detection)")
    ap.add_argument("--npz", nargs="+", help="population matrices to select from directly")
    ap.add_argument("--antxr2-floor", type=float, default=0.07)
    ap.add_argument("--n-anchors", type=int, default=N_ANCHORS,
                    help="more anchors give an empirical percentile within a single population")
    a = ap.parse_args()
    if bool(a.pops) == bool(a.npz):
        ap.error("give exactly one of --pops (with --npz-dirs) or --npz")
    rng = np.random.default_rng(SEED)
    rows = []
    jobs = a.npz if a.npz else a.pops
    for p in jobs:
        if a.npz:
            pop, tm, td, cand = candidates_from_npz(p)
            if tm < a.antxr2_floor:
                print(f"SKIP {pop}: ANTXR2 raw mean {tm:.3f} < floor {a.antxr2_floor}")
                continue
        else:
            pop, tm, cand = candidates(p)
            td = antxr2_detection(pop, a.npz_dirs)
        for rel, dt in WIDEN:
            ok = cand[(cand.raw_mean.between(tm * (1 - rel), tm * (1 + rel)))
                      & (cand.det_rate.between(td - dt, td + dt))]
            if len(ok) >= a.n_anchors:
                break
        chosen = ok.sample(min(a.n_anchors, len(ok)), random_state=rng.integers(1 << 31))
        for g, r in chosen.iterrows():
            rows.append(dict(population=pop, anchor=g, raw_mean=r.raw_mean, det_rate=r.det_rate,
                             antxr2_raw_mean=tm, antxr2_det_rate=td, rel_tol=rel, det_tol=dt,
                             n_candidates=len(ok)))
        if len(chosen) < a.n_anchors:
            print(f"WARN {pop}: only {len(chosen)} anchors at widest tolerance")
    out = pd.DataFrame(rows)
    out.to_csv(a.out, index=False)
    s = out.groupby("population").agg(n=("anchor", "size"), rel_tol=("rel_tol", "first"),
                                       n_cand=("n_candidates", "first"))
    print(s.rel_tol.value_counts().to_string())
    print(f"{s.n.sum()} anchors over {len(s)} populations; min candidates {s.n_cand.min()}")


def antxr2_detection(pop, npz_dirs):
    """ANTXR2 detection rate from the population matrix, matching the scan's
    det_rate definition (fraction of all cells in the population with a count)."""
    name = pop[len("fibroblast_"):] if not any(
        os.path.exists(os.path.join(d, pop + ".npz")) for d in npz_dirs) else pop
    for d in npz_dirs:
        f = os.path.join(d, name + ".npz")
        if os.path.exists(f):
            z = np.load(f, allow_pickle=False)
            j = list(z["genes"].astype(str)).index("ANTXR2")
            indices, n = z["indices"], int(z["shape"][0])
            return float(np.unique(np.searchsorted(z["indptr"], np.flatnonzero(indices == j),
                                                   side="right")).size / n)
    raise FileNotFoundError(pop)


if __name__ == "__main__":
    main()
