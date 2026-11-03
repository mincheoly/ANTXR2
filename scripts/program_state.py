"""Pre-specified tests for the inflamed-vs-uninflamed state contrasts
(planning_summaries/program_coexpression.md, "State contrast plan"), from
program_scores.py runs whose population A is `inflamed`. PC1 weights.

    H1  crypt       bracq:WNT_output   excess > 0, donor-bootstrap CI excludes 0
    H2  crypt       bracq:YAP_fetal    reported, no direction
    H3  enterocyte  discovery:IEG_core excess > 0, CI excludes 0; and
                    DISCOVERY_SET mean excess z > 0 vs random same-size sets

Every fixed module is reported in every run given; a hypothesis is judged
only in its own compartment's primary (vs uninflamed) run.

Writes OUT_DIR/state_contrasts.csv and prints the table.

Usage:
    python program_state.py <run> [<run> ...]
"""
import argparse

import pandas as pd

from program_config import DISCOVERY_SET, OUT_DIR
from program_replicate import library, load, set_vs_random

MODULES = ["bracq:WNT_output", "bracq:YAP_fetal", "discovery:IEG_core"]
COLS = ["n_genes", "rho_a", "rho_b", "delta", "excess", "excess_ci_lo", "excess_ci_hi",
        "excess_z", "reliability_a", "reliability_b"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    a = ap.parse_args()
    rows = []
    for run in a.runs:
        d = load(run)
        for m in MODULES:
            r = dict(run=run, unit=m)
            if m in d.index:
                r.update(d.loc[m, COLS].to_dict())
                r["ci_excludes_0"] = bool(d.loc[m, "excess_ci_lo"] > 0 or d.loc[m, "excess_ci_hi"] < 0)
            else:
                r["note"] = "not comparable"
            rows.append(r)
        lib = library(d)
        S = [u for u in DISCOVERY_SET if u in lib.index]
        obs, p, null = set_vs_random(lib, S, "greater")
        rows.append(dict(run=run, unit=f"set S ({len(S)} of {len(DISCOVERY_SET)})", excess_z=obs,
                         note=f"one-sided p(greater) = {p:.4f}; random-set mean "
                              f"{null.mean():+.2f} sd {null.std():.2f}; "
                              f"{int((lib.loc[S, 'excess'] > 0).sum())}/{len(S)} positive"))
    t = pd.DataFrame(rows)
    t.to_csv(OUT_DIR / "state_contrasts.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
