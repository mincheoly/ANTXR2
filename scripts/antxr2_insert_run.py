"""Re-run the Wnt pathway x output one-sample pipeline with ANTXR2 and three
detection-matched non-Wnt control genes inserted as if they were pathway genes.

Controls are matched PER CELL TYPE to ANTXR2's raw mean (within +/-20%,
relative) and detection rate (within 2 percentage points), per the project's
per-compartment anchor-matching rule. A single global control list would
re-create the whole-tissue matching bug that forced the earlier retraction.
"""
import os, time, pickle
import numpy as np, pandas as pd

OUT = "antxr2_insert"
os.makedirs(OUT, exist_ok=True)


def run(E, run_cts, CTRL, PATHWAY, OUTPUT, fn):
    t0 = time.time()
    for c in run_cts:
        f = f"{OUT}/{c.replace(' ', '_').replace('+', 'p')}.csv"
        if os.path.exists(f):
            continue
        sub = E[np.asarray(E.obs["Cluster"]).astype(str) == c].copy()
        genes = list(PATHWAY) + ["ANTXR2"] + list(CTRL[c])
        ht, pg, ng = fn(sub, genes, OUTPUT)
        if ht is None:
            continue
        ht["cell_type"] = c
        ht["n_groups"] = ng
        ht.to_csv(f, index=False)
        pg["cell_type"] = c
        pg.to_csv(f.replace(".csv", "_perdonor.csv"), index=False)
        print(f"{c}: {len(ht)} pairs, {ng} groups, {time.time()-t0:.0f}s", flush=True)
    return sorted(os.listdir(OUT))
