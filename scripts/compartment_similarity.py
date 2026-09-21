"""Per-compartment and global ANTXR2 co-expression similarity.

Loads every per-population profile under coexpr_out (fibroblasts) and
coexpr_comp (epithelial / endothelial / immune), computes the same family of
similarity matrices used for the fibroblast arm, applies the reliability
(attenuation) correction, and writes per-compartment plus pooled matrices.

Population keys are "<compartment>_<tissue>"; fibroblast outputs carry no
prefix and are relabelled "fibroblast_<tissue>" here so the pooled matrix has
one consistent namespace.
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

import similarity as S

COMPARTMENTS = ["fibroblast", "epithelial", "endothelial", "immune"]
REL_FLOOR = 0.20


def load_all(fibdir, compdir):
    prof, meta = {}, []
    for d, pref in [(fibdir, "fibroblast"), (compdir, None)]:
        if not d or not os.path.isdir(d):
            continue
        p = S.load_profiles(d)
        for k, v in p.items():
            key = f"{pref}_{k}" if pref else k
            prof[key] = v
        for f in sorted(os.listdir(d)):
            if not f.endswith(".json"):
                continue
            m = json.load(open(os.path.join(d, f)))
            m["population"] = (f"{pref}_{m['population']}" if pref
                               else m["population"])
            meta.append(m)
    return prof, pd.DataFrame(meta)


def compartment_of(pop):
    for c in COMPARTMENTS:
        if pop.startswith(c + "_"):
            return c
    return "other"


def tissue_of(pop):
    return pop.split("_", 1)[1] if "_" in pop else pop


def write_block(prof, outdir, tag):
    os.makedirs(outdir, exist_ok=True)
    if len(prof) < 3:
        print(f"skip {tag}: only {len(prof)} populations", flush=True)
        return None
    M = S.similarity(prof, value="corr_coef")
    D = S.similarity(prof, value="det_rate")
    R, K = S.disattenuated(prof)
    rel = pd.Series({p: S.reliability(prof[p]) for p in prof}, name="reliability")
    for k, m in M.items():
        m.to_csv(os.path.join(outdir, f"{tag}_sim_{k}.csv"))
    D["pearson_all"].to_csv(os.path.join(outdir, f"{tag}_det_pearson_all.csv"))
    R.to_csv(os.path.join(outdir, f"{tag}_sim_disattenuated.csv"))
    rel.to_csv(os.path.join(outdir, f"{tag}_profile_reliability.csv"))
    print(f"{tag}: {len(prof)} populations, "
          f"{(rel >= REL_FLOOR).sum()} above reliability floor", flush=True)
    return dict(sim=M, det=D, dis=R, rel=rel)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fibdir", default=r"C:\Data\ANTXR2_workspace\coexpr_out")
    ap.add_argument("--compdir", default=r"C:\Data\ANTXR2_workspace\coexpr_comp")
    ap.add_argument("--outdir", default=r"C:\Data\ANTXR2_workspace\results")
    a = ap.parse_args()

    prof, meta = load_all(a.fibdir, a.compdir)
    meta["compartment"] = meta.population.map(compartment_of)
    meta["tissue"] = meta.population.map(tissue_of)
    os.makedirs(a.outdir, exist_ok=True)
    meta.to_csv(os.path.join(a.outdir, "all_run_summary.csv"), index=False)
    print("loaded", len(prof), "profiles:",
          meta.compartment.value_counts().to_dict(), flush=True)

    for c in COMPARTMENTS:
        sub = {k: v for k, v in prof.items() if compartment_of(k) == c}
        write_block(sub, a.outdir, c)
    write_block(prof, a.outdir, "global")


if __name__ == "__main__":
    main()
