"""Assemble the gene-program library: Hallmark, PROGENy, CollecTRI, CytoSig,
Matrisome -> one long table <GENESET_DIR>/programs_long.csv with columns

    collection, program, gene, sign, weight

`sign` is +1/-1 (Decision 5: published signs enter the equal-weight baseline)
and `weight` is the resource's own signed weight where one exists (PROGENy
footprint weight, CytoSig centroid log-fold-change), else NaN. CollecTRI's
published weight is just its sign, so its third scheme would duplicate the
signed baseline and is left NaN on purpose.

Prioritises transcriptionally coherent resources over GO/KEGG. Matrisome uses
the eight Naba categories only: the two umbrella sets (NABA_MATRISOME, 1,694
genes; NABA_MATRISOME_ASSOCIATED, 1,212) are too broad to be one program, and
the tumour-derived NABA_* sets are not gut biology.

ANTXR2 is removed from every program here, before any filtering, so it can
never contribute to its own score.

Mapping to the dataset's genes (and the expression floor) happens in
program_scores.py, which has the var table.

Usage:
    python program_genesets.py
"""
import numpy as np
import pandas as pd
import decoupler as dc

from program_config import CYTOSIG_TOP, GENESET_DIR, PROGENY_TOP, TARGET

CYTOSIG_FILE = GENESET_DIR / "raw_cytosig" / "signature.centroid"
CYTOSIG_URL = ("https://raw.githubusercontent.com/data2intelligence/CytoSig/"
               "master/CytoSig/signature.centroid")
MATRISOME = ["NABA_BASEMENT_MEMBRANES", "NABA_COLLAGENS", "NABA_CORE_MATRISOME",
             "NABA_ECM_AFFILIATED", "NABA_ECM_GLYCOPROTEINS", "NABA_ECM_REGULATORS",
             "NABA_PROTEOGLYCANS", "NABA_SECRETED_FACTORS"]


def hallmark():
    h = dc.op.hallmark(organism="human")
    return pd.DataFrame({"collection": "hallmark", "program": h.source,
                         "gene": h.target, "sign": 1.0, "weight": np.nan})


def progeny():
    p = dc.op.progeny(organism="human", top=PROGENY_TOP)
    return pd.DataFrame({"collection": "progeny", "program": p.source,
                         "gene": p.target, "sign": np.sign(p.weight),
                         "weight": p.weight})


def collectri():
    c = dc.op.collectri(organism="human")
    return pd.DataFrame({"collection": "collectri", "program": c.source,
                         "gene": c.target, "sign": np.sign(c.weight),
                         "weight": np.nan})


def cytosig():
    if not CYTOSIG_FILE.exists():
        raise FileNotFoundError(f"{CYTOSIG_FILE} missing -- fetch {CYTOSIG_URL}")
    m = pd.read_csv(CYTOSIG_FILE, sep="\t", index_col=0)
    rows = []
    for cyt in m.columns:
        w = m[cyt].dropna()
        w = w[w != 0]
        top = w.loc[w.abs().sort_values(ascending=False).index[:CYTOSIG_TOP]]
        rows.append(pd.DataFrame({"collection": "cytosig", "program": cyt,
                                  "gene": top.index, "sign": np.sign(top.values),
                                  "weight": top.values}))
    return pd.concat(rows, ignore_index=True)


def matrisome():
    m = dc.op.resource("MSigDB", organism="human")
    m = m[m.geneset.isin(MATRISOME)]
    return pd.DataFrame({"collection": "matrisome", "program": m.geneset,
                         "gene": m.genesymbol, "sign": 1.0, "weight": np.nan})


def main():
    lib = pd.concat([hallmark(), progeny(), collectri(), cytosig(), matrisome()],
                    ignore_index=True)
    lib = lib.drop_duplicates(["collection", "program", "gene"])
    n_target = int((lib.gene == TARGET).sum())
    lib = lib[lib.gene != TARGET]
    lib["program_id"] = lib.collection + ":" + lib.program
    out = GENESET_DIR / "programs_long.csv"
    lib.to_csv(out, index=False)
    summ = lib.groupby("collection").agg(n_programs=("program", "nunique"),
                                         n_rows=("gene", "size"),
                                         median_size=("program", lambda s: s.value_counts().median()))
    print(summ.to_string())
    print(f"{TARGET} removed from {n_target} program memberships -> {out}")


if __name__ == "__main__":
    main()
