"""Genomic span (end - start + 1, bp) per gene from GENCODE v44 basic, keyed
by gene symbol and by unversioned Ensembl ID -> program_config.GENE_SPAN.
Used to match anchors to ANTXR2 on length (planning note, "Correction:
a gene-length axis").

Usage:
    python gene_span.py      # expects ROOT/annotation/gencode.v44.basic.annotation.gtf.gz
"""
import gzip
import re

import pandas as pd

from program_config import GENE_SPAN

GTF = GENE_SPAN.parent / "gencode.v44.basic.annotation.gtf.gz"


def main():
    rows = []
    with gzip.open(GTF, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t")
            if p[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^".]+)', p[8]).group(1)
            name = re.search(r'gene_name "([^"]+)', p[8]).group(1)
            rows.append((gid, name, int(p[4]) - int(p[3]) + 1))
    pd.DataFrame(rows, columns=["gene_id", "gene_name", "span"]).to_csv(GENE_SPAN, index=False)
    print(f"wrote {GENE_SPAN}: {len(rows)} genes")


if __name__ == "__main__":
    main()
