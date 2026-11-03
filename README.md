# ANTXR2

Single-cell co-expression analyses of ANTXR2 (CMG2) in gut and other tissues.

- `ANALYSIS_SUMMARY.md`: earlier phases and the project's standing corrections.
- `planning_summaries/program_coexpression.md`: the gene-program co-expression
  pilot and the IBDverse co-expression landscape: decisions, results, how to
  reproduce them on a new host ("Running it") and next steps.
- `genesets/`: that pilot's program library (`programs_long.csv`) and its
  CytoSig input; the library every result used, kept in the repo so reruns
  use the same one.
- `scripts/program_*.py`: that pilot's pipeline (`program_config.py` holds the
  settings; set `ANTXR2_ROOT` if the data root is not `/data/antxr2`).
