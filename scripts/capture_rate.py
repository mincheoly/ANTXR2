"""Capture-rate (memento `q`) derivation. Pure functions, no I/O, so the math is
testable in isolation from the h5ad plumbing.

Two independent estimates are produced, deliberately kept separate rather than
averaged -- if they disagree that is information about the assumptions, not
noise to be smoothed away:

1. `chem_capture_rate` -- literature/chemistry-based: memento's own broad
   droplet rate, corrected for how thoroughly the library was sequenced.
2. `pbmc_scaled_capture_rate` -- empirical: transfers a reference dataset's q
   by comparing UMI depth in matched immune populations.
"""
import numpy as np
from scipy.optimize import brentq


def saturation_to_detected_fraction(saturation):
    """Convert 10x 'sequencing saturation' to the fraction of library molecules
    actually observed.

    These are NOT the same quantity, and multiplying a capture rate by
    saturation directly is wrong. 10x defines sequencing saturation as the
    fraction of reads that are duplicates of an already-seen UMI -- a statement
    about read redundancy, not about molecules lost. Under Poisson sampling of
    reads over molecules with mean lambda reads/molecule:

        saturation        = 1 - (1 - e^-lambda) / lambda
        detected_fraction = 1 - e^-lambda

    So we invert the first for lambda, then evaluate the second. At saturation
    0.85 this gives lambda ~ 6.65 and detected ~ 0.9987 -- i.e. at high
    saturation essentially every library molecule has been seen, and the
    correction is nearly inert. The naive product (q * 0.85) would understate q
    by ~15%.

    Returns (lambda, detected_fraction).
    """
    if not (0.0 < saturation < 1.0):
        raise ValueError(f"saturation must be in (0, 1), got {saturation}")

    def f(lam):
        return (1.0 - (1.0 - np.exp(-lam)) / lam) - saturation

    # f is monotonically increasing in lambda: ~0 as lambda -> 0, -> 1 as
    # lambda -> inf. Bracket wide enough to cover any realistic saturation.
    lam = brentq(f, 1e-9, 500.0, xtol=1e-12)
    return lam, float(1.0 - np.exp(-lam))


def chem_capture_rate(assay, saturation, q_chem_by_assay, q_chem_default):
    """q = (chemistry capture efficiency) x (fraction of library molecules seen).

    Returns (q, lambda, detected_fraction) so the derivation stays auditable
    rather than collapsing to one opaque number.
    """
    q_chem = q_chem_by_assay.get(assay, q_chem_default)
    lam, detected = saturation_to_detected_fraction(saturation)
    return q_chem * detected, lam, detected


def gate_masks(counts_csr, var_names, gates):
    """Marker-positive gating, applied identically to our data and the PBMC
    reference so the population definitions cannot drift between them.

    Deliberately crude (>=1 count of each positive marker, 0 of each negative)
    -- the goal is a population whose UMI depth is comparable across datasets,
    not a defensible annotation. Genes absent from a dataset's var are treated
    as unsatisfiable and that gate is skipped, rather than silently passing.

    Returns {gate_name: boolean mask over cells}. Gates that cannot be
    evaluated are omitted from the result.
    """
    import scipy.sparse as sp

    idx = {g: i for i, g in enumerate(var_names)}
    out = {}
    for name, spec in gates.items():
        needed = list(spec["pos"]) + list(spec["neg"])
        if any(g not in idx for g in needed):
            continue
        mask = np.ones(counts_csr.shape[0], dtype=bool)
        for g in spec["pos"]:
            col = counts_csr[:, idx[g]]
            col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
            mask &= col > 0
        for g in spec["neg"]:
            col = counts_csr[:, idx[g]]
            col = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
            mask &= col == 0
        out[name] = mask
    return out


def median_umi_by_gate(total_umi, masks, min_cells=50):
    """Median total UMI per cell within each gate. Gates with fewer than
    `min_cells` are dropped rather than contributing a noisy median."""
    return {
        name: float(np.median(total_umi[m]))
        for name, m in masks.items()
        if m.sum() >= min_cells
    }


def pbmc_scaled_capture_rate(our_umi_by_gate, ref_umi_by_gate, q_ref):
    """Transfer a reference dataset's capture rate by UMI-depth ratio.

    For a given cell type, detected UMIs scale with capture efficiency, so
    q_ours / q_ref ~ medUMI_ours / medUMI_ref. We take the median ratio across
    matched gates so one atypical population can't dominate.

    Key limitation, which is why this is a second column and not a replacement
    for the chemistry-based one: it assumes a given immune cell type carries the
    same absolute mRNA content in gut tissue as in blood. Tissue-resident
    lymphocytes are plausibly more activated and RNA-richer than circulating
    ones, which would inflate the apparent q.

    Returns (q, ratios_by_gate). Raises if no gate is shared.
    """
    shared = sorted(set(our_umi_by_gate) & set(ref_umi_by_gate))
    if not shared:
        raise ValueError(
            f"no shared gates between ours {sorted(our_umi_by_gate)} "
            f"and reference {sorted(ref_umi_by_gate)}"
        )
    ratios = {g: our_umi_by_gate[g] / ref_umi_by_gate[g] for g in shared}
    return q_ref * float(np.median(list(ratios.values()))), ratios
