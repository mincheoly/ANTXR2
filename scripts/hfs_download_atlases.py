"""Download the CELLxGENE atlases used by the HFS cross-tissue arm.

Usage
-----
    python scripts/hfs_download_atlases.py              # all datasets
    python scripts/hfs_download_atlases.py gut skin     # a subset

Downloads go to ``HFS_RAW_DIR`` as ``<key>.h5ad``. Existing files are skipped on
size match and otherwise re-fetched, so an interrupted run is resumable.

After each download the file's own ``uns["citation"]`` is parsed and the
dataset-version UUID inside it is checked against the registry in
``hfs_config.DATASETS``. A completed transfer is not evidence of correct
content; this check is what makes the registry and the bytes on disk agree.
"""
import sys
import re
import urllib.request

import h5py

from hfs_config import DATASETS, CXG_BASE, HFS_RAW_DIR

UUID_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")


def citation_version_id(path):
    """Dataset-version UUID declared inside the file, or None."""
    with h5py.File(path, "r") as f:
        uns = f.get("uns")
        if uns is None or "citation" not in uns:
            return None
        v = uns["citation"][()]
        cit = v.decode() if isinstance(v, bytes) else str(v)
    ids = UUID_RE.findall(cit)
    return ids[0] if ids else None


def download(key, spec):
    dst = HFS_RAW_DIR / f"{key}.h5ad"
    url = CXG_BASE.format(version_id=spec["version_id"])
    if dst.exists() and dst.stat().st_size > 0.5e9:
        print(f"skip     {key:12s} {dst.stat().st_size/1e9:5.2f} GB already on disk", flush=True)
    else:
        print(f"fetch    {key:12s} {spec['gb']:5.2f} GB  {spec['title'][:60]}", flush=True)
        tmp = dst.with_suffix(".h5ad.part")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(dst)
        print(f"  wrote  {dst.stat().st_size/1e9:5.2f} GB", flush=True)

    got = citation_version_id(dst)
    if got is None:
        print(f"  WARN   {key}: no citation in uns, cannot verify identity")
    elif got != spec["version_id"]:
        raise SystemExit(
            f"{key}: file declares dataset version {got} but registry expects "
            f"{spec['version_id']} — refusing to proceed on mismatched data")
    else:
        print(f"  verify {key}: dataset version matches registry")


def main(keys=None):
    keys = keys or list(DATASETS)
    unknown = [k for k in keys if k not in DATASETS]
    if unknown:
        raise SystemExit(f"unknown dataset keys: {unknown}\nknown: {list(DATASETS)}")
    for k in keys:
        download(k, DATASETS[k])


if __name__ == "__main__":
    main(sys.argv[1:] or None)
