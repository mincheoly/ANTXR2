"""Download 10x's public PBMC reference datasets, chemistry-matched to the two
assays in our Elmentaite2021 working set, for the PBMC-scaled capture-rate
column (see capture_rate.pbmc_scaled_capture_rate).

Follows the same resumable-download + manifest pattern as download_skin_data.py.
These files are small (~30-38MB each), so no chunk/resume gymnastics are really
needed, but recording provenance in a manifest matters -- the derived capture
rate is only interpretable if the reference it came from is pinned.

Usage: python download_pbmc_reference.py
"""
import csv
import hashlib
import os
import time

import requests

from config import PBMC_REFERENCE_DIR, PBMC_REFERENCES

CHUNK = 1024 * 1024 * 8

MANIFEST_FIELDS = [
    "assay", "name", "kind", "url", "filename", "filesize_bytes", "md5",
    "download_timestamp", "download_seconds", "local_path",
]


def download_file(url, dest_path):
    """Download to dest_path (skipping if already present with nonzero size);
    returns (md5, size)."""
    if os.path.exists(dest_path) and os.path.getsize(dest_path) > 0:
        print(f"  already present, verifying: {dest_path}")
    else:
        tmp = dest_path + ".part"
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            done = 0
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=CHUNK):
                    f.write(chunk)
                    done += len(chunk)
                    print(f"\r  {os.path.basename(dest_path)}: {done/1e6:.1f}"
                          f"{'/' + format(total/1e6, '.1f') if total else ''} MB",
                          end="", flush=True)
        print()
        os.rename(tmp, dest_path)

    md5 = hashlib.md5()
    with open(dest_path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            md5.update(chunk)
    return md5.hexdigest(), os.path.getsize(dest_path)


def main():
    os.makedirs(PBMC_REFERENCE_DIR, exist_ok=True)
    rows = []
    for assay, spec in PBMC_REFERENCES.items():
        dest = os.path.join(PBMC_REFERENCE_DIR, spec["filename"])
        print(f"=== {assay}: {spec['name']} ===")
        t0 = time.time()
        md5, size = download_file(spec["url"], dest)
        rows.append({
            "assay": assay, "name": spec["name"], "kind": spec["kind"],
            "url": spec["url"], "filename": spec["filename"],
            "filesize_bytes": size, "md5": md5,
            "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "download_seconds": round(time.time() - t0, 1),
            "local_path": dest,
        })
        print(f"  {size/1e6:.1f} MB, md5={md5}")

    manifest = os.path.join(PBMC_REFERENCE_DIR, "manifest.csv")
    with open(manifest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {manifest} ({len(rows)} references)")


if __name__ == "__main__":
    main()
