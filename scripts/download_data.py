"""Download selected CELLxGene Discover datasets to /data/ANTXR2/raw/.

Resolves download URLs at runtime via the Curation API (collection IDs are stable;
dataset-level asset URLs rotate), per the Phase 1 plan. See config.py for which
collections/datasets are selected and why.

Usage: python download_data.py
"""
import csv
import hashlib
import os
import sys
import time

import requests

from config import CURATION_API_BASE, COLLECTIONS, RAW_DIR

CHUNK = 1024 * 1024 * 8  # 8MB


def fetch_collection(collection_id: str) -> dict:
    url = f"{CURATION_API_BASE}/collections/{collection_id}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.json()


def download_file(url: str, dest_path: str, expected_size: int) -> str:
    """Download with resume support; returns md5 of the final file."""
    tmp_path = dest_path + ".part"
    existing = os.path.getsize(tmp_path) if os.path.exists(tmp_path) else 0

    if os.path.exists(dest_path) and os.path.getsize(dest_path) == expected_size:
        print(f"  already downloaded, verifying checksum: {dest_path}")
    else:
        headers = {"Range": f"bytes={existing}-"} if existing else {}
        mode = "ab" if existing else "wb"
        with requests.get(url, headers=headers, stream=True, timeout=120) as r:
            r.raise_for_status()
            downloaded = existing
            with open(tmp_path, mode) as f:
                for chunk in r.iter_content(chunk_size=CHUNK):
                    f.write(chunk)
                    downloaded += len(chunk)
                    print(f"\r  {dest_path}: {downloaded / 1e9:.2f} / {expected_size / 1e9:.2f} GB", end="", flush=True)
        print()
        os.rename(tmp_path, dest_path)

    md5 = hashlib.md5()
    with open(dest_path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            md5.update(chunk)
    return md5.hexdigest()


MANIFEST_FIELDS = [
    "collection_id", "collection_name", "dataset_id", "dataset_version_id", "title",
    "tissue", "cell_count", "feature_count", "assay", "disease", "raw_data_location",
    "download_url", "filesize_bytes", "md5", "download_timestamp", "download_seconds",
    "local_path",
]


def write_manifest(manifest_path, rows):
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    os.makedirs(RAW_DIR, exist_ok=True)
    manifest_path = os.path.join(RAW_DIR, "manifest.csv")
    rows = []

    for collection_id, spec in COLLECTIONS.items():
        subdir = os.path.join(RAW_DIR, spec["name"])
        os.makedirs(subdir, exist_ok=True)
        print(f"=== Collection {collection_id} ({spec['name']}) ===")
        collection = fetch_collection(collection_id)
        wanted_titles = set(spec["dataset_titles"])
        found_titles = set()

        for ds in collection.get("datasets", []):
            if ds["title"] not in wanted_titles:
                continue
            found_titles.add(ds["title"])
            h5ad_assets = [a for a in ds["assets"] if a["filetype"] == "H5AD"]
            if not h5ad_assets:
                print(f"  WARNING: no H5AD asset for dataset {ds['dataset_id']} ({ds['title']}); skipping")
                continue
            asset = h5ad_assets[0]
            dest_path = os.path.join(subdir, f"{ds['dataset_id']}.h5ad")
            print(f"Downloading {ds['title']} ({ds['dataset_id']}), {asset['filesize'] / 1e9:.2f} GB")
            t0 = time.time()
            md5 = download_file(asset["url"], dest_path, asset["filesize"])
            elapsed = time.time() - t0

            rows.append({
                "collection_id": collection_id,
                "collection_name": spec["name"],
                "dataset_id": ds["dataset_id"],
                "dataset_version_id": ds["dataset_version_id"],
                "title": ds["title"],
                "tissue": ";".join(sorted({t["label"] for t in ds.get("tissue", [])})),
                "cell_count": ds["cell_count"],
                "feature_count": ds.get("feature_count"),
                "assay": ";".join(sorted({a["label"] for a in ds.get("assay", [])})),
                "disease": ";".join(sorted({d["label"] for d in ds.get("disease", [])})),
                "raw_data_location": ds.get("raw_data_location"),
                "download_url": asset["url"],
                "filesize_bytes": asset["filesize"],
                "md5": md5,
                "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "download_seconds": round(elapsed, 1),
                "local_path": dest_path,
            })
            write_manifest(manifest_path, rows)
            print(f"  updated manifest: {manifest_path} ({len(rows)} datasets so far)")

        missing = wanted_titles - found_titles
        if missing:
            print(f"  WARNING: expected dataset titles not found in collection: {missing}", file=sys.stderr)

    print(f"\nFinal manifest: {manifest_path} ({len(rows)} datasets)")


if __name__ == "__main__":
    main()
