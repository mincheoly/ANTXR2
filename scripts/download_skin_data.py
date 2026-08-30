"""Download the skin fibroblast atlas (Steele et al., Nat. Immunol. 2025) integrated
h5ad to SKIN_FIBROBLAST_DIR.

See config.py for how SKIN_FIBROBLAST_URL was resolved (not a CELLxGene Discover
dataset -- prompts/skin_means.md flagged this as a download blocker requiring a
stable accession rather than a guessed/scraped URL; resolved via the paper's own
public web-portal CMS API, see the comment above SKIN_FIBROBLAST_URL).

Usage: python download_skin_data.py
"""
import csv
import hashlib
import os
import time

import requests

from config import (
    SKIN_FIBROBLAST_DIR, SKIN_FIBROBLAST_EXPECTED_MD5, SKIN_FIBROBLAST_EXPECTED_SIZE,
    SKIN_FIBROBLAST_URL,
)

CHUNK = 1024 * 1024 * 8  # 8MB

MANIFEST_FIELDS = [
    "source", "study_slug", "dataset_name", "download_url", "filesize_bytes",
    "expected_md5", "md5", "md5_verified", "download_timestamp", "download_seconds",
    "local_path",
]


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


def main():
    os.makedirs(SKIN_FIBROBLAST_DIR, exist_ok=True)
    dest_path = os.path.join(SKIN_FIBROBLAST_DIR, "adata_webportal.h5ad")
    manifest_path = os.path.join(SKIN_FIBROBLAST_DIR, "manifest.csv")

    print(f"Downloading skin fibroblast atlas integrated h5ad ({SKIN_FIBROBLAST_EXPECTED_SIZE / 1e9:.2f} GB)")
    print(f"  from: {SKIN_FIBROBLAST_URL}")
    t0 = time.time()
    md5 = download_file(SKIN_FIBROBLAST_URL, dest_path, SKIN_FIBROBLAST_EXPECTED_SIZE)
    elapsed = time.time() - t0

    md5_verified = md5 == SKIN_FIBROBLAST_EXPECTED_MD5
    print(f"  md5: {md5} (expected {SKIN_FIBROBLAST_EXPECTED_MD5}, verified={md5_verified})")

    row = {
        "source": "haniffalab/skin_fibroblast_atlas web portal CMS (strapi-api-dot-haniffa-lab.nw.r.appspot.com)",
        "study_slug": "skin-fibroblast",
        "dataset_name": "Integrated atlas",
        "download_url": SKIN_FIBROBLAST_URL,
        "filesize_bytes": os.path.getsize(dest_path),
        "expected_md5": SKIN_FIBROBLAST_EXPECTED_MD5,
        "md5": md5,
        "md5_verified": md5_verified,
        "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "download_seconds": round(elapsed, 1),
        "local_path": dest_path,
    }
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerow(row)

    print(f"\nWrote manifest: {manifest_path}")
    if not md5_verified:
        print("WARNING: md5 mismatch -- expected checksum was read off the GCS ETag header, "
              "which for objects uploaded as a single part equals the MD5. If this object was "
              "uploaded multipart the ETag is not an MD5 and this warning is a false alarm -- "
              "re-check manually before trusting the file.")


if __name__ == "__main__":
    main()
