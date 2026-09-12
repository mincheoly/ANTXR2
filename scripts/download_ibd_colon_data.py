"""Download the IBD colon atlas (Smillie et al., Cell 2019, SCP259) raw MTX
matrices, metadata, and pre-computed Seurat objects, from the Human Cell Atlas
Data Coordination Platform.

See config.py for why HCA is the source instead of SCP259 directly (SCP259 is
login-gated; HCA hosts the identical, genuinely open-access file set, verified
via unauthenticated round-trips of both a tiny and a >1GB file). File uuids are
resolved at runtime from the live Azul catalog (per-file uuids are stable within
a project, but which catalog is "current" rotates over time) -- same
resolve-at-runtime-off-a-stable-id pattern as download_data.py's use of
CELLxGene collection ids.

Usage: python download_ibd_colon_data.py
"""
import csv
import hashlib
import os
import time

import requests

from config import (
    AZUL_BASE, IBD_COLON_BARCODES_FILES, IBD_COLON_COMPARTMENTS, IBD_COLON_DIR,
    IBD_COLON_GENES_FILES, IBD_COLON_HCA_PROJECT_ID, IBD_COLON_MATRIX_FILES,
    IBD_COLON_META_FILE, IBD_COLON_SEURAT_RDS_FILES, IBD_COLON_SUBSETS_FILE,
    IBD_COLON_SUBSETS_URL,
)

CHUNK = 1024 * 1024 * 8  # 8MB

MANIFEST_FIELDS = [
    "source", "file_name", "hca_file_uuid", "hca_version", "expected_size_bytes",
    "filesize_bytes", "md5", "source_etag", "etag_is_md5", "download_url",
    "download_timestamp", "download_seconds", "local_path",
]


def get_default_catalog():
    resp = requests.get(f"{AZUL_BASE}/index/catalogs", timeout=30)
    resp.raise_for_status()
    return resp.json()["default_catalog"]


def list_project_files(catalog, project_id):
    """Return {file_name: {uuid, version, size}} for every file in the project."""
    filters = '{"projectId":{"is":["%s"]}}' % project_id
    resp = requests.get(
        f"{AZUL_BASE}/index/files",
        params={"catalog": catalog, "filters": filters, "size": 500},
        timeout=60,
    )
    resp.raise_for_status()
    data = resp.json()
    out = {}
    for hit in data["hits"]:
        f = hit["files"][0]
        out[f["name"]] = {"uuid": f["uuid"], "version": f["version"], "size": f["size"]}
    return out


def resolve_download_url(catalog, file_uuid, version):
    """Follow the Azul repository/files redirect to a concrete (no-auth) URL."""
    resp = requests.get(
        f"{AZUL_BASE}/repository/files/{file_uuid}",
        params={"catalog": catalog, "version": version},
        allow_redirects=False,
        timeout=60,
    )
    resp.raise_for_status()
    return resp.headers["Location"]


def download_file(url: str, dest_path: str, expected_size: int) -> tuple:
    """Download with resume support; returns (md5, source_etag)."""
    tmp_path = dest_path + ".part"
    existing = os.path.getsize(tmp_path) if os.path.exists(tmp_path) else 0

    if os.path.exists(dest_path) and (expected_size is None or os.path.getsize(dest_path) == expected_size):
        print(f"  already downloaded: {dest_path}")
        etag = None
    else:
        headers = {"Range": f"bytes={existing}-"} if existing else {}
        mode = "ab" if existing else "wb"
        with requests.get(url, headers=headers, stream=True, timeout=120) as r:
            r.raise_for_status()
            etag = r.headers.get("ETag", "").strip('"')
            downloaded = existing
            with open(tmp_path, mode) as f:
                for chunk in r.iter_content(chunk_size=CHUNK):
                    f.write(chunk)
                    downloaded += len(chunk)
                    total = f"{expected_size / 1e9:.2f}" if expected_size else "?"
                    print(f"\r  {dest_path}: {downloaded / 1e9:.2f} / {total} GB", end="", flush=True)
        print()
        os.rename(tmp_path, dest_path)

    md5 = hashlib.md5()
    with open(dest_path, "rb") as f:
        for chunk in iter(lambda: f.read(CHUNK), b""):
            md5.update(chunk)
    return md5.hexdigest(), etag


def write_manifest(manifest_path, rows):
    with open(manifest_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    os.makedirs(IBD_COLON_DIR, exist_ok=True)
    manifest_path = os.path.join(IBD_COLON_DIR, "manifest.csv")
    rows = []

    catalog = get_default_catalog()
    print(f"Azul default catalog: {catalog}")
    files = list_project_files(catalog, IBD_COLON_HCA_PROJECT_ID)
    print(f"Found {len(files)} files on HCA for project {IBD_COLON_HCA_PROJECT_ID}:")
    for name, info in sorted(files.items()):
        print(f"  {name}: {info['size'] / 1e6:.1f} MB")

    wanted_names = (
        list(IBD_COLON_MATRIX_FILES.values())
        + list(IBD_COLON_GENES_FILES.values())
        + list(IBD_COLON_BARCODES_FILES.values())
        + list(IBD_COLON_SEURAT_RDS_FILES.values())
        + [IBD_COLON_META_FILE]
    )
    missing = [n for n in wanted_names if n not in files]
    if missing:
        raise RuntimeError(f"expected files not found on HCA: {missing}")

    for name in wanted_names:
        info = files[name]
        dest_path = os.path.join(IBD_COLON_DIR, name)
        url = resolve_download_url(catalog, info["uuid"], info["version"])
        print(f"Downloading {name} ({info['size'] / 1e6:.1f} MB)")
        t0 = time.time()
        md5, etag = download_file(url, dest_path, info["size"])
        elapsed = time.time() - t0
        # S3 ETags for multipart-uploaded objects (large files) are NOT the
        # object's MD5 (they have a "-<n_parts>" suffix) -- only single-part
        # uploads' ETags equal MD5. Recorded either way for provenance, but only
        # trust etag_is_md5==True as an actual integrity check; otherwise rely on
        # the exact-filesize match already enforced by download_file().
        etag_is_md5 = bool(etag) and "-" not in etag and etag == md5
        rows.append({
            "source": "HCA DCP (Azul) -- project cd61771b-661a-4e19-b269-6e5d95350de6",
            "file_name": name,
            "hca_file_uuid": info["uuid"],
            "hca_version": info["version"],
            "expected_size_bytes": info["size"],
            "filesize_bytes": os.path.getsize(dest_path),
            "md5": md5,
            "source_etag": etag,
            "etag_is_md5": etag_is_md5,
            "download_url": url.split("?")[0],  # strip the signed query string
            "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "download_seconds": round(elapsed, 1),
            "local_path": dest_path,
        })
        write_manifest(manifest_path, rows)

    # cell_subsets.txt: small, public, from the authors' own GitHub repo, not HCA.
    dest_path = os.path.join(IBD_COLON_DIR, IBD_COLON_SUBSETS_FILE)
    print(f"Downloading {IBD_COLON_SUBSETS_FILE} (from GitHub)")
    t0 = time.time()
    md5, etag = download_file(IBD_COLON_SUBSETS_URL, dest_path, None)
    rows.append({
        "source": "github.com/cssmillie/ulcerative_colitis (authors' analysis repo)",
        "file_name": IBD_COLON_SUBSETS_FILE,
        "hca_file_uuid": "",
        "hca_version": "",
        "expected_size_bytes": "",
        "filesize_bytes": os.path.getsize(dest_path),
        "md5": md5,
        "source_etag": etag,
        "etag_is_md5": False,
        "download_url": IBD_COLON_SUBSETS_URL,
        "download_timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "download_seconds": round(time.time() - t0, 1),
        "local_path": dest_path,
    })
    write_manifest(manifest_path, rows)

    print(f"\nDone. Wrote manifest: {manifest_path} ({len(rows)} files)")


if __name__ == "__main__":
    main()
