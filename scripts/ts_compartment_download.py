"""Download the Tabula Sapiens Epithelium / Immune / Endothelium H5ADs.

Siblings of the stromal object already on disk (same collection, same schema),
fetched from the CELLxGENE curation API. Resumable: an existing file whose size
already matches the API record is skipped.
"""
import json
import os
import sys
import time
import urllib.request

DEST = r"C:\Data\ANTXR2_workspace\data"
COLLECTION = "e5f58829-1a66-40b5-a624-9046778e74f5"
WANT = {"Tabula Sapiens - Epithelium": "ts_epithelium.h5ad",
        "Tabula Sapiens - Endothelium": "ts_endothelium.h5ad",
        "Tabula Sapiens - Immune": "ts_immune.h5ad"}


def main():
    col = json.load(urllib.request.urlopen(
        f"https://api.cellxgene.cziscience.com/curation/v1/collections/{COLLECTION}", timeout=120))
    todo = []
    for d in col["datasets"]:
        if d["title"] in WANT:
            a = [x for x in d["assets"] if x["filetype"] == "H5AD"][0]
            todo.append((d["title"], WANT[d["title"]], a["url"], int(a["filesize"])))
    todo.sort(key=lambda x: x[3])
    for title, fn, url, size in todo:
        out = os.path.join(DEST, fn)
        if os.path.exists(out) and os.path.getsize(out) == size:
            print(f"SKIP {fn} ({size/1e9:.2f} GB already present)", flush=True)
            continue
        t0 = time.time()
        tmp = out + ".part"
        with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as f:
            got = 0
            last = 0
            while True:
                chunk = r.read(1 << 22)
                if not chunk:
                    break
                f.write(chunk)
                got += len(chunk)
                if got - last > 2e9:
                    last = got
                    print(f"  {fn} {got/1e9:.1f}/{size/1e9:.1f} GB "
                          f"{got/1e6/(time.time()-t0):.0f} MB/s", flush=True)
        os.replace(tmp, out)
        ok = os.path.getsize(out) == size
        print(f"{'DONE' if ok else 'SIZE MISMATCH'} {fn} {os.path.getsize(out)/1e9:.2f} GB "
              f"in {(time.time()-t0)/60:.1f} min", flush=True)
        if not ok:
            sys.exit(1)


if __name__ == "__main__":
    main()
