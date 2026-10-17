#!/usr/bin/env python
"""Download the Perez et al. 2022 lupus PBMC h5ad from CELLxGENE (~12 GB).

* Asks the CELLxGENE Curation API for the current file URL and size, falling
  back to the URL pinned in config.yaml if the API is unreachable.
* Checks free disk space first.
* Resumes a partial download (HTTP Range) instead of starting over.
* Verifies the final file size against the size the API reports.

Usage:  python scripts/00_download.py [--config config/config.yaml]
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sle_pbmc.config import load_config  # noqa: E402

API = "https://api.cellxgene.cziscience.com/curation/v1/collections/{cid}"
CHUNK = 8 * 1024 * 1024


def resolve_asset(cfg: dict) -> tuple[str, int | None]:
    cid, did = cfg["data"]["collection_id"], cfg["data"]["dataset_id"]
    try:
        r = requests.get(API.format(cid=cid), timeout=30)
        r.raise_for_status()
        for ds in r.json().get("datasets", []):
            if ds.get("dataset_id") == did:
                for a in ds.get("assets", []):
                    if a.get("filetype") == "H5AD":
                        return a["url"], int(a["filesize"])
        print(f"[warn] dataset {did} not found in collection; using fallback URL")
    except Exception as e:  # network/API problems -> fallback
        print(f"[warn] CELLxGENE API unavailable ({e}); using fallback URL")
    return cfg["data"]["fallback_url"], None


def download(url: str, dest: Path, expected: int | None) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")

    if dest.exists() and (expected is None or dest.stat().st_size == expected):
        print(f"Already downloaded: {dest} ({dest.stat().st_size / 1e9:.2f} GB)")
        return

    free = shutil.disk_usage(dest.parent).free
    need = (expected or 12.3e9) - (part.stat().st_size if part.exists() else 0)
    if free < need * 1.05:
        sys.exit(f"Not enough disk space: need ~{need / 1e9:.1f} GB, have {free / 1e9:.1f} GB.")

    done = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={done}-"} if done else {}
    with requests.get(url, headers=headers, stream=True, timeout=60) as r:
        if done and r.status_code != 206:
            print("[warn] server ignored resume request; restarting download")
            done = 0
            part.unlink(missing_ok=True)
        r.raise_for_status()
        total = expected or (int(r.headers.get("Content-Length", 0)) + done) or None
        mode = "ab" if done else "wb"
        with open(part, mode) as fh, tqdm(
            total=total, initial=done, unit="B", unit_scale=True, desc=dest.name
        ) as bar:
            for chunk in r.iter_content(CHUNK):
                fh.write(chunk)
                bar.update(len(chunk))

    size = part.stat().st_size
    if expected is not None and size != expected:
        sys.exit(
            f"Size mismatch: got {size} bytes, expected {expected}. "
            "Re-run to resume, or delete the .part file to start over."
        )
    part.replace(dest)  # replace() also overwrites an old copy on Windows
    print(f"Saved {dest} ({size / 1e9:.2f} GB)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=None)
    args = ap.parse_args()
    cfg = load_config(args.config)
    url, size = resolve_asset(cfg)
    print(f"Source: {url}" + (f"  ({size / 1e9:.2f} GB)" if size else ""))
    download(url, Path(cfg["data"]["raw_h5ad"]), size)


if __name__ == "__main__":
    main()
