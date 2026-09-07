#!/usr/bin/env python3
"""
01_fetch_clinvar.py — PHASE 1. Pull all ClinVar records for FBN1 (Gene ID 2200).

Uses the bulk `variant_summary.txt.gz` rather than E-utilities: one reproducible artifact,
complete coverage, and every field we need (review status, conditions, dates) in one place.
An E-utilities esearch is still issued as a small live cross-check on the row count.

Keeps EVERYTHING — VUS, conflicting classifications, all molecular consequences. Nothing is
filtered on clinical significance here; labelling happens in Phase 2. The only row filter is
GeneSymbol/GeneID == FBN1.

Idempotent: skips the download if the cached archive matches the server's size, and skips the
filter if the extracted TSV is already present and its checksum is recorded.

Writes:
  data/raw/variant_summary_<release>.txt.gz     unmodified upstream archive (gitignored)
  data/raw/clinvar_fbn1_<release>.tsv           FBN1 rows, unmodified columns
  manifest/clinvar.md                           provenance
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GENE_ID, GENE_SYMBOL = "2200", "FBN1"
UA = {"User-Agent": "fbn1-marfan-pipeline/0.1 (Phase 1 ClinVar)"}

# ClinVar review status -> gold stars. Source: ClinVar review status documentation.
# Any value not in this map is a hard STOP rather than a silent 0 — an unmapped status
# would quietly corrupt the Phase 2 star cutoff.
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the single variant": 0,
    # ClinVar writes "-" when a record carries no germline classification at all
    # (ClinicalSignificance is "-" on these rows too). Zero stars, and Phase 2 will
    # exclude them from the labelled sets for lack of a classification — but they are
    # retained here and labelled rather than dropped.
    "-": 0,
}


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def head_release() -> tuple[int, str]:
    """Server size + Last-Modified, used both as the release id and the cache guard."""
    r = requests.head(URL, headers=UA, timeout=60, allow_redirects=True)
    r.raise_for_status()
    size = int(r.headers["Content-Length"])
    lm = parsedate_to_datetime(r.headers["Last-Modified"]).date().isoformat()
    return size, lm


def download(dest: Path, expect_size: int) -> None:
    if dest.is_file() and dest.stat().st_size == expect_size:
        print(f"  cached archive OK ({dest.stat().st_size:,} bytes) — skipping download")
        return
    print(f"  downloading {expect_size / 1e6:.0f} MB -> {dest.relative_to(ROOT)}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(URL, headers=UA, timeout=1800, stream=True) as r:
        r.raise_for_status()
        got = 0
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(1 << 20):
                fh.write(chunk)
                got += len(chunk)
                if got % (50 << 20) < (1 << 20):
                    print(f"    {got / 1e6:.0f} MB", flush=True)
    if tmp.stat().st_size != expect_size:
        tmp.unlink()
        raise SystemExit(f"STOP: downloaded {tmp.stat().st_size} bytes, expected {expect_size}")
    tmp.rename(dest)
    print(f"  done ({dest.stat().st_size:,} bytes)")


def filter_fbn1(archive: Path, out: Path) -> dict:
    """Stream the archive, keeping only FBN1 rows. Columns are passed through unmodified."""
    kept = 0
    total = 0
    header: list[str] = []
    unmapped: set[str] = set()
    stats = {"assemblies": {}, "review_status": {}, "clinsig": {}, "type": {}}

    with gzip.open(archive, "rt", encoding="utf-8", errors="replace") as fh, \
            out.open("w", encoding="utf-8") as fo:
        header_line = fh.readline()
        header = header_line.lstrip("#").rstrip("\n").split("\t")
        # Write the header WITHOUT ClinVar's leading '#'. That marker makes the first column
        # parse as "#AlleleID" in any standard TSV reader. This is the only modification made
        # to the upstream data and it is recorded in manifest/clinvar.md.
        fo.write("\t".join(header) + "\n")
        idx = {c: i for i, c in enumerate(header)}
        i_gid, i_sym = idx["GeneID"], idx["GeneSymbol"]

        for line in fh:
            total += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            if f[i_gid] != GENE_ID and f[i_sym] != GENE_SYMBOL:
                continue
            fo.write(line)
            kept += 1
            rs = f[idx["ReviewStatus"]].strip()
            if rs.lower() not in STARS:
                unmapped.add(rs)
            for key, col in (("assemblies", "Assembly"), ("review_status", "ReviewStatus"),
                             ("clinsig", "ClinicalSignificance"), ("type", "Type")):
                v = f[idx[col]]
                stats[key][v] = stats[key].get(v, 0) + 1

    if unmapped:
        raise SystemExit(
            "STOP: ClinVar ReviewStatus values not in the star map — refusing to guess:\n  "
            + "\n  ".join(sorted(unmapped))
        )
    return {"kept": kept, "scanned": total, "header": header, "stats": stats}


def eutils_crosscheck() -> int | None:
    """Independent count from a different access path; a large mismatch is worth knowing.

    Must use FBN1[gene]. The field "<id>[gene_id]" does NOT exist in the ClinVar E-utilities
    index: NCBI silently rewrites it to a free-text "[All Fields]" match, so "2200[gene_id]"
    returns ~13,500 records including COL2A1 c.2200T>C and other unrelated genes. The
    querytranslation is asserted below so a silent re-translation can never be trusted again.
    """
    try:
        from dotenv import load_dotenv
        load_dotenv(ROOT / ".env", override=False)
        params = {"db": "clinvar", "term": f"{GENE_SYMBOL}[gene]", "retmode": "json", "retmax": "1"}
        if os.environ.get("NCBI_API_KEY"):
            params["api_key"] = os.environ["NCBI_API_KEY"]
        if os.environ.get("NCBI_EMAIL"):
            params["email"] = os.environ["NCBI_EMAIL"]
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                         params=params, headers=UA, timeout=60)
        r.raise_for_status()
        js = r.json()["esearchresult"]
        if "[All Fields]" in js.get("querytranslation", ""):
            print(f"  WARN: esearch degraded to free text ({js['querytranslation']!r}) — "
                  "cross-check discarded as meaningless")
            return None
        return int(js["count"])
    except Exception as exc:  # noqa: BLE001 — cross-check only, never fatal
        print(f"  WARN: E-utilities cross-check failed ({type(exc).__name__}) — not fatal")
        return None


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    print("ClinVar FBN1 pull\n")

    size, release = head_release()
    print(f"  release (Last-Modified): {release}   size: {size:,} bytes")

    archive = RAW / f"variant_summary_{release}.txt.gz"
    out = RAW / f"clinvar_fbn1_{release}.tsv"

    download(archive, size)

    if out.is_file():
        print(f"  {out.relative_to(ROOT)} already present — re-filtering to verify")
    print("  filtering to FBN1 …")
    res = filter_fbn1(archive, out)
    print(f"  kept {res['kept']:,} FBN1 rows from {res['scanned']:,} scanned")

    n_eutils = eutils_crosscheck()
    grch38 = res["stats"]["assemblies"].get("GRCh38", 0)
    if n_eutils:
        print(f"  E-utilities reports {n_eutils:,} ClinVar records for FBN1; "
              f"bulk GRCh38 rows = {grch38:,}")

    digest = sha256(out)
    (ROOT / "manifest/clinvar.md").write_text(
        "# manifest/clinvar.md — ClinVar FBN1 pull\n\n"
        f"Generated by `scripts/01_fetch_clinvar.py` at {datetime.now(timezone.utc).isoformat()}.\n\n"
        f"- Source: `{URL}`\n"
        f"- Release (Last-Modified): **{release}**, {size:,} bytes\n"
        f"- Filter: `GeneID == {GENE_ID}` OR `GeneSymbol == {GENE_SYMBOL}` — no other row filter; "
        "VUS and conflicting classifications retained\n"
        f"- Output: `data/raw/{out.name}` — {res['kept']:,} rows, sha256 `{digest}`\n"
        f"- Columns: {len(res['header'])}, values passed through unmodified. **Sole "
        "modification:** the leading `#` was stripped from the header line so the first column "
        "reads `AlleleID` rather than `#AlleleID`.\n"
        f"- E-utilities cross-check (`{GENE_ID}[gene_id]`): "
        f"{n_eutils:,} records\n" if n_eutils else "- E-utilities cross-check: unavailable\n"
        "\n## Row breakdown\n\n"
        "### Assembly\n\n"
        + "\n".join(f"- {k or '(blank)'}: {v:,}" for k, v in sorted(res["stats"]["assemblies"].items()))
        + "\n\n### Review status (gold stars in parentheses)\n\n"
        + "\n".join(f"- {k or '(blank)'} ({STARS.get(k.lower(), '?')}★): {v:,}"
                    for k, v in sorted(res["stats"]["review_status"].items(), key=lambda x: -x[1]))
        + "\n\n### Variant type\n\n"
        + "\n".join(f"- {k}: {v:,}" for k, v in sorted(res["stats"]["type"].items(), key=lambda x: -x[1]))
        + "\n\n### Clinical significance (raw strings, top 25)\n\n"
        + "\n".join(f"- {k or '(blank)'}: {v:,}"
                    for k, v in sorted(res["stats"]["clinsig"].items(), key=lambda x: -x[1])[:25])
        + "\n",
        encoding="utf-8")

    (RAW / f"clinvar_fbn1_{release}.stats.json").write_text(
        json.dumps({"release": release, "rows": res["kept"], "sha256": digest,
                    "eutils_count": n_eutils, **res["stats"]}, indent=2), encoding="utf-8")

    print(f"\n  wrote {out.relative_to(ROOT)} (sha256 {digest[:16]}…)")
    print("  wrote manifest/clinvar.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
