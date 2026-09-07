#!/usr/bin/env python3
"""
02_fetch_uniprot.py — PHASE 1. Cache the UniProt P35555 record (JSON + GFF).

This is the source of truth for domain boundaries, disulfide bonds and the reference sequence
that every later phase checks residue identity against. Stored verbatim; parsing happens later.

Also emits a small derived summary (domain table + cbEGF1..43 index + disulfide pairs) for
convenience, clearly separated from the raw record.

Writes:
  data/raw/uniprot_P35555_<release>.json    verbatim REST JSON
  data/raw/uniprot_P35555_<release>.gff     verbatim GFF features
  data/raw/uniprot_P35555_<release>.fasta   reference sequence
  data/interim/uniprot_domains.tsv          derived domain table (cbEGF index, boundaries)
  data/interim/uniprot_disulfides.tsv       derived disulfide pair table
  manifest/uniprot.md                       provenance
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
RAW, INTERIM = ROOT / "data/raw", ROOT / "data/interim"
ACC = "P35555"
EXPECTED_LEN = 2871
UA = {"User-Agent": "fbn1-marfan-pipeline/0.1 (Phase 1 UniProt)"}
BASE = "https://rest.uniprot.org/uniprotkb"


def fetch(fmt: str) -> str:
    r = requests.get(f"{BASE}/{ACC}.{fmt}", headers=UA, timeout=120)
    r.raise_for_status()
    return r.text


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    INTERIM.mkdir(parents=True, exist_ok=True)
    print(f"UniProt {ACC} pull\n")

    js_text = fetch("json")
    record = json.loads(js_text)

    seq = record["sequence"]["value"]
    if record["primaryAccession"] != ACC or len(seq) != EXPECTED_LEN:
        raise SystemExit(f"STOP: got {record['primaryAccession']} length {len(seq)}, "
                         f"expected {ACC} length {EXPECTED_LEN}")

    # Release identity: use the entry's own version + last sequence update, not today's date,
    # so a re-run on an unchanged entry reuses the same filenames.
    audit = record.get("entryAudit", {})
    ver = audit.get("entryVersion", "NA")
    seq_date = audit.get("lastSequenceUpdateDate", "NA")
    tag = f"v{ver}"

    paths = {
        "json": RAW / f"uniprot_{ACC}_{tag}.json",
        "gff": RAW / f"uniprot_{ACC}_{tag}.gff",
        "fasta": RAW / f"uniprot_{ACC}_{tag}.fasta",
    }
    paths["json"].write_text(js_text, encoding="utf-8")
    for fmt in ("gff", "fasta"):
        paths[fmt].write_text(fetch(fmt), encoding="utf-8")
    for k, p in paths.items():
        print(f"  wrote {p.relative_to(ROOT)} ({p.stat().st_size:,} bytes)")

    feats = record["features"]
    doms = sorted((f for f in feats if f["type"] == "Domain"),
                  key=lambda f: f["location"]["start"]["value"])
    disulf = [f for f in feats if f["type"] == "Disulfide bond"]

    # cbEGF index 1..43 in sequence order — the numbering the FBN1 literature uses.
    rows, cb = [], 0
    for f in doms:
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        desc = f["description"]
        is_cb = desc.startswith("EGF-like") and "calcium-binding" in desc
        if is_cb:
            cb += 1
        kind = ("cbEGF" if is_cb else
                "EGF" if desc.startswith("EGF-like") else
                "TB" if desc.startswith("TB") else "other")
        rows.append({
            "uniprot_description": desc, "kind": kind,
            "cbegf_index": cb if is_cb else "",
            "start": s, "end": e, "length": e - s + 1,
            "n_cys": seq[s - 1:e].count("C"),
            "sequence": seq[s - 1:e],
        })
    if cb != 43:
        raise SystemExit(f"STOP: found {cb} calcium-binding EGF-like domains, expected 43 — "
                         "UniProt's domain model changed; do not proceed silently")

    n_egf = sum(1 for r in rows if r["kind"] in ("EGF", "cbEGF"))
    n_tb = sum(1 for r in rows if r["kind"] == "TB")
    print(f"\n  domains: {len(rows)}  (EGF-like {n_egf}, of which cbEGF {cb}; TB {n_tb})")
    print(f"  disulfide bonds: {len(disulf)}")

    with (INTERIM / "uniprot_domains.tsv").open("w", encoding="utf-8") as fh:
        cols = ["uniprot_description", "kind", "cbegf_index", "start", "end", "length",
                "n_cys", "sequence"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # Disulfide pairs, annotated with which domain each partner falls in.
    def domain_of(pos: int) -> str:
        for r in rows:
            if r["start"] <= pos <= r["end"]:
                return (f"cbEGF{r['cbegf_index']}" if r["kind"] == "cbEGF"
                        else r["uniprot_description"])
        return "outside annotated domains"

    ds_rows = []
    for f in disulf:
        loc = f["location"]
        a, b = loc["start"]["value"], loc["end"]["value"]
        if not isinstance(a, int) or not isinstance(b, int):
            continue
        ds_rows.append({
            "cys1": a, "cys2": b,
            "res1": seq[a - 1], "res2": seq[b - 1],
            "domain1": domain_of(a), "domain2": domain_of(b),
            "intra_domain": domain_of(a) == domain_of(b),
        })
    bad = [d for d in ds_rows if d["res1"] != "C" or d["res2"] != "C"]
    if bad:
        raise SystemExit(f"STOP: {len(bad)} annotated disulfide endpoints are not Cys: {bad[:5]}")
    print(f"  all {len(ds_rows)} disulfide endpoints verified as Cys")

    with (INTERIM / "uniprot_disulfides.tsv").open("w", encoding="utf-8") as fh:
        cols = ["cys1", "cys2", "res1", "res2", "domain1", "domain2", "intra_domain"]
        fh.write("\t".join(cols) + "\n")
        for d in ds_rows:
            fh.write("\t".join(str(d[c]) for c in cols) + "\n")

    binding = [f for f in feats if f["type"] == "Binding site"]
    ftypes = {}
    for f in feats:
        ftypes[f["type"]] = ftypes.get(f["type"], 0) + 1

    (ROOT / "manifest/uniprot.md").write_text(
        f"# manifest/uniprot.md — UniProt {ACC}\n\n"
        f"Generated by `scripts/02_fetch_uniprot.py` at {datetime.now(timezone.utc).isoformat()}.\n\n"
        f"- Endpoint: `{BASE}/{ACC}.{{json,gff,fasta}}`\n"
        f"- Entry version: **{ver}**, last sequence update {seq_date}\n"
        f"- Sequence length {len(seq)}, MD5 `{hashlib.md5(seq.encode()).hexdigest()}`, "
        f"SHA256 `{hashlib.sha256(seq.encode()).hexdigest()}`\n"
        + "".join(f"- `data/raw/{p.name}` — {p.stat().st_size:,} bytes, "
                  f"sha256 `{hashlib.sha256(p.read_bytes()).hexdigest()[:32]}…`\n"
                  for p in paths.values())
        + f"\n## Feature counts\n\n"
        + "\n".join(f"- {k}: {v}" for k, v in sorted(ftypes.items(), key=lambda x: -x[1]))
        + f"\n\n## Domain model\n\n"
        f"- EGF-like domains: **{n_egf}** (calcium-binding **{cb}**)\n"
        f"- TB domains as annotated by UniProt: **{n_tb}**\n"
        f"- Disulfide bonds: **{len(ds_rows)}**, every endpoint verified as Cys\n"
        f"- Binding-site features: **{len(binding)}**\n\n"
        "**Caveat carried into Phase 3:** UniProt annotates no `Binding site` features for this\n"
        "entry, so there are zero explicit calcium-site annotations to validate a computed\n"
        "coordinating-residue set against. The Phase 3 gate must use the consensus-motif\n"
        "derivation checked against disulfide annotations, the experimental Ca sites in\n"
        "`resources/reference/structures.md`, and the literature consensus instead.\n\n"
        "**Caveat on TB counts:** UniProt's TB set includes the 2 hybrid domains (its 'TB 4' is\n"
        "the domain the structural literature calls hybrid-2), so this count is not directly\n"
        "comparable with the literature's '7 TB domains'.\n",
        encoding="utf-8")

    print(f"  Binding site features: {len(binding)} "
          f"({'none — see manifest caveat' if not binding else ''})")
    print("\n  wrote data/interim/uniprot_domains.tsv, data/interim/uniprot_disulfides.tsv")
    print("  wrote manifest/uniprot.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
