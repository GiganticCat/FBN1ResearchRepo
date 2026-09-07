#!/usr/bin/env python3
"""
04_prepare_af3_requests.py — build the AlphaFold-Server (AF3) job packet for the mentor.

PHASE 4 prep, run early at the mentor's request (they have AF-Server access and generating
models has long wall-clock latency, so this is staged ahead of Phases 1-3).

Scope: the CALIBRATION batch only — the four fragments that have experimental Ca-bound
structures. CLAUDE.md Phase 4 requires calibrating AF3+Ca against known geometry BEFORE
relying on it for un-crystallized domains, and this batch is fully determined now without
needing the variant set. The production batch (which of the 43 cbEGF domains actually carry
variants) is deliberately NOT generated here — it depends on Phase 3.

Grounding (all verified, none assumed):
  * Domain boundaries and the cbEGF1..43 index come from UniProt P35555 features, cached to
    data/raw/. UniProt reports 47 EGF-like (43 calcium-binding) + 9 TB, and its "TB 4" is the
    domain the structural literature calls hybrid-2.
  * Construct ranges match the deposited experimental constructs (SIFTS UniProt ranges), so
    each model superposes residue-for-residue onto its experimental counterpart.
  * Ca2+ stoichiometry was counted from the deposited coordinates, not assumed: exactly
    ONE Ca2+ per cbEGF domain, and ZERO for TB/hybrid domains.
  * AF3 job JSON follows the documented alphafoldserver dialect (v1); "CA" is an allowed ion.

Writes:
  data/raw/uniprot_P35555_<date>.json          cached source record
  handoff/structures/inbox/af3_jobs/*.json     one uploadable job file per construct
  handoff/structures/inbox/af3_batch_all.json  all four jobs in one uploadable file
  handoff/structures/inbox/constructs.fasta    the same sequences as FASTA (reference)
  handoff/structures/inbox/manifest.csv        construct -> range -> Ca count -> destination
  manifest/af3_requests.md                     provenance record
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import date, datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
TODAY = date.today().isoformat()
UNIPROT_ACC = "P35555"

INBOX = ROOT / "handoff/structures/inbox"
JOBS = INBOX / "af3_jobs"
RAW = ROOT / "data/raw"
DEST = "structures/alphafold/"

# Calibration constructs. Ranges are the deposited experimental constructs (SIFTS UniProt
# ranges) so the AF3 model and the experimental structure cover the same residues.
# n_ca was COUNTED from the deposited coordinates (per model / per chain), not assumed.
CONSTRUCTS = [
    {
        "id": "cbEGF9-hyb2-cbEGF10",
        "start": 807, "end": 951, "n_ca": 2,
        "ref_pdb": "2W86", "ref_method": "X-ray 1.80 A",
        "note": "Ca-complete. UniProt calls the middle domain 'TB 4'; literature calls it hybrid-2.",
    },
    {
        "id": "cbEGF12-13",
        "start": 1069, "end": 1154, "n_ca": 2,
        "ref_pdb": "1LMJ", "ref_method": "NMR, 25 models",
        "note": "Neonatal-severe cluster. 2 Ca per model in the deposited ensemble.",
    },
    {
        "id": "cbEGF22-TB4-cbEGF23",
        "start": 1486, "end": 1647, "n_ca": 2,
        "ref_pdb": "1UZJ", "ref_method": "X-ray 2.25 A",
        "note": ("Use 1UZJ, NOT the higher-resolution 1UZK: 1UZJ carries 2 Ca per chain, "
                 "1UZK only 1, and 1UZP/1UZQ are apo."),
    },
    {
        "id": "cbEGF32-33",
        "start": 2127, "end": 2205, "n_ca": 2,
        "ref_pdb": "1EMN", "ref_method": "NMR (1EMO = 22-model ensemble)",
        "note": "Classic Downing 1996 rigid pair with an interdomain Ca site.",
    },
]


def fetch_uniprot() -> dict:
    """Cache the raw UniProt record; analyse off the cache (CLAUDE.md operating principle 4)."""
    RAW.mkdir(parents=True, exist_ok=True)
    cached = RAW / f"uniprot_{UNIPROT_ACC}_{TODAY}.json"
    if cached.is_file():
        print(f"  using cached {cached.relative_to(ROOT)}")
        return json.loads(cached.read_text(encoding="utf-8"))
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ACC}.json",
                     headers={"User-Agent": "fbn1-marfan-pipeline/0.1"}, timeout=60)
    r.raise_for_status()
    cached.write_text(r.text, encoding="utf-8")
    print(f"  fetched and cached {cached.relative_to(ROOT)}")
    return r.json()


def cbegf_index(record: dict) -> dict[tuple[int, int], int]:
    """Number the calcium-binding EGF-like domains 1..43 in sequence order."""
    doms = sorted((f for f in record["features"] if f["type"] == "Domain"),
                  key=lambda f: f["location"]["start"]["value"])
    idx, n = {}, 0
    for f in doms:
        if f["description"].startswith("EGF-like") and "calcium-binding" in f["description"]:
            n += 1
            idx[(f["location"]["start"]["value"], f["location"]["end"]["value"])] = n
    if n != 43:
        raise SystemExit(f"STOP: expected 43 cbEGF domains, UniProt gave {n} — "
                         "domain model changed, do not proceed on this assumption")
    return idx


def domains_in(record: dict, start: int, end: int) -> list[str]:
    doms = sorted((f for f in record["features"] if f["type"] == "Domain"),
                  key=lambda f: f["location"]["start"]["value"])
    idx = cbegf_index(record)
    out = []
    for f in doms:
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        if s >= start - 3 and e <= end + 3:
            label = f["description"]
            if (s, e) in idx:
                label += f" [cbEGF{idx[(s, e)]}]"
            out.append(f"{label} ({s}-{e})")
    return out


def main() -> int:
    print("Preparing AlphaFold-Server (AF3) calibration requests\n")
    record = fetch_uniprot()
    seq = record["sequence"]["value"]
    if len(seq) != 2871:
        raise SystemExit(f"STOP: P35555 length {len(seq)} != 2871")

    JOBS.mkdir(parents=True, exist_ok=True)
    rows, fasta, batch = [], [], []

    for c in CONSTRUCTS:
        s, e = c["start"], c["end"]
        sub = seq[s - 1:e]  # UniProt numbering is 1-based inclusive
        job_name = f"FBN1_{c['id']}_Ca{c['n_ca']}"

        # Documented alphafoldserver dialect v1. Empty modelSeeds lets the server choose;
        # the returned seed is recorded in the output and must be logged for reproducibility.
        job = {
            "name": job_name,
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": sub, "count": 1}},
                {"ion": {"ion": "CA", "count": c["n_ca"]}},
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        }
        (JOBS / f"{job_name}.json").write_text(json.dumps([job], indent=2), encoding="utf-8")
        batch.append(job)

        fasta.append(f">{job_name} P35555|{s}-{e}|{len(sub)}aa|ref={c['ref_pdb']}\n"
                     + "\n".join(sub[i:i + 60] for i in range(0, len(sub), 60)))

        rows.append({
            "job_name": job_name,
            "construct": c["id"],
            "uniprot_range": f"{s}-{e}",
            "length_aa": len(sub),
            "n_ca_ions": c["n_ca"],
            "reference_pdb": c["ref_pdb"],
            "reference_method": c["ref_method"],
            "domains_covered": "; ".join(domains_in(record, s, e)),
            "destination": f"{DEST}{job_name}/",
            "sha256_sequence": hashlib.sha256(sub.encode()).hexdigest()[:16],
            "note": c["note"],
        })
        print(f"  {job_name:34} {s}-{e}  {len(sub):3} aa  {c['n_ca']} Ca  ref {c['ref_pdb']}")

    (INBOX / "af3_batch_all.json").write_text(json.dumps(batch, indent=2), encoding="utf-8")
    (INBOX / "constructs.fasta").write_text("\n".join(fasta) + "\n", encoding="utf-8")

    with (INBOX / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    entry = (record.get("entryAudit") or {}).get("lastSequenceUpdateDate", "n/a")
    (ROOT / "manifest/af3_requests.md").write_text(
        "# manifest/af3_requests.md — AlphaFold-Server calibration batch\n\n"
        f"Generated by `scripts/04_prepare_af3_requests.py` at "
        f"{datetime.now(timezone.utc).isoformat()}.\n\n"
        f"- Source: UniProt {UNIPROT_ACC}, cached `data/raw/uniprot_{UNIPROT_ACC}_{TODAY}.json`\n"
        f"- Sequence length 2871, MD5 `{hashlib.md5(seq.encode()).hexdigest()}`, "
        f"last sequence update {entry}\n"
        f"- Constructs: {len(rows)} (calibration only; production batch awaits Phase 3)\n"
        f"- Ca2+ counts were counted from deposited coordinates, not assumed\n"
        f"- AF3 dialect: alphafoldserver v1; ion CCD code `CA`\n\n"
        "| job | range | aa | Ca | reference |\n|---|---|---|---|---|\n"
        + "\n".join(f"| {r['job_name']} | {r['uniprot_range']} | {r['length_aa']} | "
                    f"{r['n_ca_ions']} | {r['reference_pdb']} |" for r in rows) + "\n",
        encoding="utf-8")

    print(f"\nWrote {len(rows)} job files to {JOBS.relative_to(ROOT)}/")
    print(f"Wrote {(INBOX / 'af3_batch_all.json').relative_to(ROOT)} (all four in one upload)")
    print(f"Wrote {(INBOX / 'manifest.csv').relative_to(ROOT)}, constructs.fasta, "
          f"manifest/af3_requests.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
