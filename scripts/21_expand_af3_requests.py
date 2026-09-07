#!/usr/bin/env python3
"""
21_expand_af3_requests.py — AlphaFold-Server (AF3) requests, Batch 2 and Batch 3.

PHASE 4 extension, authorised by the mentor on 2026-08-08.

Batch 1 (scripts/04) covered only the four fragments with experimental Ca-bound structures.
That calibration passed, so this script asks for the two things it unlocked:

  Batch 2 — STRUCTURAL COVERAGE FOR ALL 43 cbEGF DOMAINS.
    Only 524 of the 2,897 variants that sit inside a cbEGF domain currently have any
    structure. Coverage is extended with tandem constructs of three consecutive cbEGF
    domains, stepping by two, so every domain except cbEGF1 appears at least once WITH its
    N-terminal neighbour. That matters: cbEGF calcium affinity is modulated by N-terminal
    domain linkage, so an isolated domain is the wrong construct. Whatever lies between two
    consecutive cbEGF domains (a TB or hybrid domain) is included, because the sequence is
    taken as one contiguous span — the same construct logic as the calibrated 2W86 and 1UZJ
    fragments.

  Batch 3 — MUTANT COFOLDING.
    Wild-type models cannot tell us whether a substitution destroys the metal site. These
    jobs ask AF3 to fold the mutant sequence with the same Ca2+ stoichiometry: if the site is
    destroyed, the ion should move, lose ligands, or fail to be placed. Every mutant is built
    on a CALIBRATED construct, so each has both an experimental reference and an AF3 wild-type
    counterpart. N2144S is included specifically because its calcium affinity has been
    measured experimentally (reduced ~5-fold in isolated cbEGF32, ~9-fold in the TB6-cbEGF32
    pair), which gives the whole analysis one external anchor.

Nothing here is fabricated: sequences are sliced from the cached UniProt P35555 record, the
cbEGF index is recomputed from UniProt features, and every mutant's wild-type residue is
verified against the reference sequence before the job is written.

Writes:
  handoff/structures/inbox/batch2_domains/     job JSONs + batch file + fasta + manifest
  handoff/structures/inbox/batch3_mutants/     same, for the mutant jobs
  handoff/structures/REQUESTS.md               appended (batch 1 section left intact)
  manifest/af3_requests_batch2_3.md            provenance
  logs/21_af3_requests_<stamp>.log
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/structures/inbox"
B2 = INBOX / "batch2_domains"
B3 = INBOX / "batch3_mutants"
DEST = "structures/alphafold/"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

# Batch 2 geometry: three consecutive cbEGF domains per construct, stepping by two.
DOMAINS_PER_CONSTRUCT = 3
STEP = 2

# Batch 3. Every entry is (construct id, range, mutation, why it is in the set).
# Ranges are IDENTICAL to the Batch 1 calibration constructs so the mutant model superposes
# residue-for-residue on both the AF3 wild type and the experimental structure.
CALIBRATED = {
    "cbEGF9-hyb2-cbEGF10": (807, 951, 2, "2W86"),
    "cbEGF12-13": (1069, 1154, 2, "1LMJ"),
    "cbEGF22-TB4-cbEGF23": (1486, 1647, 2, "1UZJ"),
    "cbEGF32-33": (2127, 2205, 2, "1EMN"),
}

MUTANTS = [
    ("cbEGF32-33", 2144, "N", "S",
     "EXPERIMENTAL ANCHOR. Ca affinity measured: ~5x lower in isolated cbEGF32, ~9x in the "
     "TB6-cbEGF32 pair. Pathogenic in ClinVar. The one variant where a modelled answer can be "
     "checked against a measured one."),
    ("cbEGF32-33", 2183, "N", "S",
     "Second literature mutant (the equivalent position in cbEGF33), characterised in the same "
     "domain-context work. Not in our ClinVar set — included as a modelling control."),
    ("cbEGF32-33", 2130, "E", "K",
     "Pathogenic ca_EH ligand, buried (RSA 7.3%). Charge reversal next to the ion."),
    ("cbEGF12-13", 1113, "D", "G",
     "Pathogenic ca_D1 ligand. Removes the carboxylate outright — the maximal ligand deletion."),
    ("cbEGF12-13", 1073, "E", "K",
     "Pathogenic ca_EH ligand, buried (RSA 9.0%), neonatal-severe region."),
    ("cbEGF22-TB4-cbEGF23", 1487, "D", "G",
     "Pathogenic ca_D1 ligand in the TB4-flanked construct — tests a different domain context."),
    ("cbEGF9-hyb2-cbEGF10", 913, "E", "K",
     "Pathogenic ca_EH ligand, the most buried in the set (RSA 2.4%), X-ray reference."),
    ("cbEGF12-13", 1138, "C", "F",
     "POSITIVE CONTROL for the other mechanism: removes disulfide C4. Prediction — the fold is "
     "damaged but the calcium site should survive. If the ion also moves, our two mechanisms "
     "are not as separable as the manuscript claims."),
    ("cbEGF12-13", 1141, "P", "L",
     "NEGATIVE CONTROL: non-critical residue, AlphaMissense 0.23. Prediction — nothing moves. "
     "Bounds how much apparent ion displacement is just AF3 run-to-run noise."),
]


def log_line(lines: list[str], m: str = "") -> None:
    print(m, flush=True)
    lines.append(m)


def load_record() -> tuple[dict, str]:
    """Analyse off the cached UniProt pull; never refetch silently (operating principle 4)."""
    cached = sorted((ROOT / "data/raw").glob("uniprot_P35555_v*.json"))
    if not cached:
        cached = sorted((ROOT / "data/raw").glob("uniprot_P35555_2*.json"))
    if not cached:
        raise SystemExit("STOP: no cached UniProt P35555 record in data/raw/ — run scripts/02 first")
    rec = json.loads(cached[-1].read_text(encoding="utf-8"))
    return rec, cached[-1].name


def cbegf_domains(record: dict) -> list[tuple[int, int, int]]:
    """Return [(cbegf_index, start, end)] for the 43 calcium-binding EGF-like domains."""
    doms = sorted((f for f in record["features"] if f["type"] == "Domain"),
                  key=lambda f: f["location"]["start"]["value"])
    out = []
    for f in doms:
        d = f["description"]
        if d.startswith("EGF-like") and "calcium-binding" in d:
            out.append((len(out) + 1,
                        f["location"]["start"]["value"], f["location"]["end"]["value"]))
    if len(out) != 43:
        raise SystemExit(f"STOP: expected 43 cbEGF domains, UniProt gave {len(out)} — "
                         "the domain model changed; do not proceed on this assumption")
    return out


def features_in(record: dict, start: int, end: int, idx: dict[tuple[int, int], int]) -> str:
    doms = sorted((f for f in record["features"] if f["type"] == "Domain"),
                  key=lambda f: f["location"]["start"]["value"])
    out = []
    for f in doms:
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        if s >= start - 3 and e <= end + 3:
            label = f["description"]
            if (s, e) in idx:
                label += f" [cbEGF{idx[(s, e)]}]"
            out.append(f"{label} ({s}-{e})")
    return "; ".join(out)


def job_json(name: str, seq: str, n_ca: int) -> dict:
    """Documented alphafoldserver dialect v1. Empty modelSeeds: the server picks and reports
    the seed, which must be recorded on return for reproducibility."""
    return {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": 1}},
            {"ion": {"ion": "CA", "count": n_ca}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def write_batch(outdir: Path, jobs: list[dict], rows: list[dict], fasta: list[str],
                batch_name: str) -> None:
    (outdir / "af3_jobs").mkdir(parents=True, exist_ok=True)
    for j in jobs:
        (outdir / "af3_jobs" / f"{j['name']}.json").write_text(
            json.dumps([j], indent=2), encoding="utf-8")
    (outdir / f"{batch_name}.json").write_text(json.dumps(jobs, indent=2), encoding="utf-8")
    (outdir / "constructs.fasta").write_text("\n".join(fasta) + "\n", encoding="utf-8")
    with (outdir / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    lines: list[str] = []
    log = lambda m="": log_line(lines, m)
    log(f"AF3 requests, batches 2 and 3 — {datetime.now(timezone.utc).isoformat()}")

    record, src = load_record()
    seq = record["sequence"]["value"]
    if len(seq) != 2871:
        raise SystemExit(f"STOP: P35555 length {len(seq)} != 2871")
    doms = cbegf_domains(record)
    idx = {(s, e): i for i, s, e in doms}
    log(f"  source {src}; 43 cbEGF domains; sequence MD5 {hashlib.md5(seq.encode()).hexdigest()}\n")

    # ---------------- Batch 2: coverage for all 43 cbEGF domains ----------------
    log("Batch 2 — tandem constructs covering all 43 cbEGF domains")
    b2_jobs, b2_rows, b2_fasta = [], [], []
    starts = list(range(0, len(doms) - DOMAINS_PER_CONSTRUCT + 1, STEP))
    if starts[-1] + DOMAINS_PER_CONSTRUCT < len(doms):
        starts.append(len(doms) - DOMAINS_PER_CONSTRUCT)
    for i0 in starts:
        group = doms[i0:i0 + DOMAINS_PER_CONSTRUCT]
        first, last = group[0], group[-1]
        s, e = first[1], last[2]
        sub = seq[s - 1:e]
        n_ca = len(group)                      # exactly one Ca2+ per cbEGF domain
        cid = f"cbEGF{first[0]}-{last[0]}"
        name = f"FBN1_{cid}_Ca{n_ca}"
        b2_jobs.append(job_json(name, sub, n_ca))
        b2_fasta.append(f">{name} P35555|{s}-{e}|{len(sub)}aa\n"
                        + "\n".join(sub[k:k + 60] for k in range(0, len(sub), 60)))
        b2_rows.append({
            "job_name": name, "construct": cid, "uniprot_range": f"{s}-{e}",
            "length_aa": len(sub), "n_ca_ions": n_ca,
            "cbegf_domains": ",".join(str(g[0]) for g in group),
            "domains_covered": features_in(record, s, e, idx),
            "destination": f"{DEST}{name}/",
            "sha256_sequence": hashlib.sha256(sub.encode()).hexdigest()[:16],
        })
        log(f"  {name:28} {s:>5}-{e:<5} {len(sub):>4} aa  {n_ca} Ca  cbEGF "
            f"{','.join(str(g[0]) for g in group)}")
    longest = max(b2_rows, key=lambda r: r["length_aa"])
    log(f"  {len(b2_jobs)} jobs; longest {longest['length_aa']} aa ({longest['job_name']})")
    covered = {int(x) for r in b2_rows for x in r["cbegf_domains"].split(",")}
    if covered != set(range(1, 44)):
        raise SystemExit(f"STOP: tiling misses cbEGF domains {sorted(set(range(1, 44)) - covered)}")
    log("  all 43 cbEGF domains covered; every domain except cbEGF1 appears with its "
        "N-terminal neighbour\n")

    # ---------------- Batch 3: mutant cofolding ----------------
    log("Batch 3 — mutant cofolding on calibrated constructs")
    b3_jobs, b3_rows, b3_fasta = [], [], []
    for cid, pos, wt, mut, why in MUTANTS:
        s, e, n_ca, ref = CALIBRATED[cid]
        if seq[pos - 1] != wt:
            raise SystemExit(f"STOP: wild-type mismatch at {pos}: reference has "
                             f"{seq[pos - 1]}, request says {wt}")
        if not (s <= pos <= e):
            raise SystemExit(f"STOP: position {pos} lies outside construct {cid} ({s}-{e})")
        sub = list(seq[s - 1:e])
        sub[pos - s] = mut
        sub = "".join(sub)
        name = f"FBN1_{cid}_{wt}{pos}{mut}_Ca{n_ca}"
        b3_jobs.append(job_json(name, sub, n_ca))
        b3_fasta.append(f">{name} P35555|{s}-{e}|{wt}{pos}{mut}|ref={ref}\n"
                        + "\n".join(sub[k:k + 60] for k in range(0, len(sub), 60)))
        b3_rows.append({
            "job_name": name, "construct": cid, "uniprot_range": f"{s}-{e}",
            "length_aa": len(sub), "n_ca_ions": n_ca,
            "mutation": f"{wt}{pos}{mut}", "wt_verified_against": "P35555 reference sequence",
            "reference_pdb": ref, "wildtype_counterpart": f"FBN1_{cid}_Ca{n_ca}",
            "destination": f"{DEST}{name}/",
            "sha256_sequence": hashlib.sha256(sub.encode()).hexdigest()[:16],
            "rationale": why,
        })
        log(f"  {name:38} {wt}{pos}{mut}  ref {ref}")
    log(f"  {len(b3_jobs)} mutant jobs\n")

    write_batch(B2, b2_jobs, b2_rows, b2_fasta, "af3_batch2_all")
    write_batch(B3, b3_jobs, b3_rows, b3_fasta, "af3_batch3_all")

    # ---------------- provenance ----------------
    (ROOT / "manifest/af3_requests_batch2_3.md").write_text(
        "# manifest/af3_requests_batch2_3.md — AlphaFold-Server batches 2 and 3\n\n"
        f"Generated by `scripts/21_expand_af3_requests.py` at "
        f"{datetime.now(timezone.utc).isoformat()}.\n\n"
        f"- Source: cached `data/raw/{src}`; sequence length {len(seq)}, "
        f"MD5 `{hashlib.md5(seq.encode()).hexdigest()}`\n"
        f"- Batch 2: {len(b2_jobs)} constructs, {DOMAINS_PER_CONSTRUCT} consecutive cbEGF "
        f"domains each, step {STEP}; one Ca2+ per cbEGF domain\n"
        f"- Batch 3: {len(b3_jobs)} mutant jobs on the four calibrated constructs; every "
        f"wild-type residue verified against the reference sequence\n"
        f"- AF3 dialect: alphafoldserver v1; ion CCD code `CA`\n\n"
        "## Batch 2\n\n| job | range | aa | Ca | cbEGF |\n|---|---|---|---|---|\n"
        + "\n".join(f"| {r['job_name']} | {r['uniprot_range']} | {r['length_aa']} | "
                    f"{r['n_ca_ions']} | {r['cbegf_domains']} |" for r in b2_rows)
        + "\n\n## Batch 3\n\n| job | mutation | construct | reference |\n|---|---|---|---|\n"
        + "\n".join(f"| {r['job_name']} | {r['mutation']} | {r['construct']} | "
                    f"{r['reference_pdb']} |" for r in b3_rows) + "\n",
        encoding="utf-8")

    (ROOT / f"logs/21_af3_requests_{STAMP}.log").write_text("\n".join(lines) + "\n",
                                                            encoding="utf-8")
    log(f"wrote {B2.relative_to(ROOT)}/ and {B3.relative_to(ROOT)}/ "
        f"({len(b2_jobs) + len(b3_jobs)} jobs total)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
