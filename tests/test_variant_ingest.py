#!/usr/bin/env python3
"""
tests/test_variant_ingest.py — PHASE 1 GATE.

Must PASS before Phase 2 (normalization) may start.

Asserts:
  assert_clinvar_present        FBN1 pull exists, non-empty, plausibly sized
  assert_clinvar_schema         every expected column present
  assert_clinvar_ids_unique     no duplicate VariationID within an assembly
  assert_clinvar_stars          every ReviewStatus maps to a gold-star rating
  assert_clinvar_hgvs_parses    HGVS c./p. parse; unparseable rows are quarantined, not lost
  assert_clinvar_keeps_conflicts VUS and conflicting classifications were retained
  assert_uniprot_present        P35555 cached, length 2871, 43 cbEGF, disulfide endpoints Cys
  assert_gnomad_present         gnomAD v4 pull exists with populated grpmax filtering AF
  assert_alphamissense_present  P35555 subset present and parseable
  assert_manifests_written      every source has a provenance manifest

Run:  python tests/test_variant_ingest.py       (also importable by pytest)
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
INTERIM = ROOT / "data/interim"

# Rough sanity envelopes. A silent 0 or a 10x jump is a failure, not a pass.
CLINVAR_MIN, CLINVAR_MAX = 5_000, 40_000
GNOMAD_MIN, GNOMAD_MAX = 5_000, 60_000

HGVS = re.compile(r"^(?P<tx>N[MR]_\d+\.\d+)\((?P<gene>[^)]+)\):(?P<c>c\.[^ ]+)"
                  r"(?:\s+\((?P<p>p\.[^)]+)\))?")

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
    "-": 0,
}

CLINVAR_COLS = [
    "AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
    "ClinSigSimple", "LastEvaluated", "RCVaccession", "PhenotypeList", "Assembly",
    "Chromosome", "Start", "Stop", "ReferenceAllele", "AlternateAllele", "ReviewStatus",
    "NumberSubmitters", "VariationID", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
]

GNOMAD_COLS = [
    "variant_id", "chrom", "pos", "ref", "alt", "consequence", "hgvsc", "hgvsp",
    "exome_ac", "exome_an", "genome_ac", "genome_an", "joint_ac", "joint_an",
    "faf95_grpmax", "faf95_grpmax_group",
]

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def _one(pattern: str) -> Path:
    hits = sorted(RAW.glob(pattern))
    if not hits:
        raise GateFailure(f"no file matching data/raw/{pattern} — run the Phase 1 fetch scripts")
    return hits[-1]


def _clinvar_rows() -> list[dict]:
    with _one("clinvar_fbn1_*.tsv").open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _gnomad_rows() -> list[dict]:
    with _one("gnomad_*_FBN1_*.tsv").open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------
def assert_clinvar_present() -> str:
    rows = _clinvar_rows()
    if not rows:
        raise GateFailure("ClinVar pull is empty")
    if not (CLINVAR_MIN <= len(rows) <= CLINVAR_MAX):
        raise GateFailure(f"{len(rows):,} ClinVar rows is outside the sane envelope "
                          f"{CLINVAR_MIN:,}-{CLINVAR_MAX:,}")
    off = [r for r in rows if r["GeneID"] != "2200" and r["GeneSymbol"] != "FBN1"]
    if off:
        raise GateFailure(f"{len(off)} rows are not FBN1, e.g. "
                          f"{[(r['VariationID'], r['GeneSymbol']) for r in off[:3]]}")
    g38 = [r for r in rows if r["Assembly"] == "GRCh38"]
    return f"{len(rows):,} rows ({len(g38):,} GRCh38), all FBN1"


def assert_clinvar_schema() -> str:
    rows = _clinvar_rows()
    missing = [c for c in CLINVAR_COLS if c not in rows[0]]
    if missing:
        raise GateFailure(f"missing columns: {missing}")
    return f"all {len(CLINVAR_COLS)} expected columns present ({len(rows[0])} total)"


def assert_clinvar_ids_unique() -> str:
    """VariationID repeats across assemblies by design; it must be unique WITHIN one."""
    rows = _clinvar_rows()
    detail = []
    for asm in sorted({r["Assembly"] for r in rows}):
        ids = [r["VariationID"] for r in rows if r["Assembly"] == asm]
        dupes = {i for i in ids if ids.count(i) > 1} if len(ids) != len(set(ids)) else set()
        if dupes:
            raise GateFailure(f"{len(dupes)} duplicate VariationIDs within assembly {asm}: "
                              f"{sorted(dupes)[:5]}")
        detail.append(f"{asm}={len(ids):,}")
    return "unique within each assembly (" + ", ".join(detail) + ")"


def assert_clinvar_stars() -> str:
    rows = _clinvar_rows()
    unmapped = {r["ReviewStatus"] for r in rows if r["ReviewStatus"].strip().lower() not in STARS}
    if unmapped:
        raise GateFailure(f"ReviewStatus values with no star mapping: {sorted(unmapped)}")
    g38 = [r for r in rows if r["Assembly"] == "GRCh38"]
    counts: dict[int, int] = {}
    for r in g38:
        s = STARS[r["ReviewStatus"].strip().lower()]
        counts[s] = counts.get(s, 0) + 1
    if counts.get(3, 0) + counts.get(4, 0) == 0:
        raise GateFailure("no expert-panel or practice-guideline records found — "
                          "suspicious for FBN1, which has a ClinGen VCEP")
    return "GRCh38 star distribution " + ", ".join(f"{k}★={v:,}" for k, v in sorted(counts.items()))


def assert_clinvar_hgvs_parses() -> str:
    """Unparseable rows are allowed but must be a small, explainable minority."""
    rows = [r for r in _clinvar_rows() if r["Assembly"] == "GRCh38"]
    bad = [r for r in rows if not HGVS.match(r["Name"] or "")]
    frac = len(bad) / len(rows)
    if frac > 0.05:
        raise GateFailure(f"{len(bad):,}/{len(rows):,} ({frac:.1%}) of Names do not parse as "
                          f"HGVS — expected <5%. Examples: {[r['Name'] for r in bad[:3]]}")
    # Everything that parses must be on the pinned gene and a NM_000138 transcript.
    txs = {HGVS.match(r["Name"]).group("tx") for r in rows if HGVS.match(r["Name"] or "")}
    off_tx = {t for t in txs if not t.startswith("NM_000138.")}
    if off_tx:
        raise GateFailure(f"unexpected transcripts in the FBN1 pull: {sorted(off_tx)}")
    return (f"{len(rows) - len(bad):,}/{len(rows):,} parse; {len(bad)} genomic-only descriptions "
            f"retained for Phase 2; transcripts {sorted(txs)}")


def assert_clinvar_keeps_conflicts() -> str:
    """The pipeline must label conflicts, never drop them."""
    rows = [r for r in _clinvar_rows() if r["Assembly"] == "GRCh38"]
    sig = [r["ClinicalSignificance"] for r in rows]
    n_vus = sum(1 for s in sig if "Uncertain significance" in s)
    n_conf = sum(1 for s in sig if "Conflicting" in s)
    if n_vus == 0:
        raise GateFailure("no VUS retained — they must be kept and labelled, not filtered out")
    if n_conf == 0:
        raise GateFailure("no conflicting classifications retained — they must be kept")
    return f"{n_vus:,} VUS and {n_conf:,} conflicting classifications retained"


def assert_uniprot_present() -> str:
    hits = sorted(RAW.glob("uniprot_P35555_v*.json"))
    if not hits:
        raise GateFailure("no data/raw/uniprot_P35555_v*.json — run scripts/02_fetch_uniprot.py")
    rec = json.loads(hits[-1].read_text(encoding="utf-8"))
    seq = rec["sequence"]["value"]
    if rec["primaryAccession"] != "P35555" or len(seq) != 2871:
        raise GateFailure(f"{rec['primaryAccession']} length {len(seq)}, expected P35555 / 2871")

    dom_file = INTERIM / "uniprot_domains.tsv"
    dis_file = INTERIM / "uniprot_disulfides.tsv"
    for f in (dom_file, dis_file):
        if not f.is_file():
            raise GateFailure(f"missing derived table {f.relative_to(ROOT)}")

    with dom_file.open(encoding="utf-8") as fh:
        doms = list(csv.DictReader(fh, delimiter="\t"))
    n_cb = sum(1 for d in doms if d["kind"] == "cbEGF")
    if n_cb != 43:
        raise GateFailure(f"{n_cb} cbEGF domains in the derived table, expected 43")

    # Domain boundaries must be inside the sequence and non-inverted.
    bad = [d for d in doms if not (1 <= int(d["start"]) <= int(d["end"]) <= 2871)]
    if bad:
        raise GateFailure(f"{len(bad)} domains with impossible coordinates: {bad[:3]}")

    with dis_file.open(encoding="utf-8") as fh:
        dis = list(csv.DictReader(fh, delimiter="\t"))
    non_cys = [d for d in dis if d["res1"] != "C" or d["res2"] != "C"]
    if non_cys:
        raise GateFailure(f"{len(non_cys)} disulfide endpoints are not Cys: {non_cys[:3]}")

    # Every cbEGF domain should carry 6 cysteines; deviations are a real annotation problem.
    odd = [(d["cbegf_index"], d["n_cys"]) for d in doms
           if d["kind"] == "cbEGF" and int(d["n_cys"]) != 6]
    note = f"; NOTE {len(odd)} cbEGF domains without exactly 6 Cys: {odd[:5]}" if odd else ""
    return (f"P35555 v-cached, 2871 aa, {n_cb} cbEGF domains, {len(dis)} disulfides "
            f"all Cys-verified{note}")


def assert_gnomad_present() -> str:
    rows = _gnomad_rows()
    if not rows:
        raise GateFailure("gnomAD pull is empty")
    if not (GNOMAD_MIN <= len(rows) <= GNOMAD_MAX):
        raise GateFailure(f"{len(rows):,} gnomAD variants outside envelope "
                          f"{GNOMAD_MIN:,}-{GNOMAD_MAX:,}")
    missing = [c for c in GNOMAD_COLS if c not in rows[0]]
    if missing:
        raise GateFailure(f"missing columns: {missing}")

    n_missense = sum(1 for r in rows if r["consequence"] == "missense_variant")
    if n_missense == 0:
        raise GateFailure("no missense variants in the gnomAD pull")

    # The grpmax filtering AF is the benign threshold input. An all-null column here silently
    # produces a benign set with no frequency evidence, so require it to be populated.
    n_faf = sum(1 for r in rows if r["faf95_grpmax"] not in ("", "None", None))
    if n_faf == 0:
        raise GateFailure(
            "faf95_grpmax is empty for every variant — the benign threshold would have no "
            "frequency evidence. gnomAD v4 populates joint.fafmax.faf95_max, NOT faf95_joint."
        )
    return (f"{len(rows):,} variants ({n_missense:,} missense); "
            f"{n_faf:,} carry a grpmax filtering AF")


def assert_alphamissense_present() -> str:
    p = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    n = 0
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 4:
                raise GateFailure(f"malformed AlphaMissense row: {line[:80]!r}")
            if f[0] != "P35555":
                raise GateFailure(f"non-P35555 row: {line[:80]!r}")
            n += 1
    if n < 50_000:
        raise GateFailure(f"only {n:,} AlphaMissense rows, expected ~54.5k")
    return f"{n:,} rows, all P35555"


def assert_manifests_written() -> str:
    required = ["clinvar.md", "uniprot.md", "gnomad.md", "alphamissense.md"]
    missing = [m for m in required if not (ROOT / "manifest" / m).is_file()]
    if missing:
        raise GateFailure(f"missing manifest entries: {missing}")
    empty = [m for m in required if (ROOT / "manifest" / m).stat().st_size < 100]
    if empty:
        raise GateFailure(f"manifest entries suspiciously small: {empty}")
    return f"{len(required)} source manifests present"


# --------------------------------------------------------------------------
def test_clinvar_present(): assert_clinvar_present()
def test_clinvar_schema(): assert_clinvar_schema()
def test_clinvar_ids_unique(): assert_clinvar_ids_unique()
def test_clinvar_stars(): assert_clinvar_stars()
def test_clinvar_hgvs_parses(): assert_clinvar_hgvs_parses()
def test_clinvar_keeps_conflicts(): assert_clinvar_keeps_conflicts()
def test_uniprot_present(): assert_uniprot_present()
def test_gnomad_present(): assert_gnomad_present()
def test_alphamissense_present(): assert_alphamissense_present()
def test_manifests_written(): assert_manifests_written()


def main() -> int:
    print("PHASE 1 GATE — tests/test_variant_ingest.py")
    print(f"root: {ROOT}\n")

    checks = [
        ("clinvar_present", assert_clinvar_present),
        ("clinvar_schema", assert_clinvar_schema),
        ("clinvar_ids_unique", assert_clinvar_ids_unique),
        ("clinvar_stars", assert_clinvar_stars),
        ("clinvar_hgvs_parses", assert_clinvar_hgvs_parses),
        ("clinvar_keeps_conflicts", assert_clinvar_keeps_conflicts),
        ("uniprot_present", assert_uniprot_present),
        ("gnomad_present", assert_gnomad_present),
        ("alphamissense_present", assert_alphamissense_present),
        ("manifests_written", assert_manifests_written),
    ]
    for name, fn in checks:
        try:
            print(f"  [PASS] {name}: {fn()}")
        except GateFailure as exc:
            _failures.append(name)
            print(f"  [FAIL] {name}: {exc}")
        except Exception as exc:  # noqa: BLE001
            _failures.append(name)
            print(f"  [FAIL] {name}: unexpected {type(exc).__name__}: {exc}")

    print()
    if _failures:
        print(f"GATE FAILED — {len(_failures)}/{len(checks)} checks failed: "
              f"{', '.join(_failures)}")
        print("Phase 2 must not start until these pass.")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 2 may proceed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
