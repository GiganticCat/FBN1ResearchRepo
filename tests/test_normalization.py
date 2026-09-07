#!/usr/bin/env python3
"""
tests/test_normalization.py — PHASE 2 GATE.

Must PASS before Phase 3 (domain & site annotation) may start.

Asserts:
  assert_outputs_exist          all four output tables + the master table were written
  assert_wt_identity            every kept variant's WT residue == P35555 at that position
  assert_positions_in_range     1 <= position <= 2871, integer, no NaNs
  assert_all_missense           nothing non-missense leaked into the normalized set
  assert_counts_reconcile       normalized + excluded + mismatch + unparsed == ClinVar GRCh38 in
  assert_sets_disjoint          pathogenic ∩ benign == ∅ in both tiers
  assert_no_contradictory_labels no variant is both P/LP and gnomAD-common without being flagged
  assert_no_silent_nans         required columns are fully populated; optional ones accounted for
  assert_alphamissense_identity the AM join used (position, WT, ALT), not position alone
  assert_criteria_recorded      inclusion_criteria.md exists and matches the thresholds used

Run:  python tests/test_normalization.py       (also importable by pytest)
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW, INTERIM, PROC = ROOT / "data/raw", ROOT / "data/interim", ROOT / "data/processed"

FAF_BENIGN = 0.001
AA = set("ACDEFGHIKLMNPQRSTVWY")

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def _seq() -> str:
    hits = sorted(RAW.glob("uniprot_P35555_v*.json"))
    if not hits:
        raise GateFailure("no cached UniProt record")
    return json.loads(hits[-1].read_text(encoding="utf-8"))["sequence"]["value"]


def _read(path: Path) -> list[dict]:
    if not path.is_file():
        raise GateFailure(f"missing {path.relative_to(ROOT)}")
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        return []
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def _norm() -> list[dict]:
    return _read(INTERIM / "variants_normalized.tsv")


def fnum(v):
    return float(v) if v not in ("", "None", None) else None


# --------------------------------------------------------------------------
def assert_outputs_exist() -> str:
    paths = [
        INTERIM / "variants_normalized.tsv",
        INTERIM / "excluded_from_structural.tsv",
        INTERIM / "wt_mismatch.tsv",
        INTERIM / "unparsed.tsv",
        PROC / "variant_master.tsv",
    ]
    missing = [p.relative_to(ROOT) for p in paths if not p.is_file()]
    if missing:
        raise GateFailure(f"missing outputs: {missing}")
    n = len(_norm())
    if n == 0:
        raise GateFailure("variants_normalized.tsv is empty")
    if len(_read(PROC / "variant_master.tsv")) != n:
        raise GateFailure("variant_master.tsv row count differs from variants_normalized.tsv")
    return f"5 tables present; normalized set = {n:,} variants"


def assert_wt_identity() -> str:
    """The single most important check in the pipeline."""
    seq = _seq()
    rows = _norm()
    bad = []
    for r in rows:
        p = int(r["position"])
        if seq[p - 1] != r["wt_aa"]:
            bad.append((r["variation_id"], r["hgvs_p"], r["wt_aa"], seq[p - 1], p))
    if bad:
        raise GateFailure(
            f"{len(bad)} variants whose WT residue disagrees with P35555 "
            f"(id, hgvs_p, claimed, reference, pos): {bad[:5]}"
        )
    # The recorded reference_aa column must agree too — catches a stale/derived column.
    stale = [r["variation_id"] for r in rows if r.get("reference_aa") != r["wt_aa"]]
    if stale:
        raise GateFailure(f"{len(stale)} rows where reference_aa != wt_aa: {stale[:5]}")
    return f"all {len(rows):,} variants: WT residue matches P35555"


def assert_positions_in_range() -> str:
    rows = _norm()
    bad = []
    for r in rows:
        try:
            p = int(r["position"])
        except (TypeError, ValueError):
            bad.append((r["variation_id"], r["position"], "not an integer"))
            continue
        if not (1 <= p <= 2871):
            bad.append((r["variation_id"], p, "outside 1-2871"))
    if bad:
        raise GateFailure(f"{len(bad)} bad positions: {bad[:5]}")
    ps = [int(r["position"]) for r in rows]
    return f"all positions integer and within 1-2871 (observed {min(ps)}-{max(ps)})"


def assert_all_missense() -> str:
    rows = _norm()
    bad = []
    for r in rows:
        wt, mut = r["wt_aa"], r["mut_aa"]
        if wt not in AA or mut not in AA:
            bad.append((r["variation_id"], wt, mut, "non-standard residue"))
        elif wt == mut:
            bad.append((r["variation_id"], wt, mut, "synonymous leaked in"))
    if bad:
        raise GateFailure(f"{len(bad)} non-missense rows in the normalized set: {bad[:5]}")
    return f"all {len(rows):,} rows are genuine missense (standard AA, wt != mut)"


def assert_counts_reconcile() -> str:
    cv = sorted(RAW.glob("clinvar_fbn1_*.tsv"))[-1]
    with cv.open(encoding="utf-8") as fh:
        n_in = sum(1 for r in csv.DictReader(fh, delimiter="\t") if r["Assembly"] == "GRCh38")
    parts = {
        "normalized": len(_norm()),
        "excluded": len(_read(INTERIM / "excluded_from_structural.tsv")),
        "wt_mismatch": len(_read(INTERIM / "wt_mismatch.tsv")),
        "unparsed": len(_read(INTERIM / "unparsed.tsv")),
    }
    total = sum(parts.values())
    if total != n_in:
        raise GateFailure(f"{total:,} rows out != {n_in:,} ClinVar GRCh38 rows in — "
                          f"{parts}. Rows were dropped or duplicated.")
    return (f"{n_in:,} in == {total:,} out ("
            + ", ".join(f"{k}={v:,}" for k, v in parts.items()) + ")")


def assert_sets_disjoint() -> str:
    rows = _norm()
    detail = []
    for tier in ("primary", "sensitivity"):
        col = f"set_{tier}"
        if col not in rows[0]:
            raise GateFailure(f"missing column {col}")
        path = {r["variation_id"] for r in rows if r[col] == "pathogenic"}
        ben = {r["variation_id"] for r in rows if r[col] == "benign"}
        overlap = path & ben
        if overlap:
            raise GateFailure(f"{tier}: {len(overlap)} variants in BOTH pathogenic and benign: "
                              f"{sorted(overlap)[:5]}")
        detail.append(f"{tier} P={len(path):,} B={len(ben):,}")
    return "disjoint in both tiers (" + "; ".join(detail) + ")"


def assert_no_contradictory_labels() -> str:
    """A ClinVar P/LP call that is gnomAD-common must be surfaced, never silently kept as P."""
    rows = _norm()
    hidden = []
    for r in rows:
        faf = fnum(r.get("faf95_grpmax"))
        if faf is not None and faf > FAF_BENIGN:
            for tier in ("primary", "sensitivity"):
                if r[f"set_{tier}"] == "pathogenic":
                    hidden.append((r["variation_id"], r["hgvs_p"], faf, tier))
    if hidden:
        raise GateFailure(
            f"{len(hidden)} variants labelled pathogenic despite grpmax FAF > {FAF_BENIGN} "
            f"without a contradiction flag: {hidden[:5]}"
        )
    flagged = sum(1 for r in rows if "contradiction" in (r["set_primary"], r["set_sensitivity"]))
    return f"no hidden contradictions ({flagged} explicitly flagged)"


def assert_no_silent_nans() -> str:
    rows = _norm()
    required = ["variation_id", "position", "wt_aa", "mut_aa", "reference_aa",
                "clinvar_class", "review_status", "stars", "hgvs_c", "hgvs_p"]
    empties = {c: sum(1 for r in rows if r.get(c) in ("", None)) for c in required}
    bad = {c: n for c, n in empties.items() if n}
    if bad:
        raise GateFailure(f"required columns with empty values: {bad}")
    # Optional columns are allowed to be empty, but the counts must be reported, not hidden.
    optional = ["gnomad_variant_id", "faf95_grpmax", "am_pathogenicity", "am_class"]
    filled = {c: sum(1 for r in rows if r.get(c) not in ("", None)) for c in optional}
    return (f"{len(required)} required columns fully populated; optional coverage: "
            + ", ".join(f"{c}={n:,}/{len(rows):,}" for c, n in filled.items()))


def assert_alphamissense_identity() -> str:
    """Re-derive the join independently and confirm identity matching was used."""
    seq = _seq()
    am_path = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
    am: dict[tuple[int, str, str], str] = {}
    conflicts = 0
    with am_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, var, score, cls = line.rstrip("\n").split("\t")
            w, p, a = var[0], int(var[1:-1]), var[-1]
            if seq[p - 1] != w:
                conflicts += 1
                continue
            am[(p, w, a)] = score
    if conflicts == 0:
        raise GateFailure("no AlphaMissense rows were rejected on WT identity — expected the "
                          "position-472 C-allele rows to be filtered; the join may be using "
                          "position alone")

    rows = _norm()
    wrong = []
    for r in rows:
        key = (int(r["position"]), r["wt_aa"], r["mut_aa"])
        expect = am.get(key, "")
        got = r.get("am_pathogenicity", "")
        if expect and got and abs(float(expect) - float(got)) > 1e-9:
            wrong.append((r["variation_id"], expect, got))
        if expect and not got:
            wrong.append((r["variation_id"], expect, "MISSING"))
    if wrong:
        raise GateFailure(f"{len(wrong)} AlphaMissense score mismatches: {wrong[:5]}")

    # No variant may carry a score for a position where AM's WT disagrees with the reference.
    at472 = [r for r in rows if int(r["position"]) == 472 and r.get("am_pathogenicity")]
    bad472 = [r["variation_id"] for r in at472 if r["wt_aa"] != seq[471]]
    if bad472:
        raise GateFailure(f"position-472 rows scored against the wrong WT: {bad472[:5]}")
    n = sum(1 for r in rows if r.get("am_pathogenicity"))
    return (f"{conflicts} AM rows rejected on WT identity; {n:,}/{len(rows):,} variants scored, "
            "all values re-derived independently")


def assert_criteria_recorded() -> str:
    p = ROOT / "resources/reference/inclusion_criteria.md"
    if not p.is_file():
        raise GateFailure("resources/reference/inclusion_criteria.md missing — thresholds must "
                          "be recorded before the sets are built")
    text = p.read_text(encoding="utf-8")
    needed = ["NM_000138.5", "P35555", "0.001", "2 stars", "conflicting"]
    missing = [n for n in needed if n not in text]
    if missing:
        raise GateFailure(f"inclusion_criteria.md does not document: {missing}")
    return "thresholds recorded and consistent with the values used"


# --------------------------------------------------------------------------
def test_outputs_exist(): assert_outputs_exist()
def test_wt_identity(): assert_wt_identity()
def test_positions_in_range(): assert_positions_in_range()
def test_all_missense(): assert_all_missense()
def test_counts_reconcile(): assert_counts_reconcile()
def test_sets_disjoint(): assert_sets_disjoint()
def test_no_contradictory_labels(): assert_no_contradictory_labels()
def test_no_silent_nans(): assert_no_silent_nans()
def test_alphamissense_identity(): assert_alphamissense_identity()
def test_criteria_recorded(): assert_criteria_recorded()


def main() -> int:
    print("PHASE 2 GATE — tests/test_normalization.py")
    print(f"root: {ROOT}\n")
    checks = [
        ("outputs_exist", assert_outputs_exist),
        ("wt_identity", assert_wt_identity),
        ("positions_in_range", assert_positions_in_range),
        ("all_missense", assert_all_missense),
        ("counts_reconcile", assert_counts_reconcile),
        ("sets_disjoint", assert_sets_disjoint),
        ("no_contradictory_labels", assert_no_contradictory_labels),
        ("no_silent_nans", assert_no_silent_nans),
        ("alphamissense_identity", assert_alphamissense_identity),
        ("criteria_recorded", assert_criteria_recorded),
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
        print("Phase 3 must not start until these pass.")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 3 may proceed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
