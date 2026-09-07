#!/usr/bin/env python3
"""
tests/test_annotation.py — PHASE 3 GATE.

Must PASS before Phase 4 (structure sourcing) may start.

Asserts:
  assert_outputs_exist          annotated table + derived site tables written
  assert_counts_reconcile       annotation preserved every Phase 2 row, no dupes
  assert_cysteine_annotations   every position annotated as a cysteine site really is Cys
  assert_cbegf_cysteine_count   exactly 258 cbEGF cysteines (the VCEP's stated figure)
  assert_disulfide_pattern      every cbEGF disulfide is C1-C3, C2-C4 or C5-C6
  assert_consensus_motif        derived consensus positions satisfy the VCEP motif letters
  assert_consensus_vs_structure derived Ca sites reproduce the ligands OBSERVED in 2W86/1LMJ/
                                1UZJ/1EMN — the replacement for the unavailable UniProt check
  assert_crosstab_sane          class x feature cross-tab is non-empty and internally consistent
  assert_no_silent_nans         annotation columns populated; domain coverage accounted for

NOTE on the CLAUDE.md Phase 3 gate spec: it asks that the coordinating-residue set be validated
against "UniProt's explicit Ca-binding annotations". UniProt P35555 carries ZERO `Binding site`
features, so that check is impossible. It is replaced here by a stronger one — agreement with
calcium ligands measured directly from experimental coordinates.

Run:  python tests/test_annotation.py
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW, INTERIM, PROC = ROOT / "data/raw", ROOT / "data/interim", ROOT / "data/processed"

# Calcium ligands measured from the deposited structures (UniProt numbering).
# sidechain-coordinating residues only — these are what the consensus must reproduce.
OBSERVED_SIDECHAIN_LIGANDS = {
    807: "D", 810: "E", 823: "N",        # cbEGF9  (2W86)
    910: "D", 913: "E", 928: "N",        # cbEGF10 (2W86)
    1070: "D", 1073: "E", 1088: "N",     # cbEGF12 (1LMJ)
    1113: "D", 1116: "E", 1131: "N",     # cbEGF13 (1LMJ)
    1487: "D", 1490: "E", 1504: "N",     # cbEGF22 (1UZJ)
    1606: "D", 1609: "E", 1624: "N",     # cbEGF23 (1UZJ)
    2127: "D", 2130: "E", 2144: "N",     # cbEGF32 (1EMN)
    2166: "D", 2169: "E", 2183: "N",     # cbEGF33 (1EMN)
}

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def _seq() -> str:
    return json.loads(sorted(RAW.glob("uniprot_P35555_v*.json"))[-1]
                      .read_text(encoding="utf-8"))["sequence"]["value"]


def _read(p: Path) -> list[dict]:
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    t = p.read_text(encoding="utf-8")
    return list(csv.DictReader(t.splitlines(), delimiter="\t")) if t.strip() else []


def _ann() -> list[dict]:
    return _read(PROC / "variant_annotated.tsv")


def _true(v) -> bool:
    return str(v).strip().lower() == "true"


# --------------------------------------------------------------------------
def assert_outputs_exist() -> str:
    for p in (PROC / "variant_annotated.tsv", INTERIM / "cbegf_sites.tsv",
              INTERIM / "disulfide_pairing.tsv", ROOT / "results/phase3_crosstab.md"):
        if not p.is_file():
            raise GateFailure(f"missing {p.relative_to(ROOT)}")
    n = len(_ann())
    if n == 0:
        raise GateFailure("variant_annotated.tsv is empty")
    return f"4 outputs present; {n:,} annotated variants"


def assert_counts_reconcile() -> str:
    before = len(_read(PROC / "variant_master.tsv"))
    after = _ann()
    if len(after) != before:
        raise GateFailure(f"{len(after):,} annotated != {before:,} from Phase 2 — "
                          "annotation must not add or drop rows")
    ids = [v["variation_id"] for v in after]
    if len(set(ids)) != len(ids):
        dup = [i for i, c in Counter(ids).items() if c > 1]
        raise GateFailure(f"{len(dup)} duplicate variation_ids: {dup[:5]}")
    return f"{len(after):,} in == {before:,} out, ids unique"


def assert_cysteine_annotations() -> str:
    seq = _seq()
    rows = _ann()
    bad = [(v["variation_id"], v["position"], v["wt_aa"], v["cys_role"])
           for v in rows if v["cys_role"] and v["wt_aa"] != "C"]
    if bad:
        raise GateFailure(f"{len(bad)} variants carry a cysteine role but WT is not Cys: {bad[:5]}")
    bad2 = [(v["variation_id"], v["position"]) for v in rows
            if v["cys_role"] and seq[int(v["position"]) - 1] != "C"]
    if bad2:
        raise GateFailure(f"{len(bad2)} cysteine-role positions are not Cys in P35555: {bad2[:5]}")
    # cys_removing / cys_creating must be mutually exclusive and consistent with the residues
    wrong = [v["variation_id"] for v in rows
             if _true(v["cys_removing"]) and (v["wt_aa"] != "C" or v["mut_aa"] == "C")]
    wrong += [v["variation_id"] for v in rows
              if _true(v["cys_creating"]) and (v["mut_aa"] != "C" or v["wt_aa"] == "C")]
    if wrong:
        raise GateFailure(f"{len(wrong)} inconsistent cys_removing/cys_creating flags: {wrong[:5]}")
    n = sum(1 for v in rows if v["cys_role"])
    return f"all {n:,} cysteine-site variants verified as Cys; removing/creating flags consistent"


def assert_cbegf_cysteine_count() -> str:
    sites = _read(INTERIM / "cbegf_sites.tsv")
    if len(sites) != 43:
        raise GateFailure(f"{len(sites)} cbEGF domains in the site table, expected 43")
    cys = {int(p) for r in sites for p in r["cys"].split(";")}
    if len(cys) != 258:
        raise GateFailure(f"{len(cys)} cbEGF cysteines, the ClinGen VCEP states 258")
    seq = _seq()
    bad = [p for p in cys if seq[p - 1] != "C"]
    if bad:
        raise GateFailure(f"{len(bad)} listed cbEGF cysteines are not Cys: {bad[:5]}")
    return "43 cbEGF domains carrying exactly 258 cysteines, matching the VCEP"


def assert_disulfide_pattern() -> str:
    pairs = _read(INTERIM / "disulfide_pairing.tsv")
    cb = [p for p in pairs if _true(p["in_cbegf"])]
    if not cb:
        raise GateFailure("no cbEGF disulfides classified")
    counts = Counter(p["pairing"] for p in cb)
    unexpected = set(counts) - {"C1-C3", "C2-C4", "C5-C6"}
    if unexpected:
        raise GateFailure(f"cbEGF disulfides violating the 1-3/2-4/5-6 rule: {sorted(unexpected)}")
    if not (counts["C1-C3"] == counts["C2-C4"] == counts["C5-C6"] == 43):
        raise GateFailure(f"expected 43 of each cbEGF disulfide type, got {dict(counts)}")
    return f"{len(cb)} cbEGF disulfides: 43 x C1-C3, 43 x C2-C4, 43 x C5-C6"


def assert_consensus_motif() -> str:
    seq = _seq()
    sites = _read(INTERIM / "cbegf_sites.tsv")
    rules = {"ca_D1": "D", "ca_DN2": "DN", "ca_EH": "EH", "ca_DN_bOH": "DN", "ca_YF": "YF"}
    bad = []
    for r in sites:
        for col, allowed in rules.items():
            pos, aa = r[col].split(":")
            if seq[int(pos) - 1] != aa:
                bad.append((col, pos, aa, "table/sequence disagree"))
            elif aa not in allowed:
                bad.append((col, pos, aa, f"not in [{allowed}]"))
    if bad:
        raise GateFailure(f"{len(bad)} consensus positions violate the VCEP motif: {bad[:5]}")
    return f"all {len(sites) * 5} derived consensus positions satisfy [D]/[D/N]/[E/H]/[D/N]/[Y/F]"


def assert_consensus_vs_structure() -> str:
    """Replacement for the unavailable UniProt Ca-annotation cross-check."""
    seq = _seq()
    sites = _read(INTERIM / "cbegf_sites.tsv")
    derived = set()
    for r in sites:
        for col in ("ca_D1", "ca_EH", "ca_DN_bOH"):   # the sidechain ligands
            derived.add(int(r[col].split(":")[0]))
    missed = {p: aa for p, aa in OBSERVED_SIDECHAIN_LIGANDS.items() if p not in derived}
    if missed:
        raise GateFailure(
            f"{len(missed)} calcium ligands OBSERVED in experimental structures are not in the "
            f"derived consensus set: {missed}"
        )
    wrong_res = {p: (aa, seq[p - 1]) for p, aa in OBSERVED_SIDECHAIN_LIGANDS.items()
                 if seq[p - 1] != aa}
    if wrong_res:
        raise GateFailure(f"observed ligand residues disagree with P35555: {wrong_res}")
    return (f"all {len(OBSERVED_SIDECHAIN_LIGANDS)} sidechain Ca ligands measured in 2W86/1LMJ/"
            "1UZJ/1EMN are recovered by the derived consensus (8 domains, 16 sites)")


def assert_crosstab_sane() -> str:
    rows = _ann()
    details = []
    for tier in ("primary", "sensitivity"):
        col = f"set_{tier}"
        path = [v for v in rows if v[col] == "pathogenic"]
        ben = [v for v in rows if v[col] == "benign"]
        if not path or not ben:
            raise GateFailure(f"{tier}: pathogenic={len(path)} benign={len(ben)} — "
                              "cross-tab needs both classes non-empty")
        p_crit = sum(1 for v in path if _true(v["is_critical_residue"]))
        b_crit = sum(1 for v in ben if _true(v["is_critical_residue"]))
        if p_crit == 0:
            raise GateFailure(f"{tier}: no pathogenic variant hits a critical residue — "
                              "the annotation is almost certainly broken")
        # Direction check: the study hypothesis is enrichment, so a REVERSED direction here
        # would signal a labelling or annotation inversion rather than a discovery.
        if (p_crit / len(path)) < (b_crit / len(ben)):
            raise GateFailure(
                f"{tier}: critical-residue rate is HIGHER in benign ({b_crit}/{len(ben)}) than "
                f"pathogenic ({p_crit}/{len(path)}) — suspect inverted labels or annotation")
        details.append(f"{tier}: P {p_crit}/{len(path)} vs B {b_crit}/{len(ben)}")
    # A variant cannot be both cysteine-removing and cysteine-creating.
    both = [v["variation_id"] for v in rows
            if _true(v["cys_removing"]) and _true(v["cys_creating"])]
    if both:
        raise GateFailure(f"{len(both)} variants flagged both cys_removing and cys_creating")
    return "; ".join(details)


def assert_no_silent_nans() -> str:
    rows = _ann()
    required = ["domain_kind", "in_cbegf", "in_neonatal_region", "is_ca_consensus",
                "is_cysteine_site", "cys_removing", "cys_creating", "is_critical_residue"]
    empty = {c: sum(1 for v in rows if v.get(c) in ("", None)) for c in required}
    bad = {c: n for c, n in empty.items() if n}
    if bad:
        raise GateFailure(f"annotation columns with empty values: {bad}")
    kinds = Counter(v["domain_kind"] for v in rows)
    if "cbEGF" not in kinds:
        raise GateFailure("no variants mapped into cbEGF domains")
    return ("all annotation columns populated; domain coverage "
            + ", ".join(f"{k}={v:,}" for k, v in kinds.most_common()))


# --------------------------------------------------------------------------
def test_outputs_exist(): assert_outputs_exist()
def test_counts_reconcile(): assert_counts_reconcile()
def test_cysteine_annotations(): assert_cysteine_annotations()
def test_cbegf_cysteine_count(): assert_cbegf_cysteine_count()
def test_disulfide_pattern(): assert_disulfide_pattern()
def test_consensus_motif(): assert_consensus_motif()
def test_consensus_vs_structure(): assert_consensus_vs_structure()
def test_crosstab_sane(): assert_crosstab_sane()
def test_no_silent_nans(): assert_no_silent_nans()


def main() -> int:
    print("PHASE 3 GATE — tests/test_annotation.py")
    print(f"root: {ROOT}\n")
    checks = [
        ("outputs_exist", assert_outputs_exist),
        ("counts_reconcile", assert_counts_reconcile),
        ("cysteine_annotations", assert_cysteine_annotations),
        ("cbegf_cysteine_count", assert_cbegf_cysteine_count),
        ("disulfide_pattern", assert_disulfide_pattern),
        ("consensus_motif", assert_consensus_motif),
        ("consensus_vs_structure", assert_consensus_vs_structure),
        ("crosstab_sane", assert_crosstab_sane),
        ("no_silent_nans", assert_no_silent_nans),
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
        print(f"GATE FAILED — {len(_failures)}/{len(checks)} failed: {', '.join(_failures)}")
        print("Phase 4 must not start until these pass.")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 4 may proceed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
