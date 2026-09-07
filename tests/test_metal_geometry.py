#!/usr/bin/env python3
"""
tests/test_metal_geometry.py — PHASE 5 REVISION GATE.

Guards the three new metal-aware outputs. Must PASS before any of them reaches the manuscript.

Asserts:
  assert_bvs_tracks_cmm        our bond-valence implementation reproduces CheckMyMetal's VALENCE
                               on the 16 sites the server already scored (Pearson r >= 0.95)
  assert_bvs_recomputes        BVS recomputed from the coordinates matches the stored table
  assert_ligand_accounting     ligand atom lists are chemically possible for the wild-type residue
  assert_backbone_invariant    backbone-only ligand positions lose exactly zero bond valence
  assert_sidechain_loses       side-chain donor deletions lose a non-trivial amount of it
  assert_counts_reconcile      metal table covers exactly the variant set of the upstream table
  assert_no_silent_nans        every row carries a status; blanks only where status explains them
  assert_ablation_complete     if the FoldX ablation has run, every covered variant is accounted
  assert_ablation_isolates     the ablation compares equal protein coordinates (ion removed only)
  assert_am3d_calibrated       if AllMetal3D has run, it recovered the experimental sites first

Run:  python tests/test_metal_geometry.py
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import math
import re
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
INBOX = ROOT / "handoff/cmm/inbox"
RES = ROOT / "results"

R0 = {"O": 1.967, "N": 2.14, "S": 2.45}
B = 0.37
CUTOFF = 3.2
BACKBONE = {"N", "CA", "C", "O", "OXT"}
# Which atoms each residue can offer a metal from its side chain.
SIDE_DONORS = {"D": {"OD1", "OD2"}, "E": {"OE1", "OE2"}, "N": {"OD1", "ND2"},
               "Q": {"OE1", "NE2"}, "S": {"OG"}, "T": {"OG1"}, "Y": {"OH"},
               "H": {"ND1", "NE2"}, "C": {"SG"}, "M": {"SD"}, "K": {"NZ"},
               "R": {"NE", "NH1", "NH2"}, "W": {"NE1"}}

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def read(p: Path) -> list[dict]:
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    return list(csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"))


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def pdb_atoms(path: Path) -> list[dict]:
    out = []
    for l in path.read_text(encoding="utf-8").splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        out.append({"het": l.startswith("HETATM"), "name": l[12:16].strip(),
                    "resname": l[17:20].strip(), "chain": l[21], "resseq": l[22:26].strip(),
                    "el": (l[76:78].strip() or l[12:16].strip()[0]).upper(),
                    "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54]))})
    return out


def bvs_of(atoms, ion, skip_sidechain_of=None) -> tuple[float, int]:
    v, n = 0.0, 0
    for a in atoms:
        if a is ion or a["el"] not in R0:
            continue
        if (skip_sidechain_of is not None and a["resseq"] == skip_sidechain_of
                and not a["het"] and a["name"] not in BACKBONE):
            continue
        d = math.dist(a["xyz"], ion["xyz"])
        if d < CUTOFF:
            v += math.exp((R0[a["el"]] - d) / B)
            n += 1
    return v, n


# --------------------------------------------------------------------------
def assert_bvs_tracks_cmm() -> str:
    doc = RES / "phase5_bvs_validation.md"
    if not doc.is_file():
        raise GateFailure("results/phase5_bvs_validation.md absent — run scripts/23")
    m = re.search(r"Pearson r = ([0-9.]+)", doc.read_text(encoding="utf-8"))
    if not m:
        raise GateFailure("validation document records no correlation with CheckMyMetal")
    r = float(m.group(1))
    if r < 0.95:
        raise GateFailure(f"bond valence tracks CheckMyMetal only at r = {r:.3f}; "
                          "the local implementation may not be measuring the same thing")
    return f"Pearson r = {r:.3f} against CheckMyMetal VALENCE on the scored sites"


def assert_bvs_recomputes() -> str:
    rows = [r for r in read(PROC / "metal_geometry.tsv") if r["metal_status"] == "ok"]
    if not rows:
        raise GateFailure("no rows with metal_status=ok")
    templates = {p.stem.split("_")[1]: pdb_atoms(p) for p in INBOX.glob("EXP_*.pdb")}
    bad, checked = [], 0
    for r in rows:
        atoms = templates.get(r["template_pdb"])
        if atoms is None:
            continue
        ions = [a for a in atoms if a["het"] and a["resname"] == "CA"]
        # The ion this variant touches is the one its own atoms are closest to.
        best, ion = None, None
        for cand in ions:
            d = min((math.dist(a["xyz"], cand["xyz"]) for a in atoms
                     if a["resseq"] == r["position"] and not a["het"]), default=None)
            if d is not None and (best is None or d < best):
                best, ion = d, cand
        if ion is None:
            continue
        wt, _ = bvs_of(atoms, ion)
        mut, _ = bvs_of(atoms, ion, skip_sidechain_of=r["position"])
        checked += 1
        if abs(wt - f(r["bvs_wt"])) > 0.01 or abs(mut - f(r["bvs_after_loss"])) > 0.01:
            bad.append((r["variation_id"], round(wt, 3), r["bvs_wt"]))
    if bad:
        raise GateFailure(f"{len(bad)} rows do not recompute from coordinates: {bad[:3]}")
    return f"{checked} rows recomputed from the deposited coordinates and matched"


def assert_ligand_accounting() -> str:
    rows = [r for r in read(PROC / "metal_geometry.tsv") if r["metal_status"] == "ok"]
    bad = []
    for r in rows:
        atoms = [a for a in r["ligand_atom_types"].split("+") if a]
        allowed = BACKBONE | SIDE_DONORS.get(r["wt_aa"], set())
        extra = set(atoms) - allowed
        if extra:
            bad.append((r["variation_id"], r["wt_aa"], sorted(extra)))
        n = int(r["n_ligand_atoms"])
        if n != len(atoms):
            bad.append((r["variation_id"], "count mismatch", n))
        lost = int(r["sidechain_atoms_lost"])
        if lost != sum(1 for a in atoms if a not in BACKBONE):
            bad.append((r["variation_id"], "sidechain count mismatch", lost))
    if bad:
        raise GateFailure(f"chemically impossible ligand assignments: {bad[:4]}")
    return f"{len(rows)} ligand assignments consistent with the wild-type residue's chemistry"


def assert_backbone_invariant() -> str:
    rows = [r for r in read(PROC / "metal_geometry.tsv")
            if r["metal_status"] == "ok" and r["ligand_source"] == "backbone"]
    if not rows:
        raise GateFailure("no backbone-only ligand positions found — the taxonomy is not "
                          "being applied")
    # Test for None explicitly: 0.0 is falsy, and `or` would rewrite the very value we expect.
    bad = [r["variation_id"] for r in rows
           if f(r["delta_bvs"]) is None or abs(f(r["delta_bvs"])) > 1e-9]
    if bad:
        raise GateFailure(f"backbone-only positions lost bond valence, which no substitution "
                          f"can do: {bad[:5]}")
    return (f"{len(rows)} backbone-carbonyl ligand positions, all with zero valence loss "
            "as chemistry requires")


def assert_sidechain_loses() -> str:
    rows = [r for r in read(PROC / "metal_geometry.tsv")
            if r["metal_status"] == "ok" and r["ligand_source"] != "backbone"]
    if not rows:
        raise GateFailure("no side-chain ligand positions found")
    zero = [r["variation_id"] for r in rows
            if f(r["delta_bvs"]) is None or abs(f(r["delta_bvs"])) < 1e-9]
    if zero:
        raise GateFailure(f"side-chain donors deleted with no valence change: {zero[:5]}")
    med = st.median([f(r["delta_bvs_pct"]) for r in rows])
    return f"{len(rows)} side-chain donor deletions, median {med:.1f}% of site valence lost"


def assert_counts_reconcile() -> str:
    up = {r["variation_id"] for r in read(PROC / "structural_metrics.tsv")}
    mg = read(PROC / "metal_geometry.tsv")
    ids = {r["variation_id"] for r in mg}
    if ids != up:
        raise GateFailure(f"metal_geometry covers {len(ids)} variants, upstream has {len(up)}")
    covered = {r["variation_id"] for r in read(PROC / "structural_metrics.tsv")
               if r["structural_coverage"] == "True"}
    scored = {r["variation_id"] for r in mg if r["metal_status"] in ("ok", "not_a_ligand")}
    if scored != covered:
        raise GateFailure(f"{len(covered - scored)} structurally covered variants have no "
                          "metal verdict")
    return f"{len(ids)} variants, {len(covered)} covered and all given a metal verdict"


def assert_no_silent_nans() -> str:
    rows = read(PROC / "metal_geometry.tsv")
    numeric = ("bvs_wt", "bvs_after_loss", "delta_bvs", "delta_bvs_pct")
    bad = []
    for r in rows:
        blank = [c for c in numeric if r[c] == ""]
        if r["metal_status"] == "ok" and blank:
            bad.append((r["variation_id"], "ok but blank", blank))
        if r["metal_status"] != "ok" and len(blank) != len(numeric):
            bad.append((r["variation_id"], "not ok but populated", r["metal_status"]))
    if bad:
        raise GateFailure(f"values present or absent without a matching status: {bad[:4]}")
    return f"{len(rows)} rows: every value is either present with status ok or labelled"


def assert_ablation_complete() -> str:
    p = PROC / "foldx_ion_ablation.tsv"
    if not p.is_file():
        return "not run yet — skipped (scripts/22 has not produced output)"
    rows = read(p)
    covered = [r for r in rows if r["structural_coverage"] == "True"]
    missing = [r for r in covered if r["ablation_status"] not in ("ok", "failed")]
    if missing:
        raise GateFailure(f"{len(missing)} covered variants with no ablation verdict")
    ok = [r for r in covered if r["ablation_status"] == "ok"]
    if not ok:
        raise GateFailure("no successful ablation results")
    blank = [r["variation_id"] for r in ok if r["ddG_metal_attributable"] == ""]
    if blank:
        raise GateFailure(f"ok rows without a metal-attributable value: {blank[:5]}")
    return (f"{len(ok)}/{len(covered)} covered variants have an ablation result, "
            f"{len(covered) - len(ok)} explicitly failed")


def assert_ablation_isolates() -> str:
    """The control is only meaningful if the two conditions share protein coordinates."""
    src = (ROOT / "scripts/22_foldx_ion_ablation.py").read_text(encoding="utf-8")
    if "RepairPDB" not in src:
        raise GateFailure("ablation script no longer repairs the structure")
    if src.count('foldx(["--command=RepairPDB"') > 1 or "in_Repair.pdb" not in src:
        raise GateFailure("ablation appears to repair each condition separately, which would "
                          "confound the ion with a different repair path")
    if 'l[17:20].strip() == "CA"' not in src:
        raise GateFailure("ablation no longer strips calcium specifically")
    return "ion-free condition is derived from the repaired coordinates, not a second repair"


def assert_am3d_calibrated() -> str:
    doc = RES / "phase5_allmetal3d_calibration.md"
    if not doc.is_file():
        return "not run yet — skipped (scripts/24 has not produced output)"
    txt = doc.read_text(encoding="utf-8")
    m = re.search(r"\*\*(\d+) of (\d+) sites recovered", txt)
    if not m:
        raise GateFailure("calibration document records no recovery count")
    found, total = int(m.group(1)), int(m.group(2))
    if found < total * 0.75:
        raise GateFailure(f"AllMetal3D recovered only {found}/{total} experimental sites; "
                          "it is not calibrated for this domain family and must not be used "
                          "to judge mutants")
    return f"{found}/{total} experimental calcium sites recovered before any mutant was scored"


def main() -> int:
    print("PHASE 5 REVISION GATE — tests/test_metal_geometry.py\n")
    checks = [
        ("bvs_tracks_cmm", assert_bvs_tracks_cmm),
        ("bvs_recomputes", assert_bvs_recomputes),
        ("ligand_accounting", assert_ligand_accounting),
        ("backbone_invariant", assert_backbone_invariant),
        ("sidechain_loses", assert_sidechain_loses),
        ("counts_reconcile", assert_counts_reconcile),
        ("no_silent_nans", assert_no_silent_nans),
        ("ablation_complete", assert_ablation_complete),
        ("ablation_isolates", assert_ablation_isolates),
        ("am3d_calibrated", assert_am3d_calibrated),
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
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
