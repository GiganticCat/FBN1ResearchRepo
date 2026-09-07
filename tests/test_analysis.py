#!/usr/bin/env python3
"""
tests/test_analysis.py — PHASE 5 GATE.

Must PASS before Phase 6 (figures and manuscript tables).

Asserts:
  assert_cmm_complete          every site in the inbox manifest has a parsed result or an
                               explicit MISSING flag — never silently absent
  assert_cmm_parsed_sanely     parsed CMM parameters are in physically plausible ranges
  assert_counts_reconcile      structural/ddG tables cover exactly the Phase 3 variant set
  assert_no_silent_nans        missing values are labelled by status, not dropped
  assert_structural_identity   every structurally-measured variant's WT matches P35555
  assert_ddg_accounted         every covered variant has a ddG or a recorded reason
  assert_ddg_plausible         ddG values are in a physically sensible range
  assert_stats_present         statistics ran, with effect sizes, CIs and BH correction
  assert_seeds_fixed           the analysis scripts pin a seed
  assert_reproducible          re-running the statistics reproduces identical numbers

Run:  python tests/test_analysis.py
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
INBOX = ROOT / "handoff/cmm/inbox"

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


# --------------------------------------------------------------------------
def assert_cmm_complete() -> str:
    man = list(csv.DictReader((INBOX / "manifest.csv").read_text(encoding="utf-8").splitlines()))
    sites = read(PROC / "cmm_sites.tsv")
    submitted = {r["filename"][:-4] for r in man}
    seen = {r["file"] for r in sites}
    unaccounted = submitted - seen
    if unaccounted:
        raise GateFailure(f"submitted files with no row at all in cmm_sites.tsv: "
                          f"{sorted(unaccounted)}")
    expected_sites = sum(int(r["n_ca_sites"]) for r in man)
    parsed = [r for r in sites if r["cmm_status"] == "parsed"]
    missing = [r for r in sites if r["cmm_status"] == "MISSING"]
    if len(parsed) + len(missing) < len(submitted):
        raise GateFailure("fewer result rows than submitted files")
    if len(parsed) != expected_sites and not missing:
        raise GateFailure(f"{len(parsed)} parsed sites but {expected_sites} were submitted, "
                          "and nothing is flagged MISSING — sites vanished silently")
    return (f"{len(parsed)}/{expected_sites} sites parsed"
            + (f", {len(missing)} explicitly flagged MISSING" if missing else ", none missing"))


def assert_cmm_parsed_sanely() -> str:
    sites = [r for r in read(PROC / "cmm_sites.tsv") if r["cmm_status"] == "parsed"]
    bad = []
    for r in sites:
        cn, val, nv, gr = (f(r["coordination_number"]), f(r["valence"]),
                           f(r["nvecsum"]), f(r["grmsd"]))
        if cn is None or not (3 <= cn <= 12):
            bad.append((r["file"], r["site_id"], f"coordination {cn}"))
        if val is None or not (0 < val < 4):
            bad.append((r["file"], r["site_id"], f"valence {val}"))
        if nv is None or not (0 <= nv <= 2):
            bad.append((r["file"], r["site_id"], f"nVECSUM {nv}"))
        if gr is None or not (0 <= gr <= 90):
            bad.append((r["file"], r["site_id"], f"gRMSD {gr}"))
        if r["metal"] != "Ca":
            bad.append((r["file"], r["site_id"], f"metal {r['metal']!r} — expected Ca"))
    if bad:
        raise GateFailure(f"{len(bad)} implausible CMM values: {bad[:5]}")
    return f"all {len(sites)} sites have plausible coordination/valence/nVECSUM/gRMSD, metal=Ca"


def assert_counts_reconcile() -> str:
    ann = read(PROC / "variant_annotated.tsv")
    smet = read(PROC / "structural_metrics.tsv")
    if len(smet) != len(ann):
        raise GateFailure(f"structural_metrics has {len(smet):,} rows vs {len(ann):,} annotated")
    if {r["variation_id"] for r in smet} != {r["variation_id"] for r in ann}:
        raise GateFailure("structural_metrics variation_ids do not match the annotated set")
    detail = f"structural {len(smet):,}"
    if (PROC / "foldx_ddg.tsv").is_file():
        dd = read(PROC / "foldx_ddg.tsv")
        if len(dd) != len(ann):
            raise GateFailure(f"foldx_ddg has {len(dd):,} rows vs {len(ann):,} annotated")
        detail += f", ddG {len(dd):,}"
    return f"all tables cover exactly the {len(ann):,} Phase 3 variants ({detail})"


def assert_no_silent_nans() -> str:
    smet = read(PROC / "structural_metrics.tsv")
    cov = [r for r in smet if r["structural_coverage"] == "True"]
    if not cov:
        raise GateFailure("no variants have structural coverage")
    for col in ("rsa_pct", "dist_to_nearest_ca_A", "template_pdb"):
        empty = [r["variation_id"] for r in cov if r.get(col) in ("", None)]
        if empty:
            raise GateFailure(f"{len(empty)} covered variants have empty {col}: {empty[:5]}")
    uncov = [r for r in smet if r["structural_coverage"] != "True"]
    stray = [r["variation_id"] for r in uncov if r.get("rsa_pct") not in ("", None)]
    if stray:
        raise GateFailure(f"{len(stray)} uncovered variants carry metrics they should not have")
    return (f"{len(cov):,} covered variants fully populated; "
            f"{len(uncov):,} uncovered carry no metrics")


def assert_structural_identity() -> str:
    seq = json.loads(sorted((ROOT / "data/raw").glob("uniprot_P35555_v*.json"))[-1]
                     .read_text(encoding="utf-8"))["sequence"]["value"]
    cov = [r for r in read(PROC / "structural_metrics.tsv") if r["structural_coverage"] == "True"]
    bad = [(r["variation_id"], r["position"], r["residue"], seq[int(r["position"]) - 1])
           for r in cov if r.get("residue") and r["residue"] != seq[int(r["position"]) - 1]]
    if bad:
        raise GateFailure(f"{len(bad)} structural rows whose residue disagrees with P35555: {bad[:5]}")
    bad2 = [r["variation_id"] for r in cov if r.get("residue") != r["wt_aa"]]
    if bad2:
        raise GateFailure(f"{len(bad2)} rows where the template residue != variant wt_aa: {bad2[:5]}")
    return f"all {len(cov):,} structurally-measured variants verified against P35555"


def assert_ddg_accounted() -> str:
    p = PROC / "foldx_ddg.tsv"
    if not p.is_file():
        raise GateFailure("data/processed/foldx_ddg.tsv missing — FoldX stage did not complete")
    rows = read(p)
    cov = [r for r in rows if r["structural_coverage"] == "True"]
    unexplained = [r["variation_id"] for r in cov
                   if r["foldx_status"] not in ("ok", "parse_failed", "MISSING")]
    if unexplained:
        raise GateFailure(f"{len(unexplained)} covered variants with no ddG status")
    ok = [r for r in cov if r["foldx_status"] == "ok"]
    failed = [r for r in cov if r["foldx_status"] != "ok"]
    if not ok:
        raise GateFailure("no ddG values were computed at all")
    blank = [r["variation_id"] for r in ok if f(r["ddG_kcal_mol"]) is None]
    if blank:
        raise GateFailure(f"{len(blank)} rows marked ok but carry no ddG value")
    return (f"{len(ok):,}/{len(cov):,} covered variants have ddG; "
            f"{len(failed)} explicitly flagged as failed/missing")


def assert_ddg_plausible() -> str:
    rows = [r for r in read(PROC / "foldx_ddg.tsv") if r["foldx_status"] == "ok"]
    vals = [f(r["ddG_kcal_mol"]) for r in rows]
    out = [(r["variation_id"], v) for r, v in zip(rows, vals) if v is None or not (-15 <= v <= 60)]
    if out:
        raise GateFailure(f"{len(out)} ddG values outside the plausible -15..60 kcal/mol range: "
                          f"{out[:5]}")
    # Cysteine removal must, on average, cost stability — a reversed sign means the disulfide
    # term was misparsed or the mutation strings were built wrong.
    cys = [f(r["ddG_kcal_mol"]) for r in rows if r["cys_removing"] == "True"]
    oth = [f(r["ddG_kcal_mol"]) for r in rows if r["cys_removing"] != "True"]
    if cys and oth:
        import statistics as st
        if st.median(cys) <= st.median(oth):
            raise GateFailure(
                f"median ddG for cysteine-removing variants ({st.median(cys):.2f}) is not above "
                f"other variants ({st.median(oth):.2f}) — suspect misbuilt mutation strings")
    return (f"{len(vals):,} ddG values in range; cysteine-removing median "
            f"{__import__('statistics').median(cys):.2f} vs other "
            f"{__import__('statistics').median(oth):.2f} kcal/mol")


def assert_stats_present() -> str:
    rows = read(ROOT / "results/phase5_statistics.tsv")
    if not rows:
        raise GateFailure("results/phase5_statistics.tsv is empty")
    for col in ("family", "test", "effect", "ci_low", "ci_high", "p_raw", "q_bh"):
        if col not in rows[0]:
            raise GateFailure(f"statistics table missing column {col}")
    noq = [r["test"] for r in rows if r.get("q_bh") in ("", None)]
    if noq:
        raise GateFailure(f"{len(noq)} tests without a BH-corrected q value: {noq[:3]}")
    nociv = [r["test"] for r in rows
             if f(r["ci_low"]) is None and "cmm" not in r["family"].lower()]
    if nociv:
        raise GateFailure(f"{len(nociv)} tests report no confidence interval: {nociv[:3]}")
    fams = {r["family"] for r in rows}
    if not any(x.startswith("A_") for x in fams):
        raise GateFailure("no Part A (biophysics) tests — these are the headline analysis")
    return f"{len(rows)} tests across {len(fams)} families, all with effect size, CI and q"


def assert_seeds_fixed() -> str:
    hits = []
    for s in ("12_statistics.py", "05_normalize_variants.py"):
        t = (ROOT / "scripts" / s).read_text(encoding="utf-8")
        if not re.search(r"seed\s*[=(]\s*\d{4,}", t, re.I):
            raise GateFailure(f"scripts/{s} does not pin a numeric seed")
        hits.append(s)
    return "seeds pinned in " + ", ".join(hits)


def assert_reproducible() -> str:
    """Re-run the statistics and confirm identical output."""
    before = (ROOT / "results/phase5_statistics.tsv").read_text(encoding="utf-8")
    r = subprocess.run([sys.executable, "scripts/12_statistics.py"], cwd=ROOT,
                       capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        raise GateFailure(f"re-running 12_statistics.py failed: {r.stdout[-400:]}")
    after = (ROOT / "results/phase5_statistics.tsv").read_text(encoding="utf-8")
    if before != after:
        raise GateFailure("statistics are not reproducible — a re-run produced different numbers "
                          "despite a fixed seed")
    return "re-running the statistics reproduced byte-identical results"


# --------------------------------------------------------------------------
def test_cmm_complete(): assert_cmm_complete()
def test_cmm_parsed_sanely(): assert_cmm_parsed_sanely()
def test_counts_reconcile(): assert_counts_reconcile()
def test_no_silent_nans(): assert_no_silent_nans()
def test_structural_identity(): assert_structural_identity()
def test_ddg_accounted(): assert_ddg_accounted()
def test_ddg_plausible(): assert_ddg_plausible()
def test_stats_present(): assert_stats_present()
def test_seeds_fixed(): assert_seeds_fixed()
def test_reproducible(): assert_reproducible()


def main() -> int:
    print("PHASE 5 GATE — tests/test_analysis.py")
    print(f"root: {ROOT}\n")
    checks = [
        ("cmm_complete", assert_cmm_complete),
        ("cmm_parsed_sanely", assert_cmm_parsed_sanely),
        ("counts_reconcile", assert_counts_reconcile),
        ("no_silent_nans", assert_no_silent_nans),
        ("structural_identity", assert_structural_identity),
        ("ddg_accounted", assert_ddg_accounted),
        ("ddg_plausible", assert_ddg_plausible),
        ("stats_present", assert_stats_present),
        ("seeds_fixed", assert_seeds_fixed),
        ("reproducible", assert_reproducible),
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
        print("Phase 6 must not start until these pass.")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 6 may proceed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
