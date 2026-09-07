#!/usr/bin/env python3
"""
tests/test_env.py — PHASE 0 GATE.

Must PASS before Phase 1 (variant collection) may start.

Asserts:
  assert_scaffold_complete       every directory in the CLAUDE.md layout exists
  assert_imports_succeed         every Python package the pipeline needs imports
  assert_native_tools_callable   pymol / mkdssp / foldx resolve via shutil.which AND execute
  assert_foldx_molecules_present FoldX 5 molecules/ folder sits beside the binary
  assert_secrets_present         .env supplies NCBI creds and is git-ignored (never printed)
  assert_alphamissense_usable    P35555 subset present, parseable, positions in range
  assert_databases_reachable     each REQUIRED database answered a live test query
  assert_lit_notes_complete      lit_notes.md has an entry for every PDF in resources/papers/

Deliberate non-failures (documented findings, not broken environment):
  * AlphaFold DB has no model for P35555 (>2700 aa cap). The gate requires that this be
    *confirmed against a control*, not merely unreachable.
  * AlphaFill has no FBN1 entry (it derives from AlphaFold DB). Phase 4 Tier-3 only.
  * VariantValidator's REST host is unreachable here; the gate requires at least ONE working
    HGVS validator (Mutalyzer or VariantValidator), since Phase 2 needs exactly one.

Run:  python tests/test_env.py          (also importable by pytest)
Exit: 0 all assertions pass, 1 otherwise.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ENV_JSON = ROOT / "logs/00_env_check_latest.json"

# Databases Phase 1-3 cannot proceed without.
REQUIRED_DBS = ["uniprot", "ncbi_clinvar", "pdbe_sifts", "gnomad"]
# Phase 2 needs at least one of these.
HGVS_VALIDATORS = ["mutalyzer", "variantvalidator"]

EXPECTED_DIRS = [
    "resources/papers", "resources/databases", "resources/reference",
    "data/raw", "data/interim", "data/processed",
    "structures/pdb", "structures/alphafold", "structures/ca_transplanted",
    "scripts", "tests", "figures", "results", "logs", "manifest",
    "handoff/cmm/inbox", "handoff/cmm/outbox", "handoff/structures",
]

REQUIRED_MODULES = [
    "pandas", "numpy", "scipy", "requests", "Bio", "biotite", "prody",
    "freesasa", "matplotlib", "seaborn", "statsmodels", "dotenv", "pypdf",
]

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def _load_env_report() -> dict:
    """The env-check JSON carries the live database probes; the gate refuses stale evidence."""
    if not ENV_JSON.is_file():
        raise GateFailure(
            f"{ENV_JSON.relative_to(ROOT)} not found — run `python scripts/00_env_check.py` first"
        )
    report = json.loads(ENV_JSON.read_text(encoding="utf-8"))
    ts = datetime.fromisoformat(report["timestamp_utc"])
    age = datetime.now(timezone.utc) - ts
    if age > timedelta(days=7):
        raise GateFailure(
            f"env-check evidence is {age.days} days old ({ts.isoformat()}); "
            "re-run scripts/00_env_check.py so database reachability reflects today"
        )
    return report


# --------------------------------------------------------------------------
# Assertions
# --------------------------------------------------------------------------
def assert_scaffold_complete() -> str:
    missing = [d for d in EXPECTED_DIRS if not (ROOT / d).is_dir()]
    if missing:
        raise GateFailure(f"missing directories: {missing}")
    return f"{len(EXPECTED_DIRS)} directories present"


def assert_imports_succeed() -> str:
    broken = {}
    for mod in REQUIRED_MODULES:
        try:
            import_module(mod)
        except Exception as exc:  # noqa: BLE001
            broken[mod] = f"{type(exc).__name__}: {exc}"
    if broken:
        raise GateFailure("imports failed:\n" + "\n".join(f"    {m}: {e}" for m, e in broken.items()))
    return f"{len(REQUIRED_MODULES)} modules import cleanly"


def assert_native_tools_callable() -> str:
    """which() alone is not enough — a dangling symlink resolves but cannot run."""
    checks = {
        # tool: (argv, predicate on combined output)
        "pymol": (["pymol", "-cq", "-d", "print(cmd.get_version()[0])"],
                  lambda o: any(c.isdigit() for c in o)),
        "mkdssp": (["mkdssp", "--version"], lambda o: "dssp" in o.lower()),
        # FoldX 5 has no --version flag: it prints its banner and complains. A non-zero
        # exit is expected; the banner is the proof of life.
        "foldx": (["foldx", "--version"], lambda o: "foldx" in o.lower()),
    }
    problems, found = [], []
    for tool, (argv, ok) in checks.items():
        path = shutil.which(tool)
        if path is None:
            problems.append(f"{tool}: not on PATH")
            continue
        try:
            proc = subprocess.run(argv, capture_output=True, text=True, timeout=120, cwd=ROOT)
        except Exception as exc:  # noqa: BLE001
            problems.append(f"{tool}: on PATH at {path} but failed to execute ({exc})")
            continue
        out = (proc.stdout or "") + (proc.stderr or "")
        if not ok(out):
            problems.append(f"{tool}: executed but output unrecognised: {out.strip()[:160]!r}")
        else:
            found.append(f"{tool}@{path}")
    if problems:
        raise GateFailure("native tools not usable:\n" + "\n".join(f"    {p}" for p in problems))
    return ", ".join(found)


def assert_foldx_molecules_present() -> str:
    path = shutil.which("foldx")
    if path is None:
        raise GateFailure("foldx not on PATH")
    molecules = Path(path).resolve().parent / "molecules"
    if not molecules.is_dir():
        raise GateFailure(
            f"{molecules} missing — FoldX 5 needs it (rotabase.txt replacement); "
            "pass --rotabaseLocation= when running from another directory"
        )
    n = len(list(molecules.iterdir()))
    if n == 0:
        raise GateFailure(f"{molecules} exists but is empty")
    return f"{molecules} ({n} files)"


def assert_secrets_present() -> str:
    """Presence and shape only. Values are never read into output."""
    env_file = ROOT / ".env"
    if not env_file.is_file():
        raise GateFailure(".env not found at project root")

    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    if ".env" not in gitignore:
        raise GateFailure(".env is not listed in .gitignore — the NCBI key could be committed")

    keys = {}
    for line in env_file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.strip()

    missing = [k for k in ("NCBI_API_KEY", "NCBI_EMAIL") if not keys.get(k)]
    if missing:
        raise GateFailure(f"missing or empty in .env: {missing}")
    return "NCBI_API_KEY and NCBI_EMAIL present, .env git-ignored (values not shown)"


def assert_alphamissense_usable() -> str:
    import pandas as pd

    subset = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
    if not subset.is_file():
        raise GateFailure(f"missing {subset.relative_to(ROOT)}")

    # The subset was built as "header + rows matching ^P35555", which preserved only the
    # copyright comment, not the real column header — so name the columns explicitly.
    df = pd.read_csv(subset, sep="\t", comment="#", header=None,
                     names=["uniprot_id", "protein_variant", "am_pathogenicity", "am_class"])

    problems = []
    off = df.loc[df["uniprot_id"] != "P35555"]
    if len(off):
        problems.append(f"{len(off)} rows are not P35555, e.g. {off.head(3).to_dict('records')}")

    pos = pd.to_numeric(df["protein_variant"].str[1:-1], errors="coerce")
    bad = df.loc[pos.isna()]
    if len(bad):
        problems.append(f"{len(bad)} unparseable protein_variant values, e.g. "
                        f"{bad['protein_variant'].head(3).tolist()}")
    else:
        oor = df.loc[(pos < 1) | (pos > 2871)]
        if len(oor):
            problems.append(f"{len(oor)} positions outside 1-2871, e.g. "
                            f"{oor['protein_variant'].head(3).tolist()}")

    dup = df.loc[df["protein_variant"].duplicated(keep=False)]
    if len(dup):
        problems.append(f"{len(dup)} duplicate protein_variant rows, e.g. "
                        f"{dup['protein_variant'].head(3).tolist()}")

    if not df["am_class"].isin(["benign", "pathogenic", "ambiguous"]).all():
        problems.append(f"unexpected am_class values: "
                        f"{sorted(set(df['am_class']) - {'benign', 'pathogenic', 'ambiguous'})}")

    nan_counts = df.isna().sum()
    if nan_counts.any():
        problems.append(f"NaNs present: {nan_counts[nan_counts > 0].to_dict()}")

    if problems:
        raise GateFailure("AlphaMissense subset problems:\n" + "\n".join(f"    {p}" for p in problems))

    # Not a failure, but must be visible: a position carrying two wild-type residues would
    # double-count in the Phase 2 join if merged on position alone.
    multi = df.assign(pos=pos, wt=df["protein_variant"].str[0]).groupby("pos")["wt"].nunique()
    multi = multi[multi > 1]
    note = ""
    if len(multi):
        note = (f"; NOTE {len(multi)} multi-WT position(s) {list(multi.index.astype(int))} "
                "— Phase 2 must filter these by WT identity against P35555")
    return f"{len(df):,} rows, positions 1-2871, no dups/NaNs{note}"


def assert_databases_reachable(report: dict) -> str:
    dbs = report.get("databases", {})
    down = [d for d in REQUIRED_DBS if not dbs.get(d)]
    if down:
        detail = {r["check"]: r["detail"] for r in report["results"] if r["check"].startswith("db.")}
        raise GateFailure(
            "required databases unreachable: " + ", ".join(down) + "\n"
            + "\n".join(f"    {k}: {v}" for k, v in detail.items() if k.split(".", 1)[1] in down)
        )

    if not any(dbs.get(v) for v in HGVS_VALIDATORS):
        raise GateFailure(
            "no working HGVS validator: both Mutalyzer and VariantValidator failed. "
            "Phase 2 cannot validate c./p. descriptions without one."
        )

    # AlphaFold DB absence is acceptable ONLY as a confirmed finding, never as a silent failure.
    if not dbs.get("alphafold") and not dbs.get("alphafold_absent_confirmed"):
        raise GateFailure(
            "AlphaFold DB neither returned a model for P35555 nor confirmed its absence "
            "against a control accession — cannot distinguish 'no entry' from 'endpoint broken'"
        )

    ok = [d for d in REQUIRED_DBS if dbs.get(d)]
    validator = next(v for v in HGVS_VALIDATORS if dbs.get(v))
    af = "AFDB absent-confirmed" if dbs.get("alphafold_absent_confirmed") else "AFDB present"
    return f"{', '.join(ok)} reachable; HGVS via {validator}; {af}"


def assert_lit_notes_complete() -> str:
    notes_path = ROOT / "resources/reference/lit_notes.md"
    if not notes_path.is_file():
        raise GateFailure("resources/reference/lit_notes.md not found")
    notes = notes_path.read_text(encoding="utf-8")

    pdfs = sorted(p.name for p in (ROOT / "resources/papers").glob("*.pdf"))
    if not pdfs:
        raise GateFailure("no PDFs found in resources/papers/")

    # An entry counts only if the note names the file (stem or full filename).
    unreferenced = [p for p in pdfs if p not in notes and Path(p).stem not in notes]
    if unreferenced:
        raise GateFailure(
            f"{len(unreferenced)} PDF(s) have no entry in lit_notes.md:\n"
            + "\n".join(f"    {p}" for p in unreferenced)
        )
    return f"all {len(pdfs)} PDFs have an entry"


# --------------------------------------------------------------------------
# pytest-compatible wrappers
# --------------------------------------------------------------------------
def test_scaffold_complete(): assert_scaffold_complete()
def test_imports_succeed(): assert_imports_succeed()
def test_native_tools_callable(): assert_native_tools_callable()
def test_foldx_molecules_present(): assert_foldx_molecules_present()
def test_secrets_present(): assert_secrets_present()
def test_alphamissense_usable(): assert_alphamissense_usable()
def test_databases_reachable(): assert_databases_reachable(_load_env_report())
def test_lit_notes_complete(): assert_lit_notes_complete()


def main() -> int:
    print("PHASE 0 GATE — tests/test_env.py")
    print(f"root: {ROOT}\n")

    try:
        report = _load_env_report()
    except GateFailure as exc:
        print(f"  [FAIL] env-check evidence: {exc}")
        print("\nGATE FAILED — Phase 1 must not start.")
        return 1

    print(f"env-check evidence: {ENV_JSON.relative_to(ROOT)} "
          f"({report['timestamp_utc']}, "
          f"{report['summary']['pass']} pass / {report['summary']['warn']} warn / "
          f"{report['summary']['fail']} fail)\n")

    checks = [
        ("scaffold_complete", assert_scaffold_complete, ()),
        ("imports_succeed", assert_imports_succeed, ()),
        ("native_tools_callable", assert_native_tools_callable, ()),
        ("foldx_molecules_present", assert_foldx_molecules_present, ()),
        ("secrets_present", assert_secrets_present, ()),
        ("alphamissense_usable", assert_alphamissense_usable, ()),
        ("databases_reachable", assert_databases_reachable, (report,)),
        ("lit_notes_complete", assert_lit_notes_complete, ()),
    ]

    for name, fn, args in checks:
        try:
            detail = fn(*args)
            print(f"  [PASS] {name}: {detail}")
        except GateFailure as exc:
            _failures.append(name)
            print(f"  [FAIL] {name}: {exc}")
        except Exception as exc:  # noqa: BLE001 — an unexpected error is still a gate failure
            _failures.append(name)
            print(f"  [FAIL] {name}: unexpected {type(exc).__name__}: {exc}")

    print()
    if _failures:
        print(f"GATE FAILED — {len(_failures)}/{len(checks)} checks failed: "
              f"{', '.join(_failures)}")
        print("Phase 1 must not start until these pass.")
        return 1

    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 1 may proceed on approval.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
