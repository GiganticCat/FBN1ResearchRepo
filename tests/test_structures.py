#!/usr/bin/env python3
"""
tests/test_structures.py — PHASE 4 GATE.

Must PASS before the CheckMyMetal packet is released to the mentor, and before Phase 5.

Asserts:
  assert_sifts_roundtrip        UniProt position -> author number -> residue IDENTITY, for every
                                experimental structure, using its own offset. Catches the exact
                                failure mode CLAUDE.md warns about: 2W86/1LMJ carry local
                                numbering while 1UZJ/1EMN match UniProt, so an offset applied to
                                the wrong entry works for half of them.
  assert_packet_files_exist     every file named in manifest.csv is on disk and non-trivial
  assert_every_packet_has_ca    no metal-free structure may reach CheckMyMetal
  assert_packet_numbering       packet residues really are in UniProt numbering and match P35555
  assert_calibration_ran        calibration + transplant tables exist and are populated
  assert_af3_calibration_passes AF3 Ca placement is within tolerance of experimental ground truth
  assert_transplant_rmsd_sane   leave-one-out transplant errors are below a sane threshold
  assert_provenance_recorded    ca_transplanted/provenance.json documents donors and method
  assert_instructions_written   INSTRUCTIONS.md exists and names the URL, metal and destination

Run:  python tests/test_structures.py
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"

# Structure -> (author->UniProt offset, chain). Verified in Phase 4, recorded in
# resources/reference/structures.md.
OFFSETS = {"2w86": (804, "A"), "1lmj": (1066, "A"), "1uzj": (0, "A"), "1emn": (0, "A")}

# Tolerances. AF3 is judged against crystallographic references, where sub-angstrom agreement is
# achievable; NMR references carry their own ~2 A coordinate spread so they get more room.
AF3_CA_TOL_XRAY = 1.0
AF3_CA_TOL_NMR = 2.5
TRANSPLANT_MAX_A = 5.0
XRAY_REFS = {"2W86", "1UZJ"}

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def _seq() -> str:
    return json.loads(sorted((ROOT / "data/raw").glob("uniprot_P35555_v*.json"))[-1]
                      .read_text(encoding="utf-8"))["sequence"]["value"]


def _read(p: Path) -> list[dict]:
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    return list(csv.DictReader(p.read_text(encoding="utf-8").splitlines(), delimiter="\t"))


AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _cif_residues(path: Path, chain: str) -> dict[int, str]:
    """author_seq_id -> one-letter, first model only."""
    cols, inloop, out = [], False, {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        st = line.strip()
        if st.startswith("_atom_site."):
            cols.append(st.split(".")[1]); inloop = True; continue
        if inloop and st.startswith("ATOM"):
            f = st.split()
            if len(f) < len(cols):
                continue
            d = dict(zip(cols, f))
            if d.get("pdbx_PDB_model_num", "1") != "1":
                continue
            if d.get("auth_asym_id") != chain:
                continue
            comp, num = d.get("label_comp_id"), d.get("auth_seq_id")
            if comp in AA3 and num not in (None, "."):
                try:
                    out[int(num)] = AA3[comp]
                except ValueError:
                    pass
        elif inloop and st.startswith("#") and out:
            break
    return out


def _pdb_residues(path: Path) -> dict[int, str]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM"):
            comp, num = line[17:20].strip(), line[22:26].strip()
            if comp in AA3 and num:
                try:
                    out[int(num)] = AA3[comp]
                except ValueError:
                    pass
    return out


def _pdb_ca_count(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines()
               if line.startswith("HETATM") and line[17:20].strip() == "CA")


# --------------------------------------------------------------------------
def assert_sifts_roundtrip() -> str:
    seq = _seq()
    details = []
    for pdb, (off, chain) in OFFSETS.items():
        f = ROOT / f"structures/pdb/{pdb}.cif"
        if not f.is_file():
            raise GateFailure(f"missing {f.relative_to(ROOT)}")
        obs = _cif_residues(f, chain)
        if not obs:
            raise GateFailure(f"{pdb.upper()}: no residues parsed from chain {chain}")
        checked = mism = 0
        examples = []
        for auth, aa in obs.items():
            up = auth + off
            if not (1 <= up <= len(seq)):
                continue
            checked += 1
            if seq[up - 1] != aa:
                mism += 1
                if len(examples) < 3:
                    examples.append((auth, up, aa, seq[up - 1]))
        if checked < 20:
            raise GateFailure(f"{pdb.upper()}: only {checked} residues in range")
        # 2W86 carries two non-native N-terminal cloning residues; everything else must be exact.
        allowed = 2 if pdb == "2w86" else 0
        if mism > allowed:
            raise GateFailure(
                f"{pdb.upper()}: {mism}/{checked} residues fail the identity round-trip with "
                f"offset {off:+d} (auth, uniprot, pdb_aa, ref_aa): {examples}")
        details.append(f"{pdb.upper()}{off:+d} {checked - mism}/{checked}")
    return "identity round-trip OK: " + ", ".join(details)


def assert_packet_files_exist() -> str:
    man = INBOX / "manifest.csv"
    if not man.is_file():
        raise GateFailure("handoff/cmm/inbox/manifest.csv missing")
    rows = list(csv.DictReader(man.read_text(encoding="utf-8").splitlines()))
    if not rows:
        raise GateFailure("manifest.csv is empty")
    required = {"filename", "source", "construct", "uniprot_range", "n_ca_sites", "numbering"}
    missing_cols = required - set(rows[0])
    if missing_cols:
        raise GateFailure(f"manifest.csv missing columns: {sorted(missing_cols)}")
    missing = [r["filename"] for r in rows if not (INBOX / r["filename"]).is_file()]
    if missing:
        raise GateFailure(f"manifest lists files that do not exist: {missing}")
    tiny = [r["filename"] for r in rows if (INBOX / r["filename"]).stat().st_size < 2000]
    if tiny:
        raise GateFailure(f"suspiciously small packet files: {tiny}")
    orphans = [p.name for p in INBOX.glob("*.pdb") if p.name not in {r["filename"] for r in rows}]
    if orphans:
        raise GateFailure(f"PDBs in inbox not listed in manifest.csv: {orphans}")
    return f"{len(rows)} packet files, all present and listed"


def assert_every_packet_has_ca() -> str:
    rows = list(csv.DictReader((INBOX / "manifest.csv").read_text(encoding="utf-8").splitlines()))
    bad, total = [], 0
    for r in rows:
        n = _pdb_ca_count(INBOX / r["filename"])
        total += n
        if n == 0:
            bad.append(r["filename"])
        elif n != int(r["n_ca_sites"]):
            bad.append(f"{r['filename']} (manifest says {r['n_ca_sites']}, file has {n})")
    if bad:
        raise GateFailure("a metal-free or miscounted model must never reach CheckMyMetal: " +
                          str(bad))
    return f"all {len(rows)} files carry Ca2+; {total} sites total"


def assert_packet_numbering() -> str:
    """The packet claims UniProt numbering — verify against the reference sequence."""
    seq = _seq()
    rows = list(csv.DictReader((INBOX / "manifest.csv").read_text(encoding="utf-8").splitlines()))
    problems = []
    for r in rows:
        if r["numbering"] != "UniProt P35555":
            problems.append(f"{r['filename']}: numbering claimed as {r['numbering']!r}")
            continue
        res = _pdb_residues(INBOX / r["filename"])
        lo, hi = (int(x) for x in r["uniprot_range"].split("-"))
        checked = mism = 0
        for pos, aa in res.items():
            if not (1 <= pos <= len(seq)):
                problems.append(f"{r['filename']}: residue {pos} outside 1-{len(seq)}")
                break
            checked += 1
            if seq[pos - 1] != aa:
                mism += 1
        allowed = 2 if "2W86" in r["filename"] else 0
        if mism > allowed:
            problems.append(f"{r['filename']}: {mism}/{checked} residues disagree with P35555")
        if checked and not (lo - 3 <= min(res) and max(res) <= hi + 3):
            problems.append(f"{r['filename']}: residues span {min(res)}-{max(res)}, "
                            f"manifest says {lo}-{hi}")
    if problems:
        raise GateFailure("packet numbering problems:\n  " + "\n  ".join(problems))
    return f"all {len(rows)} files verified as genuine UniProt P35555 numbering"


def assert_calibration_ran() -> str:
    cal = _read(ROOT / "results/phase4_calibration.tsv")
    loo = _read(ROOT / "results/phase4_transplant_loo.tsv")
    if len(cal) < 8:
        raise GateFailure(f"only {len(cal)} calibration rows — expected one per domain per model")
    if len(loo) < 8:
        raise GateFailure(f"only {len(loo)} transplant rows")
    return f"{len(cal)} calibration comparisons, {len(loo)} leave-one-out transplants"


def assert_af3_calibration_passes() -> str:
    cal = _read(ROOT / "results/phase4_calibration.tsv")
    bad = []
    for r in cal:
        dev = float(r["ca_deviation_A"])
        tol = AF3_CA_TOL_XRAY if r["reference_pdb"] in XRAY_REFS else AF3_CA_TOL_NMR
        if dev > tol:
            bad.append((r["domain"], r["reference_pdb"], r["model"], dev, tol))
    if bad:
        raise GateFailure(
            f"{len(bad)} AF3 models place Ca2+ beyond tolerance — AF3 must NOT be used for "
            f"un-crystallized domains until this is resolved: {bad[:5]}")
    xr = [float(r["ca_deviation_A"]) for r in cal if r["reference_pdb"] in XRAY_REFS]
    nm = [float(r["ca_deviation_A"]) for r in cal if r["reference_pdb"] not in XRAY_REFS]
    med = lambda x: sorted(x)[len(x) // 2]
    return (f"X-ray refs median {med(xr):.2f} A (tol {AF3_CA_TOL_XRAY}), "
            f"NMR refs median {med(nm):.2f} A (tol {AF3_CA_TOL_NMR})")


def assert_transplant_rmsd_sane() -> str:
    loo = _read(ROOT / "results/phase4_transplant_loo.tsv")
    bad = [(r["donor"], r["target"], float(r["ca_transplant_error_A"])) for r in loo
           if float(r["ca_transplant_error_A"]) > TRANSPLANT_MAX_A]
    if bad:
        raise GateFailure(f"{len(bad)} transplants exceed {TRANSPLANT_MAX_A} A: {bad[:5]}")
    errs = sorted(float(r["ca_transplant_error_A"]) for r in loo)
    return (f"{len(loo)} transplants, median {errs[len(errs)//2]:.2f} A, max {errs[-1]:.2f} A "
            f"(threshold {TRANSPLANT_MAX_A} A)")


def assert_provenance_recorded() -> str:
    p = ROOT / "structures/ca_transplanted/provenance.json"
    if not p.is_file():
        raise GateFailure("structures/ca_transplanted/provenance.json missing")
    js = json.loads(p.read_text(encoding="utf-8"))
    for k in ("method", "purpose", "donors", "ca_transplant_error_A"):
        if k not in js:
            raise GateFailure(f"provenance.json missing '{k}'")
    if not js["donors"]:
        raise GateFailure("provenance.json records no donor structures")
    for d in js["donors"]:
        if not {"label", "pdb", "uniprot_span"} <= set(d):
            raise GateFailure(f"donor entry incomplete: {d}")
    return f"{len(js['donors'])} donors recorded with accession and span"


def assert_instructions_written() -> str:
    p = ROOT / "handoff/cmm/INSTRUCTIONS.md"
    if not p.is_file():
        raise GateFailure("handoff/cmm/INSTRUCTIONS.md missing")
    t = p.read_text(encoding="utf-8")
    need = ["csgid.org", "calcium", "outbox", "UniProt P35555"]
    missing = [n for n in need if n not in t]
    if missing:
        raise GateFailure(f"INSTRUCTIONS.md does not mention: {missing}")
    rows = list(csv.DictReader((INBOX / "manifest.csv").read_text(encoding="utf-8").splitlines()))
    unlisted = [r["filename"] for r in rows if r["filename"] not in t]
    if unlisted:
        raise GateFailure(f"packet files not named in INSTRUCTIONS.md: {unlisted}")
    return "instructions name the URL, metal, destination and every file"


# --------------------------------------------------------------------------
def test_sifts_roundtrip(): assert_sifts_roundtrip()
def test_packet_files_exist(): assert_packet_files_exist()
def test_every_packet_has_ca(): assert_every_packet_has_ca()
def test_packet_numbering(): assert_packet_numbering()
def test_calibration_ran(): assert_calibration_ran()
def test_af3_calibration_passes(): assert_af3_calibration_passes()
def test_transplant_rmsd_sane(): assert_transplant_rmsd_sane()
def test_provenance_recorded(): assert_provenance_recorded()
def test_instructions_written(): assert_instructions_written()


def main() -> int:
    print("PHASE 4 GATE — tests/test_structures.py")
    print(f"root: {ROOT}\n")
    checks = [
        ("sifts_roundtrip", assert_sifts_roundtrip),
        ("packet_files_exist", assert_packet_files_exist),
        ("every_packet_has_ca", assert_every_packet_has_ca),
        ("packet_numbering", assert_packet_numbering),
        ("calibration_ran", assert_calibration_ran),
        ("af3_calibration_passes", assert_af3_calibration_passes),
        ("transplant_rmsd_sane", assert_transplant_rmsd_sane),
        ("provenance_recorded", assert_provenance_recorded),
        ("instructions_written", assert_instructions_written),
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
        print("Do NOT release the CheckMyMetal packet until these pass.")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Packet may be released.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
