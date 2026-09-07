#!/usr/bin/env python3
"""
67_restore_alphamissense.py — put the two AlphaMissense columns back into a released copy.

The published tables ship without `am_pathogenicity` and `am_class`. Those are Cheng et al.'s
predictions, released under CC BY-NC-SA 4.0, and not redistributing them keeps everything else in
this release cleanly CC BY 4.0. Nothing about the analysis changes -- the values are one join
away, and this script performs it.

  1. Download AlphaMissense_aa_substitutions.tsv.gz from Zenodo record 8208688
     (https://doi.org/10.5281/zenodo.8208688), which is 1.2 GB.
  2. Cut the P35555 rows out of it, which is what this project ever used:

         zcat AlphaMissense_aa_substitutions.tsv.gz \\
           | awk 'NR<=4 || /^P35555/' \\
           > resources/databases/AlphaMissense_FBN1_P35555.tsv

  3. Run this script. It rejoins on (position, wild-type, substituted) exactly as
     `05_normalize_variants.py` did, and it asserts the wild-type residue matches the reference
     at every position, which is the same identity check every other join in this pipeline
     carries. A position that agrees by number but not by residue is a transcript mismatch, and
     it is refused rather than silently accepted.

After this, `tests/test_v2.py` passes in full. Without it, the checks that read an AlphaMissense
value are the only ones that cannot run.

Writes: the columns back into data/interim/ and data/processed/, in place.
Run:    .venv/bin/python scripts/67_restore_alphamissense.py [--check]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AM = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
FASTA = ROOT / "data/raw/uniprot_P35555_v264.fasta"

TARGETS = [
    ROOT / "data/interim/variants_normalized.tsv",
    ROOT / "data/processed/structural_metrics.tsv",
    ROOT / "data/processed/structural_metrics_af3.tsv",
    ROOT / "data/processed/variant_annotated.tsv",
    ROOT / "data/processed/variant_master.tsv",
]

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def load_am() -> dict[tuple[int, str, str], tuple[str, str]]:
    if not AM.is_file():
        sys.exit(f"STOP: {AM.relative_to(ROOT)} not found.\n"
                 f"      See the header of this script for the two commands that create it.")
    out = {}
    with AM.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or not f[0].startswith("P35555"):
                continue
            var = f[1] if len(f[1]) >= 3 and f[1][0].isalpha() else None
            if not var:
                continue
            wt, pos, mut = var[0], var[1:-1], var[-1]
            if not pos.isdigit():
                continue
            out[(int(pos), wt, mut)] = (f[-2], f[-1])
    if not out:
        sys.exit(f"STOP: no P35555 rows parsed from {AM.relative_to(ROOT)}")
    return out


def reference() -> str:
    if not FASTA.is_file():
        return ""
    return "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))


def main() -> int:
    am = load_am()
    ref = reference()
    print(f"AlphaMissense: {len(am):,} P35555 substitutions loaded")
    if ref:
        print(f"reference P35555: {len(ref)} aa, wild-type identity will be checked")
    else:
        print("reference FASTA absent, wild-type identity CANNOT be checked -- "
              "run scripts/02_fetch_uniprot.py first for the full check")

    check_only = "--check" in sys.argv
    for path in TARGETS:
        if not path.is_file():
            print(f"  {path.name:34} absent, skipped")
            continue
        rows = list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))
        if not rows:
            continue
        if "am_pathogenicity" in rows[0]:
            print(f"  {path.name:34} already has the columns, skipped")
            continue

        pos_c = next((c for c in ("position", "uniprot_pos") if c in rows[0]), None)
        wt_c = next((c for c in ("wt_aa", "wt", "ref_aa") if c in rows[0]), None)
        mut_c = next((c for c in ("mut_aa", "alt_aa", "mut") if c in rows[0]), None)
        if not (pos_c and wt_c and mut_c):
            print(f"  {path.name:34} no position/wt/mut columns, skipped")
            continue

        n_hit = n_bad = 0
        for r in rows:
            try:
                pos = int(r[pos_c])
            except (ValueError, TypeError):
                r["am_pathogenicity"] = r["am_class"] = ""
                continue
            wt, mut = (r[wt_c] or "").strip(), (r[mut_c] or "").strip()
            wt, mut = AA3.get(wt.upper(), wt), AA3.get(mut.upper(), mut)
            if ref and 1 <= pos <= len(ref) and wt and ref[pos - 1] != wt:
                n_bad += 1
            hit = am.get((pos, wt, mut))
            r["am_pathogenicity"], r["am_class"] = hit if hit else ("", "")
            n_hit += bool(hit)

        if n_bad:
            sys.exit(f"STOP: {n_bad} rows in {path.name} have a wild-type residue that does not "
                     f"match P35555 at that position. Refusing to write a join built on a "
                     f"transcript mismatch.")

        print(f"  {path.name:34} {n_hit:,}/{len(rows):,} matched"
              + ("  (check only, not written)" if check_only else ""))
        if check_only:
            continue
        cols = list(rows[0].keys())
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    if check_only:
        print("\n--check given, nothing written")
    else:
        print("\nRestored. tests/test_v2.py should now pass in full.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
