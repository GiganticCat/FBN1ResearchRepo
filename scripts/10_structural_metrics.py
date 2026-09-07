#!/usr/bin/env python3
"""
10_structural_metrics.py — PHASE 5. Per-variant structural measurements.

Computed from the experimental structures for every variant that falls inside one, because these
are the measurements the manuscript headlines: they come from coordinates and owe nothing to any
clinical classification.

Per variant:
  * relative solvent accessibility (freesasa) — how buried the residue is
  * secondary structure (DSSP)
  * minimum distance from the residue to the nearest Ca2+ ion, and whether it is a direct ligand
  * disulfide partner and S-S distance, for cysteines
  * which structure supplied the measurement

Templates are the water-free packet PDBs, already renumbered into UniProt P35555 numbering, so
no offset arithmetic is needed anywhere in this script.

Writes:
  data/processed/structural_metrics.tsv
  logs/10_structural_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
SCRATCH = Path("/tmp/claude-1000/-home-harvey-research-fbn1-marfan/"
               "604537d3-2ba9-4704-b90a-04b2c233decc/scratchpad/struct")
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

TEMPLATES = [
    ("EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb", "2W86", (807, 951)),
    ("EXP_1LMJ_cbEGF12-13.pdb", "1LMJ", (1069, 1154)),
    ("EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb", "1UZJ", (1486, 1647)),
    ("EXP_1EMN_cbEGF32-33.pdb", "1EMN", (2127, 2205)),
]

# Tien et al. (2013) theoretical maximum solvent accessibility, used to turn absolute SASA
# into a relative value comparable across residue types.
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225, "E": 223, "G": 104,
          "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155,
          "T": 172, "W": 285, "Y": 263, "V": 174}

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m)
    log_lines.append(m)


def main() -> int:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    log(f"Phase 5 — structural metrics — {datetime.now(timezone.utc).isoformat()}\n")

    variants = list(csv.DictReader(
        (ROOT / "data/processed/variant_annotated.tsv").open(encoding="utf-8"), delimiter="\t"))
    log(f"{len(variants):,} annotated variants")

    import freesasa
    parser = PDBParser(QUIET=True)
    metrics: dict[int, dict] = {}

    for fname, pdb, span in TEMPLATES:
        path = INBOX / fname
        model = parser.get_structure(pdb, str(path))[0]
        chain = next(iter(model))

        # --- calcium positions --------------------------------------------
        cas = [r["CA"].coord.astype(float) for ch in model for r in ch
               if r.id[0].startswith("H_") and r.get_resname().strip() == "CA"]

        # --- SASA -----------------------------------------------------------
        s = freesasa.Structure(str(path))
        res_area = freesasa.calc(s).residueAreas()

        # --- DSSP -------------------------------------------------------------
        ss_map: dict[int, str] = {}
        out = SCRATCH / f"{pdb}.dssp"
        try:
            subprocess.run(["mkdssp", "--output-format", "dssp", str(path), str(out)],
                           capture_output=True, timeout=120, check=True)
            started = False
            for line in out.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("  #  RESIDUE"):
                    started = True
                    continue
                if started and len(line) > 16 and line[13] != "!":
                    try:
                        ss_map[int(line[5:10])] = line[16].strip() or "-"
                    except ValueError:
                        pass
        except Exception as exc:  # noqa: BLE001
            log(f"  WARN {pdb}: DSSP failed ({type(exc).__name__}) — secondary structure blank")

        n = 0
        for res in chain:
            if res.id[0] != " ":
                continue
            pos = res.id[1]
            # Stay inside the declared native span. 2W86 carries two non-native N-terminal
            # cloning residues (renumbered 805-806) that are not FBN1 sequence at all; measuring
            # them would attach real numbers to residues the protein does not have.
            if not (span[0] <= pos <= span[1]):
                continue
            aa = None
            try:
                from Bio.PDB.Polypeptide import protein_letters_3to1 as p31
                aa = p31.get(res.get_resname().capitalize()) or p31.get(res.get_resname())
            except Exception:  # noqa: BLE001
                pass
            if aa is None:
                aa = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
                      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
                      "TYR": "Y", "VAL": "V"}.get(res.get_resname().strip())
            if aa is None:
                continue

            # distance to nearest Ca, over ALL atoms of the residue
            dmin = min((float(np.linalg.norm(a.coord - c)) for a in res for c in cas),
                       default=float("nan"))
            # is it a direct ligand? oxygen/nitrogen within 3.2 A
            lig = any(a.element in ("O", "N") and
                      min(float(np.linalg.norm(a.coord - c)) for c in cas) <= 3.2
                      for a in res) if cas else False

            area = None
            try:
                ra = res_area[chain.id][str(pos)]
                area = ra.total
            except Exception:  # noqa: BLE001
                pass
            rsa = (area / MAXASA[aa] * 100) if (area is not None and aa in MAXASA) else None

            metrics[pos] = {
                "template_pdb": pdb, "template_file": fname,
                "residue": aa,
                "sasa_abs_A2": None if area is None else round(area, 1),
                "rsa_pct": None if rsa is None else round(rsa, 1),
                "buried": None if rsa is None else bool(rsa < 20.0),
                "secondary_structure": ss_map.get(pos, ""),
                "dist_to_nearest_ca_A": round(dmin, 2) if dmin == dmin else None,
                "is_direct_ca_ligand": lig,
            }
            n += 1
        log(f"  {pdb}: {n} residues measured, {len(cas)} Ca ions, "
            f"{sum(1 for p in metrics.values() if p['template_pdb'] == pdb and p['is_direct_ca_ligand'])} direct ligands")

    # --- disulfide geometry -------------------------------------------------
    dis = list(csv.DictReader((ROOT / "data/interim/uniprot_disulfides.tsv").open(encoding="utf-8"),
                              delimiter="\t"))
    partner = {}
    for r in dis:
        a, b = int(r["cys1"]), int(r["cys2"])
        partner[a], partner[b] = b, a

    ss_dist: dict[int, float] = {}
    for fname, pdb, span in TEMPLATES:
        model = parser.get_structure(pdb, str(INBOX / fname))[0]
        chain = next(iter(model))
        sg = {r.id[1]: r["SG"].coord.astype(float) for r in chain
              if r.id[0] == " " and r.get_resname().strip() == "CYS" and "SG" in r}
        for p, c in sg.items():
            q = partner.get(p)
            if q in sg:
                ss_dist[p] = round(float(np.linalg.norm(c - sg[q])), 2)

    # --- join onto variants --------------------------------------------------
    rows, n_cov = [], 0
    for v in variants:
        pos = int(v["position"])
        m = metrics.get(pos)
        row = {"variation_id": v["variation_id"], "position": pos,
               "wt_aa": v["wt_aa"], "mut_aa": v["mut_aa"],
               "set_primary": v["set_primary"], "set_sensitivity": v["set_sensitivity"],
               "cbegf_index": v["cbegf_index"], "site_label": v["site_label"],
               "is_ca_consensus": v["is_ca_consensus"], "cys_role": v["cys_role"],
               "cys_removing": v["cys_removing"], "is_critical_residue": v["is_critical_residue"],
               "am_pathogenicity": v["am_pathogenicity"], "am_class": v["am_class"],
               "structural_coverage": bool(m)}
        if m:
            n_cov += 1
            # residue identity must agree with the reference before any metric is trusted
            if m["residue"] != v["wt_aa"]:
                raise SystemExit(
                    f"STOP: position {pos} is {m['residue']} in {m['template_pdb']} but the "
                    f"variant table says {v['wt_aa']} — numbering error, refusing to continue")
            row |= m
            row["disulfide_partner"] = partner.get(pos, "")
            row["ss_bond_dist_A"] = ss_dist.get(pos, "")
        rows.append(row)

    cols = list({k: None for r in rows for k in r})
    with (ROOT / "data/processed/structural_metrics.tsv").open("w", newline="",
                                                               encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    log(f"\n  structural coverage: {n_cov:,}/{len(variants):,} variants")
    covered = [r for r in rows if r["structural_coverage"]]
    lig = sum(1 for r in covered if r.get("is_direct_ca_ligand"))
    bur = sum(1 for r in covered if r.get("buried"))
    log(f"  of covered: {lig} at direct Ca ligands, {bur} buried (RSA < 20%)")
    log(f"  residue identity verified against P35555 for all {n_cov:,} covered variants")
    log("\nwrote data/processed/structural_metrics.tsv")
    (ROOT / f"logs/10_structural_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                          encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
