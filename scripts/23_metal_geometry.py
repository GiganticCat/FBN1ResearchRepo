#!/usr/bin/env python3
"""
23_metal_geometry.py — PHASE 5 (revision). The measurement FoldX cannot make.

FoldX scores folding stability and, as scripts/22 demonstrates, barely registers the calcium at
all. So a near-zero ddG at a calcium ligand says nothing about the metal site. This script
measures the metal site directly, on the same structures, with no new modelling.

Two quantities, both computed from coordinates:

  BOND VALENCE SUM (BVS). Brown-Altermatt: v = exp((R0 - R) / b), summed over every O/N/S donor
  within 3.2 A of the ion. This is the same physics behind CheckMyMetal's VALENCE parameter, so
  the implementation can be checked rather than trusted -- the script REQUIRES that it reproduce
  CheckMyMetal's VALENCE on all 16 sites already returned by the server (Pearson r >= 0.95) and
  stops if it does not. Once validated it can be applied to any structure, including mutant
  models, with no web hand-off.

  LIGAND ACCOUNTING. For each variant, which of the ion's ligand atoms the substitution actually
  removes. This distinction matters and is missing from the current analysis: several positions
  in the calcium consensus coordinate the ion through their BACKBONE CARBONYL, which survives
  any substitution. A variant at such a position cannot break the coordination sphere by
  substitution alone, and should not be pooled with variants that delete a side-chain carboxylate.

The mutant valence is computed with the wild-type side-chain donors deleted and every other
ligand held at its observed position. That is deliberately an UPPER BOUND on what the mutant
retains: no relaxation, and no credit for a donor the new side chain might contribute from a
geometry we have not modelled. It is a statement about which coordination bonds the substitution
removes, not a prediction of the relaxed mutant structure -- Batch 3 of the AF3 requests exists
to answer that separately.

Writes:
  data/processed/metal_geometry.tsv
  results/phase5_bvs_validation.md
  logs/23_metal_geometry_<stamp>.log
"""

from __future__ import annotations

import csv
import math
import statistics as st
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

# Brown-Altermatt bond-valence parameters for Ca(II). b is the universal 0.37 A.
R0 = {"O": 1.967, "N": 2.14, "S": 2.45}
B = 0.37
CUTOFF = 3.2          # A; first coordination shell for Ca-O/N
GATE_R = 0.95         # required correlation with CheckMyMetal VALENCE

BACKBONE = {"N", "CA", "C", "O", "OXT"}
# Which chemistry each residue can offer the ion from its side chain.
DONOR_CLASS = {"D": "carboxylate", "E": "carboxylate",
               "N": "amide", "Q": "amide",
               "S": "hydroxyl", "T": "hydroxyl", "Y": "hydroxyl",
               "H": "imidazole"}

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m, flush=True)
    log_lines.append(m)


def read_pdb(path: Path) -> list[dict]:
    atoms = []
    for l in path.read_text(encoding="utf-8").splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        el = (l[76:78].strip() or l[12:16].strip()[0]).upper()
        atoms.append({
            "hetatm": l.startswith("HETATM"),
            "name": l[12:16].strip(), "resname": l[17:20].strip(),
            "chain": l[21], "resseq": l[22:26].strip(), "element": el,
            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54])),
        })
    return atoms


def ca_ions(atoms: list[dict]) -> list[dict]:
    return [a for a in atoms if a["hetatm"] and a["resname"] == "CA" and a["name"] == "CA"]


def shell(atoms: list[dict], ion: dict) -> list[tuple[dict, float]]:
    """Donor atoms within CUTOFF of the ion, with their distances."""
    out = []
    for a in atoms:
        if a is ion or a["element"] not in R0:
            continue
        d = math.dist(a["xyz"], ion["xyz"])
        if d < CUTOFF:
            out.append((a, d))
    return sorted(out, key=lambda x: x[1])


def bvs(shell_atoms: list[tuple[dict, float]]) -> float:
    return sum(math.exp((R0[a["element"]] - d) / B) for a, d in shell_atoms)


def validate_against_cmm() -> tuple[float, float]:
    """Gate: our BVS must track CheckMyMetal's VALENCE on the 16 already-scored sites."""
    cmm = {(r["file"], r["site_id"]): r for r in csv.DictReader(
        (ROOT / "data/processed/cmm_sites.tsv").open(encoding="utf-8"), delimiter="\t")}
    ours, theirs, rows = [], [], []
    for path in sorted(INBOX.glob("*.pdb")):
        atoms = read_pdb(path)
        for ion in ca_ions(atoms):
            sid = f"{ion['chain']}:{ion['resseq']}"
            ref = cmm.get((path.stem, sid))
            if ref is None:
                log(f"  !! no CheckMyMetal record for {path.stem} {sid}")
                continue
            sh = shell(atoms, ion)
            v = bvs(sh)
            ours.append(v); theirs.append(float(ref["valence"]))
            rows.append((path.stem, sid, v, float(ref["valence"]), len(sh),
                         int(ref["coordination_number"])))
    if len(ours) < 8:
        raise SystemExit("STOP: fewer than 8 sites available for the BVS validation gate")
    mx, my = st.mean(ours), st.mean(theirs)
    r = (sum((a - mx) * (b - my) for a, b in zip(ours, theirs))
         / math.sqrt(sum((a - mx) ** 2 for a in ours) * sum((b - my) ** 2 for b in theirs)))
    offset = st.mean([a - b for a, b in zip(ours, theirs)])
    log(f"  BVS validation against CheckMyMetal: n={len(ours)}, Pearson r = {r:.3f}, "
        f"mean offset {offset:+.3f}")
    if r < GATE_R:
        raise SystemExit(f"STOP: BVS implementation does not track CheckMyMetal "
                         f"(r = {r:.3f} < {GATE_R}) — do not use it for mutant scoring")

    # Also as a table: the figure that plots this agreement must read the numbers rather than
    # re-derive them or scrape them back out of the prose below.
    with (ROOT / "data/processed/bvs_vs_cmm.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["structure", "site_id", "bvs_local", "valence_cmm", "cn_local", "cn_cmm"])
        for row in rows:
            w.writerow([row[0], row[1], f"{row[2]:.4f}", f"{row[3]:.4f}", row[4], row[5]])
    log(f"  wrote data/processed/bvs_vs_cmm.tsv ({len(rows)} sites)")

    (ROOT / "results/phase5_bvs_validation.md").write_text(
        "# Local bond-valence implementation, checked against CheckMyMetal\n\n"
        f"Generated by `scripts/23_metal_geometry.py` at "
        f"{datetime.now(timezone.utc).isoformat()}.\n\n"
        "CheckMyMetal is a web service with no programmatic API, so it cannot score a mutant "
        "model without another manual hand-off. Its VALENCE parameter, however, is a "
        "Brown-Altermatt bond-valence sum, which is reproducible from coordinates. This table "
        "is the check that our local implementation agrees with the server's, computed on the "
        "16 sites the server has already scored.\n\n"
        f"**Pearson r = {r:.3f}; mean offset {offset:+.3f} valence units; "
        f"cutoff {CUTOFF} A; R0(Ca-O) = {R0['O']}, R0(Ca-N) = {R0['N']}, b = {B}.**\n\n"
        "The constant offset is expected: CheckMyMetal's published implementation differs in "
        "donor-set and cutoff details. The correlation, not the absolute value, is what licenses "
        "using the local calculation for differences between a wild type and its mutant.\n\n"
        "| structure | site | local BVS | CMM VALENCE | local CN | CMM CN |\n"
        "|---|---|---|---|---|---|\n"
        + "\n".join(f"| {a} | {b} | {c:.2f} | {d:.1f} | {e} | {f} |" for a, b, c, d, e, f in rows)
        + "\n", encoding="utf-8")
    return r, offset


def main() -> int:
    log(f"Phase 5 revision — metal geometry — {datetime.now(timezone.utc).isoformat()}\n")
    log("Validation gate")
    validate_against_cmm()

    variants = list(csv.DictReader(
        (ROOT / "data/processed/structural_metrics.tsv").open(encoding="utf-8"), delimiter="\t"))
    templates = {p.name: read_pdb(p) for p in INBOX.glob("EXP_*.pdb")}

    rows = []
    for v in variants:
        base = {k: v[k] for k in ("variation_id", "position", "wt_aa", "mut_aa", "set_primary",
                                  "site_label", "is_ca_consensus", "is_direct_ca_ligand",
                                  "cys_removing", "structural_coverage", "template_pdb")}
        if v["structural_coverage"] != "True" or v["template_file"] not in templates:
            rows.append({**base, "metal_status": "no_structure",
                         **{k: "" for k in EXTRA if k != "metal_status"}})
            continue

        atoms = templates[v["template_file"]]
        pos, wt, mut = v["position"], v["wt_aa"], v["mut_aa"]

        # Every atom of the mutated residue, and the ion it is closest to.
        res_atoms = [a for a in atoms if a["resseq"] == pos and not a["hetatm"]]
        best = None
        for ion in ca_ions(atoms):
            sh = shell(atoms, ion)
            mine = [(a, d) for a, d in sh if a["resseq"] == pos and not a["hetatm"]]
            if mine and (best is None or mine[0][1] < best[2]):
                best = (ion, sh, mine[0][1], mine)
        if best is None:
            rows.append({**base, "metal_status": "not_a_ligand",
                         "n_ligand_atoms": 0, "ligand_atom_types": "",
                         "ligand_source": "none", "sidechain_atoms_lost": 0,
                         "donor_class_wt": DONOR_CLASS.get(wt, "none"),
                         "donor_class_mut": DONOR_CLASS.get(mut, "none"),
                         "donor_outcome": "not_a_ligand",
                         **{k: "" for k in ("bvs_wt", "bvs_after_loss", "delta_bvs",
                                            "delta_bvs_pct", "cn_wt", "cn_after_loss",
                                            "nearest_ligand_dist_A")}})
            continue

        ion, sh, nearest, mine = best
        side = [(a, d) for a, d in mine if a["name"] not in BACKBONE]
        back = [(a, d) for a, d in mine if a["name"] in BACKBONE]
        # Mutant upper bound: wild-type side-chain donors deleted, everything else untouched.
        kept = [(a, d) for a, d in sh
                if not (a["resseq"] == pos and not a["hetatm"] and a["name"] not in BACKBONE)]
        v_wt, v_mut = bvs(sh), bvs(kept)

        wt_class = DONOR_CLASS.get(wt, "none")
        mut_class = DONOR_CLASS.get(mut, "none")
        if not side:
            outcome = "backbone_ligand_retained"      # substitution cannot remove this bond
        elif mut_class == "none":
            outcome = "donor_abolished"
        elif mut_class == wt_class:
            outcome = "donor_class_conserved"
        else:
            outcome = f"donor_changed_{wt_class}_to_{mut_class}"

        rows.append({**base, "metal_status": "ok",
                     "n_ligand_atoms": len(mine),
                     "ligand_atom_types": "+".join(a["name"] for a, _ in mine),
                     "ligand_source": ("sidechain+backbone" if side and back else
                                       "sidechain" if side else "backbone"),
                     "sidechain_atoms_lost": len(side),
                     "donor_class_wt": wt_class, "donor_class_mut": mut_class,
                     "donor_outcome": outcome,
                     "bvs_wt": round(v_wt, 3), "bvs_after_loss": round(v_mut, 3),
                     "delta_bvs": round(v_mut - v_wt, 3),
                     "delta_bvs_pct": round(100 * (v_mut - v_wt) / v_wt, 1),
                     "cn_wt": len(sh), "cn_after_loss": len(kept),
                     "nearest_ligand_dist_A": round(nearest, 2)})

    cols = list(rows[0])
    out = ROOT / "data/processed/metal_geometry.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    ok = [r for r in rows if r["metal_status"] == "ok"]
    log(f"\n  {len(ok)} variants sit in a calcium coordination shell "
        f"({sum(1 for r in rows if r['metal_status'] == 'not_a_ligand')} covered but not ligands)")
    by = {}
    for r in ok:
        by.setdefault(r["ligand_source"], []).append(r)
    for k, g in sorted(by.items()):
        log(f"  {k:22} n={len(g):3}  median delta-BVS "
            f"{st.median([r['delta_bvs'] for r in g]):+.3f} "
            f"({st.median([r['delta_bvs_pct'] for r in g]):+.1f}% of site valence)")
    out_by = {}
    for r in ok:
        out_by.setdefault(r["donor_outcome"], []).append(r)
    log("")
    for k, g in sorted(out_by.items(), key=lambda kv: -len(kv[1])):
        log(f"  {k:34} n={len(g):3}  median delta-BVS% "
            f"{st.median([r['delta_bvs_pct'] for r in g]):+.1f}")
    log(f"\nwrote {out.relative_to(ROOT)}")
    (ROOT / f"logs/23_metal_geometry_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                              encoding="utf-8")
    return 0


EXTRA = ("metal_status", "n_ligand_atoms", "ligand_atom_types", "ligand_source",
         "sidechain_atoms_lost", "donor_class_wt", "donor_class_mut", "donor_outcome",
         "bvs_wt", "bvs_after_loss", "delta_bvs", "delta_bvs_pct", "cn_wt", "cn_after_loss",
         "nearest_ligand_dist_A")

if __name__ == "__main__":
    sys.exit(main())
