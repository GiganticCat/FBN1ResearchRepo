#!/usr/bin/env python3
"""
07_calibrate_af3.py — PHASE 4. Calibrate AF3-placed calcium against experimental structures.

CLAUDE.md requires this BEFORE AF3-placed calcium may be trusted for un-crystallized domains:
run AF3+Ca2+ on domains that already have experimental structures and confirm it reproduces the
known geometry. If it fails here, the whole Tier-2 plan changes.

Two measurements per domain pair:

  A. CALIBRATION — superpose each AF3 model onto its experimental counterpart on backbone CA
     atoms, then measure how far the AF3-placed Ca2+ sits from the experimentally observed Ca2+.
     Backbone RMSD says "is the fold right"; the Ca displacement says "is the metal right",
     which is the number that actually matters here.

  B. LEAVE-ONE-OUT TRANSPLANT — the Tier-3 cross-check, made rigorous. AlphaFill is unavailable
     for FBN1 (no AlphaFold DB entry, 2871 aa > the 2700 aa cap), so homology transplant is done
     directly: superpose a DONOR experimental cbEGF domain onto a TARGET experimental domain,
     carry the donor's Ca2+ across, and measure how far it lands from the target's real Ca2+.
     Because both are experimental, the answer is a ground-truth error bar for what transplant
     accuracy would be on a domain that has no structure of its own.

Writes:
  results/phase4_calibration.tsv        per-model calibration numbers
  results/phase4_transplant_loo.tsv     leave-one-out transplant errors
  structures/ca_transplanted/provenance.json
  logs/07_calibrate_<stamp>.log
"""

from __future__ import annotations

import json
import sys
import warnings
from datetime import datetime, timezone
from itertools import permutations
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, Superimposer

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

# Experimental reference for each AF3 job: (pdb, author->UniProt offset, chain, construct range)
PAIRS = {
    "fold_fbn1_cbegf9_hyb2_cbegf10_ca2":  ("2w86", 804, "A", (807, 951),  "cbEGF9-hyb2-cbEGF10"),
    "fold_fbn1_cbegf12_13_ca2":           ("1lmj", 1066, "A", (1069, 1154), "cbEGF12-13"),
    "fold_fbn1_cbegf22_tb4_cbegf23_ca2":  ("1uzj", 0,   "A", (1486, 1647), "cbEGF22-TB4-cbEGF23"),
    "fold_fbn1_cbegf32_33_ca2":           ("1emn", 0,   "A", (2127, 2205), "cbEGF32-33"),
}

# Individual experimental cbEGF domains usable as transplant donors/targets:
# (pdb, chain, offset, UniProt domain span, label)
CBEGF_UNITS = [
    ("2w86", "A", 804, (807, 846), "cbEGF9"),
    ("2w86", "A", 804, (910, 951), "cbEGF10"),
    ("1lmj", "A", 1066, (1070, 1112), "cbEGF12"),
    ("1lmj", "A", 1066, (1113, 1154), "cbEGF13"),
    ("1uzj", "A", 0, (1486, 1527), "cbEGF22"),
    ("1uzj", "A", 0, (1605, 1647), "cbEGF23"),
    ("1emn", "A", 0, (2127, 2165), "cbEGF32"),
    ("1emn", "A", 0, (2166, 2205), "cbEGF33"),
]

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m)
    log_lines.append(m)


parser = MMCIFParser(QUIET=True)


def load(path: Path):
    return parser.get_structure(path.stem, str(path))[0]  # first model


def backbone(model, chain: str, offset: int, span=None) -> dict[int, np.ndarray]:
    """UniProt position -> alpha-carbon coordinate. Skips hetero residues (incl. Ca ions)."""
    out = {}
    if chain not in model:
        return out
    for res in model[chain]:
        if res.id[0] != " ":          # hetero flag set => ion/ligand, not a residue
            continue
        if "CA" not in res:
            continue
        up = res.id[1] + offset
        if span and not (span[0] <= up <= span[1]):
            continue
        out[up] = res["CA"].coord.astype(float)
    return out


def calcium(model, span=None, offset: int = 0) -> list[np.ndarray]:
    """Coordinates of every Ca2+ ion; optionally only those nearest a residue span."""
    out = []
    for ch in model:
        for res in ch:
            if res.id[0].startswith("H_") and res.get_resname().strip() == "CA":
                out.append(res["CA"].coord.astype(float))
    return out


def kabsch(mob: np.ndarray, ref: np.ndarray):
    """Return (rotation, translation, rmsd) mapping mob onto ref."""
    mc, rc = mob.mean(0), ref.mean(0)
    M, R = mob - mc, ref - rc
    U, S, Vt = np.linalg.svd(M.T @ R)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    rot = Vt.T @ D @ U.T
    rmsd = float(np.sqrt(((R - (rot @ M.T).T) ** 2).sum() / len(M)))
    return rot, rc - rot @ mc, rmsd


def main() -> int:
    log(f"Phase 4 calibration — {datetime.now(timezone.utc).isoformat()}\n")

    # ---------------- A. AF3 vs experimental --------------------------------
    # Superposition is done PER cbEGF DOMAIN, not per construct. A whole-construct fit mixes two
    # different quantities: how well AF3 placed the metal inside a domain, and how well it guessed
    # the hinge angle between domains. The interdomain term dominates and would make a correct Ca
    # placement look several angstroms wrong. Both are reported, but the per-domain number is the
    # one that answers "can AF3 place calcium?".
    rows = []
    log("A. AF3 calibration against experimental structures (per-domain superposition)")
    log(f"   {'cbEGF':8} {'ref':6} {'model':>5} {'n_res':>5} {'bbRMSD':>7} {'Ca dev':>7} {'pLDDT':>6}")

    for job, (pdb, off, chain, span, label) in PAIRS.items():
        exp = load(ROOT / f"structures/pdb/{pdb}.cif")
        exp_all = calcium(exp)
        # the cbEGF sub-domains contained in this construct
        subs = [u for u in CBEGF_UNITS if u[0] == pdb and span[0] <= u[3][0] and u[3][1] <= span[1]]

        for cif in sorted((ROOT / "structures/alphafold" / job).glob("*_model_*.cif")):
            mdl = load(cif)
            n = int(cif.stem.rsplit("_", 1)[-1])
            af_all = {span[0] + res.id[1] - 1: res for res in mdl["A"]
                      if res.id[0] == " " and "CA" in res}
            af_ca_all = calcium(mdl)

            # whole-construct fit, reported for context only
            exp_bb_full = backbone(exp, chain, off, span)
            common_full = sorted(set(af_all) & set(exp_bb_full))
            _, _, rmsd_full = kabsch(
                np.array([af_all[p]["CA"].coord.astype(float) for p in common_full]),
                np.array([exp_bb_full[p] for p in common_full]))

            for _, ch, doff, dspan, dlabel in subs:
                exp_bb = backbone(exp, ch, doff, dspan)
                af_bb = {p: r["CA"].coord.astype(float) for p, r in af_all.items()
                         if dspan[0] <= p <= dspan[1]}
                common = sorted(set(af_bb) & set(exp_bb))
                if len(common) < 20:
                    continue
                rot, tr, rmsd = kabsch(np.array([af_bb[p] for p in common]),
                                       np.array([exp_bb[p] for p in common]))
                # this domain's own experimental ion = nearest to its centroid
                cen = np.mean(list(exp_bb.values()), axis=0)
                exp_ion = min(exp_all, key=lambda c: float(np.linalg.norm(c - cen)))
                # nearest AF3 ion after transforming into the experimental frame
                moved = [rot @ c + tr for c in af_ca_all]
                dev = min(float(np.linalg.norm(m - exp_ion)) for m in moved)

                site_res = [r for p, r in af_all.items() if dspan[0] <= p <= dspan[1]]
                plddt = [a.bfactor for r in site_res for a in r
                         if min(np.linalg.norm(a.coord - c) for c in af_ca_all) <= 4.0]
                mpl = float(np.mean(plddt)) if plddt else float("nan")

                rows.append({"domain": dlabel, "construct": label, "job": job,
                             "reference_pdb": pdb.upper(), "model": n,
                             "n_residues": len(common),
                             "domain_backbone_rmsd_A": round(rmsd, 3),
                             "construct_backbone_rmsd_A": round(rmsd_full, 3),
                             "ca_deviation_A": round(dev, 3),
                             "mean_plddt_ca_site": round(mpl, 1)})
                log(f"   {dlabel:8} {pdb.upper():6} {n:>5} {len(common):5} {rmsd:7.3f} "
                    f"{dev:7.2f} {mpl:6.1f}")

    hdr = list(rows[0])
    (ROOT / "results/phase4_calibration.tsv").write_text(
        "\t".join(hdr) + "\n" + "\n".join("\t".join(str(r[c]) for c in hdr) for r in rows) + "\n",
        encoding="utf-8")

    devs_all = np.array([r["ca_deviation_A"] for r in rows])
    dom_rmsd = np.array([r["domain_backbone_rmsd_A"] for r in rows])
    con_rmsd = np.array([r["construct_backbone_rmsd_A"] for r in rows])
    log(f"\n   {len(rows)} domain-model comparisons across {len({r['domain'] for r in rows})} cbEGF domains")
    log(f"   per-domain backbone RMSD : median {np.median(dom_rmsd):.2f} A, max {dom_rmsd.max():.2f} A")
    log(f"   whole-construct RMSD     : median {np.median(con_rmsd):.2f} A "
        f"(inflated by interdomain hinge, not metal error)")
    log(f"   Ca deviation             : median {np.median(devs_all):.2f} A, "
        f"mean {devs_all.mean():.2f} A, max {devs_all.max():.2f} A")
    best = {}
    for r in rows:
        if r["domain"] not in best or r["ca_deviation_A"] < best[r["domain"]]:
            best[r["domain"]] = r["ca_deviation_A"]
    log("   best model per domain: " + ", ".join(f"{k} {v:.2f}" for k, v in sorted(best.items())))
    allmax = devs_all

    # ---------------- B. leave-one-out transplant ---------------------------
    log("\nB. Leave-one-out homology transplant (Tier-3 cross-check)")
    log("   Superpose a donor cbEGF domain onto a target, carry the donor Ca2+ across,")
    log("   and measure the error against the target's own experimental Ca2+.")

    units = []
    for pdb, chain, off, span, label in CBEGF_UNITS:
        m = load(ROOT / f"structures/pdb/{pdb}.cif")
        bb = backbone(m, chain, off, span)
        # the Ca ion belonging to THIS domain = nearest ion to the domain centroid
        cen = np.mean(list(bb.values()), axis=0)
        ions = calcium(m)
        if not ions or len(bb) < 25:
            continue
        own = min(ions, key=lambda c: float(np.linalg.norm(c - cen)))
        units.append({"label": label, "pdb": pdb.upper(), "bb": bb, "ca": own, "span": span})
    log(f"   usable experimental cbEGF units: {len(units)}")

    loo = []
    for tgt in units:
        for don in units:
            if don["label"] == tgt["label"]:
                continue
            # align by position WITHIN the domain, since different domains have different numbering
            dk, tk = sorted(don["bb"]), sorted(tgt["bb"])
            n = min(len(dk), len(tk))
            if n < 25:
                continue
            dm = np.array([don["bb"][p] for p in dk[:n]])
            tm = np.array([tgt["bb"][p] for p in tk[:n]])
            rot, tr, rmsd = kabsch(dm, tm)
            moved = rot @ don["ca"] + tr
            err = float(np.linalg.norm(moved - tgt["ca"]))
            loo.append({"target": tgt["label"], "target_pdb": tgt["pdb"],
                        "donor": don["label"], "donor_pdb": don["pdb"],
                        "n_aligned": n, "backbone_rmsd_A": round(rmsd, 3),
                        "ca_transplant_error_A": round(err, 3)})

    hdr2 = list(loo[0])
    (ROOT / "results/phase4_transplant_loo.tsv").write_text(
        "\t".join(hdr2) + "\n" + "\n".join("\t".join(str(r[c]) for c in hdr2) for r in loo) + "\n",
        encoding="utf-8")

    errs = np.array([r["ca_transplant_error_A"] for r in loo])
    log(f"   {len(loo)} donor->target transplants")
    log(f"   Ca transplant error: median {np.median(errs):.2f} A, "
        f"mean {errs.mean():.2f} A, 90th pct {np.percentile(errs, 90):.2f} A, max {errs.max():.2f} A")
    bestper = {}
    for r in loo:
        if r["target"] not in bestper or r["ca_transplant_error_A"] < bestper[r["target"]][0]:
            bestper[r["target"]] = (r["ca_transplant_error_A"], r["donor"])
    log("   best donor per target: "
        + ", ".join(f"{t}<-{d} {e:.2f}A" for t, (e, d) in sorted(bestper.items())))

    (ROOT / "structures/ca_transplanted/provenance.json").write_text(json.dumps({
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "method": "leave-one-out homology transplant between experimental cbEGF domains",
        "purpose": ("Tier-3 cross-check. AlphaFill is unavailable for FBN1 because there is no "
                    "AlphaFold DB entry (2871 aa exceeds the 2700 aa cap), so transplant accuracy "
                    "is characterised directly against experimental ground truth."),
        "donors": [{"label": u["label"], "pdb": u["pdb"], "uniprot_span": u["span"]} for u in units],
        "n_transplants": len(loo),
        "ca_transplant_error_A": {"median": round(float(np.median(errs)), 3),
                                  "mean": round(float(errs.mean()), 3),
                                  "p90": round(float(np.percentile(errs, 90)), 3),
                                  "max": round(float(errs.max()), 3)},
        "af3_worst_case_ca_deviation_A": round(float(max(allmax)), 3),
    }, indent=2), encoding="utf-8")

    log("\nwrote results/phase4_calibration.tsv, results/phase4_transplant_loo.tsv")
    log("wrote structures/ca_transplanted/provenance.json")
    (ROOT / f"logs/07_calibrate_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                         encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
