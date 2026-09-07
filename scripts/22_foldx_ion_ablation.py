#!/usr/bin/env python3
"""
22_foldx_ion_ablation.py — PHASE 5 (revision). How much of the calcium does FoldX actually see?

The manuscript's central claim rests on FoldX reporting near-zero ddG at calcium-coordinating
residues. A reviewer will ask the obvious question: is the fold genuinely undisturbed, or is the
force field blind to the interaction that was destroyed? This script answers it with a control
instead of an argument.

Design. For each template, BuildModel is run twice on the SAME repaired coordinates:

  with-ion   the repaired structure as used in scripts/11 (Ca2+ HETATMs present)
  ion-free   byte-identical protein coordinates, Ca2+ HETATM records deleted

The only difference is the ion, so `ddG(with-ion) - ddG(ion-free)` is the entire contribution the
calcium makes to FoldX's answer for that mutation. Re-repairing without the ion would confound
this with a different repair path, which is why the strip is done AFTER RepairPDB.

Both conditions use numberOfRuns=5, which also answers the second methodological objection:
FoldX is stochastic and is normally reported as a mean over replicate builds, not a single build.
The mean and standard deviation are recorded per mutation and per condition.

`partial covalent bonds` is FoldX's only metal-coordination term and is captured explicitly, so
the fraction of calcium ligands at which FoldX registers any metal interaction at all becomes a
reported number rather than an assumption.

Prototype result on 1EMN (16 ->Ala mutations, 2026-08-08): the metal term fired at 2 of 4 direct
side-chain ligands and 0 of 6 other consensus positions; removing the entire ion shifted ddG by a
median of +0.15 kcal/mol at direct ligands versus +0.06 at controls. This script establishes
whether that holds across all four templates and all 751 covered variants.

Writes:
  data/processed/foldx_ion_ablation.tsv
  logs/22_foldx_ablation_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import statistics as st
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
MOLECULES = Path.home() / "tools/foldx/molecules/"
WORK = Path(os.environ.get("FBN1_FOLDX_WORK", "/tmp/fbn1-foldx-ablation"))
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
N_RUNS = 5

TEMPLATES = [
    ("EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb", "2W86"),
    ("EXP_1LMJ_cbEGF12-13.pdb", "1LMJ"),
    ("EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb", "1UZJ"),
    ("EXP_1EMN_cbEGF32-33.pdb", "1EMN"),
]
CHAIN = "A"
TERMS = ("total energy", "partial covalent bonds", "Electrostatics", "Solvation Polar",
         "Van der Waals", "disulfide")

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m, flush=True)
    log_lines.append(m)


def foldx(args: list[str], cwd: Path, timeout: int = 36000) -> subprocess.CompletedProcess:
    return subprocess.run(["foldx", *args, f"--rotabaseLocation={MOLECULES}"],
                          cwd=cwd, capture_output=True, text=True, timeout=timeout)


def parse_dif(path: Path, n_muts: int, n_runs: int) -> list[dict[str, list[float]]] | None:
    """Return one dict of term -> [value per replicate] per mutation, in list order."""
    text = path.read_text(encoding="utf-8")
    hdr = next((l.split("\t") for l in text.splitlines() if l.startswith("Pdb\t")), None)
    rows = [l.split("\t") for l in text.splitlines() if l.startswith("in_Repair_")]
    if hdr is None or len(rows) != n_muts * n_runs:
        return None
    idx = {c: i for i, c in enumerate(hdr)}
    # FoldX names each row in_Repair_<mutation index>_<replicate index>.pdb. Group on the
    # mutation index rather than assuming row order.
    by_mut: dict[int, list[list[str]]] = {}
    for r in rows:
        try:
            m = int(r[0].removeprefix("in_Repair_").split("_")[0])
        except ValueError:
            return None
        by_mut.setdefault(m, []).append(r)
    if sorted(by_mut) != list(range(1, n_muts + 1)) or any(
            len(v) != n_runs for v in by_mut.values()):
        return None
    return [{t: [float(r[idx[t]]) for r in by_mut[m]] for t in TERMS if t in idx}
            for m in range(1, n_muts + 1)]


def build(wd: Path, muts: list[str], n_muts: int) -> list[dict[str, list[float]]] | None:
    dif = wd / "Dif_in_Repair.fxout"
    if not dif.is_file():                       # idempotent: skip work already done
        (wd / "individual_list.txt").write_text("\n".join(muts) + "\n", encoding="utf-8")
        t0 = time.time()
        r = foldx(["--command=BuildModel", "--pdb=in_Repair.pdb",
                   "--mutant-file=individual_list.txt", f"--numberOfRuns={N_RUNS}"], wd)
        if r.returncode != 0 or not dif.is_file():
            log(f"    !! BuildModel failed (rc={r.returncode}) — {r.stdout[-300:]}")
            return None
        log(f"    {n_muts} mutations x {N_RUNS} runs in {time.time() - t0:.0f}s")
    else:
        log("    using cached BuildModel output")
    return parse_dif(dif, n_muts, N_RUNS)


def main() -> int:
    if shutil.which("foldx") is None:
        raise SystemExit("STOP: foldx not on PATH — ask before omitting ddG (CLAUDE.md 5)")
    WORK.mkdir(parents=True, exist_ok=True)
    log(f"Phase 5 revision — FoldX calcium ablation — {datetime.now(timezone.utc).isoformat()}")
    log(f"  work dir {WORK}; {N_RUNS} replicate builds per condition\n")

    seq = json.loads(sorted((ROOT / "data/raw").glob("uniprot_P35555_v*.json"))[-1]
                     .read_text(encoding="utf-8"))["sequence"]["value"]
    variants = list(csv.DictReader(
        (ROOT / "data/processed/structural_metrics.tsv").open(encoding="utf-8"), delimiter="\t"))
    covered = [v for v in variants if v["structural_coverage"] == "True"]
    log(f"{len(covered):,} variants with structural coverage")

    results: dict[str, dict] = {}

    for fname, pdb in TEMPLATES:
        sub = [v for v in covered if v["template_pdb"] == pdb]
        if not sub:
            continue
        log(f"\n{pdb} — {len(sub)} mutations")
        base = WORK / pdb
        (base / "withca").mkdir(parents=True, exist_ok=True)
        (base / "noion").mkdir(parents=True, exist_ok=True)

        # --- repair once, with the ion present (as in scripts/11) --------------------
        wc = base / "withca"
        if not (wc / "in_Repair.pdb").is_file():
            shutil.copy(INBOX / fname, wc / "in.pdb")
            t0 = time.time()
            r = foldx(["--command=RepairPDB", "--pdb=in.pdb"], wc)
            if r.returncode != 0 or not (wc / "in_Repair.pdb").is_file():
                log(f"    !! RepairPDB failed — {r.stdout[-300:]}")
                continue
            log(f"    RepairPDB {time.time() - t0:.0f}s")
        n_ca = sum(1 for l in (wc / "in_Repair.pdb").read_text(encoding="utf-8").splitlines()
                   if l.startswith("HETATM") and l[17:20].strip() == "CA")
        if n_ca == 0:
            raise SystemExit(f"STOP: {pdb} repaired structure lost its Ca2+ — the ablation "
                             "control is meaningless without it")
        log(f"    repaired structure retains {n_ca} Ca2+ ion(s)")

        # --- identical protein coordinates, ion deleted ------------------------------
        ni = base / "noion"
        if not (ni / "in_Repair.pdb").is_file():
            keep = [l for l in (wc / "in_Repair.pdb").read_text(encoding="utf-8").splitlines()
                    if not (l.startswith("HETATM") and l[17:20].strip() == "CA")]
            (ni / "in_Repair.pdb").write_text("\n".join(keep) + "\n", encoding="utf-8")

        # --- mutation strings, wild type re-verified against P35555 ------------------
        muts, keep_v, bad = [], [], []
        for v in sub:
            pos, wt, mut = int(v["position"]), v["wt_aa"], v["mut_aa"]
            if seq[pos - 1] != wt:
                bad.append((v["variation_id"], pos, wt, seq[pos - 1]))
                continue
            muts.append(f"{wt}{CHAIN}{pos}{mut};")
            keep_v.append(v)
        if bad:
            raise SystemExit(f"STOP: wild-type mismatch before FoldX: {bad[:5]}")

        log("  with ion:")
        a = build(wc, muts, len(keep_v))
        log("  ion removed:")
        b = build(ni, muts, len(keep_v))
        if a is None or b is None:
            for v in keep_v:
                results[v["variation_id"]] = {"ablation_status": "failed"}
            continue

        for v, ra, rb in zip(keep_v, a, b):
            tot_a, tot_b = ra["total energy"], rb["total energy"]
            results[v["variation_id"]] = {
                "ablation_status": "ok",
                "foldx_template": pdb,
                "ddG_withca_mean": round(st.mean(tot_a), 3),
                "ddG_withca_sd": round(st.pstdev(tot_a), 3),
                "ddG_withca_run1": round(tot_a[0], 3),
                "ddG_noion_mean": round(st.mean(tot_b), 3),
                "ddG_noion_sd": round(st.pstdev(tot_b), 3),
                "ddG_metal_attributable": round(st.mean(tot_a) - st.mean(tot_b), 3),
                "partcov_withca": round(st.mean(ra.get("partial covalent bonds", [0.0])), 3),
                "partcov_fired": abs(st.mean(ra.get("partial covalent bonds", [0.0]))) > 0.01,
                "elec_withca": round(st.mean(ra.get("Electrostatics", [0.0])), 3),
                "elec_noion": round(st.mean(rb.get("Electrostatics", [0.0])), 3),
            }

    # --- write ------------------------------------------------------------------
    cols = ["variation_id", "position", "wt_aa", "mut_aa", "set_primary", "site_label",
            "is_ca_consensus", "is_direct_ca_ligand", "cys_removing", "structural_coverage",
            "ablation_status", "foldx_template", "ddG_withca_mean", "ddG_withca_sd",
            "ddG_withca_run1", "ddG_noion_mean", "ddG_noion_sd", "ddG_metal_attributable",
            "partcov_withca", "partcov_fired", "elec_withca", "elec_noion"]
    rows = []
    for v in variants:
        r = results.get(v["variation_id"], {})
        rows.append({
            **{k: v[k] for k in ("variation_id", "position", "wt_aa", "mut_aa", "set_primary",
                                 "site_label", "is_ca_consensus", "is_direct_ca_ligand",
                                 "cys_removing", "structural_coverage")},
            "ablation_status": r.get("ablation_status",
                                     "no_structure" if v["structural_coverage"] != "True"
                                     else "MISSING"),
            **{k: r.get(k, "") for k in cols[11:]},
        })
    out = ROOT / "data/processed/foldx_ion_ablation.tsv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    # --- headline summary --------------------------------------------------------
    ok = [r for r in rows if r["ablation_status"] == "ok"]
    log(f"\n  {len(ok):,} variants with an ablation result; "
        f"{sum(1 for r in rows if r['ablation_status'] == 'MISSING')} covered-but-missing")

    def group(pred):
        return [r for r in ok if pred(r)]

    for name, sel in (("direct Ca ligand", lambda r: r["is_direct_ca_ligand"] == "True"),
                      ("Ca consensus (any)", lambda r: r["is_ca_consensus"] == "True"),
                      ("cysteine-removing", lambda r: r["cys_removing"] == "True"),
                      ("all other covered", lambda r: r["is_ca_consensus"] != "True"
                       and r["cys_removing"] != "True")):
        g = group(sel)
        if not g:
            continue
        log(f"  {name:20} n={len(g):4}  median metal-attributable ddG "
            f"{st.median([r['ddG_metal_attributable'] for r in g]):+.3f} kcal/mol; "
            f"metal term fired {sum(1 for r in g if r['partcov_fired'])}/{len(g)}; "
            f"median replicate SD {st.median([r['ddG_withca_sd'] for r in g]):.3f}")

    log(f"\nwrote {out.relative_to(ROOT)}")
    (ROOT / f"logs/22_foldx_ablation_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                              encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
