#!/usr/bin/env python3
"""
24_allmetal3d.py — PHASE 5 (revision). A second, metal-aware opinion that owes nothing to FoldX.

AllMetal3D (Durr et al., LCBC-EPFL) is a 3D convolutional network that predicts where metal ions
bind in a protein, which metal it is, and the coordination geometry — from the protein alone, with
no ion in the input. That makes it the right instrument for the question FoldX cannot answer: given
the mutant protein, does a calcium site still exist there?

Two stages.

  CALIBRATION. Run on the four experimental templates with their ions stripped, and measure how
  close the predicted metal comes to the crystallographic one. On 1EMN this recovers both sites at
  0.49 A / 1.36 A with probability 1.00 / 0.96. Without this step a prediction on a mutant would be
  uninterpretable; with it we know what "found the site" looks like in this domain family.

  MUTANTS. Run on FoldX-built mutant models (ions stripped) and record, for each reference calcium
  position, the best site probability nearby and how far the prediction moved. The comparison that
  matters is between variants that delete a side-chain calcium donor and variants at positions
  that coordinate the ion only through their backbone carbonyl, which no substitution can remove.

Caveat recorded rather than hidden: FoldX mutant models are locally repacked, not relaxed with a
metal present, so a lost site here means "the local environment no longer looks like a calcium
site", not "the mutant protein cannot bind calcium". AF3 mutant cofolding (Batch 3) tests the
stronger version separately.

Requires the separate Python 3.11 environment `.venv-metal3d` (AllMetal3D pins numpy/scipy versions
that do not exist for the project interpreter). Recorded in manifest/tools.md.

Writes:
  data/processed/allmetal3d_sites.tsv
  results/phase5_allmetal3d_calibration.md
  logs/24_allmetal3d_<stamp>.log
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import shutil
import statistics as st
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
FOLDX_WORK = Path(os.environ.get("FBN1_FOLDX_WORK", "/tmp/fbn1-foldx-ablation"))
WORK = Path(os.environ.get("FBN1_AM3D_WORK", "/tmp/fbn1-allmetal3d"))
AM3D = ROOT / ".venv-metal3d/bin/allmetal3d"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

TEMPLATES = {"2W86": "EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb",
             "1LMJ": "EXP_1LMJ_cbEGF12-13.pdb",
             "1UZJ": "EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb",
             "1EMN": "EXP_1EMN_cbEGF32-33.pdb"}
N_CONTROLS = 20          # per control class, sampled with a fixed seed
SEED = 20260809
NEAR = 3.0               # A; a prediction this close counts as "the same site"

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m, flush=True)
    log_lines.append(m)


def strip_hetatm(src: Path, dst: Path) -> None:
    dst.write_text("\n".join(l for l in src.read_text(encoding="utf-8").splitlines()
                             if not l.startswith("HETATM")) + "\n", encoding="utf-8")


def ref_calciums(template: str) -> list[tuple[str, tuple[float, float, float]]]:
    out = []
    for l in (INBOX / TEMPLATES[template]).read_text(encoding="utf-8").splitlines():
        if l.startswith("HETATM") and l[17:20].strip() == "CA":
            out.append((l[22:26].strip(),
                        (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def run_am3d(pdb: Path, tag: str) -> list[tuple[tuple[float, float, float], float]] | None:
    """Predicted metal positions with their probabilities. Idempotent: cached by tag."""
    out = WORK / tag
    metals = out / f"{pdb.stem}_metals.pdb"
    if not metals.is_file():
        out.mkdir(parents=True, exist_ok=True)
        r = subprocess.run([str(AM3D), "-i", pdb.name, "-a", "allmetal3d", "-m", "fast",
                            "-o", str(out.resolve())],
                           cwd=pdb.parent, capture_output=True, text=True, timeout=3600)
        if r.returncode != 0 or not metals.is_file():
            log(f"    !! AllMetal3D failed for {tag}: {r.stderr.strip()[-200:]}")
            return None
    return [((float(l[30:38]), float(l[38:46]), float(l[46:54])), float(l[54:60]))
            for l in metals.read_text(encoding="utf-8").splitlines() if l.startswith("HETATM")]


def score(preds, refs) -> list[dict]:
    rows = []
    for rid, xyz in refs:
        if preds:
            d, p = min(((math.dist(xyz, q), pr) for q, pr in preds), key=lambda t: t[0])
        else:
            d, p = float("inf"), 0.0
        rows.append({"site": rid, "dist_A": round(d, 2) if d != float("inf") else "",
                     "probability": round(p, 2), "site_found": d <= NEAR})
    return rows


def foldx_models(pdb_code: str) -> dict[str, Path]:
    """Map 'WT<pos><MUT>' to the FoldX mutant model built with the ion present."""
    wd = FOLDX_WORK / pdb_code / "withca"
    lst = wd / "individual_list.txt"
    if not lst.is_file():
        return {}
    out = {}
    for i, line in enumerate((l for l in lst.read_text(encoding="utf-8").splitlines() if l.strip()),
                             start=1):
        m = line.strip().rstrip(";")           # e.g. DA2127A -> chain letter at index 1
        key = m[0] + m[2:]
        model = wd / f"in_Repair_{i}_0.pdb"    # replicate 0 = the deterministic first build
        if model.is_file():
            out[key] = model
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calibrate-only", action="store_true",
                    help="run stage 1 only; skip the mutant panel, which needs FoldX models")
    args = ap.parse_args()

    if not AM3D.is_file():
        raise SystemExit(f"STOP: {AM3D} not found — create .venv-metal3d first (see docstring)")
    WORK.mkdir(parents=True, exist_ok=True)
    log(f"Phase 5 revision — AllMetal3D — {datetime.now(timezone.utc).isoformat()}\n")

    # ---------------- stage 1: calibration on the experimental templates ----------------
    log("Calibration — ions stripped from the experimental structures")
    cal_rows = []
    for code, fname in TEMPLATES.items():
        apo = WORK / f"{code}_apo.pdb"
        strip_hetatm(INBOX / fname, apo)
        preds = run_am3d(apo, f"cal_{code}")
        if preds is None:
            continue
        for r in score(preds, ref_calciums(code)):
            cal_rows.append({"template": code, **r})
            log(f"  {code} Ca {r['site']:>6}: nearest prediction {r['dist_A']} A, "
                f"p = {r['probability']:.2f}")
    if not cal_rows:
        raise SystemExit("STOP: calibration produced no results")
    found = [r for r in cal_rows if r["site_found"]]
    med = st.median([r["dist_A"] for r in found]) if found else float("nan")
    log(f"  {len(found)}/{len(cal_rows)} crystallographic sites recovered within {NEAR} A; "
        f"median displacement {med:.2f} A\n")

    (ROOT / "results/phase5_allmetal3d_calibration.md").write_text(
        "# AllMetal3D calibration on the experimental cbEGF structures\n\n"
        f"Generated by `scripts/24_allmetal3d.py` at {datetime.now(timezone.utc).isoformat()}.\n\n"
        "Each experimental structure was stripped of its calcium ions and handed to AllMetal3D, "
        "which predicts metal sites from the protein alone. The question is whether the network "
        "puts the ion back where the crystallographer found it.\n\n"
        f"**{len(found)} of {len(cal_rows)} sites recovered within {NEAR} A; median displacement "
        f"{med:.2f} A.**\n\n"
        "| template | site | displacement (A) | probability | recovered |\n|---|---|---|---|---|\n"
        + "\n".join(f"| {r['template']} | {r['site']} | {r['dist_A']} | {r['probability']:.2f} | "
                    f"{'yes' if r['site_found'] else 'no'} |" for r in cal_rows)
        + "\n\nThis is what a surviving calcium site looks like in this domain family, and it is "
          "the baseline every mutant model is scored against.\n", encoding="utf-8")

    if args.calibrate_only:
        log("  --calibrate-only: stopping before the mutant panel.")
        log("  Stage 2 as written scores FoldX-repacked mutant models, which are locally repacked "
            "rather than re-folded with the metal present. Those models lived in "
            f"{FOLDX_WORK} and did not survive the reboot. They should be replaced by the Boltz "
            "mutant models rather than regenerated — same question, stronger structure.")
        return 0

    # ---------------- stage 2: mutant models ----------------
    mg = {r["variation_id"]: r for r in csv.DictReader(
        (ROOT / "data/processed/metal_geometry.tsv").open(encoding="utf-8"), delimiter="\t")}
    sm = {r["variation_id"]: r for r in csv.DictReader(
        (ROOT / "data/processed/structural_metrics.tsv").open(encoding="utf-8"), delimiter="\t")}

    ligand = [v for v, r in mg.items() if r["metal_status"] == "ok"]
    others = [v for v, r in mg.items()
              if r["metal_status"] == "not_a_ligand" and sm[v]["cys_removing"] != "True"]
    cys = [v for v, r in mg.items()
           if r["metal_status"] == "not_a_ligand" and sm[v]["cys_removing"] == "True"]
    rng = random.Random(SEED)
    panel = sorted(ligand) + sorted(rng.sample(sorted(others), min(N_CONTROLS, len(others)))) \
        + sorted(rng.sample(sorted(cys), min(N_CONTROLS, len(cys))))
    log(f"Mutant panel: {len(ligand)} ligand-shell + {min(N_CONTROLS, len(others))} non-ligand "
        f"+ {min(N_CONTROLS, len(cys))} cysteine-removing = {len(panel)} models")

    models = {code: foldx_models(code) for code in TEMPLATES}
    n_avail = sum(len(v) for v in models.values())
    log(f"  {n_avail} FoldX mutant models available in {FOLDX_WORK}")
    if n_avail == 0:
        raise SystemExit(f"STOP: no FoldX mutant models under {FOLDX_WORK} — "
                         "run scripts/22_foldx_ion_ablation.py first")

    rows, missing = [], 0
    for i, vid in enumerate(panel, 1):
        v, g = sm[vid], mg[vid]
        code = v["template_pdb"]
        key = f"{v['wt_aa']}{v['position']}{v['mut_aa']}"
        model = models.get(code, {}).get(key)
        if model is None:
            missing += 1
            rows.append({"variation_id": vid, "mutation": key, "template": code,
                         "am3d_status": "no_model", **{k: "" for k in SITE_COLS}})
            continue
        apo = WORK / f"{code}_{key}_apo.pdb"
        strip_hetatm(model, apo)
        preds = run_am3d(apo, f"mut_{code}_{key}")
        if preds is None:
            rows.append({"variation_id": vid, "mutation": key, "template": code,
                         "am3d_status": "failed", **{k: "" for k in SITE_COLS}})
            continue
        # Score the calcium site this variant actually touches; for controls, the nearest one.
        best = None
        for r in score(preds, ref_calciums(code)):
            if best is None or r["dist_A"] == "" or (best["dist_A"] != "" and
                                                     r["dist_A"] < best["dist_A"]):
                best = r
        rows.append({"variation_id": vid, "mutation": key, "template": code,
                     "am3d_status": "ok",
                     "class": ("sidechain_donor_lost" if g["ligand_source"] in
                               ("sidechain", "sidechain+backbone")
                               else "backbone_only" if g["ligand_source"] == "backbone"
                               else "cys_removing" if v["cys_removing"] == "True" else "control"),
                     "site": best["site"], "displacement_A": best["dist_A"],
                     "probability": best["probability"], "site_found": best["site_found"],
                     "delta_bvs_pct": g.get("delta_bvs_pct", ""),
                     "set_primary": v["set_primary"],
                     "am_pathogenicity": v["am_pathogenicity"]})
        if i % 20 == 0:
            log(f"  {i}/{len(panel)} models scored")

    out = ROOT / "data/processed/allmetal3d_sites.tsv"
    cols = ["variation_id", "mutation", "template", "am3d_status", *SITE_COLS]
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows([{c: r.get(c, "") for c in cols} for r in rows])

    ok = [r for r in rows if r["am3d_status"] == "ok"]
    log(f"\n  {len(ok)} scored, {missing} without a FoldX model")
    by = {}
    for r in ok:
        by.setdefault(r["class"], []).append(r)
    log(f"\n  {'class':22}{'n':>4}{'site kept':>11}{'med displ':>11}{'med p':>8}")
    for k, g in sorted(by.items()):
        kept = sum(1 for r in g if r["site_found"])
        dd = [r["displacement_A"] for r in g if r["displacement_A"] != ""]
        log(f"  {k:22}{len(g):4}{100 * kept / len(g):10.0f}%"
            f"{st.median(dd) if dd else float('nan'):11.2f}"
            f"{st.median([r['probability'] for r in g]):8.2f}")
    log(f"\nwrote {out.relative_to(ROOT)}")
    (ROOT / f"logs/24_allmetal3d_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                          encoding="utf-8")
    return 0


SITE_COLS = ["class", "site", "displacement_A", "probability", "site_found", "delta_bvs_pct",
             "set_primary", "am_pathogenicity"]

if __name__ == "__main__":
    sys.exit(main())
