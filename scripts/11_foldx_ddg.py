#!/usr/bin/env python3
"""
11_foldx_ddg.py — PHASE 5. FoldX ddG for every variant with structural coverage.

ddG is the predicted change in folding free energy on mutation (kcal/mol); positive means
destabilising. It is computed from coordinates alone and is independent of any clinical
classification, which is why it is one of the manuscript's headline measurements.

Mutation strings use the template's own residue numbering. That is normally the classic FBN1
trap — PDB author numbering is not UniProt numbering — but the packet PDBs were renumbered into
UniProt P35555 numbering in Phase 4, and the wild-type residue in every mutation string is
re-verified against the reference sequence before FoldX is invoked.

Templates are the four experimental structures. Each is repaired once (cached), then all of its
mutations are submitted in a single BuildModel run.

Writes:
  data/processed/foldx_ddg.tsv
  logs/11_foldx_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
WORK = Path("/tmp/claude-1000/-home-harvey-research-fbn1-marfan/"
            "604537d3-2ba9-4704-b90a-04b2c233decc/scratchpad/foldx")
MOLECULES = Path.home() / "tools/foldx/molecules/"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

TEMPLATES = [
    ("EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb", "2W86", (807, 951)),
    ("EXP_1LMJ_cbEGF12-13.pdb", "1LMJ", (1069, 1154)),
    ("EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb", "1UZJ", (1486, 1647)),
    ("EXP_1EMN_cbEGF32-33.pdb", "1EMN", (2127, 2205)),
]
CHAIN = "A"

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m, flush=True)
    log_lines.append(m)


def foldx(args: list[str], cwd: Path, timeout: int = 7200) -> subprocess.CompletedProcess:
    return subprocess.run(["foldx", *args, f"--rotabaseLocation={MOLECULES}"],
                          cwd=cwd, capture_output=True, text=True, timeout=timeout)


def main() -> int:
    if shutil.which("foldx") is None:
        raise SystemExit("STOP: foldx not on PATH — ask before omitting ddG (CLAUDE.md 5)")
    WORK.mkdir(parents=True, exist_ok=True)
    log(f"Phase 5 — FoldX ddG — {datetime.now(timezone.utc).isoformat()}")

    seq = json.loads(sorted((ROOT / "data/raw").glob("uniprot_P35555_v*.json"))[-1]
                     .read_text(encoding="utf-8"))["sequence"]["value"]
    variants = list(csv.DictReader(
        (ROOT / "data/processed/structural_metrics.tsv").open(encoding="utf-8"), delimiter="\t"))
    covered = [v for v in variants if v["structural_coverage"] == "True"]
    log(f"{len(covered):,} variants with structural coverage\n")

    results: dict[str, dict] = {}

    for fname, pdb, span in TEMPLATES:
        sub = [v for v in covered if v["template_pdb"] == pdb]
        if not sub:
            continue
        wd = WORK / pdb
        wd.mkdir(parents=True, exist_ok=True)
        shutil.copy(INBOX / fname, wd / "in.pdb")

        repaired = wd / "in_Repair.pdb"
        if not repaired.is_file():
            t0 = time.time()
            r = foldx(["--command=RepairPDB", "--pdb=in.pdb"], wd)
            if r.returncode != 0 or not repaired.is_file():
                log(f"  !! {pdb}: RepairPDB failed — {r.stdout[-300:]}")
                continue
            log(f"  {pdb}: RepairPDB {time.time() - t0:.0f}s")
        else:
            log(f"  {pdb}: using cached repair")

        # Re-verify every wild-type residue against P35555 before writing the mutation list.
        muts, keep, bad = [], [], []
        for v in sub:
            pos, wt, mut = int(v["position"]), v["wt_aa"], v["mut_aa"]
            if seq[pos - 1] != wt:
                bad.append((v["variation_id"], pos, wt, seq[pos - 1]))
                continue
            muts.append(f"{wt}{CHAIN}{pos}{mut};")
            keep.append(v)
        if bad:
            raise SystemExit(f"STOP: wild-type mismatch before FoldX: {bad[:5]}")

        (wd / "individual_list.txt").write_text("\n".join(muts) + "\n", encoding="utf-8")
        t0 = time.time()
        r = foldx(["--command=BuildModel", "--pdb=in_Repair.pdb",
                   "--mutant-file=individual_list.txt", "--numberOfRuns=1"], wd)
        dt = time.time() - t0

        dif = wd / "Dif_in_Repair.fxout"
        if r.returncode != 0 or not dif.is_file():
            log(f"  !! {pdb}: BuildModel failed (rc={r.returncode}) — {r.stdout[-300:]}")
            continue

        lines = [ln for ln in dif.read_text(encoding="utf-8").splitlines() if ln.startswith("in_Repair_")]
        hdr = None
        for ln in dif.read_text(encoding="utf-8").splitlines():
            if ln.startswith("Pdb\t"):
                hdr = ln.split("\t")
                break
        if hdr is None or len(lines) != len(keep):
            log(f"  !! {pdb}: expected {len(keep)} ddG rows, parsed {len(lines)} — recording as failed")
            for v in keep:
                results[v["variation_id"]] = {"foldx_status": "parse_failed"}
            continue

        idx = {c: i for i, c in enumerate(hdr)}
        for v, ln in zip(keep, lines):
            f = ln.split("\t")
            g = lambda c: float(f[idx[c]]) if c in idx and f[idx[c]] not in ("", None) else None
            results[v["variation_id"]] = {
                "foldx_status": "ok",
                "foldx_template": pdb,
                "ddG_kcal_mol": round(g("total energy"), 3),
                "ddG_vdw_clash": round(g("Van der Waals clashes") or 0, 3),
                "ddG_disulfide": round(g("disulfide") or 0, 3),
                "ddG_electrostatics": round(g("Electrostatics") or 0, 3),
                "ddG_solvation_polar": round(g("Solvation Polar") or 0, 3),
                "ddG_backbone_hbond": round(g("Backbone Hbond") or 0, 3),
                "ddG_sidechain_hbond": round(g("Sidechain Hbond") or 0, 3),
            }
        log(f"  {pdb}: {len(keep)} mutations in {dt:.0f}s ({dt/max(len(keep),1):.1f}s each)")

    # --- write ---------------------------------------------------------------
    rows = []
    for v in variants:
        r = results.get(v["variation_id"], {})
        rows.append({
            "variation_id": v["variation_id"], "position": v["position"],
            "wt_aa": v["wt_aa"], "mut_aa": v["mut_aa"],
            "set_primary": v["set_primary"], "set_sensitivity": v["set_sensitivity"],
            "site_label": v["site_label"], "cys_role": v["cys_role"],
            "cys_removing": v["cys_removing"], "is_ca_consensus": v["is_ca_consensus"],
            "structural_coverage": v["structural_coverage"],
            "foldx_status": r.get("foldx_status", "no_structure" if v["structural_coverage"] != "True" else "MISSING"),
            **{k: r.get(k, "") for k in ("foldx_template", "ddG_kcal_mol", "ddG_vdw_clash",
                                         "ddG_disulfide", "ddG_electrostatics",
                                         "ddG_solvation_polar", "ddG_backbone_hbond",
                                         "ddG_sidechain_hbond")},
        })

    cols = list(rows[0])
    with (ROOT / "data/processed/foldx_ddg.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    ok = [r for r in rows if r["foldx_status"] == "ok"]
    miss = [r for r in rows if r["foldx_status"] == "MISSING"]
    log(f"\n  ddG computed for {len(ok):,} variants; {len(miss)} covered-but-missing; "
        f"{sum(1 for r in rows if r['foldx_status'] == 'no_structure'):,} lack a structure")
    if ok:
        import statistics as st
        vals = [r["ddG_kcal_mol"] for r in ok]
        log(f"  ddG range {min(vals):.2f} .. {max(vals):.2f}, median {st.median(vals):.2f} kcal/mol")
    log("\nwrote data/processed/foldx_ddg.tsv")
    (ROOT / f"logs/11_foldx_{STAMP}.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
