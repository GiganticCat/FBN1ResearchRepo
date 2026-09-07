#!/usr/bin/env python3
"""
31_boltz_setup.py — PHASE 5 (revision). Local calcium cofolding with Boltz-2.

AlphaFold Server answered the mutant question for nine variants, one manual upload at a time.
That budget is why the batch-3 seed spread rests on five models and why 262 of the 271 side-chain
calcium-ligand variants have no mutant model at all. Boltz-2 runs the same class of calculation
on the local RTX 4060 under an MIT licence, with no quota and no hand-off, which turns both
limits into compute questions.

Two details make this cheap rather than a second infrastructure project.

  MSAs ARE ALREADY PAID FOR. Boltz needs an alignment, and generating one locally would mean
  hundreds of gigabytes of sequence databases. AlphaFold Server already returned one per
  construct -- the cbEGF32-33 unpaired alignment is 13,756 sequences deep -- so those files are
  reused directly and no MSA server is contacted.

  CALCIUM IS A LIGAND, NOT A GUESS. Boltz takes ions by CCD code, so a construct is specified as
  its sequence plus N copies of CCD `CA`, exactly as the AF3 jobs were.

Validation before use, not after. `--validate` folds a construct that already has both an AF3
model and an experimental structure, so Boltz's calcium placement can be compared against a known
answer in this domain family before any mutant result is believed. A tool that reproduces the
1EMN geometry has earned the right to be run on domains with no structure; one that does not is
reporting its own artifacts.

Usage:
  python scripts/31_boltz_setup.py --validate                  # one WT construct, compare to AF3
  python scripts/31_boltz_setup.py --write-inputs [--seeds N]  # YAMLs for WT + ligand mutants

Writes:
  structures/boltz/inputs/*.yaml
  structures/boltz/predictions/...        (created by the boltz CLI)
  logs/31_boltz_setup_<stamp>.log
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
SITES = ROOT / "data" / "processed" / "af3_ca_sites.tsv"
SM = ROOT / "data" / "processed" / "structural_metrics_af3.tsv"
FASTA = ROOT / "data" / "raw" / "uniprot_P35555_v264.fasta"
BOLTZ = ROOT / "structures" / "boltz"
INPUTS = BOLTZ / "inputs"
PREDS = BOLTZ / "predictions"
BOLTZ_BIN = ROOT / ".venv-boltz" / "bin" / "boltz"

VALIDATE_CONSTRUCT = "FBN1_cbEGF31-33_Ca3"   # contains cbEGF32-33, solved as 1EMN/1EMO

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"31_boltz_setup_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


def reference() -> str:
    return "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))


def job_dir(job_name: str) -> Path | None:
    want = re.sub(r"[^a-z0-9]+", "_", job_name.lower()).strip("_")
    for d in AFDIR.iterdir():
        if d.is_dir():
            have = re.sub(r"[^a-z0-9]+", "_", d.name.lower()).strip("_")
            if have in (want, f"fold_{want}"):
                return d
    return None


def msa_for(job_name: str) -> Path | None:
    """The AF3 unpaired alignment for a construct, reused verbatim."""
    d = job_dir(job_name)
    if d is None:
        return None
    hits = [p for p in (d / "msas").glob("*unpaired*.a3m")
            if not p.name.endswith("Zone.Identifier")]
    return hits[0] if hits else None


def write_yaml(path: Path, sequence: str, n_ions: int, msa: Path) -> None:
    """Boltz input: one protein chain plus N calcium ions given by CCD code."""
    ids = "BCDEFGHIJK"
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {sequence}", f"      msa: {msa}"]
    for i in range(n_ions):
        lines += ["  - ligand:", f"      id: {ids[i]}", "      ccd: CA"]
    path.write_text("\n".join(lines) + "\n")


def constructs() -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in csv.DictReader(open(SITES), delimiter="\t"):
        lo, hi = (int(x) for x in r["uniprot_range"].split("-"))
        d = out.setdefault(r["construct"], dict(lo=lo, hi=hi, n_ions=0,
                                                domains=r["cbegf_domains"]))
        d["n_ions"] += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--write-inputs", action="store_true")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--devices", type=int, default=1)
    ap.add_argument("--no-kernels", action="store_true",
                    help="disable the cuEquivariance triangular-multiplication kernels; "
                         "needed unless cuequivariance-torch is installed, and slower")
    args = ap.parse_args()

    log(f"Phase 5 revision — Boltz-2 setup — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    ref = reference()
    cons = constructs()
    log(f"  {len(cons)} constructs; reference {len(ref)} aa")
    INPUTS.mkdir(parents=True, exist_ok=True)
    PREDS.mkdir(parents=True, exist_ok=True)

    if args.validate:
        name = VALIDATE_CONSTRUCT
        c = cons[name]
        msa = msa_for(name)
        if msa is None:
            log(f"FAIL: no AF3 alignment found for {name}")
            return 1
        seq = ref[c["lo"] - 1:c["hi"]]
        y = INPUTS / f"{name}_wt.yaml"
        write_yaml(y, seq, c["n_ions"], msa)
        log(f"  validating on {name}: {len(seq)} aa, {c['n_ions']} Ca, "
            f"MSA {msa.name} ({sum(1 for l in msa.open() if l.startswith('>'))} seqs)")
        log(f"  wrote {y.relative_to(ROOT)}")

        cmd = [str(BOLTZ_BIN), "predict", str(y), "--out_dir", str(PREDS),
               "--diffusion_samples", "1", "--output_format", "mmcif",
               "--override"]
        if args.no_kernels:
            cmd.append("--no_kernels")
        log(f"  running: {' '.join(cmd)}")
        t = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        el = time.time() - t
        tail = (proc.stdout or "").strip().splitlines()[-12:]
        for line in tail:
            log(f"    {line}")
        if proc.returncode != 0:
            log(f"  FAILED after {el:.0f}s (exit {proc.returncode})")
            for line in (proc.stderr or "").strip().splitlines()[-15:]:
                log(f"    stderr: {line}")
            return 1
        log(f"  completed in {el:.0f}s")
        made = sorted(PREDS.rglob("*.cif"))
        log(f"  {len(made)} structure file(s) produced")
        for p in made[:5]:
            log(f"    {p.relative_to(ROOT)}")
        log()
        log("  Next: scripts/26-style site inventory on this model, compared against the AF3 "
            "model and 1EMN, before any mutant is folded.")
        return 0

    if args.write_inputs:
        rows = list(csv.DictReader(open(SM), delimiter="\t"))
        want = [r for r in rows if r["ca_ligand_sidechain"] == "True"]
        by_con = defaultdict(list)
        for r in want:
            by_con[r["construct"]].append(r)
        n = 0
        for name, c in sorted(cons.items()):
            msa = msa_for(name)
            if msa is None:
                log(f"  {name}: no alignment — skipped")
                continue
            seq = ref[c["lo"] - 1:c["hi"]]
            write_yaml(INPUTS / f"{name}_wt.yaml", seq, c["n_ions"], msa)
            n += 1
            for r in by_con.get(name, []):
                pos, wt, mut = int(r["position"]), r["wt_aa"], r["mut_aa"]
                i = pos - c["lo"]
                if not (0 <= i < len(seq)) or seq[i] != wt:
                    log(f"  {name}: {wt}{pos}{mut} does not match the construct — skipped")
                    continue
                mseq = seq[:i] + mut + seq[i + 1:]
                write_yaml(INPUTS / f"{name}_{wt}{pos}{mut}.yaml", mseq, c["n_ions"], msa)
                n += 1
        log(f"  wrote {n} Boltz inputs to {INPUTS.relative_to(ROOT)}")
        log(f"  run with: {BOLTZ_BIN} predict {INPUTS} --out_dir {PREDS} "
            f"--diffusion_samples {args.seeds} --output_format mmcif")
        return 0

    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
