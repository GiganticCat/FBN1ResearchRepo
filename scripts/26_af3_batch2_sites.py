#!/usr/bin/env python3
"""
26_af3_batch2_sites.py — PHASE 4/5. Validate AF3 batch 2 and inventory every calcium site.

Batch 2 is 21 tandem constructs covering all 43 cbEGF domains with three Ca2+ each. This script
turns those downloads into an analysis-ready site inventory, and refuses to do so quietly if the
models are not what was requested.

VALIDATION, before any measurement. For each construct: five models present; the modelled
sequence matches P35555 over the requested UniProt range exactly (identity, not just length);
the requested number of calcium ions is present. A construct failing any of these is reported
and excluded rather than silently carried forward -- a one-residue frame error here would
mis-assign every downstream ligand.

NUMBERING. AF3 numbers each construct from 1. The offset to UniProt is fixed per construct by
the manifest's uniprot_range and then *verified* by sequence identity, so every coordinating
residue below is reported in P35555 numbering and can be joined to the variant table directly.

CALIBRATION FOR FREE. Four batch-2 constructs (cbEGF9-11, 11-13, 21-23, 31-33) overlap domains
with experimental structures that were not used to build them. Their sites are flagged so the
placement accuracy measured on the batch-1 calibration set can be re-tested on constructs that
did not contribute to it.

Every quantity is a median over the five models, with the seed spread reported alongside --
the batch-3 analysis showed seed spread reaching 19% of site valence, which is larger than most
effects of interest, so a site's noise is part of its description here rather than an
afterthought.

Writes:
  data/processed/af3_ca_sites.tsv        one row per (construct, site)
  data/processed/af3_ca_ligands.tsv      one row per (construct, site, coordinating residue)
  results/phase5_af3_batch2.md
  logs/26_af3_batch2_sites_<stamp>.log
"""

from __future__ import annotations

import csv
import datetime as dt
import math
import re
import statistics as st
import sys
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
MANIFEST = ROOT / "handoff" / "structures" / "inbox" / "batch2_domains" / "manifest.csv"
FASTA = ROOT / "data" / "raw" / "uniprot_P35555_v264.fasta"
OUT_SITES = ROOT / "data" / "processed" / "af3_ca_sites.tsv"
OUT_LIG = ROOT / "data" / "processed" / "af3_ca_ligands.tsv"
OUT_MD = ROOT / "results" / "phase5_af3_batch2.md"

R0 = {"O": 1.967, "N": 2.14, "S": 2.45}
B = 0.37
CUTOFF = 3.2

# Constructs overlapping a domain with an experimental structure, usable as an out-of-sample
# check on AF3 calcium placement.
CALIBRATION_OVERLAP = {
    "FBN1_cbEGF9-11_Ca3": "2W86 (cbEGF9-hyb2-cbEGF10)",
    "FBN1_cbEGF11-13_Ca3": "1LMJ (cbEGF12-13)",
    "FBN1_cbEGF21-23_Ca3": "1UZJ (cbEGF22-TB4-cbEGF23)",
    "FBN1_cbEGF31-33_Ca3": "1EMN/1EMO (cbEGF32-33)",
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"26_af3_batch2_sites_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def read_reference() -> str:
    seq = []
    for line in FASTA.read_text().splitlines():
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def job_dir(job_name: str) -> Path | None:
    """AlphaFold Server lowercases job names and strips punctuation in the folder it returns."""
    want = re.sub(r"[^a-z0-9]+", "_", job_name.lower()).strip("_")
    for d in AFDIR.iterdir():
        if not d.is_dir():
            continue
        have = re.sub(r"[^a-z0-9]+", "_", d.name.lower()).strip("_")
        if have == want or have == f"fold_{want}":
            return d
    return None


def model_paths(job: Path) -> list[Path]:
    out = []
    for i in range(5):
        hits = [h for h in sorted(job.glob(f"*_model_{i}.cif"))
                if not h.name.endswith("Zone.Identifier")]
        if len(hits) == 1:
            out.append(hits[0])
    return out


def read_model(path: Path):
    stx = gemmi.read_structure(str(path))
    stx.setup_entities()
    mdl = stx[0]
    residues, ions, donors = [], [], []
    for ch in mdl:
        for res in ch:
            name = res.name.strip()
            if name == "CA" and len(res) == 1:
                ions.append(res[0].pos)
                continue
            if name not in AA3:
                continue
            residues.append((res.seqid.num, AA3[name], res))
            for at in res:
                el = at.element.name.upper()
                if el in R0:
                    donors.append((res.seqid.num, at.name, el, at.pos))
    residues.sort(key=lambda r: r[0])
    return residues, ions, donors


def shell_of(donors, ion_pos):
    out = []
    for seqid, aname, el, pos in donors:
        d = ion_pos.dist(pos)
        if d < CUTOFF:
            out.append((seqid, aname, el, d))
    return sorted(out, key=lambda x: x[3])


def bvs_of(shell) -> float:
    return sum(math.exp((R0[el] - d) / B) for _, _, el, d in shell)


def med(xs):
    return st.median(xs) if xs else float("nan")


def spread(xs):
    if len(xs) < 2:
        return float("nan")
    m = st.median(xs)
    return 1.4826 * st.median([abs(x - m) for x in xs])


def is_backbone(atom_name: str) -> bool:
    return atom_name.strip() in {"O", "N", "OXT"}


def main() -> int:
    log(f"Phase 4/5 — AF3 batch 2 site inventory — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  {CUTOFF} A shell; R0(Ca-O)={R0['O']}, R0(Ca-N)={R0['N']}, b={B}")
    log()

    ref = read_reference()
    log(f"reference P35555: {len(ref)} aa")
    if len(ref) != 2871:
        log(f"  *** expected 2871 aa, got {len(ref)} — refusing to map ***")
        return 1

    rows = csv.DictReader(open(MANIFEST))
    site_rows, lig_rows = [], []
    n_ok = n_fail = 0

    log()
    log("Validation and site inventory")
    for m in rows:
        job_name = m["job_name"]
        lo, hi = (int(x) for x in m["uniprot_range"].split("-"))
        want_len = int(m["length_aa"])
        want_ions = int(m["n_ca_ions"])
        cbegfs = m["cbegf_domains"].strip('"')

        d = job_dir(job_name)
        if d is None:
            log(f"  {job_name:26s} MISSING — no downloaded folder matches")
            n_fail += 1
            continue

        paths = model_paths(d)
        if len(paths) != 5:
            log(f"  {job_name:26s} FAIL — {len(paths)} models found, expected 5")
            n_fail += 1
            continue

        # ---- validate model 0 against the reference before measuring anything
        residues, ions, donors = read_model(paths[0])
        seq = "".join(r[1] for r in residues)
        expect = ref[lo - 1:hi]
        if len(seq) != want_len:
            log(f"  {job_name:26s} FAIL — modelled {len(seq)} aa, manifest says {want_len}")
            n_fail += 1
            continue
        if seq != expect:
            nmis = sum(1 for a, b in zip(seq, expect) if a != b)
            log(f"  {job_name:26s} FAIL — sequence differs from P35555 {lo}-{hi} "
                f"at {nmis} position(s)")
            n_fail += 1
            continue
        if len(ions) != want_ions:
            log(f"  {job_name:26s} FAIL — {len(ions)} Ca ions, expected {want_ions}")
            n_fail += 1
            continue

        offset = lo - residues[0][0]   # local seqid + offset = UniProt position

        # ---- reference sites from model 0, then follow them across the other four
        ref_shells = [shell_of(donors, p) for p in ions]
        ref_sigs = [{s for s, _, _, _ in sh} for sh in ref_shells]

        per_site: dict[int, list[dict]] = {i: [] for i in range(len(ions))}
        for mi, p in enumerate(paths):
            res_i, ions_i, don_i = read_model(p)
            shells_i = [shell_of(don_i, q) for q in ions_i]
            used = set()
            for k, sh in enumerate(shells_i):
                sig = {s for s, _, _, _ in sh}
                best, best_j = 0.0, None
                for j, rs in enumerate(ref_sigs):
                    if j in used or not (sig | rs):
                        continue
                    jac = len(sig & rs) / len(sig | rs)
                    if jac > best:
                        best, best_j = jac, j
                if best_j is not None and best >= 0.35:
                    used.add(best_j)
                    per_site[best_j].append(dict(model=mi, shell=sh, bvs=bvs_of(sh), cn=len(sh)))

        plddt = med([float(np.mean([a.b_iso for _, _, r in read_model(p)[0] for a in r]))
                     for p in paths])

        for j in range(len(ions)):
            obs = per_site[j]
            b = [o["bvs"] for o in obs]
            cn = [o["cn"] for o in obs]
            sd = spread(b)
            # ligand set: residues coordinating in at least 3 of the models it was found in
            counts: dict[tuple[int, str], int] = {}
            for o in obs:
                for s, aname, el, dist in o["shell"]:
                    counts[(s, aname)] = counts.get((s, aname), 0) + 1
            consensus = {k: v for k, v in counts.items() if v >= max(3, len(obs) - 1)}

            uni_res = sorted({s + offset for (s, _) in consensus})
            site_rows.append(dict(
                construct=job_name, site=j, cbegf_domains=cbegfs,
                uniprot_range=m["uniprot_range"],
                n_models_present=len(obs), plddt=round(plddt, 1),
                bvs_median=round(med(b), 3) if b else "",
                bvs_seed_sd=round(sd, 3) if not math.isnan(sd) else "",
                bvs_seed_sd_pct=round(100 * sd / med(b), 1) if b and med(b) and not math.isnan(sd) else "",
                cn_median=round(med(cn), 1) if cn else "",
                n_ligand_residues=len(uni_res),
                ligand_positions_uniprot=";".join(str(x) for x in uni_res),
                calibration_overlap=CALIBRATION_OVERLAP.get(job_name, ""),
            ))
            for (s, aname), c in sorted(consensus.items()):
                up = s + offset
                wt = ref[up - 1]
                dists = [dd for o in obs for (ss, an, _, dd) in o["shell"]
                         if ss == s and an == aname]
                lig_rows.append(dict(
                    construct=job_name, site=j, cbegf_domains=cbegfs,
                    uniprot_pos=up, wt_aa=wt, atom=aname,
                    atom_class="backbone" if is_backbone(aname) else "sidechain",
                    n_models=c, dist_median_A=round(med(dists), 2),
                ))

        flag = f"  [calibration overlap: {CALIBRATION_OVERLAP[job_name]}]" if job_name in CALIBRATION_OVERLAP else ""
        log(f"  {job_name:26s} OK  {want_len:3d} aa  cbEGF {cbegfs:8s} "
            f"pLDDT {plddt:.1f}  sites {len(ions)}{flag}")
        n_ok += 1

    log()
    log(f"{n_ok} constructs validated, {n_fail} failed")
    if n_ok == 0:
        log("FAIL: nothing to write")
        return 1

    OUT_SITES.parent.mkdir(parents=True, exist_ok=True)
    for path, data in ((OUT_SITES, site_rows), (OUT_LIG, lig_rows)):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(data)
        log(f"wrote {path.relative_to(ROOT)}  ({len(data)} rows)")

    # ---- summary
    sds = [r["bvs_seed_sd_pct"] for r in site_rows if r["bvs_seed_sd_pct"] != ""]
    bvss = [r["bvs_median"] for r in site_rows if r["bvs_median"] != ""]
    sidechain = sum(1 for r in lig_rows if r["atom_class"] == "sidechain")
    log()
    log(f"  {len(site_rows)} calcium sites across {n_ok} constructs")
    log(f"  BVS median {med(bvss):.2f} (range {min(bvss):.2f}-{max(bvss):.2f})")
    log(f"  seed spread median {med(sds):.1f}% of valence "
        f"(range {min(sds):.1f}-{max(sds):.1f}%)")
    log(f"  {len(lig_rows)} coordinating atoms; {sidechain} side-chain, "
        f"{len(lig_rows)-sidechain} backbone")
    log(f"  {len({r['uniprot_pos'] for r in lig_rows})} distinct P35555 positions coordinate a calcium")

    with open(OUT_MD, "w") as fh:
        fh.write("# AF3 batch 2 — calcium sites across all cbEGF domains\n\n")
        fh.write(f"Generated by `scripts/26_af3_batch2_sites.py` at "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}.\n\n")
        fh.write(f"{n_ok} of {n_ok+n_fail} constructs passed validation (five models, sequence "
                 f"identical to P35555 over the requested range, requested ion count present). "
                 f"{len(site_rows)} calcium sites; "
                 f"{len({r['uniprot_pos'] for r in lig_rows})} distinct P35555 positions "
                 f"coordinate a calcium in at least one model.\n\n")
        fh.write("Values are medians over five AlphaFold Server models. `seed SD` is the spread "
                 "across those five and is reported for every site, because it reaches a size "
                 "comparable to the mutant effects this project measures.\n\n")
        fh.write("| construct | site | cbEGF | pLDDT | BVS | seed SD | CN | ligands | calibration |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in site_rows:
            fh.write(f"| {r['construct']} | {r['site']} | {r['cbegf_domains']} | {r['plddt']} | "
                     f"{r['bvs_median']} | {r['bvs_seed_sd_pct']}% | {r['cn_median']} | "
                     f"{r['n_ligand_residues']} | {r['calibration_overlap'] or '—'} |\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
