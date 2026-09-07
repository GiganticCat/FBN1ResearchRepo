#!/usr/bin/env python3
"""
25_af3_mutant_geometry.py — PHASE 5 (revision). AF3 batch-3 mutants vs their wild-type folds.

Answers the two halves of the thesis with measurements that do not pass through FoldX:

  FOLD INTACT?   Ca RMSD and pLDDT, each mutant model against each wild-type model.
  SITE INTACT?   Bond-valence sum on the AF3-placed calcium, mutant against wild type.

Three things distinguish this from the exploratory probe it replaces.

  ALL FIVE SEEDS. AlphaFold Server returns five models per job. Using model_0 alone gives a
  number with no error bar, and the differences at stake here (a few percent of site valence)
  are the size of seed noise. Every quantity below is a median over five models, and the
  wild-type seed spread is what a mutant effect has to clear to count.

  IONS MATCHED BY COORDINATION SHELL, NOT BY DISTANCE. A cbEGF tandem has one calcium per
  domain, and the naive "nearest ion after superposition" match silently pairs a wild-type site
  with the wrong mutant ion when a site is genuinely vacated -- which is exactly the case the
  study cares about. Sites are identified by which residues coordinate them (Jaccard overlap
  against the wild-type reference shells, assigned with the Hungarian algorithm), so a vacated
  site is reported as vacated instead of being matched to a spurious partner 30 A away.

  BVS NEEDS NO SUPERPOSITION. Valence is intrinsic to each structure's own coordinates.
  Superposition enters only for RMSD and for reporting ion displacement.

Caveat recorded rather than hidden: AF3 is a structure predictor trained largely on holo
structures, so it is biased toward placing an ion somewhere. The benign control (P1141L) and the
disulfide control (C1138F) define what "no calcium effect" looks like under that bias; nothing
here should be read without them.

Writes:
  data/processed/af3_mutant_geometry.tsv
  results/phase5_af3_mutants.md
  logs/25_af3_mutant_geometry_<stamp>.log
"""

from __future__ import annotations

import csv
import datetime as dt
import math
import os
import re
import statistics as st
import sys
from pathlib import Path

import gemmi
import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
OUT_TSV = ROOT / "data" / "processed" / "af3_mutant_geometry.tsv"
OUT_MD = ROOT / "results" / "phase5_af3_mutants.md"

# Brown-Altermatt, identical parameters to scripts/23_metal_geometry.py so the two are comparable.
R0 = {"O": 1.967, "N": 2.14, "S": 2.45}
B = 0.37
CUTOFF = 3.2

# A shell must overlap a reference shell by at least this Jaccard index to be called the same
# site. Losing one ligand of seven drops Jaccard to ~0.86; losing three drops it to ~0.57. Below
# 0.35 the ion is somewhere else entirely.
JACCARD_MIN = 0.35

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"25_af3_mutant_geometry_{STAMP}.log"
_log_fh = None


def log(msg: str = "") -> None:
    global _log_fh
    if _log_fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _log_fh = open(LOGPATH, "w")
    print(msg)
    _log_fh.write(msg + "\n")
    _log_fh.flush()


# ---------------------------------------------------------------- structure reading


def model_paths(job: Path) -> list[Path]:
    """The five AF3 models for a job, ordered by model index."""
    out = []
    for i in range(5):
        hits = sorted(job.glob(f"*_model_{i}.cif"))
        hits = [h for h in hits if not h.name.endswith("Zone.Identifier")]
        if len(hits) == 1:
            out.append(hits[0])
    return out


def read_model(path: Path):
    """Return (residues, ions, donor_atoms) for one CIF.

    residues    (chain, seqid, one-letter, gemmi residue)
    ions        (chain, seqid, gemmi.Position)
    donors      (seqid, atom_name, element, gemmi.Position) for every protein O/N/S
    """
    st_ = gemmi.read_structure(str(path))
    st_.setup_entities()
    mdl = st_[0]
    residues, ions, donors = [], [], []
    for ch in mdl:
        for res in ch:
            name = res.name.strip()
            # AF3 writes a bare calcium as a single-atom residue named CA. A protein residue
            # also contains an atom named CA, hence the length test rather than a name test.
            if name == "CA" and len(res) == 1:
                ions.append((ch.name, res.seqid.num, res[0].pos))
                continue
            if name not in AA3:
                continue
            residues.append((ch.name, res.seqid.num, AA3[name], res))
            for at in res:
                el = at.element.name.upper()
                if el in R0:
                    donors.append((res.seqid.num, at.name, el, at.pos))
    return residues, ions, donors


def shell_of(donors, ion_pos) -> list[tuple[int, str, str, float]]:
    out = []
    for seqid, aname, el, pos in donors:
        d = ion_pos.dist(pos)
        if d < CUTOFF:
            out.append((seqid, aname, el, d))
    return sorted(out, key=lambda x: x[3])


def bvs_of(shell) -> float:
    return sum(math.exp((R0[el] - d) / B) for _, _, el, d in shell)


def sequence(residues) -> str:
    return "".join(r[2] for r in residues)


def calpha_array(residues) -> np.ndarray:
    pts = []
    for _, _, _, res in residues:
        a = res.find_atom("CA", "*")
        if a is None:
            return np.empty((0, 3))
        pts.append([a.pos.x, a.pos.y, a.pos.z])
    return np.asarray(pts)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    """Superpose P onto Q. Returns (rmsd, rotation, translation)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Qc - R @ Pc
    Pt = (R @ P.T).T + t
    rmsd = float(np.sqrt(((Pt - Q) ** 2).sum(1).mean()))
    return rmsd, R, t


# ---------------------------------------------------------------- site identity


def site_signature(shell) -> set[int]:
    """Residue numbers coordinating this ion. Identity of a site, robust to it moving."""
    return {seqid for seqid, _, _, _ in shell}


def jaccard(a: set[int], b: set[int]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def assign_sites(shells, ref_sigs: list[set[int]]) -> list[int | None]:
    """Map each ion in `shells` to a reference site index, or None if it matches nothing.

    Hungarian assignment on 1 - Jaccard, so two ions can never claim the same site.
    """
    n_ion, n_ref = len(shells), len(ref_sigs)
    if n_ion == 0 or n_ref == 0:
        return [None] * n_ion
    cost = np.ones((n_ion, n_ref))
    for i, sh in enumerate(shells):
        sig = site_signature(sh)
        for j, ref in enumerate(ref_sigs):
            cost[i, j] = 1.0 - jaccard(sig, ref)
    rows, cols = linear_sum_assignment(cost)
    out: list[int | None] = [None] * n_ion
    for i, j in zip(rows, cols):
        if 1.0 - cost[i, j] >= JACCARD_MIN:
            out[i] = int(j)
    return out


# ---------------------------------------------------------------- per-job aggregation


def job_sites(job: Path, ref_sigs: list[set[int]] | None):
    """Measure every model of a job. Returns (per_site, ref_sigs, seq, plddt_list, residues0).

    per_site maps reference-site index -> list of dicts, one per model in which that site was
    found. A model in which the site was not found contributes nothing, and the shortfall is
    what `n_models_present` reports.
    """
    paths = model_paths(job)
    if not paths:
        return None
    per_site: dict[int, list[dict]] = {}
    plddts, seqs, residues0, vacated = [], [], None, 0

    for mi, p in enumerate(paths):
        residues, ions, donors = read_model(p)
        if residues0 is None:
            residues0 = residues
        seqs.append(sequence(residues))
        plddts.append(float(np.mean([at.b_iso for _, _, _, res in residues for at in res])))

        shells = [shell_of(donors, pos) for _, _, pos in ions]

        if ref_sigs is None:
            # First model of the wild type defines the reference sites.
            ref_sigs = [site_signature(sh) for sh in shells]

        assigned = assign_sites(shells, ref_sigs)
        for k, site in enumerate(assigned):
            if site is None:
                vacated += 1
                continue
            per_site.setdefault(site, []).append(
                dict(
                    model=mi,
                    bvs=bvs_of(shells[k]),
                    cn=len(shells[k]),
                    ion_pos=ions[k][2],
                    shell=shells[k],
                    jac=jaccard(site_signature(shells[k]), ref_sigs[site]),
                )
            )
    return per_site, ref_sigs, seqs, plddts, residues0, vacated, paths


def med(xs):
    return st.median(xs) if xs else float("nan")


def spread(xs):
    """Median absolute deviation scaled to a standard-deviation equivalent."""
    if len(xs) < 2:
        return float("nan")
    m = st.median(xs)
    return 1.4826 * st.median([abs(x - m) for x in xs])


def main() -> int:
    log(f"Phase 5 revision — AF3 mutant geometry — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  {CUTOFF} A shell; R0(Ca-O)={R0['O']}, R0(Ca-N)={R0['N']}, b={B}; "
        f"site match Jaccard >= {JACCARD_MIN}")
    log()

    if not AFDIR.is_dir():
        log(f"FAIL: {AFDIR} does not exist")
        return 1

    jobs = sorted(d for d in AFDIR.iterdir() if d.is_dir())
    wt_jobs = {d.name: d for d in jobs if d.name.startswith("fold_")}
    mut_jobs = [d for d in jobs if not d.name.startswith("fold_")]

    if not wt_jobs:
        log("FAIL: no wild-type (fold_*) jobs found")
        return 1
    log(f"{len(wt_jobs)} wild-type jobs, {len(mut_jobs)} mutant jobs")
    log()

    # ---- wild types first: they define reference sites and the seed-noise floor
    wt_data = {}
    log("Wild-type seed spread (the noise a mutant effect must clear)")
    for name, path in sorted(wt_jobs.items()):
        r = job_sites(path, None)
        if r is None:
            log(f"  {name}: NO MODELS FOUND — skipped")
            continue
        per_site, ref_sigs, seqs, plddts, residues0, vacated, paths = r
        base = re.sub(r"^fold_", "", re.sub(r"_ca2$", "", name))
        wt_data[base] = dict(per_site=per_site, ref_sigs=ref_sigs, seq=seqs[0],
                             plddt=med(plddts), residues=residues0, n_models=len(paths))
        for site in sorted(per_site):
            b = [m["bvs"] for m in per_site[site]]
            log(f"  {base:34s} site {site}  n={len(b)}/{len(paths)} models  "
                f"BVS median {med(b):.3f}  seed SD {spread(b):.3f}  "
                f"({100*spread(b)/med(b):.1f}% of valence)")
        if vacated:
            log(f"  {base:34s} {vacated} ion-model(s) matched no reference site")
    log()

    rows = []
    log("Mutants")
    for path in mut_jobs:
        name = path.name
        m = re.match(r"^(.*?)_([a-z]\d+[a-z])_ca2$", name)
        if not m:
            log(f"  {name}: cannot parse a mutation from the job name — skipped")
            continue
        base, muttok = m.group(1), m.group(2).upper()
        if base not in wt_data:
            log(f"  {name}: no wild-type counterpart for '{base}' — skipped")
            continue
        wt = wt_data[base]

        r = job_sites(path, wt["ref_sigs"])
        if r is None:
            log(f"  {name}: NO MODELS FOUND — skipped")
            continue
        per_site, _, seqs, plddts, residues0, vacated, paths = r

        # verify the mutation actually present, by sequence difference against the wild type
        if len(seqs[0]) != len(wt["seq"]):
            log(f"  {name}: LENGTH MISMATCH vs wild type — skipped")
            continue
        diffs = [i for i, (a, b) in enumerate(zip(wt["seq"], seqs[0])) if a != b]
        if len(diffs) != 1:
            log(f"  {name}: expected exactly 1 sequence difference, found {len(diffs)} — skipped")
            continue
        i = diffs[0]
        local_num = wt["residues"][i][1]
        observed = f"{wt['seq'][i]}{local_num}{seqs[0][i]}"
        expected_wt, expected_mut = muttok[0], muttok[-1]
        ident_ok = (wt["seq"][i] == expected_wt and seqs[0][i] == expected_mut)

        # fold: every mutant model against wild-type model 0
        wt_ca = calpha_array(wt["residues"])
        rmsds = []
        for p in paths:
            res_m, _, _ = read_model(p)
            mu_ca = calpha_array(res_m)
            if mu_ca.shape == wt_ca.shape and mu_ca.size:
                rmsds.append(kabsch_rmsd(mu_ca, wt_ca)[0])

        log(f"  {name}")
        log(f"      mutation {observed} (job name says {muttok}) "
            f"{'OK' if ident_ok else '*** MISMATCH ***'}")
        log(f"      fold: Ca RMSD median {med(rmsds):.2f} A over {len(rmsds)} models; "
            f"pLDDT {med(plddts):.1f} vs wild type {wt['plddt']:.1f}")

        for site in sorted(wt["per_site"]):
            wt_b = [x["bvs"] for x in wt["per_site"][site]]
            wt_cn = [x["cn"] for x in wt["per_site"][site]]
            noise = spread(wt_b)
            present = per_site.get(site, [])
            n_present = len(present)

            if n_present == 0:
                log(f"      site {site}: VACATED in all {len(paths)} models "
                    f"(no ion matched the wild-type shell)")
                rows.append(dict(
                    job=name, base=base, mutation=observed, mutation_expected=muttok,
                    identity_ok=ident_ok, site=site, status="vacated",
                    n_models_present=0, n_models=len(paths),
                    rmsd_ca_A=round(med(rmsds), 3), plddt_mut=round(med(plddts), 1),
                    plddt_wt=round(wt["plddt"], 1),
                    bvs_wt=round(med(wt_b), 3), bvs_mut="", delta_bvs="", delta_bvs_pct="",
                    wt_seed_sd=round(noise, 3), z_vs_seed_noise="",
                    cn_wt=round(med(wt_cn), 1), cn_mut="", ion_shift_A="",
                    jaccard_median="",
                ))
                continue

            mu_b = [x["bvs"] for x in present]
            mu_cn = [x["cn"] for x in present]
            jac = [x["jac"] for x in present]
            delta = med(mu_b) - med(wt_b)
            pct = 100.0 * delta / med(wt_b)
            z = delta / noise if noise and not math.isnan(noise) and noise > 0 else float("nan")

            # ion displacement: superpose each mutant model on the wild type, model 0 ion as anchor
            shift = float("nan")
            wt_ion0 = wt["per_site"][site][0]["ion_pos"]
            shifts = []
            for x in present:
                res_m, _, _ = read_model(paths[x["model"]])
                mu_ca = calpha_array(res_m)
                if mu_ca.shape != wt_ca.shape or not mu_ca.size:
                    continue
                _, R, t = kabsch_rmsd(mu_ca, wt_ca)
                v = R @ np.array([x["ion_pos"].x, x["ion_pos"].y, x["ion_pos"].z]) + t
                shifts.append(float(np.linalg.norm(v - np.array([wt_ion0.x, wt_ion0.y, wt_ion0.z]))))
            shift = med(shifts)

            verdict = "below seed noise"
            if not math.isnan(z) and abs(z) >= 2.0:
                verdict = "LOSS beyond seed noise" if delta < 0 else "gain beyond seed noise"
            if n_present < len(paths):
                verdict += f"; site absent in {len(paths)-n_present}/{len(paths)} models"

            log(f"      site {site}: BVS {med(wt_b):.3f} -> {med(mu_b):.3f} "
                f"({pct:+.1f}%, z={z:+.1f} vs seed noise) CN {med(wt_cn):.1f}->{med(mu_cn):.1f} "
                f"ion moved {shift:.2f} A  [{verdict}]")

            rows.append(dict(
                job=name, base=base, mutation=observed, mutation_expected=muttok,
                identity_ok=ident_ok, site=site, status="present",
                n_models_present=n_present, n_models=len(paths),
                rmsd_ca_A=round(med(rmsds), 3), plddt_mut=round(med(plddts), 1),
                plddt_wt=round(wt["plddt"], 1),
                bvs_wt=round(med(wt_b), 3), bvs_mut=round(med(mu_b), 3),
                delta_bvs=round(delta, 3), delta_bvs_pct=round(pct, 1),
                wt_seed_sd=round(noise, 3), z_vs_seed_noise=round(z, 2),
                cn_wt=round(med(wt_cn), 1), cn_mut=round(med(mu_cn), 1),
                ion_shift_A=round(shift, 2), jaccard_median=round(med(jac), 2),
            ))
        if vacated:
            log(f"      {vacated} ion-model(s) matched no wild-type site")
    log()

    if not rows:
        log("FAIL: no mutant/wild-type pairs measured")
        return 1

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {OUT_TSV.relative_to(ROOT)}  ({len(rows)} site rows)")

    # ---- report
    with open(OUT_MD, "w") as fh:
        fh.write("# AF3 mutant cofolding — fold and calcium site\n\n")
        fh.write(f"Generated by `scripts/25_af3_mutant_geometry.py` at "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}.\n\n")
        fh.write("Each value is a median over the five AlphaFold Server models. "
                 "`wt_seed_sd` is the wild-type spread across those five seeds; "
                 "`z_vs_seed_noise` expresses the mutant shift in units of it. "
                 "Sites are matched between wild type and mutant by coordinating-residue "
                 "overlap, not by proximity, so a genuinely vacated site reports as vacated.\n\n")
        fh.write("| job | mutation | site | Ca RMSD | pLDDT mut/wt | BVS wt -> mut | "
                 "delta % | z | CN | ion shift | status |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write(
                f"| {r['job']} | {r['mutation']} | {r['site']} | {r['rmsd_ca_A']} A | "
                f"{r['plddt_mut']}/{r['plddt_wt']} | "
                f"{r['bvs_wt']} -> {r['bvs_mut'] or '—'} | "
                f"{r['delta_bvs_pct'] if r['delta_bvs_pct'] != '' else '—'} | "
                f"{r['z_vs_seed_noise'] if r['z_vs_seed_noise'] != '' else '—'} | "
                f"{r['cn_wt']}->{r['cn_mut'] or '—'} | "
                f"{r['ion_shift_A'] if r['ion_shift_A'] != '' else '—'} | {r['status']} |\n")
        fh.write("\n**Controls.** `P1141L` is the benign negative control and `C1138F` removes a "
                 "disulfide without touching a calcium ligand; both should sit within seed "
                 "noise. Any calcium-ligand variant that does not clear the noise these two "
                 "define has not been shown to disrupt the site.\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
