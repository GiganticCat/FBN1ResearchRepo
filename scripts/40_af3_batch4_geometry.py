#!/usr/bin/env python3
"""
40_af3_batch4_geometry.py — PHASE 5. Batch-4 mutant cofolding: the geometric result as a RATE.

Figure 4 of the draft rests on nine variants and 18 site measurements. It shows that a
calcium-site substitution can leave the fold intact while loosening the ion's coordination, and
it recovers the one experimentally measured case (N2144S). What it cannot say is how often that
happens, which is the first question a reviewer asks. Batch 4 answers it: 60 mutants folded with
their calcium ions against the 21 wild-type tandem constructs already folded in batch 2.

Three design decisions distinguish this from `25_af3_mutant_geometry.py`, which it extends
rather than replaces.

  COGNATE SITE, NOT ALL SITES. Each construct carries three calcium sites, but a substitution
  can only plausibly act on one of them. Reporting all three dilutes a real effect by two
  bystanders and inflates the multiple-testing burden. The cognate site is defined here as the
  site whose wild-type coordination shell contains the mutated position; for control variants,
  which coordinate nothing by construction, it is the site whose ion lies nearest the mutated
  residue. The other two sites are retained and reported as WITHIN-CONSTRUCT controls: the same
  fold, the same five seeds, a site the mutation cannot reach.

  AN EMPIRICAL NOISE FLOOR, NOT A PER-SITE z. `25` scored each site against its own wild-type
  seed spread. That denominator is estimated from five values and is sometimes near zero, which
  turns an irrelevant shift into a large z. Batch 4 was designed with 14 negative controls
  precisely so the null can be measured instead of assumed: the decision threshold below is the
  95th percentile of |ΔBVS%| over the negative controls' own cognate sites. The per-site z is
  still reported, with the denominator floored at the median wild-type spread, so the two
  criteria can be compared.

  LOCAL SUPERPOSITION, NOT GLOBAL. A three-domain construct hinges between domains, and AF3
  models that hinge differently between seeds give a whole-construct Ca RMSD of 25 A while every
  domain and every calcium site inside them is unchanged. Fold disruption and ion displacement
  are therefore measured on the site's own neighbourhood -- the residues within 12 A of the
  wild-type ion -- and the wild-type models are superposed on each other the same way to give
  the displacement its own noise floor. The global RMSD is retained in the table, because a
  reader will want to know a construct hinged, but it is not the fold criterion.

  A QUALITY FILTER FIXED ON THE WILD TYPE ALONE. One construct (cbEGF1-3, pLDDT 74) folds
  inconsistently enough that its sites move between seeds; one of its sites has a seed spread of
  39% of its own valence, and a mutant there returned a 237% BVS change that is model noise, not
  biology. Sites are admitted on wild-type evidence only -- seed spread, site occupancy and
  local confidence, all measured before any mutant is read -- so the filter cannot be tuned,
  even unconsciously, on the answer.

  DIRECTION MATTERS. The hypothesis is loss of coordination. A gain beyond the floor is recorded
  and reported but is not counted as a disruption; it is one of the ways this analysis can
  embarrass itself, and hiding it would defeat the point of folding controls.

Validation before measurement, as elsewhere in this pipeline: every mutant sequence must differ
from its wild-type construct at exactly one position, that position must be the one the manifest
claims, in UniProt numbering, and the residues must be the ones named. A job failing any of
these is excluded and named, never silently carried.

Writes:
  data/processed/af3_batch4_geometry.tsv   one row per (job, site)
  results/phase5_batch4_geometry.md
  logs/40_af3_batch4_geometry_<stamp>.log
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
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
B4DIR = AFDIR / "batch4"
B4MAN = ROOT / "handoff/structures/inbox/batch4_mutants/manifest.csv"
B2MAN = ROOT / "handoff/structures/inbox/batch2_domains/manifest.csv"
FASTA = ROOT / "data/raw/uniprot_P35555_v264.fasta"
OUT_TSV = ROOT / "data/processed/af3_batch4_geometry.tsv"
OUT_MD = ROOT / "results/phase5_batch4_geometry.md"

# Brown-Altermatt bond valence, identical to scripts/23, 25 and 26 so every BVS in this project
# is the same quantity.
R0 = {"O": 1.967, "N": 2.14, "S": 2.45}
B_BV = 0.37
CUTOFF = 3.2
JACCARD_MIN = 0.35

# Radius defining a site's neighbourhood for local superposition. 12 A reaches the whole calcium
# site and the turn either side of it without crossing into the next domain.
LOCAL_R = 12.0

# Wild-type-only admission criteria for a site. Fixed before any mutant was read.
MAX_WT_SD_PCT = 15.0     # seed spread as a fraction of the site's own valence
MIN_WT_MODELS = 5        # the site must be found in every wild-type model
MIN_WT_PLDDT = 70.0      # mean confidence over the coordinating residues

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"40_af3_batch4_geometry_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def reference_seq() -> str:
    return "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))


# ---------------------------------------------------------------- structure reading


def model_paths(job: Path) -> list[Path]:
    out = []
    for i in range(5):
        hits = [h for h in sorted(job.glob(f"*_model_{i}.cif"))
                if not h.name.endswith("Zone.Identifier")]
        if len(hits) == 1:
            out.append(hits[0])
    return out


def read_model(path: Path):
    """(residues, ions, donors) for one CIF. residues=(seqid, one-letter, gemmi residue)."""
    stx = gemmi.read_structure(str(path))
    stx.setup_entities()
    residues, ions, donors = [], [], []
    for ch in stx[0]:
        for res in ch:
            name = res.name.strip()
            if name == "CA" and len(res) == 1:      # a bare calcium, not a protein CA atom
                ions.append(res[0].pos)
                continue
            if name not in AA3:
                continue
            residues.append((res.seqid.num, AA3[name], res))
            for at in res:
                el = at.element.name.upper()
                if el in R0:
                    donors.append((res.seqid.num, at.name, el, at.pos))
    return residues, ions, donors


def shell_of(donors, ion_pos):
    out = [(sid, an, el, ion_pos.dist(p)) for sid, an, el, p in donors if ion_pos.dist(p) < CUTOFF]
    return sorted(out, key=lambda x: x[3])


def bvs_of(shell) -> float:
    return sum(math.exp((R0[el] - d) / B_BV) for _, _, el, d in shell)


def sig_of(shell) -> set[int]:
    return {sid for sid, _, _, _ in shell}


def jaccard(a: set[int], b: set[int]) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0


def assign(shells, refs):
    n_i, n_r = len(shells), len(refs)
    if not n_i or not n_r:
        return [None] * n_i
    cost = np.ones((n_i, n_r))
    for i, sh in enumerate(shells):
        s = sig_of(sh)
        for j, r in enumerate(refs):
            cost[i, j] = 1.0 - jaccard(s, r)
    rows, cols = linear_sum_assignment(cost)
    out = [None] * n_i
    for i, j in zip(rows, cols):
        if 1.0 - cost[i, j] >= JACCARD_MIN:
            out[i] = int(j)
    return out


def ca_array(residues) -> np.ndarray:
    pts = []
    for _, _, res in residues:
        a = res.find_atom("CA", "*")
        if a is None:
            return np.empty((0, 3))
        pts.append([a.pos.x, a.pos.y, a.pos.z])
    return np.asarray(pts)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    U, S, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = Qc - R @ Pc
    return float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean())), R, t


def med(xs):
    return st.median(xs) if xs else float("nan")


def spread(xs):
    """MAD scaled to a standard-deviation equivalent."""
    if len(xs) < 2:
        return float("nan")
    m = st.median(xs)
    return 1.4826 * st.median([abs(x - m) for x in xs])


def pos_of(residues, seqid):
    for sid, aa, res in residues:
        if sid == seqid:
            return res
    return None


# ---------------------------------------------------------------- wild-type constructs


def local_set(residues, ion_pos) -> set[int]:
    """Residue numbers with any atom within LOCAL_R of the ion — the site's neighbourhood."""
    out = set()
    for sid, _, res in residues:
        for at in res:
            if ion_pos.dist(at.pos) < LOCAL_R:
                out.add(sid)
                break
    return out


def ca_subset(residues, keep: set[int]) -> np.ndarray:
    pts = []
    for sid, _, res in residues:
        if sid in keep:
            a = res.find_atom("CA", "*")
            if a is None:
                return np.empty((0, 3))
            pts.append([a.pos.x, a.pos.y, a.pos.z])
    return np.asarray(pts)


def mean_plddt(residues, keep: set[int]) -> float:
    vals = [a.b_iso for sid, _, res in residues if sid in keep for a in res]
    return float(np.mean(vals)) if vals else float("nan")


def measure_wt(path: Path, offset: int):
    """Per-site wild-type description: valence, seed spread, neighbourhood, self-consistency."""
    paths = model_paths(path)
    if len(paths) != 5:
        return None
    models = [read_model(p) for p in paths]
    residues0, ions0, donors0 = models[0]
    plddt_global = med([float(np.mean([a.b_iso for _, _, r in res for a in r]))
                        for res, _, _ in models])

    shells0 = [shell_of(donors0, ip) for ip in ions0]
    refs = [sig_of(sh) for sh in shells0]
    sites = {j: dict(bvs=[], cn=[], sig=refs[j], ion0=ions0[j],
                     local=local_set(residues0, ions0[j]), plddt=[], self_shift=[])
             for j in range(len(refs))}

    for mi, (residues, ions, donors) in enumerate(models):
        shells = [shell_of(donors, ip) for ip in ions]
        for k, j in enumerate(assign(shells, refs)):
            if j is None:
                continue
            S = sites[j]
            S["bvs"].append(bvs_of(shells[k]))
            S["cn"].append(len(shells[k]))
            S["plddt"].append(mean_plddt(residues, S["sig"]))
            if mi == 0:
                continue
            # how far does the ion move between wild-type seeds, after local superposition?
            P, Q = ca_subset(residues, S["local"]), ca_subset(residues0, S["local"])
            if P.shape == Q.shape and P.size:
                _, R, t = kabsch(P, Q)
                v = R @ np.array([ions[k].x, ions[k].y, ions[k].z]) + t
                i0 = S["ion0"]
                S["self_shift"].append(
                    float(np.linalg.norm(v - np.array([i0.x, i0.y, i0.z]))))

    for j, S in sites.items():
        sd = spread(S["bvs"])
        S["sd"] = sd
        S["ok"] = (len(S["bvs"]) >= MIN_WT_MODELS
                   and not math.isnan(sd)
                   and 100.0 * sd / med(S["bvs"]) <= MAX_WT_SD_PCT
                   and med(S["plddt"]) >= MIN_WT_PLDDT)

    return dict(paths=paths, residues=residues0, plddt=plddt_global, sites=sites,
                offset=offset, seq="".join(r[1] for r in residues0))


def main() -> int:
    log(f"Batch-4 mutant geometry — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  BVS: {CUTOFF} A shell, R0(Ca-O)={R0['O']} R0(Ca-N)={R0['N']} b={B_BV}; "
        f"site match Jaccard >= {JACCARD_MIN}")
    log()

    ref = reference_seq()
    b4 = list(csv.DictReader(open(B4MAN)))
    b2 = {r["job_name"]: r for r in csv.DictReader(open(B2MAN))}
    log(f"{len(b4)} batch-4 jobs in the manifest; {len(b2)} wild-type constructs available")

    # index every folder once, both the batch-4 tree and the top level (one job was unzipped to
    # the top level as well as into batch4/ -- byte-identical, so either copy will do)
    dirs: dict[str, Path] = {}
    for d in list(B4DIR.iterdir()) + list(AFDIR.iterdir()):
        if d.is_dir():
            dirs.setdefault(norm(d.name), d)

    # ---- wild types
    wt: dict[str, dict] = {}
    log()
    log("Wild-type constructs (reference sites and seed spread)")
    for cname in sorted({r["construct"] for r in b4}):
        meta = b2.get(cname)
        if meta is None:
            log(f"  {cname}: NOT IN BATCH-2 MANIFEST — every mutant on it will be excluded")
            continue
        d = dirs.get(norm(cname)) or dirs.get(norm("fold_" + cname))
        if d is None:
            log(f"  {cname}: no folder found — excluded")
            continue
        off = int(meta["uniprot_range"].split("-")[0]) - 1
        m = measure_wt(d, off)
        if m is None:
            log(f"  {cname}: fewer than five models — excluded")
            continue
        # numbering check: construct sequence must equal P35555 over the stated range
        lo, hi = (int(x) for x in meta["uniprot_range"].split("-"))
        if m["seq"] != ref[lo - 1:hi]:
            log(f"  {cname}: SEQUENCE DOES NOT MATCH P35555 {lo}-{hi} — excluded")
            continue
        wt[cname] = m
        for j, s in sorted(m["sites"].items()):
            lig = sorted(x + off for x in s["sig"])
            log(f"  {cname:22s} site {j}  BVS {med(s['bvs']):.3f} +/- {spread(s['bvs']):.3f} "
                f"({100*spread(s['bvs'])/med(s['bvs']):4.1f}%)  CN {med(s['cn']):.1f}  "
                f"pLDDT {med(s['plddt']):.0f}  ion seed spread "
                f"{med(s['self_shift']) if s['self_shift'] else float('nan'):.2f} A  "
                f"{'admitted' if s['ok'] else 'EXCLUDED (wild-type quality)'}  "
                f"ligands {','.join(str(x) for x in lig)}")

    noise_floor = med([s["sd"] for m in wt.values() for s in m["sites"].values()
                   if not math.isnan(s["sd"])])
    log()
    log(f"Median wild-type seed spread across all sites: {noise_floor:.3f} BVS units "
        f"(used to floor the per-site z denominator)")
    log()

    # ---- mutants
    rows, excluded = [], []
    log("Mutants")
    for r in b4:
        job, cname, var = r["job_name"], r["construct"], r["variant"]
        d = dirs.get(norm(job))
        if d is None:
            excluded.append((job, "no folder")); log(f"  {job}: NO FOLDER — excluded"); continue
        if cname not in wt:
            excluded.append((job, "wild type unusable")); continue
        W = wt[cname]
        off = W["offset"]
        paths = model_paths(d)
        if len(paths) != 5:
            excluded.append((job, f"{len(paths)} models")); log(f"  {job}: {len(paths)} models — excluded"); continue

        # ---- validation: exactly one substitution, at the position the manifest names
        mut_res, _, _ = read_model(paths[0])
        mseq = "".join(x[1] for x in mut_res)
        if len(mseq) != len(W["seq"]):
            excluded.append((job, "length mismatch")); log(f"  {job}: LENGTH MISMATCH — excluded"); continue
        diffs = [i for i, (a, b) in enumerate(zip(W["seq"], mseq)) if a != b]
        if len(diffs) != 1:
            excluded.append((job, f"{len(diffs)} substitutions")); log(f"  {job}: {len(diffs)} substitutions — excluded"); continue
        i = diffs[0]
        upos = W["residues"][i][0] + off
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", var)
        want_wt, want_pos, want_mut = m.group(1), int(m.group(2)), m.group(3)
        if not (upos == want_pos and W["seq"][i] == want_wt and mseq[i] == want_mut):
            excluded.append((job, f"observed {W['seq'][i]}{upos}{mseq[i]} != {var}"))
            log(f"  {job}: observed {W['seq'][i]}{upos}{mseq[i]} but manifest says {var} — excluded")
            continue
        if ref[upos - 1] != want_wt:
            excluded.append((job, "WT residue disagrees with P35555"))
            log(f"  {job}: P35555 has {ref[upos-1]} at {upos}, not {want_wt} — excluded"); continue

        local = W["residues"][i][0]

        # ---- cognate site: the wild-type site this residue coordinates, else the nearest ion
        cognate, how = None, ""
        for j, s in W["sites"].items():
            if local in s["sig"]:
                cognate, how = j, "coordinates"
                break
        if cognate is None:
            wres = pos_of(W["residues"], local)
            anchor = wres.find_atom("CB", "*") or wres.find_atom("CA", "*")
            dists = {j: anchor.pos.dist(s["ion0"]) for j, s in W["sites"].items()}
            cognate = min(dists, key=dists.get)
            how = f"nearest ion ({dists[cognate]:.1f} A)"

        # ---- measure every model of the mutant
        models = [read_model(q) for q in paths]
        plddts = [float(np.mean([a.b_iso for _, _, rr in res for a in rr]))
                  for res, _, _ in models]
        wt_ca = ca_array(W["residues"])
        rmsds = []
        for res_m, _, _ in models:
            mu_ca = ca_array(res_m)
            if mu_ca.shape == wt_ca.shape and mu_ca.size:
                rmsds.append(kabsch(mu_ca, wt_ca)[0])

        site_keys = sorted(W["sites"])
        refs = [W["sites"][j]["sig"] for j in site_keys]
        per_site: dict[int, list[dict]] = {}
        for mi, (res_m, ions, donors) in enumerate(models):
            shells = [shell_of(donors, ip) for ip in ions]
            for k, sj in enumerate(assign(shells, refs)):
                if sj is None:
                    continue
                j = site_keys[sj]
                S = W["sites"][j]
                # superpose on this site's neighbourhood only: a hinge elsewhere in the
                # construct must not register as the ion having moved
                P, Q = ca_subset(res_m, S["local"]), ca_subset(W["residues"], S["local"])
                lrms = shift = float("nan")
                if P.shape == Q.shape and P.size:
                    lrms, R, t = kabsch(P, Q)
                    v = R @ np.array([ions[k].x, ions[k].y, ions[k].z]) + t
                    i0 = S["ion0"]
                    shift = float(np.linalg.norm(v - np.array([i0.x, i0.y, i0.z])))
                per_site.setdefault(j, []).append(dict(
                    bvs=bvs_of(shells[k]), cn=len(shells[k]), shift=shift, local_rmsd=lrms,
                    plddt=mean_plddt(res_m, S["sig"]),
                    contact=local in sig_of(shells[k]),
                    jac=jaccard(sig_of(shells[k]), S["sig"])))

        log(f"  {job}  {var}  [{r['class']}]  cognate site {cognate} ({how})  "
            f"global RMSD {med(rmsds):.2f} A  pLDDT {med(plddts):.1f} vs {W['plddt']:.1f}")

        for j in site_keys:
            S = W["sites"][j]
            wb, wc = S["bvs"], S["cn"]
            sd = S["sd"]
            denom = max(sd, noise_floor) if not math.isnan(sd) else noise_floor
            present = per_site.get(j, [])
            is_lig = local in S["sig"]
            base = dict(
                job=job, variant=var, variation_id=r["variation_id"], klass=r["class"],
                clinical=r["clinical"], construct=cname, uniprot_pos=upos,
                wt_aa=want_wt, mut_aa=want_mut, site=j,
                is_cognate=(j == cognate), cognate_by=how if j == cognate else "",
                wt_site_ok=S["ok"], is_wt_ligand_of_site=is_lig,
                n_models_present=len(present), n_models=len(paths),
                rmsd_ca_global_A=round(med(rmsds), 3),
                plddt_mut=round(med(plddts), 1), plddt_wt=round(W["plddt"], 1),
                plddt_site_wt=round(med(S["plddt"]), 1),
                bvs_wt=round(med(wb), 3), wt_seed_sd=round(sd, 3),
                wt_seed_sd_pct=round(100.0 * sd / med(wb), 1),
                cn_wt=round(med(wc), 1),
                ion_shift_wt_sd_A=(round(med(S["self_shift"]), 2) if S["self_shift"] else ""),
            )
            if not present:
                rows.append({**base, "status": "vacated", "bvs_mut": "", "delta_bvs": "",
                             "delta_bvs_pct": "", "z_vs_noise": "", "cn_mut": "",
                             "ion_shift_A": "", "local_rmsd_A": "", "plddt_site_mut": "",
                             "jaccard_median": "", "still_contacts": ""})
                if j == cognate:
                    log(f"      site {j} (cognate): VACATED in all five models")
                continue
            mb = [x["bvs"] for x in present]
            delta = med(mb) - med(wb)
            pct = 100.0 * delta / med(wb)
            rows.append({**base, "status": "present", "bvs_mut": round(med(mb), 3),
                         "delta_bvs": round(delta, 3), "delta_bvs_pct": round(pct, 1),
                         "z_vs_noise": round(delta / denom, 2),
                         "cn_mut": round(med([x["cn"] for x in present]), 1),
                         "ion_shift_A": round(med([x["shift"] for x in present]), 2),
                         "local_rmsd_A": round(med([x["local_rmsd"] for x in present]), 2),
                         "plddt_site_mut": round(med([x["plddt"] for x in present]), 1),
                         "jaccard_median": round(med([x["jac"] for x in present]), 2),
                         "still_contacts": (sum(x["contact"] for x in present) if is_lig else "")})
            if j == cognate:
                log(f"      site {j} (cognate): BVS {med(wb):.3f} -> {med(mb):.3f} "
                    f"({pct:+.1f}%, z={delta/denom:+.1f})  CN {med(wc):.1f}->"
                    f"{med([x['cn'] for x in present]):.1f}  "
                    f"ion {med([x['shift'] for x in present]):.2f} A "
                    f"(wt seed {base['ion_shift_wt_sd_A']} A)  "
                    f"local RMSD {med([x['local_rmsd'] for x in present]):.2f} A"
                    + (f"  contact kept {sum(x['contact'] for x in present)}/5" if is_lig else "")
                    + ("" if S["ok"] else "   [SITE FAILS WILD-TYPE QUALITY FILTER]"))

    log()
    if excluded:
        log(f"{len(excluded)} job(s) excluded:")
        for j, why in excluded:
            log(f"  {j}: {why}")
    else:
        log("No jobs excluded — all 60 validated and measured.")

    if not rows:
        log("FAIL: nothing measured")
        return 1

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    log(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(rows)} site rows, "
        f"{len({r['job'] for r in rows})} jobs)")

    # ---- report
    cog = [r for r in rows if r["is_cognate"]]
    with open(OUT_MD, "w") as fh:
        fh.write("# Batch-4 mutant cofolding — calcium coordination at the cognate site\n\n")
        fh.write(f"Generated by `scripts/40_af3_batch4_geometry.py`, "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}.\n\n")
        fh.write(f"{len({r['job'] for r in rows})} mutants, each folded five times with its "
                 "calcium ions and compared to the five-model wild type of the same construct. "
                 "The **cognate site** is the calcium site whose wild-type shell contains the "
                 "mutated residue; for control variants, which coordinate nothing, it is the "
                 "site with the nearest ion. The other two sites of each construct are "
                 "within-construct controls and are in the TSV.\n\n")
        fh.write("| job | variant | class | site | BVS wt → mut | Δ% | z | CN | ion shift | "
                 "contact kept |\n|---|---|---|---|---|---|---|---|---|---|\n")
        for r in sorted(cog, key=lambda x: (x["klass"], x["delta_bvs_pct"] if x["delta_bvs_pct"] != "" else 0)):
            fh.write(f"| {r['job']} | {r['variant']} | {r['klass']} | {r['site']} | "
                     f"{r['bvs_wt']} → {r['bvs_mut'] or '—'} | {r['delta_bvs_pct']} | "
                     f"{r['z_vs_noise']} | {r['cn_wt']}→{r['cn_mut'] or '—'} | "
                     f"{r['ion_shift_A']} | {r['still_contacts']}/5 |\n")
        fh.write("\nRates, thresholds and the comparison against controls are computed in "
                 "`scripts/41_batch4_rate.py` → `results/phase5_batch4_rate.md`.\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
