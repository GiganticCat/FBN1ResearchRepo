#!/usr/bin/env python3
"""
27_structural_metrics_af3.py — PHASE 5 (revision). Structural metrics over all 43 cbEGF domains.

Replaces the four-template structural table (751 variants) with one built on the validated AF3
batch-2 constructs (3,831 variants, 91% of the pathogenic set). Two definitions change here, and
both were wrong before rather than merely narrower.

CALCIUM LIGANDS ARE NOW STRUCTURE-DERIVED. UniProt annotates no calcium site anywhere in P35555,
so the previous `is_ca_consensus` flag came from a sequence motif and was never validated against
anything. Measured against the batch-2 structures it over-calls by 43%: of 197 consensus
positions, 85 coordinate no calcium in any model. This script uses the observed coordination
shells (data/processed/af3_ca_ligands.tsv) as the primary definition and carries the consensus
flag alongside, labelled, for a sensitivity analysis.

CONSTRUCT CHOICE IS EXPLICIT. The tandem constructs overlap by one domain, so most positions are
modelled more than once. Metrics are taken from the construct in which the residue sits furthest
from either terminus, because a residue near a construct edge has artificially high solvent
exposure and an artificially truncated environment. The chosen construct and the margin are
recorded per variant so the choice is auditable rather than implicit.

Every geometric quantity is a median over the five models. Seed spread is carried for the
calcium distance, since the batch-3 work showed model-to-model variation large enough to matter.

Writes:
  data/processed/structural_metrics_af3.tsv
  results/phase5_structural_af3.md
  logs/27_structural_metrics_af3_<stamp>.log
"""

from __future__ import annotations

import csv
import datetime as dt
import math
import re
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

import freesasa
import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
SITES = ROOT / "data" / "processed" / "af3_ca_sites.tsv"
LIGANDS = ROOT / "data" / "processed" / "af3_ca_ligands.tsv"
VARIANTS = ROOT / "data" / "processed" / "variant_annotated.tsv"
GFF = ROOT / "data" / "raw" / "uniprot_P35555_v264.gff"
FASTA = ROOT / "data" / "raw" / "uniprot_P35555_v264.fasta"
OUT = ROOT / "data" / "processed" / "structural_metrics_af3.tsv"
OUT_MD = ROOT / "results" / "phase5_structural_af3.md"

# Tien et al. 2013, theoretical maximum accessible surface area (A^2), used for RSA.
MAXASA = {
    "A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225, "G": 104,
    "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240, "P": 159, "S": 155,
    "T": 172, "W": 285, "Y": 263, "V": 174,
}
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BURIED_RSA = 20.0

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"27_structural_metrics_af3_{STAMP}.log"
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
    return "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))


def job_dir(job_name: str) -> Path | None:
    want = re.sub(r"[^a-z0-9]+", "_", job_name.lower()).strip("_")
    for d in AFDIR.iterdir():
        if d.is_dir():
            have = re.sub(r"[^a-z0-9]+", "_", d.name.lower()).strip("_")
            if have in (want, f"fold_{want}"):
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


def disulfides_from_gff() -> dict[int, int]:
    """UniProt disulfide bonds -> {pos: partner}. UniProt does annotate these (155 of them)."""
    pairs = {}
    for line in GFF.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) > 4 and f[2] == "Disulfide bond":
            a, b = int(f[3]), int(f[4])
            if a != b:
                pairs[a] = b
                pairs[b] = a
    return pairs


def measure_construct(job_name: str, lo: int, hi: int, ref: str):
    """Per-residue metrics for one construct, medians over its five models.

    Returns {uniprot_pos: dict} or None if the construct cannot be read.
    """
    d = job_dir(job_name)
    if d is None:
        return None
    paths = model_paths(d)
    if len(paths) != 5:
        return None

    acc: dict[int, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    aa_at: dict[int, str] = {}

    for p in paths:
        stx = gemmi.read_structure(str(p))
        stx.setup_entities()
        stx.remove_alternative_conformations()
        mdl = stx[0]

        ions, residues = [], []
        for ch in mdl:
            for res in ch:
                nm = res.name.strip()
                if nm == "CA" and len(res) == 1:
                    ions.append(res[0].pos)
                elif nm in AA3:
                    residues.append(res)
        residues.sort(key=lambda r: r.seqid.num)
        if not residues:
            return None
        offset = lo - residues[0].seqid.num

        # SASA on the protein alone: the ion would occlude its own ligands and confound
        # burial with coordination, which are separate variables in the analysis.
        tmp = gemmi.Structure()
        tmp.add_model(gemmi.Model("1"))
        chain = gemmi.Chain("A")
        for res in residues:
            chain.add_residue(res)
        tmp[0].add_chain(chain)
        tmp.setup_entities()
        pdb_text = tmp.make_pdb_string()

        sasa_by_res: dict[int, float] = {}
        try:
            struct = freesasa.Structure()
            for line in pdb_text.splitlines():
                if line.startswith("ATOM"):
                    struct.addAtom(line[12:16], line[17:20], line[22:26], line[21],
                                   float(line[30:38]), float(line[38:46]), float(line[46:54]))
            result = freesasa.calc(struct)
            for i in range(struct.nAtoms()):
                rn = int(struct.residueNumber(i))
                sasa_by_res[rn] = sasa_by_res.get(rn, 0.0) + result.atomArea(i)
        except Exception as exc:                       # pragma: no cover - reported, not hidden
            log(f"    SASA failed on {p.name}: {exc}")
            return None

        for res in residues:
            up = res.seqid.num + offset
            one = AA3[res.name.strip()]
            aa_at[up] = one
            sa = sasa_by_res.get(res.seqid.num, 0.0)
            acc[up]["sasa"].append(sa)
            acc[up]["rsa"].append(100.0 * sa / MAXASA[one])
            acc[up]["plddt"].append(float(np.mean([a.b_iso for a in res])))
            if ions:
                dmin = min(min(ion.dist(a.pos) for a in res) for ion in [i for i in ions])
                acc[up]["dca"].append(dmin)

    out = {}
    for up, vals in acc.items():
        rsa = st.median(vals["rsa"])
        dca = st.median(vals["dca"]) if vals["dca"] else float("nan")
        out[up] = dict(
            aa=aa_at[up],
            sasa=round(st.median(vals["sasa"]), 1),
            rsa=round(rsa, 1),
            buried=rsa < BURIED_RSA,
            plddt=round(st.median(vals["plddt"]), 1),
            dist_ca=round(dca, 2) if not math.isnan(dca) else "",
            dist_ca_sd=round(st.pstdev(vals["dca"]), 2) if len(vals["dca"]) > 1 else "",
            margin=min(up - lo, hi - up),
        )
    return out


def main() -> int:
    log(f"Phase 5 revision — structural metrics on AF3 batch 2 — "
        f"{dt.datetime.now(dt.timezone.utc).isoformat()}")
    ref = read_reference()
    ss_pairs = disulfides_from_gff()
    log(f"  reference {len(ref)} aa; {len(ss_pairs)//2} disulfide bonds from UniProt")

    # ---- ligand definitions
    lig_rows = list(csv.DictReader(open(LIGANDS), delimiter="\t"))
    ligand_atoms: dict[int, set[str]] = defaultdict(set)
    for r in lig_rows:
        ligand_atoms[int(r["uniprot_pos"])].add(r["atom_class"])
    sidechain_ligands = {p for p, c in ligand_atoms.items() if "sidechain" in c}
    backbone_only = {p for p, c in ligand_atoms.items() if c == {"backbone"}}
    log(f"  {len(sidechain_ligands)} side-chain calcium ligands, "
        f"{len(backbone_only)} backbone-only, from {len(lig_rows)} coordinating atoms")

    # ---- constructs
    site_rows = list(csv.DictReader(open(SITES), delimiter="\t"))
    constructs = {}
    for r in site_rows:
        lo, hi = (int(x) for x in r["uniprot_range"].split("-"))
        constructs[r["construct"]] = (lo, hi, r["cbegf_domains"])

    log()
    log(f"Measuring {len(constructs)} constructs, 5 models each")
    measured = {}
    for name, (lo, hi, dom) in sorted(constructs.items()):
        m = measure_construct(name, lo, hi, ref)
        if m is None:
            log(f"  {name:26s} FAILED — excluded")
            continue
        # verify identity against the reference before trusting any position
        bad = [p for p, v in m.items() if ref[p - 1] != v["aa"]]
        if bad:
            log(f"  {name:26s} FAILED — {len(bad)} residue identity mismatches vs P35555")
            continue
        measured[name] = m
        log(f"  {name:26s} OK  {len(m)} residues")

    if not measured:
        log("FAIL: no constructs measured")
        return 1

    # ---- pick, for each position, the construct where it sits furthest from an edge
    best: dict[int, tuple[str, dict]] = {}
    for name, m in measured.items():
        for up, v in m.items():
            if up not in best or v["margin"] > best[up][1]["margin"]:
                best[up] = (name, v)
    log()
    log(f"{len(best)} distinct P35555 positions covered by at least one construct")

    # ---- join to variants
    variants = list(csv.DictReader(open(VARIANTS), delimiter="\t"))
    rows, n_cov, n_mismatch = [], 0, 0
    for v in variants:
        if not v["position"]:
            continue
        pos = int(v["position"])
        wt = v["wt_aa"]
        if pos not in best:
            continue
        name, mm = best[pos]
        if mm["aa"] != wt:
            n_mismatch += 1
            continue
        n_cov += 1
        lo, hi, dom = constructs[name]
        is_sc_lig = pos in sidechain_ligands
        is_bb_lig = pos in backbone_only
        partner = ss_pairs.get(pos, "")
        rows.append(dict(
            variation_id=v["variation_id"], position=pos, wt_aa=wt, mut_aa=v["mut_aa"],
            set_primary=v["set_primary"], am_pathogenicity=v.get("am_pathogenicity", ""),
            construct=name, construct_range=f"{lo}-{hi}", cbegf_domains=dom,
            edge_margin=mm["margin"], plddt=mm["plddt"],
            sasa_abs_A2=mm["sasa"], rsa_pct=mm["rsa"], buried=mm["buried"],
            dist_to_nearest_ca_A=mm["dist_ca"], dist_ca_seed_sd=mm["dist_ca_sd"],
            ca_ligand_sidechain=is_sc_lig,
            ca_ligand_backbone_only=is_bb_lig,
            ca_ligand_any=is_sc_lig or is_bb_lig,
            is_ca_consensus_legacy=v.get("is_ca_consensus", ""),
            consensus_but_not_ligand=(v.get("is_ca_consensus") == "True"
                                      and not (is_sc_lig or is_bb_lig)),
            is_cysteine=wt == "C",
            cys_removing=(wt == "C" and v["mut_aa"] != "C"),
            cys_introducing=(wt != "C" and v["mut_aa"] == "C"),
            disulfide_partner=partner,
        ))

    log(f"{n_cov} variants given structural metrics; {n_mismatch} dropped on WT identity mismatch")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"wrote {OUT.relative_to(ROOT)}  ({len(rows)} rows)")

    # ---- summary
    from collections import Counter
    cls = Counter(r["set_primary"] or "unlabelled" for r in rows)
    lig = [r for r in rows if r["ca_ligand_sidechain"]]
    cys = [r for r in rows if r["cys_removing"]]
    disagree = sum(1 for r in rows if r["consensus_but_not_ligand"])
    log()
    log(f"  by class: {dict(cls)}")
    log(f"  at a side-chain calcium ligand : {len(lig)} "
        f"({Counter(r['set_primary'] or 'unlabelled' for r in lig)})")
    log(f"  cysteine-removing              : {len(cys)}")
    log(f"  consensus-flagged but not a ligand in any model: {disagree}")
    med_rsa_lig = st.median([r["rsa_pct"] for r in lig]) if lig else float("nan")
    log(f"  median RSA at calcium ligands  : {med_rsa_lig:.1f}%")

    with open(OUT_MD, "w") as fh:
        fh.write("# Structural metrics across all 43 cbEGF domains\n\n")
        fh.write(f"Generated by `scripts/27_structural_metrics_af3.py` at "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}.\n\n")
        fh.write(f"{len(rows)} variants have structural metrics, up from 751 on the four "
                 f"experimental templates. Values are medians over the five AlphaFold Server "
                 f"models of the chosen construct; where tandem constructs overlap, the one "
                 f"placing the residue furthest from a terminus is used and recorded.\n\n")
        fh.write("Calcium ligands are defined by observed coordination in the models, not by "
                 "sequence motif: UniProt annotates no calcium site in P35555. "
                 f"{disagree} variants sit at positions the motif flags as calcium-consensus but "
                 "which coordinate no calcium in any model; they are labelled "
                 "`consensus_but_not_ligand` rather than dropped.\n\n")
        fh.write("| group | n |\n|---|---|\n")
        for k, n in cls.most_common():
            fh.write(f"| {k} | {n} |\n")
        fh.write(f"| at side-chain Ca ligand | {len(lig)} |\n")
        fh.write(f"| cysteine-removing | {len(cys)} |\n")
        fh.write(f"| consensus but not a ligand | {disagree} |\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
