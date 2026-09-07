#!/usr/bin/env python3
"""
29_rosetta_ddg.py — PHASE 5 (revision). Metal-aware ddG in a second force field.

FoldX reports essentially nothing at calcium ligands: the metal term fires at 39 of 78 direct
ligands and the whole ion contributes a median +0.247 kcal/mol. That leaves the "fold is intact"
half of the thesis resting on a tool that cannot see the interaction being destroyed. Rosetta is
the second opinion, chosen because `-in:auto_setup_metals` detects the ion, identifies its
coordinating residues, and adds distance constraints plus a `metalbinding_constraint` score term
that is actually evaluated.

The two force fields have opposite biases and that is the point. FoldX under-credits metal
coordination; Rosetta's metal term is a harmonic restraint whose depth is a chosen weight, so it
over-penalises ligand loss by construction. Neither is the truth. Reported together they bracket
it, and the bracket is the honest deliverable -- swapping FoldX's artifact for Rosetta's would
not be progress.

Protocol, per mutation: load the AF3 construct with its ion, bond the ion to its ligands, relax,
then score wild type and mutant after local repacking and minimisation within a shell around the
mutated position. Repeated `--repeats` times from different random seeds; the reported value is
the mean, with the spread carried so a reader can see when a number is unstable.

Torsion space, and metal bonds without metal constraints. Neither is a preference -- both are
forced, and both were found by measurement rather than assumed (see init_worker and
_setup_metals). Cartesian minimisation cannot be used because ref2015_cart's cart_bonded term
scores the metal-ligand bonds against ideal parameters that do not exist for calcium, putting
7693 of 12336 raw REU of fictitious strain on the ions and their ligands -- the exact positions
this study mutates. Harmonic metal constraints cannot be used because the packer cannot see them,
so they detonate to +17364 REU across the repack step. Validated on controls: a disulfide-removing
substitution costs +22.4 REU, the benign control +2.4 REU.

Deliberately not done here: modelling the mutant's calcium affinity. Rosetta reports a folding
free-energy difference for a fixed metal geometry. Whether the mutant still binds calcium is the
question AF3 cofolding and the bond-valence analysis answer.

Usage:
  python scripts/29_rosetta_ddg.py --calibrate       # time a handful, write nothing
  python scripts/29_rosetta_ddg.py [--repeats 3] [--jobs 4] [--limit N]

Writes:
  data/processed/rosetta_ddg.tsv
  logs/29_rosetta_ddg_<stamp>.log
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import os
import re
import statistics as st
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
SM = ROOT / "data" / "processed" / "structural_metrics_af3.tsv"
WORK = Path("/tmp/fbn1-rosetta")
OUT = ROOT / "data" / "processed" / "rosetta_ddg.tsv"

SHELL_A = 8.0          # repack/minimise radius around the mutated residue
CST_SD_A = 0.5         # harmonic width on CA restraints during minimisation
BASE_SEED = 20260809

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"29_rosetta_ddg_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg, flush=True)
    _fh.write(msg + "\n")
    _fh.flush()


AA1TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU", "Q": "GLN",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def job_dir(job_name: str) -> Path | None:
    want = re.sub(r"[^a-z0-9]+", "_", job_name.lower()).strip("_")
    for d in AFDIR.iterdir():
        if d.is_dir():
            have = re.sub(r"[^a-z0-9]+", "_", d.name.lower()).strip("_")
            if have in (want, f"fold_{want}"):
                return d
    return None


def prepare_pdb(job_name: str, lo: int) -> tuple[Path, int] | None:
    """Write model_0 of a construct as PDB with its calcium ions, renumbered to UniProt.

    Returns (path, offset_applied). Rosetta reads author numbering, so putting UniProt numbers
    in the file removes the SIFTS-style offset bookkeeping that FoldX mutation strings needed.
    """
    import gemmi

    d = job_dir(job_name)
    if d is None:
        return None
    hits = [h for h in sorted(d.glob("*_model_0.cif")) if not h.name.endswith("Zone.Identifier")]
    if not hits:
        return None
    stx = gemmi.read_structure(str(hits[0]))
    stx.setup_entities()
    stx.remove_alternative_conformations()
    stx.remove_hydrogens()
    mdl = stx[0]

    prot, ions = [], []
    for ch in mdl:
        for res in ch:
            nm = res.name.strip()
            if nm == "CA" and len(res) == 1:
                ions.append(res)
            elif nm in AA1TO3.values():
                prot.append(res)
    prot.sort(key=lambda r: r.seqid.num)
    if not prot:
        return None
    offset = lo - prot[0].seqid.num

    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    chain = gemmi.Chain("A")
    for res in prot:
        r2 = res.clone()
        r2.seqid.num = res.seqid.num + offset
        r2.het_flag = "A"
        chain.add_residue(r2)
    out[0].add_chain(chain)
    ionchain = gemmi.Chain("B")
    for i, res in enumerate(ions, start=1):
        r2 = res.clone()
        r2.seqid.num = i
        r2.het_flag = "H"
        r2.name = " CA"
        ionchain.add_residue(r2)
    out[0].add_chain(ionchain)
    out.setup_entities()

    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / f"{job_name}.pdb"
    path.write_text(out.make_pdb_string())
    return path, offset


# ------------------------------------------------------------------ Rosetta, per worker

_PYR = {}


def init_worker(repeats: int, protocol: str = "torsion", metals: bool = True):
    """Start Rosetta in one worker.

    `protocol` selects the minimisation space, and it is not a free choice -- see the note on
    cart_bonded below. `metals` toggles -in:auto_setup_metals, which is what emits the
    metal-ligand constraints that make this force field worth running at all.
    """
    import pyrosetta
    from pyrosetta.rosetta.core.scoring import ScoreType

    # Deliberately NOT -in:auto_setup_metals. That flag bonds the ion to its ligands at load time,
    # and those bonds cannot be torn down afterwards: removing the ligand's SIDECHAIN_CONJUGATION
    # variant leaves the ion's own connection records pointing at the old partner atom, so
    # rebuilding the constraints on a mutated pose dies on a garbage atom index. Loading clean and
    # calling the setup explicitly per scoring pass (see _score_variant) means the bonds are only
    # ever built from coordinates that currently exist.
    pyrosetta.init(
        "-mute all -ignore_unrecognized_res false -ignore_zero_occupancy false "
        "-detect_disulf true",
        silent=True,
    )

    # ref2015_cart carries cart_bonded at weight 0.5. auto_setup_metals declares real chemical
    # bonds between Ca2+ and its coordinating atoms, and cart_bonded then scores those bonds
    # against ideal parameters that do not exist for calcium. Measured on FBN1_cbEGF1-3_Ca3:
    # the three ions alone contribute 7693 of 12336 raw cart_bonded, and the next worst residues
    # are the calcium ligands themselves (ASP245 407, ASP43 439, GLU4 495) against a protein
    # median of 1.10. That strain sits precisely on the positions this study mutates, so the
    # Cartesian protocol cannot be used together with metal setup. Torsion space has no
    # cart_bonded term and is the default here for that reason.
    sf = pyrosetta.create_score_function(
        "ref2015_cart" if protocol == "cart" else "ref2015")
    # auto_setup_metals emits the constraints; they only bite if the term is weighted.
    sf.set_weight(ScoreType.metalbinding_constraint, 1.0)
    sf.set_weight(ScoreType.atom_pair_constraint, 1.0)
    sf.set_weight(ScoreType.angle_constraint, 1.0)
    sf.set_weight(ScoreType.dihedral_constraint, 1.0)
    sf.set_weight(ScoreType.coordinate_constraint, 1.0)

    # Reporting score function: identical, minus the coordinate restraint. The restraint is a
    # modelling device that keeps the shell from drifting; it is not part of the protein's
    # energy, and leaving it in the reported value inflates exactly the substitutions that need
    # the most backbone movement. Measured on C1876W, 58.9 of its 139.3 REU was restraint.
    sf_report = sf.clone()
    sf_report.set_weight(ScoreType.coordinate_constraint, 0.0)

    _PYR["pyrosetta"] = pyrosetta
    _PYR["sf"] = sf
    _PYR["sf_report"] = sf_report
    _PYR["repeats"] = repeats
    _PYR["protocol"] = protocol
    _PYR["metals"] = metals
    _PYR["cache"] = {}


def relaxed_path(pdb_path: str) -> Path:
    # The protocol is in the filename: a Cartesian-relaxed pose is not a valid starting point for
    # a torsion-space run, and silently reusing one across protocols would compare two different
    # reference states.
    tag = _PYR.get("protocol", "torsion")
    metal = "metal" if _PYR.get("metals", True) else "nometal"
    return Path(pdb_path).with_name(f"{Path(pdb_path).stem}_relaxed_{tag}_{metal}.pdb")


def prerelax(pdb_path: str) -> str:
    """Cartesian-relax a construct once and cache it to disk.

    Relaxation is the expensive step and is identical for every mutation in a construct, so it
    must happen once overall rather than once per worker -- otherwise 21 constructs across 4
    workers pay for it up to 84 times.
    """
    out = relaxed_path(pdb_path)
    if out.exists():
        return str(out)
    pyrosetta = _PYR["pyrosetta"]
    from pyrosetta.rosetta.protocols.relax import FastRelax

    pose = pyrosetta.pose_from_pdb(pdb_path)
    if _PYR["metals"]:
        # Loading no longer sets metals up, so the relax must be told about them explicitly --
        # otherwise the reference state is relaxed with the ion floating free of its ligands.
        _setup_metals(pose, _PYR["sf"])
    fr = FastRelax(_PYR["sf"], 3)
    fr.cartesian(_PYR["protocol"] == "cart")
    fr.max_iter(200)
    fr.apply(pose)
    pose.dump_pdb(str(out))
    return str(out)


def _relaxed_pose(pdb_path: str):
    """Load the pre-relaxed construct, relaxing it here only if the cache is missing."""
    if pdb_path in _PYR["cache"]:
        return _PYR["cache"][pdb_path]
    pyrosetta = _PYR["pyrosetta"]
    src = relaxed_path(pdb_path)
    if not src.exists():
        src = Path(prerelax(pdb_path))
    pose = pyrosetta.pose_from_pdb(str(src))
    _PYR["cache"][pdb_path] = pose
    return pose


def _setup_metals(pose, sf) -> None:
    """Bond each ion to the ligands present in the pose's current coordinates.

    Bonds only. `auto_setup_all_metal_constraints` additionally pins every metal-ligand distance
    and angle to its input value with a harmonic, and that cannot be combined with repacking:
    constraint terms are not pairwise-decomposable, so the packer optimises blind to them and the
    rescore then sees the violation. Measured on cbEGF11-13 with a 12-residue shell,
    metalbinding_constraint goes 0.0 -> +17364 REU across the pack step alone, which swamped every
    real term and produced ddG values of -8022 REU at calcium ligands.

    Note that the constraint helper also re-enables its own score-function weights when called, so
    zeroing the weights is not a way to disable it -- it has to not be called.

    With bonds only, the controls behave: C1138F (removes a disulfide) +22.4 REU, P1141L (benign)
    +2.4 REU, and absolute scores stay near -180 REU for a 130-residue construct.
    """
    from pyrosetta.rosetta.core.util import auto_setup_all_metal_bonds
    auto_setup_all_metal_bonds(pose, 1.0, True)


def _score_variant(pose, seqpos: int, target_aa: str | None, seed: int) -> float:
    """Repack + Cartesian-minimise a shell around seqpos, optionally mutating first."""
    pyrosetta = _PYR["pyrosetta"]
    sf = _PYR["sf"]
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
    from pyrosetta.rosetta.core.kinematics import MoveMap
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    # Repacking draws from the global RNG, so the repeats only differ if it is reseeded here.
    pyrosetta.rosetta.basic.random.init_random_generators(seed, "mt19937")

    work = pose.clone()
    if target_aa is not None:
        # NOT pyrosetta.toolbox.mutate_residue: it routes the pose through
        # packed_pose.to_pose(), which returns a copy, so the mutation lands on a throwaway and
        # the caller's pose is returned unchanged. That is silent -- it raises nothing and
        # reports success -- and it made every ddG in the first calibration exactly 0.000.
        MutateResidue(target=seqpos, new_res=AA1TO3[target_aa]).apply(work)
        if work.residue(seqpos).name1() != target_aa:
            raise RuntimeError(
                f"mutation did not take at {seqpos}: still {work.residue(seqpos).name3()}")

    if _PYR["metals"]:
        # Built here, after the substitution, from the coordinates as they now stand: the ion
        # coordinates whatever ligands actually remain, which is the physical question. The
        # wild-type pass runs this identically, so the two scores differ by the substitution
        # rather than by how their metal setup was constructed.
        _setup_metals(work, sf)

    centre = work.residue(seqpos).nbr_atom_xyz()
    shell = [i for i in range(1, work.total_residue() + 1)
             if work.residue(i).nbr_atom_xyz().distance(centre) <= SHELL_A]

    # Restrain every CA to where it started, so the shell below can open up for a bulky side
    # chain without the pose drifting somewhere else entirely. Applied identically to the wild
    # type and the mutant, so the restraint is a common offset, and removed from the reported
    # score afterwards.
    from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
    from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
    from pyrosetta.rosetta.core.id import AtomID
    anchor, fn = AtomID(1, 1), HarmonicFunc(0.0, CST_SD_A)
    for i in range(1, work.total_residue() + 1):
        if work.residue(i).is_protein():
            aid = AtomID(work.residue(i).atom_index("CA"), i)
            work.add_constraint(CoordinateConstraint(aid, anchor, work.residue(i).xyz("CA"), fn))

    tf = TaskFactory()
    task = tf.create_task_and_apply_taskoperations(work)
    task.restrict_to_repacking()
    for i in range(1, work.total_residue() + 1):
        if i not in shell:
            task.nonconst_residue_task(i).prevent_repacking()
    PackRotamersMover(sf, task).apply(work)

    # Backbone AND side chains, but only inside the shell and only against the restraints above.
    #
    # A fixed backbone cannot accommodate a large residue replacing a small buried one. With the
    # backbone frozen, C1511R (Arg for a disulfide-bonded Cys) reported 3257 REU: the packer had
    # no rotamer that missed the backbone, so it returned an unrelieved collision rather than an
    # energy. 21.5% of the first full run exceeded 100 REU that way, almost all of them Y/F/R/W
    # at cysteines. Freeing the shell backbone brings C1511R to 28.7 and C2522Y to 23.5 while
    # leaving the controls where they were (disulfide loss 22.4 -> 16.7, benign 2.4 -> 3.9,
    # Ca ligands unchanged).
    #
    # An earlier attempt at backbone freedom in Cartesian space, WITHOUT restraints, diverged
    # globally and produced 43 REU for a single aspartate-to-valine. The restraints are what
    # make the freedom safe.
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    for i in shell:
        mm.set_chi(i, True)
        mm.set_bb(i, True)
    mover = MinMover(mm, sf, "lbfgs_armijo_nonmonotone", 1e-4, True)
    mover.cartesian(_PYR["protocol"] == "cart")
    mover.apply(work)
    return float(_PYR["sf_report"](work))


def run_one(task: dict) -> dict:
    """One variant: repeats x (wild type, mutant) -> mean ddG."""
    t0 = time.time()
    out = dict(task)
    try:
        pose = _relaxed_pose(task["pdb"])
        pos = task["position"]
        if pose.pdb_info() is None:
            out.update(status="no_pdb_info", ddg_mean="", ddg_sd="", n_repeats=0)
            return out
        # `pos` is a P35555 number and the pose carries UniProt numbering in its PDB info, so
        # pdb2pose is the only valid lookup -- comparing it against total_residue() is a
        # category error that rejected every position past the construct's length.
        seqpos = pose.pdb_info().pdb2pose("A", pos)
        if seqpos == 0:
            out.update(status="position_not_in_pose", ddg_mean="", ddg_sd="", n_repeats=0)
            return out
        obs = pose.residue(seqpos).name1()
        if obs != task["wt_aa"]:
            out.update(status=f"wt_mismatch_{obs}", ddg_mean="", ddg_sd="", n_repeats=0)
            return out

        vals = []
        for k in range(_PYR["repeats"]):
            seed = BASE_SEED + k
            wt = _score_variant(pose, seqpos, None, seed)
            mu = _score_variant(pose, seqpos, task["mut_aa"], seed)
            vals.append(mu - wt)
        out.update(
            status="ok",
            ddg_mean=round(st.mean(vals), 3),
            # Park et al. 2016 report Rosetta Cartesian ddG in REU and convert to kcal/mol by
            # dividing by 2.94. This protocol is not identical to theirs, so the scaled column is
            # labelled as a convention rather than a measurement; REU is the primary value.
            ddg_kcal_scaled=round(st.mean(vals) / 2.94, 3),
            ddg_sd=round(st.pstdev(vals), 3) if len(vals) > 1 else 0.0,
            n_repeats=len(vals),
            seconds=round(time.time() - t0, 1),
        )
    except Exception as exc:                                # reported per variant, never guessed
        out.update(status=f"error:{type(exc).__name__}:{str(exc)[:120]}",
                   ddg_mean="", ddg_sd="", n_repeats=0)
    return out


def build_tasks(limit: int | None) -> list[dict]:
    rows = list(csv.DictReader(open(SM), delimiter="\t"))
    T = lambda r, k: r[k] == "True"                                        # noqa: E731

    # Priority: the calcium question, its composition control, and the cysteine contrast.
    DEN = {"D", "E", "N"}
    keep = []
    for r in rows:
        why = None
        if T(r, "ca_ligand_sidechain"):
            why = "ca_ligand_sidechain"
        elif T(r, "ca_ligand_backbone_only"):
            why = "ca_ligand_backbone"
        elif T(r, "consensus_but_not_ligand"):
            why = "motif_overcalled"
        elif r["wt_aa"] in DEN and not T(r, "cys_removing"):
            why = "composition_control_DEN"
        elif T(r, "cys_removing"):
            why = "cys_removing"
        if why:
            keep.append((why, r))

    order = {"ca_ligand_sidechain": 0, "ca_ligand_backbone": 1, "motif_overcalled": 2,
             "composition_control_DEN": 3, "cys_removing": 4}
    keep.sort(key=lambda x: (order[x[0]], x[1]["construct"], int(x[1]["position"])))
    if limit:
        keep = keep[:limit]

    # one prepared PDB per construct
    ranges = {}
    for _, r in keep:
        lo = int(r["construct_range"].split("-")[0])
        ranges[r["construct"]] = lo
    prepared = {}
    for name, lo in sorted(ranges.items()):
        p = prepare_pdb(name, lo)
        if p is None:
            log(f"  could not prepare {name} — its variants are skipped")
            continue
        prepared[name] = p[0]

    tasks = []
    for why, r in keep:
        if r["construct"] not in prepared:
            continue
        tasks.append(dict(
            variation_id=r["variation_id"], position=int(r["position"]),
            wt_aa=r["wt_aa"], mut_aa=r["mut_aa"], set_primary=r["set_primary"],
            group=why, construct=r["construct"], pdb=str(prepared[r["construct"]]),
        ))
    return tasks


def main() -> int:
    ap = argparse.ArgumentParser()
    # Default 1, not 3. The repacking shell is small enough that simulated annealing reaches the
    # same optimum from every seed: over five seeds the ddG is identical to the third decimal,
    # while the underlying RNG streams verifiably differ (0.398/0.716/0.581 as first draws). The
    # extra repeats therefore buy nothing, and an sd of 0.000 computed from them would read as a
    # stability result when it only restates that the protocol is convergent.
    ap.add_argument("--repeats", type=int, default=1,
                    help="repeat count; the protocol converges, so >1 mainly costs time")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--calibrate", action="store_true")
    ap.add_argument("--protocol", choices=["torsion", "cart"], default="torsion",
                    help="minimisation space; cart is unsafe with metals (see init_worker)")
    ap.add_argument("--no-metals", dest="metals", action="store_false",
                    help="disable -in:auto_setup_metals (control for the metal terms)")
    args = ap.parse_args()

    log(f"Phase 5 revision — Rosetta ddG — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  {'ref2015_cart' if args.protocol == 'cart' else 'ref2015'}; "
        f"metal bonds {'on' if args.metals else 'off'} (rebuilt per pass, no harmonic "
        f"constraints); {args.protocol} space; {SHELL_A} A repack/min shell; "
        f"{args.repeats} repeats; {args.jobs} workers")

    # relaxed_path() reads these, and it is called in this process as well as in the workers.
    _PYR["protocol"], _PYR["metals"] = args.protocol, args.metals

    tasks = build_tasks(8 if args.calibrate else args.limit)
    if not tasks:
        log("FAIL: no tasks built")
        return 1
    counts = defaultdict(int)
    for t in tasks:
        counts[t["group"]] += 1
    log(f"  {len(tasks)} mutations: {dict(counts)}")
    log(f"  {len({t['pdb'] for t in tasks})} prepared structures in {WORK}")
    log()

    # Relax each construct once, up front, so the per-mutation cost is only repack + minimise.
    pdbs = sorted({t["pdb"] for t in tasks})
    todo = [p for p in pdbs if not relaxed_path(p).exists()]
    if todo:
        log(f"  pre-relaxing {len(todo)} construct(s) — one-off, cached to {WORK}")
        tr = time.time()
        with Pool(min(args.jobs, len(todo)), initializer=init_worker,
                  initargs=(args.repeats, args.protocol, args.metals)) as pool:
            for p in pool.imap_unordered(prerelax, todo):
                log(f"    relaxed {Path(p).name}")
        log(f"  pre-relaxation took {time.time()-tr:.0f}s")
    else:
        log(f"  all {len(pdbs)} construct(s) already relaxed in {WORK}")
    log()

    t0 = time.time()
    results = []
    with Pool(args.jobs, initializer=init_worker,
              initargs=(args.repeats, args.protocol, args.metals)) as pool:
        for i, res in enumerate(pool.imap_unordered(run_one, tasks), start=1):
            results.append(res)
            if args.calibrate or i % 25 == 0:
                log(f"  [{i}/{len(tasks)}] {res['wt_aa']}{res['position']}{res['mut_aa']} "
                    f"{res['group']:24s} {res['status']:20s} "
                    f"ddG {res.get('ddg_mean','')} ({res.get('seconds','?')}s)")

    el = time.time() - t0
    ok = [r for r in results if r["status"] == "ok"]
    log()
    log(f"  {len(ok)}/{len(results)} succeeded in {el:.0f}s "
        f"({el/max(len(results),1):.1f}s per mutation wall-clock at {args.jobs} workers)")
    bad = defaultdict(int)
    for r in results:
        if r["status"] != "ok":
            bad[r["status"].split(":")[0]] += 1
    if bad:
        log(f"  failures: {dict(bad)}")

    if args.calibrate:
        per = st.mean([r["seconds"] for r in ok]) if ok else float("nan")
        log()
        log(f"  CALIBRATION ONLY — nothing written.")
        log(f"  mean {per:.1f}s per mutation per worker; at {args.jobs} workers that is "
            f"{per/args.jobs:.1f}s each.")
        full = len(build_tasks(None))
        log(f"  full set is {full} mutations -> approximately "
            f"{full*per/args.jobs/3600:.1f} h")
        return 0

    if ok:
        for r in ok:
            grp = defaultdict(list)
        by = defaultdict(list)
        for r in ok:
            by[r["group"]].append(r["ddg_mean"])
        log()
        for g, v in sorted(by.items()):
            log(f"  {g:26s} n={len(v):4d}  median ddG {st.median(v):+.3f}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["variation_id", "position", "wt_aa", "mut_aa", "set_primary", "group", "construct",
            "status", "ddg_mean", "ddg_kcal_scaled", "ddg_sd", "n_repeats", "seconds"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    log(f"\nwrote {OUT.relative_to(ROOT)}  ({len(results)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
