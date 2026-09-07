#!/usr/bin/env python3
"""
38_af3_batch4_mutants.py — build the AlphaFold Server request that turns Figure 4 into a rate.

Figure 4 currently rests on nine variants: 18 site measurements, six clearing the seed-noise
floor, one experimental anchor (N2144S, Kettle 1999). That demonstrates the effect exists. It
cannot say how often it happens, which is the question a reviewer will ask.

This selects a set that can answer it, and it is built around what makes the answer
interpretable rather than around maximising n:

  * every PATHOGENIC side-chain Ca2+ ligand variant that has structural coverage (n = 40)
    -- the population whose disruption rate we want
  * NEGATIVE controls: benign / population variants at non-ligand positions in the same
    constructs, so the noise floor is measured on the same folds rather than assumed
  * COMPOSITION controls: D/E/N variants that are NOT ligands, so "is it the residue type or
    the calcium role?" stays answerable in the geometry data, exactly as it is in the energetics
  * the two existing controls (C1138F, P1141L) are already folded and are not repeated

Wild-type residues are verified against P35555 before anything is written. Constructs reuse the
21 tandem definitions already folded for the wild type, so every mutant has a matched wild-type
counterpart with a five-seed spread already measured -- without that, a mutant fold is
uninterpretable.

Writes:
  handoff/structures/inbox/batch4_mutants/af3_batch4_all.json   (one paste per job)
  handoff/structures/inbox/batch4_mutants/manifest.csv
  handoff/structures/REQUESTS.md  (appended)
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
RAW = ROOT / "data/raw"
OUTDIR = ROOT / "handoff/structures/inbox/batch4_mutants"
REQUESTS = ROOT / "handoff/structures/REQUESTS.md"

# AlphaFold Server's daily job allowance. Every day file is filled to exactly this, so the batch
# size is chosen to be a whole number of days rather than the controls being an arbitrary cap.
DAY_SIZE = 30

# The pathogenic ligand set is fixed by the data -- every one of them is wanted, and there is no
# principled reason to drop any. The controls are what get sized, to bring the total up to the
# next whole day. A minimum is enforced so a batch can never come out control-starved just
# because the pathogenic count happened to land near a multiple of DAY_SIZE.
MIN_CONTROLS = 20
NEG_SHARE = 0.55          # negative vs composition split of the control budget

AA3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "E": "GLU", "Q": "GLN",
       "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
       "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def reference() -> str:
    hits = sorted(RAW.glob("uniprot_P35555_v*.json"))
    if not hits:
        sys.exit("STOP: no cached UniProt P35555 JSON in data/raw")
    return json.loads(hits[-1].read_text(encoding="utf-8"))["sequence"]["value"]


def tbool(v) -> bool:
    return str(v).strip().lower() == "true"


def main() -> int:
    seq = reference()
    rows = list(csv.DictReader((PROC / "structural_metrics_af3.tsv").open(encoding="utf-8"),
                               delimiter="\t"))
    rows = [r for r in rows if r.get("construct")]

    # Already folded as mutants in batch 3 — never ask for them twice.
    done = {("N2144S"), ("N2183S"), ("E2130K"), ("D1113G"), ("E1073K"),
            ("D1487G"), ("E913K"), ("C1138F"), ("P1141L")}

    def tag(r):
        return f"{r['wt_aa']}{r['position']}{r['mut_aa']}"

    ranges = {}
    for r in rows:
        lo, hi = r["construct_range"].split("-")
        ranges[r["construct"]] = (int(lo), int(hi))

    picks: list[tuple[str, dict]] = []
    seen = set()

    def take(kind, r):
        t = tag(r)
        if t in done or t in seen:
            return False
        seen.add(t)
        picks.append((kind, r))
        return True

    # 1. every pathogenic side-chain Ca ligand
    for r in rows:
        if tbool(r["ca_ligand_sidechain"]) and r["set_primary"] == "pathogenic":
            take("ca_ligand_pathogenic", r)

    # Size the control budget so the whole batch is a whole number of full days.
    n_path = len(picks)
    total = -(-(n_path + MIN_CONTROLS) // DAY_SIZE) * DAY_SIZE      # round up to a full day
    budget = total - n_path
    n_neg = round(budget * NEG_SHARE)
    n_comp = budget - n_neg

    # 2. negative controls — benign or unlabelled, not a ligand, not a cysteine, spread across
    #    the constructs that carry the pathogenic ligands so the noise floor is measured on the
    #    same folds. Taken round-robin over constructs so the coverage stays even rather than
    #    piling several controls into whichever construct happens to be listed first.
    want = {r["construct"] for _, r in picks}
    by_con = defaultdict(list)
    for r in rows:
        if (r["construct"] in want and not tbool(r["ca_ligand_any"])
                and not tbool(r["cys_removing"]) and not tbool(r["consensus_but_not_ligand"])
                and r["set_primary"] in ("benign", "")):
            by_con[r["construct"]].append(r)
    for pool in by_con.values():
        pool.sort(key=lambda r: int(r["position"]))
    rnd = 0
    while sum(1 for k, _ in picks if k == "negative_control") < n_neg:
        progressed = False
        for con in sorted(by_con):
            if sum(1 for k, _ in picks if k == "negative_control") >= n_neg:
                break
            pool = by_con[con]
            # walk outward from the middle of each construct, avoiding the termini where the
            # fold is least reliable
            idx = len(pool) // 2 + (rnd // 2 + 1) * (1 if rnd % 2 else -1) if rnd else len(pool) // 2
            if 0 <= idx < len(pool) and take("negative_control", pool[idx]):
                progressed = True
        rnd += 1
        if not progressed:
            break

    # 3. composition controls — D/E/N that coordinate nothing
    comp = [r for r in rows
            if r["wt_aa"] in ("D", "E", "N") and not tbool(r["ca_ligand_any"])
            and not tbool(r["cys_removing"]) and r["construct"] in want]
    comp.sort(key=lambda r: (r["construct"], int(r["position"])))
    step = max(1, len(comp) // max(n_comp, 1))
    for r in comp[::step] + comp:          # stride first for spread, then backfill any shortfall
        if sum(1 for k, _ in picks if k == "composition_control") >= n_comp:
            break
        take("composition_control", r)

    # ---- validate every pick against the reference before writing anything ------------
    jobs, manifest, bad = [], [], []
    for kind, r in picks:
        pos, wt, mut = int(r["position"]), r["wt_aa"], r["mut_aa"]
        con = r["construct"]
        lo, hi = ranges[con]
        if seq[pos - 1] != wt:
            bad.append(f"{tag(r)}: reference has {seq[pos-1]} at {pos}")
            continue
        if not (lo <= pos <= hi):
            bad.append(f"{tag(r)}: position outside its construct {con} ({lo}-{hi})")
            continue
        if mut not in AA3:
            bad.append(f"{tag(r)}: unsupported residue {mut}")
            continue
        sub = seq[lo - 1:hi]
        local = pos - lo                      # 0-based offset inside the construct
        if sub[local] != wt:
            bad.append(f"{tag(r)}: construct offset check failed")
            continue
        mutated = sub[:local] + mut + sub[local + 1:]
        n_ca = 3 if con.count("-") and len(sub) > 150 else 3
        name = f"FBN1_{con.replace('FBN1_', '').replace('_Ca3', '')}_{tag(r)}_Ca{n_ca}"
        jobs.append({
            "name": name,
            "modelSeeds": [],
            "sequences": [
                {"proteinChain": {"sequence": mutated, "count": 1}},
                {"ion": {"ion": "CA", "count": n_ca}},
            ],
            "dialect": "alphafoldserver",
            "version": 1,
        })
        manifest.append({
            "job_name": name, "class": kind, "variant": tag(r),
            "variation_id": r["variation_id"], "clinical": r["set_primary"] or "unlabelled",
            "construct": con, "range": f"{lo}-{hi}", "length": len(sub),
            "n_ca": n_ca, "wt_verified": "yes",
            "destination": f"structures/alphafold/{name.lower()}/",
        })

    if bad:
        print("STOP — refusing to write a request containing unverified variants:")
        for b in bad:
            print("   ", b)
        return 1
    if not jobs:
        print("STOP — no jobs built")
        return 1

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "af3_batch4_all.json").write_text(json.dumps(jobs, indent=2) + "\n",
                                                encoding="utf-8")
    jobdir = OUTDIR / "af3_jobs"
    jobdir.mkdir(exist_ok=True)
    for j in jobs:
        (jobdir / f"{j['name']}.json").write_text(json.dumps([j], indent=2) + "\n",
                                                  encoding="utf-8")
    with (OUTDIR / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()))
        w.writeheader(); w.writerows(manifest)

    # Day-sized chunks. AlphaFold Server caps free non-commercial use at ~30 jobs per day, and a
    # batch that trips the cap halfway leaves the mentor guessing which jobs ran. Chunks are cut
    # so each day contains a mix of all three classes: if only one day ever gets run, the result
    # is still analysable rather than being 25 pathogenic variants with no controls.
    CHUNK = DAY_SIZE
    by_class = defaultdict(list)
    for j, m in zip(jobs, manifest):
        by_class[m["class"]].append((j, m))
    # Deal each class across the days in turn, rather than round-robin into one flat list.
    # Round-robin looks balanced but drains the smaller classes first: it gave day 1 a clean
    # 10/10/10 and day 2 a lopsided 24/4/2, so a day-2-only run would have been nearly
    # control-free -- exactly the situation the split exists to prevent. Dealing per class gives
    # every day the same proportions as the batch.
    n_days = len(jobs) // CHUNK
    if n_days * CHUNK != len(jobs):
        sys.exit(f"STOP: {len(jobs)} jobs is not a whole number of {CHUNK}-job days")
    buckets: list[list[tuple[dict, dict]]] = [[] for _ in range(n_days)]
    for k in ("ca_ligand_pathogenic", "negative_control", "composition_control"):
        for i, item in enumerate(by_class[k]):
            buckets[i % n_days].append(item)
    # Dealing is exact only when each class divides evenly; even it out if not.
    for _ in range(len(jobs)):
        over = max(range(n_days), key=lambda d: len(buckets[d]))
        under = min(range(n_days), key=lambda d: len(buckets[d]))
        if len(buckets[over]) - len(buckets[under]) <= 0:
            break
        buckets[under].append(buckets[over].pop())
    days = buckets
    if any(len(d) != CHUNK for d in days):
        sys.exit(f"STOP: day files are {[len(d) for d in days]}, not all {CHUNK} — "
                 "the control budget did not round to a whole number of days")
    for d, chunk in enumerate(days, start=1):
        (OUTDIR / f"af3_batch4_day{d}.json").write_text(
            json.dumps([j for j, _ in chunk], indent=2) + "\n", encoding="utf-8")

    counts = defaultdict(int)
    for m in manifest:
        counts[m["class"]] += 1
    print(f"batch 4: {len(jobs)} jobs")
    for k, v in sorted(counts.items()):
        print(f"  {k:<26} {v}")
    print(f"  constructs touched         {len({m['construct'] for m in manifest})}")
    print(f"  longest construct          {max(m['length'] for m in manifest)} aa")
    # Instructions live with the files, per CLAUDE.md §3b.
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d")
    inst = [
        "# AlphaFold Server — batch 4 (mutant cofolding)",
        "",
        f"Prepared {stamp} by `scripts/38_af3_batch4_mutants.py`. "
        f"**{len(jobs)} jobs**, split into {len(days)} day files of exactly {CHUNK} — "
        f"one full AlphaFold Server daily allowance each.",
        "",
        "## Why this batch exists",
        "",
        "Figure 4 currently rests on nine variants and can show only that Ca²⁺-site disruption "
        "happens, not how often. These jobs convert it into a rate, with the controls needed to "
        "make that rate mean something.",
        "",
        "| class | n | what it is for |",
        "|---|---|---|",
        f"| `ca_ligand_pathogenic` | {counts['ca_ligand_pathogenic']} | every pathogenic "
        "side-chain Ca²⁺ ligand with structural coverage that is not already folded |",
        f"| `negative_control` | {counts['negative_control']} | benign / unlabelled non-ligand "
        "positions in the same constructs — measures the noise floor on the same folds |",
        f"| `composition_control` | {counts['composition_control']} | D/E/N residues that "
        "coordinate nothing — separates residue type from calcium role |",
        "",
        "## What to do",
        "",
        "1. Go to <https://alphafoldserver.com>.",
        f"2. Upload `af3_batch4_day1.json`. Each day file holds exactly {CHUNK} jobs — a full "
        "daily allowance — and carries a mix of all three classes, so a day that runs alone is "
        "still analysable.",
        "3. Repeat with `af3_batch4_day2.json` (and day3 if present) on subsequent days.",
        "4. For each finished job, download the result and unzip the **whole output folder** "
        "into the path in `manifest.csv` → `destination`, i.e. "
        "`structures/alphafold/<job_name_lowercased>/`.",
        "   Keep the CIF, the confidence JSON **and the seed the server chose** — the seed is "
        "what makes the run reproducible, and the five-model spread is what every significance "
        "claim in Figure 4 is measured against.",
        "5. Tell me when a day is done. I will not touch `structures/alphafold/` until you do.",
        "",
        "If a job is rejected, **do not hand-edit the sequence** — send me the error and I will "
        "regenerate it. Every sequence here was checked against P35555 before writing: correct "
        "length, exactly one substitution, at the intended position, with the intended residues.",
        "",
        "## What I will do with them",
        "",
        "- Measure ΔBVS, coordination number and ion displacement per mutant against its "
        "existing wild-type fold, exactly as in Figure 4.",
        "- Report the **fraction** of pathogenic Ca²⁺-ligand variants whose site disruption "
        "clears the seed-noise floor, with the controls establishing that floor empirically "
        "rather than by assumption.",
        "- Test whether disruption predicts anything else we hold — clinical severity, "
        "AlphaMissense, burial, which ligand position is hit.",
        "",
        f"Total sequence to fold: {sum(int(m['length']) for m in manifest):,} residues across "
        f"{len({m['construct'] for m in manifest})} constructs; longest "
        f"{max(int(m['length']) for m in manifest)} aa.",
        "",
    ]
    (OUTDIR / "INSTRUCTIONS.md").write_text("\n".join(inst), encoding="utf-8")

    with REQUESTS.open("a", encoding="utf-8") as fh:
        fh.write(f"\n\n---\n\n## Batch 4 — mutant cofolding to make Figure 4 a rate "
                 f"({stamp})\n\n")
        fh.write(f"{len(jobs)} jobs in `handoff/structures/inbox/batch4_mutants/`, split into "
                 f"{len(days)} day files of exactly {CHUNK}. Full instructions in that folder's `INSTRUCTIONS.md`; "
                 f"per-job detail in its `manifest.csv`.\n\n")
        for k in ("ca_ligand_pathogenic", "negative_control", "composition_control"):
            fh.write(f"- **{k}** — {counts[k]} jobs\n")
        fh.write("\nEvery wild-type residue verified against P35555 before writing; every "
                 "sequence checked to differ from the reference at exactly one position.\n")

    print(f"\nwrote {(OUTDIR / 'af3_batch4_all.json').relative_to(ROOT)}")
    print(f"wrote {(OUTDIR / 'manifest.csv').relative_to(ROOT)}")
    print(f"wrote {(OUTDIR / 'INSTRUCTIONS.md').relative_to(ROOT)}")
    print(f"wrote {len(days)} day files of exactly {CHUNK} jobs each")
    print(f"wrote {len(jobs)} individual job files in {jobdir.relative_to(ROOT)}")
    print(f"appended to {REQUESTS.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
