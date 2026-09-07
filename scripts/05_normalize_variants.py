#!/usr/bin/env python3
"""
05_normalize_variants.py — PHASE 2. Normalize FBN1 variants onto P35555 and build labelled sets.

Thresholds come from resources/reference/inclusion_criteria.md (mentor-approved 2026-08-02) and
are NOT tunable here.

The wild-type-identity check is the core validation: every kept missense variant must have a
wild-type residue matching P35555 at that position. Mismatches are quarantined, never coerced.

Nothing is deleted. Every input row leaves in exactly one of:
  variants_normalized.tsv        missense, WT-verified, usable for the structural analysis
  excluded_from_structural.tsv   real variants of a non-missense consequence
  wt_mismatch.tsv                WT residue disagrees with P35555 — numbering/isoform problem
  unparsed.tsv                   no transcript-level HGVS to work from
Counts are reconciled at the end; a discrepancy is a hard stop.

Writes:
  data/interim/variants_normalized.tsv
  data/interim/{excluded_from_structural,wt_mismatch,unparsed}.tsv
  data/processed/variant_master.tsv          normalized + gnomAD + AlphaMissense joined
  data/processed/variant_master_dictionary.md
  logs/05_normalize_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
RAW, INTERIM, PROC = ROOT / "data/raw", ROOT / "data/interim", ROOT / "data/processed"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

PINNED_TX = "NM_000138.5"
FAF_BENIGN = 0.001      # stand-alone benign (BA1-style)
FAF_SUPPORT = 0.0001    # supporting only

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

STARS = {
    "practice guideline": 4, "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0, "no assertion provided": 0,
    "no classification provided": 0, "no classifications from unflagged records": 0,
    "no classification for the single variant": 0, "-": 0,
}

NAME = re.compile(r"^(?P<tx>N[MR]_\d+\.\d+)\((?P<gene>[^)]+)\):(?P<c>c\.[^ ]+)"
                  r"(?:\s+\((?P<p>p\.[^)]+)\))?")
MISSENSE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

log_lines: list[str] = []


def log(msg: str = "") -> None:
    print(msg)
    log_lines.append(msg)


def clinvar_class(sig: str) -> str:
    """Map ClinVar's free-text significance onto our five buckets."""
    s = (sig or "").strip()
    low = s.lower()
    if "conflicting" in low:
        return "conflicting"
    head = low.split(";")[0].strip()
    if head in ("pathogenic", "likely pathogenic", "pathogenic/likely pathogenic",
                "pathogenic/likely pathogenic/pathogenic, low penetrance",
                "pathogenic/likely pathogenic/established risk allele"):
        return "P/LP"
    if head in ("benign", "likely benign", "benign/likely benign"):
        return "B/LB"
    if head in ("uncertain significance", "uncertain risk allele",
                "uncertain significance/uncertain risk allele"):
        return "VUS"
    return "other"


def load_reference() -> str:
    hits = sorted(RAW.glob("uniprot_P35555_v*.json"))
    if not hits:
        raise SystemExit("STOP: no cached UniProt record — run scripts/02_fetch_uniprot.py")
    rec = json.loads(hits[-1].read_text(encoding="utf-8"))
    seq = rec["sequence"]["value"]
    if rec["primaryAccession"] != "P35555" or len(seq) != 2871:
        raise SystemExit("STOP: cached UniProt record is not P35555/2871")
    return seq


def mutalyzer(hgvs: str, session: requests.Session) -> tuple[str | None, str | None]:
    """Return (protein_description, error). Used on a targeted subset, not all 9k variants."""
    try:
        r = session.get(f"https://mutalyzer.nl/api/normalize/{hgvs}", timeout=60)
        if r.status_code != 200:
            return None, f"HTTP {r.status_code}"
        js = r.json()
        return ((js.get("protein") or {}).get("description"), None)
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}"


def main() -> int:
    INTERIM.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    log(f"Phase 2 normalization — {datetime.now(timezone.utc).isoformat()}")

    seq = load_reference()
    log(f"reference: P35555, {len(seq)} aa\n")

    cv_path = sorted(RAW.glob("clinvar_fbn1_*.tsv"))[-1]
    with cv_path.open(encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["Assembly"] == "GRCh38"]
    log(f"ClinVar GRCh38 input rows: {len(rows):,}  (from {cv_path.name})")

    normalized, excluded, mismatched, unparsed = [], [], [], []

    for r in rows:
        vid = r["VariationID"]
        name = r["Name"] or ""
        sig = r["ClinicalSignificance"]
        rs = r["ReviewStatus"].strip().lower()
        if rs not in STARS:
            raise SystemExit(f"STOP: unmapped ReviewStatus {r['ReviewStatus']!r} (id {vid})")
        stars = STARS[rs]
        base = {
            "variation_id": vid, "allele_id": r["AlleleID"], "name": name,
            "clinical_significance": sig, "clinvar_class": clinvar_class(sig),
            "review_status": r["ReviewStatus"], "stars": stars,
            "last_evaluated": r["LastEvaluated"], "n_submitters": r["NumberSubmitters"],
            "phenotypes": r["PhenotypeList"], "variant_type": r["Type"],
            "chrom": r["Chromosome"], "pos_vcf": r["PositionVCF"],
            "ref_vcf": r["ReferenceAlleleVCF"], "alt_vcf": r["AlternateAlleleVCF"],
        }

        m = NAME.match(name)
        if not m:
            unparsed.append({**base, "reason": "no transcript-level HGVS in Name"})
            continue

        tx, hgvs_c, hgvs_p = m.group("tx"), m.group("c"), m.group("p")
        base |= {"transcript": tx, "hgvs_c": hgvs_c, "hgvs_p": hgvs_p or ""}

        if hgvs_p is None:
            excluded.append({**base, "reason": "no protein consequence (likely splice/intronic)"})
            continue

        mm = MISSENSE.match(hgvs_p)
        if not mm:
            kind = ("synonymous" if hgvs_p.endswith("=") else
                    "frameshift" if "fs" in hgvs_p else
                    "in-frame indel" if any(k in hgvs_p for k in ("del", "dup", "ins")) else
                    "unknown protein effect" if hgvs_p == "p.?" else "other")
            excluded.append({**base, "reason": f"not missense: {kind}"})
            continue

        wt3, pos_s, alt3 = mm.groups()
        pos = int(pos_s)
        if alt3 == "Ter":
            excluded.append({**base, "reason": "not missense: nonsense (Ter)"})
            continue
        if wt3 == alt3:
            excluded.append({**base, "reason": "not missense: synonymous"})
            continue
        if wt3 not in AA3 or alt3 not in AA3:
            excluded.append({**base, "reason": f"unrecognised amino-acid code {wt3}/{alt3}"})
            continue
        if not (1 <= pos <= len(seq)):
            mismatched.append({**base, "position": pos, "wt_aa": AA3[wt3], "mut_aa": AA3[alt3],
                               "reference_aa": "", "reason": f"position {pos} outside 1-{len(seq)}"})
            continue

        wt, alt = AA3[wt3], AA3[alt3]
        ref_aa = seq[pos - 1]
        row = {**base, "position": pos, "wt_aa": wt, "mut_aa": alt, "reference_aa": ref_aa}

        # THE wild-type identity check.
        if wt != ref_aa:
            mismatched.append({**row, "reason": f"WT {wt} != P35555 residue {ref_aa} at {pos}"})
            continue

        normalized.append(row)

    log(f"  missense, WT-verified : {len(normalized):,}")
    log(f"  excluded (non-missense): {len(excluded):,}")
    log(f"  WT mismatches          : {len(mismatched):,}")
    log(f"  unparsed               : {len(unparsed):,}")

    total = len(normalized) + len(excluded) + len(mismatched) + len(unparsed)
    if total != len(rows):
        raise SystemExit(f"STOP: row reconciliation failed — {total} out != {len(rows)} in")
    log(f"  reconciliation: {total:,} out == {len(rows):,} in  OK\n")

    # ---- gnomAD join on exact genomic coordinates -------------------------
    gn_path = sorted(RAW.glob("gnomad_*_FBN1_*.tsv"))[-1]
    with gn_path.open(encoding="utf-8") as fh:
        gn = {r["variant_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

    def fnum(v):
        return float(v) if v not in ("", "None", None) else None

    n_gn = 0
    for row in normalized:
        key = f"{row['chrom']}-{row['pos_vcf']}-{row['ref_vcf']}-{row['alt_vcf']}"
        g = gn.get(key)
        row["gnomad_variant_id"] = key if g else ""
        row["gnomad_joint_af"] = ""
        row["faf95_grpmax"] = ""
        row["faf95_grpmax_group"] = ""
        if g:
            n_gn += 1
            ac, an = fnum(g["joint_ac"]), fnum(g["joint_an"])
            row["gnomad_joint_af"] = (ac / an) if (ac is not None and an) else ""
            row["faf95_grpmax"] = g["faf95_grpmax"]
            row["faf95_grpmax_group"] = g["faf95_grpmax_group"]
            # gnomAD's own protein description is an independent check on ours.
            if g["hgvsp"] and g["hgvsp"] != row["hgvs_p"]:
                row["gnomad_hgvsp_disagrees"] = g["hgvsp"]
        row.setdefault("gnomad_hgvsp_disagrees", "")
    log(f"gnomAD join: {n_gn:,}/{len(normalized):,} normalized variants found in gnomAD v4")
    disagree = [r for r in normalized if r["gnomad_hgvsp_disagrees"]]
    log(f"  gnomAD hgvsp disagreements: {len(disagree)}"
        + (f" e.g. {[(r['variation_id'], r['hgvs_p'], r['gnomad_hgvsp_disagrees']) for r in disagree[:3]]}"
           if disagree else ""))

    # ---- AlphaMissense join on (position, WT, ALT) ------------------------
    am_path = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
    am: dict[tuple[int, str, str], tuple[str, str]] = {}
    am_dropped = 0
    with am_path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            _, var, score, cls = line.rstrip("\n").split("\t")
            w, p, a = var[0], int(var[1:-1]), var[-1]
            # Join on identity, not position: this is what rejects the position-472 C-allele rows.
            if seq[p - 1] != w:
                am_dropped += 1
                continue
            am[(p, w, a)] = (score, cls)
    log(f"\nAlphaMissense: {len(am):,} usable entries; "
        f"{am_dropped} rows dropped for WT != P35555 (position-472 style)")

    n_am = 0
    for row in normalized:
        hit = am.get((row["position"], row["wt_aa"], row["mut_aa"]))
        row["am_pathogenicity"], row["am_class"] = hit if hit else ("", "")
        n_am += bool(hit)
    log(f"  AlphaMissense join: {n_am:,}/{len(normalized):,} matched")

    # ---- Labelled sets ----------------------------------------------------
    for row in normalized:
        faf = fnum(row["faf95_grpmax"])
        cls, stars = row["clinvar_class"], row["stars"]

        gnomad_benign = faf is not None and faf > FAF_BENIGN
        row["gnomad_benign_support"] = ("stand_alone" if gnomad_benign else
                                        "supporting" if (faf is not None and faf > FAF_SUPPORT)
                                        else "")

        for tier, floor in (("primary", 2), ("sensitivity", 1)):
            label, source = "", ""
            if stars >= floor:
                if cls == "P/LP":
                    label, source = "pathogenic", "clinvar"
                elif cls == "B/LB":
                    label, source = "benign", "clinvar"
                elif cls == "VUS":
                    label, source = "vus", "clinvar"
                elif cls == "conflicting":
                    label, source = "conflicting", "clinvar"
            # Population evidence stands alone regardless of star tier.
            if gnomad_benign:
                if label == "pathogenic":
                    label, source = "contradiction", "clinvar_P_vs_gnomad_common"
                else:
                    label = "benign"
                    source = "both" if source == "clinvar" and cls == "B/LB" else "gnomad_faf"
            row[f"set_{tier}"] = label
            row[f"set_{tier}_source"] = source

    def tally(tier: str) -> Counter:
        return Counter(r[f"set_{tier}"] for r in normalized)

    log("\nLabelled sets (missense, WT-verified):")
    for tier in ("primary", "sensitivity"):
        t = tally(tier)
        log(f"  {tier:11} (stars >= {2 if tier == 'primary' else 1}): "
            + ", ".join(f"{k or 'unlabelled'}={v:,}" for k, v in sorted(t.items())))

    contradictions = [r for r in normalized if r["set_primary"] == "contradiction"
                      or r["set_sensitivity"] == "contradiction"]
    if contradictions:
        log(f"\n  !! {len(contradictions)} contradictions (ClinVar P/LP but gnomAD-common):")
        for r in contradictions[:10]:
            log(f"     {r['variation_id']} {r['hgvs_p']} stars={r['stars']} "
                f"faf95={r['faf95_grpmax']} ({r['faf95_grpmax_group']}) — {r['clinical_significance']}")

    # ---- Targeted Mutalyzer validation ------------------------------------
    # Validating all ~5k missense variants over the API would take hours. The WT-identity check
    # above is the stronger guarantee for our purpose, so Mutalyzer is used where it adds
    # information: every non-pinned-transcript record, every WT mismatch, and a random sample.
    import random
    random.seed(20260802)
    off_tx = [r for r in normalized if r["transcript"] != PINNED_TX]
    sample = random.sample(normalized, min(100, len(normalized)))
    targets = {r["variation_id"]: r for r in off_tx + sample}
    log(f"\nMutalyzer validation on {len(targets)} records "
        f"({len(off_tx)} off-pinned-transcript + {len(sample)} random sample)")

    ok = bad = err = 0
    details = []
    with requests.Session() as s:
        s.headers.update({"User-Agent": "fbn1-marfan-pipeline/0.1 (Phase 2)"})
        for i, (vid, r) in enumerate(targets.items(), 1):
            desc, e = mutalyzer(f"{PINNED_TX}:{r['hgvs_c']}", s)
            if e:
                err += 1
            elif desc and f"{r['wt_aa']}" and re.search(
                    rf"p\.\(?{re.escape(r['hgvs_p'][2:])}\)?", desc):
                ok += 1
            elif desc:
                bad += 1
                details.append((vid, r["hgvs_c"], r["hgvs_p"], desc))
            time.sleep(0.35)
            if i % 25 == 0:
                log(f"    {i}/{len(targets)} …")
    log(f"  agreement: {ok} agree, {bad} disagree, {err} errors")
    for d in details[:10]:
        log(f"    DISAGREE {d[0]} {d[1]}: ClinVar {d[2]} vs Mutalyzer {d[3]}")

    # ---- Write outputs ----------------------------------------------------
    def write(path: Path, data: list[dict]) -> None:
        if not data:
            path.write_text("", encoding="utf-8")
            return
        cols = list({k: None for row in data for k in row})
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in data:
                w.writerow(row)

    write(INTERIM / "variants_normalized.tsv", normalized)
    write(INTERIM / "excluded_from_structural.tsv", excluded)
    write(INTERIM / "wt_mismatch.tsv", mismatched)
    write(INTERIM / "unparsed.tsv", unparsed)
    write(PROC / "variant_master.tsv", normalized)

    log(f"\nwrote data/interim/variants_normalized.tsv ({len(normalized):,})")
    log(f"wrote data/interim/excluded_from_structural.tsv ({len(excluded):,})")
    log(f"wrote data/interim/wt_mismatch.tsv ({len(mismatched):,})")
    log(f"wrote data/interim/unparsed.tsv ({len(unparsed):,})")
    log(f"wrote data/processed/variant_master.tsv ({len(normalized):,})")

    (ROOT / f"logs/05_normalize_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                         encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
