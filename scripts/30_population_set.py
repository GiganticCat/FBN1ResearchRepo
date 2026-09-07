#!/usr/bin/env python3
"""
30_population_set.py — PHASE 5 (revision). A population comparison set that is not n=29.

The ClinVar benign set contains 29 variants with structural coverage and zero of them fall at a
calcium ligand or remove a cysteine, so every enrichment odds ratio in the study came back as
infinity. That is not a weak result, it is an absent one.

The replacement is the standard population-comparison design: FBN1 missense variants observed in
gnomAD v4 that are not classified pathogenic or likely pathogenic. Marfan syndrome is a dominant
disorder with prevalence about 1 in 10,000, so an allele seen in a population reference and never
called pathogenic is, in aggregate, depleted of severe alleles -- not proven benign individually,
but the right denominator for asking whether pathogenic variants concentrate anywhere.

Three things this deliberately does not claim. Individual population variants are not benign, and
are never labelled as such. The set is not independent of ClinVar, because the exclusion uses
ClinVar. And it carries the ascertainment bias of gnomAD itself. The comparison is reported
alongside the ClinVar-benign version rather than replacing it, so a reader can see exactly what
the choice of denominator does to the answer.

Wild-type identity is checked against P35555 for every parsed variant, the same gate applied
everywhere else in this pipeline; a variant whose reference residue disagrees is dropped and
counted, never coerced.

Writes:
  data/processed/population_set.tsv
  results/phase5_enrichment_v2.md
  logs/30_population_set_<stamp>.log
"""

from __future__ import annotations

import csv
import datetime as dt
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
GNOMAD = ROOT / "data" / "raw" / "gnomad_gnomad_r4_FBN1_2026-08-02.tsv"
FASTA = ROOT / "data" / "raw" / "uniprot_P35555_v264.fasta"
SM = ROOT / "data" / "processed" / "structural_metrics_af3.tsv"
VAR = ROOT / "data" / "processed" / "variant_annotated.tsv"
OUT = ROOT / "data" / "processed" / "population_set.tsv"
OUT_MD = ROOT / "results" / "phase5_enrichment_v2.md"

MANE_ENST = "ENST00000316623"          # MANE Select for FBN1, pairs with NM_000138.5

THREE2ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
HGVSP = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"30_population_set_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def reference() -> str:
    return "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))


def fnum(x, default=float("nan")):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main() -> int:
    log(f"Phase 5 revision — population comparison set — "
        f"{dt.datetime.now(dt.timezone.utc).isoformat()}")
    ref = reference()

    # ---- pathogenic protein changes to exclude
    path_changes, clinvar_class = set(), {}
    for r in csv.DictReader(open(VAR), delimiter="\t"):
        if not r["position"]:
            continue
        key = (int(r["position"]), r["wt_aa"], r["mut_aa"])
        clinvar_class[key] = r["set_primary"]
        if r["set_primary"] == "pathogenic":
            path_changes.add(key)
    log(f"  {len(path_changes)} pathogenic protein changes to exclude")

    # ---- parse gnomAD missense
    rows = list(csv.DictReader(open(GNOMAD), delimiter="\t"))
    n_mis = n_badtx = n_unparsed = n_wtmismatch = n_excluded = 0
    pop = []
    for r in rows:
        if r["consequence"] != "missense_variant":
            continue
        n_mis += 1
        if r["transcript_id"] != MANE_ENST:
            n_badtx += 1
            continue
        m = HGVSP.match((r["hgvsp"] or "").strip())
        if not m:
            n_unparsed += 1
            continue
        wt3, pos_s, mut3 = m.groups()
        if wt3 not in THREE2ONE or mut3 not in THREE2ONE:
            n_unparsed += 1
            continue
        wt, mut, pos = THREE2ONE[wt3], THREE2ONE[mut3], int(pos_s)
        if pos < 1 or pos > len(ref) or ref[pos - 1] != wt:
            n_wtmismatch += 1
            continue
        key = (pos, wt, mut)
        if key in path_changes:
            n_excluded += 1
            continue
        pop.append(dict(
            variant_id=r["variant_id"], position=pos, wt_aa=wt, mut_aa=mut,
            hgvsp=r["hgvsp"], joint_af=fnum(r.get("joint_af")),
            joint_ac=r.get("joint_ac", ""), joint_an=r.get("joint_an", ""),
            clinvar_status=clinvar_class.get(key, "not_in_clinvar"),
        ))

    log(f"  {n_mis} gnomAD missense rows; {n_badtx} off-transcript, {n_unparsed} unparsed HGVS, "
        f"{n_wtmismatch} WT-identity mismatches, {n_excluded} excluded as pathogenic")
    log(f"  {len(pop)} variants in the population set")
    log(f"  ClinVar status within it: {dict(Counter(p['clinvar_status'] for p in pop))}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pop[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(pop)
    log(f"wrote {OUT.relative_to(ROOT)}  ({len(pop)} rows)")

    # ---- join to structural features by position
    sm = list(csv.DictReader(open(SM), delimiter="\t"))
    feat_by_pos: dict[int, dict] = {}
    for r in sm:
        feat_by_pos.setdefault(int(r["position"]), r)

    T = lambda r, k: r[k] == "True"                                        # noqa: E731
    features = {
        "side-chain Ca ligand": lambda r: T(r, "ca_ligand_sidechain"),
        "backbone-only Ca ligand": lambda r: T(r, "ca_ligand_backbone_only"),
        "motif-overcalled (not a ligand)": lambda r: T(r, "consensus_but_not_ligand"),
        "cysteine position": lambda r: T(r, "is_cysteine"),
    }

    pathogenic = [r for r in sm if r["set_primary"] == "pathogenic"]
    clinvar_benign = [r for r in sm if r["set_primary"] == "benign"]
    population = [feat_by_pos[p["position"]] for p in pop if p["position"] in feat_by_pos]
    # a position may host several population variants; that is correct, each is an observation
    log()
    log(f"  with structural coverage: pathogenic {len(pathogenic)}, "
        f"ClinVar-benign {len(clinvar_benign)}, population {len(population)}")

    results = []
    for label, pred in features.items():
        for denom_name, denom in (("population (gnomAD)", population),
                                  ("ClinVar benign", clinvar_benign)):
            a = sum(1 for r in pathogenic if pred(r))
            b = len(pathogenic) - a
            c = sum(1 for r in denom if pred(r))
            d = len(denom) - c
            try:
                orr, p = stats.fisher_exact([[a, b], [c, d]])
            except ValueError:
                orr, p = float("nan"), float("nan")
            # Haldane-Anscombe correction only for display when a cell is zero
            orr_disp = orr
            if 0 in (a, b, c, d):
                orr_disp = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            # Woolf interval on the Haldane-corrected table, so a zero cell yields a finite
            # interval instead of dropping the row from a forest plot. Always computed on the
            # +0.5 table -- reporting a point estimate from one table and an interval from
            # another would be incoherent.
            A, B, C, D = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            se = math.sqrt(1 / A + 1 / B + 1 / C + 1 / D)
            or_h = (A * D) / (B * C)
            ci_lo = math.exp(math.log(or_h) - 1.96 * se)
            ci_hi = math.exp(math.log(or_h) + 1.96 * se)
            results.append(dict(
                feature=label, denominator=denom_name,
                path_yes=a, path_no=b, denom_yes=c, denom_no=d,
                pct_path=round(100 * a / max(len(pathogenic), 1), 1),
                pct_denom=round(100 * c / max(len(denom), 1), 1),
                odds_ratio=round(orr, 3) if not math.isnan(orr) and not math.isinf(orr) else "inf",
                odds_ratio_corrected=round(orr_disp, 3),
                ci_lo=round(ci_lo, 4), ci_hi=round(ci_hi, 4),
                zero_cell=0 in (a, b, c, d),
                p_raw=p,
            ))

    # BH across the whole family
    ps = [r["p_raw"] for r in results]
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    m = len(ps)
    prev = 1.0
    adj = [1.0] * m
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, ps[i] * m / k)
        adj[i] = prev
    for r, p in zip(results, adj):
        r["p_adj"] = p

    log()
    log(f"  {'feature':34s} {'denominator':20s} {'path%':>6s} {'den%':>6s} {'OR':>9s} {'p_adj':>10s}")
    for r in results:
        log(f"  {r['feature']:34s} {r['denominator']:20s} {r['pct_path']:6.1f} "
            f"{r['pct_denom']:6.1f} {str(r['odds_ratio_corrected']):>9s} {r['p_adj']:10.3g}"
            + ("  [zero cell — OR is Haldane-corrected]" if r["zero_cell"] else ""))

    # The enrichment figure reads this table. Without it the only machine-readable form of these
    # odds ratios is the prose below, and a figure would have to duplicate the contingency logic.
    OUT_TSV = ROOT / "data" / "processed" / "enrichment_v2.tsv"
    cols = ["feature", "denominator", "path_yes", "path_no", "denom_yes", "denom_no",
            "pct_path", "pct_denom", "odds_ratio_corrected", "ci_lo", "ci_hi",
            "zero_cell", "p_raw", "p_adj"]
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    log(f"wrote {OUT_TSV.relative_to(ROOT)}  ({len(results)} rows)")

    with open(OUT_MD, "w") as fh:
        fh.write("# Enrichment against a population denominator\n\n")
        fh.write(f"Generated by `scripts/30_population_set.py` at "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}.\n\n")
        fh.write(f"The ClinVar benign set has {len(clinvar_benign)} variants with structural "
                 f"coverage and no representation at calcium ligands or cysteines, which is why "
                 f"every odds ratio in the first analysis was infinite. The population "
                 f"denominator is {len(population)} gnomAD v4 missense observations on "
                 f"{MANE_ENST} that are not classified pathogenic, each checked for wild-type "
                 f"identity against P35555.\n\n")
        fh.write("**Population variants are not benign.** They are a denominator, not a control "
                 "group: individually unclassified, collectively depleted of severe alleles in a "
                 "dominant disorder. Both denominators are reported so the effect of the choice "
                 "is visible.\n\n")
        fh.write("| feature | denominator | pathogenic | denominator | OR | p (BH) |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for r in results:
            fh.write(f"| {r['feature']} | {r['denominator']} | "
                     f"{r['path_yes']}/{r['path_yes']+r['path_no']} ({r['pct_path']}%) | "
                     f"{r['denom_yes']}/{r['denom_yes']+r['denom_no']} ({r['pct_denom']}%) | "
                     f"{r['odds_ratio_corrected']}{' *' if r['zero_cell'] else ''} | "
                     f"{r['p_adj']:.3g} |\n")
        fh.write("\n`*` odds ratio is Haldane-Anscombe corrected because a cell was zero.\n\n")
        fh.write("Enrichment at these positions is partly definitional: the ClinGen FBN1 VCEP "
                 "applies PM1 at calcium-binding and cysteine residues, so ClinVar pathogenic "
                 "labels are not independent of the features being tested. These numbers "
                 "describe the labelled data; they are not evidence that the features cause "
                 "pathogenicity.\n")
    log(f"\nwrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
