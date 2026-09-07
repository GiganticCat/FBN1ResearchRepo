#!/usr/bin/env python3
"""
28_statistics_v2.py — PHASE 5 (revision). Statistics on the corrected definitions.

Supersedes results/phase5_statistics.*, which had two defects that pushed effects in the same
direction and so could not cancel.

  DEFINITION. Calcium-coordinating positions came from a sequence motif that was never validated,
  because UniProt annotates no calcium site in P35555. Against the batch-2 structures the motif
  over-calls by 43%. Every test below uses observed coordination; the 143 variants at
  motif-flagged non-ligand positions become their own group, which doubles as a check on the new
  definition -- if the structure-derived split is real, that group should behave like background
  rather than like ligands.

  COMPOSITION. Calcium ligands are aspartate, asparagine and glutamate by definition, and FoldX
  scores buried acidic side chains as unusually cheap to mutate regardless of whether they touch
  a metal. Comparing ligands against all other positions therefore measures amino-acid identity
  as much as calcium chemistry. Every comparison is reported twice: unmatched, as before, and
  stratified by wild-type residue, where the effect is estimated within each residue type and
  pooled. Where the two disagree, the stratified estimate is the one to believe.

Effect sizes are Cliff's delta with a bootstrap CI (seeded). Rank tests are Mann-Whitney;
2x2 tables are Fisher exact. Benjamini-Hochberg runs within each family.

The benign set is n=29 and is retained only because the mentor has not approved a replacement.
Any test involving it is flagged UNDERPOWERED in the output and should not carry a conclusion.

Writes:
  results/phase5_statistics_v2.md
  results/phase5_statistics_v2.tsv
  logs/28_statistics_v2_<stamp>.log
"""

from __future__ import annotations

import csv
import datetime as dt
import math
import random
import statistics as st
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SM = ROOT / "data" / "processed" / "structural_metrics_af3.tsv"
ABL = ROOT / "data" / "processed" / "foldx_ion_ablation.tsv"
OUT_MD = ROOT / "results" / "phase5_statistics_v2.md"
OUT_TSV = ROOT / "results" / "phase5_statistics_v2.tsv"

SEED = 20260809
N_BOOT = 2000
MIN_STRATUM = 3          # per-residue stratum needs this many on each side to contribute

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"28_statistics_v2_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def cliffs_delta(a: list[float], b: list[float]) -> float:
    """P(a>b) - P(a<b), computed by ranking rather than the O(nm) double loop."""
    if not a or not b:
        return float("nan")
    allv = np.concatenate([np.asarray(a, float), np.asarray(b, float)])
    r = stats.rankdata(allv)
    ra = r[: len(a)].sum()
    # Mann-Whitney U from ranks, then delta = 2U/(nm) - 1
    u = ra - len(a) * (len(a) + 1) / 2.0
    return 2.0 * u / (len(a) * len(b)) - 1.0


def boot_ci(a: list[float], b: list[float], rng: random.Random) -> tuple[float, float]:
    if len(a) < 3 or len(b) < 3:
        return (float("nan"), float("nan"))
    aa, bb = np.asarray(a, float), np.asarray(b, float)
    out = []
    for _ in range(N_BOOT):
        sa = aa[np.random.randint(0, len(aa), len(aa))]
        sb = bb[np.random.randint(0, len(bb), len(bb))]
        out.append(cliffs_delta(list(sa), list(sb)))
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def stratified_delta(pairs_a, pairs_b) -> tuple[float, int, int]:
    """Cliff's delta pooled across wild-type residue strata.

    pairs_* are (wt_residue, value). Returns (delta, n_strata, n_pairs_used).
    Weighting is by n_a*n_b within a stratum, which is the number of comparisons the stratum
    contributes -- the natural weight for a statistic defined over all pairs.
    """
    by_a, by_b = defaultdict(list), defaultdict(list)
    for wt, v in pairs_a:
        by_a[wt].append(v)
    for wt, v in pairs_b:
        by_b[wt].append(v)
    num = den = 0.0
    n_strata = 0
    for wt in sorted(set(by_a) & set(by_b)):
        a, b = by_a[wt], by_b[wt]
        if len(a) < MIN_STRATUM or len(b) < MIN_STRATUM:
            continue
        w = len(a) * len(b)
        num += w * cliffs_delta(a, b)
        den += w
        n_strata += 1
    if den == 0:
        return (float("nan"), 0, 0)
    return (num / den, n_strata, int(den))


def mw_p(a: list[float], b: list[float]) -> float:
    if len(a) < 3 or len(b) < 3:
        return float("nan")
    try:
        return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def bh(pvals: list[float]) -> list[float]:
    idx = [i for i, p in enumerate(pvals) if not math.isnan(p)]
    m = len(idx)
    if m == 0:
        return pvals[:]
    order = sorted(idx, key=lambda i: pvals[i])
    out = pvals[:]
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        val = min(prev, pvals[i] * m / k)
        out[i] = val
        prev = val
    return out


def fnum(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def main() -> int:
    random.seed(SEED)
    np.random.seed(SEED)
    rng = random.Random(SEED)

    log(f"Phase 5 revision — statistics v2 — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  seed {SEED}; {N_BOOT} bootstrap resamples; stratum minimum {MIN_STRATUM}")

    rows = list(csv.DictReader(open(SM), delimiter="\t"))
    log(f"  {len(rows)} variants with structural metrics")

    # FoldX ddG, available only for the four experimental templates
    ddg = {}
    for r in csv.DictReader(open(ABL), delimiter="\t"):
        if r["ablation_status"] == "ok":
            ddg[r["variation_id"]] = fnum(r["ddG_withca_mean"])
    log(f"  {len(ddg)} of those carry a FoldX ddG (experimental templates only)")

    def group(pred):
        return [r for r in rows if pred(r)]

    T = lambda r, k: r[k] == "True"                                        # noqa: E731
    lig = group(lambda r: T(r, "ca_ligand_sidechain"))
    bb_lig = group(lambda r: T(r, "ca_ligand_backbone_only"))
    overcalled = group(lambda r: T(r, "consensus_but_not_ligand"))
    cysrem = group(lambda r: T(r, "cys_removing"))
    other = group(lambda r: not T(r, "ca_ligand_any") and not T(r, "cys_removing")
                  and not T(r, "consensus_but_not_ligand"))
    DEN = {"D", "E", "N"}
    other_den = [r for r in other if r["wt_aa"] in DEN]

    log()
    log(f"  side-chain Ca ligands {len(lig)}; backbone-only ligands {len(bb_lig)}; "
        f"motif-overcalled {len(overcalled)}; cys-removing {len(cysrem)}; other {len(other)} "
        f"(of which D/E/N {len(other_den)})")

    results = []

    def test(family, label, A, B, key, note=""):
        a = [(r["wt_aa"], fnum(r[key])) for r in A if not math.isnan(fnum(r[key]))]
        b = [(r["wt_aa"], fnum(r[key])) for r in B if not math.isnan(fnum(r[key]))]
        av, bv = [v for _, v in a], [v for _, v in b]
        d = cliffs_delta(av, bv)
        lo, hi = boot_ci(av, bv, rng)
        sd, ns, npair = stratified_delta(a, b)
        p = mw_p(av, bv)
        results.append(dict(
            family=family, test=label, n_a=len(av), n_b=len(bv),
            median_a=round(st.median(av), 4) if av else "",
            median_b=round(st.median(bv), 4) if bv else "",
            cliffs_delta=round(d, 3) if not math.isnan(d) else "",
            ci_lo=round(lo, 3) if not math.isnan(lo) else "",
            ci_hi=round(hi, 3) if not math.isnan(hi) else "",
            delta_wt_matched=round(sd, 3) if not math.isnan(sd) else "",
            n_strata=ns, p_raw=p, note=note,
        ))

    def test_ddg(family, label, A, B, note=""):
        A2 = [r for r in A if r["variation_id"] in ddg]
        B2 = [r for r in B if r["variation_id"] in ddg]
        a = [(r["wt_aa"], ddg[r["variation_id"]]) for r in A2]
        b = [(r["wt_aa"], ddg[r["variation_id"]]) for r in B2]
        av, bv = [v for _, v in a], [v for _, v in b]
        if len(av) < 3 or len(bv) < 3:
            results.append(dict(family=family, test=label, n_a=len(av), n_b=len(bv),
                                median_a="", median_b="", cliffs_delta="", ci_lo="", ci_hi="",
                                delta_wt_matched="", n_strata=0, p_raw=float("nan"),
                                note="too few with ddG"))
            return
        d = cliffs_delta(av, bv)
        lo, hi = boot_ci(av, bv, rng)
        sd, ns, _ = stratified_delta(a, b)
        results.append(dict(
            family=family, test=label, n_a=len(av), n_b=len(bv),
            median_a=round(st.median(av), 3), median_b=round(st.median(bv), 3),
            cliffs_delta=round(d, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
            delta_wt_matched=round(sd, 3) if not math.isnan(sd) else "",
            n_strata=ns, p_raw=mw_p(av, bv), note=note,
        ))

    # ---- Family A: AlphaMissense, full structural set
    test("A_alphamissense", "AlphaMissense: side-chain Ca ligands vs other positions",
         lig, other, "am_pathogenicity")
    test("A_alphamissense", "AlphaMissense: side-chain Ca ligands vs other D/E/N positions",
         lig, other_den, "am_pathogenicity", "composition control")
    test("A_alphamissense", "AlphaMissense: motif-overcalled positions vs other positions",
         overcalled, other, "am_pathogenicity",
         "should resemble background if the structural definition is right")
    test("A_alphamissense", "AlphaMissense: motif-overcalled vs true side-chain ligands",
         overcalled, lig, "am_pathogenicity", "definition check")
    test("A_alphamissense", "AlphaMissense: backbone-only ligands vs other positions",
         bb_lig, other, "am_pathogenicity",
         "no substitution can remove a backbone carbonyl")
    test("A_alphamissense", "AlphaMissense: cysteine-removing vs other positions",
         cysrem, other, "am_pathogenicity")

    # ---- Family B: FoldX ddG, experimental-template subset, the original claim re-tested
    test_ddg("B_foldx_ddg", "ddG: side-chain Ca ligands vs all other positions",
             lig, other, "the original comparison")
    test_ddg("B_foldx_ddg", "ddG: side-chain Ca ligands vs other D/E/N positions",
             lig, other_den, "composition-matched")
    test_ddg("B_foldx_ddg", "ddG: cysteine-removing vs other positions", cysrem, other)
    test_ddg("B_foldx_ddg", "ddG: buried vs exposed (RSA<20%)",
             [r for r in rows if T(r, "buried")], [r for r in rows if not T(r, "buried")])

    # ---- Family C: burial and geometry
    test("C_geometry", "RSA: side-chain Ca ligands vs other positions", lig, other, "rsa_pct")
    test("C_geometry", "distance to nearest Ca: pathogenic vs benign",
         [r for r in rows if r["set_primary"] == "pathogenic"],
         [r for r in rows if r["set_primary"] == "benign"],
         "dist_to_nearest_ca_A", "UNDERPOWERED: benign n=29")

    # BH within family
    fam_idx = defaultdict(list)
    for i, r in enumerate(results):
        fam_idx[r["family"]].append(i)
    for fam, idxs in fam_idx.items():
        adj = bh([results[i]["p_raw"] for i in idxs])
        for j, i in enumerate(idxs):
            results[i]["p_adj"] = adj[j]

    # ---- Family D: 2x2 enrichment, Fisher
    enrich = []

    def fisher(label, feat_pred, note=""):
        pat = [r for r in rows if r["set_primary"] == "pathogenic"]
        ben = [r for r in rows if r["set_primary"] == "benign"]
        a = sum(1 for r in pat if feat_pred(r))
        b = len(pat) - a
        c = sum(1 for r in ben if feat_pred(r))
        d = len(ben) - c
        try:
            orr, p = stats.fisher_exact([[a, b], [c, d]])
        except ValueError:
            orr, p = float("nan"), float("nan")
        enrich.append(dict(test=label, path_yes=a, path_no=b, benign_yes=c, benign_no=d,
                           odds_ratio=round(orr, 3) if not math.isnan(orr) else "",
                           p_raw=p, note=note))

    fisher("pathogenic vs benign at side-chain Ca ligands",
           lambda r: T(r, "ca_ligand_sidechain"), "UNDERPOWERED: benign n=29")
    fisher("pathogenic vs benign at cysteine-removing",
           lambda r: T(r, "cys_removing"), "UNDERPOWERED: benign n=29")
    fisher("pathogenic vs benign at motif-overcalled positions",
           lambda r: T(r, "consensus_but_not_ligand"), "UNDERPOWERED: benign n=29")
    adj = bh([e["p_raw"] for e in enrich])
    for e, p in zip(enrich, adj):
        e["p_adj"] = p

    # ---- output
    log()
    log("Rank comparisons (delta_matched is the wild-type-stratified estimate)")
    log(f"  {'family':17s} {'n_a/n_b':>10s} {'delta':>7s} {'matched':>8s} {'p_adj':>10s}  test")
    for r in results:
        log(f"  {r['family']:17s} {str(r['n_a'])+'/'+str(r['n_b']):>10s} "
            f"{str(r['cliffs_delta']):>7s} {str(r['delta_wt_matched']):>8s} "
            f"{r['p_adj']:>10.3g}  {r['test']}"
            + (f"   [{r['note']}]" if r["note"] else ""))
    log()
    log("Enrichment (Fisher)")
    for e in enrich:
        log(f"  OR {str(e['odds_ratio']):>8s}  p_adj {e['p_adj']:.3g}  {e['test']}"
            + (f"   [{e['note']}]" if e["note"] else ""))

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(results)
    log(f"\nwrote {OUT_TSV.relative_to(ROOT)}")

    with open(OUT_MD, "w") as fh:
        fh.write("# Statistics v2 — structure-derived definitions, composition-matched\n\n")
        fh.write(f"Generated by `scripts/28_statistics_v2.py` at "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}. Seed {SEED}, "
                 f"{N_BOOT} bootstrap resamples.\n\n")
        fh.write("Supersedes `phase5_statistics.md`. Calcium-coordinating positions are defined "
                 "by observed coordination in the AF3 models rather than by sequence motif, and "
                 "every comparison is reported both unmatched and stratified by wild-type "
                 "residue. **Where the two disagree, the stratified estimate is the one to "
                 "believe** — calcium ligands are D/E/N by definition, and that alone shifts "
                 "the unmatched comparison.\n\n")
        fh.write("| family | test | n | median A / B | Cliff's d | 95% CI | d (WT-matched) | strata | p (BH) |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in results:
            fh.write(f"| {r['family']} | {r['test']} | {r['n_a']}/{r['n_b']} | "
                     f"{r['median_a']} / {r['median_b']} | {r['cliffs_delta']} | "
                     f"[{r['ci_lo']}, {r['ci_hi']}] | **{r['delta_wt_matched']}** | "
                     f"{r['n_strata']} | {r['p_adj']:.3g} |"
                     + (f" {r['note']}" if r["note"] else "") + "\n")
        fh.write("\n## Enrichment\n\n")
        fh.write("| test | pathogenic yes/no | benign yes/no | OR | p (BH) | note |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for e in enrich:
            fh.write(f"| {e['test']} | {e['path_yes']}/{e['path_no']} | "
                     f"{e['benign_yes']}/{e['benign_no']} | {e['odds_ratio']} | "
                     f"{e['p_adj']:.3g} | {e['note']} |\n")
        fh.write("\nThe benign set is n=29. Every enrichment row above is underpowered and is "
                 "reported for completeness, not as a finding.\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
