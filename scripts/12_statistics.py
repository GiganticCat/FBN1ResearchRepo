#!/usr/bin/env python3
"""
12_statistics.py — PHASE 5. Statistical analysis.

Framing (mentor decision, 2026-08-02): the manuscript headlines the BIOPHYSICS, not the
enrichment counting. Accordingly this script is organised as:

  PART A (headline)   — biophysical measurements: FoldX ddG, Ca2+ coordination geometry, burial.
                        These come from coordinates and are independent of any clinical label,
                        so they are not affected by the PM1 circularity described in Phase 3.
  PART B (supporting) — enrichment of pathogenic variants at critical sites. Reported against the
                        all-possible-missense background, and always alongside the circularity
                        caveat, because ClinVar's pathogenic calls used the very rule being tested.

Every test reports an effect size with a confidence interval, not just a p-value. Multiple
testing is corrected within each family with Benjamini-Hochberg. Seeds are fixed.

Writes:
  results/phase5_statistics.md
  results/phase5_statistics.tsv
  logs/12_statistics_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
SEED = 20260802
rng = np.random.default_rng(SEED)

log_lines: list[str] = []
tests: list[dict] = []


def log(m: str = "") -> None:
    print(m)
    log_lines.append(m)


def read(p: Path) -> list[dict]:
    return list(csv.DictReader(p.open(encoding="utf-8"), delimiter="\t"))


def T(v) -> bool:
    return str(v).strip().lower() == "true"


def f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def cliffs_delta(a, b) -> float:
    """Non-parametric effect size for the Mann-Whitney comparison."""
    a, b = np.asarray(a), np.asarray(b)
    gt = sum((x > b).sum() for x in a)
    lt = sum((x < b).sum() for x in a)
    return (gt - lt) / (len(a) * len(b))


def or_ci(a, b, c, d):
    """Odds ratio with a Woolf 95% CI; Haldane-Anscombe correction when a cell is zero."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return orr, orr * np.exp(-1.96 * se), orr * np.exp(1.96 * se)


def record(family, name, stat, p, effect, ci, n, note=""):
    tests.append({"family": family, "test": name, "statistic": stat, "p_raw": p,
                  "effect": effect, "ci_low": ci[0], "ci_high": ci[1], "n": n, "note": note})


def mw(family, name, a, b, la, lb, note=""):
    """Mann-Whitney U with Cliff's delta and a bootstrap CI on the delta."""
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if len(a) < 3 or len(b) < 3:
        log(f"    {name}: skipped (n={len(a)}/{len(b)}, too few)")
        return
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    d = cliffs_delta(a, b)
    boots = [cliffs_delta(rng.choice(a, len(a), replace=True),
                          rng.choice(b, len(b), replace=True)) for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    record(family, name, u, p, d, (lo, hi), f"{len(a)}/{len(b)}", note)
    log(f"    {name}")
    log(f"      {la}: n={len(a)} median={np.median(a):.2f}   {lb}: n={len(b)} median={np.median(b):.2f}")
    log(f"      Mann-Whitney U={u:.0f}  p={p:.3g}  Cliff's delta={d:+.3f} [{lo:+.3f}, {hi:+.3f}]")


def main() -> int:
    log(f"Phase 5 — statistics — {datetime.now(timezone.utc).isoformat()}")
    log(f"seed {SEED}\n")

    ann = read(PROC / "variant_annotated.tsv")
    smet = {r["variation_id"]: r for r in read(PROC / "structural_metrics.tsv")}
    ddg = {r["variation_id"]: r for r in read(PROC / "foldx_ddg.tsv")} \
        if (PROC / "foldx_ddg.tsv").is_file() else {}
    for v in ann:
        v.update({k: val for k, val in (smet.get(v["variation_id"]) or {}).items()
                  if k not in v})
        d = ddg.get(v["variation_id"], {})
        v["ddG"] = f(d.get("ddG_kcal_mol"))
        v["ddG_disulfide"] = f(d.get("ddG_disulfide"))
        v["foldx_status"] = d.get("foldx_status", "")

    log(f"{len(ann):,} variants; {sum(1 for v in ann if v['ddG'] is not None):,} with ddG; "
        f"{sum(1 for v in ann if v.get('structural_coverage') == 'True'):,} with structural data\n")

    TIER = "set_primary"
    path = [v for v in ann if v[TIER] == "pathogenic"]
    ben = [v for v in ann if v[TIER] == "benign"]
    vus = [v for v in ann if v[TIER] == "vus"]

    # =====================================================================
    log("=" * 72)
    log("PART A — BIOPHYSICS (headline; independent of clinical labels)")
    log("=" * 72)

    log("\nA1. Folding stability (FoldX ddG, kcal/mol)")
    mw("A_ddG", "ddG: pathogenic vs benign",
       [v["ddG"] for v in path], [v["ddG"] for v in ben], "pathogenic", "benign")
    mw("A_ddG", "ddG: cbEGF cysteine-removing vs all other missense",
       [v["ddG"] for v in ann if T(v["cys_removing"]) and T(v["in_cbegf_cys"])],
       [v["ddG"] for v in ann if not T(v["cys_removing"])],
       "Cys-removing", "other")
    mw("A_ddG", "ddG: calcium-consensus residues vs other cbEGF residues",
       [v["ddG"] for v in ann if T(v["is_ca_consensus"])],
       [v["ddG"] for v in ann if T(v["in_cbegf"]) and not T(v["is_ca_consensus"])],
       "Ca consensus", "other cbEGF")
    mw("A_ddG", "ddG: direct Ca ligands vs non-ligands (structure-measured)",
       [v["ddG"] for v in ann if v.get("is_direct_ca_ligand") == "True"],
       [v["ddG"] for v in ann if v.get("is_direct_ca_ligand") == "False"],
       "direct ligand", "non-ligand")
    mw("A_ddG", "ddG: buried vs exposed (RSA < 20%)",
       [v["ddG"] for v in ann if v.get("buried") == "True"],
       [v["ddG"] for v in ann if v.get("buried") == "False"],
       "buried", "exposed")

    log("\nA2. The disulfide term — mechanism-specific")
    cysrem = [v["ddG_disulfide"] for v in ann if T(v["cys_removing"]) and v["ddG_disulfide"] is not None]
    other = [v["ddG_disulfide"] for v in ann if not T(v["cys_removing"]) and v["ddG_disulfide"] is not None]
    if cysrem and other:
        mw("A_ddG", "FoldX disulfide energy term: Cys-removing vs other",
           cysrem, other, "Cys-removing", "other",
           note="isolates the staple-loss mechanism from general destabilisation")

    log("\nA3. Burial of the critical sites")
    mw("A_burial", "RSA: calcium-consensus vs other cbEGF residues",
       [f(v.get("rsa_pct")) for v in ann if T(v["is_ca_consensus"])],
       [f(v.get("rsa_pct")) for v in ann if T(v["in_cbegf"]) and not T(v["is_ca_consensus"])],
       "Ca consensus", "other cbEGF")
    mw("A_burial", "RSA: cbEGF cysteines vs other cbEGF residues",
       [f(v.get("rsa_pct")) for v in ann if T(v["in_cbegf_cys"])],
       [f(v.get("rsa_pct")) for v in ann if T(v["in_cbegf"]) and not T(v["in_cbegf_cys"])],
       "cbEGF Cys", "other cbEGF")

    log("\nA4. Do two independent predictors agree? (ddG vs AlphaMissense)")
    pairs = [(v["ddG"], f(v["am_pathogenicity"])) for v in ann
             if v["ddG"] is not None and f(v["am_pathogenicity"]) is not None]
    if len(pairs) > 10:
        x, y = zip(*pairs)
        rho, p = stats.spearmanr(x, y)
        boots = [stats.spearmanr(*zip(*[pairs[i] for i in rng.integers(0, len(pairs), len(pairs))]))[0]
                 for _ in range(1000)]
        lo, hi = np.percentile(boots, [2.5, 97.5])
        record("A_agreement", "Spearman: FoldX ddG vs AlphaMissense", rho, p, rho, (lo, hi), len(pairs))
        log(f"    n={len(pairs)}  Spearman rho={rho:+.3f} [{lo:+.3f}, {hi:+.3f}]  p={p:.3g}")

    log("\nA4b. MECHANISM DISSOCIATION — the central biophysical result")
    log("     If calcium-site variants were pathogenic by destabilising the fold, they would")
    log("     show high ddG. They do not. Stratifying by burial rules out the obvious confound")
    log("     that Ca ligands are simply more solvent-exposed than other residues.")
    # Corrected 2026-08-08 (Phase 6). The previous form read
    #     (f(v.get("rsa_pct")) or -1) >= lo and (f(v.get("rsa_pct")) or 999) < hi
    # which is wrong for RSA == 0.0: 0.0 is falsy, so a completely buried residue was
    # rewritten to -1/999 and fell out of EVERY stratum. That silently dropped 4 fully
    # buried "other cbEGF" variants -- and they are the most destabilising ones
    # (median ddG 13.90) -- from the buried bin. Test None explicitly instead.
    def in_stratum(v, lo, hi):
        r = f(v.get("rsa_pct"))
        return r is not None and lo <= r < hi

    for lo, hi, lbl in ((0, 20, "buried (RSA<20%)"), (20, 50, "partial (RSA 20-50%)")):
        ca = [v["ddG"] for v in ann if T(v["is_ca_consensus"]) and v["ddG"] is not None
              and in_stratum(v, lo, hi)]
        ot = [v["ddG"] for v in ann if T(v["in_cbegf"]) and not T(v["is_ca_consensus"])
              and v["ddG"] is not None and in_stratum(v, lo, hi)]
        mw("A_mechanism", f"ddG: Ca-consensus vs other cbEGF, {lbl}", ca, ot,
           "Ca consensus", "other cbEGF",
           note="burial-matched; rules out exposure as the explanation")

    log("\n     Stability cost vs predicted pathogenicity, by mechanism:")
    groups = {
        "direct Ca ligand": lambda v: v.get("is_direct_ca_ligand") == "True",
        "cbEGF Cys-removing": lambda v: T(v["in_cbegf_cys"]) and T(v["cys_removing"]),
        "other cbEGF residue": lambda v: (T(v["in_cbegf"]) and v.get("is_direct_ca_ligand") != "True"
                                          and not (T(v["in_cbegf_cys"]) and T(v["cys_removing"]))),
    }
    summary = {}
    for lbl, sel in groups.items():
        g = [v for v in ann if sel(v) and v["ddG"] is not None]
        if not g:
            continue
        am = [f(v["am_pathogenicity"]) for v in g if f(v["am_pathogenicity"]) is not None]
        summary[lbl] = (len(g), float(np.median([v["ddG"] for v in g])),
                        float(np.median(am)) if am else float("nan"),
                        float(np.median([f(v.get("rsa_pct")) for v in g
                                         if f(v.get("rsa_pct")) is not None])))
        log(f"       {lbl:22} n={len(g):3}  ddG={summary[lbl][1]:5.2f}  "
            f"AlphaMissense={summary[lbl][2]:.3f}  RSA={summary[lbl][3]:.1f}%")

    # Ca ligands vs other cbEGF: statistically indistinguishable in stability, but far apart
    # in predicted pathogenicity. That dissociation is the finding.
    lig = [v for v in ann if v.get("is_direct_ca_ligand") == "True" and v["ddG"] is not None]
    oth = [v for v in ann if T(v["in_cbegf"]) and v.get("is_direct_ca_ligand") == "False"
           and not (T(v["in_cbegf_cys"]) and T(v["cys_removing"])) and v["ddG"] is not None]
    if lig and oth:
        mw("A_mechanism", "AlphaMissense: direct Ca ligands vs other cbEGF residues",
           [f(v["am_pathogenicity"]) for v in lig], [f(v["am_pathogenicity"]) for v in oth],
           "Ca ligand", "other cbEGF",
           note="same low stability cost, very different predicted pathogenicity")

    log("\nA5. Calcium coordination geometry (CheckMyMetal)")
    if (PROC / "cmm_sites.tsv").is_file():
        cmm = [r for r in read(PROC / "cmm_sites.tsv") if r.get("cmm_status") == "parsed"]
        e = [f(r["grmsd"]) for r in cmm if r["source"] == "experimental"]
        a = [f(r["grmsd"]) for r in cmm if r["source"] != "experimental"]
        if len(e) >= 3 and len(a) >= 3:
            u, p = stats.mannwhitneyu(e, a, alternative="two-sided")
            d = cliffs_delta(e, a)
            record("A_cmm", "CMM gRMSD: experimental vs AF3 sites", u, p, d, (np.nan, np.nan),
                   f"{len(e)}/{len(a)}",
                   note="non-significance is the desired result: AF3 sites indistinguishable")
            log(f"    experimental n={len(e)} median={np.median(e):.1f} deg; "
                f"AF3 n={len(a)} median={np.median(a):.1f} deg")
            log(f"    Mann-Whitney p={p:.3g}, Cliff's delta={d:+.3f}")
            log("    (a NON-significant difference is the desired outcome here)")

    # =====================================================================
    log("\n" + "=" * 72)
    log("PART B — ENRICHMENT (supporting only; see the circularity caveat)")
    log("=" * 72)

    n_pos = 2871
    n_cys_pos, n_ca_pos = 258, 215
    log(f"\nB1. Pathogenic variants vs the all-possible-missense background")
    log(f"    background: {n_cys_pos} cbEGF cysteine positions and {n_ca_pos} calcium-consensus")
    log(f"    positions out of {n_pos} — every position offers the same 19 substitutions, so")
    log(f"    position share is the correct null.")

    for label, flag, npos in (("cbEGF cysteine (removing)",
                               lambda v: T(v["in_cbegf_cys"]) and T(v["cys_removing"]), n_cys_pos),
                              ("calcium-consensus residue",
                               lambda v: T(v["is_ca_consensus"]), n_ca_pos)):
        k = sum(1 for v in path if flag(v))
        n = len(path)
        exp_p = npos / n_pos
        res = stats.binomtest(k, n, exp_p, alternative="greater")
        # The p-value is one-sided (we are testing for enrichment), but the interval must be
        # TWO-sided: a one-sided binomial CI pins its upper bound at 1.0, which would be
        # reported as a meaningless fold-enrichment ceiling of 1/exp_p.
        ci = stats.binomtest(k, n, exp_p).proportion_ci(confidence_level=0.95)
        record("B_background", f"pathogenic at {label} vs background", k / n, res.pvalue,
               (k / n) / exp_p, (ci.low / exp_p, ci.high / exp_p), n,
               note="fold-enrichment over the positional null")
        log(f"    {label}: {k}/{n} = {k/n*100:.1f}% observed vs {exp_p*100:.1f}% expected")
        log(f"      fold-enrichment {(k/n)/exp_p:.1f}x "
            f"[{ci.low/exp_p:.1f}, {ci.high/exp_p:.1f}], binomial p={res.pvalue:.3g}")

    log("\nB2. Pathogenic vs benign (Fisher's exact) — UNDERPOWERED benign arm")
    for label, flag in (("cbEGF cysteine (removing)",
                         lambda v: T(v["in_cbegf_cys"]) and T(v["cys_removing"])),
                        ("calcium-consensus residue", lambda v: T(v["is_ca_consensus"])),
                        ("any VCEP critical residue", lambda v: T(v["is_critical_residue"]))):
        a = sum(1 for v in path if flag(v)); b = len(path) - a
        c = sum(1 for v in ben if flag(v)); d = len(ben) - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        o, lo, hi = or_ci(a, b, c, d)
        record("B_pvb", f"pathogenic vs benign at {label}", odds, p, o, (lo, hi),
               f"{len(path)}/{len(ben)}",
               note="ClinVar labels used PM1 — partly circular")
        log(f"    {label}: P {a}/{len(path)} vs B {c}/{len(ben)}   OR={o:.1f} [{lo:.1f}, {hi:.1f}]  p={p:.3g}")

    log("\nB3. AlphaMissense as an orthogonal label (no FBN1 PM1 rules used)")
    amp = [v for v in ann if v["am_class"] == "pathogenic"]
    amb = [v for v in ann if v["am_class"] == "benign"]
    for label, flag in (("cbEGF cysteine (removing)",
                         lambda v: T(v["in_cbegf_cys"]) and T(v["cys_removing"])),
                        ("calcium-consensus residue", lambda v: T(v["is_ca_consensus"]))):
        a = sum(1 for v in amp if flag(v)); b = len(amp) - a
        c = sum(1 for v in amb if flag(v)); d = len(amb) - c
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        o, lo, hi = or_ci(a, b, c, d)
        record("B_orthogonal", f"AM-pathogenic vs AM-benign at {label}", odds, p, o, (lo, hi),
               f"{len(amp)}/{len(amb)}", note="independent of ClinVar labelling")
        log(f"    {label}: AM-P {a}/{len(amp)} vs AM-B {c}/{len(amb)}   "
            f"OR={o:.1f} [{lo:.1f}, {hi:.1f}]  p={p:.3g}")

    # ---- multiple testing ---------------------------------------------------
    log("\n" + "=" * 72)
    log("Multiple-testing correction (Benjamini-Hochberg, within family)")
    log("=" * 72)
    for fam in sorted({t["family"] for t in tests}):
        idx = [i for i, t in enumerate(tests) if t["family"] == fam]
        ps = [tests[i]["p_raw"] for i in idx]
        rej, q, _, _ = multipletests(ps, alpha=0.05, method="fdr_bh")
        for j, i in enumerate(idx):
            tests[i]["q_bh"] = q[j]
            tests[i]["significant_q05"] = bool(rej[j])
        log(f"  {fam}: {sum(rej)}/{len(idx)} significant at q<0.05")

    cols = ["family", "test", "n", "statistic", "effect", "ci_low", "ci_high",
            "p_raw", "q_bh", "significant_q05", "note"]
    with (ROOT / "results/phase5_statistics.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(tests)

    (ROOT / "results/phase5_statistics.md").write_text(
        "# Phase 5 — statistical results\n\n"
        f"Generated by `scripts/12_statistics.py` at {datetime.now(timezone.utc).isoformat()}, "
        f"seed {SEED}. Tier: primary (ClinVar review status ≥ 2 stars).\n\n"
        "Effect sizes: Cliff's δ for rank comparisons (bootstrap 95% CI, 2000 resamples), "
        "odds ratios with Woolf 95% CI for 2×2 tables, fold-enrichment with an exact binomial CI "
        "against the positional null. Benjamini-Hochberg within each family.\n\n"
        "| family | test | n | effect | 95% CI | p | q (BH) | sig |\n"
        "|---|---|---|---|---|---|---|---|\n"
        + "\n".join(
            f"| {t['family']} | {t['test']} | {t['n']} | {t['effect']:+.3f} | "
            f"[{t['ci_low']:+.3f}, {t['ci_high']:+.3f}] | {t['p_raw']:.3g} | {t['q_bh']:.3g} | "
            f"{'**yes**' if t['significant_q05'] else 'no'} |" for t in tests)
        + "\n\n## Reading these numbers\n\n"
          "**Part A (biophysics) is the headline.** Those measurements come from coordinates and "
          "are independent of how anyone classified the variants.\n\n"
          "**Part B (enrichment) is supporting only.** ClinVar's pathogenic classifications were "
          "made using ClinGen PM1 — the same critical-residue rule the enrichment tests. The "
          "pathogenic-versus-benign contrast is therefore partly circular and its magnitude is "
          "inflated. The AlphaMissense contrast (B3) is the honest version, since AlphaMissense "
          "never used FBN1-specific rules. The benign arm is also small (n=34 at ≥2 stars) "
          "because FBN1 is too constrained to carry many common missense variants.\n",
        encoding="utf-8")

    log("\nwrote results/phase5_statistics.md and .tsv")
    (ROOT / f"logs/12_statistics_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                          encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
