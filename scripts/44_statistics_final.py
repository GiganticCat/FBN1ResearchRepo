#!/usr/bin/env python3
"""
44_statistics_final.py — the single statistics table the manuscript cites.

`28_statistics_v2.py` ran on 2026-08-09 and covers AlphaMissense, FoldX and geometry. Two things
happened after it: PyRosetta ddG was computed for 1,934 variants and then corrected twice (the
metal-constraint detonation, then the fixed-backbone steric artefact), and batch 4 turned the
mutant-cofolding demonstration into a rate. Neither is in that table, so the manuscript has been
quoting Rosetta numbers out of a figure sidecar and the rate out of a second report. This
consolidates all of it into one place, recomputed from `data/processed/` in a single pass.

It does not replace `28`: that script's stratified (wild-type-matched) machinery is reused here
by import rather than re-typed, so the two cannot drift.

Families, each corrected within itself by Benjamini-Hochberg:
  A  sequence-based prediction        AlphaMissense across the four mechanism groups
  B  folding stability                FoldX (kcal/mol) and Rosetta (REU), never pooled
  C  position                         solvent accessibility, distance to the ion
  D  calcium-site geometry            batch-4 disruption, magnitude and donor chemistry
  E  enrichment                       pathogenic vs population, with the circularity noted

Writes:
  results/statistics_final.tsv
  results/statistics_final.md
  logs/44_statistics_final_<stamp>.log
"""
from __future__ import annotations

import datetime as dt
import importlib.util
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import figstyle as F  # noqa: E402

SEED = 20260829
OUT_TSV = ROOT / "results/statistics_final.tsv"
OUT_MD = ROOT / "results/statistics_final.md"


def _v2():
    spec = importlib.util.spec_from_file_location(
        "_stats_v2", ROOT / "scripts" / "28_statistics_v2.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


V2 = _v2()
cliffs_delta, boot_ci, stratified_delta = V2.cliffs_delta, V2.boot_ci, V2.stratified_delta
mw_p, bh = V2.mw_p, V2.bh

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"44_statistics_final_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def main() -> int:
    random.seed(SEED)
    np.random.seed(SEED)
    rng = random.Random(SEED)
    log(f"Final statistics — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  seed {SEED}; {V2.N_BOOT} bootstrap resamples; BH correction within each family")
    log()

    df = F.load_v2()
    df["g"] = F.mech_group_v2(df)
    rows = []

    def cmp_test(family, label, A, B, wtA=None, wtB=None, note=""):
        a = [v for v in A if pd.notna(v)]
        b = [v for v in B if pd.notna(v)]
        d = cliffs_delta(a, b)
        lo, hi = boot_ci(a, b, rng)
        sd, ns, _ = (stratified_delta(wtA, wtB) if wtA is not None else (float("nan"), 0, 0))
        rows.append(dict(family=family, test=label, n_a=len(a), n_b=len(b),
                         median_a=round(float(np.median(a)), 4) if a else np.nan,
                         median_b=round(float(np.median(b)), 4) if b else np.nan,
                         effect="Cliff's delta", estimate=round(d, 3),
                         ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                         estimate_wt_matched=round(sd, 3) if ns else np.nan,
                         n_strata=ns, p_raw=mw_p(a, b), note=note))

    def col(g, c):
        s = df[df.g == g]
        return s[c].to_numpy(), list(zip(s.wt_aa, s[c]))

    # ---------------- A  sequence-based prediction
    for g, lab in (("ca_ligand", "Ca ligands"), ("motif_overcalled", "motif non-ligands"),
                   ("cys_removing", "cysteine-removing")):
        a, wa = col(g, "am")
        b, wb = col("composition_control", "am")
        cmp_test("A_alphamissense", f"AlphaMissense: {lab} vs composition control",
                 a, b, wa, wb, "composition control = other D/E/N positions with coverage")

    # ---------------- B  stability, the two force fields kept apart
    for metric, unit in (("foldx_ddG", "kcal/mol"), ("rosetta_ddG", "REU")):
        for g, lab in (("ca_ligand", "Ca ligands"), ("motif_overcalled", "motif non-ligands"),
                       ("cys_removing", "cysteine-removing")):
            a, wa = col(g, metric)
            b, wb = col("composition_control", metric)
            cmp_test("B_stability", f"{metric} ({unit}): {lab} vs composition control",
                     a, b, wa, wb,
                     "FoldX covers only experimental templates" if metric == "foldx_ddG"
                     else "metal-aware; the better-powered of the two")

    # ---------------- C  position
    a, wa = col("ca_ligand", "rsa")
    b, wb = col("composition_control", "rsa")
    cmp_test("C_position", "relative solvent accessibility: Ca ligands vs composition control",
             a, b, wa, wb, "ligands are more buried; Figure 3 tests whether that explains B")
    for band, lo, hi in (("buried (RSA<20%)", -1, 20), ("partial (20-50%)", 20, 50)):
        s = df[(df.rsa >= lo) & (df.rsa < hi)]
        A = s[s.g == "ca_ligand"].rosetta_ddG
        B = s[s.g == "composition_control"].rosetta_ddG
        cmp_test("C_position", f"rosetta_ddG within {band}: Ca ligands vs composition control",
                 A.to_numpy(), B.to_numpy(), note="burial-stratified control")

    # ---------------- D  calcium-site geometry (batch 4)
    r = pd.read_csv(F.PROC / "batch4_rate.tsv", sep="\t")
    T = float(r.null_threshold_pct.iloc[0])
    lig = r[r.klass == "ca_ligand_pathogenic"]
    ctrl = r[r.klass != "ca_ligand_pathogenic"]
    cmp_test("D_geometry", "delta BVS %: pathogenic Ca ligands vs all controls",
             lig.delta_bvs_pct.to_numpy(), ctrl.delta_bvs_pct.to_numpy(),
             note=f"cognate site; null threshold {T:.1f}% from 123 no-effect sites")
    cmp_test("D_geometry", "local backbone RMSD: pathogenic Ca ligands vs all controls",
             lig.local_rmsd_A.to_numpy(), ctrl.local_rmsd_A.to_numpy(),
             note="the fold, measured on the same models as the row above")

    def rate_row(label, k, n, k2, n2, note=""):
        orr, p = stats.fisher_exact([[k, n - k], [k2, n2 - k2]])
        lo, hi = wilson(k, n)
        rows.append(dict(family="D_geometry", test=label, n_a=n, n_b=n2,
                         median_a=round(100 * k / n, 1), median_b=round(100 * k2 / n2, 1),
                         effect="odds ratio", estimate=round(orr, 3) if np.isfinite(orr) else np.inf,
                         ci_lo=round(100 * lo, 1), ci_hi=round(100 * hi, 1),
                         estimate_wt_matched=np.nan, n_strata=0, p_raw=p,
                         note=note + " (median columns are % disrupted; CI is Wilson on group A)"))

    kk = int(lig.disrupted.sum())
    rate_row("disruption rate: pathogenic Ca ligands vs all controls", kk, len(lig),
             int(ctrl.disrupted.sum()), len(ctrl))
    dk = lig[lig.donor_kept]
    dl = lig[~lig.donor_kept]
    rate_row("disruption rate: donor removed vs donor retained", int(dl.disrupted.sum()), len(dl),
             int(dk.disrupted.sum()), len(dk),
             note="within the pathogenic ligand set")
    # the donor rule is not uniform across the site: it holds at the N-terminal Asp and Glu and
    # fails at the beta-hydroxylation Asn, which is the position the ClinGen FBN1 VCEP treats
    # differently. Pooling the two would hide exactly the thing worth reporting.
    for posc, label in (("Asp/Glu", "at the Asp/Glu positions"),
                        ("beta-OH Asn", "at the β-OH Asn position")):
        s_ = lig[lig.site_pos == posc]
        a_, b_ = s_[~s_.donor_kept], s_[s_.donor_kept]
        if len(a_) and len(b_):
            rate_row(f"disruption rate: donor removed vs retained, {label}",
                     int(a_.disrupted.sum()), len(a_), int(b_.disrupted.sum()), len(b_),
                     note="stratified by consensus position")
    ns = lig[(lig.wt_aa == "N") & (lig.mut_aa == "S") & (lig.site_pos == "beta-OH Asn")]
    rest = lig[~lig.index.isin(ns.index)]
    if len(ns):
        rate_row("disruption rate: Asn→Ser at the β-OH Asn vs every other ligand substitution",
                 int(ns.disrupted.sum()), len(ns), int(rest.disrupted.sum()), len(rest),
                 note="the substitution the ClinGen FBN1 VCEP exempts from PM1")

    # geometry against energetics on the same variants
    m = lig.merge(pd.read_csv(F.PROC / "rosetta_ddg.tsv", sep="\t")
                  .query("status == 'ok'")[["variation_id", "ddg_mean"]],
                  on="variation_id", how="left").dropna(subset=["ddg_mean"])
    rho, p = stats.spearmanr(m.delta_bvs_pct, m.ddg_mean)
    rows.append(dict(family="D_geometry",
                     test="Spearman(delta BVS %, Rosetta ddG) across pathogenic Ca ligands",
                     n_a=len(m), n_b=np.nan, median_a=np.nan, median_b=np.nan,
                     effect="Spearman rho", estimate=round(float(rho), 3),
                     ci_lo=np.nan, ci_hi=np.nan, estimate_wt_matched=np.nan, n_strata=0,
                     p_raw=float(p),
                     note="the two axes of damage are independent within this set"))
    am = lig.merge(df[["variation_id", "am"]], on="variation_id", how="left")
    cmp_test("D_geometry", "AlphaMissense: disrupted vs intact Ca-ligand variants",
             am[am.disrupted].am.to_numpy(), am[~am.disrupted].am.to_numpy(),
             note="a sequence method cannot separate the variants that break the site")

    # ---------------- E  enrichment (partly circular, reported for completeness)
    enr = pd.read_csv(F.PROC / "enrichment_v2.tsv", sep="\t")
    for _, e in enr.iterrows():
        rows.append(dict(family="E_enrichment",
                         test=f"pathogenic vs {e.denominator} at {e.feature}",
                         n_a=int(e.path_yes + e.path_no), n_b=int(e.denom_yes + e.denom_no),
                         median_a=float(e.pct_path), median_b=float(e.pct_denom),
                         effect="odds ratio", estimate=float(e.odds_ratio_corrected),
                         ci_lo=float(e.ci_lo), ci_hi=float(e.ci_hi),
                         estimate_wt_matched=np.nan, n_strata=0, p_raw=float(e.p_raw),
                         note="PARTLY CIRCULAR: ClinGen rules assign pathogenicity partly on "
                              "position at these very features"))

    out = pd.DataFrame(rows)
    out["p_adj"] = np.nan
    for fam in out.family.unique():
        m_ = out.family == fam
        out.loc[m_, "p_adj"] = bh(list(out.loc[m_, "p_raw"]))
    out.to_csv(OUT_TSV, sep="\t", index=False)

    for fam in out.family.unique():
        log(f"--- {fam}")
        for _, x in out[out.family == fam].iterrows():
            ci = (f"[{x.ci_lo}, {x.ci_hi}]" if pd.notna(x.ci_lo) else "")
            log(f"  {x.test[:78]:78s} n={x.n_a}/{x.n_b}  {x.effect} {x.estimate} {ci}  "
                f"q={x.p_adj:.3g}")
        log()

    with open(OUT_MD, "w") as fh:
        fh.write("# Final statistics\n\n")
        fh.write(f"Generated by `scripts/44_statistics_final.py`, "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}. Seed {SEED}, "
                 f"{V2.N_BOOT} bootstrap resamples, Benjamini-Hochberg within each family.\n\n")
        fh.write("Supersedes `phase5_statistics_v2.md`, which predates both the corrected "
                 "PyRosetta ΔΔG and the batch-4 geometry. Where a wild-type-matched estimate is "
                 "given it is the one to believe: calcium ligands are D/E/N by definition, and "
                 "that alone shifts any unmatched comparison.\n\n")
        fh.write("| family | test | n A/B | median A / B | effect | estimate | 95% CI | "
                 "WT-matched | q |\n|---|---|---|---|---|---|---|---|---|\n")
        for _, x in out.iterrows():
            ci = f"[{x.ci_lo}, {x.ci_hi}]" if pd.notna(x.ci_lo) else "—"
            wm = f"{x.estimate_wt_matched} ({int(x.n_strata)})" if pd.notna(x.estimate_wt_matched) else "—"
            med = (f"{x.median_a} / {x.median_b}" if pd.notna(x.median_a) else "—")
            fh.write(f"| {x.family} | {x.test} | {x.n_a}/{x.n_b} | {med} | {x.effect} | "
                     f"{x.estimate} | {ci} | {wm} | {x.p_adj:.3g} |\n")
        fh.write("\n**Reading the families.** B never pools the two force fields: FoldX is in "
                 "kcal/mol over 751 variants on experimental templates, Rosetta in Rosetta "
                 "energy units over 1,934 variants with modelled coverage and explicit metal "
                 "bonds. D is the only family that measures the calcium site itself. E is "
                 "reported for completeness and is partly circular.\n")
    log(f"wrote {OUT_TSV.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)} ({len(out)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
