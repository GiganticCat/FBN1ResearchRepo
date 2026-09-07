#!/usr/bin/env python3
"""
41_batch4_rate.py — PHASE 5. How OFTEN does a pathogenic calcium-ligand variant break the site?

`33_v2fig2_geometry.py` could show that calcium-site disruption happens and can be measured
against seed noise, on nine variants. It could not say how often, and "how often" is the
difference between a demonstration and a result. Batch 4 folded 34 pathogenic side-chain ligand
variants, 12 composition controls (D/E/N substitutions that coordinate nothing) and 14 negative
controls, each five times, each against its own five-model wild type. This script turns those
into a rate with a measured null.

THE NULL IS MEASURED, NOT ASSUMED. Every construct carries three calcium sites and a variant can
reach only one of them. The other two are the same fold, the same five seeds, and a site the
substitution cannot touch -- 120 measurements of "nothing happened here", from the same models
that supply the signal. Pooled with the negative controls' own cognate sites, they define the
threshold: a change counts as disruption only if it exceeds the 95th percentile of the absolute
change seen where no effect is possible. Nothing about the threshold is chosen by hand.

TWO AXES, KEPT SEPARATE. The disruption call uses bond-valence sum, which is a property of the
calcium site. The fold is checked separately with the local backbone RMSD around that site. The
paper's claim is that these come apart -- the site loses grip while the fold does not move -- so
they must never be combined into one score.

DONOR CHEMISTRY IS TESTED, NOT ASSUMED. Half of these substitutions replace one oxygen donor
with another (D->N, D->E, E->Q, N->D). If the ion cares about the atom rather than the residue
label, those should behave differently from substitutions that remove the donor entirely, and
that split is a prediction the data can refute.

AND IT IS REFUTED IN ONE PLACE, WHICH IS THE INTERESTING PART. The consensus site has two kinds
of position: the aspartate and glutamate at the domain's N terminus (offsets 0 and +3), and the
asparagine at offset +17/+18 that is the beta-hydroxylation site. The donor rule holds cleanly at
the first and fails at the second, where replacing the asparagine carboxamide with a serine
hydroxyl still costs coordination. That distinction is not cosmetic: the ClinGen FBN1 expert
panel declines to apply PM1 to exactly this substitution at exactly this position, so the two
strata are reported separately rather than pooled.

Writes:
  data/processed/batch4_rate.tsv         one row per variant (cognate site)
  results/phase5_batch4_rate.md
  logs/41_batch4_rate_<stamp>.log
"""
from __future__ import annotations

import datetime as dt
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
OUT_TSV = PROC / "batch4_rate.tsv"
OUT_MD = ROOT / "results/phase5_batch4_rate.md"

SEED = 20260829
N_BOOT = 2000
NULL_PCTILE = 95.0

# Substitutions that keep a carboxylate/amide oxygen at the position. The ion is coordinated by
# an atom, not by a residue name, so this is the chemically meaningful split.
O_DONOR = set("DENQST")

STAMP = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
LOGPATH = ROOT / "logs" / f"41_batch4_rate_{STAMP}.log"
_fh = None


def log(msg: str = "") -> None:
    global _fh
    if _fh is None:
        LOGPATH.parent.mkdir(parents=True, exist_ok=True)
        _fh = open(LOGPATH, "w")
    print(msg)
    _fh.write(msg + "\n")
    _fh.flush()


def cliffs_delta(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    if not a.size or not b.size:
        return float("nan")
    gt = sum((a[:, None] > b[None, :]).sum(1))
    lt = sum((a[:, None] < b[None, :]).sum(1))
    return float((gt - lt) / (a.size * b.size))


def boot_ci(a, b, rng) -> tuple[float, float]:
    a, b = list(a), list(b)
    if len(a) < 3 or len(b) < 3:
        return float("nan"), float("nan")
    ds = []
    for _ in range(N_BOOT):
        ds.append(cliffs_delta([rng.choice(a) for _ in a], [rng.choice(b) for _ in b]))
    return float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5))


def wilson(k: int, n: int) -> tuple[float, float]:
    """Wilson score interval — behaves at 0/n and n/n, where the normal approximation does not."""
    if n == 0:
        return float("nan"), float("nan")
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return float((c - h) / d), float((c + h) / d)


def main() -> int:
    rng = random.Random(SEED)
    np.random.seed(SEED)
    log(f"Batch-4 disruption rate — {dt.datetime.now(dt.timezone.utc).isoformat()}")
    log(f"  seed {SEED}; {N_BOOT} bootstrap resamples; null threshold = "
        f"{NULL_PCTILE}th percentile of |delta BVS %| where no effect is possible")
    log()

    g = pd.read_csv(PROC / "af3_batch4_geometry.tsv", sep="\t")
    log(f"{len(g)} site rows, {g.job.nunique()} jobs")

    # ---- admission: wild-type quality filter, fixed in script 40 on wild-type evidence alone
    bad = g[~g.wt_site_ok]
    log(f"{bad.site.count()} site rows sit on wild-type sites that failed the quality filter "
        f"({bad.construct.nunique()} construct(s)); of these "
        f"{int(bad.is_cognate.sum())} are cognate sites and their variants drop out:")
    for _, r in bad[bad.is_cognate].iterrows():
        log(f"    {r.job} ({r.klass}) — {r.construct} site {r.site}, "
            f"wild-type seed spread {r.wt_seed_sd_pct}% of valence")
    g = g[g.wt_site_ok].copy()

    vacated = g[g.status == "vacated"]
    if len(vacated):
        log(f"{len(vacated)} site row(s) vacated in all five models — counted as disruption "
            "where cognate")
    g["delta_pct"] = pd.to_numeric(g.delta_bvs_pct, errors="coerce")

    cog = g[g.is_cognate].copy()
    byst = g[(~g.is_cognate) & (g.status == "present")].copy()

    # ---- the null
    neg_cog = cog[(cog.klass == "negative_control") & (cog.status == "present")]
    null = np.concatenate([byst.delta_pct.to_numpy(), neg_cog.delta_pct.to_numpy()])
    T = float(np.percentile(np.abs(null), NULL_PCTILE))
    log()
    log(f"NULL: {len(byst)} bystander sites + {len(neg_cog)} negative-control cognate sites "
        f"= {len(null)} measurements")
    log(f"  median {np.median(null):+.2f}%, IQR {np.percentile(null,25):+.2f} to "
        f"{np.percentile(null,75):+.2f}%, |delta| 95th percentile = {T:.2f}%")
    log(f"  bystanders alone: {np.percentile(np.abs(byst.delta_pct),NULL_PCTILE):.2f}%; "
        f"negative controls alone: "
        f"{np.percentile(np.abs(neg_cog.delta_pct),NULL_PCTILE):.2f}%")
    log(f"THRESHOLD: a cognate site is disrupted if its valence falls by more than {T:.2f}%, "
        "or if the site is vacated outright.")
    log()

    cog["disrupted"] = (cog.status == "vacated") | (cog.delta_pct <= -T)
    cog["gain"] = cog.delta_pct >= T
    cog["donor_kept"] = [m in O_DONOR for m in cog.mut_aa]

    # ---- rates
    log("DISRUPTION RATE at the cognate site")
    order = ["ca_ligand_pathogenic", "composition_control", "negative_control"]
    rates = {}
    for k in order:
        s = cog[cog.klass == k]
        n, d = len(s), int(s.disrupted.sum())
        lo, hi = wilson(d, n)
        rates[k] = (d, n)
        log(f"  {k:22s} {d:2d}/{n:2d} = {100*d/n:5.1f}%  (95% CI {100*lo:.1f}–{100*hi:.1f}%)"
            f"   median delta {s.delta_pct.median():+.1f}%")

    log()
    log("Fisher's exact, pathogenic ligands against each control group")
    fisher_rows = []
    for k in ("composition_control", "negative_control"):
        a, na = rates["ca_ligand_pathogenic"]
        b, nb = rates[k]
        orr, p = stats.fisher_exact([[a, na - a], [b, nb - b]])
        fisher_rows.append((k, a, na, b, nb, orr, p))
        log(f"  vs {k:22s} {a}/{na} vs {b}/{nb}   OR {orr:.2f}   p = {p:.4g}")
    a, na = rates["ca_ligand_pathogenic"]
    b = rates["composition_control"][0] + rates["negative_control"][0]
    nb = rates["composition_control"][1] + rates["negative_control"][1]
    orr_all, p_all = stats.fisher_exact([[a, na - a], [b, nb - b]])
    log(f"  vs all controls pooled     {a}/{na} vs {b}/{nb}   OR {orr_all:.2f}   p = {p_all:.4g}")

    # ---- magnitude
    log()
    log("MAGNITUDE of the valence change (cognate site, present in >=1 model)")
    pres = cog[cog.status == "present"]
    lig = pres[pres.klass == "ca_ligand_pathogenic"].delta_pct.to_numpy()
    comp = pres[pres.klass == "composition_control"].delta_pct.to_numpy()
    neg = pres[pres.klass == "negative_control"].delta_pct.to_numpy()
    mag_rows = []
    for label, A, B in (("pathogenic ligands vs composition controls", lig, comp),
                        ("pathogenic ligands vs negative controls", lig, neg)):
        d = cliffs_delta(A, B)
        lo, hi = boot_ci(A, B, rng)
        p = float(stats.mannwhitneyu(A, B, alternative="two-sided").pvalue)
        mag_rows.append((label, len(A), len(B), float(np.median(A)), float(np.median(B)), d, lo, hi, p))
        log(f"  {label:44s} {np.median(A):+6.1f}% vs {np.median(B):+5.1f}%  "
            f"Cliff's d {d:+.3f} [{lo:+.3f}, {hi:+.3f}]  p = {p:.4g}")

    # ---- donor chemistry
    log()
    log("DONOR CHEMISTRY within the pathogenic ligand set")
    L = cog[cog.klass == "ca_ligand_pathogenic"]
    dk, dl = L[L.donor_kept], L[~L.donor_kept]
    for label, s in (("substitution keeps an O donor", dk), ("donor removed", dl)):
        n, d = len(s), int(s.disrupted.sum())
        lo, hi = wilson(d, n)
        log(f"  {label:32s} {d:2d}/{n:2d} = {100*d/n:5.1f}%  (95% CI {100*lo:.1f}–{100*hi:.1f}%)"
            f"   median delta {s.delta_pct.median():+.1f}%")
    orr_d, p_d = stats.fisher_exact(
        [[int(dl.disrupted.sum()), len(dl) - int(dl.disrupted.sum())],
         [int(dk.disrupted.sum()), len(dk) - int(dk.disrupted.sum())]])
    log(f"  donor removed vs donor kept:  OR {orr_d:.2f}   p = {p_d:.4g}")
    # contact retention is the mechanism this predicts
    kept = pd.to_numeric(L.still_contacts, errors="coerce")
    log(f"  the mutated side chain still contacts the ion in all five models for "
        f"{int((kept == 5).sum())}/{len(L)} pathogenic ligand variants, in none for "
        f"{int((kept == 0).sum())}")

    # ---- donor chemistry, stratified by which consensus position is hit
    log()
    log("BY CONSENSUS POSITION within the pathogenic ligand set")
    dom = pd.read_csv(ROOT / "data/interim/uniprot_domains.tsv", sep="\t")
    cb = dom[dom["kind"] == "cbEGF"]

    def offset(pos):
        d = cb[(cb["start"] <= pos) & (cb["end"] >= pos)]
        return int(pos - d.iloc[0]["start"]) if len(d) == 1 else None

    L = L.copy()
    L["offset"] = L["uniprot_pos"].map(offset)
    # offsets 0 and +3 are the consensus Asp and Glu; +16 to +20 is the beta-hydroxylation Asn
    L["site_pos"] = L["offset"].map(
        lambda o: "Asp/Glu (offset 0, +3)" if o in (0, 3)
        else ("beta-OH Asn (offset +16..+20)" if o is not None and 16 <= o <= 20 else "other"))
    for pos_class in ("Asp/Glu (offset 0, +3)", "beta-OH Asn (offset +16..+20)"):
        s_ = L[L["site_pos"] == pos_class]
        for keep, lab in ((True, "donor kept"), (False, "donor removed")):
            t = s_[s_["donor_kept"] == keep]
            if not len(t):
                continue
            log(f"  {pos_class:30s} {lab:14s} {int(t.disrupted.sum())}/{len(t)}  "
                f"median {t.delta_pct.median():+6.1f}%  "
                f"[{', '.join(f'{r.wt_aa}{r.uniprot_pos}{r.mut_aa}' for r in t.itertuples())}]")
    ns = L[(L["wt_aa"] == "N") & (L["mut_aa"] == "S")
           & (L["site_pos"] == "beta-OH Asn (offset +16..+20)")]
    if len(ns):
        log(f"  Asn->Ser at the beta-hydroxylation position: {int(ns.disrupted.sum())}/{len(ns)} "
            f"disrupted, median {ns.delta_pct.median():+.1f}% "
            f"({', '.join(f'{r.wt_aa}{r.uniprot_pos}{r.mut_aa}' for r in ns.itertuples())}) "
            "— the substitution the ClinGen FBN1 VCEP exempts from PM1")
    cog["site_pos"] = cog["uniprot_pos"].map(
        lambda p_: ("Asp/Glu" if offset(p_) in (0, 3)
                    else ("beta-OH Asn" if offset(p_) is not None and 16 <= offset(p_) <= 20
                          else "other")))

    # ---- the fold, measured on the same models
    log()
    log("THE FOLD, measured locally around the same site")
    for k in order:
        s = cog[(cog.klass == k) & (cog.status == "present")]
        lr = pd.to_numeric(s.local_rmsd_A, errors="coerce")
        gr = pd.to_numeric(s.rmsd_ca_global_A, errors="coerce")
        log(f"  {k:22s} local RMSD median {lr.median():.2f} A (max {lr.max():.2f}); "
            f"whole-construct RMSD median {gr.median():.2f} A (max {gr.max():.2f})")
    lr_l = pd.to_numeric(cog[(cog.klass == 'ca_ligand_pathogenic') & (cog.status == 'present')].local_rmsd_A, errors='coerce').to_numpy()
    lr_c = pd.to_numeric(cog[(cog.klass != 'ca_ligand_pathogenic') & (cog.status == 'present')].local_rmsd_A, errors='coerce').to_numpy()
    d_fold = cliffs_delta(lr_l, lr_c)
    lo_f, hi_f = boot_ci(lr_l, lr_c, rng)
    p_fold = float(stats.mannwhitneyu(lr_l, lr_c, alternative="two-sided").pvalue)
    log(f"  pathogenic ligands vs all controls: Cliff's d {d_fold:+.3f} "
        f"[{lo_f:+.3f}, {hi_f:+.3f}], p = {p_fold:.4g}")

    # ---- does the geometric damage track the energetic damage?
    log()
    log("GEOMETRY AGAINST ENERGETICS on the same variants")
    ros = pd.read_csv(PROC / "rosetta_ddg.tsv", sep="\t")[["variation_id", "ddg_mean", "status"]]
    sm = pd.read_csv(PROC / "structural_metrics_af3.tsv", sep="\t")[
        ["variation_id", "am_pathogenicity", "rsa_pct"]]
    m = cog.merge(ros[ros.status == "ok"], on="variation_id", how="left").merge(sm, on="variation_id", how="left")
    sub = m[(m.klass == "ca_ligand_pathogenic") & m.ddg_mean.notna() & m.delta_pct.notna()]
    if len(sub) >= 5:
        r_s, p_s = stats.spearmanr(sub.delta_pct, sub.ddg_mean)
        log(f"  Spearman(delta BVS %, Rosetta ddG) over {len(sub)} pathogenic ligands: "
            f"rho = {r_s:+.3f}, p = {p_s:.3g}")
        log(f"    disrupted    ddG median {sub[sub.disrupted].ddg_mean.median():.2f} REU "
            f"(n={int(sub.disrupted.sum())})")
        log(f"    not disrupted ddG median {sub[~sub.disrupted].ddg_mean.median():.2f} REU "
            f"(n={int((~sub.disrupted).sum())})")
        am = m[(m.klass == 'ca_ligand_pathogenic') & m.am_pathogenicity.notna()]
        log(f"    AlphaMissense median: disrupted "
            f"{am[am.disrupted].am_pathogenicity.median():.3f}, not disrupted "
            f"{am[~am.disrupted].am_pathogenicity.median():.3f} — a sequence method cannot "
            "tell these apart")

    # ---- write per-variant table
    keep = ["job", "variant", "variation_id", "klass", "clinical", "construct", "uniprot_pos",
            "site_pos",
            "wt_aa", "mut_aa", "site", "cognate_by", "status", "bvs_wt", "bvs_mut",
            "delta_bvs_pct", "z_vs_noise", "wt_seed_sd_pct", "cn_wt", "cn_mut", "ion_shift_A",
            "ion_shift_wt_sd_A", "local_rmsd_A", "rmsd_ca_global_A", "plddt_site_wt",
            "plddt_site_mut", "still_contacts", "donor_kept", "disrupted", "gain"]
    out = cog[keep].copy()
    out.insert(len(out.columns), "null_threshold_pct", round(T, 2))
    out = out.sort_values(["klass", "delta_bvs_pct"])
    out.to_csv(OUT_TSV, sep="\t", index=False)
    log()
    log(f"wrote {OUT_TSV.relative_to(ROOT)} ({len(out)} variants)")

    # ---- report
    with open(OUT_MD, "w") as fh:
        fh.write("# How often does a pathogenic calcium-ligand variant break the site?\n\n")
        fh.write(f"Generated by `scripts/41_batch4_rate.py`, "
                 f"{dt.datetime.now(dt.timezone.utc).isoformat()}. Seed {SEED}, "
                 f"{N_BOOT} bootstrap resamples.\n\n")
        fh.write(f"**Threshold.** A cognate calcium site counts as disrupted when its "
                 f"bond-valence sum falls by more than **{T:.1f}%**, or when the ion no longer "
                 f"occupies the site in any model. That figure is the {NULL_PCTILE:.0f}th "
                 f"percentile of the absolute change measured at {len(null)} sites where no "
                 f"effect is possible — {len(byst)} bystander sites in the same constructs plus "
                 f"{len(neg_cog)} negative-control sites — not a chosen cut-off.\n\n")
        fh.write("## Rate\n\n| group | disrupted | rate | 95% CI | median ΔBVS |\n")
        fh.write("|---|---|---|---|---|\n")
        for k in order:
            s = cog[cog.klass == k]
            d, n = int(s.disrupted.sum()), len(s)
            lo, hi = wilson(d, n)
            fh.write(f"| {k} | {d}/{n} | {100*d/n:.0f}% | {100*lo:.0f}–{100*hi:.0f}% | "
                     f"{s.delta_pct.median():+.1f}% |\n")
        fh.write("\n| comparison | OR | p (Fisher) |\n|---|---|---|\n")
        for k, a, na, b, nb, orr, p in fisher_rows:
            fh.write(f"| pathogenic ligands vs {k} | {orr:.2f} | {p:.3g} |\n")
        fh.write(f"| pathogenic ligands vs all controls pooled | {orr_all:.2f} | {p_all:.3g} |\n")
        fh.write("\n## Magnitude\n\n| comparison | n | medians | Cliff's δ | 95% CI | p |\n")
        fh.write("|---|---|---|---|---|---|\n")
        for label, na_, nb_, ma, mb, d, lo, hi, p in mag_rows:
            fh.write(f"| {label} | {na_}/{nb_} | {ma:+.1f}% vs {mb:+.1f}% | {d:+.3f} | "
                     f"[{lo:+.3f}, {hi:+.3f}] | {p:.3g} |\n")
        fh.write("\n## The fold does not move\n\n")
        fh.write("Local backbone RMSD around the same calcium site, wild type to mutant:\n\n")
        fh.write("| group | median | max |\n|---|---|---|\n")
        for k in order:
            s = pd.to_numeric(cog[(cog.klass == k) & (cog.status == "present")].local_rmsd_A,
                              errors="coerce")
            fh.write(f"| {k} | {s.median():.2f} Å | {s.max():.2f} Å |\n")
        fh.write(f"\nPathogenic ligands against all controls: Cliff's δ {d_fold:+.3f} "
                 f"[{lo_f:+.3f}, {hi_f:+.3f}], p = {p_fold:.3g}. The whole-construct RMSD is "
                 "much larger and much more variable because these tandem constructs hinge "
                 "between domains between seeds; that is why the fold is judged locally.\n")
        fh.write("\n## Donor chemistry\n\n| subgroup | disrupted | rate | median ΔBVS |\n")
        fh.write("|---|---|---|---|\n")
        for label, s in (("keeps an oxygen donor", dk), ("removes the donor", dl)):
            fh.write(f"| {label} | {int(s.disrupted.sum())}/{len(s)} | "
                     f"{100*s.disrupted.mean():.0f}% | {s.delta_pct.median():+.1f}% |\n")
        fh.write(f"\nFisher's exact, donor removed against donor kept: OR {orr_d:.2f}, "
                 f"p = {p_d:.3g}.\n")
        fh.write("\n### Stratified by which consensus position is hit\n\n")
        fh.write("| consensus position | substitution | disrupted | median ΔBVS |\n")
        fh.write("|---|---|---|---|\n")
        for pos_class in ("Asp/Glu (offset 0, +3)", "beta-OH Asn (offset +16..+20)"):
            s_ = L[L["site_pos"] == pos_class]
            for keep_, lab in ((True, "keeps an O donor"), (False, "removes the donor")):
                t = s_[s_["donor_kept"] == keep_]
                if not len(t):
                    continue
                fh.write(f"| {pos_class} | {lab} | {int(t.disrupted.sum())}/{len(t)} | "
                         f"{t.delta_pct.median():+.1f}% |\n")
        fh.write("\nThe donor rule holds at the N-terminal Asp and Glu and fails at the "
                 "β-hydroxylation Asn, where a serine hydroxyl does not substitute for the "
                 "asparagine carboxamide. That position is the one for which the ClinGen FBN1 "
                 "VCEP exempts Asn→Ser from PM1.\n")
    log(f"wrote {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
