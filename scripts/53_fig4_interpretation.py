#!/usr/bin/env python3
"""
53_fig4_interpretation.py — Figure 4: the damage is invisible to the methods used to detect it.

The consequence of Figures 2 and 3 for how a variant is actually interpreted in a clinic. A
stability calculation is blind to a functional metal contact, and a sequence-based predictor
saturates: both call every calcium-ligand substitution damaging, and neither separates the ones
that measurably break the site from the ones that do not.

a  Rosetta ddG against the measured change in coordination, for the same 33 variants
b  AlphaMissense for those variants, split by whether the site was disrupted
c  enrichment of pathogenic variants at each feature against a population comparison group,
   with the circularity that qualifies it written on the panel

Writes: figures/fig4_interpretation.{png,pdf}, .numbers.json,
        logs/53_fig4_interpretation_<stamp>.log
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

NAME = "fig4_interpretation"
SCRIPT = Path(__file__).name

FEATURES = [
    ("side-chain Ca ligand", "Ca²⁺ ligand\n(side chain)", P.C["ink"]),
    ("backbone-only Ca ligand", "Ca²⁺ ligand\n(main-chain carbonyl only)", P.C["grey"]),
    ("motif-overcalled (not a ligand)", "Consensus motif,\nnon-ligand", P.C["grey"]),
    ("cysteine position", "Cysteine position", P.C["ink"]),
]


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(P.SEED)
    N = P.Numbers(NAME, SCRIPT)
    log = [f"Figure 4 — {NAME}"]

    df = P.load()
    rate = pd.read_csv(P.PROC / "batch4_rate.tsv", sep="\t")
    ros = pd.read_csv(P.PROC / "rosetta_ddg.tsv", sep="\t").query("status == 'ok'")
    L = (rate[rate.klass == "ca_ligand_pathogenic"]
         .merge(ros[["variation_id", "ddg_mean"]], on="variation_id", how="left")
         .merge(df[["variation_id", "am"]], on="variation_id", how="left"))
    T = float(rate.null_threshold_pct.iloc[0])

    fig = plt.figure(figsize=(P.W2, 2.48))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.78, 1.42], wspace=0.42,
                          left=0.062, right=0.988, top=0.82, bottom=0.29)

    # ------------------------------------------------------------------ a  ddG vs geometry
    ax = fig.add_subplot(gs[0, 0])
    sub = L.dropna(subset=["ddg_mean", "delta_bvs_pct"])
    rho, p = sstats.spearmanr(sub.delta_bvs_pct, sub.ddg_mean)
    ax.axvspan(-T, T, color=P.C["band"], zorder=0)
    for flag, face in ((True, P.C["ca"]), (False, "white")):
        s = sub[sub.disrupted == flag]
        ax.scatter(s.delta_bvs_pct, s.ddg_mean, s=13, facecolors=face,
                   edgecolors=P.C["ca"], linewidths=0.7, zorder=3)
    ax.set_xlabel("Change in Ca²⁺ valence (%)", labelpad=1.5)
    ax.set_ylabel("ΔΔG (Rosetta energy units)")
    ax.set_xlim(-44, 14)
    # headroom below the lowest variant so the three-line correlation label sits on empty axes
    ax.set_ylim(-9.4, 10.2)
    q = P.agrees(p, "Spearman(delta BVS %, Rosetta ddG) across pathogenic Ca ligands")
    # bottom right, which is the one corner of this panel with no points in it; the two-line
    # P/q label no longer fits bottom left without sitting on the most negative variant
    ax.annotate(f"Spearman ρ = {rho:+.2f}\n{P.fmt_pq(p, q)}", xy=(0.985, 0.03),
                xycoords="axes fraction", ha="right", va="bottom", fontsize=6.2,
                color=P.C["ink"], linespacing=1.3)
    P.label_at(ax, 0.03, 0.965, "filled, coordination reduced", P.C["ink"],
               transform=ax.transAxes, va="top")
    P.label_at(ax, 0.03, 0.885, "open, coordination retained", P.C["ink2"],
               transform=ax.transAxes, va="top")
    P.title(ax, "Rosetta ΔΔG versus Ca²⁺ valence")
    P.panel(ax, "a", dx=-0.24)
    N.set("spearman_rho", round(float(rho), 3))
    N.set("spearman_p", float(f"{p:.3g}"))
    N.set("spearman_n", int(len(sub)))
    N.set("spearman_q", float(f"{q:.3g}"))
    log.append(f"  a  rho={rho:.3f} p={p:.3g} n={len(sub)}")

    # ------------------------------------------------------------------ b  AlphaMissense
    ax = fig.add_subplot(gs[0, 1])
    # "Coordination reduced" / "Coordination retained" will not fit as two tick labels across a
    # third-width panel, so the noun goes on the axis label and the ticks carry the outcome.
    groups = [("Reduced", L.loc[L.disrupted, "am"].to_numpy(), P.C["ink"]),
              ("Retained", L.loc[~L.disrupted, "am"].to_numpy(), "open")]
    st = P.strip(ax, groups, "AlphaMissense score", rng, ms=7.0, alpha=0.8, width=0.40,
                 show_n=False)
    ax.set_ylim(0.55, 1.135)
    ax.axhline(0.564, color=P.C["rule"], lw=0.6, ls=(0, (2.5, 2)), zorder=1)
    a = L.loc[L.disrupted, "am"].dropna()
    b = L.loc[~L.disrupted, "am"].dropna()
    p_am = float(sstats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
    q_am = P.agrees(p_am, "AlphaMissense: disrupted vs intact Ca-ligand variants")
    N.set("b_am_q", float(f"{q_am:.3g}"))
    P.bracket(ax, 0, 1, 1.028, P.fmt_pq(p_am, q_am), dy=0.009)
    for k, v in st.items():
        N.set(f"b_am_{k.replace(chr(10),'_').lower()}_median", round(v["median"], 3))
        N.set(f"b_am_{k.replace(chr(10),'_').lower()}_n", v["n"])
    N.set("b_am_p", float(f"{p_am:.3g}"))
    ax.set_xticklabels([f"{lab}\nn = {st[lab]['n']}" for lab, _, _ in groups])
    ax.set_xlabel("Ca²⁺ coordination in the mutant model", labelpad=1.5)
    P.title(ax, "AlphaMissense score, same variants")
    P.panel(ax, "b", dx=-0.36)

    # ------------------------------------------------------------------ c  enrichment
    ax = fig.add_subplot(gs[0, 2])
    enr = pd.read_csv(P.PROC / "enrichment_v2.tsv", sep="\t")
    enr = enr[enr.denominator == "population (gnomAD)"]
    ys = np.arange(len(FEATURES))[::-1]
    for y, (feat, lab, col) in zip(ys, FEATURES):
        row = enr[enr.feature == feat]
        if row.empty:
            raise SystemExit(f"no enrichment row for {feat!r}")
        row = row.iloc[0]
        orr, lo, hi = float(row.odds_ratio_corrected), float(row.ci_lo), float(row.ci_hi)
        ax.plot([lo, hi], [y, y], color=col, lw=1.1, solid_capstyle="round", zorder=3)
        ax.plot([orr], [y], marker="o", ms=4.0, mfc=col, mec="white", mew=0.6, zorder=4)
        ax.annotate(f"{orr:.3g}" if orr < 10 else f"{orr:.0f}",
                    xy=(hi, y), xytext=(4, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=6, color=col)
        N.set(f"or_{feat.replace(' ', '_')}", round(orr, 3))
    ax.axvline(1, color=P.C["ink"], lw=0.5, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(0.012, 900)
    ax.set_yticks(ys)
    ax.set_yticklabels([f[1] for f in FEATURES])
    ax.set_ylim(-0.7, len(FEATURES) - 0.3)
    ax.set_xlabel("Odds ratio, pathogenic vs 1,901 gnomAD variants (95% CI)", labelpad=1.5)
    ax.annotate("depleted", xy=(0.06, len(FEATURES) - 0.55), ha="center", va="center",
                fontsize=5.8, color=P.C["ink2"])
    ax.annotate("enriched", xy=(20, len(FEATURES) - 0.55), ha="center", va="center",
                fontsize=5.8, color=P.C["ink2"])
    ax.annotate("Descriptive, and partly circular, because the ClinGen FBN1 rules assign\n"
                "pathogenicity partly on position at these features. It confirms the two\n"
                "features carry the variants and is not evidence of mechanism. No amino-acid\n"
                "substitution can remove a main-chain carbonyl, so that row is a negative control.",
                xy=(0.0, -0.30), xycoords="axes fraction", ha="left", va="top", fontsize=5.4,
                color=P.C["ink2"], linespacing=1.35)
    P.title(ax, "Enrichment relative to gnomAD v4")
    P.panel(ax, "c", dx=-0.20)

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("53_fig4_interpretation", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
