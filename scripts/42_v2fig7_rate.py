#!/usr/bin/env python3
"""
42_v2fig7_rate.py — v2 Figure 4 (final): how often a Ca2+-ligand variant breaks the site.

Replaces `33_v2fig2_geometry.py` as the manuscript's mechanism figure. That figure showed the
effect on nine variants against a seed-noise floor; this one measures a rate on 57, against a
null built from 123 measurements where no effect is possible.

a  every variant's change in site valence, sorted, with the measured null band behind it
b  the disruption rate per group, Wilson intervals
c  the chemistry that decides it: substitutions keeping an oxygen donor against those removing it
d  the same variants' Rosetta ddG against their valence change — the two axes do not track

The N2144 position carries the study's only experimental anchor (Kettle 1999, ~9x weaker Ca2+
binding with the fold unchanged); N2144H is folded here and marked in a.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig7_rate"
SCRIPT = Path(__file__).name

KLASS_ORDER = ["ca_ligand_pathogenic", "composition_control", "negative_control"]
KLASS_LABEL = {
    "ca_ligand_pathogenic": "Pathogenic Ca²⁺ ligand",
    "composition_control": "D/E/N control",
    "negative_control": "Negative control",
}
KLASS_COLOR = {
    "ca_ligand_pathogenic": F.PALETTE["ca_site"],
    "composition_control": F.PALETTE["other_cbegf"],
    "negative_control": F.PALETTE["af3"],
}


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(F.SEED)
    N = F.Numbers(NAME, SCRIPT)

    r = pd.read_csv(F.PROC / "batch4_rate.tsv", sep="\t")
    T = float(r.null_threshold_pct.iloc[0])
    N.set("null_threshold_pct", T)
    N.set("n_variants", int(len(r)))

    fig = plt.figure(figsize=(F.COL_2, 4.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.88], hspace=0.30, wspace=0.38,
                          left=0.070, right=0.988, top=0.95, bottom=0.13)

    # ---------------------------------------------------------------- a  waterfall
    ax = fig.add_subplot(gs[0, :])
    d = r.sort_values("delta_bvs_pct").reset_index(drop=True)
    x = np.arange(len(d))
    ax.axhspan(-T, T, color=F.PALETTE["grid"], zorder=0)
    ax.axhline(0, color=F.PALETTE["muted"], lw=0.5, zorder=1)
    for k in KLASS_ORDER:
        m = d.klass == k
        ax.bar(x[m], d.delta_bvs_pct[m], width=0.74, color=KLASS_COLOR[k],
               linewidth=0.0, zorder=2,
               label=f"{KLASS_LABEL[k]} (n={int(m.sum())})")
    ax.set_xlim(-0.8, len(d) - 0.2)
    ax.set_xticks([])
    ax.set_xlabel("Variant, ranked by change in site valence", labelpad=1)
    ax.set_ylabel("Δ bond-valence sum (%)")
    ax.annotate(f"null band, ±{T:.1f}%\n(95th pct of {123} sites where no effect is possible)",
                xy=(len(d) - 1.5, T), xytext=(-2, 4), textcoords="offset points",
                ha="right", va="bottom", fontsize=5.5, color=F.PALETTE["ink2"],
                linespacing=1.2)
    anchor = d.index[d.variant == "N2144H"]
    if len(anchor):
        i = int(anchor[0])
        ax.annotate("N2144H\n(N2144S: Kd ↑ ~9×, fold unchanged)",
                    xy=(i, d.delta_bvs_pct[i]), xytext=(6, 10), textcoords="offset points",
                    fontsize=5.5, color=F.PALETTE["ink"], linespacing=1.2,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=F.PALETTE["ink2"]))
        N.set("N2144H_delta_pct", float(d.delta_bvs_pct[i]))
    ax.legend(loc="lower right", ncol=1, fontsize=5.5)
    F.tidy(ax, "y")
    F.panel_label(ax, "a", dx=-0.055)

    # ---------------------------------------------------------------- b  rate
    ax = fig.add_subplot(gs[1, 0])
    for i, k in enumerate(KLASS_ORDER):
        s = r[r.klass == k]
        n, kk = len(s), int(s.disrupted.sum())
        lo, hi = wilson(kk, n)
        ax.bar(i, 100 * kk / n, width=0.56, color=KLASS_COLOR[k], linewidth=0.0, zorder=2)
        ax.errorbar(i, 100 * kk / n, yerr=[[100 * kk / n - 100 * lo], [100 * hi - 100 * kk / n]],
                    fmt="none", ecolor=F.PALETTE["ink"], elinewidth=0.6, capsize=1.8, zorder=3)
        ax.annotate(f"{kk}/{n}", xy=(i, 100 * hi), xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=5.5, color=F.PALETTE["ink2"])
        N.set(f"rate_{k}", round(100 * kk / n, 1))
        N.set(f"n_{k}", n)
        N.set(f"disrupted_{k}", kk)
    ax.set_xticks(range(len(KLASS_ORDER)))
    ax.set_xticklabels([KLASS_LABEL[k] for k in KLASS_ORDER], rotation=24, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel("Sites disrupted (%)")
    ax.set_ylim(0, 100)
    F.tidy(ax, "y")
    F.panel_label(ax, "b")

    # ---------------------------------------------------------------- c  donor chemistry
    ax = fig.add_subplot(gs[1, 1])
    L = r[r.klass == "ca_ligand_pathogenic"]
    groups = []
    for keep, lab in ((True, "keeps O donor"), (False, "removes donor")):
        s = L[L.donor_kept == keep]
        groups.append((lab, s.delta_bvs_pct.to_numpy(),
                       F.PALETTE["ca_site"] if not keep else F.PALETTE["other_cbegf"], s))
    ax.axhspan(-T, T, color=F.PALETTE["grid"], zorder=0)
    ax.axhline(0, color=F.PALETTE["muted"], lw=0.5, zorder=1)
    for i, (lab, vals, col, s) in enumerate(groups):
        bp = ax.boxplot([vals], positions=[i], widths=0.36, showfliers=False,
                        patch_artist=True, zorder=2)
        for b in bp["boxes"]:
            b.set(facecolor=col, alpha=0.16, edgecolor=col, linewidth=0.7)
        for part in ("whiskers", "caps"):
            for a_ in bp[part]:
                a_.set(color=col, linewidth=0.6)
        for m_ in bp["medians"]:
            m_.set(color=col, linewidth=1.4)
        ax.scatter(np.full(len(vals), i) + rng.uniform(-0.1, 0.1, len(vals)), vals,
                   s=5, color=col, alpha=0.75, linewidths=0.0, zorder=3)
        kk = int(s.disrupted.sum())
        ax.annotate(f"{kk}/{len(s)}\ndisrupted", xy=(i, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(0, 1), textcoords="offset points", ha="center", va="bottom",
                    fontsize=5.5, color=F.PALETTE["ink2"], linespacing=1.15)
        key = "donor_kept" if lab.startswith("keeps") else "donor_removed"
        N.set(f"{key}_n", int(len(s)))
        N.set(f"{key}_disrupted", kk)
        N.set(f"{key}_median_delta", round(float(np.median(vals)), 1))
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_xlim(-0.55, 1.55)
    ax.set_ylabel("Δ bond-valence sum (%)")
    ax.set_xlabel("Pathogenic Ca²⁺-ligand substitution", labelpad=2)
    F.tidy(ax, "y")
    F.panel_label(ax, "c")

    # ---------------------------------------------------------------- d  geometry vs energetics
    ros = pd.read_csv(F.PROC / "rosetta_ddg.tsv", sep="\t")
    ros = ros[ros.status == "ok"][["variation_id", "ddg_mean"]]
    m = L.merge(ros, on="variation_id", how="left").dropna(subset=["ddg_mean"])
    ax = fig.add_subplot(gs[1, 2])
    ax.axvspan(-T, T, color=F.PALETTE["grid"], zorder=0)
    for flag, face, lab in ((True, F.PALETTE["ca_site"], "site disrupted"),
                            (False, "none", "site intact")):
        sub = m[m.disrupted == flag]
        ax.scatter(sub.delta_bvs_pct, sub.ddg_mean, s=11,
                   facecolors=face, edgecolors=F.PALETTE["ca_site"],
                   linewidths=0.6, alpha=0.9, zorder=3, label=lab)
    ax.set_xlabel("Δ bond-valence sum (%)", labelpad=2)
    ax.set_ylabel("Rosetta ΔΔG (REU)")
    ax.legend(loc="upper left", fontsize=5.5, handletextpad=0.3)
    F.tidy(ax, "y")
    F.panel_label(ax, "d")

    N.set("spearman_n", int(len(m)))
    from scipy import stats as sstats
    rho, p = sstats.spearmanr(m.delta_bvs_pct, m.ddg_mean)
    N.set("spearman_rho", round(float(rho), 3))
    N.set("spearman_p", float(f"{p:.4g}"))
    N.set("ddg_median_disrupted", round(float(m[m.disrupted].ddg_mean.median()), 2))
    N.set("ddg_median_intact", round(float(m[~m.disrupted].ddg_mean.median()), 2))
    N.set("median_delta_ca_ligand", round(float(L.delta_bvs_pct.median()), 1))
    N.set("median_delta_composition", round(float(r[r.klass == "composition_control"].delta_bvs_pct.median()), 1))
    N.set("median_delta_negative", round(float(r[r.klass == "negative_control"].delta_bvs_pct.median()), 1))
    N.set("local_rmsd_median_ca_ligand", round(float(L.local_rmsd_A.median()), 2))
    N.set("local_rmsd_median_controls",
          round(float(r[r.klass != "ca_ligand_pathogenic"].local_rmsd_A.median()), 2))

    F.save(fig, NAME, N)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
