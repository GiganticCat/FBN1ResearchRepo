#!/usr/bin/env python3
"""
52_fig3_coordination.py — Figure 3: Ca2+-ligand variants loosen the ion while the fold stays put.

Replaces `42_v2fig7_rate.py`, which had the right content and the wrong presentation: a legend
box in a corner instead of labels on the bars, a rate panel whose comparison was left for the
reader to make, and no measurement of the fold beside the measurement of the site -- even though
"the site loosens and the fold does not" is the entire claim.

a  every variant's change in site valence, ranked, against the null band measured from 123 sites
   where no effect is possible
b  the disruption rate per group, with the comparison drawn
c  what decides it: whether the substitution leaves an oxygen donor behind, and at which
   consensus position — the rule holds at the N-terminal Asp and Glu and fails at the
   beta-hydroxylation Asn
d  the fold, measured on the same models and the same site neighbourhoods

Writes: figures/fig3_coordination.{png,pdf}, .numbers.json, logs/52_fig3_coordination_<stamp>.log
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

NAME = "fig3_coordination"
SCRIPT = Path(__file__).name

KORDER = ["ca_ligand_pathogenic", "composition_control", "negative_control"]
# Full names appear once, in the key on panel a. Panels b and d are narrow, so they carry short
# forms; spelling them out there collides them into each other.
KLAB = {"ca_ligand_pathogenic": "Pathogenic Ca²⁺ ligand",
        "composition_control": "Non-ligand Asp/Glu/Asn",
        "negative_control": "Non-ligand control"}
KSHORT = {"ca_ligand_pathogenic": "Pathogenic\nCa²⁺ ligand",
          "composition_control": "Non-ligand\nAsp/Glu/Asn", "negative_control": "Non-ligand\ncontrol"}
# Three series in one ranked panel do need separating, and fill does it: solid, grey, open.
# No hue is spent on identity.
KCOL = {"ca_ligand_pathogenic": P.C["ink"], "composition_control": P.C["grey"],
        "negative_control": "none"}
KEDGE = {"ca_ligand_pathogenic": P.C["ink"], "composition_control": P.C["grey"],
         "negative_control": P.C["ink"]}


def wilson(k, n):
    z, p = 1.959964, k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(P.SEED)
    N = P.Numbers(NAME, SCRIPT)
    log = [f"Figure 3 — {NAME}"]

    r = pd.read_csv(P.PROC / "batch4_rate.tsv", sep="\t")
    g = pd.read_csv(P.PROC / "af3_batch4_geometry.tsv", sep="\t")
    T = float(r.null_threshold_pct.iloc[0])
    n_null = int(((~g.is_cognate) & g.wt_site_ok & (g.status == "present")).sum()
                 + ((g.is_cognate) & g.wt_site_ok & (g.klass == "negative_control")
                    & (g.status == "present")).sum())
    N.set("null_threshold_pct", round(T, 2))
    N.set("n_null_sites", n_null)
    N.set("n_variants", int(len(r)))

    fig = plt.figure(figsize=(P.W2, 4.35))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.92], hspace=0.34, wspace=0.38,
                          width_ratios=[1.02, 1.48, 1.02],
                          left=0.062, right=0.990, top=0.935, bottom=0.125)

    # ------------------------------------------------------------------ a  ranked change
    ax = fig.add_subplot(gs[0, :])
    d = r.sort_values("delta_bvs_pct").reset_index(drop=True)
    x = np.arange(len(d))
    ax.axhspan(-T, T, color=P.C["band"], zorder=0)
    ax.axhline(0, color=P.C["ink"], lw=0.4, zorder=1)
    for k in KORDER:
        m = d.klass == k
        ax.bar(x[m], d.delta_bvs_pct[m], width=0.76, facecolor=KCOL[k],
               edgecolor=KEDGE[k], linewidth=0.35, zorder=2)
    ax.set_xlim(-0.9, len(d) - 0.1)
    ax.set_xticks([])
    ax.set_ylim(-42, 13)
    ax.set_ylabel("Change in Ca²⁺ bond-valence\nsum, mutant − wild type (%)", labelpad=2)
    ax.set_xlabel(f"{len(d)} cofolded variants, ranked by change in valence", labelpad=-6)
    ax.spines["bottom"].set_visible(False)

    # the null band named on the band itself
    ax.annotate(f"Shaded band ±{T:.1f}%, the 95th percentile of absolute change across\n"
                f"{n_null} sites at which no effect is possible (109 bystander sites and\n"
                f"14 negative controls)",
                xy=(len(d) - 1, -T), xytext=(-2, -2), textcoords="offset points",
                ha="right", va="top", fontsize=5.6, color=P.C["ink2"], linespacing=1.3)
    from matplotlib.patches import Rectangle as _R
    for k, yf in (("ca_ligand_pathogenic", 0.955), ("composition_control", 0.885),
                  ("negative_control", 0.815)):
        n_k = int((d.klass == k).sum())
        # a swatch carries the fill; the text stays in ink so it never sits grey-on-grey
        ax.add_patch(_R((0.012, yf - 0.024), 0.011, 0.048, transform=ax.transAxes,
                        facecolor=KCOL[k], edgecolor=KEDGE[k], linewidth=0.5,
                        clip_on=False, zorder=6))
        P.label_at(ax, 0.030, yf, f"{KLAB[k].replace(chr(10), ' ')}  (n = {n_k})", P.C["ink"],
                   transform=ax.transAxes, va="center")
    i = int(d.index[d.variant == "N2144H"][0])
    ax.annotate("N2144H",
                xy=(i, d.delta_bvs_pct[i]), xytext=(26, 6), textcoords="offset points",
                fontsize=5.8, color=P.C["ink"], linespacing=1.35, va="bottom", ha="left",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=P.C["ink2"],
                                shrinkA=1, shrinkB=2))
    N.set("N2144H_delta_pct", float(d.delta_bvs_pct[i]))
    P.panel(ax, "a", dx=-0.048, dy=1.01)

    # ------------------------------------------------------------------ b  rate
    # Horizontal, for the same reason Figure 2 is: three category names do not fit side by side
    # under a panel this narrow, and abbreviating them past recognition is the thing to avoid.
    ax = fig.add_subplot(gs[1, 0])
    for i, k in enumerate(KORDER):
        y = len(KORDER) - 1 - i
        s_ = r[r.klass == k]
        n, kk = len(s_), int(s_.disrupted.sum())
        lo, hi = wilson(kk, n)
        ax.barh(y, 100 * kk / n, height=0.55, facecolor=KCOL[k], edgecolor=KEDGE[k],
                linewidth=0.5, zorder=2)
        ax.errorbar(100 * kk / n, y, xerr=[[100 * (kk / n - lo)], [100 * (hi - kk / n)]],
                    fmt="none", ecolor=P.C["ink"], elinewidth=0.6, capsize=1.6, zorder=3)
        # inside the bar when there is room, outside it when there is not
        inside = 100 * kk / n > 22
        ax.annotate(f"{kk}/{n}", xy=(2.0 if inside else 100 * hi + 2.5, y),
                    ha="left", va="center", fontsize=6,
                    color="white" if inside else P.C["ink2"])
        N.set(f"rate_{k}", round(100 * kk / n, 1))
        N.set(f"n_{k}", n)
        N.set(f"disrupted_{k}", kk)
    a_, na = int(r[r.klass == KORDER[0]].disrupted.sum()), int((r.klass == KORDER[0]).sum())
    b_ = int(r[r.klass != KORDER[0]].disrupted.sum())
    nb = int((r.klass != KORDER[0]).sum())
    orr, p = sstats.fisher_exact([[a_, na - a_], [b_, nb - b_]])
    N.set("fisher_or", round(float(orr), 2))
    N.set("fisher_p", float(f"{p:.3g}"))
    ax.set_yticks(range(len(KORDER)))
    ax.set_yticklabels([KSHORT[k] for k in reversed(KORDER)])
    ax.set_xlabel("Ca²⁺ sites with reduced\ncoordination (%)", labelpad=1.5)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.62, len(KORDER) - 0.38)
    # the q is the Benjamini-Hochberg value of this same test in results/statistics_final.tsv,
    # and `agrees` stops the build if the panel's own Fisher test has drifted from that row
    q = P.agrees(p, "disruption rate: pathogenic Ca ligands vs all controls")
    N.set("fisher_q", float(f"{q:.3g}"))
    P.bracket_h(ax, 2, 0.5, 72, f"OR {orr:.0f}\n" + P.fmt_pq(p, q), dx=3.5)
    P.title(ax, "Sites with reduced coordination")
    P.panel(ax, "b", dx=-0.52)

    # ------------------------------------------------------------------ c  donor chemistry
    ax = fig.add_subplot(gs[1, 1])
    L = r[r.klass == "ca_ligand_pathogenic"].copy()
    ax.axhspan(-T, T, color=P.C["band"], zorder=0)
    ax.axhline(0, color=P.C["ink"], lw=0.4, zorder=1)
    cats = [("Asp/Glu", True, "Donor\nretained"), ("Asp/Glu", False, "Donor\nlost"),
            ("beta-OH Asn", True, "Donor\nretained"), ("beta-OH Asn", False, "Donor\nlost")]
    xs, subs = [], []
    for i, (posc, keep, lab) in enumerate(cats):
        x = i + (0.55 if i >= 2 else 0.0)
        xs.append(x)
        s_ = L[(L.site_pos == posc) & (L.donor_kept == keep)]
        subs.append(s_)
        col = P.C["grey"] if keep else P.C["ink"]
        v = s_.delta_bvs_pct.to_numpy()
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.scatter(np.full(v.size, x) + rng.uniform(-0.12, 0.12, v.size), v, s=6.5,
                   color=col, alpha=0.8, linewidths=0, zorder=3)
        ax.add_patch(P.plt_rect(x - 0.19, q1, 0.38, q3 - q1, col))
        ax.plot([x - 0.26, x + 0.26], [med, med], lw=1.6, color=col, zorder=4)
        ax.annotate(f"{int(s_.disrupted.sum())}/{len(s_)}", xy=(x, 0.005),
                    xycoords=("data", "axes fraction"), xytext=(0, 1),
                    textcoords="offset points", ha="center", va="bottom", fontsize=5.8,
                    color=col)
        key = f"c_{posc.replace('-','').replace(' ','_')}_{'kept' if keep else 'removed'}"
        N.set(f"{key}_n", int(len(s_)))
        N.set(f"{key}_disrupted", int(s_.disrupted.sum()))
        N.set(f"{key}_median_delta", round(float(med), 1))
    # the three Asn->Ser variants are the ones the ClinGen FBN1 panel exempts from PM1
    # These were ringed in red in an earlier version. The ring alone does the job: it is the
    # only ring in the panel, and the points it encircles are grey discs, so it survives
    # desaturation with more contrast than the red did.
    ns = L[(L.wt_aa == "N") & (L.mut_aa == "S") & (L.site_pos == "beta-OH Asn")]
    ax.scatter(np.full(len(ns), xs[2]) + rng.uniform(-0.12, 0.12, len(ns)),
               ns.delta_bvs_pct, s=30, facecolors="none", edgecolors=P.C["flag"],
               linewidths=1.0, zorder=5)
    # anchored clear of the jitter band: the rings are the mark, the text must not sit on them
    ax.annotate("Asn→Ser, exempted\nfrom PM1 by the\nClinGen FBN1 panel",
                xy=(xs[2] + 0.22, float(ns.delta_bvs_pct.median())),
                xytext=(xs[2] - 0.30, 15.5), textcoords="data",
                ha="left", va="top", fontsize=5.5,
                color=P.C["flag"], linespacing=1.3,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=P.C["flag"],
                                shrinkA=2, shrinkB=1))
    N.set("asn_ser_n", int(len(ns)))
    N.set("asn_ser_disrupted", int(ns.disrupted.sum()))
    N.set("asn_ser_median_delta", round(float(ns.delta_bvs_pct.median()), 1))
    orr_d, p_d = sstats.fisher_exact(
        [[int(subs[1].disrupted.sum()), len(subs[1]) - int(subs[1].disrupted.sum())],
         [int(subs[0].disrupted.sum()), len(subs[0]) - int(subs[0].disrupted.sum())]])
    N.set("donor_fisher_p", float(f"{p_d:.3g}"))
    ax.set_xticks(xs)
    ax.set_xticklabels([c[2] for c in cats])
    ax.set_xlim(-0.6, xs[-1] + 0.6)
    ax.set_ylim(-42, 17)
    ax.set_ylabel("Change in valence (%)", labelpad=2)
    for x0, x1, lab in ((xs[0], xs[1], "Asp / Glu\n(offsets 0, +3)"),
                        (xs[2], xs[3], "β-OH Asn\n(offsets +16 to +20)")):
        ax.annotate(lab, xy=((x0 + x1) / 2, -0.245), xycoords=("data", "axes fraction"),
                    ha="center", va="top", fontsize=6.2, color=P.C["ink"], linespacing=1.25)
        ax.plot([x0 - 0.32, x1 + 0.32], [-0.205, -0.205], transform=ax.get_xaxis_transform(),
                lw=0.6, color=P.C["ink2"], clip_on=False)
    q_d = P.agrees(p_d, "disruption rate: donor removed vs retained, at the Asp/Glu positions")
    N.set("donor_fisher_q", float(f"{q_d:.3g}"))
    P.bracket(ax, xs[0], xs[1], 8.5, P.fmt_pq(p_d, q_d), dy=1.8)
    P.title(ax, "Change in valence by oxygen donor and consensus position")
    P.panel(ax, "c", dx=-0.155)

    # ------------------------------------------------------------------ d  the fold
    ax = fig.add_subplot(gs[1, 2])
    groups = [(KSHORT[k], r.loc[r.klass == k, "local_rmsd_A"].to_numpy(),
               "open" if k == "negative_control" else KCOL[k]) for k in KORDER]
    st = P.striph(ax, groups, "Backbone rmsd within\n12 Å of the ion (Å)", rng,
                  ms=6.5, alpha=0.85, height=0.40, show_n=False)
    ax.set_xlim(0, 1.0)
    for k, v in st.items():
        N.set(f"d_rmsd_{k.replace(chr(10), ' ')}_median", round(v["median"], 3))
    N.set("d_rmsd_max", round(float(r.local_rmsd_A.max()), 3))
    P.title(ax, "Local backbone rmsd")
    P.panel(ax, "d", dx=-0.52)

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("52_fig3_coordination", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
