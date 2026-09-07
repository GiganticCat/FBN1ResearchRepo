#!/usr/bin/env python3
"""
17_fig5_calibration.py — PHASE 6, Figure 5. How much to trust a placed calcium.

An AlphaFold model contains no metals. Every calcium in a predicted structure here was PLACED,
by AF3 cofolding, and every structural metric that depends on the ion inherits that placement's
error. This figure is the error bar:

  A  placement accuracy   AF3 against X-ray and against NMR references, and homology transplant
  B  the superposition trap  whole-construct fitting inflates the apparent error ~4x
  C  CheckMyMetal          four site-quality parameters, experimental controls beside AF3 sites

Panel C's message is a NON-difference, which is the desired outcome — and it only means anything
because the experimental structures were submitted to CheckMyMetal as controls. Their scores land
in the same band as the models, well outside CMM's published "good" ranges, which says cbEGF
calcium sites are irregular rather than that the models are wrong. No "N% of sites are CMM
outliers" claim is made anywhere, and this panel is why.

Writes: figures/fig5_ca_placement.{png,pdf}, .numbers.json, logs/17_fig5_*.log
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig5_ca_placement"
XRAY_REFS = ("2W86", "1UZJ")
NMR_REFS = ("1LMJ", "1EMN")
CMM_PARAMS = [("grmsd", "gRMSD (°)", "{:.1f}"),
              ("valence", "Valence", "{:.2f}"),
              ("nvecsum", "nVECSUM", "{:.2f}"),
              ("coordination_number", "Coordination number", "{:.1f}")]


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(fs.SEED)
    log = [f"Phase 6 — Figure 5 — seed {fs.SEED}"]
    num = fs.Numbers(FIG, Path(__file__).name)

    calib = pd.read_csv(fs.RESULTS / "phase4_calibration.tsv", sep="\t")
    loo = pd.read_csv(fs.RESULTS / "phase4_transplant_loo.tsv", sep="\t")
    cmm = pd.read_csv(fs.PROC / "cmm_sites.tsv", sep="\t")
    stats = pd.read_csv(fs.RESULTS / "phase5_statistics.tsv", sep="\t")

    xray = calib[calib.reference_pdb.isin(XRAY_REFS)]
    nmr = calib[calib.reference_pdb.isin(NMR_REFS)]
    log.append(f"calibration: {len(calib)} model-domain comparisons "
               f"({len(xray)} X-ray, {len(nmr)} NMR); {len(loo)} leave-one-out transplants")

    fig = plt.figure(figsize=(7.4, 6.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], hspace=0.68, wspace=0.30)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    gsC = gs[1, :].subgridspec(1, 4, wspace=0.55)

    # ---- A: placement accuracy ------------------------------------------------------
    # Both AF3 groups share one colour on purpose: they are the same method, differing only
    # in the reference they are scored against, and the x labels already say which. Giving
    # them two blues would either fail the palette's chroma floor or fall under the
    # normal-vision separation floor — two steps of one hue cannot clear both.
    groups = [
        ("AF3 vs\nX-ray ref", xray.ca_deviation_A.values, fs.PALETTE["af3"]),
        ("AF3 vs\nNMR ref", nmr.ca_deviation_A.values, fs.PALETTE["af3"]),
        ("homology\ntransplant", loo.ca_transplant_error_A.values, fs.PALETTE["transplant"]),
    ]
    sA = fs.dist_panel(axA, groups, "Ca²⁺ displacement (Å)", rng, annotate_fmt="{:.2f}",
                       ylim=(0, 4.6))
    axA.set_title("How far a placed ion sits from the real one", loc="left", fontsize=8.5)
    axA.set_xlabel("each point is one cbEGF domain (A) or one donor→target pair (transplant)",
                   fontsize=6.4, color=fs.PALETTE["ink2"], labelpad=6)
    for key, lbl in (("af3_xray", "AF3 vs\nX-ray ref"), ("af3_nmr", "AF3 vs\nNMR ref"),
                     ("transplant", "homology\ntransplant")):
        num.set(f"ca_dev_median_{key}", sA[lbl]["median"])
        num.set(f"ca_dev_n_{key}", sA[lbl]["n"])
        log.append(f"  {key:11} n={sA[lbl]['n']:3} median={sA[lbl]['median']:.3f} A")
    ratio = sA["homology\ntransplant"]["median"] / sA["AF3 vs\nX-ray ref"]["median"]
    num.set("transplant_over_af3_xray_ratio", round(float(ratio), 2))
    # Placed below the n/median strip at the top of the axes, in the empty band above the
    # (very low) X-ray box.
    axA.annotate(f"AF3 is {ratio:.1f}× closer\nthan transplant", xy=(0.03, 0.74),
                 xycoords="axes fraction", fontsize=6.8, va="top", ha="left",
                 color=fs.PALETTE["af3"], fontweight="bold")

    # ---- B: the superposition trap --------------------------------------------------
    # Filled vs hollow, not two blues. The two series are the same method against different
    # reference types, and fill state is a print- and CVD-safe secondary encoding that costs
    # the palette nothing.
    for refs, marker, facecolor, lbl in (
            (XRAY_REFS, "o", fs.PALETTE["af3"], "X-ray reference"),
            (NMR_REFS, "^", "none", "NMR reference")):
        s = calib[calib.reference_pdb.isin(refs)]
        axB.scatter(s.domain_backbone_rmsd_A, s.construct_backbone_rmsd_A, s=20,
                    marker=marker, facecolors=facecolor, edgecolors=fs.PALETTE["af3"],
                    linewidths=1.0, alpha=0.85, label=lbl, zorder=3)
    lim = (0, max(calib.construct_backbone_rmsd_A.max(), calib.domain_backbone_rmsd_A.max()) * 1.1)
    axB.plot(lim, lim, ls=":", lw=0.8, color=fs.PALETTE["muted"], zorder=1)
    axB.set_xlim(*lim)
    axB.set_ylim(*lim)
    axB.set_xlabel("Backbone RMSD, per-domain fit (Å)")
    axB.set_ylabel("Backbone RMSD,\nwhole-construct fit (Å)")
    axB.set_title("Fit each domain, not the pair", loc="left", fontsize=8.5)
    axB.legend(loc="lower right", handletextpad=0.2, borderpad=0.2)
    axB.yaxis.grid(True, zorder=0)
    axB.xaxis.grid(True, zorder=0)
    axB.set_axisbelow(True)
    md, mc = float(calib.domain_backbone_rmsd_A.median()), float(
        calib.construct_backbone_rmsd_A.median())
    num.set("backbone_rmsd_median_per_domain", md)
    num.set("backbone_rmsd_median_whole_construct", mc)
    axB.annotate(f"median {md:.2f} Å per domain\nvs {mc:.2f} Å per construct —\n"
                 f"the gap is interdomain hinge motion,\nnot placement error",
                 xy=(0.04, 0.96), xycoords="axes fraction", fontsize=6.4, va="top",
                 ha="left", color=fs.PALETTE["ink2"])
    log.append(f"  backbone RMSD median: per-domain {md:.3f} A, whole-construct {mc:.3f} A")

    # ---- C: CheckMyMetal ------------------------------------------------------------
    exp = cmm[cmm.source == "experimental"]
    af3 = cmm[cmm.source != "experimental"]
    num.set("cmm_n_experimental", len(exp))
    num.set("cmm_n_af3", len(af3))
    for i, (col, label, fmt) in enumerate(CMM_PARAMS):
        ax = fig.add_subplot(gsC[0, i])
        for j, (name, sub, colour) in enumerate((("exp.", exp, fs.PALETTE["experimental"]),
                                                 ("AF3", af3, fs.PALETTE["af3"]))):
            vals = sub[col].astype(float).values
            ax.scatter(np.full(len(vals), j) + rng.uniform(-0.11, 0.11, len(vals)), vals,
                       s=16, color=colour, alpha=0.8, linewidths=0.4,
                       edgecolors=fs.PALETTE["surface"], zorder=3)
            ax.hlines(np.median(vals), j - 0.26, j + 0.26, color=colour, lw=2.0, zorder=4)
            num.set(f"cmm_{col}_median_{'experimental' if j == 0 else 'af3'}",
                    float(np.median(vals)))
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["exp.", "AF3"])
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel(label, fontsize=7)
        ax.yaxis.grid(True, zorder=0)
        ax.set_axisbelow(True)
        if i == 0:
            ax.text(-0.42, 1.16, "C", transform=ax.transAxes, fontsize=11, fontweight="bold",
                    va="bottom", ha="left", color=fs.PALETTE["ink"])
            ax.text(0.0, 1.16, "CheckMyMetal: AF3 calcium sites score like the experimental "
                               "controls", transform=ax.transAxes, fontsize=8.5,
                    va="bottom", ha="left", color=fs.PALETTE["ink"])
        log.append(f"  CMM {col:20} exp median {np.median(exp[col].astype(float)):.2f}  "
                   f"AF3 median {np.median(af3[col].astype(float)):.2f}")

    r_cmm = stats[stats.test == "CMM gRMSD: experimental vs AF3 sites"].iloc[0]
    num.set("cmm_grmsd_p", float(r_cmm.p_raw))
    fig.text(0.03, 0.045, f"gRMSD experimental vs AF3: Mann-Whitney p = {r_cmm.p_raw:.2f} — "
                          f"not significant, which is the desired outcome. "
                          f"n = {len(exp)} experimental and {len(af3)} AF3 sites.\n"
                          f"Both groups sit outside CheckMyMetal's published \"good\" bands: "
                          f"cbEGF sites are irregular and carbonyl-dominated, so the absolute "
                          f"thresholds do not transfer. This is why the controls were submitted.",
             fontsize=6.4, color=fs.PALETTE["ink2"], va="top", ha="left", linespacing=1.6)

    for ax, letter in ((axA, "A"), (axB, "B")):
        fs.panel_label(ax, letter, dx=-0.16, dy=1.04)
    fig.suptitle("Calcium placement is a modelling assumption — this is its error bar",
                 fontsize=10, x=0.02, ha="left", y=0.985)
    fig.subplots_adjust(bottom=0.17, top=0.90)

    fs.save(fig, FIG, num)
    plt.close(fig)
    print(f"  log -> {fs.write_log('17_fig5', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
