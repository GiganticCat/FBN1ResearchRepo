#!/usr/bin/env python3
"""
54_figS1_validation.py — Supplementary Figure S1: the modelled calcium sites are trustworthy.

Every geometric statement in Figure 3 depends on two things being true: that AlphaFold3 puts the
ion where an experiment would, and that the local bond-valence sum tracks the established
metal-site validation service. Both are testable, and both are tested here rather than asserted.

a  modelled ion position against the experimental one, for the eight cbEGF domains with a
   calcium-bound structure, plotted against the backbone deviation of the same domain
b  local bond-valence sum against CheckMyMetal, over the 16 sites that service has scored
c  the nine-variant pilot series that established the measurement before batch 4 scaled it

Writes: figures/figS1_validation.{png,pdf}, .numbers.json, logs/54_figS1_validation_<stamp>.log
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

NAME = "figS1_validation"
SCRIPT = Path(__file__).name

# Batch-3 job names carry construct-local numbering; these are the P35555 substitutions.
PILOT = {"N18S": ("N2144S", "Ca²⁺ ligand"), "N57S": ("N2183S", "Ca²⁺ ligand"),
         "E4K": ("E2130K", "Ca²⁺ ligand"), "D45G": ("D1113G", "Ca²⁺ ligand"),
         "E5K": ("E1073K", "Ca²⁺ ligand"), "D2G": ("D1487G", "Ca²⁺ ligand"),
         "E107K": ("E913K", "Ca²⁺ ligand"), "C70F": ("C1138F", "control, S–S loss"),
         "P73L": ("P1141L", "control, benign")}


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(P.SEED)
    N = P.Numbers(NAME, SCRIPT)
    log = [f"Supplementary Figure S1 — {NAME}"]

    fig = plt.figure(figsize=(P.W2, 2.40))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25], wspace=0.40,
                          left=0.058, right=0.988, top=0.82, bottom=0.20)

    # ------------------------------------------------------------------ a  placement
    ax = fig.add_subplot(gs[0, 0])
    cal = pd.read_csv(P.RESULTS / "phase4_calibration.tsv", sep="\t")
    c0 = cal[cal.model == 0]
    lim = float(max(c0.domain_backbone_rmsd_A.max(), c0.ca_deviation_A.max())) * 1.12
    # shade the half-plane where the ion is placed more accurately than the fold around it,
    # which is the claim being made; all eight domains fall in it
    ax.fill_between([0, lim], [0, 0], [0, lim], color=P.C["band"], zorder=0)
    ax.plot([0, lim], [0, lim], color=P.C["muted"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.scatter(c0.domain_backbone_rmsd_A, c0.ca_deviation_A, s=14, color=P.C["ca"],
               linewidths=0, zorder=3)
    for _, r in c0.iterrows():
        ax.annotate(r.domain.replace("cbEGF", ""), xy=(r.domain_backbone_rmsd_A, r.ca_deviation_A),
                    xytext=(3, 1.5), textcoords="offset points", fontsize=5.4, color=P.C["ink2"])
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Domain backbone r.m.s.d. (Å)", labelpad=1.5)
    ax.set_ylabel("Ca²⁺ position deviation (Å)")
    ax.annotate("Ca²⁺ placed more\naccurately than\nthe backbone", xy=(lim * 0.60, lim * 0.16),
                ha="center", va="center", fontsize=5.6, color=P.C["ink2"], linespacing=1.3)
    med_ca, med_bb = float(c0.ca_deviation_A.median()), float(c0.domain_backbone_rmsd_A.median())
    ax.annotate(f"median {med_ca:.2f} Å vs {med_bb:.2f} Å\n"
                f"{int((c0.ca_deviation_A < c0.domain_backbone_rmsd_A).sum())} of {len(c0)} domains",
                xy=(0.03, 0.97), xycoords="axes fraction", ha="left", va="top",
                fontsize=5.8, color=P.C["ink"], linespacing=1.3)
    N.set("a_median_ca_deviation_A", round(med_ca, 3))
    N.set("a_median_backbone_rmsd_A", round(med_bb, 3))
    N.set("a_n_domains", int(len(c0)))
    P.title(ax, "Modeled versus experimental Ca²⁺ position")
    P.panel(ax, "a", dx=-0.26)

    # ------------------------------------------------------------------ b  BVS vs CMM
    ax = fig.add_subplot(gs[0, 1])
    b = pd.read_csv(P.PROC / "bvs_vs_cmm.tsv", sep="\t")
    x, y = b.valence_cmm.to_numpy(float), b.bvs_local.to_numpy(float)
    r = float(np.corrcoef(x, y)[0, 1])
    off = float(np.mean(y - x))
    lo, hi = min(x.min(), y.min()) - 0.15, max(x.max(), y.max()) + 0.15
    ax.plot([lo, hi], [lo, hi], color=P.C["muted"], lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.scatter(x, y, s=14, color=P.C["ink"], linewidths=0, zorder=3)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("CheckMyMetal valence", labelpad=1.5)
    ax.set_ylabel("Local bond-valence sum")
    ax.annotate(f"r = {r:.3f}, n = {len(b)} sites\nconstant offset +{off:.3f}",
                xy=(0.03, 0.97), xycoords="axes fraction", ha="left", va="top",
                fontsize=5.8, color=P.C["ink"], linespacing=1.3)
    ax.annotate("A constant offset cancels when a mutant\nis subtracted from its wild type",
                xy=(0.97, 0.04), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=5.4, color=P.C["ink2"], linespacing=1.3)
    N.set("b_pearson_r", round(r, 4))
    N.set("b_n_sites", int(len(b)))
    N.set("b_mean_offset", round(off, 4))
    P.title(ax, "Bond-valence sum versus CheckMyMetal")
    P.panel(ax, "b", dx=-0.26)

    # ------------------------------------------------------------------ c  the pilot
    ax = fig.add_subplot(gs[0, 2])
    g = pd.read_csv(P.PROC / "af3_mutant_geometry.tsv", sep="\t")
    g = g[g.status == "present"].copy()
    g["variant"] = g.mutation.map(lambda m: PILOT[m][0] if m in PILOT else None)
    g["role"] = g.mutation.map(lambda m: PILOT[m][1] if m in PILOT else None)
    if g.variant.isna().any():
        raise SystemExit(f"unmapped pilot mutations: {sorted(set(g.loc[g.variant.isna(),'mutation']))}")
    order = (g.groupby("variant").delta_bvs_pct.min().sort_values().index.tolist())
    order = ([v for v in order if not g.loc[g.variant == v, "role"].iloc[0].startswith("control")]
             + [v for v in order if g.loc[g.variant == v, "role"].iloc[0].startswith("control")])
    ypos = {v: i for i, v in enumerate(reversed(order))}
    ax.axvspan(-2, 2, color=P.C["band"], zorder=0)
    ax.axvline(0, color=P.C["ink"], lw=0.4, zorder=1)
    for _, row in g.iterrows():
        yy = ypos[row.variant]
        z = float(row.z_vs_seed_noise)
        beyond = abs(z) >= 2.0
        col = P.C["grey"] if row.role.startswith("control") else P.C["ink"]
        ax.plot([z], [yy], marker="o", ms=4.0, mfc=col if beyond else "white",
                mec=col, mew=0.8, zorder=3)
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([v for v in ypos])
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("Change in valence, in units of the\nwild-type seed spread (z)", labelpad=1.5)
    ax.annotate("shaded, within seed spread\nfilled, beyond seed spread", xy=(0.97, 0.05),
                xycoords="axes fraction", ha="right", va="bottom", fontsize=5.4,
                color=P.C["ink2"], linespacing=1.3)
    n_beyond = int((g.z_vs_seed_noise.abs() >= 2.0).sum())
    N.set("c_n_site_measurements", int(len(g)))
    N.set("c_n_beyond_noise", n_beyond)
    P.title(ax, f"Pilot series, {g.variant.nunique()} variants at {len(g)} sites")
    P.panel(ax, "c", dx=-0.36)
    log.append(f"  c  {n_beyond}/{len(g)} site measurements beyond 2x seed noise")

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("54_figS1_validation", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
