#!/usr/bin/env python3
"""
14_fig2_mechanism.py — PHASE 6, Figure 2. THE HEADLINE.

The two textbook Marfan missense mechanisms — losing a cbEGF disulfide, and breaking the
calcium consensus — are always described together as "structurally critical". This figure
shows they are not the same kind of damage:

  A  FoldX ddG            cysteine loss is expensive; calcium-ligand loss is not
  B  AlphaMissense        both are called pathogenic
  C  the two together     calcium ligands occupy a quadrant that stability alone cannot reach
  D  FoldX disulfide term the cysteine cost is specifically the lost staple, not general strain

Medians and n come from `data/processed/` and are recomputed here. Effect sizes and BH-corrected
q values are READ from `results/phase5_statistics.tsv` rather than recomputed, so the figure and
the statistics table cannot drift apart; the figure never runs its own hypothesis test.

Wording note carried into `results/phase6_report.md`: the Phase 5 prose calls calcium-ligand ddG
"statistically indistinguishable" from other cbEGF residues. The table it cites does not say
that — it reports a small but significant effect in the OPPOSITE direction to destabilisation
(Cliff's delta -0.23, q = 0.001; ligands are slightly *cheaper* than non-ligands). This figure
annotates the measured effect, not the prose.

Writes: figures/fig2_mechanism_dissociation.{png,pdf}, .numbers.json, logs/14_fig2_*.log
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig2_mechanism_dissociation"
STATS = fs.RESULTS / "phase5_statistics.tsv"

# ddG has a 50 kcal/mol tail; clipping keeps the boxes readable and dist_panel reports how
# many points fall off the top of each axis.
DDG_LIM = (-8, 18)


def stat_row(tbl: pd.DataFrame, test: str) -> pd.Series:
    hit = tbl[tbl.test == test]
    if len(hit) != 1:
        raise SystemExit(f"expected exactly one row for {test!r} in {STATS}, found {len(hit)}")
    return hit.iloc[0]


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(fs.SEED)
    log = [f"Phase 6 — Figure 2 — seed {fs.SEED}"]
    num = fs.Numbers(FIG, Path(__file__).name)

    df = fs.load_merged()
    df["mech"] = fs.mech_group(df)
    tbl = pd.read_csv(STATS, sep="\t")

    cov = df[df.ddG.notna()].copy()
    log.append(f"{len(cov)} variants with FoldX ddG (of {len(df)} total)")
    num.set("n_with_ddg", len(cov))

    # Short axis labels: the full names live in the legend of panel C and in the caption.
    # Long tick labels collide at this panel width and are the commonest way a four-panel
    # figure becomes unreadable in print.
    SHORT = {"ca_ligand": "Ca²⁺\nligand", "cys_removing": "cbEGF Cys\nremoved",
             "other_cbegf": "other\ncbEGF"}
    groups = [(SHORT[k], cov[cov.mech == k], fs.MECH_COLOR[k]) for k in fs.MECH_ORDER]

    fig, axs = plt.subplots(2, 2, figsize=(7.4, 7.0))
    (axA, axB), (axC, axD) = axs
    fig.subplots_adjust(hspace=0.52, wspace=0.30, top=0.90, bottom=0.09)

    # ---- A: folding stability ------------------------------------------------------
    sA = fs.dist_panel(axA, [(lbl, g.ddG.values, c) for lbl, g, c in groups],
                       "FoldX ΔΔG (kcal/mol)", rng, ylim=DDG_LIM)
    axA.axhline(0, color=fs.PALETTE["muted"], lw=0.7, ls=":", zorder=1)
    axA.set_title("Predicted cost to fold stability", loc="left", fontsize=8.5)
    for key, (lbl, _, _) in zip(fs.MECH_ORDER, groups):
        num.set(f"ddG_median_{key}", sA[lbl]["median"])
        num.set(f"ddG_n_{key}", sA[lbl]["n"])
        log.append(f"  ddG {key:12} n={sA[lbl]['n']:3}  median={sA[lbl]['median']:.4f}")

    r_cys = stat_row(tbl, "ddG: cbEGF cysteine-removing vs all other missense")
    r_lig = stat_row(tbl, "ddG: direct Ca ligands vs non-ligands (structure-measured)")
    axA.set_xlabel(f"Cys-removing vs rest  δ = {r_cys.effect:+.2f}, q = {r_cys.q_bh:.0e}\n"
                   f"Ca²⁺ ligand vs non-ligand  δ = {r_lig.effect:+.2f}, "
                   f"q = {r_lig.q_bh:.1e}  (lower, not higher)",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    num.set("stat_ddG_cys_vs_rest_delta", float(r_cys.effect))
    num.set("stat_ddG_ligand_vs_nonligand_delta", float(r_lig.effect))

    # ---- B: predicted pathogenicity ------------------------------------------------
    sB = fs.dist_panel(axB, [(lbl, g.am.values, c) for lbl, g, c in groups],
                       "AlphaMissense pathogenicity", rng, annotate_fmt="{:.3f}",
                       ylim=(0, 1.18))
    axB.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.axhline(0.564, color=fs.PALETTE["muted"], lw=0.7, ls="--", zorder=1)
    axB.text(2.45, 0.564, "AM pathogenic\nthreshold", fontsize=6, va="center",
             ha="left", color=fs.PALETTE["muted"])
    axB.set_title("Predicted pathogenicity (orthogonal to FoldX)", loc="left", fontsize=8.5)
    for key, (lbl, _, _) in zip(fs.MECH_ORDER, groups):
        num.set(f"am_median_{key}", sB[lbl]["median"])
        num.set(f"am_n_{key}", sB[lbl]["n"])
        log.append(f"  AM  {key:12} n={sB[lbl]['n']:3}  median={sB[lbl]['median']:.4f}")
    r_am = stat_row(tbl, "AlphaMissense: direct Ca ligands vs other cbEGF residues")
    axB.set_xlabel(f"Ca²⁺ ligand vs other cbEGF  δ = {r_am.effect:+.2f}, q = {r_am.q_bh:.0e}",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    num.set("stat_am_ligand_vs_other_delta", float(r_am.effect))

    # ---- C: the dissociation, jointly ----------------------------------------------
    # Shade the quadrant the finding lives in: predicted damaging, yet cheap to fold.
    axC.axhspan(0.564, 1.05, xmin=0, xmax=(2.0 - DDG_LIM[0]) / (DDG_LIM[1] - DDG_LIM[0]),
                color=fs.PALETTE["ca_site"], alpha=0.07, lw=0, zorder=0)
    for key in fs.MECH_ORDER:
        g = cov[(cov.mech == key) & cov.am.notna()]
        axC.scatter(g.ddG.clip(*DDG_LIM), g.am, s=11, color=fs.MECH_COLOR[key],
                    alpha=0.55, linewidths=0.4, edgecolors=fs.PALETTE["surface"],
                    label=SHORT[key].replace("\n", " "), zorder=3)
    for key in fs.MECH_ORDER:
        g = cov[(cov.mech == key) & cov.am.notna()]
        axC.plot(np.median(g.ddG), np.median(g.am), marker="D", ms=6,
                 color=fs.MECH_COLOR[key], mec=fs.PALETTE["surface"], mew=1.0, zorder=5)
    axC.axhline(0.564, color=fs.PALETTE["muted"], lw=0.7, ls="--", zorder=1)
    axC.axvline(2.0, color=fs.PALETTE["muted"], lw=0.7, ls=":", zorder=1)
    axC.set_xlim(*DDG_LIM)
    axC.set_ylim(0, 1.05)
    axC.set_xlabel("FoldX ΔΔG (kcal/mol)")
    axC.set_ylabel("AlphaMissense pathogenicity")
    axC.set_title("Damaging but stable: the calcium quadrant", loc="left", fontsize=8.5)
    axC.yaxis.grid(True, zorder=0)
    axC.set_axisbelow(True)
    axC.legend(loc="lower right", markerscale=1.4, handletextpad=0.2, borderpad=0.2,
               labelspacing=0.3, fontsize=6.8)

    rho = stat_row(tbl, "Spearman: FoldX ddG vs AlphaMissense")
    axC.text(0.02, 0.03, f"Spearman ρ = {rho.effect:+.2f}\n"
                         f"[{rho.ci_low:+.2f}, {rho.ci_high:+.2f}], n = {int(rho.n)}",
             transform=axC.transAxes, fontsize=6.8, va="bottom", ha="left",
             color=fs.PALETTE["ink2"])
    num.set("stat_spearman_rho", float(rho.effect))

    n_quadrant = int(((cov.mech == "ca_ligand") & (cov.am > 0.564) & (cov.ddG < 2.0)).sum())
    n_lig_am = int(((cov.mech == "ca_ligand") & cov.am.notna()).sum())
    num.set("ca_ligand_damaging_but_stable", n_quadrant)
    num.set("ca_ligand_with_am", n_lig_am)
    axC.annotate(f"{n_quadrant}/{n_lig_am} Ca²⁺ ligands sit here:\n"
                 f"AM-pathogenic, ΔΔG < 2 kcal/mol",
                 xy=(0.985, 0.70), xycoords="axes fraction", fontsize=6.8,
                 color=fs.PALETTE["ca_site"], ha="right", va="top", fontweight="bold")
    log.append(f"  {n_quadrant}/{n_lig_am} Ca ligands are AM-pathogenic with ddG < 2")

    # ---- D: mechanism specificity --------------------------------------------------
    sD = fs.dist_panel(axD, [(lbl, g.ddG_disulfide.values, c) for lbl, g, c in groups],
                       "FoldX disulfide term (kcal/mol)", rng, ylim=(-0.6, 4.6))
    axD.set_title("The cysteine cost is the lost staple", loc="left", fontsize=8.5)
    for key, (lbl, _, _) in zip(fs.MECH_ORDER, groups):
        num.set(f"disulfide_median_{key}", sD[lbl]["median"])
    r_dis = stat_row(tbl, "FoldX disulfide energy term: Cys-removing vs other")
    axD.set_xlabel(f"Cys-removing vs other  δ = {r_dis.effect:+.2f}, p = {r_dis.p_raw:.0e}",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    num.set("stat_disulfide_delta", float(r_dis.effect))

    for ax, letter in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D")):
        fs.panel_label(ax, letter, dx=-0.14, dy=1.03)

    fig.suptitle("Two Marfan mechanisms, two different kinds of damage",
                 fontsize=10, x=0.02, ha="left", y=0.985)

    fs.save(fig, FIG, num)
    plt.close(fig)
    print(f"  log -> {fs.write_log('14_fig2', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
