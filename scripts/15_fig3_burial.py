#!/usr/bin/env python3
"""
15_fig3_burial.py — PHASE 6, Figure 3. Ruling out the obvious objection.

Figure 2 shows calcium-ligand variants are predicted damaging yet cost almost nothing in fold
stability. The obvious objection is that calcium ligands are simply surface residues, and
surface mutations are always cheap. This figure kills that explanation:

  A  where each class sits    cysteines are buried, calcium ligands only partly so
  B  distance to the ion      the calcium class is defined by geometry, not by position alone
  C  burial-matched ddG       even among BURIED residues the calcium class is ~4x cheaper

Panel C follows the Phase 5 test exactly, including its class definition: `is_ca_consensus`
(annotation-level) against cbEGF residues that are not consensus positions. Figure 2 used the
structure-measured `is_direct_ca_ligand` instead; the two are different populations and are
never mixed within a panel.

The RSA >= 50% stratum in panel C carries no test because it contains ZERO calcium-consensus
variants — no calcium-consensus residue in the four solved constructs is that exposed. That
absence is drawn rather than hidden, since it is itself evidence against the exposure objection.

Writes: figures/fig3_burial_control.{png,pdf}, .numbers.json, logs/15_fig3_*.log
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig3_burial_control"
STATS = fs.RESULTS / "phase5_statistics.tsv"
STRATA = [(0.0, 20.0, "buried\nRSA < 20%"), (20.0, 50.0, "partial\nRSA 20–50%"),
          (50.0, np.inf, "exposed\nRSA ≥ 50%")]
DDG_LIM = (-6, 16)


def stat_row(tbl: pd.DataFrame, test: str) -> pd.Series:
    hit = tbl[tbl.test == test]
    if len(hit) != 1:
        raise SystemExit(f"expected one row for {test!r}, found {len(hit)}")
    return hit.iloc[0]


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(fs.SEED)
    log = [f"Phase 6 — Figure 3 — seed {fs.SEED}"]
    num = fs.Numbers(FIG, Path(__file__).name)

    df = fs.load_merged()
    df["mech"] = fs.mech_group(df)
    tbl = pd.read_csv(STATS, sep="\t")
    cov = df[df.ddG.notna()].copy()

    SHORT = {"ca_ligand": "Ca²⁺\nligand", "cys_removing": "cbEGF Cys\nremoved",
             "other_cbegf": "other\ncbEGF"}

    fig = plt.figure(figsize=(7.4, 5.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], hspace=0.62, wspace=0.28)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, :])

    # ---- A: solvent accessibility --------------------------------------------------
    sA = fs.dist_panel(axA, [(SHORT[k], cov[cov.mech == k].rsa.values, fs.MECH_COLOR[k])
                             for k in fs.MECH_ORDER],
                       "Relative solvent accessibility (%)", rng,
                       annotate_fmt="{:.1f}", ylim=(-4, 118))
    axA.set_yticks([0, 25, 50, 75, 100])
    axA.axhspan(-4, 20, color=fs.PALETTE["grid"], alpha=0.6, lw=0, zorder=0)
    axA.text(2.48, 8, "buried", fontsize=6, color=fs.PALETTE["muted"], ha="left", va="center")
    axA.set_title("Burial differs between the classes", loc="left", fontsize=8.5)
    r_rsa = stat_row(tbl, "RSA: cbEGF cysteines vs other cbEGF residues")
    axA.set_xlabel(f"cbEGF Cys vs other cbEGF  δ = {r_rsa.effect:+.2f}, q = {r_rsa.q_bh:.0e}",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    for k in fs.MECH_ORDER:
        num.set(f"rsa_median_{k}", sA[SHORT[k]]["median"])
        num.set(f"rsa_n_{k}", sA[SHORT[k]]["n"])
        log.append(f"  RSA {k:12} n={sA[SHORT[k]]['n']:3} median={sA[SHORT[k]]['median']:.2f}%")

    # ---- B: distance to the modelled ion -------------------------------------------
    sB = fs.dist_panel(axB, [(SHORT[k], cov[cov.mech == k].dist_to_nearest_ca_A.values,
                              fs.MECH_COLOR[k]) for k in fs.MECH_ORDER],
                       "Distance to nearest Ca²⁺ (Å)", rng, annotate_fmt="{:.1f}",
                       ylim=(0, 33))
    axB.axhline(3.0, color=fs.PALETTE["muted"], lw=0.7, ls="--", zorder=1)
    axB.text(2.48, 3.0, "coordination\nshell", fontsize=6, color=fs.PALETTE["muted"],
             ha="left", va="center")
    axB.set_title("The calcium class is defined by geometry", loc="left", fontsize=8.5)
    axB.set_xlabel("distance measured in the Ca²⁺-placed models\n"
                   "(experimental where available, else AF3-cofolded)",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    for k in fs.MECH_ORDER:
        num.set(f"dist_ca_median_{k}", sB[SHORT[k]]["median"])
        log.append(f"  dCa {k:12} median={sB[SHORT[k]]['median']:.2f} A")

    # ---- C: burial-matched stability cost ------------------------------------------
    ca = cov[cov.is_ca_consensus]
    other = cov[cov.in_cbegf & ~cov.is_ca_consensus]
    width = 0.34
    for si, (lo, hi, lbl) in enumerate(STRATA):
        for j, (name, sub, colour) in enumerate((("Ca-consensus", ca, fs.PALETTE["ca_site"]),
                                                 ("other cbEGF", other,
                                                  fs.PALETTE["other_cbegf"]))):
            vals = sub[(sub.rsa >= lo) & (sub.rsa < hi)].ddG.values
            key = f"panelC_{lbl.split(chr(10))[0]}_{name.replace(' ', '_')}"
            num.set(f"{key}_n", int(len(vals)))
            pos = si + (j - 0.5) * width * 1.15
            if len(vals) == 0:
                axC.text(pos, 0.6, "no\nvariants", ha="center", va="bottom", fontsize=6,
                         color=fs.PALETTE["muted"], style="italic")
                log.append(f"  panel C {lbl.split(chr(10))[0]:8} {name:12} n=0")
                continue
            num.set(f"{key}_median", float(np.median(vals)))
            bp = axC.boxplot([vals], positions=[pos], widths=width, showfliers=False,
                             patch_artist=True, zorder=2)
            for box in bp["boxes"]:
                box.set(facecolor=colour, alpha=0.20, edgecolor=colour, linewidth=1.1)
            for part in ("whiskers", "caps"):
                for art in bp[part]:
                    art.set(color=colour, linewidth=1.0)
            for med in bp["medians"]:
                med.set(color=colour, linewidth=2.0)
            axC.scatter(np.full(len(vals), pos) + rng.uniform(-0.09, 0.09, len(vals)), vals,
                        s=5.5, color=colour, alpha=0.45, linewidths=0.3,
                        edgecolors=fs.PALETTE["surface"], zorder=3)
            axC.annotate(f"n={len(vals)}\n{np.median(vals):.2f}", xy=(pos, DDG_LIM[1]),
                         xytext=(0, -2), textcoords="offset points", ha="center", va="top",
                         fontsize=6.4, color=fs.PALETTE["ink2"])
            n_above = int((vals > DDG_LIM[1]).sum())
            if n_above:
                axC.plot([pos], [DDG_LIM[1]], marker="^", ms=4, color=colour,
                         clip_on=False, zorder=4)
            log.append(f"  panel C {lbl.split(chr(10))[0]:8} {name:12} n={len(vals):3} "
                       f"median={np.median(vals):.2f} ({n_above} off scale)")

    axC.set_xticks(range(len(STRATA)))
    axC.set_xticklabels([s[2] for s in STRATA])
    axC.set_xlim(-0.6, len(STRATA) - 0.4)
    axC.set_ylim(*DDG_LIM)
    axC.axhline(0, color=fs.PALETTE["muted"], lw=0.7, ls=":", zorder=1)
    axC.set_ylabel("FoldX ΔΔG (kcal/mol)")
    axC.yaxis.grid(True, zorder=0)
    axC.set_axisbelow(True)
    axC.set_title("Matched for burial, the calcium class is still the cheaper one",
                  loc="left", fontsize=8.5)
    handles = [Line2D([0], [0], marker="s", ls="none", ms=6, mfc=c, mec=c,
                      label=n) for n, c in (("Calcium-consensus residue", fs.PALETTE["ca_site"]),
                                            ("Other cbEGF residue", fs.PALETTE["other_cbegf"]))]
    # Anchored below the top strip: that strip carries the per-box n/median labels, and the
    # empty upper-right corner (exposed stratum, high ddG) is the only clear space.
    # The legend is not optional here — panel C's green means calcium-CONSENSUS, whereas
    # panels A and B use green for the narrower structure-measured Ca ligand class.
    axC.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.0, 0.90), ncol=1,
               handletextpad=0.3, borderpad=0.2, labelspacing=0.3)

    r_bur = stat_row(tbl, "ddG: Ca-consensus vs other cbEGF, buried (RSA<20%)")
    r_par = stat_row(tbl, "ddG: Ca-consensus vs other cbEGF, partial (RSA 20-50%)")
    axC.set_xlabel(f"buried  δ = {r_bur.effect:+.2f}, q = {r_bur.q_bh:.1e}          "
                   f"partial  δ = {r_par.effect:+.2f}, q = {r_par.q_bh:.2f} (n.s.)          "
                   f"exposed  no calcium-consensus variants exist",
                   fontsize=6.6, color=fs.PALETTE["ink2"], labelpad=6)
    num.set("stat_buried_delta", float(r_bur.effect))
    num.set("stat_partial_delta", float(r_par.effect))

    for ax, letter, dx in ((axA, "A", -0.15), (axB, "B", -0.15), (axC, "C", -0.068)):
        fs.panel_label(ax, letter, dx=dx, dy=1.04)

    fig.suptitle("The calcium result is not a solvent-exposure artifact",
                 fontsize=10, x=0.02, ha="left", y=0.985)

    fs.save(fig, FIG, num)
    plt.close(fig)
    print(f"  log -> {fs.write_log('15_fig3', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
