#!/usr/bin/env python3
"""
37_v2fig6_enrichment.py — v2 Figure 6: enrichment, reported as supporting evidence only.

This figure is deliberately last and deliberately caveated. The ClinGen FBN1 VCEP applies PM1 at
calcium-binding and cysteine residues, so ClinVar pathogenic labels are NOT independent of the
features being tested here. The odds ratios describe the labelled data; they are not evidence
that these features cause pathogenicity. That circularity is drawn onto the figure rather than
left to a caption, because a forest plot of 140-fold enrichment is exactly the kind of panel that
gets lifted out of a paper and quoted on its own.

a  odds ratios with 95% Woolf intervals, each feature against both denominators
b  the two denominators side by side, showing why the choice matters

The ClinVar benign set has 29 structurally covered variants and zero at calcium ligands or
cysteines, which is what made every odds ratio in the first analysis infinite. The population
denominator is 1,901 gnomAD observations that are not classified pathogenic -- a denominator,
not a control group.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig6_enrichment"
SCRIPT = Path(__file__).name

SHORT = {
    "side-chain Ca ligand": "Ca²⁺ ligand (side chain)",
    "backbone-only Ca ligand": "Ca²⁺ ligand (backbone only)",
    "motif-overcalled (not a ligand)": "Motif, non-ligand",
    "cysteine position": "cbEGF cysteine",
}
DEN_COLOR = {"population (gnomAD)": F.PALETTE["af3"],
             "ClinVar benign": F.PALETTE["other_cbegf"]}
DEN_SHORT = {"population (gnomAD)": "gnomAD population", "ClinVar benign": "ClinVar benign"}


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    df = pd.read_csv(F.PROC / "enrichment_v2.tsv", sep="\t")
    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: enrichment forest"]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(F.COL_2, 2.5),
                                   gridspec_kw={"width_ratios": [1.65, 1.0]})

    # ---- a: forest -------------------------------------------------------------------
    feats = [f for f in SHORT if f in set(df["feature"])]
    ypos, labels = [], []
    for i, feat in enumerate(feats):
        base = len(feats) - 1 - i
        for j, den in enumerate(("population (gnomAD)", "ClinVar benign")):
            row = df[(df["feature"] == feat) & (df["denominator"] == den)]
            if row.empty:
                continue
            r = row.iloc[0]
            y = base + (0.17 if j == 0 else -0.17)
            orr, lo, hi = float(r["odds_ratio_corrected"]), float(r["ci_lo"]), float(r["ci_hi"])
            axA.plot([lo, hi], [y, y], color=DEN_COLOR[den], lw=0.8, zorder=2,
                     solid_capstyle="butt")
            axA.scatter([orr], [y], s=11, color=DEN_COLOR[den],
                        marker="D" if r["zero_cell"] else "o", linewidths=0, zorder=3)
            for k, v in (("or", orr), ("ci_lo", lo), ("ci_hi", hi),
                         ("zero_cell", bool(r["zero_cell"]))):
                nums.set(f"a_{SHORT[feat]}_{DEN_SHORT[den]}_{k}", v)
            log.append(f"  a  {SHORT[feat]:<28} {DEN_SHORT[den]:<18} "
                       f"OR={orr:>9.3f} [{lo:.3f}, {hi:.3f}] zero_cell={bool(r['zero_cell'])}")
        ypos.append(base); labels.append(SHORT[feat])

    axA.axvline(1.0, color=F.PALETTE["ink"], lw=0.6, zorder=1)
    axA.set_xscale("log")
    # Wide enough for the ClinVar-benign cysteine interval (upper bound 2330). A
    # clipped whisker on a log axis silently understates uncertainty.
    axA.set_xlim(5e-3, 1.2e4)
    axA.set_yticks(ypos); axA.set_yticklabels(labels)
    axA.set_ylim(-0.6, len(feats) - 0.4)
    axA.set_xlabel("odds ratio, pathogenic vs denominator (95% CI)")
    axA.annotate("depleted", xy=(7e-3, len(feats) - 0.52), fontsize=5,
                 color=F.PALETTE["ink2"], ha="left", va="top")
    axA.annotate("enriched", xy=(9e3, len(feats) - 0.52), fontsize=5,
                 color=F.PALETTE["ink2"], ha="right", va="top")
    axA.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=3, color=DEN_COLOR["population (gnomAD)"],
               label="gnomAD population (n = 1,901)"),
        Line2D([], [], marker="o", ls="none", ms=3, color=DEN_COLOR["ClinVar benign"],
               label="ClinVar benign (n = 29)"),
        Line2D([], [], marker="D", ls="none", ms=3, color=F.PALETTE["ink2"],
               label="zero cell, Haldane-corrected"),
    ], loc="lower left", fontsize=5.2, handlelength=1.0)
    F.tidy(axA, "x")
    F.panel_label(axA, "a", dx=-0.34)

    # ---- b: why the denominator matters -----------------------------------------------
    pop = df[df["denominator"] == "population (gnomAD)"]
    ben = df[df["denominator"] == "ClinVar benign"]
    width = 0.34
    xs = np.arange(len(feats))
    for j, (sub, den) in enumerate(((pop, "population (gnomAD)"), (ben, "ClinVar benign"))):
        vals = []
        for feat in feats:
            row = sub[sub["feature"] == feat]
            vals.append(float(row.iloc[0]["pct_denom"]) if not row.empty else 0.0)
        axB.bar(xs + (j - 0.5) * width, vals, width * 0.88, color=DEN_COLOR[den],
                alpha=0.9, linewidth=0, zorder=2,
                label=f"{DEN_SHORT[den]}")
        for x, v in zip(xs + (j - 0.5) * width, vals):
            axB.annotate(f"{v:.1f}", xy=(x, v), xytext=(0, 1.5), textcoords="offset points",
                         ha="center", va="bottom", fontsize=4.8, color=F.PALETTE["ink2"])
        for feat, v in zip(feats, vals):
            nums.set(f"b_{SHORT[feat]}_{DEN_SHORT[den]}_pct", v)
    axB.set_xticks(xs)
    axB.set_xticklabels([SHORT[f].replace(" (", "\n(") for f in feats],
                        rotation=32, ha="right", rotation_mode="anchor", fontsize=5)
    axB.set_ylabel("% of denominator at the feature")
    axB.legend(loc="upper right", fontsize=5.2)
    F.tidy(axB)
    F.panel_label(axB, "b", dx=-0.20)

    # The caveat is part of the figure, not the caption.
    fig.text(0.5, 0.012,
             "Partly circular: ClinGen PM1 assigns pathogenicity at Ca²⁺-binding and cysteine "
             "positions, so the labels are not independent of the tested features.",
             ha="center", va="bottom", fontsize=5.2, style="italic",
             color=F.PALETTE["ink2"])

    fig.subplots_adjust(left=0.175, right=0.995, top=0.965, bottom=0.30, wspace=0.55)
    F.save(fig, NAME, nums)
    F.write_log("37_v2fig6_enrichment", log)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
