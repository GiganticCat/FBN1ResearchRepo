#!/usr/bin/env python3
"""
34_v2fig3_landscape.py — v2 figure 3: where the variants are, and how much we can now measure.

Panel A  every classified variant along the 2871-residue protein, with the domain architecture
Panel B  structural coverage before and after the AlphaFold3 batch (751 -> 3831)
Panel C  what the covered variants are: the four mechanism groups

Panel B is the reason the rest of the v2 analysis exists. The first pass could only measure
variants that happened to fall in one of four experimentally solved fragments -- the domains
crystallographers chose, not a random sample of the gene. Folding all 43 calcium-binding
domains removes that selection.

Writes figures/v2_fig3_landscape.{png,pdf} + numbers sidecar.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig3_landscape"
SCRIPT = Path(__file__).name
PROT_LEN = 2871

PLAIN = {
    "pathogenic": "Pathogenic /\nlikely pathogenic",
    "vus": "Uncertain\nsignificance",
    "benign": "Benign /\nlikely benign",
}


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    df = F.load_v2()
    df["mech"] = F.mech_group_v2(df)
    dom = F.domains()
    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: v2 landscape and coverage"]

    fig = plt.figure(figsize=(F.COL_2, 3.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.0, 1.0],
                          hspace=0.85, wspace=0.42)
    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])

    # ---- A: the landscape -----------------------------------------------------------
    lanes = {"pathogenic": 2.0, "vus": 1.0, "benign": 0.0}
    for cls, y in lanes.items():
        sub = df[df["set_primary"] == cls]
        axA.vlines(sub["position"], y + 0.06, y + 0.72, color=F.PALETTE[cls],
                   lw=0.45, alpha=0.55)
        axA.annotate(f"{PLAIN[cls]}\nn = {len(sub)}", xy=(-40, y + 0.39), ha="right",
                     va="center", fontsize=5.4, color=F.PALETTE[cls], linespacing=1.25)
        nums.set(f"A_n_{cls}", int(len(sub)))
        log.append(f"  A {cls:<12} n={len(sub)}")

    # domain architecture track under the lanes
    track_y, track_h = -0.62, 0.32
    for _, d in dom.iterrows():
        colour = {"cbEGF": F.PALETTE["dom_cbEGF"], "EGF": F.PALETTE["dom_EGF"],
                  "TB": F.PALETTE["dom_TB"]}.get(d["kind"], F.PALETTE["dom_other"])
        axA.add_patch(Rectangle((d["start"], track_y), d["end"] - d["start"], track_h,
                                facecolor=colour, edgecolor="none"))
    axA.add_patch(Rectangle((1, track_y + track_h / 2 - 0.015), PROT_LEN, 0.03,
                            facecolor=F.PALETTE["dom_other"], edgecolor="none", zorder=0))
    axA.annotate("domain\narchitecture", xy=(-40, track_y + track_h / 2), ha="right",
                 va="center", fontsize=5.4, color=F.PALETTE["ink2"], linespacing=1.25)
    axA.annotate("dark blocks, 43 cbEGF domains", xy=(30, track_y - 0.20), ha="left",
                 va="top", fontsize=5.2, color=F.PALETTE["ink2"])

    axA.set_xlim(-780, PROT_LEN + 30)
    axA.set_ylim(-1.05, 2.95)
    axA.set_yticks([])
    axA.set_xticks([1, 500, 1000, 1500, 2000, 2500, PROT_LEN])
    axA.set_xlabel("residue (P35555)")
    # The title states the measured fraction rather than characterising it. An earlier draft
    # said "overwhelmingly inside the calcium-binding domains"; the actual split is 74% of
    # pathogenic versus 62% of benign, which the eye cannot read off this panel anyway, and
    # with 34 benign variants the contrast is not something this figure establishes.
    pat_cb = int(df.loc[df["set_primary"] == "pathogenic", "in_cbegf"].sum())
    n_pat = int((df["set_primary"] == "pathogenic").sum())
    nums.set("A_pathogenic_in_cbegf", pat_cb)
    nums.set("A_pathogenic_in_cbegf_pct", round(100 * pat_cb / n_pat, 1))
    F.panel_label(axA, "a", dx=-0.075)
    for s in ("left", "right", "top"):
        axA.spines[s].set_visible(False)

    # ---- B: coverage before/after ---------------------------------------------------
    before = int(df["foldx_ddG"].notna().sum())
    after = int(df["has_structure"].sum())
    total = int(len(df))
    bars = [("experimental\ntemplates (4)", before, F.PALETTE["other_cbegf"]),
            ("AF3, all 43\ncbEGF domains", after, F.PALETTE["af3"])]
    for i, (lbl, v, c) in enumerate(bars):
        axB.bar(i, v, 0.55, color=c, alpha=0.9, linewidth=0)
        axB.annotate(f"{v:,}\n{100*v/total:.0f}%", xy=(i, v), xytext=(0, 2),
                     textcoords="offset points", ha="center", va="bottom",
                     fontsize=5.4, color=F.PALETTE["ink2"], linespacing=1.2)
    axB.axhline(total, color=F.PALETTE["ink2"], lw=0.9, ls=(0, (4, 3)))
    axB.annotate(f"all {total:,}", xy=(-0.45, total), xytext=(0, 2),
                 textcoords="offset points", ha="left", va="bottom", fontsize=5.2,
                 color=F.PALETTE["ink2"])
    axB.set_xticks([0, 1]); axB.set_xticklabels([b[0] for b in bars])
    axB.set_ylabel("variants with structural metrics")
    axB.set_ylim(0, total * 1.30)
    F.tidy(axB)
    F.panel_label(axB, "b")
    nums.set("B_covered_before", before)
    nums.set("B_covered_after", after)
    nums.set("B_total_variants", total)
    log.append(f"  B coverage {before} -> {after} of {total}")

    # ---- C: what the covered variants are -------------------------------------------
    counts = [(k, int((df["mech"] == k).sum())) for k in F.MECH2_ORDER]
    ys = np.arange(len(counts))[::-1]
    for y, (k, n) in zip(ys, counts):
        axC.barh(y, n, 0.6, linewidth=0, color=F.MECH2_COLOR[k], alpha=0.9)
        axC.annotate(f"{n:,}", xy=(n, y), xytext=(2.5, 0), textcoords="offset points",
                     va="center", fontsize=5.4, color=F.PALETTE["ink2"])
        nums.set(f"C_n_{k}", n)
        log.append(f"  C {k:<22} n={n}")
    axC.set_yticks(ys)
    axC.set_yticklabels([F.MECH2_LABEL[k].replace("\n", " ") for k, _ in counts])
    axC.set_xlabel("variants")
    axC.set_xlim(0, max(n for _, n in counts) * 1.18)
    F.tidy(axC, "x")
    F.panel_label(axC, "c", dx=-0.52)

    F.save(fig, NAME, nums)
    F.write_log("34_v2fig3_landscape", log)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
