#!/usr/bin/env python3
"""
36_v2fig5_validation.py — v2 Figure 5: is the modelled calcium site trustworthy?

Every geometric claim in this study rests on two things that must be shown rather than asserted:
that AlphaFold3 places Ca2+ where the experimental structures put it, and that the local
bond-valence implementation reproduces what CheckMyMetal returns. Neither is assumed anywhere
else in the paper.

a  AF3 vs experimental Ca2+ position, per domain, against the backbone agreement of the same
   domain -- the ion is placed at least as accurately as the fold it sits in
b  local BVS vs CheckMyMetal VALENCE over the 16 server-scored sites
c  coordination number, local vs CheckMyMetal, on the same sites

Panel b is why no second CheckMyMetal hand-off was needed for the mutants: the constant offset
is a documented difference in donor-set and cutoff conventions, and it cancels in a wild-type
minus mutant comparison. The correlation, not the absolute value, is what licenses that.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig5_validation"
SCRIPT = Path(__file__).name


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    cal = pd.read_csv(F.RESULTS / "phase4_calibration.tsv", sep="\t")
    bvs = pd.read_csv(F.PROC / "bvs_vs_cmm.tsv", sep="\t")
    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: placement and bond-valence validation"]

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(F.COL_2, 2.25),
                                        gridspec_kw={"width_ratios": [1.15, 1.0, 0.85]})

    # ---- a: Ca placement vs backbone agreement ---------------------------------------
    # Model 0 only: the five AF3 models of a construct are not independent measurements of the
    # same thing, and averaging them here would hide the seed spread reported in Figure 4.
    c0 = cal[cal["model"] == 0].copy()
    axA.scatter(c0["domain_backbone_rmsd_A"], c0["ca_deviation_A"], s=9,
                color=F.PALETTE["af3"], linewidths=0, zorder=3)
    # Limits from the data, not from a guess: the NMR references (1LMJ, 1EMN) reach 2.5 A
    # backbone r.m.s.d., and an axis capped at 1.05 silently dropped half the domains.
    lim = float(max(c0["domain_backbone_rmsd_A"].max(), c0["ca_deviation_A"].max())) * 1.10
    axA.plot([0, lim], [0, lim], color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    axA.annotate("y = x", xy=(lim * 0.82, lim * 0.86), fontsize=5, color=F.PALETTE["ink2"],
                 rotation=38, ha="center", va="center")
    for _, r in c0.iterrows():
        axA.annotate(r["domain"].replace("cbEGF", ""), xy=(r["domain_backbone_rmsd_A"],
                     r["ca_deviation_A"]), xytext=(2.2, 1.4), textcoords="offset points",
                     fontsize=4.6, color=F.PALETTE["ink2"])
    axA.set_xlabel("domain backbone r.m.s.d. (Å)")
    axA.set_ylabel("Ca²⁺ position deviation (Å)")
    axA.set_xlim(0, lim); axA.set_ylim(0, lim)
    axA.set_aspect("equal", adjustable="box")
    F.tidy(axA)
    F.panel_label(axA, "a")
    med_ca = float(c0["ca_deviation_A"].median())
    med_bb = float(c0["domain_backbone_rmsd_A"].median())
    below = int((c0["ca_deviation_A"] < c0["domain_backbone_rmsd_A"]).sum())
    nums.set("a_median_ca_deviation_A", med_ca)
    nums.set("a_median_backbone_rmsd_A", med_bb)
    nums.set("a_n_domains", int(len(c0)))
    nums.set("a_n_ca_better_than_backbone", below)
    log.append(f"  a  n={len(c0)} domains; median Ca deviation {med_ca:.3f} A, "
               f"median backbone rmsd {med_bb:.3f} A; {below}/{len(c0)} below the diagonal")

    # ---- b: local BVS vs CheckMyMetal -------------------------------------------------
    x, y = bvs["valence_cmm"].to_numpy(float), bvs["bvs_local"].to_numpy(float)
    r = float(np.corrcoef(x, y)[0, 1])
    offset = float(np.mean(y - x))
    axB.scatter(x, y, s=9, color=F.PALETTE["experimental"], linewidths=0, zorder=3)
    lo, hi = 1.15, 2.15
    axB.plot([lo, hi], [lo, hi], color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    axB.plot([lo, hi], [lo + offset, hi + offset], color=F.PALETTE["transplant"],
             lw=0.8, zorder=2)
    axB.set_xlabel("CheckMyMetal VALENCE")
    axB.set_ylabel("local bond-valence sum")
    axB.set_xlim(lo, hi); axB.set_ylim(lo, hi)
    axB.annotate(f"r = {r:.3f}\noffset {offset:+.3f}", xy=(0.04, 0.96),
                 xycoords="axes fraction", fontsize=5.5, va="top",
                 color=F.PALETTE["ink2"], linespacing=1.3)
    axB.legend(handles=[
        Line2D([], [], color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), label="y = x"),
        Line2D([], [], color=F.PALETTE["transplant"], lw=0.8, label="constant offset"),
    ], loc="lower right", fontsize=5.5, handlelength=1.4)
    F.tidy(axB)
    F.panel_label(axB, "b")
    nums.set("b_pearson_r", round(r, 4))
    nums.set("b_mean_offset", round(offset, 4))
    nums.set("b_n_sites", int(len(bvs)))
    log.append(f"  b  n={len(bvs)} sites; Pearson r={r:.4f}; mean offset {offset:+.4f}")

    # ---- c: coordination number agreement ---------------------------------------------
    # Jittered because several sites share the same integer pair and would overplot to a single
    # dot, which would misrepresent how many observations support each cell.
    rng = np.random.default_rng(F.SEED)
    jx = bvs["cn_cmm"].to_numpy(float) + rng.uniform(-0.11, 0.11, len(bvs))
    jy = bvs["cn_local"].to_numpy(float) + rng.uniform(-0.11, 0.11, len(bvs))
    axC.scatter(jx, jy, s=9, color=F.PALETTE["experimental"], linewidths=0, zorder=3)
    axC.plot([4.5, 7.5], [4.5, 7.5], color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    axC.set_xlabel("CheckMyMetal CN")
    axC.set_ylabel("local CN")
    axC.set_xlim(4.5, 7.5); axC.set_ylim(4.5, 7.5)
    axC.set_xticks([5, 6, 7]); axC.set_yticks([5, 6, 7])
    agree = int((bvs["cn_local"] == bvs["cn_cmm"]).sum())
    axC.annotate(f"{agree}/{len(bvs)} exact", xy=(0.04, 0.96), xycoords="axes fraction",
                 fontsize=5.5, va="top", color=F.PALETTE["ink2"])
    F.tidy(axC)
    F.panel_label(axC, "c")
    nums.set("c_n_exact_agreement", agree)
    log.append(f"  c  {agree}/{len(bvs)} sites agree exactly on coordination number")

    fig.subplots_adjust(left=0.075, right=0.995, top=0.955, bottom=0.19, wspace=0.42)
    F.save(fig, NAME, nums)
    F.write_log("36_v2fig5_validation", log)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
