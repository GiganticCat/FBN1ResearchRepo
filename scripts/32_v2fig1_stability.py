#!/usr/bin/env python3
"""
32_v2fig1_stability.py — v2 Figure 2: folding stability by site class, in two force fields.

The claim under test: variants at Ca2+-coordinating positions do not destabilise the cbEGF fold,
while cysteine-removing variants do. That claim previously rested on FoldX alone, and FoldX is
measurably blind to the metal (ion ablation shifts ddG by a median of +0.15 kcal/mol at direct
ligands). Rosetta forms explicit bonds between the ion and its ligands and is the independent
check.

a  FoldX ddG by site class
b  Rosetta ddG by site class
c  Cliff's delta against the composition-matched control, both force fields, both mechanisms
d  AlphaMissense, a sequence-only predictor that never saw a structure

Double-column width. Panel letters are lowercase and no panel carries a title -- the reading
belongs in the caption (results/v2_figure_captions.md).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig1_stability"
SCRIPT = Path(__file__).name

# Rosetta's tail is steric, not thermodynamic: 21.5% of values exceed 100 REU and they are
# almost entirely bulky substitutions (Y/F/R/W) at cysteines, where a fixed backbone and an 8 A
# repack shell cannot absorb the new side chain. Medians are unaffected -- the calcium groups
# contain no such value at all (max 66.9) -- but the axis is clipped so the boxes stay readable,
# and dist_panel marks how many points sit above it.
ROS_CLIP = 60.0
N_BOOT = 2000          # matches 28_statistics_v2.py


def cliffs_delta_ci(a: np.ndarray, b: np.ndarray, rng, n_boot: int = N_BOOT):
    """Cliff's delta of a vs b, with a percentile bootstrap CI. Same estimator as Phase 5 v2."""
    def delta(x, y):
        # P(x > y) - P(x < y), computed by rank rather than by the O(n*m) pair loop
        allv = np.concatenate([x, y])
        order = allv.argsort(kind="mergesort")
        ranks = np.empty(len(allv), dtype=float)
        ranks[order] = np.arange(1, len(allv) + 1)
        # average ranks over ties so exact equalities contribute zero, not a spurious direction
        _, inv, cnt = np.unique(allv, return_inverse=True, return_counts=True)
        sums = np.zeros(len(cnt)); np.add.at(sums, inv, ranks)
        ranks = (sums / cnt)[inv]
        rx = ranks[:len(x)].sum()
        u = rx - len(x) * (len(x) + 1) / 2.0
        return 2.0 * u / (len(x) * len(y)) - 1.0

    d = float(delta(a, b))
    boots = np.empty(n_boot)
    for i in range(n_boot):
        boots[i] = delta(rng.choice(a, len(a), replace=True),
                         rng.choice(b, len(b), replace=True))
    return d, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(F.SEED)
    df = F.load_v2()
    df["mech"] = F.mech_group_v2(df)
    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: v2 stability figure", f"  seed {F.SEED}"]

    fig, axes = plt.subplots(1, 4, figsize=(F.COL_2, 2.15))
    axA, axB, axC, axD = axes

    keys = F.MECH2_ORDER

    def groups(col):
        return [(F.MECH2_LABEL[k],
                 df.loc[df["mech"] == k, col].dropna().to_numpy(dtype=float),
                 F.MECH2_COLOR[k]) for k in keys]

    # ---- a: FoldX -------------------------------------------------------------------
    sA = F.dist_panel(axA, groups("foldx_ddG"), "FoldX ΔΔG (kcal mol⁻¹)", rng,
                      annotate_fmt="{:.2f}", ylim=(-3, 22))
    F.panel_label(axA, "a")

    # ---- b: Rosetta -----------------------------------------------------------------
    sB = F.dist_panel(axB, groups("rosetta_ddG"), "Rosetta ΔΔG (REU)", rng,
                      annotate_fmt="{:.2f}", ylim=(-12, ROS_CLIP))
    F.panel_label(axB, "b")

    for ax in (axA, axB):
        ax.axhline(0, color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), zorder=1)

    # ---- c: composition-matched effect sizes ----------------------------------------
    # NOT a ratio of medians: FoldX's control median is 0.24 kcal/mol, and dividing by a
    # near-zero denominator renders a 0.8 kcal/mol difference as "4.25x", an artifact of the
    # denominator. NOT two force fields on a shared energy axis either -- kcal/mol and REU are
    # different quantities. Cliff's delta is unit-free and carries an interval, so "no
    # difference" is shown rather than inferred from an unlabelled gap.
    combos = [("ca_ligand", "foldx_ddG"), ("ca_ligand", "rosetta_ddG"),
              ("cys_removing", "foldx_ddG"), ("cys_removing", "rosetta_ddG")]
    xs = [0, 0.5, 1.5, 2.0]
    for x, (mech, col) in zip(xs, combos):
        tgt = df.loc[df["mech"] == mech, col].dropna().to_numpy(dtype=float)
        ctl = df.loc[df["mech"] == "composition_control", col].dropna().to_numpy(dtype=float)
        d, lo, hi = cliffs_delta_ci(tgt, ctl, rng)
        crosses = lo <= 0.0 <= hi
        colour = F.MECH2_COLOR[mech]
        filled = col == "rosetta_ddG"
        axC.errorbar([x], [d], yerr=[[d - lo], [hi - d]], fmt="o" if filled else "s",
                     ms=3.4, mfc=colour if filled else F.PALETTE["surface"], mec=colour,
                     color=colour, ecolor=colour, elinewidth=0.8, capsize=1.8,
                     capthick=0.6, mew=0.8, zorder=3)
        for key, val in (("delta", d), ("ci_lo", lo), ("ci_hi", hi), ("n", len(tgt)),
                         ("n_control", len(ctl)), ("ci_crosses_zero", crosses)):
            nums.set(f"c_{mech}_{col}_{key}", val)
        log.append(f"  c  {mech:<14} {col:<12} d={d:+.3f} [{lo:+.3f}, {hi:+.3f}] "
                   f"n={len(tgt)} vs {len(ctl)}  crosses_zero={crosses}")

    axC.axhline(0.0, color=F.PALETTE["ink"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    axC.set_xticks([0.25, 1.75])
    axC.set_xticklabels(["Ca²⁺ ligand", "cbEGF Cys loss"], rotation=32, ha="right",
                        rotation_mode="anchor")
    axC.set_xlim(-0.4, 2.4)
    axC.set_ylim(-0.2, 0.95)
    axC.set_ylabel("Cliff's δ vs D/E/N control")
    axC.legend(handles=[
        Line2D([], [], marker="s", ls="none", ms=3.2, mfc=F.PALETTE["surface"],
               mec=F.PALETTE["ink2"], mew=0.8, label="FoldX"),
        Line2D([], [], marker="o", ls="none", ms=3.2, color=F.PALETTE["ink2"], label="Rosetta"),
    ], loc="upper left", fontsize=5.5, handlelength=1.0)
    F.tidy(axC)
    F.panel_label(axC, "c")

    # ---- d: AlphaMissense -----------------------------------------------------------
    sD = F.dist_panel(axD, groups("am"), "AlphaMissense score", rng,
                      annotate_fmt="{:.3f}", ylim=(-0.05, 1.15))
    F.panel_label(axD, "d")

    for k in keys:
        for tag, s in (("a_foldx", sA), ("b_rosetta", sB), ("d_am", sD)):
            st = s[F.MECH2_LABEL[k]]
            nums.set(f"{tag}_{k}_median", st["median"])
            nums.set(f"{tag}_{k}_n", st["n"])
            log.append(f"  {tag:<10} {k:<22} n={st['n']:<5} median={st['median']:.4f}")

    for ax in axes:
        ax.tick_params(axis="x", pad=1.5)
    fig.subplots_adjust(left=0.075, right=0.995, top=0.93, bottom=0.34, wspace=0.55)
    F.save(fig, NAME, nums)
    F.write_log("32_v2fig1_stability", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
