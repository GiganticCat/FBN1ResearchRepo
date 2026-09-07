#!/usr/bin/env python3
"""
35_v2fig4_burial.py — v2 Figure 3: the Ca2+ result is not a solvent-exposure artifact.

The obvious objection to Figure 2 is confounding by burial. Ca2+ ligands are more buried than an
average residue (median RSA 21.7% vs 38.4%), and buried positions are generally more costly to
mutate. If the calcium groups looked cheap simply because of where they sit, that would be an
artifact rather than a finding.

a  relative solvent accessibility by site class -- establishes that the confound is real
b  Rosetta ddG against RSA, Ca2+ ligands over the D/E/N control cloud
c  ddG stratified by burial: within each stratum, Ca2+ ligands vs the composition control

Panel c is the actual control. If the calcium effect were an exposure artifact, the difference
would vanish once burial is held fixed. It does not appear in any stratum, which is the point.

Note on strata: `(rsa or -1)` style guards are avoided deliberately. In Phase 6 that idiom
dropped every residue with RSA exactly 0.0 -- falsy in Python -- and those are the most buried,
most destabilising positions in the set.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig4_burial"
SCRIPT = Path(__file__).name

MIN_N_BAR = 10   # below this, show the raw points; a median+IQR from n<10 overstates precision

STRATA = [("buried\n<20%", 0.0, 20.0), ("partial\n20–50%", 20.0, 50.0),
          ("exposed\n≥50%", 50.0, 1e9)]


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(F.SEED)
    df = F.load_v2()
    df["mech"] = F.mech_group_v2(df)
    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: burial control", f"  seed {F.SEED}"]

    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(F.COL_2, 2.3),
                                        gridspec_kw={"width_ratios": [1.0, 1.05, 1.25]})

    keys = F.MECH2_ORDER

    # ---- a: RSA by class ------------------------------------------------------------
    gA = [(F.MECH2_LABEL[k], df.loc[df["mech"] == k, "rsa"].dropna().to_numpy(float),
           F.MECH2_COLOR[k]) for k in keys]
    sA = F.dist_panel(axA, gA, "relative solvent accessibility (%)", rng,
                      annotate_fmt="{:.1f}", ylim=(-4, 105))
    for k in keys:
        st = sA[F.MECH2_LABEL[k]]
        nums.set(f"a_rsa_{k}_median", st["median"]); nums.set(f"a_rsa_{k}_n", st["n"])
        log.append(f"  a  RSA {k:<22} n={st['n']:<5} median={st['median']:.2f}")
    F.panel_label(axA, "a")

    # ---- b: ddG vs RSA --------------------------------------------------------------
    ctl = df[(df["mech"] == "composition_control") & df["rosetta_ddG"].notna() & df["rsa"].notna()]
    lig = df[(df["mech"] == "ca_ligand") & df["rosetta_ddG"].notna() & df["rsa"].notna()]
    axB.scatter(ctl["rsa"], ctl["rosetta_ddG"], s=2.2, color=F.PALETTE["other_cbegf"],
                alpha=0.45, linewidths=0, zorder=2, rasterized=True)
    axB.scatter(lig["rsa"], lig["rosetta_ddG"], s=3.4, color=F.PALETTE["ca_site"],
                alpha=0.8, linewidths=0, zorder=3, rasterized=True)
    axB.set_xlabel("relative solvent accessibility (%)")
    axB.set_ylabel("Rosetta ΔΔG (REU)")
    axB.set_xlim(-3, 100); axB.set_ylim(-12, 45)
    axB.axhline(0, color=F.PALETTE["ink2"], lw=0.5, ls=(0, (3, 2)), zorder=1)
    axB.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=2.6, color=F.PALETTE["ca_site"],
               label=f"Ca²⁺ ligand ({len(lig)})"),
        Line2D([], [], marker="o", ls="none", ms=2.6, color=F.PALETTE["other_cbegf"],
               label=f"D/E/N control ({len(ctl)})"),
    ], loc="upper right", fontsize=5.5, handlelength=1.0)
    F.tidy(axB)
    F.panel_label(axB, "b")
    nums.set("b_n_ca_ligand", int(len(lig))); nums.set("b_n_control", int(len(ctl)))

    # ---- c: stratified by burial ----------------------------------------------------
    width = 0.34
    for si, (slabel, lo, hi) in enumerate(STRATA):
        for j, k in enumerate(("ca_ligand", "composition_control")):
            sub = df[(df["mech"] == k) & df["rsa"].notna() & df["rosetta_ddG"].notna()]
            sel = sub[(sub["rsa"] >= lo) & (sub["rsa"] < hi)]["rosetta_ddG"].to_numpy(float)
            x = si + (j - 0.5) * width
            nums.set(f"c_{slabel.split(chr(10))[0]}_{k}_n", int(len(sel)))
            if len(sel) == 0:
                # An empty stratum is evidence, not a gap: no Ca ligand in this study is
                # exposed. Marking it beats leaving a blank the reader must interpret.
                axC.annotate("none", xy=(x, 0.6), fontsize=5, rotation=90, ha="center",
                             va="bottom", color=F.PALETTE["ink2"])
                nums.set(f"c_{slabel.split(chr(10))[0]}_{k}_median", None)
                continue
            med = float(np.median(sel))
            if len(sel) < MIN_N_BAR:
                # A bar with an interquartile range drawn from one observation is a lie about
                # precision: the single exposed Ca ligand in the whole study (ΔΔG 16.2 REU) set
                # the y-axis and made the exposed stratum look like a finding. Below the
                # threshold the raw points are drawn instead, with no summary statistic.
                axC.scatter(np.full(len(sel), x), sel, s=6, color=F.MECH2_COLOR[k],
                            linewidths=0, zorder=3)
                axC.annotate(f"n={len(sel)}", xy=(x, max(sel)), xytext=(0, 3),
                             textcoords="offset points", ha="center", va="bottom",
                             fontsize=5, color=F.PALETTE["ink2"])
                nums.set(f"c_{slabel.split(chr(10))[0]}_{k}_median", med)
                log.append(f"  c  {slabel.replace(chr(10),' '):<16} {k:<22} "
                           f"n={len(sel):<5} median={med:.3f}  [below bar threshold]")
                continue
            q1, q3 = np.percentile(sel, [25, 75])
            axC.bar(x, med, width * 0.88, color=F.MECH2_COLOR[k], alpha=0.85,
                    zorder=2, linewidth=0)
            axC.errorbar([x], [med], yerr=[[med - q1], [q3 - med]], fmt="none",
                         ecolor=F.PALETTE["ink"], elinewidth=0.6, capsize=1.6,
                         capthick=0.5, zorder=3)
            axC.annotate(f"{len(sel)}", xy=(x, 0), xytext=(0, -1.5),
                         textcoords="offset points", ha="center", va="top", fontsize=5,
                         color=F.PALETTE["ink2"])
            nums.set(f"c_{slabel.split(chr(10))[0]}_{k}_median", med)
            log.append(f"  c  {slabel.replace(chr(10),' '):<16} {k:<22} "
                       f"n={len(sel):<5} median={med:.3f}")
    axC.set_xticks(range(len(STRATA)))
    axC.set_xticklabels([s[0] for s in STRATA])
    axC.set_ylabel("Rosetta ΔΔG (REU), median ± IQR")
    axC.set_ylim(0, 9.5)
    axC.set_xlim(-0.55, len(STRATA) - 0.45)
    axC.legend(handles=[
        Line2D([], [], marker="s", ls="none", ms=3.4, color=F.PALETTE["ca_site"],
               label="Ca²⁺ ligand"),
        Line2D([], [], marker="s", ls="none", ms=3.4, color=F.PALETTE["other_cbegf"],
               label="D/E/N control"),
    ], loc="upper right", fontsize=5.5, handlelength=1.0)
    F.tidy(axC)
    F.panel_label(axC, "c")

    fig.subplots_adjust(left=0.075, right=0.995, top=0.955, bottom=0.30, wspace=0.42)
    F.save(fig, NAME, nums)
    F.write_log("35_v2fig4_burial", log)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
