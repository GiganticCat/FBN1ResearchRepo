#!/usr/bin/env python3
"""
33_v2fig2_geometry.py — v2 Figure 4: Ca2+ coordination in cofolded mutant models.

If Ca2+-ligand variants do not destabilise the fold (Figure 2), the damage must be elsewhere.
This measures the ion site directly in AlphaFold3 models of each mutant against the wild type of
the same construct.

The quantity is the bond-valence sum (BVS) at the ion. The local implementation tracks
CheckMyMetal's VALENCE at r = 0.961 over the 16 server-scored sites (Figure 5b), so wild-type to
mutant differences are trustworthy without a further manual hand-off.

a  change in BVS per variant, both sites
b  the same change in units of the wild-type seed spread (5 folds per construct)
c  coordination number, wild type -> mutant

C1138F (disulfide loss, no Ca2+ ligand) and P1141L (benign) are controls: neither should perturb
the site, and they set the floor for believing the rest.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as F  # noqa: E402

NAME = "v2_fig2_geometry"
SCRIPT = Path(__file__).name

ROLE = {
    "N2144S": ("ligand", "Kd 5–9× ↑ (Kettle 1999)"),
    "N2183S": ("ligand", ""), "E2130K": ("ligand", ""), "D1113G": ("ligand", ""),
    "E1073K": ("ligand", ""), "D1487G": ("ligand", ""), "E913K": ("ligand", ""),
    "C1138F": ("control", "S–S loss"), "P1141L": ("control", "benign"),
}
UNIPROT = {"N18S": "N2144S", "N57S": "N2183S", "E4K": "E2130K", "D45G": "D1113G",
           "E5K": "E1073K", "D2G": "D1487G", "E107K": "E913K", "C70F": "C1138F",
           "P73L": "P1141L"}
NOISE_Z = 2.0     # a shift below twice the wild-type seed spread is not a finding


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    df = pd.read_csv(F.PROC / "af3_mutant_geometry.tsv", sep="\t")
    df = df[df["status"] == "present"].copy()
    df["variant"] = df["mutation"].map(UNIPROT)
    if df["variant"].isna().any():
        raise SystemExit(f"unmapped mutations: {sorted(df.loc[df['variant'].isna(),'mutation'])}")
    df["role"] = df["variant"].map(lambda v: ROLE[v][0])

    worst = df.groupby("variant")["delta_bvs_pct"].min().sort_values()
    order = ([v for v in worst.index if ROLE[v][0] == "ligand"]
             + [v for v in worst.index if ROLE[v][0] == "control"])
    ypos = {v: i for i, v in enumerate(reversed(order))}

    nums = F.Numbers(NAME, SCRIPT)
    log = [f"{SCRIPT}: Ca2+ site geometry", f"  {len(df)} site rows, {len(order)} variants"]

    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=(F.COL_2, 2.35), gridspec_kw={"width_ratios": [1.22, 1.0, 0.72]})

    def colour(v):
        return F.PALETTE["ca_site"] if ROLE[v][0] == "ligand" else F.PALETTE["other_cbegf"]

    def yof(r):
        return ypos[r["variant"]] + (0.17 if r["site"] == 0 else -0.17)

    # ---- a: BVS change --------------------------------------------------------------
    for _, r in df.iterrows():
        y = yof(r)
        axA.plot([0, r["delta_bvs_pct"]], [y, y], color=colour(r["variant"]),
                 lw=0.8, alpha=0.6, zorder=2, solid_capstyle="butt")
        axA.scatter([r["delta_bvs_pct"]], [y], s=8, color=colour(r["variant"]),
                    linewidths=0, zorder=3)
        nums.set(f"a_{r['variant']}_site{int(r['site'])}_delta_bvs_pct",
                 float(r["delta_bvs_pct"]))
    axA.axvline(0, color=F.PALETTE["ink"], lw=0.5)
    axA.set_yticks(range(len(order)))
    axA.set_yticklabels(list(reversed(order)))
    axA.set_xlabel("ΔBVS (%)")
    axA.set_xlim(-33, 26)
    F.tidy(axA, "x")
    F.panel_label(axA, "a", dx=-0.30)
    for v in order:
        if ROLE[v][1]:
            axA.annotate(ROLE[v][1], xy=(25, ypos[v]), fontsize=4.8, style="italic",
                         color=F.PALETTE["ink2"], ha="right", va="center")

    # ---- b: against seed noise ------------------------------------------------------
    axB.axvspan(-NOISE_Z, NOISE_Z, color=F.PALETTE["grid"], zorder=0)
    for _, r in df.iterrows():
        beyond = abs(r["z_vs_seed_noise"]) > NOISE_Z
        axB.scatter([r["z_vs_seed_noise"]], [yof(r)], s=9 if beyond else 7,
                    color=colour(r["variant"]) if beyond else F.PALETTE["surface"],
                    edgecolors=colour(r["variant"]), linewidths=0.6, zorder=3)
        nums.set(f"b_{r['variant']}_site{int(r['site'])}_z", float(r["z_vs_seed_noise"]))
    axB.axvline(0, color=F.PALETTE["ink"], lw=0.5)
    axB.set_yticks(range(len(order))); axB.set_yticklabels([])
    axB.set_xlabel("ΔBVS / wild-type seed s.d.")
    axB.set_xlim(-7.5, 7.5)
    F.tidy(axB, "x")
    F.panel_label(axB, "b", dx=-0.06)
    axB.annotate("seed noise", xy=(0, len(order) - 0.4), fontsize=5,
                 color=F.PALETTE["ink2"], ha="center", va="center")

    # ---- c: coordination number -----------------------------------------------------
    for _, r in df.iterrows():
        y = yof(r)
        if abs(r["cn_mut"] - r["cn_wt"]) < 1e-9:
            # Unchanged is a result: the ion kept every contact. A zero-length arrow renders
            # nothing, which would read as missing data.
            axC.scatter([r["cn_wt"]], [y], s=7, facecolors=F.PALETTE["surface"],
                        edgecolors=colour(r["variant"]), linewidths=0.6, zorder=3)
        else:
            axC.annotate("", xy=(r["cn_mut"], y), xytext=(r["cn_wt"], y),
                         arrowprops=dict(arrowstyle="-|>,head_width=0.12,head_length=0.28",
                                         color=colour(r["variant"]), lw=0.7,
                                         shrinkA=0, shrinkB=0))
        nums.set(f"c_{r['variant']}_site{int(r['site'])}_cn_wt", float(r["cn_wt"]))
        nums.set(f"c_{r['variant']}_site{int(r['site'])}_cn_mut", float(r["cn_mut"]))
    axC.set_yticks(range(len(order))); axC.set_yticklabels([])
    axC.set_xlabel("coordination no.")
    axC.set_xlim(3.5, 8.5); axC.set_xticks([4, 5, 6, 7, 8])
    F.tidy(axC, "x")
    F.panel_label(axC, "c", dx=-0.09)

    for ax in (axA, axB, axC):
        ax.set_ylim(-0.65, len(order) - 0.35)
        ax.spines["left"].set_visible(ax is axA)
        ax.tick_params(axis="y", length=0 if ax is not axA else 2.2)

    fig.legend(handles=[
        Line2D([], [], marker="o", ls="none", ms=3, color=F.PALETTE["ca_site"],
               label="Ca²⁺ ligand variant"),
        Line2D([], [], marker="o", ls="none", ms=3, color=F.PALETTE["other_cbegf"],
               label="control"),
        Line2D([], [], marker="o", ls="none", ms=3, mfc="white",
               mec=F.PALETTE["ink2"], mew=0.6, label="within seed noise / unchanged"),
    ], loc="lower center", ncol=3, fontsize=5.5, handlelength=1.0,
        bbox_to_anchor=(0.55, -0.035))

    n_beyond = int((df["z_vs_seed_noise"].abs() > NOISE_Z).sum())
    lig_beyond = int(((df["role"] == "ligand") &
                      (df["z_vs_seed_noise"].abs() > NOISE_Z)).sum())
    nums.set("n_site_measurements", int(len(df)))
    nums.set("n_beyond_noise", n_beyond)
    nums.set("n_ligand_sites_beyond_noise", lig_beyond)
    log.append(f"  {n_beyond}/{len(df)} site measurements exceed ±{NOISE_Z} seed s.d. "
               f"({lig_beyond} at Ca ligand variants)")

    fig.subplots_adjust(left=0.115, right=0.99, top=0.965, bottom=0.235, wspace=0.10)
    F.save(fig, NAME, nums)
    F.write_log("33_v2fig2_geometry", log)
    print("\n".join(log))
    return 0


if __name__ == "__main__":
    sys.exit(main())
