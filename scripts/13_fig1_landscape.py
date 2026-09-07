#!/usr/bin/env python3
"""
13_fig1_landscape.py — PHASE 6, Figure 1. Where the variants sit.

Three needle lanes (pathogenic / VUS / benign, primary >=2-star tier) along the 2,871 residues
of P35555, a domain-architecture track beneath them, and a per-cbEGF-domain breakdown of the
pathogenic variants by mechanism.

Design decisions worth stating:

* Separate lanes, not one mirrored track. Stacking two classes into a single stem makes the
  height ambiguous; three baselines make every count readable directly. The benign lane looks
  almost empty and that is the point — 34 benign variants at >=2 stars is the real ceiling
  (LIMITATIONS.md #2), so the figure should show it rather than hide it behind a shared axis.
* The architecture track uses lightness steps, not hues. It is annotation, not a data series,
  so it must not compete with the palette that carries meaning in Figures 2-6.
* Panel C uses the ANNOTATION-level calcium class (`is_ca_consensus`, n=414 across the gene),
  not the structure-measured `is_direct_ca_ligand` used in Figures 2-3 (n=78, limited to the
  four solved constructs). Mixing the two would be the exact "domain-resident vs
  calcium-coordinating" conflation CLAUDE.md warns about, so both the panel title and the
  recorded numbers name which definition is in play.

Writes: figures/fig1_variant_landscape.{png,pdf}, .numbers.json, logs/13_fig1_*.log
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig1_variant_landscape"
PROTEIN_LEN = 2871
NEONATAL = (659, 1362)   # exon 24-32 severe/neonatal cluster, as flagged in Phase 3


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    log: list[str] = [f"Phase 6 — Figure 1 — seed {fs.SEED}"]
    num = fs.Numbers(FIG, Path(__file__).name)

    df = fs.load_merged()
    dom = fs.domains()
    log.append(f"{len(df):,} variants; {len(dom)} annotated domains")

    if int(df.position.max()) > PROTEIN_LEN or int(df.position.min()) < 1:
        raise SystemExit("variant positions fall outside P35555 1-2871")

    num.set("protein_length", PROTEIN_LEN)
    num.set("n_variants_total", len(df))
    num.set("neonatal_region", list(NEONATAL))

    # 6 rows: three needle lanes, the architecture track, an empty spacer (row 4) that keeps
    # panel C's title clear of the track's axis labels, then the per-domain breakdown.
    fig = plt.figure(figsize=(7.4, 7.6))
    gs = fig.add_gridspec(6, 1, height_ratios=[2.0, 2.0, 1.0, 0.95, 0.5, 2.4], hspace=0.55)

    lanes = [
        ("pathogenic", "Pathogenic /\nlikely pathogenic"),
        ("vus", "Uncertain\nsignificance"),
        ("benign", "Benign /\nlikely benign"),
    ]
    axes = []
    for row, (cls, ylab) in enumerate(lanes):
        ax = fig.add_subplot(gs[row, 0])
        axes.append(ax)
        sub = df[df.set_primary == cls]
        counts = sub.groupby("position").size()
        num.set(f"n_{cls}", len(sub))
        num.set(f"n_positions_{cls}", len(counts))
        num.set(f"max_stack_{cls}", int(counts.max()) if len(counts) else 0)
        log.append(f"  {cls:11} n={len(sub):4}  distinct positions={len(counts):4}  "
                   f"tallest stem={int(counts.max()) if len(counts) else 0}")

        colour = fs.PALETTE[cls]
        ax.axvspan(*NEONATAL, color=fs.PALETTE["grid"], alpha=0.55, lw=0, zorder=0)
        if len(counts):
            ax.vlines(counts.index, 0, counts.values, color=colour, lw=0.8, alpha=0.85, zorder=2)
            ax.scatter(counts.index, counts.values, s=7, color=colour, zorder=3,
                       linewidths=0.3, edgecolors=fs.PALETTE["surface"])
        ax.set_xlim(0, PROTEIN_LEN)
        top = max(3, (int(counts.max()) if len(counts) else 0) + 1)
        ax.set_ylim(0, top)
        ax.set_ylabel(ylab, fontsize=7.5, color=fs.PALETTE["ink"])
        ax.set_yticks(range(0, top + 1, max(1, top // 3)))
        ax.tick_params(labelbottom=False)
        ax.spines["bottom"].set_visible(row == len(lanes) - 1)
        # direct label — relief for the low-contrast VUS slot, and saves a legend
        ax.text(0.995, 0.92, f"n = {len(sub)}", transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5, color=colour, fontweight="bold")

    axes[0].set_title("Pathogenic, uncertain and benign FBN1 missense variants along "
                      "fibrillin-1 (P35555)", fontsize=9, pad=14, loc="left")
    axes[0].annotate("neonatal / severe region\n(residues 659–1362)",
                     xy=(np.mean(NEONATAL), axes[0].get_ylim()[1]), xytext=(0, -1),
                     textcoords="offset points", ha="center", va="top",
                     fontsize=6.5, color=fs.PALETTE["ink2"])

    # ---- domain architecture track -------------------------------------------------
    axd = fig.add_subplot(gs[3, 0])
    axd.set_xlim(0, PROTEIN_LEN)
    axd.set_ylim(0, 1)
    axd.axhline(0.5, color=fs.PALETTE["dom_other"], lw=3.5, zorder=1,
                solid_capstyle="butt")
    kind_colour = {"cbEGF": fs.PALETTE["dom_cbEGF"], "EGF": fs.PALETTE["dom_EGF"],
                   "TB": fs.PALETTE["dom_TB"]}
    for _, d in dom.iterrows():
        axd.add_patch(Rectangle((d.start, 0.27), d.end - d.start, 0.46,
                                facecolor=kind_colour.get(d.kind, fs.PALETTE["dom_other"]),
                                edgecolor=fs.PALETTE["surface"], linewidth=0.35, zorder=2))
    for _, d in dom[dom.kind == "cbEGF"].iterrows():
        if int(d.cbegf_index) % 5 == 0 or int(d.cbegf_index) == 1:
            axd.text((d.start + d.end) / 2, 0.06, str(int(d.cbegf_index)),
                     ha="center", va="bottom", fontsize=5.8, color=fs.PALETTE["ink2"])
    axd.text(0.0, 0.92, "domain architecture  (numbers = cbEGF index)", transform=axd.transAxes,
             fontsize=6.8, color=fs.PALETTE["ink2"], va="bottom")
    axd.set_yticks([])
    for s in ("left", "top", "right"):
        axd.spines[s].set_visible(False)
    axd.set_xlabel("Residue position in P35555")
    counts_kind = dom.kind.value_counts().to_dict()
    handles = [Rectangle((0, 0), 1, 1, facecolor=kind_colour[k], edgecolor="none")
               for k in ("cbEGF", "EGF", "TB")]
    axd.legend(handles, [f"cbEGF ({counts_kind.get('cbEGF', 0)})",
                         f"EGF ({counts_kind.get('EGF', 0)})",
                         f"TB / hybrid ({counts_kind.get('TB', 0)})"],
               loc="lower right", bbox_to_anchor=(1.0, 0.86), ncol=3, handlelength=1.0,
               handleheight=0.8, columnspacing=1.0, fontsize=6.8, borderpad=0.15)
    for k, v in counts_kind.items():
        num.set(f"n_domains_{k}", int(v))

    # ---- panel C: pathogenic variants per cbEGF domain, by mechanism ----------------
    axc = fig.add_subplot(gs[5, 0])
    path = df[df.set_primary == "pathogenic"]
    cys = path.in_cbegf_cys & path.cys_removing
    cls_of = np.where(cys, "cys", np.where(path.is_ca_consensus, "ca", "other"))
    path = path.assign(_cls=cls_of)

    idx = np.arange(1, 44)
    series = {
        "cys": ("cbEGF cysteine-removing", fs.PALETTE["cbegf_cys"]),
        "ca": ("Calcium-consensus residue", fs.PALETTE["ca_site"]),
        "other": ("Other cbEGF residue", fs.PALETTE["other_cbegf"]),
    }
    bottom = np.zeros(len(idx))
    for key, (label, colour) in series.items():
        vals = np.array([int(((path.cbegf_index == i) & (path._cls == key)).sum())
                         for i in idx], dtype=float)
        axc.bar(idx, vals, bottom=bottom, width=0.72, color=colour,
                edgecolor=fs.PALETTE["surface"], linewidth=0.6, label=label, zorder=2)
        bottom += vals
        num.set(f"panelC_total_{key}", int(vals.sum()))
        log.append(f"  panel C {key:6}: {int(vals.sum())} pathogenic variants")

    num.set("panelC_pathogenic_in_cbegf", int(path.in_cbegf.sum()))
    num.set("panelC_pathogenic_total", len(path))
    axc.set_xlim(0.3, 43.7)
    axc.set_xticks([1] + list(range(5, 44, 5)))
    axc.set_xlabel("cbEGF domain")
    axc.set_ylabel("Pathogenic variants")
    axc.yaxis.grid(True, zorder=0)
    axc.set_axisbelow(True)
    axc.legend(loc="upper left", ncol=3, handlelength=1.0, handleheight=0.8,
               columnspacing=1.2, borderpad=0.15)
    axc.set_title(f"Pathogenic variants per cbEGF domain, by affected site  "
                  f"({int(path.in_cbegf.sum())} of {len(path)} lie in a cbEGF domain)",
                  fontsize=8, pad=16, loc="left")
    axc.set_ylim(0, bottom.max() * 1.42)

    fs.panel_label(axes[0], "A", dx=-0.075, dy=1.36)
    fs.panel_label(axd, "B", dx=-0.075, dy=1.10)
    fs.panel_label(axc, "C", dx=-0.075, dy=1.10)

    fs.save(fig, FIG, num)
    plt.close(fig)
    log.append(f"wrote figures/{FIG}.png")
    print(f"  log -> {fs.write_log('13_fig1', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
