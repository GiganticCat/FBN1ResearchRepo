#!/usr/bin/env python3
"""
50_fig1_module.py — Figure 1: the cbEGF module, and where variants fall on it.

The figure the earlier draft did not have. Its Figure 1 spent two of three panels on numbers that
belong in a sentence ("751 versus 3,831 variants"; four group sizes as a bar chart) and gave the
third to three rows of tick marks in which the uncertain-significance row was a solid block of
ink. A reader learned the protein is long and that variants are everywhere.

What a reader actually needs before anything else is the module: six cysteines in a fixed
disulfide pattern, and a small set of residues at the N-terminal end that grip a calcium ion.
Panel b draws that, and — this is the point — the positions are not a cartoon. They are the
modal offsets measured across all 43 domains, and the histogram beneath is where pathogenic
variants actually fall on the same axis. The two peaks are the paper.

a  variants along fibrillin-1, pathogenic above the domain architecture and the gnomAD
   population comparison below it, on a shared axis
b  the consensus cbEGF module, with pathogenic variant counts aligned to the same positions
c  the calcium site of cbEGF9 in the 1.8 A crystal structure of the cbEGF9-hyb2-cbEGF10 fragment

Writes: figures/fig1_module.{png,pdf}, .numbers.json, logs/50_fig1_module_<stamp>.log
"""
from __future__ import annotations

import collections
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

NAME = "fig1_module"
SCRIPT = Path(__file__).name
CASITE_PANEL = P.FIGURES / "panels" / "v2_fig8_structure_b_casite.png"


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc, Rectangle

    rng = np.random.default_rng(P.SEED)
    N = P.Numbers(NAME, SCRIPT)
    log = [f"Figure 1 — {NAME}"]

    df = P.load()
    dom = P.domains()
    pop = pd.read_csv(P.PROC / "population_set.tsv", sep="\t")
    covered = set(df.loc[df.has_structure, "position"])
    pop = pop[pop.position.isin(covered)]
    path = df[df.set_primary == "pathogenic"]
    N.set("n_pathogenic", int(len(path)))
    N.set("n_population", int(len(pop)))
    N.set("n_cbegf_domains", int((dom.kind == "cbEGF").sum()))

    fig = plt.figure(figsize=(P.W2, 4.35))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.28], width_ratios=[1.42, 1.0],
                          hspace=0.36, wspace=0.20,
                          left=0.062, right=0.988, top=0.925, bottom=0.075)

    # ================================================================== a  along the protein
    ax = fig.add_subplot(gs[0, :])
    L = 2871

    def needles(positions, up: bool, color: str):
        c = collections.Counter(int(p) for p in positions)
        xs = np.array(sorted(c))
        ys = np.array([c[x] for x in xs], float) * (1 if up else -1)
        ax.vlines(xs, 0, ys, color=color, lw=0.45, alpha=0.85, zorder=3)
        return c

    cpath = needles(path.position, True, P.C["ink"])
    cpop = needles(pop.position, False, P.C["grey"])
    top = max(cpath.values())
    bot = -max(cpop.values())

    # domain architecture ribbon, on the zero line
    hh = top * 0.085
    kind_color = {"cbEGF": P.C["dom_cbegf"], "EGF": P.C["dom_egf"], "TB": P.C["dom_tb"]}
    ax.add_patch(Rectangle((1, -hh), L, 2 * hh, facecolor=P.C["dom_other"],
                           edgecolor="none", zorder=2))
    for r in dom.itertuples():
        ax.add_patch(Rectangle((r.start, -hh), r.end - r.start, 2 * hh,
                               facecolor=kind_color.get(r.kind, P.C["dom_other"]),
                               edgecolor="none", zorder=3))

    ax.set_xlim(-30, L + 30)
    ax.set_ylim(bot * 1.22, top * 1.22)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 2))
    ax.set_xticks([1, 500, 1000, 1500, 2000, 2500, 2871])
    ax.set_xlabel("Residue in fibrillin-1 (UniProt P35555)", labelpad=1.5)

    P.label_at(ax, 30, top * 1.13, f"ClinVar pathogenic / likely pathogenic, n = {len(path)}",
               P.C["ink"], va="center", fontweight="bold")
    P.label_at(ax, 30, bot * 1.13, f"gnomAD v4 population comparison, n = {len(pop):,}",
               P.C["grey"], va="center", fontweight="bold")
    # architecture key: swatches on the ribbon's own line, at the free right-hand end
    for i, (lab, col) in enumerate((("cbEGF (43)", P.C["dom_cbegf"]), ("EGF", P.C["dom_egf"]),
                                    ("TB / hybrid", P.C["dom_tb"]))):
        x0 = 2000 + i * 300
        ax.add_patch(Rectangle((x0, bot * 0.72), 62, top * 0.085,
                               facecolor=col, edgecolor="none", zorder=5))
        ax.annotate(lab, xy=(x0 + 78, bot * 0.72 + top * 0.04), ha="left", va="center",
                    fontsize=5.6, color=P.C["ink2"], zorder=5)
    ax.axhline(0, color=P.C["ink"], lw=0.4, zorder=4)
    P.panel(ax, "a", dx=-0.048, dy=1.00)
    log.append(f"  a  {len(path)} pathogenic, {len(pop)} population, {(dom.kind=='cbEGF').sum()} cbEGF")

    # ================================================================== b  the module
    ax = fig.add_subplot(gs[1, 0])
    cb = dom[dom.kind == "cbEGF"]
    sites = pd.read_csv(P.INTERIM / "cbegf_sites.tsv", sep="\t")
    lig = pd.read_csv(P.PROC / "af3_ca_ligands.tsv", sep="\t")

    # cysteine offsets: median over domains, by rank in sequence (C1..C6)
    ranks = collections.defaultdict(list)
    for r in sites.itertuples():
        offs = sorted(int(c) - r.start for c in str(r.cys).split(";"))
        for i, o in enumerate(offs):
            ranks[i].append(o)
    cys_off = [float(np.median(ranks[i])) for i in range(6)]
    N.set("cys_offsets", [round(x, 1) for x in cys_off])

    # calcium donor offsets: measured in the AF3 models, mapped into each domain
    rows = []
    for r in lig.itertuples():
        d = cb[(cb.start <= r.uniprot_pos) & (cb.end >= r.uniprot_pos)]
        if len(d) == 1:
            rows.append((int(d.iloc[0].cbegf_index), r.uniprot_pos - d.iloc[0].start,
                         r.wt_aa, r.atom_class))
    # one row per (domain, offset, class): a residue contributing two oxygens to the same ion
    # must not count twice, and the tally is "in how many of the 43 domains", not "how many atoms"
    Ld = pd.DataFrame(rows, columns=["idx", "off", "aa", "cls"]).drop_duplicates(
        ["idx", "off", "cls"])
    side = collections.Counter(Ld.loc[Ld.cls == "sidechain", "off"])
    back = collections.Counter(Ld.loc[Ld.cls == "backbone", "off"])
    side_off = sorted(o for o, n in side.items() if n >= 10)
    back_off = sorted(o for o, n in back.items() if n >= 10)
    N.set("sidechain_donor_offsets", side_off)
    N.set("backbone_donor_offsets", back_off)
    log.append(f"  b  Cys offsets {cys_off}; side-chain donors {side_off}; "
               f"backbone donors {back_off}")

    # pathogenic variants per offset within cbEGF domains
    poff = []
    for r in path.itertuples():
        d = cb[(cb.start <= r.position) & (cb.end >= r.position)]
        if len(d) == 1:
            poff.append(r.position - d.iloc[0].start)
    hist = collections.Counter(poff)
    xmax = 43
    counts = np.array([hist.get(i, 0) for i in range(xmax + 1)], float)
    N.set("n_pathogenic_in_cbegf", int(sum(hist.values())))
    N.set("peak_offset_counts", {str(o): int(hist.get(o, 0)) for o in
                                 sorted(set(int(x) for x in cys_off) | set(side_off))})

    hmax = counts.max()
    base = hmax * 1.55          # the schematic sits above the histogram
    ax.bar(range(xmax + 1), counts, width=0.78, color=P.C["grey_lt"], linewidth=0, zorder=2)
    # Three categories share one bar series, so they separate by lightness: cysteines black,
    # calcium donors mid grey, everything else light grey. The schematic directly above aligns
    # a marker over each highlighted bar, so the reader never has to match a shade to a key.
    for o in [int(round(x)) for x in cys_off]:
        ax.bar([o], [counts[o]], width=0.78, color=P.C["ink"], linewidth=0, zorder=3)
    for o in side_off:
        ax.bar([o], [counts[o]], width=0.78, color=P.C["ca_mid"], linewidth=0, zorder=3)

    # the module backbone
    ax.plot([0, xmax], [base, base], color=P.C["ink2"], lw=0.9, zorder=3,
            solid_capstyle="round")
    # disulfides 1-3, 2-4, 5-6 as arcs beneath the line
    for a, b in ((0, 2), (1, 3), (4, 5)):
        x1, x2 = cys_off[a], cys_off[b]
        w = x2 - x1
        ax.add_patch(Arc(((x1 + x2) / 2, base), w, hmax * 0.62, theta1=180, theta2=360,
                         edgecolor=P.C["ink"], lw=0.8, zorder=4))
    for i, o in enumerate(cys_off):
        ax.plot([o], [base], marker="o", ms=3.4, mfc=P.C["ink"], mec="white", mew=0.5, zorder=6)
        ax.annotate(f"C{i+1}", xy=(o, base), xytext=(0, 4.5), textcoords="offset points",
                    ha="center", va="bottom", fontsize=5.8, color=P.C["ink"])
    # Calcium donors are DIAMONDS and cysteines are CIRCLES. Both are black, so shape is what
    # separates them; filled versus open then carries side-chain versus backbone donor. Hue did
    # this job in an earlier version and nothing else has to change to replace it.
    for o in side_off:
        ax.plot([o], [base], marker="D", ms=2.9, mfc=P.C["ca"], mec="white", mew=0.4, zorder=6)
    for o in back_off:
        ax.plot([o], [base], marker="D", ms=2.7, mfc="white", mec=P.C["ca"], mew=0.75, zorder=6)

    # name the three side-chain donors by their consensus residue
    consensus = {}
    for o in side_off:
        sub = Ld[(Ld.cls == "sidechain") & (Ld.off == o)]
        consensus[o] = collections.Counter(sub.aa).most_common(1)[0][0]
    N.set("sidechain_donor_residues", {str(k): v for k, v in consensus.items()})
    # adjacent donors of the same identity share one label rather than printing over each other
    aa3 = {"D": "Asp", "E": "Glu", "N": "Asn"}
    runs, cur = [], [side_off[0]]
    for o in side_off[1:]:
        if o - cur[-1] <= 1 and consensus[o] == consensus[cur[-1]]:
            cur.append(o)
        else:
            runs.append(cur); cur = [o]
    runs.append(cur)
    for run in runs:
        ax.annotate(aa3.get(consensus[run[0]], consensus[run[0]]),
                    xy=(float(np.mean(run)), base), xytext=(0, -5.5),
                    textcoords="offset points", ha="center", va="top", fontsize=5.8,
                    color=P.C["ca"])
    ax.annotate("Ca²⁺ site", xy=(float(np.mean(side_off + back_off)), base), xytext=(0, 13),
                textcoords="offset points", ha="center", va="bottom", fontsize=6.6,
                color=P.C["ca"], fontweight="bold")
    ax.annotate("Disulfide bonds\nC1–C3, C2–C4, C5–C6", xy=(31, base - hmax * 0.36),
                xytext=(0, 0), textcoords="offset points", ha="center", va="top",
                fontsize=5.8, color=P.C["ink"], linespacing=1.25)
    ax.annotate("Ca²⁺ ligands   ◆ side-chain donor     ◇ main-chain carbonyl donor",
                xy=(xmax + 1.4, base + hmax * 0.42), ha="right", va="bottom", fontsize=5.4,
                color=P.C["ink2"])

    ax.set_xlim(-1.6, xmax + 1.6)
    ax.set_ylim(0, base * 1.30)
    ax.set_yticks(np.arange(0, hmax + 1, 20 if hmax > 60 else 10))
    ax.set_ylim(0, base * 1.30)
    ax.spines["left"].set_bounds(0, hmax)
    ax.set_ylabel("Pathogenic\nvariants (n)", labelpad=2)
    ax.set_xlabel("Residue offset from the cbEGF domain N terminus", labelpad=1.5)
    ax.set_xticks([0, 10, 20, 30, 40])
    P.panel(ax, "b", dx=-0.135, dy=1.00)

    # ================================================================== c  the real site
    ax = fig.add_subplot(gs[1, 1])
    if not CASITE_PANEL.is_file():
        raise SystemExit(f"missing render {CASITE_PANEL} — run scripts/46_v2fig8_structure.py")
    ax.imshow(plt.imread(CASITE_PANEL))
    ax.set_axis_off()
    ax.annotate("cbEGF9 Ca²⁺ site, PDB 2W86 (1.8 Å)", xy=(0.5, -0.01),
                xycoords="axes fraction", ha="center", va="top", fontsize=6.4,
                color=P.C["ink"])
    P.panel(ax, "c", dx=0.0, dy=0.995)

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("50_fig1_module", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
