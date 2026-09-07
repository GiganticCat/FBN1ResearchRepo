#!/usr/bin/env python3
"""
51_fig2_stability.py — Figure 2: cysteine loss destabilises the domain; calcium-ligand loss does not.

Laid out to the conventions on the pages in `resources/papers/`. Three things about it are
deliberate and were wrong in earlier versions of this figure.

  CATEGORIES ARE ON THE Y AXIS. Four descriptive category names do not fit side by side beneath a
  half-column panel; the earlier version collided them into "no Ca2+ contactno Ca2+ contact".
  Cheng 2023 (Science) puts thirteen method names on a vertical axis for exactly this reason.

  THE CATEGORIES ARE NAMED AS THE LITERATURE NAMES THEM. Jensen 2009 calls a coordinating
  residue a "Ca2+ ligand", so the rows read "Ca2+ ligand (side chain)", "Non-ligand Asp/Glu/Asn",
  "Consensus motif, non-ligand" and "Cysteine substitution". Earlier versions labelled the rows
  with a description of how each group was computed, which is not a name a reader recognises.

  COLOUR IS NOT USED TO NAME A GROUP. The rows already name themselves. Fill separates the two
  structural classes (solid) from the two controls (grey, open), which is what the reference
  pages do and what survives greyscale printing.

Effect sizes and corrected q values live in `results/statistics_final.tsv`; the panel carries the
comparison the reader needs to see, drawn on the data.

Writes: figures/fig2_stability.{png,pdf}, .numbers.json, logs/51_fig2_stability_<stamp>.log
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

NAME = "fig2_stability"
SCRIPT = Path(__file__).name

ORDER = ["ca_ligand", "composition_control", "motif_overcalled", "cys_removing"]
FILL = {"ca_ligand": P.C["ink"], "composition_control": P.C["grey"],
        "motif_overcalled": "open", "cys_removing": P.C["ink"]}


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(P.SEED)
    N = P.Numbers(NAME, SCRIPT)
    log = [f"Figure 2 — {NAME}"]

    df = P.load()

    def groups_for(col):
        return [(P.GLABEL[g], df.loc[df.group == g, col].to_numpy(), FILL[g]) for g in ORDER]

    fig = plt.figure(figsize=(P.W2, 3.55))
    gs = fig.add_gridspec(2, 2, wspace=0.30, hspace=0.62,
                          left=0.185, right=0.988, top=0.915, bottom=0.095)

    # ------------------------------------------------------------------ a  Rosetta
    ax = fig.add_subplot(gs[0, 0])
    g = groups_for("rosetta_ddG")
    stats = P.striph(ax, g, "ΔΔG (Rosetta energy units)", rng)
    ax.set_xlim(-6, 40)
    P.clip_note_h(ax, g, 34)
    ax.axvline(0, color=P.C["rule"], lw=0.5, zorder=1)
    P.bracket_h(ax, 3, 2, 20.0,
                P.pq("rosetta_ddG (REU): Ca ligands vs composition control"), dx=1.1)
    # dropped off the midpoint, which is the consensus-motif row and carries its clipped count
    P.bracket_h(ax, 2, 0, 20.0,
                P.pq("rosetta_ddG (REU): cysteine-removing vs composition control"),
                dx=1.1, ly=0.62)
    P.title(ax, "Rosetta ref2015, explicit Ca²⁺ bonds")
    P.panel(ax, "a", dx=-0.42)
    for k, v in stats.items():
        N.set(f"a_rosetta_{k.replace(chr(10), ' ')}_median", round(v["median"], 3))
        N.set(f"a_rosetta_{k.replace(chr(10), ' ')}_n", v["n"])
    log.append(f"  a  {stats}")

    # ------------------------------------------------------------------ b  FoldX
    ax = fig.add_subplot(gs[0, 1])
    g = groups_for("foldx_ddG")
    stats = P.striph(ax, g, "ΔΔG (kcal mol⁻¹)", rng, ms=3.0, alpha=0.55)
    ax.set_xlim(-3, 17.5)
    P.clip_note_h(ax, g, 15)
    ax.axvline(0, color=P.C["rule"], lw=0.5, zorder=1)
    P.bracket_h(ax, 3, 2, 8.2,
                P.pq("foldx_ddG (kcal/mol): Ca ligands vs composition control"), dx=0.45)
    P.bracket_h(ax, 2, 0, 8.2,
                P.pq("foldx_ddG (kcal/mol): cysteine-removing vs composition control"), dx=0.45)
    P.title(ax, "FoldX 5.1, experimental templates")
    P.panel(ax, "b", dx=-0.42)
    for k, v in stats.items():
        N.set(f"b_foldx_{k.replace(chr(10), ' ')}_median", round(v["median"], 3))
        N.set(f"b_foldx_{k.replace(chr(10), ' ')}_n", v["n"])

    # ------------------------------------------------------------------ c  burial control
    ax = fig.add_subplot(gs[1, 0])
    # `key` is what the sidecar and the gate have always called each band; `shown` is what the
    # panel prints. Relative solvent accessibility is abbreviated RSA, as in the Methods.
    bands = [("buried", "RSA < 20%", -1, 20), ("partly buried", "RSA 20–50%", 20, 50)]
    rows = []
    for key, shown, lo, hi in bands:
        sub = df[(df.rsa >= lo) & (df.rsa < hi)]
        for grp, gl in (("ca_ligand", "Ca²⁺ ligand"),
                        ("composition_control", "Non-ligand Asp/Glu/Asn")):
            v = sub.loc[sub.group == grp, "rosetta_ddG"].dropna().to_numpy()
            rows.append((f"{gl},\n{shown}", v,
                         P.C["ink"] if grp == "ca_ligand" else P.C["grey"]))
            N.set(f"c_{key}_{grp}_median", round(float(np.median(v)), 3))
            N.set(f"c_{key}_{grp}_n", int(v.size))
    stats = P.striph(ax, rows, "ΔΔG (Rosetta energy units)", rng, ms=2.2)
    ax.set_xlim(-6, 38)
    ax.axvline(0, color=P.C["rule"], lw=0.5, zorder=1)
    P.bracket_h(ax, 3, 2, 22.0, P.pq(
        "rosetta_ddG within buried (RSA<20%): Ca ligands vs composition control"), dx=1.1)
    P.bracket_h(ax, 1, 0, 22.0, P.pq(
        "rosetta_ddG within partial (20-50%): Ca ligands vs composition control"), dx=1.1)
    P.title(ax, "Rosetta ΔΔG stratified by RSA")
    P.panel(ax, "c", dx=-0.42)

    # ------------------------------------------------------------------ d  AlphaMissense
    ax = fig.add_subplot(gs[1, 1])
    g = groups_for("am")
    stats = P.striph(ax, g, "AlphaMissense score", rng)
    ax.set_xlim(-0.03, 1.09)
    ax.axvline(0.564, color=P.C["rule"], lw=0.6, ls=(0, (2.5, 2)), zorder=1)
    ax.annotate("pathogenic\nthreshold", xy=(0.564, -0.60), xytext=(-2, 0),
                textcoords="offset points", ha="right", va="bottom", fontsize=5.4,
                color=P.C["muted"], linespacing=1.25)
    P.bracket_h(ax, 3, 2, 1.045,
                P.pq("AlphaMissense: Ca ligands vs composition control"), dx=0.028)
    P.title(ax, "AlphaMissense pathogenicity score")
    P.panel(ax, "d", dx=-0.42)
    for k, v in stats.items():
        N.set(f"d_am_{k.replace(chr(10), ' ')}_median", round(v["median"], 4))
        N.set(f"d_am_{k.replace(chr(10), ' ')}_n", v["n"])
    # The cysteine row is 1,132 here and 1,038 in a and c. AlphaMissense is a sequence method and
    # scores every variant, while a ddG needs a modelled construct, and 94 cysteine substitutions
    # lie outside the span of the 21 constructs. Recorded so the gate ties the two together.
    n_cys_uncovered = int(((df.group == "cys_removing") & df.rosetta_ddG.isna()).sum())
    N.set("cys_without_structural_coverage", n_cys_uncovered)
    log.append(f"  cysteine variants without a modelled construct: {n_cys_uncovered}")

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("51_fig2_stability", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
