#!/usr/bin/env python3
"""
55_figS2_structure.py — Supplementary Figure S2: both mechanisms on experimental coordinates.

The calcium site itself is Figure 1c. This shows the module it sits in and where the pathogenic
variants of this fragment fall, on deposited coordinates rather than on a model: PDB 2W86,
cbEGF9-hyb2-cbEGF10 at 1.8 A, crystallised in 20 mM CaCl2 so both calcium sites are occupied.

The PyMOL renders are produced by `46_v2fig8_structure.py`, which is imported for its headless
rendering machinery rather than copied. This script restyles the montage to the conventions in
`paperstyle.py` and drops the calcium-site panel, which was promoted to Figure 1.

a  the fragment: two cbEGF domains either side of the hybrid domain, six disulfides, two Ca2+
b  pathogenic ClinVar variants on the same fragment, coloured by mechanism

Writes: figures/figS2_structure.{png,pdf}, .numbers.json, logs/55_figS2_structure_<stamp>.log
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paperstyle as P  # noqa: E402

ROOT = P.ROOT
NAME = "figS2_structure"
SCRIPT = Path(__file__).name
PANELS = P.FIGURES / "panels"


def main() -> int:
    P.apply()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    N = P.Numbers(NAME, SCRIPT)
    log = [f"Supplementary Figure S2 — {NAME}"]

    # the renders come from script 46; run it if its panels are not on disk
    a_png = PANELS / "v2_fig8_structure_a_fragment.png"
    c_png = PANELS / "v2_fig8_structure_c_variants.png"
    if not (a_png.is_file() and c_png.is_file()):
        spec = importlib.util.spec_from_file_location(
            "_fig8", ROOT / "scripts" / "46_v2fig8_structure.py")
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        m.main()

    df = P.load()
    lo, hi = 807, 951
    on = df[(df.position >= lo) & (df.position <= hi) & (df.set_primary == "pathogenic")]
    n_cys = int(on.loc[on.group == "cys_removing", "position"].nunique())
    n_lig = int(on.loc[on.group == "ca_ligand", "position"].nunique())
    n_oth = int(on.loc[~on.group.isin(["cys_removing", "ca_ligand"]), "position"].nunique())
    N.set("n_pathogenic_variants", int(len(on)))
    N.set("n_positions_cys_removing", n_cys)
    N.set("n_positions_ca_ligand", n_lig)
    N.set("n_positions_other", n_oth)
    log.append(f"  {len(on)} pathogenic variants; {n_lig} Ca-ligand, {n_cys} cysteine, "
               f"{n_oth} other positions")

    fig = plt.figure(figsize=(P.W15, 1.95))
    gs = fig.add_gridspec(1, 2, wspace=0.02, left=0.005, right=0.995, top=0.985, bottom=0.20)
    for i, (png, letter, cap) in enumerate((
            (a_png, "a", "cbEGF9–hybrid 2–cbEGF10, PDB 2W86"),
            (c_png, "b", f"{len(on)} pathogenic variants, {on.position.nunique()} positions"))):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(plt.imread(png))
        ax.set_axis_off()
        ax.annotate(cap, xy=(0.5, -0.005), xycoords="axes fraction", ha="center", va="top",
                    fontsize=6.2, color=P.C["ink"])
        P.panel(ax, letter, dx=0.0, dy=0.965)

    # Swatches are read straight off the render palette in script 46 rather than repeated as
    # literals, so the key and the picture cannot drift apart. Squares are cartoon, circles are
    # atoms — the same shape-carries-class rule the plots use.
    PAL = {"af3": "#c2c2c2", "other_cbegf": "#e2e2e2", "ca_site": "#000000",
           "cbegf_cys": "#585858", "outside_cbegf": "#9a9a9a"}
    handles = [
        Line2D([], [], marker="s", ls="", ms=3.4, mfc=PAL["af3"], mec=P.C["ink2"], mew=0.3,
               label="cbEGF9 / cbEGF10"),
        Line2D([], [], marker="s", ls="", ms=3.4, mfc=PAL["other_cbegf"], mec=P.C["ink2"],
               mew=0.3, label="hybrid 2 (UniProt “TB 4”)"),
        Line2D([], [], marker="o", ls="", ms=3.4, mfc=PAL["ca_site"], mec="none",
               label=f"Ca²⁺ and its ligands ({n_lig} position)"),
        Line2D([], [], marker="o", ls="", ms=3.4, mfc=PAL["cbegf_cys"], mec="none",
               label=f"disulfide / Cys loss ({n_cys})"),
        Line2D([], [], marker="o", ls="", ms=3.4, mfc=PAL["outside_cbegf"], mec="none",
               label=f"other pathogenic ({n_oth})"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=5.8,
               bbox_to_anchor=(0.5, 0.005), handletextpad=0.3, columnspacing=1.2)

    P.save(fig, NAME, N)
    plt.close(fig)
    P.write_log("55_figS2_structure", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
