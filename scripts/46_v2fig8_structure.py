#!/usr/bin/env python3
"""
46_v2fig8_structure.py — Supplementary Figure S2: the two mechanisms on experimental coordinates.

Everything else in this study is a number. This is the picture the numbers are about, drawn from
deposited coordinates rather than from a model: PDB 2W86, cbEGF9-hyb2-cbEGF10 at 1.8 A,
crystallised in 20 mM CaCl2 so both calcium sites are occupied (Jensen 2009).

It rebuilds `16_fig4_structure.py` for v2, and the rebuild is not cosmetic:

  THE CYSTEINE SET WAS UNDERCOUNTED. v1 coloured only the 258 conserved *cbEGF* cysteines, so on
  this fragment it drew 12 pathogenic cysteine positions. The hybrid domain's cysteines are just
  as real and just as disulfide-bonded; v2's definition is any cysteine-removing substitution,
  and the correct figure for this fragment is 23 positions.

  THE AF3 OVERLAY PANEL IS GONE. It duplicated Figure 5, which now carries model-versus-
  experiment validation across all eight calibration domains rather than one.

  JOURNAL FORMAT. No panel titles, lowercase panel letters, double-column width — the
  conventions the rest of the v2 set follows and `tests/test_v2.py` F3 enforces.

Panel b is also a check on the analysis: the six residues drawn coordinating cbEGF9's ion are the
six the crystallographers report, and they are the six the AF3 models return independently.

The PyMOL machinery -- headless rendering of committed .pml scripts, deterministic cameras, white
-margin cropping -- is imported from `16_fig4_structure.py` rather than copied, so the two figures
cannot drift in style.

Writes: figures/v2_fig8_structure.{png,pdf}, .numbers.json, figures/pml/v2_fig8_*.pml,
        figures/panels/v2_fig8_*.png, logs/46_v2fig8_structure_<stamp>.log
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import figstyle as F  # noqa: E402
import paperstyle as P  # noqa: E402

spec = importlib.util.spec_from_file_location("_fig4_v1", ROOT / "scripts" / "16_fig4_structure.py")
V1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(V1)

# v1 carries the backbone donor positions but not their residue names; naming them is what makes
# panel b checkable against the crystallographers' own list (Jensen 2009).
V1.CA9_BACKBONE_LABEL = {808: "Ile808", 824: "Ser824", 827: "Ser827"}

# The renders inherit `figstyle.PALETTE` through script 16. Override the entries they use with
# the monochrome set in `paperstyle`, so the structure panels match the plots.
#
# These are four steps on one lightness ramp, ordered so that salience tracks importance: the
# calcium and its ligands are black, the disulfides mid-dark, unassigned pathogenic positions
# mid, and the protein itself pale enough that every marked atom sits clearly on top of it.
# Panel b recolours the whole cartoon grey85, so all three sphere classes there are darker than
# their background by construction. The old scheme used one blue and the conventional sulfur
# yellow; both are gone, and nothing was encoded by them that a lightness step cannot carry.
V1.fs.PALETTE.update({
    "ca_site": "#000000",          # calcium and its ligands — the darkest thing in the panel
    "cbegf_cys": "#585858",        # disulfides, and cysteine-removing variant positions
    "af3": "#c2c2c2",              # cbEGF domains
    "other_cbegf": "#e2e2e2",      # hybrid domain
    "outside_cbegf": "#9a9a9a",    # pathogenic variants of no assigned mechanism
    "experimental": "#3d3d3d",
})

NAME = "v2_fig8_structure"
SCRIPT = Path(__file__).name


def main() -> int:
    F.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    log = [f"Supplementary Figure S2 — {NAME}", f"PyMOL {F.pymol_version()}"]
    N = F.Numbers(NAME, SCRIPT)
    V1.PML_DIR.mkdir(parents=True, exist_ok=True)
    V1.PANEL_DIR.mkdir(parents=True, exist_ok=True)
    if not V1.EXP.is_file():
        raise SystemExit(f"missing {V1.EXP}")

    # ---- which variants sit on this fragment, under the v2 definitions
    df = F.load_v2()
    df["mech"] = F.mech_group_v2(df)
    lo, hi = V1.CONSTRUCT
    sub = df[(df.position >= lo) & (df.position <= hi)]
    path = sub[sub.set_primary == "pathogenic"]
    pos_lig = sorted(int(x) for x in path[path.mech == "ca_ligand"].position.unique())
    pos_cys = sorted(int(x) for x in path[path.mech == "cys_removing"].position.unique())
    pos_oth = sorted(int(x) for x in path[~path.mech.isin(["ca_ligand", "cys_removing"])]
                     .position.unique())
    N.set("template_pdb", "2W86")
    N.set("construct_range", f"{lo}-{hi}")
    N.set("n_pathogenic_variants", int(len(path)))
    N.set("n_positions_ca_ligand", len(pos_lig))
    N.set("n_positions_cys_removing", len(pos_cys))
    N.set("n_positions_other", len(pos_oth))
    N.set("positions_ca_ligand", pos_lig)
    log.append(f"2W86 fragment {lo}-{hi}: {len(path)} pathogenic variants at "
               f"{path.position.nunique()} positions "
               f"({len(pos_lig)} Ca ligand, {len(pos_cys)} cysteine, {len(pos_oth)} other)")

    ss_sel = "+".join(str(r) for pair in V1.SS_BONDS for r in pair)
    common = f"""
load {V1.EXP}, exp
remove exp and resi {V1.CLONING_ARTEFACT[0]}-{V1.CLONING_ARTEFACT[1]}
remove exp and solvent
hide everything
select exp9,  exp and resi {V1.CBEGF9[0]}-{V1.CBEGF9[1]}
select exphy, exp and resi {V1.HYBRID2[0]}-{V1.HYBRID2[1]}
select exp10, exp and resi {V1.CBEGF10[0]}-{V1.CBEGF10[1]}
"""

    # ---- a: the fragment
    pA = V1.write_and_run(f"{NAME}_a_fragment", common + f"""
show cartoon, exp and polymer
# the inter-domain linkers belong to no selection below; colour everything first so they can
# never come out in PyMOL's default green
color grey70,    exp and polymer
color pal_af3,   exp9
color pal_other, exphy
color pal_af3,   exp10
show sticks, exp and resi {ss_sel} and not name N+C+O
set stick_radius, 0.17, exp and resi {ss_sel}
color pal_cys, exp and resi {ss_sel} and elem S
color pal_cys, exp and resi {ss_sel} and elem C
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
orient exp and polymer
turn x, -10
zoom exp and polymer, 1.0
""", log)

    # ---- b: the cbEGF9 calcium site
    side = "+".join(str(r) for r in V1.CA9_SIDECHAIN)
    back = "+".join(str(r) for r in V1.CA9_BACKBONE)
    labels = "\n".join(
        f'label exp and resi {r} and name CB, "{n}"' for r, n in V1.CA9_SIDECHAIN.items())
    back_labels = "\n".join(
        f'label exp and resi {r} and name O, "{n}"' for r, n in V1.CA9_BACKBONE_LABEL.items())
    pB = V1.write_and_run(f"{NAME}_b_casite", common + f"""
select ca9, exp and elem Ca within 6 of (exp9 and name CA)
# only the site's own neighbourhood: the rest of the domain is a slab that crowds the labels
select site9, byres (exp9 within 9.5 of ca9)
show cartoon, site9
set cartoon_transparency, 0.55, site9
color grey80, site9
show sticks, exp and resi {side} and sidechain
color pal_ca, exp and resi {side} and elem O
color grey50, exp and resi {side} and elem C
show sticks, exp and resi {back} and name C+O
color grey60, exp and resi {back} and name C+O
show spheres, ca9
color pal_ca, ca9
set sphere_scale, 0.22, ca9
distance d_side, ca9, (exp and resi {side} and sidechain and elem O), 3.2
distance d_back, ca9, (exp and resi {back} and name O), 3.2
hide labels, d_side
hide labels, d_back
{labels}
{back_labels}
# Asn823 and Ser827 sit almost on top of each other from this angle; nudge them apart rather
# than choosing a camera that hides one of the six ligands
set label_position, (-1.6, -1.9, 2.2), exp and resi 823
set label_position, ( 1.9,  1.2, 2.2), exp and resi 827
set label_size, 15
orient site9
turn y, 15
turn x, 8
center ca9
zoom ca9, 6.2
""", log)

    # ---- c: pathogenic variants by mechanism
    sel = lambda p: "+".join(str(x) for x in p) if p else "0"  # noqa: E731
    pC = V1.write_and_run(f"{NAME}_c_variants", common + f"""
show cartoon, exp and polymer
color grey85, exp and polymer
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
set sphere_scale, 0.35, elem Ca
show spheres, exp and resi {sel(pos_oth)} and name CA
color pal_outside, exp and resi {sel(pos_oth)} and name CA
show spheres, exp and resi {sel(pos_cys)} and name CA
color pal_cys, exp and resi {sel(pos_cys)} and name CA
show spheres, exp and resi {sel(pos_lig)} and name CA
color pal_ca, exp and resi {sel(pos_lig)} and name CA
set sphere_scale, 0.85, name CA
orient exp and polymer
turn x, -10
zoom exp and polymer, 1.0
""", log)

    # ---- montage
    fig, axs = plt.subplots(1, 3, figsize=(F.COL_2, 2.62))
    fig.subplots_adjust(left=0.005, right=0.995, top=0.955, bottom=0.115, wspace=0.02)
    for ax, png, letter in zip(axs, (pA, pB, pC), "abc"):
        ax.imshow(plt.imread(png))
        ax.set_axis_off()
        F.panel_label(ax, letter, dx=0.0, dy=0.985)

    handles = [
        Line2D([], [], marker="o", ls="", ms=4, mfc=P.C["ca"],
               mec="none", label="Ca²⁺ and its ligands"),
        Line2D([], [], marker="o", ls="", ms=4, mfc="#b8860b",
               mec="none", label="disulfide / Cys-removing"),
        Line2D([], [], marker="o", ls="", ms=4, mfc="#8a8a8a",
               mec="none", label="other pathogenic"),
        Line2D([], [], marker="s", ls="", ms=4, mfc="#6f6f6f",
               mec="none", label="cbEGF9 / cbEGF10"),
        Line2D([], [], marker="s", ls="", ms=4, mfc="#b5b5b5",
               mec="none", label="hybrid 2 (UniProt “TB 4”)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=5.8,
               bbox_to_anchor=(0.5, 0.002), handletextpad=0.35, columnspacing=1.4)

    F.save(fig, NAME, N)
    plt.close(fig)
    F.write_log("46_v2fig8_structure", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
