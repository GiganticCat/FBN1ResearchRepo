#!/usr/bin/env python3
"""
16_fig4_structure.py — PHASE 6, Figure 4. The structure itself.

Renders four PyMOL panels of the cbEGF9–hybrid2–cbEGF10 fragment (PDB 2W86, X-ray 1.8 A,
crystallised with calcium) and assembles them into one figure:

  A  the fragment          two cbEGF domains either side of a hybrid domain, 2 Ca2+, 6 disulfides
  B  the calcium site      cbEGF9's ion with its coordinating side chains and backbone carbonyls
  C  AF3 vs experiment     the cofolded model superposed on the crystal, both ions shown
  D  variants              ClinVar pathogenic positions, coloured by mechanism

Why 2W86 and not another entry: it is the highest-quality calcium-bound FBN1 fragment available
(X-ray, crystallised in Ca2+ rather than titrated), it carries the largest share of structurally
covered variants (254), and it is where AF3 cofolding calibrated best (0.13 A on cbEGF9's ion).
Using an NMR reference for panel C would show the reference's own coordinate spread rather than
the method's accuracy.

Reproducibility: PyMOL runs headless via subprocess on generated `.pml` scripts, which are kept
in `figures/pml/` alongside the panels they produce. Camera state is set with `orient` on a named
selection plus fixed `turn`/`zoom` commands — deterministic across machines, unlike a
hand-posed view. Each panel writes its resulting view matrix into the run log.

Numbering: the packet PDBs were renumbered into UniProt P35555 numbering in Phase 4, so every
residue label in this figure is directly comparable with the tables. 2W86's two non-native
cloning residues (renumbered 805-806) are removed before rendering rather than drawn as if they
were fibrillin.

Naming note: the construct is called cbEGF9-hyb2-cbEGF10 in the literature (Jensen 2009), but
UniProt annotates residues 851-902 as "TB 4", not as a hybrid domain. The panel labels carry
both names rather than silently picking one.

Writes: figures/fig4_structure.{png,pdf}, .numbers.json, figures/pml/*.pml,
        figures/panels/fig4_*.png, logs/16_fig4_*.log
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig4_structure"
PML_DIR = fs.FIGURES / "pml"
PANEL_DIR = fs.FIGURES / "panels"
EXP = fs.ROOT / "handoff/cmm/inbox/EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb"
AF3 = fs.ROOT / "handoff/cmm/inbox/AF3_cbEGF9-hyb2-cbEGF10_model0.pdb"

CONSTRUCT = (807, 951)
CLONING_ARTEFACT = (805, 806)          # non-native residues in 2W86; not fibrillin
CBEGF9, HYBRID2, CBEGF10 = (807, 846), (851, 902), (910, 951)
SS_BONDS = [(811, 821), (816, 830), (832, 845), (914, 926), (921, 935), (937, 950)]
# cbEGF9 calcium site, from data/interim/cbegf_sites.tsv (side chains) and
# structural_metrics.tsv (backbone carbonyls measured within coordinating distance)
CA9_SIDECHAIN = {807: "Asp807", 810: "Glu810", 823: "Asn823"}
CA9_BACKBONE = [808, 824, 827]
RENDER = (1500, 1350)

STYLE = """
# ---- shared conventions for every Phase 6 structural panel ----------------------------
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.25
set cartoon_transparency, 0.0
set cartoon_loop_radius, 0.20
set cartoon_tube_radius, 0.28
set stick_radius, 0.11
set sphere_scale, 0.45, elem Ca
set dash_gap, 0.28
set dash_radius, 0.035
set dash_color, grey40
set label_size, 17
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.25
set label_position, (0, 0, 2.2)
set label_font_id, 7
set_color pal_ca,     [%(ca_r).3f, %(ca_g).3f, %(ca_b).3f]
set_color pal_cys,    [%(cys_r).3f, %(cys_g).3f, %(cys_b).3f]
set_color pal_other,  [%(oth_r).3f, %(oth_g).3f, %(oth_b).3f]
set_color pal_af3,    [%(af3_r).3f, %(af3_g).3f, %(af3_b).3f]
set_color pal_exp,    [%(exp_r).3f, %(exp_g).3f, %(exp_b).3f]
set_color pal_outside,[%(out_r).3f, %(out_g).3f, %(out_b).3f]
"""


def hex_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def style_block() -> str:
    c = {}
    for key, name in (("ca_site", "ca"), ("cbegf_cys", "cys"), ("other_cbegf", "oth"),
                      ("af3", "af3"), ("experimental", "exp"), ("outside_cbegf", "out")):
        r, g, b = hex_rgb(fs.PALETTE[key])
        c[f"{name}_r"], c[f"{name}_g"], c[f"{name}_b"] = r, g, b
    return STYLE % c


def crop_margins(png: Path, pad: int = 14) -> tuple[int, int]:
    """
    Trim the white border PyMOL leaves around the molecule.

    Without this the montage is mostly empty page: PyMOL renders into a fixed canvas and the
    subject occupies whatever fraction of it the camera happens to give. Cropping to the
    rendered content means each panel fills its cell regardless of molecule shape.
    """
    from PIL import Image
    import numpy as np

    im = Image.open(png).convert("RGB")
    a = np.asarray(im)
    ink = (a < 248).any(axis=2)
    if not ink.any():
        raise SystemExit(f"{png.name} rendered blank — nothing was drawn")
    ys, xs = np.where(ink)
    box = (max(int(xs.min()) - pad, 0), max(int(ys.min()) - pad, 0),
           min(int(xs.max()) + pad + 1, a.shape[1]), min(int(ys.max()) + pad + 1, a.shape[0]))
    im.crop(box).save(png)
    return box[2] - box[0], box[3] - box[1]


def run_pymol(pml: Path, log: list[str]) -> None:
    out = subprocess.run(["pymol", "-cq", str(pml)], capture_output=True, text=True, timeout=900)
    for line in (out.stdout + out.stderr).splitlines():
        if line.strip():
            log.append(f"    [pymol] {line.rstrip()}")
    if out.returncode != 0:
        raise SystemExit(f"PyMOL failed on {pml.name} (exit {out.returncode}):\n{out.stderr}")


def write_and_run(name: str, body: str, log: list[str]) -> Path:
    png = PANEL_DIR / f"{name}.png"
    pml = PML_DIR / f"{name}.pml"
    pml.write_text(style_block() + body + f"""
ray {RENDER[0]}, {RENDER[1]}
png {png}, dpi=300
print("VIEW {name}", cmd.get_view())
""", encoding="utf-8")
    log.append(f"  rendering {name}")
    run_pymol(pml, log)
    if not png.is_file():
        raise SystemExit(f"{png} was not produced")
    w, h = crop_margins(png)
    log.append(f"    cropped to {w}x{h} px")
    return png


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    log = [f"Phase 6 — Figure 4 — seed {fs.SEED}", f"PyMOL {fs.pymol_version()}"]
    num = fs.Numbers(FIG, Path(__file__).name)
    PML_DIR.mkdir(parents=True, exist_ok=True)
    PANEL_DIR.mkdir(parents=True, exist_ok=True)

    for p in (EXP, AF3):
        if not p.is_file():
            raise SystemExit(f"missing input structure {p}")

    df = fs.load_merged()
    df["mech"] = fs.mech_group(df)
    sub = df[df.template_pdb == "2W86"]
    path = sub[sub.set_primary == "pathogenic"]
    by_mech = {m: sorted(int(x) for x in path[path.mech == m].position.unique())
               for m in fs.MECH_ORDER}
    other_path = sorted(int(x) for x in path[path.mech.isna()].position.unique())
    num.set("template_pdb", "2W86")
    num.set("n_variants_on_template", len(sub))
    num.set("n_pathogenic_on_template", len(path))
    for m in fs.MECH_ORDER:
        num.set(f"pathogenic_positions_{m}", by_mech[m])
    num.set("pathogenic_positions_outside_cbegf", other_path)
    log.append(f"2W86 carries {len(sub)} covered variants, {len(path)} pathogenic")
    log.append(f"  by mechanism: { {m: len(v) for m, v in by_mech.items()} }, "
               f"outside cbEGF: {len(other_path)}")

    calib = pd.read_csv(fs.RESULTS / "phase4_calibration.tsv", sep="\t")
    c9 = calib[(calib.domain == "cbEGF9") & (calib.model == 0)].iloc[0]
    num.set("af3_ca_deviation_cbEGF9_A", float(c9.ca_deviation_A))
    num.set("af3_backbone_rmsd_cbEGF9_A", float(c9.domain_backbone_rmsd_A))
    log.append(f"cbEGF9 AF3 model 0: Ca deviation {c9.ca_deviation_A} A, "
               f"backbone RMSD {c9.domain_backbone_rmsd_A} A")

    ss_sel = "+".join(str(r) for pair in SS_BONDS for r in pair)
    common = f"""
load {EXP}, exp
remove exp and resi {CLONING_ARTEFACT[0]}-{CLONING_ARTEFACT[1]}
hide everything
"""

    # ---- A: the whole fragment ------------------------------------------------------
    # The two cbEGF domains take the same light tone and the hybrid a darker one, so the
    # eye reads "calcium-binding module, spacer, calcium-binding module" at a glance.
    pA = write_and_run("fig4_A_fragment", common + f"""
show cartoon, exp and polymer
color grey85, exp and polymer
color skyblue, exp and resi {CBEGF9[0]}-{CBEGF9[1]}
color grey50,  exp and resi {HYBRID2[0]}-{HYBRID2[1]}
color skyblue, exp and resi {CBEGF10[0]}-{CBEGF10[1]}
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
set sphere_scale, 0.55, elem Ca
show sticks, exp and resi {ss_sel} and (sidechain or name CA)
color pal_cys, exp and resi {ss_sel} and elem C
color yellow,  exp and resi {ss_sel} and elem S
set stick_radius, 0.18, exp and resi {ss_sel}
orient exp and polymer
turn x, -12
zoom exp and polymer, 1.0
""", log)

    # ---- B: the calcium site --------------------------------------------------------
    sc_sel = "+".join(str(r) for r in CA9_SIDECHAIN)
    bb_sel = "+".join(str(r) for r in CA9_BACKBONE)
    dist_cmds = "\n".join(
        f"distance d{r}, (exp and resi {r} and sidechain and (elem O+N)), (ca9), 3.2, 0"
        for r in CA9_SIDECHAIN) + "\n" + "\n".join(
        f"distance b{r}, (exp and resi {r} and name O), (ca9), 3.4, 0"
        for r in CA9_BACKBONE)
    labels = "\n".join(f'label exp and resi {r} and name CB, "{n}"'
                       for r, n in CA9_SIDECHAIN.items())
    # The ion is drawn small on purpose. At its ionic radius it hides the very side chains
    # the panel exists to show, and an oversized sphere would also misrepresent how much
    # room the metal occupies in the site.
    pB = write_and_run("fig4_B_casite", common + f"""
select ca9, exp and elem Ca within 8 of (exp and resi {CBEGF9[0]}-{CBEGF9[1]})
show cartoon, exp and resi {CBEGF9[0]}-{CBEGF9[1]}
color grey85, exp and resi {CBEGF9[0]}-{CBEGF9[1]}
set cartoon_transparency, 0.6, exp
show spheres, ca9
color pal_ca, ca9
set sphere_scale, 0.28, ca9
show sticks, exp and resi {sc_sel} and (sidechain or name CA)
color pal_ca, exp and resi {sc_sel} and elem C
show sticks, exp and resi {bb_sel} and name C+O+CA
color grey45, exp and resi {bb_sel} and elem C
{dist_cmds}
hide labels, d* b*
{labels}
orient (exp and resi {sc_sel}+{bb_sel}) or ca9
turn y, 20
zoom ((exp and resi {sc_sel}+{bb_sel}) or ca9), 2.2
""", log)

    # ---- C: AF3 against the crystal -------------------------------------------------
    pC = write_and_run("fig4_C_af3_overlay", common + f"""
load {AF3}, af3
hide everything
select exp9, exp and resi {CBEGF9[0]}-{CBEGF9[1]}
select af39, af3 and resi {CBEGF9[0]}-{CBEGF9[1]}
align af39 and name CA+C+N+O, exp9 and name CA+C+N+O
show cartoon, exp9
show cartoon, af39
color pal_exp, exp9
color pal_af3, af39
select caE, exp and elem Ca within 8 of exp9
select caA, af3 and elem Ca within 8 of af39
show spheres, caE
show spheres, caA
color pal_exp, caE
color pal_af3, caA
# Deliberately unequal radii, and the experimental ion is translucent. The two ions are
# 0.13 A apart, so equal opaque spheres would occlude one another and the panel would look
# like it contains a single ion. A small solid blue ball inside a larger translucent grey
# shell reads as "the model put it in the same place".
set sphere_scale, 0.45, caE
set sphere_scale, 0.18, caA
set sphere_transparency, 0.55, caE
print("ATOMS caE", cmd.count_atoms("caE"), "caA", cmd.count_atoms("caA"))
orient exp9
turn y, 20
zoom exp9, 1.5
""", log)

    # ---- D: pathogenic variants -----------------------------------------------------
    def sel(positions):
        return "+".join(str(p) for p in positions) if positions else "0"

    pD = write_and_run("fig4_D_variants", common + f"""
show cartoon, exp and polymer
color grey85, exp and polymer
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
set sphere_scale, 0.35, elem Ca
show spheres, exp and resi {sel(by_mech['cys_removing'])} and name CA
color pal_cys, exp and resi {sel(by_mech['cys_removing'])} and name CA
show spheres, exp and resi {sel(by_mech['ca_ligand'])} and name CA
color pal_ca,  exp and resi {sel(by_mech['ca_ligand'])} and name CA
show spheres, exp and resi {sel(other_path)} and name CA
color pal_outside, exp and resi {sel(other_path)} and name CA
set sphere_scale, 0.85, name CA
orient exp and polymer
turn x, -12
zoom exp and polymer, 1.0
""", log)

    # ---- montage --------------------------------------------------------------------
    n_path_pos = int(path.position.nunique())
    num.set("n_pathogenic_positions_on_template", n_path_pos)

    fig, axs = plt.subplots(2, 2, figsize=(7.4, 7.2))
    fig.subplots_adjust(hspace=0.24, wspace=0.06, top=0.91, bottom=0.10,
                        left=0.02, right=0.99)
    titles = [
        (pA, "A", "cbEGF9 – hybrid 2 – cbEGF10 (PDB 2W86, 1.8 Å)",
         "blue = the two cbEGF domains, grey = hybrid 2 (UniProt “TB 4”);\n"
         "purple/yellow sticks = the six cbEGF disulfides; green = Ca²⁺"),
        (pB, "B", "The cbEGF9 calcium site",
         "three side-chain ligands (labelled) and three backbone carbonyls (grey),\n"
         "dashed to the ion at 2.3–2.6 Å; ion drawn under-sized so the ligands stay visible"),
        (pC, "C", "AF3 cofolded model on the crystal structure",
         f"cbEGF9 only; backbone RMSD {c9.domain_backbone_rmsd_A:.2f} Å and the two Ca²⁺\n"
         f"ions differ by {c9.ca_deviation_A:.2f} Å — they overlap at this scale"),
        (pD, "D", f"{len(path)} pathogenic variants at {n_path_pos} positions",
         "Cα drawn as a sphere, coloured by mechanism; positions, not variants:\n"
         f"{len(by_mech['cys_removing'])} cysteine-removing, "
         f"{len(by_mech['ca_ligand'])} calcium-ligand, {len(other_path)} outside a cbEGF"),
    ]
    for ax, (png, letter, title, sub_t) in zip(axs.ravel(), titles):
        ax.imshow(plt.imread(png))
        ax.set_axis_off()
        # The panel letter rides in the title string. Placing it in axes coordinates
        # collides with the title, because an image axis is cropped tight to its content.
        ax.set_title(rf"$\bf{{{letter}}}$   {title}", fontsize=8.2, loc="left", pad=4)
        ax.text(0, -0.02, sub_t, transform=ax.transAxes, fontsize=6.4, va="top",
                ha="left", color=fs.PALETTE["ink2"], linespacing=1.5)

    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc=fs.PALETTE["ca_site"],
               mec="none", label="Ca²⁺ ion / calcium-ligand variant"),
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc=fs.PALETTE["cbegf_cys"],
               mec="none", label="cbEGF disulfide / cysteine-removing variant"),
        Line2D([0], [0], marker="o", ls="none", ms=6, mfc=fs.PALETTE["outside_cbegf"],
               mec="none", label="pathogenic variant outside a cbEGF domain"),
        Patch(facecolor=fs.PALETTE["af3"], label="AF3 cofolded model (panel C)"),
        Patch(facecolor=fs.PALETTE["experimental"], label="experimental structure (panel C)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=6.8,
               handletextpad=0.4, columnspacing=1.4, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle("Calcium, disulfides and pathogenic variants in a solved FBN1 fragment",
                 fontsize=10, x=0.02, ha="left", y=0.985)

    fs.save(fig, FIG, num)
    plt.close(fig)
    print(f"  log -> {fs.write_log('16_fig4', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
