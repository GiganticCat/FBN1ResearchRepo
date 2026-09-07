#!/usr/bin/env python3
"""
paperstyle.py — figure conventions taken from the journals this work would be submitted to.

Replaces `figstyle.py` for the final figure set. The earlier module encoded three rules that the
reference papers in `resources/papers/` do not follow, and one habit none of them has.

  PANEL TITLES ARE NORMAL. `figstyle.apply_style()` banned them outright and the v2 gate enforced
  the ban. Cheng 2023 (Science) titles nearly every panel -- "ClinVar (Class-balanced 18924
  variants)", "ACMG genes"; Godwin 2023 (Nat Struct Mol Biol) titles panels with the measured
  constant, "PF3 Kd = 17.690 +/- 5.638 nM". What those titles never do is state a conclusion.
  A title says what is plotted; the caption says what it means. That is the rule here.

  COMPARISONS ARE ANNOTATED ON THE DATA. Godwin 2023 draws a bracket and "****" directly over the
  two groups being compared. The old Figure 2 instead gave a whole panel to Cliff's delta with
  bootstrap intervals -- a quantity no reader of a structural-biology paper reads, occupying a
  quarter of the figure to restate what the two panels beside it already showed. Effect sizes
  belong in the statistics table. Brackets and P values belong on the plot.

  SERIES ARE LABELLED WHERE THEY ARE DRAWN. Godwin labels its two traces "Fibrillin" and
  "Fibrillin + LTBP-1" inside the axes, in the colour of the trace. A legend box in a corner is a
  lookup table the reader has to hold in their head.

  NOTHING IS ROTATED. Rotated 30-degree category labels are a matplotlib default and appear in
  none of the reference papers. Categories get short horizontal labels, wrapped to two lines when
  they must be.

Widths are Nature's: 89 mm single column, 120 mm 1.5, 183 mm double. Type is 7 pt for labels and
6 pt for ticks, which is what those pages measure at.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
INTERIM = ROOT / "data/interim"
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"

SEED = 20260830

# Nature column widths, inches.
W1 = 3.50      #  89 mm
W15 = 4.72     # 120 mm
W2 = 7.20      # 183 mm

# ---------------------------------------------------------------------------------------
# Colour
#
# Measured off the pages in `resources/papers/`, not chosen from a palette generator. The data
# panels in Godwin 2023 (Nat Struct Mol Biol) and Jensen 2009 (Structure) are almost entirely
# black, grey and white; colour appears only where it carries information a position or a shape
# cannot -- a domain type in a schematic, one highlighted series. An earlier version of this
# module gave every group its own saturated hue, which is what a plotting library does by
# default and what no journal page does.
#
# The rule here: in a panel where groups already occupy separate x positions, colour is
# redundant and is not used. Groups are distinguished by fill -- solid black, mid grey, open --
# in that order.
#
# THE FIGURES ARE MONOCHROME. An earlier version kept two accents, blue for calcium and red for
# the highlighted group. They are gone: the whole set is now black, grey and white, which is what
# the data panels of Godwin 2023 and Jensen 2009 actually are, and what survives a photocopier,
# a greyscale proof and every form of colour-vision deficiency without a separate check.
#
# Removing hue costs nothing here because nothing was ever encoded by hue alone. Where two marks
# must be told apart they are told apart by a property that survives desaturation, in this order
# of preference:
#
#   1. POSITION   -- separate x or y, which most panels already use
#   2. SHAPE      -- circle vs diamond; this is how calcium donors separate from cysteines in
#                    Fig. 1b, where both would otherwise be black marks on one schematic
#   3. FILL       -- solid vs open, which also carries side-chain vs backbone donor
#   4. LIGHTNESS  -- black / mid grey / light grey, in bar fills where three categories share
#                    one series and shape is not available
#
# `ca` and `flag` survive as names because they carry meaning the figures still need -- calcium
# is black everywhere it appears, and the flagged group is ringed rather than recoloured -- but
# both now resolve to ink.
# ---------------------------------------------------------------------------------------
C = {
    "ink":     "#000000",
    "ink2":    "#3d3d3d",
    "grey":    "#8a8a8a",   # second series
    "grey_lt": "#c4c4c4",   # third series / box fills
    "rule":    "#b8b8b8",
    "band":    "#e8e8e8",   # "no effect possible" shading
    "muted":   "#6e6e6e",
    "surface": "#ffffff",
    # calcium, and the highlighted group -- both ink, separated by shape and fill
    "ca":      "#000000",   # calcium, everywhere it appears
    "ca_mid":  "#767676",   # calcium where it shares a bar series with black cysteine bars
    "flag":    "#000000",   # the one thing to look at: ringed, not recoloured
    # domain architecture: lightness steps, as in the reference papers' schematics
    "dom_cbegf": "#2f2f2f",
    "dom_egf":   "#8f8f8f",
    "dom_tb":    "#c9c9c9",
    "dom_other": "#eaeaea",
}

# Group order, fills and labels.
#
# The labels name each class in the vocabulary the cbEGF literature already uses for it. Jensen
# 2009 calls a coordinating residue a "Ca2+ ligand" and describes a site as three side-chain
# oxygens and three main-chain carbonyls, so that is the term used here. An earlier version
# labelled the groups by how they were computed -- "Matches motif, no Ca2+ contact",
# "Coordinates Ca2+" -- which is a description of the pipeline rather than a name for the class,
# and read as jargon to a reader who knows the field. See results/style_reference.md.
GROUPS = ["ca_ligand", "composition_control", "motif_overcalled", "cys_removing"]
GFILL = {"ca_ligand": C["ink"], "composition_control": C["grey"],
         "motif_overcalled": C["grey_lt"], "cys_removing": C["ink"]}
GEDGE = {"ca_ligand": C["ink"], "composition_control": C["grey"],
         "motif_overcalled": C["ink2"], "cys_removing": C["ink"]}
GLABEL = {
    "ca_ligand": "Ca²⁺ ligand\n(side chain)",
    "composition_control": "Non-ligand\nAsp/Glu/Asn",
    "motif_overcalled": "Consensus motif,\nnon-ligand",
    "cys_removing": "Cysteine\nsubstitution",
}
# kept for callers that still ask for a single colour
GCOLOR = GFILL


def apply() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.015,
        "savefig.facecolor": C["surface"],
        "figure.facecolor": C["surface"],
        "axes.facecolor": C["surface"],
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.2,
        "legend.fontsize": 6,
        "axes.edgecolor": C["ink"],
        "axes.labelcolor": C["ink"],
        "text.color": C["ink"],
        "xtick.color": C["ink"],
        "ytick.color": C["ink"],
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.titlelocation": "center",
        "axes.titlepad": 3.5,
        "legend.frameon": False,
        "lines.linewidth": 1.0,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax, letter: str, dx: float = -0.16, dy: float = 1.06) -> None:
    """Bold lowercase letter at the panel's top-left, Nature style."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=C["ink"])


def title(ax, text: str) -> None:
    """A panel title states what is plotted. It never states a conclusion."""
    ax.set_title(text, fontsize=6.8, color=C["ink"], pad=3.5)


def pstar(p: float) -> str:
    """The convention the reference papers use: symbols for the eye, the value in the caption."""
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# ---------------------------------------------------------------------------------------
# P values, and the q value beside them
#
# Every headline panel now prints the Benjamini-Hochberg q as well as the raw P, because a
# reviewer cannot check a correction that lives only in a supplementary file. The correction is
# applied within a family of tests in `44_statistics_final.py`, and `results/statistics_final.tsv`
# is the single place both values are stored, so the figures read them from there instead of
# recomputing them. Two comparisons change status under the correction -- FoldX at the calcium
# ligands (P = 0.032, q = 0.060) and the Rosetta motif non-ligand row -- and printing only the
# raw value would hide that.
# ---------------------------------------------------------------------------------------

_STATS: pd.DataFrame | None = None

SUP = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def stat(test: str) -> tuple[float, float]:
    """The raw P and the Benjamini-Hochberg q of one named test, from the final table."""
    global _STATS
    if _STATS is None:
        _STATS = pd.read_csv(RESULTS / "statistics_final.tsv", sep="\t")
    row = _STATS[_STATS.test == test]
    if row.empty:
        raise SystemExit(f"no statistics row for {test!r}")
    return float(row.iloc[0].p_raw), float(row.iloc[0].p_adj)


def _one(sym: str, v: float) -> str:
    if v >= 0.05:
        return f"{sym} = {v:.2f}"
    if v >= 1e-4:
        return f"{sym} = {v:.1g}"
    e = int(np.floor(np.log10(v)))
    m = v / 10.0 ** e
    return f"{sym} = {m:.0f}\u00d710{str(e).translate(SUP)}"


def fmt_pq(p: float, q: float, ns: bool = True) -> str:
    """P over its q, the two-line form every headline comparison carries."""
    tail = ", n.s." if (ns and q >= 0.05) else ""
    return f"{_one('P', p)}\n{_one('q', q)}{tail}"


def pq(test: str, ns: bool = True) -> str:
    """`fmt_pq` straight off the named row of results/statistics_final.tsv."""
    return fmt_pq(*stat(test), ns=ns)


def agrees(p: float, test: str, rtol: float = 0.06) -> float:
    """
    Guard for a panel that recomputes a P value it could have read.

    Three panels compute their own Fisher or Mann-Whitney P from the geometry tables. Those
    values must be the ones the statistics table corrected, or the q printed beside them belongs
    to a different test. This returns the table's q and stops if the two P values have drifted.
    """
    p_ref, q_ref = stat(test)
    if not np.isclose(p, p_ref, rtol=rtol):
        raise SystemExit(f"panel P {p:.4g} does not match {test!r} at {p_ref:.4g}")
    return q_ref


def bracket(ax, x1: float, x2: float, y: float, label: str, dy: float = 0.0,
            lw: float = 0.6, fontsize: float = 6.2, color: str | None = None) -> None:
    """A significance bracket over two groups, with the tick marks turned down at the ends."""
    color = color or C["ink"]
    h = dy or (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.022
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color=color,
            clip_on=False, solid_joinstyle="miter")
    ax.text((x1 + x2) / 2, y + h, label, ha="center", va="bottom",
            fontsize=fontsize, color=color, clip_on=False)


def strip(ax, groups, ylabel, rng, log=False, ms=1.9, alpha=0.42,
          show_n=True, width=0.30):
    """
    Individual points with a median bar and an interquartile box.

    Every point is drawn, as Godwin 2023 does for its mass-mapping data, because a box alone
    hides the shape of a distribution and these distributions are the argument. Returns the
    medians so the caller can hand the same numbers to the sidecar rather than recompute them.
    """
    out = {}
    for i, (label, vals, color) in enumerate(groups):
        v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
        out[label] = {"n": int(v.size), "median": float(np.median(v)) if v.size else float("nan")}
        if not v.size:
            continue
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        # an open group is drawn as a ring; a filled one as a disc. Fill, not hue, separates
        # series here, exactly as it does on the reference pages.
        face = "none" if color == "open" else color
        edge = C["ink"] if color == "open" else color
        ax.scatter(np.full(v.size, i) + rng.uniform(-width * 0.42, width * 0.42, v.size), v,
                   s=ms, facecolors=face, edgecolors=edge,
                   linewidths=0.0 if face != "none" else 0.35,
                   alpha=alpha, zorder=2, rasterized=True)
        ax.add_patch(plt_rect(i - width / 2, q1, width, q3 - q1, edge))
        ax.plot([i - width * 0.62, i + width * 0.62], [med, med], lw=1.5, color=edge, zorder=4,
                solid_capstyle="butt")
        if show_n:
            # inside the axes: above it the label fights the panel title, which the earlier
            # version of this figure lost
            ax.annotate(f"n = {v.size:,}", xy=(i, 0.995), xycoords=("data", "axes fraction"),
                        xytext=(0, -1), textcoords="offset points", ha="center", va="top",
                        fontsize=5.8, color=C["ink2"])
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlim(-0.62, len(groups) - 0.38)
    return out


def striph(ax, groups, xlabel, rng, ms=1.9, alpha=0.42, height=0.30, show_n=True):
    """
    The same distribution panel, laid on its side.

    Categories go on the y axis. This is not a stylistic preference: four descriptive category
    names do not fit side by side under a half-column panel, and the alternative -- rotating them
    30 degrees, or abbreviating them into jargon -- is what the reference papers avoid by putting
    the categories on the vertical axis instead (Cheng 2023, Fig. 2a-c). Labels then read
    horizontally with as much room as they need.
    """
    out = {}
    n_g = len(groups)
    for i, (label, vals, color) in enumerate(groups):
        y = n_g - 1 - i                      # first group at the top, as a reader expects
        v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
        out[label] = {"n": int(v.size), "median": float(np.median(v)) if v.size else float("nan")}
        if not v.size:
            continue
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        face = "none" if color == "open" else color
        edge = C["ink"] if color == "open" else color
        ax.scatter(v, np.full(v.size, y) + rng.uniform(-height * 0.42, height * 0.42, v.size),
                   s=ms, facecolors=face, edgecolors=edge,
                   linewidths=0.0 if face != "none" else 0.35,
                   alpha=alpha, zorder=2, rasterized=True)
        ax.add_patch(plt_rect(q1, y - height / 2, q3 - q1, height, edge))
        ax.plot([med, med], [y - height * 0.62, y + height * 0.62], lw=1.5, color=edge,
                zorder=4, solid_capstyle="butt")
    ax.set_yticks(range(n_g))
    # n goes in the tick label, where the reference pages put it. Anywhere inside the axes it
    # collides with either the data or the significance brackets.
    # n is appended to the label's last line rather than added as a third line: at these row
    # spacings a third line runs into the row above
    ax.set_yticklabels([f"{g[0]} (n = {out[g[0]]['n']:,})" if show_n else g[0]
                        for g in reversed(groups)])
    ax.set_xlabel(xlabel)
    ax.set_ylim(-0.62, n_g - 0.38)
    return out


def clip_note_h(ax, groups, hi):
    """How many points a clipped x axis is hiding, per row."""
    n_g = len(groups)
    for i, (_, vals, color) in enumerate(groups):
        v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
        n = int((v > hi).sum())
        if n:
            ax.annotate(f"{n} →", xy=(hi, n_g - 1 - i), xytext=(-1, -6.5),
                        textcoords="offset points", ha="right", va="top",
                        fontsize=5.6, color=C["ink2"])


def bracket_h(ax, y1, y2, x, label, dx=0.0, lw=0.6, fontsize=6.2, color=None, ly=None):
    """
    A significance bracket beside two rows.

    `ly` moves the label off the midpoint. The two-line P/q labels are tall enough that a
    midpoint landing on a category row collides with that row's clipped-point count, and nudging
    the text is less disruptive than moving either the bracket or the count.
    """
    color = color or C["ink"]
    w = dx or (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.022
    ax.plot([x, x + w, x + w, x], [y1, y1, y2, y2], lw=lw, color=color, clip_on=False)
    ax.text(x + w, (y1 + y2) / 2 if ly is None else ly, " " + label, ha="left", va="center",
            fontsize=fontsize, color=color, clip_on=False)


def plt_rect(x, y, w, h, color):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, facecolor="none", edgecolor=color,
                     linewidth=0.55, zorder=3)


def clip_note(ax, groups, hi):
    """Say how many points a clipped axis is hiding, per group. Never hide them silently.

    Placed one line below the n, because both want the top of the axes and the n wins.
    """
    for i, (_, vals, color) in enumerate(groups):
        v = np.asarray([x for x in vals if x is not None and np.isfinite(x)], float)
        n = int((v > hi).sum())
        if n:
            ax.annotate(f"{n} ↑", xy=(i, 0.995), xycoords=("data", "axes fraction"),
                        xytext=(0, -9), textcoords="offset points",
                        ha="center", va="top", fontsize=5.6, color=color)


def label_at(ax, x, y, text, color, **kw):
    """Series label in the plot area, in the colour of the series it names."""
    ax.text(x, y, text, color=color, fontsize=6.4, **kw)


# ---------------------------------------------------------------------------------------
# Data — the joins are figstyle's, reused rather than reimplemented so the two cannot disagree
# ---------------------------------------------------------------------------------------

def _fs():
    import importlib.util
    spec = importlib.util.spec_from_file_location("_figstyle", ROOT / "scripts" / "figstyle.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load() -> pd.DataFrame:
    fs = _fs()
    df = fs.load_v2()
    df["group"] = fs.mech_group_v2(df)
    return df


def domains() -> pd.DataFrame:
    return pd.read_csv(INTERIM / "uniprot_domains.tsv", sep="\t")


class Numbers:
    """Every value a figure prints, written beside it so the gate can re-derive each one."""

    def __init__(self, figure: str, script: str):
        self.figure, self.script, self.values = figure, script, {}

    def set(self, key, value):
        if isinstance(value, (np.integer,)):
            value = int(value)
        elif isinstance(value, (np.floating,)):
            value = float(value)
        elif isinstance(value, (np.bool_,)):
            value = bool(value)
        self.values[key] = value
        return value

    def write(self) -> Path:
        p = FIGURES / f"{self.figure}.numbers.json"
        p.write_text(json.dumps({
            "figure": self.figure, "script": self.script,
            "generated": datetime.now(timezone.utc).isoformat(), "seed": SEED,
            "values": self.values}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return p


def save(fig, name: str, numbers: Numbers | None = None) -> None:
    FIGURES.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{ext}")
    print(f"  wrote figures/{name}.png ({(FIGURES / f'{name}.png').stat().st_size/1024:.0f} kB)"
          f" and .pdf")
    if numbers is not None:
        print(f"  wrote {numbers.write().relative_to(ROOT)}")


def write_log(stem: str, lines: list[str]) -> Path:
    LOGS.mkdir(exist_ok=True)
    p = LOGS / f"{stem}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.log"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p
