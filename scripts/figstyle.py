#!/usr/bin/env python3
"""
figstyle.py — PHASE 6. Shared figure conventions.

Imported by every `1[3-8]_fig*.py` script so that all panels in the manuscript speak one
visual language. Nothing here reads or writes analysis data on its own; it supplies

  * PALETTE      — the fixed entity -> colour map, validated for colour-vision deficiency
  * apply_style  — matplotlib rcParams for print
  * load_merged  — the Phase 2/3/5 tables joined the *same* way `12_statistics.py` joins them
  * mech_group   — the three mechanism classes, defined byte-for-byte as in Phase 5
  * Numbers      — records every value a figure draws, for the Phase 6 gate to re-derive

WHY the Numbers recorder exists: the Phase 6 gate has to check that "numbers in figures match
data/processed/". A PNG is opaque, so each figure script declares the values it annotated into
`figures/<name>.numbers.json`, and `tests/test_figures.py` recomputes them from the processed
tables and asserts equality. A figure that quietly drifts from the data fails the gate.

Palette provenance — each set was checked with the data-viz validator
(`validate_palette.js --mode light --pairs all`) against the print surface:

  clinical  #e34948 / #eda100 / #2a78d6   ALL PASS (worst CVD dE 15.3, normal 20.8)
  site      #1baf7a / #4a3aa7             ALL PASS (worst CVD dE 31.1, normal 35.8)
  method    #2a78d6 / #eb6834             ALL PASS (worst CVD dE 24.7, normal 33.6)

`#eda100` (VUS) and `#1baf7a` (calcium) fall below 3:1 contrast on white, which the validator
flags as needing relief. Relief is supplied everywhere they appear: every distribution panel
carries a direct n= and median label, so identity is never colour-alone.

`other_cbegf` is deliberately neutral grey. It is the reference series the other two are read
against, not a fourth identity — the same role a baseline plays. It is therefore excluded from
the categorical validation above rather than forced to carry a hue.
"""

from __future__ import annotations

import json
import subprocess
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

SEED = 20260802

# --------------------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------------------

PALETTE = {
    # clinical class (ClinVar, primary >=2-star tier)
    "pathogenic": "#e34948",
    "vus": "#eda100",
    "benign": "#2a78d6",
    # mechanism / site class
    "ca_site": "#1baf7a",
    "cbegf_cys": "#4a3aa7",
    "other_cbegf": "#7a7a74",
    # structure provenance
    "af3": "#2a78d6",
    "transplant": "#eb6834",
    "experimental": "#52514e",
    # pathogenic variants that fall outside any cbEGF domain — a residual bucket, not an
    # identity, so it stays neutral. Exactly PyMOL's grey40, which panel D of Figure 4
    # renders with, so the legend swatch and the render cannot drift apart.
    "outside_cbegf": "#666666",
    # ink and furniture
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#8a8a85",
    "grid": "#e3e2dd",
    "surface": "#ffffff",
    # domain architecture track: lightness steps, not hues, so the track never competes
    # with the data palette for attention
    "dom_cbEGF": "#4a4a46",
    "dom_EGF": "#9a9a94",
    "dom_TB": "#c9c8c2",
    "dom_other": "#ecebe6",
}

CLINICAL_ORDER = ["pathogenic", "vus", "benign"]
CLINICAL_LABEL = {
    "pathogenic": "Pathogenic / likely pathogenic",
    "vus": "Uncertain significance",
    "benign": "Benign / likely benign",
}

MECH_ORDER = ["ca_ligand", "cys_removing", "other_cbegf"]
MECH_COLOR = {
    "ca_ligand": PALETTE["ca_site"],
    "cys_removing": PALETTE["cbegf_cys"],
    "other_cbegf": PALETTE["other_cbegf"],
}
MECH_LABEL = {
    "ca_ligand": "Direct Ca²⁺ ligand",
    "cys_removing": "cbEGF cysteine-removing",
    "other_cbegf": "Other cbEGF residue",
}


def apply_style() -> None:
    """
    Journal-plate defaults, matched to the conventions in `resources/papers/`.

    Checked against Godwin 2023 (Nat Struct Mol Biol 30:608) and Jensen 2009 (Structure 17:759),
    both of which sit in `resources/papers/`. What those pages actually do, and what this
    reproduces:

      * Type is SMALL and uniform -- roughly 7 pt for labels, 6 pt for ticks. Screen-sized type
        is the single clearest sign of a figure that has not been prepared for print.
      * Panels carry a lowercase bold letter in the top-left corner and NO title. The
        interpretation lives in the caption, never on the axes. An earlier version of these
        figures put sentences like "A force field that barely sees the calcium" above each panel;
        no journal does that.
      * Furniture is hairline (0.5 pt spines and ticks), ticks point outward and are short.
      * Axis labels are terse and carry units in parentheses.
      * Panels are packed tightly; whitespace is a layout failure, not breathing room.

    Nothing here is decorative: every value is set so a 90 mm single-column or 180 mm
    double-column placement stays legible at print size.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,          # journals ask for 600 dpi on line/combination art
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "savefig.facecolor": PALETTE["surface"],
        "figure.facecolor": PALETTE["surface"],
        "axes.facecolor": PALETTE["surface"],
        # DejaVu Sans is the only sans face present in this venv (checked in Phase 6);
        # Arial/Helvetica are absent, so pinning them would silently fall back and make
        # renders machine-dependent. Kept explicit for reproducibility.
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 7,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.edgecolor": PALETTE["ink"],
        "axes.labelcolor": PALETTE["ink"],
        "text.color": PALETTE["ink"],
        "xtick.color": PALETTE["ink"],
        "ytick.color": PALETTE["ink"],
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.size": 2.2,
        "ytick.major.size": 2.2,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": PALETTE["grid"],
        "grid.linewidth": 0.4,
        "legend.frameon": False,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.35,
        "legend.borderpad": 0.2,
        "lines.linewidth": 0.9,
        "patch.linewidth": 0.5,
        "pdf.fonttype": 42,   # editable text in the vector output, not outlines
        "ps.fonttype": 42,
    })


# Journal column widths in inches. Figures are sized to these rather than to whatever looked
# right on screen, because a figure scaled down to fit a column is a figure with unreadable type.
COL_1 = 3.46      #  88 mm, single column
COL_15 = 5.51     # 140 mm, 1.5 column
COL_2 = 7.09      # 180 mm, double column


def panel_label(ax, letter: str, dx: float = -0.115, dy: float = 1.02) -> None:
    """Lowercase bold letter, top-left, Nature style. No trailing punctuation."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=PALETTE["ink"])


def tidy(ax, grid_axis: str | None = "y") -> None:
    """Hairline grid behind the data, on one axis only. Called by every panel."""
    if grid_axis:
        getattr(ax, f"{grid_axis}axis").grid(True, zorder=0)
        ax.set_axisbelow(True)


# --------------------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------------------

def tbool(s: pd.Series) -> pd.Series:
    """TSV booleans arrive as a mix of bool, 'True'/'False' strings and NaN. Normalise."""
    return s.map(lambda v: str(v).strip().lower() == "true")


def load_merged() -> pd.DataFrame:
    """
    variant_annotated + structural_metrics + foldx_ddg, joined on variation_id.

    Mirrors the join in `12_statistics.py`: the annotation table is authoritative and the
    structural/FoldX tables contribute only columns it does not already carry.
    """
    ann = pd.read_csv(PROC / "variant_annotated.tsv", sep="\t", low_memory=False)
    smet = pd.read_csv(PROC / "structural_metrics.tsv", sep="\t", low_memory=False)
    ddg = pd.read_csv(PROC / "foldx_ddg.tsv", sep="\t", low_memory=False)

    smet = smet[["variation_id"] + [c for c in smet.columns if c not in ann.columns]]
    ddg = ddg[["variation_id"] + [c for c in ddg.columns if c not in ann.columns
                                  and c not in smet.columns]]
    df = ann.merge(smet, on="variation_id", how="left").merge(ddg, on="variation_id", how="left")

    if len(df) != len(ann):
        raise SystemExit(f"join changed row count {len(ann)} -> {len(df)}; duplicate keys?")

    for c in ("in_cbegf", "in_cbegf_cys", "cys_removing", "cys_creating", "is_ca_consensus",
              "is_ca_sidechain_ligand", "is_critical_residue", "is_cysteine_site",
              "is_direct_ca_ligand", "buried", "structural_coverage", "in_neonatal_region"):
        if c in df.columns:
            df[c] = tbool(df[c])

    df["ddG"] = pd.to_numeric(df.get("ddG_kcal_mol"), errors="coerce")
    df["am"] = pd.to_numeric(df.get("am_pathogenicity"), errors="coerce")
    df["rsa"] = pd.to_numeric(df.get("rsa_pct"), errors="coerce")
    return df


def mech_group(df: pd.DataFrame) -> pd.Series:
    """
    The three mechanism classes, defined exactly as in `12_statistics.py` A4b.

    `ca_ligand` is structure-measured (a side chain actually within coordinating distance of
    the modelled ion), NOT merely "sits at a consensus position" — the distinction CLAUDE.md
    insists on. `other_cbegf` is every remaining cbEGF residue and acts as the control.
    """
    g = pd.Series(pd.NA, index=df.index, dtype="object")
    cys = df["in_cbegf_cys"] & df["cys_removing"]
    lig = df["is_direct_ca_ligand"]
    g[df["in_cbegf"] & ~lig & ~cys] = "other_cbegf"
    g[cys] = "cys_removing"
    g[lig] = "ca_ligand"
    return g


def load_v2() -> pd.DataFrame:
    """
    The v2 join: annotation + AF3-wide structural metrics + both force fields.

    Differs from `load_merged` in what defines a calcium ligand. v1 used the sequence motif,
    which over-calls: 143 positions the motif flags coordinate no calcium in any model, and
    UniProt annotates no calcium site in P35555 at all, so the motif had nothing to check it
    against. v2 uses observed coordination in the AF3 models (`ca_ligand_sidechain`), and
    structural coverage rises from 751 variants on four experimental templates to 3831 across
    all 43 cbEGF domains.

    Rosetta ddG is present for the 1934 variants in the mechanism groups; FoldX only for the
    751 with experimental-template coverage. Both are left as NaN elsewhere rather than filled.
    """
    ann = pd.read_csv(PROC / "variant_annotated.tsv", sep="\t", low_memory=False)
    smet = pd.read_csv(PROC / "structural_metrics_af3.tsv", sep="\t", low_memory=False)
    fx = pd.read_csv(PROC / "foldx_ddg.tsv", sep="\t", low_memory=False)
    ros = pd.read_csv(PROC / "rosetta_ddg.tsv", sep="\t", low_memory=False)

    smet = smet[["variation_id"] + [c for c in smet.columns if c not in ann.columns]]
    fx = fx[["variation_id", "ddG_kcal_mol", "ddG_disulfide"]].rename(
        columns={"ddG_kcal_mol": "foldx_ddG"})
    ros = ros[ros["status"] == "ok"][["variation_id", "ddg_mean", "group"]].rename(
        columns={"ddg_mean": "rosetta_ddG", "group": "rosetta_group"})

    df = (ann.merge(smet, on="variation_id", how="left")
             .merge(fx, on="variation_id", how="left")
             .merge(ros, on="variation_id", how="left"))
    if len(df) != len(ann):
        raise SystemExit(f"join changed row count {len(ann)} -> {len(df)}; duplicate keys?")

    for c in ("ca_ligand_sidechain", "ca_ligand_backbone_only", "ca_ligand_any",
              "consensus_but_not_ligand", "is_cysteine", "cys_removing", "buried",
              "in_cbegf", "in_cbegf_cys"):
        if c in df.columns:
            df[c] = tbool(df[c])

    df["foldx_ddG"] = pd.to_numeric(df["foldx_ddG"], errors="coerce")
    df["rosetta_ddG"] = pd.to_numeric(df["rosetta_ddG"], errors="coerce")
    df["am"] = pd.to_numeric(df.get("am_pathogenicity"), errors="coerce")
    df["rsa"] = pd.to_numeric(df.get("rsa_pct"), errors="coerce")
    df["has_structure"] = df["construct"].notna()
    return df


def mech_group_v2(df: pd.DataFrame) -> pd.Series:
    """
    Four groups, matching `28_statistics_v2.py`.

    The composition control exists because calcium ligands are D/E/N by definition, and that
    alone shifts any comparison against a mixed background. Comparing ligands to *other* D/E/N
    residues is what isolates the calcium role from the amino-acid identity.
    """
    g = pd.Series(pd.NA, index=df.index, dtype="object")
    den = df["wt_aa"].isin(["D", "E", "N"])
    g[df["has_structure"] & den & ~df["cys_removing"]] = "composition_control"
    g[df["consensus_but_not_ligand"]] = "motif_overcalled"
    g[df["cys_removing"]] = "cys_removing"
    g[df["ca_ligand_sidechain"]] = "ca_ligand"
    return g


MECH2_ORDER = ["ca_ligand", "composition_control", "motif_overcalled", "cys_removing"]
MECH2_COLOR = {
    "ca_ligand": PALETTE["ca_site"],
    "composition_control": PALETTE["other_cbegf"],
    "motif_overcalled": PALETTE["transplant"],
    "cys_removing": PALETTE["cbegf_cys"],
}
# Two label sets, deliberately. The figures carry standard nomenclature, because a journal
# figure that says "Touches the calcium" on its axis does not read as a journal figure. The
# plain-language wording moved to the accompanying document, where it belongs and where there is
# room to define terms. MECH2_PLAIN is what `results/v2_report.html` uses.
MECH2_LABEL = {
    "ca_ligand": "Ca²⁺ ligand",
    "composition_control": "D/E/N control",
    "motif_overcalled": "Motif non-ligand",
    "cys_removing": "cbEGF Cys loss",
}
MECH2_PLAIN = {
    "ca_ligand": "Touches the calcium",
    "composition_control": "Same letter, no calcium",
    "motif_overcalled": "Looks right, isn't",
    "cys_removing": "Breaks a cysteine bridge",
}


def domains() -> pd.DataFrame:
    return pd.read_csv(INTERIM / "uniprot_domains.tsv", sep="\t")


# --------------------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------------------

class Numbers:
    """
    Collects every value a figure annotates, so the gate can re-derive it from the
    processed tables. Keys are free-form but must match what `tests/test_figures.py` checks.
    """

    def __init__(self, figure: str, script: str):
        self.figure = figure
        self.script = script
        self.values: dict = {}

    def set(self, key: str, value):
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
            "figure": self.figure,
            "script": self.script,
            "generated": datetime.now(timezone.utc).isoformat(),
            "seed": SEED,
            "inputs": sorted(str(q.relative_to(ROOT)) for q in PROC.glob("*.tsv")),
            "values": self.values,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return p


def save(fig, name: str, numbers: Numbers | None = None) -> None:
    """Write the figure as PNG (raster, 300 dpi) and PDF (vector), plus its numbers sidecar."""
    FIGURES.mkdir(exist_ok=True)
    png = FIGURES / f"{name}.png"
    pdf = FIGURES / f"{name}.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    print(f"  wrote {png.relative_to(ROOT)}  ({png.stat().st_size/1024:.0f} kB)")
    print(f"  wrote {pdf.relative_to(ROOT)}")
    if numbers is not None:
        print(f"  wrote {numbers.write().relative_to(ROOT)}")


def write_log(stem: str, lines: list[str]) -> Path:
    LOGS.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    p = LOGS / f"{stem}_{stamp}.log"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def pymol_version() -> str:
    try:
        out = subprocess.run(["pymol", "-cq", "-d", "print(cmd.get_version()[0])"],
                             capture_output=True, text=True, timeout=120)
        return out.stdout.strip().splitlines()[-1] if out.stdout.strip() else "unknown"
    except Exception as exc:  # noqa: BLE001 - reported, never silently swallowed
        return f"unavailable ({exc})"


# --------------------------------------------------------------------------------------
# Distribution panel — one implementation, used by every comparison figure
# --------------------------------------------------------------------------------------

def dist_panel(ax, groups: list[tuple[str, np.ndarray, str]], ylabel: str,
               rng: np.random.Generator, log_y: bool = False,
               annotate_fmt: str = "{:.2f}", ylim: tuple | None = None) -> dict:
    """
    Box + jittered strip for a small number of groups.

    Returns the medians it drew, so the caller can hand them to the Numbers recorder —
    the annotation and the recorded value come from the same computation, never two.

    Direct n= and median labels are mandatory here, not decorative: they are the relief
    that lets the low-contrast palette slots stay legible (see module docstring).
    """
    stats_out: dict = {}
    positions = np.arange(len(groups))
    for i, (label, vals, color) in enumerate(groups):
        vals = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
        stats_out[label] = {"n": int(len(vals)), "median": float(np.median(vals))}
        if len(vals) == 0:
            continue
        bp = ax.boxplot([vals], positions=[i], widths=0.38, showfliers=False,
                        patch_artist=True, zorder=2)
        for box in bp["boxes"]:
            box.set(facecolor=color, alpha=0.16, edgecolor=color, linewidth=0.7)
        for part in ("whiskers", "caps"):
            for art in bp[part]:
                art.set(color=color, linewidth=0.6)
        for med in bp["medians"]:
            med.set(color=color, linewidth=1.4)
        jitter = rng.uniform(-0.115, 0.115, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=1.6, color=color,
                   alpha=0.38, linewidths=0.0, zorder=3, rasterized=True)

    ax.set_xticks(positions)
    ax.set_xticklabels([g[0] for g in groups], rotation=32, ha="right",
                       rotation_mode="anchor")
    ax.set_ylabel(ylabel)
    if log_y:
        ax.set_yscale("symlog", linthresh=1.0)
    if ylim:
        ax.set_ylim(*ylim)
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)

    # A clipped axis hides points. Say how many, per group, rather than letting the reader
    # assume the strip shows everything — FoldX ddG has a long tail (max 50 kcal/mol) that
    # would otherwise flatten every box into a line.
    lo_lim, hi_lim = ax.get_ylim()
    top = hi_lim
    for i, (label, vals, color) in enumerate(groups):
        s = stats_out[label]
        arr = np.asarray([v for v in vals if v is not None and np.isfinite(v)], dtype=float)
        n_above = int((arr > hi_lim).sum())
        n_below = int((arr < lo_lim).sum())
        s["n_above_axis"], s["n_below_axis"] = n_above, n_below
        note = f"{s['n']}"
        if n_above:
            note += f"\n({n_above}↑)"
        ax.annotate(note, xy=(i, top), xytext=(0, -1.5), textcoords="offset points",
                    ha="center", va="top", fontsize=5.5, color=PALETTE["ink2"],
                    linespacing=1.1)
        if n_above:
            ax.plot([i], [hi_lim], marker="^", ms=2.6, color=color, clip_on=False, zorder=4)
    return stats_out
