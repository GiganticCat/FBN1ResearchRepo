#!/usr/bin/env python3
"""
18_fig6_enrichment.py — PHASE 6, Figure 6. Enrichment, reported as supporting context.

This figure is deliberately NOT the headline. ClinVar's pathogenic classifications were made
using ClinGen's PM1 critical-residue rule — the very rule tested here — so the ClinVar-labelled
contrasts partly measure the classification rule rather than biology (LIMITATIONS.md #1). The
figure is drawn so a reader cannot take the large odds ratios at face value:

  A  against the all-possible-missense background — position share is the correct null, since
     every position offers the same 19 substitutions
  B  odds ratios, with the circular ClinVar contrasts and the orthogonal AlphaMissense
     contrasts on the same axis and visibly distinguished

Effect sizes and intervals are READ from `results/phase5_statistics.tsv`; this script runs no
test of its own. Both panels use a log x-axis because these are ratio measures whose intervals
are multiplicative — a linear axis would make an OR of 692 with a 41–11450 interval unplottable
next to an OR of 3.9.

Writes: figures/fig6_enrichment.{png,pdf}, .numbers.json, logs/18_fig6_*.log
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

FIG = "fig6_enrichment"

BACKGROUND = [
    ("pathogenic at cbEGF cysteine (removing) vs background",
     "cbEGF cysteine-removing"),
    ("pathogenic at calcium-consensus residue vs background",
     "Calcium-consensus residue"),
]

ODDS = [
    ("pathogenic vs benign at any VCEP critical residue",
     "Any VCEP critical residue", "clinvar"),
    ("pathogenic vs benign at cbEGF cysteine (removing)",
     "cbEGF cysteine-removing", "clinvar"),
    ("pathogenic vs benign at calcium-consensus residue",
     "Calcium-consensus residue", "clinvar"),
    ("AM-pathogenic vs AM-benign at cbEGF cysteine (removing)",
     "cbEGF cysteine-removing", "orthogonal"),
    ("AM-pathogenic vs AM-benign at calcium-consensus residue",
     "Calcium-consensus residue", "orthogonal"),
]

SOURCE_COLOR = {"clinvar": fs.PALETTE["pathogenic"], "orthogonal": fs.PALETTE["benign"]}
SOURCE_LABEL = {
    "clinvar": "ClinVar pathogenic vs benign — CIRCULAR: these labels were assigned using PM1",
    "orthogonal": "AlphaMissense pathogenic vs benign — never used FBN1-specific rules",
}


def rows(tbl: pd.DataFrame, spec) -> list[dict]:
    out = []
    for entry in spec:
        test, label = entry[0], entry[1]
        hit = tbl[tbl.test == test]
        if len(hit) != 1:
            raise SystemExit(f"expected one row for {test!r}, found {len(hit)}")
        r = hit.iloc[0]
        out.append({"label": label, "effect": float(r.effect), "lo": float(r.ci_low),
                    "hi": float(r.ci_high), "p": float(r.p_raw), "q": float(r.q_bh),
                    "n": str(r.n), "source": entry[2] if len(entry) > 2 else None})
    return out


def forest(ax, data: list[dict], colours, xlabel: str) -> None:
    y = np.arange(len(data))[::-1]
    for yy, d in zip(y, data):
        c = colours(d)
        ax.plot([d["lo"], d["hi"]], [yy, yy], color=c, lw=1.6, solid_capstyle="round",
                zorder=3)
        for edge in ("lo", "hi"):
            ax.plot([d[edge]], [yy], marker="|", ms=6, color=c, zorder=3)
        ax.plot([d["effect"]], [yy], marker="o", ms=6.5, color=c,
                mec=fs.PALETTE["surface"], mew=0.9, zorder=4)
    ax.axvline(1.0, color=fs.PALETTE["muted"], lw=0.9, ls="--", zorder=1)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([d["label"] for d in data], fontsize=7.5)
    ax.set_ylim(-0.7, len(data) - 0.3)
    ax.set_xlabel(xlabel)
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


def main() -> int:
    fs.apply_style()
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    log = [f"Phase 6 — Figure 6 — seed {fs.SEED}"]
    num = fs.Numbers(FIG, Path(__file__).name)
    tbl = pd.read_csv(fs.RESULTS / "phase5_statistics.tsv", sep="\t")

    bg = rows(tbl, BACKGROUND)
    od = rows(tbl, ODDS)

    fig = plt.figure(figsize=(7.4, 5.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.9], hspace=0.95)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[1, 0])

    # ---- A: against the all-possible-missense background ---------------------------
    forest(axA, bg, lambda d: fs.PALETTE["cbegf_cys"] if "cysteine" in d["label"]
           else fs.PALETTE["ca_site"],
           "Fold enrichment of pathogenic variants  (log scale)")
    axA.set_title("Against the background of every possible missense change",
                  loc="left", fontsize=8.5, pad=8)
    for d, yy in zip(bg, np.arange(len(bg))[::-1]):
        axA.annotate(f"{d['effect']:.1f}× [{d['lo']:.1f}, {d['hi']:.1f}]   p = {d['p']:.0e}",
                     xy=(d["hi"], yy), xytext=(7, 0), textcoords="offset points",
                     va="center", ha="left", fontsize=6.8, color=fs.PALETTE["ink2"])
        key = "cys" if "cysteine" in d["label"] else "ca"
        num.set(f"background_fold_{key}", d["effect"])
        num.set(f"background_ci_{key}", [d["lo"], d["hi"]])
        log.append(f"  background {key}: {d['effect']:.3f} [{d['lo']:.3f}, {d['hi']:.3f}]")
    axA.set_xlim(0.8, 30)

    # ---- B: odds ratios -------------------------------------------------------------
    forest(axB, od, lambda d: SOURCE_COLOR[d["source"]],
           "Odds ratio, variants at the site class vs elsewhere  (log scale)")
    axB.set_title("Pathogenic vs benign — and the same contrast without ClinVar's labels",
                  loc="left", fontsize=8.5, pad=8)
    for d, yy in zip(od, np.arange(len(od))[::-1]):
        txt = f"OR {d['effect']:.1f} [{d['lo']:.0f}, {d['hi']:.0f}]   n = {d['n']}"
        axB.annotate(txt, xy=(d["hi"], yy), xytext=(7, 0), textcoords="offset points",
                     va="center", ha="left", fontsize=6.8, color=fs.PALETTE["ink2"])
        num.set(f"or_{d['source']}_{d['label'][:14].strip().replace(' ', '_')}", d["effect"])
        log.append(f"  OR {d['source']:10} {d['label']:28} {d['effect']:.2f} "
                   f"[{d['lo']:.2f}, {d['hi']:.2f}]  n={d['n']}")
    axB.set_xlim(0.4, 4e6)
    handles = [Line2D([0], [0], marker="o", ls="-", lw=1.6, ms=6, color=SOURCE_COLOR[k],
                      mec=fs.PALETTE["surface"], label=SOURCE_LABEL[k])
               for k in ("clinvar", "orthogonal")]
    axB.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.42, -0.62), ncol=1,
               fontsize=6.8, handletextpad=0.5, labelspacing=0.45)

    for ax, letter in ((axA, "A"), (axB, "B")):
        fs.panel_label(ax, letter, dx=-0.27, dy=1.06)

    fig.suptitle("Enrichment at critical residues — supporting evidence, not the headline",
                 fontsize=10, x=0.02, ha="left", y=0.985)
    fig.text(0.02, 0.055,
             "The red intervals are inflated by construction and are not an independent "
             "discovery: ClinVar's pathogenic calls used ClinGen PM1, the same critical-residue "
             "rule\nbeing tested here. The blue intervals reproduce the direction from "
             "AlphaMissense scores alone, which owe nothing to FBN1-specific rules. Panel A's "
             "asymmetry —\ncysteines ~6×, calcium residues only ~1.4× — is consistent with "
             "calcium-site variants being under-ascertained because they look thermodynamically "
             "benign (Fig. 2).",
             fontsize=6.4, color=fs.PALETTE["ink2"], va="top", ha="left", linespacing=1.7)
    fig.subplots_adjust(left=0.245, right=0.80, top=0.90, bottom=0.235)

    fs.save(fig, FIG, num)
    plt.close(fig)
    print(f"  log -> {fs.write_log('18_fig6', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
