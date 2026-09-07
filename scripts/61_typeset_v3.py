#!/usr/bin/env python3
"""
61_typeset_v3.py -- set the citation-free manuscript as a two-column Word file and PDF.

`results/manuscript_v4.md` is the rewritten draft. It carries no [@key] citations, because the
mentor is inserting references from Mendeley, and its Results headings were reworded, so the two
existing typesetters need two small pieces of help.

  THE DOCX comes from `59_typeset_docx.py`, which now takes `--source`, `--out` and
  `--no-references`. Nothing else about it changes, so the two documents are the same page.

  THE PDF comes from `56_typeset.py`, which floats each figure to the top of a page and finds
  the page by matching a Results heading. Those headings were reworded, so PLACEMENT is
  overridden here rather than edited in 56, which the keyed pipeline still uses unchanged.

Writes: results/FBN1_manuscript_v4.docx, results/FBN1_manuscript_v4.pdf
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/manuscript_v4.md"
DOCX = ROOT / "results/FBN1_manuscript_v4.docx"
PDF = ROOT / "results/FBN1_manuscript_v4.pdf"
STAMP = ROOT / "results/.manuscript_v4.build.json"

# Results headings of manuscript_v4.md, in the order the figures are cited. Two moves on
# 2026-09-06: the geometry validation went ahead of the claims that depend on it and enrichment
# became its own section, and then the whole set was rewritten as neutral noun phrases for a
# traditional journal, so the heading a figure hangs from changed twice in one day.
PLACEMENT_V3 = [
    ("Distribution of pathogenic variants across the cbEGF array", 1, "fig1_module"),
    ("Folding stability at cysteine and calcium-ligand positions", 2, "fig2_stability"),
    ("Calcium coordination in the mutant models", 3, "fig3_coordination"),
    ("Comparison with stability prediction and AlphaMissense", 4, "fig4_interpretation"),
]


def front_matter(md: str) -> tuple[str, str]:
    """
    The author block, taken from the source instead of retyped here.

    56 drops every paragraph above the abstract and sets its own subtitle and standfirst, which
    was fine while the source carried nothing there but a note. The draft now carries the author
    line, the affiliation and the corresponding author, and a submission cannot lose them, so
    they are read out of the markdown and handed to 56 as those two fields.
    """
    body = md.split("## Abstract", 1)[0]
    paras, buf = [], []
    for line in body.splitlines()[1:]:          # skip the title line
        if line.strip() in ("", "---"):
            if buf:
                paras.append(" ".join(buf))
                buf = []
            continue
        buf.append(line.strip())
    if buf:
        paras.append(" ".join(buf))
    if len(paras) < 2:
        sys.exit("STOP: no author block found above the abstract")
    return paras[0], "<br/>".join(paras[1:])


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else ""


def guard_hand_edits(force: bool) -> None:
    """
    Refuse to overwrite a Word file somebody has been editing by hand.

    This script regenerates the DOCX from the markdown, so a hand edit made in Word is destroyed
    silently the next time anything rebuilds. A modification time cannot detect that, because
    this script writes the DOCX a second after the markdown and the DOCX is therefore always the
    newer file. So the hash of what was written is recorded instead, and a mismatch means someone
    else wrote it.

    `--force` overwrites anyway, and is the right answer only once the edits are back in
    `results/manuscript_v4.md`, which is the source of truth.
    """
    if not STAMP.is_file():
        return
    rec = json.loads(STAMP.read_text(encoding="utf-8"))
    edited = [p.name for p, k in ((DOCX, "docx"), (PDF, "pdf"))
              if rec.get(k) and p.is_file() and _sha(p) != rec[k]]
    if not edited:
        return
    if force:
        print(f"WARNING: overwriting hand-edited {', '.join(edited)} because --force was given")
        return
    sys.exit(
        f"STOP: {' and '.join(edited)} changed since this script last wrote it, so it is being\n"
        f"      edited by hand. Rebuilding would destroy those edits.\n\n"
        f"      The markdown is the source of truth. Either put the edits into\n"
        f"      results/manuscript_v4.md and re-run, or re-run with --force to discard them.\n"
        f"      To keep a copy first:  cp results/{edited[0]} results/{edited[0]}.handedit")


def _load(stem: str):
    spec = importlib.util.spec_from_file_location(f"_{stem}", ROOT / "scripts" / f"{stem}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main() -> int:
    if not SRC.is_file():
        sys.exit(f"STOP: {SRC.relative_to(ROOT)} does not exist")
    guard_hand_edits("--force" in sys.argv)
    md = SRC.read_text(encoding="utf-8")
    if "[@" in md:
        sys.exit("STOP: the v3 source is meant to be citation-free but carries [@key] markers")
    for h, _, _ in PLACEMENT_V3:
        if f"### {h}" not in md:
            sys.exit(f"STOP: PLACEMENT_V3 heading not found in the manuscript -- {h!r}")

    print("--- Word ---")
    t59 = _load("59_typeset_docx")
    t59.SRC, t59.OUT, t59.NO_REFS = SRC, DOCX, True
    # 59 keys figure placement on Results headings too, and the v3 draft renamed all of them.
    t59.PLACEMENT = {h: stem for h, _, stem in PLACEMENT_V3}
    rc = t59.main()
    if rc:
        return rc

    print("\n--- PDF ---")
    t56 = _load("56_typeset")
    t56.SRC, t56.OUT = SRC, PDF
    t56.PLACEMENT = PLACEMENT_V3
    t56.HEADER_RIGHT = "Preprint — 6 September 2026"
    authors, standfirst = front_matter(md)
    t56.SUBTITLE = t56.inline(authors)
    t56.STANDFIRST = "<br/>".join(t56.inline(x) for x in standfirst.split("<br/>"))
    rc = t56.main()
    if rc == 0:
        STAMP.write_text(json.dumps({"docx": _sha(DOCX), "pdf": _sha(PDF),
                                     "source": _sha(SRC)}, indent=2) + "\n", encoding="utf-8")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
