#!/usr/bin/env python3
"""
45_manuscript_final_docx.py — render the final manuscript to a styled Word document.

`results/manuscript_final.md` is the source of truth: edit the markdown, re-run this, and the
document follows. The rendering machinery is `39_manuscript_v2_docx.py`, imported rather than
copied, so the two documents cannot drift in styling; only what changed between drafts is
written here — the source file, the output name, and which figure follows which section.

Each figure is placed after the Results section that first discusses it. The two supplementary
figures are not placed inline; their legends stay in the legends section at the end.

Writes: results/FBN1_manuscript_final.docx
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

spec = importlib.util.spec_from_file_location(
    "_docx_v2", ROOT / "scripts" / "39_manuscript_v2_docx.py")
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)

R.SRC = ROOT / "results" / "manuscript_final.md"
R.OUT = ROOT / "results" / "FBN1_manuscript_final.docx"
R.PLACEMENT = [
    ("Variants concentrate at the two features of the module", 1, "fig1_module"),
    ("Cysteine loss destabilises the domain; calcium-ligand loss does not", 2, "fig2_stability"),
    ("Donor chemistry decides the outcome", 3, "fig3_coordination"),
    ("Neither instrument in routine use sees the difference", 4, "fig4_interpretation"),
]

if __name__ == "__main__":
    sys.exit(R.main())
