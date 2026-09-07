#!/usr/bin/env python3
"""
test_manuscript.py — gate on the typeset manuscript, its references and its figures.

`tests/test_v2.py` guards the analysis and the figure NUMBERS. This guards the document built
from them: that every reference is real and cited, that the two renderings agree on every
number they print, and that the figures are actually monochrome.

Each check is one way the manuscript could rot silently:

  M1  a citation key with no reference entry, or an entry nothing cites
  M2  a DOI that does not resolve, or the same paper entered twice
  M3  the numbers in the text disagreeing with the numbers in the list
  M4  the Word file losing its two-column geometry
  M5  a figure missing from the document
  M6  the DOCX and the PDF drifting apart
  M7  an unresolved citation marker printed as "[@key]"
  M8  colour creeping back into a figure set that has to survive greyscale

Run: .venv/bin/python tests/test_manuscript.py
"""
from __future__ import annotations

import re
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/manuscript_paper.md"
REFS = ROOT / "data/processed/references_final.tsv"
VERIFY = ROOT / "data/processed/references_verified.tsv"   # all 63; the
# 43-entry references.tsv is the earlier pass over the superseded manuscript
DOCX = ROOT / "results/FBN1_manuscript.docx"
PDF = ROOT / "results/FBN1_manuscript.pdf"
FIGS = ["fig1_module", "fig2_stability", "fig3_coordination", "fig4_interpretation",
        "figS1_validation", "figS2_structure"]

MIN_REFERENCES = 40          # the floor this manuscript was asked to clear

results: list[tuple[bool, str, str]] = []


def check(ok: bool, name: str, detail: str) -> None:
    results.append((bool(ok), name, detail))


def source_text() -> str:
    return re.sub(r"<!--.*?-->", "", SRC.read_text(encoding="utf-8"), flags=re.S)


def first_appearance(src: str) -> list[str]:
    order: list[str] = []
    for m in re.finditer(r"\[@([a-z0-9,]+)\]", src):
        for k in m.group(1).split(","):
            if k not in order:
                order.append(k)
    return order


def docx_paragraphs(path: Path) -> list[str]:
    from docx import Document
    return [p.text for p in Document(str(path)).paragraphs]


def main() -> int:
    src = source_text()
    refs = pd.read_csv(REFS, sep="\t").fillna("")
    order = first_appearance(src)

    # ---- M1 every key resolves, every reference is used ---------------------------------
    defined, used = set(refs.key), set(order)
    check(not (used - defined) and not (defined - used) and len(refs) >= MIN_REFERENCES,
          "M1 references cited",
          f"{len(refs)} references (floor {MIN_REFERENCES}), "
          f"{len(used - defined)} unresolved keys, {len(defined - used)} uncited entries")

    # ---- M2 DOIs are unique and were resolved against Crossref/Europe PMC ---------------
    dup = refs.doi[refs.doi.duplicated()].tolist()
    bad = [d for d in refs.doi if not re.match(r"^10\.\d{4,9}/\S+$", str(d))]
    resolved, unresolved_ok = "no verification file — run scripts/57_references.py", False
    if VERIFY.is_file():
        v = pd.read_csv(VERIFY, sep="\t")
        found = v.found_crossref | v.found_epmc
        resolved = f"{int(found.sum())}/{len(v)} resolved against Crossref or Europe PMC"
        # the verification must cover the list being shipped, not an older one
        unresolved_ok = bool(found.all()) and len(v) == len(refs)
    check(not dup and not bad and unresolved_ok, "M2 DOIs resolve",
          f"{len(refs)} DOIs, {len(dup)} duplicated, {len(bad)} malformed; {resolved}")

    # ---- M3 in-text numbers match the list ------------------------------------------------
    paras = docx_paragraphs(DOCX)
    listed = [m.group(2) for p in paras if (m := re.match(r"^(\d+)\.\s+(\S+)", p))]
    nums = [int(m.group(1)) for p in paras if (m := re.match(r"^(\d+)\.\s", p))]
    expect_first = [refs.set_index("key").formatted[k].split()[0] for k in order]
    check(nums == list(range(1, len(order) + 1)) and listed == expect_first,
          "M3 numbering by first appearance",
          f"{len(nums)} entries numbered 1..{max(nums) if nums else 0}, "
          f"order matches the text: {listed == expect_first}")

    # ---- M4 two-column geometry survived --------------------------------------------------
    from docx import Document
    from docx.oxml.ns import qn
    doc = Document(str(DOCX))
    cols = []
    for s in doc.sections:
        c = s._sectPr.find(qn("w:cols"))
        cols.append(int(c.get(qn("w:num")) or 1) if c is not None else 1)
    pw = doc.sections[0].page_width.pt
    ph = doc.sections[0].page_height.pt
    check(cols.count(2) >= 5 and abs(pw - 595.28) < 0.5 and abs(ph - 790.87) < 0.5,
          "M4 two-column geometry",
          f"{len(doc.sections)} sections, {cols.count(2)} two-column, "
          f"page {pw:.1f}x{ph:.1f} pt")

    # ---- M5 every figure is embedded -------------------------------------------------------
    with zipfile.ZipFile(DOCX) as z:
        media = [n for n in z.namelist() if n.startswith("word/media/")]
    check(len(media) == len(FIGS), "M5 figures embedded",
          f"{len(media)} images in the document, {len(FIGS)} expected")

    # ---- M6 the DOCX and the PDF agree ------------------------------------------------------
    pdf_ok, pdf_detail = True, "PDF not built"
    if PDF.is_file():
        import pymupdf
        with pymupdf.open(PDF) as d:
            ptext = "".join(pg.get_text() for pg in d)
            npages = d.page_count
        plist = re.findall(r"^\s*(\d+)\.\s+(\S+)", ptext, re.M)[:len(order)]
        pdf_ok = [x[1] for x in plist] == listed
        pdf_detail = (f"{npages} pages, {len(plist)} references, "
                      f"same order as the DOCX: {pdf_ok}")
    check(pdf_ok, "M6 DOCX and PDF agree", pdf_detail)

    # ---- M7 no unresolved markers reached either artifact -----------------------------------
    leaked_docx = [p for p in paras if "[@" in p]
    leaked_pdf = "[@" in ptext if PDF.is_file() else False
    check(not leaked_docx and not leaked_pdf, "M7 no unresolved citations",
          f"{len(leaked_docx)} in the DOCX, {'yes' if leaked_pdf else 'none'} in the PDF")

    # ---- M8 the figures are monochrome -------------------------------------------------------
    worst, worst_name = 0, ""
    for stem in FIGS:
        png = ROOT / "figures" / f"{stem}.png"
        if not png.is_file():
            worst, worst_name = 999, f"{stem} missing"
            break
        a = np.asarray(Image.open(png).convert("RGB"), dtype=np.int16)
        sat = int((a.max(2) - a.min(2)).max())
        if sat > worst:
            worst, worst_name = sat, stem
    check(worst == 0, "M8 figures monochrome",
          f"largest RGB spread {worst} across {len(FIGS)} figures"
          + (f" (worst: {worst_name})" if worst else " — pixel-exact greyscale"))

    print("manuscript gate — the document, its references and its figures\n")
    for ok, name, detail in results:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<34} {detail}")
    n_ok = sum(1 for ok, _, _ in results if ok)
    print(f"\n{n_ok}/{len(results)} checks passed")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
