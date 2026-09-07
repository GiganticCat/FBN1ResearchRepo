#!/usr/bin/env python3
"""
59_typeset_docx.py -- set the manuscript as a two-column Word document.

Source is `results/manuscript_paper.md` (citations as [@key]) plus
`data/processed/references_final.tsv`. Reference numbers are assigned HERE, by order of first
appearance, so inserting a citation anywhere renumbers the paper correctly and no number in the
text can disagree with the list.

Geometry is the same measurement `56_typeset.py` took off
`resources/papers/Godwin_2023_NSMB_microfibril_cryoEM.pdf` with PyMuPDF, so the Word file and the
PDF are the same page:

    page          595.28 x 790.87 pt
    text block    x 39.7 to 563.4, y 46.6 to 745.7
    columns       two, 256.9 pt wide, 9.9 pt gutter
    body          serif, 8.2 pt, justified
    captions      serif, 7.0 pt, opening with a bold "Fig. N | ..." run

Two things Word does that reportlab does not, and that this file has to work around:

  COLUMNS ARE A PROPERTY OF A SECTION, not of a frame. A full-width figure inside two-column
  text is therefore three sections -- two columns, then a continuous break to one column for the
  figure, then a continuous break back. python-docx exposes no column API at all, so `_columns`
  writes the w:cols element directly.

  SUPERSCRIPT IS A RUN PROPERTY, so a citation cannot be part of the surrounding text run. Every
  paragraph is tokenised into runs and the citation numbers get their own.

Writes: results/FBN1_manuscript.docx, logs/59_typeset_docx_<stamp>.log
"""
from __future__ import annotations

import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt, RGBColor
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/manuscript_paper.md"
REFS = ROOT / "data/processed/references_final.tsv"
OUT = ROOT / "results/FBN1_manuscript.docx"

# `--source X --out Y --no-references` sets the same document from a citation-free markdown.
# The mentor is inserting citations from Mendeley, so `results/manuscript_v3.md` carries no
# [@key] markers at all and must not gain a numbered reference list. Defaults are unchanged, so
# a bare run still builds the keyed paper exactly as before.
NO_REFS = False
FIGDIR = ROOT / "figures"
LOGS = ROOT / "logs"

# --- geometry, measured from the reference page ------------------------------------------
PW, PH = 595.28, 790.87
ML, MR, MT, MB = 39.7, 31.9, 46.6, 45.2
TEXT_W = PW - ML - MR                    # 523.7
GUTTER = 9.9
COL_W = (TEXT_W - GUTTER) / 2            # 256.9

SERIF = "Times New Roman"
BODY_PT = 8.2
CAPTION_PT = 7.0
INK = RGBColor(0x10, 0x10, 0x10)

# Which figure follows which Results subsection. Same order as `56_typeset.py`.
PLACEMENT = {
    "Variants concentrate at the two features of the module": "fig1_module",
    "Cysteine loss destabilises the domain; calcium-ligand loss does not": "fig2_stability",
    "The donor rule holds at the Asp and Glu positions and fails at the β-hydroxylation Asn":
        "fig3_coordination",
    "Neither instrument in routine use sees the difference": "fig4_interpretation",
}
SUPP_FIGS = ["figS1_validation", "figS2_structure"]

TITLE = ("Calcium-ligand and cysteine variants in fibrillin-1 damage the cbEGF module "
         "in different ways")

# Sections that are set full width rather than in columns.
FRONT_MATTER = {"Abstract"}
# Dropped: captions are placed under their figures instead, as on a journal page.
DROP_SECTIONS = {"Figure legends"}


# ------------------------------------------------------------------------------------------
# Word plumbing python-docx does not expose
# ------------------------------------------------------------------------------------------

def _columns(section, num: int, space_pt: float = GUTTER) -> None:
    """Set the column count on a section. python-docx has no API for this."""
    sectPr = section._sectPr
    cols = sectPr.find(qn("w:cols"))
    if cols is None:
        cols = OxmlElement("w:cols")
        sectPr.append(cols)
    cols.set(qn("w:num"), str(num))
    cols.set(qn("w:space"), str(int(round(space_pt * 20))))   # twips
    cols.set(qn("w:equalWidth"), "1")


def _page(section) -> None:
    section.page_width, section.page_height = Pt(PW), Pt(PH)
    section.left_margin, section.right_margin = Pt(ML), Pt(MR)
    section.top_margin, section.bottom_margin = Pt(MT), Pt(MB)


def _new_section(doc, ncols: int):
    s = doc.add_section(WD_SECTION.CONTINUOUS)
    _page(s)
    _columns(s, ncols)
    return s


def _style(doc) -> None:
    n = doc.styles["Normal"]
    n.font.name = SERIF
    n.font.size = Pt(BODY_PT)
    n.font.color.rgb = INK
    # East-Asian font mapping, or Word substitutes something else for the run
    n.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)
    pf = n.paragraph_format
    pf.space_before, pf.space_after = Pt(0), Pt(0)
    pf.line_spacing = 1.0
    pf.widow_control = True


# ------------------------------------------------------------------------------------------
# Markdown -> runs
# ------------------------------------------------------------------------------------------

TOKEN = re.compile(
    r"(\*\*.+?\*\*)"          # bold
    r"|(?<!\*)(\*[^*]+?\*)(?!\*)"   # italic
    r"|(`[^`]+?`)"            # code
    r"|(\[@[a-z0-9,]+\])"     # citation, superscript
    r"|(\[#[a-z0-9,]+\])"     # citation read as prose: "(refs. 8,9)"
)

# "of which there are 43 (refs. [@a,b])" is read aloud, so those numbers sit on the baseline
# rather than in superscript. This is a pure text rewrite that assigns no numbers, so it cannot
# disturb the first-appearance ordering the way an early numbering pass would.
PROSE_CITE = re.compile(r"\(refs?\.\s*\[@([a-z0-9,]+)\]\)")


class Numbering:
    """Reference numbers, assigned in order of first appearance."""

    def __init__(self, keys: set[str]):
        self.known, self.order, self.n = keys, [], {}

    def cite(self, key: str) -> int:
        if key not in self.known:
            raise KeyError(key)
        if key not in self.n:
            self.order.append(key)
            self.n[key] = len(self.order)
        return self.n[key]

    def render(self, group: str) -> str:
        """"[@a,b,c]" -> "1,2,3", collapsing runs of three or more into a range."""
        nums = sorted(self.cite(k) for k in group.split(","))
        out, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2
                       else ",".join(str(x) for x in nums[i:j + 1]))
            i = j + 1
        return ",".join(out)


def add_runs(p, text: str, num: Numbering, size: float = BODY_PT,
             bold: bool = False, italic: bool = False) -> None:
    """Tokenise one paragraph of markdown into Word runs."""
    text = PROSE_CITE.sub(r"(refs. [#\1])", text)
    pos = 0
    for m in TOKEN.finditer(text):
        if m.start() > pos:
            _run(p, text[pos:m.start()], size, bold, italic)
        b, i, c, cite, prose = m.groups()
        if b:
            _run(p, b[2:-2], size, True, italic)
        elif i:
            _run(p, i[1:-1], size, bold, True)
        elif c:
            r = _run(p, c[1:-1], size - 0.6, bold, italic)
            r.font.name = "Consolas"
            r.element.rPr.rFonts.set(qn("w:eastAsia"), "Consolas")
        elif cite:
            r = _run(p, num.render(cite[2:-1]), size, False, False)
            r.font.superscript = True
        elif prose:
            _run(p, num.render(prose[2:-1]), size, False, False)
        pos = m.end()
    if pos < len(text):
        _run(p, text[pos:], size, bold, italic)


def _run(p, text: str, size: float, bold: bool, italic: bool):
    r = p.add_run(text)
    r.font.name, r.font.size = SERIF, Pt(size)
    r.font.bold, r.font.italic = bold, italic
    r.font.color.rgb = INK
    r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)
    return r


def para(doc, text: str, num: Numbering, *, size=BODY_PT, justify=True, first_indent=0.0,
         space_before=0.0, space_after=0.0, bold=False, italic=False, align=None):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.alignment = align if align is not None else (
        WD_ALIGN_PARAGRAPH.JUSTIFY if justify else WD_ALIGN_PARAGRAPH.LEFT)
    pf.first_line_indent = Pt(first_indent)
    pf.space_before, pf.space_after = Pt(space_before), Pt(space_after)
    pf.line_spacing = 1.02
    add_runs(p, text, num, size=size, bold=bold, italic=italic)
    return p


# ------------------------------------------------------------------------------------------
# Source parsing
# ------------------------------------------------------------------------------------------

def parse(md: str) -> tuple[list[dict], dict[str, str]]:
    """Flatten the markdown into a block list, and pull the figure captions out separately."""
    md = re.sub(r"<!--.*?-->", "", md, flags=re.S)
    blocks, captions = [], {}
    section = None
    buf: list[str] = []

    def flush():
        if buf:
            blocks.append({"kind": "p", "section": section, "text": " ".join(buf).strip()})
            buf.clear()

    for raw in md.splitlines():
        line = raw.rstrip()
        if line.strip() == "---":
            flush()
            continue
        if line.startswith("# "):
            flush()
            section = None
            continue                       # the title is set from TITLE
        if line.startswith("## "):
            flush()
            section = line[3:].strip()
            blocks.append({"kind": "h2", "section": section, "text": section})
            continue
        if line.startswith("### "):
            flush()
            blocks.append({"kind": "h3", "section": section, "text": line[4:].strip()})
            continue
        if not line.strip():
            flush()
            continue
        buf.append(line.strip())
    flush()

    # captions: "**Fig. 1 | title.** body..." inside the Figure legends section
    cap_blocks = [b for b in blocks if b["section"] == "Figure legends" and b["kind"] == "p"]
    for b in cap_blocks:
        m = re.match(r"\*\*(Supplementary Fig\.|Fig\.)\s*(S?\d+)\s*\|", b["text"])
        if m:
            captions[m.group(2)] = b["text"]
    return blocks, captions


def figure_key(name: str) -> str:
    """fig3_coordination -> '3';  figS1_validation -> 'S1'."""
    m = re.match(r"fig(S?\d+)_", name)
    return m.group(1) if m else name


# ------------------------------------------------------------------------------------------

def add_figure(doc, num: Numbering, name: str, captions: dict, log: list) -> None:
    """A full-width figure and its caption, between continuous section breaks."""
    png = FIGDIR / f"{name}.png"
    if not png.is_file():
        sys.exit(f"STOP: missing {png}")
    _new_section(doc, 1)

    with Image.open(png) as im:
        w, h = im.size
    width_pt = TEXT_W
    height_pt = width_pt * h / w
    p = doc.add_paragraph()
    p.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(2.5)
    p.add_run().add_picture(str(png), width=Pt(width_pt))

    key = figure_key(name)
    cap = captions.get(key)
    if cap is None:
        sys.exit(f"STOP: no caption for figure {key} ({name})")
    cp = para(doc, cap, num, size=CAPTION_PT, justify=True, space_after=7)
    cp.paragraph_format.line_spacing = 1.0
    log.append(f"  figure {key:>2}  {name}  {w}x{h}px -> {width_pt:.0f}x{height_pt:.0f}pt")
    _new_section(doc, 2)


def main() -> int:
    md = SRC.read_text(encoding="utf-8")
    refs = pd.read_csv(REFS, sep="\t").fillna("")
    blocks, captions = parse(md)
    num = Numbering(set(refs.key))
    title = TITLE
    for line in md.splitlines():
        if line.startswith("# "):
            title = line[2:].strip()
            break
    log = [f"59_typeset_docx {datetime.now(timezone.utc).isoformat()}"]

    placed: set[str] = set()

    doc = Document()
    _style(doc)
    _page(doc.sections[0])
    _columns(doc.sections[0], 1)

    # ---------------------------------------------------------------- title block, full width
    t = doc.add_paragraph()
    t.paragraph_format.space_after = Pt(5)
    r = t.add_run(title)
    r.font.name, r.font.size, r.font.bold = SERIF, Pt(15), True
    r.font.color.rgb = INK
    r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)

    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(3)
    rr = p.add_run("Structural bioinformatics of pathogenic FBN1 missense variation "
                   "in Marfan syndrome")
    rr.font.name, rr.font.size, rr.font.italic = SERIF, Pt(8.6), True
    rr.font.color.rgb = INK
    rr.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)

    # the draft line and provenance note, carried from the source rather than retyped
    strap = [b for b in blocks if b["section"] is None and b["kind"] == "p"]
    for b in strap:
        sp = para(doc, b["text"], num, size=6.8, justify=False, space_after=9)
        for r in sp.runs:
            r.font.color.rgb = RGBColor(0x55, 0x55, 0x55)

    p = doc.add_paragraph()
    p.paragraph_format.space_after = Pt(2)
    ra = p.add_run("ABSTRACT")
    ra.font.name, ra.font.size, ra.font.bold = SERIF, Pt(8.2), True
    ra.font.color.rgb = INK
    ra.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)

    # ---------------------------------------------------------------- abstract, full width
    abstract = [b for b in blocks if b["section"] == "Abstract" and b["kind"] == "p"]
    for i, b in enumerate(abstract):
        para(doc, b["text"], num, size=8.6, justify=True,
             space_after=4 if i < len(abstract) - 1 else 9)

    # ---------------------------------------------------------------- body, two columns
    _new_section(doc, 2)

    pending_fig: str | None = None
    body = [b for b in blocks
            if b["section"] not in FRONT_MATTER | DROP_SECTIONS | {None, "References"}]

    for b in body:
        if b["kind"] == "h2":
            if pending_fig:
                add_figure(doc, num, pending_fig, captions, log)
                placed.add(pending_fig)
                pending_fig = None
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.space_after = Pt(2.5)
            r = p.add_run(b["text"])
            r.font.name, r.font.size, r.font.bold = SERIF, Pt(10.5), True
            r.font.color.rgb = INK
            r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)
        elif b["kind"] == "h3":
            if pending_fig:
                add_figure(doc, num, pending_fig, captions, log)
                placed.add(pending_fig)
                pending_fig = None
            p = doc.add_paragraph()
            p.paragraph_format.space_before = Pt(6)
            p.paragraph_format.space_after = Pt(1.5)
            r = p.add_run(b["text"])
            r.font.name, r.font.size, r.font.bold = SERIF, Pt(8.6), True
            r.font.color.rgb = INK
            r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)
            pending_fig = PLACEMENT.get(b["text"])
        else:
            para(doc, b["text"], num, first_indent=0, space_after=3.2)

    if pending_fig:
        add_figure(doc, num, pending_fig, captions, log)
        placed.add(pending_fig)

    # PLACEMENT is keyed on Results headings, so a reworded heading silently drops its figure.
    # That happened once, when the v3 draft renamed every Results section.
    never = [f for f in PLACEMENT.values() if f not in placed]
    if never:
        sys.exit("STOP: figures never placed: " + ", ".join(never)
                 + " — a heading in PLACEMENT does not match the manuscript")

    # ---------------------------------------------------------------- supplementary figures
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.space_after = Pt(2.5)
    r = p.add_run("Supplementary figures")
    r.font.name, r.font.size, r.font.bold = SERIF, Pt(10.5), True
    r.font.color.rgb = INK
    r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)
    for name in SUPP_FIGS:
        add_figure(doc, num, name, captions, log)

    # ---------------------------------------------------------------- references
    if NO_REFS:
        if "[@" in md:
            sys.exit("STOP: --no-references given but the source still carries [@key] markers")
        LOGS.mkdir(exist_ok=True)
        doc.save(OUT)
        log.append("  no reference list (--no-references)")
        log.append(f"  {len(body)} body blocks, {len(captions)} captions")
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        (LOGS / f"59_typeset_docx_{stamp}.log").write_text("\n".join(log) + "\n",
                                                           encoding="utf-8")
        print("\n".join(log[1:]))
        shown = OUT.relative_to(ROOT) if OUT.is_relative_to(ROOT) else OUT
        print(f"wrote {shown} ({OUT.stat().st_size/1024:.0f} kB)")
        return 0

    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(8)
    p.paragraph_format.space_after = Pt(2.5)
    r = p.add_run("References")
    r.font.name, r.font.size, r.font.bold = SERIF, Pt(10.5), True
    r.font.color.rgb = INK
    r.element.rPr.rFonts.set(qn("w:eastAsia"), SERIF)

    by_key = {row.key: row for row in refs.itertuples()}
    uncited = [k for k in refs.key if k not in num.n]
    if uncited:
        sys.exit("STOP: references defined but never cited: " + ", ".join(uncited))

    for i, key in enumerate(num.order, start=1):
        row = by_key[key]
        rp = doc.add_paragraph()
        pf = rp.paragraph_format
        pf.space_after, pf.line_spacing = Pt(1.6), 1.0
        pf.left_indent, pf.first_line_indent = Pt(11), Pt(-11)
        pf.alignment = WD_ALIGN_PARAGRAPH.LEFT
        _run(rp, f"{i}. ", CAPTION_PT, False, False)
        add_runs(rp, row.formatted + ".", num, size=CAPTION_PT)
        if row.doi:
            _run(rp, f" doi:{row.doi}", CAPTION_PT, False, False)

    LOGS.mkdir(exist_ok=True)
    doc.save(OUT)

    log.append(f"  {len(num.order)} references, numbered by first appearance")
    log.append(f"  {len(body)} body blocks, {len(captions)} captions")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (LOGS / f"59_typeset_docx_{stamp}.log").write_text("\n".join(log) + "\n", encoding="utf-8")

    print("\n".join(log[1:]))
    shown = OUT.relative_to(ROOT) if OUT.is_relative_to(ROOT) else OUT
    print(f"wrote {shown} ({OUT.stat().st_size/1024:.0f} kB)")
    return 0


def _cli() -> int:
    global SRC, OUT, NO_REFS
    args = sys.argv[1:]
    while args:
        a = args.pop(0)
        if a == "--source":
            SRC = (ROOT / args.pop(0)).resolve()
        elif a == "--out":
            OUT = (ROOT / args.pop(0)).resolve()
        elif a == "--no-references":
            NO_REFS = True
        else:
            sys.exit(f"STOP: unknown argument {a}")
    return main()


if __name__ == "__main__":
    raise SystemExit(_cli())
