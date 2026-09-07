#!/usr/bin/env python3
"""
60_typeset_pdf.py -- render the keyed manuscript to the two-column PDF.

`59_typeset_docx.py` produces the Word file; this produces the PDF proof of the same paper, from
the same source and with the same reference numbering, so the two cannot disagree.

It does not reimplement the typography. `56_typeset.py` already carries the page geometry
measured off Godwin 2023 and every workaround the base-14 Times family needs; this script
flattens `results/manuscript_paper.md` into the shape that script expects -- citation keys
resolved to superscript numbers, the reference list appended in first-appearance order -- and
then hands it over.

Writes: results/FBN1_manuscript.pdf, logs/60_typeset_pdf_<stamp>.log
"""
from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

SRC = ROOT / "results/manuscript_paper.md"
REFS = ROOT / "data/processed/references_final.tsv"
OUT = ROOT / "results/FBN1_manuscript.pdf"
FLAT = ROOT / "data/interim/manuscript_flattened.md"

SUPER = {str(i): c for i, c in enumerate("⁰¹²³⁴⁵⁶⁷⁸⁹")}


def _load56():
    spec = importlib.util.spec_from_file_location("_t56", ROOT / "scripts" / "56_typeset.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def sup(s: str) -> str:
    return "".join(SUPER.get(c, "⁻" if c == "-" else c) for c in s)


def main() -> int:
    md = SRC.read_text(encoding="utf-8")
    md = re.sub(r"<!--.*?-->", "", md, flags=re.S).lstrip()
    refs = pd.read_csv(REFS, sep="\t").fillna("")
    known = set(refs.key)

    order: list[str] = []

    def number(key: str) -> int:
        if key not in known:
            sys.exit(f"STOP: citation [@{key}] has no entry in {REFS.name}")
        if key not in order:
            order.append(key)
        return order.index(key) + 1

    def label(keys: str) -> str:
        """Resolve a key group to its number string, collapsing runs of three or more."""
        nums = sorted(number(k) for k in keys.split(","))
        out, i = [], 0
        while i < len(nums):
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append(f"{nums[i]}-{nums[j]}" if j - i >= 2
                       else ",".join(str(x) for x in nums[i:j + 1]))
            i = j + 1
        return ",".join(out)

    # ONE left-to-right pass over both citation forms. Numbering is assigned in the order the
    # substitution visits them, so running a separate earlier pass for the "(refs. ...)" form --
    # which is what this did at first -- hands numbers 1 and 2 to whichever references sit in
    # that one construction, and every number in the paper is then wrong. The alternation is
    # what keeps it to a single pass.
    both = re.compile(r"\(refs?\.\s*\[@([a-z0-9,]+)\]\)|\[@([a-z0-9,]+)\]")

    def one(m: re.Match) -> str:
        if m.group(1) is not None:
            return f"(refs. {label(m.group(1))})"      # prose, on the baseline
        return sup(label(m.group(2)))                  # superscript

    body = both.sub(one, md)

    if "[@" in body:
        sys.exit("STOP: unresolved citation markers remain: "
                 + ", ".join(sorted(set(re.findall(r"\[@[^\]]*\]", body)))[:5]))

    uncited = [k for k in refs.key if k not in order]
    if uncited:
        sys.exit("STOP: references never cited: " + ", ".join(uncited))

    by_key = {r.key: r for r in refs.itertuples()}
    lines = [body.rstrip(), "", "---", "", "## References", ""]
    for i, key in enumerate(order, start=1):
        r = by_key[key]
        doi = f" doi:{r.doi}" if r.doi else ""
        lines.append(f"{i}. {r.formatted}.{doi}")
    FLAT.parent.mkdir(parents=True, exist_ok=True)
    FLAT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    t56 = _load56()
    t56.SRC, t56.OUT = FLAT, OUT
    rc = t56.main()
    print(f"\n{len(order)} references, numbered by first appearance")
    print(f"flattened source kept at {FLAT.relative_to(ROOT)}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
