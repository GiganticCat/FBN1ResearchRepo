#!/usr/bin/env python3
"""
62_prose_diagnostic.py — measure the two things the mentor says are wrong with the prose.

The mentor rejected two drafts for the same two reasons, sentences that carry too much and
vocabulary a reader cannot decode. This applies the diagnostic from the academic-paper-composer
skill (`Issue 3: Unclear Prose` and `Mistake 4: Jargon Overload` in
`.claude/skills/academic-paper-composer/references/`) as a measurement rather than an opinion,
so a revision can be checked instead of argued about.

  LOAD    how much one sentence carries. Counted as clauses, from finite-verb and connective
          markers, plus parenthetical asides and numbers. A sentence over the threshold is
          asking the reader to hold too many facts at once.

  JARGON  every domain term, and whether it was defined in plain words before first use.
          The skill's test is whether a graduate student in an adjacent field could follow.

Writes: results/prose_diagnostic.md
Run:    .venv/bin/python scripts/62_prose_diagnostic.py [source.md]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT = ROOT / "results/manuscript_v4.md"
OUT = ROOT / "results/prose_diagnostic.md"

# Domain terms a non-specialist cannot decode without help. `plain` is the wording that would
# count as having defined it.
TERMS = {
    "bond-valence sum":     "a geometric score of how well an ion is held",
    "valence":              "the same score, used as shorthand",
    "cognate site":         "the one calcium site the mutated residue belongs to",
    "bystander":            "the other sites in the same construct",
    "coordination":         "the set of contacts holding the ion",
    "ligand":               "a residue that contacts the ion",
    "donor":                "the atom that makes the contact",
    "ΔΔG":                  "the change in folding stability",
    "Rosetta energy units": "the arbitrary units of that stability score",
    "construct":            "the two- or three-domain piece that was modelled",
    "cofolded":             "modelled together with the ion",
    "seed":                 "one run of the structure predictor",
    "pLDDT":                "the predictor's own confidence score",
    "r.m.s.d.":             "how far atoms moved",
    "relative solvent accessibility": "how exposed a residue is",
    "consensus":            "the recurring sequence pattern",
    "PM1":                  "a specific evidence code in the classification guidelines",
    "odds ratio":           "how many times more likely",
    "Cliff's":              "an effect-size statistic",
    "Jaccard":              "an overlap score",
    "Wilson":               "a kind of confidence interval",
}

CLAUSE_MARKERS = re.compile(
    r"\b(although|whereas|because|since|while|which|whose|that|when|where|if|so that|"
    r"rather than|instead of|compared with|against|given that|such that|in which|"
    r"at which|by which|and|but|yet|however|therefore|nonetheless)\b", re.I)

LOAD_LIMIT = 6          # clause markers
WORD_LIMIT = 34         # words


def sentences(md: str) -> list[tuple[str, str]]:
    """(section, sentence) for prose only, skipping headings, tables and code."""
    out, section = [], "front matter"
    for para in md.split("\n\n"):
        p = para.strip()
        if not p or p.startswith(("|", "---", "<!--")):
            continue
        if p.startswith("#"):
            section = p.lstrip("# ").strip()
            continue
        p = re.sub(r"`[^`]*`", "CODE", p)
        p = re.sub(r"\s+", " ", p)
        for s in re.split(r"(?<=[.!?])\s+(?=[A-Z*(“])", p):
            s = s.strip()
            if len(s.split()) >= 4:
                out.append((section, s))
    return out


def load_of(s: str) -> tuple[int, int, int]:
    words = len(s.split())
    clauses = len(CLAUSE_MARKERS.findall(s))
    asides = s.count("(") + s.count(",")
    return words, clauses, asides


def main() -> int:
    src = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT
    md = src.read_text(encoding="utf-8")
    sents = sentences(md)

    # ---- LOAD ----------------------------------------------------------------------------
    scored = [(w, c, a, sec, s) for sec, s in sents for (w, c, a) in [load_of(s)]]
    overloaded = sorted([x for x in scored if x[1] > LOAD_LIMIT or x[0] > WORD_LIMIT],
                        key=lambda x: (-x[1], -x[0]))
    med_w = sorted(x[0] for x in scored)[len(scored) // 2]
    med_c = sorted(x[1] for x in scored)[len(scored) // 2]

    # ---- JARGON --------------------------------------------------------------------------
    flat = re.sub(r"\s+", " ", md)
    undefined = []
    for term, plain in sorted(TERMS.items()):
        m = re.search(re.escape(term), flat, re.I)
        if not m:
            continue
        before = flat[:m.start()]
        # counted as defined if the plain-language gloss, or an explicit definition verb,
        # appears within the 300 characters around first use
        window = flat[max(0, m.start() - 300):m.end() + 300]
        defined = bool(re.search(r"\b(defined as|which is|that is|meaning|a measure of|"
                                 r"we call|termed|i\.e\.)\b", window, re.I))
        n = len(re.findall(re.escape(term), flat, re.I))
        undefined.append((term, plain, n, defined,
                          len(before.split())))          # word position of first use

    never = [u for u in undefined if not u[3]]

    # ---- report --------------------------------------------------------------------------
    # a source outside the tree is legitimate -- an earlier draft is compared from a scratch
    # copy -- so the header falls back to the bare name instead of raising
    try:
        shown = src.relative_to(ROOT)
    except ValueError:
        shown = src.name
    L = [f"# Prose diagnostic — `{shown}`", "",
         "Measured with `scripts/62_prose_diagnostic.py`, applying the clarity diagnostic from "
         "the `academic-paper-composer` skill. Not an opinion about the writing, a count.", "",
         "## 1. Sentence load", "",
         f"- {len(scored)} sentences, median {med_w} words and {med_c} clause markers",
         f"- **{len(overloaded)} sentences exceed the limit** "
         f"(more than {LOAD_LIMIT} clause markers or more than {WORD_LIMIT} words)",
         f"- that is {100*len(overloaded)/len(scored):.0f}% of the paper", "",
         "### The worst offenders, heaviest first", ""]
    for w, c, a, sec, s in overloaded[:18]:
        L += [f"**{c} clauses, {w} words, {a} asides** — *{sec}*", "", f"> {s}", ""]

    L += ["## 2. Undefined vocabulary", "",
          f"- {len(undefined)} domain terms used, **{len(never)} never explained in plain "
          "words before first use**", "",
          "| term | uses | first use at word | plain wording it needs |",
          "|---|---|---|---|"]
    for term, plain, n, defined, pos in sorted(never, key=lambda x: x[4]):
        L += [f"| {term} | {n} | {pos} | {plain} |"]

    L += ["", "## 3. What this means", "",
          "The skill's test for jargon is whether a graduate student in an adjacent field "
          "could follow the text. On the count above the paper fails that test at "
          f"{len(never)} terms, and it asks the reader to hold more than {LOAD_LIMIT} clauses "
          f"in {len(overloaded)} sentences.", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"{len(scored)} sentences, median {med_w} words / {med_c} clauses")
    print(f"{len(overloaded)} overloaded ({100*len(overloaded)/len(scored):.0f}%)")
    print(f"{len(never)} of {len(undefined)} domain terms never defined before first use")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
