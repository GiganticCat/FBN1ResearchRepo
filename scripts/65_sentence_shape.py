#!/usr/bin/env python3
"""
65_sentence_shape.py — how a long sentence is BUILT, corpus against this draft.

`63_style_extract.py` counts how long our sentences are and how many clause markers they carry.
Both of those can sit inside the corpus range while the sentences are still built the wrong way,
which is what happened on 2026-09-06. The mentor read the draft and said the multi-clause
sentences still read as machine-written even though every aggregate matched.

Reading the five papers shows what the difference is. A published long sentence has **one main
clause** and gets its length from noun phrases, parentheses and a trailing participial phrase.

  Jensen 2009  "Superpositions of the individual domains cbEGF9, hyb2 and cbEGF10 onto the
                corresponding domains in cbEGF22-TB4-cbEGF23 (Figure 3B) showed higher degrees
                of similarity (rmsd = 0.829, 1.844, and 1.668, respectively), consistent with
                the conserved nature of the cbEGF fold and similarities between TB and hyb
                domain sequences."

One subject, one verb, two parentheses, one trailing adjectival phrase. What it never does is
bolt a second independent clause on with ", and the model is ...", which is the construction
this draft reaches for. So the measurements here are

  COORDINATION   ", and/but/so/yet/while" followed by a new subject and finite verb -- a second
                 independent clause welded to the first. The AI tell, counted directly.
  BALANCED PAIR  "X does A, and Y does B" and "not X but Y" -- the symmetrical two-part sentence.
  PARENTHESES    how much length comes from asides rather than from more clauses.
  TAIL           a trailing participial or adjectival phrase (", suggesting that", ", consistent
                 with", ", demonstrating") -- how the corpus attaches an inference cheaply.

Everything is reported for long sentences (>= 30 words) separately, because that is where the
two registers actually diverge; short sentences look alike either way.

Writes: results/sentence_shape.md
Run:    .venv/bin/python scripts/65_sentence_shape.py [source.md]
"""
from __future__ import annotations

import importlib.util
import re
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "resources/papers"
DEFAULT = ROOT / "results/manuscript_v4.md"
OUT = ROOT / "results/sentence_shape.md"

LONG = 30          # words, the band where the two registers diverge

# A second independent clause welded on. The subject slot allows a determiner, a pronoun or a
# capitalised or lower-case noun, and the verb slot the auxiliaries and reporting verbs these
# papers actually use. Deliberately conservative -- it under-counts rather than over-counts.
# The subject slot allows up to three words, because "a sequence-based predictor scores" and
# "the control arm of 24 produced" both slipped through a one-word subject and were the exact
# construction being counted.
COORD = re.compile(
    r",\s+(and|but|so|yet|while|whereas)\s+"
    r"(the|this|that|these|those|its|their|our|it|they|we|a|an|no|each|every|both|neither|"
    r"[A-Za-z][A-Za-z-]+)"
    r"(?:\s+[A-Za-z][A-Za-z0-9-]*){0,3}\s+"
    r"(is|are|was|were|has|have|had|does|do|did|would|will|can|could|may|might|should|must|"
    r"[a-z]+s|[a-z]+ed)\b", re.I)

# The "not X but Y" seesaw only. An earlier version also matched "and", which fired on
# "positions that are not ligands (n = 143), and the fourth is ..." -- an ordinary list, not a
# seesaw -- so the "and" arm is gone.
BALANCED = re.compile(r"\b(is|are|was|were|does|do|did)\s+not\s+[^.,]{0,60}\bbut\b"
                      r"|\bnot\s+whether\b[^.]{0,80}\bbut\s+whether\b", re.I)

# How the corpus attaches an inference without starting a new clause.
TAIL = re.compile(r",\s+(suggesting|indicating|demonstrating|showing|consistent with|"
                  r"reflecting|implying|supporting|in contrast to|compared with|"
                  r"presumably|possibly|respectively)\b", re.I)

SUBORD_OPEN = re.compile(r"^(although|whereas|while|because|since|when|where|if|given|"
                         r"in the absence|in contrast|to (understand|confirm|overcome|test))\b",
                         re.I)

# ---------------------------------------------------------------------------------------
# Clause load, counted ONE SENTENCE AT A TIME
#
# This is the measure the distribution checks hid. `63_style_extract.py` reported 12% of our
# sentences carrying 3+ clauses against 12% in the corpus, and the abstract still contained a
# four-clause sentence. Counting markers per sentence with the same instrument on both sides
# gives a ceiling instead of an average, and the ceiling is unambiguous:
#
#   1,164 published sentences -- 79.0% carry ONE marker, 18.3% TWO, 2.7% THREE, and NONE FOUR.
#
# So four is not "rare in the corpus", it is absent from it. A distribution that matches on
# average can still contain sentences no published paper contains, which is exactly what
# happened. "of which" is one marker and not two, since the preposition belongs to the relative.
# ---------------------------------------------------------------------------------------

LOAD_SUB = re.compile(r"\b(although|though|whereas|while|because|since|unless|when|where|"
                      r"whether|so that|such that|until)\b", re.I)
LOAD_REL = re.compile(r"\b(which|who|whom|whose)\b|\bthat\b(?=\s+\w+(s|ed|ing)\b)", re.I)


def clause_load(s: str) -> int:
    """Clause markers in one sentence, plus the main clause itself."""
    return (1 + len(LOAD_SUB.findall(s)) + len(LOAD_REL.findall(s))
            + len(COORD.findall(s)))


def shape(s: str) -> dict:
    return {
        "words": len(s.split()),
        "coord": len(COORD.findall(s)),
        "balanced": 1 if BALANCED.search(s) else 0,
        "parens": s.count("("),
        "commas": s.count(","),
        "tail": 1 if TAIL.search(s) else 0,
        "subord_open": 1 if SUBORD_OPEN.match(s) else 0,
        "load": clause_load(s),
    }


def profile(ss: list[str]) -> dict:
    long_ = [s for s in ss if len(s.split()) >= LONG]
    sh = [shape(s) for s in ss]
    lsh = [shape(s) for s in long_]
    n, nl = len(sh), max(len(lsh), 1)
    return {
        "n": n, "n_long": len(lsh),
        "pct_coord": 100 * sum(1 for x in sh if x["coord"]) / n,
        "pct_coord_long": 100 * sum(1 for x in lsh if x["coord"]) / nl,
        "coord_per_long": sum(x["coord"] for x in lsh) / nl,
        "pct_balanced": 100 * sum(x["balanced"] for x in sh) / n,
        "med_parens_long": st.median([x["parens"] for x in lsh]) if lsh else 0,
        "med_commas_long": st.median([x["commas"] for x in lsh]) if lsh else 0,
        "pct_tail": 100 * sum(x["tail"] for x in sh) / n,
        "pct_subord_open": 100 * sum(x["subord_open"] for x in sh) / n,
        "pct_load1": 100 * sum(1 for x in sh if x["load"] == 1) / n,
        "pct_load3": 100 * sum(1 for x in sh if x["load"] == 3) / n,
        "pct_load4": 100 * sum(1 for x in sh if x["load"] >= 4) / n,
        "max_load": max(x["load"] for x in sh),
    }


ROWS = [
    ("sentences", "n", "{:.0f}"),
    ("long sentences (30+ words)", "n_long", "{:.0f}"),
    ("% of ALL sentences welding a second clause", "pct_coord", "{:.1f}%"),
    ("% of LONG sentences welding a second clause", "pct_coord_long", "{:.1f}%"),
    ("welded clauses per long sentence", "coord_per_long", "{:.2f}"),
    ("% balanced pair / not-X-but-Y", "pct_balanced", "{:.1f}%"),
    ("median parentheses in a long sentence", "med_parens_long", "{:.1f}"),
    ("median commas in a long sentence", "med_commas_long", "{:.1f}"),
    ("% with a participial tail", "pct_tail", "{:.1f}%"),
    ("% opening on a subordinate clause", "pct_subord_open", "{:.1f}%"),
    ("% carrying ONE clause marker", "pct_load1", "{:.1f}%"),
    ("% carrying THREE", "pct_load3", "{:.1f}%"),
    ("% carrying FOUR or more", "pct_load4", "{:.1f}%"),
    ("heaviest sentence in the text", "max_load", "{:.0f}"),
]


def main() -> int:
    import pymupdf
    spec = importlib.util.spec_from_file_location("_se", ROOT / "scripts/63_style_extract.py")
    se = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(se)

    src = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT
    rows, pooled = {}, []
    for label, fn in se.CORPUS.items():
        d = pymupdf.open(PAPERS / fn)
        ss = se.sents(se.clean_pdf("\n".join(p.get_text() for p in d)))
        rows[label] = profile(ss)
        pooled += ss
    rows["— pooled —"] = profile(pooled)

    md = src.read_text(encoding="utf-8")
    body = "\n\n".join(b for b in md.split("\n\n")
                       if not b.strip().startswith(("#", "|", "---", "<!--", "**Fig", "**Supp")))
    body = re.sub(r"`[^`]*`", "CODE", body)
    body = re.sub(r"\*+", "", body)
    body = re.sub(r"\s+", " ", body)
    ours = se.sents(body)
    rows["THIS PAPER"] = profile(ours)

    L = [f"# Sentence shape — five published papers against `{src.name}`", "",
         "Produced by `scripts/65_sentence_shape.py`. `63_style_extract.py` measures how long a "
         "sentence is; this measures how it was built. A published long sentence carries one "
         "main clause and gets its length from noun phrases, parentheses and a trailing "
         "participial phrase. Welding a second independent clause on with `, and the model is "
         "...` is the construction that reads as machine-written.", "",
         "| measure | " + " | ".join(rows) + " |", "|---" * (len(rows) + 1) + "|"]
    for lab, k, fmt in ROWS:
        L.append(f"| {lab} | " + " | ".join(fmt.format(rows[r][k]) for r in rows) + " |")

    c, o = rows["— pooled —"], rows["THIS PAPER"]
    L += ["", "## Where this draft sits", ""]
    for lab, k, fmt in ROWS[2:]:
        if c[k] == 0:
            verdict = "**OVER**" if o[k] > 0 else "ok"
        else:
            r = o[k] / c[k]
            verdict = "**OVER**" if r > 1.35 else ("under" if r < 0.65 else "ok")
        L.append(f"- {lab} — corpus {fmt.format(c[k])}, ours {fmt.format(o[k])} — {verdict}")

    heavy = sorted(((shape(s), s) for s in ours if shape(s)["load"] >= 3),
                   key=lambda x: (-x[0]["load"], -x[0]["words"]))
    L += ["", f"## The {len(heavy)} sentences carrying three or more clause markers", "",
          "The corpus carries three in 2.7% of sentences and four in none of 1,164. Anything "
          "listed at four has to go.", ""]
    for sh_, s_ in heavy:
        L.append(f"**{sh_['load']} markers, {sh_['words']} words**")
        L += ["", f"> {s_}", ""]

    offenders = sorted(((shape(s), s) for s in ours
                        if shape(s)["coord"] or shape(s)["balanced"]),
                       key=lambda x: (-x[0]["coord"], -x[0]["words"]))
    L += ["", f"## The {len(offenders)} sentences that weld clauses, heaviest first", ""]
    for sh_, s in offenders:
        L.append(f"**{sh_['coord']} welded, {sh_['words']} words, {sh_['parens']} parens**")
        L += ["", f"> {s}", ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"corpus long sentences weld a clause {c['pct_coord_long']:.1f}% of the time, "
          f"ours {o['pct_coord_long']:.1f}%")
    print(f"over all sentences: corpus {c['pct_coord']:.1f}%, ours {o['pct_coord']:.1f}%")
    print(f"median parentheses in a long sentence: corpus {c['med_parens_long']:.1f}, "
          f"ours {o['med_parens_long']:.1f}")
    print(f"clause load: corpus {c['pct_load1']:.1f}/{c['pct_load3']:.1f}/{c['pct_load4']:.1f}"
          f" (one/three/four+), ours {o['pct_load1']:.1f}/{o['pct_load3']:.1f}/{o['pct_load4']:.1f}")
    print(f"heaviest sentence: corpus {c['max_load']:.0f} markers, ours {o['max_load']:.0f}")
    print(f"{len(offenders)} welded sentences listed in {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
