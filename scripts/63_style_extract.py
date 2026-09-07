#!/usr/bin/env python3
"""
63_style_extract.py — measure the sentence structure of the five reference papers and of ours.

The mentor's objection is that the prose "sounds like AI" and that some of its vocabulary is
invented rather than standard. Both are measurable against the papers themselves, so this script
does the `academic-paper-strategist` Phase 1 style extraction ("First-person vs passive voice
usage", "Use of technical terminology") as a count over the corpus in `resources/papers/` and
reports the delta against `results/manuscript_v3.md`.

Two questions it answers:

  STRUCTURE   how the reference papers build a sentence -- length, clause load, how a sentence
              opens, how often the authors say "we", how often they use the passive. The openers
              matter most. Starting a sentence with an adverbial or subordinate phrase instead of
              its subject is the habit that reads as machine-written when it is done constantly.

  TERMS       whether a term this paper uses appears anywhere in the five published papers. A
              term that appears nowhere in the corpus was invented here, and either has to be
              replaced with the standard term or defined at first use.

Writes: results/style_extraction.md
"""
from __future__ import annotations

import re
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "resources/papers"
MS = ROOT / "results/manuscript_v4.md"
OUT = ROOT / "results/style_extraction.md"

CORPUS = {
    "Jensen 2009 (Structure)":  "1-s2.0-S0969212609001610-main.pdf",
    "Kettle 1999 (JMB)":        "kettle1999.pdf",
    "Downing 1996 (Cell)":      "Downing_1996_Cell_cbEGF_pair.pdf",
    "Handford 1995 (JBC)":      "Handford_1995_JBC_cbEGF.pdf",
    "Godwin 2023 (NSMB)":       "Godwin_2023_NSMB_microfibril_cryoEM.pdf",
}

# Terms this manuscript leans on. Checked against the corpus, not against intuition.
OUR_TERMS = [
    "cbEGF", "calcium ligand", "Ca2+ ligand", "coordination", "bond-valence", "valence",
    "cognate site", "bystander", "composition-matched", "consensus motif", "motif",
    "donor", "side-chain oxygen", "main-chain carbonyl", "backbone carbonyl",
    "disulphide", "disulfide", "r.m.s.d.", "rmsd", "seed", "construct", "cofold",
    "pathogenic", "missense", "substitution", "affinity", "stability", "buried",
    "solvent accessibility", "beta-hydroxylation", "hydroxylation", "over-call",
    "non-ligand", "negative control", "threshold", "odds ratio",
]

SUBORD = ("although", "while", "whereas", "because", "since", "if", "when", "where",
          "after", "before", "given", "having", "using", "applied", "across", "among",
          "within", "under", "despite", "unlike", "in", "on", "at", "for", "with",
          "by", "to", "from", "as", "here", "thus", "therefore", "however", "moreover",
          "furthermore", "finally", "first", "second", "additionally", "notably",
          "importantly", "overall", "taken", "together", "based", "consistent")
CLAUSE = re.compile(r"\b(although|whereas|because|since|while|which|whose|that|when|where|"
                    r"if|so that|rather than|instead of|compared with|and|but|yet|however|"
                    r"therefore|nonetheless|in which|at which|by which)\b", re.I)
PASSIVE = re.compile(r"\b(was|were|is|are|been|being|be)\s+\w+(ed|en)\b", re.I)
FIRSTP = re.compile(r"\b(we|our|us)\b", re.I)


def clean_pdf(txt: str) -> str:
    """Strip the wreckage PDF text extraction leaves behind."""
    txt = txt.replace("­", "").replace("ﬁ", "fi").replace("ﬂ", "fl").replace("®", "fi")
    lines = []
    for ln in txt.splitlines():
        s = ln.strip()
        if not s or len(s) < 3:
            continue
        if re.match(r"^\W*\d+\W*$", s):                      # page numbers
            continue
        if re.search(r"(Structure 17|Cell 85|ª\d{4}|©\s?\d{4}|Elsevier|Academic Press|"
                     r"All rights reserved|doi:|DOI |www\.|http|Nature Structural)", s, re.I):
            continue
        if s.isupper() and len(s.split()) < 8:               # section headers
            continue
        if re.match(r"^(Figure|Fig\.|Table|Supplementary)\s", s):
            continue
        lines.append(s)
    t = " ".join(lines)
    t = re.sub(r"-\s+(?=[a-z])", "", t)                      # rejoin hyphenated line breaks
    t = re.sub(r"\([^)]{0,80}et al\.,? \d{4}[^)]{0,120}\)", "", t)   # author-year citations
    t = re.sub(r"\s+", " ", t)
    return t


def sents(t: str) -> list[str]:
    out = []
    for s in re.split(r"(?<=[.!?])\s+(?=[A-Z])", t):
        s = s.strip()
        w = s.split()
        # keep real prose only: sane length, mostly words, ends in a full stop
        if 6 <= len(w) <= 80 and s.endswith(".") and \
           sum(c.isalpha() or c.isspace() for c in s) / len(s) > 0.82:
            out.append(s)
    return out


def profile(ss: list[str]) -> dict:
    lens = [len(s.split()) for s in ss]
    cl = [len(CLAUSE.findall(s)) for s in ss]
    opener_sub = sum(1 for s in ss if s.split()[0].lower().strip(",") in SUBORD)
    opener_comma = sum(1 for s in ss if "," in s[:38])
    return {
        "n": len(ss),
        "median_words": st.median(lens),
        "mean_words": sum(lens) / len(lens),
        "p90_words": sorted(lens)[int(0.9 * len(lens)) - 1],
        "pct_over_34": 100 * sum(1 for x in lens if x > 34) / len(lens),
        "pct_under_15": 100 * sum(1 for x in lens if x < 15) / len(lens),
        "median_clauses": st.median(cl),
        "pct_3plus_clauses": 100 * sum(1 for x in cl if x >= 3) / len(cl),
        "pct_subord_opener": 100 * opener_sub / len(ss),
        "pct_early_comma": 100 * opener_comma / len(ss),
        "pct_passive": 100 * sum(1 for s in ss if PASSIVE.search(s)) / len(ss),
        "pct_first_person": 100 * sum(1 for s in ss if FIRSTP.search(s)) / len(ss),
    }


def main() -> int:
    import pymupdf

    rows, corpus_sents, corpus_text = {}, [], []
    for label, fn in CORPUS.items():
        d = pymupdf.open(PAPERS / fn)
        t = clean_pdf("\n".join(p.get_text() for p in d))
        ss = sents(t)
        rows[label] = profile(ss)
        corpus_sents += ss
        corpus_text.append(t.lower())
    rows["— all five pooled —"] = profile(corpus_sents)

    md = MS.read_text(encoding="utf-8")
    body = "\n\n".join(b for b in md.split("\n\n")
                       if not b.strip().startswith(("#", "|", "---", "<!--", "**Fig", "**Supp")))
    body = re.sub(r"`[^`]*`", "CODE", body)
    body = re.sub(r"\*+", "", body)
    body = re.sub(r"\s+", " ", body)
    ours = profile(sents(body))
    rows["THIS PAPER"] = ours

    keys = [("median_words", "median words", "{:.0f}"), ("mean_words", "mean words", "{:.1f}"),
            ("p90_words", "90th pct words", "{:.0f}"),
            ("pct_over_34", "% over 34 words", "{:.0f}%"),
            ("pct_under_15", "% under 15 words", "{:.0f}%"),
            ("median_clauses", "median clauses", "{:.0f}"),
            ("pct_3plus_clauses", "% with 3+ clauses", "{:.0f}%"),
            ("pct_subord_opener", "% not starting on the subject", "{:.0f}%"),
            ("pct_early_comma", "% with a comma in the first 38 chars", "{:.0f}%"),
            ("pct_passive", "% passive", "{:.0f}%"),
            ("pct_first_person", "% saying we/our", "{:.0f}%")]

    L = ["# Style extraction — five published papers against this manuscript", "",
         "Produced by `scripts/63_style_extract.py`. Prose sentences only, with captions, "
         "headers, tables and author-year citations stripped. This is the "
         "`academic-paper-strategist` Phase 1 style extraction done as a count.", "",
         "## 1. Sentence structure", "",
         "| measure | " + " | ".join(rows) + " |",
         "|---" * (len(rows) + 1) + "|"]
    for k, lab, fmt in keys:
        L.append(f"| {lab} | " + " | ".join(fmt.format(rows[r][k]) for r in rows) + " |")
    L += ["", f"Sentences sampled: " +
          ", ".join(f"{r} {rows[r]['n']}" for r in rows), ""]

    # ---- rhetorical habits ------------------------------------------------------------------
    # The length and opener counts above came out close to the corpus, so what reads as
    # machine-written is not sentence length. It is a handful of constructions used at many
    # times the rate the published papers use them.
    PATTERNS = {
        "appositive tail \", which \"":        r", which ",
        "sentence-final \", which is ...\"":   r", which is [^.]{0,60}\.$",
        "\"rather than\"":                     r"\brather than\b",
        "\"namely\"":                          r"\bnamely\b",
        "meta-commentary about the paper":  r"\b(what is new|matters more|the point is|"
                                            r"worth noting|we note|it is worth|"
                                            r"the relevant question|is itself informative|"
                                            r"should be read as)\b",
        "hedge stack (\"indicate that ... whether\")":
                                            r"\b(indicate|suggest)s? that\b[^.]{0,80}\bwhether\b",
        "\"is not X, it is Y\"":               r"\bis not\b[^.]{0,60}\b(but|rather|it is)\b",
    }
    L += ["## 2. Rhetorical habits", "",
          "Sentence length, openers, passive rate and first-person rate all sit inside the "
          "corpus range. These do not.", "",
          "| construction | five papers | this paper | ratio |", "|---|---|---|---|"]
    for lab, pat in PATTERNS.items():
        c = 100 * sum(1 for s_ in corpus_sents if re.search(pat, s_, re.I)) / len(corpus_sents)
        o = 100 * sum(1 for s_ in sents(body) if re.search(pat, s_, re.I)) / len(sents(body))
        ratio = "—" if c == 0 else f"{o / c:.0f}x"
        mark = " **" if (c == 0 and o > 0) or (c and o / c > 1.8 and o > 1.5) else ""
        L.append(f"| {lab} | {c:.1f}% | {o:.1f}%{mark} | {ratio if c else 'never in corpus'} |")
    L += ["", f"Also, {rows['THIS PAPER']['pct_3plus_clauses']:.0f}% of our sentences carry "
          f"three or more clauses against "
          f"{rows['— all five pooled —']['pct_3plus_clauses']:.0f}% in the corpus.", ""]

    # ---- terminology ----------------------------------------------------------------------
    joined = " ".join(corpus_text)
    ours_l = body.lower()
    L += ["## 3. Terminology — does the corpus actually use our words?", "",
          "| term | uses here | appears in the five papers |", "|---|---|---|"]
    invented = []
    for term in OUR_TERMS:
        n_ours = len(re.findall(re.escape(term.lower()), ours_l))
        if not n_ours:
            continue
        n_corpus = len(re.findall(re.escape(term.lower()), joined))
        L.append(f"| {term} | {n_ours} | " +
                 (f"yes, {n_corpus}" if n_corpus else "**NO — invented here**") + " |")
        if not n_corpus:
            invented.append((term, n_ours))
    L += ["", f"**{len(invented)} terms appear nowhere in the corpus**: "
          + ", ".join(f"{t} ({n})" for t, n in sorted(invented, key=lambda x: -x[1])), ""]

    OUT.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"corpus {rows['— all five pooled —']['n']} sentences, ours {ours['n']}")
    for k, lab, fmt in keys:
        print(f"  {lab:<38} corpus {fmt.format(rows['— all five pooled —'][k]):>7}"
              f"   ours {fmt.format(ours[k]):>7}")
    print(f"\ninvented terms: {[t for t, _ in invented]}")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
