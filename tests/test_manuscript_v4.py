#!/usr/bin/env python3
"""
test_manuscript_v4.py — gate on the rewritten, citation-free draft.

`results/manuscript_v4.md` was rewritten on 2026-09-05 against the register of five published
papers (`results/style_reference.md`) after the mentor rejected the previous draft's vocabulary
and prose. This gate holds the four constraints that rewrite was made under, so that none of them
can be reintroduced by a later edit, and it re-checks that the rewrite changed no number.

  V1  no colons or semicolons anywhere in the draft, URLs and DOIs excepted
  V2  no citation markers and no reference list
  V3  every claim-bearing value in the keyed manuscript survives verbatim in the rewrite
  V4  no short dramatic sentence of the form the mentor rejected
  V5  the figure vocabulary in the scripts matches the agreed terms
  V6  the built DOCX and PDF exist, are newer than the source, and carry the six figures
  V7  sentence load and length stay inside the range measured from the five reference papers
  V8  no term absent from those five papers unless it is on the allowed list
  V9  the title and the Results headings are neutral noun phrases
  V10 sentences are BUILT the way the five papers build them, and none exceeds the heaviest
      sentence in the 1,164 published ones (a ceiling, because an average hides one bad sentence)
  V11 American spelling throughout, including anything drawn into a figure
  V12 no sentence announces a result instead of stating it

Run: .venv/bin/python tests/test_manuscript_v3.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/manuscript_v4.md"
KEYED = ROOT / "results/manuscript_paper.md"
DOCX = ROOT / "results/FBN1_manuscript_v4.docx"
PDF = ROOT / "results/FBN1_manuscript_v4.pdf"

results: list[tuple[bool, str, str]] = []


def check(ok: bool, name: str, detail: str) -> None:
    results.append((ok, name, detail))


# Values that carry a claim. If the rewrite dropped or altered one of these, it changed the
# science rather than the prose, which is the one thing it was not allowed to do.
CLAIM_VALUES = [
    "4,451", "578", "823", "425", "3,831", "751", "17%", "86%", "143", "271", "323", "1,132",
    "10.81", "1.77", "1.84", "0.42", "21.7", "38.4", "1.72", "4.14", "1.92", "2.01", "0.54",
    "0.989", "0.998", "0.260", "9.1%", "123", "109", "45%", "−5.7", "−1.7", "−38.0", "28.5",
    "0.29", "0.24", "0.80", "29.8", "12 of 16", "3 of 17", "33.0", "0.47", "11.9", "−1.0",
    "−12.2", "−14.6", "−17.2", "N589S", "N1046S", "N2624S", "18 of 66", "0.13", "−0.09", "0.62",
    "0.982", "0.61", "3.41", "2.21", "5.27", "0.049", "0.017", "0.212", "1.49", "0.90", "2.51",
    "0.45", "1.35", "0.55", "0.961", "0.147", "31.6", "1.44", "24.76", "21.5", "0.75",
    "1.6 mM", "2,871", "264", "1.967", "2.14", "0.37", "0.35", "2,000", "129", "2026-07-28",
    "NM_000138.5", "P35555", "2W86", "1LMJ", "1UZJ", "1EMN", "1EMO", "1,901",
    "9 of 13",
]

# Claims the rewrite may state in words instead of digits. Either form passes.
CLAIM_ALTERNATIVES = {"0 of 12": ["none of 12"]}

FIGURE_TERMS = {
    "scripts/paperstyle.py": ["Ca²⁺ ligand\\n(side chain)", "Non-ligand\\nAsp/Glu/Asn",
                              "Consensus motif,\\nnon-ligand", "Cysteine\\nsubstitution"],
    "scripts/52_fig3_coordination.py": ["Donor\\nretained", "Donor\\nlost"],
    "scripts/53_fig4_interpretation.py": ["coordination reduced", "coordination retained"],
}
BANNED_FIGURE_TERMS = ["Matches motif", "Coordinates\\nCa²⁺", "Cysteine\\nremoved",
                       "Keeps an\\nO donor", "Site disrupted", "The fold, same models"]


def main() -> int:
    print("manuscript v4 gate — the rewritten, citation-free draft\n")
    src = SRC.read_text(encoding="utf-8")

    # ---------------------------------------------------------------- V1 punctuation
    # The house rule is about prose punctuation, so a web address or a DOI in the availability
    # statement is not an offence. Those are stripped before the line is checked, which lets the
    # mentor paste the real repository link in without the gate turning red.
    def strip_links(ln: str) -> str:
        return re.sub(r"(https?://\S+|\bdoi\.org/\S+|\bdoi\s*10\.\S+)", "LINK", ln)

    offenders = [(i, ln) for i, ln in enumerate(src.splitlines(), 1)
                 if ":" in strip_links(ln) or ";" in strip_links(ln)]
    check(not offenders, "V1 no colons or semicolons",
          "clean" if not offenders
          else f"{len(offenders)} lines, first at {offenders[0][0]}: {offenders[0][1][:60]}")

    # ---------------------------------------------------------------- V2 no citations
    markers = re.findall(r"\[@[^\]]*\]", src)
    has_reflist = bool(re.search(r"^##+\s*References", src, re.M))
    check(not markers and not has_reflist, "V2 no citations",
          f"{len(markers)} [@key] markers, reference list present: {has_reflist}")

    # ---------------------------------------------------------------- V3 numbers preserved
    keyed = KEYED.read_text(encoding="utf-8")
    flat_src = re.sub(r"\s+", " ", src)
    flat_keyed = re.sub(r"\s+", " ", re.sub(r"\[@[^\]]*\]", " ", keyed))
    def present(v: str) -> bool:
        return v in flat_src or any(a in flat_src for a in CLAIM_ALTERNATIVES.get(v, []))

    lost = [v for v in CLAIM_VALUES + list(CLAIM_ALTERNATIVES)
            if v in flat_keyed and not present(v)]
    absent = [v for v in CLAIM_VALUES if v not in flat_keyed]
    check(not lost and not absent, "V3 claim values preserved",
          f"{len(CLAIM_VALUES) + len(CLAIM_ALTERNATIVES)} values checked against the keyed source, "
          f"{len(lost)} lost, {len(absent)} not in the keyed source"
          + (f" — lost {lost[:4]}" if lost else "")
          + (f" — absent {absent[:4]}" if absent else ""))

    # ---------------------------------------------------------------- V4 no dramatic fragments
    # The rejected pattern is a very short sentence that reverses the one before it. Guard the
    # reversal verbs directly rather than sentence length, which has legitimate short cases.
    body = re.sub(r"`[^`]*`", "X", src)
    body = re.sub(r"\s+", " ", body)
    banned = [r"\.\s+It is not\.", r"\.\s+It does not\.", r"\.\s+They are not\.",
              r"\.\s+Neither does\.", r"\.\s+It cannot\.", r"\.\s+They do not\.",
              r"\.\s+We found that they aren", r"\.\s+Not so\."]
    hits = [p for p in banned if re.search(p, body)]
    check(not hits, "V4 no dramatic reversals", "clean" if not hits else f"found {hits}")

    # ---------------------------------------------------------------- V5 figure vocabulary
    missing, present_banned = [], []
    for rel, terms in FIGURE_TERMS.items():
        t = (ROOT / rel).read_text(encoding="utf-8")
        missing += [f"{rel}:{x}" for x in terms if x not in t]
    for rel in list(FIGURE_TERMS) + ["scripts/51_fig2_stability.py"]:
        # the rejected labels survive in the comments that record why they were replaced, and
        # that is the point of those comments, so only code lines are checked
        code = "\n".join(ln for ln in (ROOT / rel).read_text(encoding="utf-8").splitlines()
                          if not ln.lstrip().startswith("#"))
        present_banned += [f"{rel}:{x}" for x in BANNED_FIGURE_TERMS if f'"{x}' in code]
    ok = not missing and not present_banned
    check(ok, "V5 figure vocabulary",
          "agreed terms present, none of the rejected labels drawn" if ok
          else f"missing {missing[:3]} banned {present_banned[:3]}")

    # ---------------------------------------------------------------- V6 built documents
    det = []
    ok = True
    for p in (DOCX, PDF):
        if not p.is_file():
            ok = False
            det.append(f"{p.name} missing")
        elif p.stat().st_mtime < SRC.stat().st_mtime:
            ok = False
            det.append(f"{p.name} older than the source")
        else:
            det.append(f"{p.name} {p.stat().st_size/1024:.0f} kB")
    if PDF.is_file():
        try:
            import pymupdf
            d = pymupdf.open(PDF)
            n_img = sum(len(d[i].get_images()) for i in range(d.page_count))
            det.append(f"{d.page_count} pages, {n_img} images")
            if n_img < 6:
                ok = False
        except ImportError:
            det.append("pymupdf unavailable, image count not checked")
    check(ok, "V6 documents built", "; ".join(det).replace(";", ","))

    # ---------------------------------------------------------------- V7 corpus style
    # Targets measured from the five reference papers by scripts/63_style_extract.py. The
    # tolerances are deliberately one-sided where a lower value is not a fault: shorter
    # sentences and fewer clauses than the corpus are safe, longer and more are the failure
    # the mentor rejected twice.
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "scripts/63_style_extract.py")],
                       capture_output=True, text=True, cwd=ROOT)
    out = r.stdout
    def grab(label: str) -> float:
        m = re.search(rf"{re.escape(label)}\s+corpus\s+([\d.]+)%?\s+ours\s+([\d.]+)%?", out)
        return (float(m.group(1)), float(m.group(2))) if m else (float("nan"),) * 2

    checks, style_ok = [], True
    for label, cap in [("% with 3+ clauses", 1.4), ("% over 34 words", None)]:
        c, o = grab(label)
        if cap and o > c * cap:
            style_ok = False
            checks.append(f"{label} {o:.0f}% vs corpus {c:.0f}% OVER")
        else:
            checks.append(f"{label} {o:.0f}% vs {c:.0f}%")
    med_c, med_o = grab("median words")
    if med_o > med_c * 1.15:
        style_ok = False
    checks.append(f"median words {med_o:.0f} vs {med_c:.0f}")
    check(style_ok, "V7 style within corpus range", ", ".join(checks))

    # ---------------------------------------------------------------- V8 no invented terms
    # Terms absent from all five reference papers. Six are allowed, because they come from
    # crystallography, statistics and structure prediction and simply do not arise in these
    # papers. Anything else absent from the corpus was coined here.
    ALLOWED = {"bond-valence", "valence", "seed", "cofold", "negative control", "odds ratio"}
    m = re.search(r"invented terms: \[(.*?)\]", out)
    found = {x.strip().strip("'\"") for x in m.group(1).split(",")} if m and m.group(1) else set()
    coined = found - ALLOWED
    check(not coined, "V8 no coined terminology",
          f"{len(found)} terms absent from the corpus, all expected"
          if not coined else f"coined here: {sorted(coined)}")

    # ---------------------------------------------------------------- V9 heading register
    # The mentor asked on 2026-09-06 for the traditional-journal register, where a heading names
    # its subject and does not state its finding. Jensen 2009 heads sections "Structure of the
    # cbEGF domain pair"; the set this replaced read "Cysteine substitutions destabilise the
    # domain and calcium-ligand substitutions do not". A parser is overkill for eight strings, so
    # this is a tripwire on the finite verbs the rejected set used, plus the give-away openers.
    FINITE = ("concentrate", "reproduce", "destabilise", "destabilises", "depends", "holds",
              "distinguishes", "is reduced", "are enriched", "shows", "show", "reveals",
              "reveal", "does not", "do not", "cannot", "loosen", "loosens", "preserve",
              "preserves", "fails", "fail", "detects", "detect")
    OPENERS = ("neither ", "the outcome ", "pathogenic variants ", "modelled ")
    heads = [ln[2:].strip() for ln in src.splitlines() if ln.startswith("# ")]
    heads += [ln[4:].strip() for ln in src.splitlines() if ln.startswith("### ")]
    offending = sorted({f"{h!r} ({w})" for h in heads for w in FINITE
                        if re.search(rf"\b{re.escape(w)}\b", h, re.I)}
                       | {f"{h!r} (opens {o.strip()!r})" for h in heads for o in OPENERS
                          if h.lower().startswith(o)})
    check(not offending, "V9 headings are noun phrases",
          f"{len(heads)} headings, all noun phrases" if not offending
          else f"{len(offending)} state a finding: {offending[:2]}")

    # ---------------------------------------------------------------- V10 sentence shape
    # V7 checks how long a sentence is. This checks how it was built, which is the thing V7
    # cannot see: on 2026-09-06 every V7 aggregate was inside the corpus range while 28.6% of
    # sentences still bolted a second independent clause on with ", and the model is ...",
    # against 4.0% in the five papers. `65_sentence_shape.py` measures that construction, plus
    # the three devices the corpus uses instead -- a parenthesis, a participial tail and a
    # subordinate opener. Ceilings only where more is the fault; floors where less is.
    r10 = subprocess.run([sys.executable, str(ROOT / "scripts/65_sentence_shape.py")],
                         capture_output=True, text=True, cwd=ROOT)
    o10 = r10.stdout

    def grab10(pat: str) -> tuple[float, float]:
        m = re.search(pat, o10)
        return (float(m.group(1)), float(m.group(2))) if m else (float("nan"),) * 2

    shape_ok, notes = True, []
    c_all, o_all = grab10(r"over all sentences: corpus ([\d.]+)%, ours ([\d.]+)%")
    if o_all > c_all * 1.35:
        shape_ok = False
    notes.append(f"welded {o_all:.1f}% vs corpus {c_all:.1f}%")
    c_long, o_long = grab10(r"long sentences weld a clause ([\d.]+)% of the time, ours ([\d.]+)%")
    if o_long > c_long * 1.35:
        shape_ok = False
    notes.append(f"welded in long {o_long:.1f}% vs {c_long:.1f}%")
    c_par, o_par = grab10(r"median parentheses in a long sentence: corpus ([\d.]+), ours ([\d.]+)")
    if o_par < c_par:
        shape_ok = False
    notes.append(f"parens per long sentence {o_par:.0f} vs {c_par:.0f}")

    # A CEILING, not an average. The distribution checks above and V7 both passed on a draft
    # whose abstract carried a four-clause sentence, because an average cannot see a single bad
    # sentence. No sentence in the 1,164 published ones carries four clause markers, so no
    # sentence here may either.
    m3 = re.search(r"clause load: corpus ([\d.]+)/([\d.]+)/([\d.]+) \(one/three/four\+\), "
                   r"ours ([\d.]+)/([\d.]+)/([\d.]+)", o10)
    mx = re.search(r"heaviest sentence: corpus (\d+) markers, ours (\d+)", o10)
    if m3 and mx:
        c1, c3, c4 = (float(m3.group(i)) for i in (1, 2, 3))
        n1, n3, n4 = (float(m3.group(i)) for i in (4, 5, 6))
        c_max, o_max = int(mx.group(1)), int(mx.group(2))
        if o_max > c_max or n4 > c4:
            shape_ok = False
        if n3 > c3 * 1.35:
            shape_ok = False
        notes.append(f"clause load one/three/four+ {n1:.0f}/{n3:.1f}/{n4:.1f}% vs "
                     f"{c1:.0f}/{c3:.1f}/{c4:.1f}%, heaviest {o_max} vs {c_max}")
    else:
        shape_ok = False
        notes.append("clause load not reported")
    check(shape_ok, "V10 sentences built like the corpus", ", ".join(notes))

    # ---------------------------------------------------------------- V11 spelling register
    # British until 2026-09-06, American from then on by the mentor's call. The reference papers
    # are mostly British journals, so the pull is towards -ise and disulphide every time a
    # sentence is rewritten against them. `disulfide` is safe for V8 because Jensen 2009 and
    # Handford 1995 both use it. `grey` in the figure scripts is the name of a colour-key entry,
    # `P.C["grey"]`, and is never drawn, so only strings that reach a figure are checked.
    BRITISH = [
        "modelled", "modelling", "labelled", "signalling", "characterised", "crystallised",
        "destabilise", "destabilised", "destabilises", "destabilising", "normalised",
        "minimisation", "summarising", "recognises", "optimises", "disulphide", "disulphides",
        "sulphur", "haemophilia", "colour", "colours", "behaviour", "neighbours",
        "neighbouring", "neighbourhood", "artefact", "artefacts", "grey", "greyscale",
        "towards", "judgement", "centre", "fibre", "whilst", "amongst", "programme",
        "analyse", "analysed", "analysing", "utilise", "organised", "standardised",
        "stabilise", "stabilised", "hypothesised",
    ]
    found = sorted({w for w in BRITISH if re.search(rf"\b{w}\b", src, re.I)})

    # plus anything a figure script actually draws
    drawn = re.compile(r"(set_xlabel|set_ylabel|set_title|annotate|P\.title|label_at|"
                       r"ticklabels)\s*\([^)]{0,400}", re.S)
    fig_hits = []
    for f in sorted((ROOT / "scripts").glob("*.py")):
        t = f.read_text(encoding="utf-8")
        for blk in drawn.finditer(t):
            for lit in re.findall(r'"([^"\n]{3,160})"', blk.group(0)):
                fig_hits += [f"{f.name}:{w}" for w in BRITISH
                             if re.search(rf"\b{w}\b", lit, re.I)]
    fig_hits = sorted(set(fig_hits))
    check(not found and not fig_hits, "V11 American spelling",
          "no British forms in the draft or in any drawn figure label"
          if not found and not fig_hits
          else f"draft {found[:4]}, figures {fig_hits[:3]}")

    # ---------------------------------------------------------------- V12 no scaffolding
    # The academic-paper-composer checklist reduces to "specific, not vague", and the five papers
    # obey it almost absolutely. A sentence that names a category without giving its content --
    # "The outcome was decided by the replacing residue." -- occurs twice in 1,164 published
    # sentences (0.2%), and both of those are real claims. This draft carried nine (3.6%), each
    # followed immediately by the sentence that said the thing.
    #
    # A sentence is caught when it is short, carries no number and no named entity, and its verb
    # is one of the empty ones. The bolded lead of a Limitations paragraph is exempt, because
    # that section states a claim and then evidences it, which is its deliberate structure.
    VAGUE_V = re.compile(
        r"\b(decided|determined|depends?|concerns?|lies?|matters?|differ|differs|is testable|"
        r"is systematic|are systematic|is a floor|varies|gives a different|is not equivalent|"
        r"is thin|is small|is modeled|rests? on|bears? on|follows?|applies|implication)\b", re.I)
    NAMED = re.compile(
        r"\d|\b(cbEGF|ClinVar|gnomAD|FoldX|Rosetta|AlphaMissense|AlphaFold|UniProt|CheckMyMetal|"
        r"PDB|SIFTS|Ca|Asp|Glu|Asn|Cys|Ser|His|Tyr|PM1|REVEL|ClinGen|FBN1|Mutalyzer|FreeSASA|"
        r"DSSP|PyMOL|Å|MFS|TGF|rmsd|ΔΔG)\b")

    import importlib.util as _il
    _sp = _il.spec_from_file_location("_se12", ROOT / "scripts/63_style_extract.py")
    _se = _il.module_from_spec(_sp)
    _sp.loader.exec_module(_se)

    stub, n_sent, section = [], 0, "front"
    for para in src.split("\n\n"):
        t = para.strip()
        if not t or t.startswith(("|", "---", "<!--")):
            continue
        if t.startswith("#"):
            section = t.lstrip("# ").strip()
            continue
        if section in ("Figure legends", "Declarations", "Data and code availability"):
            continue
        bold_lead = t.startswith("**")
        flat = re.sub(r"\s+", " ", re.sub(r"\*+", "", re.sub(r"`[^`]*`", "CODE", t)))
        for i, sent in enumerate(_se.sents(flat)):
            n_sent += 1
            if bold_lead and i == 0:
                continue
            if len(sent.split()) <= 16 and not NAMED.search(sent) and VAGUE_V.search(sent):
                stub.append(sent)
    rate = 100 * len(stub) / max(n_sent, 1)
    # 0.2% in the corpus. Two are allowed here because two of ours are genuine claims about the
    # literature rather than throat-clears, which is what the corpus's own two are.
    check(len(stub) <= 2, "V12 no announce-without-stating",
          f"{len(stub)}/{n_sent} = {rate:.1f}% against 0.2% in the five papers"
          + ("" if len(stub) <= 2 else f" — {stub[:2]}"))

    # ---------------------------------------------------------------- report
    for good, name, detail in results:
        print(f"  {'PASS' if good else 'FAIL'}  {name:<32} {detail}")
    n_ok = sum(1 for g, _, _ in results if g)
    print(f"\n{n_ok}/{len(results)} checks passed")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
