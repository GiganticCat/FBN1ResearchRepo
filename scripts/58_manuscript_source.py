#!/usr/bin/env python3
"""
58_manuscript_source.py -- build the keyed manuscript source and the final reference database.

`results/manuscript_final.md` carries its citations as literal superscript digits tied to a
hand-numbered list. That is fine until a reference is inserted, at which point every superscript
after the insertion point is wrong and nothing detects it. This script converts that file once
into a source that cites by KEY -- `[@kettle1999]` -- and emits the references as data. The docx
builder then assigns numbers by order of first appearance, so a citation added anywhere renumbers
the whole paper correctly and silently.

It also folds in the twenty references added after `57_references.py` verified the original list,
and applies the two title corrections that verification turned up:

  * Godwin 2023 -- the cited title was truncated before "affecting a key regulatory latent
    TGF-beta-binding site", which is the half of the title that names what the paper is about.
  * Delhomme 2022 -- the cited title was a paraphrase, not the published title.

Outputs
  results/manuscript_paper.md          the keyed source, prose byte-identical to the original
  data/processed/references_final.tsv  63 references, every field from Crossref or Europe PMC

Usage
  .venv/bin/python scripts/58_manuscript_source.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/manuscript_final.md"
OUT_MD = ROOT / "results/manuscript_paper.md"
OUT_TSV = ROOT / "data/processed/references_final.tsv"

# Keys for the original 43, in the order they appear in that file's reference list.
KEYS_43 = [
    "sakai1986", "dietz1991", "handford1995", "jensen2012", "downing1996", "knott1996",
    "godwin2023", "baudhuin2019", "drackley2024", "handford1991", "kettle1999", "reinhardt2000",
    "mcgettrick2000", "whiteman2003", "whiteman2007", "haller2020", "rao1995", "lee2004",
    "smallridge2003", "jensen2009", "landrum2018", "karczewski2020", "chen2024", "abramson2024",
    "schymkowitz2005", "alford2017", "cheng2023", "brown1985", "zheng2017", "pejaver2022",
    "dana2019", "uniprot2023", "mitternacht2016", "kabsch1983", "richards2015", "loeys2010",
    "milewicz2021", "jumper2021", "varadi2024", "hekkelman2023", "asano2022", "delhomme2022",
    "vallat2026",
]

# Corrections from `results/reference_verification.md`. Everything else in the original list
# checked out; the year flags there are all the NAR/Nature "online year vs issue year" split,
# where the issue year is the correct citation and is what the file already had.
TITLE_FIXES = {
    "godwin2023": (
        "Fibrillin microfibril structure identifies long-range effects of inherited "
        "pathogenic mutations.",
        "Fibrillin microfibril structure identifies long-range effects of inherited "
        "pathogenic mutations affecting a key regulatory latent TGFβ-binding site."),
    "delhomme2022": (
        "Mitral valve disease in Marfan syndrome: genotype–phenotype correlations.",
        "Genotype–mitral valve phenotype correlations in Marfan syndrome with *FBN1* "
        "pathogenic variants."),
}

# The twenty added here. Every field was resolved through Crossref and Europe PMC by
# `57_references.py`; see `data/interim/reference_candidates_verified.tsv`. Where the two
# services disagree on year the journal issue year is used, as it is for the original list.
NEW_REFS = [
    ("ioannidis2016", "Ioannidis NM, *et al.*",
     "REVEL: an ensemble method for predicting the pathogenicity of rare missense variants",
     "Am J Hum Genet", 2016, "99", "877–885", "10.1016/j.ajhg.2016.08.016"),
    ("benjamini1995", "Benjamini Y, Hochberg Y",
     "Controlling the false discovery rate: a practical and powerful approach to multiple testing",
     "J R Stat Soc B", 1995, "57", "289–300", "10.1111/j.2517-6161.1995.tb02031.x"),
    ("mann1947", "Mann HB, Whitney DR",
     "On a test of whether one of two random variables is stochastically larger than the other",
     "Ann Math Stat", 1947, "18", "50–60", "10.1214/aoms/1177730491"),
    ("wilson1927", "Wilson EB",
     "Probable inference, the law of succession, and statistical inference",
     "J Am Stat Assoc", 1927, "22", "209–212", "10.1080/01621459.1927.10502953"),
    ("cliff1993", "Cliff N",
     "Dominance statistics: ordinal analyses to answer ordinal questions",
     "Psychol Bull", 1993, "114", "494–509", "10.1037/0033-2909.114.3.494"),
    ("harris2020", "Harris CR, *et al.*", "Array programming with NumPy",
     "Nature", 2020, "585", "357–362", "10.1038/s41586-020-2649-2"),
    ("virtanen2020", "Virtanen P, *et al.*",
     "SciPy 1.0: fundamental algorithms for scientific computing in Python",
     "Nat Methods", 2020, "17", "261–272", "10.1038/s41592-019-0686-2"),
    ("hunter2007", "Hunter JD", "Matplotlib: a 2D graphics environment",
     "Comput Sci Eng", 2007, "9", "90–95", "10.1109/MCSE.2007.55"),
    ("touw2015", "Touw WG, *et al.*", "A series of PDB-related databanks for everyday needs",
     "Nucleic Acids Res", 2015, "43", "D364–D368", "10.1093/nar/gku1028"),
    ("nivon2013", "Nivón LG, Moretti R, Baker D",
     "A Pareto-optimal refinement method for protein design scaffolds",
     "PLoS One", 2013, "8", "e59004", "10.1371/journal.pone.0059004"),
    ("harding2006", "Harding MM",
     "Small revisions to predicted distances around metal sites in proteins",
     "Acta Crystallogr D", 2006, "62", "678–682", "10.1107/S0907444906014594"),
    ("morales2022", "Morales J, *et al.*",
     "A joint NCBI and EMBL-EBI transcript set for clinical genomics and research",
     "Nature", 2022, "604", "310–315", "10.1038/s41586-022-04558-8"),
    ("dendunnen2016", "den Dunnen JT, *et al.*",
     "HGVS recommendations for the description of sequence variants: 2016 update",
     "Hum Mutat", 2016, "37", "564–569", "10.1002/humu.22981"),
    ("wildeman2008", "Wildeman M, van Ophuizen E, den Dunnen JT, Taschner PEM",
     "Improving sequence variant descriptions in mutation databases and literature using the "
     "Mutalyzer sequence variation nomenclature checker",
     "Hum Mutat", 2008, "29", "6–13", "10.1002/humu.20654"),
    ("rehm2015", "Rehm HL, *et al.*", "ClinGen — the Clinical Genome Resource",
     "N Engl J Med", 2015, "372", "2235–2242", "10.1056/NEJMsr1406261"),
    ("rees1988", "Rees DJG, *et al.*",
     "The role of β-hydroxyaspartate and adjacent carboxylate residues in the first EGF domain "
     "of human factor IX",
     "EMBO J", 1988, "7", "2053–2061", "10.1002/j.1460-2075.1988.tb03045.x"),
    ("neptune2003", "Neptune ER, *et al.*",
     "Dysregulation of TGF-β activation contributes to pathogenesis in Marfan syndrome",
     "Nat Genet", 2003, "33", "407–411", "10.1038/ng1116"),
    ("ramirez2009", "Ramirez F, Sakai LY", "Biogenesis and function of fibrillin assemblies",
     "Cell Tissue Res", 2009, "339", "71–82", "10.1007/s00441-009-0822-x"),
    ("collodberoud2003", "Collod-Béroud G, *et al.*",
     "Update of the UMD-FBN1 mutation database and creation of an FBN1 polymorphism database",
     "Hum Mutat", 2003, "22", "199–208", "10.1002/humu.10249"),
    ("wojdyr2022", "Wojdyr M", "GEMMI: a library for structural biology",
     "J Open Source Softw", 2022, "7", "4200", "10.21105/joss.04200"),
]

SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
       "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}
# Three real corruptions were found by inspecting the output of this rule, and each guard below
# fixes one of them. Superscripts in this file mean two different things and the pattern has to
# tell them apart:
#
#   "P < 10⁻⁴⁷", "2 × 10⁻⁴", "kcal mol⁻¹"  -- exponents. A run that starts straight after a
#       superscript minus, or after another superscript digit, is part of an exponent, never a
#       citation. Without the lookbehind, 10⁻⁴ became citation 4 in six P values, and 10⁻⁴⁷ lost
#       its ⁷ to citation 7.
#   "Ca²⁺"                                 -- a charge. A run followed by a superscript sign is
#       a charge, never a citation. Without the lookahead, every Ca²⁺ in Methods became a
#       citation of Dietz 1991.
#   "gnomAD v4²³", "AlphaMissense v1.0.0²⁷" -- genuine citations that happen to follow an ASCII
#       digit, so ASCII digits must NOT be in the lookbehind or these are silently dropped.
SUP_RUN = re.compile(r"(?<![⁻⁰¹²³⁴⁵⁶⁷⁸⁹])[⁰¹²³⁴⁵⁶⁷⁸⁹][⁰¹²³⁴⁵⁶⁷⁸⁹⁻, ]*(?![⁺⁻])")


def parse_original(text: str) -> list[dict]:
    """Split the existing reference list into (authors, title, journal, year, vol, page, doi)."""
    start = text.index("## References")
    tail = text[start + len("## References"):]
    end = tail.find("\n**Note on sources.")
    block = tail[: end if end != -1 else len(tail)]

    entries, cur = [], None
    for raw in block.splitlines():
        line = raw.strip()
        if not line:
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", line)
        if m and (cur is None or int(m.group(1)) == cur[0] + 1):
            if cur:
                entries.append(cur)
            cur = [int(m.group(1)), m.group(2)]
        elif cur:
            cur[1] += " " + line
    if cur:
        entries.append(cur)

    # The formatted string is carried through verbatim rather than split into fields and
    # reassembled. Splitting it is where a reference list gets quietly mangled: author lists
    # contain "*et al.*", titles contain colons, gene names are italicised inside titles, and
    # every one of those breaks a field regex in a way that looks plausible in the output.
    # These strings are already correct and already verified against Crossref by script 57.
    out = []
    for n, body in entries:
        doi = ""
        md = re.search(r"doi:\s*(10\.\S+?)\s*$", body)
        if md:
            doi = md.group(1).rstrip(".")
            body = body[: md.start()].strip()
        out.append({"n": n, "formatted": body.rstrip().rstrip("."), "doi": doi})
    return out


def superscript_to_keys(body: str, num2key: dict[int, str]) -> tuple[str, int]:
    """Replace every superscript citation run with a bracketed key list."""
    n_sub = 0

    def expand(run: str) -> list[int]:
        run = run.replace(" ", "")
        nums = []
        for part in run.split(","):
            part = part.strip()
            if not part:
                continue
            if "⁻" in part:
                a, b = part.split("⁻")
                a = int("".join(SUP[c] for c in a))
                b = int("".join(SUP[c] for c in b))
                nums.extend(range(a, b + 1))
            else:
                nums.append(int("".join(SUP[c] for c in part)))
        return nums

    def repl(m: re.Match) -> str:
        nonlocal n_sub
        try:
            nums = expand(m.group(0))
        except (KeyError, ValueError):
            return m.group(0)
        if not nums or any(x not in num2key for x in nums):
            return m.group(0)
        n_sub += 1
        return "[@" + ",".join(num2key[x] for x in nums) + "]"

    body = SUP_RUN.sub(repl, body)

    # the one prose-form citation: "(refs. 3,4)"
    def refs_repl(m: re.Match) -> str:
        nonlocal n_sub
        nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
        if any(x not in num2key for x in nums):
            return m.group(0)
        n_sub += 1
        return "(refs. [@" + ",".join(num2key[x] for x in nums) + "])"

    body = re.sub(r"\(refs?\.\s*([\d,\s–-]+)\)", refs_repl, body)
    return body, n_sub


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    orig = parse_original(text)
    if len(orig) != len(KEYS_43):
        print(f"expected {len(KEYS_43)} references, parsed {len(orig)}", file=sys.stderr)
        return 1

    rows, n_fixed = [], 0
    for key, r in zip(KEYS_43, orig):
        formatted = r["formatted"]
        if key in TITLE_FIXES:
            old, new = TITLE_FIXES[key]
            if old not in formatted:
                print(f"title fix for {key} did not match:\n  {formatted}", file=sys.stderr)
                return 1
            formatted = formatted.replace(old, new)
            n_fixed += 1
        rows.append({"key": key, "formatted": formatted, "doi": r["doi"]})

    for key, authors, title, journal, year, vol, page, doi in NEW_REFS:
        # Nature style, matching the strings already in the file: authors, title, italic
        # journal, volume, pages, year.
        # "Hunter JD" needs the sentence period the original list has; "Ioannidis NM, *et al.*"
        # already ends in one, inside the italic markers.
        stop = "" if authors.rstrip("*").endswith(".") else "."
        rows.append({"key": key,
                     "formatted": f"{authors}{stop} {title}. *{journal}* {vol}, {page} ({year})",
                     "doi": doi})

    refs = pd.DataFrame(rows)
    if refs.key.duplicated().any():
        print("duplicate keys: " + ", ".join(refs.key[refs.key.duplicated()]), file=sys.stderr)
        return 1
    blank = refs[(refs.formatted.str.len() < 30) | (refs.doi == "")]
    if len(blank):
        print("incomplete entries:\n" + blank.to_string(), file=sys.stderr)
        return 1
    dup_doi = refs[refs.doi.duplicated(keep=False)]
    if len(dup_doi):
        print("duplicate DOIs — the same paper cited twice:\n" + dup_doi.to_string(),
              file=sys.stderr)
        return 1
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    refs.to_csv(OUT_TSV, sep="\t", index=False)

    # body = everything before the reference list; the list itself is now data
    body = text[: text.index("## References")].rstrip() + "\n"
    num2key = {i + 1: k for i, k in enumerate(KEYS_43)}
    body, n_sub = superscript_to_keys(body, num2key)

    leftover = SUP_RUN.findall(body)
    if leftover:
        print(f"WARNING: {len(leftover)} superscript runs not converted: {leftover[:8]}",
              file=sys.stderr)

    # This script is a ONE-SHOT migration. Once the keyed file exists it is the hand-maintained
    # source of the paper -- citations get inserted into it -- so regenerating would silently
    # discard that editing. Refuse, and write beside it instead.
    if OUT_MD.exists() and "--force" not in sys.argv:
        alt = OUT_MD.with_suffix(".regenerated.md")
        writing_to = alt
        print(f"NOTE: {OUT_MD.relative_to(ROOT)} already exists and is the maintained source.\n"
              f"      Writing the regenerated conversion to {alt.relative_to(ROOT)} instead.\n"
              f"      Pass --force to overwrite the maintained file.")
    else:
        writing_to = OUT_MD

    header = (
        "<!--\n"
        "  Keyed manuscript source. Citations are bracketed keys; the reference list is\n"
        "  data/processed/references_final.tsv and numbering is assigned at build time by\n"
        "  scripts/59_typeset_docx.py, in order of first appearance. Do not hand-number.\n"
        "  Regenerate this file's citation scheme with scripts/58_manuscript_source.py.\n"
        "-->\n\n"
    )
    writing_to.write_text(header + body, encoding="utf-8")

    print(f"references: {len(refs)} ({len(KEYS_43)} original + {len(NEW_REFS)} added)")
    print(f"citations converted: {n_sub} runs")
    print(f"title corrections applied: {n_fixed} ({', '.join(TITLE_FIXES)})")
    print(f"wrote {writing_to.relative_to(ROOT)} and {OUT_TSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
