#!/usr/bin/env python3
"""
57_references.py -- verify every manuscript reference against Crossref and Europe PMC.

The reference list in `results/manuscript_final.md` was assembled by hand over several sessions.
Hand-assembled reference lists carry three kinds of error that a reader can check and an author
usually does not: a DOI that resolves to a different paper, a volume/page/year that disagrees with
the record, and an author list that has drifted. This script checks all three against two
independent sources and writes the result as data.

It never edits the manuscript. It reports, and a human decides -- a wrong DOI might be a typo in
our file or a genuinely wrong citation, and those need different fixes.

Outputs
  data/processed/references.tsv        one row per reference, with the authoritative metadata
  results/reference_verification.md    the report, including every disagreement

Usage
  .venv/bin/python scripts/57_references.py
"""
from __future__ import annotations

import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANUSCRIPT = ROOT / "results/manuscript_final.md"
OUT_TSV = ROOT / "data/processed/references.tsv"
OUT_MD = ROOT / "results/reference_verification.md"
LOGS = ROOT / "logs"

MAILTO = "crabonc@proton.me"          # Crossref asks for a contact; it buys the polite pool
UA = f"fbn1-marfan-refcheck/1.0 (mailto:{MAILTO})"
CROSSREF = "https://api.crossref.org/works/"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


# ---------------------------------------------------------------------------------------
# Parsing the manuscript's own list
# ---------------------------------------------------------------------------------------

REF_LINE = re.compile(r"^(\d+)\.\s+(.*)$")
DOI_RE = re.compile(r"doi:\s*(10\.\d{4,9}/\S+?)\s*$", re.I)


def parse_manuscript_refs(path: Path) -> list[dict]:
    """Pull the numbered list out of the `## References` section."""
    text = path.read_text(encoding="utf-8")
    start = text.index("## References")
    # the list ends at the next h2, or at the source note
    tail = text[start + len("## References"):]
    end = tail.find("\n**Note on sources.")
    if end == -1:
        end = len(tail)
    block = tail[:end]

    refs, cur = [], None
    for raw in block.splitlines():
        line = raw.rstrip()
        if not line.strip():
            continue
        m = REF_LINE.match(line.strip())
        if m and (cur is None or int(m.group(1)) == cur["n"] + 1):
            if cur:
                refs.append(cur)
            cur = {"n": int(m.group(1)), "text": m.group(2).strip()}
        elif cur is not None:
            cur["text"] += " " + line.strip()
    if cur:
        refs.append(cur)

    for r in refs:
        m = DOI_RE.search(r["text"])
        r["doi"] = m.group(1).rstrip(".") if m else ""
        r["cited_text"] = DOI_RE.sub("", r["text"]).strip()
    return refs


# ---------------------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------------------

def _get(url: str, tries: int = 4) -> dict | None:
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA,
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as fh:
                return json.loads(fh.read().decode("utf-8"))
        except Exception as exc:                      # noqa: BLE001 - report, never guess
            if attempt == tries - 1:
                print(f"    ! {type(exc).__name__}: {exc}", file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


def crossref(doi: str) -> dict | None:
    if not doi:
        return None
    url = CROSSREF + urllib.parse.quote(doi, safe="") + f"?mailto={MAILTO}"
    js = _get(url)
    if not js or js.get("status") != "ok":
        return None
    m = js["message"]

    def first_page(p):
        return (p or "").split("-")[0]

    issued = m.get("issued", {}).get("date-parts", [[None]])[0]
    authors = []
    for a in m.get("author", []) or []:
        fam, giv = a.get("family", ""), a.get("given", "")
        if fam:
            initials = "".join(x[0] for x in re.split(r"[\s.-]+", giv) if x)
            authors.append(f"{fam} {initials}".strip())
    return {
        "source": "crossref",
        "title": (m.get("title") or [""])[0],
        "journal": (m.get("short-container-title") or m.get("container-title") or [""])[0],
        "journal_full": (m.get("container-title") or [""])[0],
        "year": issued[0] if issued else None,
        "volume": m.get("volume", ""),
        "page": first_page(m.get("page", "")),
        "article_number": m.get("article-number", ""),
        "authors": authors,
        "type": m.get("type", ""),
    }


def epmc(doi: str) -> dict | None:
    if not doi:
        return None
    q = urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json",
                                "resultType": "core", "pageSize": 1})
    js = _get(f"{EPMC}?{q}")
    if not js:
        return None
    hits = js.get("resultList", {}).get("result", [])
    if not hits:
        return None
    r = hits[0]
    return {
        "source": "epmc",
        "title": (r.get("title") or "").rstrip("."),
        "journal": (r.get("journalInfo", {}).get("journal", {}) or {}).get("medlineAbbreviation")
                   or (r.get("journalInfo", {}).get("journal", {}) or {}).get("title", ""),
        "year": int(r["pubYear"]) if r.get("pubYear", "").isdigit() else None,
        "volume": (r.get("journalInfo", {}) or {}).get("volume", ""),
        "page": (r.get("pageInfo") or "").split("-")[0],
        "authors": [a.strip() for a in (r.get("authorString") or "").rstrip(".").split(",")],
        "pmid": r.get("pmid", ""),
        "pmcid": r.get("pmcid", ""),
    }


# ---------------------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------------------

def norm(s: str) -> str:
    """Fold to something comparable: no accents, no punctuation, no case, no runs of space."""
    s = unicodedata.normalize("NFKD", s or "")
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.replace("–", "-").replace("—", "-").replace("’", "'")
    s = re.sub(r"[^a-z0-9 ]+", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def title_matches(cited: str, authoritative: str) -> bool:
    """
    The cited text is a whole reference; the authoritative title is a substring of it, modulo
    punctuation. Fall back to a token-overlap test for the cases where a subtitle was trimmed.
    """
    c, a = norm(cited), norm(authoritative)
    if not a:
        return False
    if a in c:
        return True
    at, ct = set(a.split()), set(c.split())
    if len(at) < 4:
        return False
    return len(at & ct) / len(at) >= 0.80


def check(ref: dict, cr: dict | None, ep: dict | None) -> list[str]:
    """Return a list of disagreements. An empty list means the citation checks out."""
    problems = []
    best = cr or ep
    if best is None:
        problems.append("no record found at Crossref or Europe PMC")
        return problems

    cited = ref["cited_text"]
    if not title_matches(cited, best["title"]):
        problems.append(f'title: DOI resolves to "{best["title"]}"')

    if best.get("year"):
        if not re.search(rf"\({best['year']}\)", cited):
            years = re.findall(r"\((\d{4})\)", cited)
            problems.append(f"year: record says {best['year']}, citation says "
                            f"{years[0] if years else 'none'}")

    vol = str(best.get("volume") or "")
    if vol and not re.search(rf"(?<!\d){re.escape(vol)}(?!\d)", cited):
        problems.append(f"volume: record says {vol}")

    page = str(best.get("page") or "") or str(best.get("article_number") or "")
    if page and not re.search(rf"(?<!\d){re.escape(page)}(?!\d)", cited):
        problems.append(f"page/article: record says {page}")

    # first author surname must appear
    if best.get("authors"):
        surname = best["authors"][0].split()[0]
        if norm(surname) and norm(surname) not in norm(cited):
            problems.append(f'first author: record says {best["authors"][0]}')

    # where both sources answered, make them agree with each other too
    if cr and ep:
        if cr.get("year") and ep.get("year") and cr["year"] != ep["year"]:
            problems.append(f"Crossref year {cr['year']} != Europe PMC year {ep['year']}")
    return problems


# ---------------------------------------------------------------------------------------

def main() -> int:
    if not MANUSCRIPT.exists():
        print(f"missing {MANUSCRIPT}", file=sys.stderr)
        return 1

    refs = parse_manuscript_refs(MANUSCRIPT)
    print(f"parsed {len(refs)} references from {MANUSCRIPT.relative_to(ROOT)}")
    if not refs:
        return 1

    rows, log = [], [f"reference verification {datetime.now(timezone.utc).isoformat()}"]
    n_clean = n_dirty = n_missing = 0

    for ref in refs:
        print(f"  [{ref['n']:>2}] {ref['doi'] or '(no doi)'}", flush=True)
        cr = crossref(ref["doi"])
        ep = epmc(ref["doi"])
        time.sleep(0.15)                       # both services ask for restraint; this is plenty
        problems = check(ref, cr, ep)
        best = cr or ep or {}

        if not ref["doi"]:
            problems.insert(0, "no DOI in the citation")
        if "no record found" in " ".join(problems):
            n_missing += 1
        elif problems:
            n_dirty += 1
        else:
            n_clean += 1

        rows.append({
            "n": ref["n"],
            "doi": ref["doi"],
            "cited_text": ref["cited_text"],
            "found_crossref": bool(cr),
            "found_epmc": bool(ep),
            "title": best.get("title", ""),
            "journal": best.get("journal", ""),
            "year": best.get("year") or "",
            "volume": best.get("volume", ""),
            "page": best.get("page", "") or best.get("article_number", ""),
            "authors": "; ".join(best.get("authors", [])),
            "pmid": (ep or {}).get("pmid", ""),
            "problems": " | ".join(problems),
        })
        for p in problems:
            log.append(f"[{ref['n']}] {ref['doi']}: {p}")

    import pandas as pd
    df = pd.DataFrame(rows)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TSV, sep="\t", index=False)

    lines = [
        "# Reference verification",
        "",
        f"Generated by `scripts/57_references.py`, {datetime.now(timezone.utc).isoformat()}.",
        "Every DOI in `results/manuscript_final.md` was resolved against the Crossref REST API and",
        "Europe PMC independently, and the returned title, year, volume, first page and first",
        "author compared with the citation as written.",
        "",
        f"- references checked: **{len(refs)}**",
        f"- clean (all fields agree): **{n_clean}**",
        f"- disagreements: **{n_dirty}**",
        f"- no record found: **{n_missing}**",
        f"- resolved at Crossref: {int(df.found_crossref.sum())} / Europe PMC: "
        f"{int(df.found_epmc.sum())}",
        "",
    ]
    bad = df[df.problems.astype(bool)]
    if len(bad):
        lines += ["## Disagreements", "", "| # | DOI | problem |", "|---|---|---|"]
        for _, r in bad.iterrows():
            lines.append(f"| {r.n} | `{r.doi}` | {r.problems.replace('|', '/')} |")
        lines.append("")
    else:
        lines += ["Every reference checks out against both services.", ""]

    lines += ["## Verified list", "",
              "| # | authors | title | journal | year | vol | page |",
              "|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        a = r.authors.split("; ")
        who = a[0] + (" et al." if len(a) > 1 else "") if r.authors else "—"
        lines.append(f"| {r.n} | {who} | {r.title} | {r.journal} | {r.year} | {r.volume} "
                     f"| {r.page} |")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    LOGS.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (LOGS / f"57_references_{stamp}.log").write_text("\n".join(log) + "\n", encoding="utf-8")

    print(f"\n{n_clean} clean, {n_dirty} with disagreements, {n_missing} not found")
    print(f"wrote {OUT_TSV.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
