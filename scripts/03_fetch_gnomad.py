#!/usr/bin/env python3
"""
03_fetch_gnomad.py — PHASE 1. Pull gnomAD v4 variant frequencies for FBN1.

Schema was verified live against the API before writing this query (the `Variant` type exposes
`faf95_joint { popmax popmax_population }`, and `joint` has no `af` field — allele frequency is
computed as ac/an). Nothing here is filtered on frequency: the benign-set threshold is a Phase 2
decision that needs mentor sign-off, so this stage only retrieves and records.

`faf95_joint.popmax` is the **grpmax filtering allele frequency** that CLAUDE.md specifies for the
benign threshold — captured per variant along with the genetic ancestry group it came from.

Idempotent: skips the pull if a cached JSON for the same dataset exists.

Writes:
  data/raw/gnomad_<dataset>_FBN1_<date>.json   verbatim API response
  data/raw/gnomad_<dataset>_FBN1_<date>.tsv    flattened one row per variant
  manifest/gnomad.md                           provenance
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from datetime import date, datetime, timezone
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"
GENE, REF = "FBN1", "GRCh38"
TODAY = date.today().isoformat()
UA = {"Content-Type": "application/json", "User-Agent": "fbn1-marfan-pipeline/0.1 (Phase 1)"}

QUERY = """
query($sym:String!, $ref:ReferenceGenomeId!, $ds:DatasetId!){
  gene(gene_symbol:$sym, reference_genome:$ref){
    gene_id symbol chrom start stop canonical_transcript_id
    mane_select_transcript{ ensembl_id refseq_id }
    variants(dataset:$ds){
      variant_id chrom pos ref alt rsids flags
      consequence hgvsc hgvsp transcript_id gene_symbol
      exome{ ac an homozygote_count filters }
      genome{ ac an homozygote_count filters }
      joint{
        ac an homozygote_count filters
        fafmax{ faf95_max faf95_max_gen_anc faf99_max faf99_max_gen_anc }
      }
    }
  }
}"""
# NOTE: the top-level `faf95_joint { popmax }` field exists in the schema but is null for every
# FBN1 variant, including ones that certainly have a defined FAF (p.Pro1148Ala, joint AF 1.2%,
# faf95 0.254 in eas). The populated location in gnomAD v4 is `joint.fafmax.faf95_max`, verified
# live before this query was written. Using the empty field would have silently produced a
# benign set with no frequency evidence at all.


def post(query: str, variables: dict, tries: int = 4) -> dict:
    """gnomAD rate-limits; back off rather than hammering."""
    delay = 5.0
    for attempt in range(1, tries + 1):
        r = requests.post(API, json={"query": query, "variables": variables},
                          headers=UA, timeout=600)
        if r.status_code == 429:
            print(f"  rate-limited, sleeping {delay:.0f}s (attempt {attempt}/{tries})")
            time.sleep(delay)
            delay *= 2
            continue
        r.raise_for_status()
        js = r.json()
        if js.get("errors"):
            raise SystemExit("STOP: gnomAD GraphQL errors:\n  "
                             + json.dumps(js["errors"], indent=2)[:1500])
        return js
    raise SystemExit("STOP: gnomAD rate limit not cleared after retries")


def af(ac, an):
    return (ac / an) if (ac is not None and an) else None


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    print(f"gnomAD {DATASET} pull for {GENE} ({REF})\n")

    raw_path = RAW / f"gnomad_{DATASET}_{GENE}_{TODAY}.json"
    if raw_path.is_file():
        print(f"  using cached {raw_path.relative_to(ROOT)}")
        js = json.loads(raw_path.read_text(encoding="utf-8"))
    else:
        js = post(QUERY, {"sym": GENE, "ref": REF, "ds": DATASET})
        raw_path.write_text(json.dumps(js, indent=2), encoding="utf-8")
        print(f"  wrote {raw_path.relative_to(ROOT)} ({raw_path.stat().st_size:,} bytes)")

    gene = js["data"]["gene"]
    variants = gene["variants"]
    mane = gene.get("mane_select_transcript") or {}
    print(f"  gene {gene['gene_id']} {gene['symbol']} chr{gene['chrom']}:"
          f"{gene['start']}-{gene['stop']}")
    print(f"  MANE Select: {mane.get('refseq_id')} / {mane.get('ensembl_id')}")
    print(f"  variants returned: {len(variants):,}")

    # MANE cross-check: CLAUDE.md pins NM_000138.5. Report rather than assert equality,
    # because gnomAD may carry a different minor version than the transcript we pin.
    refseq = (mane.get("refseq_id") or "")
    if not refseq.startswith("NM_000138"):
        print(f"  WARN: gnomAD MANE Select is {refseq!r}, expected NM_000138.x — recorded, "
              "Phase 2 must reconcile")

    rows = []
    for v in variants:
        ex, ge, jt = v.get("exome") or {}, v.get("genome") or {}, v.get("joint") or {}
        fm = jt.get("fafmax") or {}
        rows.append({
            "variant_id": v["variant_id"], "chrom": v["chrom"], "pos": v["pos"],
            "ref": v["ref"], "alt": v["alt"],
            "rsids": ";".join(v.get("rsids") or []),
            "consequence": v.get("consequence"),
            "hgvsc": v.get("hgvsc"), "hgvsp": v.get("hgvsp"),
            "transcript_id": v.get("transcript_id"),
            "flags": ";".join(v.get("flags") or []),
            "exome_ac": ex.get("ac"), "exome_an": ex.get("an"),
            "exome_af": af(ex.get("ac"), ex.get("an")),
            "exome_nhomalt": ex.get("homozygote_count"),
            "exome_filters": ";".join(ex.get("filters") or []),
            "genome_ac": ge.get("ac"), "genome_an": ge.get("an"),
            "genome_af": af(ge.get("ac"), ge.get("an")),
            "genome_nhomalt": ge.get("homozygote_count"),
            "genome_filters": ";".join(ge.get("filters") or []),
            "joint_ac": jt.get("ac"), "joint_an": jt.get("an"),
            "joint_af": af(jt.get("ac"), jt.get("an")),
            "joint_nhomalt": jt.get("homozygote_count"),
            # grpmax filtering AF — the value CLAUDE.md specifies for the benign threshold
            "faf95_grpmax": fm.get("faf95_max"),
            "faf95_grpmax_group": fm.get("faf95_max_gen_anc"),
            "faf99_grpmax": fm.get("faf99_max"),
            "faf99_grpmax_group": fm.get("faf99_max_gen_anc"),
        })

    tsv = RAW / f"gnomad_{DATASET}_{GENE}_{TODAY}.tsv"
    with tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    import collections
    cons = collections.Counter(r["consequence"] for r in rows)
    n_faf = sum(1 for r in rows if r["faf95_grpmax"] not in (None, 0))
    n_missense = cons.get("missense_variant", 0)
    print(f"  missense: {n_missense:,}   with non-zero grpmax FAF95: {n_faf:,}")

    digest = hashlib.sha256(tsv.read_bytes()).hexdigest()
    (ROOT / "manifest/gnomad.md").write_text(
        f"# manifest/gnomad.md — gnomAD {DATASET} FBN1\n\n"
        f"Generated by `scripts/03_fetch_gnomad.py` at {datetime.now(timezone.utc).isoformat()}.\n\n"
        f"- Endpoint: `{API}` (GraphQL)\n"
        f"- Dataset: **{DATASET}**, reference genome **{REF}**, gene_symbol `{GENE}`\n"
        f"- Retrieved: {TODAY}\n"
        f"- gnomAD gene id: `{gene['gene_id']}`, locus chr{gene['chrom']}:"
        f"{gene['start']}-{gene['stop']}\n"
        f"- gnomAD MANE Select: `{mane.get('refseq_id')}` / `{mane.get('ensembl_id')}` "
        f"(project pins NM_000138.5)\n"
        f"- Raw response: `data/raw/{raw_path.name}`\n"
        f"- Flattened: `data/raw/{tsv.name}` — **{len(rows):,} variants**, sha256 `{digest}`\n"
        f"- No frequency filter applied at this stage; the benign AF cutoff is a Phase 2 "
        "decision requiring mentor approval.\n\n"
        "## Frequency fields captured\n\n"
        "`exome_*`, `genome_*` and `joint_*` allele counts (af computed as ac/an, since the\n"
        "joint type exposes no `af` field), plus **`faf95_grpmax`** / `faf99_grpmax` — the\n"
        "filtering allele frequency at the genetic-ancestry group maximum, with the group named\n"
        "in `faf95_grpmax_group`. This is the value CLAUDE.md specifies for the benign threshold.\n\n"
        f"- variants with a non-zero grpmax FAF95: **{n_faf:,}**\n\n"
        "## Consequence breakdown\n\n"
        + "\n".join(f"- {k or '(none)'}: {v:,}" for k, v in cons.most_common()) + "\n",
        encoding="utf-8")

    print(f"\n  wrote {tsv.relative_to(ROOT)} (sha256 {digest[:16]}…)")
    print("  wrote manifest/gnomad.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
