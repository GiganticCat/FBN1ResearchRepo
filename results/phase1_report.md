# Phase 1 report — Variant collection

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED (10/10) · **Awaiting approval for Phase 2.**

All four sources are cached in `data/raw/` with provenance manifests. Nothing was filtered on
clinical significance or frequency — labelling and thresholds are Phase 2 decisions that need
your sign-off.

---

## 1. What was collected

| Source | Artifact | Rows | Version / release |
|---|---|---|---|
| ClinVar | `data/raw/clinvar_fbn1_2026-07-28.tsv` | **18,770** (9,327 GRCh38) | bulk `variant_summary`, release **2026-07-28** |
| UniProt | `data/raw/uniprot_P35555_v264.{json,gff,fasta}` | 2871 aa | entry version **264** |
| gnomAD | `data/raw/gnomad_gnomad_r4_FBN1_2026-08-02.tsv` | **12,107** | dataset **gnomad_r4**, GRCh38 |
| AlphaMissense | `resources/databases/AlphaMissense_FBN1_P35555.tsv` | **54,568** | Zenodo 8208688 v1.0.0 |

Scripts: `01_fetch_clinvar.py`, `02_fetch_uniprot.py`, `03_fetch_gnomad.py`.
Manifests: `manifest/{clinvar,uniprot,gnomad,alphamissense}.md`. All three fetchers are
idempotent — a re-run reuses the cached archive rather than re-downloading.

### ClinVar composition (GRCh38, 9,327 records)

| Review status | Stars | n |
|---|---|---|
| reviewed by expert panel | 3★ | **120** |
| criteria provided, multiple submitters, no conflicts | 2★ | 2,648 |
| criteria provided, single submitter | 1★ | 5,437 |
| criteria provided, conflicting classifications | 1★ | 704 |
| no assertion criteria provided | 0★ | 396 |
| no classification / `-` | 0★ | 22 |

Classifications: 2,740 VUS · 2,009 Pathogenic · 1,969 Likely benign · 1,119 Likely pathogenic ·
704 Conflicting · 484 P/LP · 161 Benign · 119 B/LB. **VUS and conflicts were retained and
labelled, never dropped** — verified by the gate.

### UniProt domain model

47 EGF-like domains, **43 calcium-binding (cbEGF)**, 9 TB, 155 disulfide bonds. Every annotated
disulfide endpoint was verified to be an actual cysteine. Derived tables:
`data/interim/uniprot_domains.tsv` (with a cbEGF1..43 index in sequence order) and
`data/interim/uniprot_disulfides.tsv` (with the domain each partner sits in).

### gnomAD v4

12,107 variants, **2,783 missense**; 4,847 carry a grpmax filtering allele frequency. Captured
exome/genome/joint allele counts plus `faf95_grpmax` / `faf99_grpmax` and the ancestry group each
came from.

---

## 2. Two bugs caught — one of which would have corrupted the study

### 2.1 🔴 The Phase 0 ClinVar count was wrong (now corrected)

Phase 0 reported *"13,501 records for Gene ID 2200"*. **That number was an artifact of a broken
query.** `[gene_id]` is not a field in the ClinVar E-utilities index, and NCBI **silently
rewrote** `2200[gene_id]` into `2200[All Fields]` — a free-text search for the string "2200".
It was matching unrelated genes: COL2A1 `c.2276G>A`, ALK `c.2200T>C`, and so on.

I found this because the bulk pull returned 9,327 GRCh38 rows against E-utilities' 13,501, and
the ~4,000 "missing" records turned out on inspection to be **other genes entirely**.

- Correct query: **`FBN1[gene]` → 9,475 records.**
- Our bulk pull holds 9,446 distinct VariationIDs; the diff is **29 records present in the live
  index but not the 2026-07-28 snapshot**, and **0 records in our file that are absent upstream**
  — i.e. a clean subset, with the 29 being additions in the 5 days since the release.
- `scripts/00_env_check.py` now uses `FBN1[gene]` **and fails if the returned `querytranslation`
  degrades to `[All Fields]`**, so this class of silent mistranslation cannot recur.
- The bad number never propagated — it appeared only in the Phase 0 reachability table, which is
  corrected with a visible note.

### 2.2 🔴 The gnomAD filtering-AF field was silently empty

The obvious field, `faf95_joint { popmax }`, **exists in the schema and returns `null` for every
one of the 12,107 FBN1 variants** — including p.Pro1148Ala, which has a joint AF of 1.2% and
certainly has a defined FAF. Had I trusted it, the benign set would have been built with **no
frequency evidence whatsoever**, and it would have looked like a legitimate "no common variants"
result rather than a bug.

The populated location in gnomAD v4 is **`joint.fafmax.faf95_max`** (with
`faf95_max_gen_anc` naming the ancestry group). Verified live: p.Pro1148Ala →
**faf95_max = 0.254 in `eas`**. Now 4,847 variants carry a filtering AF.

The Phase 1 gate asserts this column is non-empty, so an all-null FAF can never pass silently.

---

## 3. Things Phase 2 must handle

1. **Two transcript versions.** 9,260 records are on `NM_000138.5` (the pinned MANE Select) and
   **16 on `NM_000138.4`**. These must be reconciled, not mixed silently.
2. **51 genomic-only descriptions** (0.5%) have no transcript-level HGVS, e.g.
   `NC_000015.9:g.48793164_49095744del`. Retained; they need genomic→protein mapping or explicit
   exclusion with a reason.
3. **AlphaMissense position 472 has two wild-type residues** (`C` and `Y`). UniProt says **`Y`**;
   the `C472*` rows must be dropped by the WT-identity check. Joining on position alone
   double-counts.
4. **gnomAD reports MANE Select as `NM_000138`** without a version suffix — compatible with our
   pin, but the version cannot be confirmed from the API and should not be assumed.
5. **No UniProt calcium-site annotations exist** (0 `Binding site` features). The Phase 3 gate
   as written in CLAUDE.md wants to validate computed coordinating residues against UniProt's
   explicit Ca-binding annotations; that cross-check is unavailable and must be replaced (see
   `resources/reference/structures.md`).

---

## 4. Decisions I need from you before Phase 2 finalises the labelled sets

CLAUDE.md requires your sign-off on these, and they change which variants enter the analysis:

1. **ClinVar review-status cutoff.** Options: ≥2★ only (2,768 records, high confidence, but drops
   the large 1★ single-submitter tier); ≥1★ (8,909); or ≥1★ with 3★ expert-panel calls given
   precedence on conflicts. **My recommendation: ≥2★ for the primary analysis, with a ≥1★
   sensitivity analysis** — FBN1 has a ClinGen VCEP whose 120 expert-panel calls are the highest
   quality evidence available, and Baudhuin 2019 showed 1★ single-submitter FBN1 calls are
   frequently overclassified.
2. **gnomAD benign AF threshold.** The ClinGen/ACMG convention for a dominant disorder is a
   filtering-AF cutoff (BA1/BS1). **My recommendation: `faf95_grpmax` > 0.001 (0.1%) as
   stand-alone benign (BA1-like), > 0.0001 (0.01%) as supporting** — but Marfan's population
   prevalence (~1 in 10,000, from Handford 1995) argues for a strict cutoff, so I'd like your
   call rather than mine.
3. **Conflicting classifications** (704 records): label and analyse separately, or exclude from
   the primary comparison? **My recommendation: retain, labelled, excluded from the primary
   P vs B contrast but reported.**

I'll write whatever you decide into `resources/reference/inclusion_criteria.md` before building
the labelled sets.

---

## 5. Gate

`tests/test_variant_ingest.py` — **10/10 pass**. Self-tested against simulated failures (empty
FAF column, unmapped review status); both were caught and the files restored. The gate checks
row-count envelopes, schema, per-assembly VariationID uniqueness, star mapping completeness, HGVS
parse rate, retention of VUS/conflicts, UniProt domain and disulfide integrity, gnomAD FAF
population, and manifest presence.

**No Phase 2 work has started.**
