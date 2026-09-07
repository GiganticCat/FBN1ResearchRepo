# Phase 2 report — Variant normalization & validation

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED (10/10) · **Awaiting approval for Phase 3.**

Thresholds you approved are recorded in `resources/reference/inclusion_criteria.md` and were
fixed **before** any set was built, so they could not be tuned to a result.

---

## 1. Headline: the reference mapping is clean

**0 wild-type mismatches out of 4,451 missense variants.** Every variant's wild-type residue,
parsed from ClinVar's protein description, matches the residue actually present at that position
in UniProt P35555. This is the check that catches transcript/isoform misalignment, and it came
back perfect — which is the outcome that lets Phase 3 trust its positions.

Three further independent confirmations:

- **Mutalyzer:** 100 randomly sampled variants (seed 20260802) re-validated against
  `NM_000138.5`. **100/100 agree**, 0 disagreements, 0 errors.
- **gnomAD:** of the 1,523 variants also present in gnomAD v4, **0 have a protein description
  that differs from ours**.
- **Row reconciliation:** 4,451 normalized + 4,825 excluded + 0 mismatched + 51 unparsed =
  **9,327**, exactly the ClinVar GRCh38 input. Nothing was dropped or duplicated.

All 9,327 ClinVar records on the pinned transcript turned out to be `NM_000138.5`; the 16
`NM_000138.4` records flagged in Phase 1 were on GRCh37 rows, so the GRCh38 working set is
transcript-homogeneous. No mixing occurred.

---

## 2. What the 9,327 records became

| Outcome | n | Notes |
|---|---|---|
| **Missense, WT-verified** | **4,451** | the structural analysis set |
| No protein consequence (splice/intronic) | 1,899 | retained, labelled |
| Synonymous | 1,295 | retained, labelled |
| Frameshift | 975 | retained, labelled |
| Nonsense (Ter) | 491 | retained, labelled |
| In-frame indel | 162 | retained, labelled |
| Unknown/other protein effect | 3 | retained, labelled |
| Genomic-only description (no transcript HGVS) | 51 | retained, labelled |
| **Wild-type mismatch** | **0** | — |

Nothing was deleted. Every excluded row carries an explicit reason.

## 3. Labelled sets

| Set | Primary (≥2★) | Sensitivity (≥1★) |
|---|---|---|
| pathogenic | **578** | 1,409 |
| benign | **34** | 76 |
| VUS | 823 | 2,295 |
| conflicting | 0 | 492 |
| unlabelled (below star floor) | 3,016 | 179 |

`pathogenic ∩ benign = ∅` in both tiers. **0 contradictions** — no ClinVar P/LP call exceeds the
gnomAD stand-alone benign threshold.

### Independent agreement with AlphaMissense

AlphaMissense played no part in building the sets, so this is a genuine external check:

| Our label (primary) | AM pathogenic | AM ambiguous | AM benign |
|---|---|---|---|
| pathogenic (578) | **562 (97.2%)** | 4 | 12 |
| benign (34) | 3 | 1 | **30 (88.2%)** |
| VUS (823) | 265 | 95 | 463 |

That concordance is strong but not suspiciously perfect, which is what I'd expect from correct
labelling. The 12 pathogenic-but-AM-benign and 3 benign-but-AM-pathogenic variants are worth a
look in Phase 5 as potential mechanism outliers rather than errors.

---

## 4. ⚠️ The one thing that limits the study: the benign set is small

**34 benign missense variants at ≥2★, 76 at ≥1★.** I checked whether gnomAD could enlarge it:

> Of the 2,783 gnomAD missense variants, only **24** exceed the stand-alone benign filtering-AF
> threshold (faf95_grpmax > 0.001) — and **all 24 are already in our ClinVar-derived set**.
> There are exactly **0** gnomAD-only variants available to add.

So this is **real biology, not a pipeline gap**: FBN1 is under strong selective constraint and
simply has very few common missense variants. It does mean a direct pathogenic-vs-benign
enrichment test at cbEGF calcium sites and cysteines will be **underpowered on the benign arm**,
especially once split by site class.

**My recommendation for Phase 5**, which I'll implement unless you object:

1. **Primary comparator: the all-possible-missense background.** AlphaMissense covers all 54,549
   possible substitutions on P35555, so "are pathogenic variants enriched at coordinating
   residues / cysteines relative to what a random missense change would hit?" is well-powered
   and is exactly the comparison the study goal names. This becomes the headline test.
2. **Secondary comparator: observed benign variants**, reported with exact confidence intervals
   and an explicit power statement rather than a bare p-value.
3. **Report both star tiers**, so the benign arm can borrow power from ≥1★ (76) while the
   pathogenic arm stays conservative.

This changes emphasis, not scope — CLAUDE.md already names all three comparators. Flagging it
because it affects how the headline result is framed.

---

## 5. Notes carried into Phase 3

- **No UniProt calcium-site annotations exist** (0 `Binding site` features on P35555). The Phase 3
  gate in CLAUDE.md wants computed coordinating residues validated against UniProt's explicit
  Ca-binding annotations; that cross-check is unavailable. Replacement: derive from the cbEGF
  consensus, then validate against (a) the 155 UniProt disulfide annotations, (b) the
  experimental Ca sites in `resources/reference/structures.md`, and (c) the consensus in
  Handford 1995 / Rao 1995 / Baudhuin 2019.
- **Now that the ClinGen FBN1 VCEP PDF is correct**, Phase 3's critical-residue definition will
  follow the VCEP rules, with Baudhuin 2019's taxonomy as corroboration.
- The position-472 AlphaMissense anomaly is resolved: 19 rows rejected on WT identity, and the
  gate fails if the join ever reverts to position-only matching.

---

## 6. Outputs

| Path | Contents |
|---|---|
| `data/processed/variant_master.tsv` | 4,451 variants × 34 columns |
| `data/processed/variant_master_dictionary.md` | full column documentation |
| `data/interim/variants_normalized.tsv` | same set, interim copy |
| `data/interim/excluded_from_structural.tsv` | 4,825 non-missense, with reasons |
| `data/interim/wt_mismatch.tsv` | 0 rows |
| `data/interim/unparsed.tsv` | 51 genomic-only |
| `resources/reference/inclusion_criteria.md` | approved thresholds |
| `logs/05_normalize_*.log` | run transcript |

**Gate:** `tests/test_normalization.py` — 10/10 pass, self-tested against four simulated
failures (corrupted WT residue, pathogenic/benign overlap, hidden contradiction, dropped rows);
all four were caught and the data restored.

**No Phase 3 work has started.**
