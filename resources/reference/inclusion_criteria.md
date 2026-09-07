# Inclusion criteria — FBN1 / Marfan variant sets

**Approved by the mentor on 2026-08-02** (they endorsed the Phase 1 report's recommendations,
explicitly retaining oversight). These thresholds are fixed BEFORE the labelled sets are built,
so they cannot be tuned to a result. Any later change must be recorded here with a date and a
reason.

## Reference pins

- Transcript: **NM_000138.5** (MANE Select). ClinVar records on `NM_000138.4` are re-validated,
  not assumed equivalent.
- Protein: **UniProt P35555**, 2871 aa, entry version **264**, sequence MD5
  `f8f84ac9dc30f3ff640ad4f691fee2b1`. Signal peptide 1–27.
- Every protein position in this project is a P35555 position.

## 1. ClinVar review status (gold stars)

| Tier | Rule | Use |
|---|---|---|
| **Primary** | **≥ 2 stars** | the headline analysis |
| **Sensitivity** | **≥ 1 star** | reported alongside; conclusions must survive it |
| Excluded from both | 0 stars | retained in the master table, flagged `star_excluded` |

Star mapping: practice guideline 4★ · reviewed by expert panel 3★ · criteria provided, multiple
submitters, no conflicts 2★ · criteria provided, single submitter 1★ · criteria provided,
conflicting classifications 1★ · no assertion criteria / no classification / `-` 0★.

*Rationale:* FBN1 has a ClinGen VCEP supplying 120 expert-panel calls, and Baudhuin 2019 showed
1★ single-submitter FBN1 missense calls are frequently overclassified (31.6% of non-critical-
residue calls judged overclassified once MAF was considered). A 2★ floor buys precision; the 1★
sensitivity tier guards against the resulting loss of power.

**Expert-panel precedence:** where a 3★/4★ classification exists for a variant, it takes
precedence over lower-starred submissions for that variant.

## 2. Pathogenic set (P/LP)

ClinVar `ClinicalSignificance` in {`Pathogenic`, `Likely pathogenic`, `Pathogenic/Likely
pathogenic`}, at the star tier in force.

## 3. Benign set (B/LB)

A variant qualifies as benign if **either**:

- **ClinVar evidence:** `ClinicalSignificance` in {`Benign`, `Likely benign`, `Benign/Likely
  benign`} at the star tier in force; **or**
- **Population evidence (gnomAD v4.1, `gnomad_r4`):** `faf95_grpmax` **> 0.001** (0.1%) —
  a stand-alone BA1-style benign call.

Supporting-only tier, recorded but **not** sufficient alone: `faf95_grpmax` **> 0.0001** (0.01%).

The source supporting each benign call is recorded per variant (`benign_source` =
`clinvar` / `gnomad_faf` / `both`).

*Rationale:* Marfan syndrome is autosomal dominant with population prevalence ≈ 1 in 10,000
(Handford 1995), so a filtering allele frequency above 0.1% is incompatible with a fully
penetrant dominant disease allele. Using `faf95_grpmax` (the 95% CI lower bound of the allele
frequency in the highest genetic-ancestry group) rather than a global AF is what makes this
robust to ancestry-specific polymorphisms — e.g. p.Pro1148Ala reaches faf95 0.254 in East Asians
while its global joint AF is only 0.012.

## 4. VUS

`Uncertain significance` at the star tier in force. Retained as a labelled set, scored with
AlphaMissense, and reported — but never merged into either the pathogenic or benign set.

## 5. Conflicting classifications

`Conflicting classifications of pathogenicity` (704 records). **Retained and labelled
`conflicting`; excluded from the primary pathogenic-vs-benign contrast; reported separately.**
Never silently dropped or coerced to a class.

## 6. Structural-analysis subset

Only **missense** variants enter the structural/enrichment analysis. Synonymous, nonsense,
frameshift, splice, in-frame indels and genomic-only descriptions are retained in a labelled
`excluded_from_structural` table with the reason, never deleted.

## 7. Hard disqualifiers (quarantine, do not coerce)

- **Wild-type identity mismatch** against P35555 → `wt_mismatch` table. A mismatch means wrong
  numbering, wrong isoform, or wrong HGVS; it is never fixed by overwriting the residue.
- Protein position outside 1–2871.
- AlphaMissense rows whose wild-type residue disagrees with P35555 at that position (notably
  **position 472**, where AlphaMissense carries both `C` and `Y`; P35555 says **`Y`**).

## 8. Disjointness

`pathogenic ∩ benign = ∅` is asserted by the Phase 2 gate. A variant appearing in both is a
contradiction to be surfaced, not resolved by precedence.
