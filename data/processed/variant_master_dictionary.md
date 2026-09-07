# Data dictionary — `data/processed/variant_master.tsv`

Produced by `scripts/05_normalize_variants.py` (Phase 2). One row per **ClinVar GRCh38 missense
variant whose wild-type residue was verified against UniProt P35555**. 4,451 rows, 34 columns.

Criteria: `resources/reference/inclusion_criteria.md` (mentor-approved 2026-08-02).
Gate: `tests/test_normalization.py`.

## Provenance / identity

| Column | Type | Meaning |
|---|---|---|
| `variation_id` | int | ClinVar VariationID. Unique within this table. |
| `allele_id` | int | ClinVar AlleleID. |
| `name` | str | ClinVar `Name`, verbatim — the source of `hgvs_c` / `hgvs_p`. |
| `variant_type` | str | ClinVar `Type` (all rows here are `single nucleotide variant`). |
| `transcript` | str | Transcript in the ClinVar name. All rows: `NM_000138.5`. |
| `hgvs_c` | str | cDNA description on that transcript. |
| `hgvs_p` | str | Protein description, 3-letter (e.g. `p.Cys1039Tyr`). |

## Position and residues (P35555)

| Column | Type | Meaning |
|---|---|---|
| `position` | int | Protein position, **1–2871, P35555 numbering**. |
| `wt_aa` | char | Wild-type residue, 1-letter, parsed from `hgvs_p`. |
| `mut_aa` | char | Substituted residue, 1-letter. Always ≠ `wt_aa`. |
| `reference_aa` | char | Residue actually at `position` in P35555. **Guaranteed == `wt_aa`** — rows failing this are quarantined in `data/interim/wt_mismatch.tsv`, never included here. |

## Genomic coordinates (GRCh38)

`chrom`, `pos_vcf`, `ref_vcf`, `alt_vcf` — ClinVar's VCF-style coordinates, used to join gnomAD.

## ClinVar classification

| Column | Type | Meaning |
|---|---|---|
| `clinical_significance` | str | ClinVar's raw string, unmodified. |
| `clinvar_class` | enum | Bucketed: `P/LP`, `B/LB`, `VUS`, `conflicting`, `other`. |
| `review_status` | str | ClinVar review status, verbatim. |
| `stars` | int 0–4 | Gold stars derived from `review_status`. Every value is mapped; an unmapped status is a hard stop, never a silent 0. |
| `last_evaluated` | date | ClinVar last-evaluated date. |
| `n_submitters` | int | Number of submitters. |
| `phenotypes` | str | ClinVar `PhenotypeList`. |

## gnomAD v4 (`gnomad_r4`, GRCh38)

| Column | Type | Meaning |
|---|---|---|
| `gnomad_variant_id` | str | `chrom-pos-ref-alt`; empty if the variant is absent from gnomAD (2,928 of 4,451 — expected for rare pathogenic alleles). |
| `gnomad_joint_af` | float | Joint allele frequency, computed as `joint.ac / joint.an` (the API exposes no `af` on the joint type). |
| `faf95_grpmax` | float | **Filtering allele frequency**, 95% CI lower bound at the maximum genetic-ancestry group. Read from `joint.fafmax.faf95_max`. Empty where gnomAD reports none. |
| `faf95_grpmax_group` | str | Ancestry group that maximum came from (e.g. `eas`). |
| `gnomad_benign_support` | enum | `stand_alone` (faf95 > 0.001), `supporting` (> 0.0001), or empty. |
| `gnomad_hgvsp_disagrees` | str | gnomAD's protein description when it differs from ours. **Empty for all 4,451 rows** — the two sources agree completely. |

## AlphaMissense

| Column | Type | Meaning |
|---|---|---|
| `am_pathogenicity` | float | AlphaMissense score. Joined on **(position, wild-type, alternate)**, never position alone. |
| `am_class` | enum | `benign` / `ambiguous` / `pathogenic`, as shipped. Numeric cutoffs are not reproduced because the PDF in `resources/papers/` does not state them. |

Coverage is 4,451/4,451 (100%). 19 AlphaMissense rows were rejected because their wild-type
residue disagrees with P35555 — all at **position 472**, where AlphaMissense carries both `C` and
`Y` and the reference is `Y`.

## Labelled sets

Two tiers are provided; **both are computed for every row** so the sensitivity analysis needs no
re-derivation.

| Column | Meaning |
|---|---|
| `set_primary` | Label at the **≥2 star** cutoff: `pathogenic`, `benign`, `vus`, `conflicting`, `contradiction`, or empty (below the star floor). |
| `set_primary_source` | `clinvar`, `gnomad_faf`, or `both`. |
| `set_sensitivity` | Same at the **≥1 star** cutoff. |
| `set_sensitivity_source` | As above. |

`contradiction` marks a ClinVar P/LP call whose grpmax filtering AF exceeds the stand-alone
benign threshold. **There are currently 0 of these.** The gate fails if such a variant is ever
labelled `pathogenic` without the flag.

### Counts

| Set | Primary (≥2★) | Sensitivity (≥1★) |
|---|---|---|
| pathogenic | 578 | 1,409 |
| benign | 34 | 76 |
| vus | 823 | 2,295 |
| conflicting | 0 (1★ by definition) | 492 |
| unlabelled | 3,016 | 179 |

`pathogenic ∩ benign = ∅` in both tiers, asserted by the gate.

## Companion quarantine tables (nothing is deleted)

| File | Rows | Contents |
|---|---|---|
| `data/interim/excluded_from_structural.tsv` | 4,825 | Real variants of a non-missense consequence, each with a `reason`. |
| `data/interim/wt_mismatch.tsv` | **0** | Would hold variants whose WT disagrees with P35555. |
| `data/interim/unparsed.tsv` | 51 | Genomic-only descriptions with no transcript-level HGVS. |

4,451 + 4,825 + 0 + 51 = **9,327**, exactly the ClinVar GRCh38 input count.
