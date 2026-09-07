# manifest/alphamissense.md — AlphaMissense P35555 subset

- Source release: `AlphaMissense_aa_substitutions.tsv.gz`, Zenodo record **8208688**
  (v1.0.0, DOI 10.5281/zenodo.8208688), placed by the mentor 2026-08-01
- Full archive sha256: `bc18d954a9aa6b9efab8e40d6d1912e32b774906cd3846a1f51066588e9c176d`
  (1,207,278,510 bytes)
- Subset: `resources/databases/AlphaMissense_FBN1_P35555.tsv` — **54,568 data rows**, 1,649,641 bytes
- Subset sha256: `925b3d1a2f96fdf9e8ed1ecb40b243d1946b88c111ec62d522e3e897e5f00b69`
- Derivation: copyright comment line + rows matching `^P35555` from the full archive
- Licence: CC BY-NC-SA 4.0, academic use. Cite Cheng et al. (2023) Science 381:eadg7492,
  DOI 10.1126/science.adg7492 (PDF present in `resources/papers/`)
- Verified by `tests/test_variant_ingest.py` (Phase 1 gate) and `scripts/00_env_check.py`

## Parsing notes

The subset retained only the copyright comment line, **not** the real column header. Columns must
be named explicitly when reading:
`uniprot_id`, `protein_variant`, `am_pathogenicity`, `am_class`.

## Known anomaly carried into Phase 2

Position **472** carries two wild-type residues (`C` and `Y`), giving 38 rows instead of 19 and
explaining the total of 54,568 vs the expected 2871 x 19 = 54,549. **UniProt P35555 residue 472 is
`Y`**, so the `C472*` rows correspond to a different reference allele and must be dropped by the
wild-type identity check. Joining on position alone would double-count this position.

## Class boundaries

Use the shipped `am_class` column (`benign` / `ambiguous` / `pathogenic`). The Science paper in
`resources/papers/` is the Research Article Summary and does **not** state the numeric score
cutoffs, so they are not reproduced here rather than cited from memory.
