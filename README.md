# Calcium coordination and fold stability of pathogenic *FBN1* missense variants

Analysis code and derived data for *Distinct mechanisms of calcium-ligand and cysteine
substitutions in the fibrillin-1 cbEGF module*.

All 43 calcium-binding EGF-like (cbEGF) domains of fibrillin-1 were modeled with their calcium
ions as 21 overlapping tandem constructs, together with 60 mutant constructs, and the calcium
ligands were identified from observed contact with the modeled ion instead of from the sequence
consensus. Cysteine substitutions destabilize the domain in two independent energy functions.
Calcium-ligand substitutions do not, yet 45% of the pathogenic ones lose calcium coordination
while the backbone around the ion stays within half a bond length of wild type.

---

## What is here

| path | contents |
|---|---|
| `scripts/` | the pipeline, one numbered stage per file, `01_` through `65_` |
| `tests/` | the validation gates; every one must pass before a stage's output is used |
| `data/processed/` | analysis-ready tables, including the variant master table and its data dictionary |
| `data/interim/` | cleaned intermediates, including the HGVS failures that were excluded |
| `results/` | statistics, the manuscript and its supplementary workbook, per-phase reports |
| `figures/` | the six manuscript figures, each with the script that made it and a sidecar listing every value it prints |
| `manifest/` | provenance for every download, with date, query, row count and checksum |
| `handoff/` | the record of the two manual steps, CheckMyMetal and AlphaFold Server |
| `resources/reference/` | domain definitions, transcript identifiers, inclusion criteria, literature notes |
| `structures/alphafold/` | the AlphaFold Server job submissions and per-model confidence summaries |



```bash
# after downloading Zenodo record 8208688 and cutting the P35555 rows, see the script header
.venv/bin/python scripts/67_restore_alphamissense.py
```

The restored values are identical to those the paper was computed from. Until then, the checks in
`tests/test_v2.py` that read an AlphaMissense value are the only ones that cannot run. See
NOTICE.

**Journal PDFs.** `resources/papers/` holds publisher-copyright articles. The citations are in
`resources/reference/lit_notes.md` and `data/processed/references_final.tsv`.

**FoldX.** Licensed per user from [foldxsuite.crg.eu](https://foldxsuite.crg.eu). The binary and
its `molecules/` directory are not redistributable.

---

## Reproducing the analysis

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/00_env_check.py     # tool versions and one test query per database
.venv/bin/python tests/test_v2.py            # 15 checks on the analysis tables and the figures
.venv/bin/python tests/test_manuscript_v4.py # 12 checks on the manuscript
```

`tests/test_v2.py` re-runs every figure script and fails if any printed value moves, so it takes
about a minute. Both gates run on what is in this repository, with no download and no structure
prediction.

Three things need more than `requirements.txt`. **PyRosetta** is in a separate Python 3.11
environment. **FoldX 5** and **mkdssp** are external binaries. **PyMOL** is called headlessly for
one supplementary figure. Versions and paths are recorded in `manifest/tools.md`, and every stage
that needs one says so in its module docstring.

## Regenerating from raw data

Stages `01_` to `05_` re-download from ClinVar, UniProt, gnomAD and the PDB. They need an NCBI API
key and email in a `.env` file, which is git-ignored and must never be committed. Stages `20_`
onward run off the cached tables and need no network.

The two manual steps cannot be automated. `handoff/structures/REQUESTS.md` records what was
submitted to AlphaFold Server and `handoff/cmm/INSTRUCTIONS.md` what was submitted to
CheckMyMetal, in both cases with the settings used and the outputs returned.

---

## Where each number in the paper comes from

Every figure is produced by exactly one script and writes a `.numbers.json` sidecar listing every
value it prints. `tests/test_v2.py` re-derives those values from `data/processed/` and fails if
one has drifted. Supplementary Table S6 in `results/supplementary_tables_v3.xlsx` carries all 29
statistical tests with effect sizes and Benjamini-Hochberg corrected *q* values, and
`results/statistics_final.tsv` is the same table in plain text.

## Licence

Code is MIT. Data and figures are CC BY 4.0, with no exceptions — see NOTICE for why the
AlphaMissense columns are absent and how to restore them.

## Citation

Please cite the paper. If you use the models, cite the Zenodo archive as well.
