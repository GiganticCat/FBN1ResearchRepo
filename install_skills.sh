#!/usr/bin/env bash
# install_skills.sh — create the 7 FBN1/Marfan project skills under .claude/skills/
# Run from your project root (e.g. ~/research/fbn1-marfan). Safe to re-run (overwrites the 7 SKILL.md files).
set -euo pipefail
echo "Installing skills into: $PWD/.claude/skills"

mkdir -p ".claude/skills/data-validation"
cat > ".claude/skills/data-validation/SKILL.md" << '__SKILL_MD_EOF__'
---
name: data-validation
description: >
  Write and run the validation-gate scripts that guard every phase of the pipeline — schema
  checks, deduplication, disjoint-set checks, NaN accounting, count reconciliation between
  stages, and reproducibility on re-run. Use this skill WHENEVER a phase produces data that the
  next phase depends on — "validate the output", "write the test/gate", "QC this table", "check
  the data before continuing", or any Phase 0-6 gate. Use it proactively at the END of every
  data-producing step, because the whole pipeline's trustworthiness rests on catching silent
  errors (bad numbering, dropped rows, leaking labels) before they propagate.
---

# Validation gates (the checks that make results trustworthy)

Goal: after each data-producing stage, a matching `tests/` script must PASS before the next
stage runs. A gate failure stops the pipeline, logs why, and surfaces it — never proceed on
broken data.

## Universal checks (apply to every stage)
- **Schema:** expected columns present, correct dtypes, no unexpected columns.
- **Non-empty & sane size:** row counts within expected order of magnitude; a suspicious zero or
  a 10x jump is a failure, not a pass.
- **Dedup:** no duplicate primary keys (e.g. VariationID; domain/site/variant).
- **NaN accounting:** count and report NaNs per column; assert none were dropped silently —
  dropped rows must land in a labeled quarantine table, not vanish.
- **Count reconciliation:** the row counts entering a stage reconcile with those leaving the
  previous stage (in + quarantined + newly-excluded == out).
- **Reproducibility:** re-running the stage on the cached inputs yields identical outputs
  (checksum-stable); seeds are fixed.

## Phase-specific gates
- **Phase 0 (`test_env`):** imports succeed; `pymol`, `mkdssp`, `foldx` callable
  (`shutil.which`); one test query per database reachable; `lit_notes.md` has one entry per PDF.
- **Phase 1 (`test_variant_ingest`):** non-empty pulls; HGVS strings parse; no duplicate
  VariationIDs; review-status/stars captured; manifests written.
- **Phase 2 (`test_normalization`):** every kept variant maps to a real residue whose wild-type
  identity matches the P35555 reference; positions in 1-2871; no contradictory labels;
  benign ∩ pathogenic == ∅.
- **Phase 3 (`test_annotation`):** annotated cysteine positions are actually Cys; coordinating-
  residue set matches UniProt annotation within tolerance; class × domain-feature cross-tab is
  sane and non-empty.
- **Phase 4 (`test_structures`):** SIFTS numbering round-trips (position → author num → residue
  identity matches); every CMM-packet PDB contains a Ca ion + source/provenance; transplants
  have RMSD below threshold; `manifest.csv` complete and all referenced files exist.
- **Phase 5 (`test_analysis`):** every CMM site in the inbox manifest has a parsed result or an
  explicit failed/missing flag; counts reconcile with Phase 3; statistical assumptions checked;
  seeds fixed; reproducible.
- **Phase 6 (`test_figures`):** every figure has a regenerating script + caption stub; numbers
  in figures match `data/processed/`.

## Writing style for gates
Prefer small, deterministic scripts with clear assertion names that read well in output (e.g.
`assert_no_wt_mismatch`, `assert_sets_disjoint`). On failure, print the offending rows/values,
not just a boolean. A gate that can't explain WHY it failed is only half a gate.
__SKILL_MD_EOF__

mkdir -p ".claude/skills/hgvs-normalization"
cat > ".claude/skills/hgvs-normalization/SKILL.md" << '__SKILL_MD_EOF__'
---
name: hgvs-normalization
description: >
  Normalize and validate FBN1 variant nomenclature onto the canonical transcript NM_000138.5 and
  protein P35555, reconcile cDNA (c.) with protein (p.) descriptions, and enforce the wild-type
  residue-identity check that catches transcript/isoform misalignment. Use this skill WHENEVER
  variants must be turned into trustworthy protein positions — "normalize the variants",
  "validate HGVS", "map c. to p.", "check the reference residue", building the labeled
  pathogenic/benign/VUS sets, joining AlphaMissense scores, or Phase 2 of the pipeline. This is
  the single most common place silent errors enter variant work, so use it even when mapping
  "looks trivial" — a position that is off by one isoform poisons every downstream metric.
---

# HGVS normalization & wild-type-identity validation

Goal: convert raw ClinVar/gnomAD variant strings into validated (position, wild-type AA,
mutant AA) tuples on **P35555**, so structural mapping and scoring are trustworthy.

## Pin the reference first
- Canonical transcript: **NM_000138.5** (MANE Select). Protein: **UniProt P35555** (2871 aa,
  signal peptide 1-27). Record the exact versions used in `manifest/`.
- All protein positions in the project are P35555 positions. If a source uses a different
  transcript/isoform, re-map it — never assume equivalence.

## Validate every HGVS string
- Run each variant through **VariantValidator** (REST: `rest.variantvalidator.org`) or
  **Mutalyzer** against NM_000138.5. Confirm the c. and p. descriptions agree; if a source gives
  only c., derive p. from the validated result (don't hand-translate).
- Normalize to a canonical form (e.g., 3'-shifted per HGVS) so duplicates collapse consistently.
- Missense only for the structural analysis: separate out synonymous, nonsense, frameshift,
  splice, and in-frame indels into a labeled "excluded_from_structural" table (kept, not deleted).

## The wild-type-identity check (do not skip)
For every kept missense variant, assert that the wild-type amino acid named in the p. description
equals the residue at that position in the P35555 reference sequence.
- Mismatch → the variant is on a different numbering/isoform, or the HGVS is wrong. Flag it,
  quarantine it in a `wt_mismatch` table, and report counts. Do NOT coerce it.
- Apply the SAME check when joining AlphaMissense (its rows are P35555 substitutions): join on
  (position, wild-type AA), not position alone.

## Build the labeled sets
Produce disjoint sets with criteria written to `resources/reference/inclusion_criteria.md`
(mentor-approved thresholds — ask before finalizing):
- **Pathogenic/likely-pathogenic**, **Benign/likely-benign**, **VUS**.
- ClinVar review-status cutoff (e.g. ≥2 stars); handle conflicts by labeling
  (`conflicting`), never dropping.
- Benign set may combine ClinVar B/LB with gnomAD common variants above the agreed AF cutoff;
  record which source supports each benign call.
- Assert benign ∩ pathogenic = ∅ before proceeding.

## Output
`data/interim/variants_normalized.tsv` with columns: variation_id, hgvs_c, hgvs_p, position,
wt_aa, mut_aa, class, review_status, sources, validation_status. Plus quarantine tables
(`wt_mismatch`, `excluded_from_structural`, `conflicting`). Hand to data-validation for the
Phase 2 gate (real residues, WT matches reference, positions in 1-2871, disjoint sets, no
contradictory labels).
__SKILL_MD_EOF__

mkdir -p ".claude/skills/metal-site-analysis"
cat > ".claude/skills/metal-site-analysis/SKILL.md" << '__SKILL_MD_EOF__'
---
name: metal-site-analysis
description: >
  Place calcium into cbEGF domain structures using the tiered strategy (experimental → AF3
  cofolding → AlphaFill/transplant cross-check), calibrate predicted models against experimental
  Ca geometry, and run the CheckMyMetal validation via the manual mentor hand-off. Use this skill
  WHENEVER the task involves calcium coordination, metal sites, CheckMyMetal, "add Ca to the
  model", transplanting or cofolding calcium, validating coordination geometry, or Phase 4/4.5 of
  the pipeline. Use it even if someone just says "check the calcium binding", because raw
  AlphaFold models contain NO metals and running metal validation on a metal-free model produces
  confident-looking nonsense — this skill exists to prevent exactly that.
---

# Calcium placement & CheckMyMetal validation (tiered + hand-off)

Goal: produce Ca-bearing cbEGF structures whose coordination geometry has been validated, so
variant effects on calcium binding are measured, not guessed. Work per domain / domain pair.

## Core hazard
AF2 models have no ions. A metal-free model sent to CheckMyMetal (CMM) is meaningless. Every
structure entering CMM must contain a Ca ion whose origin is documented as a modeling assumption.

## Tiered Ca sourcing (prefer the cheapest that works)
1. **Experimental (Tier 1).** For domains with deposited Ca-bound structures (Appendix B), use
   the experimental PDB directly — gold standard AND Ca donor.
2. **AF3 cofolding (Tier 2, primary for un-crystallized domains).** Request the mentor run
   AlphaFold-Server/AF3 (or Chai-1/Boltz) on the domain sequence **with one Ca2+ ligand**
   (handoff/structures/REQUESTS.md, §3b). Model wild-type; where feasible also the variant, to
   read predicted coordination change.
3. **AlphaFill / manual transplant (Tier 3, CROSS-CHECK only).** Superpose an AF2 domain model
   onto an experimental Ca site; copy the Ca ion in; record donor accession + superposition RMSD
   to a provenance JSON in `structures/ca_transplanted/`. Compare its Ca position to the
   AF3-placed one — agreement builds confidence, disagreement flags an artifact.

## Calibrate before trusting predictions
Before relying on AF3 for domains without structures, run AF3+Ca2+ on the experimental domains
(cbEGF12-13, 32-33, 9-10) and confirm it reproduces the known Ca geometry (via CMM). Report the
calibration agreement (coordination number match, Ca-position RMSD) as a methods result. Only
then extend AF3 to un-crystallized domains.

## Pick ONE primary pipeline, uniformly
Choose experimental-where-available, else AF3-cofolded, and apply it to ALL domains, recording
per-domain which method/source was used. Mixing methods silently biases the pathogenic-vs-benign
enrichment. The transplant stays a secondary cross-check.

## CheckMyMetal is a MANUAL hand-off (never automate/fake it)
CMM is web-only. Prepare, then stop:
1. Assemble the **submission packet** in `handoff/cmm/inbox/`: one prepared PDB per Ca site
   (experimental + Ca-placed models), descriptive filenames, plus `manifest.csv`
   (filename → domain → site → method/source → intended variant(s)).
2. Write `handoff/cmm/INSTRUCTIONS.md`: URL https://csgid.org/csgid/metal_sites, metal = **Ca**,
   settings, processing order, and which outputs to save into `handoff/cmm/outbox/`.
3. STOP and announce per CLAUDE.md §3b. Resume only when the mentor confirms outbox is populated.

## Parse returned CMM results
Validate every submitted site has a result. Parse coordination number, geometry classification,
bond lengths, valence/nVECSUM, and flags into a tidy table keyed to domain/site/variant. Report
failed/missing submissions explicitly — never fill gaps with guesses.

## Output
`structures/ca_transplanted/*` (+ provenance JSON), the CMM results table, and a per-site
comparison of predicted vs experimental/transplanted Ca geometry, feeding the Phase 5 stats.
__SKILL_MD_EOF__

mkdir -p ".claude/skills/pymol-figures"
cat > ".claude/skills/pymol-figures/SKILL.md" << '__SKILL_MD_EOF__'
---
name: pymol-figures
description: >
  Generate reproducible, publication-quality PyMOL figures of cbEGF domains showing mutated
  residues, calcium ions, and disulfide bonds, by shelling out to headless PyMOL with committed
  .pml scripts. Use this skill WHENEVER the task involves rendering, drawing, or producing a
  structural figure — "make the PyMOL figure", "render the domain", "show the mutation on the
  structure", "figure of the calcium site", or Phase 6 figure generation. Use it even for a
  quick one-off render, because ad-hoc GUI clicking is not reproducible and this project requires
  every figure to regenerate from a script with matching numbers.
---

# Publication PyMOL figures (headless, reproducible)

Goal: every structural figure is produced by a committed script, regenerates identically, and
uses consistent conventions — no GUI, no manual posing.

## Run headless, shell out
- Do NOT `import pymol2` inside the venv. Call open-source PyMOL as a subprocess:
  `pymol -cq figure_<name>.pml`. `-c` = no GUI, `-q` = quiet.
- Keep each figure's `.pml` (or a PyMOL-driving Python script) in `scripts/` or next to its
  output; the figure and its script are delivered together.
- Set a fixed random/seed-independent view: define the camera explicitly with `set_view`
  (capture it once, paste it in) so renders are deterministic across machines.

## Standard cbEGF figure recipe
Within the `.pml`:
- `bg_color white`; `hide everything`; `show cartoon`; thin cartoon for the domain backbone.
- **Calcium:** `show spheres` for `resn CA` (Ca ion), colored a consistent accent; scale
  `set sphere_scale, 0.5, resn CA`.
- **Mutated residue(s):** `show sticks` + label; color by class using a fixed palette
  (pathogenic / benign / VUS) so panels are comparable across figures.
- **Disulfides:** `show sticks, resn CYS and name CA+CB+SG`; optionally `set dash` bonds between
  paired cysteines to highlight the 1-3/2-4/5-6 pattern.
- **Coordinating residues:** show sticks for the Ca-coordinating side chains; optionally draw
  distance objects from those atoms to the Ca ion.
- Ray-trace at print resolution: `set ray_opaque_background, 0`; `ray 2000, 1500`;
  `png figures/<name>.png, dpi=300`.

## Consistency rules (so figures are comparable)
- One shared color palette and one atom/representation convention across ALL figures — define it
  once in a `references/style.pml` snippet and `@`-include it.
- Use SIFTS-mapped numbering when labeling residues in P35555 terms; if the structure uses author
  numbering, translate labels so the figure speaks the manuscript's numbering.
- Every figure needs a caption stub in `results/figure_captions.md` and the numbers/residues in
  it must match `data/processed/` (data-validation gate for Phase 6).

## Output
`figures/<name>.png` (+ the `.pml`/script that made it) and a caption stub. Never leave a figure
without its regenerating script.
__SKILL_MD_EOF__

mkdir -p ".claude/skills/reproducibility"
cat > ".claude/skills/reproducibility/SKILL.md" << '__SKILL_MD_EOF__'
---
name: reproducibility
description: >
  Enforce reproducibility and provenance across the pipeline — pinned environments, fixed seeds,
  per-run logging of tool versions and database access dates, checksums, and a manifest that
  records where every artifact came from (including the non-pip tools PyMOL, DSSP, and FoldX that
  requirements.txt cannot capture). Use this skill WHENEVER setting up runs, writing logs,
  recording provenance, pinning versions, or preparing results for a manuscript — "make this
  reproducible", "record the versions", "write the manifest", "log this run", or the methods/
  provenance parts of any phase. Use it proactively, because a result whose origin isn't recorded
  can't be defended to a reviewer and can't be regenerated.
---

# Reproducibility & provenance

Goal: any result can be traced to its inputs, tools, and code, and regenerated identically. This
matters more here than usual because the environment mixes a Python venv with external binaries.

## Environment pinning
- Python layer: `pip freeze --local > requirements.txt` (the venv uses
  `include-system-site-packages = true`, so `--local` keeps system packages out of the pin).
- **Non-pip tools can't be captured by pip** — record them explicitly in `manifest/tools.md`:
  PyMOL (`pymol -cq -d "print(cmd.get_version()[0])"`), DSSP (`mkdssp --version`), and FoldX
  (`foldx --version`, plus the binary path `~/tools/foldx/foldx_20261231` and its `molecules/`
  location). Capture their versions AND paths.
- Record database/tool versions used for data: ClinVar release date, gnomAD v4.1.0, UniProt
  release, AlphaFold DB version, AlphaMissense Zenodo record 8208688 (v1.0.0), and the
  CheckMyMetal run date.

## Provenance manifest (every artifact)
For each downloaded or generated file write a `manifest/` entry: source (URL/endpoint/query or
generating script), accession/version, retrieval or generation date, row/atom count, and
**sha256**. Derived tables additionally record the input files + script that produced them. Never
overwrite a raw file in place — raw is immutable once written.

## Determinism
- Fix all random seeds (numpy, python `random`, any stochastic predictor). State the seed in the
  script and the log.
- Cache raw network pulls in `data/raw/` and run analysis off the cache so results don't depend
  on live network state.
- Make stages idempotent: skip work when a valid output + matching checksum already exists.

## Per-run logging
Each script writes a timestamped log to `logs/`: start/end time, tool versions, input checksums,
row counts in/out, seed, and any rows quarantined. A run you can't audit later is a run you
can't trust.

## Secrets discipline
The NCBI API key and any credentials live in `.env` (git-ignored). They must never appear in
logs, manifests, committed files, or printed output. When recording provenance for an
API-sourced file, record the endpoint and query — never the key.

## Manuscript readiness
Maintain `results/LIMITATIONS.md` capturing modeling assumptions that affect interpretation:
the Ca-placement method per domain (experimental vs AF3 vs transplant), AlphaFold confidence,
ClinVar label quality, benign-set ascertainment bias, and absence of glycosylation/
β-hydroxylation (not required for Ca binding). Reviewers will ask; have it ready.
__SKILL_MD_EOF__

mkdir -p ".claude/skills/structure-mapping"
cat > ".claude/skills/structure-mapping/SKILL.md" << '__SKILL_MD_EOF__'
---
name: structure-mapping
description: >
  Retrieve FBN1 structures and map residues between UniProt P35555 numbering and PDB author
  numbering using SIFTS, with round-trip validation of residue identity. Use this skill WHENEVER
  the task touches structures — "find the FBN1 PDB structures", "which domains have structures",
  "map my variant onto the structure", "get the AlphaFold model", "SIFTS mapping", handling
  pLDDT, or the structure-sourcing part of Phase 4. Use it even when a residue index "obviously"
  lines up, because PDB fragments carry local numbering that rarely equals UniProt numbering, and
  a silent off-by-N mapping mislabels every structural metric and every FoldX mutation.
---

# Structure retrieval & SIFTS residue mapping

Goal: for each cbEGF domain of interest, obtain the right structure(s) and a validated mapping
from P35555 positions to structure positions, so distances, coordination, and ΔΔG are computed on
the correct atoms.

## Enumerate what exists (don't trust memory)
- Query the **PDBe SIFTS / Graph API** (or RCSB) for all PDB entries mapped to **P35555**. This
  returns the current list with per-entry UniProt residue ranges. FBN1 entries are small
  fragments (e.g. cbEGF pairs / TB regions) — most of the 43 cbEGF domains have NO structure.
- Record each entry's covered domain range, method, resolution, and whether it is Ca-bound.
  Known Ca-bound anchors: 1EMN/1EMO (cbEGF32-33), 1LMJ (cbEGF12-13), 2W86 (cbEGF9-hyb2-cbEGF10),
  and a cbEGF22-TB4-cbEGF23 fragment — see CLAUDE.md Appendix B.

## SIFTS mapping (the core service)
- Download the **SIFTS** mapping for each chosen PDB (`ebi.ac.uk/pdbe/api/mappings/P35555` or the
  per-entry SIFTS XML). Build a lookup: (PDB id, chain, author_resnum) ↔ (P35555 position).
- **Round-trip validate:** take a known P35555 position → map to author numbering → read the
  residue in the coordinate file → assert its three-letter code matches the expected amino acid.
  If identity fails, the mapping is wrong — stop and report; do not compute on it.
- Store the mapping table per structure in `structures/pdb/<id>_sifts.tsv`.

## AlphaFold / predicted models
- AF2 full-length AF-P35555 exists but inter-domain geometry is low-confidence and it has NO
  metals. Use **per-domain** slices, not the whole chain, for geometry.
- Carry **per-residue pLDDT** alongside coordinates; flag low-confidence residues so downstream
  metrics can be filtered or caveated.
- For predicted models the "author numbering" typically equals P35555 numbering, but STILL run
  the identity check — don't assume.

## Handoff for structures you can't fetch
If a needed structure isn't programmatically available (AF3 Ca-cofolded model, AlphaFill model,
a manual PDB pick), do NOT fabricate it. Append an itemized request to
`handoff/structures/REQUESTS.md` (domain, residue range, tool + exact settings, destination
folder) and STOP for the mentor per CLAUDE.md §3b.

## Output
Chosen structures in `structures/pdb/` and `structures/alphafold/`, each with a validated
`_sifts.tsv` mapping and a manifest entry. Feed these to the metal-site-analysis skill (Ca
placement + CheckMyMetal) and to Phase 5 (SASA, distances, FoldX numbering).
__SKILL_MD_EOF__

mkdir -p ".claude/skills/variant-collection"
cat > ".claude/skills/variant-collection/SKILL.md" << '__SKILL_MD_EOF__'
---
name: variant-collection
description: >
  Collect FBN1 (UniProt P35555) variants and their annotations from ClinVar, UniProt, gnomAD v4,
  and the AlphaMissense subset for the Marfan structural-bioinformatics pipeline. Use this
  skill WHENEVER the task involves pulling, downloading, refreshing, or assembling variant data
  from any of these sources — including phrases like "get the ClinVar variants", "pull gnomAD
  frequencies", "load AlphaMissense", "collect the FBN1 mutations", or Phase 1 of the pipeline.
  Use it even when the user just says "start collecting data" in this project, because the
  retrieval details (rate limits, review-status capture, keeping VUS/conflicts, manifest
  writing) are easy to get wrong and corrupt everything downstream.
---

# Variant collection (ClinVar / UniProt / gnomAD v4 / AlphaMissense)

Goal: land raw, fully provenanced variant data in `data/raw/` so later phases work off a stable
cache. Never analyze in this phase — only retrieve, validate shape, and record provenance.

## Golden rules
- **Cache-first.** Write every raw pull to `data/raw/` unchanged, then work off the cache. A
  re-run must not re-download if a valid cached file + checksum exists (idempotent).
- **Keep everything, label nothing away.** Retain VUS and conflicting classifications — they are
  data. Never silently drop rows; add a status column instead.
- **Provenance is mandatory.** For each source write a `manifest/<source>.md` entry: endpoint/URL,
  exact query, accession/version, retrieval date, row count, sha256.
- **Secrets stay in `.env`.** Load the NCBI key via `python-dotenv`; never print or log it.

## ClinVar (FBN1, Gene ID 2200)
Two access paths — pick based on need, prefer the bulk file for completeness:
- **E-utilities** (`Bio.Entrez` esearch→efetch). Set `Entrez.email` and `Entrez.api_key` from
  `.env` (raises limit to ~10 req/s). Batch with `WebEnv`/`query_key`; sleep between requests;
  retry on HTTP 429/500 with backoff.
- **Bulk** `variant_summary.txt.gz` (NCBI ClinVar FTP), filtered to `GeneSymbol == FBN1`.
Capture per variant: VariationID, HGVS **c.** and **p.**, ClinicalSignificance,
**ReviewStatus (gold stars)**, Condition(s), Assembly, LastEvaluated. Verify the endpoint with
one tiny test query before bulk-pulling.

## UniProt P35555
REST (`https://rest.uniprot.org/uniprotkb/P35555`) in JSON + GFF. Extract: sequence, and domain
features (EGF / cbEGF / TB / hybrid boundaries), disulfide-bond annotations, and calcium-binding
site annotations. This is the source of truth for domain mapping downstream — store it verbatim.

## gnomAD v4
GraphQL API (`https://gnomad.broadinstitute.org/api`). Query FBN1 variants for allele frequency
(use the joint/grpmax filtering AF). Record dataset version (v4.1.0). Used to build/vet the
benign set and flag common variants.

## AlphaMissense
Read the pre-placed P35555 subset in `resources/databases/` (header + `^P35555` rows). Do NOT
re-download the multi-GB release. Columns: UniProt ID, protein substitution, am_pathogenicity,
am_class. Defer the wild-type-identity join check to the hgvs-normalization / data-validation
skills — but note here that positions are on canonical P35555.

## Output of this phase
`data/raw/clinvar_fbn1.tsv`, `data/raw/uniprot_P35555.json` (+ features), `data/raw/gnomad_fbn1.tsv`,
and a confirmed-present `resources/databases/AlphaMissense_FBN1_P35555.tsv`, each with a manifest
entry. Then hand to the data-validation skill for the Phase 1 gate (non-empty, expected columns,
HGVS parses, no duplicate VariationIDs, stars captured).

## Retry/adapter pattern (apply to every source)
Wrap each source in a small `scripts/` adapter with: timeout, exponential-backoff retry, a
"skip if cached + checksum matches" guard, and a manifest write on success. If a source is
unreachable, record the gap and continue — do not fabricate rows to fill it.
__SKILL_MD_EOF__

echo "Done. Created:"
find .claude/skills -name SKILL.md | sort
