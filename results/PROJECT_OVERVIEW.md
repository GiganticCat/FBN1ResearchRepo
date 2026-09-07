# FBN1 / Marfan structural bioinformatics — project overview

**State as of 2026-08-02:** Phases 0–3 complete, all four gates passing. Phase 4 not started.

The scientific question: **are pathogenic FBN1 missense variants enriched at structurally
critical positions in cbEGF domains** — the calcium-coordinating residues and the six conserved
cysteines — relative to benign variants and to the background of all possible missense changes?
And how much do those variants perturb Ca²⁺ geometry, disulfide bonding and stability?

---

## 1. Directory map

```
fbn1-marfan/
├── CLAUDE.md              the operating brief (project constitution)
├── requirements.md        human setup instructions
├── requirements.txt       pinned Python environment (33 packages)
├── .env                   NCBI key + email — git-ignored, never logged
│
├── resources/             READ-ONLY human input
│   ├── papers/            19 literature PDFs
│   ├── databases/         AlphaMissense P35555 subset (+ 1.2 GB source archive)
│   └── reference/         lit_notes.md · structures.md · inclusion_criteria.md
│
├── data/
│   ├── raw/               immutable downloads (447 MB, git-ignored)
│   ├── interim/           cleaned/derived intermediates
│   └── processed/         analysis-ready tables + data dictionary
│
├── structures/
│   ├── pdb/               11 experimental FBN1 entries (17 MB)
│   ├── alphafold/         4 AF3 Ca²⁺-cofolded models (63 MB)
│   └── ca_transplanted/   empty — Phase 4
│
├── scripts/               7 numbered pipeline stages
├── tests/                 4 validation gates
├── handoff/               structures/ + cmm/ — mentor hand-off packets
├── manifest/              provenance for every artifact
├── logs/                  timestamped run transcripts
├── results/               phase reports + cross-tabs
└── figures/               empty — Phase 6
```

---

## 2. Scripts (`scripts/`)

| Script | Phase | What it does |
|---|---|---|
| `00_env_check.py` | 0 | Verifies scaffold, imports, PyMOL/DSSP/FoldX callability, `.env` presence, AlphaMissense subset, and one live shape-checked query per database. Writes `manifest/tools.md` + a machine-readable JSON the gate consumes. |
| `01_fetch_clinvar.py` | 1 | Downloads ClinVar's bulk `variant_summary.txt.gz` (442 MB) and filters to FBN1. Refuses to run if any ReviewStatus lacks a star mapping. Cache-guarded. |
| `02_fetch_uniprot.py` | 1 | Caches P35555 as JSON/GFF/FASTA; derives the domain table (with the cbEGF 1–43 index) and the disulfide table, verifying every disulfide endpoint really is Cys. |
| `03_fetch_gnomad.py` | 1 | gnomAD v4 GraphQL pull; captures exome/genome/joint counts and the **grpmax filtering AF**. |
| `04_prepare_af3_requests.py` | 4 (early) | Built the AlphaFold-Server calibration packet: sliced domain sequences from the cached UniProt record and emitted uploadable AF3 job JSON. |
| `05_normalize_variants.py` | 2 | Parses HGVS, classifies consequence, runs the **wild-type-identity check** against P35555, joins gnomAD + AlphaMissense, and builds the labelled sets. |
| `06_annotate_domains.py` | 3 | Derives cbEGF calcium sites from cysteine anchors, classifies disulfides, and annotates every variant with domain/site/VCEP-PM1 flags. |

## 3. Gates (`tests/`) — all must pass to advance

| Gate | Checks | Status |
|---|---|---|
| `test_env.py` | scaffold, imports, tool callability, secrets, DB reachability, one lit-note per PDF | **8/8** |
| `test_variant_ingest.py` | row envelopes, schema, ID uniqueness, star mapping, HGVS parse rate, VUS/conflict retention, FAF populated | **10/10** |
| `test_normalization.py` | WT identity, position range, missense purity, count reconciliation, set disjointness, hidden contradictions, AM join identity | **10/10** |
| `test_annotation.py` | cysteine annotations, 258-cysteine count, disulfide pattern, consensus motif, **agreement with measured Ca ligands**, cross-tab sanity | **9/9** |

Every gate was **self-tested against simulated failures** — corrupted residues, inverted labels,
dropped rows, off-by-one calcium ligands — to confirm it actually fails when it should.

---

## 4. Data lineage

```
ClinVar bulk (2026-07-28)  ──┐
UniProt P35555 v264        ──┤
gnomAD v4 (gnomad_r4)      ──┼─→ 05_normalize ─→ variant_master.tsv     (4,451 × 34)
AlphaMissense (Zenodo)     ──┘                          │
                                                        ↓
11 PDB entries + 4 AF3 models ─→ 06_annotate ─→ variant_annotated.tsv (4,451 × 56)
```

### From 9,327 ClinVar GRCh38 records

| Outcome | n |
|---|---|
| **Missense, WT-verified** | **4,451** |
| splice/intronic (no protein consequence) | 1,899 |
| synonymous | 1,295 |
| frameshift | 975 |
| nonsense | 491 |
| in-frame indel | 162 |
| other / unknown protein effect | 3 |
| genomic-only description | 51 |
| **wild-type mismatch** | **0** |

Nothing was deleted; every excluded row carries a reason and the counts reconcile exactly.

### Labelled sets

| Set | Primary (≥2★) | Sensitivity (≥1★) |
|---|---|---|
| pathogenic | 578 | 1,409 |
| benign | 34 | 76 |
| VUS | 823 | 2,295 |
| conflicting | 0 | 492 |

---

## 5. Key scientific findings so far

1. **Reference mapping is clean.** 0/4,451 wild-type mismatches; 100/100 Mutalyzer agreement;
   0 gnomAD protein-description disagreements.
2. **The calcium consensus was derived, not assumed.** Anchored on cysteines, calibrated against
   Ca²⁺ ligands measured in 2W86/1LMJ/1UZJ/1EMN, and validated against the VCEP motif (43/43
   domains, zero exceptions) and 9/9 published critical positions from Baudhuin 2019.
3. **258 cbEGF cysteines** — exactly the ClinGen VCEP's stated figure.
4. **All 129 cbEGF disulfides** follow C1–C3 / C2–C4 / C5–C6 with no exceptions.
5. **Ca²⁺ stoichiometry is 1 per cbEGF domain, 0 for TB/hybrid** — counted from coordinates.
6. **The benign set is small for real biological reasons.** Only 24 gnomAD missense variants
   clear the benign filtering-AF threshold, and all 24 are already in the ClinVar set. FBN1 is
   too constrained to have many common missense variants.

### Bugs caught that would have corrupted the study

- **`2200[gene_id]` is not a ClinVar field.** NCBI silently rewrote it to a free-text
  `2200[All Fields]` search that matched COL2A1 and ALK records. Corrected to `FBN1[gene]`; the
  env check now fails if a query translation degrades to `[All Fields]`.
- **gnomAD's `faf95_joint` is null for every FBN1 variant.** The populated field is
  `joint.fafmax.faf95_max`. Trusting the empty one would have produced a benign set with no
  frequency evidence at all, looking like a clean result rather than a bug.
- **AlphaMissense position 472 carries two wild-type residues** (C and Y; P35555 says Y). The
  join uses (position, WT, ALT) so the 19 spurious rows are rejected.

---

## 6. ⚠️ The circularity caveat

**91% of pathogenic variants sit at a VCEP critical residue vs 0% of benign.** This is partly an
artifact: ClinVar's pathogenic calls were themselves made using PM1 — the very critical-residue
rule being tested.

An orthogonal check (AlphaMissense, which never used FBN1 PM1 rules) reproduces the *direction*
cleanly — 33.1% of AM-pathogenic at cbEGF cysteines vs **0.0%** of AM-benign — so the biology is
real. But the *magnitude* from ClinVar labels is inflated and must not be headlined as an
independent discovery.

**Consequence for Phase 5:** the headline result should be the **physics** — Ca²⁺ coordination
geometry (CheckMyMetal), FoldX ΔΔG, solvent accessibility, disulfide disruption — which are
computed from structures and are fully independent of any classification scheme.

---

## 7. Open items

| # | Item | Needs |
|---|---|---|
| 1 | **Phase 4 structure plan** — which cbEGF domains to model | mentor decision (standing rule: ask before building requests) |
| 2 | **Tier-3 cross-check replacement** — AlphaFill is unavailable for FBN1 (no AlphaFold DB entry, 2871 aa > 2700 cap) | mentor decision |
| 3 | **VCEP Supplementary Table 1** — would resolve the 655 vs 633 PM1 position count | optional download |
| 4 | Primary Jensen 2009 (2W86) paper | optional |
| 5 | `data/raw/uniprot_P35555_2026-08-02.json` duplicates `uniprot_P35555_v264.json` | cosmetic |

## 8. Reproducibility

- `requirements.txt` pins the Python layer; `manifest/tools.md` records PyMOL 3.1.0, mkdssp
  4.5.7, FoldX 5.1 with paths and versions.
- Every download carries a manifest entry with source URL, version/release, date, row count and
  checksum. All fetchers are idempotent.
- Seeds fixed (Mutalyzer sampling: 20260802). AF3 model seeds recorded.
- The NCBI key has never appeared in any log, manifest, report or data file — audited by grep
  after every phase.
