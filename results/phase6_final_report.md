# Phase 5/6 final report — 2026-08-29

Closes the pipeline. Batch 4 arrived, the geometric result became a rate, and the manuscript,
supplement, statistics and figure set were rebuilt on the final tables.

---

## What was produced

| artifact | script | rows / size | what it is |
|---|---|---|---|
| `data/processed/af3_batch4_geometry.tsv` | `40_af3_batch4_geometry.py` | 180 rows (60 jobs × 3 sites) | wild-type vs mutant coordination at every site of every batch-4 construct |
| `data/processed/batch4_rate.tsv` | `41_batch4_rate.py` | 57 variants | one row per variant at its cognate site, with the disruption call |
| `results/phase5_batch4_geometry.md` | `40_…` | — | per-variant table, cognate sites |
| `results/phase5_batch4_rate.md` | `41_…` | — | the rate, its threshold, magnitude, donor chemistry |
| `results/statistics_final.tsv` / `.md` | `44_statistics_final.py` | 26 tests, 5 families | the single statistics table the manuscript cites |
| `figures/v2_fig7_rate.*` | `42_v2fig7_rate.py` | 4 panels | **Figure 4** — the rate |
| `figures/v2_fig8_structure.*` | `46_v2fig8_structure.py` | 3 panels | **Supplementary Figure S2** — the mechanisms on 2W86 |
| `results/supplementary_tables_v2.xlsx` | `43_supplementary_v2.py` | 10 sheets, 1.25 MB | S1–S9 + README, every column defined |
| `results/manuscript_final.md` | written | 478 lines | full draft with methods and 16 references |
| `results/FBN1_manuscript_final.docx` | `45_manuscript_final_docx.py` | 2.0 MB | the same, figures placed inline |

Updated in place: `results/LIMITATIONS.md` (sections 1, 3, 5, 7 revised; 13–15 added),
`results/v2_figure_captions.md`, `results/SESSION_CONTEXT.md`, `tests/test_v2.py`.

## Counts, and what was dropped

- **60 of 60** batch-4 jobs found, validated and measured. None excluded for a technical reason:
  every model set was complete (5/5), every sequence differed from its wild type at exactly one
  position, and every substitution matched the manifest and P35555.
- **3 of 60 variants dropped from the rate** — D490Y, N261S, E487A — because their cognate sites
  are on the cbEGF1–3 construct, which fails the wild-type quality filter (one site has a seed
  spread of 39% of its own valence; another is occupied in only 4 of 5 wild-type models). The
  filter is computed on wild-type evidence alone, before any mutant is read.
- **0 sites vacated.** In every mutant, an ion still occupied every wild-type site.
- Final groups at the cognate site: 33 pathogenic Ca²⁺ ligands, 10 composition controls, 14
  negative controls.
- Null: **123** measurements (109 bystander sites + 14 negative-control cognate sites), giving a
  threshold of **9.14%** of site valence.

## Headline results

| | value |
|---|---|
| pathogenic Ca²⁺-ligand variants disrupting their site | **15/33 = 45%** (95% CI 30–62%) |
| controls disrupting | 1/24 |
| odds ratio | 19.2, Fisher p = 7.0 × 10⁻⁴ |
| magnitude, ligands vs controls | −5.7% vs −1.7% / +0.1%; Cliff's δ = −0.473 [−0.717, −0.192] |
| donor removed vs donor retained | 12/16 (75%) vs 3/17 (18%); OR 14.0, p = 0.002 |
| local fold, ligands vs controls | 0.29 Å vs 0.24 Å backbone RMSD |
| Rosetta ΔΔG vs coordination change | Spearman ρ = −0.089, q = 0.62 |
| AlphaMissense, disrupted vs intact | 0.982 vs 0.982 |

Two of these are new claims the earlier drafts could not make: the **rate**, and the fact that
**donor chemistry predicts it while neither ΔΔG nor AlphaMissense does**.

## Corrections made to earlier text

- The draft said "six of 43 domains have experimental structures". It is **eight** (cbEGF9, 10,
  12, 13, 22, 23, 32, 33). Corrected in the manuscript and `LIMITATIONS.md`.
- The Figure 5 calibration claim "the ion was placed more accurately than the protein around it,
  in all eight cases" holds on the top-ranked model of each construct, which is what the figure
  plots. Summarising each domain by its median over five models gives **seven of eight** (0.55 Å
  against 1.35 Å), the exception being cbEGF12 where the two are equal to within 0.02 Å. Both are
  now stated.
- `phase5_statistics_v2.md` predates the corrected PyRosetta ΔΔG and is superseded by
  `statistics_final.tsv`. The Rosetta numbers the manuscript quotes had been living only in a
  figure sidecar; they now have a statistics table.

## Gate

`tests/test_v2.py`, **14/14 passing**. Four checks added: A7 batch-4 completeness, A8 the null
threshold re-derives from the data, A9 every disruption call follows from that threshold, A10 the
statistics table is coherent and never compares FoldX to Rosetta. F1–F4 now cover eight figures.

## Open questions for the mentor

1. **Eighteen pathogenic Ca²⁺-ligand variants keep their calcium**, and they are overwhelmingly
   the conservative substitutions (D→N, D→E). Is the right reading β-hydroxylation loss, an
   affinity change below what a modelled valence can show, or ClinVar overclassification? The
   manuscript states all three and picks none. This is the one place a reviewer will push.
2. **Reference 16 (ClinGen FBN1 VCEP)** is cited from its DOI because the supplied PDF is a
   different paper. The AlphaFold 3 methods paper is missing from `resources/papers/` entirely and
   must be added before submission.
3. **Is Boltz-2 still wanted?** It is the only outstanding independent check on AF3's calcium
   placement, and the environment is still broken.
4. **Should v1 be retained as a supplementary sensitivity analysis**, as the earlier draft
   proposed, or dropped? Nothing currently cites it.
