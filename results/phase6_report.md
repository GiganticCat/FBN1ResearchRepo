# Phase 6 report — Figures & manuscript tables

**Date:** 2026-08-08 (UTC) · **Status:** complete, gate PASSED (9/9) · **All seven phases done.**

Six publication figures, each from one seeded script, plus the supplementary workbook. Framing
follows your 2026-08-02 direction: Figures 2–4 carry the biophysics, Figure 6 reports enrichment
as supporting evidence with its circularity caveat drawn into the figure rather than left to the
caption.

---

## 1. What was produced

| Figure | Content | Script |
|---|---|---|
| 1 | Variant landscape along P35555 by clinical class + domain architecture + per-cbEGF breakdown | `13_fig1_landscape.py` |
| 2 | **Headline.** ΔΔG, AlphaMissense, their joint distribution, and the FoldX disulfide term | `14_fig2_mechanism.py` |
| 3 | Burial control — the calcium result is not a solvent-exposure artifact | `15_fig3_burial.py` |
| 4 | PyMOL renders of 2W86: fragment, calcium site, AF3 overlay, pathogenic positions | `16_fig4_structure.py` |
| 5 | Calcium-placement accuracy and the CheckMyMetal comparison | `17_fig5_calibration.py` |
| 6 | Enrichment forest plots, supporting only | `18_fig6_enrichment.py` |

Each figure is written as a 300 dpi PNG **and** a vector PDF with editable text (`pdf.fonttype
42`), so a journal can typeset from it. Captions are in `results/figure_captions.md`.

`results/supplementary_tables.xlsx` (`19_supplementary_tables.py`) holds eight sheets: a README
carrying the three governing caveats, the 4,451-variant master table (S1), the 751 structurally
covered variants (S2), the 16 CheckMyMetal sites (S3), the 20 statistical tests (S4), the AF3
calibration (S5), the leave-one-out transplants (S6), and a data dictionary (S7).

## 2. Two defects found and fixed

### 2.1 A stratification bug in Phase 5 — corrected, conclusion unchanged

Building Figure 3 I could not reproduce the Phase 5 buried-stratum counts. The cause was in
`12_statistics.py`:

```python
(f(v.get("rsa_pct")) or -1) >= lo and (f(v.get("rsa_pct")) or 999) < hi
```

`0.0` is falsy in Python, so a residue with RSA **exactly 0.0** — completely buried — was
rewritten to `-1`/`999` and fell out of *every* stratum. Four fully buried "other cbEGF"
variants were silently dropped from the buried bin, and they are the most destabilising ones
(median ΔΔG 13.90).

Fixed to test for `None` explicitly. Re-running Phase 5: **one row of twenty changed**.

| | before | after |
|---|---|---|
| buried, other cbEGF | n = 170, median 3.54 | n = 174, median **3.62** |
| Cliff's δ | −0.557 | **−0.562** |
| q (BH) | 2.5×10⁻⁶ | 2.5×10⁻⁶ |

The effect moved *away* from zero, so the bug had been making the result slightly weaker, not
stronger. No calcium-consensus variant has RSA 0.0, so that arm was never affected. The Phase 5
gate still passes 10/10 including byte-identical reproducibility. `phase5_report.md` and
`FULL_PROJECT_REPORT.md` are updated and carry a note recording the correction.

### 2.2 A wording claim the statistics do not support

`phase5_report.md` §1 says calcium-ligand ΔΔG is "statistically indistinguishable" from other
cbEGF residues. The table it cites does not say that. It reports a **small but significant**
effect in the *opposite* direction to destabilisation: direct Ca²⁺ ligands versus non-ligands
gives Cliff's δ = −0.23, q = 1.1×10⁻³, i.e. ligands are slightly *cheaper* to mutate, not
equal-cost.

I have not rewritten that sentence, because it is your prose and the fix is a judgment call
about emphasis. Figure 2 annotates the measured effect, with the words "lower, not higher".
**The finding is unaffected and arguably strengthened** — calcium-ligand loss does not
destabilise the fold — but "indistinguishable" should become something like "no more
destabilising than an ordinary cbEGF residue, and by a small margin less" before submission.

## 3. Judgment calls worth your review

- **Structure for Figure 4: 2W86.** X-ray at 1.8 Å, crystallised in calcium rather than
  titrated, carries the most structurally covered variants (254), and is where AF3 calibrated
  best (0.13 Å on cbEGF9). An NMR construct would have shown the reference's own coordinate
  spread rather than the method's accuracy. 2W86's two non-native cloning residues are removed
  before rendering.
- **Two different calcium classes appear in the figures, and they are not interchangeable.**
  Figures 2–3 (panels A–B) use the *structure-measured* `is_direct_ca_ligand` — an O or N within
  3.2 Å of the modelled ion, n = 78. Figure 1C and Figure 3C use the *annotation-level*
  `is_ca_consensus`, defined for all 4,451 variants, n = 414. Figure 3C uses the annotation-level
  class because that is what the Phase 5 test used. Every panel title and the captions say which
  is in play; conflating them would be exactly the "domain-resident vs calcium-coordinating"
  error CLAUDE.md warns about.
- **UniProt calls the hybrid domain "TB 4".** The construct is cbEGF9–hyb2–cbEGF10 in the
  literature, but UniProt annotates residues 851–902 as a TB feature. Figure 4 carries both
  names rather than silently picking one.
- **Figure 3C's exposed stratum is drawn empty.** No calcium-consensus variant has RSA ≥ 50%,
  so that cell contains only the comparison group. Showing the absence is evidence against the
  exposure objection; hiding the stratum would have looked tidier and said less.
- **Figure 1's benign lane is nearly empty.** That is the 34-variant ceiling from
  `LIMITATIONS.md` §2, drawn rather than hidden behind a shared axis.

## 4. Reproducibility and provenance

- **One seed**, `20260802`, defined in `scripts/figstyle.py` and used by every figure script.
  The gate refuses a script that does not draw from it.
- **One palette.** Every colour comes from `figstyle.PALETTE`; the gate scans the figure scripts
  and fails on any hex literal outside it. Each palette group was checked against the data-viz
  validator for colour-vision deficiency, lightness banding and contrast on the print surface;
  the two slots that fall below 3:1 contrast (`#eda100` uncertain-significance, `#1baf7a`
  calcium) carry direct n/median labels everywhere they appear, so identity is never
  colour-alone.
- **PyMOL runs headless** via subprocess on generated `.pml` scripts kept in `figures/pml/`.
  Camera state is set with `orient` on a named selection plus fixed `turn`/`zoom` — deterministic
  across machines — and each panel writes its resulting view matrix into the run log.
- **`requirements.txt` re-pinned** with `pip freeze --local`. It gained `openpyxl==3.1.5` and
  `et_xmlfile==2.0.0`, installed this phase to write XLSX. The re-pin also captured
  `gemmi==0.7.5`, which was present in the venv but had been missing from the pin; no script
  imports it directly.
- **Every figure carries a `numbers.json` sidecar** listing each value it annotated, with the
  script, the seed and the input files. This is what makes "the numbers in the figure match the
  data" checkable rather than asserted.

## 5. Gate

`tests/test_figures.py` — **9/9**:

| Check | What it proves |
|---|---|
| `figures_exist` | every figure has a PNG and a vector PDF, both non-trivial |
| `scripts_exist` | each figure names the script that made it, and it is present |
| `captions_present` | one caption per figure, naming its file and its script |
| `numbers_match_data` | **104 recorded values re-derive** from `data/processed/` and `results/` |
| `pymol_panels` | the four `.pml` scripts and panels exist, load the UniProt-renumbered packet, strip the cloning residues; PyMOL is callable |
| `supplementary` | eight sheets, row counts match source, no undocumented column, 200 sampled rows agree with `data/processed/` |
| `shared_palette` | no figure script hard-codes a colour outside `figstyle.PALETTE` |
| `seeds_fixed` | every figure script uses the shared seed |
| `reproducible` | all six scripts re-ran and reproduced identical recorded numbers |

**The gate was proven able to fail.** Five simulated defects were each caught: a corrupted
recorded median (4.19 → 4.99), a deleted caption section, a missing `.pml`, a supplementary
column stripped from the dictionary, and an off-palette colour introduced into a figure script.
`shared_palette` also caught three *real* violations on its first run — two ad-hoc blues and one
grey I had written by hand. Rather than widen the check I gave the grey a named palette entry
that PyMOL and the legend now both read from, and replaced the second blue with a filled/hollow
marker distinction: two steps of one hue cannot clear both the palette's chroma floor and its
normal-vision separation floor, so the encoding, not the rule, had to change.

## 6. Open items for you

1. **The "indistinguishable" wording in `phase5_report.md` §1** — §2.2 above. Yours to decide.
2. **Structural coverage stays at 17%** per your instruction; no AF3 requests were prepared.
3. **The VCEP PM1 count discrepancy (655 vs 633)** remains open per your instruction; it is not
   referenced in any figure.
4. **The Jensen 2009 (2W86) primary paper** is still absent from `resources/papers/`. Figure 4
   is built entirely on the deposited coordinates and the PDB entry metadata, so nothing in it
   depends on that PDF, but the methods section will want the citation.

## 7. Outputs

| Path | Contents |
|---|---|
| `figures/fig1…fig6.{png,pdf}` | six publication figures, 300 dpi raster + vector |
| `figures/*.numbers.json` | every value each figure annotates, for the gate |
| `figures/pml/*.pml`, `figures/panels/*.png` | the PyMOL scripts and rendered panels behind Figure 4 |
| `results/figure_captions.md` | manuscript-ready captions + the supplementary sheet index |
| `results/supplementary_tables.xlsx` | eight sheets, 944 kB |
| `scripts/figstyle.py` | shared palette, matplotlib defaults, data joins, the numbers recorder |
| `scripts/13_…19_*.py` | one script per figure, plus the workbook builder |
| `tests/test_figures.py` | the Phase 6 gate |
| `logs/13_…19_*.log` | run transcripts, including PyMOL view matrices and versions |
