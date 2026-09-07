# Figure captions

Manuscript-ready captions for the Phase 6 figures. Every number quoted here is recorded in the
matching `figures/<name>.numbers.json` and re-derived from `data/processed/` by
`tests/test_figures.py`; if a caption and the data ever disagree, the gate fails.

Each figure is produced by exactly one script and regenerates identically from it (seed
`20260802`). Figures are written as 300 dpi PNG and as vector PDF with editable text.

| Figure | File | Script |
|---|---|---|
| 1 | `figures/fig1_variant_landscape.png` | `scripts/13_fig1_landscape.py` |
| 2 | `figures/fig2_mechanism_dissociation.png` | `scripts/14_fig2_mechanism.py` |
| 3 | `figures/fig3_burial_control.png` | `scripts/15_fig3_burial.py` |
| 4 | `figures/fig4_structure.png` | `scripts/16_fig4_structure.py` |
| 5 | `figures/fig5_ca_placement.png` | `scripts/17_fig5_calibration.py` |
| 6 | `figures/fig6_enrichment.png` | `scripts/18_fig6_enrichment.py` |

---

## Figure 1 — Pathogenic, uncertain and benign FBN1 missense variants along fibrillin-1

**(A)** ClinVar missense variants mapped to UniProt P35555 (2,871 residues) on transcript
NM_000138.5, shown in three lanes by clinical class at the primary ≥2-star review tier:
pathogenic/likely pathogenic (n = 578, at 398 distinct positions), uncertain significance
(n = 823, 712 positions) and benign/likely benign (n = 34, 32 positions). Stem height is the number of
variants at that position (tallest stem: 7 pathogenic, 4 uncertain, 2 benign). The shaded band
marks residues 659–1362, the exon 24–32
neonatal/severe cluster. The benign lane is sparse because it cannot be otherwise: FBN1 is
under strong constraint and only 34 benign missense variants reach ≥2 stars, with no
gnomAD-only additions available (`results/LIMITATIONS.md` §2).
**(B)** Domain architecture from UniProt: 43 calcium-binding EGF-like (cbEGF) domains, 4
non-calcium-binding EGF-like domains, and 9 TB/hybrid domains. Numerals give the cbEGF index.
The track is drawn in lightness steps rather than colour so it does not compete with the
palette that carries meaning in Figures 2–6.
**(C)** The 425 pathogenic variants that fall inside a cbEGF domain, per domain, split by the
site they affect: cbEGF cysteine-removing (304), calcium-consensus residue (61), other cbEGF
residue (60). Cysteine substitutions dominate in nearly every domain. Panel C uses the
*annotation-level* calcium class (`is_ca_consensus`), which is defined for all 4,451 variants;
Figures 2–3 use the narrower structure-measured class and are not interchangeable with it.

## Figure 2 — Two Marfan mechanisms, two different kinds of damage

The 751 variants with structural coverage, grouped into direct Ca²⁺ ligands (n = 78; a side
chain or backbone atom within 3.2 Å of the modelled ion), cbEGF cysteine-removing variants
(n = 146) and all other cbEGF residues (n = 300). Boxes show median and interquartile range;
points are individual variants; the count of points above each axis limit is stated on the panel.

**(A)** FoldX ΔΔG. Cysteine removal is expensive (median 4.19 kcal/mol; Cliff's δ = +0.64
versus all other missense, q = 9×10⁻³²), driven by the loss of a disulfide. Calcium-ligand loss
is not (median 0.76 versus 0.56 for other cbEGF residues); the measured effect versus
non-ligands is small and in the *opposite* direction to destabilisation (δ = −0.23,
q = 1.1×10⁻³).
**(B)** AlphaMissense pathogenicity for the same variants. Calcium ligands score 0.975 against
0.314 for other cbEGF residues (δ = +0.55, q = 2×10⁻¹³) — damaging by a predictor that never
saw a force field.
**(C)** The two together. 47 of 78 calcium-ligand variants fall in the shaded quadrant:
AlphaMissense-pathogenic yet costing under 2 kcal/mol to fold. FoldX and AlphaMissense agree
overall (Spearman ρ = +0.61 [+0.56, +0.65], n = 751) but disagree systematically here.
**(D)** The FoldX disulfide term alone. It is 2.83 kcal/mol for cysteine-removing variants and
0.00 for both other classes (δ = +0.79, p = 3×10⁻¹¹²), confirming the cysteine cost is the lost
staple rather than general strain.

*Interpretation.* Calcium-site variants are pathogenic functionally, not thermodynamically. A
variant predictor that ranks FBN1 missense changes by predicted destabilisation alone will
systematically under-call one of the two principal disease mechanisms.

## Figure 3 — The calcium result is not a solvent-exposure artifact

**(A)** Relative solvent accessibility by class. cbEGF cysteines are buried (median 12.4%),
consistent with stapling the core, while other cbEGF residues sit at 45.4%
(δ = −0.71, q = 5×10⁻³⁶); calcium ligands are intermediate at 26.9%.
**(B)** Distance from each variant position to the nearest Ca²⁺ ion in the calcium-placed
structure (experimental where available, otherwise AF3-cofolded). The calcium class sits at
2.5 Å — inside the coordination shell — against ~11 Å for both other classes, confirming the
class is defined by measured geometry rather than by sequence position alone.
**(C)** FoldX ΔΔG for calcium-consensus positions versus other cbEGF residues, within burial
strata. Among buried residues (RSA < 20%) the calcium class costs 0.87 kcal/mol against 3.62
for its neighbours (δ = −0.56, q = 1.9×10⁻⁶); in the partial stratum the two are
indistinguishable (δ = −0.05, q = 0.58). No calcium-consensus variant is highly exposed
(RSA ≥ 50%), so that stratum contains only the comparison group — itself evidence against
exposure as the explanation. Panel C follows the Phase 5 test exactly and therefore uses the
annotation-level calcium class, unlike panels A and B.

## Figure 4 — Calcium, disulfides and pathogenic variants in a solved FBN1 fragment

All panels show PDB **2W86** (cbEGF9–hybrid 2–cbEGF10, X-ray 1.8 Å, crystallised with calcium),
renumbered into P35555 numbering via SIFTS; the two non-native cloning residues present in the
entry are removed rather than drawn as fibrillin. UniProt annotates residues 851–902 as "TB 4";
the literature name for the same region is hybrid domain 2.

**(A)** The fragment: two cbEGF domains (blue) flanking the hybrid domain (grey), each cbEGF
carrying one Ca²⁺ ion (green) and three disulfides in the invariant C1–C3, C2–C4, C5–C6 pattern
(sticks).
**(B)** The cbEGF9 calcium site. Three side-chain ligands — Asp807, Glu810 and Asn823, the
β-hydroxylation site — together with three backbone carbonyls (grey), dashed to the ion at
2.3–2.6 Å. The ion is drawn below its ionic radius so the ligands remain visible.
**(C)** The AF3 cofolded model (blue) superposed on the crystal structure (grey) over cbEGF9:
backbone RMSD 0.34 Å, and the cofolded ion sits 0.13 Å from the crystallographic one — drawn as
a small solid sphere inside a translucent shell because at this scale they coincide.
**(D)** The 30 ClinVar pathogenic variants at 20 distinct positions on this fragment, Cα shown
as a sphere: 8 cysteine-removing positions, 1 calcium-ligand position, and 11 positions outside
a cbEGF domain.

## Figure 5 — Calcium placement is a modelling assumption; this is its error bar

**(A)** Displacement between a placed Ca²⁺ ion and the experimentally observed one. AF3
cofolding reproduces X-ray references to a median of 0.32 Å (n = 20 domain comparisons) and NMR
references to 1.25 Å (n = 20) — the larger NMR figure tracking those references' own ~2 Å
coordinate spread. Homology transplant, evaluated leave-one-out across all donor→target cbEGF
pairs, is 4.6× worse (median 1.44 Å, n = 56). AF3 is therefore the primary method and transplant
is a cross-check.
**(B)** Why superposition frame matters. Fitting the whole two-domain construct gives a median
backbone RMSD of 2.30 Å against 1.26 Å per domain; the gap is interdomain hinge motion, not
placement error. An earlier whole-construct fit produced an apparent 4.46 Å worst case and would
have led us to discard a method accurate to 0.32 Å.
**(C)** CheckMyMetal parameters for all 16 submitted calcium sites — 8 from experimental
structures, submitted as controls, and 8 from AF3 models. gRMSD does not differ
(Mann-Whitney p = 0.25); non-significance is the desired outcome. Both groups fall outside
CheckMyMetal's published "good" bands, because cbEGF sites are irregular and carbonyl-dominated,
so the absolute thresholds do not transfer. Submitting the experimental controls is what
revealed this; no "N% of sites are outliers" statement is made anywhere in this work.

## Figure 6 — Enrichment at critical residues (supporting evidence, not the headline)

**(A)** Enrichment of the 578 pathogenic variants against the correct null — the share of
positions in each class, since every position offers the same 19 substitutions. cbEGF
cysteine-removing changes are enriched 5.9× [5.4, 6.3] (p = 8×10⁻¹⁵⁸); calcium-consensus
residues only 1.4× [1.1, 1.8] (p = 5×10⁻³). The asymmetry is consistent with Figure 2: calcium-
site variants are under-ascertained precisely because they look thermodynamically benign.
**(B)** Odds ratios on a log axis. The red intervals contrast ClinVar pathogenic against ClinVar
benign and are **inflated by construction** — ClinVar's pathogenic calls used ClinGen PM1, the
same critical-residue rule under test — so they are not an independent discovery. The blue
intervals repeat the contrast using AlphaMissense classes only, which never used FBN1-specific
rules, and reproduce the direction independently (calcium-consensus OR 8.9 [6.1, 12.9]).

Note the benign arm carries 34 variants, which is why the red intervals are wide; see
`results/LIMITATIONS.md` §2.

---

## Supplementary tables

`results/supplementary_tables.xlsx`, produced by `scripts/19_supplementary_tables.py`.

| Sheet | Contents |
|---|---|
| README | Provenance, the three governing caveats, and a sheet index |
| S1_variants | Master table, one row per ClinVar GRCh38 missense variant (4,451 × 47) |
| S2_structural | The 751 structurally covered variants with measured geometry (751 × 25) |
| S3_cmm_sites | 16 CheckMyMetal calcium sites: 8 experimental controls, 8 AF3 models |
| S4_statistics | 20 statistical tests with effect size, 95% CI and BH-corrected q |
| S5_calibration | AF3 model versus experimental reference, per cbEGF domain (40 rows) |
| S6_transplant | Leave-one-out Ca²⁺ transplant, donor → target (56 rows) |
| S7_dictionary | Definition of every column in S1 and S2 (53 rows) |
