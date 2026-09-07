# Limitations

Written to be read before any result in this project is believed. Ordered by how much each one
could change a conclusion.

> **Updated 2026-08-30 for the final analysis.** Sections 1, 3, 5 and 7 have been revised where
> the v2 work changed the facts, and sections 13–15 were added for the batch-4 mutant cofolding.
> Where a number below differs from an older report, this file is the current one.

---

## 1. The enrichment result is partly circular

**ClinVar's pathogenic classifications were made using ClinGen PM1 — the same critical-residue
rule this project tests.** Measuring enrichment of ClinVar-pathogenic variants at ClinGen-defined
critical residues therefore partly measures the classification rule rather than biology. The
headline figure from that comparison (91% of pathogenic vs 0% of benign at critical residues, ≥2
stars) is inflated and is **not** an independent discovery.

**What survives the caveat.** AlphaMissense, which never used FBN1-specific PM1 rules, reproduces
the direction independently: 33.1% of its predicted-pathogenic variants remove a cbEGF cysteine
versus 0.0% of predicted-benign. The biology is real; the magnitude from ClinVar labels is not
trustworthy.

**Consequence for the manuscript.** The biophysical measurements (ΔΔG in two force fields, Ca²⁺
coordination geometry, burial) are the headline, because they are computed from coordinates and
owe nothing to any classification scheme. Enrichment is Figure 6, last, with the caveat printed
inside the figure.

**Updated figure (2026-08-30).** The enrichment analysis now uses a 1,901-variant gnomAD
population denominator rather than the 34 ClinVar benign variants, which gave infinite odds
ratios. Side-chain Ca²⁺ ligands: OR 3.41 (2.21–5.27). cbEGF cysteines: OR 140 (93–205). The
backbone-only ligands are *depleted* (OR 0.049), which is the negative control the analysis
passes, and the motif non-ligands show no enrichment (1.49, 0.90–2.51).

## 2. The benign set is small, and this is not fixable

**34 benign missense variants at ≥2 stars; 76 at ≥1 star.** This is a property of the gene, not a
pipeline failure: of 2,783 gnomAD missense variants only 24 exceed the stand-alone benign
filtering-AF threshold, and all 24 were already captured from ClinVar. There are **zero**
gnomAD-only benign variants available to add. FBN1 is under strong selective constraint and simply
does not carry many common missense variants.

Any pathogenic-versus-benign contrast is therefore underpowered on the benign arm, especially once
split by site class. Comparisons against the all-possible-missense background (54,549
substitutions) are used where power matters, and benign comparisons are reported with confidence
intervals rather than bare p-values.

## 3. Calcium placement in predicted models is a modelling assumption

AlphaFold and AlphaFold-Server models contain no metals unless cofolded with them. Every calcium
ion in a predicted model here was **placed by AF3 cofolding**, not observed.

Calibration against experimental structures (Phase 4) gives the error bar: median **0.32 Å**
displacement against X-ray references, 1.25 Å against NMR references (whose own backbone spread is
~2 Å); pooled across all eight calibration domains the median is **0.45 Å**, against a median
backbone deviation of 1.35 Å for the same domains. Homology transplant, the alternative, is
roughly 4× worse (median 1.44 Å, 56 leave-one-out transplants). CheckMyMetal cannot distinguish
AF3 sites from experimental ones.

The ligand *identities* were checked independently of the geometry: for cbEGF9, where the 1.8 Å
crystal structure names its coordinating residues (Jensen 2009), the models return exactly the
same six residues with the same side-chain/backbone split.

This validates AF3 for this domain family. It does **not** make a predicted calcium an observed
one, and any structural metric derived from a predicted site inherits this uncertainty.

## 4. CheckMyMetal's absolute thresholds do not apply to cbEGF sites

Deposited, refined structures (2W86 at 1.80 Å, 1UZJ at 2.25 Å) return gRMSD ≈ 25–27° and valence
≈ 1.3–1.8 — values CMM's published bands call *outlier* and *below the formal Ca²⁺ valence of
2.0*. cbEGF calcium sites are intrinsically irregular: six or seven ligands dominated by backbone
carbonyls, in a low-symmetry arrangement that does not fit CMM's idealised geometric templates.

**Only the experimental-versus-predicted comparison is valid here.** No statement of the form
"N% of sites are CMM outliers" is made anywhere in this project, and none should be.

## 5. Structural coverage is partial — but much less so than in v1

Experimental structures cover 8 of the 43 cbEGF domains (cbEGF9, 10, 12, 13, 22, 23, 32, 33).
In the v1 analysis **751 of 4,451 missense variants (17%)** fell inside one of the four calibrated
constructs, and that set was not a random sample of the gene: it was the set of domains
crystallographers chose to solve.

**v2 removed that selection.** All 43 cbEGF domains were cofolded with calcium as 21 overlapping
tandem constructs, raising coverage to **3,831 of 4,451 (86%)**. What remains uncovered is the
non-cbEGF portion of the protein: the TB domains, the proline-rich region and the termini.

FoldX ΔΔG is still limited to the 751 on experimental templates — that is why the FoldX and
Rosetta comparisons in Figure 2 are computed on different subsets, and why the metal-aware,
better-powered Rosetta result is the one reported as primary.

## 6. Waters were removed before CheckMyMetal submission

AF3 models contain no water, so waters were stripped from the experimental structures too, for a
like-for-like comparison. Two experimental sites genuinely have a water ligand (2W86 site 2,
1UZJ site 2) and were therefore scored with one fewer ligand than as-deposited. 1LMJ and 1EMN are
NMR structures with no waters and are unaffected.

The effect makes the experimental *control* marginally worse, not the AF3 models better, so it
cannot inflate the calibration result. An as-deposited control would require re-submitting those
two files with waters retained.

## 7. ΔΔG is a prediction, not a measurement

FoldX ΔΔG is an empirical force-field estimate with typical reported error around 0.5 kcal/mol,
and it was computed on a single repaired structure per template rather than an ensemble. It
captures disulfide loss and packing well; it handles long-range electrostatics and large
conformational rearrangements poorly. Values are used comparatively (which variants are worse),
not as absolute thermodynamic quantities.

The same applies to the PyRosetta ΔΔG added in v2, with two protocol caveats specific to it.
First, a **steric artefact** was found and corrected: with the backbone held rigid, bulky
substitutions at buried cysteines produced unrelieved collisions rather than energies (21.5% of
values above 100 REU, worst 3,257). Allowing restrained backbone relaxation inside the repack
shell reduced that to 1.0% (worst 633) and changed the cysteine median from 24.76 to 10.81 REU.
The direction of every conclusion was unchanged; **the pre-correction magnitudes must not be
quoted**, and `tests/test_v2.py` A2 fails if the regression returns. Second, Rosetta energy units
are not kcal/mol; the two force fields are never pooled or plotted on a shared axis.

## 8. ClinVar label quality varies

Review status is captured per variant and a ≥2-star floor is applied for the primary analysis,
with a ≥1-star sensitivity tier. Even so, 5,437 of 9,327 GRCh38 records are single-submitter, and
Baudhuin 2019 found 31.6% of non-critical-residue FBN1 classifications overclassified once allele
frequency was considered. 704 records carry conflicting classifications; these are retained,
labelled, and excluded from the primary contrast rather than resolved by fiat.

## 9. Benign-set ascertainment bias

Benign variants reach ClinVar largely because someone sequenced a patient and needed the variant
adjudicated. They are not a random sample of tolerated variation, and they are enriched for
positions that look suspicious enough to be tested. This biases the benign set toward the same
domains as the pathogenic set, which if anything makes the enrichment contrast conservative.

## 10. Post-translational modifications are absent from all models

Fibrillin-1 is glycosylated, and the calcium-binding consensus includes a β-hydroxylation site
(the C3+2 Asn). Neither glycans nor β-hydroxylation are present in any structure or model used
here. β-hydroxylation is **not required** for calcium binding, so this does not invalidate the
coordination analysis, but predicted stabilities omit any contribution from these modifications.

## 11. Reference and version pinning

Results are tied to specific releases: ClinVar bulk **2026-07-28**, UniProt P35555 entry version
**264**, gnomAD **v4 (gnomad_r4)**, AlphaMissense Zenodo **8208688**. ClinVar in particular is
re-classified continuously — 29 FBN1 records already existed in the live index but not in the
snapshot at the time of analysis. Re-running later will not reproduce identical counts, which is
why every source is checksummed and versioned in `manifest/`.

## 12. Two literature items remain unresolved

- The VCEP's Supplementary Table 1 (the authoritative PM1 position list) was not available. My
  derivation reproduces **655** PM1 positions against the paper's stated **633** — a 3.4%
  discrepancy in the auxiliary categories (which glycines, which non-cbEGF cysteines). The two
  primary feature classes are exact: 258 cbEGF cysteines and 215 calcium-consensus positions.
  I did not tune parameters to force the totals to agree.
- The primary Jensen et al. (2009) paper for PDB 2W86 is absent; the file under that name is the
  2012 review by the same group, which supplied the PDB↔domain mapping but is a secondary source.

---

## 13. The disruption rate rests on 33 variants, and on a modelled ion

The headline rate — 45% of pathogenic side-chain Ca²⁺-ligand variants disrupt their site — comes
from 33 variants that passed the wild-type site-quality filter, out of 34 folded. The Wilson
interval is **30–62%**, which is wide, and the point estimate should be read with it. The control
arm is 24 variants and produced a single disruption, so the odds ratio (19.2) is even less
precisely determined than the rate.

The measured quantity is a **bond-valence sum over a modelled ion**, not a binding constant. It is
validated against CheckMyMetal on 16 sites (r = 0.961) but not against any measured affinity. The
single case where both exist — N2144 — agrees in direction and rank, which is the strongest
statement the data support.

AlphaFold3 is trained largely on holo structures and is biased toward placing an ion somewhere.
Every number here is therefore a **difference** against a wild-type model of the same construct,
folded with the same number of ions, and the decision threshold comes from controls folded under
the same bias. An absolute statement about any single site would not be defensible; a difference
against a matched wild type is.

## 14. Three variants were excluded, and the filter that excluded them is stated in advance

The construct spanning cbEGF1–3 folds inconsistently (global pLDDT 74; one site has a seed spread
of 39% of its own valence, and one site is occupied in only four of five wild-type models). Three
batch-4 variants fall on its sites and were dropped: D490Y, N261S and E487A.

The filter — seed spread ≤ 15% of valence, site occupied in all five wild-type models, local
pLDDT ≥ 70 — is computed **on the wild type alone, before any mutant is read**, so it cannot be
tuned on the answer. Since 2026-09-06 its outcome for every site is recorded rather than logged:
`scripts/64_site_qc.py` re-applies it to all 63 sites of the 21 wild-type constructs and writes
`data/processed/af3_site_qc.tsv`, which is sheet S10 of the supplementary workbook. **60 of 63
sites are admitted.** Two of the three failures are the cbEGF1–3 sites above; the third is
cbEGF5–7 site 1 (spread 16.2%), which carries no variant in this study. Only **one** admitted
site — cbEGF19–21 site 2 — sits within 20% of any limit, so the filter is not doing quiet work at
the margin. `tests/test_v2.py` A11 re-derives every admission call from the three criteria. Without it, one of those variants returned a +237% valence change that is
model noise: the ion had moved 25 Å between folds.

Related: whole-construct backbone RMSD reaches 29.8 Å in these tandem models because the
constructs hinge between domains from seed to seed. That is inter-domain motion, not a broken
fold, and it is why fold disruption and ion displacement are measured after superposing only the
residues within 12 Å of the ion. Anyone re-analysing multi-domain predicted structures should
expect this trap.

## 15. What the rate does not explain

Eighteen of 33 pathogenic Ca²⁺-ligand variants did **not** disrupt their site, and they are
overwhelmingly the substitutions that leave an oxygen donor behind (D→N, D→E, N→D, E→Q). Those
variants are classified pathogenic in ClinVar, and this analysis does not account for why. Three
possibilities it cannot distinguish: loss of β-hydroxylation at the consensus asparagine (not
modelled here); an affinity change too small to register as a change in modelled valence; or
overclassification in ClinVar, which Baudhuin 2019 documents for FBN1 at a rate of 31.6% outside
critical residues.

This is the clearest open question the project leaves, and it is testable: the prediction is that
D→N and D→E substitutions at these positions have milder calcium-binding defects than D→G or D→Y
at the same position.

## 16. The load-bearing assumption, stated plainly

Everything geometric in this project rests on one assumption that is not itself tested here:
**that an AlphaFold Server model of a mutant construct reports a real change in Ca²⁺ coordination
rather than a modelling artefact.** No affinity was measured in this work.

The experimental anchoring is thin, and it is worth naming exactly what it consists of.

| variant | what was measured experimentally | agreement with the models |
|---|---|---|
| N2144S | *K*d, 1.6 mM → ~14 mM (≈9-fold loss), fold unchanged | quantitative direction and rank; −28.5% valence |
| N548I  | proteolytic susceptibility of cbEGF4, local structural change | qualitative; −33.0%, no contact retained |
| E1073K | proteolytic susceptibility, same study | qualitative; −11.9% |

Three variants are a weak anchor for a class of 271. Until a direct binding measurement exists for
a modelled variant other than N2144S, the 45% is best described as a property of the models that
is consistent with the one measured case, not as an established property of fibrillin-1.

**The rate is a floor, and for a second reason beyond the donor chemistry.** A structure predictor
given a sequence that differs from the wild type at one position is drawn toward the wild-type
fold, so its error on a single substitution runs toward *no change*. That biases the disruption
call toward false negatives. Both the donor-chemistry effect and the predictor's conservatism push
the measured rate down and neither pushes it up, so the true proportion of pathogenic Ca²⁺-ligand
variants that weaken their site is **at least** 45%.

**Bond-valence sum is a geometric proxy.** It summarises how well the donor atoms surround the ion
and is not a binding free energy, so it supports statements of direction and rank and not
statements of affinity. This is stated in §13 and repeated here because it is the single most
common way a reader could over-read the result.
