<!--
  Keyed manuscript source. Citations are bracketed keys; the reference list is
  data/processed/references_final.tsv and numbering is assigned at build time by
  scripts/59_typeset_docx.py, in order of first appearance. Do not hand-number.
  Regenerate this file's citation scheme with scripts/58_manuscript_source.py.
-->

# Calcium-ligand and cysteine variants in fibrillin-1 damage the cbEGF module in different ways

**Draft, 2 September 2026.** Every number traces to a table in `data/processed/`, a row of
`results/statistics_final.tsv`, or a figure sidecar in `figures/*.numbers.json`.

---

## Abstract

Most Marfan syndrome missense variants fall on the calcium-binding EGF-like (cbEGF) module of
fibrillin-1, at one of two features: the six cysteines that form the module's three disulfide
bridges, or the small set of residues that coordinate a calcium ion. Variant-interpretation
guidelines treat the two as equivalent evidence. We tested whether they are equivalent in
mechanism.

We cofolded all 43 cbEGF domains with their calcium ions as 21 overlapping tandem constructs,
raising the fraction of 4,451 ClinVar missense variants with structural measurements from 17% to
86%, and defined coordinating residues from observed contact with the ion rather than from a
sequence motif. Removing a cysteine destabilised the domain in both an empirical and a
metal-aware energy function (median ΔΔG 10.8 versus 1.8 Rosetta energy units against a
composition-matched control, *P* < 10⁻⁴⁷). Removing a calcium ligand did not: it was
indistinguishable from D/E/N positions that coordinate nothing (*P* = 0.42), and within buried
positions it cost less, not more (*P* = 2 × 10⁻⁴).

We then cofolded 60 mutant constructs with calcium. Fifteen of 33 pathogenic calcium-ligand
variants (45%, 95% CI 30–62%) reduced the ion's coordination beyond a threshold measured from 123
sites where no effect was possible, against 1 of 24 controls (odds ratio 19, *P* = 7 × 10⁻⁴),
while the backbone within 12 Å of the ion moved a median of 0.29 Å. What decided the outcome was
donor chemistry, and where in the consensus the substitution fell: at the N-terminal Asp and Glu,
substitutions removing the coordinating oxygen disrupted the site in 9 of 13 cases and those
replacing it with another oxygen donor in 0 of 12 (*P* = 5 × 10⁻⁴). At the β-hydroxylation Asn the
rule fails — all three Asn→Ser substitutions we modelled lost coordination (median −17.2%), and
this is precisely the substitution the ClinGen FBN1 expert panel exempts from its critical-residue
criterion. Neither the stability calculation (Spearman ρ = −0.09) nor AlphaMissense (median
0.982 in both groups) separated the two outcomes.

Calcium-ligand and cysteine variants are therefore two mechanistically distinct classes. Stability
prediction detects the second and is blind to the first, and the property that distinguishes
damaging from tolerated calcium-site substitutions — whether an oxygen donor survives — is visible
to neither of the methods currently used as computational evidence.

---

## Introduction

Fibrillin-1 is the principal structural component of the 10–12 nm extracellular microfibrils that
give connective tissue its elasticity[@sakai1986], and it directs both their assembly and the
regulatory functions they carry in the matrix[@ramirez2009,asano2022]. Variants in *FBN1* cause
Marfan syndrome[@dietz1991], a dominantly inherited disorder defined clinically by the revised
Ghent nosology[@loeys2010] in which structural deficiency of the microfibril is compounded by
dysregulated TGF-β signalling[@neptune2003,milewicz2021]. The mature
protein is 2,871 residues built largely from one repeating module: the calcium-binding EGF-like
domain, of which there are 43 (refs. [@handford1995,jensen2012]). Each is a small β-hairpin unit closed by three
disulfide bridges in a fixed C1–C3, C2–C4, C5–C6 pattern, and each binds a calcium ion at the
junction with the preceding domain. Solution structures of consecutive pairs showed that the ion
rigidifies that junction into a rod-like arrangement[@downing1996,knott1996] and cryo-EM of intact microfibrils has
since placed those rods within the assembled fibril[@godwin2023].

Missense variants concentrate at two features of this module (**Fig. 1**), a pattern apparent
since the first curated *FBN1* mutation databases[@collodberoud2003] and one whose clinical
correlates depend in part on where in the protein the substitution falls[@delhomme2022]. Replacing one of the six
cysteines leaves an unpaired thiol and removes a bridge from the domain core. Replacing one of the
residues that coordinate the ion — in the consensus these are an aspartate, a glutamate and an
asparagine, joined by three backbone carbonyls — removes a contact with an ion that is not part of
the protein. The ClinGen *FBN1* variant curation expert panel[@rehm2015], applying the ACMG/AMP
classification framework[@richards2015], treats both as critical-residue (PM1) evidence,
though not at the same strength: cysteine-removing variants in any of the 43 cbEGF domains are
raised to PM1_Strong, while variants at the consensus calcium-binding residues
[D]-X-[D/N]-[E/H]-X*m*-[D/N]-X*n*-[Y/F] receive PM1 at moderate strength, with one exception —
Asn→Ser at the second D/N position, which the panel exempts on the basis of its frequency in
gnomAD[@drackley2024]. That ranking is justified in the specification by the prevalence of cysteine variants
among reported cases, not by a mechanistic comparison, and the calcium rule is applied to a
sequence consensus rather than to observed coordination.

There is reason to ask whether that structure survives a mechanistic test, and the reason is old. In coagulation factor IX, whose
first EGF-like domain carries the same consensus, Handford *et al.* showed in 1991 that removing a
coordinating aspartate reduced calcium affinity more than a thousandfold and caused haemophilia
B[@handford1991]. In fibrillin-1 itself, Kettle *et al.* characterised N2144S — a substitution that removes a
coordinating asparagine of cbEGF32 — and found that it does not alter the native fold of cbEGF32
or its neighbours while reducing calcium affinity roughly ninefold, from a *K*d of 1.6 mM to about
14 mM[@kettle1999]. Related work found that calcium-ligand substitutions increase susceptibility to
proteolysis[@reinhardt2000] that their molecular consequences depend on the domain preceding the
site[@mcgettrick2000] and that some cbEGF substitutions impair secretion rather than folding[@whiteman2003,whiteman2007]. Prior
simulation has examined the reverse direction — how cysteine mutations perturb calcium binding
under mechanical load[@haller2020] — but not whether calcium-ligand substitutions destabilise the domain.

Whether N2144S is representative of its class has not been established, largely because so few
fibrillin-1 domains have experimental structures: eight of 43, from four crystallographic and NMR
fragments[@rao1995,lee2004,smallridge2003,jensen2009]. We addressed this by generating calcium-bound models of all 43 cbEGF domains and
of 60 mutants, defining coordinating residues from observed contact with the ion, and comparing
the energetic and geometric consequences of the two variant classes on a common footing.

---

## Results

### Variants concentrate at the two features of the module

We assembled 4,451 *FBN1* missense variants from ClinVar[@landrum2018] normalised to the MANE Select transcript
NM_000138.5 and UniProt P35555, comprising 578 pathogenic or likely pathogenic, 823 of uncertain
significance and 34 benign or likely benign at a two-star review floor. Pathogenic variants fall
across the cbEGF array while the gnomAD[@karczewski2020,chen2024] population comparison group is spread more evenly
(**Fig. 1a**).

Aligning every cbEGF domain on its own N terminus makes the concentration explicit (**Fig. 1b**).
Of the 425 pathogenic variants inside a cbEGF domain, the tallest peaks sit exactly on the six
cysteines and on the calcium donors. Those positions are not assumed: across all 43 domains the
side-chain donors are at offsets 0 (Asp, 43/43 domains), +3 (Glu, 34/43) and +17/+18 (Asn), the
backbone carbonyls at +1 and +19 to +22, and the cysteines at median offsets +4, +11, +16, +25,
+27 and +40. The site itself is visible in the 1.8 Å crystal structure of the
cbEGF9–hybrid 2–cbEGF10 fragment[@jensen2009] (**Fig. 1c**), where the ion is held by Asp807, Glu810 and
Asn823 through their side chains and by the carbonyls of Ile808, Ser824 and Ser827.

Experimental structures cover only eight of the 43 domains, and the variants inside them are not a
representative sample — they are whichever positions were solved decades ago. We therefore folded
the full array as 21 overlapping two- and three-domain constructs with AlphaFold Server[@abramson2024] each
cofolded with three Ca²⁺ ions and five models per construct. This raised structural coverage from
751 variants (17%) to 3,831 (86%). Coordinating residues were then defined by observed contact
with the modelled ion. That definition matters: the conventional sequence motif over-calls,
flagging 143 positions that coordinate nothing in any model, and these are kept throughout as a
separate group that tests whether the structural definition is doing work.

Four groups are compared below: variants at side-chain calcium ligands (n = 271); variants at
aspartate, glutamate and asparagine positions that coordinate nothing (n = 323, the
composition-matched control, which exists because calcium ligands are D/E/N by definition);
motif positions that are not ligands (n = 143); and variants removing a cbEGF cysteine
(n = 1,132).

### Cysteine loss destabilises the domain; calcium-ligand loss does not

We computed folding stability changes with two unrelated energy functions: FoldX[@schymkowitz2005] an empirical
force field, and PyRosetta's ref2015[@alford2017] which forms explicit bonds between the ion and its
coordinating atoms and is therefore the calculation that could show a penalty at these positions.

Neither does (**Fig. 2a,b**). In the metal-aware calculation cysteine-removing variants cost a
median 10.81 Rosetta energy units against 1.77 for the composition-matched control
(*P* = 3 × 10⁻⁴⁷), while calcium ligands cost 1.84 — indistinguishable from the same control
(*P* = 0.42). FoldX gives the same ordering, with a small residual difference at calcium ligands
(*P* = 0.03) computed on a fifth as many variants, because it is limited to the domains with
experimental templates.

Calcium ligands are more buried than the control (median relative solvent accessibility 21.7%
versus 38.4%), so the comparison could in principle reflect location rather than role. It does
not (**Fig. 2c**). Among buried positions calcium ligands cost *less* than buried controls, 1.72
against 4.14 Rosetta energy units (*P* = 2 × 10⁻⁴); among partly buried positions the two are
level (1.92 against 2.01, *P* = 0.54). Only one of the 271 calcium ligands is substantially
exposed, so that band cannot be assessed.

The contrast with sequence-based prediction is stark (**Fig. 2d**). AlphaMissense[@cheng2023] which never
sees a structure, scores calcium-ligand variants at a median 0.989 out of 1 — as damaging as
cysteine variants (0.998) and far above the composition control (0.260). These variants are
therefore recognised as damaging by an orthogonal method while being, by direct calculation,
energetically unremarkable.

### Nearly half of pathogenic calcium-ligand variants loosen the ion

If the fold survives, the consequence must lie in the site. We cofolded 60 further constructs with
calcium — 34 pathogenic side-chain ligand variants, 12 composition controls and 14 negative
controls — five models each, against the wild-type models of the same constructs.

Each construct carries three calcium sites and a substitution can reach only one. We measured the
*cognate* site, defined as the site whose wild-type coordination shell contains the mutated
residue, and used the other two as within-construct controls: same fold, same five seeds, a site
the mutation cannot touch. Those 109 bystander measurements, pooled with the 14 negative-control
sites, give 123 measurements of a change where none is possible. The 95th percentile of their
absolute change is 9.1% of site valence, and that is the threshold used below — measured, not
chosen. Three variants fell on sites of one poorly converged construct (cbEGF1–3, pLDDT 74) that
failed a quality filter fixed on wild-type evidence before any mutant was read, leaving 57.

Fifteen of 33 pathogenic calcium-ligand variants (45%, 95% CI 30–62%) lost coordination beyond the
threshold, against 1 of 24 controls (odds ratio 19, Fisher's exact *P* = 7 × 10⁻⁴)
(**Fig. 3a,b**). Median change was −5.7% for pathogenic ligands against −1.7% and +0.1% for the
two control groups. The largest loss, −38.0%, is N2144H at the position whose serine substitution
Kettle *et al.* measured; N2144S itself, folded in a pilot series, loses 28.5% (**Supplementary
Fig. S1c**).

Meanwhile the fold does not move (**Fig. 3d**). Backbone r.m.s.d. over residues within 12 Å of the
same ion is a median 0.29 Å for pathogenic ligands against 0.24 Å for controls, with a maximum of
0.80 Å across all 57 variants — a Cα displacement smaller than a bond length. Whole-construct
r.m.s.d. is far larger and far more variable, up to 29.8 Å, because these tandem constructs hinge
between domains from seed to seed; that is why the fold is judged locally, and it is a trap for
anyone measuring domain displacement on multi-domain predicted structures.

### Donor chemistry decides the outcome

The rate is not 100%, and the exceptions are systematic. Substitutions that remove the
coordinating oxygen entirely disrupted the site in 12 of 16 cases (75%); substitutions that
replace it with another oxygen donor — D→N, D→E, N→D, E→Q — did so in 3 of 17 (18%), odds ratio
14, *P* = 0.002. In 18 of 33 pathogenic ligand variants the substituted side chain still contacted
the ion in all five models. The ion is coordinated by an atom, not by a residue label, and a
conservative substitution that preserves the donor largely preserves the site — with one
systematic exception, described next.

Two of these variants have published experimental characterisation, and both agree. N548I, which
Reinhardt *et al.* showed renders cbEGF4 susceptible to proteolysis with local structural change
close to the mutation[@reinhardt2000] loses 33.0% of its coordination here, retains no contact with the ion in
any model, and moves the local backbone by 0.47 Å. E1073K, characterised in the same work, loses
11.9% at its cognate site — above the 9.1% threshold, although it did not clear the per-site
noise criterion used in the smaller pilot series, which is one reason the pooled null replaced it.

### The donor rule holds at the Asp and Glu positions and fails at the β-hydroxylation Asn

The consensus site has two kinds of position, and they do not behave alike (**Fig. 3c**). At the
aspartate and glutamate at the domain's N terminus (offsets 0 and +3), the donor rule is clean:
none of 12 donor-preserving substitutions disrupted the site (median −1.0%) and 9 of 13
donor-removing substitutions did (median −12.2%; Fisher's exact *P* = 5 × 10⁻⁴). At the
asparagine at offset +16 to +20 — the β-hydroxylation site[@rees1988] — it fails. Three of five
donor-preserving substitutions disrupted the site, at a median of −14.6%, and all three of the
Asn→Ser substitutions we modelled (N589S, N1046S, N2624S) lost coordination, median −17.2%. A
serine hydroxyl is a shorter and weaker donor than an asparagine carboxamide, and at this position
it does not reach.

Asn→Ser at this position is the one substitution the ClinGen *FBN1* panel exempts from PM1[@drackley2024]. It is
also the most common substitution at that position in our cohort (18 of 66 variants), and four of
those 18 are two-star pathogenic or likely pathogenic in ClinVar — including N2144S, the single
*FBN1* variant whose calcium affinity has been measured, which falls roughly ninefold[@kettle1999]. Our
models place N2144S at −28.5% (**Supplementary Fig. S1c**) and N2144H at −38.0%. With three
modelled Asn→Ser variants the comparison against all other ligand substitutions is not itself
significant (*P* = 0.13), so this is an indication rather than a demonstration; but the direction
is consistent across every line of evidence available, and the exemption rests on a
population-frequency argument that these variants — all singletons or absent in gnomAD — do not
obviously support.

### Neither instrument in routine use sees the difference

Across the same 33 variants, Rosetta ΔΔG is uncorrelated with the change in coordination
(Spearman ρ = −0.09, *P* = 0.62) (**Fig. 4a**), and AlphaMissense assigns a median score of 0.982
to the variants that disrupt the site and 0.982 to those that do not (*P* = 0.61) (**Fig. 4b**). A
stability calculation cannot see the damage; a sequence-based predictor sees damage everywhere and
cannot resolve which variants cause it. The *FBN1* panel specifies REVEL[@ioannidis2016] rather
than AlphaMissense for PP3 (applied at REVEL ≥ 0.75)[@drackley2024] and we did not test REVEL; the saturation we observe is a
property of the predictor we measured, and whether REVEL separates these outcomes is an open and
easily answered question.

Pathogenic variants are enriched at both features relative to the population comparison group
(**Fig. 4c**): side-chain calcium ligands, odds ratio 3.41 (2.21–5.27); cysteine positions, 140
(93–205). Two features of this analysis matter more than the headline numbers. It is partly
circular, because the guidelines used to assign the clinical labels apply a rule at precisely
these positions[@drackley2024]. And the backbone-only calcium ligands are strongly *depleted* among pathogenic
variants (odds ratio 0.049, 0.017–0.212) — the expected result, since no amino-acid substitution
can remove a backbone carbonyl, and the negative control that the analysis passes. Motif positions
that are not ligands show no significant enrichment (1.49, 0.90–2.51), the corresponding check on
the structural definition.

### The modelled calcium sites reproduce experimental ones

Every geometric statement above depends on the models placing calcium correctly. Across the eight
domains with a calcium-bound experimental reference, the modelled ion sat a median 0.45 Å from its
experimental position against a median backbone deviation of 1.35 Å for the same domains — placed
more accurately than the protein around it, in all eight (**Supplementary Fig. S1a**). Summarising
each domain by the median over its five models instead gives 0.55 Å against 1.35 Å and seven of
eight, the exception being cbEGF12, where the two are equal to within 0.02 Å.

Ligand identity agrees independently of geometry. Jensen *et al.* report the cbEGF9 site as three
side-chain oxygens (Asp807, Glu810, Asn823), three main-chain carbonyls (Ile808, Ser824, Ser827)
and a water[@jensen2009]; our models, which never saw that assignment, return exactly those six residues with
the same side-chain/backbone split. We also verified the coordination measure — a Brown–Altermatt
bond-valence sum[@brown1985] — against CheckMyMetal[@zheng2017] on the 16 sites that service has scored: *r* = 0.961
with a constant offset of +0.147 (**Supplementary Fig. S1b**). A constant offset cancels when a
mutant is subtracted from its wild type, which is the only operation performed here.

---

## Discussion

The two classes of fibrillin-1 variant that carry critical-residue evidence are not equivalent in
mechanism, and the difference is not the one the evidence weighting encodes. Removing a cysteine
destabilises the domain, and any stability-based method detects it. Removing a calcium-coordinating side chain does not destabilise the domain — by
direct calculation in two independent energy functions, with composition and burial controlled —
yet in nearly half of cases it measurably loosens the ion while the backbone around the site stays
within half a bond length of where it started. That is the behaviour Kettle *et al.* measured for
N2144S at the bench in 1999[@kettle1999] and the behaviour Handford *et al.* inferred for the equivalent
substitution in factor IX eight years earlier[@handford1991]; what is new here is that it is the majority
behaviour of the class, and that its exceptions have a chemical explanation.

The 45% figure is a floor, and the shortfall is the informative part. The variants that spare the
site are overwhelmingly those that leave an oxygen donor behind, and a domain that keeps its
ligand keeps most of its coordination. Those variants are nonetheless classified pathogenic, and
this analysis does not account for why. Three possibilities it cannot separate: loss of
β-hydroxylation at the consensus asparagine, which is not modelled here; an affinity change too
small to register as a change in modelled valence; or overclassification, which Baudhuin *et al.*
document for *FBN1* at 31.6% outside critical residues[@baudhuin2019]. The distinction is testable, and the
prediction is specific: among pathogenic calcium-ligand variants, D→N and D→E substitutions should
show milder calcium-binding defects than D→G or D→Y at the same position.

Two further mechanisms sit outside what we measured and would compound the effect we did. Calcium
occupancy protects cbEGF domains from proteolysis, and calcium-ligand substitutions increase
proteolytic susceptibility[@reinhardt2000]; a site that binds its ion more weakly spends more time unprotected.
And the site is mechanosensitive — steered molecular dynamics shows cbEGF calcium binding weakens
under strain, with disulfide mutations altering that response[@haller2020] — so a weakened site in a
load-bearing microfibril may fail at forces a wild-type site tolerates. Our measurements are made
at zero load on a static model and cannot address either; they establish the resting-state defect
that both would amplify.

Two of our results bear directly on the *FBN1* PM1 rule as written. First, that rule identifies
calcium residues from the sequence consensus, and the consensus over-calls: 143 positions it flags
coordinate no calcium in any of our models, and those positions show no enrichment among
pathogenic variants (odds ratio 1.49, 0.90–2.51) while true side-chain ligands show 3.41
(2.21–5.27). Defining the ligand set structurally would concentrate the same evidence on fewer
positions. Second, the published exemption for Asn→Ser at the β-hydroxylation asparagine is the
one place our donor rule breaks: every Asn→Ser variant we modelled there lost coordination, the
position's best-characterised variant loses ninefold affinity at the bench[@kettle1999] and four of the
eighteen such variants in our cohort are two-star pathogenic. We are not in a position to overturn
a rule built on curated case data with three modelled variants, and we do not claim to; we are in
a position to say that the exemption deserves re-examination with the evidence that now exists.

There is a practical consequence for variant interpretation more broadly. Computational evidence now enters
classification at calibrated strengths[@pejaver2022] and stability prediction is a reasonable proxy for one
class of damage. Applied to calcium-site variants it returns values indistinguishable from benign
substitutions, and applied alone it would misclassify them. Where a variant sits at a
metal-coordinating position the relevant question is not whether the protein still folds but
whether the site still binds — and our results indicate the answer turns on whether the
substitution preserves an oxygen donor, a property neither a stability calculation nor a
sequence-based score currently resolves. A donor-aware rule is simple enough to state and cheap
enough to apply, and our data suggest it should be position-aware too: conservative at the
N-terminal Asp and Glu, not at the β-hydroxylation Asn.

Finally, the cbEGF module is not unique to fibrillin-1. The same fold and the same consensus carry
disease variants in coagulation factor IX[@handford1991] and in the fibrillin and LTBP superfamily more
broadly[@jensen2012]. If calcium-site substitutions behave this way here, the same logic should apply wherever
a small domain uses a bound ion to rigidify an interface, and the same blind spot in
stability-based interpretation would follow.

---

## Limitations

**Structural coverage is modelled.** Eight of 43 cbEGF domains have a calcium-bound experimental
structure; the rest are predictions. They reproduce the solved cases closely, including the ion's
position and the identity of its ligands, but that is a reason for confidence, not a substitute
for measurement.

**AlphaFold3 is trained largely on holo structures** and is biased toward placing an ion
somewhere. Every measurement here is therefore a difference against a wild-type model of the same
construct folded with the same number of ions, and the threshold comes from controls folded under
the same bias.

**Bond-valence sum is a geometric proxy for affinity, not a measurement of it.** It is validated
against CheckMyMetal but not against a binding constant. The cases where both exist — N2144, and
N548I and E1073K through proteolytic susceptibility — agree in direction and rank, which is the
strongest statement the data support.

**The rate rests on 33 variants** with a 95% interval of 30–62%, and the control arm of 24
produced a single disruption, so the odds ratio is less precisely determined than the rate.

**The position-stratified result rests on small strata.** The β-hydroxylation Asn stratum contains
eight modelled variants and the Asn→Ser subset three. The Asp/Glu stratum (25 variants) is what
carries statistical weight; the Asn stratum is a signal to be tested, not a result to be acted on
alone.

**The benign comparison group is small**: 34 variants, 29 with structural coverage, none at a
calcium ligand or a cysteine. This is why a population comparison group was substituted for the
enrichment analysis.

**Enrichment is partly circular**, as described above.

**The two energy functions differ slightly at calcium ligands** (*P* = 0.03 for FoldX, 0.42 for
Rosetta), computed on 50 and 271 variants respectively. We regard the better-powered, metal-aware
result as the more reliable and report both.

**A steric artefact was identified and corrected.** With the backbone held rigid, large residues
replacing buried cysteines produced unrelieved collisions rather than energies: 21.5% of values
exceeded 100 Rosetta energy units. Allowing local backbone relaxation under positional restraints
reduced this to 1.0% and changed the cysteine median from 24.76 to 10.81. Directions were
unaffected; the earlier magnitudes should not be quoted.

**Post-translational modification is absent from the models.** cbEGF domains carry
β-hydroxylation at the consensus asparagine and O-linked glycosylation; neither is required for
calcium binding, and neither is modelled.

---

## Methods

**Variants.** *FBN1* (Gene ID 2200) missense variants were retrieved from ClinVar[@landrum2018] (release
2026-07-28) via NCBI E-utilities and normalised to the MANE Select transcript
NM_000138.5[@morales2022] and UniProt P35555 (entry version
264; sequence MD5 f8f84ac9dc30f3ff640ad4f691fee2b1). HGVS descriptions[@dendunnen2016] were
validated with Mutalyzer[@wildeman2008], and the wild-type residue was checked against the
reference sequence at every position;
failures are listed in `data/interim/wt_mismatch.tsv` and were excluded. Clinical classes were
assigned from ClinVar review status with a two-star floor for the primary analysis and a one-star
sensitivity set; conflicting classifications are labelled, never dropped. Population frequencies
are from gnomAD v4[@chen2024] and pathogenicity scores from AlphaMissense v1.0.0[@cheng2023] joined with the same
wild-type identity check.

**Structures.** Experimental structures from the PDB[@vallat2026] were used where they exist, as
gold standard and as
calibration targets: 2W86 (cbEGF9–hyb2–cbEGF10)[@jensen2009] 1LMJ (cbEGF12–13)[@smallridge2003] 1UZJ
(cbEGF22–TB4–cbEGF23)[@lee2004] and 1EMN/1EMO (cbEGF32–33)[@downing1996]. Residue numbering was mapped through SIFTS[@dana2019]
and verified by residue identity, not index. All 43 cbEGF domains were additionally modelled as 21
overlapping two- and three-domain constructs on AlphaFold Server[@abramson2024] each cofolded with three Ca²⁺
ions, five models per construct; construct sequences were verified against P35555 over the stated
range before any measurement. Sixty mutant constructs were folded the same way, each verified to
differ from its wild type at exactly one position, at the position claimed, with the residues
claimed. Cofolding was chosen over the two alternatives for placing the ion. The AlphaFold
Protein Structure Database[@varadi2024] holds no model for P35555, whose length exceeds its cap,
and a monomeric AlphaFold2 model[@jumper2021] would in any case contain no metal. Transplanting
the ion by homology[@hekkelman2023] was evaluated directly: over 56 leave-one-out transplants
between the calibration domains it placed calcium a median 1.44 Å from the crystallographic
position, against 0.45 Å for cofolding.

**Site definitions.** A residue is a calcium ligand if one of its oxygen, nitrogen or sulfur atoms
lies within 3.2 Å of a modelled Ca²⁺, separating side-chain from backbone-only donors. Cysteine
roles and disulfide partners were taken from UniProt[@uniprot2023] and confirmed against the models; cbEGF
domains pair C1–C3, C2–C4 and C5–C6 without exception across all 129 cbEGF disulfides. UniProt
annotates no calcium site in P35555, so there is no annotation against which to validate the site
definition; the sequence-motif definition was retained as a separate comparison group instead.

**Energetics.** FoldX 5.1 ΔΔG[@schymkowitz2005] was computed with RepairPDB followed by BuildModel, mutation
strings written in author numbering via SIFTS. PyRosetta ref2015 ΔΔG[@alford2017] was computed with explicit
metal bonds (`auto_setup_all_metal_bonds`; metal *constraints* were not used, because the packer
optimises blind to them and the rescore then sees the violation), minimisation in torsion space,
and backbone relaxation inside an 8 Å repack shell under harmonic Cα restraints (σ = 0.5 Å,
after the restrained-relaxation protocol of Nivón *et al.*[@nivon2013])
applied during minimisation and removed before scoring. The two are different quantities in
different units and are never pooled.

**Geometry.** Coordination was quantified as a Brown–Altermatt bond-valence sum[@brown1985] over donor atoms
within 3.2 Å (R₀ = 1.967 Å for Ca–O, 2.14 Å for Ca–N, b = 0.37), validated against
CheckMyMetal[@zheng2017]. The 3.2 Å cutoff follows the surveyed distance distributions for
calcium sites in protein structures[@harding2006].
Sites were matched between wild type and mutant by coordinating-residue overlap (Jaccard ≥ 0.35,
assigned by the Hungarian algorithm) rather than by proximity, so a vacated site reports as
vacated rather than being matched to a spurious partner. Every quantity is a median over five
models. Ion displacement and backbone r.m.s.d. were computed after superposing the site's own
neighbourhood — residues with any atom within 12 Å of the wild-type ion — because these tandem
constructs hinge between domains and a global superposition converts a hinge into an apparent ion
displacement. Sites were admitted on wild-type evidence alone, before any mutant was read: seed
spread ≤ 15% of the site's own valence, site occupied in all five wild-type models, and mean
pLDDT ≥ 70 over the coordinating residues. Solvent accessibility was computed with FreeSASA[@mitternacht2016] and
secondary structure with DSSP[@kabsch1983] as implemented in mkdssp[@touw2015].

**Statistics.** Group comparisons use Mann–Whitney *U* tests[@mann1947] with Cliff's
δ[@cliff1993] and 2,000-resample
bootstrap confidence intervals, additionally reported stratified by wild-type residue where a
composition effect could confound. Rates use Fisher's exact test with Wilson score
intervals[@wilson1927].
Benjamini–Hochberg correction[@benjamini1995] is applied within each family of tests; the figures show
uncorrected *P* values and `results/statistics_final.tsv` carries both. All analyses are seeded
and reproduce identically on re-run, enforced by an automated gate.

**Software.** Python 3.14.4 (analysis) and 3.11.15 (PyRosetta); gemmi[@wojdyr2022], FreeSASA
2.2.1, mkdssp 4.5.7, PyMOL 3.1.0, NumPy 2.5.1[@harris2020], SciPy 1.18.0[@virtanen2020],
pandas 3.0.5, Matplotlib 3.11.1[@hunter2007]. Versions and paths
are in `manifest/tools.md`; every download is logged in `manifest/` with its date, query and row
count.

---

## Data and code availability

Processed tables are in `data/processed/`, statistical output in `results/statistics_final.tsv`,
and the supplementary workbook in `results/supplementary_tables_v2.xlsx`. Each figure is generated
by exactly one script in `scripts/` and writes a `.numbers.json` sidecar listing every value it
prints, which `tests/test_v2.py` re-derives from the processed tables. AlphaFold Server outputs,
including the server-assigned seeds, are in `structures/alphafold/`; their provenance, with a
SHA-256 over all returned coordinate files, is in `manifest/`.

---

## Figure legends

**Fig. 1 | Fibrillin-1 variants concentrate on the two features of the cbEGF module.**
**a**, ClinVar pathogenic and likely pathogenic variants (above, black) and the gnomAD v4 population
comparison group (below, grey) along fibrillin-1, as counts per residue, either side of the domain
architecture. **b**, Every cbEGF domain aligned on its own N terminus. Bars are pathogenic variant
counts per position (425 variants in cbEGF domains); the schematic above shows the consensus
module, with the six cysteines (circles, C1–C6) and their three disulfides drawn as arcs, and the
calcium donors (diamonds; filled, side-chain; open, backbone carbonyl) at the offsets measured
across all 43 domains. Bars are shaded to match: black at cysteine positions, mid grey at
side-chain donors, light grey elsewhere. **c**,
The cbEGF9 calcium site in the 1.8 Å crystal structure of cbEGF9–hybrid 2–cbEGF10 (PDB 2W86),
with the three side-chain and three backbone donors labelled.

**Fig. 2 | Cysteine loss destabilises the domain; calcium-ligand loss does not.**
**a**, ΔΔG from PyRosetta ref2015, which forms explicit bonds between the ion and its coordinating
atoms, for the four site classes. **b**, ΔΔG from FoldX, restricted to the domains with
experimental templates. **c**, The Rosetta ΔΔG of **a** stratified by relative solvent
accessibility, comparing calcium ligands with the composition-matched control within each band.
**d**, AlphaMissense score for the same four groups; the dashed line is the developers' pathogenic
threshold. Points are individual variants; bars are medians and boxes the interquartile range,
with *n* in each axis label and off-scale points counted at the axis. Groups are separated by
fill — solid, grey, open — rather than by colour; every figure in this paper is monochrome, so
none of them depends on colour reproduction. *P* values are two-sided
Mann–Whitney; effect sizes and corrected *q* values are in `results/statistics_final.tsv`. FoldX
(kcal mol⁻¹) and Rosetta (Rosetta energy units) are different quantities and are never plotted on
a shared axis.

**Fig. 3 | Nearly half of pathogenic calcium-ligand variants loosen the ion while the fold stays
put.** **a**, Change in bond-valence sum at each variant's cognate calcium site, mutant minus wild
type, ranked, for 57 cofolded variants. The shaded band is ±9.1%, the 95th percentile of absolute
change across 123 sites where no effect is possible (109 bystander sites in the same constructs
and 14 negative-control sites). **b**, Fraction of each group whose cognate site fell below the
threshold, with Wilson 95% intervals; odds ratio and Fisher's exact *P* compare pathogenic ligands
with both control groups pooled. **c**, The same change for pathogenic calcium-ligand variants
only, split by whether the substituted residue can still donate an oxygen and by which consensus
position is hit: the N-terminal Asp and Glu (offsets 0 and +3) or the β-hydroxylation Asn (offsets
+16 to +20). Open rings mark the three Asn→Ser substitutions that the ClinGen *FBN1* panel exempts
from PM1. Counts below each group are the number disrupted. **d**, Backbone
r.m.s.d. over residues within 12 Å of the same ion, mutant against wild type, on the same models.

**Fig. 4 | The damage is invisible to the methods used to detect it.**
**a**, Rosetta ΔΔG against the change in calcium valence for the 33 pathogenic calcium-ligand
variants; filled points are the variants whose site was disrupted. **b**, AlphaMissense score for
the same variants, split by whether the site was disrupted. **c**, Odds ratios with 95% intervals
for pathogenic variants at each feature against the 1,901-variant gnomAD population comparison
group, on a logarithmic scale.

**Supplementary Fig. S1 | The modelled calcium sites are trustworthy.**
**a**, Deviation of the modelled ion from its experimental position against the backbone deviation
of the same domain, for the eight domains with a calcium-bound structure; the shaded half-plane is
where the ion is placed more accurately than the fold. **b**, Local bond-valence sum against
CheckMyMetal valence over the 16 sites that service scored. **c**, The nine-variant pilot series,
scored against per-site wild-type seed spread; filled points lie beyond twice that spread.

**Supplementary Fig. S2 | Both mechanisms on experimental coordinates.**
**a**, The cbEGF9–hybrid 2–cbEGF10 fragment (PDB 2W86, 1.8 Å, crystallised in 20 mM CaCl₂): two
cbEGF domains either side of the hybrid domain, the six cbEGF disulfides as sticks, two Ca²⁺ ions.
**b**, Pathogenic ClinVar variants on the same fragment, Cα drawn as a sphere: black at the
calcium ligand, as the ion itself is drawn; mid grey at cysteines; light grey elsewhere.

---
