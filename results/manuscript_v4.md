# Distinct mechanisms of calcium-ligand and cysteine substitutions in the fibrillin-1 cbEGF module

[Author name]¹

¹ [Department, Institution, City, Postcode, Country]

**Corresponding author.** [Author name], [email address]

**Draft, 6 September 2026.** Citations have been removed from this version and will be inserted
from Mendeley. The keyed source carrying all 63 verified references remains at
`results/manuscript_paper.md`.

---

## Abstract

Most missense variants that cause Marfan syndrome (MFS) fall at one of two features of the
calcium-binding epidermal growth factor-like (cbEGF) module of fibrillin-1, either the six
cysteines that form its three disulfide bridges or the small group of residues that bind the
calcium ion. Current interpretation guidelines treat the two features as equivalent evidence of
pathogenicity. We asked whether they are equivalent in mechanism. All 43 cbEGF domains were
modeled with their calcium ions as 21 overlapping tandem constructs, together with 60 mutants, the
ligands being identified from observed contact with the ion instead of from the sequence
consensus. Substitution of a cysteine destabilized the domain in two independent energy functions,
whereas substitution of a calcium ligand did not, its computed folding penalty being
indistinguishable from that of aspartate, glutamate and asparagine positions that bind no ion. In
the mutant models 45% of pathogenic calcium-ligand variants lost calcium coordination (95% CI 30
to 62%), against 1 of 24 controls, while the backbone around the ion moved by a median of 0.29 Å.
Substitutions that removed the contacting oxygen atom broke the site, whereas those supplying
another oxygen atom largely spared it. The exception was the β-hydroxylated asparagine, at which
every Asn→Ser substitution modeled lost coordination. That is the substitution the ClinGen *FBN1*
expert panel exempts from its critical-residue criterion, its best-characterized example losing
calcium affinity approximately ninefold at the bench. Neither the stability calculation nor
AlphaMissense separated the variants that lost coordination from those that kept it.
Calcium-ligand and cysteine substitutions are therefore two mechanistically distinct classes of
*FBN1* missense variant. The property that separates a damaging calcium-site substitution from a
tolerated one is resolved by neither computational method now used as evidence.

---

## Introduction

Fibrillin-1 is the principal structural component of the 10–12 nm extracellular microfibrils
responsible for the elasticity of connective tissue, directing both their assembly and the
regulatory functions they carry in the matrix. Variants in *FBN1* cause Marfan syndrome, a
dominantly inherited disorder defined clinically by the revised Ghent nosology. In MFS the
microfibril is structurally deficient and TGF-β signaling is dysregulated. The mature protein is
2,871 residues and is built largely from one repeating module, the calcium-binding EGF-like
(cbEGF) domain, of which there are 43. Each cbEGF domain is a small β-hairpin unit closed by three
disulfide bridges in a fixed pattern (C1–C3, C2–C4 and C5–C6), binding a calcium ion at the
junction with the preceding domain. Solution structures of consecutive pairs showed the ion
rigidifying this junction into a rod-like arrangement. Cryo-EM of intact microfibrils has since
placed these rods within the assembled fibril.

Missense variants concentrate at two features of this module (**Fig. 1**), a pattern apparent
since the first curated *FBN1* mutation databases, with clinical correlates depending in part on
the position of the substitution within the protein. Substitution of one of the six cysteines
leaves an unpaired thiol and removes a bridge from the domain core. Substitution of a residue that
binds the ion removes a contact with an ion that is not part of the protein. In the consensus
these residues are an aspartate, a glutamate and an asparagine, together with three main-chain
carbonyls. The ClinGen *FBN1* variant curation expert panel applies the ACMG/AMP classification
framework and treats both features as critical-residue (PM1) evidence, though not at the same
strength. Cysteine-removing variants in any of the 43 cbEGF domains are raised to PM1_Strong,
while variants at the consensus calcium-binding residues [D]-X-[D/N]-[E/H]-X*m*-[D/N]-X*n*-[Y/F]
receive PM1 at moderate strength, with one exemption for Asn→Ser at the second D/N position (on
the grounds of its frequency in gnomAD). The specification justifies this ranking by the
prevalence of cysteine variants among reported cases, not by a mechanistic comparison, applying
the calcium rule to a sequence consensus instead of to observed coordination.

The two features have not been compared directly, though there is longstanding evidence that they
differ. In coagulation factor IX, whose first EGF-like domain carries the same consensus, removal
of a coordinating aspartate reduced calcium affinity more than a thousandfold and caused
hemophilia B. In fibrillin-1 itself, the N2144S substitution removes a coordinating asparagine of
cbEGF32 without altering the native fold of that domain or of its neighbors, reducing calcium
affinity approximately ninefold (from a *K*d of 1.6 mM to about 14 mM). Related work has shown
that calcium-ligand substitutions increase susceptibility to proteolysis, that their molecular
consequences depend on the domain preceding the site, and that some cbEGF substitutions impair
secretion instead of folding. Earlier simulation examined how cysteine mutations perturb calcium
binding under mechanical load, but not whether calcium-ligand substitutions destabilize the
domain.

Whether N2144S is representative of its class is not established, largely because few fibrillin-1
domains have experimental structures. Only eight of the 43 cbEGF domains are covered, by four
crystallographic and NMR fragments. We have therefore generated calcium-bound models of all 43
cbEGF domains and of 60 mutants, identifying the ligands from observed contact with the ion, for a
comparison of the energetic and geometric consequences of the two variant classes on a common
footing.

---

## Results

### Distribution of pathogenic variants across the cbEGF array

A total of 4,451 *FBN1* missense variants were assembled from ClinVar and normalized to the MANE
Select transcript NM_000138.5 and to UniProt P35555. At a two-star review floor these comprise 578
pathogenic or likely pathogenic, 823 of uncertain significance and 34 benign or likely benign
variants. Pathogenic variants fall across the cbEGF array, whereas the gnomAD population
comparison group is distributed more evenly (**Fig. 1a**).

When every cbEGF domain is aligned on its own N terminus the concentration becomes explicit
(**Fig. 1b**), and of the 425 pathogenic variants inside a cbEGF domain the tallest peaks coincide
with the six cysteines and with the calcium ligands. These positions were derived from the models
and were not assumed. Across all 43 domains the side-chain ligands were found at offsets 0 (Asp,
43 of 43 domains), +3 (Glu, 34 of 43) and +17 or +18 (Asn), the main-chain carbonyls at +1 and at
+19 to +22, and the cysteines at median offsets +4, +11, +16, +25, +27 and +40. The site is
visible in the 1.8 Å crystal structure of the cbEGF9-hyb2-cbEGF10 fragment (**Fig. 1c**), in which
the ion is held by the side chains of Asp807, Glu810 and Asn823 and by the main-chain carbonyls of
Ile808, Ser824 and Ser827.

Experimental structures cover only eight of the 43 domains, leaving the variants inside them an
unrepresentative sample, whichever positions earlier work happened to solve. The full array was
therefore folded as 21 overlapping two- and three-domain constructs on AlphaFold Server. Each
construct was cofolded with three Ca²⁺ ions (that is, protein and ions were modeled together),
five models being returned for each. Structural coverage rose from 751 variants (17%) to 3,831
(86%). Ligands were then identified by observed contact with the modeled ion, a definition not
equivalent to the sequence consensus, which flags 143 positions that bind no ion in any model.
Those 143 positions are kept throughout as a separate group and test whether the structural
definition does any work.

Four groups are compared below. The first is variants at side-chain calcium ligands (n = 271) and
the second variants at aspartate, glutamate and asparagine positions that bind no ion (n = 323,
the control set), a set needed because calcium ligands are aspartate, glutamate or asparagine by
definition. The third is variants at consensus positions that are not ligands (n = 143). The
fourth is variants that substitute a cysteine (n = 1,132), of which 851 lie in a cbEGF domain and
the remainder in the TB domains and the non-calcium EGF domains.

### Validation of the modeled calcium sites

Across the eight domains with a calcium-bound experimental reference (2W86, 1LMJ, 1UZJ and 1EMN
together with 1EMO), the modeled ion was placed a median of 0.45 Å from its experimental position.
The median backbone deviation for the same domains was 1.35 Å, placing the ion more accurately
than the protein around it in all eight (**Supplementary Fig. S1a**). Summarizing each domain by
the median over its five models gives 0.55 Å against 1.35 Å and seven of eight. The exception is
cbEGF12, where the two values are equal to within 0.02 Å.

Ligand identity agrees independently of geometry. The cbEGF9 site has been reported as three
side-chain oxygens (Asp807, Glu810 and Asn823), three main-chain carbonyls (Ile808, Ser824 and
Ser827) and a single water molecule, an assignment made from the 1.8 Å crystal structure. Our
models never saw that assignment and return exactly those six residues with the same division
between side chain and main chain, indicating that the coordination shell is recovered and not
imposed. Coordination was quantified throughout as a Brown-Altermatt bond-valence sum, a geometric
summary of how well the donor atoms surround the ion. The sum describes geometry and not binding
free energy, and was verified against CheckMyMetal on the 16 sites that service has scored, the
two agreeing at *r* = 0.961 with a constant offset of +0.147 (**Supplementary Fig. S1b**). A
constant offset cancels when a mutant is subtracted from its wild type, that subtraction being the
only operation performed here.

A site was admitted to the analysis only if its five wild-type models agreed with one another, on
criteria fixed before any mutant was read. Of the 63 calcium sites in the 21 constructs, 60 were
admitted (**Supplementary Table S10**). Two sites of the cbEGF1-3 construct failed, one on
seed-to-seed spread of valence (39.2%) and one on occupancy (four of five models). A third
failure, at one site of cbEGF5-7 (spread 16.2%), carries no variant in this study. Only one
admitted site (site 2 of cbEGF19-21) lies within 20% of any of the three limits, indicating that
the filter separates a small number of poorly converged sites from a well-behaved majority.

### Folding stability at cysteine and calcium-ligand positions

Folding stability changes were computed with two unrelated energy functions, an empirical force
field (FoldX 5.1) and a metal-aware one (PyRosetta ref2015). The Rosetta function forms explicit
bonds between the ion and its contacting atoms, making it the calculation able to register a
penalty at these positions.

Neither function registers one at the calcium ligands (**Fig. 2a,b**). In the metal-aware
calculation, cysteine-removing variants were found to cost a median of 10.81 Rosetta energy units
against 1.77 for the control set (*P* = 3 × 10⁻⁴⁷). Calcium ligands cost 1.84 and are
indistinguishable from the same control (*P* = 0.42). FoldX gives the same ordering, with a small
residual difference at calcium ligands (*P* = 0.03) computed on a fifth as many variants, because
it is restricted to the domains with experimental templates. That residual difference does not
survive correction for multiple testing (*q* = 0.06). A folding calculation needs a modeled
structure, restricting the Rosetta comparison to 1,038 of the 1,132 cysteine substitutions. The 94
excluded lie outside the span of the 21 constructs, in the N-terminal EGF-like domains, in TB1 and
in the linkers at either end of the protein. AlphaMissense works from sequence alone and so covers
all 1,132.

Calcium ligands are more buried than the control set (median relative solvent accessibility 21.7%
against 38.4%), leaving open the possibility that the comparison reflects location and not role.
That explanation is not supported when the comparison is stratified by burial (**Fig. 2c**). Among
buried positions calcium ligands were found to cost less than buried controls, 1.72 against 4.14
Rosetta energy units (*P* = 2 × 10⁻⁴). Among partly buried positions the two are level (1.92
against 2.01, *P* = 0.54). Only one of the 271 calcium ligands is substantially exposed, leaving
that band unassessed.

Sequence-based prediction gives a different result on the same variants (**Fig. 2d**).
AlphaMissense, which does not use a structure, scores calcium-ligand variants at a median of 0.989
out of 1, as damaging as cysteine variants (0.998) and far above the control set (0.260). An
orthogonal method therefore recognizes these variants as damaging while direct calculation finds
them energetically unremarkable.

### Calcium coordination in the mutant models

Since the fold survives the substitution, any consequence must lie in the site itself. A further
60 constructs were therefore cofolded with calcium, comprising 34 pathogenic side-chain ligand
variants, 12 controls from the Asp, Glu and Asn set, and 14 negative controls at positions that
touch no ion (five models each). Every mutant was compared with wild-type models of the same
construct.

Each construct binds three Ca²⁺ ions, of which a substitution can affect only the one its residue
contacts. That ion was the one scored, the other two in the same construct serving as controls,
since they lie in the same fold, come from the same five prediction runs and cannot be reached by
the substitution. Those 109 control ions, together with the 14 negative-control sites, gave 123
measurements in which no change should occur. The 95th percentile of their absolute change, 9.1%
of site valence, was the measured threshold used throughout. Three variants sat at the two sites
of the cbEGF1-3 construct that the admission filter excluded and were dropped, leaving 57.

Fifteen of 33 pathogenic calcium-ligand variants lost coordination beyond the threshold, against 1
of 24 controls (45%, 95% CI 30 to 62%, odds ratio 19, Fisher's exact *P* = 7 × 10⁻⁴) (**Fig.
3a,b**). The median change was −5.7% for pathogenic ligands, against −1.7% and +0.1% for the two
control groups. The largest loss is N2144H at −38.0%, falling at the position whose serine
substitution has been characterized experimentally, and N2144S itself was folded in a pilot series
and loses 28.5% (**Supplementary Fig. S1c**).

The fold does not move over the same models (**Fig. 3d**). Root-mean-square deviation (rmsd) of
the backbone, measured over residues within 12 Å of the ion, was a median of 0.29 Å for pathogenic
ligands and 0.24 Å for controls. The maximum across all 57 variants is 0.80 Å, a Cα displacement
smaller than a bond length. Whole-construct rmsd is far larger and far more variable, reaching
29.8 Å, a consequence of hinging between domains from seed to seed in these tandem constructs. The
fold was therefore judged locally, because on a multi-domain predicted structure a global
superposition turns such a hinge into an apparent ion displacement.

### Dependence of the outcome on the substituting residue

A rate of 45% means that some pathogenic ligand substitutions leave the site intact. These
exceptions are systematic. Substitutions that removed the contacting oxygen atom entirely reduced
coordination in 12 of 16 cases (75%). Substitutions that replaced it with another oxygen donor,
Asp→Asn, Asp→Glu, Asn→Asp and Glu→Gln, did so in 3 of 17 (18%), an odds ratio of 14 (*P* = 0.002).
In 18 of the 33 pathogenic ligand variants the substituted side chain still contacted the ion in
all five models. An atom binds the ion, not a residue identity. A conservative substitution that
keeps the atom largely keeps the site, with one systematic exception described below.

Two of these variants have been characterized experimentally, both results agreeing with the
models. N548I (cbEGF4), which renders the domain susceptible to proteolysis with local structural
change close to the mutation, loses 33.0% of its coordination here, retains no contact with the
ion in any model and moves the local backbone by 0.47 Å. E1073K was characterized in the same work
and loses 11.9% at its own site. That value is above the 9.1% threshold, though it did not clear
the per-site noise criterion used in the smaller pilot series. The pooled null replaced that
criterion for this reason.

### The donor relationship at the two kinds of consensus position

The consensus site contains two kinds of position and they do not behave alike (**Fig. 3c**). At
the aspartate and glutamate near the domain N terminus (offsets 0 and +3) the relationship is
clean. Coordination fell at none of 12 substitutions that kept an oxygen atom (median −1.0%) and
at 9 of 13 that removed it (median −12.2%, Fisher's exact *P* = 5 × 10⁻⁴). The relationship breaks
down at the asparagine at offsets +16 to +20 (the β-hydroxylation site), where three of five
substitutions that kept an oxygen atom reduced coordination, at a median of −14.6%. All three of
the Asn→Ser substitutions modeled here (N589S, N1046S and N2624S) lost coordination at a median of
−17.2%, consistent with the ninefold affinity loss measured for N2144S at the same position. A
serine hydroxyl is a shorter and weaker donor than an asparagine carboxamide, which at this
position does not reach the ion.

Asn→Ser at this position is the substitution the ClinGen *FBN1* panel exempts from PM1. It is also
the most common substitution at that position in our cohort (18 of 66 variants, four of them
classified two-star pathogenic or likely pathogenic in ClinVar). One of the four is N2144S, the
single *FBN1* variant with a measured calcium affinity, falling approximately ninefold. Our models
place it at −28.5% (**Supplementary Fig. S1c**), with N2144H at −38.0%. With only three modeled
Asn→Ser variants the comparison against all other ligand substitutions does not reach significance
(*P* = 0.13), indicating a direction without demonstrating it.

### Comparison with stability prediction and AlphaMissense

Across the same 33 variants, Rosetta ΔΔG shows no correlation with the change in coordination
(Spearman ρ = −0.09, *P* = 0.62) (**Fig. 4a**). AlphaMissense assigns a median score of 0.982 to
the variants that lose coordination and 0.982 to those that keep it (*P* = 0.61), showing no
separation between the two outcomes (**Fig. 4b**). A stability calculation therefore does not
register the damage, while a sequence-based predictor scores these variants as damaging without
separating the two outcomes. The *FBN1* panel specifies REVEL for PP3 at a threshold of 0.75 and
does not use AlphaMissense. We did not test REVEL. The saturation reported here is a property of
the predictor we measured. Whether REVEL separates these outcomes is open and easily answered.

### Enrichment of pathogenic variants at the two features

Pathogenic variants are enriched at both features relative to the population comparison group
(**Fig. 4c**). The odds ratio is 3.41 (2.21 to 5.27) at side-chain calcium ligands and 140 (93 to
205) at cysteine positions. Because the guidelines used to assign the clinical labels apply a rule
at these same positions, the comparison is partly circular and is reported for completeness (a
descriptive result and not a mechanistic one). Ligands that contribute only a main-chain carbonyl
are strongly depleted among pathogenic variants (odds ratio 0.049, 0.017 to 0.212), as expected,
since no amino-acid substitution can remove a main-chain carbonyl. That row is the negative
control the analysis passes. Consensus positions that are not ligands show no significant
enrichment (1.49, 0.90 to 2.51), the corresponding check on the structural definition.

---

## Discussion

Here we report calcium-bound models of all 43 cbEGF domains of fibrillin-1 and of 60 mutants,
together with a comparison of the two classes of missense variant that carry critical-residue
evidence in *FBN1* (calcium ligands and cbEGF cysteines). The two classes differ in kind, and not
in the degree of evidence the current weighting assigns them. Substitution of a cysteine
destabilizes the domain, which any stability-based method detects. Substitution of a
calcium-binding side chain does not destabilize the domain, by direct calculation in two
independent energy functions with composition and burial controlled. In nearly half of cases it
nonetheless loosens the ion, with the backbone around the site staying within half a bond length
of its wild-type position. This is the behavior measured for N2144S at the bench in 1999 and
inferred for the equivalent substitution in factor IX eight years earlier. It is shown here to be
the majority behavior of the class, with a chemical explanation for its exceptions.

The whole geometric argument rests on one assumption, that an AlphaFold Server model of a mutant
construct reports a real change in coordination and not a modeling artifact. The experimental
anchoring for that assumption is thin, amounting to N2144S, whose calcium affinity was measured
directly, together with qualitative agreement at N548I and E1073K, characterized by proteolytic
susceptibility and not by affinity. Three variants are a weak anchor for a class of 271. The
assumption should be tested by direct binding measurements before the rate reported here is
treated as a property of fibrillin-1 and not of the models.

The figure of 45% is a floor, for two reasons. The variants that spare the site are overwhelmingly
those that leave an oxygen atom in contact with the ion (Asp→Asn, Asp→Glu, Asn→Asp and Glu→Gln),
leaving most of the coordination intact. A structure predictor is also conservative when a single
side chain is changed. The mutant sequence differs from the wild type at one position, drawing the
model toward the wild-type fold. A real but modest loss of contact can therefore be returned as no
change. Both effects push the measured rate down and neither pushes it up.

Those variants are nonetheless classified as pathogenic, for a reason this analysis does not
identify. Three possibilities cannot be separated here (loss of β-hydroxylation at the consensus
asparagine, an affinity change too small to register as a change in modeled valence, and
overclassification, documented for *FBN1* at 31.6% outside critical residues). β-hydroxylation is
not modeled in this work. The distinction is testable, the specific prediction being that among
pathogenic calcium-ligand variants Asp→Asn and Asp→Glu substitutions should show milder
calcium-binding defects than Asp→Gly or Asp→Tyr at the same position.

Two mechanisms we did not measure, proteolytic protection and mechanical load, would compound the
effect we did. Calcium occupancy protects cbEGF domains from proteolysis, an effect that
calcium-ligand substitutions weaken by increasing proteolytic susceptibility. A site that binds
its ion more weakly therefore spends more time unprotected. The site is also mechanosensitive,
steered molecular dynamics showing cbEGF calcium binding to weaken under strain, with disulfide
mutations altering that response. A weakened site in a load-bearing microfibril may therefore fail
at forces a wild-type site tolerates. Our measurements were made at zero load on a static model
and address neither mechanism, though they establish the resting-state defect that both would
amplify.

Two of these results bear on the *FBN1* PM1 rule as written, both stated here as computational
hypotheses with experimentally testable predictions attached. The first concerns how the calcium
residues are identified, the rule reading them off the sequence consensus, which flags 143
positions that bind no calcium in any of our models. Those positions show no enrichment among
pathogenic variants (odds ratio 1.49, 0.90 to 2.51), while true side-chain ligands show 3.41 (2.21
to 5.27). The prediction is that an independent *FBN1* cohort (curated without reference to
position) will show the same split, with pathogenic variants concentrated at the structurally
defined ligands and absent from the 143. Should it hold, the same evidence would be concentrated
on fewer positions than the rule currently covers.

The second concerns the published exemption for Asn→Ser at the β-hydroxylation asparagine, the one
position at which the donor relationship breaks down. Every Asn→Ser variant modeled there lost
coordination, the best-characterized variant at the position loses ninefold affinity at the bench,
and four of the eighteen such variants in our cohort are two-star pathogenic. The prediction is
that a direct binding measurement on an isolated cbEGF pair carrying N589S, N1046S or N2624S will
show a reduction in calcium affinity comparable to the ninefold reported for N2144S, whereas
Asp→Asn or Asp→Glu at the N-terminal ligands of the same domains will not. Three modeled variants
at *P* = 0.13 cannot overturn a rule built on curated case data, nor is any such claim made here.
They are enough for the exemption to warrant re-examination once that measurement exists.

Computational evidence is now admitted to classification at calibrated strengths, with stability
prediction a reasonable proxy for one class of damage. Applied to calcium-site variants, however,
it returns values indistinguishable from those of benign substitutions. Used on its own it would
misclassify them. At a metal-binding position the property under test is whether the site still
binds, not whether the protein still folds. Our results place that property in the substituting
residue, specifically in whether it keeps an oxygen atom in contact with the ion. Neither a
stability calculation nor a sequence-based score resolves it. A donor-aware rule would be simple
to state and cheap to apply. Our data suggest it would also need to be position-aware,
conservative at the N-terminal aspartate and glutamate but not at the β-hydroxylation asparagine.
Whether such a rule classifies better than the one now in use is an empirical question, answerable
by applying both to a curated set held back from the models.

Finally, the cbEGF module is not unique to fibrillin-1. The same fold and the same consensus carry
disease variants in coagulation factor IX and, more broadly, throughout the fibrillin and LTBP
superfamily. If calcium-site substitutions behave this way here, the same reasoning should apply
wherever a small domain uses a bound ion to rigidify an interface, and the same limitation of
stability-based interpretation would follow.

---

## Limitations

**The load-bearing assumption is that mutant models report real coordination changes.** Every
geometric result here is a difference between two AlphaFold Server predictions, no calcium
affinity having been measured in this work, and the direct experimental support amounts to one
variant (N2144S) plus qualitative agreement at N548I and E1073K. That assumption is the first
thing a reader should weigh, and the first thing an experiment should test.

**The 45% disruption rate is a floor.** A structure predictor changed at one position is drawn
toward the wild-type fold, making its error on a single substitution an error toward no change.
The rate is therefore biased toward false negatives, and the true proportion of pathogenic
calcium-ligand variants that weaken their site is at least the figure reported.

**Structural coverage is modeled.** Eight of the 43 cbEGF domains have a calcium-bound
experimental structure and the remainder are predictions. Those predictions reproduce the solved
cases closely, including the position of the ion and the identity of its ligands. That is a reason
for confidence and not a substitute for measurement.

**AlphaFold3 is trained largely on holo structures** and is biased toward placing an ion
somewhere. Every measurement reported here is a difference against a wild-type model of the same
construct folded with the same number of ions, with the threshold taken from controls folded under
the same bias.

**Bond-valence sum is a geometric proxy for affinity and not a measurement of it**, validated
against CheckMyMetal though not against a binding constant. Both exist for three variants (N2144
directly, and N548I and E1073K through proteolytic susceptibility), all three agreeing in
direction and in rank, the strongest statement the data support.

**The rate rests on 33 variants**, with a 95% interval of 30 to 62%. The control arm of 24
produced a single disruption, leaving the odds ratio less precisely determined than the rate.

**The position-stratified result rests on small strata.** The β-hydroxylation asparagine stratum
contains only eight modeled variants, of which the Asn→Ser subset is three. The statistical weight
therefore rests on the aspartate and glutamate stratum (25 variants), leaving the asparagine
stratum as a signal to be tested rather than a result to act on alone.

**The benign comparison group is small.** It is made up of 34 variants (29 with structural
coverage), none of which lies at a calcium ligand or at a cysteine. A population comparison group
was substituted for the enrichment analysis for that reason.

**The enrichment analysis is descriptive and partly circular**, as described above, with the
mechanistic argument resting on the coordination measurements instead.

**The two energy functions differ slightly at calcium ligands** (*P* = 0.03 for FoldX against 0.42
for Rosetta), computed on 50 and 271 variants respectively. The FoldX difference does not survive
Benjamini-Hochberg correction, and the better-powered, metal-aware result is regarded here as the
more reliable, though both are reported.

**A steric artifact was identified and corrected.** With the backbone held rigid, large residues
replacing buried cysteines produced unrelieved collisions instead of energies, and 21.5% of values
exceeded 100 Rosetta energy units. Allowing local backbone relaxation under positional restraints
reduced this to 1.0% and changed the cysteine median from 24.76 to 10.81. The directions were
unaffected, though the earlier magnitudes should not be quoted.

**Post-translational modification is absent from the models.** cbEGF domains carry β-hydroxylation
at the consensus asparagine and O-linked glycosylation, neither required for calcium binding,
neither modeled here.

---

## Methods

**Variants.** *FBN1* (Gene ID 2200) missense variants were retrieved from ClinVar (release
2026-07-28) through the NCBI E-utilities and normalized to the MANE Select transcript NM_000138.5
and to UniProt P35555 (entry version 264, sequence MD5 f8f84ac9dc30f3ff640ad4f691fee2b1). HGVS
descriptions were validated with Mutalyzer, with the wild-type residue checked against the
reference sequence at every position. Failures are listed in `data/interim/wt_mismatch.tsv` and
were excluded. Clinical classes were assigned from ClinVar review status, using a two-star floor
for the primary analysis and a one-star sensitivity set. Conflicting classifications are labeled
and are never dropped. Population frequencies are from gnomAD v4 and pathogenicity scores from
AlphaMissense v1.0.0, joined with the same wild-type identity check.

**Structures.** Experimental structures from the PDB were used where they exist, both as gold
standard and as calibration targets. These are 2W86 (cbEGF9-hyb2-cbEGF10), 1LMJ (cbEGF12-13), 1UZJ
(cbEGF22-TB4-cbEGF23) and 1EMN together with 1EMO (cbEGF32-33). Residue numbering was mapped
through SIFTS and verified by residue identity, not by index. All 43 cbEGF domains were
additionally modeled as 21 overlapping two- and three-domain constructs on AlphaFold Server, each
cofolded with three Ca²⁺ ions and with five models per construct. Construct sequences were
verified against P35555 over the stated range before any measurement. Sixty mutant constructs were
folded the same way, each verified to differ from its wild type at exactly one position, at the
position claimed and with the residues claimed. Cofolding was chosen over the two alternatives for
placing the ion. The AlphaFold Protein Structure Database holds no model for P35555, whose length
exceeds its cap. A monomeric AlphaFold2 model would in any case contain no metal. Transplanting
the ion by homology was evaluated directly. Over 56 leave-one-out transplants between the
calibration domains it placed calcium a median of 1.44 Å from the crystallographic position,
against 0.45 Å for cofolding.

**Site definitions.** A residue was taken to be a calcium ligand if one of its oxygen, nitrogen or
sulfur atoms lies within 3.2 Å of a modeled Ca²⁺. Side-chain donors were distinguished from
main-chain donors. Cysteine roles and disulfide partners were taken from UniProt and confirmed
against the models. cbEGF domains pair C1–C3, C2–C4 and C5–C6 without exception across all 129
cbEGF disulfides. UniProt annotates no calcium site in P35555, leaving no annotation against which
the site definition could be validated. The sequence-motif definition was retained as a separate
comparison group instead.

**Energetics.** FoldX 5.1 ΔΔG was computed with RepairPDB followed by BuildModel, with mutation
strings written in author numbering by way of SIFTS. PyRosetta ref2015 ΔΔG was computed with
explicit metal bonds (`auto_setup_all_metal_bonds`), with minimization in torsion space, and with
backbone relaxation inside an 8 Å repack shell under harmonic Cα restraints (σ = 0.5 Å, following
the restrained-relaxation protocol of Nivón *et al.*) applied during minimization and removed
before scoring. Metal constraints were not used, because the packer optimizes blind to them and
the rescore then registers the violation. The two energy functions are expressed as different
quantities in different units and were never pooled. Because a ΔΔG requires a structure, the
Rosetta calculation was run on the 1,934 variants that lie in a modeled construct and belong to
one of the site classes compared here, the FoldX calculation on the 751 that fall on an
experimental template.

**Geometry.** Coordination was quantified as a Brown-Altermatt bond-valence sum over donor atoms
within 3.2 Å (R₀ = 1.967 Å for Ca–O and 2.14 Å for Ca–N, b = 0.37) and was validated against
CheckMyMetal. The sum is a geometric description of the donor shell and not a binding free energy,
supporting statements of direction and rank but not of affinity. The 3.2 Å cutoff follows the
surveyed distance distributions for calcium sites in protein structures. Sites were matched
between wild type and mutant by overlap of their coordinating residues (Jaccard ≥ 0.35, assigned
by the Hungarian algorithm) and not by proximity. A vacated site is therefore reported as vacated
and never matched to a spurious partner. Every quantity reported is a median over five models. Ion
displacement and backbone rmsd were computed after superposing the neighborhood of the site
itself, taken as residues with any atom within 12 Å of the wild-type ion. These tandem constructs
hinge between domains, a motion that a global superposition converts into an apparent ion
displacement. Sites were admitted on wild-type evidence alone, before any mutant was read.
Admission required a seed spread no greater than 15% of the site's own valence, occupancy in all
five wild-type models and a mean pLDDT of at least 70 over the coordinating residues (60 of the 63
sites qualified). Every site is listed against those three criteria in Supplementary Table S10.
Solvent accessibility was computed with FreeSASA and secondary structure with DSSP as implemented
in mkdssp.

**Statistics.** Group comparisons were made with Mann-Whitney *U* tests, with Cliff's δ and
2,000-resample bootstrap confidence intervals. Where a composition effect could confound a
comparison, the result is also reported stratified by wild-type residue. Rates were compared with
Fisher's exact test and are reported with Wilson score intervals. Benjamini-Hochberg correction is
applied within each family of tests, with every headline comparison carrying its corrected *q*
value beside the raw *P* value in the figure. Supplementary Table S6 carries both for all 29
tests. All analyses are seeded and reproduce identically on re-run, which an automated gate
enforces.

**Software.** Python 3.14.4 was used for the analysis and 3.11.15 for PyRosetta, with gemmi,
FreeSASA 2.2.1, mkdssp 4.5.7, PyMOL 3.1.0, NumPy 2.5.1, SciPy 1.18.0, pandas 3.0.5 and Matplotlib
3.11.1. Versions and paths are recorded in `manifest/tools.md`, with every download logged in
`manifest/` under its date, query and row count.

---

## Declarations

**Funding.** The author received no specific funding for this work.

**Competing interests.** The author has declared that no competing interests exist.

**Ethics.** This study used only variant records and allele frequencies already in the public
domain, from ClinVar and from the aggregated gnomAD v4 release. No identifiable individual-level
data were accessed and no human participants were recruited. The work therefore required neither
ethical approval nor informed consent.

**Author contributions.** [Author name] designed the study, wrote the analysis code, performed the
analysis and wrote the manuscript.

---

## Data and code availability

All data underlying the findings are publicly available, together with the code that produced
them.

**Repository.** The complete analysis pipeline is at [GitHub repository URL] and is archived with
a permanent identifier at Zenodo under [DOI]. The archive contains every numbered stage script,
the validation gates, the figure scripts and the manifests recording the date, query and row count
of every download.

**Source data.** ClinVar variant records were retrieved through the NCBI E-utilities for Gene ID
2200 (release 2026-07-28). Population frequencies are from gnomAD v4.1.0 and pathogenicity scores
from AlphaMissense v1.0.0 (Zenodo record 8208688). The reference sequence is UniProt P35555 entry
version 264 and the transcript is MANE Select NM_000138.5. Experimental coordinates are PDB
entries 2W86, 1LMJ, 1UZJ, 1EMN and 1EMO.

**Structural models.** All 96 AlphaFold Server jobs are deposited in the archive at five models
each, comprising the 21 wild-type array constructs, the 60 mutants and 15 calibration, pilot and
repeat constructs. Each carries its server-assigned seed, and a SHA-256 over every one of the 480
returned coordinate files is in the archive manifest.

**Underlying numbers.** These are in the supplementary workbook, one sheet at a time.
Supplementary Table S1 carries all 4,451 variants with their labels, annotation and metrics. S2
holds the 3,831 with structural coverage and S3 every modeled calcium site. S4 and S5 hold the
mutant geometry and the disruption calls. S6 lists every statistical test with its effect size and
its corrected *q*. S7 is the enrichment analysis and S8 the placement and CheckMyMetal validation.
S9 defines every column used above and S10 is the per-site quality and admission record. Each
figure is generated by exactly one script and writes a sidecar listing every value it prints,
which an automated gate re-derives from the processed tables.

---

## Figure legends

**Fig. 1 | Fibrillin-1 variants concentrate at two features of the cbEGF module.** **a**, ClinVar
pathogenic and likely pathogenic variants (above, black) and the gnomAD v4 population comparison
group (below, gray) along fibrillin-1, as counts per residue, drawn either side of the domain
architecture. **b**, Every cbEGF domain aligned on its own N terminus. Bars are pathogenic variant
counts per position (425 variants in cbEGF domains). The schematic above shows the consensus
module, with the six cysteines (circles, C1 to C6) and their three disulfides drawn as arcs, and
the calcium ligands (diamonds, filled for side-chain donors and open for main-chain carbonyl
donors) at the offsets measured across all 43 domains. Bars are shaded to match, black at cysteine
positions, mid gray at side-chain ligands and light gray elsewhere. **c**, The cbEGF9 calcium site
in the 1.8 Å crystal structure of cbEGF9-hyb2-cbEGF10 (PDB 2W86), with the three side-chain and
three main-chain donors labeled.

**Fig. 2 | Cysteine substitutions destabilize the domain and calcium-ligand substitutions do
not.** **a**, ΔΔG from PyRosetta ref2015, which forms explicit bonds between the ion and the atoms
that contact it, for the four site classes. **b**, ΔΔG from FoldX 5.1, restricted to the domains
with experimental templates. **c**, The Rosetta ΔΔG of **a** stratified by relative solvent
accessibility (RSA), comparing calcium ligands with the control set within each band. **d**,
AlphaMissense score for the same four groups, with the dashed line marking the developers'
pathogenic threshold. The cysteine row is 1,132 variants in **d** and 1,038 in **a** and **c**,
because AlphaMissense scores every variant from sequence while a ΔΔG needs a modeled construct.
The 94 cysteine substitutions outside the span of the 21 constructs account for the difference.
Points are individual variants, bars are medians and boxes the interquartile range, with n in each
axis label and off-scale points counted at the axis. Groups are separated by fill, solid, gray and
open, and not by color. Every figure in this paper is monochrome and none depends on color
reproduction. *P* values are two-sided Mann-Whitney tests and *q* values are their
Benjamini-Hochberg corrections within the family of tests each belongs to. Effect sizes are in
Supplementary Table S6. FoldX (kcal mol⁻¹) and Rosetta (Rosetta energy units) are different
quantities and are never plotted on a shared axis.

**Fig. 3 | Pathogenic calcium-ligand substitutions reduce coordination while the fold is
preserved.** **a**, Change in bond-valence sum at the calcium site each variant contacts, mutant
minus wild type, ranked, for 57 cofolded variants. The shaded band is ±9.1%, the 95th percentile
of absolute change across 123 sites at which no effect is possible (109 control ions in the same
constructs and 14 negative-control sites). N2144H is the largest loss. The serine substitution at
the same position is the one *FBN1* variant whose calcium affinity has been measured, falling
approximately ninefold with the fold unchanged. **b**, Proportion of each group whose site fell
beyond the threshold, with Wilson 95% intervals. The odds ratio and Fisher's exact *P* compare
pathogenic ligands with both control groups pooled, and *q* is the Benjamini-Hochberg correction
within the geometry family. **c**, The same change for pathogenic calcium-ligand variants only,
split by whether the substituted residue can still donate an oxygen atom and by which consensus
position is affected, either the N-terminal aspartate and glutamate (offsets 0 and +3) or the
β-hydroxylation asparagine (offsets +16 to +20). Open rings mark the three Asn→Ser substitutions
the ClinGen *FBN1* panel exempts from PM1. Counts below each group are the number that fell beyond
the threshold. **d**, Backbone rmsd over residues within 12 Å of the same ion, mutant against wild
type, on the same models.

**Fig. 4 | The reduction in coordination is not detected by the methods used as computational
evidence.** **a**, Rosetta ΔΔG against the change in calcium valence for the 33 pathogenic
calcium-ligand variants. Filled points are the variants whose coordination fell beyond the
threshold. **b**, AlphaMissense score for the same variants, split by whether coordination was
reduced or retained. **c**, Odds ratios with 95% intervals for pathogenic variants at each feature
against the 1,901-variant gnomAD population comparison group, on a logarithmic scale. Panel **c**
is descriptive and partly circular, as the note beneath it records, with the mechanistic argument
resting on **a** and **b** and on Fig. 3. Cysteine positions here are every cysteine inside a
modeled construct, 304 of the 370 pathogenic ones lying in a cbEGF domain.

**Supplementary Fig. S1 | Validation of the modeled calcium sites.** **a**, Deviation of the
modeled ion from its experimental position against the backbone deviation of the same domain, for
the eight domains with a calcium-bound structure. In the shaded half-plane the ion is placed more
accurately than the fold. **b**, Local bond-valence sum against CheckMyMetal valence over the 16
sites that service has scored. **c**, The nine-variant pilot series, scored against the per-site
wild-type seed spread. Filled points lie beyond twice that spread.

**Supplementary Fig. S2 | Both mechanisms drawn on experimental coordinates.** **a**, The
cbEGF9-hyb2-cbEGF10 fragment (PDB 2W86, 1.8 Å, crystallized in 20 mM CaCl₂), showing two cbEGF
domains either side of the hybrid domain, the six cbEGF disulfides as sticks and two Ca²⁺ ions.
**b**, Pathogenic ClinVar variants on the same fragment, with Cα drawn as a sphere, black at the
calcium ligand, as the ion itself is drawn, mid gray at cysteines and light gray elsewhere.

---
