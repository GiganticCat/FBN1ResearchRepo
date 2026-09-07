# Prose diagnostic — `results/manuscript_v4.md`

Measured with `scripts/62_prose_diagnostic.py`, applying the clarity diagnostic from the `academic-paper-composer` skill. Not an opinion about the writing, a count.

## 1. Sentence load

- 315 sentences, median 21 words and 1 clause markers
- **33 sentences exceed the limit** (more than 6 clause markers or more than 34 words)
- that is 10% of the paper

### The worst offenders, heaviest first

**5 clauses, 40 words, 4 asides** — *Figure legends*

> The cysteine row is 1,132 variants in **d** and 1,038 in **a** and **c**, because AlphaMissense scores every variant from sequence while a ΔΔG needs a modelled construct, and 94 cysteine substitutions lie outside the span of the 21 constructs.

**4 clauses, 50 words, 12 asides** — *Distribution of pathogenic variants across the cbEGF array*

> Across all 43 domains the side-chain ligands were found at offsets 0 (Asp, 43 of 43 domains), +3 (Glu, 34 of 43) and +17 or +18 (Asn), the main-chain carbonyls at +1 and at +19 to +22, and the cysteines at median offsets +4, +11, +16, +25, +27 and +40.

**4 clauses, 46 words, 5 asides** — *Figure legends*

> **c**, The same change for pathogenic calcium-ligand variants only, split by whether the substituted residue can still donate an oxygen atom and by which consensus position is affected, either the N-terminal aspartate and glutamate (offsets 0 and +3) or the β-hydroxylation asparagine (offsets +16 to +20).

**4 clauses, 44 words, 6 asides** — *Distribution of pathogenic variants across the cbEGF array*

> The first is variants at side-chain calcium ligands (n = 271) and the second variants at aspartate, glutamate and asparagine positions that bind no ion (n = 323, the control set), a set needed because calcium ligands are aspartate, glutamate or asparagine by definition.

**4 clauses, 43 words, 3 asides** — *Methods*

> Because a ΔΔG requires a structure, the Rosetta calculation was run on the 1,934 variants that lie in a modelled construct and belong to one of the site classes compared here, the FoldX calculation on the 751 that fall on an experimental template.

**4 clauses, 41 words, 4 asides** — *Distribution of pathogenic variants across the cbEGF array*

> The site is visible in the 1.8 Å crystal structure of the cbEGF9-hyb2-cbEGF10 fragment (**Fig. 1c**), in which the ion is held by the side chains of Asp807, Glu810 and Asn823 and by the main-chain carbonyls of Ile808, Ser824 and Ser827.

**4 clauses, 38 words, 3 asides** — *Abstract*

> Substitution of a cysteine destabilised the domain in two independent energy functions, whereas substitution of a calcium ligand did not, its computed folding penalty being indistinguishable from that of aspartate, glutamate and asparagine positions that bind no ion.

**4 clauses, 37 words, 3 asides** — *Validation of the modelled calcium sites*

> The sum describes geometry and not binding free energy, and was verified against CheckMyMetal on the 16 sites that service has scored, the two agreeing at *r* = 0.961 with a constant offset of +0.147 (**Supplementary Fig.

**4 clauses, 37 words, 5 asides** — *Calcium coordination in the mutant models*

> A further 60 constructs were therefore cofolded with calcium, comprising 34 pathogenic side-chain ligand variants, 12 controls from the Asp, Glu and Asn set, and 14 negative controls at positions that touch no ion (five models each).

**4 clauses, 36 words, 3 asides** — *Figure legends*

> **Fig. 2 | Cysteine substitutions destabilise the domain and calcium-ligand substitutions do not.** **a**, ΔΔG from PyRosetta ref2015, which forms explicit bonds between the ion and the atoms that contact it, for the four site classes.

**3 clauses, 49 words, 3 asides** — *Limitations*

> **The load-bearing assumption is that mutant models report real coordination changes.** Every geometric result here is a difference between two AlphaFold Server predictions, no calcium affinity having been measured in this work, and the direct experimental support amounts to one variant (N2144S) plus qualitative agreement at N548I and E1073K.

**3 clauses, 45 words, 3 asides** — *Abstract*

> Most missense variants that cause Marfan syndrome (MFS) fall at one of two features of the calcium-binding epidermal growth factor-like (cbEGF) module of fibrillin-1, either the six cysteines that form its three disulphide bridges or the small group of residues that bind the calcium ion.

**3 clauses, 45 words, 6 asides** — *Figure legends*

> The schematic above shows the consensus module, with the six cysteines (circles, C1 to C6) and their three disulphides drawn as arcs, and the calcium ligands (diamonds, filled for side-chain donors and open for main-chain carbonyl donors) at the offsets measured across all 43 domains.

**3 clauses, 42 words, 2 asides** — *Introduction*

> We have therefore generated calcium-bound models of all 43 cbEGF domains and of 60 mutants, identifying the ligands from observed contact with the ion, for a comparison of the energetic and geometric consequences of the two variant classes on a common footing.

**3 clauses, 40 words, 2 asides** — *Distribution of pathogenic variants across the cbEGF array*

> When every cbEGF domain is aligned on its own N terminus the concentration becomes explicit (**Fig. 1b**), and of the 425 pathogenic variants inside a cbEGF domain the tallest peaks coincide with the six cysteines and with the calcium ligands.

**3 clauses, 39 words, 3 asides** — *Methods*

> Cofolding was chosen over the two alternatives for placing the ion, since the AlphaFold Protein Structure Database holds no model for P35555, whose length exceeds its cap, and a monomeric AlphaFold2 model would in any case contain no metal.

**3 clauses, 38 words, 2 asides** — *Discussion*

> Here we report calcium-bound models of all 43 cbEGF domains of fibrillin-1 and of 60 mutants, together with a comparison of the two classes of missense variant that carry critical-residue evidence in *FBN1* (calcium ligands and cbEGF cysteines).

**3 clauses, 38 words, 3 asides** — *Limitations*

> **Structural coverage is modelled.** Although eight of the 43 cbEGF domains have a calcium-bound experimental structure, the remainder are predictions, which reproduce the solved cases closely (including the position of the ion and the identity of its ligands).

## 2. Undefined vocabulary

- 18 domain terms used, **16 never explained in plain words before first use**

| term | uses | first use at word | plain wording it needs |
|---|---|---|---|
| ligand | 53 | 5 | a residue that contacts the ion |
| construct | 25 | 142 | the two- or three-domain piece that was modelled |
| consensus | 16 | 162 | the recurring sequence pattern |
| coordination | 25 | 212 | the set of contacts holding the ion |
| valence | 12 | 667 | the same score, used as shorthand |
| bond-valence sum | 5 | 1540 | a geometric score of how well an ion is held |
| donor | 13 | 1549 | the atom that makes the contact |
| seed | 8 | 1662 | one run of the structure predictor |
| Rosetta energy units | 4 | 1803 | the arbitrary units of that stability score |
| relative solvent accessibility | 2 | 1943 | how exposed a residue is |
| odds ratio | 8 | 2294 | how many times more likely |
| ΔΔG | 9 | 2954 | the change in folding stability |
| Jaccard | 1 | 5465 | an overlap score |
| pLDDT | 1 | 5581 | the predictor's own confidence score |
| Cliff's | 1 | 5633 | an effect-size statistic |
| Wilson | 2 | 5668 | a kind of confidence interval |

## 3. What this means

The skill's test for jargon is whether a graduate student in an adjacent field could follow the text. On the count above the paper fails that test at 16 terms, and it asks the reader to hold more than 6 clauses in 33 sentences.

