# House style, derived from five published papers

Written 2026-09-05, in response to the criticism that the manuscript and figures used vague,
non-standard vocabulary and a rhetorical register that no journal prints.

The five papers below were read in full from `resources/papers/` and the rules in §2 and §3 were
taken from them. They were chosen because all five are structural or biophysical studies of the
same domain in the same protein, so their vocabulary is the vocabulary this paper must use.

## 1. The reference set

| Paper | Journal | What was taken from it |
|---|---|---|
| Jensen SA *et al.* (2009) Structure and interdomain interactions of a hybrid domain | *Structure* 17:759 | Results-section register, how a coordination shell is described in prose, r.m.s.d. reporting |
| Kettle S *et al.* (1999) Defective calcium binding to fibrillin-1 | *J Mol Biol* 285:1277 | How a single substitution's consequences are stated and hedged, affinity reporting |
| Downing AK *et al.* (1996) Solution structure of a pair of cbEGF domains | *Cell* 85:597 | Handling of an unexplained observation, sequence-analysis prose |
| Handford PA *et al.* (1995) The calcium binding properties and molecular organization of EGF-like domains | *J Biol Chem* 270:6751 | Mutagenesis-result phrasing, consensus-motif nomenclature |
| Godwin ARF *et al.* (2023) Fibrillin microfibril structure | *Nat Struct Mol Biol* 30:608 | Modern Discussion structure, how alternatives and limitations are conceded |

## 2. Register rules taken from those papers

1. **Sentences are of even, moderate length**, roughly 20 to 35 words. None of the five papers
   contains a short dramatic sentence in Results or Discussion. Constructions of the form
   "Standard rules treat these as equivalent. They are not." do not occur and are not used here.
2. **Conclusions are attributed to the data and hedged in the standard way.** The verbs the five
   papers use are *suggest*, *indicate*, *are consistent with*, *may reflect*, *is likely due to*,
   *we conclude*, *taken together*. A result is never announced, it is reported and then
   interpreted.
3. **Contrast is carried by conjunctions**, not by juxtaposed fragments. The five papers use
   *whereas*, *although*, *in contrast to*, *by comparison*. Kettle 1999 writes "This is in direct
   contrast to those interactions between cbEGF32 and 33 that result in a rigid linkage."
4. **Uncertainty is stated plainly rather than dramatised.** Downing 1996 writes "The reason for
   the approximate 25-fold difference in calcium affinity observed for the binding sites in the N-
   and C-terminal domains is also unclear." Jensen 2009 writes "The functional significance of the
   rearrangement of disulphides in a hyb domain compared with a TB domain is unclear."
5. **Section headings name their subject as a neutral noun phrase.** The corpus carries both
   registers — Jensen 2009 heads sections "Structure of the cbEGF domain pair" and "Dynamics of
   the cbEGF domain pair", and also "cbEGF9-hyb2-cbEGF10 and cbEGF22-TB4-cbEGF23 constructs of
   fibrillin-1 have a similar shape". **The mentor chose the noun-phrase register on 2026-09-06**,
   for the traditional-journal convention, and it applies to the paper title as well. A heading
   names what the section is about; the finding belongs in the section. `test_manuscript_v4.py`
   V9 fails if a heading picks up a finite verb again.
6. **Abbreviations are defined once and then used without further comment** — cbEGF, TB, MFS,
   ECM, r.m.s.d., *K*d, RSA, BVS, ΔΔG.
7. **Numbers carry their units and their n**, in parentheses, at the point of the claim.
8. **The Discussion opens by restating what the study produced**, in one sentence, before any
   interpretation. Godwin 2023 opens "Here we present the cryo-EM structure of a native fibrillin
   microfibril from mammalian tissue." Jensen 2009 opens "Our study has provided the first
   high-resolution structure and possible evolutionary origin of a hybrid domain."

## 3. Three local rules, imposed by the mentor rather than taken from the papers

9. **No colons and no semicolons.** The reference papers use both, Kettle 1999 heavily. They are
   nonetheless not used anywhere in this manuscript or in any figure. A list becomes a sentence
   with *and*, an apposition becomes a parenthesis or a separate sentence.
10a. **American spelling, from 2026-09-06.** `-ise` becomes `-ize` (destabilize, characterize,
    normalize, minimization, summarizing, recognizes, optimizes), `-our` becomes `-or` (color,
    behavior, neighbor), a doubled `l` before a suffix becomes single (modeled, modeling,
    labeled, signaling), and it is disulfide, sulfur, hemophilia, gray, artifact, toward.
    **Four of the five reference papers are British journals**, so every rewrite made against
    them pulls the spelling back — `disulphide` appears 31 times in the corpus and `disulfide`
    only 3. Watch for it. `disulfide` is nonetheless safe for the terminology gate because
    Jensen 2009 and Handford 1995 both use it. `grey` in the figure scripts is the name of a
    colour-key entry, `P.C["grey"]`, and is never drawn, so it stays.
    `test_manuscript_v4.py` **V11** checks the draft and every string a figure script draws.

10. **No citations in the manuscript source.** They will be inserted from Mendeley. The keyed
    source with all 63 references is preserved at `results/manuscript_paper.md` and the verified
    reference list at `data/processed/references_final.tsv`, so nothing is lost.

## 3b. How a long sentence is built — measured 2026-09-06, after the counts had all passed

The mentor read the draft with every §2 and `63_style_extract.py` aggregate inside the corpus
range and said the multi-clause sentences still read as machine-written. That was correct, and
the reason is that sentence *length* and sentence *shape* are different things.
`scripts/65_sentence_shape.py` measures the shape. What it found on the draft that had passed
every other check:

| | five papers | that draft |
|---|---|---|
| sentences welding a second independent clause | 4.0% | **28.6%** |
| long sentences (30+ words) doing it | 7.5% | **48.2%** |
| median parentheses in a long sentence | 1 | **0** |
| sentences with a participial tail | 2.7% | 1.0% |
| sentences opening on a subordinate clause | 3.0% | 1.0% |

**A published long sentence has one main clause.** Its length comes from noun phrases, from one
or two parentheses and from a trailing participial or adjectival phrase. Jensen 2009 —
"Superpositions of the individual domains cbEGF9, hyb2 and cbEGF10 onto the corresponding domains
in cbEGF22-TB4-cbEGF23 (Figure 3B) showed higher degrees of similarity (rmsd = 0.829, 1.844, and
1.668, respectively), consistent with the conserved nature of the cbEGF fold and similarities
between TB and hyb domain sequences." One subject, one verb, two parentheses, one tail.

**What it never does is bolt a second independent clause on with `, and the model is ...`.** That
is the construction to avoid, and it is the one an unedited draft reaches for. When two facts
belong in one sentence, use one of the four devices the corpus uses, in this order.

1. **A relative clause** — `, which flags 143 positions that bind no calcium`.
2. **A subordinator** — `Because every geometric statement depends on ...`, `Since the fold
   survives ...`, `whereas substitution of a calcium ligand did not`.
3. **A participial or absolute phrase** — `, leaving the odds ratio less precisely determined`,
   `, the other two ions serving as controls`, `, consistent with the ninefold affinity loss`.
4. **A parenthesis** — demote an appositive into an aside. `Only one admitted site (site 2 of
   cbEGF19-21) lies within 20% of any limit`.

**Do not fix this by splitting everything into short sentences.** That was tried first and it
overshot in the other direction — median 17 words against the corpus 23, and 37% of sentences
under 15 words against 21%. Choppy reads as machine-written too. Four passes were needed, and
passes one and three overshot in opposite directions, so change this by measuring.

### The ceiling, which is the rule that actually matters

**Count clause markers one sentence at a time. An average cannot see one bad sentence.** Every
distribution check above passed on a draft whose abstract still read

> Substitutions that removed the contacting oxygen atom broke the site, whereas those supplying
> another oxygen atom largely spared it, with one exception at the β-hydroxylated asparagine,
> where every Asn→Ser substitution modelled lost coordination.

Four clauses, in the abstract. The same instrument run over the five papers one sentence at a
time gives the ceiling the averages hid.

| clause markers in one sentence | five papers (n = 1,164) |
|---|---|
| one | 79.0% |
| two | 18.3% |
| three | 2.7% |
| **four or more** | **0.0% — none, in 1,164 sentences** |

So four is not rare in this literature, it is absent from it. **No sentence in this manuscript may
carry more than three, and three should stay near 3%.** Count "of which" as one marker, not two.

`test_manuscript_v4.py` **V10** enforces all of it — welding at or below the corpus rate, a
parenthesis in the median long sentence, three-marker sentences at or below 3.5%, and **no
sentence heavier than the heaviest published one**. The last of those is a ceiling rather than an
average, and it is the check that would have caught the abstract.

## 3c. No sentence may announce a result instead of stating it

The `academic-paper-composer` checklist reduces to one rule, **specific, not vague**, and the five
papers obey it almost absolutely. A sentence that names a category without giving its content
occurs **twice in 1,164 published sentences (0.2%)**, and both of those are real claims. This
draft carried **nine (3.6%)**, each one followed immediately by the sentence that said the thing.

| deleted | because the next sentence already said |
|---|---|
| "The outcome was decided by the replacing residue." | which substitutions broke the site and which spared it |
| "Every geometric statement below depends on the models placing calcium correctly. The models were therefore measured first against the domains for which the answer is known." | the 0.45 Å measurement, under a heading that already reads *Validation of the modeled calcium sites* |
| "Its value lies in two internal checks." | both checks, each of which names itself |
| "There is a broader implication for variant interpretation." | the implication |

Two were rewritten to say the thing rather than point at it. "The two classes differ in mechanism,
though not in the way the current evidence weighting encodes" became "The two classes differ in
kind, and not in the degree of evidence the current weighting assigns them". "Two further
mechanisms lie outside what we measured" became "Two mechanisms we did not measure, proteolytic
protection and mechanical load, ...".

**Two are kept, and the distinction matters.** A gap statement ("The two features have not been
compared directly") and an inference from the previous section ("Since the fold survives the
substitution, any consequence must lie in the site itself") are claims. A throat-clear is not.
The bolded lead of a Limitations paragraph is also exempt, since stating a claim and then
evidencing it is that section's deliberate structure.

`test_manuscript_v4.py` **V12** fails if more than two survive.

## 4. Figure vocabulary

Category labels name the class in the terms the literature uses for it, not in a description of
how the class was computed. Jensen 2009 calls a coordinating residue a "Ca2+ ligand" and this
paper does the same.

| Replaced | Used now |
|---|---|
| Coordinates Ca²⁺ | Ca²⁺ ligand (side chain) |
| Asp/Glu/Asn, no Ca²⁺ contact | Non-ligand Asp/Glu/Asn |
| Matches motif, no Ca²⁺ contact | Consensus motif, non-ligand |
| Cysteine removed | Cysteine substitution |
| Keeps an O donor / Removes the donor | Oxygen donor retained / Oxygen donor lost |
| Site disrupted / Site intact | Coordination reduced / Coordination retained |
| The fold, same models | Local backbone r.m.s.d. |
| The same ΔΔG, stratified by burial | Rosetta ΔΔG by relative solvent accessibility |
| Stability against coordination | Rosetta ΔΔG versus change in Ca²⁺ valence |
| AlphaMissense, sequence only | AlphaMissense score |
| Enrichment against a population comparison | Enrichment relative to gnomAD v4 |
| Buried, <20% exposed / Partly buried, 20–50% | RSA < 20% / RSA 20–50% |
