# Sentence shape — five published papers against `manuscript_v4.md`

Produced by `scripts/65_sentence_shape.py`. `63_style_extract.py` measures how long a sentence is; this measures how it was built. A published long sentence carries one main clause and gets its length from noun phrases, parentheses and a trailing participial phrase. Welding a second independent clause on with `, and the model is ...` is the construction that reads as machine-written.

| measure | Jensen 2009 (Structure) | Kettle 1999 (JMB) | Downing 1996 (Cell) | Handford 1995 (JBC) | Godwin 2023 (NSMB) | — pooled — | THIS PAPER |
|---|---|---|---|---|---|---|---|
| sentences | 223 | 219 | 203 | 142 | 377 | 1164 | 270 |
| long sentences (30+ words) | 91 | 55 | 74 | 56 | 96 | 372 | 51 |
| % of ALL sentences welding a second clause | 5.4% | 3.2% | 7.4% | 3.5% | 1.9% | 4.0% | 2.2% |
| % of LONG sentences welding a second clause | 7.7% | 9.1% | 12.2% | 3.6% | 5.2% | 7.5% | 7.8% |
| welded clauses per long sentence | 0.09 | 0.11 | 0.12 | 0.04 | 0.06 | 0.08 | 0.08 |
| % balanced pair / not-X-but-Y | 0.0% | 0.5% | 0.0% | 0.7% | 0.3% | 0.3% | 0.0% |
| median parentheses in a long sentence | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| median commas in a long sentence | 2.0 | 1.0 | 3.0 | 1.0 | 2.0 | 2.0 | 2.0 |
| % with a participial tail | 3.1% | 2.7% | 2.0% | 2.8% | 2.7% | 2.7% | 2.2% |
| % opening on a subordinate clause | 3.1% | 1.8% | 1.0% | 6.3% | 3.4% | 3.0% | 2.2% |
| % carrying ONE clause marker | 78.9% | 76.7% | 81.3% | 76.1% | 80.1% | 79.0% | 71.9% |
| % carrying THREE | 3.1% | 4.6% | 2.0% | 1.4% | 2.7% | 2.8% | 2.6% |
| % carrying FOUR or more | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| heaviest sentence in the text | 3 | 3 | 3 | 3 | 3 | 3 | 3 |

## Where this draft sits

- % of ALL sentences welding a second clause — corpus 4.0%, ours 2.2% — under
- % of LONG sentences welding a second clause — corpus 7.5%, ours 7.8% — ok
- welded clauses per long sentence — corpus 0.08, ours 0.08 — ok
- % balanced pair / not-X-but-Y — corpus 0.3%, ours 0.0% — under
- median parentheses in a long sentence — corpus 1.0, ours 1.0 — ok
- median commas in a long sentence — corpus 2.0, ours 2.0 — ok
- % with a participial tail — corpus 2.7%, ours 2.2% — ok
- % opening on a subordinate clause — corpus 3.0%, ours 2.2% — ok
- % carrying ONE clause marker — corpus 79.0%, ours 71.9% — ok
- % carrying THREE — corpus 2.8%, ours 2.6% — ok
- % carrying FOUR or more — corpus 0.0%, ours 0.0% — ok
- heaviest sentence in the text — corpus 3, ours 3 — ok

## The 7 sentences carrying three or more clause markers

The corpus carries three in 2.7% of sentences and four in none of 1,164. Anything listed at four has to go.

**3 markers, 44 words**

> Cysteine-removing variants in any of the 43 cbEGF domains are raised to PM1_Strong, while variants at the consensus calcium-binding residues [D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F] receive PM1 at moderate strength, with one exemption for Asn→Ser at the second D/N position (on the grounds of its frequency in gnomAD).

**3 markers, 23 words**

> A stability calculation therefore does not register the damage, while a sequence-based predictor scores these variants as damaging without separating the two outcomes.

**3 markers, 21 words**

> Substitution of a residue that binds the ion removes a contact with an ion that is not part of the protein.

**3 markers, 20 words**

> Earlier simulation examined how cysteine mutations perturb calcium binding under mechanical load, but not whether calcium-ligand substitutions destabilize the domain.

**3 markers, 20 words**

> At a metal-binding position the property under test is whether the site still binds, not whether the protein still folds.

**3 markers, 19 words**

> Substitutions that removed the contacting oxygen atom broke the site, whereas those supplying another oxygen atom largely spared it.

**3 markers, 18 words**

> Whether N2144S is representative of its class is not established, largely because few fibrillin-1 domains have experimental structures.


## The 6 sentences that weld clauses, heaviest first

**1 welded, 44 words, 1 parens**

> Cysteine-removing variants in any of the 43 cbEGF domains are raised to PM1_Strong, while variants at the consensus calcium-binding residues [D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F] receive PM1 at moderate strength, with one exemption for Asn→Ser at the second D/N position (on the grounds of its frequency in gnomAD).

**1 welded, 38 words, 1 parens**

> Every geometric result here is a difference between two AlphaFold Server predictions, no calcium affinity having been measured in this work, and the direct experimental support amounts to one variant (N2144S) plus qualitative agreement at N548I and E1073K.

**1 welded, 37 words, 1 parens**

> The sum describes geometry and not binding free energy, and was verified against CheckMyMetal on the 16 sites that service has scored, the two agreeing at r = 0.961 with a constant offset of +0.147 (Supplementary Fig.

**1 welded, 32 words, 0 parens**

> Related work has shown that calcium-ligand substitutions increase susceptibility to proteolysis, that their molecular consequences depend on the domain preceding the site, and that some cbEGF substitutions impair secretion instead of folding.

**1 welded, 23 words, 0 parens**

> A stability calculation therefore does not register the damage, while a sequence-based predictor scores these variants as damaging without separating the two outcomes.

**1 welded, 20 words, 0 parens**

> Earlier simulation examined how cysteine mutations perturb calcium binding under mechanical load, but not whether calcium-ligand substitutions destabilize the domain.

