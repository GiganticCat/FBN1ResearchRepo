# Reference verification — final list

Generated 2026-09-02T07:15:22.876790+00:00 by `scripts/57_references.py` helpers over
`data/processed/references_final.tsv`. Every DOI was resolved independently against the
Crossref REST API and Europe PMC, and the returned title and first author compared with the
citation as written in the manuscript.

- references: **63**
- resolved at Crossref: **63**, at Europe PMC: **56**
- unresolved at both: **0**
- flagged for review: **2**

## Flagged

| key | DOI | note |
|---|---|---|
| uniprot2023 | `10.1093/nar/gkac1052` | first author: record says Bateman A |
| loeys2010 | `10.1136/jmg.2009.072785` | title mismatch: DOI resolves to "The revised Ghent nosology for the Marfan syndrome: Table 1" |

**Both flags were adjudicated and neither is an error.**

* `uniprot2023` — Crossref lists Bateman A as first author of the UniProt update. The paper is
  conventionally cited, and asks to be cited, as *The UniProt Consortium*. Left as written.
* `loeys2010` — Crossref's indexed title carries a trailing ": Table 1", which is an artifact of
  how the table's own DOI was registered, not part of the article title. Left as written.

Year disagreements between Crossref and Europe PMC were also reviewed and are not errors: they
are the online-publication-year versus issue-year split that *Nucleic Acids Research* database
issues and *Nature* papers routinely carry. A citation carries the issue year, which is what the
manuscript uses.

Two genuine errors were found and corrected during preparation: the Godwin 2023 title had been
truncated before "affecting a key regulatory latent TGFβ-binding site", and the Delhomme 2022
title was a paraphrase rather than the published title. A third was caught before it entered the
manuscript at all — the DOI first recorded for Mutalyzer resolved to a paper on muscarinic
acetylcholine receptor mutations, and was replaced with Wildeman 2008.

## Resolved list

| # | key | title | journal | year |
|---|---|---|---|---|
| 1 | sakai1986 | Fibrillin, a new 350-kD glycoprotein, is a component of extracellular microfibrils. | The Journal of cell biology | 1986 |
| 2 | dietz1991 | Marfan syndrome caused by a recurrent de novo missense mutation  in the fibrillin gene | Nature | 1991 |
| 3 | handford1995 | The Calcium Binding Properties and Molecular Organization of Epidermal Growth Factor-like Domains in Human Fibrillin-1 | Journal of Biological Chemistry | 1995 |
| 4 | jensen2012 | Dissecting the Fibrillin Microfibril: Structural Insights into Organization and Function | Structure | 2012 |
| 5 | downing1996 | Solution Structure of a Pair of Calcium-Binding Epidermal Growth Factor-like Domains: Implications for the Marfan Syndrome and Other Genetic Disorders | Cell | 1996 |
| 6 | knott1996 | Calcium Binding Properties of an Epidermal Growth Factor-like Domain Pair from Human Fibrillin-1 | Journal of Molecular Biology | 1996 |
| 7 | godwin2023 | Fibrillin microfibril structure identifies long-range effects of inherited pathogenic mutations affecting a key regulatory latent TGFβ-binding site | Nat Struct Mol Biol | 2023 |
| 8 | baudhuin2019 | Variability in gene-based knowledge impacts variant classification: an analysis of FBN1 missense variants in ClinVar | Eur J Hum Genet | 2019 |
| 9 | drackley2024 | Interpretation and classification of FBN1 variants associated with Marfan syndrome: consensus recommendations from the Clinical Genome Resource’s FBN1 variant curation expert panel | Genome Med | 2024 |
| 10 | handford1991 | Key residues involved in calcium-binding motifs in EGF-like domains | Nature | 1991 |
| 11 | kettle1999 | Defective Calcium Binding to Fibrillin-1: Consequence of an N2144S Change for Fibrillin-1 Structure and Function | Journal of Molecular Biology | 1999 |
| 12 | reinhardt2000 | Mutations in Calcium-binding Epidermal Growth Factor Modules Render Fibrillin-1 Susceptible to Proteolysis | Journal of Biological Chemistry | 2000 |
| 13 | mcgettrick2000 | Molecular effects of calcium binding mutations in Marfan syndrome depend on domain context | Human Molecular Genetics | 2000 |
| 14 | whiteman2003 | Defective secretion of recombinant fragments of fibrillin-1: implications of protein misfolding for the pathogenesis of Marfan syndrome and related disorders | Human Molecular Genetics | 2003 |
| 15 | whiteman2007 | Cellular and molecular studies of Marfan syndrome mutations identify co-operative protein folding in the cbEGF12–13 region of fibrillin-1 | Human Molecular Genetics | 2007 |
| 16 | haller2020 | Steered molecular dynamic simulations reveal Marfan syndrome mutations disrupt fibrillin-1 cbEGF domain mechanosensitive calcium binding | Sci Rep | 2020 |
| 17 | rao1995 | The structure of a Ca2+-binding epidermal growth factor-like domain: Its role in protein-protein interactions | Cell | 1995 |
| 18 | lee2004 | Structure of the Integrin Binding Fragment from Fibrillin-1 Gives New Insights into Microfibril Organization | Structure | 2004 |
| 19 | smallridge2003 | Solution Structure and Dynamics of a Calcium Binding Epidermal Growth Factor-like Domain Pair from the Neonatal Region of Human Fibrillin-1 | Journal of Biological Chemistry | 2003 |
| 20 | jensen2009 | Structure and Interdomain Interactions of a Hybrid Domain: A Disulphide-Rich Module of the Fibrillin/LTBP Superfamily of Matrix Proteins | Structure | 2009 |
| 21 | landrum2018 | ClinVar: improving access to variant interpretations and supporting evidence | Nucleic Acids Research | 2017 |
| 22 | karczewski2020 | The mutational constraint spectrum quantified from variation in 141,456 humans | Nature | 2020 |
| 23 | chen2024 | A genomic mutational constraint map using variation in 76,156 human genomes | Nature | 2023 |
| 24 | abramson2024 | Accurate structure prediction of biomolecular interactions with AlphaFold 3 | Nature | 2024 |
| 25 | schymkowitz2005 | The FoldX web server: an online force field | Nucleic Acids Research | 2005 |
| 26 | alford2017 | The Rosetta All-Atom Energy Function for Macromolecular Modeling and Design | J. Chem. Theory Comput. | 2017 |
| 27 | cheng2023 | Accurate proteome-wide missense variant effect prediction with AlphaMissense | Science | 2023 |
| 28 | brown1985 | Bond-valence parameters obtained from a systematic analysis of the Inorganic Crystal Structure Database | Acta Crystallogr B Struct Sci | 1985 |
| 29 | zheng2017 | <i>CheckMyMetal</i>
                    : a macromolecular metal-binding validation tool | Acta Crystallogr D Struct Biol | 2017 |
| 30 | pejaver2022 | Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria | The American Journal of Human Genetics | 2022 |
| 31 | dana2019 | SIFTS: updated Structure Integration with Function, Taxonomy and Sequences resource allows 40-fold increase in coverage of structure-based annotations for proteins | Nucleic Acids Research | 2018 |
| 32 | uniprot2023 | UniProt: the Universal Protein Knowledgebase in 2023 | Nucleic Acids Research | 2022 |
| 33 | mitternacht2016 | FreeSASA: An open source C library for solvent accessible surface area calculations | F1000Res | 2016 |
| 34 | kabsch1983 | Dictionary of protein secondary structure: Pattern recognition of hydrogen‐bonded and geometrical features | Biopolymers | 1983 |
| 35 | richards2015 | Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology | Genetics in Medicine | 2015 |
| 36 | loeys2010 | The revised Ghent nosology for the Marfan syndrome: Table 1 | J Med Genet | 2010 |
| 37 | milewicz2021 | Marfan syndrome | Nat Rev Dis Primers | 2021 |
| 38 | jumper2021 | Highly accurate protein structure prediction with AlphaFold | Nature | 2021 |
| 39 | varadi2024 | AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences | Nucleic Acids Research | 2023 |
| 40 | hekkelman2023 | AlphaFill: enriching AlphaFold models with ligands and cofactors | Nat Methods | 2022 |
| 41 | asano2022 | The Multiple Functions of Fibrillin-1 Microfibrils in Organismal Physiology | IJMS | 2022 |
| 42 | delhomme2022 | Genotype-Mitral Valve Phenotype Correlations in Marfan Syndrome With FBN1 Pathogenic Variants | JACC: Advances | 2022 |
| 43 | vallat2026 | RCSB Protein Data Bank: Delivering integrative structures alongside experimental structures and computed structure models | Nucleic Acids Research | 2025 |
| 44 | ioannidis2016 | REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants | The American Journal of Human Genetics | 2016 |
| 45 | benjamini1995 | Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing | Journal of the Royal Statistical Society Series B: Statistical Methodology | 1995 |
| 46 | mann1947 | On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other | Ann. Math. Statist. | 1947 |
| 47 | wilson1927 | Probable Inference, the Law of Succession, and Statistical Inference | Journal of the American Statistical Association | 1927 |
| 48 | cliff1993 | Dominance statistics: Ordinal analyses to answer ordinal questions. | Psychological Bulletin | 1993 |
| 49 | harris2020 | Array programming with NumPy | Nature | 2020 |
| 50 | virtanen2020 | SciPy 1.0: fundamental algorithms for scientific computing in Python | Nat Methods | 2020 |
| 51 | hunter2007 | Matplotlib: A 2D Graphics Environment | Comput. Sci. Eng. | 2007 |
| 52 | touw2015 | A series of PDB-related databanks for everyday needs | Nucleic Acids Research | 2014 |
| 53 | nivon2013 | A Pareto-Optimal Refinement Method for Protein Design Scaffolds | PLoS ONE | 2013 |
| 54 | harding2006 | Small revisions to predicted distances around metal sites in proteins | Acta Crystallogr D Biol Crystallogr | 2006 |
| 55 | morales2022 | A joint NCBI and EMBL-EBI transcript set for clinical genomics and research | Nature | 2022 |
| 56 | dendunnen2016 | HGVS Recommendations for the Description of Sequence Variants: 2016 Update | Human Mutation | 2016 |
| 57 | wildeman2008 | Improving sequence variant descriptions in mutation databases and literature using the Mutalyzer sequence variation nomenclature checker | Hum. Mutat. | 2007 |
| 58 | rehm2015 | ClinGen — The Clinical Genome Resource | N Engl J Med | 2015 |
| 59 | rees1988 | The role of beta‐hydroxyaspartate and adjacent carboxylate residues in the first EGF domain of human factor IX. | EMBO J | 1988 |
| 60 | neptune2003 | Dysregulation of TGF-β activation contributes to pathogenesis in Marfan syndrome | Nat Genet | 2003 |
| 61 | ramirez2009 | Biogenesis and function of fibrillin assemblies | Cell Tissue Res | 2009 |
| 62 | collodberoud2003 | Update of the UMD-<i>FBN1</i>mutation database and creation of an<i>FBN1</i>polymorphism database | Hum. Mutat. | 2003 |
| 63 | wojdyr2022 | GEMMI: A library for structural biology | JOSS | 2022 |
