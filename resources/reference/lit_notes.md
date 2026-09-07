# Literature notes — FBN1 / Marfan structural bioinformatics

Phase 0 literature intake. One entry per PDF actually present in `resources/papers/`.
Every citation and every number below was read out of the PDF in this folder; DOIs were
verified against the text of the file itself, not from memory. Where a value the project
wants is *not* in the PDF, that is stated explicitly rather than filled in.

- Text extracted with `pypdf` 6.14.2 on 2026-08-02.
- Cross-check against CLAUDE.md Appendix A is in §3 below and in `results/phase0_report.md`.

---

## 1. FBN1 / Marfan biology and variant interpretation

### Baudhuin LM, Kluge ML, Kotzer KE, Lagerstedt SA (2019)
*Variability in gene-based knowledge impacts variant classification: an analysis of FBN1
missense variants in ClinVar.* Eur J Hum Genet 27:1550–1560. DOI 10.1038/s41431-019-0440-3
→ `FBN1_ClinVar_variability_2019.pdf`

**The single most important rulebook now present in this folder**, and the practical stand-in
for the missing ClinGen FBN1 VCEP paper (see §3). It defines FBN1 "critical residue" categories
operationally and reanalyses 674 FBN1 missense variants from 18 ClinVar submitters. Its critical
residue taxonomy is exactly the Phase 3 annotation schema: (1) creates/destroys a Cys in a
cbEGF-like domain; (2) alters any of the five residues of the calcium-binding consensus
**D/N-x-D/N-E/Q-Xm-D/N-Xn-Y/F** in a cbEGF domain; (3) alters the critical Gly/Ala between Cys2
and Cys3 in a subset of cbEGF domains; (4) alters interdomain packing residues; (5)(6)(7) Cys
created/altered in a TB, EGF, or hybrid2 domain (Cys in hybrid1 was explicitly *not* counted
toward PM1); plus fibrillin FUN N-terminus.

**Numbers to reuse.** Conflicting classifications in 30.7% of multi-submitter missense variants.
451 classifications of 361 critical-residue variants, 80.0% (L)P. Cysteine variants 91.3%
(283/310) (L)P vs non-cysteine critical residues 55.3% (78/141) — this cysteine-vs-Ca-site
asymmetry is a direct prior for our enrichment test. Per-category counts (their Table 1, as
classification calls): cbEGF Cys 231 total (90.0% (L)P, 0 (L)B); cbEGF Ca-binding consensus 110
total (60.0% (L)P, 37.2% VUS, 2.7% (L)B); cbEGF C2–C3 Gly/Ala 21 (42.9% (L)P); interdomain
packing 9 (28.6% (L)P); TB Cys 37 (91.9% (L)P); hybrid2 Cys 9 (100% (L)P); EGF Cys 26 (96.2%
(L)P); FUN 8 (100% (L)P). 31.6% (165/522) of non-critical-residue classifications were judged
overclassified, especially once MAF was considered.

**Tolerated-substitution exceptions we must encode, not treat as pathogenic** (they kept these at
VUS/unknown-effect): Asn→Ser at the *second* D/N of the consensus — p.Asn615Ser, p.Asn1030Ser,
p.Asn1282Ser, p.Asn1489Ser, p.Asn2526Ser, p.Asn2650Ser; and D↔N swaps at a D/N position —
p.Asp1197Asn, p.Asp1240Asn, p.Asn1907Asp. 30.6% of (L)P Ca-binding calls were at the *first* D/N,
which they highlight as under-recognised.

### Zhang M, Chen Z, Chen T, Sun X, Jiang Y
*Cysteine substitution and calcium-binding mutations in FBN1 cbEGF-like domains are associated
with severe ocular involvement in patients with congenital ectopia lentis.*
Front Cell Dev Biol. DOI 10.3389/fcell.2021.816397
→ `FBN1_cbEGF_ectopia_lentis_ocular.pdf`

Retrospective cohort of 68 congenital ectopia lentis probands with FBN1 cbEGF-like mutations,
split by genotype into a **cysteine group (n=43)**, a **Ca²⁺-binding group (n=13)**, and others
(n=12) — the same two-mechanism split this project tests structurally, and evidence that the two
classes are phenotypically distinguishable rather than interchangeable. States the domain rules
we rely on: FBN1 has **43 cbEGF-like domains**; the six conserved cysteines form disulfides in a
**C1–C3, C2–C4, C5–C6** pattern (consistent with CLAUDE.md); the Ca²⁺-binding sequence comprises
an N-terminal loop and a C-terminal β-hairpin and is required for FBN1-to-microfibril assembly.

**Numbers to reuse.** Age-adjusted, both cysteine substitutions and Ca²⁺-binding mutations
contributed to axial length elongation (standardized coefficients 0.410, p=0.008 and 0.367,
p=0.017). Cataract more frequent with Ca²⁺-binding mutations (observed 3 vs expected 1.0,
p=0.036). Cysteine-group patients had the poorest preoperative visual acuity (p=0.012).

### Delhomme C, et al. (2022)
*Genotype–mitral valve phenotype correlations in Marfan syndrome with FBN1 pathogenic variants.*
JACC Adv 1:100149. DOI 10.1016/j.jacadv.2022.100149
→ `Delhomme_2022_JACCAdv_MV_genotype.pdf`

Research letter, 250 MFS patients (2015–2017), the first genotype–phenotype correlation for
mitral valve prolapse (MVP) and mitral annular disjunction (MAD). Its value here is the
**cysteine-content framing**: in-frame variants are split by whether they remove a cysteine
(−Cys), add one (+Cys), or leave cysteine count unchanged (noCys) — a variant annotation we can
compute directly in Phase 3 and a justification for treating Cys-removing and Cys-adding
variants separately rather than lumping "cysteine variants".

**Numbers to reuse.** 99 (39.6%) PTC vs 151 (60.4%) in-frame variants; among in-frame, 27 +Cys,
49 −Cys, 75 noCys. Overall MVP 93/250 (37.2%), similar for PTC vs IFV (37.4% vs 37.1%, P=0.96).
Within IFV, MVP was 51% (25/49) −Cys vs 34.7% (26/75) noCys vs 18.5% (5/27) +Cys, P=0.016; MAD
32.6% vs 19.1% vs 7.4%, P=0.033. Background note: MVP prevalence in MFS reported as high as 40%.

### Li L, Huang J, Liu Y (2023)
*The extracellular matrix glycoprotein fibrillin-1 in health and disease.*
Front Cell Dev Biol 11:1302285. DOI 10.3389/fcell.2023.1302285
→ `Fibrillin1_health_disease_2023_FCELL.pdf`

Broad review; used here to corroborate the domain architecture this pipeline pins in Phase 3.
Confirms fibrillin-1 is a **2871-amino-acid** cysteine-rich, calcium-binding ECM glycoprotein,
and that **of the 47 EGF domains, 43 are calcium-binding (cbEGF)**, alongside **2 hybrid
domains** — matching CLAUDE.md §1 and giving an independent check on the UniProt feature pull.
Also covers FBN1's interactions with microfibril-associated proteins, growth factors and cell
surface receptors, and the disease spectrum beyond Marfan. Note the domain counts here are
review-level; **UniProt P35555 features remain the source of truth** for boundaries.

### Asano K, Cantalupo A, Sedes L, Ramirez F (2022)
*The multiple functions of fibrillin-1 microfibrils in organismal physiology.*
Int J Mol Sci 23:1892. DOI 10.3390/ijms23031892
→ `Asano_2022_IJMS_fibrillin_review.pdf`

Review of microfibril biogenesis and function. Provides the mechanistic context for why loss of
Ca²⁺ binding or a disulfide matters at the tissue level: fibrillin-1 microfibrils both bear load
and control **TGFβ bioavailability**, tethering the latent TGFβ complex (TGFβ pro-peptides bound
to LTBP-1/3/4) to the ECM as a spatially organised, rapidly mobilisable reservoir. Reports the
cardinal observation that MFS patient connective tissue shows a significant decrease in
microfibrils. Useful for the manuscript introduction/discussion; contributes no numeric
thresholds to the analysis.

---

## 1b. Primary structural literature on cbEGF domains

*(These five PDFs were added by the mentor on 2026-08-02 while Phase 0 was running — the gate
caught them and they are noted here. Four are correct; see §3 for the one filename discrepancy.)*

### Handford PA, Downing AK, Rao Z, Hewett DR, Sykes BC, Kielty CM (1995)
*The calcium binding properties and molecular organization of epidermal growth factor-like
domains in human fibrillin-1.* J Biol Chem 270(12):6751–6756. DOI 10.1074/jbc.270.12.6751
→ `Handford_1995_JBC_cbEGF.pdf`

**This PDF is a scanned page image with no extractable text layer** — it was read visually
(embedded TIFF → PNG) rather than parsed, so notes here come from the printed page.
The foundational measurement behind this whole project: it establishes that fibrillin-1 cbEGF
domains actually bind Ca²⁺ and that consensus mutations reduce that binding. Confirms the
architecture we pin in Phase 3 — a **350-kDa** glycoprotein of 10-nm microfibrils with **47
EGF-like domains, 43 of which carry the calcium-binding consensus**, plus **7 TGFβ1-binding
(8-cysteine) domains**, **2 hybrid domains**, a proline-rich region, and 2 unique terminal regions.

**Consensus sequence as originally written:** **Asp–Asp/Asn–Glu/Gln–Asp\*/Asn\*–Tyr/Phe**, where
`*` marks a **β-hydroxylated** residue. (Compatible with Baudhuin 2019's
D/N-x-D/N-E/Q-Xm-D/N-Xn-Y/F; the two differ in notation, not content.)

**Numbers to reuse.** A synthesized wild-type cbEGF (residues 2126–2165) vs the Marfan mutant
**Asn-2144→Ser** showed **>5-fold reduction in calcium affinity** by ¹H-NMR, while rotary
shadowing showed the mutation does **not** prevent microfibril assembly but alters the interbead
region — i.e. a Ca-site variant can be pathogenic without abolishing assembly. Marfan incidence
quoted as at least **1 in 10,000**; microfibril beaded periodicity **50–55 nm**.

**⚠️ Numbering trap:** residue numbers in this paper follow **Pereira et al. (1993)**, not
UniProt P35555. Do **not** copy `Asn-2144` (or any 1990s-era position) into our tables without
re-mapping and re-checking WT identity against P35555. Same caution applies to Rao 1995 and
Downing 1996 below.

### Rao Z, Handford P, Mayhew M, Knott V, Brownlee GG, Stuart D (1995)
*The structure of a Ca²⁺-binding epidermal growth factor-like domain: its role in protein-protein
interactions.* Cell 82:131–141. DOI 10.1016/0092-8674(95)90059-4
→ `Rao_1995_Cell_cbEGF_structure.pdf`

The original crystal structure of a cbEGF domain (from human clotting factor IX) — the structural
basis for treating the five consensus residues as a genuine coordination site rather than an
alignment artifact. States the consensus as **Asp/Asn, Asp/Asn, Gln/Glu, Asp\*/Asn\*, Tyr/Phe**,
with `*` a possible β-hydroxylated residue, and notes that the importance of each consensus
residue was established by extensive mutagenesis on the factor IX cbEGF domain. Establishes the
disease logic we extend to FBN1: mutations at consensus residues cause hemophilia B in factor IX,
just as they cause Marfan syndrome in fibrillin-1. Useful for the Phase 4 rationale and the
manuscript's methods framing of "calcium-coordinating residue".

### Downing AK, Knott V, Werner JM, Cardy CM, Campbell ID, Handford PA (1996)
*Solution structure of a pair of calcium-binding epidermal growth factor-like domains:
implications for the Marfan syndrome and other genetic disorders.* Cell 85:597–605.
DOI 10.1016/S0092-8674(00)81259-3
→ `Downing_1996_Cell_cbEGF_pair.pdf`

The **fib32–33 (cbEGF32-33) pair** — the structure behind PDB 1EMN/1EMO and the reason this
project analyses domain *pairs*, not isolated domains. Shows the two domains adopt a **rigid,
rod-like arrangement stabilized by interdomain calcium binding**, packing at roughly **130°** in
the crystal asymmetric unit, with a calcium ligand **donated from one domain to the other** and
the interface further stabilized by interdomain hydrophobic packing. Reports a difference in
calcium affinity between the N-terminal and C-terminal sites, with the N-terminal site's affinity
enhanced by the interdomain interface.

**Direct consequence for our analysis:** because a Ca²⁺ site can be completed by residues from
the *neighbouring* domain, a variant's coordination environment is not necessarily contained
within its own domain. Phase 4 should model **domain pairs** where the site is interdomain — which
is what CLAUDE.md already implies with "per cbEGF domain (or domain pair)". Residues at the
interdomain packing interface are also a critical-residue category in Baudhuin 2019.

### Jensen SA, Robertson IB, Handford PA (2012)
*Dissecting the fibrillin microfibril: structural insights into organization and function.*
Structure 20(2):215–225. DOI 10.1016/j.str.2011.12.008
→ `Jensen_2009_Structure_2W86.pdf` — **filename says 2009; the file is the 2012 review** (§3)

A review of fibrillin microfibril structure by the group that solved most of the FBN1 fragments.
Its single most useful contribution here is an **authoritative PDB ↔ domain mapping**, which
independently confirms the Phase 4 structure inventory and resolves the entry CLAUDE.md
Appendix B left unnamed:

| Fragment | PDB |
|---|---|
| cbEGF9–hyb2–cbEGF10 | **2W86** |
| cbEGF12–cbEGF13 | **1LMJ** |
| cbEGF22–TB4–cbEGF23 | **1UZJ** |
| isolated TB6 | **1APJ** |
| cbEGF32–cbEGF33 | **1EMN** |

Also defines the **"neonatal region" as TB3 → cbEGF18**, where substitutions causing the severe
neonatal form of Marfan syndrome cluster — a region flag worth carrying in the Phase 3 annotation
table and a natural stratum for the enrichment analysis. Notes that pairwise domain interactions
in these structures suggest a near-linear organization stabilized by calcium, and that calcium
binding rigidifies the domain linkage.

*Caveat:* this is a secondary source (a review). The primary 2W86 paper named in Appendix A
(Jensen et al. 2009, Structure 17:759) is still not present.

---

## 2. Methods, databases and tools

### Landrum MJ, et al. (2018)
*ClinVar: improving access to variant interpretations and supporting evidence.*
Nucleic Acids Res 46:D1062–D1067. DOI 10.1093/nar/gkx1153
→ `ClinVar_2018_NAR.pdf`

The citation for the Phase 1 variant source and, more usefully, for the **review-status concept
that drives our star cutoff**. Review status is derived from the submissions behind a variant —
how many submitters agree, whether assertion criteria were provided, and whether an **expert
panel or practice guideline** submitted the interpretation. Submissions are SCV records
aggregated into a variant-level interpretation. Confirms the design decision in CLAUDE.md to
record review status per variant and prefer ≥2-star/expert-panel while labelling rather than
dropping conflicts. The paper describes filters for higher review statuses but does not
enumerate the star mapping numerically — take the star definitions from ClinVar's own
documentation at pull time and record the accessed date.

### The UniProt Consortium (2023)
*UniProt: the Universal Protein Knowledgebase in 2023.*
Nucleic Acids Res 51:D523–D531. DOI 10.1093/nar/gkac1052
→ `UniProt_2023_NAR.pdf`

Citation for P35555 as the source of truth for sequence, domain features (EGF/cbEGF/TB/hybrid),
disulfide bonds and calcium-binding annotations used in Phase 3. Documents the Swiss-Prot
(manually reviewed) vs TrEMBL distinction that justifies treating the P35555 feature table as
curated evidence. Record the UniProt release used at pull time — this Phase 0 run retrieved
P35555 at length 2871 with sequence MD5 `f8f84ac9dc30f3ff640ad4f691fee2b1`.

### Jumper J, et al. (2021)
*Highly accurate protein structure prediction with AlphaFold.* Nature 596:583–589.
DOI 10.1038/s41586-021-03819-2
→ `Jumper_2021_AlphaFold_Nature.pdf`

The AlphaFold2 method citation. Relevant caveat for this project: AF2 predicts **protein
coordinates only** — the method has no notion of bound metals or ligands, which is the formal
basis for CLAUDE.md's rule that a raw AlphaFold model cannot be handed to CheckMyMetal. Cite for
methods; no reusable numeric thresholds beyond pLDDT, whose interpretation bands are given in the
AlphaFold DB paper below.

### Varadi M, et al. (2024)
*AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million
protein sequences.* Nucleic Acids Res 52:D368–D375. DOI 10.1093/nar/gkad1011
→ `AlphaFoldDB_2024_NAR.pdf`

**Directly explains why FBN1 has no AlphaFold DB model** — a Phase 0 finding (see
`results/phase0_report.md`). The paper states which sequences are *not* covered: those under 16
amino acids, or **>2700 residues for Swiss-Prot/proteome sequences** (1280 for other UniProt
sequences), those with non-standard amino acids, those absent from the UniProt "one sequence per
gene" FASTA, and viral proteins. **P35555 is 2871 aa and therefore falls outside the >2700
cutoff**, which is consistent with the live API returning 404 for P35555 while returning 200 for
a 2477-aa control, and with UniProt carrying no AlphaFoldDB cross-reference for P35555.

**Numbers to reuse — pLDDT interpretation bands** (stored in the B-factor field of AFDB PDB
files): >90 modelled to high accuracy; 70–90 generally well modelled; 50–70 low confidence;
<50 often "spaghetti-like" / likely disordered. PAE is provided separately for inter-domain
positioning.

### Hekkelman ML, de Vries I, Joosten RP, Perrakis A (2023)
*AlphaFill: enriching AlphaFold models with ligands and cofactors.* Nat Methods 20:205–213.
DOI 10.1038/s41592-022-01685-y
→ `AlphaFill_2023_NatMethods.pdf`

The Phase 4 Tier-3 cross-check method. States the problem this project must solve explicitly —
AlphaFold models "all lack coordinates for small molecules", with **zinc-finger motifs lacking
the zinc ions essential for structural integrity** given as a direct analogue of our Ca²⁺ case.
AlphaFill transplants ligands/ions from experimentally determined homologous structures.

**Method parameters to reuse.** Transplants are made from structures with **sequence identity
>25% over an aligned stretch of at least 85 residues**; identity levels of 25/30/50/70% are
reported as indicative tiers; transplant quality is reported as a **local RMSD** of the
environment around the transplanted compound, which is the metric to record in our
`structures/ca_transplanted/` provenance JSON. Scale: 12,029,789 transplants across 995,411
AlphaFold models, covering 2,694 distinct compounds. **Caveat for us:** AlphaFill is built on
AlphaFold DB entries, so with no AFDB model for P35555 there is no AlphaFill entry for FBN1
(confirmed live: status `unknown`); a Tier-3 cross-check would have to be run on a per-domain
model rather than fetched.

### Zheng H, Cooper DR, Porebski PJ, Shabalin IG, Handing KB, Minor W (2017)
*CheckMyMetal: a macromolecular metal-binding validation tool.* Acta Cryst D73:223–233.
DOI 10.1107/S2059798317001061
→ `CheckMyMetal_2017_ActaD.pdf`

The Phase 4.5 hand-off tool and the Phase 5 parsing spec. Server URL confirmed in the paper as
**http://csgid.org/csgid/metal_sites** (CLAUDE.md gives the https form). Motivation: metals are
modelled in roughly 40% of PDB macromolecular structures and a significant fraction of those
sites are poorly modelled — i.e. this validation step is not a formality. CMM reports **eight
parameters**, each flagged as acceptable / borderline / outlier, including ligand identity,
coordination number, geometry class, **gRMSD**, **overall valence** (bond-valence model),
**nVECSUM**, **vacancy**, and B-factor comparison.

**Thresholds to reuse when parsing `handoff/cmm/outbox/`** (borderline / outlier): nVECSUM
**>0.10 / >0.23**; gRMSD **>13.5° / >21.5°**; vacancy **>10% / >25%**. The valence parameter is
judged against the expected formal charge for the metal — **+2 for Ca²⁺** — with borderline and
outlier zones set relative to the tabulated range. The B-factor parameter compares the metal's B
factor with the bond-valence-weighted mean B of its ligand atoms.

### Cheng J, Novati G, Pan J, et al. (2023)
*Accurate proteome-wide missense variant effect prediction with AlphaMissense.*
Science 381:eadg7492. DOI 10.1126/science.adg7492
→ `AlphaMissense_2023_Science.pdf` — **present, though Appendix A expected it missing**

Citation for the orthogonal pathogenicity comparator loaded in Phase 1 and joined in Phase 2.
Classifies variants into three ACMG-like categories — **likely pathogenic, likely benign, and
ambiguous** — with cutoffs chosen so that both confident classes carry **90% expected precision
estimated from ClinVar**. Reports classifying 92.9% of ClinVar test variants at 90% precision vs
67.1% for EVE (+25.8 percentage points). Proteome-wide, ~32% of missense variants were called
likely pathogenic and 57% likely benign at that cutoff.

**Important for us:** the PDF in this folder is the Research Article Summary and **does not state
the numeric score cutoffs** (the commonly quoted 0.34 / 0.564 boundaries are not in this text).
We therefore use the `am_class` column shipped with the data rather than re-deriving classes from
`am_pathogenicity`, and we do not cite numeric cutoffs we cannot source. License is
CC BY-NC-SA 4.0 (academic use), Zenodo record 8208688.

### Chen S, Francioli LC, Goodrich JK, et al. (2024)
*A genomic mutational constraint map using variation in 76,156 human genomes.*
Nature 625:92–100. DOI 10.1038/s41586-023-06045-0
→ `gnomAD_v4_constraint_2024_Nature.pdf`

The gnomAD constraint citation for the benign-set work in Phase 2. Aggregates **76,156 human
genomes** and builds a genome-wide constraint map (**Gnocchi**) that extends constraint modelling
beyond protein-coding regions, where classical gene-level metrics do not apply. Relevant to us as
the methodological citation for using gnomAD as a population reference and for the principle that
depletion of variation marks functionally important sequence.

**Caveat for our use:** this paper is about **genome-wide constraint**, not the gnomAD v4 variant
frequency release itself. Our benign threshold depends on the **Grpmax filtering allele frequency
from gnomAD v4.1.0**, which comes from the gnomAD API/release notes, not from this PDF. Record
the dataset version at pull time and do not cite this paper for an AF cutoff it does not state.

---

## 3. Cross-check against CLAUDE.md Appendix A

**19 PDFs** are now present in `resources/papers/`. The mentor added the five Appendix-A
"missing" papers on 2026-08-02 *during* this Phase 0 run; the Phase 0 gate detected them and they
have been read and noted above.

Seventeen match their filenames. **Two do not** — the files exist but contain entirely different
papers, verified by full-text search returning **zero** occurrences of FBN1, fibrillin, Marfan or
cbEGF in either. **These two were not fixed by the 2026-08-02 top-up and are still wrong.**
One further file is the right topic and authors but the wrong year/article.

| File | Appendix A expects | PDF actually contains | Status |
|---|---|---|---|
| `ClinGen_FBN1_VCEP_2024_GenomeMed.pdf` | ClinGen FBN1 VCEP (2024), Genome Med 16:141, DOI 10.1186/s13073-024-01423-3 — the **critical-residue rulebook** | A PNAS 2024 paper on amino acids modulating intracellular liquid–liquid phase separation and stress granules, DOI **10.1073/pnas.2407633121** | **WRONG PAPER** |
| `Godwin_2023_NSMB_microfibril_cryoEM.pdf` | Godwin A, et al. (2023) fibrillin microfibril cryo-EM, Nat Struct Mol Biol 30:608, DOI 10.1038/s41594-023-00950-8 | Jia BB, Jussila A, Kern C, Zhu Q, Ren B. *A spatial genome aligner for resolving chromatin architectures from multiplexed DNA FISH.* Nat Biotechnol 41:1004–1017 (2023), DOI **10.1038/s41587-022-01568-9** | **WRONG PAPER** |
| `Jensen_2009_Structure_2W86.pdf` | Jensen SA, et al. (2009) cbEGF9-hyb2-cbEGF10 primary paper, Structure 17:759 | Jensen SA, Robertson IB, Handford PA (2012) *Dissecting the fibrillin microfibril*, Structure 20:215–225, DOI **10.1016/j.str.2011.12.008** — a **review**, same group, right topic, wrong article/year | **PARTIAL — usable** |

**Present and correct (16):** `FBN1_ClinVar_variability_2019`, `FBN1_cbEGF_ectopia_lentis_ocular`,
`Delhomme_2022_JACCAdv_MV_genotype`, `Fibrillin1_health_disease_2023_FCELL`,
`Asano_2022_IJMS_fibrillin_review`, `ClinVar_2018_NAR`, `UniProt_2023_NAR`,
`Jumper_2021_AlphaFold_Nature`, `AlphaFoldDB_2024_NAR`, `AlphaFill_2023_NatMethods`,
`CheckMyMetal_2017_ActaD`, `AlphaMissense_2023_Science`, `Rao_1995_Cell_cbEGF_structure`,
`Downing_1996_Cell_cbEGF_pair`, `Handford_1995_JBC_cbEGF`, `gnomAD_v4_constraint_2024_Nature`.

**Resolved on 2026-08-02** (were Appendix-A gaps, now present and read): Rao 1995, Downing 1996,
Handford 1995, gnomAD v4 constraint 2024, plus AlphaMissense 2023 (added earlier).

**Still missing / unusable (3) — recorded as gaps, contents NOT inferred or invented:**

1. **ClinGen FBN1 VCEP (2024)** — Genome Med 16:141, DOI 10.1186/s13073-024-01423-3. The
   nominally present file is the wrong paper. This is the one CLAUDE.md flags as the
   critical-residue rulebook, and Phase 3 requires aligning our "critical residue" definition to
   it. *Mitigation:* Baudhuin 2019 supplies an operational FBN1 critical-residue taxonomy, but it
   is **not** the VCEP specification and its PM1/PP3 calibration will differ.
2. **Godwin A, et al. (2023)** — fibrillin microfibril cryo-EM, NSMB 30:608,
   DOI 10.1038/s41594-023-00950-8. The nominally present file is the wrong paper. *Mitigation:*
   Jensen 2012 (above) covers microfibril organization, though it predates the cryo-EM structure.
3. **Jensen SA, et al. (2009)** — the primary cbEGF9-hyb2-cbEGF10 / PDB 2W86 paper,
   Structure 17:759. The file under this name is the 2012 review. *Mitigation:* the 2012 review
   supplies the PDB↔domain mapping we needed; only the primary structural citation is missing.

**Extraction note:** `Handford_1995_JBC_cbEGF.pdf` is a **scanned image PDF with no text layer**.
It cannot be parsed with `pypdf` and was read visually. Any future automated pass over
`resources/papers/` must not treat its empty text extraction as an empty paper.

---

## Added 2026-08-16 — two of the gaps above are now closed

The mentor supplied two previously paywalled papers. Both were read from the PDFs now in
`resources/papers/`; gap items 3 and (partly) the N2144S anchor are resolved.

### Kettle S, Yuan X, Grundy G, Knott V, Downing AK, Handford PA (1999)
*Defective calcium binding to fibrillin-1: consequence of an N2144S change for fibrillin-1
structure and function.* **J Mol Biol** 285:1277–1291.
File: `kettle1999.pdf` (11 pp).

**This is the most directly load-bearing paper in the whole collection**, because it is the one
experimental result this project's central claim can be checked against.

NMR study of wild-type and N2144S TB6–cbEGF32 and cbEGF32–33 domain pairs. N2144S removes one
of the calcium-binding ligands of cbEGF32. Findings we reuse:

1. **The substitution does not alter the native fold** of TB6, cbEGF32 or cbEGF33. Stated
   explicitly in the abstract and supported by the NMR assignments.
2. **Calcium affinity falls ~9-fold** in the TB6–cbEGF32 pair: Kd 1.6 mM wild-type → ~14 mM
   mutant. Simulations varying the assumed protein concentration gave 1.8–1.1 mM (wild-type) and
   14.2–13.6 mM (mutant), so the ninefold figure is robust to that uncertainty. The isolated
   cbEGF32 domain shows a ~5-fold reduction, which is where the "5–9×" range in our Figure 4
   annotation comes from.
3. **The defect is local.** Reduced cbEGF32 affinity does not lower cbEGF33 affinity, and the
   authors attribute the effect to the flexible TB6–cbEGF32 interdomain linkage.
4. **Biosynthesis is unaffected**: ³⁵S pulse-chase in N2144S patient fibroblasts showed no
   detectable change in fibrillin-1 synthesis, secretion rate, processing, or deposition of
   reducible fibrillin-1 into the matrix.
5. Context Kd values reused for scale: isolated cbEGF32 ~4.3 mM; cbEGF33 in the cbEGF32–33 pair
   ~350 µM; free ionised calcium physiologically 1.0–1.2 mM, so the wild-type TB6–cbEGF32 site is
   only **partially saturated in vivo**.

**Why this matters for us:** points 1 and 2 together are precisely the thesis of this study —
fold intact, calcium binding broken — demonstrated experimentally for one variant, twenty-seven
years before we modelled it. It is the anchor for Figure 4 and the strongest single justification
for trusting the AF3 mutant geometry on variants nobody has measured.

**Unit caution:** the PDF's text layer mangles µ/m and renders `fi`/`fl` ligatures as `®`/`¯`.
Quantities above were reconciled against internal consistency (a 1.6 mM Kd against 1.0–1.2 mM
free calcium giving partial saturation) rather than taken from the extraction verbatim. Check any
further number against the rendered page, not the text layer.

### Jensen SA, Iqbal S, Lowe ED, Redfield C, Handford PA (2009)
*Structure and interdomain interactions of a hybrid domain: a disulphide-rich module of the
fibrillin/LTBP superfamily of matrix proteins.* **Structure** 17:759–768.
DOI 10.1016/j.str.2009.03.014. File: `1-s2.0-S0969212609001610-main.pdf` (10 pp).

The primary paper for **PDB 2W86**, closing gap item 3 above. Crystal structure of the
fibrillin-1 cbEGF9–hyb2–cbEGF10 fragment at **1.8 Å**, crystallised in **20 mM CaCl₂** to
saturate both cbEGF calcium sites — which is why we treat it as a gold-standard Ca²⁺ donor and an
AF3 calibration target rather than a Ca-free model.

Reusable specifics:

- **Ca²⁺ coordination, cbEGF9**: three side-chain oxygens (Asp807 OD1, Glu810 OE1, Asn823 OD1),
  three main-chain carbonyls (Ile808, Ser824, Ser827) and one water — approximately pentagonal
  bipyramidal.
- **Ca²⁺ coordination, cbEGF10**: side chains Asp910, Glu913, Asn928; main-chain carbonyls
  Ile911, Thr929, Ser932; one water.
- These are exactly the donor sets our coordination-number and bond-valence code should recover
  on 2W86, and they confirm the **three side-chain + three backbone + one water** pattern that
  motivates separating side-chain from backbone-only ligands in the v2 grouping.
- The hyb2 domain adopts a TB-like fold with rearranged disulphides (C1-3, 2-5, 4-6, 7-8 in hyb
  versus C1-3, 2-6, 4-7, 5-8 in TB), and the **hyb2–cbEGF10 interface is Ca²⁺-dependent**.
- Superposition onto cbEGF22–TB4–cbEGF23 gave 2.96 Å overall, with individual domains at
  0.829 / 1.844 / 1.668 Å.

**Why this matters for us:** it supplies the experimentally determined ligand identity for the
domains AF3 was calibrated on, and its explicit backbone-carbonyl contribution is the literature
basis for the "backbone-only ligand" class, which the enrichment analysis finds *depleted* of
pathogenic variants (OR 0.049) — as expected, since no substitution can remove a backbone
carbonyl.

**Remaining gaps after this update:** the ClinGen FBN1 VCEP 2024 paper and Godwin 2023 are still
the wrong files under the right names (items 1 and 2 above). Whiteman & Handford 2000, Reinhardt
2000 are still absent and would add further measured calcium affinities.

---

## Added 2026-08-30 — two "missing" papers were never missing, and 12 references added from the literature

### Correction: `ClinGen_FBN1_VCEP_2024_GenomeMed.pdf` and `Godwin_2023_NSMB_microfibril_cryoEM.pdf` are correct

§3 above records both as **WRONG PAPER**. They are not, and were probably replaced by the mentor
after that note was written. Verified 2026-08-30 by reading the PDFs:

- `ClinGen_FBN1_VCEP_2024_GenomeMed.pdf` → Drackley A, Somerville C, Arnaud P, *et al.*
  *Interpretation and classification of FBN1 variants associated with Marfan syndrome: consensus
  recommendations from the Clinical Genome Resource's FBN1 variant curation expert panel.*
  Genome Medicine 16:154 (2024). DOI 10.1186/s13073-024-01423-3.
  **Note the article number is 154, not 141 as CLAUDE.md Appendix A states.**
- `Godwin_2023_NSMB_microfibril_cryoEM.pdf` → Godwin ARF, *et al.* *Fibrillin microfibril structure
  identifies long-range effects of inherited pathogenic mutations affecting a key regulatory latent
  TGFβ-binding site.* Nat Struct Mol Biol 30:608–618 (2023). DOI 10.1038/s41594-023-00950-8.

**What the VCEP actually specifies** — this is the rulebook the manuscript argues with, so the
rules are recorded verbatim in substance:

- **PM1_Strong** for cysteine-removing variants in any of the 43 cbEGF domains (258 cysteines).
- **PM1 (moderate)** for cysteine substitutions outside cbEGF domains, cysteine-*creating* variants
  in cbEGF domains, interdomain-packing glycines, β-hydroxylation sites, and variants altering the
  conserved residues of the consensus calcium-binding sequence
  **[D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F]**.
- **Exemption:** Asn→Ser at the *second* [D/N] of that consensus "might be tolerated based on the
  frequency of this type of missense variant in gnomAD", and PM1 should not be applied there.
- PM1 applies at moderate or strong strength to 22.0% (633/2871) of positions, across 375 amino
  acid positions for the non-cysteine categories.
- **PP3 uses REVEL ≥ 0.75**; BP4 uses REVEL ≤ 0.326. Not AlphaMissense.

Three consequences for this project, all now in the manuscript. The two features are **not** given
equal weight, so the earlier framing was wrong. The calcium rule keys on the **sequence consensus**,
which our models show over-calls by 143 positions. And the Asn→Ser exemption falls exactly where
our donor-chemistry rule breaks.

### References added on 2026-08-30, verified against Europe PMC or the publisher record

Not present as PDFs in this folder; cited from the verified bibliographic record only, and used
for claims that do not depend on reading the full text beyond the abstract-level facts noted.

| Citation | Why it is cited |
|---|---|
| Sakai LY, Keene DR, Engvall E. *J Cell Biol* 103:2499 (1986). doi:10.1083/jcb.103.6.2499 | the original identification of fibrillin as a 350-kD microfibril component |
| Dietz HC, *et al.* *Nature* 352:337 (1991). doi:10.1038/352337a0 | the first FBN1 missense mutation shown to cause Marfan syndrome |
| **Handford PA, *et al.* *Nature* 351:164 (1991). doi:10.1038/351164a0** | **the precedent for the whole thesis**: in factor IX, removing a coordinating Asp reduced Ca²⁺ affinity >1000-fold and caused haemophilia B; Asp47Glu reduced it >4-fold |
| Knott V, Downing AK, Cardy CM, Handford PA. *J Mol Biol* 255:22 (1996). doi:10.1006/jmbi.1996.0003 | calcium affinities of the cbEGF32–33 pair |
| **Reinhardt DP, *et al.* *J Biol Chem* 275:12339 (2000). doi:10.1074/jbc.275.16.12339** | **a second experimental anchor**: N548I and E1073K — both in our cofolding set — render fibrillin-1 susceptible to proteolysis, with cleavage sites mapping close to the mutation, i.e. local structural change. Our models give N548I −33.0% coordination and E1073K −11.9% |
| McGettrick AJ, Knott V, Willis A, Handford PA. *Hum Mol Genet* 9:1987 (2000). doi:10.1093/hmg/9.13.1987 | calcium-ligand mutations behave differently depending on whether the preceding domain is a TB or a cbEGF. **Note: CLAUDE.md and earlier notes attribute this to "Whiteman & Handford" — the first author is McGettrick.** |
| Whiteman P, Handford PA. *Hum Mol Genet* 12:727 (2003). doi:10.1093/hmg/ddg081 | some cbEGF substitutions impair secretion rather than folding |
| Whiteman P, *et al.* *Hum Mol Genet* 16:907 (2007). doi:10.1093/hmg/ddm035 | co-operative folding across the cbEGF12–13 region |
| Smallridge RS, *et al.* *J Biol Chem* 278:12199 (2003). doi:10.1074/jbc.M208266200 | the cbEGF12–13 (neonatal region) solution structure, PDB 1LMJ |
| Lee SS, *et al.* *Structure* 12:717 (2004). doi:10.1016/j.str.2004.02.023 | the cbEGF22–TB4–cbEGF23 fragment, PDB 1UZJ |
| **Haller SJ, Roitberg AE, Dudley AT. *Sci Rep* 10:16844 (2020). doi:10.1038/s41598-020-73969-2** | **the closest prior computational work.** Steered MD on cbEGF12–13 and TB4–cbEGF23 with four Cys→Arg mutations (C1086R, C1117R, C1111R, C1138R): cbEGF calcium binding weakens under mechanical strain, and disulfide mutations disrupt that response. Complementary to ours — they ask how cysteine loss affects calcium; we ask how calcium-ligand loss affects stability |
| Loeys BL, *et al.* *J Med Genet* 47:476 (2010). doi:10.1136/jmg.2009.072785 | the revised Ghent nosology |
| Milewicz DM, *et al.* *Nat Rev Dis Primers* 7:64 (2021). doi:10.1038/s41572-021-00298-7 | current clinical review of Marfan syndrome |
| Richards S, *et al.* *Genet Med* 17:405 (2015). doi:10.1038/gim.2015.30 | the ACMG/AMP framework the VCEP specialises |
| Pejaver V, *et al.* *Am J Hum Genet* 109:2163 (2022). doi:10.1016/j.ajhg.2022.10.013 | calibration of computational evidence to PP3/BP4 strengths |
| Abramson J, *et al.* *Nature* 630:493 (2024). doi:10.1038/s41586-024-07487-w | AlphaFold 3, the model behind the cofolding used here — **this was the reference the earlier draft was missing** |
| Alford RF, *et al.* *J Chem Theory Comput* 13:3031 (2017). doi:10.1021/acs.jctc.7b00125 | the Rosetta ref2015 energy function |
| Schymkowitz J, *et al.* *Nucleic Acids Res* 33:W382 (2005). doi:10.1093/nar/gki387 | FoldX |
| Brown ID, Altermatt D. *Acta Crystallogr B* 41:244 (1985). doi:10.1107/S0108768185002063 | the bond-valence parameters used for every coordination measurement here |
| Kabsch W, Sander C. *Biopolymers* 22:2577 (1983). doi:10.1002/bip.360221211 | DSSP |
| Mitternacht S. *F1000Res* 5:189 (2016). doi:10.12688/f1000research.7931.1 | FreeSASA |
| Dana JM, *et al.* *Nucleic Acids Res* 47:D482 (2019). doi:10.1093/nar/gky1114 | SIFTS, the UniProt↔PDB residue mapping |
| Karczewski KJ, *et al.* *Nature* 581:434 (2020). doi:10.1038/s41586-020-2308-7 | gnomAD constraint |
| Vallat B, *et al.* *Nucleic Acids Res* 54:D489 (2026). doi:10.1093/nar/gkaf1187 | the current RCSB PDB reference |

**Still not present and still not fabricated:** no full text was consulted for the rows above
beyond publisher abstracts and, where quoted, the passages reproduced in this file. Claims in the
manuscript that require detailed content — Kettle 1999's Kd values, Jensen 2009's ligand list,
Reinhardt 2000's mutation set, Haller 2020's simulated mutations, the VCEP's PM1 wording — were
each taken from a source read in full or from a verified abstract, and are marked as such above.
