# Session context — resume point

## 0. What changed on 2026-09-06 — READ THIS FIRST

**The mentor returned a submission review of `manuscript_v4.md` and every item is now addressed
except citations, which were explicitly deferred.** The draft is still citation-free and still
built by `61_typeset_v3.py`. Target venue is now **PLOS ONE**, which fixes the abstract format.

**The title and every Results heading are neutral noun phrases**, which is the mentor's second
call of the day and supersedes the declarative set that replaced the first one. The title is now
*Distinct mechanisms of calcium-ligand and cysteine substitutions in the fibrillin-1 cbEGF
module*. `style_reference.md` rule 5 records the choice and `test_manuscript_v4.py` **V9** fails
if a heading picks up a finite verb again. Note that `61_typeset_v3.PLACEMENT_V3` hangs each
figure off a Results heading, so it has to move whenever the headings do. British spelling
throughout, and the one American spelling in the body (`sulfur`) is now `sulphur`.

**THE REAL FIX WAS SENTENCE SHAPE, NOT SENTENCE LENGTH. Read `style_reference.md` §3b before
editing prose.** The mentor read a draft in which every `63_style_extract.py` aggregate sat inside
the corpus range and said the multi-clause sentences still read as machine-written. That was
correct. `scripts/65_sentence_shape.py` was written to measure what 63 cannot see, which is how a
long sentence is *built*, and it found the draft welding a second independent clause on with
", and the model is ..." in **28.6%** of sentences against **4.0%** in the five papers, and
**48.2%** of long sentences against 7.5%. The median long sentence carried **0** parentheses
against the corpus 1. A published long sentence has one main clause and takes its length from
noun phrases, parentheses and a trailing participial phrase.

**121 sentences were reshaped over four passes**, using a relative clause, a subordinator, a
participial or absolute phrase, or a parenthesis in place of every weld. All eight shape measures
now sit at the corpus rate or below. **Passes one and three overshot in opposite directions** —
pass one split so much that the median fell to 17 words against the corpus 23 and 37% of sentences
were under 15 words, which reads as choppy and is its own tell. Do not fix welding by splitting.
**Then the mentor found a four-clause sentence in the abstract, with every one of those checks
passing.** That was the real lesson of the day and it is now the first rule in
`style_reference.md` §3b. **An average cannot see one bad sentence.** Counting clause markers one
sentence at a time, with the same instrument on both sides, gives a ceiling instead of a mean, and
the ceiling is absolute — of 1,164 published sentences **79.0% carry one marker, 18.3% two, 2.7%
three and NOT ONE carries four**. Three of ours did. Eleven more sentences were fixed, the heaviest
sentence in the paper is now 3 markers against the corpus maximum of 3, and four-marker sentences
are at 0.0%.

`test_manuscript_v4.py` **V10** now enforces the ceiling as well as the distribution — no sentence
heavier than the heaviest published one, three-marker sentences at or below 3.5%, welding at or
below the corpus rate, and a parenthesis in the median long sentence. Long sentences are still 8%
against the corpus 21%, which is the deliberate short bias inherited from 2026-09-05.

**Do not report a prose pass as finished on aggregate numbers alone.** Read
`results/sentence_shape.md`, which now lists every sentence at three or more markers by name.

**Nor by trading one construction for another that scores well and reads no better.** The first
replacement for that abstract sentence stacked two participles off an absolute — "substitutions
that removed the contacting oxygen atom breaking the site and those supplying another oxygen atom
largely sparing it". Two clause markers, and still three things to hold. It now reads as a plain
`whereas` contrast at three markers, which is what the corpus does.

**A sentence that announces a result instead of stating it is the third prose fault found today,
and it is the one the installed skill actually names.** `academic-paper-composer`'s checklist is
"specific, not vague". The mentor caught "The outcome was decided by the replacing residue." in
the abstract, immediately followed by the sentence saying which substitutions broke the site.
Measured against the corpus, that construction appears **twice in 1,164 published sentences
(0.2%)** and both are real claims. This draft had **nine (3.6%)**. Four were deleted, two
rewritten to say the thing, and two kept because they are claims and not throat-clears — see
`style_reference.md` §3c for which and why. `test_manuscript_v4.py` **V12** now fails if more
than two survive.

**SPELLING IS AMERICAN FROM 2026-09-06 and this reverses the earlier decision.** 75 substitutions
across 25 words in the manuscript, plus one figure panel title (`figS1` "Modelled versus
experimental"). The trap is that four of the five reference papers are British journals — the
corpus writes `disulphide` 31 times and `disulfide` 3 — so any rewrite made by reading them pulls
the spelling straight back. `test_manuscript_v4.py` **V11** checks the draft and every string a
figure script draws. `grey` inside the figure scripts is a colour-key name, `P.C["grey"]`, never
drawn, and is deliberately left alone.

**Both older prose diagnostics were also run, and all seven corpus tells are now at or
below the published rate** — including three that were not, and are the cleanest they have
measured. `63_style_extract.py` caught a 3+ clause overshoot at 17% against the corpus 12% after
the first editing pass, fixed to 15% by splitting four sentences. `62_prose_diagnostic.py`, which
had never been run on this draft, then caught a hedge stack, an "is not X, it is Y" and a
sentence-final ", which is …" — all three now 0.0%.

**The two diagnostics disagree with each other and that is not a bug.**
`62_prose_diagnostic.py` counts a term as defined if "which is", "that is" or "i.e." appears
within 300 characters of first use, and `63_style_extract.py` penalises exactly those
constructions because the five reference papers barely use them. The jargon count has read
**16 of 18 terms undefined** since 2026-09-05 and did not move in either direction this session.
Do not "fix" it by scattering `, which is …` glosses, which is what the checker literally rewards
and what the corpus check would then flag. If the jargon count is to be addressed, it needs
appositives and separate sentences, and it is a mentor decision because it changes the register.

**The 1,038 / 1,132 discrepancy was real and is now explained rather than patched.** The cysteine
group is every cysteine substitution in *FBN1* (n = 1,132, of which 851 are in a cbEGF domain and
281 are in TB or non-calcium EGF domains, so the old wording "cbEGF cysteine" was wrong and is
gone). A ΔΔG needs a modelled construct and 94 of them lie outside the span of the 21 constructs,
which is why Fig. 2a and 2c carry 1,038 and Fig. 2d carries 1,132. The figure caption says so, the
figure sidecar records the 94, and `test_v2.py` F2 checks it.

**Every headline panel now prints its Benjamini-Hochberg *q* beside the raw *P*.**
`paperstyle.stat/pq/fmt_pq/agrees` read `results/statistics_final.tsv` so a figure can never print
a *q* from a different test — `agrees()` stops the build if a panel's own recomputed *P* has
drifted from the table row whose *q* it is about to print. One status changes under correction and
is now stated in the text, FoldX at the calcium ligands (*P* = 0.032, *q* = 0.060).

**The per-site admission table exists.** `scripts/64_site_qc.py` re-runs the three admission
criteria over all 63 sites of the 21 wild-type constructs, using script 40's own functions so the
numbers cannot drift, and cross-checks against `af3_ca_sites.tsv` and `af3_batch4_geometry.tsv`.
**60 of 63 admitted.** The manuscript's old wording ("cbEGF1-3 converged poorly at pLDDT 74") was
imprecise — the failing criteria were seed spread (39.2%) at one site and occupancy (4/5) at
another, plus cbEGF5-7 site 1 on spread, which carries no variant. Outputs are
`data/processed/af3_site_qc.tsv` and `results/site_qc.md`, and it is sheet **S10** of the new
workbook.

**The supplementary workbook is now `supplementary_tables_v3.xlsx`** (11 sheets), built by the
same `43_supplementary_v2.py`. `supplementary_tables_v2.xlsx` is unchanged on disk and superseded.

**Structure changes to the manuscript.** Geometry validation moved from last to second in Results,
ahead of the claims that depend on it. Enrichment became its own short section and is framed as
descriptive and confirmatory, in the text and in the note printed under Fig. 4c. The abstract is
one unstructured 290-word block leading with the finding, roughly a third fewer numbers than
before. Results headings are one consistent declarative set. `61_typeset_v3.PLACEMENT_V3` moved
with them and now also lifts the author block out of the source for the PDF, which previously
dropped everything above the abstract.

**Administrative sections added** — author line, affiliation, corresponding author, Declarations
(funding, competing interests, ethics, contributions) and a real Data and code availability
section naming sheets rather than file paths. **Five distinct bracketed placeholders remain and must be
filled before submission**, `[Author name]`, `[Department, Institution, City, Postcode, Country]`,
`[email address]`, `[GitHub repository URL]` and `[DOI]`. Nothing was invented. `test_manuscript_v4.py`
V1 now strips URLs and DOIs before looking for colons, so pasting the real link will not fail the
gate.

**Scientific scoping.** The regulatory message is now two computational hypotheses each with a
named experiment attached, and the ClinGen Asn→Ser point reads "warrants re-examination once that
measurement exists". A new Discussion paragraph and Limitations §16 state the load-bearing
assumption (AF Server mutant models report real coordination changes), name the three variants
that anchor it, and say the 45% is a floor partly because a predictor is conservative on a single
substitution. Bond-valence sum is called a geometric proxy in Results, Methods and Limitations.

**Gates: 31 checks across three files, all passing.** `test_v2.py` is now 15/15 with **A11**,
which re-derives every admission call from the three criteria and checks the QC table against the
mutant geometry. `test_manuscript_v4.py` 8/8, `test_manuscript.py` 8/8 on the keyed paper.

**Still not fixed and still inherited.** The PDF wastes pages 2 and 4, because `56_typeset.py`
forces a page break to float a full-width figure. It is 11 pages now rather than 10. The DOCX is
unaffected and is the submission format.

## 0a. What changed on 2026-09-05 — READ THIS BEFORE TOUCHING THE MANUSCRIPT OR FIGURES

**Current manuscript is `results/manuscript_v4.md`.** `manuscript_v3.md` was an intermediate
rewrite the mentor also rejected and is superseded. `manuscript_paper.md` (63 keyed citations)
remains the citation source of record and its own gate still passes.

**Style is now measured against the five reference papers, not asserted.**
`scripts/63_style_extract.py` parses `resources/papers/` into 1,164 published prose sentences and
reports our numbers beside theirs. `results/style_extraction.md` is the report.
`scripts/62_prose_diagnostic.py` does the sentence-load and jargon count from the
`academic-paper-composer` skill. Two skills from `lishix520/academic-paper-skills` are installed
in `.claude/skills/`.

**What that measurement found, after two rewrites had failed for the wrong reasons.** Sentence
length, openers, passive rate and first-person rate were never the problem — they were already
inside the corpus range. What read as machine-written was six constructions used at many times
the published rate. "rather than" at 17x, sentence-final ", which is …" at 7x, ", which" tails at
2x, and "namely", "meta-commentary about the paper" and stacked hedges which appear **zero times
in 1,164 published sentences**. Do not reintroduce them. The corpus never lets a paper talk about
itself.

**How the corpus builds a long sentence, which is the thing to copy.** Its long sentences carry
a median of 2 clauses, 1 parenthetical and 2 commas. Ours carried 3 clauses, 0 parentheses and 3
commas. Length comes from parentheticals and noun phrases, never from chaining independent
clauses with ", and". Six passes were needed to land this, and passes 3 and 5 overshot in
opposite directions, so change it by measuring, not by eye.

**Terminology is checked against the corpus.** Twelve terms appeared nowhere in the five papers.
Five were coined here and are gone (`cognate site`, `bystander`, `composition-matched`,
`over-call`, `matches motif`). One was the wrong variant of a real term — the corpus writes
`main-chain carbonyl` 5 times and `backbone carbonyl` never, and `rmsd` not `r.m.s.d.`, both now
fixed in the text and in the figures. Six are allowed and listed in the gate, being crystallography
and statistics terms these structure papers had no occasion to use.

**Where it landed.** 3+ clause sentences 14% against 12%, median words 21 against 23, short
sentences 21% against 21%, passive 34% against 38%, all six rhetorical tells at or below corpus
rate. One target was not met and was left deliberately — long sentences are 10% against the
corpus 21%. Going back up re-created the packed sentences the mentor rejected, so the draft errs
short. That is a judgement call, not an oversight.

**Gate: `tests/test_manuscript_v4.py`, 8 checks.** V7 fails if clause load or median length drifts
above corpus, V8 fails if a term absent from the corpus is not on the allowed list, V3 re-derives
95 claim values from the keyed source so an editing pass cannot move a number.

---


**The manuscript source moved again, and the reason is editorial rather than analytical.** The
mentor rejected the 2026-09-02 draft's prose and the figure labels. Five specific complaints, all
now addressed, none of which changed a number.

| | superseded | current |
|---|---|---|
| Manuscript source | `manuscript_paper.md` (63 `[@key]` citations) | **`manuscript_v3.md`** (no citations) |
| Word file | `FBN1_manuscript.docx` | **`FBN1_manuscript_v3.docx`** |
| PDF | `FBN1_manuscript.pdf` | **`FBN1_manuscript_v3.pdf`** |
| Build | `59_typeset_docx.py` + `60_typeset_pdf.py` | **`61_typeset_v3.py`** (drives 59 and 56) |
| Gate | `tests/test_manuscript.py` | **plus `tests/test_manuscript_v3.py`** |

**Citations were removed on purpose.** The mentor is inserting them from Mendeley.
`manuscript_paper.md` and `data/processed/references_final.tsv` are intact and
`tests/test_manuscript.py` still passes on them, so nothing was lost. Do not re-add `[@key]`
markers to the v3 source. If citations are ever wanted back in the prose, they go into the keyed
source and the v3 wording is ported across, not the other way round.

**The register was taken from five published papers**, read in full and recorded in
**`results/style_reference.md`** with the rule each one contributed — Jensen 2009 (*Structure*),
Kettle 1999 (*JMB*), Downing 1996 (*Cell*), Handford 1995 (*JBC*), Godwin 2023 (*NSMB*). Read that
file before editing the prose. The four rules that matter most are moderate even sentence length,
conclusions hedged and attributed to the data, contrast carried by conjunctions rather than by
juxtaposed fragments, and no colons or semicolons anywhere. The last is a house rule, not
something the reference papers do.

**The figure vocabulary changed and the old labels must not come back.** "Matches motif, no Ca²⁺
contact", "Coordinates Ca²⁺", "Cysteine removed", "Keeps an O donor", "Site disrupted" and "The
fold, same models" were replaced by the terms the cbEGF literature already uses. The mapping is
the table at the end of `style_reference.md`, the labels live in `paperstyle.GLABEL`, and
`tests/test_manuscript_v3.py` V5 fails if a rejected label is drawn again.

**Two real bugs were found while rebuilding and are fixed.**
`56_typeset.py` never placed the supplementary figures in *any* PDF, so every PDF the project has
produced was missing Supp. Figs. S1 and S2 while the text cited them throughout. It also printed
an empty "Figure legends" heading over nothing. Both PDFs are now 10 pages with 6 figures.
`59_typeset_docx.py` keys figure placement on Results headings and dropped all four main figures
without complaining when the v3 draft renamed them, so it now stops if a figure is never placed.

**Still inherited and still unfixed.** The PDF wastes most of page 2 and page 5, because
`56_typeset.py` reaches a full-width figure by forcing a page break, which abandons the rest of
the current page. It predates this session and affects the keyed PDF identically. Fixing it needs
a real paginator in 56, not a placement tweak, because the waste depends on how much text falls
between figure anchors. The DOCX is unaffected, since Word reflows continuously.

**Gates: 28 checks across three files, all passing.** `test_v2.py` 14/14, `test_manuscript.py`
8/8 on the keyed paper, `test_manuscript_v3.py` 6/6 on the rewrite. V3 of that last gate
re-derives 95 claim-bearing values from the keyed source and fails if the rewrite moved one, which
is the check that the editing pass changed only prose.

---

**Updated 2026-09-02.** Read this first in a new chat session, then `CLAUDE.md`, then
`results/manuscript_paper.md`. This file records what state the project is in and what a fresh
session needs to know that is not obvious from the files themselves. It replaces the 2026-08-16
version; that session's content is preserved below under §6.

---

## 0. What changed on 2026-09-02 — READ THIS BEFORE TOUCHING THE MANUSCRIPT

**The manuscript source moved.** `results/manuscript_final.md` is now superseded by
**`results/manuscript_paper.md`**, and the difference is not cosmetic:

| | superseded | current |
|---|---|---|
| Manuscript source | `manuscript_final.md` (43 refs, hand-numbered superscripts) | **`manuscript_paper.md`** (63 refs, `[@key]` citations) |
| Reference list | inline in the markdown | **`data/processed/references_final.tsv`** |
| Word file | `FBN1_manuscript_final.docx` | **`FBN1_manuscript.docx`** (true two-column) |
| PDF | `FBN1_manuscript_final.pdf` | **`FBN1_manuscript.pdf`** |
| Figures | near-monochrome, blue + red accents | **pixel-exact greyscale** |

**Citations are keys, never numbers.** Reference numbers are assigned at build time by
`scripts/59_typeset_docx.py` and `60_typeset_pdf.py`, in order of first appearance. Do not
hand-number anything; add a citation as `[@key]` and rebuild. `scripts/58_manuscript_source.py`
was the one-shot migration and now refuses to overwrite the maintained file without `--force`.

**Three real corruptions were found and fixed during that migration**, all silent, all in the
superscript→key conversion, and all documented in the comment above `SUP_RUN` in script 58:
`10⁻⁴` became citation 4 in six *P* values; `10⁻⁴⁷` lost its `⁷` to citation 7; every `Ca²⁺`
became a citation of Dietz 1991. A fourth bug lived in the renderers: running a separate early
pass for the `(refs. …)` citation form handed numbers 1–2 to the wrong references and shifted
every number in the paper. All four are now guarded, and `tests/test_manuscript.py` M3/M7 fail
if numbering or marker resolution regresses.

**References are verified against two services, not asserted.** `scripts/57_references.py`
resolves every DOI through Crossref *and* Europe PMC and compares title, year, volume, first page
and first author with the citation as written; the report is
`results/reference_verification.md`. It caught a DOI that resolved to a muscarinic-receptor
paper instead of Mutalyzer, a truncated Godwin 2023 title and a paraphrased Delhomme 2022 title.
The year flags it raises for NAR and Nature entries are the online-year/issue-year split and are
**not** errors — the issue year is what a citation carries.

**The figures are monochrome by construction, not by convention.** `paperstyle.py` no longer has
accents; where two marks must be told apart they use position, then shape (circles vs diamonds in
Fig. 1b), then fill, then lightness — in that order. `tests/test_manuscript.py` M8 fails if any
figure's largest RGB spread is not exactly 0.

**Second gate.** `tests/test_manuscript.py` — 8 checks, all passing, on the document rather than
the analysis. `tests/test_v2.py` still passes 14/14 with the monochrome figures.

**Known cosmetic issue, inherited:** page 2 of the PDF is sparse because `56_typeset.py` floats
each figure to the top of a page and breaks to get it there. `FBN1_manuscript_final.pdf` has the
same gap, so this predates the rebuild. The DOCX is unaffected — Word reflows continuously.

---

## 1. The project is complete through Phase 6

Every phase gate passes, the analysis is closed, the figures are final and the manuscript is
written. What is outstanding is listed in §5 and is optional: nothing there changes a conclusion.

**The single most important thing to know is still that superseded files remain on disk.** There
are now three layers, and only the newest is current:

| | v1 (superseded) | v2 (superseded in part) | final (current) |
|---|---|---|---|
| Statistics | `phase5_statistics.md` | `phase5_statistics_v2.md` | **`statistics_final.tsv` / `.md`** |
| Manuscript | `FBN1_manuscript.docx` | `manuscript_v2_draft.md` | **`manuscript_final.md`** + `FBN1_manuscript_final.docx` |
| Supplement | — | `supplementary_tables.xlsx` | **`supplementary_tables_v2.xlsx`** |
| Figures | `fig1`–`fig6` | `v2_fig1`–`v2_fig8` | **`fig1_module`, `fig2_stability`, `fig3_coordination`, `fig4_interpretation`, `figS1_validation`, `figS2_structure`** |
| Figure style module | `figstyle.py` | `figstyle.py` | **`paperstyle.py`** |
| Ca sites defined by | sequence motif | observed AF3 coordination | observed AF3 coordination |
| Variants with structure | 751 | 3,831 | 3,831 |

`phase5_statistics_v2.md` is superseded because it ran on 2026-08-09, before the PyRosetta ΔΔG was
corrected and before batch 4 existed. Do not quote it.

**The manuscript is now typeset as a two-column PDF**, `results/FBN1_manuscript_final.pdf`, by
`scripts/56_typeset.py` (reportlab). Its geometry is measured off Godwin 2023 with PyMuPDF, not
invented: 595.28 × 790.87 pt page, text block x 39.7–563.4 / y 46.6–745.7, two 256.9 pt columns
with a 9.9 pt gutter, 8.2 pt serif body, 7 pt captions, figures floated full width to the top of
a page. Two traps: reportlab's base-14 Times is Latin-1 only, so Greek and arrows are set in an
embedded DejaVu Serif face (the classic Symbol-font trick renders as boxes because base-14 fonts
are not embedded); and the front-matter frame height must include each style's spaceBefore and
spaceAfter or the last abstract paragraph spills into column one. The markdown remains the source
of truth — edit it and re-run 56 and 45.

**The whole figure set was rebuilt twice on 2026-08-30.** The second rebuild answered a specific
criticism: the figures still looked like plotting-library output. Three fixes, all now in
`paperstyle.py` with the reasoning: **near-monochrome** (black / grey / open, with exactly two
accents — blue for calcium wherever it appears, red for the one thing the reader is being asked
to look at), because the data panels on the reference pages are almost entirely black and white;
**categories on the y axis** wherever a panel is narrow, because four descriptive names do not fit
side by side and the alternatives are rotation or jargon; and **labels that say what a group is**
— "Matches motif, no Ca²⁺ contact", not "Motif, not ligand", which meant nothing without the
Methods open beside it. Groups are separated by fill, not hue, so the panels survive greyscale.

The first rebuild on 2026-08-30 after the mentor rejected the previous one as
confusing and un-publishable. What was wrong, so it is not repeated: a whole panel was spent on
Cliff's δ with bootstrap intervals instead of putting brackets and *P* values on the data;
comparisons were exiled to a separate panel; category labels were rotated 30° (a matplotlib
default that appears in none of the reference papers); half of Figure 1 was two bar charts of
numbers that belong in a sentence; and `figstyle.apply_style()` **banned panel titles**, which was
simply wrong — Cheng 2023 (Science) and Godwin 2023 (NSMB) title nearly every panel. The rule that
replaced it: a panel title says what is plotted, never what it means. `scripts/paperstyle.py`
carries the conventions and the reasoning; `tests/test_v2.py` F3 enforces them.

---

## 2. What the 2026-08-29 session did

**Batch 4 came back.** The mentor ran all 60 AlphaFold Server mutant jobs; they are in
`structures/alphafold/batch4/`. Every one validated: right length, exactly one substitution, at
the position and with the residues the manifest claimed, five models each.

*(Two housekeeping notes. `fbn1_cbegf5_7_d723a_ca3` was unzipped both into `batch4/` and at the
top level of `structures/alphafold/` — byte-identical, harmless. `fbn1_cbegf25_27_r1832s_ca3_2`
is a genuine second run of one job with a different server seed; it is not in the manifest and is
not used, but it is a free replicate if anyone wants one.)*

**Figure 4 became a rate.** `scripts/40_af3_batch4_geometry.py` measures every mutant against its
wild-type construct; `41_batch4_rate.py` turns that into the rate and its statistics.

- **45% (15/33, 95% CI 30–62%)** of pathogenic side-chain Ca²⁺-ligand variants disrupt their
  calcium site, against **1 of 24** controls. OR 19.2, p = 7 × 10⁻⁴.
- **Donor chemistry decides it**: substitutions that remove the coordinating oxygen disrupt in
  12/16 (75%); those that leave another oxygen donor (D→N, D→E, N→D, E→Q) in 3/17 (18%). OR 14.0,
  p = 0.002. **This is the session's genuinely new mechanistic finding.**
- **Neither instrument in routine use sees it.** Rosetta ΔΔG is uncorrelated with the coordination
  change across those variants (ρ = −0.089), and AlphaMissense scores 0.982 for both the variants
  that break the site and those that do not.
- The local fold does not move: 0.29 Å backbone RMSD around the site, against 0.24 Å for controls.

**Three methodological decisions worth not re-litigating:**

1. **The cognate site, not all three.** A construct has three calcium sites and a substitution can
   reach one. Measuring all three dilutes a real effect by two bystanders. The other two are kept
   as within-construct controls, which is where most of the null comes from.
2. **The null is measured, not assumed.** 109 bystander sites + 14 negative-control sites = 123
   measurements where no effect is possible; the threshold is the 95th percentile of their
   absolute change, **9.14%**. `tests/test_v2.py` A8 re-derives it and fails if it drifts.
3. **Local superposition, not global.** Whole-construct RMSD reaches 29.8 Å because these tandem
   models hinge between domains between seeds — that is inter-domain motion, not a broken fold,
   and a global superposition turns it into a fake ion displacement. Everything is superposed on
   the residues within 12 Å of the ion.

**A quality filter fixed on wild-type evidence alone** (seed spread ≤ 15% of valence, site
occupied in 5/5 wild-type models, local pLDDT ≥ 70) excluded the cbEGF1–3 construct's sites and
with them three variants: D490Y, N261S, E487A. Without it, one of those returned a +237% valence
change that was the ion having moved 25 Å between folds. The filter is computed before any mutant
is read so it cannot be tuned on the answer.

**Two corrections to the project's own records** (2026-08-30). `resources/reference/lit_notes.md`
recorded the ClinGen VCEP and Godwin PDFs as "WRONG PAPER"; both are present and correct — the
mentor must have replaced them. The ClinGen article number is **154, not 141** as CLAUDE.md
Appendix A says, and HMG 9:1987 is **McGettrick** *et al.*, not Whiteman.

**Reading the VCEP changed two claims in the paper.** It does **not** weight the two features
equally: cysteine-removing in cbEGF is PM1_**Strong**, the calcium consensus is PM1 moderate. It
keys the calcium rule to the sequence consensus — which our models show over-calls by 143
positions that show no enrichment. And it **exempts Asn→Ser at the second [D/N]** from PM1, which
is exactly where our donor-chemistry rule breaks: all three Asn→Ser variants we folded there lost
15–18% of coordination, and N2144S at that position is the one FBN1 variant with a measured
ninefold affinity loss. That is now a Results subsection and the sharpest practical claim in the
paper. It also uses **REVEL ≥ 0.75** for PP3, not AlphaMissense — so our saturation observation is
about a predictor the guideline does not use, and testing REVEL is an easy open question.

**Also produced:** the consolidated statistics table (`44_statistics_final.py`), the v2
supplementary workbook (`43_supplementary_v2.py`, 10 sheets, every column documented — it refuses
to build if one is not), the final manuscript in markdown and Word (`45_...docx.py`), and a
supplementary structural figure (`46_v2fig8_structure.py`) that rebuilds the old PyMOL render with
v2 definitions. That rebuild mattered: v1 coloured only conserved cbEGF cysteines and so drew 12 pathogenic
cysteine variants at 8 positions on 2W86; under the v2 definition it is 23 variants at 15
positions.

**A validation worth keeping.** Jensen 2009 names cbEGF9's ligands: Asp807, Glu810, Asn823 side
chains plus Ile808, Ser824, Ser827 carbonyls. The AF3 models return exactly those six with the
same split, having never seen that assignment. It is in the manuscript and drawn in Supp. Fig. S2b.

---

## 3. Decisions the mentor has made (do not re-ask)

- **v2 is canonical.** (2026-08-15)
- **No second CheckMyMetal hand-off.** Local bond-valence matches the server at r = 0.961 with a
  constant offset, which cancels in wild-type-minus-mutant differences.
- **PyRosetta is worth doing** — it answers the sharpest reviewer objection to the headline.
- **Enrichment is demoted**, to Figure 6, with the circularity caveat printed inside the figure.
- **Git does not matter for now.** The repo still has **zero commits** — there is no safety net,
  so prefer writing new files over overwriting good ones. That is why v1 and v2 artifacts are all
  still on disk. *A first commit would be cheap insurance and has never been made.*
- Standing rule from memory: **ask the mentor before any step needing an external structure
  service** (AlphaFold Server, CheckMyMetal, AlphaFill) — before preparing the request, not just
  before running it.

---

## 4. The gate

`tests/test_v2.py` — **14 checks, all passing**, 10 on the analysis tables and 4 on the eight
figures. Run it with `.venv/bin/python tests/test_v2.py`; F4 re-runs every figure script, which
now includes a PyMOL render, so it takes about a minute.

The four added on 2026-08-29 are worth knowing because each guards a specific way this analysis
could rot: **A7** every requested mutant was measured and carries the substitution requested;
**A8** the null threshold re-derives from the bystander sites; **A9** every disruption call
follows from that threshold rather than being set by hand; **A10** the statistics table is
internally coherent and never compares FoldX to Rosetta.

The v1 gates are left in place and still pass — they validate v1 artifacts, which still exist.
A green v1 gate says nothing about the final analysis.

---

## 5. Open work — all optional

1. ~~References~~ — **done**. 43, every one verified against Europe PMC, the publisher record or
   the PDF. Twelve were added from the literature on 2026-08-30, several of which changed the
   paper: Handford 1991 (*Nature* 351:164, the factor IX precedent), Reinhardt 2000 (*JBC*
   275:12339 — **a second experimental anchor: N548I and E1073K are both in our cofolding set**),
   Haller 2020 (*Sci Rep* 10:16844, the closest prior computational work) and the AlphaFold 3
   paper that was previously missing. All are tabulated in `lit_notes.md`.
2. **Boltz-2** as an independent check that AF3's calcium placement is not a single-model
   artifact. `.venv-boltz` still fails with `ImportError: cannot import name
   'is_fx_symbolic_tracing'`; the error text says install `triton==3.3.0`. Setup is
   `scripts/31_boltz_setup.py`, input YAML already written. Mentor called it worth one timeboxed
   attempt; drop it if the environment fight drags.
3. **The 18 unexplained variants** (§15 of `LIMITATIONS.md`): pathogenic Ca²⁺-ligand variants that
   keep their calcium. β-hydroxylation, an effect too small to model, or ClinVar
   overclassification — this analysis cannot separate them. The testable prediction is that D→N
   and D→E at these positions have milder binding defects than D→G or D→Y at the same position.
4. **The "indistinguishable" wording** in `results/phase5_report.md` §1 — still the v1 wording,
   still the mentor's call. The final manuscript does not use it.
5. **A first git commit.**

---

## 6. Environment facts that cost time to rediscover

- Three virtualenvs. `.venv` (Python 3.14) is the project default and has pandas, matplotlib,
  gemmi. **`.venv-rosetta` (3.11) is the only one with PyRosetta**, `.venv-boltz` (3.11) is for
  Boltz. `python` is not on PATH — invoke `.venv/bin/python` explicitly.
- Join v2 data through `figstyle.load_v2()` / `mech_group_v2()`. That module is the only place the
  four mechanism-group definitions are written down.
- Rosetta scratch lives in `/tmp/fbn1-rosetta` and **`/tmp` gets cleared** — pre-relaxed constructs
  are silently rebuilt (~80 s each) after a reboot. The relaxed-pose filename encodes protocol and
  metal setup, so a Cartesian-relaxed pose can never be reused for a torsion run.
- 24 cores; `--jobs 8` is comfortable. GPU is an RTX A6000.
- **The four PyRosetta traps found on 2026-08-16**, all fixed and commented in
  `scripts/29_rosetta_ddg.py`, all silent if reintroduced: `pyrosetta.toolbox.mutate_residue` is a
  no-op that returns a copy; Cartesian minimisation cannot be combined with metals
  (`cart_bonded` scores calcium bonds against parameters that do not exist — 7,693 of 12,336 REU
  from three ions); harmonic metal constraints detonate during repacking (0 → +17,364 REU) and the
  constraint helper re-enables its own score-function weights when called, so it must simply not
  be called; and `--repeats` bought nothing because a 12-residue repack shell converges to one
  optimum. `MutateResidue` does clean up disulfides correctly — that one was checked and is fine.
