# FBN1 / Marfan project — briefing for a new team member

*Written for someone with introductory-biology background and no prior research experience.
No jargon is used without being defined.*

**Status: 2026-08-02.** Four of seven phases complete.

---

## 1. The biology you need

**Fibrillin-1** is a very large protein (2,871 amino acids) that cells secrete into the space
between them. Copies of it assemble end-to-end into **microfibrils** — think of tiny ropes woven
through the tissue that gives your blood vessels, skin, ligaments and the suspension of your eye
lens their stretch and strength. The gene that codes for it is called **FBN1**, on chromosome 15.

When FBN1 is broken, you get **Marfan syndrome** — a dominant genetic disorder (one bad copy is
enough) affecting roughly **1 in 10,000 people**. The dangerous part is the aorta: the wall
weakens and can tear.

Fibrillin-1 is built from repeating modules called **domains**, like beads on a string. The most
common bead type — **43 of them** — is the **calcium-binding EGF-like domain**, abbreviated
**cbEGF**. Each cbEGF bead has two features that matter here:

- **One calcium ion (Ca²⁺) clamped into it.** The calcium is held by a handful of specific amino
  acids whose side chains point inward and grip it. The calcium acts like a splint: it stiffens
  the joint between neighbouring beads so the rope stays straight and rigid.
- **Six cysteines forming three disulfide bonds.** Cysteine is the one amino acid that can form a
  covalent bridge to another cysteine. In each cbEGF the six cysteines pair up in a fixed
  pattern — 1st with 3rd, 2nd with 4th, 5th with 6th — stapling the bead into its correct shape.

So each bead is held together by three internal staples and stiffened by one calcium clamp.

A **missense variant** is a single DNA change that swaps one amino acid for a different one. If
that swap deletes a cysteine, one staple can't form. If it deletes a calcium-gripping residue,
the clamp loosens. Either way the rope gets floppy — which is the leading explanation for how
Marfan syndrome happens.

---

## 2. The question this project asks

> **Do the disease-causing FBN1 variants pile up specifically at the calcium-gripping residues
> and the six cysteines — more than you would expect by chance, and more than harmless variants
> do? And when they do, how badly do they actually distort the calcium clamp and the protein's
> stability?**

The first half is a *counting* question. The second half is a *physics* question, and it is the
part that would be genuinely new.

---

## 3. What has been done so far

### Phase 0 — Set up and read the literature
Checked that all the software works (PyMOL for 3D graphics, DSSP and FoldX for structural
calculations), confirmed we can reach every database, and read **19 research papers**, writing a
summary of each.

Two of the supplied PDFs turned out to be **completely the wrong papers** — the file named for
the official variant-classification rulebook actually contained a paper about protein droplets.
Caught by searching each PDF for the word "FBN1" and finding zero hits. Both were replaced.

### Phase 1 — Collect the data
Downloaded from four sources:

| Source | What it is | What we got |
|---|---|---|
| **ClinVar** | Public archive where labs deposit variants and say whether they cause disease | **9,327** FBN1 variants |
| **UniProt** | Reference database for protein sequence and structure annotation | The 2,871-amino-acid sequence, all 56 domain boundaries, 155 disulfide bonds |
| **gnomAD** | Catalogue of variants found in ~800,000 healthy people — if a variant is common here, it probably isn't causing a rare disease | **12,107** FBN1 variants with frequencies |
| **AlphaMissense** | An AI model that scores *every* possible amino-acid swap for harmfulness | **54,568** scores (every possible change to this protein) |

### Phase 2 — Clean and label the variants
Converted every variant into "position + original amino acid + new amino acid," then ran the
single most important check in the project: **does the original amino acid the database claims is
at that position actually match the real protein sequence there?** If a database is quietly using
a different version of the gene, everything downstream is wrong by a few positions.

**Result: 4,451 out of 4,451 matched. Zero errors.** Independently confirmed by re-checking 100
random variants through a separate validation service (100/100 agreed).

Of the 9,327 ClinVar records, **4,451 were missense** (the type we can analyse structurally). The
other 4,876 — variants that truncate the protein, don't change the protein, or affect splicing —
were kept and labelled, never deleted.

Then the variants were sorted into groups using thresholds agreed in advance:

| Group | Count | Meaning |
|---|---|---|
| **Pathogenic** | **578** | confidently disease-causing |
| **Benign** | **34** | confidently harmless |
| **Uncertain (VUS)** | **823** | genuinely unknown — the clinical problem |

**Why so few benign?** Not a mistake. FBN1 is under such strong evolutionary pressure that
harmless changes to it are rare — only 24 FBN1 missense variants are common enough in healthy
people to call benign on frequency alone, and all 24 were already in our set. This is a real
constraint on the study, and it shapes how the statistics must be done.

### Phase 3 — Find the calcium sites and the cysteines
Rather than copy the calcium-binding pattern out of a textbook, it was **derived from actual 3D
structures**. Eleven experimental structures of FBN1 fragments were downloaded, and the
distances from each calcium ion to every nearby atom were measured. Across **16 calcium sites**
in 8 different domains, the exact same set of gripping residues appeared every single time.

That measured rule was then applied to all 43 cbEGF domains and checked three ways:

- It matches the official clinical rulebook's pattern in **43 out of 43** domains, no exceptions.
- It predicts **258 cysteines** in cbEGF domains — the rulebook independently states 258.
- All **129** cbEGF disulfide bonds follow the expected 1–3 / 2–4 / 5–6 pattern, no exceptions.
- **9 out of 9** specific positions published in an independent 2019 paper land exactly where
  the rule predicts.

---

## 4. What the data shows — and one important warning

Among the confidently-labelled variants:

| Feature | Pathogenic (578) | Benign (34) |
|---|---|---|
| Removes a cbEGF cysteine | **304 (52.6%)** | **0 (0%)** |
| Hits a calcium-gripping residue | 61 (10.6%) | 1 (2.9%) |
| Hits any "critical" position | **526 (91%)** | **0 (0%)** |

That looks like a spectacular confirmation. **It is partly circular, and must not be reported as
a discovery.** The reason: the experts who labelled these variants as pathogenic *used the rule
that they sit at critical residues* as part of their evidence. So we are partly measuring the
labelling rule, not nature.

To test this, the same comparison was run using the AI model **AlphaMissense**, which was built
without any knowledge of the FBN1-specific clinical rules. It reproduces the pattern
independently: **33.1%** of its predicted-harmful variants remove a cbEGF cysteine versus
**0.0%** of its predicted-harmless ones. So the underlying biology is real — but the 91%-vs-0%
figure is inflated.

**Consequence:** the headline of this study should be the *physics* — how far the calcium ion is
pushed out of position, how much stability is lost — because those are measured from structures
and owe nothing to anyone's classification decision.

### Bugs caught that would have silently ruined the results

1. A database query used a field name that doesn't exist. Instead of erroring, the NCBI server
   quietly reinterpreted it as a plain text search for "2200" and returned thousands of records
   **from completely different genes**.
2. The obvious gnomAD frequency field returns *empty* for every FBN1 variant. Using it would have
   produced a benign set with no frequency evidence at all — which would have looked like a
   clean result rather than a bug.
3. AlphaMissense lists **two different** "original" amino acids at position 472. Matching on
   position alone would have double-counted it; matching on identity rejects the 19 bad rows.

---

## 5. Deliverables produced so far

| File | What it contains |
|---|---|
| `data/processed/variant_annotated.tsv` | **The main product.** 4,451 variants × 56 columns: position, amino acids, clinical label, frequency, AI score, which domain, whether it hits a calcium site or cysteine, which disulfide it breaks |
| `data/processed/variant_master_dictionary.md` | Plain-English definition of every column |
| `data/interim/cbegf_sites.tsv` | The derived calcium-site positions for all 43 cbEGF domains |
| `data/interim/disulfide_pairing.tsv` | All 155 disulfide bonds, classified by pairing pattern |
| `resources/reference/structures.md` | The 11 experimental structures, their calcium content, and a numbering warning |
| `resources/reference/inclusion_criteria.md` | The thresholds, fixed in advance so they can't be tuned to a result |
| `results/phase0–3_report.md` | Full record of each phase |
| `manifest/` | Where every file came from — source, version, date, checksum |
| `tests/` (4 scripts) | Automatic checks that must pass before the next phase runs |
| `structures/` | 11 experimental structures + 4 AI-predicted models with calcium |

**On the test scripts:** each one was deliberately fed corrupted data — wrong amino acids,
swapped labels, deleted rows, a calcium site shifted by one position — to confirm it actually
fails. A check that can't fail is worthless.

---

## 6. What happens next

- **Phase 4** — build 3D models of the cbEGF domains with calcium in place, and have them
  validated by CheckMyMetal, a web service that judges whether a metal site is geometrically
  believable. *(Requires the supervisor to run these.)*
- **Phase 5** — the real measurements: how much each variant distorts the calcium geometry
  (CheckMyMetal), how much stability it costs (FoldX), how buried the residue is, and the
  statistics.
- **Phase 6** — publication figures and tables.

**One thing worth understanding about how this project is run:** nothing is ever guessed. If a
value can't be computed or a file is missing, it gets recorded as a gap. A missing number is a
result; an invented one is misconduct.
