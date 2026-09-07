# Phase 3 report — Domain & site annotation

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED (9/9) · **Awaiting approval for Phase 4.**

The calcium-binding consensus was **derived from the data, not copied from a paper**, then
validated three independent ways before use. One important caveat about circularity is in §4 —
please read it before the enrichment numbers are treated as a result.

---

## 1. How the calcium sites were derived

Each cbEGF domain is anchored on its six cysteines; consensus positions are defined as offsets
from those anchors. The offsets came from **measuring actual Ca²⁺ coordination in the deposited
structures**, not from assumption:

| Derived position | Residue observed | Role |
|---|---|---|
| **C1−4** | Asp (43/43) | sidechain Ca²⁺ ligand |
| C1−3 | variable | backbone carbonyl ligand |
| **C1−2** | Asp or Asn only (25/18) | consensus, not a direct ligand |
| **C1−1** | Glu (43/43) | sidechain Ca²⁺ ligand |
| **C3+2** | Asn (42) / Asp (1) | sidechain ligand — the β-hydroxylation site |
| C3+3, C3+6 | variable | backbone carbonyl ligands |
| **C4−2** | Tyr or Phe only (23/20) | aromatic core |

### Three independent validations

1. **Structure (ground truth).** Ca²⁺ ligands were measured directly from 2W86, 1LMJ, 1UZJ and
   1EMN — 8 cbEGF domains, 16 calcium sites. Every site shows the identical 5–6 ligand pattern,
   and **all 24 sidechain ligands are recovered by the derived rule**. The gate re-checks this
   and fails on a one-residue shift.
2. **Sequence.** Across all 43 cbEGF domains the derived positions carry exactly the residues the
   ClinGen VCEP consensus `[D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F]` specifies — **zero exceptions**.
   Note this settles a discrepancy between sources: the VCEP's `[D]` (Asp only) and `[E/H]` are
   correct for FBN1, where Baudhuin 2019's looser `D/N` and `E/Q` would have been over-permissive.
3. **Literature.** All six of Baudhuin 2019's published tolerated Asn>Ser positions (615, 1030,
   1282, 1489, 2526, 2650) land **exactly** on the derived **C1−2**, and their three D↔N
   positions (1197, 1240, 1907) land on derived consensus positions. **9/9.** This also resolves
   which position the VCEP's "second D/N" exception refers to.

### Independent count checks

- **258 cbEGF cysteines** — matches the VCEP's stated figure **exactly**.
- **All 129 cbEGF disulfides follow C1–C3, C2–C4, C5–C6** (43 of each), with no exceptions,
  confirming the canonical pattern against UniProt's own annotations.
- 215 calcium-consensus positions (43 × 5); 87 interdomain-packing glycines.

**One discrepancy, reported rather than tuned away.** My full VCEP PM1 position set totals
**655**; the VCEP paper states **633** (22.0% of 2,871). The 3.4% gap is in the auxiliary
categories (which glycines and which non-cbEGF cysteines qualify), whose exact membership lives
in the paper's Supplementary Table 1, which we don't have. **I stopped adjusting parameters
rather than reverse-engineer a match** — tuning until a number agrees would manufacture false
precision. The study's two primary feature classes are unaffected and exactly validated:
cbEGF cysteines (258, exact) and calcium-consensus residues (215, structurally confirmed).

---

## 2. Replacement for an impossible gate check

CLAUDE.md's Phase 3 gate asks that the coordinating-residue set be validated against "UniProt's
explicit Ca-binding annotations". **UniProt P35555 has zero `Binding site` features**, so that
check cannot be performed. It is replaced by a stronger one: agreement with calcium ligands
**measured from experimental coordinates** (§1.1). The gate documents this substitution in its
own docstring.

---

## 3. Cross-tabulation (primary tier, ≥2★)

| feature | pathogenic (n=578) | benign (n=34) | VUS (n=823) |
|---|---|---|---|
| in cbEGF domain | 425 (73.5%) | 21 (61.8%) | 486 (59.1%) |
| Ca consensus residue | 61 (10.6%) | 1 (2.9%) | 48 (5.8%) |
| — of which sidechain ligand | 42 (7.3%) | **0 (0.0%)** | 26 (3.2%) |
| cbEGF cysteine (removing) | **304 (52.6%)** | **0 (0.0%)** | 0 (0.0%) |
| cysteine-creating | 52 (9.0%) | 0 (0.0%) | 24 (2.9%) |
| interdomain packing glycine | 22 (3.8%) | 0 (0.0%) | 21 (2.6%) |
| VCEP critical residue | 526 (91.0%) | 0 (0.0%) | 85 (10.3%) |
| neonatal region (TB3–cbEGF18) | 134 (23.2%) | 14 (41.2%) | 185 (22.5%) |

Full tables including the ≥1★ tier: `results/phase3_crosstab.md`.

---

## 4. ⚠️ The result is partly circular — this must shape how it is reported

**91% of pathogenic variants sit at a VCEP critical residue, versus 0% of benign.** That looks
spectacular, and CLAUDE.md says to sanity-check anything that looks too strong before
celebrating. It is partly an artifact:

> **ClinVar's pathogenic classifications were themselves made using PM1 — the very
> critical-residue definition being tested.** Testing enrichment of ClinVar-pathogenic variants
> at ClinGen-defined critical residues is, in part, measuring the classification rule rather
> than biology.

I ran three probes:

| Probe | Result |
|---|---|
| **AlphaMissense as an orthogonal label** (never used FBN1 PM1 rules) | AM-pathogenic: 33.1% at cbEGF cysteines, 14.8% at Ca consensus. AM-benign: **0.0%** and 1.9%. **Direction and separation reproduce independently.** |
| ClinVar P/LP evaluated **before 2019** | 16/16 (100%) at critical residues — but n=16, too small to lean on |
| ClinVar P/LP evaluated **2024+** | 325/362 (89.8%) — i.e. the modern, VCEP-era calls |

**Interpretation.** The *direction* of enrichment is real and independently corroborated: an
orthogonal predictor that never saw the VCEP rules separates the classes just as cleanly at
cysteines (33.1% vs 0.0%). But the *magnitude* from ClinVar labels (91% vs 0%) is inflated by
circularity and must not be presented as an independent discovery.

**How I propose to handle it in Phase 5 and the manuscript:**

1. **Do not headline the ClinVar-vs-ClinVar enrichment.** Report it as a consistency check on the
   annotation, with the circularity stated plainly in the same sentence.
2. **Headline the physics.** Phase 5's biophysical measures — Ca²⁺ coordination geometry from
   CheckMyMetal, FoldX ΔΔG, solvent accessibility, disulfide disruption — are computed from
   structures and are **completely independent** of any classification scheme. That is the
   study's genuine contribution: not "pathogenic variants are at critical residues" (known), but
   *how much* they perturb coordination geometry and stability, and whether that predicts
   severity.
3. **Use the all-possible-missense background** (AlphaMissense's 54,549 substitutions) as the
   denominator, so the question becomes "do pathogenic variants concentrate at these sites more
   than chance?" rather than "do two label sets differ?".
4. **State it in `results/LIMITATIONS.md`** as a primary caveat.

This is a change in framing, not scope. Flagging it now because it determines what Phase 5
computes and what the paper can claim.

---

## 5. Outputs

| Path | Contents |
|---|---|
| `data/processed/variant_annotated.tsv` | 4,451 variants × 56 columns (Phase 2 table + 22 annotation columns) |
| `data/interim/cbegf_sites.tsv` | per-domain derived site positions, 43 rows |
| `data/interim/disulfide_pairing.tsv` | 155 disulfides classified C1-C3 / C2-C4 / C5-C6 |
| `results/phase3_crosstab.md` | cross-tabs for both star tiers |
| `logs/06_annotate_*.log` | run transcript |

**Gate:** `tests/test_annotation.py` — 9/9. Self-tested against five simulated failures
(cysteine role on a non-Cys residue, inverted enrichment direction, contradictory cys flags,
motif violation, and a one-residue shift of a structurally-observed calcium ligand); all five
were caught.

**No Phase 4 work has started** — and per your standing instruction I will not prepare any
structure or CheckMyMetal request until you approve it.
