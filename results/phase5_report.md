# Phase 5 report — Structural & biophysical analysis

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED (10/10) · **Phase 6 may proceed.**

Framed per your instruction: the biophysics is the headline; enrichment is supporting context
carrying its circularity caveat.

---

## 1. The central result — the two Marfan mechanisms are biophysically dissociable

The hallmark FBN1 missense mechanisms are cysteine substitution (losing a disulfide) and
disruption of the calcium-binding consensus. They have always been described together as
"structurally critical." **They are not the same kind of damage at all.**

| Variant class | n | Median ΔΔG (kcal/mol) | Median AlphaMissense | Median RSA |
|---|---|---|---|---|
| **Direct Ca²⁺ ligand** | 78 | **0.76** | **0.975** | 26.9% |
| **cbEGF cysteine-removing** | 146 | **4.19** | 0.998 | 12.4% |
| Other cbEGF residue | 300 | 0.56 | 0.314 | 45.4% |

Read the first two rows against the third:

- **Cysteine removal destabilises the fold.** ΔΔG 4.19 vs 0.56 for control residues
  (Cliff's δ = +0.636, q < 0.001), driven by the FoldX disulfide term specifically
  (2.87 vs 0.00, δ = +0.791, p = 3×10⁻¹¹²). The staple is lost and the domain pays for it.
- **Calcium-ligand loss does not destabilise the fold at all.** ΔΔG 0.76 — statistically
  indistinguishable from ordinary cbEGF residues at 0.56 — **yet AlphaMissense scores them 0.975
  versus 0.314** (δ = +0.551, p = 6×10⁻¹⁴). Same negligible stability cost, completely different
  predicted consequence.

So calcium-site variants are pathogenic **functionally, not thermodynamically**: they abolish the
metal clamp while leaving the fold intact.

### It is not a burial artifact

The obvious objection is that calcium ligands are surface residues, and surface mutations are
always cheap. Stratifying by solvent accessibility kills that explanation:

| Burial stratum | Ca-consensus ΔΔG | Other cbEGF ΔΔG | p |
|---|---|---|---|
| Buried (RSA < 20%) | 0.87 (n=29) | **3.62** (n=174) | 1.3×10⁻⁶ |
| Partial (RSA 20–50%) | 0.78 (n=47) | 0.90 (n=132) | 0.58 |

**Even among buried residues**, calcium-consensus positions cost ~4× less folding energy than
their neighbours (δ = −0.562). The effect is a property of the site, not of its exposure.

> **Corrected 2026-08-08 during Phase 6.** The buried row originally read 3.54 (n=170,
> δ = −0.557). The stratification used `(rsa_pct or -1) >= lo`, and because `0.0` is falsy in
> Python a residue with RSA exactly 0.0 was rewritten to −1 and fell out of *every* stratum.
> That silently dropped 4 completely buried "other cbEGF" variants — the most destabilising
> ones, median ΔΔG 13.90 — from the comparison. `12_statistics.py` now tests for `None`
> explicitly. Only this one row of the 20 changed, the effect moved away from zero rather than
> toward it, and q is unchanged at 2.5×10⁻⁶; the conclusion is unaffected. No Ca-consensus
> variant has RSA 0.0, so the Ca arm was never affected.

### Why this matters beyond FBN1

**A ΔΔG-based variant predictor would systematically miss calcium-site pathogenic variants in
this gene.** They look thermodynamically benign because they are. Any pipeline that ranks FBN1
missense variants by predicted destabilisation alone will under-call one of the two principal
disease mechanisms. That is a concrete, testable methodological warning and it falls directly out
of the measurements.

---

## 2. Supporting biophysics

| Comparison | n | Median | Effect (Cliff's δ) | q (BH) |
|---|---|---|---|---|
| ΔΔG pathogenic vs benign | 98 / 5 | 3.94 vs 0.69 | +0.580 [+0.278, +0.833] | 0.027 |
| ΔΔG buried vs exposed | 323 / 428 | 3.20 vs 0.63 | +0.480 [+0.400, +0.555] | <0.001 |
| RSA cbEGF Cys vs other cbEGF | 146 / 378 | 12.5% vs 36.8% | −0.709 [−0.768, −0.647] | <0.001 |
| Spearman ΔΔG vs AlphaMissense | 751 | ρ = +0.609 [+0.561, +0.651] | — | <0.001 |

The cysteines are buried (12.5% RSA) — consistent with their role stapling the hydrophobic core —
while calcium ligands sit at the surface where a metal can reach them. Two entirely independent
predictors, a physics force field and a deep-learning model, correlate at ρ = 0.61: enough to be
mutually corroborating, not so much as to be redundant.

## 3. CheckMyMetal — AF3 calcium is indistinguishable from experimental calcium

| Parameter | Experimental (n=8) | AF3 (n=8) |
|---|---|---|
| gRMSD (°) | 25.5 | 24.6 |
| nVECSUM | 0.22 | 0.26 |
| Valence | 1.45 | 1.65 |
| Coordination number | 6.5 | 6.0 |

Mann-Whitney on gRMSD: **p = 0.248, not significant** — and here non-significance is the desired
outcome. Combined with the 0.32 Å median displacement from Phase 4, AF3 cofolding is validated
for this domain family.

**⚠️ CMM's absolute thresholds are meaningless for cbEGF sites.** The experimental controls score
as "badly" as the models: deposited structures at 1.80 Å and 2.25 Å return gRMSD ≈ 25–27° and
valence ≈ 1.3–1.8, which CMM's published bands call outliers. cbEGF calcium sites are irregular,
carbonyl-dominated and low-symmetry; the gRMSD measures template mismatch, not site quality. This
is exactly why the experimental controls were submitted alongside the models — without them we
would have reported the AF3 sites as failing validation. **No "N% of sites are CMM outliers"
statement is made anywhere in this project.**

## 4. Enrichment (supporting only)

Against the all-possible-missense background (position share is the correct null, since every
position offers the same 19 substitutions):

| Feature | Observed | Expected | Fold-enrichment | p |
|---|---|---|---|---|
| cbEGF cysteine (removing) | 52.6% | 9.0% | **5.9× [5.4, 6.3]** | 8×10⁻¹⁵⁸ |
| Calcium-consensus residue | 10.6% | 7.5% | **1.4× [1.1, 1.8]** | 0.005 |

Note the asymmetry — cysteines are enriched ~6-fold, calcium residues only 1.4-fold. That is
consistent with §1: cysteine variants are both more damaging *and* more readily recognised as
pathogenic, whereas calcium-site variants are under-ascertained precisely because they look benign
to stability-based reasoning.

**Circularity caveat.** The pathogenic-versus-benign contrast (cbEGF cysteine OR = 76.5, any
critical residue OR = 692) is inflated: ClinVar's pathogenic calls used ClinGen PM1, the same rule
being tested. The orthogonal AlphaMissense contrast — which never used FBN1-specific rules — gives
OR = 8.9 [6.1, 12.9] at calcium-consensus residues and confirms the direction independently.

## 5. Gate and reproducibility

`tests/test_analysis.py` — **10/10**, including:

- All 16 CheckMyMetal sites accounted for, none silently missing; parameters within physically
  plausible ranges.
- All three processed tables cover exactly the 4,451 Phase 3 variants.
- 751/751 structurally-covered variants have a ΔΔG; **0 failures**, and all residue identities
  re-verified against P35555 before FoldX was invoked.
- A directional sanity check: cysteine-removing ΔΔG must exceed other variants, or the mutation
  strings were built wrong.
- **Re-running the statistics reproduces byte-identical output** (seed 20260802).

### One statistical bug caught and fixed

The fold-enrichment confidence intervals initially read `5.9× [5.5, 11.1]` and `1.4× [1.1, 13.4]`.
Those upper bounds were artifacts: I had taken the interval from a **one-sided** binomial test,
which pins its upper bound at 1.0 and therefore reports a meaningless ceiling of 1/expected. The
p-value is legitimately one-sided (we test for enrichment), but the interval must be two-sided.
Corrected to `[5.4, 6.3]` and `[1.1, 1.8]`.

---

## 6. Outputs

| Path | Contents |
|---|---|
| `data/processed/cmm_sites.tsv` | 16 calcium sites × all CMM parameters + threshold bands |
| `data/processed/structural_metrics.tsv` | 4,451 rows; RSA, DSSP, Ca distance, ligand status, disulfide geometry for the 751 covered |
| `data/processed/foldx_ddg.tsv` | 4,451 rows; ΔΔG + energy decomposition for the 751 covered |
| `results/phase5_statistics.md` / `.tsv` | 17 tests with effect sizes, CIs and BH-corrected q |
| `results/phase5_cmm_calibration.md` | the CMM verdict and its caveats |
| `results/LIMITATIONS.md` | twelve caveats, ordered by how much each could change a conclusion |
| `logs/09–12_*.log` | run transcripts |

## 7. Scope note

751 of 4,451 variants (17%) have structural coverage — those inside the four calibrated
constructs. The remaining 3,700 carry annotation and AlphaMissense scores but no biophysics, and
are labelled `no_structure` rather than left blank. Extending coverage would need AF3 models for
further domains, which is a request I will not prepare without asking you first.
