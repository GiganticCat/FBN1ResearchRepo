# AlphaFold Server — batch 4 (mutant cofolding)

Prepared 2026-08-16 by `scripts/38_af3_batch4_mutants.py`. **60 jobs**, split into 2 day files of exactly 30 — one full AlphaFold Server daily allowance each.

## Why this batch exists

Figure 4 currently rests on nine variants and can show only that Ca²⁺-site disruption happens, not how often. These jobs convert it into a rate, with the controls needed to make that rate mean something.

| class | n | what it is for |
|---|---|---|
| `ca_ligand_pathogenic` | 34 | every pathogenic side-chain Ca²⁺ ligand with structural coverage that is not already folded |
| `negative_control` | 14 | benign / unlabelled non-ligand positions in the same constructs — measures the noise floor on the same folds |
| `composition_control` | 12 | D/E/N residues that coordinate nothing — separates residue type from calcium role |

## What to do

1. Go to <https://alphafoldserver.com>.
2. Upload `af3_batch4_day1.json`. Each day file holds exactly 30 jobs — a full daily allowance — and carries a mix of all three classes, so a day that runs alone is still analysable.
3. Repeat with `af3_batch4_day2.json` (and day3 if present) on subsequent days.
4. For each finished job, download the result and unzip the **whole output folder** into the path in `manifest.csv` → `destination`, i.e. `structures/alphafold/<job_name_lowercased>/`.
   Keep the CIF, the confidence JSON **and the seed the server chose** — the seed is what makes the run reproducible, and the five-model spread is what every significance claim in Figure 4 is measured against.
5. Tell me when a day is done. I will not touch `structures/alphafold/` until you do.

If a job is rejected, **do not hand-edit the sequence** — send me the error and I will regenerate it. Every sequence here was checked against P35555 before writing: correct length, exactly one substitution, at the intended position, with the intended residues.

## What I will do with them

- Measure ΔBVS, coordination number and ion displacement per mutant against its existing wild-type fold, exactly as in Figure 4.
- Report the **fraction** of pathogenic Ca²⁺-ligand variants whose site disruption clears the seed-noise floor, with the controls establishing that floor empirically rather than by assumption.
- Test whether disruption predicts anything else we hold — clinical severity, AlphaMissense, burial, which ligand position is hit.

Total sequence to fold: 9,657 residues across 17 constructs; longest 284 aa.
