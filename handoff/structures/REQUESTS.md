# Structure requests — AlphaFold-Server (AF3) · Batch 1: CALIBRATION

**Status:** ready for you · **Prepared:** 2026-08-02 · **Jobs: 4** · **Est. runtime: minutes each**

---

## Read this first: why only 4 jobs, not 43

You offered to generate "the missing domains". I've deliberately **not** asked for those yet,
for two reasons:

1. **We don't yet know which domains matter.** Which cbEGF domains carry pathogenic variants is
   the output of Phases 1–3 (variant collection → normalization → domain annotation). Generating
   all 43 now would mean ~42 pair jobs, most of which we might never use.
2. **CLAUDE.md Phase 4 requires calibration first.** Before AF3-placed calcium can be trusted for
   an un-crystallized domain, we must show AF3 reproduces *known* Ca²⁺ geometry on domains where
   an experimental answer exists. If it fails calibration, the whole Tier-2 plan changes — so
   running production models first would risk wasting your effort.

**These 4 jobs are the calibration set.** They're fully determined right now, independent of the
variant data, and they're the gate that decides whether Batch 2 is worth running. Once Phase 3
tells us which domains carry variants, I'll prepare Batch 2 as a targeted list.

---

## What to run

Everything is prepared in `handoff/structures/inbox/`:

```
inbox/af3_batch_all.json      ← all 4 jobs in one file (recommended: single upload)
inbox/af3_jobs/*.json         ← the same 4 jobs individually, if you prefer one at a time
inbox/constructs.fasta        ← the same sequences as FASTA, for reference/manual entry
inbox/manifest.csv            ← construct → range → Ca count → destination folder
```

| # | Job name | UniProt range | Length | Ca²⁺ | Compare against |
|---|---|---|---|---|---|
| 1 | `FBN1_cbEGF9-hyb2-cbEGF10_Ca2` | 807–951 | 145 aa | 2 | **2W86** (X-ray 1.80 Å) |
| 2 | `FBN1_cbEGF12-13_Ca2` | 1069–1154 | 86 aa | 2 | **1LMJ** (NMR, 25 models) |
| 3 | `FBN1_cbEGF22-TB4-cbEGF23_Ca2` | 1486–1647 | 162 aa | 2 | **1UZJ** (X-ray 2.25 Å) |
| 4 | `FBN1_cbEGF32-33_Ca2` | 2127–2205 | 79 aa | 2 | **1EMN** (NMR) |

**Why 2 calcium ions on every job:** I counted the ions in the deposited structures rather than
assuming — it is exactly **one Ca²⁺ per cbEGF domain, and zero for TB/hybrid domains**. Each
construct contains two cbEGF domains (jobs 1 and 3 also contain a hybrid/TB domain, which
contributes none). Please **do not** change these counts.

---

## Step-by-step

1. Go to **https://alphafoldserver.com** and sign in (free Google account; non-commercial /
   academic use only, which is what this is).

2. **Preferred — upload the batch file.** Use the JSON upload option and give it
   `handoff/structures/inbox/af3_batch_all.json`. That queues all four jobs at once with the
   correct sequences and ion counts, and removes any chance of a transcription error.

   *If the UI has changed and JSON upload isn't obvious,* fall back to manual entry:
   - Add **one protein chain**, paste the sequence from `inbox/constructs.fasta` (copy the
     sequence lines only, not the `>` header).
   - Add an **ion** entity, choose **Ca²⁺ (CA)**, set **count = 2**.
   - Name the job exactly as in the table above (the names are how I match results to requests).
   - Repeat for each of the four.

3. **Leave the seed field empty / default.** The server will pick one and report it — that's
   fine, but see step 5: I need it recorded.

4. When each job finishes, **download the results zip**.

5. **Where to put them** — create one folder per job under `structures/alphafold/`, named
   exactly after the job, and unzip into it:

   ```
   structures/alphafold/FBN1_cbEGF9-hyb2-cbEGF10_Ca2/
   structures/alphafold/FBN1_cbEGF12-13_Ca2/
   structures/alphafold/FBN1_cbEGF22-TB4-cbEGF23_Ca2/
   structures/alphafold/FBN1_cbEGF32-33_Ca2/
   ```

   Keep **everything** in the zip — the `.cif` model files, the confidence/summary JSONs, and any
   `terms_of_use` or job-definition file. The confidence JSON carries pLDDT and PAE, and the job
   definition carries **the model seed**, which I need for the reproducibility record. Don't
   rename the files inside.

   From bash, that's roughly:
   ```bash
   cd ~/research/fbn1-marfan
   mkdir -p structures/alphafold/FBN1_cbEGF12-13_Ca2
   unzip ~/Downloads/<downloaded>.zip -d structures/alphafold/FBN1_cbEGF12-13_Ca2/
   ```

6. **Tell me "AF3 results are in"** and I'll take it from there.

---

## What I'll do when they land

1. Verify each model actually contains 2 Ca²⁺ ions (a metal-free model cannot be validated, and
   CLAUDE.md forbids sending one to CheckMyMetal).
2. Superpose each AF3 model onto its experimental counterpart and report backbone RMSD **and**
   the Ca²⁺-position deviation — the number that actually decides whether AF3 places calcium
   correctly.
3. Build the CheckMyMetal packet (`handoff/cmm/inbox/`) containing both the experimental sites and
   the AF3-placed sites, so CMM scores them side by side. That's the calibration result, and it
   goes in the methods.
4. Report agreement to you before any un-crystallized domain is modelled.

---

## Notes / things worth knowing

- **AF-Server has a daily job limit.** I could not verify the current number (the FAQ is a
  JavaScript page I can't read reliably), but 4 jobs is well inside any published limit. If you
  hit a cap, run them in the order listed — job 2 (`cbEGF12-13`) and job 4 (`cbEGF32-33`) are the
  most informative, being clean 2-domain pairs with unambiguous experimental Ca sites.
- **Sequences are verified, not transcribed.** They were sliced from the cached UniProt P35555
  record by `scripts/04_prepare_af3_requests.py` and checked residue-by-residue against the
  deposited coordinates (1LMJ 86/86, 1UZJ 162/162, 1EMN 79/79, 2W86 145/145 native residues).
- **Modelled as wild-type.** Variant models come later, once Phase 3 identifies the variants worth
  modelling.
- **AF3 output is a modelling assumption, not a measurement.** Calcium placed by AF3 will be
  documented as such in `results/LIMITATIONS.md`, distinct from experimentally observed calcium.
- If a job fails or the server rejects an input, **don't hand-edit the sequence** — tell me what
  it said and I'll regenerate the file.

---

## Not requested (recorded so nothing is silently dropped)

- **Batch 2 — production models** for cbEGF domains carrying variants: blocked on Phase 3.
- **AlphaFill cross-check:** unavailable for FBN1. AlphaFill derives from AlphaFold DB, and FBN1
  has no AlphaFold DB entry (2871 aa, above the 2700 aa cutoff). Proposed replacement is a manual
  superposition transplant from an experimental Ca²⁺ site — awaiting your decision (Phase 0
  report §6, question 2).
- **2M74 (45–178) and 5MS9 (113–287):** present in the PDB, not yet assigned to named domains, and
  2M74 contains no calcium. No action needed now.


---
---

# Structure requests — AlphaFold-Server (AF3) · Batches 2 and 3

**Status:** ready for you · **Prepared:** 2026-08-09 · **Jobs: 30** (21 + 9) · **Authorised by you on 2026-08-08**

Calibration (Batch 1) passed: AF3-placed calcium sits a median 0.32 Å from the X-ray position and
CheckMyMetal cannot distinguish the modelled sites from the experimental ones. That unlocks the two
things this batch asks for.

## Why these two batches

**Batch 2 closes the coverage gap.** 2,897 variants sit inside a cbEGF domain; only 524 currently
have any structure, because only 8 of the 43 domains were ever solved. Everything biophysical in
the manuscript — ΔΔG, solvent accessibility, calcium distance, coordination — is computed on that
17% and extrapolated in prose to the rest. These 21 constructs put all 43 domains on the same
footing.

**Batch 3 measures the mechanism instead of inferring it.** The manuscript's central claim is that
calcium-ligand variants damage the metal site while leaving the fold intact. Wild-type models cannot
test that. These 9 jobs fold the *mutant* sequence with the same calcium stoichiometry and ask
whether AF3 can still place the ion, and where. Two are deliberate controls whose outcome would
falsify the claim if they behave unexpectedly.

## Batch 2 — coverage for all 43 cbEGF domains (21 jobs)

Constructs are three consecutive cbEGF domains, stepping by two, so **every domain except cbEGF1
appears at least once with its N-terminal neighbour** — calcium affinity in this family is modulated
by N-terminal domain linkage, so an isolated domain is the wrong construct. Any TB or hybrid domain
lying between two cbEGF domains is included, exactly as in the calibrated 2W86 and 1UZJ fragments.
One Ca²⁺ per cbEGF domain, so three ions per job.

Files: `inbox/batch2_domains/af3_batch2_all.json` (single upload) or
`inbox/batch2_domains/af3_jobs/*.json` (one at a time).

| # | Job name | UniProt range | Length | Ca²⁺ | Domains |
|---|---|---|---|---|---|
| 1 | `FBN1_cbEGF1-3_Ca3` | 246-529 | 284 aa | 3 | cbEGF 1,2,3 |
| 2 | `FBN1_cbEGF3-5_Ca3` | 490-612 | 123 aa | 3 | cbEGF 3,4,5 |
| 3 | `FBN1_cbEGF5-7_Ca3` | 572-764 | 193 aa | 3 | cbEGF 5,6,7 |
| 4 | `FBN1_cbEGF7-9_Ca3` | 723-846 | 124 aa | 3 | cbEGF 7,8,9 |
| 5 | `FBN1_cbEGF9-11_Ca3` | 807-1069 | 263 aa | 3 | cbEGF 9,10,11 |
| 6 | `FBN1_cbEGF11-13_Ca3` | 1028-1154 | 127 aa | 3 | cbEGF 11,12,13 |
| 7 | `FBN1_cbEGF13-15_Ca3` | 1113-1237 | 125 aa | 3 | cbEGF 13,14,15 |
| 8 | `FBN1_cbEGF15-17_Ca3` | 1197-1321 | 125 aa | 3 | cbEGF 15,16,17 |
| 9 | `FBN1_cbEGF17-19_Ca3` | 1280-1403 | 124 aa | 3 | cbEGF 17,18,19 |
| 10 | `FBN1_cbEGF19-21_Ca3` | 1363-1486 | 124 aa | 3 | cbEGF 19,20,21 |
| 11 | `FBN1_cbEGF21-23_Ca3` | 1446-1647 | 202 aa | 3 | cbEGF 21,22,23 |
| 12 | `FBN1_cbEGF23-25_Ca3` | 1606-1807 | 202 aa | 3 | cbEGF 23,24,25 |
| 13 | `FBN1_cbEGF25-27_Ca3` | 1766-1890 | 125 aa | 3 | cbEGF 25,26,27 |
| 14 | `FBN1_cbEGF27-29_Ca3` | 1849-1972 | 124 aa | 3 | cbEGF 27,28,29 |
| 15 | `FBN1_cbEGF29-31_Ca3` | 1930-2054 | 125 aa | 3 | cbEGF 29,30,31 |
| 16 | `FBN1_cbEGF31-33_Ca3` | 2013-2205 | 193 aa | 3 | cbEGF 31,32,33 |
| 17 | `FBN1_cbEGF33-35_Ca3` | 2166-2290 | 125 aa | 3 | cbEGF 33,34,35 |
| 18 | `FBN1_cbEGF35-37_Ca3` | 2247-2443 | 197 aa | 3 | cbEGF 35,36,37 |
| 19 | `FBN1_cbEGF37-39_Ca3` | 2402-2523 | 122 aa | 3 | cbEGF 37,38,39 |
| 20 | `FBN1_cbEGF39-41_Ca3` | 2485-2606 | 122 aa | 3 | cbEGF 39,40,41 |
| 21 | `FBN1_cbEGF41-43_Ca3` | 2567-2687 | 121 aa | 3 | cbEGF 41,42,43 |

## Batch 3 — mutant cofolding (9 jobs)

Every mutant is built on a **calibrated** construct, so each has both an experimental reference and
an AF3 wild-type counterpart already in hand. Wild-type residues were verified against the P35555
reference before the files were written.

Files: `inbox/batch3_mutants/af3_batch3_all.json` or `inbox/batch3_mutants/af3_jobs/*.json`.

| # | Job name | Mutation | Reference | Why it is in the set |
|---|---|---|---|---|
| 1 | `FBN1_cbEGF32-33_N2144S_Ca2` | N2144S | 1EMN | EXPERIMENTAL ANCHOR. |
| 2 | `FBN1_cbEGF32-33_N2183S_Ca2` | N2183S | 1EMN | Second literature mutant (the equivalent position in cbEGF33), characterised in the same domain-context work. |
| 3 | `FBN1_cbEGF32-33_E2130K_Ca2` | E2130K | 1EMN | Pathogenic ca_EH ligand, buried (RSA 7. |
| 4 | `FBN1_cbEGF12-13_D1113G_Ca2` | D1113G | 1LMJ | Pathogenic ca_D1 ligand. |
| 5 | `FBN1_cbEGF12-13_E1073K_Ca2` | E1073K | 1LMJ | Pathogenic ca_EH ligand, buried (RSA 9. |
| 6 | `FBN1_cbEGF22-TB4-cbEGF23_D1487G_Ca2` | D1487G | 1UZJ | Pathogenic ca_D1 ligand in the TB4-flanked construct — tests a different domain context. |
| 7 | `FBN1_cbEGF9-hyb2-cbEGF10_E913K_Ca2` | E913K | 2W86 | Pathogenic ca_EH ligand, the most buried in the set (RSA 2. |
| 8 | `FBN1_cbEGF12-13_C1138F_Ca2` | C1138F | 1LMJ | POSITIVE CONTROL for the other mechanism: removes disulfide C4. |
| 9 | `FBN1_cbEGF12-13_P1141L_Ca2` | P1141L | 1LMJ | NEGATIVE CONTROL: non-critical residue, AlphaMissense 0. |

**`N2144S` is the one that matters most.** Its calcium affinity has been measured — roughly 5-fold
reduced in isolated cbEGF32 and 9-fold in the TB6–cbEGF32 pair. It is the only variant in this study
where a modelled answer can be checked against an experimental one.

## What to do

1. Upload `af3_batch2_all.json` at <https://alphafoldserver.com>. If the server caps you below 21
   jobs a day, run them in listed order and finish the rest the next day.
2. Upload `af3_batch3_all.json` the same way.
3. For each finished job, download the result and unzip it into the folder named in the batch's
   `manifest.csv` `destination` column, i.e. `structures/alphafold/<job_name>/`. Keep the whole
   output folder — I need the CIF, the confidence JSON **and the seed the server chose**, which is
   what makes the run reproducible.
4. Tell me when a batch is done. I will not touch anything in `structures/alphafold/` until you do.

If a job is rejected, **do not hand-edit the sequence** — send me the error and I will regenerate.

## What I will do with them

- Re-run the placement calibration on the constructs that overlap the solved domains (cbEGF9–11,
  11–13, 21–23, 31–33 all contain calibration domains) — free extra validation of the 0.32 Å figure
  on constructs that were not used to derive it.
- Extend every structural metric to all 43 domains, raising structural coverage from 524 to
  approximately 2,897 cbEGF variants.
- For Batch 3, measure ion displacement, ligand retention and bond-valence loss in each mutant
  against its wild-type counterpart, and assemble a CheckMyMetal packet of wild-type/mutant pairs.

## Still outstanding from me to you (not AF3)

- **Four paywalled PDFs** for `resources/papers/`: Kettle *et al.* 1999 *J Mol Biol* (N2144S calcium
  binding); Whiteman & Handford 2000 *Hum Mol Genet* 9:1987 (domain context); Reinhardt *et al.*
  2000 *JBC* (proteolytic susceptibility); Handford *et al.* 1995 *JBC* (cbEGF calcium affinity).
  These carry the measured Kd values that anchor the whole analysis to experiment.
- **A Rosetta or PyRosetta academic licence** (free, registration only), if you want the ΔΔG
  repeated in a force field that actually models metal coordination.


---

## Batch 4 — mutant cofolding to make Figure 4 a rate (2026-08-16)

60 jobs in `handoff/structures/inbox/batch4_mutants/`, split into 2 day files of exactly 30. Full instructions in that folder's `INSTRUCTIONS.md`; per-job detail in its `manifest.csv`.

- **ca_ligand_pathogenic** — 34 jobs
- **negative_control** — 14 jobs
- **composition_control** — 12 jobs

Every wild-type residue verified against P35555 before writing; every sequence checked to differ from the reference at exactly one position.
