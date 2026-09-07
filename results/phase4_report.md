# Phase 4 report — Structure sourcing & metal setup

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED (9/9) ·
**STOPPED for the CheckMyMetal hand-off.**

Phase 4 ran entirely on data already in hand — the 11 experimental structures and the 4 AF3
models you generated. **No new external jobs were needed**, which is why I proceeded without
asking: the calibration-first plan I recommended requires nothing from you until the
CheckMyMetal step below.

---

## 1. Headline: AF3 passes calibration

CLAUDE.md requires proving AF3 reproduces known calcium geometry *before* it may be trusted on
domains with no experimental structure. It does.

| Reference type | n | Ca²⁺ deviation (median) | Max | Backbone RMSD (median) |
|---|---|---|---|---|
| **X-ray** (2W86, 1UZJ) | 20 | **0.32 Å** | 0.59 Å | 0.52 Å |
| **NMR** (1LMJ, 1EMN) | 20 | 1.25 Å | 2.24 Å | 2.19 Å |

Against crystallographic ground truth AF3 places calcium **sub-ångström** — comfortably inside
the resolution of the structures it is being compared to. The larger NMR figure is not an AF3
failure: those references have ~2 Å backbone spread of their own, and the calcium deviation
tracks that spread rather than exceeding it.

Best single models: cbEGF9 **0.09 Å**, cbEGF23 0.24 Å, cbEGF10 0.30 Å, cbEGF22 0.36 Å.

### A methodological trap worth recording

My first calibration superposed each AF3 model onto its experimental counterpart **across the
whole construct**, and produced an alarming worst case of **4.46 Å**. That number was an
artifact. These constructs contain two cbEGF domains joined by a flexible hinge, so a global fit
mixes two different quantities — how well AF3 placed the metal *inside* a domain, and how well it
guessed the *angle between* domains. The interdomain term dominates and makes correct calcium
placement look badly wrong.

Superposing **per cbEGF domain** separates them: the same models go from 4.46 Å worst-case to a
0.54 Å median. Both numbers are reported in `results/phase4_calibration.tsv`
(`construct_backbone_rmsd_A` vs `domain_backbone_rmsd_A`), because the gap between them is itself
a finding — it quantifies how much interdomain flexibility CLAUDE.md was right to warn about.

---

## 2. Tier-3 cross-check, rebuilt

AlphaFill is unavailable for FBN1 (no AlphaFold DB entry — 2,871 aa exceeds the 2,700 aa cap), so
the homology-transplant cross-check was done directly, and made stricter than AlphaFill would
have been:

> **Leave-one-out transplant.** Superpose a donor cbEGF domain onto a target domain, carry the
> donor's Ca²⁺ across, and measure how far it lands from the target's *own* experimental calcium.
> Because both structures are experimental, this is a ground-truth error bar for what transplant
> accuracy would be on a domain that has no structure of its own.

**56 donor→target transplants across 8 experimental cbEGF domains:**
median **1.44 Å**, mean 1.61 Å, 90th percentile 2.79 Å, max 3.42 Å.
Best donors: cbEGF22↔cbEGF23 0.46 Å, cbEGF10↔cbEGF13 0.70 Å.

**Conclusion: AF3 (0.32 Å vs X-ray) is roughly 4× more accurate than homology transplant
(1.44 Å).** That justifies AF3 as the primary Tier-2 method and relegates transplant to what
CLAUDE.md intended — a cross-check, not a source. It also gives us a defensible uncertainty to
quote for any domain we later model.

---

## 3. Structure selection decisions

- **1UZJ, not 1UZK,** for cbEGF22–TB4–cbEGF23. 1UZK has better resolution (1.35 Å) but only one
  of two calcium sites occupied; 1UZP and 1UZQ are apo. Resolution is the wrong criterion when
  the metal *is* the measurement.
- **1EMN** (single model) as the cbEGF32–33 representative; 1EMO is the 22-model ensemble.
- **2M74 and 5MS9** remain unassigned to named domains; 2M74 contains no calcium. Neither is used.
- All packet files were renumbered into **UniProt P35555 numbering**. This is not cosmetic:
  2W86 and 1LMJ carry local author numbering (offsets **+804** and **+1066**) while 1UZJ and 1EMN
  already coincide with UniProt. Leaving that mixture in place would have silently mislabelled
  half the CheckMyMetal output — and, being right for half the entries, would have been easy to
  miss.

---

## 4. The CheckMyMetal packet — ready for you

`handoff/cmm/inbox/` · **8 files, 16 calcium sites** · instructions in `handoff/cmm/INSTRUCTIONS.md`

| File | Role |
|---|---|
| `EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb` | ground truth / control |
| `EXP_1LMJ_cbEGF12-13.pdb` | ground truth / control |
| `EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb` | ground truth / control |
| `EXP_1EMN_cbEGF32-33.pdb` | ground truth / control |
| `AF3_cbEGF9-hyb2-cbEGF10_model0.pdb` | predicted site, compare vs 2W86 |
| `AF3_cbEGF12-13_model0.pdb` | predicted site, compare vs 1LMJ |
| `AF3_cbEGF22-TB4-cbEGF23_model0.pdb` | predicted site, compare vs 1UZJ |
| `AF3_cbEGF32-33_model0.pdb` | predicted site, compare vs 1EMN |

Experimental and predicted structures are submitted **together and in interleaved order** so CMM
judges both by its own criteria. The four `EXP_*` files are the control: if CMM flags *those*,
the problem is our file preparation, not the models.

Waters and non-calcium heteroatoms removed; hydrogens dropped; NMR ensembles collapsed to model 1.
Every file was asserted to contain Ca²⁺ before release — the gate refuses to let a metal-free
model reach CMM, since that produces confident-looking nonsense.

---

## 5. Gate

`tests/test_structures.py` — **9/9 pass**, including:

- **SIFTS identity round-trip** for all four structures using each one's own offset:
  2W86 145/147 (the 2 exceptions are non-native cloning residues), 1LMJ 86/86, 1UZJ 162/162,
  1EMN 82/82.
- Every packet file present, listed, non-trivial, and free of orphans.
- Every packet file carries Ca²⁺, with counts matching the manifest.
- Packet residues verified as **genuine** UniProt numbering against P35555, not merely claimed.
- AF3 calibration within tolerance (1.0 Å vs X-ray, 2.5 Å vs NMR).
- Transplant errors below 5.0 Å.

Self-tested against four simulated failures — a stripped calcium ion, a wrong SIFTS offset
applied to 1UZJ, an out-of-tolerance AF3 deviation, and an unlisted stray PDB in the inbox — all
four were caught.

---

## 6. 🛑 MANUAL STEP REQUIRED

I've prepared **8 files for CheckMyMetal** in `handoff/cmm/inbox/`. Please follow
`handoff/cmm/INSTRUCTIONS.md`: submit each at **https://csgid.org/csgid/metal_sites** with metal
type **calcium**, save the full result pages into `handoff/cmm/outbox/` using the same base names,
and tell me to resume.

Please keep the **full per-site parameter table**, not just the pass/fail colour — Phase 5 parses
coordination number, geometry class, bond valence, nVECSUM, gRMSD and vacancy per site.

**Phase 5 will not start until you confirm results are in the outbox.** If any file is rejected
or reports zero metal sites, tell me which and what it said rather than editing it — I would
rather regenerate the file than have a patched one enter the analysis.

---

## 7. Outputs

| Path | Contents |
|---|---|
| `results/phase4_calibration.tsv` | 40 domain-model comparisons: backbone RMSD (per-domain and per-construct), Ca deviation, pLDDT at the site |
| `results/phase4_transplant_loo.tsv` | 56 leave-one-out transplants with errors |
| `structures/ca_transplanted/provenance.json` | donors, accessions, spans, method, error summary |
| `handoff/cmm/inbox/` | 8 prepared PDBs + `manifest.csv` |
| `handoff/cmm/INSTRUCTIONS.md` | submission instructions |
| `logs/07_calibrate_*.log` | run transcript |
