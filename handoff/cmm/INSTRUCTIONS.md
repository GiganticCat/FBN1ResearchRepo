# CheckMyMetal — instructions

**Prepared:** 2026-08-02 ·
**8 files · 16 calcium sites total**

Everything is in `handoff/cmm/inbox/`. Please put results in `handoff/cmm/outbox/`.

---

## What this is for

We need CheckMyMetal's own verdict on whether the calcium sites are geometrically believable —
both in the **experimental** structures (which should pass, and act as our control) and in the
**AF3-predicted** models (which is what we actually want to validate).

Distance measurements already say AF3 places calcium a median of **0.32 Å** from the
crystallographic position. CMM tests something different and stricter: coordination number,
geometry class, bond valence, and whether the site looks like a real Ca²⁺ site at all.

---

## Steps

1. Go to **https://csgid.org/csgid/metal_sites** — select metal type **calcium (CA)**.

2. Upload the files **in the order below**, so each AF3 model is scored right after its
   experimental counterpart and the two are easy to compare:

| # | File | What it is |
|---|---|---|
| 1 | `AF3_cbEGF12-13_model0.pdb` | predicted site, compare against 1LMJ |
| 2 | `EXP_1LMJ_cbEGF12-13.pdb` | ground truth for calibration |
| 3 | `AF3_cbEGF22-TB4-cbEGF23_model0.pdb` | predicted site, compare against 1UZJ |
| 4 | `EXP_1UZJ_cbEGF22-TB4-cbEGF23.pdb` | ground truth for calibration |
| 5 | `AF3_cbEGF32-33_model0.pdb` | predicted site, compare against 1EMN |
| 6 | `EXP_1EMN_cbEGF32-33.pdb` | ground truth for calibration |
| 7 | `AF3_cbEGF9-hyb2-cbEGF10_model0.pdb` | predicted site, compare against 2W86 |
| 8 | `EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb` | ground truth for calibration |

3. For each file, save the CMM output into `handoff/cmm/outbox/` using **the same base name**
   with the CMM result suffix — e.g. `EXP_2W86_cbEGF9-hyb2-cbEGF10_cmm.html` (or `.csv`/`.txt`,
   whichever the server offers). Keep the full result page, not a screenshot: I need the
   per-site parameter table, not just the pass/fail colour.

4. Tell me **"CMM results are in"** and I will parse them.

---

## Things worth knowing

- **All residues are numbered in UniProt P35555 numbering**, not the PDB's own. I renumbered them
  deliberately: 2W86 and 1LMJ use local numbering (offset +804 and +1066) while 1UZJ and 1EMN
  already match UniProt. Leaving that mixture in place would have silently mislabelled half the
  results. So a residue CMM calls 1070 really is P35555 residue 1070.
- Waters and non-calcium heteroatoms were removed; hydrogens dropped. Only protein + Ca²⁺ remain.
- NMR entries (1LMJ, 1EMN) were collapsed to **model 1**. CMM expects a single model.
- If CMM rejects a file or reports **zero** metal sites, don't try to fix the file — tell me
  which one and what it said. A file that silently loses its calcium would corrupt the
  calibration, so I would rather regenerate it than have it patched.
- Every submitted file is guaranteed to contain at least one Ca²⁺ ion; this was asserted by the
  Phase 4 gate before the packet was released.

## Expected outcome

The four `EXP_*` files should come back clean — they are deposited, refined structures. They are
the control: if CMM flags **those**, the problem is our file preparation, not the models. The
four `AF3_*` files are the real test.
