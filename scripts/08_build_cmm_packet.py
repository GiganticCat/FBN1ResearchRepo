#!/usr/bin/env python3
"""
08_build_cmm_packet.py — PHASE 4. Assemble the CheckMyMetal submission packet.

CheckMyMetal (CMM) is a web service with no usable programmatic API here, so this script only
PREPARES the inputs. The mentor runs them. Per CLAUDE.md 3b the pipeline stops after this.

Packet contents: 8 PDB files — the 4 experimental Ca-bound structures (ground truth) and the
4 AF3 Ca-cofolded models (top-ranked model_0), so CMM scores predicted and experimental sites
side by side and the calibration is judged by CMM's own criteria rather than by distance alone.

Every file is renumbered into **UniProt P35555 numbering** before writing, so CMM's output maps
straight onto the variant table with no offset arithmetic. This matters: 2W86 and 1LMJ use local
author numbering (offsets +804 and +1066) while 1UZJ and 1EMN happen to match UniProt already —
a mixture that would otherwise silently mislabel half the results.

Writes:
  handoff/cmm/inbox/*.pdb          8 prepared structures
  handoff/cmm/inbox/manifest.csv   filename -> domain -> site -> source/method
  handoff/cmm/INSTRUCTIONS.md      numbered instructions for the mentor
"""

from __future__ import annotations

import csv
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
INBOX = ROOT / "handoff/cmm/inbox"
CMM_URL = "https://csgid.org/csgid/metal_sites"

# (key, pdb, author->UniProt offset, chain, construct span, label, cbEGF sub-domains)
EXPERIMENTAL = [
    ("2w86", 804, "A", (807, 951), "cbEGF9-hyb2-cbEGF10", ["cbEGF9", "cbEGF10"]),
    ("1lmj", 1066, "A", (1069, 1154), "cbEGF12-13", ["cbEGF12", "cbEGF13"]),
    ("1uzj", 0, "A", (1486, 1647), "cbEGF22-TB4-cbEGF23", ["cbEGF22", "cbEGF23"]),
    ("1emn", 0, "A", (2127, 2205), "cbEGF32-33", ["cbEGF32", "cbEGF33"]),
]
AF3 = {
    "fold_fbn1_cbegf9_hyb2_cbegf10_ca2": ((807, 951), "cbEGF9-hyb2-cbEGF10", "2W86"),
    "fold_fbn1_cbegf12_13_ca2":          ((1069, 1154), "cbEGF12-13", "1LMJ"),
    "fold_fbn1_cbegf22_tb4_cbegf23_ca2": ((1486, 1647), "cbEGF22-TB4-cbEGF23", "1UZJ"),
    "fold_fbn1_cbegf32_33_ca2":          ((2127, 2205), "cbEGF32-33", "1EMN"),
}

parser = MMCIFParser(QUIET=True)
io = PDBIO()


class KeepProteinAndCa(Select):
    """Drop waters and any non-calcium heteroatoms; CMM should see the metal site, not solvent."""
    def accept_residue(self, residue):
        het = residue.id[0]
        if het == " ":
            return 1
        return 1 if residue.get_resname().strip() == "CA" else 0

    def accept_atom(self, atom):
        return 0 if atom.element == "H" else 1


def count_ca(structure) -> int:
    return sum(1 for ch in structure for r in ch
               if r.id[0].startswith("H_") and r.get_resname().strip() == "CA")


def main() -> int:
    INBOX.mkdir(parents=True, exist_ok=True)
    for old in INBOX.glob("*.pdb"):
        old.unlink()
    rows = []
    print("Building CheckMyMetal packet\n")

    # ---- experimental ------------------------------------------------------
    for pdb, off, chain, span, label, subs in EXPERIMENTAL:
        st = parser.get_structure(pdb, str(ROOT / f"structures/pdb/{pdb}.cif"))
        model = st[0]                      # first model (NMR ensembles collapse to model 1)
        for ch in list(model):
            if ch.id != chain:
                model.detach_child(ch.id)
        # renumber protein residues into UniProt numbering
        for res in model[chain]:
            if res.id[0] == " ":
                res.id = (" ", res.id[1] + off, " ")
        n_ca = count_ca(model)
        out = INBOX / f"EXP_{pdb.upper()}_{label}.pdb"
        io.set_structure(model)
        io.save(str(out), select=KeepProteinAndCa())
        rows.append({"filename": out.name, "source": "experimental", "method": "X-ray/NMR",
                     "pdb_or_job": pdb.upper(), "construct": label,
                     "domains": "+".join(subs), "uniprot_range": f"{span[0]}-{span[1]}",
                     "n_ca_sites": n_ca, "numbering": "UniProt P35555",
                     "author_offset_applied": off,
                     "role": "ground truth for calibration",
                     "intended_variants": "n/a (wild-type reference)"})
        print(f"  {out.name:52} {n_ca} Ca  (offset {off:+d})")

    # ---- AF3 top-ranked models --------------------------------------------
    for job, (span, label, ref) in AF3.items():
        cif = ROOT / "structures/alphafold" / job / f"{job}_model_0.cif"
        st = parser.get_structure(job, str(cif))
        model = st[0]
        for res in model["A"]:
            if res.id[0] == " ":
                res.id = (" ", span[0] + res.id[1] - 1, " ")
        n_ca = count_ca(model)
        out = INBOX / f"AF3_{label}_model0.pdb"
        io.set_structure(model)
        io.save(str(out), select=KeepProteinAndCa())
        rows.append({"filename": out.name, "source": "AF3 (AlphaFold-Server)",
                     "method": "cofolded with 2 Ca2+", "pdb_or_job": job, "construct": label,
                     "domains": label, "uniprot_range": f"{span[0]}-{span[1]}",
                     "n_ca_sites": n_ca, "numbering": "UniProt P35555",
                     "author_offset_applied": 0,
                     "role": f"predicted site, compare against {ref}",
                     "intended_variants": "n/a (wild-type model)"})
        print(f"  {out.name:52} {n_ca} Ca  (ref {ref})")

    with (INBOX / "manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    total_sites = sum(r["n_ca_sites"] for r in rows)
    (ROOT / "handoff/cmm/INSTRUCTIONS.md").write_text(f"""# CheckMyMetal — instructions

**Prepared:** {datetime.now(timezone.utc).date().isoformat()} ·
**{len(rows)} files · {total_sites} calcium sites total**

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

1. Go to **{CMM_URL}** — select metal type **calcium (CA)**.

2. Upload the files **in the order below**, so each AF3 model is scored right after its
   experimental counterpart and the two are easy to compare:

| # | File | What it is |
|---|---|---|
{chr(10).join(f"| {i} | `{r['filename']}` | {r['role']} |" for i, r in enumerate(sorted(rows, key=lambda x: (x['construct'], x['source'])), 1))}

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
""", encoding="utf-8")

    print(f"\n  {len(rows)} files, {total_sites} Ca sites")
    print("  wrote handoff/cmm/inbox/manifest.csv and handoff/cmm/INSTRUCTIONS.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
