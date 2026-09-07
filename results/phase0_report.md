# Phase 0 report — Environment, scaffold & literature intake

**Date:** 2026-08-02 (UTC) · **Status:** complete, gate PASSED · **Awaiting your approval before Phase 1.**

Phase 0 gate `tests/test_env.py`: **8/8 checks pass**.
Env check `scripts/00_env_check.py`: **38 PASS / 5 WARN / 0 FAIL**.

Nothing here is blocking Phase 1. Four findings need your attention, and **two of them change
what later phases can do** — please read §4 before approving.

---

## 1. What was produced

| Path | Contents |
|---|---|
| `scripts/00_env_check.py` | Environment/API verification; writes the log, JSON and tools manifest |
| `tests/test_env.py` | **Phase 0 gate** — 8 assertions, self-tested against 4 simulated failures |
| `resources/reference/lit_notes.md` | One entry per PDF (**19**) + Appendix A cross-check |
| `manifest/tools.md` | Tool versions, paths, DB access date, P35555 sequence MD5 |
| `logs/00_env_check_*.log`, `logs/00_env_check_latest.json` | Run transcript + machine-readable evidence |
| `requirements.txt` | Re-pinned via `pip freeze --local` (added `pypdf==6.14.2`) |
| `handoff/structures/` | Created — was missing from the scaffold |

---

## 2. Tool versions

| Tool | Version | Path |
|---|---|---|
| Python | 3.14.4 | `.venv/bin/python` (venv active) |
| PyMOL | 3.1.0 | `/usr/bin/pymol` |
| DSSP (`mkdssp`) | 4.5.7 | `/usr/bin/mkdssp` |
| FoldX | 5.1 | `/home/harvey/.local/bin/foldx` → `~/tools/foldx/foldx_20261231` |
| FoldX `molecules/` | present (15 files) | `~/tools/foldx/molecules` |

All three native tools resolve via `shutil.which` **and execute** — the gate runs each one, since
a dangling symlink resolves but cannot run. FoldX has no `--version` flag; its banner is parsed
instead and its non-zero exit is expected.

Python packages verified by import: pandas 3.0.5, numpy 2.5.1, scipy 1.18.0, requests 2.34.2,
biopython 1.87, biotite 1.7.1, ProDy 2.6.1, freesasa 2.2.1, matplotlib 3.11.1, seaborn 0.13.2,
statsmodels 0.14.6, python-dotenv 1.2.2, pypdf 6.14.2.

`.env` supplies `NCBI_API_KEY` and `NCBI_EMAIL` (presence and length only — values were never
read into any log, and `.env` is confirmed git-ignored).

---

## 3. Database reachability (live test queries, 2026-08-02)

| Source | Status | Evidence returned |
|---|---|---|
| UniProt P35555 | ✅ | 2871 aa as expected; sequence MD5 `f8f84ac9dc30f3ff640ad4f691fee2b1` |
| ClinVar (NCBI E-utilities) | ✅ | **9,475** records for `FBN1[gene]` — *corrected in Phase 1, see note* |
| PDBe SIFTS | ✅ | 11 PDB entries mapped to P35555 |
| gnomAD GraphQL | ✅ | FBN1 = ENSG00000166147 |
| Mutalyzer | ✅ | `NM_000138.5:c.364C>T` → `p.(Arg122Cys)` |
| AlphaFold DB | ⚠️ **no entry for FBN1** | see §4.2 |
| AlphaFill | ⚠️ no entry for FBN1 | status `unknown`; depends on AlphaFold DB |
| VariantValidator | ⚠️ unreachable | connection refused; Mutalyzer covers Phase 2 |

Every check confirms the **response shape**, not just HTTP 200 — e.g. UniProt is checked for
accession *and* length 2871, gnomAD for the returned gene symbol.

> **Correction (made during Phase 1).** This table originally reported *13,501 records for Gene
> ID 2200*. **That number was wrong.** The query used `2200[gene_id]`, and `[gene_id]` is not a
> field in the ClinVar E-utilities index — NCBI silently rewrote it to `2200[All Fields]`, a
> free-text match that returned records for unrelated genes (e.g. COL2A1 `c.2200T>C`, ALK
> `c.2200T>C`). The correct query is `FBN1[gene]`, giving **9,475**, which reconciles with the
> bulk pull. `scripts/00_env_check.py` now uses `FBN1[gene]` and **fails the check if the
> returned `querytranslation` degrades to `[All Fields]`**, so this class of silent
> mistranslation cannot recur. No downstream artifact used the bad number — it appeared only in
> this reachability table.

---

## 4. Findings that need your attention

### 4.1 ✅ RESOLVED — literature gap closed

**Superseded by your uploads on 2026-08-02.** Recorded here because it changed the Phase 3 plan
mid-flight.

You topped the folder up twice during this phase: five Appendix-A gaps at 09:25 (Rao 1995,
Downing 1996, Handford 1995, Jensen, gnomAD v4 constraint), then the two wrong files at 09:35.
Both replacements are **verified correct** — `ClinGen_FBN1_VCEP_2024_GenomeMed.pdf` now contains
80 FBN1 / 51 VCEP mentions, and `Godwin_2023_NSMB_microfibril_cryoEM.pdf` now contains 191
fibrillin / 222 microfibril / 23 cbEGF mentions. The gate caught both changes on re-run.

**Phase 3 will therefore be aligned to the actual ClinGen FBN1 VCEP specification**, with
Baudhuin 2019 as the corroborating operational taxonomy rather than the substitute.

Two residual items, neither blocking:

- `Jensen_2009_Structure_2W86.pdf` is the **2012 review** (Jensen, Robertson & Handford,
  *Structure* 20:215–225), not the 2009 primary 2W86 paper. It is genuinely useful — it supplied
  the authoritative PDB↔domain mapping in §5 — so only the primary structural citation is missing.
- `Handford_1995_JBC_cbEGF.pdf` is a **scanned image PDF with no text layer**. It was read
  visually rather than parsed. Any future automated pass must not treat its empty text extraction
  as an empty paper.

**⚠️ Numbering trap in the 1990s papers.** Handford 1995 numbers residues per Pereira et al.
(1993), not UniProt — e.g. its headline mutation "Asn-2144→Ser". Positions quoted in Rao 1995,
Downing 1996 and Handford 1995 must be **re-mapped and WT-identity-checked against P35555** before
entering any table. I will not copy them across directly.

Two Windows `:Zone.Identifier` artifacts (25-byte ADS markers from the `/mnt/c` copy) were left
alongside the new PDFs; I confirmed their contents and removed them.

**⚠️ Numbering trap in the newly added 1990s papers.** Handford 1995 numbers residues per Pereira
et al. (1993), not UniProt — e.g. its headline mutation "Asn-2144→Ser". Positions quoted in Rao
1995, Downing 1996 and Handford 1995 must be **re-mapped and WT-identity-checked against P35555**
before entering any table. I will not copy them across directly.

### 4.2 🔴 FBN1 has **no AlphaFold DB model** — this contradicts CLAUDE.md §4

`https://alphafold.ebi.ac.uk/api/prediction/P35555` returns **404**. This is not a transient
outage or a moved endpoint — I verified it three ways:

- A control accession under the size cap (**P02751**, fibronectin, 2477 aa) returns **200** from
  the same endpoint, so the API is live and correct.
- **UniProt carries no AlphaFoldDB cross-reference for P35555** (the control does).
- The AlphaFold DB 2024 paper in `resources/papers/` states the exclusion rule outright:
  sequences **>2700 residues** for Swiss-Prot/proteome entries are not covered. **FBN1 is 2871 aa.**
  Other long human proteins behave identically (P98160 4391 aa → 404; Q8WZ42 titin → 404).

**Consequences for Phase 4:**
- **Tier 1 (experimental) is unaffected** — 11 PDB entries are available (§5).
- **Tier 2 (AF3 cofolding with Ca²⁺) is unaffected** — it works on *domain sequences* you submit
  to AlphaFold-Server, not on an AFDB lookup. This remains the primary model source.
- **Tier 3 (AlphaFill cross-check) cannot be fetched as written.** AlphaFill is built on top of
  AlphaFold DB entries, so with no AFDB model there is no AlphaFill entry for FBN1 (confirmed
  live: status `unknown`). A cross-check would need either a per-domain model uploaded to
  AlphaFill (I have **not** verified that it accepts user uploads — I won't assume) or the manual
  superposition transplant from an experimental Ca site, which CLAUDE.md already describes and
  which we can do locally with the donors in §5.
- We also lose AFDB per-residue **pLDDT** for FBN1; confidence will have to come from whatever
  AF3/AF-Server returns per domain.

Since CLAUDE.md already forbids using the full-length model for geometry ("work per cbEGF domain,
never the full-length model"), the practical loss is limited to the Tier-3 cross-check and pLDDT.
**Decision needed** — see §6.

### 4.3 🟡 AlphaMissense has two wild-type residues at position 472

The subset has 54,568 rows, 19 more than the expected 2871 × 19 = 54,549. Cause: **position 472
carries both `C` and `Y` as wild type** (38 rows instead of 19). No duplicate variant strings, no
NaNs, all positions in 1–2871, all rows are P35555.

**UniProt P35555 residue 472 is `Y`.** The `C472*` rows correspond to a different reference
allele and must be dropped. Phase 2's WT-identity check will catch this **only if the join
asserts residue identity rather than matching on position alone** — which is exactly what
CLAUDE.md §5 requires. Flagged so it cannot silently double-count.

Also worth noting for parsing: the subset kept only the copyright comment line, **not** the real
column header, so columns must be named explicitly as
`uniprot_id, protein_variant, am_pathogenicity, am_class`.

### 4.4 🟡 VariantValidator unreachable; Mutalyzer is the Phase 2 validator

`rest.variantvalidator.org` refuses connections from this host (all paths, repeated attempts).
**Mutalyzer works** and correctly returns `NM_000138.5(NP_000129.3):p.(Arg122Cys)` for
`NM_000138.5:c.364C>T`. CLAUDE.md §4 names either tool, so Phase 2 will use Mutalyzer as primary.
The gate requires **at least one** working validator and fails if both are down.

---

## 5. Structure inventory (groundwork for Phase 4)

Enumerated live from the PDBe SIFTS API rather than a static list, as CLAUDE.md requires.
**11 entries** map to P35555:

| PDB | Method | Resolution | UniProt range | Domain assignment |
|---|---|---|---|---|
| 2M74 | NMR | — | 45–178 | not in Appendix B — new candidate |
| 5MS9 | NMR | — | 113–287 | not in Appendix B — new candidate |
| 2W86 | X-ray | 1.80 Å | 807–951 | ✅ cbEGF9–hyb2–cbEGF10 |
| 1LMJ | NMR | — | 1069–1154 | ✅ cbEGF12–13 |
| 1UZJ | X-ray | 2.25 Å | 1486–1647 | ✅ **cbEGF22–TB4–cbEGF23** (confirmed, Jensen 2012) |
| 1UZK | X-ray | 1.35 Å | 1486–1647 | same fragment — **highest resolution of the four** |
| 1UZP | X-ray | 1.78 Å | 1486–1647 | same fragment |
| 1UZQ | X-ray | 2.40 Å | 1486–1647 | same fragment |
| 1APJ | NMR | — | 2052–2125 | ✅ **isolated TB6** (confirmed, Jensen 2012) |
| 1EMN | NMR | — | 2124–2205 | ✅ cbEGF32–33 |
| 1EMO | NMR | — | 2124–2205 | ✅ cbEGF32–33 |

The Jensen 2012 review you supplied resolves what was previously an inference: it states the
mapping explicitly, confirming **1UZJ = cbEGF22–TB4–cbEGF23** and identifying **1APJ as TB6**
(an entry Appendix B did not name). All of Appendix B's anchors are present.

**Not yet verified:** which of these actually contain Ca²⁺ ions. Appendix B asserts it for the
anchors; I have not downloaded coordinates, so I am not repeating the claim as fact. Phase 4 will
check each file directly before treating any entry as a Ca²⁺ donor.

Coverage remains the expected problem: these span roughly 8 of 43 cbEGF domains, so most variants
will fall in domains needing Tier 2 models.

**Two design consequences from the newly added structural papers:**

- **Model domain *pairs*, not isolated domains, wherever the Ca²⁺ site is interdomain.** Downing
  1996 shows the cbEGF32–33 pair is a rigid rod stabilized by interdomain calcium, with a
  **calcium ligand donated from one domain to the other**. A variant's coordination shell can
  therefore lie partly in the neighbouring domain, and an isolated-domain model would
  systematically misrepresent it. This sharpens CLAUDE.md's "domain (or domain pair)" into a rule.
- **A useful stratum for the analysis:** Jensen 2012 defines the **"neonatal region" as TB3 →
  cbEGF18**, where severe neonatal Marfan substitutions cluster. Worth carrying as a Phase 3 flag.

---

## 6. Open questions for you

1. ~~Re-download the two wrong PDFs~~ — **done, verified (§4.1).** Phase 3 will follow the actual
   VCEP rules. Only the primary Jensen 2009 (2W86) citation remains outstanding, and the 2012
   review covers what Phase 4 needs.
2. **Phase 4 Tier-3 cross-check, given no AlphaFold DB / AlphaFill entry.** My recommendation:
   drop the AlphaFill route and use the **manual superposition transplant** from an experimental
   Ca²⁺ site (2W86 / 1LMJ / 1EMN) as the cross-check, recording donor accession and superposition
   RMSD as CLAUDE.md already specifies. Confirm, or tell me to investigate AlphaFill uploads.
3. **Phase 2 thresholds** (not needed to approve Phase 1, but flagging early since CLAUDE.md
   requires your sign-off before the labeled sets are finalised): ClinVar review-status star
   cutoff, and the gnomAD v4 Grpmax filtering-AF cutoff for the benign set. I will write these to
   `resources/reference/inclusion_criteria.md` once you decide.
4. **Modelling domain pairs (§5).** I plan to treat interdomain Ca²⁺ sites as requiring a
   **domain-pair** model rather than an isolated domain. This increases the number of AF3 jobs
   you will be asked to run in Phase 4.5. Flagging now since it affects your workload.

Note also that `resources/reference/structures.md` and `inclusion_criteria.md` referenced in
CLAUDE.md §2 do not exist yet — §5 above is the raw material for the former.

---

## 7. If approved, Phase 1 will

Pull FBN1 variants from ClinVar via E-utilities (capturing HGVS c./p., clinical significance,
**review status/stars**, condition, variation ID and dates — keeping VUS and conflicts labelled,
never dropped), load the AlphaMissense P35555 subset, and pull gnomAD v4 frequencies; land
everything immutably in `data/raw/` with manifest entries (source, query, version, date, row
count, checksum); then run the Phase 1 gate `tests/test_variant_ingest.py`.

**No Phase 1 work has started.**
