# FBN1 / P35555 experimental structure inventory

Enumerated live from the PDBe SIFTS API on 2026-08-02 (not a static list), then verified
against the deposited coordinates. Every number below was measured from the downloaded mmCIF
files in `structures/pdb/`, not taken from a paper or assumed.

## Inventory

| PDB | Method | Res. | UniProt range | Domain assignment | Models | Ca²⁺ per model | Author numbering |
|---|---|---|---|---|---|---|---|
| 2M74 | NMR | — | 45–178 | not yet assigned | 20 | **0** | — |
| 5MS9 | NMR | — | 113–287 | not yet assigned | 20 | 1 | — |
| **2W86** | X-ray | 1.80 Å | 807–951 | cbEGF9–hyb2–cbEGF10 | 1 | **2** | **local 1–147** (offset +804) |
| **1LMJ** | NMR | — | 1069–1154 | cbEGF12–13 | 25 | **2** | **local 3–88** (offset +1066) |
| **1UZJ** | X-ray | 2.25 Å | 1486–1647 | cbEGF22–TB4–cbEGF23 | 1 (3 chains A/B/C) | **2 per chain** | **= UniProt** |
| 1UZK | X-ray | 1.35 Å | 1486–1647 | cbEGF22–TB4–cbEGF23 | 1 | 1 (partial) | = UniProt |
| 1UZP | X-ray | 1.78 Å | 1486–1647 | cbEGF22–TB4–cbEGF23 | 1 | **0 (apo)** | = UniProt |
| 1UZQ | X-ray | 2.40 Å | 1486–1647 | cbEGF22–TB4–cbEGF23 | 1 | **0 (apo)** | = UniProt |
| 1APJ | NMR | — | 2052–2125 | isolated TB6 | 21 | **0** | — |
| **1EMN** | NMR | — | 2124–2205 | cbEGF32–33 | 1 | **2** | **= UniProt** |
| 1EMO | NMR | — | 2124–2205 | cbEGF32–33 | 22 | 2 | — |

Domain assignments for 2W86 / 1LMJ / 1UZJ / 1APJ / 1EMN are stated explicitly in Jensen,
Robertson & Handford (2012), *Structure* 20:215–225, and independently agree with the UniProt
domain boundaries and the cbEGF1..43 sequence-order index.

## ⚠️ Numbering — the trap this project must not fall into

**Author numbering differs between entries covering the same protein.** Verified residue-by-residue
against the P35555 sequence:

| PDB | Author → UniProt | Evidence |
|---|---|---|
| 2W86 | **+804** | 145/147 observed residues match; the 2 exceptions are non-native N-terminal cloning residues (auth 1–2) |
| 1LMJ | **+1066** | 86/86 observed residues match |
| 1UZJ | **+0** | 162/162 observed residues match |
| 1EMN | **+0** | 79/79 observed residues match |

**Consequences.** FoldX mutation strings and any per-residue structural metric must be written in
the *entry's own* author numbering. A P35555 position used unmodified against 2W86 or 1LMJ is
wrong by ~800 / ~1066 residues; used against 1UZJ or 1EMN it happens to be right — which makes
this failure mode especially dangerous, because it works for half the entries. **Always map via
SIFTS and re-check residue identity, never index.**

## Ca²⁺ stoichiometry (counted, not assumed)

**Exactly one Ca²⁺ per cbEGF domain; zero for TB and hybrid domains.** Confirmed by:
2W86 (2 cbEGF + 1 hybrid → 2 Ca), 1UZJ (2 cbEGF + 1 TB → 2 Ca per chain), 1EMN and 1LMJ
(2 cbEGF → 2 Ca), 1APJ (TB6 alone → 0 Ca).

Cysteine counts corroborate the architecture: 6 Cys per cbEGF domain, 8 per TB/hybrid domain
(cbEGF12–13 = 12 Cys; cbEGF32–33 = 12; cbEGF9–hyb2–cbEGF10 = 20; cbEGF22–TB4–cbEGF23 = 20).

## Structure selection guidance for Phase 4

- **cbEGF22–TB4–cbEGF23: use 1UZJ, not 1UZK.** 1UZK has the best resolution (1.35 Å) but only one
  of the two Ca²⁺ sites is occupied; 1UZP and 1UZQ are apo. Resolution is the wrong criterion when
  the metal is the measurement.
- **cbEGF32–33: 1EMN** is the single-model representative; **1EMO** is the 22-model ensemble and is
  the better choice if ensemble spread is wanted.
- **2M74 and 5MS9** were not named in CLAUDE.md Appendix B and their domain assignment has not yet
  been made. 2M74 contains no Ca²⁺.
- Coverage: these entries span roughly 8 of the 43 cbEGF domains, so most variants will need
  Tier-2 (AF3-cofolded) models.

## UniProt annotation caveat for Phase 3

UniProt P35555 annotates **155 disulfide bonds** but has **no `Binding site` features at all** —
there are **zero explicit calcium-site annotations**. CLAUDE.md Phase 3 asks us to validate the
computed coordinating-residue set against UniProt's explicit Ca-binding annotations; that
cross-check is **not available** and must be replaced by the consensus-motif derivation validated
against the disulfide annotations, the experimental Ca sites above, and the literature consensus
(Handford 1995 / Rao 1995 / Baudhuin 2019).

UniProt reports **47 EGF-like domains (43 calcium-binding)** and **9 TB domains**. Note the count
discrepancy with the literature's "7 TB": UniProt's TB set includes the 2 hybrid domains — e.g.
UniProt's **"TB 4" (851–902) is the domain the structural literature calls hybrid-2** (2W86).
Do not double-count these when tabulating domain classes.
