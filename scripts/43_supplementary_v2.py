#!/usr/bin/env python3
"""
43_supplementary_v2.py — PHASE 6. The v2 supplementary workbook.

Replaces `19_supplementary_tables.py`, which reads the v1 tables: sequence-motif calcium sites,
751 structurally covered variants, an enrichment denominator of 34 benign variants. Those are
superseded (see `results/SESSION_CONTEXT.md`), so a workbook built on them would contradict the
manuscript sheet by sheet.

Sheets:
  README             what each sheet is, where it came from, and the caveats to read first
  S1_variants        all 4,451 variants: labels, annotation, AF3 structural metrics, both ddGs
  S2_structural      the 3,831 with AF3 structural coverage, with their measured geometry
  S3_ca_sites        every calcium site in the 21 cofolded constructs, with its ligands
  S4_mutant_geometry the 60 batch-4 mutants x 3 sites, wild-type vs mutant coordination
  S5_disruption_rate one row per batch-4 variant at its cognate site, with the disruption call
  S6_statistics      every group comparison, effect size, CI and BH-corrected q
  S7_enrichment      pathogenic enrichment against the population set, circularity noted
  S8_validation      AF3 vs experimental calcium placement, and local BVS vs CheckMyMetal
  S9_dictionary      definition of every column used above
  S10_site_qc        every calcium site in those constructs against the three admission criteria

Every value is read from `data/processed/` or `results/`; nothing is recomputed here. A
supplementary table that recalculates its own numbers is a second chance to disagree with the
paper.

S10 was added on 2026-09-06. The geometry result rests on an admission filter applied to the
wild-type models, and until now the only place its outcome appeared was a run log, so a reviewer
could see that one construct was dropped but not how close anything else came. The sheet is built
by `scripts/64_site_qc.py`, which re-runs the filter with the functions script 40 used.

Writes: results/supplementary_tables_v3.xlsx, logs/43_supplementary_v2_<stamp>.log
"""
from __future__ import annotations

import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

OUT = fs.RESULTS / "supplementary_tables_v3.xlsx"

S1_COLUMNS = [
    "variation_id", "hgvs_c", "hgvs_p", "position", "wt_aa", "mut_aa",
    "clinvar_class", "review_status", "stars", "set_primary", "set_sensitivity",
    "gnomad_joint_af", "faf95_grpmax", "am_pathogenicity", "am_class",
    "domain", "domain_kind", "cbegf_index", "in_cbegf", "in_neonatal_region",
    "construct", "construct_range", "cbegf_domains", "plddt", "edge_margin",
    "ca_ligand_sidechain", "ca_ligand_backbone_only", "ca_ligand_any",
    "consensus_but_not_ligand", "is_cysteine", "cys_removing", "cys_introducing",
    "disulfide_partner", "mechanism_class",
    "rsa_pct", "buried", "dist_to_nearest_ca_A", "dist_ca_seed_sd",
    "rosetta_ddG", "foldx_ddG", "is_critical_residue", "vcep_pm1_category",
]

S2_COLUMNS = [
    "variation_id", "position", "wt_aa", "mut_aa", "set_primary", "mechanism_class",
    "construct", "cbegf_domains", "plddt", "edge_margin", "sasa_abs_A2", "rsa_pct", "buried",
    "dist_to_nearest_ca_A", "dist_ca_seed_sd", "ca_ligand_sidechain",
    "ca_ligand_backbone_only", "consensus_but_not_ligand", "cys_removing",
    "disulfide_partner", "rosetta_ddG", "foldx_ddG", "am_pathogenicity",
]

def _v1_definitions() -> list[tuple]:
    """Reuse the definitions written for the v1 workbook rather than paraphrasing them.

    The annotation columns (domain, cysteine role, VCEP category, labelled sets) are unchanged
    between v1 and v2, and a second wording of the same definition is a second thing to keep in
    step. Only the structural and energetic columns are redefined below, because those are what
    v2 changed.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_supp_v1", Path(__file__).resolve().parent / "19_supplementary_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    superseded = {"structural_coverage", "template_pdb", "residue", "sasa_abs_A2", "rsa_pct",
                  "buried", "dist_to_nearest_ca_A", "is_direct_ca_ligand", "mechanism_class"}
    return [d for d in mod.EXTRA_DEFINITIONS if d[1] not in superseded]


# Columns this workbook introduces that neither the Phase 2 dictionary nor the v1 workbook can
# document, because they did not exist when those were written.
EXTRA_DEFINITIONS = [
    ("Batch-4 identifiers", "job", "str", "AlphaFold Server job name for this mutant fold"),
    ("Batch-4 identifiers", "variant", "str", "the substitution in P35555 numbering"),
    ("Batch-4 identifiers", "klass", "str", "ca_ligand_pathogenic | composition_control | negative_control, fixed when the batch was designed"),
    ("Batch-4 identifiers", "clinical", "str", "ClinVar class of the variant at the time the batch was built"),
    ("Batch-4 identifiers", "uniprot_pos", "int", "position in P35555, re-derived from the model sequence rather than trusted from the request"),
    ("Batch-4 identifiers", "site", "int", "index of the calcium site within the construct, 0-2"),
    ("Batch-4 identifiers", "site_pos", "str", "which consensus position the variant hits: 'Asp/Glu' (offsets 0 and +3 from the domain start) or 'beta-OH Asn' (offsets +16 to +20). The donor-chemistry rule behaves differently in the two, and the ClinGen FBN1 VCEP treats Asn->Ser in the second differently again"),
    ("Batch-4 identifiers", "cognate_by", "str", "how the cognate site was identified: 'coordinates' (the residue is a wild-type ligand) or 'nearest ion'"),
    ("Batch-4 identifiers", "status", "str", "present | vacated — whether an ion still occupies the wild-type site in the mutant models"),
    ("v2 structure", "construct", "str", "AF3 cofolded tandem construct this variant's position falls in"),
    ("v2 structure", "construct_range", "str", "UniProt residue range of that construct"),
    ("v2 structure", "cbegf_domains", "str", "cbEGF domains the construct spans"),
    ("v2 structure", "plddt", "float", "AF3 per-residue confidence at this position, median over five models"),
    ("v2 structure", "edge_margin", "int", "residues between this position and the nearest construct terminus; small values are less reliable"),
    ("v2 structure", "sasa_abs_A2", "float", "absolute solvent-accessible surface area (freesasa) in the AF3 model"),
    ("v2 structure", "rsa_pct", "float", "relative solvent accessibility, % of the Tien et al. maximum for that residue"),
    ("v2 structure", "buried", "bool", "RSA below 20%"),
    ("v2 structure", "dist_to_nearest_ca_A", "float", "distance from any side-chain atom to the nearest modelled Ca2+, median over five models"),
    ("v2 structure", "dist_ca_seed_sd", "float", "seed-to-seed spread of that distance"),
    ("v2 site class", "is_cysteine", "bool", "wild-type residue is cysteine (the v2 name for the v1 is_cysteine_site)"),
    ("v2 site class", "cys_introducing", "bool", "the substitution introduces a cysteine where there was none, leaving an unpaired thiol"),
    ("v2 site class", "ca_ligand_sidechain", "bool", "a side-chain atom of this residue coordinates Ca2+ (within 3.2 A) in the AF3 models"),
    ("v2 site class", "ca_ligand_backbone_only", "bool", "only this residue's backbone carbonyl coordinates Ca2+; no substitution can remove it"),
    ("v2 site class", "ca_ligand_any", "bool", "either of the above"),
    ("v2 site class", "consensus_but_not_ligand", "bool", "the cbEGF sequence motif flags this position but it coordinates no Ca2+ in any model"),
    ("v2 site class", "mechanism_class", "str", "ca_ligand / composition_control / motif_overcalled / cys_removing / not_in_cbEGF"),
    ("v2 energetics", "rosetta_ddG", "float", "PyRosetta ref2015 ddG in Rosetta energy units, metal bonds formed explicitly, torsion-space minimisation with restrained backbone relaxation in the repack shell"),
    ("v2 energetics", "foldx_ddG", "float", "FoldX 5 ddG in kcal/mol; available only for the 751 variants on experimental templates"),
    ("v2 geometry", "bvs_wt", "float", "bond-valence sum at the calcium site in the wild-type construct, median over five AF3 models"),
    ("v2 geometry", "bvs_mut", "float", "the same quantity in the mutant construct"),
    ("v2 geometry", "delta_bvs_pct", "float", "percentage change in bond-valence sum, mutant relative to wild type"),
    ("v2 geometry", "z_vs_noise", "float", "that change in units of the wild-type seed spread, denominator floored at the median site spread"),
    ("v2 geometry", "wt_seed_sd_pct", "float", "wild-type seed spread as a percentage of the site's own valence"),
    ("v2 geometry", "cn_wt", "float", "coordination number of the wild-type site"),
    ("v2 geometry", "cn_mut", "float", "coordination number of the mutant site"),
    ("v2 geometry", "ion_shift_A", "float", "displacement of the ion after superposing the site's own neighbourhood"),
    ("v2 geometry", "ion_shift_wt_sd_A", "float", "the same displacement between wild-type seeds; the noise floor for the column above"),
    ("v2 geometry", "local_rmsd_A", "float", "backbone RMSD over residues within 12 A of the ion, mutant vs wild type"),
    ("v2 geometry", "rmsd_ca_global_A", "float", "whole-construct backbone RMSD; large values usually mean inter-domain hinging, not a broken fold"),
    ("v2 geometry", "plddt_site_wt", "float", "AF3 confidence over the coordinating residues, wild type"),
    ("v2 geometry", "plddt_site_mut", "float", "the same, mutant"),
    ("v2 geometry", "still_contacts", "int", "of five models, how many retain contact between the substituted side chain and the ion; blank where the residue was never a ligand"),
    ("v2 geometry", "donor_kept", "bool", "the substituted residue can still donate an oxygen (D/E/N/Q/S/T)"),
    ("v2 geometry", "disrupted", "bool", "valence loss exceeds the measured null threshold, or the site is vacated"),
    ("v2 geometry", "gain", "bool", "valence rose beyond the same threshold"),
    ("v2 geometry", "null_threshold_pct", "float", "95th percentile of absolute valence change over 123 sites where no effect is possible"),
    ("v2 geometry", "is_cognate", "bool", "this is the site the substitution can reach; the other two sites of the construct are within-construct controls"),
    ("v2 geometry", "wt_site_ok", "bool", "the wild-type site passed the quality filter (seed spread <=15% of valence, found in 5/5 models, local pLDDT >=70)"),
    ("v2 site QC", "job_name", "str", "the AlphaFold Server job that produced this construct's five wild-type models"),
    ("v2 site QC", "uniprot_range", "str", "the P35555 residue range the construct spans, verified against the reference sequence"),
    ("v2 site QC", "n_models_present", "int", "of five wild-type models, how many place an ion at this site; the second admission criterion, limit 5"),
    ("v2 site QC", "ligand_positions_uniprot", "str", "P35555 positions of the residues coordinating this site in the first wild-type model"),
    ("v2 site QC", "plddt_site", "float", "mean AlphaFold confidence over the residues coordinating this site, median over five wild-type models"),
    ("v2 site QC", "plddt_construct", "float", "mean AlphaFold confidence over the whole construct, median over five wild-type models"),
    ("v2 site QC", "seed_sd_bvs", "float", "seed-to-seed spread of the site's bond-valence sum, as a robust standard deviation"),
    ("v2 site QC", "seed_sd_pct", "float", "that spread as a percentage of the site's own valence; the first admission criterion, limit 15%"),
    ("v2 site QC", "ion_seed_shift_A", "float", "how far the ion moves between wild-type seeds after superposing the site's own neighbourhood"),
    ("v2 site QC", "pass_seed_spread", "bool", "seed spread is at most 15% of the site's valence"),
    ("v2 site QC", "pass_occupancy", "bool", "the site is found in all five wild-type models"),
    ("v2 site QC", "pass_plddt", "bool", "local pLDDT over the coordinating residues is at least 70"),
    ("v2 site QC", "admitted", "bool", "all three criteria are met, so the site is measurable in a mutant"),
    ("v2 site QC", "margin_to_nearest_limit", "float", "distance to the nearest of the three limits as a fraction of that limit; small positive values are the borderline sites"),
    ("v2 site QC", "n_variants_at_site", "int", "how many of the 60 cofolded variants have this site as their cognate site"),
]

README_ROWS = [
    ("Study", "Pathogenic FBN1 missense variants at calcium-binding positions impair ion coordination without destabilising the cbEGF fold"),
    ("Protein / transcript", "UniProt P35555 (fibrillin-1, 2,871 aa); MANE Select NM_000138.5"),
    ("Generated", datetime.now(timezone.utc).isoformat()),
    ("Built by", "scripts/43_supplementary_v2.py, from data/processed/ and results/"),
    ("VERSION", "This is the v3 workbook. It is the v2 workbook plus S10_site_qc, and both supersede supplementary_tables.xlsx, which was built on the sequence-motif definition of a calcium site and on 751 structurally covered variants."),
    ("Calcium sites", "Defined by observed coordination in AlphaFold3 models cofolded with Ca2+, not by the sequence motif. The motif over-calls: 143 positions it flags coordinate no calcium in any model, and UniProt annotates no calcium site in P35555 at all."),
    ("Structural coverage", "All 43 cbEGF domains modelled as 21 overlapping tandem constructs with three Ca2+ each, five models per construct. 3,831 of 4,451 variants (86%) have structural measurements."),
    ("Ca2+ is modelled", "AlphaFold2 models contain no metals. Every calcium position here was either taken from an experimental structure or placed by AF3 cofolding. Placement was calibrated against eight domains with experimental structures (median ion deviation 0.45 A). It remains a modelling assumption, not a measurement."),
    ("Two energy functions", "FoldX 5 (kcal/mol) and PyRosetta ref2015 (Rosetta energy units). These are different quantities and must never be pooled or plotted on a shared axis."),
    ("Enrichment caveat", "The enrichment sheet is partly circular: ClinGen's FBN1 rules assign pathogenicity partly on the basis of position at a critical residue, which is what the analysis tests. It is reported for completeness, not as evidence of causation."),
    ("Benign set", "34 benign/likely benign variants, 29 with structural coverage, none at a calcium ligand or a cysteine. Comparisons against it are underpowered and are reported with intervals rather than p-values."),
    ("Guideline rules referenced", "ClinGen FBN1 VCEP (Genome Med 16:154, 2024): PM1_Strong for cysteine-removing variants in cbEGF domains; PM1 (moderate) for the consensus calcium residues [D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F]; Asn>Ser at the second [D/N] exempted; PP3 keyed to REVEL >= 0.75."),
    ("Site admission", "A calcium site is measurable in a mutant only if its wild-type models agree with each other. S10_site_qc lists all 63 sites against the three criteria, which were fixed before any mutant was read: seed spread at most 15% of the site's own valence, occupancy in all five wild-type models, and local pLDDT at least 70. 60 sites are admitted. Two sites of cbEGF1-3 fail, one on seed spread and one on occupancy, and with them three variants; one site of cbEGF5-7 fails on seed spread and carries no variant. One admitted site, cbEGF19-21 site 2, lies within 20% of a limit."),
    ("Contact", "See results/manuscript_v4.md for the accompanying text and figure list."),
]


def parse_dictionary(md: Path) -> pd.DataFrame:
    rows, section = [], ""
    for line in md.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            section = line[3:].strip()
        m = re.match(r"^\|\s*`([^`]+)`\s*\|\s*([^|]*?)\s*\|\s*(.*?)\s*\|\s*$", line)
        if m:
            rows.append({"section": section, "column": m.group(1),
                         "type": m.group(2), "meaning": m.group(3)})
    return pd.DataFrame(rows)


def autoformat(writer, sheet: str, df: pd.DataFrame, wrap_cols: set[str] = frozenset()) -> None:
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter

    ws = writer.sheets[sheet]
    head_fill = PatternFill("solid", fgColor="EFEEE9")
    for j, col in enumerate(df.columns, start=1):
        cell = ws.cell(row=1, column=j)
        cell.font = Font(bold=True)
        cell.fill = head_fill
        cell.alignment = Alignment(vertical="top", wrap_text=True)
        body = df[col].astype(str)
        width = max(len(str(col)), int(body.str.len().quantile(0.95)) if len(body) else 0) + 2
        if col in wrap_cols:
            width = 70
            for i in range(2, len(df) + 2):
                ws.cell(row=i, column=j).alignment = Alignment(vertical="top", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = min(max(width, 9), 70)
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions


def main() -> int:
    log = [f"Phase 6 — v2 supplementary workbook — {datetime.now(timezone.utc).isoformat()}"]

    df = fs.load_v2()
    df["mechanism_class"] = fs.mech_group_v2(df).fillna("not_in_cbEGF")

    missing = [c for c in S1_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"S1 columns absent from the processed tables: {missing}")

    s1 = df[S1_COLUMNS].sort_values("position").reset_index(drop=True)
    s2 = (df[df.has_structure][S2_COLUMNS]
          .sort_values(["construct", "position"]).reset_index(drop=True))
    s3 = pd.read_csv(fs.PROC / "af3_ca_sites.tsv", sep="\t")
    lig = pd.read_csv(fs.PROC / "af3_ca_ligands.tsv", sep="\t")
    s4 = pd.read_csv(fs.PROC / "af3_batch4_geometry.tsv", sep="\t")
    s5 = pd.read_csv(fs.PROC / "batch4_rate.tsv", sep="\t")
    # statistics_final supersedes phase5_statistics_v2: the latter predates both the
    # corrected PyRosetta ddG and batch 4, and has no Rosetta family at all
    s6 = pd.read_csv(fs.RESULTS / "statistics_final.tsv", sep="\t")
    s7 = pd.read_csv(fs.PROC / "enrichment_v2.tsv", sep="\t")
    cal = pd.read_csv(fs.RESULTS / "phase4_calibration.tsv", sep="\t")
    cmm = pd.read_csv(fs.PROC / "bvs_vs_cmm.tsv", sep="\t")
    s10 = pd.read_csv(fs.PROC / "af3_site_qc.tsv", sep="\t")
    s8 = pd.concat([cal.assign(comparison="AF3 vs experimental placement"),
                    cmm.assign(comparison="local BVS vs CheckMyMetal")], ignore_index=True)

    docs = pd.concat([parse_dictionary(fs.PROC / "variant_master_dictionary.md"),
                      pd.DataFrame(_v1_definitions() + EXTRA_DEFINITIONS,
                                   columns=["section", "column", "type", "meaning"])],
                     ignore_index=True)
    used = set(S1_COLUMNS) | set(S2_COLUMNS) | set(s4.columns) | set(s5.columns) | set(s10.columns)
    s9 = docs[docs.column.isin(used)].drop_duplicates("column").reset_index(drop=True)
    undocumented = sorted(c for c in set(S1_COLUMNS) | set(S2_COLUMNS) | set(s5.columns)
                          | set(s10.columns) if c not in set(s9.column))
    if undocumented:
        raise SystemExit(f"columns in the supplementary tables with no definition: "
                         f"{undocumented}")

    if len(s1) != 4451:
        raise SystemExit(f"S1 has {len(s1)} rows, expected 4,451")
    if len(s2) != int(df.has_structure.sum()):
        raise SystemExit("S2 row count does not match the structurally covered set")
    if len(s5) != s5.variation_id.nunique():
        raise SystemExit("S5 has duplicate variants")
    if len(s10) != 63 or int(s10.admitted.sum()) != 60:
        raise SystemExit(f"S10 has {len(s10)} sites and {int(s10.admitted.sum())} admitted, "
                         "expected 63 and 60")

    # ligand detail belongs beside the sites it explains
    s3 = s3.merge(lig.groupby(["construct", "site"]).apply(
        lambda g: "; ".join(f"{r.wt_aa}{r.uniprot_pos}:{r.atom}({r.atom_class[:2]})"
                            for r in g.itertuples()), include_groups=False)
        .rename("ligand_atoms").reset_index(), on=["construct", "site"], how="left")

    readme = pd.DataFrame(README_ROWS, columns=["Item", "Detail"])

    with pd.ExcelWriter(OUT, engine="openpyxl") as w:
        sheets = [
            ("README", readme, {"Detail"}),
            ("S1_variants", s1, set()),
            ("S2_structural", s2, set()),
            ("S3_ca_sites", s3, {"ligand_atoms", "ligand_positions_uniprot"}),
            ("S4_mutant_geometry", s4, set()),
            ("S5_disruption_rate", s5, set()),
            ("S6_statistics", s6, {"note", "test"}),
            ("S7_enrichment", s7, set()),
            ("S8_validation", s8, set()),
            ("S9_dictionary", s9, {"meaning"}),
            ("S10_site_qc", s10, {"cbegf_domains", "ligand_positions_uniprot"}),
        ]
        for name, frame, wrap in sheets:
            frame.to_excel(w, sheet_name=name, index=False)
            autoformat(w, name, frame, wrap)
            log.append(f"  {name:20} {len(frame):5} rows x {len(frame.columns):2} cols")
            print(f"  {name:20} {len(frame):5} rows × {len(frame.columns):2} cols")

    log.append(f"wrote {OUT.relative_to(fs.ROOT)} ({OUT.stat().st_size/1024:.0f} kB)")
    print(f"  wrote {OUT.relative_to(fs.ROOT)} ({OUT.stat().st_size/1024:.0f} kB)")
    fs.write_log("43_supplementary_v2", log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
