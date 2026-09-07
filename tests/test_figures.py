#!/usr/bin/env python3
"""
tests/test_figures.py — PHASE 6 GATE.

The Phase 6 requirement is that every figure has a regenerating script and a caption stub, and
that the numbers in the figures match `data/processed/`. A PNG is opaque, so each figure script
records every value it annotated into `figures/<name>.numbers.json`; this gate recomputes those
values from the processed tables and the Phase 4/5 result files and asserts they agree.

Asserts:
  assert_figures_exist       every figure has a PNG and a vector PDF, both non-trivial
  assert_scripts_exist       each figure names the script that made it, and that script is present
  assert_captions_present    every figure has a caption naming its file and its script
  assert_numbers_match_data  every recorded value re-derives from data/processed or results/
  assert_pymol_panels        the structural figure's .pml scripts and panels exist; PyMOL runs
  assert_supplementary       the workbook has every sheet, right row counts, no undocumented column
  assert_shared_palette      no figure script hard-codes a colour outside figstyle.PALETTE
  assert_seeds_fixed         every figure script draws its seed from the shared module
  assert_reproducible        re-running each figure script reproduces identical recorded numbers

Run:  python tests/test_figures.py
Exit: 0 all pass, 1 otherwise.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIGURES = ROOT / "figures"
RESULTS = ROOT / "results"
PROC = ROOT / "data/processed"
sys.path.insert(0, str(ROOT / "scripts"))

import figstyle as fs  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

FIGS = {
    "fig1_variant_landscape": "13_fig1_landscape.py",
    "fig2_mechanism_dissociation": "14_fig2_mechanism.py",
    "fig3_burial_control": "15_fig3_burial.py",
    "fig4_structure": "16_fig4_structure.py",
    "fig5_ca_placement": "17_fig5_calibration.py",
    "fig6_enrichment": "18_fig6_enrichment.py",
}
TOL = 5e-4

_failures: list[str] = []


class GateFailure(AssertionError):
    pass


def numbers(fig: str) -> dict:
    p = FIGURES / f"{fig}.numbers.json"
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    return json.loads(p.read_text(encoding="utf-8"))


def close(a, b) -> bool:
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        return (len(a) == len(b)) and all(close(x, y) for x, y in zip(a, b))
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) or isinstance(b, str):
        return a == b
    return abs(float(a) - float(b)) <= TOL * max(1.0, abs(float(b)))


def check(found: dict, key: str, expected, mismatches: list) -> None:
    if key not in found:
        mismatches.append(f"{key}: not recorded by the figure")
        return
    if not close(found[key], expected):
        mismatches.append(f"{key}: figure says {found[key]!r}, data says {expected!r}")


def merged() -> pd.DataFrame:
    df = fs.load_merged()
    df["mech"] = fs.mech_group(df)
    return df


def stats_table() -> pd.DataFrame:
    return pd.read_csv(RESULTS / "phase5_statistics.tsv", sep="\t")


def effect(tbl: pd.DataFrame, test: str) -> float:
    hit = tbl[tbl.test == test]
    if len(hit) != 1:
        raise GateFailure(f"phase5_statistics.tsv has {len(hit)} rows for {test!r}")
    return float(hit.iloc[0].effect)


# --------------------------------------------------------------------------
def assert_figures_exist() -> str:
    for fig in FIGS:
        for ext, floor in (("png", 40_000), ("pdf", 5_000)):
            p = FIGURES / f"{fig}.{ext}"
            if not p.is_file():
                raise GateFailure(f"missing {p.relative_to(ROOT)}")
            if p.stat().st_size < floor:
                raise GateFailure(f"{p.relative_to(ROOT)} is only {p.stat().st_size} bytes — "
                                  f"probably an empty or failed render")
    return f"{len(FIGS)} figures, each with a 300 dpi PNG and a vector PDF"


def assert_scripts_exist() -> str:
    for fig, script in FIGS.items():
        p = ROOT / "scripts" / script
        if not p.is_file():
            raise GateFailure(f"{fig} has no regenerating script at scripts/{script}")
        declared = numbers(fig)["script"]
        if declared != script:
            raise GateFailure(f"{fig} was written by {declared}, but the gate expects {script}")
    return f"all {len(FIGS)} figures regenerate from a committed script"


def assert_captions_present() -> str:
    cap = RESULTS / "figure_captions.md"
    if not cap.is_file():
        raise GateFailure("results/figure_captions.md is missing")
    text = cap.read_text(encoding="utf-8")
    for fig, script in FIGS.items():
        if f"{fig}.png" not in text:
            raise GateFailure(f"figure_captions.md never mentions {fig}.png")
        if script not in text:
            raise GateFailure(f"figure_captions.md does not name scripts/{script}")
    headings = re.findall(r"^## Figure (\d+)", text, flags=re.M)
    if sorted(headings) != [str(i) for i in range(1, len(FIGS) + 1)]:
        raise GateFailure(f"expected one '## Figure N' caption per figure, found {headings}")
    return f"{len(headings)} captions, each naming its output file and its script"


def assert_numbers_match_data() -> str:
    df = merged()
    tbl = stats_table()
    cov = df[df.ddG.notna()]
    mism: list[str] = []
    checked = 0

    # ---- Figure 1 -------------------------------------------------------------------
    v = numbers("fig1_variant_landscape")["values"]
    check(v, "n_variants_total", len(df), mism)
    for cls in ("pathogenic", "vus", "benign"):
        sub = df[df.set_primary == cls]
        check(v, f"n_{cls}", len(sub), mism)
        check(v, f"n_positions_{cls}", int(sub.position.nunique()), mism)
        check(v, f"max_stack_{cls}", int(sub.groupby("position").size().max()), mism)
    dom = fs.domains()
    for kind in ("cbEGF", "EGF", "TB"):
        check(v, f"n_domains_{kind}", int((dom.kind == kind).sum()), mism)
    path = df[df.set_primary == "pathogenic"]
    pcys = path.in_cbegf_cys & path.cys_removing
    check(v, "panelC_pathogenic_total", len(path), mism)
    check(v, "panelC_pathogenic_in_cbegf", int(path.in_cbegf.sum()), mism)
    check(v, "panelC_total_cys", int((pcys & path.cbegf_index.notna()).sum()), mism)
    check(v, "panelC_total_ca", int((~pcys & path.is_ca_consensus
                                    & path.cbegf_index.notna()).sum()), mism)
    check(v, "panelC_total_other", int((~pcys & ~path.is_ca_consensus
                                        & path.cbegf_index.notna()).sum()), mism)
    checked += len(v)

    # ---- Figure 2 -------------------------------------------------------------------
    v = numbers("fig2_mechanism_dissociation")["values"]
    check(v, "n_with_ddg", len(cov), mism)
    for key in ("ca_ligand", "cys_removing", "other_cbegf"):
        g = cov[cov.mech == key]
        check(v, f"ddG_n_{key}", len(g), mism)
        check(v, f"ddG_median_{key}", float(np.median(g.ddG)), mism)
        check(v, f"am_median_{key}", float(np.median(g.am.dropna())), mism)
        check(v, f"disulfide_median_{key}",
              float(np.median(g.ddG_disulfide.dropna())), mism)
    lig = cov[(cov.mech == "ca_ligand") & cov.am.notna()]
    check(v, "ca_ligand_with_am", len(lig), mism)
    check(v, "ca_ligand_damaging_but_stable",
          int(((lig.am > 0.564) & (lig.ddG < 2.0)).sum()), mism)
    check(v, "stat_ddG_cys_vs_rest_delta",
          effect(tbl, "ddG: cbEGF cysteine-removing vs all other missense"), mism)
    check(v, "stat_ddG_ligand_vs_nonligand_delta",
          effect(tbl, "ddG: direct Ca ligands vs non-ligands (structure-measured)"), mism)
    check(v, "stat_am_ligand_vs_other_delta",
          effect(tbl, "AlphaMissense: direct Ca ligands vs other cbEGF residues"), mism)
    check(v, "stat_disulfide_delta",
          effect(tbl, "FoldX disulfide energy term: Cys-removing vs other"), mism)
    check(v, "stat_spearman_rho",
          effect(tbl, "Spearman: FoldX ddG vs AlphaMissense"), mism)
    checked += len(v)

    # ---- Figure 3 -------------------------------------------------------------------
    v = numbers("fig3_burial_control")["values"]
    for key in ("ca_ligand", "cys_removing", "other_cbegf"):
        g = cov[cov.mech == key]
        check(v, f"rsa_n_{key}", int(g.rsa.notna().sum()), mism)
        check(v, f"rsa_median_{key}", float(np.median(g.rsa.dropna())), mism)
        check(v, f"dist_ca_median_{key}",
              float(np.median(g.dist_to_nearest_ca_A.dropna())), mism)
    ca = cov[cov.is_ca_consensus]
    other = cov[cov.in_cbegf & ~cov.is_ca_consensus]
    for lo, hi, lbl in ((0.0, 20.0, "buried"), (20.0, 50.0, "partial"),
                        (50.0, np.inf, "exposed")):
        for name, sub in (("Ca-consensus", ca), ("other_cbEGF", other)):
            vals = sub[(sub.rsa >= lo) & (sub.rsa < hi)].ddG.values
            check(v, f"panelC_{lbl}_{name}_n", int(len(vals)), mism)
            if len(vals):
                check(v, f"panelC_{lbl}_{name}_median", float(np.median(vals)), mism)
    check(v, "stat_buried_delta",
          effect(tbl, "ddG: Ca-consensus vs other cbEGF, buried (RSA<20%)"), mism)
    check(v, "stat_partial_delta",
          effect(tbl, "ddG: Ca-consensus vs other cbEGF, partial (RSA 20-50%)"), mism)
    checked += len(v)

    # ---- Figure 4 -------------------------------------------------------------------
    v = numbers("fig4_structure")["values"]
    sub = df[df.template_pdb == "2W86"]
    p2 = sub[sub.set_primary == "pathogenic"]
    check(v, "n_variants_on_template", len(sub), mism)
    check(v, "n_pathogenic_on_template", len(p2), mism)
    check(v, "n_pathogenic_positions_on_template", int(p2.position.nunique()), mism)
    for key in ("ca_ligand", "cys_removing", "other_cbegf"):
        check(v, f"pathogenic_positions_{key}",
              sorted(int(x) for x in p2[p2.mech == key].position.unique()), mism)
    check(v, "pathogenic_positions_outside_cbegf",
          sorted(int(x) for x in p2[p2.mech.isna()].position.unique()), mism)
    calib = pd.read_csv(RESULTS / "phase4_calibration.tsv", sep="\t")
    c9 = calib[(calib.domain == "cbEGF9") & (calib.model == 0)].iloc[0]
    check(v, "af3_ca_deviation_cbEGF9_A", float(c9.ca_deviation_A), mism)
    check(v, "af3_backbone_rmsd_cbEGF9_A", float(c9.domain_backbone_rmsd_A), mism)
    checked += len(v)

    # ---- Figure 5 -------------------------------------------------------------------
    v = numbers("fig5_ca_placement")["values"]
    loo = pd.read_csv(RESULTS / "phase4_transplant_loo.tsv", sep="\t")
    xray = calib[calib.reference_pdb.isin(("2W86", "1UZJ"))]
    nmr = calib[calib.reference_pdb.isin(("1LMJ", "1EMN"))]
    check(v, "ca_dev_n_af3_xray", len(xray), mism)
    check(v, "ca_dev_n_af3_nmr", len(nmr), mism)
    check(v, "ca_dev_n_transplant", len(loo), mism)
    check(v, "ca_dev_median_af3_xray", float(np.median(xray.ca_deviation_A)), mism)
    check(v, "ca_dev_median_af3_nmr", float(np.median(nmr.ca_deviation_A)), mism)
    check(v, "ca_dev_median_transplant",
          float(np.median(loo.ca_transplant_error_A)), mism)
    check(v, "backbone_rmsd_median_per_domain",
          float(np.median(calib.domain_backbone_rmsd_A)), mism)
    check(v, "backbone_rmsd_median_whole_construct",
          float(np.median(calib.construct_backbone_rmsd_A)), mism)
    cmm = pd.read_csv(PROC / "cmm_sites.tsv", sep="\t")
    exp, af3 = cmm[cmm.source == "experimental"], cmm[cmm.source != "experimental"]
    check(v, "cmm_n_experimental", len(exp), mism)
    check(v, "cmm_n_af3", len(af3), mism)
    for col in ("grmsd", "valence", "nvecsum", "coordination_number"):
        check(v, f"cmm_{col}_median_experimental", float(np.median(exp[col])), mism)
        check(v, f"cmm_{col}_median_af3", float(np.median(af3[col])), mism)
    checked += len(v)

    # ---- Figure 6 -------------------------------------------------------------------
    v = numbers("fig6_enrichment")["values"]
    for key, test in (("cys", "pathogenic at cbEGF cysteine (removing) vs background"),
                      ("ca", "pathogenic at calcium-consensus residue vs background")):
        row = tbl[tbl.test == test].iloc[0]
        check(v, f"background_fold_{key}", float(row.effect), mism)
        check(v, f"background_ci_{key}", [float(row.ci_low), float(row.ci_high)], mism)
    for test in ("pathogenic vs benign at any VCEP critical residue",
                 "pathogenic vs benign at cbEGF cysteine (removing)",
                 "pathogenic vs benign at calcium-consensus residue",
                 "AM-pathogenic vs AM-benign at cbEGF cysteine (removing)",
                 "AM-pathogenic vs AM-benign at calcium-consensus residue"):
        row = tbl[tbl.test == test].iloc[0]
        src = "clinvar" if test.startswith("pathogenic") else "orthogonal"
        label = ("Any VCEP critical residue" if "VCEP" in test else
                 "cbEGF cysteine-removing" if "cysteine" in test else
                 "Calcium-consensus residue")
        check(v, f"or_{src}_{label[:14].strip().replace(' ', '_')}", float(row.effect), mism)
    checked += len(v)

    if mism:
        raise GateFailure(f"{len(mism)} figure value(s) disagree with the data: "
                          + "; ".join(mism[:6]) + (" ..." if len(mism) > 6 else ""))
    return f"{checked} recorded figure values all re-derive from data/processed and results/"


def assert_pymol_panels() -> str:
    pml_dir, panel_dir = FIGURES / "pml", FIGURES / "panels"
    expected = ["fig4_A_fragment", "fig4_B_casite", "fig4_C_af3_overlay", "fig4_D_variants"]
    for name in expected:
        for d, ext in ((pml_dir, "pml"), (panel_dir, "png")):
            p = d / f"{name}.{ext}"
            if not p.is_file():
                raise GateFailure(f"missing {p.relative_to(ROOT)}")
    packet = ROOT / "handoff/cmm/inbox/EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb"
    for name in expected:
        body = (pml_dir / f"{name}.pml").read_text(encoding="utf-8")
        if str(packet) not in body:
            raise GateFailure(f"{name}.pml does not load the UniProt-renumbered packet PDB — "
                              f"labels would be in PDB author numbering")
        if "remove exp and resi 805-806" not in body:
            raise GateFailure(f"{name}.pml does not strip 2W86's non-native cloning residues")
    out = subprocess.run(["pymol", "-cq", "-d", "print(cmd.get_version()[0])"],
                         capture_output=True, text=True, timeout=180)
    if out.returncode != 0:
        raise GateFailure("PyMOL is not callable, so Figure 4 cannot be regenerated")
    version = out.stdout.strip().splitlines()[-1] if out.stdout.strip() else "?"
    return f"{len(expected)} PyMOL panels with their .pml scripts; PyMOL {version} callable"


def assert_supplementary() -> str:
    xlsx = RESULTS / "supplementary_tables.xlsx"
    if not xlsx.is_file():
        raise GateFailure("results/supplementary_tables.xlsx is missing")
    sheets = pd.read_excel(xlsx, sheet_name=None)
    want = ["README", "S1_variants", "S2_structural", "S3_cmm_sites", "S4_statistics",
            "S5_calibration", "S6_transplant", "S7_dictionary"]
    missing = [s for s in want if s not in sheets]
    if missing:
        raise GateFailure(f"workbook is missing sheet(s): {missing}")

    df = merged()
    expect_rows = {
        "S1_variants": len(df),
        "S2_structural": int(df.structural_coverage.sum()),
        "S3_cmm_sites": len(pd.read_csv(PROC / "cmm_sites.tsv", sep="\t")),
        "S4_statistics": len(stats_table()),
        "S5_calibration": len(pd.read_csv(RESULTS / "phase4_calibration.tsv", sep="\t")),
        "S6_transplant": len(pd.read_csv(RESULTS / "phase4_transplant_loo.tsv", sep="\t")),
    }
    for sheet, n in expect_rows.items():
        if len(sheets[sheet]) != n:
            raise GateFailure(f"{sheet} has {len(sheets[sheet])} rows, source has {n}")

    documented = set(sheets["S7_dictionary"]["column"])
    undocumented = sorted((set(sheets["S1_variants"].columns)
                           | set(sheets["S2_structural"].columns)) - documented)
    if undocumented:
        raise GateFailure(f"supplementary columns with no dictionary entry: {undocumented}")

    # Spot-check that the workbook is a view of the data, not a re-derivation of it.
    s1 = sheets["S1_variants"]
    src = df.set_index("variation_id")
    sample = s1.sample(min(200, len(s1)), random_state=fs.SEED)
    for _, row in sample.iterrows():
        ref = src.loc[row.variation_id]
        if int(ref.position) != int(row.position) or ref.wt_aa != row.wt_aa:
            raise GateFailure(f"S1 row for variation {row.variation_id} disagrees with "
                              f"data/processed")
    return (f"{len(want)} sheets; row counts match source; "
            f"{len(documented)} documented columns, none missing")


def assert_shared_palette() -> str:
    """
    One palette across all figures. A stray hex literal in a figure script is how a set of
    panels quietly stops being comparable, so the gate refuses them outside figstyle.py.
    """
    allowed = {c.lower() for c in fs.PALETTE.values()}
    offenders = []
    for script in FIGS.values():
        text = (ROOT / "scripts" / script).read_text(encoding="utf-8")
        for m in re.finditer(r"#[0-9a-fA-F]{6}\b", text):
            if m.group(0).lower() not in allowed:
                line = text[:m.start()].count("\n") + 1
                offenders.append(f"{script}:{line} {m.group(0)}")
    if offenders:
        raise GateFailure("colours outside figstyle.PALETTE: " + ", ".join(offenders))
    return f"all {len(FIGS)} figure scripts draw only from figstyle.PALETTE"


def assert_seeds_fixed() -> str:
    for script in FIGS.values():
        text = (ROOT / "scripts" / script).read_text(encoding="utf-8")
        if "fs.SEED" not in text:
            raise GateFailure(f"scripts/{script} does not use the shared seed")
    if not re.search(r"^SEED = \d+", (ROOT / "scripts/figstyle.py").read_text(encoding="utf-8"),
                     flags=re.M):
        raise GateFailure("figstyle.py does not pin a seed")
    return f"seed {fs.SEED} pinned in figstyle.py and used by all {len(FIGS)} figure scripts"


def assert_reproducible() -> str:
    before = {fig: numbers(fig)["values"] for fig in FIGS}
    for fig, script in FIGS.items():
        out = subprocess.run([sys.executable, str(ROOT / "scripts" / script)],
                             capture_output=True, text=True, timeout=1800, cwd=ROOT)
        if out.returncode != 0:
            raise GateFailure(f"re-running scripts/{script} failed:\n{out.stderr[-800:]}")
    for fig in FIGS:
        after = numbers(fig)["values"]
        if before[fig].keys() != after.keys():
            raise GateFailure(f"{fig} recorded a different set of values on re-run")
        drift = [k for k in after if not close(after[k], before[fig][k])]
        if drift:
            raise GateFailure(f"{fig} changed on re-run: {drift[:5]}")
    return f"all {len(FIGS)} figure scripts re-ran and reproduced identical recorded numbers"


def main() -> int:
    print("PHASE 6 GATE — tests/test_figures.py")
    print(f"root: {ROOT}\n")
    checks = [
        ("figures_exist", assert_figures_exist),
        ("scripts_exist", assert_scripts_exist),
        ("captions_present", assert_captions_present),
        ("numbers_match_data", assert_numbers_match_data),
        ("pymol_panels", assert_pymol_panels),
        ("supplementary", assert_supplementary),
        ("shared_palette", assert_shared_palette),
        ("seeds_fixed", assert_seeds_fixed),
        ("reproducible", assert_reproducible),
    ]
    for name, fn in checks:
        try:
            print(f"  [PASS] {name}: {fn()}")
        except GateFailure as exc:
            _failures.append(name)
            print(f"  [FAIL] {name}: {exc}")
        except Exception as exc:  # noqa: BLE001
            _failures.append(name)
            print(f"  [FAIL] {name}: unexpected {type(exc).__name__}: {exc}")
    print()
    if _failures:
        print(f"GATE FAILED — {len(_failures)}/{len(checks)} failed: {', '.join(_failures)}")
        return 1
    print(f"GATE PASSED — {len(checks)}/{len(checks)} checks. Phase 6 is complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
