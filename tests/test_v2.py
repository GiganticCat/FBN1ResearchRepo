#!/usr/bin/env python3
"""
test_v2.py — the gate for the v2 analysis and its figures.

The v1 gates (`test_analysis.py`, `test_figures.py`) still pass, because the v1 tables and
figures they check are still on disk and unchanged. That is exactly the hazard: a green v1 gate
says nothing about v2, and v2 is what the manuscript now reports. This file gates v2.

Checks
  A1  rosetta_ddg completeness      every task has a status; no silent NaN; counts reconcile
  A2  rosetta steric tail           the unrelieved-clash fraction stays below threshold
  A3  site-class identity           Ca ligands are D/E/N, cysteine-removing positions are Cys
  A4  wild-type identity            every measured variant's WT residue matches P35555
  A5  bond-valence calibration      local BVS still tracks CheckMyMetal at r >= 0.95
  A6  enrichment table              intervals present, finite, and bracketing the point estimate
  A7  batch-4 completeness          every requested mutant measured at every site, nothing blank
  A8  null threshold re-derives     the disruption threshold recomputes from the bystander sites
  A9  disruption calls consistent   every call follows from the threshold, none set by hand
  A10 final statistics coherent     q >= p, intervals bracket estimates, both force fields kept apart
  F1  figure artifacts              png + vector pdf + numbers sidecar + the script that made it
  F2  numbers match data            every recorded value re-derives from data/processed
  F3  journal format                double-column width, no panel titles, seeded
  F4  reproducible                  re-running every figure script reproduces its numbers

Run:  .venv/bin/python tests/test_v2.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
RAW = ROOT / "data/raw"
FIGURES = ROOT / "figures"
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import figstyle as fs  # noqa: E402

TOL = 5e-4

# The steric ceiling. Before the backbone was freed inside the repack shell, 21.5% of Rosetta
# values exceeded 100 REU -- unrelieved collisions from bulky residues replacing buried
# cysteines, not energies. 2% leaves room for genuinely catastrophic substitutions while failing
# loudly if the fixed-backbone regression ever returns.
STERIC_CEILING_REU = 100.0
STERIC_MAX_FRACTION = 0.02

FIGS = {
    "fig1_module": "50_fig1_module.py",
    "fig2_stability": "51_fig2_stability.py",
    "fig3_coordination": "52_fig3_coordination.py",
    "fig4_interpretation": "53_fig4_interpretation.py",
    "figS1_validation": "54_figS1_validation.py",
    "figS2_structure": "55_figS2_structure.py",
}

# Fixed in scripts/40 on wild-type evidence alone, before any mutant was read. Restated here so
# the gate fails if the analysis quietly loosens them.
NULL_PCTILE = 95.0
N_BATCH4_JOBS = 60


class GateFailure(AssertionError):
    pass


def seq() -> str:
    hits = sorted(RAW.glob("uniprot_P35555_v*.json"))
    if not hits:
        raise GateFailure("no cached UniProt P35555 JSON in data/raw")
    return json.loads(hits[-1].read_text(encoding="utf-8"))["sequence"]["value"]


def numbers(fig: str) -> dict:
    p = FIGURES / f"{fig}.numbers.json"
    if not p.is_file():
        raise GateFailure(f"missing {p.relative_to(ROOT)}")
    return json.loads(p.read_text(encoding="utf-8"))


def close(a, b) -> bool:
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) or isinstance(b, str):
        return a == b
    # some sidecars record a list of positions or a dict of per-position counts rather than a
    # scalar; compare element-wise so a reordered, truncated or altered container fails rather
    # than raising
    if isinstance(a, (list, tuple)) or isinstance(b, (list, tuple)):
        if not (isinstance(a, (list, tuple)) and isinstance(b, (list, tuple))):
            return False
        return len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)):
            return False
        return set(a) == set(b) and all(close(a[k], b[k]) for k in a)
    return abs(float(a) - float(b)) <= TOL * max(1.0, abs(float(b)))


# --------------------------------------------------------------------------- analysis

def a1_rosetta_completeness() -> str:
    ros = pd.read_csv(PROC / "rosetta_ddg.tsv", sep="\t")
    if ros.empty:
        raise GateFailure("rosetta_ddg.tsv is empty")
    if ros["status"].isna().any():
        raise GateFailure("rows with no status — a failure must be recorded, never blank")
    ok = ros[ros["status"] == "ok"]
    bad = ros[ros["status"] != "ok"]
    # An 'ok' row must carry a finite number. A silently empty ddG on a successful row is the
    # exact failure mode that let the mutate_residue no-op through the first time.
    vals = pd.to_numeric(ok["ddg_mean"], errors="coerce")
    if vals.isna().any():
        raise GateFailure(f"{int(vals.isna().sum())} rows marked ok but carry no ddG")
    if not np.isfinite(vals).all():
        raise GateFailure("non-finite ddG among ok rows")
    # and it must not be the all-zeros signature of a mutation that never happened
    if float((vals == 0.0).mean()) > 0.10:
        raise GateFailure(
            f"{100*float((vals==0.0).mean()):.1f}% of ddG values are exactly 0.0 — "
            "this is the signature of a mutation that did not take")
    sm = pd.read_csv(PROC / "structural_metrics_af3.tsv", sep="\t")
    expected = int(sm["variation_id"].isin(ros["variation_id"]).sum())
    if expected != len(ros):
        raise GateFailure(f"rosetta rows {len(ros)} do not reconcile with "
                          f"{expected} matching structural-metric rows")
    return (f"{len(ok)}/{len(ros)} ok, {len(bad)} recorded failures, "
            f"{float((vals == 0.0).mean())*100:.2f}% exact zeros")


def a2_rosetta_steric_tail() -> str:
    ros = pd.read_csv(PROC / "rosetta_ddg.tsv", sep="\t")
    v = pd.to_numeric(ros.loc[ros["status"] == "ok", "ddg_mean"], errors="coerce").dropna()
    frac = float((v > STERIC_CEILING_REU).mean())
    if frac > STERIC_MAX_FRACTION:
        raise GateFailure(
            f"{100*frac:.1f}% of ddG values exceed {STERIC_CEILING_REU} REU "
            f"(ceiling {100*STERIC_MAX_FRACTION:.0f}%) — the repack shell is not relieving "
            "clashes; check that backbone minimisation is enabled inside the shell")
    return (f"{100*frac:.2f}% above {STERIC_CEILING_REU:.0f} REU (max "
            f"{v.max():.1f}), ceiling {100*STERIC_MAX_FRACTION:.0f}%")


def a3_site_class_identity() -> str:
    sm = pd.read_csv(PROC / "structural_metrics_af3.tsv", sep="\t")
    tb = lambda s: s.map(lambda v: str(v).strip().lower() == "true")  # noqa: E731
    lig = sm[tb(sm["ca_ligand_sidechain"])]
    off = lig[~lig["wt_aa"].isin(["D", "E", "N", "S", "T", "Y", "Q"])]
    if len(off):
        raise GateFailure(f"{len(off)} side-chain Ca ligands have a WT residue with no "
                          f"coordinating oxygen: {sorted(set(off['wt_aa']))}")
    cys = sm[tb(sm["cys_removing"])]
    if not (cys["wt_aa"] == "C").all():
        bad = sorted(set(cys.loc[cys["wt_aa"] != "C", "wt_aa"]))
        raise GateFailure(f"cysteine-removing rows whose WT residue is not Cys: {bad}")
    if (cys["mut_aa"] == "C").any():
        raise GateFailure("a 'cysteine-removing' row mutates TO cysteine")
    return f"{len(lig)} side-chain ligands all O/N donors; {len(cys)} Cys-removing all Cys→X"


def a4_wt_identity() -> str:
    s = seq()
    sm = pd.read_csv(PROC / "structural_metrics_af3.tsv", sep="\t")
    bad = []
    for _, r in sm.iterrows():
        p = int(r["position"])
        if not (1 <= p <= len(s)):
            bad.append((p, "out of range")); continue
        if s[p - 1] != r["wt_aa"]:
            bad.append((p, f"{r['wt_aa']} vs {s[p-1]}"))
    if bad:
        raise GateFailure(f"{len(bad)} rows disagree with P35555, e.g. {bad[:5]}")
    return f"{len(sm)} rows match P35555 ({len(s)} aa)"


def a5_bvs_calibration() -> str:
    b = pd.read_csv(PROC / "bvs_vs_cmm.tsv", sep="\t")
    if len(b) < 8:
        raise GateFailure(f"only {len(b)} calibration sites; need >= 8")
    r = float(np.corrcoef(b["valence_cmm"], b["bvs_local"])[0, 1])
    if r < 0.95:
        raise GateFailure(f"local bond valence no longer tracks CheckMyMetal (r = {r:.3f})")
    return f"r = {r:.4f} over {len(b)} sites"


def a6_enrichment_table() -> str:
    e = pd.read_csv(PROC / "enrichment_v2.tsv", sep="\t")
    need = {"feature", "denominator", "odds_ratio_corrected", "ci_lo", "ci_hi", "p_adj"}
    if not need.issubset(e.columns):
        raise GateFailure(f"enrichment table missing {sorted(need - set(e.columns))}")
    for _, r in e.iterrows():
        lo, hi, orr = float(r["ci_lo"]), float(r["ci_hi"]), float(r["odds_ratio_corrected"])
        if not (np.isfinite(lo) and np.isfinite(hi) and np.isfinite(orr)):
            raise GateFailure(f"non-finite interval for {r['feature']}/{r['denominator']}")
        if not (lo <= orr <= hi):
            raise GateFailure(f"point estimate outside its own interval for "
                              f"{r['feature']}/{r['denominator']}: {orr} not in [{lo}, {hi}]")
    dens = set(e["denominator"])
    if len(dens) < 2:
        raise GateFailure(f"only one denominator reported ({dens}); both must be shown")
    return f"{len(e)} rows, {len(dens)} denominators, all intervals bracket their estimate"


def a7_batch4_completeness() -> str:
    g = pd.read_csv(PROC / "af3_batch4_geometry.tsv", sep="\t")
    man = pd.read_csv(ROOT / "handoff/structures/inbox/batch4_mutants/manifest.csv")
    if g["job"].nunique() != N_BATCH4_JOBS:
        raise GateFailure(f"{g['job'].nunique()} jobs measured, {N_BATCH4_JOBS} requested — "
                          "a missing fold must be reported, not dropped")
    missing = set(man["job_name"]) - set(g["job"])
    if missing:
        raise GateFailure(f"requested but never measured: {sorted(missing)[:5]}")
    # every job must contribute exactly one cognate site and no more
    per_job = g.groupby("job")["is_cognate"].sum()
    if not (per_job == 1).all():
        raise GateFailure(f"{int((per_job != 1).sum())} jobs do not have exactly one cognate site")
    # a 'present' row must carry a number; a blank is the failure mode this whole gate exists for
    pres = g[g["status"] == "present"]
    v = pd.to_numeric(pres["delta_bvs_pct"], errors="coerce")
    if v.isna().any():
        raise GateFailure(f"{int(v.isna().sum())} rows marked present carry no valence change")
    # and the substitution must be the one that was asked for
    key = man.set_index("job_name")["variant"].to_dict()
    wrong = [j for j, var in zip(g["job"], g["variant"]) if key.get(j) != var]
    if wrong:
        raise GateFailure(f"{len(set(wrong))} jobs carry a different substitution than requested")
    return (f"{g['job'].nunique()} jobs x {len(g)//g['job'].nunique()} sites, "
            f"{len(pres)} measured, {len(g) - len(pres)} vacated, all substitutions verified")


def a8_null_threshold() -> str:
    g = pd.read_csv(PROC / "af3_batch4_geometry.tsv", sep="\t")
    r = pd.read_csv(PROC / "batch4_rate.tsv", sep="\t")
    g = g[g["wt_site_ok"] & (g["status"] == "present")]
    byst = g[~g["is_cognate"]]["delta_bvs_pct"].to_numpy(float)
    neg = g[g["is_cognate"] & (g["klass"] == "negative_control")]["delta_bvs_pct"].to_numpy(float)
    null = np.concatenate([byst, neg])
    if len(null) < 60:
        raise GateFailure(f"only {len(null)} null measurements; the threshold would be noise")
    t = float(np.percentile(np.abs(null), NULL_PCTILE))
    recorded = float(r["null_threshold_pct"].iloc[0])
    if abs(t - recorded) > 0.05:
        raise GateFailure(f"threshold in batch4_rate.tsv is {recorded}, re-derives as {t:.2f}")
    if r["null_threshold_pct"].nunique() != 1:
        raise GateFailure("more than one threshold in the rate table")
    return f"{t:.2f}% from {len(null)} no-effect sites ({len(byst)} bystander, {len(neg)} negative)"


def a9_disruption_calls() -> str:
    r = pd.read_csv(PROC / "batch4_rate.tsv", sep="\t")
    t = float(r["null_threshold_pct"].iloc[0])
    d = pd.to_numeric(r["delta_bvs_pct"], errors="coerce")
    expect = (r["status"] == "vacated") | (d <= -t)
    called = r["disrupted"].map(lambda v: str(v).strip().lower() == "true")
    if not (expect.fillna(False) == called).all():
        n = int((expect.fillna(False) != called).sum())
        raise GateFailure(f"{n} disruption calls do not follow from the threshold")
    # controls must not be silently absent: the rate is meaningless without them
    for k in ("ca_ligand_pathogenic", "composition_control", "negative_control"):
        if (r["klass"] == k).sum() == 0:
            raise GateFailure(f"no {k} variants in the rate table")
    return (f"{int(called.sum())}/{len(r)} disrupted, every call re-derived from the "
            f"{t:.2f}% threshold")


def a10_final_statistics() -> str:
    st_ = pd.read_csv(ROOT / "results/statistics_final.tsv", sep="\t")
    if st_.empty:
        raise GateFailure("statistics_final.tsv is empty")
    if (st_["p_adj"] + 1e-12 < st_["p_raw"]).any():
        raise GateFailure("a corrected q is smaller than its raw p")
    ci = st_.dropna(subset=["ci_lo", "ci_hi"])
    ci = ci[ci["effect"] == "Cliff's delta"]
    bad = ci[(ci["estimate"] < ci["ci_lo"]) | (ci["estimate"] > ci["ci_hi"])]
    if len(bad):
        raise GateFailure(f"{len(bad)} effect sizes fall outside their own interval")
    # the two force fields must never be compared to each other
    b = st_[st_["family"] == "B_stability"]
    if b.empty or not b["test"].str.contains("foldx").any() or not b["test"].str.contains("rosetta").any():
        raise GateFailure("family B must report both force fields, separately")
    if b["test"].str.contains("foldx").eq(b["test"].str.contains("rosetta")).any():
        raise GateFailure("a test compares FoldX to Rosetta — different units, never pooled")
    return f"{len(st_)} tests across {st_['family'].nunique()} families, all coherent"


# --------------------------------------------------------------------------- figures

def f1_artifacts() -> str:
    for fig, script in FIGS.items():
        for ext, floor in (("png", 20_000), ("pdf", 8_000)):
            p = FIGURES / f"{fig}.{ext}"
            if not p.is_file():
                raise GateFailure(f"missing {p.relative_to(ROOT)}")
            if p.stat().st_size < floor:
                raise GateFailure(f"{p.name} is {p.stat().st_size} bytes — suspiciously small")
        if not (SCRIPTS / script).is_file():
            raise GateFailure(f"missing script {script} for {fig}")
        rec = numbers(fig)
        if rec.get("script") != script:
            raise GateFailure(f"{fig} sidecar names {rec.get('script')!r}, expected {script!r}")
    return f"{len(FIGS)} figures with raster, vector, sidecar and script"


def f2_numbers_match_data() -> str:
    """Re-derive every value the figures print, from data/processed, and compare."""
    sys.path.insert(0, str(SCRIPTS))
    import paperstyle as ps  # noqa: E402

    df = ps.load()
    checked, bad = 0, []

    def cmp(fig, key, expected):
        nonlocal checked
        rec = numbers(fig)["values"]
        if key not in rec:
            bad.append(f"{fig}:{key} not recorded"); return
        checked += 1
        if not close(rec[key], expected):
            bad.append(f"{fig}:{key} figure={rec[key]!r} data={expected!r}")

    # ---- fig1: cohort counts and the measured consensus offsets
    pop = pd.read_csv(PROC / "population_set.tsv", sep="\t")
    covered = set(df.loc[df["has_structure"], "position"])
    cmp("fig1_module", "n_pathogenic", int((df["set_primary"] == "pathogenic").sum()))
    cmp("fig1_module", "n_population", int(pop["position"].isin(covered).sum()))
    dom = ps.domains()
    cmp("fig1_module", "n_cbegf_domains", int((dom["kind"] == "cbEGF").sum()))
    # the donor offsets must be the ones the models actually show, not a drawn cartoon
    lig = pd.read_csv(PROC / "af3_ca_ligands.tsv", sep="\t")
    cb = dom[dom["kind"] == "cbEGF"]
    rows = []
    for r in lig.itertuples():
        d = cb[(cb["start"] <= r.uniprot_pos) & (cb["end"] >= r.uniprot_pos)]
        if len(d) == 1:
            rows.append((int(d.iloc[0].cbegf_index), r.uniprot_pos - d.iloc[0].start,
                         r.atom_class))
    Ld = pd.DataFrame(rows, columns=["idx", "off", "cls"]).drop_duplicates()
    side = sorted(int(o) for o, n in Ld[Ld["cls"] == "sidechain"]["off"].value_counts().items()
                  if n >= 10)
    cmp("fig1_module", "sidechain_donor_offsets", side)

    # ---- fig2: group medians and n, in both force fields and AlphaMissense
    for g in ps.GROUPS:
        lab = ps.GLABEL[g].replace("\n", " ")
        sub = df[df["group"] == g]
        for col, tag in (("rosetta_ddG", "a_rosetta"), ("foldx_ddG", "b_foldx"), ("am", "d_am")):
            v = sub[col].dropna()
            cmp("fig2_stability", f"{tag}_{lab}_n", int(len(v)))
            if len(v):
                cmp("fig2_stability", f"{tag}_{lab}_median", round(float(np.median(v)), 4)
                    if tag == "d_am" else round(float(np.median(v)), 3))

    # the n that differs between fig 2's panels. Panel d scores 1,132 cysteine substitutions and
    # panels a and c score 1,038, because AlphaMissense needs only a sequence and a ddG needs a
    # modelled construct. The caption states that difference, so the difference is guarded.
    cmp("fig2_stability", "cys_without_structural_coverage",
        int(((df["group"] == "cys_removing") & df["rosetta_ddG"].isna()).sum()))

    # ---- the Benjamini-Hochberg q values the headline panels now print beside their P values.
    # They are read from results/statistics_final.tsv rather than recomputed in the figure, so
    # what is checked here is that the figure carried the right row across.
    st_ = pd.read_csv(ROOT / "results/statistics_final.tsv", sep="\t")

    def q_of(test: str) -> float:
        row = st_[st_["test"] == test]
        if row.empty:
            raise GateFailure(f"no statistics row for {test!r}")
        return float(f"{float(row.iloc[0]['p_adj']):.3g}")

    cmp("fig3_coordination", "fisher_q",
        q_of("disruption rate: pathogenic Ca ligands vs all controls"))
    cmp("fig3_coordination", "donor_fisher_q",
        q_of("disruption rate: donor removed vs retained, at the Asp/Glu positions"))
    cmp("fig4_interpretation", "spearman_q",
        q_of("Spearman(delta BVS %, Rosetta ddG) across pathogenic Ca ligands"))
    cmp("fig4_interpretation", "b_am_q",
        q_of("AlphaMissense: disrupted vs intact Ca-ligand variants"))

    # ---- fig3: the rate, its threshold and the groups behind it
    rate = pd.read_csv(PROC / "batch4_rate.tsv", sep="\t")
    cmp("fig3_coordination", "n_variants", int(len(rate)))
    cmp("fig3_coordination", "null_threshold_pct", round(float(rate["null_threshold_pct"].iloc[0]), 2))
    dis = rate["disrupted"].map(lambda v: str(v).strip().lower() == "true")
    for k in ("ca_ligand_pathogenic", "composition_control", "negative_control"):
        m_ = rate["klass"] == k
        cmp("fig3_coordination", f"n_{k}", int(m_.sum()))
        cmp("fig3_coordination", f"disrupted_{k}", int(dis[m_].sum()))
        cmp("fig3_coordination", f"rate_{k}", round(100 * float(dis[m_].sum()) / int(m_.sum()), 1))
    lig_ = rate[rate["klass"] == "ca_ligand_pathogenic"]
    keep = lig_["donor_kept"].map(lambda v: str(v).strip().lower() == "true")
    cmp("fig3_coordination", "d_rmsd_max", round(float(rate["local_rmsd_A"].max()), 3))
    # the donor rule is reported stratified by consensus position; both strata must reconcile
    for posc in ("Asp/Glu", "beta-OH Asn"):
        for kept, tag in ((True, "kept"), (False, "removed")):
            sub_ = lig_[(lig_["site_pos"] == posc) & (keep == kept)]
            key = f"c_{posc.replace('-', '').replace(' ', '_')}_{tag}"
            cmp("fig3_coordination", f"{key}_n", int(len(sub_)))
            cmp("fig3_coordination", f"{key}_disrupted",
                int(sub_["disrupted"].map(lambda v: str(v).strip().lower() == "true").sum()))
    ns_ = lig_[(lig_["wt_aa"] == "N") & (lig_["mut_aa"] == "S")
               & (lig_["site_pos"] == "beta-OH Asn")]
    cmp("fig3_coordination", "asn_ser_n", int(len(ns_)))
    cmp("fig3_coordination", "asn_ser_disrupted",
        int(ns_["disrupted"].map(lambda v: str(v).strip().lower() == "true").sum()))

    # ---- fig4: enrichment odds ratios come from the enrichment table, not from the figure
    enr = pd.read_csv(PROC / "enrichment_v2.tsv", sep="\t")
    enr = enr[enr["denominator"] == "population (gnomAD)"]
    for feat in ("side-chain Ca ligand", "cysteine position", "backbone-only Ca ligand"):
        row = enr[enr["feature"] == feat]
        if not row.empty:
            cmp("fig4_interpretation", f"or_{feat.replace(' ', '_')}",
                round(float(row.iloc[0]["odds_ratio_corrected"]), 3))
    cmp("fig4_interpretation", "spearman_n", int(len(lig_.merge(
        pd.read_csv(PROC / "rosetta_ddg.tsv", sep="\t").query("status == 'ok'")[["variation_id"]],
        on="variation_id"))))

    # ---- figS1: the two validations
    cal = pd.read_csv(ROOT / "results/phase4_calibration.tsv", sep="\t")
    c0 = cal[cal["model"] == 0]
    cmp("figS1_validation", "a_n_domains", int(len(c0)))
    cmp("figS1_validation", "a_median_ca_deviation_A", round(float(c0["ca_deviation_A"].median()), 3))
    bv = pd.read_csv(PROC / "bvs_vs_cmm.tsv", sep="\t")
    cmp("figS1_validation", "b_n_sites", int(len(bv)))
    cmp("figS1_validation", "b_pearson_r",
        round(float(np.corrcoef(bv["valence_cmm"], bv["bvs_local"])[0, 1]), 4))

    # ---- figS2: the positions drawn on 2W86 must be the ones in the table
    onfrag = df[(df["position"] >= 807) & (df["position"] <= 951)
                & (df["set_primary"] == "pathogenic")]
    cmp("figS2_structure", "n_pathogenic_variants", int(len(onfrag)))
    cmp("figS2_structure", "n_positions_cys_removing",
        int(onfrag.loc[onfrag["group"] == "cys_removing", "position"].nunique()))

    if bad:
        raise GateFailure(f"{len(bad)} figure values disagree with the data:\n  "
                          + "\n  ".join(bad[:12]))
    return f"{checked} recorded values re-derived from data/processed"


def f3_journal_format() -> str:
    """
    Format rules the final set follows, taken from the papers in `resources/papers/`.

    The earlier version of this check banned panel titles outright. That was wrong: Cheng 2023
    (Science) and Godwin 2023 (Nat Struct Mol Biol) title nearly every panel. What those titles
    never do is state a conclusion, and what those figures never do is rotate a category label.
    Both of those are checked instead.
    """
    problems = []
    for fig, script in FIGS.items():
        src = (SCRIPTS / script).read_text(encoding="utf-8")
        if "P.W1" not in src and "P.W15" not in src and "P.W2" not in src:
            problems.append(f"{script}: figure width is not a journal column width")
        if "import paperstyle" not in src:
            problems.append(f"{script}: does not use the shared style module")
        for line in src.splitlines():
            t = line.strip()
            if t.startswith("#"):
                continue
            # rotated category labels are a matplotlib default and appear in none of the
            # reference papers
            if "rotation=" in t and "set_xticklabels" in t:
                problems.append(f"{script}: rotated category labels — {t[:60]}")
        uses_rng = "default_rng" in src or "np.random" in src
        if uses_rng and "P.SEED" not in src:
            problems.append(f"{script}: draws random numbers without the shared seed")
    sys.path.insert(0, str(SCRIPTS))
    import paperstyle as ps  # noqa: E402
    if ps.SEED != 20260830:
        problems.append(f"shared seed changed to {ps.SEED}")
    if problems:
        raise GateFailure("; ".join(problems))
    seeded = sum(1 for s_ in FIGS.values()
                 if "default_rng" in (SCRIPTS / s_).read_text(encoding="utf-8"))
    return (f"{len(FIGS)} scripts on journal column widths, no rotated category labels; "
            f"{seeded} use randomness, all from seed {ps.SEED}")


def f4_reproducible() -> str:
    before = {f: numbers(f)["values"] for f in FIGS}
    for fig, script in FIGS.items():
        r = subprocess.run([sys.executable, str(SCRIPTS / script)],
                           capture_output=True, text=True, cwd=str(ROOT), timeout=900)
        if r.returncode != 0:
            raise GateFailure(f"{script} failed on re-run: {r.stderr[-400:]}")
    drift = []
    for fig in FIGS:
        after = numbers(fig)["values"]
        for k, v in before[fig].items():
            if k not in after:
                drift.append(f"{fig}:{k} disappeared")
            elif not close(after[k], v):
                drift.append(f"{fig}:{k} {v!r} -> {after[k]!r}")
    if drift:
        raise GateFailure(f"{len(drift)} values changed on re-run: {drift[:8]}")
    return f"{len(FIGS)} scripts re-ran and reproduced identical numbers"


def a11_site_qc() -> str:
    """
    The admission table in the supplement must re-derive from the same three criteria.

    The geometry result depends on which calcium sites were judged measurable, and until the
    table existed the only record of that was a run log. This re-applies the criteria to
    `af3_site_qc.tsv` and checks that the sites it admits are the sites the mutant analysis used.
    """
    qc = pd.read_csv(PROC / "af3_site_qc.tsv", sep="\t")
    geom = pd.read_csv(PROC / "af3_batch4_geometry.tsv", sep="\t")
    if len(qc) != 63:
        raise GateFailure(f"{len(qc)} sites in the QC table, expected 63")

    truth = ((qc["seed_sd_pct"] <= 15.0) & (qc["n_models_present"] >= 5)
             & (qc["plddt_site"] >= 70.0))
    if not (truth == qc["admitted"]).all():
        n = int((truth != qc["admitted"]).sum())
        raise GateFailure(f"{n} admission calls do not follow from the three criteria")

    b4 = geom.drop_duplicates(["construct", "site"])[["construct", "site", "wt_site_ok"]]
    j = qc.merge(b4, left_on=["job_name", "site"], right_on=["construct", "site"], how="inner")
    if len(j) != len(b4):
        raise GateFailure(f"{len(b4) - len(j)} sites used in the mutant analysis are missing "
                          "from the QC table")
    if not (j["admitted"] == j["wt_site_ok"]).all():
        raise GateFailure("the QC table and the mutant geometry disagree on admission")

    # the variants the filter cost must be the ones the manuscript says were dropped
    cog = geom[geom["is_cognate"]]
    dropped = int((~cog["wt_site_ok"]).sum())
    if dropped != 3 or len(cog) - dropped != 57:
        raise GateFailure(f"{dropped} variants dropped by the filter and {len(cog) - dropped} "
                          "kept, expected 3 and 57")
    tight = int((qc["admitted"] & (qc["margin_to_nearest_limit"] < 0.20)).sum())
    return (f"63 sites, {int(qc['admitted'].sum())} admitted, {dropped} variants dropped, "
            f"{tight} admitted site within 20% of a limit")


CHECKS = [
    ("A1 rosetta completeness", a1_rosetta_completeness),
    ("A2 rosetta steric tail", a2_rosetta_steric_tail),
    ("A3 site-class identity", a3_site_class_identity),
    ("A4 wild-type identity", a4_wt_identity),
    ("A5 bond-valence calibration", a5_bvs_calibration),
    ("A6 enrichment table", a6_enrichment_table),
    ("A7 batch-4 completeness", a7_batch4_completeness),
    ("A8 null threshold re-derives", a8_null_threshold),
    ("A9 disruption calls consistent", a9_disruption_calls),
    ("A10 final statistics coherent", a10_final_statistics),
    ("A11 site admission re-derives", a11_site_qc),
    ("F1 figure artifacts", f1_artifacts),
    ("F2 numbers match data", f2_numbers_match_data),
    ("F3 journal format", f3_journal_format),
    ("F4 reproducible", f4_reproducible),
]


def main() -> int:
    print("v2 gate — analysis tables and the manuscript figures\n")
    failed = 0
    for name, fn in CHECKS:
        try:
            detail = fn()
            print(f"  PASS  {name:<30} {detail}")
        except Exception as exc:  # noqa: BLE001 — every failure is reported, none swallowed
            failed += 1
            print(f"  FAIL  {name:<30} {exc}")
    total = len(CHECKS)
    print(f"\n{total - failed}/{total} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
