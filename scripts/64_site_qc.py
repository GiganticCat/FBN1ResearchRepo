#!/usr/bin/env python3
"""
64_site_qc.py — the per-construct calcium-site quality and admission table.

A reviewer reading the geometry result is entitled to see how the admission filter behaved
across the whole array, not only that one construct failed it. `40_af3_batch4_geometry.py`
applies the filter to the 17 constructs that carry a mutant and logs the outcome, and nothing
until now reported it for the other four or put any of it in the supplement.

This script runs the same three criteria over all 21 wild-type constructs, using the same
functions from script 40 rather than reimplementing them, so the numbers cannot drift from the
ones the analysis used.

  seed spread    the seed-to-seed spread of the site's bond-valence sum, as a percentage of the
                 site's own valence, must be at most 15%
  occupancy      the site must be found in all five wild-type models
  local pLDDT    mean AlphaFold confidence over the coordinating residues, at least 70

The margin column is the distance to the nearest of those three limits, expressed as a fraction
of the limit, so a site that scraped in is visible as a small positive number rather than only
as a pass.

Writes: data/processed/af3_site_qc.tsv, results/site_qc.md, logs/64_site_qc_<stamp>.log
"""
from __future__ import annotations

import csv
import importlib.util
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
AFDIR = ROOT / "structures" / "alphafold"
B2MAN = ROOT / "handoff/structures/inbox/batch2_domains/manifest.csv"
SITES = ROOT / "data/processed/af3_ca_sites.tsv"
GEOM = ROOT / "data/processed/af3_batch4_geometry.tsv"
OUT_TSV = ROOT / "data/processed/af3_site_qc.tsv"
OUT_MD = ROOT / "results/site_qc.md"


def _load(stem: str):
    spec = importlib.util.spec_from_file_location(f"_{stem}", ROOT / "scripts" / f"{stem}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main() -> int:
    g40 = _load("40_af3_batch4_geometry")
    ref = g40.reference_seq()

    dirs = {g40.norm(d.name): d for d in AFDIR.iterdir() if d.is_dir()}
    b2 = list(csv.DictReader(open(B2MAN)))
    sites_tab = pd.read_csv(SITES, sep="\t")
    geom = pd.read_csv(GEOM, sep="\t")
    n_var = (geom[geom.is_cognate]
             .groupby(["construct", "site"]).size().rename("n_variants_at_site"))

    log = [f"Per-site quality and admission — {datetime.now(timezone.utc).isoformat()}",
           f"  criteria  seed spread <= {g40.MAX_WT_SD_PCT}% of valence, "
           f"occupancy {g40.MIN_WT_MODELS}/5 models, local pLDDT >= {g40.MIN_WT_PLDDT}", ""]

    rows = []
    for meta in b2:
        cname = meta["job_name"]
        d = dirs.get(g40.norm(cname)) or dirs.get(g40.norm("fold_" + cname))
        if d is None:
            log.append(f"  {cname}: no downloaded folder — skipped")
            continue
        lo, hi = (int(x) for x in meta["uniprot_range"].split("-"))
        m = g40.measure_wt(d, lo - 1)
        if m is None:
            log.append(f"  {cname}: fewer than five models — skipped")
            continue
        if m["seq"] != ref[lo - 1:hi]:
            raise SystemExit(f"{cname}: construct sequence does not match P35555 {lo}-{hi}")

        for j, s in sorted(m["sites"].items()):
            occ = len(s["bvs"])
            bvs = g40.med(s["bvs"])
            sd_pct = 100.0 * s["sd"] / bvs if bvs else float("nan")
            plddt = g40.med(s["plddt"])
            shift = g40.med(s["self_shift"]) if s["self_shift"] else float("nan")
            # how much room each criterion had left, as a fraction of its own limit
            margin = min((g40.MAX_WT_SD_PCT - sd_pct) / g40.MAX_WT_SD_PCT,
                         (occ - g40.MIN_WT_MODELS + 1) / g40.MIN_WT_MODELS if occ >= 5
                         else (occ - g40.MIN_WT_MODELS) / g40.MIN_WT_MODELS,
                         (plddt - g40.MIN_WT_PLDDT) / g40.MIN_WT_PLDDT)
            rows.append(dict(
                construct=meta["construct"], job_name=cname,
                cbegf_domains=meta["cbegf_domains"].strip('"'),
                uniprot_range=meta["uniprot_range"], site=j,
                ligand_positions_uniprot=";".join(str(x + lo - 1) for x in sorted(s["sig"])),
                n_models_present=occ, plddt_site=round(plddt, 1),
                plddt_construct=round(m["plddt"], 1),
                bvs_wt=round(bvs, 3), cn_wt=round(g40.med(s["cn"]), 1),
                seed_sd_bvs=round(s["sd"], 3), seed_sd_pct=round(sd_pct, 1),
                ion_seed_shift_A=round(shift, 2) if shift == shift else "",
                pass_seed_spread=bool(sd_pct <= g40.MAX_WT_SD_PCT),
                pass_occupancy=bool(occ >= g40.MIN_WT_MODELS),
                pass_plddt=bool(plddt >= g40.MIN_WT_PLDDT),
                admitted=bool(s["ok"]), margin_to_nearest_limit=round(margin, 3)))

    qc = pd.DataFrame(rows).sort_values(["construct", "site"]).reset_index(drop=True)

    # cross-check against the two tables that already carry some of these numbers
    chk = qc.merge(sites_tab[["construct", "site", "bvs_median", "bvs_seed_sd_pct"]],
                   left_on=["job_name", "site"], right_on=["construct", "site"], how="left",
                   suffixes=("", "_tab"))
    bad = chk[chk.bvs_median.notna()
              & (((chk.bvs_wt - chk.bvs_median).abs() > 0.002)
                 | ((chk.seed_sd_pct - chk.bvs_seed_sd_pct).abs() > 0.15))]
    if len(bad):
        raise SystemExit(f"{len(bad)} sites disagree with data/processed/af3_ca_sites.tsv")
    b4 = geom.drop_duplicates(["construct", "site"])[["construct", "site", "wt_site_ok",
                                                      "plddt_site_wt"]]
    chk2 = qc.merge(b4, left_on=["job_name", "site"], right_on=["construct", "site"],
                    how="inner", suffixes=("", "_b4"))
    if not (chk2.admitted == chk2.wt_site_ok).all():
        raise SystemExit("admission calls disagree with af3_batch4_geometry.tsv")
    if (chk2.plddt_site - chk2.plddt_site_wt).abs().max() > 0.05:
        raise SystemExit("site pLDDT disagrees with af3_batch4_geometry.tsv")
    log.append(f"cross-check  {len(chk2)} sites agree with the batch-4 wild-type measurements, "
               f"{int(chk.bvs_median.notna().sum())} agree with af3_ca_sites.tsv")

    qc = qc.merge(n_var.reset_index().rename(columns={"construct": "job_name"}),
                  on=["job_name", "site"], how="left")
    qc["n_variants_at_site"] = qc["n_variants_at_site"].fillna(0).astype(int)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(OUT_TSV, sep="\t", index=False)

    n, ok = len(qc), int(qc.admitted.sum())
    tight = qc[qc.admitted & (qc.margin_to_nearest_limit < 0.20)]
    log += ["", f"{n} calcium sites across {qc.construct.nunique()} constructs, "
                f"{ok} admitted and {n - ok} excluded",
            f"{len(tight)} admitted sites sit within 20% of a limit"]
    for _, r in qc[~qc.admitted].iterrows():
        why = [k.replace("pass_", "") for k in ("pass_seed_spread", "pass_occupancy",
                                                "pass_plddt") if not r[k]]
        log.append(f"  EXCLUDED  {r.construct} site {r.site}  on {', '.join(why)}  "
                   f"(spread {r.seed_sd_pct}%, {r.n_models_present}/5 models, "
                   f"pLDDT {r.plddt_site})  {r.n_variants_at_site} variant(s) dropped")
    for _, r in tight.iterrows():
        log.append(f"  borderline  {r.construct} site {r.site}  margin "
                   f"{r.margin_to_nearest_limit:.2f}  (spread {r.seed_sd_pct}%, "
                   f"pLDDT {r.plddt_site})")

    hdr = ["construct", "cbegf_domains", "site", "n_models_present", "plddt_site",
           "seed_sd_pct", "bvs_wt", "cn_wt", "ion_seed_shift_A", "admitted",
           "margin_to_nearest_limit", "n_variants_at_site"]
    md = ["# Calcium-site quality and admission, all 21 wild-type constructs", "",
          f"Produced by `scripts/64_site_qc.py` on {datetime.now(timezone.utc).date()}. "
          "A site is admitted on wild-type evidence alone, fixed before any mutant was read. "
          f"Admission requires a seed spread no greater than {g40.MAX_WT_SD_PCT:.0f}% of the "
          f"site's own valence, occupancy in all {g40.MIN_WT_MODELS} wild-type models, and a "
          f"mean pLDDT of at least {g40.MIN_WT_PLDDT:.0f} over the coordinating residues.", "",
          f"{ok} of {n} sites are admitted. {len(tight)} of the admitted sites lie within 20% "
          "of one of the three limits and are listed as borderline below.", "",
          "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for _, r in qc.iterrows():
        md.append("| " + " | ".join(str(r[c]) for c in hdr) + " |")
    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    (ROOT / "logs" / f"64_site_qc_{stamp}.log").write_text("\n".join(log) + "\n",
                                                           encoding="utf-8")
    print("\n".join(log))
    print(f"\n  wrote {OUT_TSV.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
