#!/usr/bin/env python3
"""
06_annotate_domains.py — PHASE 3. Map variants onto domains and structural site classes.

The cbEGF calcium-binding consensus is DERIVED, not hard-coded from a paper. Each cbEGF domain
is anchored on its six cysteines, and the consensus positions are defined by offsets from those
anchors. The derivation was validated three independent ways before being used:

  1. STRUCTURE (ground truth). Ca2+ ligands were measured directly from the deposited
     coordinates of 2W86, 1LMJ, 1UZJ and 1EMN (8 cbEGF domains, 16 Ca sites). Every site shows
     the same pattern: sidechain ligands at C1-4 (Asp), C1-1 (Glu) and C3+2 (Asn), plus backbone
     carbonyls at C1-3, C3+3 and C3+6.
  2. SEQUENCE. Across all 43 cbEGF domains the derived positions carry exactly the residues the
     ClinGen FBN1 VCEP consensus [D]-X-[D/N]-[E/H]-Xm-[D/N]-Xn-[Y/F] specifies:
     C1-4 = D (43/43), C1-2 = D or N only, C1-1 = E (43/43), C3+2 = N (42) or D (1),
     C4-2 = Y or F only. No exceptions.
  3. LITERATURE. All six of Baudhuin 2019's published "tolerated Asn>Ser" positions
     (615, 1030, 1282, 1489, 2526, 2650) land exactly on the derived C1-2, and their three
     D<->N positions (1197, 1240, 1907) land on derived consensus positions. 9/9.

Rules follow the ClinGen FBN1 VCEP (Drackley et al. 2024, Genome Med 16:154):
  * cysteine-REMOVING in any cbEGF domain -> PM1_Strong (258 such cysteines)
  * cysteine-creating in cbEGF, or any cysteine change in EGF/TB/hybrid -> PM1 moderate
  * calcium-consensus residues, beta-hydroxylation site, and interdomain-packing glycines
    (between C2-C3, or between C3-C4 where an upstream cbEGF exists) -> PM1 moderate
  * Asn>Ser at the second D/N position is TOLERATED -> PM1 not applied

Writes:
  data/interim/cbegf_sites.tsv          per-domain derived site positions
  data/interim/disulfide_pairing.tsv    UniProt disulfides classified as C1-C3 / C2-C4 / C5-C6
  data/processed/variant_annotated.tsv  the master table + domain/site annotation
  results/phase3_crosstab.md            class x feature cross-tabulations
  logs/06_annotate_<stamp>.log
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW, INTERIM, PROC = ROOT / "data/raw", ROOT / "data/interim", ROOT / "data/processed"
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")

# Offsets from the cysteine anchors. See the module docstring for how these were established.
CA_CONSENSUS = {          # offset spec -> (label, VCEP motif element, ligand role)
    ("C1", -4): ("ca_D1", "[D]", "sidechain_ligand"),
    ("C1", -2): ("ca_DN2", "[D/N]", "consensus_only"),
    ("C1", -1): ("ca_EH", "[E/H]", "sidechain_ligand"),
    ("C3", +2): ("ca_DN_bOH", "[D/N] beta-OH", "sidechain_ligand"),
    ("C4", -2): ("ca_YF", "[Y/F]", "aromatic_core"),
}
CA_BACKBONE = {("C1", -3): "ca_backbone_1", ("C3", +3): "ca_backbone_2", ("C3", +6): "ca_backbone_3"}

log_lines: list[str] = []


def log(m: str = "") -> None:
    print(m)
    log_lines.append(m)


def main() -> int:
    seq = json.loads(sorted(RAW.glob("uniprot_P35555_v*.json"))[-1]
                     .read_text(encoding="utf-8"))["sequence"]["value"]
    doms = sorted(csv.DictReader((INTERIM / "uniprot_domains.tsv").open(encoding="utf-8"),
                                 delimiter="\t"), key=lambda d: int(d["start"]))
    log(f"Phase 3 annotation — {datetime.now(timezone.utc).isoformat()}")
    log(f"reference P35555 {len(seq)} aa; {len(doms)} annotated domains\n")

    # ---- derive per-domain site positions ---------------------------------
    site_of: dict[int, dict] = {}        # position -> annotation
    domain_of: dict[int, dict] = {}      # position -> domain record
    cys_role: dict[int, str] = {}        # position -> C1..C6 within its cbEGF
    site_rows = []
    cbegf_cys: set[int] = set()
    packing_gly: set[int] = set()

    for di, d in enumerate(doms):
        st, en = int(d["start"]), int(d["end"])
        for p in range(st, en + 1):
            domain_of[p] = d
        if d["kind"] != "cbEGF":
            continue
        s = d["sequence"]
        C = [i for i, c in enumerate(s) if c == "C"]
        if len(C) != 6:
            raise SystemExit(f"STOP: cbEGF{d['cbegf_index']} has {len(C)} cysteines, expected 6")
        anchors = {"C1": C[0], "C2": C[1], "C3": C[2], "C4": C[3], "C5": C[4], "C6": C[5]}

        for i, ci in enumerate(C, 1):
            cbegf_cys.add(st + ci)
            cys_role[st + ci] = f"C{i}"

        rec = {"cbegf_index": int(d["cbegf_index"]), "start": st, "end": en,
               "cys": ";".join(str(st + c) for c in C)}
        for (anchor, off), (label, motif, role) in CA_CONSENSUS.items():
            p = st + anchors[anchor] + off
            site_of[p] = {"site_class": "ca_consensus", "site_label": label,
                          "vcep_motif": motif, "ca_ligand_role": role,
                          "cbegf_index": int(d["cbegf_index"])}
            rec[label] = f"{p}:{seq[p-1]}"
        for (anchor, off), label in CA_BACKBONE.items():
            p = st + anchors[anchor] + off
            site_of.setdefault(p, {"site_class": "ca_backbone", "site_label": label,
                                   "vcep_motif": "", "ca_ligand_role": "backbone_ligand",
                                   "cbegf_index": int(d["cbegf_index"])})
        # Interdomain-packing glycines (VCEP): between C2-C3 always; between C3-C4 only when
        # an upstream cbEGF domain exists, since the contact is with the preceding domain.
        for i in range(C[1] + 1, C[2]):
            if s[i] == "G":
                packing_gly.add(st + i)
        if di > 0 and doms[di - 1]["kind"] == "cbEGF":
            for i in range(C[2] + 1, C[3]):
                if s[i] == "G":
                    packing_gly.add(st + i)
        site_rows.append(rec)

    log(f"derived cbEGF sites: {len(cbegf_cys)} cysteines, "
        f"{sum(1 for v in site_of.values() if v['site_class'] == 'ca_consensus')} consensus "
        f"positions, {len(packing_gly)} packing glycines")
    if len(cbegf_cys) != 258:
        raise SystemExit(f"STOP: {len(cbegf_cys)} cbEGF cysteines, VCEP states 258")
    log("  cbEGF cysteine count matches the VCEP's stated 258 exactly")

    # ---- disulfide pairing validation -------------------------------------
    dis = list(csv.DictReader((INTERIM / "uniprot_disulfides.tsv").open(encoding="utf-8"),
                              delimiter="\t"))
    partner: dict[int, int] = {}
    pair_rows, pattern = [], Counter()
    for r in dis:
        a, b = int(r["cys1"]), int(r["cys2"])
        partner[a], partner[b] = b, a
        ra, rb = cys_role.get(a, ""), cys_role.get(b, "")
        both_cb = a in cbegf_cys and b in cbegf_cys
        label = f"{ra}-{rb}" if both_cb else "non-cbEGF"
        pattern[label] += 1
        pair_rows.append({"cys1": a, "cys2": b, "role1": ra, "role2": rb,
                          "in_cbegf": both_cb, "pairing": label})
    log(f"\ndisulfide pairing (UniProt annotations, n={len(dis)}):")
    for k, v in pattern.most_common():
        log(f"    {k:12} {v}")
    expected = {"C1-C3", "C2-C4", "C5-C6"}
    cb_labels = {k for k in pattern if k != "non-cbEGF"}
    unexpected = cb_labels - expected
    if unexpected:
        log(f"  !! cbEGF disulfides not following the 1-3/2-4/5-6 rule: {sorted(unexpected)}")
    else:
        log("  all cbEGF disulfides follow the canonical C1-C3, C2-C4, C5-C6 pattern")

    with (INTERIM / "disulfide_pairing.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(pair_rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(pair_rows)

    cols = ["cbegf_index", "start", "end", "cys"] + [v[0] for v in CA_CONSENSUS.values()]
    with (INTERIM / "cbegf_sites.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(site_rows)

    # ---- neonatal region (Jensen 2012: TB3 -> cbEGF18) --------------------
    tb3 = next(d for d in doms if d["uniprot_description"].strip() == "TB 3")
    cb18 = next(d for d in doms if d["kind"] == "cbEGF" and int(d["cbegf_index"]) == 18)
    neo = (int(tb3["start"]), int(cb18["end"]))
    log(f"\nneonatal region (TB3 -> cbEGF18): {neo[0]}-{neo[1]}")

    # ---- annotate variants ------------------------------------------------
    variants = list(csv.DictReader((PROC / "variant_master.tsv").open(encoding="utf-8"),
                                   delimiter="\t"))
    log(f"\nannotating {len(variants):,} missense variants")

    for v in variants:
        p, wt, mut = int(v["position"]), v["wt_aa"], v["mut_aa"]
        d = domain_of.get(p)
        v["domain"] = d["uniprot_description"] if d else ""
        v["domain_kind"] = d["kind"] if d else "outside_domain"
        v["cbegf_index"] = d["cbegf_index"] if d and d["kind"] == "cbEGF" else ""
        v["in_cbegf"] = d is not None and d["kind"] == "cbEGF"
        v["in_neonatal_region"] = neo[0] <= p <= neo[1]

        site = site_of.get(p, {})
        v["site_class"] = site.get("site_class", "")
        v["site_label"] = site.get("site_label", "")
        v["vcep_motif"] = site.get("vcep_motif", "")
        v["ca_ligand_role"] = site.get("ca_ligand_role", "")
        v["is_ca_consensus"] = site.get("site_class") == "ca_consensus"
        v["is_ca_sidechain_ligand"] = site.get("ca_ligand_role") == "sidechain_ligand"

        # cysteine handling
        v["is_cysteine_site"] = wt == "C"
        v["cys_role"] = cys_role.get(p, "")
        v["disulfide_partner"] = partner.get(p, "") if wt == "C" else ""
        v["cys_removing"] = wt == "C" and mut != "C"
        v["cys_creating"] = wt != "C" and mut == "C"
        v["in_cbegf_cys"] = p in cbegf_cys
        v["is_packing_glycine"] = p in packing_gly and wt == "G"

        # VCEP PM1
        tolerated = (site.get("site_label") == "ca_DN2" and wt == "N" and mut == "S")
        v["vcep_tolerated_NtoS"] = tolerated
        if tolerated:
            cat, strength = "tolerated_N>S_at_second_D/N", "not_applied"
        elif v["cys_removing"] and v["in_cbegf_cys"]:
            cat, strength = "cbEGF_cysteine_removing", "PM1_Strong"
        elif v["cys_removing"]:
            cat, strength = "non_cbEGF_cysteine_removing", "PM1"
        elif v["cys_creating"] and v["in_cbegf"]:
            cat, strength = "cbEGF_cysteine_creating", "PM1"
        elif v["is_ca_consensus"]:
            cat, strength = "calcium_consensus_residue", "PM1"
        elif v["is_packing_glycine"]:
            cat, strength = "interdomain_packing_glycine", "PM1"
        elif v["cys_creating"]:
            cat, strength = "cysteine_creating_other", "PM1"
        else:
            cat, strength = "", ""
        v["vcep_pm1_category"] = cat
        v["vcep_pm1_strength"] = strength
        v["is_critical_residue"] = bool(strength) and strength != "not_applied"

    # ---- validation --------------------------------------------------------
    bad_cys = [v["variation_id"] for v in variants if v["cys_role"] and v["wt_aa"] != "C"]
    if bad_cys:
        raise SystemExit(f"STOP: positions annotated as cbEGF cysteines are not Cys: {bad_cys[:5]}")
    log("  every position annotated with a cysteine role has Cys as its reference residue")

    bad_motif = []
    for p, s in site_of.items():
        if s["site_class"] != "ca_consensus":
            continue
        aa, m = seq[p - 1], s["vcep_motif"]
        ok = ((m == "[D]" and aa == "D") or (m == "[D/N]" and aa in "DN")
              or (m == "[E/H]" and aa in "EH") or (m == "[D/N] beta-OH" and aa in "DN")
              or (m == "[Y/F]" and aa in "YF"))
        if not ok:
            bad_motif.append((p, aa, m))
    if bad_motif:
        raise SystemExit(f"STOP: consensus positions violating the VCEP motif: {bad_motif[:5]}")
    log("  all derived consensus positions satisfy the VCEP motif letters")

    # ---- cross-tabulations -------------------------------------------------
    def xtab(setcol: str) -> str:
        feats = [("in cbEGF domain", lambda v: v["in_cbegf"]),
                 ("Ca consensus residue", lambda v: v["is_ca_consensus"]),
                 ("  of which sidechain ligand", lambda v: v["is_ca_sidechain_ligand"]),
                 ("cbEGF cysteine (removing)",
                  lambda v: v["in_cbegf_cys"] and v["cys_removing"]),
                 ("cysteine-creating", lambda v: v["cys_creating"]),
                 ("packing glycine", lambda v: v["is_packing_glycine"]),
                 ("VCEP critical residue", lambda v: v["is_critical_residue"]),
                 ("neonatal region", lambda v: v["in_neonatal_region"])]
        sets = ["pathogenic", "benign", "vus", "conflicting"]
        n = {s: sum(1 for v in variants if v[setcol] == s) for s in sets}
        out = [f"| feature | " + " | ".join(f"{s} (n={n[s]})" for s in sets) + " |",
               "|---|" + "---|" * len(sets)]
        for label, fn in feats:
            cells = []
            for s in sets:
                sub = [v for v in variants if v[setcol] == s]
                k = sum(1 for v in sub if fn(v))
                cells.append(f"{k} ({k/len(sub)*100:.1f}%)" if sub else "0 (–)")
            out.append(f"| {label} | " + " | ".join(cells) + " |")
        return "\n".join(out)

    log("\nprimary tier (>=2 stars) cross-tab:")
    log(xtab("set_primary"))

    (ROOT / "results/phase3_crosstab.md").write_text(
        "# Phase 3 — class x structural-feature cross-tabulation\n\n"
        f"Generated by `scripts/06_annotate_domains.py` at "
        f"{datetime.now(timezone.utc).isoformat()}.\n"
        "Percentages are within-column (share of that class carrying the feature).\n\n"
        "## Primary tier (ClinVar review status >= 2 stars)\n\n" + xtab("set_primary") +
        "\n\n## Sensitivity tier (>= 1 star)\n\n" + xtab("set_sensitivity") + "\n",
        encoding="utf-8")

    # ---- write ------------------------------------------------------------
    cols = list({k: None for v in variants for k in v})
    with (PROC / "variant_annotated.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader(); w.writerows(variants)
    log(f"\nwrote data/processed/variant_annotated.tsv ({len(variants):,} x {len(cols)})")
    log("wrote data/interim/cbegf_sites.tsv, data/interim/disulfide_pairing.tsv")
    log("wrote results/phase3_crosstab.md")

    (ROOT / f"logs/06_annotate_{STAMP}.log").write_text("\n".join(log_lines) + "\n",
                                                        encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
