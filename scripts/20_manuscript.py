#!/usr/bin/env python3
"""
20_manuscript.py — Builds the project manuscript as a Word document.

Written as a script rather than by hand for the same reason every other deliverable here is:
so it regenerates, so the figures it embeds are the ones on disk, and so the numbers can be
checked against `data/processed/`. Numbers quoted in the text are pulled from the figure
`numbers.json` sidecars and the Phase 4/5 result tables wherever practical, and are marked in
this source with the constant they came from.

Every literature citation corresponds to a PDF actually present in `resources/papers/` and
verified in `resources/reference/lit_notes.md`. Tool citations that are NOT in that folder
(FoldX, DSSP, freesasa, PyMOL, SIFTS, Mutalyzer) are deliberately listed as software rather
than dressed up as read sources, and the gap is stated in the Limitations section.

Writes: results/FBN1_manuscript.docx, logs/20_manuscript_*.log
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs  # noqa: E402

OUT = fs.RESULTS / "FBN1_manuscript.docx"


def load_numbers() -> dict:
    """Every figure's recorded values, so the prose and the figures cannot disagree."""
    n = {}
    for p in sorted(fs.FIGURES.glob("*.numbers.json")):
        d = json.loads(p.read_text(encoding="utf-8"))
        n[d["figure"]] = d["values"]
    return n


def main() -> int:
    from docx import Document
    from docx.enum.section import WD_SECTION
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Inches, Pt, RGBColor

    log = ["Phase 6+ — manuscript build"]
    N = load_numbers()
    f1, f2, f3, f4, f5, f6 = (N["fig1_variant_landscape"], N["fig2_mechanism_dissociation"],
                              N["fig3_burial_control"], N["fig4_structure"],
                              N["fig5_ca_placement"], N["fig6_enrichment"])
    stats = pd.read_csv(fs.RESULTS / "phase5_statistics.tsv", sep="\t")

    def eff(test: str) -> pd.Series:
        return stats[stats.test == test].iloc[0]

    doc = Document()

    # ---- base styles ---------------------------------------------------------------
    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(11)
    normal.paragraph_format.space_after = Pt(8)
    normal.paragraph_format.line_spacing = 1.15

    for name, size, colour, before in (("Heading 1", 15, "1F3864", 16),
                                       ("Heading 2", 12.5, "2E5496", 12),
                                       ("Heading 3", 11.5, "2E5496", 10)):
        st = doc.styles[name]
        st.font.name = "Calibri"
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor.from_string(colour)
        st.paragraph_format.space_before = Pt(before)
        st.paragraph_format.space_after = Pt(4)

    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(0.9)
        s.left_margin = s.right_margin = Inches(1.0)

    # ---- helpers -------------------------------------------------------------------
    def h(text, level=1):
        doc.add_heading(text, level=level)

    def p(text, italic=False, size=None, align=None, space_after=None):
        par = doc.add_paragraph()
        run = par.add_run(text)
        run.italic = italic
        if size:
            run.font.size = Pt(size)
        par.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.JUSTIFY
        if space_after is not None:
            par.paragraph_format.space_after = Pt(space_after)
        return par

    def rich(*chunks):
        """Paragraph from (text, bold, italic) triples, for inline emphasis."""
        par = doc.add_paragraph()
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        for text, bold, italic in chunks:
            run = par.add_run(text)
            run.bold, run.italic = bold, italic
        return par

    def bullet(text):
        par = doc.add_paragraph(text, style="List Bullet")
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        par.paragraph_format.space_after = Pt(4)
        return par

    def numbered(text):
        par = doc.add_paragraph(text, style="List Number")
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        par.paragraph_format.space_after = Pt(4)
        return par

    def table(headers, rows, widths=None):
        t = doc.add_table(rows=1, cols=len(headers))
        t.style = "Light Grid Accent 1"
        for i, head in enumerate(headers):
            cell = t.rows[0].cells[i]
            cell.text = ""
            run = cell.paragraphs[0].add_run(head)
            run.bold = True
            run.font.size = Pt(9.5)
        for row in rows:
            cells = t.add_row().cells
            for i, val in enumerate(row):
                cells[i].text = ""
                run = cells[i].paragraphs[0].add_run(str(val))
                run.font.size = Pt(9.5)
        if widths:
            for r in t.rows:
                for i, w in enumerate(widths):
                    r.cells[i].width = Inches(w)
        doc.add_paragraph().paragraph_format.space_after = Pt(4)
        return t

    def figure(stem, caption, width=6.4):
        png = fs.FIGURES / f"{stem}.png"
        if not png.is_file():
            raise SystemExit(f"missing figure {png}")
        doc.add_picture(str(png), width=Inches(width))
        doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        cap = doc.add_paragraph()
        cap.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        run = cap.add_run(caption)
        run.font.size = Pt(9)
        run.italic = True
        cap.paragraph_format.space_after = Pt(14)
        log.append(f"  embedded {png.name}")

    def callout(title, body):
        t = doc.add_table(rows=1, cols=1)
        t.style = "Light Shading Accent 1"
        cell = t.rows[0].cells[0]
        cell.text = ""
        head = cell.paragraphs[0]
        r = head.add_run(title)
        r.bold = True
        r.font.size = Pt(10.5)
        for para in body:
            q = cell.add_paragraph()
            q.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
            rr = q.add_run(para)
            rr.font.size = Pt(10)
        doc.add_paragraph().paragraph_format.space_after = Pt(6)

    # ================================================================================
    # TITLE
    # ================================================================================
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    tr = title.add_run("Two Marfan mechanisms, two kinds of damage: calcium-site and "
                       "cysteine variants in fibrillin-1 are biophysically dissociable")
    tr.bold = True
    tr.font.size = Pt(17)
    tr.font.color.rgb = RGBColor.from_string("1F3864")
    title.paragraph_format.space_after = Pt(6)

    sub = doc.add_paragraph()
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    sr = sub.add_run("A computational structural-genomics analysis of 4,451 FBN1 missense "
                     "variants")
    sr.italic = True
    sr.font.size = Pt(12)
    sub.paragraph_format.space_after = Pt(14)

    meta = doc.add_paragraph()
    meta.alignment = WD_ALIGN_PARAGRAPH.CENTER
    mr = meta.add_run("Working manuscript · 8 August 2026\n"
                      "Analysis pipeline, data and figures: this repository\n"
                      "Canonical references: UniProt P35555 (fibrillin-1, 2,871 aa); "
                      "MANE Select transcript NM_000138.5")
    mr.font.size = Pt(9.5)
    mr.font.color.rgb = RGBColor.from_string("444444")
    meta.paragraph_format.space_after = Pt(18)

    # ================================================================================
    h("Abstract")
    rich(("Background. ", True, False),
         ("Marfan syndrome is a dominantly inherited disorder of connective tissue caused by "
          "variants in FBN1, the gene encoding fibrillin-1. Most of the protein consists of "
          "43 calcium-binding epidermal-growth-factor-like (cbEGF) domains, each of which grips "
          "a single calcium ion and is stapled together by three disulfide bonds. Two kinds of "
          "missense change are described as the hallmark disease mechanisms: substitutions that "
          "remove one of the six conserved cysteines, and substitutions that disturb the "
          "calcium-binding consensus. They are conventionally grouped together as damage to "
          "“structurally critical” residues.", False, False))

    rich(("Methods. ", True, False),
         (f"We assembled {f1['n_variants_total']:,} ClinVar missense variants mapped to "
          f"UniProt P35555 and verified residue-by-residue against the reference sequence, "
          f"annotated every position for domain membership, calcium-coordination role and "
          f"disulfide participation, and measured the biophysical consequence of each variant "
          f"that falls inside a calcium-bound structure "
          f"(n = {f2['n_with_ddg']}). Measurements comprised solvent accessibility, distance to "
          f"the bound ion, disulfide geometry and predicted folding free-energy change "
          f"(ΔΔG, FoldX), alongside AlphaMissense pathogenicity scores. Where no "
          f"experimental structure exists we did not model calcium silently: ion placement by "
          f"AlphaFold 3 cofolding was first calibrated against the four experimental "
          f"calcium-bound fragments and validated with CheckMyMetal.", False, False))

    rich(("Results. ", True, False),
         (f"The two mechanisms are not the same kind of damage. Cysteine-removing variants "
          f"destabilise the fold substantially (median ΔΔG "
          f"{f2['ddG_median_cys_removing']:.2f} kcal/mol versus "
          f"{f2['ddG_median_other_cbegf']:.2f} for other cbEGF residues; Cliff's "
          f"δ = {f2['stat_ddG_cys_vs_rest_delta']:+.2f}), and the cost is attributable "
          f"specifically to the lost disulfide term. Variants at residues that directly "
          f"coordinate calcium do not destabilise the fold at all "
          f"(median ΔΔG {f2['ddG_median_ca_ligand']:.2f} kcal/mol; if anything "
          f"marginally lower than non-ligand positions, δ = "
          f"{f2['stat_ddG_ligand_vs_nonligand_delta']:+.2f}), yet AlphaMissense scores them as "
          f"strongly pathogenic ({f2['am_median_ca_ligand']:.3f} versus "
          f"{f2['am_median_other_cbegf']:.3f}). "
          f"{f2['ca_ligand_damaging_but_stable']} of {f2['ca_ligand_with_am']} calcium-ligand "
          f"variants are predicted damaging while costing less than 2 kcal/mol to fold. The "
          f"dissociation survives matching for burial: among buried residues the "
          f"calcium-consensus class still costs roughly four-fold less folding energy than its "
          f"neighbours.", False, False))

    rich(("Conclusions. ", True, False),
         ("Calcium-site variants in fibrillin-1 appear to be pathogenic functionally rather "
          "than thermodynamically: they abolish the metal clamp while leaving the domain fold "
          "intact. The practical consequence is that any variant-effect predictor which ranks "
          "FBN1 missense changes by predicted destabilisation alone will systematically "
          "under-call one of the two principal disease mechanisms. We report this as the "
          "study's substantive contribution and treat the more familiar enrichment statistics "
          "as supporting context, because those statistics are partly circular by construction.",
          False, False))

    # ================================================================================
    h("Plain-language summary")
    callout("What this study is about, without the jargon", [
        "Fibrillin-1 is a long, thread-like protein that the body uses to build the elastic "
        "scaffolding in blood vessels, skin, ligaments and the fibres that hold the lens of the "
        "eye in place. When the gene for it is faulty, the result is Marfan syndrome, whose most "
        "dangerous feature is a weakened aorta.",

        "Most of the protein is made of 43 repeating units. Picture each unit as a small parcel "
        "tied with three pieces of string (chemical bonds called disulfides) and holding one "
        "calcium ion like a clasp. Two kinds of spelling mistake in the gene are known to cause "
        "disease: those that cut one of the strings, and those that break the clasp. Textbooks "
        "list them side by side, as though they were the same sort of harm.",

        "We measured both, computationally, across every well-documented fibrillin-1 spelling "
        "mistake we could map onto a three-dimensional structure. Cutting a string does what "
        "you would expect: the parcel becomes much less stable. Breaking the clasp does not "
        "destabilise the parcel at all — and yet an independent artificial-intelligence "
        "predictor, which knows nothing about our physics calculations, flags those same changes "
        "as clearly harmful.",

        "The lesson is practical. Clinical laboratories increasingly use software that ranks "
        "genetic variants by how much they are predicted to destabilise a protein. For this "
        "gene, such software would miss an entire class of genuinely harmful variants, because "
        "those variants do not destabilise anything — they simply stop the protein from "
        "holding on to calcium. Different kinds of damage need different kinds of test."])

    # ================================================================================
    h("1. Introduction")
    p("Fibrillin-1 is a 2,871-residue extracellular glycoprotein encoded by FBN1 on chromosome "
      "15. Secreted copies assemble into microfibrils, the load-bearing filaments of elastic "
      "and non-elastic connective tissue — the aorta, the skin, the ligaments, and the "
      "zonule that suspends the lens of the eye (Godwin et al., 2023; Asano et al., 2022). "
      "Pathogenic variants cause Marfan syndrome, an autosomal dominant disorder with a "
      "population prevalence of roughly 1 in 10,000 whose principal danger is aortic "
      "dissection.")

    p("The protein is strikingly modular. UniProt annotates 47 epidermal-growth-factor-like "
      "(EGF) domains, of which 43 are calcium-binding (cbEGF), together with transforming "
      "growth factor-β-binding (TB) and hybrid domains. Each cbEGF domain has two "
      "structural features that matter here. First, it binds a single calcium ion, coordinated "
      "by the side chains of an aspartate, a glutamate and an asparagine, plus three backbone "
      "carbonyl oxygens; the ion rigidifies the junction to the neighbouring domain and thereby "
      "sets the shape of the whole filament (Handford et al., 1995; Rao et al., 1995; Downing "
      "et al., 1996). Second, it contains six cysteines that form three disulfide bonds in an "
      "invariant 1–3, 2–4, 5–6 pattern.")

    p("Two missense mechanisms dominate the Marfan literature: substitution of one of the six "
      "conserved cysteines, which removes a disulfide, and substitution within the "
      "calcium-binding consensus, which is presumed to abolish ion binding. The ClinGen FBN1 "
      "Variant Curation Expert Panel formalised both as evidence for pathogenicity under the "
      "PM1 criterion (Drackley et al., 2024). What the literature does not generally do is ask "
      "whether the two mechanisms damage the protein in the same physical way. That is the "
      "question this study addresses.")

    rich(("We deliberately did not organise the study around the more common question of "
          "whether pathogenic variants are ", False, False),
         ("enriched", False, True),
         (" at these residues. That question is partly circular: the clinical classifications "
          "in ClinVar were themselves assigned using the ClinGen critical-residue rule, so "
          "measuring enrichment of ClinVar-pathogenic variants at ClinGen-defined critical "
          "residues partly measures the classification rule rather than the biology. Enrichment "
          "results are reported here, but as supporting context with that caveat attached in "
          "the same breath.", False, False))

    # ================================================================================
    h("2. Methods")

    h("2.1 Reference sequences and versions", 2)
    p("All protein positions in this work are positions in UniProt P35555 (entry version 264, "
      "2,871 residues). All cDNA descriptions are on the MANE Select transcript NM_000138.5. "
      "Pinning both, and recording the versions, matters more than it may appear: FBN1 has "
      "several transcript records in circulation, and a position that is off by one isoform "
      "silently corrupts every downstream measurement.")

    table(["Source", "Version pinned", "Records retrieved"],
          [["ClinVar bulk variant_summary", "release 2026-07-28", "18,770 FBN1 rows (9,327 GRCh38)"],
           ["UniProt P35555", "entry version 264", "2,871 aa; 56 domain features; 155 disulfides"],
           ["gnomAD", "v4 (gnomad_r4), GRCh38", "12,107 variants (2,783 missense)"],
           ["AlphaMissense", "Zenodo 8208688 v1.0.0", "54,568 substitution scores"],
           ["PDB / PDBe SIFTS", "queried 2026-08-02", "11 entries mapped to P35555"]],
          widths=[2.2, 1.9, 2.4])

    h("2.2 Variant normalisation and the wild-type identity check", 2)
    rich(("Every ClinVar record was converted to a (position, wild-type residue, mutant "
          "residue) triple on P35555 and then subjected to a check that is easy to skip and "
          "expensive to omit: does the residue the database claims match the residue actually "
          "present in the reference sequence at that position? ", False, False),
         (f"Across {f1['n_variants_total']:,} missense variants there were zero mismatches.",
          True, False),
         (" The result was corroborated three ways — an independent HGVS validator agreed "
          "on 100 of 100 randomly sampled variants; gnomAD's independently derived protein "
          "descriptions disagreed on none of 1,523 shared variants; and the row counts "
          "reconcile exactly against the ClinVar input, so nothing was dropped or duplicated.",
          False, False))

    p("Labelled sets were built using thresholds fixed in advance and recorded before the sets "
      "existed, so that they could not be tuned toward a result. The primary tier requires a "
      "ClinVar review status of at least two gold stars; a one-star tier is retained for "
      "sensitivity analysis. Benign status could be met either by ClinVar classification or by "
      "a gnomAD filtering allele frequency above 0.1%, which is incompatible with a fully "
      "penetrant dominant disease allele. Variants of uncertain significance and records with "
      "conflicting classifications were retained and labelled throughout, never silently "
      "discarded.")

    table(["Labelled set", "Primary tier (≥ 2 stars)", "Sensitivity tier (≥ 1 star)"],
          [["Pathogenic / likely pathogenic", f"{f1['n_pathogenic']}", "1,409"],
           ["Benign / likely benign", f"{f1['n_benign']}", "76"],
           ["Uncertain significance", f"{f1['n_vus']}", "2,295"],
           ["Conflicting classifications", "0", "492"]],
          widths=[2.6, 2.0, 1.9])

    p("As an external check on the labelling, AlphaMissense — which played no part in "
      "building the sets — independently classifies 97.2% of the pathogenic set as "
      "pathogenic and 88.2% of the benign set as benign.")

    h("2.3 Deriving the calcium consensus from data rather than from a paper", 2)
    p("Rather than copying consensus positions from the literature, we derived them. Each cbEGF "
      "domain was anchored on its six cysteines; candidate consensus positions were expressed "
      "as offsets from those anchors; and the offsets themselves were obtained by measuring "
      "which residues actually coordinate calcium in the deposited coordinates of the four "
      "experimental structures (eight domains, sixteen calcium sites, an identical ligand "
      "pattern at every one).")

    p("The derived rule was then checked against independent facts. It recovers all 24 "
      "side-chain calcium ligands measured in the experimental structures; it satisfies the "
      "ClinGen consensus motif in 43 of 43 domains with no exceptions; it predicts 258 cbEGF "
      "cysteines, exactly the figure the expert panel states; and all 129 cbEGF disulfides it "
      "identifies follow the 1–3, 2–4, 5–6 pattern, 43 of each, without "
      "exception.")

    rich(("A distinction is enforced throughout and is worth stating explicitly, because "
          "conflating the two is the commonest error in this literature: a ", False, False),
         ("calcium-coordinating residue", False, True),
         (" is not the same thing as a ", False, False),
         ("residue that happens to lie inside a cbEGF domain", False, True),
         (". This work uses two nested definitions and never mixes them within a single "
          "comparison. The annotation-level class (a consensus position of the calcium motif, "
          "215 positions genome-wide) is available for every variant. The structure-measured "
          "class (an oxygen or nitrogen within 3.2 Å of the modelled ion) is stricter, is "
          "available only where a structure exists, and is the class used for the biophysical "
          "comparisons.", False, False))

    h("2.4 Structures, and the problem of absent calcium", 2)
    rich(("AlphaFold models contain no metals. ", True, False),
         ("A raw AlphaFold model of a cbEGF domain has no calcium ion in it, so any tool asked "
          "to validate its metal site has nothing to score and will either fail or, worse, "
          "return confident-looking output about a site that does not exist. Calcium therefore "
          "had to be placed deliberately, and its placement treated as a modelling assumption "
          "with an error bar rather than as an observation.", False, False))

    p("A further complication emerged during setup: FBN1 has no AlphaFold Database entry at "
      "all. The database excludes sequences longer than 2,700 residues, and fibrillin-1 is "
      "2,871. This was diagnosed rather than assumed — a shorter control protein returns a "
      "model from the same endpoint, UniProt carries no AlphaFold cross-reference for P35555, "
      "and the exclusion rule is stated in the database's own publication (Varadi et al., "
      "2024). Because AlphaFill operates by transplanting ligands into AlphaFold Database "
      "models (Hekkelman et al., 2023), its absence also removed the intended cross-check, "
      "which had to be rebuilt from first principles.")

    p("Structures were therefore sourced in tiers. Where an experimental calcium-bound "
      "structure exists it was used directly: PDB 2W86 (cbEGF9–hybrid2–cbEGF10, X-ray "
      "1.8 Å), 1LMJ (cbEGF12–13, NMR), 1UZJ (cbEGF22–TB4–cbEGF23, X-ray) "
      "and 1EMN (cbEGF32–33, NMR). Elsewhere, AlphaFold 3 cofolding with calcium supplied "
      "the model. Crucially, the cofolding method was calibrated on the domains where the "
      "answer is already known before being trusted anywhere else, and a homology-transplant "
      "alternative was characterised in parallel by leave-one-out cross-validation.")

    p("All structure files were renumbered from PDB author numbering into P35555 numbering via "
      "SIFTS before any measurement was taken. This is not cosmetic: two of the four entries "
      "carry local author numbering offset by 804 and 1,066 residues respectively while the "
      "other two happen to coincide with UniProt, a mixture that would have silently "
      "mislabelled half of the results while being correct for the other half.")

    h("2.5 Biophysical measurement and statistics", 2)
    p("For each variant with structural coverage we computed relative solvent accessibility, "
      "secondary structure, distance to the nearest calcium ion, direct-ligand status, "
      "disulfide partner and sulfur–sulfur distance, and the predicted change in folding "
      "free energy (ΔΔG) using FoldX 5.1 (repair once per template, then build each "
      "mutant). Every wild-type residue was re-verified against P35555 immediately before the "
      "mutation string was issued. AlphaMissense scores were attached by matching on both "
      "position and wild-type identity, not position alone.")

    p("Group comparisons use the Mann–Whitney U test reported with Cliff's delta as an "
      "effect size and a bootstrap confidence interval; categorical contrasts use odds ratios "
      "with Woolf intervals; multiple testing is corrected within families by the "
      "Benjamini–Hochberg procedure. All random seeds are fixed, and re-running the "
      "analysis reproduces byte-identical output.")

    # ================================================================================
    h("3. Results")

    h("3.1 Where the variants are", 2)
    rich((f"Of {f1['n_variants_total']:,} verified missense variants, "
          f"{f1['n_pathogenic']} are pathogenic or likely pathogenic at the two-star tier, "
          f"{f1['n_vus']} are of uncertain significance, and {f1['n_benign']} are benign or "
          f"likely benign. ", False, False),
         (f"{f1['panelC_pathogenic_in_cbegf']} of the {f1['panelC_pathogenic_total']} "
          f"pathogenic variants fall inside a cbEGF domain", True, False),
         (f", and within those domains the distribution is heavily skewed toward cysteines: "
          f"{f1['panelC_total_cys']} remove a conserved cysteine, {f1['panelC_total_ca']} sit "
          f"at a calcium-consensus position, and {f1['panelC_total_other']} are elsewhere in "
          f"the domain.", False, False))

    p("The sparseness of the benign lane in Figure 1 is a real feature of the gene rather than "
      "a gap in the data collection, and it constrains what the study can conclude; this is "
      "discussed in Section 5.")

    figure("fig1_variant_landscape",
           "Figure 1. Pathogenic, uncertain and benign FBN1 missense variants along "
           "fibrillin-1. (A) Variants mapped to UniProt P35555, in three lanes by clinical "
           "class; stem height is the number of variants at that position; the shaded band "
           "marks the neonatal/severe region. (B) Domain architecture, drawn in lightness steps "
           "so it does not compete with the data palette. (C) Pathogenic variants per cbEGF "
           "domain, split by the site affected. Panel C uses the annotation-level calcium "
           "class, which is defined for all variants.")

    h("3.2 The central result: the two mechanisms are physically different", 2)
    rich(("Restricting attention to the ", False, False),
         (f"{f2['n_with_ddg']}", True, False),
         (" variants that lie inside a calcium-bound structure and therefore carry real "
          "measurements, the two hallmark mechanisms separate cleanly.", False, False))

    table(["Variant class", "n", "Median ΔΔG (kcal/mol)", "Median AlphaMissense",
           "Median RSA (%)"],
          [["Direct Ca²⁺ ligand", f2["ddG_n_ca_ligand"],
            f"{f2['ddG_median_ca_ligand']:.2f}", f"{f2['am_median_ca_ligand']:.3f}",
            f"{f3['rsa_median_ca_ligand']:.1f}"],
           ["cbEGF cysteine-removing", f2["ddG_n_cys_removing"],
            f"{f2['ddG_median_cys_removing']:.2f}", f"{f2['am_median_cys_removing']:.3f}",
            f"{f3['rsa_median_cys_removing']:.1f}"],
           ["Other cbEGF residue", f2["ddG_n_other_cbegf"],
            f"{f2['ddG_median_other_cbegf']:.2f}", f"{f2['am_median_other_cbegf']:.3f}",
            f"{f3['rsa_median_other_cbegf']:.1f}"]],
          widths=[1.9, 0.6, 1.5, 1.4, 1.1])

    r_cys, r_lig = (eff("ddG: cbEGF cysteine-removing vs all other missense"),
                    eff("ddG: direct Ca ligands vs non-ligands (structure-measured)"))
    r_am = eff("AlphaMissense: direct Ca ligands vs other cbEGF residues")
    r_dis = eff("FoldX disulfide energy term: Cys-removing vs other")

    rich(("Cysteine removal destabilises the fold. ", True, False),
         (f"Median ΔΔG is {f2['ddG_median_cys_removing']:.2f} kcal/mol against "
          f"{f2['ddG_median_other_cbegf']:.2f} for ordinary cbEGF residues (Cliff's "
          f"δ = {r_cys.effect:+.2f}, q = {r_cys.q_bh:.0e}). Decomposing the FoldX energy "
          f"shows the cost is not general strain but the specific loss of the disulfide: the "
          f"disulfide term is {f2['disulfide_median_cys_removing']:.2f} kcal/mol for these "
          f"variants and 0.00 for both other classes "
          f"(δ = {r_dis.effect:+.2f}, p = {r_dis.p_raw:.0e}).", False, False))

    rich(("Calcium-ligand loss does not destabilise the fold. ", True, False),
         (f"Median ΔΔG is {f2['ddG_median_ca_ligand']:.2f} kcal/mol, close to the "
          f"{f2['ddG_median_other_cbegf']:.2f} of ordinary cbEGF residues. The measured effect "
          f"against non-ligand positions is small and runs in the direction opposite to "
          f"destabilisation (δ = {r_lig.effect:+.2f}, q = {r_lig.q_bh:.1e}) — that "
          f"is, calcium ligands are if anything marginally ", False, False),
         ("cheaper", False, True),
         (f" to mutate, not equally costly. Yet AlphaMissense, which owes nothing to any force "
          f"field, scores the same variants at {f2['am_median_ca_ligand']:.3f} against "
          f"{f2['am_median_other_cbegf']:.3f} (δ = {r_am.effect:+.2f}, "
          f"q = {r_am.q_bh:.0e}). "
          f"{f2['ca_ligand_damaging_but_stable']} of {f2['ca_ligand_with_am']} calcium-ligand "
          f"variants are simultaneously predicted damaging and cost under 2 kcal/mol to fold.",
          False, False))

    p(f"The two predictors are not in general disagreement — across all "
      f"{f2['n_with_ddg']} measured variants they correlate at Spearman "
      f"ρ = {f2['stat_spearman_rho']:+.2f}. They disagree specifically, and reproducibly, "
      f"at the calcium sites.")

    figure("fig2_mechanism_dissociation",
           "Figure 2. Two Marfan mechanisms, two different kinds of damage. (A) Predicted "
           "folding cost. (B) Predicted pathogenicity from an orthogonal method. (C) The two "
           "together; the shaded quadrant contains variants that are predicted damaging yet "
           "cheap to fold. (D) The FoldX disulfide term alone, showing the cysteine cost is the "
           "lost staple rather than general strain. Boxes show median and interquartile range; "
           "counts above each axis limit are stated on the panel.")

    h("3.3 The obvious objection, and why it fails", 2)
    p("The natural counter-argument is that calcium ligands are simply surface residues, and "
      "mutations at the surface are always energetically cheap. Three lines of evidence rule "
      "this out.")

    r_bur = eff("ddG: Ca-consensus vs other cbEGF, buried (RSA<20%)")
    r_par = eff("ddG: Ca-consensus vs other cbEGF, partial (RSA 20-50%)")

    bullet(f"Calcium ligands are not especially exposed. Their median relative solvent "
           f"accessibility is {f3['rsa_median_ca_ligand']:.1f}%, between the buried cysteines "
           f"at {f3['rsa_median_cys_removing']:.1f}% and other cbEGF residues at "
           f"{f3['rsa_median_other_cbegf']:.1f}%.")
    bullet(f"The class is defined by measured geometry, not by position in the sequence. "
           f"Calcium ligands sit {f3['dist_ca_median_ca_ligand']:.1f} Å from the ion "
           f"— inside the coordination shell — against roughly "
           f"{f3['dist_ca_median_other_cbegf']:.0f} Å for other cbEGF residues.")
    bullet(f"Matching for burial does not remove the effect. Among buried residues "
           f"(accessibility below 20%) the calcium-consensus class costs "
           f"{f3['panelC_buried_Ca-consensus_median']:.2f} kcal/mol against "
           f"{f3['panelC_buried_other_cbEGF_median']:.2f} for its neighbours "
           f"(δ = {r_bur.effect:+.2f}, q = {r_bur.q_bh:.1e}); in the partially exposed "
           f"stratum the two are indistinguishable (δ = {r_par.effect:+.2f}, "
           f"q = {r_par.q_bh:.2f}).")

    p("A fourth observation is negative but informative: no calcium-consensus variant in the "
      "covered set has an accessibility above 50%. The class simply does not occur at fully "
      "exposed positions, which is what one would expect of residues whose job is to hold a "
      "metal against the protein body, and which is itself an argument against the exposure "
      "explanation.")

    figure("fig3_burial_control",
           "Figure 3. The calcium result is not a solvent-exposure artifact. (A) Burial by "
           "class. (B) Distance from each variant position to the nearest calcium ion in the "
           "calcium-placed structure. (C) Folding cost within burial strata; the exposed "
           "stratum is drawn with its calcium arm empty because no such variant exists in the "
           "covered set. Panel C follows the corresponding statistical test and therefore uses "
           "the annotation-level calcium class, unlike panels A and B.")

    h("3.4 What the structures show directly", 2)
    p(f"Figure 4 renders the highest-quality calcium-bound fibrillin-1 fragment available. The "
      f"cbEGF9 calcium site is coordinated by three side chains — Asp807, Glu810 and "
      f"Asn823, the last of which is also the β-hydroxylation site — together with "
      f"three backbone carbonyls, at 2.3–2.6 Å. This fragment alone carries "
      f"{f4['n_pathogenic_on_template']} pathogenic variants at "
      f"{f4['n_pathogenic_positions_on_template']} distinct positions.")

    figure("fig4_structure",
           "Figure 4. Calcium, disulfides and pathogenic variants in a solved FBN1 fragment "
           "(PDB 2W86, X-ray 1.8 Å), renumbered into P35555 numbering via SIFTS. (A) The "
           "fragment: two cbEGF domains flanking a hybrid domain. (B) The cbEGF9 calcium site "
           "with its ligands; the ion is drawn below its ionic radius so the ligands remain "
           "visible. (C) The AlphaFold 3 cofolded model superposed on the crystal structure. "
           "(D) Pathogenic positions, coloured by mechanism.")

    h("3.5 How far a modelled calcium ion can be trusted", 2)
    rich(("Because most cbEGF domains have never been solved, any extension of this analysis "
          "depends on placed rather than observed calcium. We therefore measured the placement "
          "error directly. ", False, False),
         (f"Against X-ray references, AlphaFold 3 cofolding puts the ion a median "
          f"{f5['ca_dev_median_af3_xray']:.2f} Å from its crystallographic position",
          True, False),
         (f"; against NMR references the figure is {f5['ca_dev_median_af3_nmr']:.2f} Å, "
          f"which tracks those references' own coordinate spread rather than exceeding it. "
          f"Homology transplant, evaluated leave-one-out across all donor-to-target pairs, is "
          f"{f5['transplant_over_af3_xray_ratio']:.1f} times worse "
          f"({f5['ca_dev_median_transplant']:.2f} Å). Cofolding is therefore the primary "
          f"method and transplant a cross-check, not the reverse.", False, False))

    p("One methodological trap is worth recording because it nearly produced the opposite "
      "conclusion. Superposing the whole two-domain construct rather than each domain "
      "separately inflated the apparent error roughly four-fold. These constructs are two rigid "
      "domains on a flexible hinge, so a global fit conflates metal placement with interdomain "
      "angle. Both numbers are retained in the output, because the gap between them quantifies "
      "the hinge motion rather than any modelling failure.")

    rich(("Sixteen calcium sites were submitted to CheckMyMetal (Zheng et al., 2017) — "
          "eight from AlphaFold 3 models and, deliberately, eight from the experimental "
          "structures as controls. The two groups are statistically indistinguishable "
          f"(gRMSD Mann–Whitney p = {f5['cmm_grmsd_p']:.2f}), which is the desired "
          "outcome. ", False, False),
         ("Submitting the controls also revealed something that would otherwise have been "
          "reported as a failure: the experimental structures score as “badly” as the "
          "models, well outside CheckMyMetal's published acceptable bands.", True, False),
         (" cbEGF calcium sites are irregular and dominated by backbone carbonyls, so the "
          "tool's absolute thresholds, calibrated on typical metalloprotein sites, do not "
          "transfer. No statement of the form “N% of sites are outliers” appears "
          "anywhere in this work, and none should.", False, False))

    figure("fig5_ca_placement",
           "Figure 5. Calcium placement is a modelling assumption; this is its error bar. "
           "(A) Displacement between a placed ion and the observed one, for cofolding against "
           "each reference type and for homology transplant. (B) Why the superposition frame "
           "matters. (C) CheckMyMetal parameters for all sixteen sites, experimental controls "
           "beside models.")

    h("3.6 Enrichment (reported as supporting context only)", 2)
    p("Against the correct null — the share of positions in each class, since every "
      "position offers the same nineteen possible substitutions — pathogenic variants are "
      f"enriched {f6['background_fold_cys']:.1f}-fold "
      f"[{f6['background_ci_cys'][0]:.1f}, {f6['background_ci_cys'][1]:.1f}] at "
      f"cysteine-removing changes, but only {f6['background_fold_ca']:.1f}-fold "
      f"[{f6['background_ci_ca'][0]:.1f}, {f6['background_ci_ca'][1]:.1f}] at "
      f"calcium-consensus residues.")

    rich(("That asymmetry is itself consistent with the central result. ", False, False),
         ("If calcium-site variants look thermodynamically benign, they are precisely the "
          "variants a clinical laboratory is least likely to have confidently classified as "
          "pathogenic", False, True),
         (" — so they are under-represented among the labels, and their apparent "
          "enrichment is correspondingly modest.", False, False))

    p("The pathogenic-versus-benign odds ratios in Figure 6 are much larger, and much less "
      "trustworthy, for the reason given in the Introduction: ClinVar's pathogenic calls used "
      "the same critical-residue rule the test is evaluating. The AlphaMissense contrast, which "
      "never used FBN1-specific rules, reproduces the direction independently at a far more "
      "modest magnitude, and is the version that should be believed.")

    figure("fig6_enrichment",
           "Figure 6. Enrichment at critical residues, supporting evidence only. (A) Against "
           "the background of every possible missense change. (B) Odds ratios on a log axis; "
           "the ClinVar-labelled contrasts are inflated by construction because those labels "
           "were assigned using the rule under test, while the AlphaMissense contrasts are "
           "independent of it.")

    # ================================================================================
    h("4. Discussion")
    p("The finding of this study is a dissociation rather than an association. Two mechanisms "
      "that the literature groups together — disulfide loss and calcium-site disruption "
      "— turn out to damage fibrillin-1 in physically different ways. Cysteine "
      "substitution is a thermodynamic lesion: the domain is measurably less stable, and the "
      "energy decomposition attributes the cost to the specific bond that was lost. Calcium-site "
      "substitution is not a thermodynamic lesion at all. The domain folds as well as it ever "
      "did; what it can no longer do is hold the ion that rigidifies its junction with the next "
      "domain.")

    p("This is biologically coherent. A cbEGF domain's contribution to a microfibril is not "
      "merely to exist but to sit at a defined angle to its neighbour, an angle the bound "
      "calcium sets. A variant that removes the clasp without unfolding the parcel leaves a "
      "structurally intact but functionally floppy protein, which is exactly the kind of defect "
      "that a stability calculation is blind to.")

    rich(("The practical consequence is a concrete methodological warning. ", False, False),
         ("A variant-effect predictor that ranks FBN1 missense changes by predicted "
          "destabilisation alone will systematically under-call one of the two principal "
          "disease mechanisms in this gene.", True, False),
         (" That prediction is testable, falls directly out of the measurements, and does not "
          "depend on any clinical label. It also suggests a specific remedy: for genes whose "
          "function depends on ligand or metal coordination, a ligand-contact term belongs "
          "alongside a stability term in variant interpretation, and the two should be reported "
          "separately rather than blended into a single score.", False, False))

    p("A secondary contribution is methodological. Cofolding an ion with a predicted structure "
      "is increasingly casual practice, and it is easy to forget that the resulting ion is a "
      "prediction rather than an observation. Calibrating the placement against structures "
      "where the answer is known — and submitting those known structures as controls to "
      "the validation tool — cost little and twice prevented a wrong conclusion: once when "
      "an inappropriate superposition frame inflated the apparent error, and once when a "
      "validation tool's thresholds turned out not to transfer to this domain family.")

    # ================================================================================
    h("5. Limitations, and problems with this study")
    p("The following are stated at length rather than compressed into a closing sentence, "
      "because several of them materially constrain what may be claimed. They fall into three "
      "groups: limits of the evidence, limits of the methods, and defects in the conduct of the "
      "work itself.")

    h("5.1 Limits of the evidence", 2)
    rich(("The enrichment analysis is partly circular. ", True, False),
         ("ClinVar's pathogenic classifications were made using the ClinGen PM1 "
          "critical-residue rule, which is the rule the enrichment tests. Any "
          "pathogenic-versus-benign contrast at critical residues therefore partly measures the "
          "classification process rather than the biology, and the resulting odds ratios are "
          "inflated. The AlphaMissense contrast is the orthogonal check and reproduces the "
          "direction at a far smaller magnitude. This is the principal reason the biophysical "
          "measurements, which are computed from coordinates and owe nothing to any "
          "classification scheme, are presented as the study's contribution.", False, False))

    rich(("The benign comparison group is very small, and cannot be enlarged. ", True, False),
         (f"There are {f1['n_benign']} benign or likely benign missense variants at the two-star "
          f"tier. Of 2,783 gnomAD missense variants, only 24 exceed the stand-alone benign "
          f"filtering-frequency threshold, and all 24 were already present in the ClinVar set "
          f"— there are zero additions available. FBN1 is under strong selective constraint "
          f"and simply does not carry many common missense variants. Every "
          f"pathogenic-versus-benign contrast in this work is therefore underpowered on one arm, "
          f"which is why those results are reported with confidence intervals rather than bare "
          f"p-values, and why they are not the headline.", False, False))

    rich(("Structural coverage is 17% and is not a random sample. ", True, False),
         (f"Only {f2['n_with_ddg']} of {f1['n_variants_total']:,} variants fall inside one of "
          f"the four calibrated constructs and therefore carry measured biophysics. The covered "
          f"set is the set of domains crystallographers chose to solve, which skews toward "
          f"well-behaved and historically interesting domains. Structural conclusions should "
          f"not be extrapolated to the whole protein without saying so.", False, False))

    rich(("Benign variants are ascertainment-biased. ", True, False),
         ("Benign classifications reach ClinVar largely because someone sequenced a patient and "
          "needed the variant adjudicated. They are not a random sample of tolerated variation "
          "and are enriched for positions that looked suspicious enough to be tested. This "
          "biases the benign set toward the same domains as the pathogenic set, which if "
          "anything makes the enrichment contrast conservative.", False, False))

    rich(("Clinical labels are of uneven quality. ", True, False),
         ("Even with a two-star floor, 5,437 of 9,327 records are single-submitter, and an "
          "independent audit found 31.6% of non-critical-residue FBN1 classifications "
          "overclassified once allele frequency was taken into account (Baudhuin et al., 2019). "
          "704 records carry conflicting classifications; these are retained and labelled, "
          "excluded from the primary contrast rather than resolved by fiat.", False, False))

    h("5.2 Limits of the methods", 2)
    rich(("Calcium in predicted models was placed, not observed. ", True, False),
         (f"Every ion in a predicted structure here was put there by cofolding. Calibration "
          f"gives the error bar — {f5['ca_dev_median_af3_xray']:.2f} Å against X-ray "
          f"references — and CheckMyMetal cannot distinguish the modelled sites from the "
          f"experimental ones. That validates the approach for this domain family. It does not "
          f"make a predicted calcium an observed one, and every structural metric derived from "
          f"a predicted site inherits this uncertainty.", False, False))

    rich(("ΔΔG is a prediction, not a measurement. ", True, False),
         ("FoldX is an empirical force field with a typical reported error near 0.5 kcal/mol, "
          "and values here were computed on a single repaired structure per template rather "
          "than on a conformational ensemble. It handles disulfide loss and packing well and "
          "long-range electrostatics and large rearrangements poorly. Values are used "
          "comparatively — which variants are worse than which — not as absolute "
          "thermodynamic quantities. No experimental stability measurement was performed to "
          "validate any of them.", False, False))

    rich(("CheckMyMetal's absolute thresholds do not apply here, and this was only discovered "
          "by accident of good design. ", True, False),
         ("Deposited, refined structures return values the tool's published bands call "
          "outliers. Only the experimental-versus-predicted comparison is meaningful. Had the "
          "controls not been submitted, the modelled sites would have been reported as failing "
          "validation.", False, False))

    rich(("Water molecules were removed before metal validation. ", True, False),
         ("Predicted models contain no water, so waters were stripped from the experimental "
          "structures too for a like-for-like comparison. Two experimental sites genuinely have "
          "a water ligand and were therefore scored with one fewer ligand than as deposited. The "
          "effect makes the control marginally worse rather than the models better, so it "
          "cannot inflate the calibration result, but an as-deposited control was not run.",
          False, False))

    rich(("Post-translational modifications are absent from every model. ", True, False),
         ("Fibrillin-1 is glycosylated, and the calcium consensus includes a "
          "β-hydroxylation site. Neither glycans nor β-hydroxylation are present in "
          "any structure used here. β-hydroxylation is not required for calcium binding, "
          "so the coordination analysis stands, but predicted stabilities omit any contribution "
          "from these modifications.", False, False))

    rich(("Results are tied to specific database releases. ", True, False),
         ("ClinVar is re-classified continuously; 29 FBN1 records already existed in the live "
          "index but not in the analysed snapshot. Re-running later will not reproduce identical "
          "counts. Every source is checksummed and version-pinned so that the discrepancy is "
          "visible rather than silent.", False, False))

    h("5.3 Defects in the conduct of the work", 2)
    p("These are recorded because each was capable of producing a confident and wrong result, "
      "and because a study that reports only its successes gives a reader no way to calibrate "
      "its failures.")

    numbered("A database query silently returned the wrong gene. A ClinVar search written as a "
             "gene-identifier field was rewritten by the service into an unrestricted free-text "
             "search, which returned records from unrelated genes as though they were FBN1 and "
             "inflated the apparent record count by roughly 40%. The environment check now "
             "fails outright if a query translation degrades in this way.")
    numbered("A population-frequency field was null for every variant. The filtering allele "
             "frequency field first used is empty throughout for this gene; had it not been "
             "checked, the benign set would have carried no frequency evidence at all while "
             "appearing to be a clean result.")
    numbered("A stratification bug hid the most destabilising variants. Burial strata were "
             "selected with an expression that treats a solvent accessibility of exactly zero "
             "as a missing value, because zero is falsy in the implementation language. Four "
             "completely buried variants — the most destabilising in the comparison, "
             "median ΔΔG 13.90 kcal/mol — fell out of every stratum. This was "
             "found only in the final phase, while attempting to reproduce an earlier number "
             "for a figure. The correction changed one of twenty test results and moved the "
             "effect away from zero, so the published direction was never at risk, but the bug "
             "had survived a passing validation gate.")
    numbered("A confidence interval was computed from the wrong tail. Fold-enrichment intervals "
             "were initially taken from a one-sided binomial test, which pins its upper bound "
             "and reports a meaningless ceiling. The test is legitimately one-sided; the "
             "interval must not be.")
    numbered("Two of the supplied literature PDFs were entirely the wrong papers. The file named "
             "for the expert-panel rulebook contained an unrelated paper on protein phase "
             "separation, and the file named for the microfibril cryo-electron microscopy study "
             "contained a paper on chromatin imaging. Both were detected by full-text search "
             "returning zero hits for the gene, and both were subsequently replaced and "
             "verified. The record of that cross-check in the project's literature notes was "
             "not updated after replacement and is therefore stale.")
    numbered("A citation in the project's own specification is incorrect. The expert-panel "
             "guideline is cited internally as Genome Medicine 16:141; the article actually "
             "present is 16:154. The digital object identifier is correct.")
    numbered("One primary source is still missing. The original publication describing the "
             "highest-quality structure used here is absent; the file under that name is a "
             "later review by the same group. Figure 4 rests on the deposited coordinates "
             "rather than on that paper, but the citation must be obtained before submission.")
    numbered("An internal count could not be reconciled with the published one. Our derivation "
             "of the expert panel's critical-residue set yields 655 positions against the "
             "panel's stated 633, a 3.4% discrepancy confined to auxiliary categories whose "
             "membership is defined in a supplementary table we do not have. The two primary "
             "classes match exactly. We stopped adjusting parameters rather than reverse-engineer "
             "agreement, since forcing a match would manufacture false precision.")
    numbered("A claim in an internal report overstates its own test. An earlier draft described "
             "calcium-ligand ΔΔG as “statistically indistinguishable” from "
             "other cbEGF residues. The corresponding test reports a small but significant "
             "effect in the direction opposite to destabilisation. The conclusion is unaffected "
             "and arguably strengthened, but the wording asserts more than the statistic "
             "supports and has been corrected in this manuscript.")
    numbered("Software citations are incomplete. The primary publications for the force field, "
             "secondary-structure assignment, solvent-accessibility, molecular-graphics and "
             "residue-mapping tools used here are not present in the project's literature "
             "folder. They are named with versions in the Methods and Software sections but are "
             "deliberately not listed as read references, and must be added before submission.")

    h("5.4 What this study is not", 2)
    p("This is a computational study end to end. No experiment was performed. Every stability "
      "value is a force-field prediction, every pathogenicity score is a machine-learning "
      "output, and most calcium positions outside four small fragments are predictions rather "
      "than observations. The dissociation reported here is a dissociation between two "
      "computational predictors that disagree in a specific and mechanistically interpretable "
      "way. It is a hypothesis worth testing at the bench — for instance by measuring "
      "calcium affinity and thermal stability for matched pairs of calcium-site and "
      "cysteine-site variants in the same domain — not a demonstration that has already "
      "been tested there.")

    # ================================================================================
    h("6. Conclusions")
    p("Pathogenic missense variants in fibrillin-1 damage the protein in at least two physically "
      "distinct ways. Cysteine substitutions destabilise the cbEGF fold, and the energetic cost "
      "is traceable to the specific disulfide that is lost. Calcium-site substitutions leave the "
      "fold intact and are predicted damaging for an entirely different reason: they abolish the "
      "metal coordination that sets the geometry of the microfibril. The two lesions have been "
      "described together for thirty years; they are not the same lesion.")
    p("For variant interpretation this implies that stability-based prediction, used alone, has "
      "a structured blind spot in this gene, and that a coordination-aware term should be "
      "reported alongside it rather than merged into it. For anyone modelling metalloproteins "
      "with predicted structures, it implies that a placed ion needs a calibration and a control "
      "before it is allowed to carry a conclusion.")

    # ================================================================================
    h("7. Data, code and reproducibility")
    p("The full pipeline runs in seven stages, each ending in a validation gate that must pass "
      "before the next begins; all seven pass on a clean re-run. Random seeds are fixed, source "
      "downloads are cached and checksummed, and re-running the statistical analysis reproduces "
      "byte-identical output. Each figure regenerates from a single committed script and carries "
      "a machine-readable record of every number it displays, which an automated check "
      "re-derives from the analysis tables; a figure that drifts from the data fails the build. "
      "Molecular graphics are produced headlessly from committed scripts with deterministic "
      "camera settings rather than by interactive posing.")
    p("Supplementary tables accompany this manuscript as a workbook containing the full variant "
      "table, the structurally measured subset, the metal-site validation results, every "
      "statistical test with effect size and corrected significance, the calcium-placement "
      "calibration, the transplant cross-validation, and a data dictionary defining every "
      "column.")

    h("Software", 2)
    p("Python 3.14; pandas 3.0.5; NumPy 2.5.1; SciPy 1.18.0; statsmodels 0.14.6; matplotlib "
      "3.11.1; Biopython 1.87; ProDy 2.6.1; freesasa 2.2.1; DSSP (mkdssp) 4.5.7; FoldX 5.1; "
      "PyMOL 3.1.0 (open-source); openpyxl 3.1.5. Residue mapping between UniProt and PDB "
      "numbering used SIFTS via the PDBe API; nomenclature validation used Mutalyzer.")

    # ================================================================================
    h("References")
    p("Every reference below corresponds to a full text held in the project library and read "
      "during the work. Software is cited by name and version in the section above rather than "
      "here, because those primary publications are not in the library — see Section 5.3, "
      "item 10.", italic=True, size=9.5)

    refs = [
        "Asano K, Cantalupo A, Sedes L, Ramirez F (2022). The multiple functions of fibrillin-1 "
        "microfibrils in organ development and disease. Int J Mol Sci 23:1892. "
        "DOI 10.3390/ijms23031892",

        "Baudhuin LM, Kluge ML, Kotzer KE, Lagerstedt SA (2019). Variability in gene-based "
        "knowledge impacts variant classification: an analysis of FBN1 missense variants in "
        "ClinVar. Eur J Hum Genet 27:1550–1560. DOI 10.1038/s41431-019-0440-3",

        "Chen S, Francioli LC, Goodrich JK, et al. (2024). A genomic mutational constraint map "
        "using variation in 76,156 human genomes. Nature 625:92–100. "
        "DOI 10.1038/s41586-023-06045-0",

        "Cheng J, Novati G, Pan J, et al. (2023). Accurate proteome-wide missense variant effect "
        "prediction with AlphaMissense. Science 381:eadg7492. DOI 10.1126/science.adg7492",

        "Delhomme C, et al. (2022). Genotype–mitral valve phenotype correlations in "
        "Marfan syndrome with FBN1 pathogenic variants. JACC Adv 1:100149. "
        "DOI 10.1016/j.jacadv.2022.100149",

        "Downing AK, Knott V, Werner JM, Cardy CM, Campbell ID, Handford PA (1996). Solution "
        "structure of a pair of calcium-binding epidermal growth factor-like domains: "
        "implications for the Marfan syndrome and other genetic disorders. Cell 85:597–605. "
        "DOI 10.1016/S0092-8674(00)81259-3",

        "Drackley A, Somerville C, Arnaud P, Baudhuin LM, Hanna N, Kluge ML, Kotzer K, "
        "Boileau C, et al. (2024). Interpretation and classification of FBN1 variants "
        "associated with Marfan syndrome: consensus recommendations from the Clinical Genome "
        "Resource's FBN1 variant curation expert panel. Genome Med 16:154. "
        "DOI 10.1186/s13073-024-01423-3",

        "Godwin ARF, Dajani R, Zhang X, Thomson J, et al. (2023). Fibrillin microfibril "
        "structure identifies long-range effects of inherited pathogenic mutations affecting a "
        "key regulatory latent TGFβ-binding site. Nat Struct Mol Biol 30:608–618. "
        "DOI 10.1038/s41594-023-00950-8",

        "Handford PA, Downing AK, Rao Z, Hewett DR, Sykes BC, Kielty CM (1995). The calcium "
        "binding properties and molecular organization of epidermal growth factor-like domains "
        "in human fibrillin-1. J Biol Chem 270:6751–6756. DOI 10.1074/jbc.270.12.6751",

        "Hekkelman ML, de Vries I, Joosten RP, Perrakis A (2023). AlphaFill: enriching AlphaFold "
        "models with ligands and cofactors. Nat Methods 20:205–213. "
        "DOI 10.1038/s41592-022-01685-y",

        "Jensen SA, Robertson IB, Handford PA (2012). Dissecting the fibrillin microfibril: "
        "structural insights into organization and function. Structure 20:215–225. "
        "DOI 10.1016/j.str.2011.12.008",

        "Jumper J, et al. (2021). Highly accurate protein structure prediction with AlphaFold. "
        "Nature 596:583–589. DOI 10.1038/s41586-021-03819-2",

        "Landrum MJ, et al. (2018). ClinVar: improving access to variant interpretations and "
        "supporting evidence. Nucleic Acids Res 46:D1062–D1067. DOI 10.1093/nar/gkx1153",

        "Li L, Huang J, Liu Y (2023). The extracellular matrix glycoprotein fibrillin-1 in "
        "health and disease. Front Cell Dev Biol 11:1302285. DOI 10.3389/fcell.2023.1302285",

        "Rao Z, Handford P, Mayhew M, Knott V, Brownlee GG, Stuart D (1995). The structure of a "
        "Ca²⁺-binding epidermal growth factor-like domain: its role in protein-protein "
        "interactions. Cell 82:131–141. DOI 10.1016/0092-8674(95)90059-4",

        "The UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. "
        "Nucleic Acids Res 51:D523–D531. DOI 10.1093/nar/gkac1052",

        "Varadi M, et al. (2024). AlphaFold Protein Structure Database in 2024: providing "
        "structure coverage for over 214 million protein sequences. Nucleic Acids Res "
        "52:D368–D375. DOI 10.1093/nar/gkad1011",

        "Zhang M, Chen Z, Chen T, Sun X, Jiang Y. Cysteine substitution and calcium-binding "
        "mutations in FBN1 cbEGF-like domains are associated with severe ocular involvement in "
        "patients with congenital ectopia lentis. Front Cell Dev Biol. "
        "DOI 10.3389/fcell.2021.816397",

        "Zheng H, Cooper DR, Porebski PJ, Shabalin IG, Handing KB, Minor W (2017). CheckMyMetal: "
        "a macromolecular metal-binding validation tool. Acta Crystallogr D 73:223–233. "
        "DOI 10.1107/S2059798317001061",
    ]
    for r in refs:
        par = doc.add_paragraph()
        par.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        par.paragraph_format.left_indent = Inches(0.3)
        par.paragraph_format.first_line_indent = Inches(-0.3)
        par.paragraph_format.space_after = Pt(6)
        run = par.add_run(r)
        run.font.size = Pt(9.5)

    doc.save(OUT)
    size_kb = OUT.stat().st_size / 1024
    log.append(f"wrote {OUT.relative_to(fs.ROOT)} ({size_kb:.0f} kB)")
    print(f"  wrote {OUT.relative_to(fs.ROOT)}  ({size_kb:.0f} kB)")
    print(f"  log -> {fs.write_log('20_manuscript', log).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
