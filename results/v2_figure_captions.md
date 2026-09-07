# Figure notes — SUPERSEDED 2026-08-30

> The figure set described below was rebuilt from scratch on 2026-08-30 and renamed. The
> current set is `fig1_module`, `fig2_stability`, `fig3_coordination`, `fig4_interpretation`,
> `figS1_validation` and `figS2_structure`, built by `scripts/50`–`55` on
> `scripts/paperstyle.py`. Their legends are in `results/manuscript_final.md`; the reasoning
> behind the style is in the `paperstyle.py` docstring and in `results/SESSION_CONTEXT.md`.
> Everything below refers to the previous `v2_fig*` files and is kept only for its caveats.

---

# v2 figures — what each one shows, in plain language

> **FINAL FIGURE SET, 2026-08-30.** The manuscript's own legends are in
> `results/manuscript_final.md`; this file is the plain-language companion. The set is now six
> figures plus one supplementary. **Figure 4 changed**: the nine-variant geometry pilot
> (`v2_fig2_geometry`) was replaced by the 57-variant rate figure (`v2_fig7_rate`) once batch 4
> came back, and the pilot became Supplementary Figure S1.
>
> | # | File | Script | Carries |
> |---|---|---|---|
> | 1 | `v2_fig3_landscape` | `34_v2fig3_landscape.py` | coverage 751 → 3,831 |
> | 2 | `v2_fig1_stability` | `32_v2fig1_stability.py` | the fold survives (two force fields) |
> | 3 | `v2_fig4_burial` | `35_v2fig4_burial.py` | burial does not explain it |
> | 4 | `v2_fig7_rate` | `42_v2fig7_rate.py` | **how often the ion loses grip — 45%** |
> | 5 | `v2_fig5_validation` | `36_v2fig5_validation.py` | AF3 placement + BVS vs CheckMyMetal |
> | 6 | `v2_fig6_enrichment` | `37_v2fig6_enrichment.py` | enrichment, circularity attached |
> | S1 | `v2_fig2_geometry` | `33_v2fig2_geometry.py` | the nine-variant pilot, incl. N2144S |

---

## Figure 4 — `v2_fig7_rate.png` (the new one)
**Script:** `scripts/42_v2fig7_rate.py`

### What you're looking at
Sixty variants were folded with their calcium ions and compared to the normal version of the same
protein piece. Each bar in **panel a** is one variant: how much grip its calcium site lost. Green
bars are the pathogenic variants that touch the calcium; grey and blue are the two control groups.

The grey band is the honest definition of "nothing happened". Every construct has three calcium
sites and a variant can only reach one, so the other two are measurements of no effect — 109 of
them, plus 14 negative controls. The band is where 95% of those land. **A bar has to stick out
below the band to count.**

- **Panel b** — the rate: how many of each group broke their site.
- **Panel c** — the reason the rate isn't 100%. Substitutions that swap one oxygen-donating amino
  acid for another (D→N, D→E) mostly keep the calcium; ones that remove the oxygen mostly don't.
- **Panel d** — the same variants plotted against how much they destabilise the protein. The cloud
  is shapeless on purpose: the two kinds of damage are unrelated.

### The finding it supports
**Nearly half of pathogenic calcium-touching variants (15 of 33) measurably loosen the calcium
site, against 1 of 24 controls.** The backbone around the site moves by 0.29 Å — nothing. And
neither of the tools people actually use for this can tell which variants broke the site: the
stability calculation is uncorrelated with the damage (ρ = −0.09), and the sequence-based AI gives
the same score, 0.982, to both groups.

### What it does not show — read this before quoting it
- **45% is a point estimate on 33 variants.** The interval is 30–62%. Do not quote it as "about
  half" without the interval somewhere nearby.
- **The other 18 variants are unexplained.** They are pathogenic in ClinVar and they keep their
  calcium. Section 15 of `LIMITATIONS.md` lists the three things that could be — β-hydroxylation,
  an effect too small to model, or overclassification — and this analysis cannot separate them.
- **The ion is modelled, not observed.** The measure is a geometric proxy for grip, not a binding
  constant. It was checked against the one variant with laboratory data (N2144) and agrees, once.
- **Panel a's tallest bar is not the most important variant**, it is the most extreme measurement.
  N2144H is marked because its sister substitution N2144S has bench data, not because 38% is a
  meaningful threshold.

---

> **The three notes below predate the final set — see `results/v2_report.html` and the block above.**
> On 2026-08-16 the figure set was rebuilt to journal specification (7 pt type, lowercase panel
> letters, no panel titles, hairline furniture, 600 dpi, sized to 180 mm double-column) and
> expanded from three figures to six. Figures now carry standard nomenclature rather than
> plain-language axis labels — the plain reading moved into the report, where there is room to
> define terms. This file is retained for the three original figures' caveats, which remain
> accurate, but the report is the current document.
>
> **The six, in manuscript order:**
>
> | # | File | Script | Carries |
> |---|---|---|---|
> | 1 | `v2_fig3_landscape` | `34_v2fig3_landscape.py` | coverage 751 → 3,831 |
> | 2 | `v2_fig1_stability` | `32_v2fig1_stability.py` | the fold survives (two force fields) |
> | 3 | `v2_fig4_burial` | `35_v2fig4_burial.py` | burial does not explain it |
> | 4 | `v2_fig2_geometry` | `33_v2fig2_geometry.py` | the ion loses grip |
> | 5 | `v2_fig5_validation` | `36_v2fig5_validation.py` | AF3 placement + BVS vs CheckMyMetal |
> | 6 | `v2_fig6_enrichment` | `37_v2fig6_enrichment.py` | enrichment, circularity attached |

Written 2026-08-16. These figures replace `fig1`–`fig6`, which were built on the earlier
analysis (751 variants, calcium sites guessed from the sequence pattern). The old files are
still in `figures/` and are **superseded**, not deleted.

Every figure regenerates from one script and writes a `.numbers.json` listing every value it
prints on the page, so any number here can be traced back to the data table it came from.

---

## The one-sentence story these figures tell together

Marfan syndrome variants in the calcium-binding domains of fibrillin-1 break the protein in
**two different ways**, and the two look nothing alike. Variants that remove a cysteine bridge
physically destabilise the protein's shape. Variants that touch the calcium **leave the shape
intact** and instead stop the site from holding onto its calcium ion. A tool that only measures
"does the protein still fold" will call the second group harmless — and be wrong.

---

## Figure 1 — `v2_fig1_stability.png`
**Script:** `scripts/32_v2fig1_stability.py`

### What you're looking at
Four groups of variants, compared four different ways. The groups are:

| Label on the figure | What it actually means |
|---|---|
| **Touches the calcium** | The changed amino acid is one that physically grips the calcium ion |
| **Same letter, no calcium** | The same *kinds* of amino acid (D, E, N), but not gripping any calcium — the fair comparison group |
| **Looks right, isn't** | Positions the textbook sequence pattern flags as calcium sites, but which grip no calcium in any model |
| **Breaks a cysteine bridge** | Removes one of the six cysteines that staple the domain together |

**Panels A and B** both ask "how much does this change destabilise the protein's shape?" — using
two completely independent pieces of software. Higher = more destabilising. **Panel C** boils
that down to a single number per group, where **0 means "no different from an ordinary look-alike
position."** **Panel D** asks a different question entirely: a sequence-only AI predictor's
opinion on whether each variant is damaging.

### The finding it supports
**Calcium-touching variants do not destabilise the protein, but they are still damaging.**

- Breaking a cysteine bridge is expensive in both programs (FoldX median 4.32 kcal/mol vs 0.24
  for controls; Rosetta 24.76 vs 2.11).
- Touching the calcium is **not**: Rosetta puts it at 2.33 against a control of 2.11 — effect
  size +0.01, with a confidence interval straddling zero, i.e. genuinely no difference.
- Yet in Panel D the sequence-only predictor scores those same calcium variants at 0.989 out of
  1 — as damaging as anything in the study.

That gap between "doesn't destabilise" and "clearly damaging" is the whole point of the paper,
and it's what Figure 2 goes on to explain.

### What it does not show — read this before quoting it
- **The two programs disagree slightly on calcium.** FoldX gives +0.23 with an interval that
  just misses zero; Rosetta gives +0.01. They are computed on different variant sets (50 vs 271),
  because FoldX only covers the four experimentally solved fragments. The honest reading is that
  the effect is small in both and absent in the better-powered, metal-aware one — not that the
  two agree exactly.
- **The Rosetta cysteine number is inflated.** 451 of those points sit off the top of the panel.
  They're mostly large amino acids (tryptophan, tyrosine, arginine) replacing a small cysteine in
  a tight space, where the calculation can't relax the backbone to accommodate them, so it reports
  an impossible collision rather than a real energy. The *direction* is solid; **do not quote
  24.76 as a physical quantity.**

---

## Figure 2 — `v2_fig2_geometry.png`
**Script:** `scripts/33_v2fig2_geometry.py`

### What you're looking at
Nine variants, each folded by AlphaFold3 **with its calcium ion present**, then compared to the
normal version of the same protein piece. Each variant has two calcium sites, drawn as two dots.

- **Panel A** — how much grip the ion loses. Left = the ion is being held less tightly.
- **Panel B** — the same change, divided by how much the answer wobbles between repeat runs of
  the same calculation. Every protein piece was folded five times; the grey band is that natural
  wobble. **Only dots outside the grey band mean anything.** Filled dots are outside it; hollow
  dots are inside and should be read as "no detectable change."
- **Panel C** — how many contacts still hold the ion in place, normal → mutant.

The bottom two rows are **controls**: P1141L is a harmless variant, and C1138F breaks a cysteine
bridge without touching any calcium. Neither should move the calcium site, and that's what sets
the bar for believing the others.

### The finding it supports
**This is the damage that Figure 1 couldn't see.** The most important row is **N2144S** — the one
variant in this study whose real-world behaviour has been measured in a laboratory (it binds
calcium 5–9× more weakly). The model reproduces that independently: it loses **28.5%** of its grip
on the ion, more than four times the run-to-run wobble. That agreement is what licenses trusting
the same measurement on variants nobody has tested in a lab.

### What it does not show — read this before quoting it
**The effect is real but not universal.** Only **6 of 18** site measurements exceed the noise
band, 5 of them at calcium-touching variants. Several calcium variants (E1073K, N2183S) do not
clear the bar at either site, and E913K moves in *opposite* directions at its two sites. With nine
variants this is a demonstration that the effect exists and can be measured — **not** an estimate
of how often it happens. The figure's title says this out loud on purpose.

---

## Figure 3 — `v2_fig3_landscape.png`
**Script:** `scripts/34_v2fig3_landscape.py`

### What you're looking at
- **Panel A** — every classified variant as a tick mark, placed where it sits along the protein's
  2,871 amino acids. Three rows: causes Marfan, unclear, harmless. The grey bar underneath is the
  protein's architecture; dark blocks are the 43 calcium-binding domains.
- **Panel B** — how many variants we can measure structurally, before and after folding all 43
  domains.
- **Panel C** — how the measurable variants divide into the four groups used in Figure 1.

### The finding it supports
**This is a methods result, and it's the reason the rest of the analysis is trustworthy.**
Before, only 751 variants (17%) could be measured — and not a random 17%: they were whichever
variants happened to fall inside the four protein fragments that crystallographers chose to solve
decades ago. Folding all 43 domains raised it to **3,831 (86%)** and removed that bias.

Panel A also shows the honest scale of the harmless group: **34 variants**, against 578 that cause
Marfan. That near-empty row is drawn rather than hidden, because it is the reason several
comparisons in this project are reported with confidence intervals instead of p-values.

### What it does not show
Panel A's title reports that 74% of Marfan-causing variants fall in a calcium-binding domain. The
matching figure for harmless variants is 62% — a much smaller gap than it might appear, and one
this panel is not powered to establish with 34 harmless variants. It is context, not a result.

---

## Resolved (was: still to do on the figures)

1. ~~An enrichment figure~~ — **done**, `v2_fig6_enrichment`, built on the 1,901-variant gnomAD
   population set with the circularity caveat printed inside the figure.
2. **The structural render** `figures/fig4_structure.png` (the PyMOL picture of a real solved
   fragment with its calcium ions) has not been rebuilt. It is drawn from deposited experimental
   coordinates rather than from the analysis tables, so it is still correct — but its caption
   still uses the old sequence-pattern definition of a calcium site and needs rewording.
3. ~~The Rosetta cysteine collision problem~~ — **fixed** by restrained backbone relaxation inside
   the repack shell; 21.5% → 1.0% above 100 REU, cysteine median 24.76 → 10.81 REU. Figure 1's
   caveat above is stale in its numbers; the direction it describes still holds.

4. **Still open:** the PyMOL structural render (item 2 above) is the one figure the final set does
   not include. It would sit before Figure 2 as an orientation panel.
