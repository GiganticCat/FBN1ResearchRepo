
# ---- shared conventions for every Phase 6 structural panel ----------------------------
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.25
set cartoon_transparency, 0.0
set cartoon_loop_radius, 0.20
set cartoon_tube_radius, 0.28
set stick_radius, 0.11
set sphere_scale, 0.45, elem Ca
set dash_gap, 0.28
set dash_radius, 0.035
set dash_color, grey40
set label_size, 17
set label_color, black
set label_outline_color, white
set label_bg_color, white
set label_bg_transparency, 0.25
set label_position, (0, 0, 2.2)
set label_font_id, 7
set_color pal_ca,     [0.000, 0.000, 0.000]
set_color pal_cys,    [0.345, 0.345, 0.345]
set_color pal_other,  [0.886, 0.886, 0.886]
set_color pal_af3,    [0.761, 0.761, 0.761]
set_color pal_exp,    [0.239, 0.239, 0.239]
set_color pal_outside,[0.604, 0.604, 0.604]

load /home/harvey/research/fbn1-marfan/handoff/cmm/inbox/EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb, exp
remove exp and resi 805-806
remove exp and solvent
hide everything
select exp9,  exp and resi 807-846
select exphy, exp and resi 851-902
select exp10, exp and resi 910-951

select ca9, exp and elem Ca within 6 of (exp9 and name CA)
# only the site's own neighbourhood: the rest of the domain is a slab that crowds the labels
select site9, byres (exp9 within 9.5 of ca9)
show cartoon, site9
set cartoon_transparency, 0.55, site9
color grey80, site9
show sticks, exp and resi 807+810+823 and sidechain
color pal_ca, exp and resi 807+810+823 and elem O
color grey50, exp and resi 807+810+823 and elem C
show sticks, exp and resi 808+824+827 and name C+O
color grey60, exp and resi 808+824+827 and name C+O
show spheres, ca9
color pal_ca, ca9
set sphere_scale, 0.22, ca9
distance d_side, ca9, (exp and resi 807+810+823 and sidechain and elem O), 3.2
distance d_back, ca9, (exp and resi 808+824+827 and name O), 3.2
hide labels, d_side
hide labels, d_back
label exp and resi 807 and name CB, "Asp807"
label exp and resi 810 and name CB, "Glu810"
label exp and resi 823 and name CB, "Asn823"
label exp and resi 808 and name O, "Ile808"
label exp and resi 824 and name O, "Ser824"
label exp and resi 827 and name O, "Ser827"
# Asn823 and Ser827 sit almost on top of each other from this angle; nudge them apart rather
# than choosing a camera that hides one of the six ligands
set label_position, (-1.6, -1.9, 2.2), exp and resi 823
set label_position, ( 1.9,  1.2, 2.2), exp and resi 827
set label_size, 15
orient site9
turn y, 15
turn x, 8
center ca9
zoom ca9, 6.2

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/v2_fig8_structure_b_casite.png, dpi=300
print("VIEW v2_fig8_structure_b_casite", cmd.get_view())
