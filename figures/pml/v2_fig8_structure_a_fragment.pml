
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

show cartoon, exp and polymer
# the inter-domain linkers belong to no selection below; colour everything first so they can
# never come out in PyMOL's default green
color grey70,    exp and polymer
color pal_af3,   exp9
color pal_other, exphy
color pal_af3,   exp10
show sticks, exp and resi 811+821+816+830+832+845+914+926+921+935+937+950 and not name N+C+O
set stick_radius, 0.17, exp and resi 811+821+816+830+832+845+914+926+921+935+937+950
color pal_cys, exp and resi 811+821+816+830+832+845+914+926+921+935+937+950 and elem S
color pal_cys, exp and resi 811+821+816+830+832+845+914+926+921+935+937+950 and elem C
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
orient exp and polymer
turn x, -10
zoom exp and polymer, 1.0

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/v2_fig8_structure_a_fragment.png, dpi=300
print("VIEW v2_fig8_structure_a_fragment", cmd.get_view())
