
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
color grey85, exp and polymer
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
set sphere_scale, 0.35, elem Ca
show spheres, exp and resi 880+882+883+884 and name CA
color pal_outside, exp and resi 880+882+883+884 and name CA
show spheres, exp and resi 811+816+830+832+853+862+875+876+887+890+908+914+921+926+937 and name CA
color pal_cys, exp and resi 811+816+830+832+853+862+875+876+887+890+908+914+921+926+937 and name CA
show spheres, exp and resi 913 and name CA
color pal_ca, exp and resi 913 and name CA
set sphere_scale, 0.85, name CA
orient exp and polymer
turn x, -10
zoom exp and polymer, 1.0

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/v2_fig8_structure_c_variants.png, dpi=300
print("VIEW v2_fig8_structure_c_variants", cmd.get_view())
