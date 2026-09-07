
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
set_color pal_ca,     [0.106, 0.686, 0.478]
set_color pal_cys,    [0.290, 0.227, 0.655]
set_color pal_other,  [0.478, 0.478, 0.455]
set_color pal_af3,    [0.165, 0.471, 0.839]
set_color pal_exp,    [0.322, 0.318, 0.306]
set_color pal_outside,[0.400, 0.400, 0.400]

load /home/harvey/research/fbn1-marfan/handoff/cmm/inbox/EXP_2W86_cbEGF9-hyb2-cbEGF10.pdb, exp
remove exp and resi 805-806
hide everything

show cartoon, exp and polymer
color grey85, exp and polymer
show spheres, exp and elem Ca
color pal_ca, exp and elem Ca
set sphere_scale, 0.35, elem Ca
show spheres, exp and resi 811+816+830+832+914+921+926+937 and name CA
color pal_cys, exp and resi 811+816+830+832+914+921+926+937 and name CA
show spheres, exp and resi 913 and name CA
color pal_ca,  exp and resi 913 and name CA
show spheres, exp and resi 853+862+875+876+880+882+883+884+887+890+908 and name CA
color pal_outside, exp and resi 853+862+875+876+880+882+883+884+887+890+908 and name CA
set sphere_scale, 0.85, name CA
orient exp and polymer
turn x, -12
zoom exp and polymer, 1.0

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/fig4_D_variants.png, dpi=300
print("VIEW fig4_D_variants", cmd.get_view())
