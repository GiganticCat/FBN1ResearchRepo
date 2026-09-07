
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

load /home/harvey/research/fbn1-marfan/handoff/cmm/inbox/AF3_cbEGF9-hyb2-cbEGF10_model0.pdb, af3
hide everything
select exp9, exp and resi 807-846
select af39, af3 and resi 807-846
align af39 and name CA+C+N+O, exp9 and name CA+C+N+O
show cartoon, exp9
show cartoon, af39
color pal_exp, exp9
color pal_af3, af39
select caE, exp and elem Ca within 8 of exp9
select caA, af3 and elem Ca within 8 of af39
show spheres, caE
show spheres, caA
color pal_exp, caE
color pal_af3, caA
# Deliberately unequal radii, and the experimental ion is translucent. The two ions are
# 0.13 A apart, so equal opaque spheres would occlude one another and the panel would look
# like it contains a single ion. A small solid blue ball inside a larger translucent grey
# shell reads as "the model put it in the same place".
set sphere_scale, 0.45, caE
set sphere_scale, 0.18, caA
set sphere_transparency, 0.55, caE
print("ATOMS caE", cmd.count_atoms("caE"), "caA", cmd.count_atoms("caA"))
orient exp9
turn y, 20
zoom exp9, 1.5

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/fig4_C_af3_overlay.png, dpi=300
print("VIEW fig4_C_af3_overlay", cmd.get_view())
