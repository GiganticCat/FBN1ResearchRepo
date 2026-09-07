
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

select ca9, exp and elem Ca within 8 of (exp and resi 807-846)
show cartoon, exp and resi 807-846
color grey85, exp and resi 807-846
set cartoon_transparency, 0.6, exp
show spheres, ca9
color pal_ca, ca9
set sphere_scale, 0.28, ca9
show sticks, exp and resi 807+810+823 and (sidechain or name CA)
color pal_ca, exp and resi 807+810+823 and elem C
show sticks, exp and resi 808+824+827 and name C+O+CA
color grey45, exp and resi 808+824+827 and elem C
distance d807, (exp and resi 807 and sidechain and (elem O+N)), (ca9), 3.2, 0
distance d810, (exp and resi 810 and sidechain and (elem O+N)), (ca9), 3.2, 0
distance d823, (exp and resi 823 and sidechain and (elem O+N)), (ca9), 3.2, 0
distance b808, (exp and resi 808 and name O), (ca9), 3.4, 0
distance b824, (exp and resi 824 and name O), (ca9), 3.4, 0
distance b827, (exp and resi 827 and name O), (ca9), 3.4, 0
hide labels, d* b*
label exp and resi 807 and name CB, "Asp807"
label exp and resi 810 and name CB, "Glu810"
label exp and resi 823 and name CB, "Asn823"
orient (exp and resi 807+810+823+808+824+827) or ca9
turn y, 20
zoom ((exp and resi 807+810+823+808+824+827) or ca9), 2.2

ray 1500, 1350
png /home/harvey/research/fbn1-marfan/figures/panels/fig4_B_casite.png, dpi=300
print("VIEW fig4_B_casite", cmd.get_view())
