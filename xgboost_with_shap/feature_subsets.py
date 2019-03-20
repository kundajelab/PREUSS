# ----------------------------------
# region Base Subsets

sub_struct_overall = [
    "all_stem_length",
    "free_energy*",
    "sim_nor_score*",
    "probability_active_conf*"
]

# Mutation
sub_mut_overall = [
   "num_mutations*"
]

sub_mut_seq = [
   "mut_pos*",
   "mut_site_dist*",
   "mut_ref_nt*",
   "mut_nt*"
]

sub_mut_struct = [
   "mut_struct*",
   "mut_ref_struct*",
   "mut_prev_struct*",
   "mut_next_struct*",
   "mut_same_as_site*"
]

sub_mut_other = [
   "mut_exist*",
   "mut_type*"
]

# Editing Site
sub_site_seq = [
   "site_prev_nt*",
   "site_next_nt*"
]

sub_site_struct = [
   "site_struct*",
   "site_prev_struct*",
   "site_next_struct*",
   "site_1_1*",
   "site_length*",
   "site_length_stem*",
   "site_length_hairpin*",
   "site_length_bulge*",
   "site_length_internal_es*",
   "site_length_internal_ecs*",
   "site_5prm_cp_hairpin*",
   "site_5prm_cp_bulge*",
   "site_3prm_cp_bulge*",
   "site_5prm_cp_internal*",
   "site_3prm_cp_internal*"
]

# Upstream
sub_upstream = [
   "u_count*",
   "u_all_stem_length*",
   "u_hairpin_length*",

   "u1_exist*",
   "u1_distance*",
   "u1_struct*",
   "u1_length*",
   "u1_length_stem*",
   "u1_length_hairpin*",
   "u1_length_bulge*",
   "u1_length_internal_es*",
   "u1_length_internal_ecs*",
   "u1_5prm_cp_hairpin*",
   "u1_5prm_cp_bulge*",
   "u1_3prm_cp_bulge*",
   "u1_5prm_cp_internal*",
   "u1_3prm_cp_internal*",

   "u2_exist*",
   "u2_distance*",
   "u2_struct*",
   "u2_length*",
   "u2_length_stem*",
   "u2_length_hairpin*",
   "u2_length_bulge*",
   "u2_length_internal_es*",
   "u2_length_internal_ecs*",
   "u2_5prm_cp_hairpin*",
   "u2_5prm_cp_bulge*",
   "u2_3prm_cp_bulge*",
   "u2_5prm_cp_internal*",
   "u2_3prm_cp_internal*",

   "u3_exist*",
   "u3_distance*",
   "u3_struct*",
   "u3_length*",
   "u3_length_stem*",
   "u3_length_hairpin*",
   "u3_length_bulge*",
   "u3_length_internal_es*",
   "u3_length_internal_ecs*",
   "u3_5prm_cp_hairpin*",
   "u3_5prm_cp_bulge*",
   "u3_3prm_cp_bulge*",
   "u3_5prm_cp_internal*",
   "u3_3prm_cp_internal*"
]

# Downstream
sub_downstream = [
   "d_count*",
   "d_all_stem_length*",

   "d1_exist*",
   "d1_distance*",
   "d1_struct*",
   "d1_length*",
   "d1_length_stem*",
   "d1_length_hairpin*",
   "d1_length_bulge*",
   "d1_length_internal_es*",
   "d1_length_internal_ecs*",
   "d1_5prm_cp_hairpin*",
   "d1_5prm_cp_bulge*",
   "d1_3prm_cp_bulge*",
   "d1_5prm_cp_internal*",
   "d1_3prm_cp_internal*",

   "d2_exist*",
   "d2_distance*",
   "d2_struct*",
   "d2_length*",
   "d2_length_stem*",
   "d2_length_hairpin*",
   "d2_length_bulge*",
   "d2_length_internal_es*",
   "d2_length_internal_ecs*",
   "d2_5prm_cp_hairpin*",
   "d2_5prm_cp_bulge*",
   "d2_3prm_cp_bulge*",
   "d2_5prm_cp_internal*",
   "d2_3prm_cp_internal*",

   "d3_exist*",
   "d3_distance*",
   "d3_struct*",
   "d3_length*",
   "d3_length_stem*",
   "d3_length_hairpin*",
   "d3_length_bulge*",
   "d3_length_internal_es*",
   "d3_length_internal_ecs*",
   "d3_5prm_cp_hairpin*",
   "d3_5prm_cp_bulge*",
   "d3_3prm_cp_bulge*",
   "d3_5prm_cp_internal*",
   "d3_3prm_cp_internal*"
]

# endregion

# Overall Structure
overall_structure = sub_struct_overall

# "Mutation" Related
mut = sub_mut_overall + sub_mut_seq + sub_mut_struct
mut_seq = sub_mut_seq
mut_struct = sub_mut_struct

# "Editing Site" Related
site = sub_site_seq + sub_site_struct
site_struct = sub_site_struct
site_seq = sub_site_seq
      
# Upstream + Downstream
u = sub_upstream
d = sub_downstream
u_d = sub_upstream + sub_downstream

# Mutations + Editing Site
mut_site = mut + site
# "Sequence" related of "Mutations" + "Editing Site"
mut_seq_site = sub_mut_seq + sub_site_seq

# NO "Upstream"
no_up = sub_struct_overall + mut + site + d
# NO "Downstream"
no_down = sub_struct_overall + mut + site + u
# NO "Mutation"
no_mut = sub_struct_overall + site + u + d

