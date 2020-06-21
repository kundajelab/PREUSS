#!/usr/bin/env Rscript

## =========================================
# R Script - Analysis Script
#
# Date: 2018-12-19
#
# Example:
#   Rscript  gen_notebook_matrix.R --pwd /d1/d2/d3/
#

library(rmarkdown);
library(optparse)

# Load `pandoc`
# Call `Sys.getenv("RSTUDIO_PANDOC")` in "Rstudio" to get the actual path of `pandoc`
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

## ---- Prep - Decide the "folder" of THIS running script -------------------

initial.options <- commandArgs(trailingOnly = FALSE) # Set `trailingOnly` to `FALSE` to include "basic" info.
# Get the "input" of the R script, by searching for `--file=`.
script.parameter <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
# Be sure to "normailize" the path just in case it contains "relative" path
script.name = normalizePath(file.path(getwd(), script.parameter))

# If you are running inside "RStudio"
# library(rstudioapi)
# script.name = rstudioapi::getActiveDocumentContext()$path

script.dir <- dirname(script.name)

## ---- Prep - Parse the Parameters -------------------

option_list = list(
    make_option(c("-p", "--pwd"), type="character", default=NULL,
              help="Path to working folder.", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## ---- Generate Reports -------------------

render_report = function(title, working_folder, dataset_codename_str, split,
                            external_test_dataset_codename_str = '', feature_column_str = '',
                            model_file = '') {
    #
    rmarkdown::render(
        "./notebooks/ml-practice-r/rf_tunning.rmd",
        params = list(
            title = title,
            working_folder = working_folder,
            dataset_codename_str = dataset_codename_str,
            training_validation_split = split,
            external_test_dataset_codename_str = external_test_dataset_codename_str,
            feature_column_str = feature_column_str,
            model_file = model_file
        ),
        output_file = normalizePath(file.path(working_folder, paste0("_output/reports/rf_matrix.", tolower(title), ".html")))
    )
}

# Defaults
training_validation_split = 0.8
working_folder = '/Users/cowfox/Desktop/_rna_lib_analysis/__analysis/2018_12_10-ML_Random_Forest-Tuning'
if (!is.null(opt$pwd)){
  # print_help(opt_parser)
  # stop("At least one argument must be supplied (input file).n", call.=FALSE
    working_folder = opt$pwd
}


# Feature Strings
# ALL
feature_all_str = 'site_prev_nt, site_next_nt, num_mutations, mut_pos, mut_site_dist, mut_ref_nt, mut_nt, mut_prev_struct, mut_same_as_site, site_struct, site_length_interior_es, site_length_interior_ecs, site_3prm_cp_interior, d_all_stem_length, d2_struct, d1_distance, d1_length_stem, d2_distance, d3_distance, d3_length_stem, d2_length_interior_ecs, d2_length_interior_es, d_count, u1_length_stem, u3_distance, u_all_stem_length, u2_distance, u1_distance, u2_length_stem, u2_struct, u_count, d2_3prm_cp_interior'

# "Mut-excluded" Features
feature_no_mut_str = 'site_prev_nt, site_next_nt, site_struct, site_length_interior_es, site_length_interior_ecs, site_3prm_cp_interior, d_all_stem_length, d2_struct, d1_distance, d1_length_stem, d2_distance, d3_distance, d3_length_stem, d2_length_interior_ecs, d2_length_interior_es, d_count, u1_length_stem, u3_distance, u_all_stem_length, u2_distance, u1_distance, u2_length_stem, u2_struct, u_count, d2_3prm_cp_interior'

# "Site-exclude" features
feature_no_site_str = 'site_prev_nt, site_next_nt, num_mutations, mut_pos, mut_site_dist, mut_ref_nt, mut_nt, mut_prev_struct, mut_same_as_site,  d_all_stem_length, d2_struct, d1_distance, d1_length_stem, d2_distance, d3_distance, d3_length_stem, d2_length_interior_ecs, d2_length_interior_es, d_count, u1_length_stem, u3_distance, u_all_stem_length, u2_distance, u1_distance, u2_length_stem, u2_struct, u_count, d2_3prm_cp_interior'

# "Downstream-exclude" features
feature_no_downstream_str = 'site_prev_nt, site_next_nt, num_mutations, mut_pos, mut_site_dist, mut_ref_nt, mut_nt, mut_prev_struct, mut_same_as_site, site_struct, site_length_interior_es, site_length_interior_ecs, site_3prm_cp_interior, u1_length_stem, u3_distance, u_all_stem_length, u2_distance, u1_distance, u2_length_stem, u2_struct, u_count'

# "Upstream-exclude" features
feature_no_upstream_str = 'site_prev_nt, site_next_nt, num_mutations, mut_pos, mut_site_dist, mut_ref_nt, mut_nt, mut_prev_struct, mut_same_as_site, site_struct, site_length_interior_es, site_length_interior_ecs, site_3prm_cp_interior, d_all_stem_length, d2_struct, d1_distance, d1_length_stem, d2_distance, d3_distance, d3_length_stem, d2_length_interior_ecs, d2_length_interior_es, d_count, d2_3prm_cp_interior'


# Neil1
render_report("Neil1_All_2", working_folder,
dataset_codename_str = 'neil1_computational',
split = training_validation_split,
feature_column_str = feature_all_str)
#
# render_report("Neil1_No_Mut", working_folder,
# 'neil1_computational',
# training_validation_split,
# feature_column_str = feature_no_mut_str)
#
# render_report("Neil1_No_Site", working_folder,
# 'neil1_computational',
# training_validation_split,
# feature_column_str = feature_no_site_str)
#
# render_report("Neil1_No_Downstream", working_folder,
# 'neil1_computational',
# training_validation_split,
# feature_column_str = feature_no_downstream_str)
#
# render_report("Neil1_No_Upstream", working_folder,
# 'neil1_computational',
# training_validation_split,
# feature_column_str = feature_no_upstream_str)


# AJUBA
# render_report("AJUBA_BC_All", working_folder,
# 'ajuba_bc_computational',
# training_validation_split,
# feature_column_str = feature_all_str)
#
# render_report("AJUBA_BC_No_Mut", working_folder,
# 'ajuba_bc_computational',
# training_validation_split,
# feature_column_str = feature_no_mut_str)
#
# render_report("AJUBA_BC_No_Site", working_folder,
# 'ajuba_bc_computational',
# training_validation_split,
# feature_column_str = feature_no_site_str)
#
# render_report("AJUBA_BC_No_Downstream", working_folder,
# 'ajuba_bc_computational',
# training_validation_split,
# feature_column_str = feature_no_downstream_str)
#
# render_report("AJUBA_BC_No_Upstream", working_folder,
# 'ajuba_bc_computational',
# training_validation_split,
# feature_column_str = feature_no_upstream_str)


# TTYH2
# render_report("TTYH2_All", working_folder,
# 'ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_all_str)
#
# render_report("TTYH2_No_Mut", working_folder,
# 'ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_mut_str)
#
# render_report("TTYH2_No_Site", working_folder,
# 'ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_site_str)
#
# render_report("TTYH2_No_Downstream", working_folder,
# 'ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_downstream_str)
#
# render_report("TTYH2_No_Upstream", working_folder,
# 'ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_upstream_str)


# Combined
# render_report("Combined_All", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_all_str)
#
# render_report("Combined_No_Mut", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_mut_str)
#
# render_report("Combined_No_Site", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_site_str)
#
# render_report("Combined_No_Downstream", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_downstream_str)
#
# render_report("Combined_No_Upstream", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_computational',
# training_validation_split,
# feature_column_str = feature_no_upstream_str)

