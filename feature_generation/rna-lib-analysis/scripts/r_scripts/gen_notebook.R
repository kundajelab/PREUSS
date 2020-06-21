#!/usr/bin/env Rscript

## =========================================
# R Script - Generate Final Reports

library(rmarkdown);

# Load `pandoc`
# Call `Sys.getenv("RSTUDIO_PANDOC")` in "Rstudio" to get the actual path of `pandoc`
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

## ---- Prep -------------------

# Get the "dir" of THIS running script
initial.options <- commandArgs(trailingOnly = FALSE)
# Get the "input" of the R script
script.parameter <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
# Be sure to "normailize" the path just in case it contains "relative" path
script.name = normalizePath(file.path(getwd(), script.parameter))

# If you are running inside "RStudio"
# library(rstudioapi)
# script.name = rstudioapi::getActiveDocumentContext()$path

script.dir <- dirname(script.name)
root_dir <- dirname(dirname(script.dir))

## ---- Generate Reports -------------------

render_report = function(title, working_folder, dataset_codenames, split, test_dataset_codenames = '') {
    #
    rmarkdown::render(
        "./notebooks/ml-practice-r/rf_basic.rmd",
        params = list(
            title = title,
            working_folder = working_folder,
            dataset_codenames = dataset_codenames,
            training_validation_split = split,
            external_test_dataset_codename = test_dataset_codenames
        ),
        output_file = normalizePath(file.path(working_folder, paste0("_output/reports/adar_editing_prediction.", tolower(title), ".html")))
    )
}

working_folder = '/Users/cowfox/Desktop/_rna_lib_analysis/__analysis/2018_12_06-ML_Random_Forest'
training_validation_split = 0.8

# # Neil1 - Computational
# render_report("Neil1_Computational", working_folder,
# 'neil1_computational',
# training_validation_split)
#
# # Neil1 - Experimental
# render_report("Neil1_Experimental", working_folder,
# 'neil1_experimental',
# training_validation_split)
#
# # AJUBA_BC - Computational
# render_report("AJUBA_BC_Computational", working_folder,
# 'ajuba_bc_computational',
# training_validation_split)
#
# # TTYH2_BC - Computational
# render_report("TTYH2_BC_Computational", working_folder,
# 'ttyh2_bc_computational',
# training_validation_split,
# 'ttyh2_bc_degenerate_computational')
#
# # TTYH2_ECS_BC - Computational
# render_report("TTYH2_ECS_BC_Computational", working_folder,
# 'ttyh2_ecs_bc_computational',
# training_validation_split)
#
# # TTYH2_ECS_BC - Experimental
# render_report("TTYH2_ECS_BC_Experimental", working_folder,
# 'ttyh2_ecs_bc_experimental',
# training_validation_split)
#
# # TTYH2 Dataset - Computational
# render_report("TTYH2_Computational", working_folder,
# 'ttyh2_bc_computational, ttyh2_ecs_bc_computational',
# training_validation_split)
#
# # All - Computational
# render_report("Computational", working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_bc_computational, ttyh2_ecs_bc_computational, ttyh2_computational',
# training_validation_split)


# Test with "Degenerate Data"
# Neil1 - Computational
# render_report("Neil1_Computational_Test_with_Degenerate",
# working_folder,
# 'neil1_computational',
# training_validation_split,
# 'neil1_degenerate_computational')

#
# render_report("TTYH2_Computational_Test_with_Degenerate",
# working_folder,
# 'ttyh2_computational',
# training_validation_split,
# 'ttyh2_bc_degenerate_computational')

# All - Computational
# render_report("Computational_Test_with_Degenerate",
# working_folder,
# 'neil1_computational, ajuba_bc_computational, ttyh2_bc_computational, ttyh2_ecs_bc_computational, ttyh2_computational',
# training_validation_split,
# 'neil1_degenerate_computational, ttyh2_bc_degenerate_computational')