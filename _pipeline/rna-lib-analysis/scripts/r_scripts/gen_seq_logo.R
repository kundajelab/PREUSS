#!/usr/bin/env Rscript

## =========================================
# R Script - Generate "Sequence Logo" Chart
#
# Date: 2019-01-10
#
# Example:
#   Rscript  gen_seq_logo.R --pwd /d1/d2/d3/
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

render_report = function(title, working_folder,
                            dataset_source,
                            editing_value_min, editing_value_max,
                            sequence_range) {
    #
    rmarkdown::render(
        "./notebooks/seq_logo.Rmd",
        params = list(
            title = title,
            working_folder = working_folder,

            #
            dataset_source = dataset_source,
            editing_value_min = editing_value_min,
            editing_value_max = editing_value_max,
            sequence_range = sequence_range
        ),
        output_file = normalizePath(file.path(working_folder, paste0("_output/reports/seq_logo.", tolower(title), ".html")))
    )
}

# Defaults
sequence_range_default = 10
working_folder = '/Users/cowfox/Desktop/_rna_lib_analysis/__analysis/'
if (!is.null(opt$pwd)){
  # print_help(opt_parser)
  # stop("At least one argument must be supplied (input file).n", call.=FALSE
    working_folder = opt$pwd
}


# Neil1
# render_report(
#     "Neil1_Computational-lt_0.03",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.03,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "Neil1_Computational-0.03_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.03,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "Neil1_Computational-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )

# TTYH2_BC
# render_report(
#     "TTYH2_BC_Computational-lt_0.03",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.03,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_BC_Computational-0.03_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.03,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_BC_Computational-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )

# TTYH2_Combined
# render_report(
#     "TTYH2_Combined_Computational-lt_0.03",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_Combined-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.03,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_Combined_Computational-0.03_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_Combined-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.03,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_Combined_Computational-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_Combined-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )


# AJUBA_BC
# render_report(
#     "AJUBA_BC_Computational-lt_0.03",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/AJUBA_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.03,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "AJUBA_BC_Computational-0.03_0.1",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/AJUBA_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.03,
#     editing_value_max = 0.1,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "AJUBA_BC_Computational-gt_0.1",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/AJUBA_BC-rna_lib-structure_summary-indel_fixed.csv',
#     editing_value_min = 0.1,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )

# Neil1 - Degenerate
# render_report(
#     "Neil1_Degenrate-lt_0.1",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-Degenerated_Isoforms-structure_summary-fixed.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.1,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "Neil1_Degenrate-0.1_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-Degenerated_Isoforms-structure_summary-fixed.csv',
#     editing_value_min = 0.1,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "Neil1_Degenrate-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/Neil1-Degenerated_Isoforms-structure_summary-fixed.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )

# TTYH2_BC-Degenerate - Degenerate
# render_report(
#     "TTYH2_BC_Degenrate-lt_0.1",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-Degenerated_Isoforms-structure_summary.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.1,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_BC_Degenrate-0.1_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-Degenerated_Isoforms-structure_summary.csv',
#     editing_value_min = 0.1,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "TTYH2_BC_Degenrate-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/TTYH2_BC-Degenerated_Isoforms-structure_summary.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )


# Genomic - HEK293.WT1
# render_report(
#     "HEK293.WT1-lt_0.02",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-HEK293.WT1-filtered.csv',
#     editing_value_min = NA,
#     editing_value_max = 0.02,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "HEK293.WT1-0.02_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-HEK293.WT1-filtered.csv',
#     editing_value_min = 0.02,
#     editing_value_max = 0.3,
#     sequence_range = sequence_range_default
# )
#
# render_report(
#     "HEK293.WT1-gt_0.3",
#     working_folder,
#     dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-HEK293.WT1-filtered.csv',
#     editing_value_min = 0.3,
#     editing_value_max = NA,
#     sequence_range = sequence_range_default
# )

# Genomic - Hela_no_gRNA_merge
render_report(
    "Hela_no_gRNA_merge.sorted-lt_0.02",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA_merge.sorted-filtered.csv',
    editing_value_min = NA,
    editing_value_max = 0.02,
    sequence_range = sequence_range_default
)

render_report(
    "Hela_no_gRNA_merge.sorted-0.02_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA_merge.sorted-filtered.csv',
    editing_value_min = 0.02,
    editing_value_max = 0.3,
    sequence_range = sequence_range_default
)

render_report(
    "Hela_no_gRNA_merge.sorted-gt_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA_merge.sorted-filtered.csv',
    editing_value_min = 0.3,
    editing_value_max = NA,
    sequence_range = sequence_range_default
)

# Genomic - Hela_no_gRNA+IFN_merge
render_report(
    "Hela_no_gRNA+IFN_merge.sorted-lt_0.02",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA+IFN_merge.sorted-filtered.csv',
    editing_value_min = NA,
    editing_value_max = 0.02,
    sequence_range = sequence_range_default
)

render_report(
    "Hela_no_gRNA+IFN_merge.sorted-0.02_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA+IFN_merge.sorted-filtered.csv',
    editing_value_min = 0.02,
    editing_value_max = 0.3,
    sequence_range = sequence_range_default
)

render_report(
    "Hela_no_gRNA+IFN_merge.sorted-gt_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-Hela_no_gRNA+IFN_merge.sorted-filtered.csv',
    editing_value_min = 0.3,
    editing_value_max = NA,
    sequence_range = sequence_range_default
)


# Genomic - u87_ctrl_vs_kd
render_report(
    "u87_ctrl_vs_kd-lt_0.02",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-u87_ctrl_vs_kd-filtered.csv',
    editing_value_min = NA,
    editing_value_max = 0.02,
    sequence_range = sequence_range_default
)

render_report(
    "u87_ctrl_vs_kd-0.02_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-u87_ctrl_vs_kd-filtered.csv',
    editing_value_min = 0.02,
    editing_value_max = 0.3,
    sequence_range = sequence_range_default
)

render_report(
    "u87_ctrl_vs_kd-gt_0.3",
    working_folder,
    dataset_source = '/Users/cowfox/Desktop/_rna_lib_analysis/__rna_library/_structure_summary/_chromosome/chromosone-u87_ctrl_vs_kd-filtered.csv',
    editing_value_min = 0.3,
    editing_value_max = NA,
    sequence_range = sequence_range_default
)

