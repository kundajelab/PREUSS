#!/usr/bin/env bash

## =========================================
# Bash Script - Run Python Script to Generate Features.
#
# - Run `bpRNA` for each dataset - both "computational" and "experimental" if available
# - Generate feature metrics
#

# Full path of the "THIS" script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
# The "ROOT" is the parent folder of the script.
ROOT_DIR="$( dirname "${SCRIPT_DIR}" )"

# Change "PWD" to "ROOT" folder
#cd "${ROOT_DIR}"

## ---- Prep -------------------
## Unzip the data source file
#DATA_SOURCE_DIR="${ROOT_DIR}/_data_source"

# It uese `unzip` (MACOSX)
#for i in `find . -name "*.zip" -type f`; do unzip -o "$i" -d "${DATA_SOURCE_DIR}"; done

cd /Users/cowfox/Desktop/_rna_lib_analysis/__analysis/2018_12_10-ML_Random_Forest-Tuning


## ---- Features -------------------
# Neil1 - Computational
#ana_gen_ml_features.py \
#_data_source/Neil1-rna_lib-structure_summary-indel_fixed.json \
#--data_type computational \
#--wt AGCCUGCCCUCUGAUCUCUGCCUGUUCCUCUGUCCCACAGGGGGCAAAGGCUACGGGUCAGAGAGCGGGGAGGAGGAC \
#--out _output/feature_matrices/neil1_computational.features.csv \
#--out_bprna _output/bpRNA_output/neil1_computational.bprna.txt
#
## AJUBA_BC - Computational
##ana_gen_ml_features.py \
##_data_source/AJUBA_BC-rna_lib-structure_summary-indel_fixed.json \
##--data_type computational \
##--wt UUUUGGGGUUGUGGUUGAUGCAGUGUGGGAUGUCCCUGAGAGGUAGCAAGUCUAGGGUGGUGAGUUCCUGCUAGGCAACCAAAUUUGUACACUUGGUUUCCUAGUAGAAGCUCACUUGCCACCUCUCAGAGGGGUCCCGGAUUGCAUCCAUCACAAUCCCAAAAC \
##--out _output/feature_matrices/ajuba_bc_computational.features.csv \
##--out_bprna _output/bpRNA_output/ajuba_bc_computational.bprna.txt
#
## AJUBA_BC - Computational
## Sequence Sliced
#ana_gen_ml_features.py \
#_data_source/AJUBA_BC-rna_lib-structure_summary-indel_fixed-seq_sliced.json \
#--data_type computational \
#--wt UUUUGGGGUUGUGGUUGAUGCAGUGUGGGAUGUCCCUGAGAGGUAGCAAGUCUAGGGUGUUGCCACCUCUCAGAGGGGUCCCGGAUUGCAUCCAUCACAAUCCCAAAAC \
#--out _output/feature_matrices/ajuba_bc_computational.features.csv \
#--out_bprna _output/bpRNA_output/ajuba_bc_computational.bprna.txt
#
## TTYH2_BC - Computational
#ana_gen_ml_features.py \
#_data_source/TTYH2_BC-rna_lib-structure_summary-indel_fixed.json \
#--data_type computational \
#--wt CAUGCUUCAUACCCAGAGAGAAGCCCCCGGCUGCCCAGGCAUGCUUAGGCUUACACGUGCUUAGGCUUAGGCGUGCCUGGGUGACCAGGGCGCUUCUCUCUGGGUGUGAAGAACU \
#--out _output/feature_matrices/ttyh2_bc_computational.features.csv \
#--out_bprna _output/bpRNA_output/ttyh2_bc_computational.bprna.txt
#
## TTYH2_ECS_BC - Computational
#ana_gen_ml_features.py \
#_data_source/TTYH2_ECS_BC-rna_lib-structure_summary-indel_fixed.json \
#--data_type computational \
#--wt CAUGCUUCAUACCCAGAGAGAAGCCCCCGGCUGCCCAGGCAUGCUUAGGCUUACACGUGCUUAGGCUUAGGCGUGCCUGGGUGACCAGGGCGCUUCUCUCUGGGUGUGAAGAACU \
#--out _output/feature_matrices/ttyh2_ecs_bc_computational.features.csv \
#--out_bprna _output/bpRNA_output/ttyh2_ecs_bc_computational.bprna.txt
#
## ----
## Neil1 - Experimental
#ana_gen_ml_features.py \
#_data_source/Neil1-rna_lib-structure_summary-indel_fixed.json \
#--data_type experimental \
#--wt AGCCUGCCCUCUGAUCUCUGCCUGUUCCUCUGUCCCACAGGGGGCAAAGGCUACGGGUCAGAGAGCGGGGAGGAGGAC \
#--out _output/feature_matrices/neil1_experimental.features.csv \
#--out_bprna _output/bpRNA_output/neil1_experimental.bprna.txt
#
## TTYH2_ECS_BC - Experimental
#ana_gen_ml_features.py \
#_data_source/TTYH2_ECS_BC-rna_lib-structure_summary-indel_fixed.json \
#--data_type experimental \
#--wt CAUGCUUCAUACCCAGAGAGAAGCCCCCGGCUGCCCAGGCAUGCUUAGGCUUACACGUGCUUAGGCUUAGGCGUGCCUGGGUGACCAGGGCGCUUCUCUCUGGGUGUGAAGAACU  \
#--out _output/feature_matrices/ttyh2_ecs_bc_experimental.features.csv \
#--out_bprna _output/bpRNA_output/ttyh2_ecs_bc_experimental.bprna.txt
#
## ----
## Neil1 - Degenerate
#ana_gen_ml_features.py \
#_data_source/Neil1-Degenerated_Isoforms-structure_summary-fixed.json \
#--data_type computational \
#--wt AGCCUGCCCUCUGAUCUCUGCCUGUUCCUCUGUCCCACAGGGGGCAAAGGCUACGGGUCAGAGAGCGGGGAGGAGGAC \
#--out _output/feature_matrices/neil1_degenerate_computational.features.csv \
#--out_bprna _output/bpRNA_output/neil1_degenerate_computational.bprna.txt
#
#
## TTYH2_BC - Degenerate
#ana_gen_ml_features.py \
#_data_source/TTYH2_BC-Degenerated_Isoforms-structure_summary.json \
#--data_type computational \
#--wt CAUGCUUCAUACCCAGAGAGAAGCCCCCGGCUGCCCAGGCAUGCUUAGGCUUACACGUGCUUAGGCUUAGGCGUGCCUGGGUGACCAGGGCGCUUCUCUCUGGGUGUGAAGAACU \
#--out _output/feature_matrices/ttyh2_bc_degenerate_computational.features.csv \
#--out_bprna _output/bpRNA_output/ttyh2_bc_degenerate_computational.bprna.txt

# -----
ana_gen_ml_features.py \
_data_source/chromosone-HEK293.WT1-filtered.json \
--data_type computational \
--wt NO_USE \
--out _output/feature_matrices/HEK293.WT1.features.csv \
--out_bprna _output/bpRNA_output/HEK293.WT1.bprna.txt

ana_gen_ml_features.py \
_data_source/chromosone-Hela_no_gRNA_merge.sorted-filtered.json \
--data_type computational \
--wt NO_USE \
--out _output/feature_matrices/Hela_no_gRNA.features.csv \
--out_bprna _output/bpRNA_output/Hela_no_gRNA.bprna.txt

ana_gen_ml_features.py \
_data_source/chromosone-Hela_no_gRNA+IFN_merge.sorted-filtered.json \
--data_type computational \
--wt NO_USE \
--out _output/feature_matrices/Hela_no_gRNA+IFN.features.csv \
--out_bprna _output/bpRNA_output/Hela_no_gRNA+IFN.bprna.txt

ana_gen_ml_features.py \
_data_source/chromosone-u87_ctrl_vs_kd-filtered.json \
--data_type computational \
--wt NO_USE \
--out _output/feature_matrices/u87_ctrl_vs_kd.features.csv \
--out_bprna _output/bpRNA_output/u87_ctrl_vs_kd.bprna.txt


