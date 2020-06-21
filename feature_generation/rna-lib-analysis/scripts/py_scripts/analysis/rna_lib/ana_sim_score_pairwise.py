#!/usr/bin/env python

# ================================================================
# Python Script - Analysis - RNA Structure Similarity Score Pairwise
#
# It loads the "RNA Lib Structure Summary" file (in JSON) and calculates the "similarity score" for both
#   "computational" and "experimental" RNA structures.
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "json" format.
#
# ## Output
# - The feature output is in `csv` format.
#
#

import os
import sys
import argparse

# Add "py_scripts" into module path, relative to "current" script
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.realpath(__file__)))))
sys.path.append(local_module_path)

import json
import csv

from collections import defaultdict

import logging
from py_scripts import setup_logging

from typing import Tuple, List, Any

from neoRNA.sequence.sequence import Sequence

from neoRNA.util.runner.simtree_runner import SimTreeRunner


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser = argparse.ArgumentParser(description='RNA Lib Analysis Script - RNA Structure Similarity Score Pairwise.')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.json',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Parameters


# Output
arguments_parser.add_argument('--out', dest='out',
                              action='store', default='output.csv',
                              help='Filename / path of output file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)


# Validation
if not os.path.exists(rna_lib_struct_summary_file_path.strip()):
    raise ValueError('"RNA Lib Structure Summary" file does not exist. ')

#

output_file_path = args.out.strip()

# endregion

# ----------------------------------
# region Global Config

# endregion


# ----------------------------------------------------------------
# region Set Logger
#

setup_logging(logging_level=logging.INFO)
# Logger Name
logger_name = 'neo_rna.script'
logger = logging.getLogger(logger_name)

# endregion


# ----------------------------------
# region Prep

cwd = os.getcwd()
if not os.path.isabs(output_file_path):
    output_file_path = os.path.join(cwd, output_file_path)

# endregion

# ----------------------------------
# region Load "RNA Lib Structure Summary" Data

# RSample results, indexed by "RNA ID"
rna_lib_struct_summary_dict = {}

with open(rna_lib_struct_summary_file_path) as infile:
    rna_lib_struct_summary_dict = json.load(infile)

# endregion

# ----------------------------------
# region Process

# WT
wt_dot_bracket_structure_computational = None
wt_dot_bracket_structure_experimental = None

# Prepare the "secondary structure" for all RNA items, also identify "WT"
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']

    # Both "computational_structure" and "experimental_structure"
    structure_string_computational = rna_item['computational_structure']
    structure_string_experimental = rna_item['experimental_structure']

    # Check if it is "WT"
    if 'mutation_syntax' in rna_item \
            and rna_item['mutation_syntax'].strip().lower().startswith('wt'):
        #
        wt_dot_bracket_structure_computational = structure_string_computational
        wt_dot_bracket_structure_experimental = structure_string_experimental


#
# The runner for SimTree
simtree_jar_location = '/Users/cowfox/Desktop/_rna_lib_analysis/__tool/SimTree_v1.2.3.jar'
simtree_runner = SimTreeRunner(simtree_jar_location, os.getcwd())

#
output_data_list = list()
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']
    logger.info('---------- RNA Item: {} ----------'.format(rna_id))

    sequence_string = rna_item['sequence_string']
    sequence = Sequence(sequence_string)

    # Both "computational_structure" and "experimental_structure"
    structure_string_computational = rna_item['computational_structure']
    structure_string_experimental = rna_item['experimental_structure']

    #
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    # RNA Structure Similarity against WT
    simtree_normalized_score_computational = None
    simtree_flipping_modes_computational = None
    if wt_dot_bracket_structure_computational:
        simtree_normalized_score_computational, simtree_flipping_modes_computational \
            = simtree_runner.compare_rna_structure(structure_string_computational,
                                                   wt_dot_bracket_structure_computational)

    simtree_normalized_score_experimental = None
    simtree_flipping_modes_experimental = None
    if wt_dot_bracket_structure_experimental:
        simtree_normalized_score_experimental, simtree_flipping_modes_experimental \
            = simtree_runner.compare_rna_structure(structure_string_experimental,
                                                   wt_dot_bracket_structure_experimental)

    #
    data_entry = list()
    data_entry.append(rna_id)
    data_entry.append(editing_value)
    data_entry.append(simtree_normalized_score_computational)
    data_entry.append(simtree_normalized_score_experimental)

    output_data_list.append(data_entry)

# endregion


# ----------------------------------
# region Output

headers = [
    'rna_id',
    'editing_value',
    'sim_nor_score_computational',
    'sim_nor_score_experimental',
    ]

# Write the "headers" line
writer = csv.writer(open(output_file_path, 'w'))
writer.writerow(headers)

# Write the content rows
for entry in output_data_list:
    # Replace all "None" value with string "None"
    # writer.writerow([val if val is not None else "None" for val in entry])
    writer.writerow(entry)

# endregion
