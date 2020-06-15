#!/usr/bin/env python

# ================================================================
# Python Script - Analysis - Mutation Categorization
#
# It loads the "RNA Lib Structure Summary" file (in csv format) and check the "mutation types".
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "csv" format.
# - Mutation Type
#
# ## Output
# - The feature output is in `csv` format.
#
#

import os
import sys
import argparse
import subprocess

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


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser = argparse.ArgumentParser(description='RNA Lib Analysis Script - Mutation Categorization.')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.csv',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Parameters
arguments_parser.add_argument('--mut_type',
                              action="store", default='single',
                              help='Mutation type - single | double. Default: "single".')

# Output
arguments_parser.add_argument('--out', dest='out',
                              action='store', default='output.csv',
                              help='Filename / path of features output file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)


# Validation
if not os.path.exists(rna_lib_struct_summary_file_path.strip()):
    raise ValueError('"RNA Lib Structure Summary" file does not exist. ')

#
mut_type = args.mut_type.strip().lower()

output_file_path = args.out.strip()

# endregion

# ----------------------------------
# region Global Config

# The mutation categories
# two factors
#   - mutation, like `AtoG`, etc.
#   - If "mutation" breaks the original pair - "1" as Break
mut_categories = {
    'AtoG-0': 'Transition',
    'GtoA-0': 'Transition',
    'CtoU-0': 'Transition',
    'UtoC-0': 'Transition',
    'AtoC-0': 'Transversion',
    'AtoU-0': 'Transversion',
    'GtoC-0': 'Transversion',
    'GtoU-0': 'Transversion',
    'CtoA-0': 'Transversion',
    'CtoG-0': 'Transversion',
    'UtoA-0': 'Transversion',
    'UtoG-0': 'Transversion',

    'AtoG-1': 'Transition+Break',
    'GtoA-1': 'Transition+Break',
    'CtoU-1': 'Transition+Break',
    'UtoC-1': 'Transition+Break',
    'AtoC-1': 'Transversion+Break',
    'AtoU-1': 'Transversion+Break',
    'GtoC-1': 'Transversion+Break',
    'GtoU-1': 'Transversion+Break',
    'CtoA-1': 'Transversion+Break',
    'CtoG-1': 'Transversion+Break',
    'UtoA-1': 'Transversion+Break',
    'UtoG-1': 'Transversion+Break',
}

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

def parse_mut_syntax(mut_syntax_str: str):
    """
    Parse the mutation syntax, and return the following info:
    - mut position
    - mut

    Parameters
    ----------
    mut_syntax_str str The mutation syntax string.

    Returns
    -------

    """

    #
    import re

    # Regex format to match the string
    regex = "^(?P<mut_position>[0-9]+)(?P<mutation>\w+)$"
    matched = re.match(regex, mut_syntax_str.strip())

    if not matched:
        return None

    #
    return matched.groupdict()


#
rna_lib_struct_summary_list = []

# WT
wt_struct_summary = {}

with open(rna_lib_struct_summary_file_path) as infile:
    reader = csv.DictReader(infile)
    #
    for item in reader:
        #
        mut_syntax_str = item['mutation_syntax']

        # Check if it is "wt"
        if mut_syntax_str.startswith('wt'):
            wt_struct_summary = item
            continue

        if mut_syntax_str.startswith('indel'):
            continue

        # Determine it is to check "single" or "double" mutations
        mut_syntax_list = mut_syntax_str.split(',')
        if mut_type == 'single' and len(mut_syntax_list) == 1:
            #
            rna_lib_struct_summary_list.append(item)

# endregion


# ----------------------------------
# region Process Data

def determine_mut_category(mut_position, mut_str, dot_bracket_annotation):
    #
    wt_dot_bracket_annotation = wt_struct_summary['computational_dot_bracket_annotation']

    position_annotation = dot_bracket_annotation[int(mut_position) - 1]
    wt_position_annotation = wt_dot_bracket_annotation[int(mut_position) - 1]

    # Check if the "mut" breaks the "pair"
    flag_break_pair = 0
    if wt_position_annotation == 'S' and position_annotation != wt_position_annotation:
        flag_break_pair = 1

    #
    mut_category = mut_categories['{}-{}'.format(mut_str, flag_break_pair)]

    #
    return position_annotation, wt_position_annotation, mut_category


# Output Data
output_list = []

#
for rna in rna_lib_struct_summary_list:
    #
    mut_syntax_str = rna['mutation_syntax']
    dot_bracket_annotation = rna['computational_dot_bracket_annotation']

    # Parse mut syntax
    mut_syntax_list = mut_syntax_str.split(',')
    mut_syntax = mut_syntax_list[0]
    mut_syntax_parsed = parse_mut_syntax(mut_syntax)

    # Determine the "mut category"
    mut_position = mut_syntax_parsed['mut_position']
    mut_str = mut_syntax_parsed['mutation']

    wt_sequence_str = wt_struct_summary['sequence_string']
    #
    position_annotation, wt_position_annotation, mut_category \
        = determine_mut_category(mut_position, mut_str, dot_bracket_annotation)

    #
    rna_item = []
    rna_item.append(rna['rna_id'])
    rna_item.append(rna['A-to-I_editing_level'])
    rna_item.append(rna['mutation_syntax'])
    rna_item.append(mut_position)
    rna_item.append(wt_sequence_str[int(mut_position) - 1])
    rna_item.append(wt_position_annotation)
    rna_item.append(position_annotation)
    rna_item.append(mut_category)

    #
    output_list.append(rna_item)

# endregion


# ----------------------------------
# region Output

headers = [
    'rna_id',
    'editing_value',
    'mutation_syntax',
    'mut_position',
    'wt_position_nt',
    'wt_position_annotation',
    'position_annotation',
    'mut_category',
]

# Write the "headers" line
writer = csv.writer(open(output_file_path, 'w'))
writer.writerow(headers)

# Write the content rows
for entry in output_list:
    # Replace all "None" value with string "None"
    # writer.writerow([val if val is not None else "None" for val in entry])
    writer.writerow(entry)

# endregion
