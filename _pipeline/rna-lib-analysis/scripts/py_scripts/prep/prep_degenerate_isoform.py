#!/usr/bin/env python

# ================================================================
# Python Script - Data Prep - Process "Degenerate Isoforms" data
#
# The original data includes three columns:
#
# - `isoforms` - It list the “nt” on “ALL” mutation positions.
# - `replicate 1` - “Editing Value”, replicate #1.
# - `replicate 2` - “Editing Value”, replicate #2.
#
# ## Input
# - Original data file
#   - It is in "csv" format.
# - Reference Sequence
# - Mutation Positions
#   - It is a "string" format which represents a list of "positions"
# - Editing Position
# - ID Prefix
# - Is use "RNA" as sequence
#
# ## Output
# - The output includes two files:
#   - RNA Lib Structure Summary file, in "JSON" format
#   - Quick summary file, in "csv" format.
#
#

import os
import argparse
import subprocess


import logging

import json
import csv

import numpy as np

from neoRNA.sequence.sequence import Sequence
from neoRNA.util.runner.matlab_runner import MatlabRunner
from neoRNA.util.file_utils import FileUtils


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser = argparse.ArgumentParser(description='RNA Lib Prep Script - Process "Degenerate Isoforms" data.')

# Inputs
arguments_parser.add_argument('original_data_file',
                              metavar='original_data_file.csv',
                              help='The file path to the "Original Data" file.')

# Parameters
arguments_parser.add_argument('--ref_sequence', dest='reference_sequence',
                              action='store', default='',
                              help='The reference sequence.')
arguments_parser.add_argument('--mut_positions', dest='mutation_positions',
                              action='store', default='',
                              help='The mutation positions - by a string.')
arguments_parser.add_argument('--editing_position', dest='editing_position',
                              action='store', default=0,
                              help='The editing position.')
arguments_parser.add_argument('--id_prefix', dest='id_prefix',
                              action='store', default='DG',
                              help='The "prefix" used for generate the ID.')
arguments_parser.add_argument('--use_rna',
                              action='store_true',
                              help='If using RNA sequence.')

# Output
arguments_parser.add_argument('--o_summary', dest='out_csv',
                              action='store', default='output.csv',
                              help='Filename of summary file.')
arguments_parser.add_argument('--o_structure', dest='out_json',
                              action='store', default='output.json',
                              help='Filename of "RNA Structure Summary" file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
original_data_file_path = os.path.abspath(args.original_data_file)


# Validation
if not os.path.exists(original_data_file_path.strip()):
    raise ValueError('"Original Data" file does not exist. ')

#
reference_sequence_str = args.reference_sequence.strip()
mutation_positions = list(map(int, args.mutation_positions.strip().split(',')))  # Convert a string to a list of "int".
editing_position = int(args.editing_position.strip()) if args.editing_position.strip() is not None else 0
id_prefix = args.id_prefix.strip()
use_rna = args.use_rna

output_csv_file_path = args.out_csv.strip()
output_json_file_path = args.out_json.strip()

# endregion


# ----------------------------------
# region Prep

cwd = os.getcwd()
output_csv_file_path = os.path.join(cwd, output_csv_file_path)
output_json_file_path = os.path.join(cwd, output_json_file_path)

# endregion


# ----------------------------------
# region Load Original Data

isoforms = []
editing_value_1_dict = {}
editing_value_2_dict = {}

# Columns to load
isoform_column_name = 'isoforms'
editing_value_1_column_name = 'replicate 1'
editing_value_2_column_name = 'replicate 2'

#
original_data = csv.DictReader(open(original_data_file_path, 'rU'))
for line in original_data:
    isoform = line[isoform_column_name]
    isoforms.append(isoform)
    #
    value_1 = line[editing_value_1_column_name]
    editing_value_1_dict[isoform] = float(value_1) if value_1 != "NA" and value_1 != "#N/A" else None

    value_2 = line[editing_value_2_column_name]
    editing_value_2_dict[isoform] = float(value_2) if value_2 != "NA" and value_2 != "#N/A" else None

# endregion


# ----------------------------------
# region Process Isoform

def generate_isoform_sequence(reference_sequence, mutation_positions, isoform):
    r"""
    Generate the sequence based on the given "isoform" sequence.

    Parameters
    ----------
    reference_sequence: str
        The reference sequence.
    mutation_positions: str
        The list of "mutation positions" (Start with "1").

    Returns
    -------
    isoform_sequence: str
        The isoform sequence

    """

    if len(reference_sequence) == 0:
        return ''

    # The length of "mut positions" and "isoform" string should be identical.
    if len(mutation_positions) != len(isoform):
        return ''

    isoform_sequence_list = list(reference_sequence)
    isoform_index = 0
    for position in mutation_positions:
        isoform_sequence_list[position - 1] = isoform[isoform_index]
        isoform_index += 1

    return "".join(isoform_sequence_list)


#
quick_summary_items = []
rna_library_structure_summary_items = list()

#
MatlabRunner.start_engine()
#
id_count = 1
reference_sequence = Sequence(reference_sequence_str)
for isoform in isoforms:
    #
    isoform_sequence_str = generate_isoform_sequence(reference_sequence_str, mutation_positions, isoform)

    if isoform_sequence_str == '':
        continue
    isoform_sequence = Sequence(isoform_sequence_str)

    #
    isoform_id = '{0}{1:0>4}'.format(id_prefix, id_count)
    id_count += 1

    # Mut Syntax
    if use_rna:
        seq_mut_positions, seq_mut_syntax = isoform_sequence.generate_mutation_syntax(
            str(reference_sequence.get_rna_sequence()), seq_type_rna=True)
    else:
        seq_mut_positions, seq_mut_syntax = isoform_sequence.generate_mutation_syntax(
            str(reference_sequence.get_dna_sequence()))

    # QC on "syntax"
    # Neil1 - discard the record which includes "46AtoG: or "46AtoG"
    # if '46AtoG' in seq_mut_syntax or '48AtoG' in seq_mut_syntax:
    #     continue

    # Decide output sequence string, RNA or DNA
    output_isoform_sequence_str = str(isoform_sequence.get_rna_sequence()) if use_rna \
        else str(isoform_sequence.get_dna_sequence())

    #
    editing_value_avg = (editing_value_1_dict[isoform] + editing_value_2_dict[isoform]) / 2.0

    # Predict "Reference Structure"
    biers_seqpos_out = range(1, len(output_isoform_sequence_str))
    biers_offset = 0  # ALWAYS "0"
    structure_na, bpp_na, bootstrap_na \
        = MatlabRunner.rna_structure(output_isoform_sequence_str, biers_offset, biers_seqpos_out)

    # Data entry for "summary" output
    entry = list()
    entry.append(isoform_id)
    entry.append(isoform)
    entry.append(output_isoform_sequence_str)
    entry.append(','.join(seq_mut_syntax))
    entry.append(len(seq_mut_positions))
    entry.append(editing_position)
    entry.append(editing_value_1_dict[isoform])
    entry.append(editing_value_2_dict[isoform])
    entry.append(editing_value_avg)

    quick_summary_items.append(entry)

    # Data entry for "structure summary" output
    summary_item = {
        "rna_id": isoform_id,
        "sequence_string": output_isoform_sequence_str,
        "mutation_syntax": ','.join(seq_mut_syntax) if seq_mut_syntax else None,

        "computational_structure": structure_na,
        "computational_bpp": np.array(bpp_na).tolist(),

        #
        "A-to-I_editing_level": editing_value_avg,
        "A-to-I_editing_site": editing_position,
    }

    rna_library_structure_summary_items.append(summary_item)

#
MatlabRunner.stop_engine()

# endregion


# ----------------------------------
# region Output - Quick Summary file

headers = [
    'dg_id',
    'isoform',
    'isoform_sequence',
    'mut_syntax',
    'mut_count',
    'editing_position',
    'editing_value_1',
    'editing_value_1',
    'editing_value_avg',
]

writer = csv.writer(open(output_csv_file_path, 'w'))
writer.writerow(headers)

for entry in quick_summary_items:
    writer.writerow(entry)

# endregion

# ----------------------------------
# region Output - RNA Lib Structure Summary file

#
rna_library_structure_summary_dict = {
    "items": rna_library_structure_summary_items,
}
json_data = json.dumps(rna_library_structure_summary_dict,
                       sort_keys=True,
                       indent=4, ensure_ascii=True)

# Output JSON as file
FileUtils.save_file(output_json_file_path, json_data)

# endregion
