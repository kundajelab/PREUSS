#!/usr/bin/env python

# ================================================================
# "RNA Lib" Util Script - Adjust RNA Sequence
#
# It loads the "RNA Lib Structure Summary" file (in JSON) and adjust the sequence based on requests.
#
# NOTE:
# - If the "adjustment" affect the "editing position" / "mutation position", also need to update it.
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "json" format.
# - "Sequence" to be "extracted"
#
#

import os
import argparse

import numpy as np

import json

# Add "py_scripts" into module path, relative to "current" script
import sys
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.realpath(__file__))))
sys.path.append(local_module_path)

import logging
from py_scripts import setup_logging

from neoRNA.sequence.sequence import Sequence
from neoRNA.util.runner.matlab_runner import MatlabRunner
from neoRNA.util.file_utils import FileUtils

# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Util - Adjust RNA Sequence')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.json',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Params
arguments_parser.add_argument('--wt', dest='wt_sequence',
                              action="store",
                              help='The "WT sequence" (already extracted).')
arguments_parser.add_argument('--extract', dest='sequence_to_extract',
                              action="store",
                              help='The "sequence" to be extracted.')


# Output
arguments_parser.add_argument('--out',
                              action='store',
                              help='The updated "RNA Lib Structure Summary" file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)

# Validation
if not os.path.exists(rna_lib_struct_summary_file_path):
    raise ValueError('"RNA Lib" definition file does not exist.')

#
sequence_to_extract = args.sequence_to_extract.strip()
wt_sequence_str = args.wt_sequence.strip()

# Output
output_file_path = args.out

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

if not os.path.isabs(output_file_path):
    cwd = os.getcwd()
    output_folder_path = os.path.join(cwd, output_file_path)

# endregion


# ----------------------------------
# region Load "RNA Lib Structure Summary" Data

# RSample results, indexed by "RNA ID"
rna_lib_struct_summary_dict = {}

with open(rna_lib_struct_summary_file_path) as infile:
    rna_lib_struct_summary_dict = json.load(infile)

# endregion


# ----------------------------------
# region Generate

#
rna_item_list_updated = list()

#
MatlabRunner.start_engine()
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']
    logger.info('---------- RNA Item: {} ----------'.format(rna_id))

    sequence_string = rna_item['sequence_string']
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    mutation_syntax_str = rna_item['mutation_syntax']

    #
    computational_structure = rna_item['computational_structure']
    computational_bpp = rna_item['computational_bpp']

    # Determine if the RNA item includes the "target sequence"
    if len(sequence_to_extract) > 0:
        #
        found_position_start = sequence_string.find(sequence_to_extract)
        if found_position_start >= 0:
            # Found
            found_position_end = found_position_start + len(sequence_to_extract)    # "end" NOT included

            # Slice
            sequence_string_length = len(sequence_string)
            sequence_string_sliced = sequence_string[0:found_position_start] \
                                     + sequence_string[found_position_end:sequence_string_length]
            sequence_string = sequence_string_sliced

            # Check if need to update the "editing position"
            if editing_position >= found_position_end:
                editing_position = editing_position - len(sequence_to_extract)

            # Determine "mutation syntax"
            sequence = Sequence(sequence_string)
            seq_mut_positions, seq_mut_syntax = sequence.generate_mutation_syntax(wt_sequence_str, seq_type_rna=True)
            #
            if seq_mut_syntax is None:
                mutation_syntax_str = 'indel'
            elif len(seq_mut_syntax) == 0:
                # WT
                mutation_syntax_str = 'wt'
            elif len(seq_mut_syntax) > 0:
                mutation_syntax_str = ','.join(seq_mut_syntax)

            # Update RNA Structure
            seqpos_out = range(1, len(sequence_string))
            structure_na, bpp_na, bootstrap_na = MatlabRunner.rna_structure(sequence_string, 0, seqpos_out)
            computational_structure = structure_na
            computational_bpp = np.array(bpp_na).tolist()

    #
    item_data = {
        "rna_id": rna_id,
        "sequence_string": sequence_string,

        #
        "A-to-I_editing_level": editing_value,
        "A-to-I_editing_site": editing_position,
        "mutation_syntax": mutation_syntax_str,

        #
        "computational_structure": computational_structure,
        "computational_bpp": computational_bpp,
    }
    rna_item_list_updated.append(item_data)

# endregion


# ----------------------------------
# region Output

rna_lib_struct_summary_dict['items'] = rna_item_list_updated
json_data = json.dumps(rna_lib_struct_summary_dict,
                       sort_keys=True,
                       indent=4, ensure_ascii=True)

# Output JSON as file
FileUtils.save_file(output_file_path, json_data)

# endregion
