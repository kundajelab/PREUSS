#!/usr/bin/env python

# ================================================================
# Python Script - Data Prep - Convert a "RNA data sheet" to a "Structure Summary" JSON file.
#
# The original data includes the following columns:
#
# - `sequence_name` - The "name", will be used as "ID".
# - `editing_position` - The "editing position".
# - `editing_level` - The "editing level" value.
# - `sequence` - The sequence.
#
# ## Input
# - "RNA data sheet", in `csv` format.
# - ID Prefix
# - Is use "RNA" as sequence
#
# ## Output
# - RNA Lib Structure Summary file, in "JSON" format
#
#

import os
import argparse
import subprocess


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

from typing import List, Any, Dict

import json
import csv

import numpy as np

from neoRNA.sequence.sequence import Sequence
from neoRNA.util.runner.matlab_runner import MatlabRunner
from neoRNA.util.file_utils import FileUtils


# ----------------------------------
# region Global Variables - Configs

# Column Names
col_sequence_name = 'Gene'
col_editing_position = 'RNA_Editing_Site'
col_editing_level = 'editlvl'
col_sequence_str = 'Sequence_From_Gokul_WT1'

# endregion

# ----------------------------------
# region Global Variables - Basic

# 'folder path' of THIS script file.
script_dir = os.path.dirname(os.path.realpath(__file__))
# 'Root' folder
root_dir = os.path.dirname(script_dir)

# Logger
logger = None

# endregion


# ----------------------------------
# region Global Variables - Arguments

original_data_file_path = None
id_prefix = None
use_rna = None

output_file_path = None

# endregion


# ----------------------------------
# region Global Variables

original_data_dict = dict()

output_data_dict = dict()
output_data_dict['items'] = []

# endregion


# ----------------------------------
# region Parse Arguments

def parse_args():

    # region Parser
    # -----------------------------------------
    # Parser
    # -----------------------------------------
    arguments_parser \
        = argparse.ArgumentParser(description='Convert a "RNA data sheet" to a "Structure Summary" JSON file.')

    # Inputs
    arguments_parser.add_argument('original_data_file',
                                  metavar='original_data_file.csv',
                                  help='The file path to the "Original Data" file.')

    # Parameters
    arguments_parser.add_argument('--id_prefix', dest='id_prefix',
                                  action='store', default='',
                                  help='The "prefix" used for generate the ID.')
    arguments_parser.add_argument('--use_rna',
                                  action='store_true',
                                  help='If using RNA sequence.')

    # Output
    arguments_parser.add_argument('--o_structure', dest='out_json',
                                  action='store', default='output.json',
                                  help='Filename of "RNA Structure Summary" file.')

    # endregion

    #
    global original_data_file_path
    global id_prefix
    global use_rna
    global output_file_path

    args = arguments_parser.parse_args()

    #
    original_data_file_path = os.path.abspath(args.original_data_file)
    id_prefix = args.id_prefix.strip()
    use_rna = args.use_rna
    #
    output_file_path = args.out_json.strip()

    # region Validate
    # -----------------------------------------
    # Validate
    # -----------------------------------------

    validate_args()

    # endregion


#
def validate_args():
    #
    global original_data_file_path

    if not os.path.exists(original_data_file_path):
        raise ValueError('The "Original Data" file does not exist. ')

# endregion


# ----------------------------------
# region Prep

def get_logger(logger_name: str, logging_level=logging.INFO):
    setup_logging(logging_level=logging_level)
    #
    return logging.getLogger(logger_name)


def prep():
    #
    global logger
    global output_file_path

    #
    logger = get_logger('neo.script')

    #
    cwd = os.getcwd()
    if not os.path.isabs(output_file_path):
        output_file_path = os.path.join(cwd, output_file_path)

# endregion


# ----------------------------------
# region Load Data

def load_data():
    #
    global original_data_file_path

    global original_data_dict

    global col_sequence_name
    global col_editing_position
    global col_editing_level
    global col_sequence_str

    #
    with open(original_data_file_path, 'rU') as infile:
        #
        data_dict = csv.DictReader(infile)

        #
        for line_item in data_dict:
            #
            sequence_name = line_item[col_sequence_name]

            #
            original_data_dict[sequence_name] = {
                col_sequence_name: sequence_name,
                col_editing_position: int(line_item[col_editing_position]),
                col_editing_level: float(line_item[col_editing_level]),
                col_sequence_str: line_item[col_sequence_str],
            }


# endregion


# ----------------------------------
# region Process

def process():
    #
    #
    MatlabRunner.start_engine()

    #
    for rna_id in original_data_dict:
        #
        rna_item = original_data_dict[rna_id]

        #
        sequence_name = rna_item[col_sequence_name]
        editing_position = int(rna_item[col_editing_position])
        editing_level = float(rna_item[col_editing_level])
        sequence_str = rna_item[col_sequence_str]

        #
        if sequence_str == '':
            continue
        sequence = Sequence(sequence_str)

        #
        output_sequence_str = str(sequence.get_rna_sequence()) if use_rna \
            else str(sequence.get_dna_sequence())

        #
        # Predict "Reference Structure"
        biers_seqpos_out = range(1, len(output_sequence_str))
        biers_offset = 0  # ALWAYS "0"
        structure_na, bpp_na, bootstrap_na \
            = MatlabRunner.rna_structure(output_sequence_str, biers_offset, biers_seqpos_out)

        #
        summary_item = {
            "rna_id": sequence_name,
            "sequence_string": output_sequence_str,

            "computational_structure": structure_na,
            # "computational_bpp": np.array(bpp_na).tolist(),

            #
            "A-to-I_editing_level": editing_level,
            "A-to-I_editing_site": editing_position,
        }

        output_data_dict['items'].append(summary_item)

    #
    MatlabRunner.stop_engine()

# endregion


# ----------------------------------
# region Output

def output():
    """
    Output as JSON file.
    """
    #
    global output_data_dict
    global output_file_path

    #
    json_contents = json.dumps(output_data_dict,
                               sort_keys=True,
                               indent=4, ensure_ascii=False)
    #
    with open(output_file_path, 'w') as outfile:
        outfile.write(json_contents)


# endregion


# ----------------------------------
# region Main Script

def main():
    # Parse Args
    parse_args()

    # Prep
    prep()

    # Load original data
    load_data()

    # Process data & output
    process()

    # Output
    output()


#
if __name__ == "__main__":
    main()

# endregion


