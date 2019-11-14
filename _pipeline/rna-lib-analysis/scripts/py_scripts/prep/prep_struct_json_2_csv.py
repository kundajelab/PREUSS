#!/usr/bin/env python

# ================================================================
# Python Script - Data Prep - Convert RNA Lib Structure Summary File (in "JSON") to a CSV File
#
# ## Input
# - RNA Lib Structure Summary file, in "JSON" format
#
# ## Output
# - RNA Lib Structure Summary file, in "CSV" format
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

from typing import List, Any, Dict, Tuple

import json
import csv

from neoRNA import io
from neoRNA.structure.secondary_structure import SecondaryStructure

# ----------------------------------
# region Parse Arguments

def parse_args():

    # region Parser
    # -----------------------------------------
    # Parser
    # -----------------------------------------
    arguments_parser \
        = argparse.ArgumentParser(description='RNA Lib Prep - Convert RNA Lib Structure Summary File (in "JSON") to a CSV File')

    # Inputs
    arguments_parser.add_argument('--struct', dest='struct_summary',
                                  metavar='struct_summary.json',
                                  help='The file path to RNA Lib Structure Summary file, in "JSON" format.')
    # Output
    arguments_parser.add_argument('--out',
                                  action='store', default='structure.csv',
                                  help='Filename of RNA Lib Structure file, in CSV format.')

    # endregion

    # region Validate
    # -----------------------------------------
    # Validate
    # -----------------------------------------

    args = arguments_parser.parse_args()
    validate_args(args)

    # endregion

    #
    return args


#
def validate_args(args: Any):
    #
    struct_summary_file_path = os.path.abspath(args.struct_summary.strip())
    if not os.path.exists(struct_summary_file_path):
        raise ValueError('"RNA Lib Structure Summary" file does not exist. ')

# endregion


# ----------------------------------
# region Prep

def get_logger(logger_name: str, logging_level = logging.INFO):
    setup_logging(logging_level=logging_level)
    #
    return logging.getLogger(logger_name)


def prep(args: Any):
    #
    logger = get_logger('neo_rna.script')

    #
    struct_summary_file_path = os.path.abspath(args.struct_summary)

    #
    cwd = os.getcwd()
    output_file_path = args.out.strip()
    if not os.path.isabs(output_file_path):
        output_file_path = os.path.join(cwd, output_file_path)

    # return
    return logger, struct_summary_file_path, output_file_path

# endregion


# ----------------------------------
# region Load Data

def load(logger, config, struct_summary_file_path):
    #
    rna_lib_struct_summary_dict = {}

    with open(struct_summary_file_path) as infile:
        rna_lib_struct_summary_dict = json.load(infile)

    #
    return rna_lib_struct_summary_dict

# endregion


# ----------------------------------
# region Process

def retrieve_num_mut(mutation_syntax_str: str):
    #
    if mutation_syntax_str == 'wt':
        return 0
    if mutation_syntax_str == 'indel':
        return -1

    #
    parts = mutation_syntax_str.split(',')
    return len(parts)


def process(logger, config, rna_lib_struct_summary_dict):
    #
    output_entries = list()

    #
    for rna_item in rna_lib_struct_summary_dict['items']:
        #
        rna_id = rna_item['rna_id']
        logger.info('---------- RNA Item: {} ----------'.format(rna_id))

        sequence_string = rna_item['sequence_string']
        editing_value = rna_item['A-to-I_editing_level']
        editing_position = rna_item['A-to-I_editing_site']

        mutation_syntax_str = ''
        num_mut = -1
        if 'mutation_syntax' in rna_item:
            mutation_syntax_str = rna_item['mutation_syntax']
            num_mut = retrieve_num_mut(mutation_syntax_str)

        computational_structure_str = rna_item['computational_structure']
        # Get the structure annotation info
        computational_comment = ','.join(['>' + rna_id, 'computational'])
        computational_bprna_annotation, computational_secondary_structure = parse_bprna_annotation(
            '\n'.join([computational_comment, sequence_string, computational_structure_str]))
        computational_dot_bracket_annotation_str = computational_secondary_structure.dot_bracket_annotation

        experimental_structure_str = ''
        experimental_dot_bracket_annotation_str = ''
        if 'experimental_structure' in rna_item:
            experimental_structure_str = rna_item['experimental_structure']
            experimental_comment = ','.join(['>' + rna_id, 'computational'])
            experimental_bprna_annotation, experimental_secondary_structure = parse_bprna_annotation(
                '\n'.join([experimental_comment, sequence_string, computational_structure_str]))
            experimental_dot_bracket_annotation_str = experimental_secondary_structure.dot_bracket_annotation

        entry = list()
        entry.append(rna_id)
        entry.append(sequence_string)
        entry.append(str(editing_position))
        entry.append(str(editing_value))
        entry.append(mutation_syntax_str)
        entry.append(str(num_mut))
        entry.append(computational_structure_str)
        entry.append(computational_dot_bracket_annotation_str)
        entry.append(experimental_structure_str)
        entry.append(experimental_dot_bracket_annotation_str)

        output_entries.append(entry)

    #
    return output_entries

# endregion


# ----------------------------------
# region Output

def output(logger, config, output_file_path, output_entries):
    #
    headers = [
        'rna_id',
        'sequence_string',
        'A-to-I_editing_site',
        'A-to-I_editing_level',
        'mutation_syntax',
        'num_mutation',
        'computational_structure',
        'computational_dot_bracket_annotation',
        'experimental_structure',
        'experimental_dot_bracket_annotation',
    ]

    writer = csv.writer(open(output_file_path, 'w'))
    writer.writerow(headers)

    for entry in output_entries:
        writer.writerow(entry)

# endregion


# ----------------------------------
# region Utils

def parse_bprna_annotation(structure_content: str) -> Tuple[str, SecondaryStructure]:

    # use an intermediate tmp file because bpRNA needs it
    out_tmp = open('tmp.dbn', 'w')
    out_tmp.write(structure_content + '\n')
    out_tmp.close()

    #
    subprocess.call(["bpRNA.pl", "tmp.dbn"], stdout=subprocess.PIPE)

    #
    in_tmp = open("tmp.st", 'r')
    annotation = in_tmp.read()
    in_tmp.close()

    # Also parse the bpRNA contents
    with open("tmp.st") as handle:
        for secondary_structure in io.parse(handle, "bp-rna"):
            #
            return annotation, secondary_structure

# endregion


# ----------------------------------
# region Main Script

def main():
    #
    args = parse_args()

    # Global Config
    config = {}

    # Prep
    logger, struct_summary_file_path, output_file_path = prep(args)

    # Load original data
    rna_lib_struct_summary_dict = load(logger, config, struct_summary_file_path)

    # Process data
    output_entries = process(logger, config, rna_lib_struct_summary_dict)

    # Output
    output(logger, config, output_file_path, output_entries)


#
if __name__ == "__main__":
    main()

# endregion
