#!/usr/bin/env python

# ================================================================
# Python Script - Data Prep - Extract Editing Level Values from a Chromosone Dataset
#
# The original dataset includes two parts:
# - A "folder" hosts the "bpRNA" results of all "analyzed" chromosome DNAs
# - A "file" includes the "editing level" values.
#
# ## Input
# - bpRNA folder
# - Editing Level file
#
# ## Output
# - RNA Lib Structure Summary file, in "JSON" format
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

from neoRNA import io
from neoRNA.sequence.sequence import Sequence
from neoRNA.structure.secondary_structure import SecondaryStructure

from neoRNA.util.runner.matlab_runner import MatlabRunner
from neoRNA.util.file_utils import FileUtils


# ----------------------------------
# region Parse Arguments

def parse_args():

    # region Parser
    # -----------------------------------------
    # Parser
    # -----------------------------------------
    arguments_parser \
        = argparse.ArgumentParser(description='RNA Lib Prep - Extract Editing Level Values from a Chromosone Dataset')

    # Inputs
    arguments_parser.add_argument('--bprna', dest='bprna_folder',
                                  metavar='bprna_folder',
                                  help='The folder path to bpRNA folder.')
    arguments_parser.add_argument('--editing_level', dest='editing_level_file',
                                  metavar='editing_level_file',
                                  help='The file path to the Editing Level file')

    # Parameters
    arguments_parser.add_argument('--threshold',
                                  action='store', type=float, default=0.0,
                                  help='The "threshold" of editing level. Default - 0.0. ')

    # Output
    arguments_parser.add_argument('--out',
                                  action='store', default='chromosone.json',
                                  help='Filename of RNA Lib Structure file.')

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


def validate_args(args: Any):
    #
    bprna_folder_path = os.path.abspath(args.bprna_folder.strip())
    if not os.path.exists(bprna_folder_path):
        raise ValueError('"bpRNA" folder does not exist. ')
    editing_level_file_path = os.path.abspath(args.editing_level_file.strip())
    if not os.path.exists(editing_level_file_path):
        raise ValueError('"Editing Level" file does not exist. ')

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
    bprna_folder_path = os.path.abspath(args.bprna_folder)
    editing_level_file_path = os.path.abspath(args.editing_level_file)
    editing_level_threshold = args.threshold

    #
    cwd = os.getcwd()
    output_file_path = args.out.strip()
    if not os.path.isabs(output_file_path):
        output_file_path = os.path.join(cwd, output_file_path)

    # return
    return logger, bprna_folder_path, editing_level_file_path, editing_level_threshold, output_file_path

# endregion


# ----------------------------------
# region Load Data

def load(logger, config, bprna_folder_path, editing_level_file_path):
    #
    # Only load "bpRNA" file
    bprna_files = [f for f in os.listdir(bprna_folder_path)
                   if os.path.isfile(os.path.join(bprna_folder_path, f))
                   and os.path.splitext(f)[1] == '.st']

    # bpRNA structures
    bprna_structure_dict: Dict[str, SecondaryStructure] = dict()
    for file in bprna_files:
        #
        file_path = os.path.join(bprna_folder_path, file)

        with open(file_path) as handle:
            for secondary_structure in io.parse(handle, "bp-rna"):
                #
                chromosome_id = secondary_structure.comment
                bprna_structure_dict[chromosome_id] = secondary_structure
    #
    return bprna_files, bprna_structure_dict

# endregion


# ----------------------------------
# region Process

def process(logger, config, bprna_structure_dict, editing_level_file_path, editing_level_threshold):
    #
    output_entries = list()

    #
    bprna_id_list = bprna_structure_dict.keys()

    #
    with open(editing_level_file_path) as infile:
        reader = csv.DictReader(infile, dialect='excel-tab')
        for record in reader:
            #
            chromosome_id = '{}_{}'.format(record['#chrom'], record['position'])
            logger.info('---------- Chromosome Item: {} ----------'.format(chromosome_id))

            # Check if this "chromosome" has "bpRNA" results.
            # If "not", skip it
            if chromosome_id not in bprna_id_list:
                continue

            editing_level = float(record['editlevel']) if record['editlevel'] != 'N/A' else None
            logger.info('---------- Editing Level: {} ----------'.format(str(editing_level)))
            # Check if the "editing level" passes the "threshold"
            if not editing_level or editing_level < editing_level_threshold:
                continue

            # Other Info
            gene = record['gene']
            strand = record['strand']
            annotation = ' | '.join([record['annot1'], record['annot2']])
            total_reads = int(record['coverage'])
            edited_reads = int(record['editedreads'])

            #
            bprna_structure = bprna_structure_dict[chromosome_id]

            # Get the "editing position"
            sequence_str = str(bprna_structure.sequence.get_rna_sequence())
            editing_position = int((len(sequence_str) - 1) / 2)

            #
            entry = {
                'rna_id': chromosome_id,
                'gene': gene,
                'strand': strand,
                'annotation': annotation,
                'total_reads': total_reads,
                'edited_reads': edited_reads,
                #
                'sequence_string': sequence_str,
                'computational_structure': bprna_structure.dot_bracket,

                #
                'A-to-I_editing_level': editing_level,
                'A-to-I_editing_site': editing_position,
            }
            output_entries.append(entry)

    #
    return output_entries

# endregion


# ----------------------------------
# region Output

def output(logger, config, output_file_path, output_entries):
    #
    rna_library_structure_summary_dict = {
        "items": output_entries,
    }
    json_data = json.dumps(rna_library_structure_summary_dict,
                           sort_keys=True,
                           indent=4, ensure_ascii=True)

    # Output JSON as file
    FileUtils.save_file(output_file_path, json_data)

# endregion


# ----------------------------------
# region Main Script

def main():
    #
    args = parse_args()

    # Global Config
    config = {}

    # Prep
    logger, bprna_folder_path, editing_level_file_path, editing_level_threshold, output_file_path = prep(args)

    # Load original data
    bprna_files, bprna_structure_dict = load(logger, config, bprna_folder_path, editing_level_file_path)

    # Process data
    output_entries = process(logger, config, bprna_structure_dict, editing_level_file_path, editing_level_threshold)

    # Output
    output(logger, config, output_file_path, output_entries)


#
if __name__ == "__main__":
    main()

# endregion
