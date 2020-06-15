#!/usr/bin/env python

# ================================================================
# "RNA Lib" Analysis Script - QC "Mutation Mapping" in RNA Lib Items
#
# Extra Info:
#
#   - Total Read Depth for "Modified" / "Untreated" reads.
#       - It loads the "results" files from "cutadapt" process.
#
# ## Input
#
# - "RNA Lib" profiling file
#     - The "python object" file to contain all "RNA Lib" item objects which include the "profiling info".
#
# ## Output
#
# - Updated "RNA Lib" profiling filev
#

import os
import argparse

import csv

# Add "py_scripts" into module path, relative to "current" script
import sys
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.realpath(__file__)))))
sys.path.append(local_module_path)

import logging
from py_scripts import setup_logging

from neoRNA.io.library_io import LibraryIO
from neoRNA.library.rna_library import RnaLibrary


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Util - Add "extra" info to profiling')

# Inputs
arguments_parser.add_argument('--rna_lib_profiling',
                              metavar='rna_lib_profiling.rbin',
                              help='The file path to RNA Lib profiling result file.')

# Output
arguments_parser.add_argument('--out',
                              action='store', default='qc_mut_mapping.csv',
                              help='The output filename.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_profiling_file_path = os.path.abspath(args.rna_lib_profiling)

# Validation
if not os.path.exists(rna_lib_profiling_file_path):
    raise ValueError('"RNA Lib profiling" file does not exist.')

# Values

# Output
output_filename = args.out

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
# region Prep - Path

#
cwd = os.getcwd()

# endregion


# ----------------------------------
# region Prep - RNA Lib Object

rna_library: RnaLibrary = LibraryIO.as_python_object(rna_lib_profiling_file_path)

# endregion


# ----------------------------------
# region Analyze

#
results = []

#
wt_id = rna_library.wide_type_rna_id
wt_rna_sequence = rna_library.wide_type_rna_sequence

for rna_item in rna_library.rna_items:
    #
    rna_id = rna_item.rna_id
    barcode = rna_item.barcode.barcode
    modified_depth = rna_item.modified_read_depth

    #
    if rna_id == wt_id:
        continue

    #
    mut_positions, mut_syntax_list = rna_item.sequence.generate_mutation_syntax(wt_rna_sequence,
                                                                                seq_type_rna=True, trim_right=13)

    # Ignore if no "mutation"
    if not mut_positions:
        continue

    for index, position in enumerate(mut_positions):
        #
        profile = rna_item.profile_dict[str(position)]
        nt_sequence = profile.nt_sequence
        modified_position_depth = profile.modified_read_depth
        modified_position_effective_depth = profile.modified_effective_depth

        #
        item = []

        item.append(rna_id)
        item.append(barcode)

        item.append(position)
        item.append(nt_sequence)
        item.append(mut_syntax_list[index])

        item.append(modified_position_depth)
        item.append(modified_position_effective_depth)
        item.append(modified_depth)
        item.append(modified_position_depth / modified_depth)

        results.append(item)

# endregion


# ----------------------------------
# region Output

output_file_path = os.path.join(cwd, output_filename)
with open(output_file_path, 'w') as outfile:
    # Head Line
    head_line = [
        'RNA_ID',
        'Barcode',
        'Mut Position',
        'Mut Sequence',
        'Mut Syntax',
        'Mut Depth',
        'Mut Effective Depth',
        'Total Read Depth',
        'Mapping Rate'
    ]
    outfile.write(','.join(head_line) + '\n')
    csv_writer = csv.writer(outfile)

    for result in results:
        # outfile.write(','.join(data_line))
        csv_writer.writerow(result)

# endregion


