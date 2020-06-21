#!/usr/bin/env python

# ================================================================
# RNA Lib Analysis Pipeline - Quality Control on "Barcode"
#
# It does "qc" on the "barcode" sequence to help evaluate the "reads" from the experiment.
#
# NOTES:
#
#
#

import os
import sys
import subprocess
import argparse
import logging

# Add "py_scripts" into module path, relative to "current" script
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.realpath(__file__)))))
sys.path.append(local_module_path)

import csv

from py_scripts import setup_logging

from typing import List

from neoRNA import io
from neoRNA.library.library_item import LibraryItem

from neoRNA.sequence.sequence import Sequence


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Analysis - Quality Control on "Barcode"')

# Inputs
arguments_parser.add_argument('--rna_lib_def',
                              metavar='rna_lib_def_file',
                              help='The file path to RNA Lib definition file.')

#
arguments_parser.add_argument('--dataset_folder',
                              action='store',
                              help='The "dataset" folder.')
arguments_parser.add_argument('--read_file',
                              action='store',
                              help='The "read" file to check.')
arguments_parser.add_argument('--prefix_count',
                              action='store',
                              help='The "number" of nt "before" the barcode.')
arguments_parser.add_argument('--extra_nt',
                              action='store', default='C',
                              help='The "extra" nt which is behind the barcode')

# Output
arguments_parser.add_argument('--out',
                              action='store', default='qc_barcode.csv',
                              help='The output file.')

# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_def_file_path = os.path.abspath(args.rna_lib_def)
dataset_folder_path = os.path.abspath(args.dataset_folder)

# Validation
if not os.path.exists(rna_lib_def_file_path):
    raise ValueError('"RNA Lib" definition file does not exist.')
if not os.path.exists(dataset_folder_path):
    raise ValueError('"Dataset" folder does not exist.')

#
read_file_name = args.read_file
prefix = args.prefix_count
extra_nt = args.extra_nt

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
    output_file_path = os.path.join(cwd, output_file_path)

# endregion


# ----------------------------------------------------------------
# region RNA Lib Items
#

rna_lib_items: List[LibraryItem] = []
for rna_item in io.parse(rna_lib_def_file_path, "rna-lib-def"):
    rna_lib_items.append(rna_item)

# endregion


# ----------------------------------
# region QC

qc_results = []

#
os.chdir(dataset_folder_path)

# barcode = 'CGCGGTTGT'
# reverse = 'ACAACCGCG'
# cmd = 'cat {} | grep -E -c "^[ATCG]{{{}}}{}[ATCG]{{4}}{}"'.format(read_file_name, prefix, barcode, reverse)
# logger.info('----- CMD: {}'.format(cmd))
# try:
#     output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
# except:
#     output = 0

#
for rna_lib_item in rna_lib_items:
    #
    rna_id = rna_lib_item.rna_id
    barcode = rna_lib_item.barcode.barcode

    barcode_sequence = Sequence(barcode)
    reverse = barcode_sequence.get_reverse_complement()

    logger.info('----- RNA ID: {} | Barcode: {}'.format(rna_id, barcode))

    result = []
    result.append(rna_id)
    result.append(barcode)

    # Barcode
    cmd = 'cat {} | grep -E -c "^[ATCG]{{{}}}{}"'.format(read_file_name, prefix, barcode)
    logger.info('----- CMD: {}'.format(cmd))
    try:
        output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
    except:
        output = 0
    result.append(output)

    # Barcode + Reverse
    cmd = 'cat {} | grep -E -c "^[ATCG]{{{}}}{}[ATCG]{{4}}{}"'.format(read_file_name, prefix, barcode, str(reverse))
    logger.info('----- CMD: {}'.format(cmd))
    try:
        output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
    except:
        output = 0
    result.append(output)

    # Barcode (plus 1 extra nt)
    cmd = 'cat {} | grep -E -c "^[ATCG]{{{}}}{}{}"'.format(read_file_name, prefix, barcode, extra_nt)
    logger.info('----- CMD: {}'.format(cmd))
    try:
        output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
    except:
        output = 0
    result.append(output)

    # Barcode (plus 1 extra nt) + Reverse
    cmd = 'cat {} | grep -E -c "^[ATCG]{{{}}}{}{}[ATCG]{{3}}{}"'.format(read_file_name, prefix, barcode, extra_nt, str(reverse))
    logger.info('----- CMD: {}'.format(cmd))
    try:
        output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
    except:
        output = 0
    result.append(output)

    #
    qc_results.append(result)


# endregion


# ----------------------------------
# region Output

with open(output_file_path, 'w') as outfile:
    # Head Line
    head_line = [
        'RNA_ID',
        'Barcode',
        'Matched - Barcode',
        'Matched - Barcode + Revers',
        'Matched - Barcode(extra 1 nt)',
        'Matched - Barcode(extra 1 nt) + Reverse'
    ]
    outfile.write(','.join(head_line) + '\n')
    csv_writer = csv.writer(outfile)

    for result in qc_results:
        # outfile.write(','.join(data_line))
        csv_writer.writerow(result)

# endregion
