#!/usr/bin/env python

# ================================================================
# Python Script - Analysis - Isoform Ensemble Probability
#
# It loads the "RNA Lib Structure Summary" file (in JSON) and analysis the "Ensemble Probability" for each isoform.
#
# # Constraint Definition
# - . (no constraint for this base)
# - | (the corresponding base has to be paired
# - x (the base is unpaired)
# - < (base i is paired with a base j>i)
# - > (base i is paired with a base j<i)
# - () (base i pairs base j)
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "json" format.
# - WT "Constraint"
#
# ## Output
# - Plot (Editing Level vs Ensemble Probability)
#
#

import os
import argparse

# Add "py_scripts" into module path, relative to "current" script
import sys
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.realpath(__file__)))))
sys.path.append(local_module_path)

import json
import csv

import logging
from py_scripts import setup_logging

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

from neoRNA.util.runner.rnafold_runner import RnaFoldRunner

# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser = argparse.ArgumentParser(description='RNA Lib Analysis Script - Isoform Ensemble Probability.')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.json',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Parameters
arguments_parser.add_argument('--wt', dest='wt_sequence',
                              action="store",
                              help='"WT" sequence.')
arguments_parser.add_argument('--wt_structure', dest='wt_structure',
                              action="store",
                              help='"WT" secondary structure.')
arguments_parser.add_argument('--wt_pair', dest='wt_pair_structure',
                              action="store",
                              help='"WT" "Pair" structure.')
arguments_parser.add_argument('--wt_constraint',
                              action="store",
                              help='"WT" constraint string.')

# Output
arguments_parser.add_argument('--out', dest='out',
                              action='store', default='output.png',
                              help='Filename / path of output image file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)


# Validation
if not os.path.exists(rna_lib_struct_summary_file_path.strip()):
    raise ValueError('"RNA Lib Structure Summary" file does not exist. ')

#
wt_sequence = args.wt_sequence.strip()
wt_structure = args.wt_structure.strip()
wt_pair_structure = args.wt_pair_structure.strip()
wt_constraint = args.wt_constraint.strip()

if not wt_constraint:
    raise ValueError('"WT" constraint required. ')
if not wt_sequence:
    raise ValueError('"WT" sequence required. ')

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
# region Extract FE and Calc Prob

# Canonical BPs
canonical_base_pairs = [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')]


# Decode the "WT Constraint" string to find "Soft Pairs".
def decode_soft_pair(constraint_str: str) -> list:
    r"""
    Decode the "WT Constraint" string to find "Soft Pairs".

    Parameters
    ----------
    constraint_str: str

    Returns
    -------

    """

    #
    left_bracket_pos_list = list()
    soft_pair_list = list()

    #
    left_bracket_char = '('
    right_bracket_char = ')'

    #
    constraint_char_list = list(constraint_str)
    for idx in range(len(constraint_char_list)):
        #
        if constraint_char_list[idx] == left_bracket_char:
            left_bracket_pos_list.append(idx)
        elif constraint_char_list[idx] == right_bracket_char:
            # Find a "PREVIOUS" left to be a pair
            if len(left_bracket_pos_list) > 0:
                left_pos = left_bracket_pos_list.pop()
                soft_pair_list.append((left_pos, idx))

    #
    return soft_pair_list


rna_fold_runner = RnaFoldRunner()

editing_value_list = list()
fe_probability_list = list()

csv_data_list = list()

# No "Non-canonical"
editing_value_all_canonical_list = list()
fe_probability_all_canonical_list = list()

# Decode the "WT pair" structure
wt_defined_bp_list = decode_soft_pair(wt_pair_structure)
logger.info('---- WT Defined Pairs: {} '.format(wt_defined_bp_list))

# Loop
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']
    # logger.info('---------- RNA Item: {} ----------'.format(rna_id))

    #
    sequence_string = rna_item['sequence_string']
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    # Skip the RNA item if
    # - the "sequence length" is NOT the SAME as WT
    # - Editing value is not available
    if len(sequence_string) != len(wt_sequence) or not editing_value:
        continue

    #
    csv_rna_item = list()
    csv_rna_item.append(rna_id)
    csv_rna_item.append(float(editing_value))
    editing_value_list.append(float(editing_value))

    # FInd Non-canonical BPs within "All" soft pairs
    non_canonical_pair_positions = []
    for soft_bp in wt_defined_bp_list:
        #
        bp = (sequence_string[soft_bp[0]], sequence_string[soft_bp[1]])
        if bp not in canonical_base_pairs:
            #
            non_canonical_pair_positions.append(soft_bp[0])
            non_canonical_pair_positions.append(soft_bp[1])

    # Extract MFE & Ensemble FE
    mfe_wt_constraint, structure_wt_constraint = rna_fold_runner.extract_free_energy(sequence_string, wt_constraint)
    fe_ensemble, structure_ensemble = rna_fold_runner.extract_free_energy(sequence_string)

    # Probability
    probability = np.exp(-1 * mfe_wt_constraint / 0.6) / np.exp(-1 * fe_ensemble / 0.6)
    # Apply a threshold
    # probability = threshold(probability)

    # logger.info('---- Ensemble-FE: {} | MFE: {} | prob: {}'.format(str(fe_ensemble), str(mfe_wt_constraint),
    #                                                                str(probability)))

    # Apply "Non-canonical BP" penalty
    num_non_canonical = len(set(non_canonical_pair_positions))
    if num_non_canonical > 0:
        #
        probability /= num_non_canonical
        # probability /= (0.5 * num_non_canonical)
        # probability = max( 0.0, probability - num_non_canonical * k )
        # probability -= num_non_canonical * 0.01
    else:
        # Add it to a special list for "NO Non-canonical BP"
        editing_value_all_canonical_list.append(float(editing_value))
        fe_probability_all_canonical_list.append(probability)

    #
    fe_probability_list.append(probability)
    csv_rna_item.append(probability)

    #
    csv_data_list.append(csv_rna_item)

#
rna_fold_runner.cleanup()

# endregion


# ----------------------------------
# region Plot

# Examine the probability size
logger.info('========= Examine Probability ========== \n---- MAX: {} \n---- MIN: {}'
            .format(str(np.max(fe_probability_list)), str(np.min(fe_probability_list))))

corr, p_value = pearsonr(fe_probability_list, editing_value_list)
logger.info('---- Correlation coefficient: {} | P value: {}'.format(str(corr), str(p_value)))

#
plt.scatter(fe_probability_list, editing_value_list,
            c='b', s=45.0, lw=0, alpha=0.5)
plt.xlabel('Probability of WT Secondary Structure')
plt.ylabel('Editing Level')
plt.axis([-0.03, 0.4, -0.1, 1.05])

# Text
plt.text(-0.0, 0.99, wt_structure, family='monospace', fontsize=10)
plt.text(-0.0, 0.95, wt_constraint, family='monospace', fontsize=10)
plt.text(-0.0, 0.91, 'Pearson Corr = {}, P_value = {}'.format(str(corr), str(p_value)),
         family='monospace', fontsize=8)

# Save
fig = plt.gcf()
fig.set_size_inches(18.5, 10.5)
fig.savefig(output_file_path, dpi=100)
# plt.show()

# endregion

# ----------------------------------
# region Output Plot Data

with open('{}.csv'.format(output_file_path), 'w') as outfile:
    # Head Line
    head_line = [
        'rna_id',
        'editing_value',
        'probability',
    ]
    outfile.write(','.join(head_line) + '\n')
    csv_writer = csv.writer(outfile)

    for item in csv_data_list:
        # outfile.write(','.join(data_line))
        csv_writer.writerow(item)

# endregion
