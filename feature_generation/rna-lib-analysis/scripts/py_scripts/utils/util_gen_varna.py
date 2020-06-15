#!/usr/bin/env python

# ================================================================
# "RNA Lib" Util Script - Generate VARNA RNA Structure Image
#
# It loads the "RNA Lib Structure Summary" file (in JSON) and call "VARNA" to output the RNA structure image.
#
# NOTE:
# - If the "RNA Lib Structure Summary" file also include the "experimental" structure, also generate image for that.
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "json" format.
#
#

import os
import argparse

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

from neoRNA.util.runner.verna_runner import VarnaRunner

# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Util - Generate VARNA RNA Structure Image')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.json',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Params
arguments_parser.add_argument('--varna', dest='varna_location',
                              action="store",
                              help='The "path" to VARNA jar package.')


# Output
arguments_parser.add_argument('--out_folder',
                              action='store',
                              help='The output folder.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)

# Validation
if not os.path.exists(rna_lib_struct_summary_file_path):
    raise ValueError('"RNA Lib" definition file does not exist.')

#
varna_location = args.varna_location.strip()

# Output
output_folder_path = args.out_folder

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

if not os.path.isabs(output_folder_path):
    cwd = os.getcwd()
    output_folder_path = os.path.join(cwd, output_folder_path)

if not os.path.exists(output_folder_path):
    os.mkdir(output_folder_path)

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

def decode_mutation_syntax_str(mutation_syntax_str):
    r"""
    Decode "mutation syntax" string to get a list of "Mutation Position".

    Parameters
    ----------
    mutation_syntax_str

    Returns
    -------

    """

    if not mutation_syntax_str:
        return []

    if mutation_syntax_str.lower().startswith('indel') or mutation_syntax_str.lower().startswith('wt'):
        return []

    mutation_syntax_list = mutation_syntax_str.split(',')
    mutation_positions = list()
    for mutation_string in mutation_syntax_list:
        mutation_position = mutation_string[:-4]
        mutation_positions.append(mutation_position)

    return mutation_positions

# VARNA runner
varna = VarnaRunner(varna_location)

for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']
    logger.info('---------- RNA Item: {} ----------'.format(rna_id))

    sequence_string = rna_item['sequence_string']
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    # Check if there is `mutation_syntax`
    mutation_syntax_str = None
    mutation_positions = []
    if 'mutation_syntax' in rna_item:
        mutation_syntax_str = rna_item['mutation_syntax']
        mutation_positions = decode_mutation_syntax_str(mutation_syntax_str)

    #
    rna_id_str = '{}-{}'.format(rna_id, mutation_syntax_str) if mutation_syntax_str is not None \
        else '{}'.format(rna_id)

    # Generate Annotations & VARNA highlight regions
    highlight_regions = list()
    annotations = list()

    # Add Editing Position
    highlight_regions.append({
        'nt_range': '{}-{}'.format(str(editing_position), str(editing_position)),
        'fill': '#0B4F6C'
    })
    annotations.append({
        'annotation_str': 'E-{}'.format(str(editing_position)),
        'anchor': '{}'.format(str(editing_position)),
    })
    for mutation_position in mutation_positions:
        #
        highlight_region_options = {
            'nt_range': '{}-{}'.format(str(mutation_position), str(mutation_position)),
            'fill': '#C55337'    # Mut
        }
        highlight_regions.append(highlight_region_options)

        #
        annotation_options = {
            'annotation_str': 'M-{}'.format(str(mutation_position)),
            'anchor': '{}'.format(str(mutation_position)),
        }
        annotations.append(annotation_options)

    # VARNA options
    options = {
        'baseNum': '#000000',
        'highlightRegion': varna.gen_highlight_region_str(highlight_regions),
        'annotations': varna.gen_annotation_str(annotations),

        #
        'title': rna_id_str
    }

    # "computational"
    computational_structure_string = rna_item['computational_structure']
    computational_output_folder = os.path.join(output_folder_path, 'computational')
    if not os.path.exists(computational_output_folder):
        os.mkdir(computational_output_folder)
    varna.gen_struct_image(sequence_string, computational_structure_string,
                           out_file=os.path.join(computational_output_folder, '{}.png'.format(rna_id_str)),
                           options=options)

    # "experimental", Optional
    if 'experimental_structure' in rna_item:
        experimental_structure_string = rna_item['experimental_structure']
        experimental_output_folder = os.path.join(output_folder_path, 'experimental')
        if not os.path.exists(experimental_output_folder):
            os.mkdir(experimental_output_folder)
        varna.gen_struct_image(sequence_string, experimental_structure_string,
                               out_file=os.path.join(experimental_output_folder, '{}.png'.format(rna_id_str)),
                               options=options)


# endregion

