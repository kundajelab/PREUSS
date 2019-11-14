#!/usr/bin/env python

# ================================================================
# Python Script - Analysis - Generate ML Features
#
# It loads the "RNA Lib Structure Summary" file (in JSON) as well as the "bpRNA" info and generate the
#   ML features based on given platform.
#
# ## Input
# - RNA Lib Structure Summary file
#   - It is in "json" format.
# - Data Type
# - WT sequence
#
# ## Output
# - The feature output is in `csv` format.
# - Optional. If need to output the bpRNA results.
#   - It is in "txt" format
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

from neoRNA import io
from neoRNA.sequence.sequence import Sequence
from neoRNA.structure import SecondaryStructureElementType
from neoRNA.structure.secondary_structure import SecondaryStructure
from neoRNA.analysis.editing_analysis import EditingAnalysis
from neoRNA.analysis.editing_analysis_item import EditingAnalysisItem
from neoRNA.util.file_utils import FileUtils

from neoRNA.util.runner.rnafold_runner import RnaFoldRunner
from neoRNA.util.runner.simtree_runner import SimTreeRunner


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser = argparse.ArgumentParser(description='RNA Lib Analysis Script - Generate ML Features.')

# Inputs
arguments_parser.add_argument('rna_lib_struct_summary_file',
                              metavar='rna_lib_struct_summary_file.json',
                              help='The file path to the "RNA Lib Structure Summary" file.')

# Parameters
arguments_parser.add_argument('--data_type',
                              action="store", default='computational',
                              help='Data type - computational | experimental. Default: "computational".')
arguments_parser.add_argument('--wt', dest='wt_sequence',
                              action="store", default='',
                              help='"WT" sequence.')

# Output
arguments_parser.add_argument('--out', dest='out',
                              action='store', default='output.csv',
                              help='Filename / path of features output file.')
arguments_parser.add_argument('--out_bprna',
                              action='store', default=None,
                              help='Filename of output file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
rna_lib_struct_summary_file_path = os.path.abspath(args.rna_lib_struct_summary_file)


# Validation
if not os.path.exists(rna_lib_struct_summary_file_path.strip()):
    raise ValueError('"RNA Lib Structure Summary" file does not exist. ')

#
data_type = args.data_type.strip().lower()
wt_sequence = args.wt_sequence.strip()

# if not wt_sequence:
#     raise ValueError('"WT" sequence required. ')

output_file_path = args.out.strip()
bprna_output_file_path = args.out_bprna

# endregion

# ----------------------------------
# region Global Config

# Number of "structure elements" to add for both "upstream" and "downstream".
num_upstream_elements = 3
num_downstream_elements = 3

# Defaults value for "EMPTY" values
empty_value_str = None
empty_value_numeric = None

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

if bprna_output_file_path and not os.path.isabs(output_file_path):
    bprna_output_file_path = os.path.join(cwd, bprna_output_file_path)

# endregion


# ----------------------------------
# region Load "RNA Lib Structure Summary" Data

# RSample results, indexed by "RNA ID"
rna_lib_struct_summary_dict = {}

with open(rna_lib_struct_summary_file_path) as infile:
    rna_lib_struct_summary_dict = json.load(infile)

# endregion


# ----------------------------------
# region Process RNA Lib Items

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


def extract_analysis_item_features(analysis_item: EditingAnalysisItem = None,
                                   represent_sequence: Sequence = None) -> List[Any]:
    r"""
    Extract ML features from an "analysis item".

    Current features for each item:
    - struct - element type
    - length
    - length_stem
    - length_hairpin
    - length_bulge
    - length_interior_es
    - length_interior_ecs
    - 5prm_cp_hairpin
    - 5prm_cp_bulge
    - 3prm_cp_bulge
    - 5prm_cp_interior
    - 3prm_cp_interior

    Parameters
    ----------
    analysis_item: EditingAnalysisItem
    represent_sequence: Sequence

    Returns
    -------

    """
    # Defaults
    struct = analysis_item.element_type if analysis_item else empty_value_str
    length = empty_value_numeric
    length_stem = empty_value_numeric
    length_hairpin = empty_value_numeric
    length_bulge = empty_value_numeric
    length_interior_es = empty_value_numeric
    length_interior_ecs = empty_value_numeric
    five_pprm_cp_hairpin = empty_value_str
    five_pprm_cp_bulge = empty_value_str
    three_pprm_cp_bulge = empty_value_str
    five_pprm_cp_interior = empty_value_str
    three_pprm_cp_interior = empty_value_str

    if analysis_item:

        # Strand Length
        if represent_sequence:
            length = represent_sequence.length

        #
        if struct == SecondaryStructureElementType.Stem:
            length_stem = length
        elif struct == SecondaryStructureElementType.Hairpin:
            length_hairpin = length
            five_pprm_cp_hairpin = analysis_item.closing_base_pair()
        elif struct == SecondaryStructureElementType.Bulge:
            length_bulge = length
            five_pprm_cp_bulge = analysis_item.closing_base_pair()
            three_pprm_cp_bulge = analysis_item.closing_base_pair(five_prim=False)
        elif struct == SecondaryStructureElementType.Interior:
            length_interior_es = length
            five_pprm_cp_interior = analysis_item.closing_base_pair()

            # Try to get "complementary strand" sequence
            represent_sequence_ecs = analysis_item.complementary_strand_sequence(represent_sequence)
            length_interior_ecs = represent_sequence_ecs.length if represent_sequence_ecs else 0
            three_pprm_cp_interior = analysis_item.closing_base_pair(five_prim=False)

    #
    return [struct, length,
            length_stem, length_hairpin, length_bulge, length_interior_es, length_interior_ecs,
            five_pprm_cp_hairpin, five_pprm_cp_bulge, three_pprm_cp_bulge, five_pprm_cp_interior, three_pprm_cp_interior]


def extract_bprna_features(sequence: Sequence, secondary_structure: SecondaryStructure, editing_position: int) -> List[Any]:
    r"""
    Extract ML features from bpRNA structure.

    Current Features:
    - u_count
    - u_stem_length
    - u_hairpin_length
    - d_count
    - "editing site" structure element
    - "upstream" structure elements
    - "downstream" structure elements

    Parameters
    ----------
    sequence: Sequence
    secondary_structure: SecondaryStructure
    editing_position: int

    Returns
    -------

    """
    #
    # info_dict = CommentParser.parse_bp_rna(secondary_structure.comment)

    # Defaults
    all_stem_length = empty_value_numeric
    u_count = 0
    u_all_stem_length = empty_value_numeric
    u_hairpin_length = empty_value_numeric
    d_count = 0
    d_all_stem_length = empty_value_numeric

    site_structure_features = list()
    upstream_structure_feature_list = list()
    downstream_structure_feature_list = list()

    # Analyze
    editing_analysis = EditingAnalysis(secondary_structure, editing_position)
    editing_analysis.analysis()

    # Separate items, by "distance"
    editing_site_items = list()
    upstream_items = list()
    downstream_items = list()
    analysis_items_by_distance = editing_analysis.analysis_items_by_distance

    # Loop from "smallest", stack the "positive" items  - the "last" one is the FIRST "up stream" item
    for distance in sorted(analysis_items_by_distance.keys()):
        item, distance_sequence = analysis_items_by_distance[distance]
        if distance < 0:
            upstream_items.append((item, distance, distance_sequence))
        elif distance == 0:
            editing_site_items.append((item, distance_sequence))
        else:
            break
    # Loop from "largest", stach the "positive" items - the "last" one is the FIRST "down stream" item
    for distance in sorted(analysis_items_by_distance.keys(), reverse=True):
        item, distance_sequence = analysis_items_by_distance[distance]
        if distance > 0:
            downstream_items.append((item, distance, distance_sequence))
        else:
            break

    # ----------------------------------
    # region "Editing Site" element
    #

    # site_feature = secondary_structure.get_annotation(editing_position)
    # editing_site_items = editing_analysis.get_analysis_items([EditingAnalysisItemType.Contain])
    if len(editing_site_items) > 0:
        # Only get the "first" one as "editing level" item
        site_item, site_included_sequence = editing_site_items[0]
        site_structure_features = extract_analysis_item_features(site_item, site_included_sequence)

        # Special features
        site_prev_nt = sequence.get_nt(editing_position - 1)
        site_next_nt = sequence.get_nt(editing_position + 1)
        site_prev_struct = secondary_structure.get_annotation(site_included_sequence.start_position - 1)
        site_next_struct = secondary_structure.get_annotation(site_included_sequence.end_position + 1)

        # site_1_1
        site_1_1 = empty_value_str
        if site_item.element_type == SecondaryStructureElementType.Interior:
            # Get the "ECS" sequence
            site_included_sequence_ecs = site_item.complementary_strand_sequence(site_included_sequence)
            # Only get "1:1 Interior Loop"
            if site_included_sequence.length == 1 \
                    and site_included_sequence_ecs and site_included_sequence_ecs.length == 1:
                # Keep the record
                site_1_1 = '{}:{}'.format(site_included_sequence.sequence_str,
                                                   site_included_sequence_ecs.sequence_str)
        # Stem
        if site_item.element_type == SecondaryStructureElementType.Stem:
            #
            site_nt = sequence.get_nt(editing_position)
            site_ecs_position, site_ecs_nt = site_item.complementary_nt(site_included_sequence, editing_position)
            if site_ecs_nt:
                # Keep the record
                site_1_1 = '{}:{}'.format(site_nt, site_ecs_nt)

        # Insert above features into a proper place - "index -> 1"
        site_structure_features.insert(1, site_1_1)
        site_structure_features.insert(1, site_next_struct)
        site_structure_features.insert(1, site_prev_struct)
        site_structure_features.insert(1, site_next_nt)
        site_structure_features.insert(1, site_prev_nt)

        # Stem length
        if site_item.element_type == SecondaryStructureElementType.Stem:
            if site_item and site_item.loop_length() is not None:
                if all_stem_length is None:
                    all_stem_length = 0
                all_stem_length += site_item.loop_length()

    # endregion

    # ----------------------------------
    # region Upstream Items

    #
    upstreams_hairpin_found = False
    upstreams_found = 0
    while upstream_items:
        u_item, u_distance, u_distance_sequence = upstream_items.pop()
        u_structure_features = extract_analysis_item_features(u_item, u_distance_sequence)

        #
        u_structure_features.insert(0, u_distance)
        u_structure_features.insert(0, 1)  # set "exist" to "1"

        # Use "insert" so that later we can "pop" it
        upstream_structure_feature_list.insert(0, u_structure_features)

        #
        if u_item.element_type == SecondaryStructureElementType.Stem:
            if u_item and u_item.loop_length() is not None:
                if u_all_stem_length is None:
                    u_all_stem_length = 0
                u_all_stem_length += u_item.loop_length()
        if not upstreams_hairpin_found and u_item.element_type == SecondaryStructureElementType.Hairpin:
            u_hairpin_length = u_item.loop_length()
            upstreams_hairpin_found = True

        #
        upstreams_found += 1

    u_count = upstreams_found
    # Stem Length
    if u_all_stem_length is not None:
        if all_stem_length is None:
            all_stem_length = 0
        all_stem_length += u_all_stem_length

    # endregion

    # ----------------------------------
    # region Downstream

    #
    downstreams_found = 0
    while downstream_items:
        d_item, d_distance, d_distance_sequence = downstream_items.pop()
        d_structure_features = extract_analysis_item_features(d_item, d_distance_sequence)

        #
        d_structure_features.insert(0, d_distance)
        d_structure_features.insert(0, 1)  # set "exist" to "1"

        #
        downstream_structure_feature_list.insert(0, d_structure_features)

        #
        if d_item.element_type == SecondaryStructureElementType.Stem:
            if d_item and d_item.loop_length() is not None:
                if d_all_stem_length is None:
                    d_all_stem_length = 0
                d_all_stem_length += d_item.loop_length()

        #
        downstreams_found += 1
    d_count = downstreams_found
    # Stem Length
    if d_all_stem_length is not None:
        if all_stem_length is None:
            all_stem_length = 0
        all_stem_length += d_all_stem_length

    # endregion

    # Build features
    bprna_features = list()

    #
    bprna_features.append(all_stem_length)

    # site features
    bprna_features += site_structure_features

    # upstream features
    bprna_features.append(u_count)
    bprna_features.append(u_all_stem_length)
    bprna_features.append(u_hairpin_length)
    # Structure element features
    for index in range(0, num_upstream_elements):
        if len(upstream_structure_feature_list) > 0:
            u_features = upstream_structure_feature_list.pop()
            bprna_features += u_features
        else:
            # Insert "empty" element
            empty_u_features = extract_analysis_item_features()
            empty_u_features.insert(0, 0)  # set "distance" to "0"
            empty_u_features.insert(0, 0)  # set "exist" to "0"

            bprna_features += empty_u_features

    # downstream features
    bprna_features.append(d_count)
    bprna_features.append(d_all_stem_length)
    # Structure element features
    for index in range(0, num_downstream_elements):
        if len(downstream_structure_feature_list) > 0:
            d_features = downstream_structure_feature_list.pop()
            bprna_features += d_features
        else:
            # Insert "empty" element
            empty_d_features = extract_analysis_item_features()
            empty_d_features.insert(0, 0)   # set "distance" to "0"
            empty_d_features.insert(0, 0)   # set "exist" to "0"

            bprna_features += empty_d_features

    return bprna_features


# Dict to keep the bpRNA content (text-based)
# rna_lib_bprna_dict = dict()
rna_lib_bprna_dict = defaultdict(Any)
rna_lib_secondary_structure_dict = defaultdict(Any)

# WT structure
wt_secondary_structure = None

# Prepare the "secondary structure" for all RNA items, also identify "WT"
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']

    sequence_string = rna_item['sequence_string']
    sequence = Sequence(sequence_string)
    # Either "computational_structure" or "experimental_structure"
    structure_string = rna_item['_'.join([data_type, 'structure'])]

    #
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    #
    rna_lib_bprna_dict[rna_id] = dict()

    #
    comment = ','.join(['>' + rna_id, data_type])
    bprna_annotation, secondary_structure = parse_bprna_annotation('\n'.join([comment, sequence_string, structure_string]))
    rna_lib_bprna_dict[rna_id][data_type] = bprna_annotation
    rna_lib_secondary_structure_dict[rna_id] = secondary_structure

    # Check if it is "WT"
    if 'mutation_syntax' in rna_item \
            and rna_item['mutation_syntax'].strip().lower().startswith('wt'):
        #
        wt_secondary_structure = secondary_structure

    # Consider "bootstrap" structures
    # if "bootstrap_structures" in rna_item:
    #     #
    #     logger.info('  -- Process "bootstrap" data -->')
    #
    #     #
    #     bootstrap_structures = rna_item['bootstrap_structures']
    #     bootstrap_counts = dict()
    #     for key, group in groupby(bootstrap_structures):
    #         #
    #         bootstrap_counts[key] = len(list(group))
    #
    #     #
    #     bootstrap_index = 0
    #     rna_lib_bprna_dict[rna_id]['bootstraps'] = dict()
    #     for bootstrap_structure, bootstrap_count in bootstrap_counts.items():
    #         #
    #         comment = ','.join(['>' + rna_id, 'bootstrap-' + str(bootstrap_index), 'frequency-' + str(bootstrap_count)])
    #         bprna_annotation, secondary_structure = parse_bprna_annotation('\n'.join([comment, sequence_string, bootstrap_structure]))
    #         rna_lib_bprna_dict[rna_id]['bootstraps'][bootstrap_structure] = bprna_annotation
    #         bootstrap_index += 1


# Construct feature list for each RNA item.
rna_lib_features_list = list()

# The runner for getting the "free energy"
rna_fold_runner = RnaFoldRunner()

# The runner for SimTree
simtree_jar_location = '/Users/cowfox/Desktop/_rna_lib_analysis/__tool/SimTree_v1.2.3.jar'
simtree_runner = SimTreeRunner(simtree_jar_location, os.getcwd())

# os.removedirs()

#
for rna_item in rna_lib_struct_summary_dict['items']:
    #
    rna_id = rna_item['rna_id']
    logger.info('---------- RNA Item: {} ----------'.format(rna_id))

    # if rna_id != '076':
    #     continue

    sequence_string = rna_item['sequence_string']
    sequence = Sequence(sequence_string)
    # Either "computational_structure" or "experimental_structure"
    structure_string = rna_item['_'.join([data_type, 'structure'])]

    #
    editing_value = rna_item['A-to-I_editing_level']
    editing_position = rna_item['A-to-I_editing_site']

    # Free Energy
    fe_ensemble, structure_ensemble = rna_fold_runner.extract_free_energy(sequence_string)

    # RNA Structure Similarity against WT
    simtree_normalized_score = empty_value_numeric
    simtree_flipping_modes = empty_value_str
    if wt_secondary_structure:
        simtree_normalized_score, simtree_flipping_modes \
            = simtree_runner.compare_rna_structure(structure_string, wt_secondary_structure.dot_bracket)

    #
    secondary_structure = rna_lib_secondary_structure_dict[rna_id]

    # Parse the ML features
    structure_features = extract_bprna_features(sequence, secondary_structure, editing_position)
    site_struct = structure_features[0]  # The first one is for "editing site"

    # Check Other Features - "mutation" related
    num_mutations = 0
    mut_exist = 0
    mut_type = empty_value_str
    mut_pos = empty_value_numeric
    mut_ref_nt = empty_value_str
    mut_nt = empty_value_str
    mut_site_distance = empty_value_numeric
    mut_struct = empty_value_str
    mut_ref_struct = empty_value_str
    mut_prev_struct = empty_value_str
    mut_next_struct = empty_value_str
    mut_same_as_site = None
    if 'mutation_syntax' in rna_item:
        mutation_syntax = rna_item['mutation_syntax'].strip()
        if mutation_syntax.lower().startswith('indel'):
            # indel
            mut_type = 'indel'
            mut_exist = 1
            num_mutations = len(structure_string) - len(wt_sequence)

            #
            features = list()
            features.append(rna_id)
            features.append(editing_value)
            features.append(fe_ensemble)
            features.append(simtree_normalized_score)
            # features.append(simtree_flipping_modes)

            features.append(num_mutations)
            features.append(mut_exist)
            features.append(mut_type)
            features.append(mut_pos)
            features.append(mut_site_distance)
            features.append(mut_ref_nt)
            features.append(mut_nt)
            features.append(mut_struct)
            features.append(mut_ref_struct)
            features.append(mut_prev_struct)
            features.append(mut_next_struct)
            features.append(mut_same_as_site)

            #
            rna_lib_features_list.append(features + structure_features)
        elif mutation_syntax.lower().startswith('wt'):
            # WT
            mut_type = 'wt'

            #
            features = list()
            features.append(rna_id)
            features.append(editing_value)
            features.append(fe_ensemble)
            features.append(simtree_normalized_score)
            # features.append(simtree_flipping_modes)

            features.append(num_mutations)
            features.append(mut_exist)
            features.append(mut_type)
            features.append(mut_pos)
            features.append(mut_site_distance)
            features.append(mut_ref_nt)
            features.append(mut_nt)
            features.append(mut_struct)
            features.append(mut_ref_struct)
            features.append(mut_prev_struct)
            features.append(mut_next_struct)
            features.append(mut_same_as_site)

            #
            rna_lib_features_list.append(features + structure_features)
        else:
            mut_type = 'mismatch'
            mutation_syntax_list = mutation_syntax.split(',')
            num_mutations = len(mutation_syntax_list)
            mut_exist = 1

            # Has mutations - loop each mutation
            # EACH mutation will have its own record
            for mutation_string in mutation_syntax_list:
                mut_pos = mutation_string[:-4]
                mut_ref_nt = mutation_string[-4]
                mut_nt = mutation_string[-1]
                mut_site_distance = int(mut_pos) - editing_position
                mut_struct = secondary_structure.get_annotation(int(mut_pos))
                # Struct Type from "WT" from the "SAME" position
                mut_ref_struct = wt_secondary_structure.get_annotation(int(mut_pos))

                #
                mut_same_as_site = 1 if mut_struct == site_struct else 0

                #
                mut_prev_struct = None
                pos = int(mut_pos)
                while pos:
                    pos -= 1
                    struct = secondary_structure.get_annotation(pos)
                    if struct != mut_struct:
                        mut_prev_struct = struct
                        break

                mut_next_struct = None
                pos = int(mut_pos)
                while pos:
                    pos += 1
                    struct = secondary_structure.get_annotation(pos)
                    if struct != mut_struct:
                        mut_next_struct = struct
                        break

                #
                features = list()
                features.append(rna_id)
                features.append(editing_value)
                features.append(fe_ensemble)
                features.append(simtree_normalized_score)
                # features.append(simtree_flipping_modes)

                features.append(num_mutations)
                features.append(mut_exist)
                features.append(mut_type)
                features.append(mut_pos)
                features.append(mut_site_distance)
                features.append(mut_ref_nt)
                features.append(mut_nt)
                features.append(mut_struct)
                features.append(mut_ref_struct)
                features.append(mut_prev_struct)
                features.append(mut_next_struct)
                features.append(mut_same_as_site)

                #
                rna_lib_features_list.append(features + structure_features)
    else:
        # No mutation related features
        features = list()
        features.append(rna_id)
        features.append(editing_value)
        features.append(fe_ensemble)
        features.append(simtree_normalized_score)
        # features.append(simtree_flipping_modes)

        features.append(num_mutations)
        features.append(mut_exist)
        features.append(mut_type)
        features.append(mut_pos)
        features.append(mut_site_distance)
        features.append(mut_ref_nt)
        features.append(mut_nt)
        features.append(mut_struct)
        features.append(mut_ref_struct)
        features.append(mut_prev_struct)
        features.append(mut_next_struct)
        features.append(mut_same_as_site)

        #
        rna_lib_features_list.append(features + structure_features)

# endregion


# ----------------------------------
# region Output - Features

headers = [
    'rna_id',
    'editing_value',
    'free_energy',
    'sim_nor_score',
    # 'sim_flip_mode',

    'num_mutations',
    'mut_exist',
    'mut_type',
    'mut_pos',
    'mut_site_dist',
    'mut_ref_nt',
    'mut_nt',
    'mut_struct',
    'mut_ref_struct',
    'mut_prev_struct',
    'mut_next_struct',
    'mut_same_as_site',

    'all_stem_length',

    'site_struct',

    'site_prev_nt',
    'site_next_nt',
    'site_prev_struct',
    'site_next_struct',
    'site_1_1',

    'site_length',
    'site_length_stem',
    'site_length_hairpin',
    'site_length_bulge',
    'site_length_internal_es',
    'site_length_internal_ecs',
    'site_5prm_cp_hairpin',
    'site_5prm_cp_bulge',
    'site_3prm_cp_bulge',
    'site_5prm_cp_internal',
    'site_3prm_cp_internal',

    'u_count',
    'u_all_stem_length',
    'u_hairpin_length',
]

# Build "header" for upstream structure elements
feature_template = [
    'exist',
    'distance',

    'struct',
    'length',
    'length_stem',
    'length_hairpin',
    'length_bulge',
    'length_internal_es',
    'length_internal_ecs',
    '5prm_cp_hairpin',
    '5prm_cp_bulge',
    '3prm_cp_bulge',
    '5prm_cp_internal',
    '3prm_cp_internal',
]
for index in range(1, num_upstream_elements + 1):
    #
    for template in feature_template:
        headers.append('u{}_{}'.format(str(index), template))

#
headers.append('d_count')
headers.append('d_all_stem_length')
for index in range(1, num_downstream_elements + 1):
    #
    for template in feature_template:
        headers.append('d{}_{}'.format(str(index), template))

# Write the "headers" line
writer = csv.writer(open(output_file_path, 'w'))
writer.writerow(headers)

# Write the content rows
for entry in rna_lib_features_list:
    # Replace all "None" value with string "None"
    # writer.writerow([val if val is not None else "None" for val in entry])
    writer.writerow(entry)

# endregion

# ----------------------------------
# region Output - bpRNA Results

if bprna_output_file_path:
    #
    output_content = ''

    #
    for rna_id, annotation_types in rna_lib_bprna_dict.items():
        for data_type, bprna_annotation in annotation_types.items():
            title = ','.join(['>' + rna_id, data_type])
            output_content += '\n'.join([title, bprna_annotation, ''])  # Add a "newline" to each of the record.

    #
    FileUtils.save_file(bprna_output_file_path, output_content, 'a')

# endregion

