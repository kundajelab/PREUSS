# -*- coding: utf-8 -*-

"""
Process Script - Process RNA Lib "Structure" Results
================
"""

import os
import json
import csv

from typing import List, Optional

from neoRNA.library.rna_library import RnaLibrary
from neoRNA.sequence.sequence import Sequence

from neoRNA.util.file_utils import FileUtils

from neoRNA.util.json_serializable import as_python_object, PythonObjectEncoder


def generate_rna_lib_structure_results(rna_library_json: str,
                                       editing_level_file_path: str, biers_results_folder_path: str, bprna_results_folder_path: str,
                                       structure_summary_file_path,
                                       sequence_start: int = 1, sequence_end: Optional[int] = None,
                                       sequence_position_counter_offset: int = 0):

    r"""
    Generate RNA Lib "Structure" results file.

    It should be called after "Biers & bpRNA" finishes the run.

    NOTE:
        - `sequence_start` and `sequence_end` will be applied to the "original" sequence, so be sure to make them valid
        against the "original" sequence.
        - `sequence_position_counter_offset` will be "ALWAYS" applied to the "initial position" - "1".

    Parameters
    ----------
    rna_library_json: str
        The "RNA Lib" object, in `string` format
    editing_level_file_path: str
        The "file path" to editing level file.
    biers_results_folder_path: str
        The "folder path" to Biers results.
    bprna_results_folder_path: str
        The "folder path" to bpRNA results.
    structure_summary_file_path: str
        The "file path" to Structure Results file.
    sequence_start: int
        Start point of the sequence. Default to "1" (from beginning).
    sequence_end: int
        End point of the sequence. Default to "None" (till end).
    sequence_position_counter_offset: int
        The "position counter offset" of the sequence to be inferred.
        Default to "0" - no offset, start the counter from "1".

    Returns
    -------

    """

    # Parse the "Python Object"
    rna_library: RnaLibrary = json.loads(rna_library_json, object_hook=as_python_object)

    # Load Editing Level data
    # Editing Level values & positions, indexed by "RNA ID"
    editing_levels = {}
    editing_positions = {}

    # Columns to load
    rna_id_column_name = 'RNA_ID_STR'
    editing_level_column_name = 'Avg'
    editing_position_column_name = 'Editing_position'  # Based on the "FULL" sequence

    #
    editing_level_data = csv.DictReader(open(editing_level_file_path, 'rU'))
    for line in editing_level_data:
        rna_id = line[rna_id_column_name]
        #
        score = line[editing_level_column_name]
        editing_levels[rna_id] = float(score) if score != "NA" and score != "#N/A" else None
        #
        position = line[editing_position_column_name]
        editing_positions[rna_id] = position

    wt_sequence = Sequence(rna_library.wide_type_rna_sequence)

    # Load data from each RNA Item
    rna_library_structure_summary_items = []
    for rna_item in rna_library.rna_items:
        #
        rna_id = rna_item.rna_id
        sequence = rna_item.sequence
        barcode = rna_item.barcode

        # Load Biers Results
        biers_result_json_file_path = os.path.join(biers_results_folder_path, '{}.json'.format(rna_id))
        with open(biers_result_json_file_path) as infile:
            biers_result_json = json.load(infile)

        # Mut Syntax
        mut_positions, mut_syntax = sequence.generate_mutation_syntax(wt_sequence.sequence_str,
                                                                      seq_type_rna=True,
                                                                      sequence_start=sequence_start, sequence_end=sequence_end)

        # Skip the item if there is "NO" editing values
        if rna_id not in editing_levels.keys():
            continue

        # Need to check if "Biers results" include "interred" data
        if 'experimental_structure' in biers_result_json:
            summary_item = {
                "rna_id": biers_result_json['rna_id'],
                "sequence_string": biers_result_json['sequence_string'],
                "mutation_syntax": ','.join(mut_syntax) if mut_syntax else None,
                "reactivity": biers_result_json['reactivity'],

                "computational_structure": biers_result_json['computational_structure'],
                "computational_bpp": biers_result_json['computational_bpp'],
                "experimental_structure": biers_result_json['experimental_structure'],
                "experimental_bpp": biers_result_json['experimental_bpp'],
                "bootstrap_structures": biers_result_json['bootstrap_structures'],

                #
                "A-to-I_editing_level": editing_levels[rna_id],
                "A-to-I_editing_site": int(editing_positions[rna_id]) - (sequence_start - 1),
            }
        else:
            summary_item = {
                "rna_id": biers_result_json['rna_id'],
                "sequence_string": biers_result_json['sequence_string'],
                "mutation_syntax": ','.join(mut_syntax) if mut_syntax else None,

                "computational_structure": biers_result_json['computational_structure'],
                "computational_bpp": biers_result_json['computational_bpp'],

                #
                "A-to-I_editing_level": editing_levels[rna_id],
                "A-to-I_editing_site": int(editing_positions[rna_id]) - (sequence_start - 1),
            }

        rna_library_structure_summary_items.append(summary_item)

    #
    rna_library_structure_summary_dict = {
        "code": rna_library.data_source_code,
        "title": rna_library.data_source_title,
        "date": rna_library.data_source_date,
        "notes": rna_library.running_notes,
        "wt_id": rna_library.wide_type_rna_id,
        "failed_rna_id": [],
        "indel_rna_id": [],
        "items": rna_library_structure_summary_items,
    }
    json_data = json.dumps(rna_library_structure_summary_dict,
                           sort_keys=True,
                           indent=4, ensure_ascii=True)

    # Output JSON as file
    FileUtils.save_file(structure_summary_file_path, json_data)

