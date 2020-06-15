# -*- coding: utf-8 -*-

"""
Process Script - Process RNA Lib "Profiling" Results
================
"""

import os
import json

import numpy as np

from typing import List, Optional

from neoRNA.library.library_item import LibraryItem
from neoRNA.sequence import Sequence
from neoRNA.structure import DotBracketNotation
from neoRNA.util.file_utils import FileUtils

from neoRNA.util.json_serializable import as_python_object, PythonObjectEncoder
from neoRNA.util.runner.matlab_runner import MatlabRunner

#
OUTPUT_FOLDER_BIERS_RESULTS = 'biers_results'
OUTPUT_FOLDER_BIERS_INFERENCE_STRUCTURE = 'biers_inference_structure'


def biers_rna_structure(rna_items_json: str, output_folder: str,
                        max_bootstrap: int = 20,
                        override: bool = True,
                        sequence_offset: int = 0, sequence_start: int = 1, sequence_end: Optional[int] = None):
    r"""
    Run "RNA Secondary Structure Inferring" (by Biers) and output the results.

    In Biers, it utilizes "RNAStructure" package to help infer the RNA Secondary Structure.

    ## Output

    The output includes two parts:
    - RNA Secondary Structure Inferring results
    - A "dot-bracket" file to record the "inference structure"
        - This file is also used for "bpRNA" step.

    NOTE:
        - These two files will be in "different" folders.
        - For the "results" file, it will be in "JSON" format, so that it can be easily processed in the following steps.

    Parameters
    ----------
    rna_items_json: str
        The "RNA Lib Item" object, as string
    output_folder: str
        The "path" of output folder. The results will be under this folder.
    max_bootstrap: int
        The "#" of bootstrap that Biers will run. "0" means NO bootstrap to run.
    override: bool
        If need to override the results.
    sequence_offset: int
        The "position offset" of the sequence to be inferred. Default to "0" - no offset.
    sequence_start: int
        The "start position" of the sequence to be inferred. Default to "1" - from the beginning of the sequence.
    sequence_end: Optional[int]
        The "end position" of the sequence to be inferred. Default to "None" - till the "end" of the sequence.

    Returns
    -------

    """

    # Parse the "Python Object"
    rna_item: LibraryItem = json.loads(rna_items_json, object_hook=as_python_object)

    bootstrap_enabled = True if max_bootstrap > 0 else False

    #
    rna_id = rna_item.rna_id
    sequence: Sequence = rna_item.sequence

    # Prep for the inference
    if sequence_end is None:
        # Use the "entire" sequence
        seqpos_out = range(sequence_start + sequence_offset, sequence.length + sequence_offset)
    else:
        #
        seqpos_out = range(sequence_start + sequence_offset, sequence_end + sequence_offset)

    # Slice used to help extract the "sequence" and "reactivity array"
    sequence_slice = slice(sequence_start - 1, sequence_end)
    target_rna_sequence_str = str(sequence.get_rna_sequence())[sequence_slice]
    # target_reactivity = rna_item.neo_reactivity_list[sequence_slice]
    reactivity_list = rna_item.flatten_reactivity_list(reactivity_type='shape')
    target_reactivity = reactivity_list[sequence_slice]

    #
    MatlabRunner.start_engine()

    # Predict "Reference Structure"
    structure_na, bpp_na, bootstrap_na \
        = MatlabRunner.rna_structure(target_rna_sequence_str, sequence_offset, seqpos_out)

    if bootstrap_enabled:
        # Predict Structure based on reactivity data (in 1D, DMS format)
        structure_dms_1d_fold, bpp_dms_1d_fold, structure_dms_1d_fold_bootstrap \
            = MatlabRunner.rna_structure(target_rna_sequence_str, sequence_offset, seqpos_out,
                                         max_bootstrap=max_bootstrap, reactivity_1d=target_reactivity)

        # Output - Inference Structure Dot-Bracket Notation
        structure_dms = DotBracketNotation('{}'.format(rna_id),
                                           target_rna_sequence_str, structure_dms_1d_fold)
        structure_dms_file_path = os.path.join(output_folder,
                                               OUTPUT_FOLDER_BIERS_INFERENCE_STRUCTURE, '{}.dbn'.format(rna_id))
        structure_dms.to_file(structure_dms_file_path)

    #
    MatlabRunner.stop_engine()

    # Output - Biers results
    if bootstrap_enabled:
        json_data = {
            "rna_id": rna_id,
            "sequence_string": target_rna_sequence_str,
            "reactivity": target_reactivity,

            #
            "computational_structure": structure_na,
            "computational_bpp": np.array(bpp_na).tolist(),
            "experimental_structure": structure_dms_1d_fold,
            "experimental_bpp": np.array(bpp_dms_1d_fold).tolist(),
            "bootstrap_structures": structure_dms_1d_fold_bootstrap,
        }
    else:
        json_data = {
            "rna_id": rna_id,
            "sequence_string": target_rna_sequence_str,

            #
            "computational_structure": structure_na,
            "computational_bpp": np.array(bpp_na).tolist(),
        }

    #
    biers_results_file_path = os.path.join(output_folder,
                                           OUTPUT_FOLDER_BIERS_RESULTS, '{}.json'.format(rna_id))
    FileUtils.save_json_to_file(biers_results_file_path, json_data)


def biers_results_folder_path(working_folder):
    #
    folder_path = os.path.join(working_folder, OUTPUT_FOLDER_BIERS_RESULTS)

    # Check if need to build the folder
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    #
    return folder_path


def biers_inference_structure_folder_path(working_folder):
    #
    folder_path = os.path.join(working_folder, OUTPUT_FOLDER_BIERS_INFERENCE_STRUCTURE)

    # Check if need to build the folder
    if not os.path.exists(folder_path):
        os.mkdir(folder_path)

    #
    return folder_path
