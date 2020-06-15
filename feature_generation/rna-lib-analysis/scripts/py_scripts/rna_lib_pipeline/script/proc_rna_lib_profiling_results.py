# -*- coding: utf-8 -*-

"""
Process Script - Process RNA Lib "Profiling" Results
================
"""

import json

from typing import List

from neoRNA import io
from neoRNA.library.library_config import RnaLibConfig
from neoRNA.library.library_item import LibraryItem
from neoRNA.library.rna_library import RnaLibrary
from neoRNA.library.shape_mapper import get_shape2_path
from neoRNA.util.file_utils import FileUtils

from neoRNA.util.json_serializable import as_python_object, PythonObjectEncoder


def generate_rna_lib_profiling_results(configs_json: str, rna_items_json: str, output_file: str):
    r"""
    Generate RNA Lib "profiling" results file.

    It should be called after "ShapeMapper 2.x" finishes the run.

    Parameters
    ----------
    configs_json: str
        The "RNA Lib" configs, in `string` format
    rna_items_json: str
        The RNA Lib Item object
    output_file: str
        The "file path" of data output.

    Returns
    -------

    """

    # Parse the "Python Object"
    configs: RnaLibConfig = json.loads(configs_json, object_hook=as_python_object)
    rna_items: List[LibraryItem] = json.loads(rna_items_json, object_hook=as_python_object)

    #
    updated_rna_items = []
    for rna_item in rna_items:
        #
        rna_id = rna_item.rna_id
        barcode = rna_item.barcode

        # Load "Profile" data
        shape_profile = get_shape2_path('profile', configs.working_folder, rna_id, barcode.barcode)
        for nt_profile in io.parse(shape_profile, "shape-profile"):
            #
            nt_position = nt_profile.nt_position
            rna_item.profile_list.append(nt_profile)
            rna_item.profile_dict[nt_position] = nt_profile

        # Load "Shape Reactivity" data
        shape_reactivity = get_shape2_path('shape', configs.working_folder, rna_id, barcode.barcode)
        for nt_reactivity in io.parse(shape_reactivity, "shape-reactivity"):
            #
            nt_position = nt_reactivity.nt_position
            rna_item.shape_reactivity_list.append(nt_reactivity)
            rna_item.shape_reactivity_dict[nt_position] = nt_reactivity

        updated_rna_items.append(rna_item)

    # Build RNA Library object
    library = RnaLibrary()
    library.load_meta(configs)

    #
    library.rna_items = updated_rna_items

    # Output
    data_str = json.dumps(library, cls=PythonObjectEncoder)
    FileUtils.save_file(output_file, data_str)
