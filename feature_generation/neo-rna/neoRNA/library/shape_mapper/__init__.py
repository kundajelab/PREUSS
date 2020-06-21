# -*- coding: utf-8 -*-

"""
ShapeMapper Module
================

It defines helper classes for understand the process / data by ShapeMapper.
"""

import os


# ----------------------------------
# region Constants

SHAPE_RESULTS_FOLDER_NAME = 'shapemapper_results'
SHAPE_TEMP_FOLDER_NAME = 'shapemapper_temp'

# The "naming pattern" for the RNA Lib item result folder.
# Format: [rna_id]_[barcode]_out
RNA_ITEM_SHAPE_RESULT_FOLDER_NAME_PATTERN = '{}_{}_out'
RNA_ITEM_SHAPE_TEMP_FOLDER_NAME_PATTERN = '{}_{}_temp'

# Files
RNA_ITEM_SHAPE_PROFILE_FILE_NAME_PATTERN = 'Pipeline_{}_profile.txt'
RNA_ITEM_SHAPE_SHAPE_FILE_NAME_PATTERN = 'Pipeline_{}.shape'

# endregion


# ----------------------------------
# region Naming Conventions
#

def get_shape2_path(path_type: str, working_folder: str,
                    rna_id: str = None, rna_barcode: str = None):
    r"""
    Get the "path" of "file/folder" under ShapeMapper 2.x

    Parameters
    ----------
    path_type: str
        The "type" of the path to get.
        Available options:
        - results - folder, result folder
        - temp - folder, temp folder
        - profile - file, profile file
        - shape - file, shape data file
    working_folder: str
        Path to the working folder.
    rna_id: str
        RNA ID
    rna_barcode: str
        Barcode of the RNA item

    Returns
    -------
    path: str
        The "path" str.

    """

    if path_type == 'results':
        return os.path.join(working_folder, SHAPE_RESULTS_FOLDER_NAME)
    if path_type == 'temp':
        return os.path.join(working_folder, SHAPE_TEMP_FOLDER_NAME)
    if path_type == 'profile':
        return os.path.join(working_folder, SHAPE_RESULTS_FOLDER_NAME,
                            RNA_ITEM_SHAPE_RESULT_FOLDER_NAME_PATTERN.format(rna_id, rna_barcode),
                            RNA_ITEM_SHAPE_PROFILE_FILE_NAME_PATTERN.format(rna_id))
    if path_type == 'shape':
        return os.path.join(working_folder, SHAPE_RESULTS_FOLDER_NAME,
                            RNA_ITEM_SHAPE_RESULT_FOLDER_NAME_PATTERN.format(rna_id, rna_barcode),
                            RNA_ITEM_SHAPE_SHAPE_FILE_NAME_PATTERN.format(rna_id))

    return None

# endregion
