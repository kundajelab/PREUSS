# -*- coding: utf-8 -*-

"""
ShapeMapper 2.x Profile Item
================
"""

from typing import List


class ShapeProfileItem(object):
    """
    ShapeMapper 2.x Profile Item Object.

    Each item represent the "profile info" for a single "Nucleotide".

    It includes the following elements:
    - Nucleotide
    - Sequence
    - Modified_mutations
    - Modified_read_depth
    - Modified_effective_depth
    - Modified_rate
    - Untreated_mutations
    - Untreated_read_depth
    - Untreated_effective_depth
    - Untreated_rate
    - Denatured_mutations
    - Denatured_read_depth
    - Denatured_effective_depth
    - Denatured_rate
    - Reactivity_profile
    - Std_err
    - HQ_profile
    - HQ_stderr
    - Norm_profile
    - Norm_stderr

    """

    # ----------------------------------
    # region Init

    def __init__(self, data_list: List[str]):
        r"""
        Init

        The parameter is a data list (in `str`).
        While doing init, we will need to convert it to "float" for some data elements.

        NOTE:
        Ne sure to update the "list index" if the data format changes.

        Parameters
        ----------
        data_list: List[str]
            The data, in a "list" format
        """

        self.nt_position = data_list[0].strip()
        self.nt_sequence = data_list[1].strip()

        # For others, need to convert it to `float`
        self.modified_mutations = float(data_list[2].strip()) if data_list[2].strip() != 'nan' else None
        self.modified_read_depth = float(data_list[3].strip()) if data_list[3].strip() != 'nan' else None
        self.modified_effective_depth = float(data_list[4].strip()) if data_list[4].strip() != 'nan' else None
        self.modified_rate = float(data_list[5].strip()) if data_list[5].strip() != 'nan' else None

        self.untreated_mutations = float(data_list[6].strip()) if data_list[6].strip() != 'nan' else None
        self.untreated_read_depth = float(data_list[7].strip()) if data_list[7].strip() != 'nan' else None
        self.untreated_effective_depth = float(data_list[8].strip()) if data_list[8].strip() != 'nan' else None
        self.untreated_rate = float(data_list[9].strip()) if data_list[9].strip() != 'nan' else None

        self.denatured_mutations = float(data_list[10].strip()) if data_list[10].strip() != 'nan' else None
        self.denatured_read_depth = float(data_list[11].strip()) if data_list[11].strip() != 'nan' else None
        self.denatured_effective_depth = float(data_list[12].strip()) if data_list[12].strip() != 'nan' else None
        self.denatured_rate = float(data_list[13].strip()) if data_list[13].strip() != 'nan' else None

        self.reactivity_profile = float(data_list[14].strip()) if data_list[14].strip() != 'nan' else None
        self.reactivity_stderr = float(data_list[15].strip()) if data_list[15].strip() != 'nan' else None
        self.hq_profile = float(data_list[16].strip()) if data_list[16].strip() != 'nan' else None
        self.hq_stderr = float(data_list[17].strip()) if data_list[17].strip() != 'nan' else None

        # It is possible to not have "norm" value if the data source is in "low quality"
        self.in_high_quality = True
        if data_list[18]:
            self.norm_profile = float(data_list[18].strip()) if data_list[18].strip() != 'nan' else None
            self.norm_stderr = float(data_list[19].strip()) if data_list[19].strip() != 'nan' else None
        else:
            # Not in "HQ"
            self.in_high_quality = False
            self.norm_profile = None
            self.norm_profile = None

    # endregion

    # ----------------------------------
    # region Properties

    # endregion

    # ----------------------------------
    # region Methods

    # endregion

    # ----------------------------------
    # region JSON Object

    # endregion
