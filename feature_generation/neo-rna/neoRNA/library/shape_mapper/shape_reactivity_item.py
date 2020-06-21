# -*- coding: utf-8 -*-

"""
ShapeMapper 2.x - Shape Reactivity Item
================
"""

from typing import List


class ShapeReactivityItem(object):
    """
    ShapeMapper 2.x Shape Reactivity Item Object.

    Each item represent the "Shape Reactivity" for a single "Nucleotide".

    It includes the following elements:
    - Nucleotide Position
    - Shape Reactivity

    """

    # "-999" means no reactivity data.
    NUMBER_FOR_NONE = -999.0

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

        # For others, need to convert it to `float`
        self.shape_reactivity = float(data_list[1].strip()) \
            if data_list[1].strip() != self.NUMBER_FOR_NONE else None

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
