# -*- coding: utf-8 -*-

"""
RNA Library Item
================
"""

import numpy

from typing import List, Dict

from neoRNA.library.shape_mapper.shape_profile_item import ShapeProfileItem
from neoRNA.library.shape_mapper.shape_reactivity_item import ShapeReactivityItem
from neoRNA.sequence.sequence import Sequence
from neoRNA.sequence.barcode import Barcode


class LibraryItem(object):
    """
    RNA library Item Object.

    Each item should include the following elements:

    - RNA ID: The unique ID (usually a number) throughout the library.
    - RNA Barcode: The unique "barcode" sequence that represents this RNA
    - RNA sequence: The actual target sequence of this RNA
    - Notes: A string to describe this RNA, supplementary info.

    """

    # The default value for "invalid value"
    INVALID_VALUE = -999.0

    # ----------------------------------
    # region Init

    def __init__(self, rna_id: str, rna_barcode_string: str, rna_sequence_string: str,
                 notes: str = None):
        r"""
        Init

        Parameters
        ----------
        rna_id: str
            RNA ID
        rna_barcode_string: str
            The RNA Library Item barcode string
        rna_sequence_string: str
            The RNA Library Item sequence string
        notes: str
            Notes content
        """

        self.rna_id: str = rna_id
        self.barcode: Barcode = Barcode(rna_barcode_string)
        self.sequence: Sequence = Sequence(rna_sequence_string)
        self.notes: str = notes

        # ----------
        # Profile data from ShapeMapper 2.x
        # - The list of profile data for each of "nt".
        #   - The elements in the list follows the "ordering" of nt position.
        # - The `dict` is indexed by "nt position"
        self.profile_list: List[ShapeProfileItem] = []
        self.profile_dict: Dict[str, ShapeProfileItem] = {}

        # ----------
        # Shape Reactivity
        # - The value is from ShapeMapper 2.x
        # - The `dict` is indexed by "nt position"
        self.shape_reactivity_list: List[ShapeReactivityItem] = []
        self.shape_reactivity_dict: Dict[str, ShapeReactivityItem] = {}

        # ----------
        # Reactivity from "OWN" method - a diff. method from ShapeMapper 2.x
        #
        # NOTE:
        # - Based on the calculation algorithm, the elements inside this list may not
        #   include the "entire" original sequence.
        self.neo_reactivity_list: List[float] = []

        # ----------------------------------
        # Stats
        self.modified_read_depth = None
        self.untreated_read_depth = None

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def total_nt(self) -> int:
        r"""
        Get the total number of the "nt" - the length of sequence.

        Returns
        -------
        nt_length: int
            The "length" of the sequence.

        """
        return self.sequence.length

    # endregion

    # ----------------------------------
    # region Methods - Stats

    def shape_profile_list_low_quality(self, nt_a_c_only: bool = True) -> List[ShapeProfileItem]:
        r"""
        Retrieve the list of "ShapeMapper 2.x profile" which have "low quality".

        "low quality" is determined by the flag - `in_high_quality` inside the "Profile Item".

        Parameters
        ----------
        nt_a_c_only: bool
            If only do the calculation based on "A", "C" nt.

        Returns
        -------
        low_quality_list: List[ShapeProfileItem]
            The list of "low quality" profile item, by ShapeMapper 2.x.

        """
        #
        if nt_a_c_only:
            self.sequence.calculate_length()
            return [item for index, item in enumerate(self.profile_list)
                    if self.sequence.is_nt_ac(index + 1) and item.in_high_quality is False]

        #
        return [item for item in self.profile_list if item.in_high_quality is False]

    def total_nt_with_condition(self, nt_a_c_only: bool = True) -> int:
        r"""
        Get the "total #" of "nt", with condition.

        Parameters
        ----------
        nt_a_c_only: bool
            If only do the calculation based on "A", "C" nt.

        Returns
        -------
        nt_length: int
            The "length" of the sequence, with condition.
        """

        #
        if not nt_a_c_only:
            return self.total_nt

        #
        total_nt = 0
        for index, item in enumerate(self.profile_list):
            #
            if self.sequence.is_nt_ac(index + 1):
                total_nt += 1

        #
        return total_nt

    # endregion

    # ----------------------------------
    # region Methods - Reactivity List

    def flatten_reactivity_list(self, reactivity_type: str = 'own') -> List[float]:
        r"""
        Flatten the "reactivity object list" to a list of "single reactivity value".

        Parameters
        ----------
        reactivity_type: str
            Which type of "reactivity" data to get.

        Returns
        -------
        flatten_reactivity_list: List[float]
            The list of "reactivity" data, ordered by "nt position"
        """

        if reactivity_type == 'shape':
            # Convert "None"
            return [item.shape_reactivity if item.shape_reactivity is not None else ShapeReactivityItem.NUMBER_FOR_NONE
                    for item in self.shape_reactivity_list]
            # return [item.shape_reactivity for item in self.shape_reactivity_list]
        if reactivity_type == 'own':
            return self.neo_reactivity_list

    # endregion

    # ----------------------------------
    # region Methods Reactivity Calculation

    def calculate_reactivity_v1(self,
                                sequence_slice: slice = slice(0, None),
                                nt_a_c_only: bool = True):
        r"""
        Calculate "reactivity" with own method.

        Version 1 - The simplest version.

        General Algorithm:
        - For the rate of "each" nt, calculate a `adjusted rate` = `mod rate` - `non_mod rate`
            - if it is `< 0`, use `0`
        - Normalize the "rate" by dividing it against the **maximum rate** front the "targeting" list.
        - Use this `normalized rate` list as the "1D Reactivity" data.

        NOTE:
        - Since this algorithm needs to find the "max rate" from the "targeting" list,
            it needs to pass the actual "slice" to help define the "targeting" list.

        Parameters
        ----------
        sequence_slice: slice
            The target "sequence slice".
            Default to `slice(0, None)` - the entire sequence.
        nt_a_c_only: bool
            If only do the calculation based on "A", "C" nt.

        Returns
        -------

        """

        # Get the list of "rates" - "modified" and "untreated"
        # - Directly use `reactivity_profile` (= modified - untreated)
        # modified_rate = [item.modified_rate for item in self.profile_list]
        # untreated_rate = [item.untreated_rate for item in self.profile_list]
        reactivity_profile = [item.reactivity_profile for item in self.profile_list]

        # --------------
        # Adjusted profile rate
        # Rules
        #  - negative value -> 0
        #  - "None" -> INVALID
        #  - If "AC-only", INVALID for "GU" nt.
        #
        # - negative value -> 0
        reactivity_profile_adjusted \
            = [abs(rate) if rate is not None and rate < 0.0 else rate for rate in reactivity_profile]
        # - "None" -> INVALID
        reactivity_profile_adjusted \
            = [rate if rate is not None else self.INVALID_VALUE for rate in reactivity_profile_adjusted]
        if nt_a_c_only:
            # INVALID for "GU" nt.
            self.sequence.calculate_length()
            reactivity_profile_adjusted \
                = [rate if rate is not None and self.sequence.is_nt_ac(index + 1) else self.INVALID_VALUE for index, rate in enumerate(reactivity_profile_adjusted)]

        # Convert it to "numpy array"
        reactivity_profile_adjusted = numpy.array(reactivity_profile_adjusted, dtype=float)

        # Normalize the rates
        reactivity_profile_adjusted = reactivity_profile_adjusted[sequence_slice]  # Apply the "sequence slice"
        max_rate = numpy.amax(reactivity_profile_adjusted)
        self.neo_reactivity_list = \
            numpy.array([rate / max_rate if rate != self.INVALID_VALUE else rate for rate in reactivity_profile_adjusted])
        self.neo_reactivity_list = self.neo_reactivity_list.tolist()

    # endregion

    # ----------------------------------
    # region Method


    # endregion
