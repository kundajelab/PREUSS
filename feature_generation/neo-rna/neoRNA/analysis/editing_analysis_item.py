# -*- coding: utf-8 -*-

"""
RNA Editing
--------------------
"""

from typing import Optional, Tuple

from neoRNA.sequence import Sequence
from neoRNA.structure.secondary_structure_element import SecondaryStructureElement, SecondaryStructureElementType


class EditingAnalysisItemType(object):
    """
    Editing Analysis Item Types

    Current Types:
    - Contained - Contained by an element

    """

    # Types
    Contain = 'Contain'
    After = 'After'
    Before = 'Before'

    @classmethod
    def get_type_string(cls, analysis_type):
        """
        Get the "readbale" type string based on the given "type"

        :param analysis_type:
        :return:
        """
        #
        if analysis_type == EditingAnalysisItemType.Contain:
            return 'Contain the editing level'
        if analysis_type == EditingAnalysisItemType.Before:
            return 'Before the editing level'
        if analysis_type == EditingAnalysisItemType.After:
            return 'After the editing level'
        else:
            return None


class EditingAnalysisItem(object):
    """
    RNA Editing Analysis Item

    The "basic" analysis item that can help convey different types of analysis results.

    The "analysis type" may include:
    - the "editing position" is within certain secondary structure element
    - the RNA contains "5' Internal loop" - to the 3' of "editing position"
    - the RNA contains "3' Internal loop" - to the 5' of "editing position"

    """

    # ----------------------------------
    # region Init

    def __init__(self, analysis_item_type: EditingAnalysisItemType, related_element: SecondaryStructureElement):
        """
        Init

        :param analysis_item_type: The "type" of analysis result.
        :param related_element: The related secondary structure element.
        """
        #
        self.__analysis_item_type = analysis_item_type
        self.__related_element = related_element

    # endregion

    # ----------------------------------
    # region Property

    @property
    def analysis_type(self):
        return self.__analysis_item_type

    @property
    def element(self):
        return self.__related_element

    @property
    def element_type(self):
        if self.__related_element is not None:
            return self.__related_element.ele_type

        return None

    # endregion

    # ----------------------------------
    # region Methods - Analysis

    def sequence(self, five_prim: bool = True) -> Optional[Sequence]:
        r"""
        Get "sequence", from 5' or 3'.

        Parameters
        ----------
        five_prim: bool
            If it is from 5'. Only used for "Interior" type element.

        Returns
        -------
        sequence: Sequence
        """

        if self.__related_element is None:
            return None

        # Get the "length" of the element
        sequence_list = self.__related_element.sequence_list

        # For "Interior" type, use `five_prim` to determine which one to return
        if self.__related_element.ele_type == SecondaryStructureElementType.Interior:
            # `sequence_list` has at least "2" sequences
            #   - first - 5'
            #   - second - 3'
            if five_prim is True:
                return sequence_list[0]
            else:
                return sequence_list[1]
        #
        return sequence_list[0]

    def complementary_strand_sequence(self, ref_sequence: Sequence) -> Optional[Sequence]:
        r"""
        Get the "Complementary Strand Sequence" by the given "ref_sequence".

        Only the following "element type" can have "Complementary Strand Sequence":
        - Stem
        - Interior

        Parameters
        ----------
        ref_sequence: Sequence

        Returns
        -------

        """

        if self.__related_element is None:
            return None

        element_type = self.__related_element.ele_type
        if element_type not in [
            SecondaryStructureElementType.Interior,
            SecondaryStructureElementType.Stem
        ]:
            return None

        #
        sequence_list = self.__related_element.sequence_list
        for sequence in sequence_list:
            #
            if sequence.sequence_str != ref_sequence.sequence_str:
                return sequence

        #
        return None

    def complementary_nt(self, ref_sequence: Sequence, ref_position: int) -> Tuple[Optional[int], Optional[str]]:
        r"""
        Get the "Complementary nt" by the given "ref_sequence" and "ref_position".

        Only the following "element type" can have "Complementary nt":
        - Stem


        Parameters
        ----------
        ref_sequence: Sequence
        ref_position: int

        Returns
        -------

        """

        #
        if self.__related_element is None:
            return None, None

        #
        element_type = self.__related_element.ele_type
        if element_type not in [
            SecondaryStructureElementType.Stem
        ]:
            return None, None

        # Get the "Complementary Sequence"
        sequence_list = self.__related_element.sequence_list
        complementary_sequence = None
        for sequence in sequence_list:
            #
            if sequence.sequence_str != ref_sequence.sequence_str:
                complementary_sequence = sequence

        #
        if complementary_sequence is None:
            return None, None

        #
        local_ordering = ref_sequence.nt_ordering(ref_position)
        return complementary_sequence.get_nt_by_ordering(-ref_sequence.length + local_ordering - 1)

    def loop_length(self, five_prim=True):
        """
        Get "Loop Length", from 5' or 3'.

        :param five_prim: If it is from 5'. Only used for "Interior" type element.
        :return:
        """

        if self.__related_element is None:
            return None

        # Get the "length" of the element
        length = self.__related_element.length()

        # For "Interior" type, use `five_prim` to determine which one to return
        if self.__related_element.ele_type == SecondaryStructureElementType.Interior:
            # `length` is a tuple,
            #   - first - 5'
            #   - second - 3'
            if five_prim is True:
                return length[0]
            else:
                return length[1]
        elif self.__related_element.ele_type == SecondaryStructureElementType.Multiloop:
            return None

        #
        return length

    def loop_length_by_position(self, position: int) -> int:
        r"""
        Get "Loop Length", whose "sequence" includes the given "position".

        Parameters
        ----------
        position: int
            The given "position"

        Returns
        -------
        loop_length: int
        """

        if self.__related_element is None:
            return None

        # Get the "length" of the element, by the "position"
        length = self.__related_element.length(ref_position=position)

        # Get the "first" one
        return length[0]

    def closing_base_pair(self, five_prim=True):
        """
        Get "closing pair", from 5' or 3'.

        :param five_prim: If it is from 5'. Used for "Interior", "Bulge", "Unpaired" types.
        :return:
        """

        if self.__related_element is None:
            return None

        #
        closing_pair = self.__related_element.closing_base_pair()

        if not closing_pair:
            return None

        # For "Interior" type, use `five_prim` to determine which one to return
        if self.__related_element.ele_type == SecondaryStructureElementType.Bulge \
                or self.__related_element.ele_type == SecondaryStructureElementType.Unpaired \
                or self.__related_element.ele_type == SecondaryStructureElementType.Interior:
            # `closing pair` is a tuple,
            #   - first - 5'
            #   - second - 3'
            if five_prim is True:
                return closing_pair[0].get_nt_pair_string()
            else:
                return closing_pair[1].get_nt_pair_string()
        #
        return closing_pair.get_nt_pair_string()

    def mismatch(self, five_prim=True):
        """
        Get "mismatch", from 5' or 3'.

        :param five_prim:
        :return:
        """

        if self.__related_element is None:
            return None

        return None

    # endregion

    # ----------------------------------
    # region Magic

    def __str__(self):

        return '[{} | {}]'.format(self.__analysis_item_type, str(self.__related_element))

# endregion
