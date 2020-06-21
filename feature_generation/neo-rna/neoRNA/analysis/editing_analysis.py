# -*- coding: utf-8 -*-

"""
RNA Editing
--------------------
"""

from typing import List, Dict

from neoRNA.analysis.editing_analysis_item import EditingAnalysisItem, EditingAnalysisItemType
from neoRNA.sequence import Sequence
from neoRNA.structure.secondary_structure_element import SecondaryStructureElementType
from neoRNA.structure.secondary_structure import SecondaryStructure


class EditingAnalysis(object):
    """
    RNA Editing Analysis

    It helps introduce the info and related analysis processes of "RNA Editing".

    The analysis results will be represented as a list of "analysis items".
    - Each of the analysis items may be in different types.

    """

    # ----------------------------------
    # region Init

    def __init__(self, secondary_structure: SecondaryStructure, editing_position):
        """
        Init

        :param secondary_structure: The secondary structure object.
        :param editing_position: The "editing position" (a.k.a nt position)
        """
        #
        self.__secondary_structure = secondary_structure
        self.__editing_position = editing_position

        #
        self.__analysis_items = list()
        self.__analysis_items_by_distance = dict()  # Organize the "analysis item" by its "distance"

    # endregion

    # ----------------------------------
    # region Property

    @property
    def analysis_items(self) -> List[EditingAnalysisItem]:
        return self.__analysis_items

    @property
    def rna_id(self) -> str:
        return self.__secondary_structure.reference_id

    @property
    def sequence(self) -> Sequence:
        return self.__secondary_structure.sequence

    @property
    def analysis_items_by_distance(self) -> Dict[int, EditingAnalysisItem]:
        return self.__analysis_items_by_distance

    # endregion

    # ----------------------------------
    # region Methods

    def has_result(self):
        """
        Check if the analysis has results returned.

        NOTE:
        - It just check if there is any "analysis item" generated.

        :return:
        """
        return len(self.__analysis_items) > 0

    def generate_string_for_items(self):

        ret = ''
        for item in self.__analysis_items:
            ret += str(item)

        return ret

    def get_analysis_items(self, analysis_types: List[EditingAnalysisItemType] = None,
                           element_types: List[SecondaryStructureElementType] = None) -> List[EditingAnalysisItem]:
        r"""
        Get a sublist of "analysis items" based on the given "type lists".

        The "type" can be:
        - Analysis Type
        - Element Type

        Parameters
        ----------
        analysis_types: List[EditingAnalysisItemType]
            The list of "Analysis Type" needed.
        element_types: List[SecondaryStructureElementType]
            The list of "Element Type" needed.

        Returns
        -------
        analysis_items: List[EditingAnalysisItem]
            A sublist of "analysis items".

        """

        # If not set, return all
        if analysis_types is None and element_types is None:
            return self.__analysis_items

        items = []

        if element_types is None:
            # Check only on "analysis type"
            for item in self.__analysis_items:
                if item.analysis_type in analysis_types:
                    items.append(item)

            return items
        elif analysis_types is None:
            # Check only on "element type"
            for item in self.__analysis_items:
                if item.element_type in element_types:
                    items.append(item)

            return items
        else:
            # Check both "analysis type" and "element type"
            for item in self.__analysis_items:
                if item.analysis_type in analysis_types and item.element_type in element_types:
                    items.append(item)
            return items

    # endregion

    # ----------------------------------
    # region Methods - Analysis

    def analysis(self):
        """
        Analysis each element inside the list.

        For each of the element:
            - Determine the "element type" (before, contain, after)

        :return:
        """

        for element in self.__secondary_structure.elements:
            item = None

            # Check if an element "contains" the "editing position".
            # `is_contain` returns includes "two" elements - first one is the flag
            if EditingAnalysis.can_analyze_contain(element) and element.is_contain(self.__editing_position)[0]:
                item = EditingAnalysisItem(EditingAnalysisItemType.Contain, element)
            elif EditingAnalysis.can_analyze_before(element) and element.is_before(self.__editing_position):
                # Check if an element "is before" the "editing position".
                item = EditingAnalysisItem(EditingAnalysisItemType.Before, element)
            elif element.is_after(self.__editing_position):
                # Check if an element "if after" the "editing position"
                item = EditingAnalysisItem(EditingAnalysisItemType.After, element)

            #
            if item:
                # Add the "item" to analysis list
                self.__analysis_items.append(item)
                distance, distance_sequence = element.distance(self.__editing_position)
                if type(distance) is int:
                    self.__analysis_items_by_distance[distance] = item, distance_sequence

    # endregion

    # ----------------------------------
    # region Internal Methods

    @classmethod
    def can_analyze_contain(cls, element):
        """
        Check if the element can be used to analyze "contain" type.

        NOTE:
        The following types of element can analyze for "contain" type.
        - Stem
        - Bulge
        - Hairpin
        - Interior

        :param element:
        :return:
        """

        if element is None:
            return False

        return element.ele_type == SecondaryStructureElementType.Stem \
               or element.ele_type == SecondaryStructureElementType.Hairpin \
               or element.ele_type == SecondaryStructureElementType.Interior \
               or element.ele_type == SecondaryStructureElementType.Multiloop \
               or element.ele_type == SecondaryStructureElementType.Unpaired \
               or element.ele_type == SecondaryStructureElementType.Bulge \
               or element.ele_type == SecondaryStructureElementType.End

    @classmethod
    def can_analyze_before(cls, element):
        """
        Check if the element can be used to analyze "before" type.

        NOTE:
        The following types of element can analyze for "before" type.
        - Stem
        - Bulge
        - Hairpin
        - Interior

        :param element:
        :return:
        """

        if element is None:
            return False

        return element.ele_type == SecondaryStructureElementType.Stem \
               or element.ele_type == SecondaryStructureElementType.Hairpin \
               or element.ele_type == SecondaryStructureElementType.Multiloop \
               or element.ele_type == SecondaryStructureElementType.Unpaired \
               or element.ele_type == SecondaryStructureElementType.Bulge \
               or element.ele_type == SecondaryStructureElementType.End

    # endregion

    # ----------------------------------
    # region Magic

    def __str__(self):

        #
        return '[{} | {} | {}]'.format(str(self.__secondary_structure), self.__editing_position,
                                       self.generate_string_for_items())

    # endregion
