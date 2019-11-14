# -*- coding: utf-8 -*-

"""
RNA Secondary Structure
--------------------
"""

from typing import List, Optional

from neoRNA.sequence.sequence import Sequence
from neoRNA.structure.secondary_structure_element import SecondaryStructureElement, SecondaryStructureElementType
from neoRNA.util.parser.comment_parser import CommentParser


class SecondaryStructure(object):
    """
    Secondary Structure

    It defines the presentation of "Secondary Structure" info of a RNA.
    The info consists of a list of "Secondary Structure Elements" as well as some basic info.

    """

    # ----------------------------------
    # region Init

    def __init__(self, comment: str, sequence_str: str,
                 dot_bracket_str: str, dot_bracket_annotation_str: str, dot_bracket_validation_str: str,
                 elements: List[SecondaryStructureElement],
                 reference_id=None):
        r"""
        Init

        Parameters
        ----------
        comment: str
            The "comment info" of the structure.
        sequence_str: str
            The associated sequence string.
        dot_bracket_str: str
        dot_bracket_annotation_str: str
        dot_bracket_validation_str: str
        elements: [SecondaryStructureElement]
        reference_id: str
        """
        #
        self.__comment = comment
        self.__sequence = Sequence(sequence_str)
        self.__dot_bracket_str = dot_bracket_str
        self.__dot_bracket_annotation_str = dot_bracket_annotation_str
        self.__dot_bracket_validation_str = dot_bracket_validation_str
        self.__elements = elements

        #
        self.__reference_id = reference_id if reference_id else self.__parse_id_from_comment(self.__comment)

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def comment(self) -> str:
        return self.__comment

    @property
    def dot_bracket(self) -> str:
        return self.__dot_bracket_str

    @property
    def dot_bracket_annotation(self) -> str:
        return self.__dot_bracket_annotation_str

    @property
    def elements(self) -> List[SecondaryStructureElement]:
        return self.__elements

    @property
    def reference_id(self)-> str:
        return self.__reference_id

    @property
    def sequence(self) -> Sequence:
        return self.__sequence

    # endregion

    # ----------------------------------
    # region Methods

    def get_elements(self, element_type=None):
        """
        Get the generator of the elements list, based on "type" if needed

        :param element_type:
        :return:
        """

        for element in self.__elements:
            if element_type is None:
                yield element
            elif element.ele_type == element_type:
                yield element

    def get_annotation(self, position: int) -> Optional[str]:
        r"""
        Get the "annotation string" based on the given "position".

        Parameters
        ----------
        position: int
            The "position" to be checked. Start from "1".

        Returns
        -------
        annotation: str
            The "annotation string".
        """

        #
        if 0 < position <= len(self.__dot_bracket_annotation_str):
            return self.__dot_bracket_annotation_str[position - 1]

        return None

    # endregion

    # ----------------------------------
    # region Private Methods

    def __parse_id_from_comment(self, comment_str):

        if not comment_str:
            return None

        info_dict = CommentParser.parse_bp_rna(comment_str)
        if info_dict and 'id' in info_dict:
            return info_dict['id']

        return None

    # endregion

    # ----------------------------------
    # region Export

    def export_csv(self, element_type):
        """
        Export the secondary structure info in a "CSV" format.

        The actual exported data will be diff. depending on the "element type".

        :param element_type:
        :return:
        """

        data = []

        if element_type == SecondaryStructureElementType.Stem:
            # "Stem" data: RNA_ID, Total Element Count, Total bp, bp for each Stem Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_bp = 0
            element_total_bp = []
            element_raw_strings = []
            for element in element_list:
                total_bp += element.length()
                element_total_bp.append(element.length())
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_bp)
            data.append(' | '.join([str(e) for e in element_total_bp]))
            data.append(' | '.join(element_raw_strings))

            #
            return data
        elif element_type == SecondaryStructureElementType.Hairpin:
            # "Hairpin" data: RNA_ID, Total Element Count, Total nt, nt for each Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_nt = 0
            element_total_nt = []
            element_raw_strings = []
            for element in element_list:
                total_nt += element.length()
                element_total_nt.append(element.length())
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_nt)
            data.append(' | '.join([str(e) for e in element_total_nt]))
            data.append(' | '.join(element_raw_strings))

            #
            # return ','.join(data)
            return data
        elif element_type == SecondaryStructureElementType.Bulge:
            # "Bulge" data: RNA_ID, Total Element Count, Total nt, nt for each Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_nt = 0
            element_total_nt = []
            element_raw_strings = []
            for element in element_list:
                total_nt += element.length()
                element_total_nt.append(element.length())
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_nt)
            data.append(' | '.join([str(e) for e in element_total_nt]))
            data.append(' | '.join(element_raw_strings))

            #
            # return ','.join(data)
            return data
        elif element_type == SecondaryStructureElementType.Interior:
            # "Interior" data: RNA_ID, Total Element Count, Total nt, nt for each Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_nt = 0
            element_total_nt = []
            element_raw_strings = []
            for element in element_list:
                left_length, right_length = element.length()
                total_nt += left_length
                total_nt += right_length
                element_total_nt.append(element.length())
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_nt)
            data.append(' | '.join([str(e) for e in element_total_nt]))
            data.append(' | '.join(element_raw_strings))

            #
            # return ','.join(data)
            return data
        elif element_type == SecondaryStructureElementType.End:
            # "End" data: RNA_ID, Total Element Count, Total nt, nt for each Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_nt = 0
            element_total_nt = []
            element_raw_strings = []
            for element in element_list:
                total_nt += element.length()
                element_total_nt.append(element.length())
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_nt)
            data.append(' | '.join([str(e) for e in element_total_nt]))
            data.append(' | '.join(element_raw_strings))

            #
            # return ','.join(data)
            return data
        elif element_type == SecondaryStructureElementType.Segment:
            # "End" data: RNA_ID, Total Element Count, Total bp, bp for each Loop, Raw strings
            element_list = list(self.get_elements(element_type))
            element_count = len(element_list)
            total_bp = 0
            element_raw_strings = []
            for element in element_list:
                total_bp += element.base_pair_count
                element_raw_strings.append(element.raw_string)
            #
            data.append(self.__reference_id)
            data.append(element_count)
            data.append(total_bp)
            data.append(' | '.join(element_raw_strings))

            #
            # return ','.join(data)
            return data

    # endregion

    # ----------------------------------
    # region Magic

    def __str__(self):

        return '[{} | {}]'.format(self.__reference_id, str(self.__sequence))

    # endregion
