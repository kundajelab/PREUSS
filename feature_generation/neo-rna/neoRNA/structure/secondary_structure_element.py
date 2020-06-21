# -*- coding: utf-8 -*-

"""
Secondary Structure Element
--------------------
"""

from typing import Optional, List, Tuple

from neoRNA.sequence import Sequence, BasePair


class SecondaryStructureElementType(object):
    """
    Secondary Structure Element Types

    Current Types:
    - Stem Loop - S
    - Hairpin Loop - H
    - Bulge Loop - B
    - Interior Loop - I
    - Multi-Loop - M
    - End sequence - E
    - Segment - segment

    """
    Stem = 'S'
    Hairpin = 'H'
    Bulge = 'B'
    Interior = 'I'
    Multiloop = 'M'
    End = 'E'
    Unpaired = 'X'
    Segment = 'SEGMENT'

    @classmethod
    def get_type(cls, type_str):
        """
        Get the element type based on the "type string"
        :param type_str:
        :return:
        """

        if type_str == 'S':
            return SecondaryStructureElementType.Stem
        elif type_str == 'H':
            return SecondaryStructureElementType.Hairpin
        elif type_str == 'B':
            return SecondaryStructureElementType.Bulge
        elif type_str == 'I':
            return SecondaryStructureElementType.Interior
        elif type_str == 'M':
            return SecondaryStructureElementType.Multiloop
        elif type_str == 'E':
            return SecondaryStructureElementType.End
        elif type_str == 'X':
            return SecondaryStructureElementType.Unpaired
        elif type_str == 'SEGMENT':
            return SecondaryStructureElementType.Segment
        else:
            return None


class SecondaryStructureElement(object):
    """
    Secondary Structure Element

    It helps define different types of "Secondary Structure Element" of a DNA / RNA, such as:
    - Stem Loop
    - Hairpin Loop
    - Bulge Loop
    - Interior Loop
    - etc.

    Ref - http://what-when-how.com/molecular-biology/rna-structure-molecular-biology/

    ## Structure Element

    Each of the "structure element" type has its own properties to help define it.

    """

    # ----------------------------------
    # region Init

    def __init__(self, ele_type, raw_string, base_pair_count=None):
        """
        Init

        :param ele_type: The Element Type
        :param raw_string: The raw string
        :param base_pair_count: The "base pair count" of a "segment"
        """
        #
        self.__ele_type = ele_type
        self.__raw_string = raw_string
        self.__base_pair_count = base_pair_count

        # Init Properties
        self.__sequence_list: List[Sequence] = list()
        self.__base_pair_list: List[BasePair] = list()

    # endregion

    # ----------------------------------
    # region Properties - General

    @property
    def ele_type(self):
        """
        Element Type

        :return:
        """
        return self.__ele_type

    @property
    def raw_string(self):
        """
        The raw string to define the structure element.

        :return: string
        """
        return self.__raw_string

    @raw_string.setter
    def raw_string(self, raw_string):
        self.__raw_string = raw_string

    @property
    def base_pair_count(self):
        return self.__base_pair_count

    @base_pair_count.setter
    def base_pair_count(self, count):
        self.__base_pair_count = count

    # endregion

    # ----------------------------------
    # region Properties - sequence Info

    @property
    def sequence_list(self) -> List[Sequence]:
        """
        The sequence list that this element contains.

        Depending on diff. element type, the `sequence list` may have diff. numbers of sequences.
        For example:
        - Stem Loop - 2 sequences
        - Hairpin Loop - 2 sequences

        :return:
        """
        return self.__sequence_list

    @sequence_list.setter
    def sequence_list(self, sequence_list: List[Sequence]):
        self.__sequence_list = sequence_list

    @property
    def base_pair_list(self) -> List[BasePair]:
        """
        The base pair list that this element contains.

        Depending on diff. element type, the `sequence list` may have diff. numbers of sequences.
        For example:
        - Stem Loop - 0
        - Hairpin Loop - 1 pairs
        - Bulge Loop - 2 pairs

        :return:
        """
        return self.__base_pair_list

    @base_pair_list.setter
    def base_pair_list(self, base_pair_list: List[BasePair]):
        self.__base_pair_list = base_pair_list

    # endregion

    # ----------------------------------
    # region Methods - Adjust

    def add_sequence(self, sequence: Sequence):
        r"""
        Add new sequence into the element.

        Parameters
        ----------
        sequence: Sequence

        Returns
        -------

        """

        self.__sequence_list.append(sequence)

    def add_base_pair(self, base_pair: BasePair):
        r"""
        Add new base pair into the element.

        Parameters
        ----------
        base_pair: BasePair

        Returns
        -------

        """

        self.__base_pair_list.append(base_pair)

    # endregion

    # ----------------------------------
    # region Methods - Position

    def is_contain(self, nt_position: int) -> Tuple[bool, Optional[Sequence]]:
        r"""
        Check if the element contains the given "nt position".

        Parameters
        ----------
        nt_position: int
            The "nt" position

        Returns
        -------
        ret_tuple: Tuple[bool, Optional[Sequence]]
            - flag: "True" if the element contains the given "nt position".
            - sequence: The sequence which contains the given "nt position"

        """

        if nt_position is None:
            return False, None

        for sequence in self.__sequence_list:
            if sequence.contain(nt_position) is True:
                return True, sequence

        return False, None

    def is_before(self, nt_position: int) -> bool:
        r"""
        Check if the element is before the given "nt position".

        Criteria:
            - "All" sequences inside must be "before" the given "nt position".

        Parameters
        ----------
        nt_position: int
            The "nt" position

        Returns
        -------
        flag: bool
            "True" if the element is before the given "nt position".
        """

        if nt_position is None:
            return False

        for sequence in self.__sequence_list:
            # Ignore "empty" sequence
            if not sequence.is_empty() and not sequence.is_before(nt_position):
                return False

        return True

    def is_after(self, nt_position: int) -> bool:
        r"""
        Check if the element is "after" the given "nt position".

        If the element is NOT in "before" or "contain", it is "after".

        Parameters
        ----------
        nt_position: int
            The "nt" position

        Returns
        -------
        flag: bool
            "True" if the element is after the given "nt position".
        """

        if nt_position is None:
            return False

        if self.is_contain(nt_position)[0] or self.is_before(nt_position):
            return False

        return True

    def distance(self, nt_position: int) -> Tuple[Optional[int], Optional[Sequence]]:
        r"""
        Calculate the "distance" between the element and given "nt position". The "return" also includes the "sequence"
        which contributes to the "distance".

        Depending on the "relationship" of the element and the nt position, it has diff. calc way:
        - contain - give "0" as the "distance"
        - before - Use "sequence's end-nt" to calc the "distance".
            If having multiple sequences, pick the "MAX" distance. All the distances are "negative".
        - after - Use "sequence's start-nt" to calc the "distance".
            If having multiple sequences, pick the "MIN POSITIVE" distance

        NOTE:
        - "Segment" type do not need to calc the distance

        Parameters
        ----------
        nt_position: int
            The "nt" position

        Returns
        -------

        """

        if not nt_position:
            return None, None

        if self.ele_type in [SecondaryStructureElementType.Segment]:
            return None, None

        # Try if it is "contain" relationship
        is_contain, contain_sequence = self.is_contain(nt_position)
        if is_contain:
            return 0, contain_sequence

        #
        from sys import maxsize
        if self.is_before(nt_position):
            ret_distance = -maxsize - 1  # Max negative int
            ret_sequence = None
            for sequence in self.__sequence_list:
                distance = sequence.distance(nt_position)
                if distance and distance > ret_distance:
                    ret_distance = distance
                    ret_sequence = sequence

            #
            return ret_distance, ret_sequence
        else:
            # After
            ret_distance = maxsize
            ret_sequence = None
            for sequence in self.__sequence_list:
                distance = sequence.distance(nt_position, is_before=False)
                if distance and 0 < distance < ret_distance:
                    ret_distance = distance
                    ret_sequence = sequence

            #
            return ret_distance, ret_sequence

    # endregion

    # ----------------------------------
    # region Methods - Stats

    def length(self, ref_position: int = None):
        r"""
        Get the "length" of the element.

        If "ref_position" is given, get the "length" of the "sequence" which includes the "ref_position".

        Parameters
        ----------
        ref_position: int
            The given "position" as reference.
            This is only useful for "Interior", "Multiloop" elements - since this element may have "diff." sequences.

        Returns
        -------

        """

        if self.__ele_type == SecondaryStructureElementType.Stem:
            # The length of "one" sequence
            if len(self.__sequence_list) == 2:
                return self.__sequence_list[0].length
        elif self.__ele_type == SecondaryStructureElementType.Hairpin \
                or self.__ele_type == SecondaryStructureElementType.Bulge \
                or self.__ele_type == SecondaryStructureElementType.End:
            # The length of "FIRST" sequence
            if len(self.__sequence_list) == 1:
                return self.__sequence_list[0].length
        elif self.__ele_type == SecondaryStructureElementType.Interior:
            # If "nt position" not provide, the length is a tuple for "BOTH" sequences
            if len(self.__sequence_list) == 2:
                if not ref_position:
                    # Return the "tuple"
                    return self.__sequence_list[0].length, self.__sequence_list[1].length
                else:
                    for sequence in self.__sequence_list:
                        if sequence.contain(ref_position):
                            return sequence.length
        elif self.__ele_type == SecondaryStructureElementType.Multiloop:
            # If "nt position" not provide, the length is a list
            if len(self.__sequence_list) > 0:
                if not ref_position:
                    return [sequence.length for sequence in self.__sequence_list]
                else:
                    for sequence in self.__sequence_list:
                        if sequence.contain(ref_position):
                            return sequence.length
        else:
            # For other types, like `Segment`, return None
            return None

    def closing_base_pair(self, ref_position: int = None):
        r"""
        Get the "Closing Base Pair(s)" of the element.

        Depending on the "element type", it will return differently.

        Parameters
        ----------
        ref_position: int
            The given "position" as reference.
            This is only useful for the element which has "multiple" sequences.

        Returns
        -------

        """

        if self.__ele_type == SecondaryStructureElementType.Hairpin:
            # Return the "FIRST" pair
            if len(self.__base_pair_list) == 1:
                return self.__base_pair_list[0]
        elif self.__ele_type == SecondaryStructureElementType.Bulge \
                or self.__ele_type == SecondaryStructureElementType.Interior:
            # Return a tuple for "TWO" pairs
            if len(self.__base_pair_list) == 2:
                return self.__base_pair_list[0], self.__base_pair_list[1]
        else:
            # For other types, like `Stem`, `End`, `Segment`, return None
            return None

    # endregion

    # ----------------------------------
    # region Magic

    def __str__(self):

        return '[{} | {}]'.format(self.ele_type, self.raw_string)

    # endregion
