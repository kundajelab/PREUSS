# -*- coding: utf-8 -*-

"""
Dot-Bracket Notation
================

Process the "dot-bracket" notation of RNA secondary structure.
"""


import os
from Bio.Seq import Seq


class DotBracketNotation(object):
    """
    Process the "dot-bracket" notation of RNA secondary structure.

    The "dot-bracker" notation is consisting of dots '.', opening '(' and closing ')' parentheses.

    - "Dotted" position represents unpaired nucleotide,
    - The "matching parenthesized" positions represent base-pairing nucleotides.

    The purpose of processing this notation is try to identify the different types of nucleotide, such as:

    - stem - Regions of contiguous canonical Watson-Crick base-paired nucleotides.
    - interior loop - Double-stranded unpaired regions flanked by stems on either side.
    - hairpin loop - aka. Stem-loop.
    - multiloop - Single-stranded unpaired regions.
    - etc.

    ## Current Support

    To keep simple, the current "type notation" support the following types:

    - "s" - stem, both `(` and `)`
    - "h" - hairpin loop
    - "-" - Unknown type


    link: http://ultrastudio.org/en/Dot-Bracket_Notation
    """

    __SEQUENCE_FOCUS = slice(22, 63)

    def __init__(self, rna_structure, sequence=None):
        """

        :param rna_structure:
        :param sequence:
        """
        self.__rna_structure = rna_structure
        self.__sequence = sequence

        # Generate the "type" notation
        self.__type_notation = self.generate_type_string_v1(self.__rna_structure)

        # Stats
        self.__stem_length = self.__calculate_length_stem()
        self.__hairpin_length = self.__calculate_length_hairpin_loop()

        self.__stem_length_with_focus = self.__calculate_length_stem(self.__SEQUENCE_FOCUS)
        self.__hairpin_length_with_focus = self.__calculate_length_hairpin_loop(self.__SEQUENCE_FOCUS)

    # region Property

    @property
    def rna_structure(self):
        return self.__rna_structure

    @property
    def type_notation(self):
        return self.__type_notation

    @property
    def stem_length(self):
        return self.__stem_length

    @property
    def hairpin_length(self):
        return self.__hairpin_length

    @property
    def stem_length_with_focus(self):
        """
        Only count the length of stem within certain "focus range".
        :return:
        """
        return self.__stem_length_with_focus

    @property
    def hairpin_length_with_focus(self):
        """
        Only count the length of hairpin loop within certain "focus range".
        :return:
        """
        return self.__hairpin_length_with_focus

    # endregion

    # region Private Methods


    def __calculate_length_hairpin_loop(self, sequence_focus=None):
        """ Calculate the length of hairpin loops.

        If having multiple ones, sum them up.

        :return:
        """
        if sequence_focus is not None:
            return self.__type_notation[sequence_focus].count('h')

        return self.__type_notation.count('h')

    def __calculate_length_stem(self, sequence_focus=None):
        """ Calculate the length of stems.

        If having multiple ones, sum them up.

        :return:
        """
        if sequence_focus is not None:
            return self.__type_notation[sequence_focus].count('s') / 2

        return self.__type_notation.count('s') / 2

    # endregion


    # region Class Methods

    @classmethod
    def generate_type_string_v1(cls, rna_structure):
        """Generate "type" string based on the given structure notation.

        Example:
        ..(((..(((((......)))))..))).. will be notated as:
        --ssss--sssshhhhhhsssss--sss--

        :return: The "type" string, with the same length of given dot-bracket notation string.
        """

        if len(rna_structure) == 0:
            return None

        type_string = []

        # Flag to indicate if there is an "open" stem found previously.
        open_stem_marked = False
        # Keep track on the "type" of the previous position
        prev_notation = None
        number_of_pending_dot = 0

        #
        for idx, element in enumerate(rna_structure):
            type_char = None

            if element == '(':
                # In any case, a '(' should mark an open stem (if not marked)
                if prev_notation != '(':
                    open_stem_marked = True
                if prev_notation == '.' and number_of_pending_dot > 0:
                    # Solve a list of pending dots as "multiloop"
                    solved_sot_string = '-' * number_of_pending_dot
                    type_string.append(solved_sot_string)
                    number_of_pending_dot = 0
                type_char = 's'
            elif element == ')':
                if prev_notation == '.':
                    # This ')' just secure a hairpin loop
                    solved_sot_string = 'h' * number_of_pending_dot
                    type_string.append(solved_sot_string)
                    number_of_pending_dot = 0
                type_char = 's'
            elif element == '.':
                if prev_notation == ')':
                    # A '.' after a closed stem.
                    open_stem_marked = False
                    type_char = '-'
                elif prev_notation == '(':
                    # Still pending - could belong to a "hairpin loop" or a "multiloop"
                    number_of_pending_dot += 1
                elif prev_notation == '.':
                    if number_of_pending_dot > 0:
                        # continue pending for the dots
                        number_of_pending_dot += 1
                    else:
                        type_char = '-'
                else:
                    type_char = '-'

            if type_char is not None:
                type_string.append(type_char)
            prev_notation = element

        return ''.join(type_string)

    # endregion