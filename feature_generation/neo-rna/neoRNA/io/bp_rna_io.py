# -*- coding: utf-8 -*-

"""
IO - bpRNA
================

It help process the result file from the tool - 'bpRNA'.
"""


import os

from typing import Tuple, Optional

from neoRNA.sequence.sequence import Sequence
from neoRNA.sequence.base_pair import BasePair
from neoRNA.structure.secondary_structure_element import SecondaryStructureElement, SecondaryStructureElementType
from neoRNA.structure.secondary_structure import SecondaryStructure


class BpRnaIO(object):
    r"""
    Parse the "result" file of 'bpRNA'.

    The result file defines the "RNA secondary structure".

    The file is a "plaintext" file.
        - Refer to the example file `/tests/io/example_files/bprna_example.st` for its format.
        - ref: bpRNA - http://bprna.cgrb.oregonstate.edu/download.php#bpRNA

    Current Format:
    - "#" line - references, may have multiple lines
    - Line 1 - RNA sequence
    - Line 2 - Dot-Bracket string
    - Line 3 - Dot-Bracket annotation string
    - Line 4 - Dot-Bracket annotation validation string
    - Other lines - Secondary structure elements, like stem, bulge, etc.

    """

    # The "marker" for each of the "Record"
    # - A record usually starts with a "comment".
    RECORD_MARKER = '#'

    # Used for storing parsing results
    INFO = {}
    ELEMENTS = []

    # Temp location for "Interior Loop"
    interior_left_raw_string = None
    interior_left_sequence = None
    interior_left_base_pair = None

    # Helper variables for "Multi-loop"
    multiloop_current_no = 0
    multiloop_current_element: SecondaryStructureElement = None
    multiloop_current_element_raw_str_list = list()

    # ----------------------------------
    # region Iterator Generator

    @classmethod
    def parse_iterator(cls, handle):
        r"""
        Iterate over records and parse it as objects.

        Parameters
        ----------
        handle: any
            input file.

        Returns
        -------
        parsed_objects: BasePairProbabilityIO
            Parsed objects.


        Usage
        -------

        >>> with open("bp-rna.st") as handle:
        ...     for record in BpRnaIO.parse_iterator(handle):
        ...         print(record.comment)
        ...

        """

        for comment, info, elements in cls.parse(handle):
            yield SecondaryStructure(comment,
                                     info['sequence'],
                                     info['dot_bracket'],
                                     info['dot_bracket_annotation'],
                                     info['dot_bracket_validation'],
                                     elements)

    # endregion

    # ----------------------------------
    # region Parser

    @classmethod
    def parse(cls, handle):
        """
        Parse the file.

        Parameters
        ----------
        handle: handle
            input file.

        Returns
        -------
            A tuple of strings (comment, bpp_str).

        """

        # Skip any text before the first record (e.g. blank lines, comments)
        while True:
            line = handle.readline()
            if line == "":
                return
            if line[0] == cls.RECORD_MARKER:  # Find the "first line" of a record.
                break

        while True:
            #
            if line[0] != cls.RECORD_MARKER:
                raise ValueError(
                    "Records should start with '{}' character!".format(cls.RECORD_MARKER))

            # Reset
            cls.INFO = {}
            cls.ELEMENTS = []

            # Helper variables for "Interior loop"
            cls.interior_left_raw_string = None
            cls.interior_left_sequence = None
            cls.interior_left_base_pair = None

            # Helper variable for "Multiloop"
            cls.multiloop_current_no = 0
            cls.multiloop_current_element = None
            cls.multiloop_current_element_raw_str_list = []

            # Get "Reference Info" from the first line
            # Ex: #Name: 001_with_reactivity
            comment = line[7:].rstrip()  # Remove the first part - "#Name: ".
            # Skip other 2 lines
            for index in range(1, 3):
                handle.readline()

            # Parse the "Sequence Info"
            # It has "4" lines
            for index in range(1, 5):
                cls.parse_basic_info(handle.readline(), index)

            # Parse the "elements" till the next "record marker"
            line = handle.readline()
            while True:
                if not line:
                    break
                if line[0] == cls.RECORD_MARKER:
                    break

                #
                cls.parse_element(line)

                #
                line = handle.readline()

            # Check if there is any "leftover" before "return"
            if cls.multiloop_current_element:
                cls.multiloop_current_element.raw_string = ' | '.join(cls.multiloop_current_element_raw_str_list)
                cls.ELEMENTS.append(cls.multiloop_current_element)

                #
                cls.multiloop_current_element = None
                cls.multiloop_current_element_raw_str_list = []
                cls.multiloop_current_no = 0

            yield comment, cls.INFO, cls.ELEMENTS

            if not line:
                return  # StopIteration

    # endregion

    # ----------------------------------
    # region Methods - Parsing Functions

    @classmethod
    def parse_basic_info(cls, line, line_parsed):
        """
        Parse the basic info from the "read line".

        :param line:
        :param line_parsed:
        :return:
        """

        if line_parsed == 1:
            # sequence
            cls.INFO['sequence'] = line.strip()
        elif line_parsed == 2:
            # Dot-bracket string
            cls.INFO['dot_bracket'] = line.strip()
        elif line_parsed == 3:
            # Dot-bracket annotation
            cls.INFO['dot_bracket_annotation'] = line.strip()
        elif line_parsed == 4:
            # Dot-bracket validation
            cls.INFO['dot_bracket_validation'] = line.strip()
        else:
            return

    @classmethod
    def parse_element(cls, line):
        """
        Parse the "secondary structure element" info.

        :param line:
        :return:
        """

        if not line.strip():
            return

        import re

        __delimiters = '\s'
        parts = re.split(__delimiters, line.strip())

        if not len(parts) >= 2:
            raise ValueError('The line must have at least 2 parts', line)

        # Decode the element index info
        element_type, element_no, element_index = cls.decode_element_index(parts[0])

        element = SecondaryStructureElement(element_type, line.strip())
        if element_type == SecondaryStructureElementType.Stem:
            # It has "2" sequences
            if len(parts) == 5:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                # 2
                start_position, end_position = cls.determine_position_pair(parts[3])
                sequence_str = cls.determine_sequence(parts[4])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))

                #
                cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.Hairpin:
            # It has "1" sequence and "1" base pair
            if len(parts) == 5:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                # 2
                position_pair = cls.determine_position_pair(parts[3])
                nt_str = cls.determine_nt_pair(parts[4])
                element.add_base_pair(BasePair(nt_str, position_pair))

                #
                cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.Bulge:
            # It has "1" sequence and "2" base pairs
            if len(parts) == 7:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                # 2
                position_pair = cls.determine_position_pair(parts[3])
                nt_str = cls.determine_nt_pair(parts[4])
                element.add_base_pair(BasePair(nt_str, position_pair))
                # 3
                position_pair = cls.determine_position_pair(parts[5])
                nt_str = cls.determine_nt_pair(parts[6])
                element.add_base_pair(BasePair(nt_str, position_pair))

                #
                cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.Unpaired:
            # It has "1" sequence and "2" base pairs
            if len(parts) == 7:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                # 2
                position_pair = cls.determine_position_pair(parts[3])
                nt_str = cls.determine_nt_pair(parts[4])
                element.add_base_pair(BasePair(nt_str, position_pair))
                # 3
                position_pair = cls.determine_position_pair(parts[5])
                nt_str = cls.determine_nt_pair(parts[6])
                element.add_base_pair(BasePair(nt_str, position_pair))

                #
                cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.Multiloop:
            # It has "1" sequence and "2" base pairs
            if len(parts) == 7:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                # 2
                position_pair_1 = cls.determine_position_pair(parts[3])
                nt_str_1 = cls.determine_nt_pair(parts[4])
                # 3
                position_pair_2 = cls.determine_position_pair(parts[5])
                nt_str_2 = cls.determine_nt_pair(parts[6])

                # Check if it is a "new" multiloop element
                if cls.multiloop_current_no == 0 or cls.multiloop_current_no != element_no:
                    # Save the previous element
                    if cls.multiloop_current_no != 0:
                        cls.multiloop_current_element.raw_string = ' | '.join(
                            cls.multiloop_current_element_raw_str_list)
                        cls.ELEMENTS.append(cls.multiloop_current_element)
                        #
                        cls.multiloop_current_element = None
                        cls.multiloop_current_element_raw_str_list = []
                        cls.multiloop_current_no = 0

                    # A "new" multiloop
                    cls.multiloop_current_no = element_no
                    cls.multiloop_current_element = element

                    #
                    cls.multiloop_current_element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                    cls.multiloop_current_element.add_base_pair(BasePair(nt_str_1, position_pair_1))
                    cls.multiloop_current_element.add_base_pair(BasePair(nt_str_2, position_pair_2))
                    cls.multiloop_current_element_raw_str_list.append(element.raw_string)
                else:
                    # For an "exist" multiloop
                    cls.multiloop_current_element.add_sequence(
                        Sequence(sequence_str, range(start_position, end_position + 1)))
                    cls.multiloop_current_element.add_base_pair(BasePair(nt_str_1, position_pair_1))
                    cls.multiloop_current_element.add_base_pair(BasePair(nt_str_2, position_pair_2))
                    cls.multiloop_current_element_raw_str_list.append(element.raw_string)
        elif element_type == SecondaryStructureElementType.Interior:
            # Each "Interior Loop" includes "TWO" lines to present
            # Each line includes "1" sequence and "1" base pair
            if len(parts) == 5:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                sequence = Sequence(sequence_str, range(start_position, end_position + 1))
                # 2
                position_pair = cls.determine_position_pair(parts[3])
                nt_str = cls.determine_nt_pair(parts[4])
                base_pair = BasePair(nt_str, position_pair)

                # Determine if needs to add a new element or just save it as a "temp" data'
                if cls.interior_left_sequence is None:
                    cls.interior_left_raw_string = element.raw_string
                    cls.interior_left_sequence = sequence
                    cls.interior_left_base_pair = base_pair
                else:
                    # Ready to add "Interior Loop" as a new element
                    element.add_sequence(cls.interior_left_sequence)
                    element.add_sequence(sequence)
                    element.add_base_pair(cls.interior_left_base_pair)
                    element.add_base_pair(base_pair)
                    element.raw_string = '{} | {}'.format(cls.interior_left_raw_string, element.raw_string)

                    cls.interior_left_raw_string = None
                    cls.interior_left_sequence = None
                    cls.interior_left_base_pair = None
                    #
                    cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.End:
            # It has "1" sequence
            if len(parts) == 3:
                # 1
                start_position, end_position = cls.determine_position_pair(parts[1])
                sequence_str = cls.determine_sequence(parts[2])
                element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))

                #
                cls.ELEMENTS.append(element)
        elif element_type == SecondaryStructureElementType.Segment:
            # It has "1" sequence
            if len(parts) > 2 and len(parts) % 2 == 0:
                # 1
                count = cls.determine_base_pair_count(parts[1])
                element.base_pair_count = count

                # Loop till the end to add all sequences
                index = 2
                while index < len(parts):
                    start_position, end_position = cls.determine_position_pair(parts[index])
                    sequence_str = cls.determine_sequence(parts[index+1])
                    element.add_sequence(Sequence(sequence_str, range(start_position, end_position + 1)))
                    index += 2

                #
                cls.ELEMENTS.append(element)

    # endregion

    # ----------------------------------
    # region Methods - Static

    @staticmethod
    def decode_element_index(element_index_str: str) -> Tuple[Optional[str], Optional[int], Optional[int]]:
        r"""
        Decode the "element index string" to retrieve the following parts:
        - element type
        - element no#
        - element internal index (Optional)

        Example "element index string":
        - S1
        - I1.1
        - segment1

        Parameters
        ----------
        element_index_str: str

        Returns
        -------

        """

        if not element_index_str.strip():
            return None, None, None

        import re

        # regex = '([\D]+)'  # Non-digit string
        regex = '^(?P<type>[\D]+)(?P<no>[\d]+)?\D*(?P<index>[\d]+)?$'
        found = re.search(regex, element_index_str.strip())

        # Decode the results
        result_dict = found.groupdict() if found else {}
        element_type = None
        element_no = None
        element_index = None
        if 'type' in result_dict:
            element_type = SecondaryStructureElementType.get_type(result_dict['type'].upper())
        if 'no' in result_dict:
            element_no = result_dict['no']
        if 'index' in result_dict:
            element_index = result_dict['index']

        return element_type, element_no, element_index

    @staticmethod
    def determine_sequence(sequence_str: str):
        r"""
        Determine the "sequence" from a given string.

        Example format:
        - "CUU"
        - ""

        Parameters
        ----------
        sequence_str: str

        Returns
        -------

        """

        if not sequence_str.strip():
            return None

        import re
        regex = '[\W]*([a-zA-Z]+)'  # Retrieve one string
        matched = re.match(regex, sequence_str.strip())

        if not matched:
            return None
        if len(matched.groups()) != 1:
            raise ValueError('The string does not contain 1 string', sequence_str)

        return matched.group(1)

    @staticmethod
    def determine_position_pair(position_pair_str: str):
        r"""
        Determine the "position pair" ( 2 digits) from a given string.

        The pair could be "start" and "end" of a sequence, or two positions of a base pair.

        Example format:
        - 48..50
        - (47,30)

        Parameters
        ----------
        position_pair_str: str

        Returns
        -------

        """

        if not position_pair_str.strip():
            return None

        import re
        regex = '[\D]*([\d]+)[\D]+([\d]+)' # Retrieve two numbers
        matched = re.match(regex, position_pair_str.strip())

        if len(matched.groups()) != 2:
            raise ValueError('The string does not contain 2 numbers', position_pair_str)

        return int(matched.group(1)), int(matched.group(2))

    @staticmethod
    def determine_nt_pair(string):
        """
        Determine the "nt pair" ( 2 letters) from a given string.

        Example format:
        - G:C

        :param string:
        :return: A range
        """

        if not string.strip():
            return None

        import re
        regex = '[\W]*([\w]):([\w])' # Retrieve two letters
        matched = re.match(regex, string.strip())

        if len(matched.groups()) != 2:
            raise ValueError('The string does not contain 2 letters', string)

        return matched.group(1), matched.group(2)

    @staticmethod
    def determine_base_pair_count(string):
        """
        Determine the "Base Pair Count" based on the given "string".

        The format of the string may be like:
        - 25bp

        :param string:
        :return:
        """

        if not string.strip():
            return None

        import re

        regex = '([\d]+)'  # Retrieve a number
        matched = re.match(regex, string.strip())

        if len(matched.groups()) != 1:
            raise ValueError('The string does not contain 1 number', string)

        return int(matched.group(1))

    # endregion
