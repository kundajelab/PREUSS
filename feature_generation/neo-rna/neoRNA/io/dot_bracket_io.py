# -*- coding: utf-8 -*-

"""
IO - Dot Bracket Notation
================

It help process the "Dot Bracket" format file.
"""

from neoRNA.structure.dot_bracket_notation import DotBracketNotation


class DotBracketIO(object):
    """
    Parse the "dot bracket" file.

    Example:
        >001 centroid 1 655
        GGGAGCCUGCCCUCUGAUCUCUGCCUGUUCCUCUGUCCCACAGAGGGCAAAGGCUACGGGUCAGAGAGCGGGGAGGAGGAC
        .....(((.(((((((.(((((((((((.(((.(((((......))))).)))..))))).)))))).))))).)))))..

    The format of "dot bracket file":
    - First line - comment line, starting with ">". Usually it contains the "basic info" of the sequence structure.
    - Second line - Sequence string.
    - Third line - dot-bracket notation string.

    """

    # The "marker" for each of the "Record"
    # - A record usually starts with a "comment".
    RECORD_MARKER = '>'

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
        parsed_objects: DotBracketNotation
            Parsed objects.


        Usage
        -------

        >>> with open("dot-bracket.dbn") as handle:
        ...     for record in DotBracketIO.parse_iterator(handle):
        ...         print(record.comment)
        ...

        """

        for comment, sequence_string, structure_string in cls.parse(handle):
            yield DotBracketNotation(comment, sequence_string, structure_string)

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
            A tuple of strings (comment, sequence_string, structure_annotation_string).

        """

        # Skip any text before the first record (e.g. blank lines, comments)
        while True:
            line = handle.readline()
            if line == "":
                return
            if line[0] == cls.RECORD_MARKER:  # Find the "first line" of a record.
                break

        while True:
            if line[0] != cls.RECORD_MARKER:
                raise ValueError(
                    "Records should start with '{}' character!".format(cls.RECORD_MARKER))

            #
            comment = line[1:].rstrip()  # Remove the first ">".
            sequence_string = handle.readline()
            structure_annotation_string = handle.readline()

            # Skip other lines till the "next" ">"
            line = handle.readline()
            while True:
                if not line:
                    break
                if line[0] == cls.RECORD_MARKER:
                    break
                line = handle.readline()

            # Remove trailing whitespace
            yield comment, sequence_string.strip(), structure_annotation_string.strip()

            if not line:
                return  # StopIteration

    # endregion
