# -*- coding: utf-8 -*-

"""
IO - RNA Library Definition
================
"""

import re

from neoRNA.library.library_item import LibraryItem


class LibraryDefinitionIO(object):
    r"""
    IO to parse "RNA library" definition file, return as a list of RNA library Items (in entry level).

    The file follows the following format:

    - "Comment line" start with "#"
    - Each line represents "1" RNA library Items.
    - For each item, it could have multiple attributes, separated by `\t`.
    - Currently it has 4 attributes:
        - RNA ID
        - RNA Barcode
        - RNA sequence
        - Note

    Example:
        #RNA ID	BC_Seq	"RNA Sequence (-30, right)"	Note
        001	GGTGCCGGT	GGGAGCCTGCCCTCTGATCTCTGCCTGTTC  "Sample notes"
        002	GGTGCCGGT	GGGAGCCTGCCCTCTGATCTCTGCCTGTTC  "Sample notes"

    """

    # The "marker" to help indicate a "start"
    STARTER_MARKER = '#'

    # Number of attributes
    NUM_ATTRIBUTES = 4
    # The delimiter used to split the content
    DELIMITERS = '\t'

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
        parsed_objects: Library
            Parsed objects.


        Usage
        -------
        >>> with open("rna_lib.rlib") as handle:
        ...     for record in LibraryIO.parse_iterator(handle):
        ...         print(record)
        ...

        """

        for rna_id, rna_barcode_string, rna_sequence_string, notes in cls.parse(handle):
            yield LibraryItem(rna_id, rna_barcode_string, rna_sequence_string, notes)

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
            A tuple of strings (rna_id, rna_barcode, rna_sequence, notes).

        """

        # Skip any text before the first record (e.g. blank lines, comments)
        while True:
            line = handle.readline()
            if line == "":
                return
            if line[0] == cls.STARTER_MARKER:  # Find the "first line" of a record.
                break

        while True:
            if line[0] != cls.STARTER_MARKER:
                raise ValueError(
                    "Records should start with '{}' character!".format(cls.STARTER_MARKER))

            #
            header = line[1:].rstrip()  # Remove the first "marker".

            # Loop in
            line = handle.readline()
            while line:
                # Parse it
                __delimiters = cls.DELIMITERS
                __max_attributes = cls.NUM_ATTRIBUTES
                parts = re.split(__delimiters, line.strip())

                # The line has to have the first THREE attributes
                if not len(parts) >= 3:
                    raise ValueError('The line must have at least 3 attributes', line)

                # Add `None` to missing attributes
                parts += [None] * (__max_attributes - len(parts))

                rna_id, rna_barcode_string, rna_sequence_string, notes = parts

                # Remove trailing whitespace
                yield rna_id, rna_barcode_string.strip(), rna_sequence_string.strip(), \
                    notes.strip() if notes else None

                #
                line = handle.readline()

            if not line:
                return  # StopIteration

    # endregion

