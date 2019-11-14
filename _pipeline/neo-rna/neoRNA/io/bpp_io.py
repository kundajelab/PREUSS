# -*- coding: utf-8 -*-

"""
IO - Base Pair Probability
================

It help process the "Base Pair Probability" format file.
"""

from neoRNA.structure.base_pair_probability import BasePairProbability


class BasePairProbabilityIO(object):
    """
    Parse the "bpp" file.

    Example:
        >81 001 bpp 2 345
        0  	0  	0  	0  	0  	0.973913  	0.9826087  	0.9826087  	0.9884058  	0.9797101  	0.0115942  	0.9826087  	1  	1
        1  	1  	0.9971014  	0.9710145  	0.9768116  	0.8318841  	0.08695652  	0.2753623  	0.9942029  	1  	1
        0.8869565  	0.1884058  	0.09275362  	0.7797101  	0.9826087  	0.9942029  	1  	1  	1  	0.9971014  	0.2173913
        0  	0  	0  	0.2173913  	0.9971014  	1  	1  	1  	0.9942029  	0.9826087  	0.7797101  	0.1884058  	0.08985507
        0.008695652  	0.8869565  	1  	1  	0.9942029  	0.005797101  	0.2753623  	0  	0.8318841  	0.9768116  	0.9710145
        0.9971014  	1  	1  	1  	1  	0.9826087  	0  	0.9797101  	0.9884058  	0.9826087  	0.9826087  	0.973913  	0.002898551
        0.0115942  	0.0115942  	0.0115942  	0  	0  	0  	0  	0

    The format of "dot bracket file":
    - First line - comment line, starting with ">". Usually it contains the "basic info" of the sequence structure.
    - Second line - The list of "bpp", separated by "space" (or "tab")

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
        parsed_objects: BasePairProbabilityIO
            Parsed objects.


        Usage
        -------

        >>> with open("bpp-data.bpp") as handle:
        ...     for record in BasePairProbabilityIO.parse_iterator(handle):
        ...         print(record.comment)
        ...

        """

        for comment, bpp_str in cls.parse(handle):
            yield BasePairProbability(comment, bpp_str)

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
            if line[0] != cls.RECORD_MARKER:
                raise ValueError(
                    "Records should start with '{}' character!".format(cls.RECORD_MARKER))

            #
            comment = line[1:].rstrip()  # Remove the first ">".
            bpp_str = handle.readline()

            # Skip other lines till the "next" ">"
            line = handle.readline()
            while True:
                if not line:
                    break
                if line[0] == cls.RECORD_MARKER:
                    break
                line = handle.readline()

            # Remove trailing whitespace
            yield comment, bpp_str.strip()

            if not line:
                return  # StopIteration

    # endregion
