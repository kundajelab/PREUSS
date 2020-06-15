# -*- coding: utf-8 -*-

"""
Barcode Sequence
================
"""

from typing import Tuple


class Barcode(object):
    r"""
    Barcode Sequence

    Normally, a barcode is used in the "demultiplexing" process, to help split the a huge sequence file into
    several individual one - each of which contains all the sequences under a specific barcode.

    The barcode can appear at the beginning of a sequence or at the end of a sequence.
    Sometimes, it may not be at the "very" beginning / end, but just inside a part of beginning / end.
    In order to skip the "unneeded" part, use "random base indicator" - "N".

    For example, if a barcode starts from the "3rd" position and has a length of "4" (like "ATCG"),
    then the barcode will be "NNATCG"

    NOTE: When including `N` in the sequence, the "actual barcode" should be the in "ONE" part which does not include
    "N". In other word, string like `NNATGCNACG` will not be working, since it has two parts of "non-N" sequences.

    """

    # ----------------------------------
    # region Init

    def __init__(self, barcode_str: str, reference_id: str = None):
        r"""
        Init

        Parameters
        ----------
        barcode_str: str
            The barcode string.
        reference_id: str
            The reference ID, usually "RNA ID".
        """

        #
        self.__raw_str: str = barcode_str
        self.__reference_id: str = reference_id

        # Extract the "actual" barcode
        self.__actual_barcode, self.__barcode_slice = self.__extract_barcode(self.__raw_str)

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def raw_str(self) -> str:
        return self.__raw_str

    @property
    def reference_id(self) -> str:
        return self.__reference_id

    @property
    def barcode(self) -> str:
        return self.__actual_barcode

    @property
    def barcode_slice(self):
        r"""
        Return the positions of the actual barcode, in a `slice` format.
        """
        return self.__barcode_slice

    # endregion

    # ----------------------------------
    # region Private Methods

    @classmethod
    def __extract_barcode(cls, raw_str) -> Tuple[str, slice]:
        r"""
        Extract the "actual" barcode from the given raw string.

        For example, if the raw string is `NNATCG` then the actual barcode is `ATCG` (without `N`s).

        Parameters
        ----------
        raw_str: str
            The barcode raw string.

        Returns
        -------
        barcode: str
            The actual barcode string.
        barcode_slice: slice
            The "slice" of the barcode string.
        """

        import re

        actual_barcode = None
        barcode_slice = None

        # The regex used for finding the "non-N" part.
        regex_string = "[^N]+"

        # Do the matches
        matches = re.search(regex_string, raw_str)
        first_match = matches.group()

        if first_match is not '':
            actual_barcode = first_match

            # Get the "start" / "end" position of the first match
            start_position = matches.start()
            end_position = matches.end()
            barcode_slice = slice(start_position, end_position)

        return actual_barcode, barcode_slice

    # endregion

