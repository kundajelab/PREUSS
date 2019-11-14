# -*- coding: utf-8 -*-

"""
DNA/RNA sequence
--------------------
"""

from typing import Tuple, Optional, List
from Bio.Seq import Seq


class Sequence(object):
    r"""
    An Object definition of DNA/RNA sequence.

    The sequence object is a wrapper based on "Seq" object, with some additional attributes and properties, such as:
    - Position Range - It is useful when the sequence position is "NOT" from "1".

    """

    # ----------------------------------
    # region Init

    def __init__(self, sequence_str: str, position_range: range = None):
        r"""
        Init

        Parameters
        ----------
        sequence_str: str
            The raw sequence string.
        position_range: range
            The position range. Optional
        """

        #
        self.__sequence: Seq = Seq(sequence_str) if sequence_str else None
        self.__position_range: range = position_range

        # Check sequence type
        self.__is_rna_sequence: bool = True
        self.__check_sequence_type()

        # Length Info
        self.__length = None
        self.__start_position = None
        self.__end_position = None
        self.calculate_length()

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def sequence(self) -> Optional[Seq]:
        return self.__sequence

    @property
    def sequence_str(self) -> str:
        return str(self.__sequence) if self.__sequence else None

    @property
    def position_range(self) -> range:
        return self.__position_range

    @property
    def length(self) -> Optional[int]:
        return self.__length

    @property
    def start_position(self) -> Optional[int]:
        return self.__start_position

    @property
    def end_position(self) -> Optional[int]:
        return self.__end_position

    # endregion

    # ----------------------------------
    # region Methods - Sequence

    def is_empty(self) -> bool:
        r"""
        Check if the "sequence" us empty.

        Returns
        -------

        """

        #
        return self.__sequence is None

    def get_reverse_complement(self) -> Optional[Seq]:
        r"""
        Get the "reverse complement" sequence.

        Returns
        -------
        reverse_complement: Optional[Seq]
            The "reverse complement" sequence.
        """

        return self.__sequence.reverse_complement() if self.__sequence else None

    def get_rna_sequence(self) -> Optional[Seq]:
        r"""
        Get the "RNA type" of sequence.

        Depending on its current "sequence type", it will return the `sequence` itself
        or a `transcribed` version.

        Returns
        -------
        rna_sequence: Optional[Seq]
            The sequence object, in "RNA" type.
        """

        #
        if self.__is_rna_sequence:
            return self.__sequence
        else:
            # A DNA sequence
            return self.__sequence.transcribe() if self.__sequence else None

    def get_dna_sequence(self) -> Optional[Seq]:
        r"""
        Get the "DNS type" of sequence.

        Depending on its current "sequence type", it will return the `sequence` itself
        or a `transcribed` version.

        Returns
        -------
        rna_sequence: Optional[Seq]
            The sequence object, in "DNA" type.
        """

        if self.__is_rna_sequence:
            return self.__sequence.back_transcribe() if self.__sequence else None
        else:
            # A DNA sequence
            return self.__sequence

    # endregion

    # ----------------------------------
    # region Methods - nt

    def get_nt(self, nt_position: int) -> Optional[str]:
        r"""
        Get the "nt" based on the given "nt position".

        Parameters
        ----------
        nt_position: int
            The "position" of nt. Should be consistent with "Start / End" position of this sequence.

        Returns
        -------
        ret_nt: Optional[str]
            The returned "nt". "None" if not available.
        """

        #
        if self.contain(nt_position):
            #
            return self.sequence_str[nt_position - self.start_position]

        return None

    def is_nt_ac(self, nt_position: int) -> bool:
        r"""
        Check if the given "nt position" is "A" or "C".

        Parameters
        ----------
        nt_position: int
            The "position" of nt.

        Returns
        -------
        is_ac: bool
            If it is "A" or "C".

        """

        #
        if self.contain(nt_position):
            #
            return self.sequence_str[nt_position - 1] == 'A' or self.sequence_str[nt_position - 1] == 'C'

        return False

    def nt_ordering(self, nt_position: int) -> Optional[int]:
        r"""
        Determine the "local ordering" of the given "nt position".

        NOTE:
        - The ordering starts from "1".

        Parameters
        ----------
        nt_position: int
            The "position" of nt.

        Returns
        -------

        """

        #
        if self.contain(nt_position):
            #
            return nt_position - self.start_position + 1

        return None

    def get_nt_by_ordering(self, local_ordering: int) -> Tuple[Optional[int], Optional[str]]:
        r"""
        Get the "nt" based on the given "local ordering".

        NOTE:
        - Normally, an ordering should start from "1".
        - If an ordering is negative, it is from the "LAST" nt of the sequence.

        Parameters
        ----------
        local_ordering: int

        Returns
        -------

        """

        #
        if local_ordering >= 1:
            #
            position = self.start_position + local_ordering - 1
            return position, self.get_nt(position)

        if local_ordering < 0:
            #
            position = self.end_position + local_ordering + 1
            return position, self.get_nt(position)

        # No found
        return None, None

    # endregion

    # ----------------------------------
    # region Methods - Sequence-Position Comparison

    def contain(self, nt_position: int) -> bool:
        r"""
        Check if the sequence "contains" the given "nt position".

        Parameters
        ----------
        nt_position: int
            The given "nt position" to check.

        Returns
        -------
        if_contain: bool
            If the sequence "contains" the given "nt position".

        """

        #
        if not nt_position or not self.__start_position or not self.__end_position:
            return False

        return self.__start_position <= int(nt_position) <= self.__end_position

    def is_after(self, nt_position: int) -> bool:
        r"""
        Check if the sequence is "after" the given "nt position".

        Parameters
        ----------
        nt_position: int
            The given "nt position" to check.

        Returns
        -------
        if_is_after: bool
            If the sequence "is after" the given "nt position".
        """

        #
        if not nt_position or not self.__start_position or not self.__end_position:
            return False

        return self.__start_position > int(nt_position)

    def is_before(self, nt_position: int) -> bool:
        r"""
        Check if the sequence is "before" the given "nt position".

        Parameters
        ----------
        nt_position: int
            The given "nt position" to check.

        Returns
        -------
        if_is_before: bool
            If the sequence "is before" the given "nt position".

        """

        #
        if not nt_position or not self.__start_position or not self.__end_position:
            return False

        return self.__end_position < int(nt_position)

    def distance(self, nt_position: int, is_before: bool = True):
        r"""
        Calc the "distance" by a given "nt position".

        NOTE:
            - If it is in a "before" mode, use "end-nt"
            - Otherwise, use "start-nt"

        Parameters
        ----------
        nt_position: int
            The given "nt position" to check.
        is_before: bool
            If it is to calc in a "before" mode.

        Returns
        -------

        """

        #
        if not nt_position or not self.__start_position or not self.__end_position:
            return None

        #
        if is_before:
            if not self.__end_position:
                return None
            return self.__end_position - nt_position
        else:
            # After
            if not self.__start_position:
                return None
            return self.__start_position - nt_position

    # endregion

    # ----------------------------------
    # region Methods - Mutations

    def generate_mutation_syntax(self, base_sequence_string: str,
                                 seq_type_rna: bool = False,
                                 sequence_start: int = 1, sequence_end: Optional[int] = None,
                                 sequence_position_counter_offset: int = 0
                                 ) -> Tuple[Optional[List[int]], Optional[List[str]]]:
        r"""
        Generate the "mutation syntax" based on a given "base" sequence.

        The "mutation syntax" follows the format: [position][base_nt]to[mut_nt]

        Example:
        `2AtoG` - at nt position `2`, there is a mutation from "A" to "G".
        `2AtoG, 12GtoU` - two mutations happen on the sequence

        Parameters
        ----------
        base_sequence_string: str
            The "Base" sequence to be compared.
        seq_type_rna: bool
            If using "RNA" sequence.
        sequence_start: int
            The "start position" of the sequence to be inferred. Must be no less than then "1".
            Default to "1" - from the beginning of the sequence.
        sequence_end: Optional[int]
            The "end position" of the sequence to be inferred.
            It can accept a "negative" number which means to trim the "last #" nt.
            Default to "None" - till the "end" of the sequence.
        sequence_position_counter_offset: int
            The "position counter offset" of the sequence to be inferred.
            Default to "0" - no offset, start the counter from "1".

        Returns
        -------
        mutation_positions: List[int]
            The list of "positions" of mutation.
        mutation_syntax: List[str]
            The list of "syntax" of mutation.


        """
        #
        if not self.__sequence:
            return None, None

        # Get sequence
        if seq_type_rna:
            sequence_string = str(self.get_rna_sequence())
        else:
            sequence_string = str(self.get_dna_sequence())

        # Ignore the generation, if the `length` is not the same
        if len(sequence_string) != len(base_sequence_string):
            return None, None

        # Get the sub-sequence based on "trimming info"
        # sequence_slice = slice(trim_left, len(base_sequence_string) - trim_right)
        sequence_slice = slice(sequence_start - 1, sequence_end)
        sequence_string = sequence_string[sequence_slice]
        base_sequence_string_trimmed = base_sequence_string[sequence_slice]

        mut_positions = []
        mut_syntax = []

        # Loop into the "base sequence" and keep record of each "difference"
        for idx, nt in enumerate(base_sequence_string_trimmed):
            if nt != sequence_string[idx]:
                mut_syntax_str = '{}{}to{}'.format(idx + 1 + sequence_position_counter_offset,
                                                   nt, sequence_string[idx])
                mut_positions.append(idx + 1 + sequence_position_counter_offset)
                mut_syntax.append(mut_syntax_str)

        # Return all syntax as one string
        return mut_positions, mut_syntax

    # endregion

    # ----------------------------------
    # region Internal Methods

    def calculate_length(self) -> None:
        r"""
        Calculate the "length" info, including length, start position, end position

        NOTE:
        - If "position range" not manually assigned, it will need to calc it here.
        - The default "position range" starts from "1".
        - `start position` and `end position` are based on `position_range`.

        Returns
        -------

        """

        self.__length = len(self.__sequence) if self.__sequence else None

        if self.__position_range is None:
            # Use default way to determine the range, start from 1
            self.__position_range = range(1, self.__length + 1)

        # Calculate "start" and "end" position of a range
        # Python 2.7, range() returns `list`.
        if type(self.__position_range) is range and len(self.__position_range) > 0:
            self.__start_position = self.__position_range[0]
            self.__end_position = self.__start_position + self.__length - 1

    def __check_sequence_type(self) -> None:
        r"""
        Check the "sequence Type" of the property `__sequence` - DNA or RNA.

        When checking, it only needs to check if the sequence contains `U` or `T`

        - `U` -> RNA Type
        - `T` -> DNA Type

        Returns
        -------
        """
        #
        if not self.__sequence:
            self.__is_rna_sequence = False
            return

        #
        if self.__sequence.upper().find('U') >= 0:
            self.__is_rna_sequence = True
        else:
            self.__is_rna_sequence = False

    # endregion

    # ----------------------------------
    # region Methods - Magic ones

    def __str__(self):
        """
        Return a string format of the sequence

        :return:
        """
        if not self.__sequence:
            return ''

        # Example: (48,50) CAA
        return '{} ({}, {})'.format(str(self.__sequence), self.__start_position, self.__end_position)

    # endregion
