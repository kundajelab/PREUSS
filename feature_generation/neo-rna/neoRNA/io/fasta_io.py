# -*- coding: utf-8 -*-

from Bio import SeqIO


class FastaIO(object):
    """
    io interface for "FASTA" file format - Nucleotide sequence File.

    It utilizes the module "Bio.SeqIO" to operate with sequence file.

    NOTE: This parser only focus on the "Nucleotide" sequence - "A","T","G","C".

    ## Example "fa" File Format

    The ".fa" file could contain multiple "sequences". For example:

    ```
    >Seq-1
    CGTAACAAG

    >Deq-2
    CGTAACAAG
    ```

    :link: https://en.wikipedia.org/wiki/FASTA_format
    """

    # File type parameter, used by `SeqIO`
    __FILE_TYPE = 'fasta'

    # region Init
    # -----------------------------------------
    # Init
    # -----------------------------------------

    def __init__(self, filename):
        self.__filename = filename

        # Build the sequence list
        self.__sequence_list = self.__parse_sequence_file(filename)

    # endregion

    # region Methods
    # -----------------------------------------
    # Methods
    # -----------------------------------------

    def sequences(self):
        """
        Return the sequence list.

        :return: The parsed sequence list.
        """
        return self.__sequence_list

    def total(self):
        """How many sequences parsed.

        :return: The total number of sequences parsed.
        """
        return len(self.__sequence_list)

    def first_sequence(self):
        """Return the "first" sequence.

        :return: The first sequence.
        """
        return self.__sequence_list[0]

    def last_sequence(self):
        """Return the "last" sequence.

        :return: The last sequence.
        """
        return self.__sequence_list[-1]

    # endregion

    # region PrivateMethods
    # -----------------------------------------
    # Private Methods
    # -----------------------------------------

    def __parse_sequence_file(self, filename):
        """Parse the sequence file and return the "list" of sequences.

        :param filename: The filename of the sequence file
        :return: The list of sequences, in "sequence Object"
        """
        # Transfer the sequence into "Upper Case".
        return list(record.upper() for record in SeqIO.parse(filename, self.__FILE_TYPE))

    # endregion
