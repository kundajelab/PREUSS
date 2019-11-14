# -*- coding: utf-8 -*-

from Bio import SeqIO

from neoRNA.util.file_utils import FileUtils


class FastqIO(object):
    """
    io interface for "FASTQ" file format - sequence File.

    In the file, each of the sequence contains "4" lines:
    - Line 1 begins with a '@' character and is followed by a sequence identifier and an optional description
    (like a FASTA title line).
    - Line 2 is the raw sequence letters.
    - Line 3 begins with a '+' character and is optionally followed by the same sequence identifier
    (and any description) again.
    - Line 4 encodes the quality values for the sequence in Line 2, and must contain the same number
    of symbols as letters in the sequence.

    It utilizes the module "Bio.SeqIO" to operate with sequence file.

    :link: https://en.wikipedia.org/wiki/FASTQ_format
    """

    # File type parameter, used by `SeqIO`
    __FILE_TYPE = 'fastq'

    # Batch size - max # of records read from the file for each batch
    __BATCH_SIZE = 100000

    # region Init
    # -----------------------------------------
    # Init
    # -----------------------------------------
    def __init__(self, filename):
        self.__filename = filename

        self.__sequence_list = None
        self.__sequence_batch_iterator = None

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
        if not self.__sequence_list:
            # Build the sequence list
            self.__sequence_list = self.__parse_sequence_file(self.__filename)

        return self.__sequence_list

    @property
    def sequence_batch_iterator(self):
        """
        Return the sequence batch iterator.

        :return: The sequence batch iterator.
        """
        if not self.__sequence_batch_iterator:
            # Build the sequence batch iterator
            self.__sequence_batch_iterator = self.__parse_sequence_file_with_batch(self.__filename)

        return self.__sequence_batch_iterator

    def total(self):
        """
        Return the total number of the sequences in the file.

        This is only available when the file is fully parsed. For "batched" parsing, it returns `None`.

        :return: Total number of sequences, or `None` for batched parsing.
        """
        if self.__sequence_list:
            return len(self.__sequence_list)

        return None

    # endregion

    # region Private Methods
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

    def __parse_sequence_file_with_batch(self, filename, batch_size=None):
        """
        Parse the "large" file and return a list of batches.

        :param filename: The filename of the sequence file
        :return: The list of batches.
        """
        # Decide the batch size. Use default one if not provided.
        size = batch_size if batch_size else self.__BATCH_SIZE

        self.__sequence_batch_iterator = FileUtils.batch_iterator(SeqIO.parse(open(filename), self.__FILE_TYPE), size)

        return self.__sequence_batch_iterator

    # endregion
