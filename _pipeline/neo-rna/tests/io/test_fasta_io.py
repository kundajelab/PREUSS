# -*- coding: utf-8 -*-

import os
import pytest

from neoRNA.io.fasta_io import FastaIO

parametrize = pytest.mark.parametrize


class TestFastaIO(object):
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/fasta_example.fa'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_fileLoad(self):
        parser = FastaIO(self.__EXAMPLE_FILE_PATH)
        assert parser is not None

    def test_loadSequence(self):
        parser = FastaIO(self.__EXAMPLE_FILE_PATH)
        assert parser.total() == 2
        assert parser.first_sequence() is not None

    def test_readSequenceObject(self):
        parser = FastaIO(self.__EXAMPLE_FILE_PATH)
        first = parser.first_sequence()

        # "id", "name", "description" should all refer to the part after ">".
        assert first.id == 'Sequence_ID_1'
        assert first.name == 'Sequence_ID_1'
        assert first.description == 'Sequence_ID_1'

        # The sequence should be "all" upper case
        assert str(first.seq) == 'GGGAGCCTGCCCTCTGA'
