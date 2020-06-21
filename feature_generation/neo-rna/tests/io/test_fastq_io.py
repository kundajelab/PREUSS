# -*- coding: utf-8 -*-

import os
import pytest

from neoRNA.io.fastq_io import FastqIO

parametrize = pytest.mark.parametrize

class TestFastqIO(object):
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/fastq_example.fastq'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_fileLoad(self):
        parser = FastqIO(self.__EXAMPLE_FILE_PATH)
        assert parser is not None

    def test_loadSequence(self):
        parser = FastqIO(self.__EXAMPLE_FILE_PATH)
        assert parser.total() == 2
        assert parser.sequences() is not None

    def test_readSequenceObject(self):
        parser = FastqIO(self.__EXAMPLE_FILE_PATH)
        sequences = parser.sequences()

        assert len(sequences) == 2

        first = sequences[0]

        # "id", "name", "description" should all refer to the part after "@".
        assert first.id == 'NS500418:AACCAGAT'
        assert first.name == 'NS500418:AACCAGAT'
        assert first.description == 'NS500418:AACCAGAT'

        # The sequence should be "all" upper case
        assert str(first.seq) == 'CTTTANAATATAGATCTTGTTG'
