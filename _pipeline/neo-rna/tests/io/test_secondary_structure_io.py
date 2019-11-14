# -*- coding: utf-8 -*-

import os
import pytest

from Bio.Seq import Seq
import numpy as numpy

from neoRNA.io.secondary_structure_io import SecondaryStructureIO

parametrize = pytest.mark.parametrize


class TestSecondaryStructureIO(object):

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/bprna_example.st'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_init(self):
        io = SecondaryStructureIO()
        secondary_structure_info = io.parse(self.__EXAMPLE_FILE_PATH)
        assert secondary_structure_info is not None

    def test_data_parsing(self):
        io = SecondaryStructureIO()
        secondary_structure_info = io.parse(self.__EXAMPLE_FILE_PATH)

        assert len(secondary_structure_info.elements) == 11
