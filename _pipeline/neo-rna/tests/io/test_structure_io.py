# -*- coding: utf-8 -*-

import os
import pytest

from Bio.Seq import Seq
import numpy as numpy

from neoRNA.io.structure_io import StructureIO

parametrize = pytest.mark.parametrize


class TestStructureIO(object):

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/structure_example.react'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_init(self):
        structure_data = StructureIO.parse_react_to_dict(self.__EXAMPLE_FILE_PATH)
        assert structure_data is not None

    def test_dataParsing(self):
        structure_data = StructureIO.parse_react_to_dict(self.__EXAMPLE_FILE_PATH)

        # assert structure_data['date'] == '2017-02-01'
        # assert structure_data['notes'] == 'First Experiment.'
        # assert len(structure_data['items']) == 1
        #
        # rna_item = structure_data['items'][0]
        # assert rna_item['rna_id'] == '001'
        # assert rna_item['barcode_string'] == 'GGTGCCGGT'
        # assert rna_item['reference_structure'] == '.....(((.(((((((.(((((((((((.(((.(((((......))))).)))..))))).)))))).))))).)))))..'

    def test_JSONConversion(self):
        target_file = os.path.join(self.fileDir, 'structure_example.json')
        # Enable it if needed to run the test
        # StructureIO.convert_to_json(self.__EXAMPLE_FILE_PATH, target_file)

