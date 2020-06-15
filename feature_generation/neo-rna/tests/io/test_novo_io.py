# -*- coding: utf-8 -*-

import os
import pytest

from Bio.Seq import Seq
import numpy as numpy

from neoRNA.io.novo_io import NovoIO

parametrize = pytest.mark.parametrize


class TestNovoIO(object):

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_BARCODE_FILENAME = 'tests/io/example_files/novobarcode_example.result'
    __EXAMPLE_BARCODE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_BARCODE_FILENAME)

    def test_parseNovobarcodeResults(self):
        barcode_results = NovoIO.parse_novobarcode_results(self.__EXAMPLE_BARCODE_FILE_PATH)
        assert barcode_results is not None
        assert barcode_results['CCGCCGCGT'] == 4960
        assert barcode_results['NC'] == 10860
        assert barcode_results['barcode'] == 33181
        assert barcode_results['total'] == 44041
