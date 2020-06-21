# -*- coding: utf-8 -*-

import os
import pytest

from Bio.Seq import Seq
import numpy as numpy
from neoRNA.io.mut_count_io import MutCountIO

parametrize = pytest.mark.parametrize


class TestMutRateIO(object):

    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/mutation_rate_example.csv'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    __SEQUENCE = Seq('GGGAGCCTGCCCTCTGATCTCTGCCTGTTCCTCTGTCCCACAGAGGGCAAAGGCTACGGGTCAGAGAGCGGGGAGGAGGACGGTGCCGGTTTCG')

    def test_init(self):
        mut_rate_record = MutCountIO(self.__SEQUENCE)
        assert mut_rate_record is not None

    def test_initWithRatesList(self):
        mut_rate_record = MutCountIO(self.__SEQUENCE)
        rates_list = [0.1, 0.2, 0.3]
        mut_rate_record.rates = rates_list
        assert mut_rate_record is not None
        assert mut_rate_record.rates.size == 3

    def test_initWithFileLoad(self):
        mut_rate_record = MutCountIO(self.__SEQUENCE, self.__EXAMPLE_FILE_PATH)
        assert mut_rate_record is not None

    def test_ratesParsing(self):
        mut_rate_record = MutCountIO(self.__SEQUENCE, self.__EXAMPLE_FILE_PATH)
        assert len(str(self.__SEQUENCE)) == mut_rate_record.rates.size

    def test_ratesFixing(self):
        mut_rate_record = MutCountIO(self.__SEQUENCE, self.__EXAMPLE_FILE_PATH)
        rates = mut_rate_record.rates
        assert rates[0] == -999
        fixed_rates = mut_rate_record.get_fixed_rates()
        assert fixed_rates[0] == 0

    def test_ACOnlyRates(self):
        sequence_string = 'ACGTATG'
        rates_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, -999]
        mut_rate_record = MutCountIO(Seq(sequence_string))
        mut_rate_record.rates = rates_list

        assert mut_rate_record.get_ac_only_rates().tolist() == [0.1, 0.2, -999.0, -999.0, 0.5, -999.0, -999.0]
