
import os
import pytest

from neoRNA.io.mut_count_io import MutCountIO
from neoRNA.library.reactivity import Reactivity

from Bio.Seq import Seq
import numpy as numpy

parametrize = pytest.mark.parametrize

class TestReactivity(object):

    def test_1DCalculationv2(self):
        sequence = Seq('AACG')

        mod_rates = [0.0, 0.0, 0.2, 0.5]
        non_mod_rates = [0.0, 0.2, 0.0, 0.2]
        mod_rates = numpy.array(mod_rates, dtype=float)
        non_mod_rates = numpy.array(non_mod_rates, dtype=float)

        mod_rates_record = MutCountIO(sequence)
        mod_rates_record.rates = mod_rates
        non_mod_rates_record = MutCountIO(sequence)
        non_mod_rates_record.rates = non_mod_rates

        reactivity = Reactivity(sequence, mod_rates_record, non_mod_rates_record)
        reactivity.calculate_1d_reactivity(2)
        results = reactivity.reactivity_1d

        assert len(results) == 4
        assert results[0] == -999.0


