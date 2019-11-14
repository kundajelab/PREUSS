
import os
import pytest

from neoRNA.library.mut_count_old import MutCountOld

parametrize = pytest.mark.parametrize

class TestMutCountOld(object):


    def test_Init(self):
        start = 1
        end = 11
        mutation_string = '~s~||A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record is not None

    def test_InvalidRecord(self):
        start = 1
        end = 100
        mutation_string = 's~s||A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record.is_valid is False

    def test_InvalidRecordByNonOverlap(self):
        start = 1
        end = 11
        mutation_string = '~s~||A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record.is_valid is True

        mutation_string = '~s|~|A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record.is_valid is False

    def test_MutationStringAdjust(self):
        start = 1
        end = 11
        mutation_string = 's~s||A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record.is_valid is True
        assert record.start_position == 4
        assert record.end_position == 9
        assert record.mutation_string == '||A|G|'

    def test_SimpleMutationString(self):
        start = 1
        end = 11
        mutation_string = 's~s||A|G|s~'
        record = MutCountOld(start, end, mutation_string)
        assert record.is_valid is True
        assert record.start_position == 4
        assert record.end_position == 9
        assert record.simple_mutation_string() == '001010'