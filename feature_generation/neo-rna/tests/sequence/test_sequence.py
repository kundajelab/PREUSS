import os
import pytest

from neoRNA.sequence.sequence import Sequence

parametrize = pytest.mark.parametrize


class TestSequence(object):

    # Sample sequences
    sequence = 'CUCUGAUCUCU'
    sequence_range = range(12, 22)

    def test_sequence_init(self):
        sequence = Sequence(self.sequence)
        assert sequence is not None

    def test_sequence_init_with_range(self):
        sequence = Sequence(self.sequence, self.sequence_range)
        assert sequence is not None

    def test_sequence_length(self):
        sequence = Sequence(self.sequence)

        assert type(range(1, 3)) is list  # Python 2.7, range() returns `list`.
        assert type(sequence.position_range) is list
        assert sequence.length is 11
        assert sequence.start_position is 1
        assert sequence.end_position is 11

    def test_sequence_length_with_range(self):
        sequence = Sequence(self.sequence, self.sequence_range)

        assert sequence.length is 11
        assert sequence.start_position is 12
        assert sequence.end_position is 22
