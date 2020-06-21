import os
import pytest

from neoRNA.sequence.base_pair import BasePair

parametrize = pytest.mark.parametrize


class TestBasePair(object):

    # Sample sequences
    nt_pair = ('A', 'U')
    pos_pair = (56, 22)

    def test_base_pair_init(self):
        bp = BasePair(self.nt_pair, self.pos_pair)
        assert bp is not None

    def test_base_pair_position(self):
        bp = BasePair(self.nt_pair, self.pos_pair)

        assert bp.left_position is 56
        assert bp.left_nt is 'A'
        assert bp.right_position is 22
        assert bp.right_nt is 'U'
