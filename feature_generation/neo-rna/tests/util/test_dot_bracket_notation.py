# -*- coding: utf-8 -*-

import os
import pytest

from neoRNA.util.dot_bracket_notation import DotBracketNotation

parametrize = pytest.mark.parametrize


class TestDotBracketNotation(object):

    def test_TypeGenrationv1(self):

        rna_structure = '..((((((((......))))))))..'
        dbn = DotBracketNotation(rna_structure)

        type_string = dbn.generate_type_string_v1(rna_structure)
        assert type_string == '--sssssssshhhhhhssssssss--'

        rna_structure = '..(((((((((...((((((.........))))))........((((((.......))))))..)))))))))...'
        # type_string = dbn.generate_type_string_v1(rna_structure)
        dbn = DotBracketNotation(rna_structure)
        assert dbn.type_notation == '--sssssssss---sssssshhhhhhhhhssssss--------sssssshhhhhhhssssss--sssssssss---'
        assert dbn.stem_length == 21
        assert dbn.hairpin_length == 16

        rna_structure = '.....(((((.(((((((((..((((...(((((........)))))...))))...))))))))).)))))..........'
        type_string = dbn.generate_type_string_v1(rna_structure)
        assert type_string == '-----sssss-sssssssss--ssss---ssssshhhhhhhhsssss---ssss---sssssssss-sssss----------'

        rna_structure = '.....(((((.(((((((((..(((.((((((((.........)))))))))))...))))))))).))))).........'
        dbn = DotBracketNotation(rna_structure)
        assert dbn.type_notation == '-----sssss-sssssssss--sss-sssssssshhhhhhhhhsssssssssss---sssssssss-sssss---------'
        assert dbn.stem_length == 25
        assert dbn.hairpin_length == 9
        assert dbn.hairpin_length_with_focus == 9

