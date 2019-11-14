
import os
import pytest

from neoRNA.library.library_item import Library

parametrize = pytest.mark.parametrize

class TestRNALib(object):

    rna_lib_line = '001\tGGTGCCGGT\tGGGAGCCTGCCCTCTGATCTCTGCCTGTTC\tThis is a note'
    rna_lib_line_missing_attribute = '001\tGGTGCCGGT\tGGGAGCCTGCCCTCTGATCTCTGCCTGTTC'
    rna_lib_line_with_error = '001\tGGTGCCGGT\t'

    def test_RNALibItemInit(self):
        rna_lib_item = Library(self.rna_lib_line)
        assert rna_lib_item is not None

    def test_RNALibItemInitWithMissingAttributes(self):
        rna_lib_item = Library(self.rna_lib_line_missing_attribute)
        assert rna_lib_item is not None

    def test_RNALibItemInitWithError(self):
        with pytest.raises(ValueError) as excepton_info:
            Library(self.rna_lib_line_with_error)

    def test_RNALibItemProperty(self):
        rna_lib_item = Library(self.rna_lib_line)
        assert rna_lib_item.id == '001'

    def test_RNALibItemPropertyOfSequence(self):
        rna_lib_item = Library(self.rna_lib_line)
        assert str(rna_lib_item.barcode) == 'GGTGCCGGT'

    def test_SequenceTypeCheck(self):
        rna_lib_item = Library(self.rna_lib_line)
        assert rna_lib_item.get_dna_sequence() == 'GGGAGCCTGCCCTCTGATCTCTGCCTGTTC'
        assert rna_lib_item.get_rna_sequence() == 'GGGAGCCUGCCCUCUGAUCUCUGCCUGUUC'

    def test_RNASequence(self):
        new_line = '001\tGGTGCCGGT\tGGGAGCCU'
        rna_lib_item = Library(new_line)
        assert rna_lib_item.get_dna_sequence() == 'GGGAGCCT'
        assert rna_lib_item.get_rna_sequence() == 'GGGAGCCU'

    def test_MutationSyntax(self):
        new_line = '001\tGGTGCCGGT\tGGGAGCCU'
        rna_lib_item = Library(new_line)
        base_sequence_string = 'GGGCGCCU'
        assert rna_lib_item.generate_mutation_syntax(base_sequence_string) == '4CtoA'
        base_sequence_string_2 = 'GGGCGACU'
        assert rna_lib_item.generate_mutation_syntax(base_sequence_string_2) == '4CtoA,6AtoC'
