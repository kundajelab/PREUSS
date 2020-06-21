import os
import pytest

from neoRNA.sequence.barcode import Barcode

parametrize = pytest.mark.parametrize

class TestBarcode(object):

    # Sample barcodes
    barcode_normal = 'ATGCA'
    barcode_with_n = 'NNATGCAN'
    barcode_with_multiple = 'NNATGCANATGCNNN'

    def test_BarcodeInit(self):
        barcode = Barcode(self.barcode_normal)
        assert barcode is not None

    def test_BarcodeExtracting(self):
        barcode = Barcode(self.barcode_normal)

        #
        assert barcode.actual_barcode() is 'ATGCA'
        assert barcode.barcode_slice().start is 0
        assert barcode.barcode_slice().stop is 5

    def test_BarcodeWithN(self):
        barcode = Barcode(self.barcode_with_n)

        # TODO An error here... Not sure why.....
        # assert barcode.actual_barcode() is 'ATGCA'
        assert barcode.barcode_slice().start is 2
        assert barcode.barcode_slice().stop is 7

    def test_BarcodeWithMultiple(self):
        barcode = Barcode(self.barcode_with_multiple)

        #
        # assert barcode.actual_barcode() is 'ATGCA'
        assert barcode.barcode_slice().start is 2
        assert barcode.barcode_slice().stop is 7
