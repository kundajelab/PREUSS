import os

import pytest

from neoRNA.io.library_def_io import LibraryIO

parametrize = pytest.mark.parametrize


class TestLibraryIO(object):
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/io/example_files/rlib_example.rlib'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_load_file(self):
        #
        with open(self.__EXAMPLE_FILE_PATH) as handle:
            assert handle is not None

    def test_parse_record(self):
        #
        records = []
        with open(self.__EXAMPLE_FILE_PATH) as handle:
            for record in LibraryIO.parse_iterator(handle):
                records.append(record)

        #
        assert len(records) == 2
        rna_lib_item = records[0]
        assert str(rna_lib_item.barcode) == 'GGTGCCGGT'
        assert str(rna_lib_item.note) == 'This is a note'
