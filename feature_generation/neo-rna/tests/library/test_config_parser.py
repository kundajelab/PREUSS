# -*- coding: utf-8 -*-

import os
import pytest

from neoRNA.util.parser.config_parser import RnaLibConfigParser

parametrize = pytest.mark.parametrize


class TestConfigParser(object):
    fileDir = os.path.dirname(os.path.realpath('__file__'))
    __EXAMPLE_FILENAME = 'tests/util/rna-lib-config-v2-example.yml'
    __EXAMPLE_FILE_PATH = os.path.join(fileDir, __EXAMPLE_FILENAME)

    def test_load_file(self):
        parser = RnaLibConfigParser(self.__EXAMPLE_FILE_PATH)
        assert parser is not None

    def test_init_with_defaults(self):
        #
        parser = RnaLibConfigParser(self.__EXAMPLE_FILE_PATH, defaults={
            'missing_key': 'new_path'
        })
        assert parser.missing_key == 'new_path'

