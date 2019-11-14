# -*- coding: utf-8 -*-

"""
RNA Library
================
"""

import json
from typing import List

from neoRNA.library.library_config import RnaLibConfig
from neoRNA.library.library_item import LibraryItem


class RnaLibrary(object):
    """
    RNA library Object.

    It normally represents a run of "RNA Lib Items" - a collection of data for each of the items.

    """

    # ----------------------------------
    # region Init

    def __init__(self, rna_items: List[LibraryItem] = None):
        r"""
        Init

        Parameters
        ----------
        rna_items: List[LibraryItem]
            The list of RNA Lib items.
        """

        self.__rna_items = rna_items

        # Meta, `public` access
        self.data_source_code = None
        self.data_source_title = None
        self.data_source_date = None
        self.data_source_description = None

        self.running_code = None
        self.running_date = None
        self.running_notes = None

        # RNA ID of "WT"
        self.wide_type_rna_id = None
        self.wide_type_rna_sequence = None

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def rna_items(self) -> List[LibraryItem]:
        return self.__rna_items

    @rna_items.setter
    def rna_items(self, items: List[LibraryItem]) -> None:
        self.__rna_items = items

    # endregion

    # ----------------------------------
    # region Methods - Meta

    def load_meta(self, configs: RnaLibConfig) -> None:
        r"""
        Load "meta" from RNA Lib config file. 
        
        Parameters
        ----------
        configs: RnaLibConfig
            The RNA Lib config object

        Returns
        -------

        """

        #
        self.data_source_code = configs['data_source_code']
        self.data_source_title = configs['data_source_title']
        self.data_source_date = configs['data_source_date']
        self.data_source_description = configs['data_source_description']

        self.running_code = configs['running_code']
        self.running_date = configs['running_date']
        self.running_notes = configs['running_notes']

        self.wide_type_rna_id = configs['wide_type_rna_id']
        self.wide_type_rna_sequence = configs['wide_type_rna_sequence']

    # endregion

    # ----------------------------------
    # region JSON Object

    # endregion
