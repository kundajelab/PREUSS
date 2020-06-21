# -*- coding: utf-8 -*-

"""
RNA Lib Config Object
================

Parser for Config File, in "YAML" format.
"""

import os
import yaml

from typing import Optional, Dict, Any


class RnaLibConfig(object):
    r"""
    A "YAML" Config File Parser.

    NOTE:
        Check "tests/util" for an example of config file.
    """

    # ----------------------------------
    # region Constants - Config Keys

    #
    KEY_DATE_SOURCE_COMBINED = 'date_source_combined'
    KEY_DATE_SOURCE_BARCODE_AS_FOLDER_NAME = 'data_source_barcode_as_folder_name'

    #
    KEY_WORKING_FOLDER_PATH = 'working_folder_path'
    KEY_OUTPUT_FOLDER_PATH = 'output_folder_path'

    #
    DEFAULT_OUTPUT_FOLDER_NAME = '_rna_lib_results'

    # endregion

    # ----------------------------------
    # region Init

    def __init__(self, config_file: str, defaults: Dict[str, Any] = None):
        r"""
        Init the parser by passing the config file with "optional" defaults values.

        Parameters
        ----------
        config_file: str
            The file path "yaml" config file.
        defaults: Dict[str, Any]
            The default values. Optional.
        """

        #
        if defaults is None:
            defaults = {}

        # Set up the defaults.
        for key, value in defaults.items():
            setattr(self, key, value)

        with open(config_file, 'rU') as infile:
            self.__config: Dict[str, Any] = yaml.safe_load(infile)

        # Override by the config file.
        for key, value in self.__config.items():
            setattr(self, key, value)

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def config(self) -> Dict[str, Any]:
        return self.__config

    @property
    def working_folder(self) -> Optional[str]:
        if self.KEY_WORKING_FOLDER_PATH in self.__config:
            return self.__config[self.KEY_WORKING_FOLDER_PATH]

        return None

    # endregion

    # ----------------------------------
    # region Methods

    def decide_working_folder(self, cwd: str) -> None:
        r"""
        Decide the "working folder".

        In the current pipeline setting,
        - The "working folder" is decided by the attribute - "output_folder_path" inside the config.
        - If "output_folder_path" is not valid, use "cwd" instead.

        Parameters
        ----------
        cwd: str
            "Full" path pf "current working directory".

        Returns
        -------

        """

        #
        if self.__config and self.KEY_OUTPUT_FOLDER_PATH in self.__config:
            output_folder_path = self.__config[self.KEY_OUTPUT_FOLDER_PATH]
            if os.path.isabs(output_folder_path):
                # Working folder is just the "output folder"
                # Add it to "configs"
                self.__config[self.KEY_WORKING_FOLDER_PATH] = output_folder_path
                return

        # Do not have "valid" output folder path, use "cwd" instead
        if os.path.isabs(cwd):
            default_output_folder_path = os.path.join(cwd, self.DEFAULT_OUTPUT_FOLDER_NAME)
            self.__config[self.KEY_OUTPUT_FOLDER_PATH] = default_output_folder_path
            self.__config[self.KEY_WORKING_FOLDER_PATH] = default_output_folder_path

    def is_combined_data_source(self) -> bool:
        r"""
        Check if the "data source" is a combined one.

        For a "combined" dataset, the processing method will be different.

        Returns
        -------
        if_is_combined: bool
            If the "data source" is a combined one.
        """

        #
        if self.KEY_DATE_SOURCE_COMBINED in self.__config:
            return self.__config[self.KEY_DATE_SOURCE_COMBINED]

        #
        return False

    def use_barcode_as_folder_name(self) -> bool:
        r"""
        Check if the "data source" uses "barcode" as the folder name of each RNA Lib item.

        Returns
        -------
        if_use_barcode: bool
            If the "data source" uses "barcode" as the folder name of each RNA Lib item.
        """

        #
        if self.KEY_DATE_SOURCE_BARCODE_AS_FOLDER_NAME in self.__config:
            return self.__config[self.KEY_DATE_SOURCE_BARCODE_AS_FOLDER_NAME]

        #
        return False

    # endregion

    # ----------------------------------
    # region Magic

    def __getitem__(self, item: str):
        r"""
        Try to get the element from the configs.

        Parameters
        ----------
        item: str
            The "key" of the element.

        Returns
        -------
        element: Any
            The value of the element. Return `None` if key is not available.

        """

        #
        if item in self.__config:
            return self.__config[item]

        return None

    # endregion


