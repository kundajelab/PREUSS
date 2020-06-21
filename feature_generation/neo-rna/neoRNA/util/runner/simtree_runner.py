# -*- coding: utf-8 -*-

"""
SimTree Runners
================

The running wrapper for "Simtree.jar"
"""

import os
import sys
import subprocess

import re

from typing import List

from neoRNA.util.file_utils import FileUtils


class SimTreeRunner(object):
    r"""
    SimTree Runner

    This runner is based on the "jar" paackage of SimTree.

    Ref: http://bioinfo.cs.technion.ac.il/SimTree/
    """

    # CMD Template
    # 4 arguments - jar path, options, structure #1, structure #2
    BASE_CMD_STR_TEMPLATE = 'java -Xmx128m -jar {} {} -structures {} {}'
    # SIM_TREE_JAR_PATH = '/Users/cowfox/Desktop/_rna_lib_analysis/__tool/SimTree_v1.2.3.jar'

    # Default options used for "VARNA.jar"
    DEFAULT_OPTIONS = {
        'details': 'yes',
        'mapping': 'yes',
        'flip': 4,
    }

    # Default highlight region options
    DEFAULT_HIGHLIGHT_REGION_OPTIONS = {
        'fill': "#C55337",
        'outline': '#FFFFFF',
    }

    # Default annotation options
    DEFAULT_ANNOTATION_OPTIONS = {
        'type': "B",
        'anchor': "66",
        'outline': '#000000',
        'size': '6',
    }

    # ----------------------------------
    # region Init

    def __init__(self, simtree_location: str, cwd: str):
        r"""
        Init

        Parameters
        ----------
        simtree_location: str
            The "path" to the SimTree jar package.
        cwd: str
            The "path" to the "Working Folder".

        """

        self.simtree_location = simtree_location
        self.simtree_options = self.DEFAULT_OPTIONS

        #
        self.cwd = cwd

    # endregion

    # ----------------------------------
    # region Methods - RNA Structure Compare

    def compare_rna_structure(self, rna_structure_1: str, rna_structure_2: str):
        r"""
        Compare "two" RNA structures by using the SimTree method.

        Parameters
        ----------
        rna_structure_1: str
        rna_structure_2: str

        Returns
        -------

        """

        #
        temp_structure_1_file_path = os.path.join(self.cwd, 'rna_structure_1.txt')
        temp_structure_2_file_path = os.path.join(self.cwd, 'rna_structure_2.txt')
        FileUtils.save_file(temp_structure_1_file_path, '{}\n'.format(rna_structure_1))
        FileUtils.save_file(temp_structure_2_file_path, '{}\n'.format(rna_structure_2))

        #
        # Build CMD
        cmd = self.__gen_cmd_str(self.simtree_options, temp_structure_1_file_path, temp_structure_2_file_path)

        # Call
        output = None
        try:
            print(cmd)
            output = subprocess.check_output([cmd], shell=True).decode(sys.stdout.encoding).strip()
        except OSError as error:
            print('CMD called with error - ', error.args)
            return None, None

        #
        normalized_score, flipping_modes = self.retrieve_info(output)

        #
        return normalized_score, flipping_modes

    # endregion

    # ----------------------------------
    # region Methods - Utils

    def retrieve_info(self, raw_output):
        #
        # Regex format to match the string
        regex = "Normalized Score:\s*(?P<normalized_score>[\w.]+)\s*Flipping Nodes: (?P<flipping_modes>[\w.]+)"
        rx_info = re.compile(regex, re.MULTILINE)
        matched = rx_info.findall(raw_output.strip())

        if not matched and len(matched) == 0:
            return None

        #
        return matched[0]

    # endregion

    # ----------------------------------
    # region Internal Methods - CMD Building

    def __gen_cmd_str(self, params,
                      rna_structure_file_path_1: str, rna_structure_file_path_2: str):
        r"""
        Generate "CMD string" for SimTree package.

        Parameters
        ----------
        params: dict
            The parameter object.
        rna_structure_file_path_1: str
        rna_structure_file_path_2: str

        Returns
        -------
        cmd_str: str
            The CMD string.
        """

        #
        params_str = self.__gen_param_str(params, use_quote=False)

        #
        return self.BASE_CMD_STR_TEMPLATE.format(self.simtree_location, params_str,
                                                 rna_structure_file_path_1, rna_structure_file_path_2)

    def __gen_param_str(self, params: dict,
                        prefix: str = '-', equal_sign: str = ' ', separator: str = ' ',
                        use_quote: bool = True) -> str:
        r"""
        Generate parameter string.

        Each element in the "parameter" dict will be string.

        Parameters
        ----------
        params: dict
            The parameter object.
        prefix: str
            The "prefix" used when generate "each" parameter pair.
        equal_sign: str
            The "equal sign" used between the "key" and "value".
        separator: str
            The "separator string" between each parameter item.
        use_quote: bool
            If use "quote" on the "value".

        Returns
        -------
        param_str: str
            The parameter string.
        """

        if not params:
            return ''

        param_str_items = []
        for param, value in params.items():
            #
            str_template = '{}{}{}"{}"' if use_quote else '{}{}{}{}'

            item_str = str_template.format(prefix, param, equal_sign, value)
            param_str_items.append(item_str)

        #
        return separator.join(param_str_items)

    # endregion
