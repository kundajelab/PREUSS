# -*- coding: utf-8 -*-

"""
RNAfold Runners
================

The running wrapper for "RNAfold" tool.
"""

import os
import subprocess

import re
from typing import List, Tuple, Optional

from neoRNA.util.file_utils import FileUtils


class RnaFoldRunner(object):
    r"""
    RNAfold Runner

    This runner is for "RNAfold" tool.

    Ref: https://www.tbi.univie.ac.at/RNA/RNAfold.1.html
    """

    # CMD template
    BASE_CMD_STR_TEMPLATE = 'RNAfold {flags} < {input_file} > {output_file}'

    # Temp files
    DEFAULT_SEQUENCE_FILE = 'tmp_rna_fold_input.txt'
    DEFAULT_OUTPUT_FILE = 'tmp_rna_fold_output.txt'

    # ----------------------------------
    # region Init

    def __init__(self):
        r"""
        Init

        Parameters
        ----------
        """

    # endregion

    # ----------------------------------
    # region Methods - Free Energy Extraction

    def extract_free_energy(self, sequence: str, constraint: str = None,
                            flags: List[str] = None) -> Tuple[Optional[float], Optional[str]]:
        r"""
        Extract the "free energy" from the RNAfold results.

        Parameters
        ----------
        sequence: str
            The input sequence.
        constraint: str
            The "constraint" - a string containing constraints on the structure encoded.
        flags: List[str]
            The "manual" flag list.

        Returns
        -------
        fe_pair: Tuple[float, str]
            The result pair of "FE".
        """

        #
        flag_list = flags if flags else []
        content_list = [sequence]

        #
        if not constraint:
            flag_list.append('-p')  # Calculate the partition function and base pairing probability matrix.
        else:
            flag_list.append('-C')
            content_list.append(constraint)

        #
        FileUtils.save_file(self.DEFAULT_SEQUENCE_FILE, '\n'.join(content_list))

        # Call
        cmd = self.BASE_CMD_STR_TEMPLATE.format(flags=' '.join(flag_list),
                                                input_file=self.DEFAULT_SEQUENCE_FILE,
                                                output_file=self.DEFAULT_OUTPUT_FILE)
        try:
            # print(cmd)
            subprocess.call(cmd, shell=True, stdout=subprocess.PIPE)
        except OSError as error:
            print('RNAfold command running error - ', error.args)

        #
        result_content_str = FileUtils.load_file_as_str(self.DEFAULT_OUTPUT_FILE, no_newline=False)

        # Parse the value from the contents
        if not constraint:
            # Ensemble FE only
            return self.__parse_ensemble_fe(result_content_str)
        else:
            # MFE only
            return self.__parse_mfe(result_content_str)

    # endregion

    # ----------------------------------
    # region Utils

    def cleanup(self):
        #
        os.remove(self.DEFAULT_SEQUENCE_FILE) if os.path.exists(self.DEFAULT_SEQUENCE_FILE) else None
        os.remove(self.DEFAULT_OUTPUT_FILE) if os.path.exists(self.DEFAULT_SEQUENCE_FILE) else None

    # endregion

    # ----------------------------------
    # region Internal Methods

    def __parse_mfe(self, content_str: str) -> Tuple[Optional[float], Optional[str]]:
        r"""
        Parse the "MFE" (Minimum Free Energy) from the content str.

        NOTE:
            - The "MFE" is inside the "()".

        Parameters
        ----------
        content_str: str
            The content str.

        Returns
        -------
        mfe_pair: Tuple[float, str]
            The result pair of "MFE".
        """
        #
        regex_mfe = '(?P<dot_bracket>[^a-zA-Z\s]+)\s+\((?P<value>[-.0-9]+)\)'

        #
        for ret_dict in self.__regex_find(content_str, regex_mfe):
            # Only pick the "FIRST" one
            if ret_dict:
                return float(ret_dict['value']), ret_dict['dot_bracket']
            else:
                return None, None

    def __parse_ensemble_fe(self, content_str: str) -> Tuple[Optional[float], Optional[str]]:
        r"""
        Parse the "Ensemble FE" (Free Energy) from the content str.

        NOTE:
            - The "Ensemble FE" is inside the "[]".

        Parameters
        ----------
        content_str: str
            The content str.

        Returns
        -------
        fe_pair: Tuple[float, str]
            The result pair of "FE".
        """
        #
        regex_fe = '(?P<dot_bracket>[^a-zA-Z\s]+)\s+\[(?P<value>[-.0-9]+)\]'

        #
        for ret_dict in self.__regex_find(content_str, regex_fe):
            # Only pick the "FIRST" one
            if ret_dict:
                return float(ret_dict['value']), ret_dict['dot_bracket']
            else:
                return None, None

    # endregion

    # ----------------------------------
    # region CLass Methods

    @classmethod
    def __regex_find(cls, content_str: str, regex: str):
        r"""
        Find "matching" results based on given "regex".

        Parameters
        ----------
        content_str: str
        regex: str

        Returns
        -------
        """

        #
        matches = re.finditer(regex, content_str, re.MULTILINE)
        for match_num, match in enumerate(matches):
            #
            yield match.groupdict()

    # endregion
