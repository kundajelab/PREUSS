# -*- coding: utf-8 -*-

"""
Matlab Runners
================

It defines several runners that utilize MATLAB.
"""

import os
import subprocess

from typing import List

# MATLAB Engine
# Ref - https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html
import matlab.engine


class MatlabRunner(object):
    r"""
    Matlab Runners

    It utlizes "MATLAB API Engine for Python" (shipped with MATLAB app) to help invoke "MatLab" functions.

    Current Matlab Functions:
    - Biers
        - rna_structure

    ## NOTE
    - Be sure to config "MatLab" and "MatLab Engine" correctly.
    - For the "MatLab" functions, be sure to add them in the "path" of MATLAB, so that it can be run directly.

    ## Links
    - MATLAB API for Python - https://www.mathworks.com/help/matlab/matlab-engine-for-python.html
    - Biers - https://github.com/ribokit/Biers
    """

    # ----------------------------------
    # region Engine

    MATLAB_ENGINE = None

    @classmethod
    def start_engine(cls):
        r"""
        Start the engine.
        """

        if cls.MATLAB_ENGINE:
            return cls.MATLAB_ENGINE

        try:
            cls.MATLAB_ENGINE = matlab.engine.start_matlab()
        except OSError:
            raise ValueError("Can't start MATLAB Engine.\n"
                             "Please be sure to set up it correctly - "
                             "https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html")

        return cls.MATLAB_ENGINE

    @classmethod
    def stop_engine(cls):
        r"""
        Stop the engine.
        """

        if not cls.MATLAB_ENGINE:
            return

        try:
            cls.MATLAB_ENGINE.quit()
        except OSError:
            raise ValueError("Can't stop MATLAB Engine.")

    # endregion

    # ----------------------------------
    # region Biers

    @classmethod
    def rna_structure(cls, sequence_str: str, position_offset: int, position_range,
                      reactivity_1d: List[float] = None, max_bootstrap: int = 0):
        r"""
        Run `rna_structure` from Biers package (Matlab package)

        Parameters
        ----------
        sequence_str: str
            The sequence, in string.
        position_offset: int
            The offset of the sequence. Normally used when the "sequence" to be checked is a partial sequence.
        position_range: range
            The range of sequence position, used with "position_offset".
        reactivity_1d: List[float]
            1-dimensional reactivity data, in DMS format.
        max_bootstrap: int
            Max rounds of bootstrap.

        Returns
        -------
        structure: str
            Predicted structure string, in "Dot-Bracket" format.
        bpp: matrix
            2-D array for base pair probability (bpp) data.
        bootstrap: List(str)
            The list of bootstrap structure strings.
        """

        try:
            # NOTES
            # - Use `matlab.int8()` or `matlab.double()` to pass the "numeric array"
            #   instead of "cell array" if passing "list"
            # - Convert all "int" value to "float" (double) since MATLAB uses "double" here.....
            # - Add `nargout=#` to determine the # of returns values
            structure, bpp, bootstrap \
                = cls.MATLAB_ENGINE.rna_structure(sequence_str, matlab.double([]),
                                                  float(position_offset),
                                                  matlab.uint8(position_range),
                                                  matlab.double([]),
                                                  float(max_bootstrap),
                                                  0,  # Use "Fold" method
                                                  matlab.double(reactivity_1d),
                                                  nargout=3)
        except OSError:
            raise OSError("Error when predicting RNA structure for sequence - {}".format(sequence_str))

        return structure, bpp, bootstrap

    # endregion

