# -*- coding: utf-8 -*-

"""
IO - RNA Library
================
"""

import json

from neoRNA.library.rna_library import RnaLibrary
from neoRNA.util.file_utils import FileUtils
from neoRNA.util.json_serializable import as_python_object


class LibraryIO(object):
    r"""
    IO to parse "RNA library" results file, return as a list of RNA library Items ("full" version).

    Each of the RNA Lib items contains all the results info.
    """

    # ----------------------------------
    # region Parse as Python object

    @classmethod
    def as_python_object(cls, rna_lib_profiling_file: str) -> RnaLibrary:
        r"""
        Parse the RNA Lib bin file and decode it as a "Python Object".

        Parameters
        ----------
        rna_lib_profiling_file: str
            The "file path" to the RNA Lib bin file.

        Returns
        -------
        parsed_rna_lib: RnaLibrary
            The parsed RNA Lib object.

        """

        object_str = FileUtils.load_file_as_str(rna_lib_profiling_file)
        python_object: RnaLibrary = json.loads(object_str, object_hook=as_python_object)

        #
        return python_object

    # endregion
