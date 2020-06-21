# -*- coding: utf-8 -*-

"""
The io interfaces to parse and operate on different types of "files".

For each of the "file type", it provides:

- Basic "parser" and "writer" interface with the file.
- Additional utilities.

"""

from __future__ import print_function

from typing import Any

from neoRNA.util.file_utils import FileUtils

# IOes
from .fasta_io import FastaIO
from .library_def_io import LibraryDefinitionIO

from .shape_profile_io import ShapeProfileIO
from .shape_reactivity_io import ShapeReactivityIO

from .dot_bracket_io import DotBracketIO
from .bpp_io import BasePairProbabilityIO
from .bp_rna_io import BpRnaIO


# ----------------------------------
# region File Type to IO Mapping

_FileTypeToIO = {
    # Sequence
    "fasta": FastaIO,
    "dot-bracket": DotBracketIO,

    # RNA Lib
    "rna-lib-def": LibraryDefinitionIO,

    #
    "bpp": BasePairProbabilityIO,

    # Tool - ShapeMapper 2
    "shape-profile": ShapeProfileIO,
    "shape-reactivity": ShapeReactivityIO,

    # Tool - bpRNA
    "bp-rna": BpRnaIO,
}

# endregion

# ----------------------------------
# region Generic Parser


def parse(handle: Any, file_type: str):
    r"""
    Parse a file based on its "file type" and return an iterator.

    Parameters
    ----------
    handle: Any
        Handle to the file, or the filename as a string.
    file_type: str
        The "format" of the file - it loads the specific parser based on the "type".

    Returns
    -------
    iterator: Iterator
        An iterator of "parsing object".


    ## Usage
    -------

    >>> from neoRNA import io
    >>> filename = "rna_001.dbn"
    >>> for record in io.parse(filename, "rna-lib-def"):
    ...    print("Comment %s" % record.comment)

    """

    # Parameter check
    if not file_type:
        raise ValueError("File type required.")
    if file_type not in _FileTypeToIO:
        raise ValueError("File type not exist. Available options: {}".format(" | ".join(list(_FileTypeToIO.keys()))))

    # Handle the input
    mode = 'rU'
    iterator = None
    with FileUtils.as_handle(handle, mode) as fp:
        # Map the file type to a parsing iterator:
        if file_type in _FileTypeToIO:
            io = _FileTypeToIO[file_type]
            iterator = io.parse_iterator(fp)

        # Iterate
        if iterator:
            for record in iterator:
                yield record

# endregion


# ----------------------------------
# region Writer

def write(records: Any, handle: Any, file_type: str):
    r"""
    It writes a set of "records" to a file.

    Parameters
    ----------
    records: Any
        A list (or iterator) of records.
    handle: Any
        Handle to the file, or the filename as a string.
    file_type: str
        The "format" of the file - it loads the specific parser based on the "type".

    Returns
    -------
    num_records: int
        The number of records written.

    """

    # Parameter check
    if not file_type:
        raise ValueError("File type required.")
    if file_type not in _FileTypeToIO:
        raise ValueError(
            "File type not exist. Available options: {}".format(" | ".join(list(_FileTypeToIO.keys()))))

    #
    count = 0
    mode = 'w'

    with FileUtils.as_handle(handle, mode) as fp:
        # Map the file format to a writer
        if file_type in _FileTypeToIO:
            io = _FileTypeToIO[file_type]
            iterator = io.parse_iterator(fp)

    return count

# endregion




