# -*- coding: utf-8 -*-

"""
File Operators
================

Utility for File related operations.
"""

from __future__ import print_function

import codecs
import os
import sys
import contextlib
import itertools

import json
import csv

from typing import Dict, Any, List, Optional

from jinja2 import Environment, FileSystemLoader


class FileUtils(object):
    """
    Utilities for File related operations.
    """

    # ----------------------------------
    # region File Handler

    @classmethod
    @contextlib.contextmanager
    def as_handle(cls, handleish, mode='r', **kwargs):
        r"""Context manager to ensure we are using a handle.

        Context manager for arguments that can be passed for read, write,
        and parse methods: either file objects or path-like objects (strings, pathlib.Path
        instances, or more generally, anything that can be handled by the builtin 'open'
        function).

        When given a path-like object, returns an open file handle to that path, with provided
        mode, which will be closed when the manager exits.

        All other inputs are returned, and are *not* closed.

        Arguments:
         - handleish  - Either a file handle or path-like object (anything which can be
                        passed to the builtin 'open' function: str, bytes, and under
                        Python >= 3.6, pathlib.Path, os.DirEntry)
         - mode       - Mode to open handleish (used only if handleish is a string)
         - kwargs     - Further arguments to pass to open(...)

        Examples
        --------
        >>> with as_handle('seqs.fasta', 'w') as fp:
        ...     fp.write('>test\nACGT')
        >>> fp.closed
        True
        >>> handle = open('seqs.fasta', 'w')
        >>> with as_handle(handle) as fp:
        ...     fp.write('>test\nACGT')
        >>> fp.closed
        False
        >>> fp.close()

        Note that if the mode argument includes U (for universal new lines)
        this will be removed under Python 3 where is is redundant and has
        been deprecated (this happens automatically in text mode).

        """
        # If we're running under a version of Python that supports PEP 519, try
        # to convert `handleish` to a string with `os.fspath`.
        if hasattr(os, 'fspath'):
            try:
                handleish = os.fspath(handleish)
            except TypeError:
                # handleish isn't path-like, and it remains unchanged -- we'll yield it below
                pass

        if isinstance(handleish, str):
            if sys.version_info[0] >= 3 and "U" in mode:
                mode = mode.replace("U", "")
            if 'encoding' in kwargs:
                with codecs.open(handleish, mode, **kwargs) as fp:
                    yield fp
            else:
                with open(handleish, mode, **kwargs) as fp:
                    yield fp
        else:
            yield handleish

    # endregion

    # ----------------------------------
    # region Load File

    @classmethod
    def load_file_as_str(cls, file_path: str, mode: str = 'rU',
                         no_newline: bool = True) -> str:
        r"""
        Load file as a single string.

        Parameters
        ----------
        file_path: str
            The "file path" to be loaded.
        mode: str
            The "mode" when loading the file. Default to "rU".
        no_newline: bool
            If need to remove "newline" char.

        Returns
        -------
        content_str: str
            The loaded content string.

        """

        #
        with cls.as_handle(file_path, mode) as fp:
            #
            content = fp.read()

            #
            return content.replace('\n', '') if no_newline else content

    @classmethod
    def load_csv(cls, file_path: str, mode: str = 'rU',
                 delimiter: str = ',',
                 as_dict: bool = False) -> Any:
        r"""
        Load a "csv" file (also for "tsv" file) and return an object.

        Parameters
        ----------
        file_path: str
            The "file path" to be loaded.
        mode: str
            The "mode" when loading the file. Default to "rU".
        delimiter: str
            The "delimiter" for the data.
        as_dict: bool
            If returns a "dict" object. While in "dict", each row is a dict, which uses "header" as "key".

        Returns
        -------

        """

        #
        with cls.as_handle(file_path, mode) as fp:
            #
            if not as_dict:
                return csv.reader(fp, delimiter=delimiter)
            else:
                return csv.DictReader(fp, dialect='excel-tab')

    # endregion

    # ----------------------------------
    # region Save File

    @classmethod
    def save_file(cls, file_path: str, file_contents: str,
                  mode: str = 'w'):
        """
        Save contents as a file.

        Parameters
        ----------
        file_path: str
            File path to be saved.
        file_contents: str
            The contents. Could be `string`, `json dump`, etc.
        mode: str
            File "open mode", such as "w", "r", "a", etc.

        Returns
        -------

        """
        with open(file_path, mode) as outfile:
            outfile.write(file_contents)

    @classmethod
    def save_json_to_file(cls, file_path: str, json_object: Dict[str, Any]):
        """
        Save JSON object as a file.

        Parameters
        ----------
        file_path: str
            File path to be saved.
        json_object: Dict[str, Any]
            The JSON object to be saved.

        Returns
        -------

        """

        #
        json_contents = json.dumps(json_object,
                                   sort_keys=True,
                                   indent=4, ensure_ascii=True)
        cls.save_file(file_path, json_contents)

    # endregion

    # ----------------------------------
    # region Utility

    @classmethod
    def get_filename(cls, file_path, without_extension=False):
        """
        Retrieve the "filename" from a given "file path".

        You can choose to include / exclude the file extension.

        Parameters
        ----------
        file_path: str
            The given file path.
        without_extension: boolean
            With/without the extension.

        Returns
        -------
        filename: str
            The "filename".
        """
        if without_extension:
            return os.path.splitext(os.path.basename(file_path))[0]
        else:
            return os.path.basename(file_path)

    @classmethod
    def get_dirname(cls, file_path):
        """
        Retrieve the "dir path" from a given "file path".

        Parameters
        ----------
        file_path: str
            The given file path.

        Returns
        -------
        dirname: str
            The path of its "dir".

        """
        return os.path.dirname(file_path)

    @classmethod
    def create_folder(cls, folder_path):
        r"""
        Create the folder if not exists.

        Parameters
        ----------
        folder_path

        Returns
        -------

        """

        #
        if os.path.isabs(folder_path) and not os.path.exists(folder_path):
            #
            os.mkdir(folder_path)

    # endregion

    # ----------------------------------
    # region Template Rendering

    @classmethod
    def render_template(cls, source_folder, template_file, render_contents, target_folder, target_file):
        """
        Render a "template" and then save as a "file".

        The template rendering is based on the package `jinja2`.

        Parameters
        ----------
        source_folder
        template_file
        render_contents
        target_folder
        target_file

        Returns
        -------
        """

        # Load the template file
        jinja2_env = Environment(loader=FileSystemLoader(source_folder), trim_blocks=True)
        template = jinja2_env.get_template(template_file)

        # Render the contents
        output = template.render(render_contents)

        # Save the output as a file
        cls.save_file(os.path.join(target_folder, target_file), output)

    # endregion

    # ----------------------------------
    # region Batch Iterator

    @classmethod
    def batch_iterator(cls, iterator, batch_size):
        """
        A special helper to help separate the iterator into several batches (based on the `batch size`).

        It is very useful to help handle "large" file reading.

        Original Link - http://biopython.org/wiki/Split_large_file

        Parameters
        ----------
        iterator: iterator
            The given iterator.
        batch_size: int
            The size of records in each batch.

        Returns
        -------
        iterators: list
            A list of "batch" iterators.

        """

        entry = True  # Make sure we loop once
        while entry:
            batch = []
            while len(batch) < batch_size:
                try:
                    entry = iterator.next()
                except StopIteration:
                    entry = None

                #
                if entry is None:
                    # End of file
                    break

                #
                batch.append(entry)

            #
            if batch:
                yield batch

    # endregion
