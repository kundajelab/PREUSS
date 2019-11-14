# -*- coding: utf-8 -*-

"""
Project Helper
================

Helper Functions for project setup and running.
"""

from __future__ import print_function

import os
import sys
import subprocess

from distutils import spawn


# ----------------------------------
# region Constants

CODE_DIRECTORY = 'neoRNA'
DOCS_DIRECTORY = 'docs'
TESTS_DIRECTORY = 'tests'
PYTEST_FLAGS = ['--doctest-modules']

# endregion


# ----------------------------------
# region General

def read_file(filename):
    r"""
    Return the contents of a file.

    Parameters
    ----------
    filename: str
        File path.

    Returns
    -------
    file_content: str
        The file's content
    """
    with open(os.path.join(os.path.dirname(__file__), filename)) as infile:
        return infile.read()


def get_project_files():
    r"""
    Retrieve a list of project files, ignoring hidden files.

    Returns
    -------
    file_list: list
        Sorted list of project files

    """

    if is_git_project() and has_git():
        return get_git_project_files()

    project_files = []
    for top, subdirs, files in os.walk('.'):
        for subdir in subdirs:
            if subdir.startswith('.'):
                subdirs.remove(subdir)

        for f in files:
            if f.startswith('.'):
                continue
            project_files.append(os.path.join(top, f))

    return project_files

# endregion


# ----------------------------------
# region Git

def is_git_project():
    return os.path.isdir('.git')


def has_git():
    return bool(spawn.find_executable("git"))


def get_git_project_files():
    r"""
    Retrieve a list of all non-ignored files, including untracked files, excluding deleted files.

    Returns
    -------
    file_list: list
        Sorted list of git project files
    """

    cached_and_untracked_files = git_ls_files(
        '--cached',  # All files cached in the index
        '--others',  # Untracked files
        # Exclude untracked files that would be excluded by .gitignore, etc.
        '--exclude-standard')
    uncommitted_deleted_files = git_ls_files('--deleted')

    # Since sorting of files in a set is arbitrary, return a sorted list to
    # provide a well-defined order to tools like flake8, etc.
    return sorted(cached_and_untracked_files - uncommitted_deleted_files)


def git_ls_files(*cmd_args):
    r"""
    Run ``git ls-files`` in the top-level project directory. Arguments go directly to execution call.

    Parameters
    ----------
    cmd_args

    Returns
    -------

    """

    cmd = ['git', 'ls-files']
    cmd.extend(cmd_args)
    return set(subprocess.check_output(cmd).splitlines())

# endregion


# ----------------------------------
# region Tests

def print_success_message(message):
    r"""
    Print a message indicating success in green color to STDOUT.

    Parameters
    ----------
    message

    Returns
    -------

    """

    try:
        import colorama
        print(colorama.Fore.GREEN + message + colorama.Fore.RESET)
    except ImportError:
        print(message)


def print_failure_message(message):
    r"""
    Print a message indicating failure in red color to STDERR.

    Parameters
    ----------
    message

    Returns
    -------

    """

    try:
        import colorama
        print(colorama.Fore.RED + message + colorama.Fore.RESET,
              file=sys.stderr)
    except ImportError:
        print(message, file=sys.stderr)


def test():
    r"""
    Run the unit tests.

    Returns
    -------

    """

    # Make sure to import pytest in this function. For the reason, see here:
    # <http://pytest.org/latest/goodpractises.html#integration-with-setuptools-test-commands>  # NOPEP8
    import pytest

    # This runs the unit tests.
    # It also runs doctest, but only on the modules in TESTS_DIRECTORY.
    return pytest.main(PYTEST_FLAGS + [TESTS_DIRECTORY])


def lint():
    """Run lint and return an exit code."""
    # Flake8 doesn't have an easy way to run checks using a Python function, so
    # just fork off another process to do it.

    # Python 3 compat:
    # - The result of subprocess call outputs are byte strings, meaning we need
    #   to pass a byte string to endswith.
    project_python_files = [filename for filename in get_project_files()
                            if filename.endswith(b'.py')]

    return_code = subprocess.call(
        ['flake8', '--max-complexity=10'] + project_python_files)

    if return_code == 0:
        print_success_message('No style errors')

    return return_code

# endregion




