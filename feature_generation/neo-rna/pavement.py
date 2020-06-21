# -*- coding: utf-8 -*-

"""Paver Tasks

:link: https://pythonhosted.org/Paver/
"""

from __future__ import print_function

import os
import sys
import imp

from paver.easy import options, task, needs, consume_args
from paver.setuputils import install_distutils_tasks

# import project_helper as helper
helper = imp.load_source('project_helper', 'project_helper.py')

# Add "current path"
sys.path.append('.')


def print_passed():
    # generated on http://patorjk.com/software/taag/#p=display&f=Small&t=PASSED
    helper.print_success_message(r'''
 ___  _   ___ ___ ___ ___
 | _ \/_\ / __/ __| __|   \
 |  _/ _ \\__ \__ \ _|| |) |
 |_|/_/ \_\___/___/___|___/
''')


def print_failed():
    # generated on http://patorjk.com/software/taag/#p=display&f=Small&t=FAILED
    helper.print_failure_message(r'''
  ___ _   ___ _    ___ ___
 | __/_\ |_ _| |  | __|   \
 | _/ _ \ | || |__| _|| |) |
 |_/_/ \_\___|____|___|___/
''')

# ---------------------------------------------
# Tasks

@task
@needs('doc_html', 'setuptools.command.sdist')
def sdist():
    """Build the HTML docs and the tarball."""
    pass

@task
def test():
    """Run the unit tests."""
    raise SystemExit(helper.test())

@task
def lint():
    # This refuses to format properly when running `paver help' unless
    # this ugliness is used.
    ('Perform PEP8 style check, run PyFlakes, and run McCabe complexity '
     'metrics on the code.')
    raise SystemExit(helper.lint())


