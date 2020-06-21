# -*- coding: utf-8 -*-

"""
neoRNA
================

A Python Toolkit for RNA NGS and Analysis.
"""

from neoRNA import metadata

# ----------------------------------
# region Meta Info

__version__ = metadata.version
__build__ = 0x021901
__author__ = metadata.authors[0]
__license__ = metadata.license
__copyright__ = metadata.copyright

# endregion


# ----------------------------------
# region Logging

# Set default logging handler to avoid "No handler found" warnings.
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

# endregion

