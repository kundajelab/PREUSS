# -*- coding: utf-8 -*-

"""
Python Scripts for `neoRNA` package
================

"""

import os
import yaml
import logging.config

script_dir = os.path.dirname(os.path.realpath(__file__))


# ----------------------------------------------------------------
# Load "logging" config
#
def setup_logging(
    config_path='logging.yaml',
    logging_level=logging.INFO,
):
    r"""
    Setup logging configuration
    """
    path = os.path.join(script_dir, config_path)
    if os.path.exists(path):
        with open(path, 'rt') as infile:
            config = yaml.safe_load(infile.read())
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=logging_level)
