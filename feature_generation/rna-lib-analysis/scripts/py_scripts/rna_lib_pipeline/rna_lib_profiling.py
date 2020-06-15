#!/usr/bin/env python

# ================================================================
# RNA Lib Pipeline - Profiling
#
# The "Entry" script for "RNA Lib" pipeline - to run "profiling" on the RNA lib.
#
#
# ## Running Steps
# - QC, by AfterQC
# - Demultiplexing, by Novobarcode
# - Profiling, by ShapeMapper 2.x
#
# ## Inputs
# - RNA Lib "config" file
#
# ## Outputs
# - RNA Lib profiling file
#   - It is a "bin" file which encodes the "results object".
#
#

import os
import argparse

import json

# Add "py_scripts" into module path, relative to "current" script
import sys
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.realpath(__file__))))
sys.path.append(local_module_path)

import logging
from py_scripts import setup_logging

from pyppl import PyPPL, Proc, Channel

from neoRNA import io
from neoRNA.library.library_config import RnaLibConfig
from neoRNA.util.json_serializable import PythonObjectEncoder


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Pipeline - Profiling')

# Inputs
arguments_parser.add_argument('config',
                              metavar='rna_lib_config_file',
                              help='The file path to RNA Lib config file.')


# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
config_file_path = os.path.abspath(args.config)

# Validation
if not os.path.exists(config_file_path):
    raise ValueError('"RNA Lib" config file does not exist.')

# endregion


# ----------------------------------------------------------------
# region Set Logger
#

setup_logging(logging_level=logging.INFO)
# Logger Name
logger_name = 'neo_rna.script'
logger = logging.getLogger(logger_name)

# endregion


# ----------------------------------
# region Prep - Config

# Load the config
configs = RnaLibConfig(config_file_path)

# Validation
if not configs:
    raise ValueError('"RNA Lib" configs load error. Please double check. ')

# Decide the working folder
configs.decide_working_folder(os.getcwd())
if not os.path.exists(configs.working_folder):
    os.mkdir(configs.working_folder)
os.chdir(configs.working_folder)

# endregion


# ----------------------------------
# region Pipeline - Prep
#
#
# The tool "PyPPL" cannot accept "python object" as input, we will need to convert it to
#  "string" and then convert it back while using it.

# Load RNA Items
if not configs['rna_lib_file']:
    raise ValueError('"RNA Lib" file load error. Please double check. ')

rna_lib_items = []
# A list of "RNA Lib items", each of which is in "python object string" format.
# It is used to pass via "Channel".
rna_lib_item_object_json_list = []
for rna_item in io.parse(configs['rna_lib_file'], "rna-lib-def"):
    rna_lib_items.append(rna_item)
    # Convert each item object to "json string"
    rna_lib_item_object_json_list.append(json.dumps(rna_item, cls=PythonObjectEncoder))

# Convert the objects to "string" and pass it as "argument"
configs_object_json = json.dumps(configs, cls=PythonObjectEncoder)
# A "python object string" format of a list, which includes a list of "RNA Lib items".
# It is used to passed as "whole"
rna_lib_items_list_json = json.dumps(rna_lib_items, cls=PythonObjectEncoder)

# endregion

# ----------------------------------
# region Pipeline - AfterQC

# endregion

# ----------------------------------
# region Pipeline - Novobarcode

# endregion


# ----------------------------------
# region Pipeline - ShapeMapper

# ShapeMapper output folder path
shape_output_folder = os.path.join(configs.working_folder, 'shapemapper_results')

pShape = Proc(desc='Run ShapeMapper 2.x')
pShape.input = {"rna_item_json:var": Channel.create(rna_lib_item_object_json_list)}
# Define the "output" channel - the "output folder"
pShape.output = "shape_output_folder:var: {}".format(shape_output_folder)
pShape.forks = 4

#
pShape.args.configs = configs_object_json
pShape.lang = 'python'
pShape.script = """
#!/usr/bin/env python

from neoRNA.library.shape_mapper.shape_runner import ShapeMapperRunner
ShapeMapperRunner.shape_mapper_v2({{args.configs | squote}}, {{in.rna_item_json | squote}})
"""
# endregion


# ----------------------------------
# region Pipeline - Process RNA Lib Running Result

# RNA Lib output file path
rna_lib_file_path = os.path.join(configs.working_folder,
                                 'data_{}_run_{}_at_{}.rbin'.format(
                                     configs['data_source_code'],
                                     configs['running_code'],
                                     configs['running_date']))

pRnaLib = Proc(desc='Process RNA Lib Running Result from ShapeMapper 2.x')

#
if configs['run_shapemapper'] is True:
    pRnaLib.depends = pShape
    # automatically inferred from pShape.output
    # Collapse the channel into "1".
    pRnaLib.input = {"shape_output_folder:var": lambda ch: ch.collapse(col=0)}

pRnaLib.output = "rna_lib_output_file:var: {}".format(rna_lib_file_path)
pRnaLib.forks = 1

#
pRnaLib.args.configs = configs_object_json
pRnaLib.args.rna_items = rna_lib_items_list_json
pRnaLib.lang = 'python'
pRnaLib.script = """
#!/usr/bin/env python

from py_scripts.rna_lib_pipeline.script.proc_rna_lib_profiling_results import generate_rna_lib_profiling_results
generate_rna_lib_profiling_results({{args.configs | squote}}, {{args.rna_items | squote}}, {{out.rna_lib_output_file | squote}})
"""

# endregion


# ----------------------------------
# region Run Pipelines

pyppl_config = {
    'proc': {
        'echo': True,  # Output all
        'errhow': 'ignore',  # Ignore the error job
        'errntry': 3
    }
}

if configs['run_shapemapper'] is True:
    PyPPL(pyppl_config).start(pShape).run()
else:
    # Only run "RNA Lib" process
    PyPPL(pyppl_config).start(pRnaLib).run()

# endregion

