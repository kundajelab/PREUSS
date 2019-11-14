#!/usr/bin/env python

# ================================================================
# RNA Lib Pipeline - Structure Inferring
#
# The script loads the "profiling results" and performs "RNA Secondary Structure" inferring.
#
# ## Running Steps
# - Structure Inferring, by Biers
# - Structure Interpretation, by bpRNA
#
# ## Inputs
# - RNA Lib profiling file
#   - The "python object" file to contain all "RNA Lib" item objects which include the "profiling info".
# - Editing Level file
#   - A `csv` file to include the "editing info", such "editing position", "editing level".
#
# ## Outputs
# - RNA Lib structure file
#   - The "python object" file to contain all "RNA Lib" item objects which include the "structure info", including both
#       "structure inferring info" and "structure element info".
# - RNA Lib structure summary file
#   - A "human readable" JSON file to help include the "structure inferring info" by Biers.
#

import os
import argparse

import json
import csv

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

from neoRNA.io.library_io import LibraryIO
from neoRNA.library.rna_library import RnaLibrary
from neoRNA.util.json_serializable import PythonObjectEncoder

from py_scripts.rna_lib_pipeline.script.proc_biers_rna_structure \
    import biers_results_folder_path, biers_inference_structure_folder_path

# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Pipeline - Structure Inferring')

# Inputs
arguments_parser.add_argument('rna_lib_profiling',
                              metavar='rna_lib_profiling.rbin',
                              help='The file path to RNA Lib profiling result file.')
arguments_parser.add_argument('--editing_level',
                              metavar='editing_level_file',
                              help='The file path to the Editing Level file')

# Sequence Slice
arguments_parser.add_argument('--start', dest='sequence_start',
                              action='store', default=1,
                              help='Start point of the sequence. Default to "1" (from beginning).')
arguments_parser.add_argument('--end', dest='sequence_end',
                              action='store', default=None,
                              help='End point of the sequence. Default to "None" (till end).')

#
arguments_parser.add_argument('--biers_bootstrap',
                              action='store', default=20,
                              help='The total of bootstrap runs.')
arguments_parser.add_argument('--indel',
                              action="store", default='None',
                              metavar='indel_rna_id_range',
                              help='The rna id range for "indel" cases')

# Output
arguments_parser.add_argument('--structure_out',
                              action='store', default='rna_lib_structure.rbin',
                              help='Filename of RNA Lib Structure file.')
arguments_parser.add_argument('--structure_summary_out',
                              action='store', default='rna_lib-structure_summary.json',
                              help='Filename of RNA Lib Structure Summary file.')

# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
profiling_file_path = os.path.abspath(args.rna_lib_profiling)
editing_level_file_path = os.path.abspath(args.editing_level)

# Validation
if not os.path.exists(profiling_file_path):
    raise ValueError('"RNA Lib" profiling file does not exist.')
if editing_level_file_path and not os.path.exists(profiling_file_path):
    raise ValueError('"Editing Level" file does not exist. ')

# Values
sequence_start = int(args.sequence_start) if args.sequence_start is not None else 1
sequence_end = int(args.sequence_end) if args.sequence_end is not None else None
biers_max_bootstrap = int(args.biers_bootstrap) if args.biers_bootstrap is not None else 20
indel_rna_id_range = args.indel

# Output
structure_file = args.structure_out
structure_summary_file = args.structure_summary_out

# Internal Flags
# -- If need to run Biers
biers_run = True
# -- If need to override results
biers_override = True
# -- If run bpRNA
bprna_run = True
# -- If only run summary step
summary_only = False

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
# region Working Folder

# Check the "main" output file
# - If it is "absolute path", use its folder as "working folder".
# - If not, use "current folder" as "working folder"

if os.path.isabs(structure_file):
    #
    os.chdir(os.path.dirname(structure_file))
else:
    #
    cwd = os.getcwd()
    structure_file = os.path.join(cwd, structure_file)

#
cwd = os.getcwd()

#
if not os.path.isabs(structure_summary_file):
    #
    structure_summary_file = os.path.join(cwd, structure_summary_file)

# endregion


# ----------------------------------
# region Prep - RNA Lib Object

rna_library: RnaLibrary = LibraryIO.as_python_object(profiling_file_path)
#
rna_library_object_json = json.dumps(rna_library, cls=PythonObjectEncoder)

# endregion


# ----------------------------------
# region Prep - Editing Level

# Editing Level values & positions, indexed by "RNA ID"
editing_levels = {}
editing_positions = {}

# Columns to load
rna_id_column_name = 'RNA_ID_STR'
editing_level_column_name = 'Avg'
editing_position_column_name = 'Editing_position'

#
editing_level_data = csv.DictReader(open(editing_level_file_path, 'rU'))
for line in editing_level_data:
    rna_id = line[rna_id_column_name]
    #
    score = line[editing_level_column_name]
    editing_levels[rna_id] = float(score) if score != "NA" and score != "#N/A" else None
    #
    position = line[editing_position_column_name]
    editing_positions[rna_id] = position

# endregion


# ----------------------------------
# region Pipeline - Prep

# A list of "RNA Lib items", each of which is in "python object string" format.
# It is used to pass via "Channel".
rna_lib_item_object_json_list = []
rna_id_list = []
for rna_item in rna_library.rna_items:
    rna_id_list.append(rna_item.rna_id)
    # Convert it to "json string"
    rna_lib_item_object_json_list.append(json.dumps(rna_item, cls=PythonObjectEncoder))

# FOR TESTING - Only run a few
# rna_lib_item_object_json_list = rna_lib_item_object_json_list[slice(0, 2)]
# rna_id_list = rna_id_list[slice(0, 2)]

# endregion


# ----------------------------------
# region Pipeline - Biers

# Biers output folder path
biers_inference_structure_folder = biers_inference_structure_folder_path(cwd)

pBiers = Proc(desc='Run Biers RNA Structure Inferring.')
pBiers.input = {
    "rna_item_json:var": Channel.create(rna_lib_item_object_json_list),
    "rna_id:var": Channel.create(rna_id_list)
}
# Define the "output" channel - the RNA ID
pBiers.output = "rna_id:var:{{in.rna_id}}"
pBiers.forks = 1  # MatLab ONLY does "1" thread....

#
pBiers.args.local_module_path = local_module_path
pBiers.args.working_folder = cwd
pBiers.args.sequence_start = sequence_start
pBiers.args.sequence_end = sequence_end
pBiers.args.max_bootstrap = biers_max_bootstrap
pBiers.args.override = biers_override
pBiers.lang = 'python'
pBiers.script = """
#!/usr/bin/env python

import os, sys
sys.path.append({{args.local_module_path | squote}})

from py_scripts.rna_lib_pipeline.script.proc_biers_rna_structure import biers_rna_structure
biers_rna_structure({{in.rna_item_json | squote}}, {{args.working_folder | squote}}, override={{args.override}}, sequence_start={{args.sequence_start}}, sequence_end={{args.sequence_end}}, max_bootstrap={{args.max_bootstrap}})
"""
# endregion


# ----------------------------------
# region Pipeline - bpRNA

pBpRNA = Proc(desc='Run bpRNA to interpret RNA secondary structure.')
if bprna_run and biers_run:
    pBpRNA.depends = pBiers
    pBpRNA.input = "rna_id:var"
else:
    pBpRNA.input = {
        "rna_id:var": Channel.create(rna_id_list)
    }
pBpRNA.output = "rna_id:var:{{in.rna_id}}"
pBpRNA.forks = 1

#
pBpRNA.args.working_folder = cwd
pBpRNA.args.biers_inference_structure_folder_path = biers_inference_structure_folder
pBpRNA.args.bprna_results_folder_name = 'bpRNA_results'
pBpRNA.script = """
bpRNA.pl {{args.biers_inference_structure_folder_path}}/{{in.rna_id}}.dbn

mkdir -p {{args.working_folder}}/{{args.bprna_results_folder_name}}
mv {{in.rna_id}}.st {{args.working_folder}}/{{args.bprna_results_folder_name}}
"""

# endregion


# ----------------------------------
# region Pipeline - Summary

biers_results_folder = biers_results_folder_path(cwd)
bprna_results_folder = os.path.join(cwd, 'bpRNA_results')

#
pSummary = Proc(desc='Process Structure Summary Results.')

#
if summary_only is False:
    if bprna_run:
        # It has dependency
        pSummary.depends = pBpRNA
        # automatically inferred from pShape.output
        # Collapse the channel into "1".
        pSummary.input = {"rna_id:var": lambda ch: ch.collapse(col=0)}
    if biers_run:
        pSummary.depends = pBiers
        pSummary.input = {"rna_id:var": lambda ch: ch.collapse(col=0)}

pSummary.output = "summary_file:var: {}".format(structure_summary_file)
pSummary.forks = 1

#
pSummary.args.local_module_path = local_module_path

pSummary.args.rna_library = rna_library_object_json

pSummary.args.sequence_start = sequence_start
pSummary.args.sequence_end = sequence_end

pSummary.args.editing_level_file_path = editing_level_file_path
pSummary.args.biers_results_folder_path = biers_results_folder
pSummary.args.bprna_results_folder_path = bprna_results_folder

pSummary.args.structure_summary_file_path = structure_summary_file
pSummary.lang = 'python'
pSummary.script = """
#!/usr/bin/env python

import os, sys
sys.path.append({{args.local_module_path | squote}})

from py_scripts.rna_lib_pipeline.script.proc_rna_lib_structure_results import generate_rna_lib_structure_results
generate_rna_lib_structure_results({{args.rna_library | squote}}, {{args.editing_level_file_path | squote}}, {{args.biers_results_folder_path | squote}}, {{args.bprna_results_folder_path | squote}}, {{args.structure_summary_file_path | squote}}, {{args.sequence_start}}, {{args.sequence_end}})
"""

# endregion


# ----------------------------------
# region Pipeline - Run

# Config
pyppl_config = {
    'proc': {
        'echo': 'stderr',  # Output all stderr
        'errhow': 'ignore',  # Ignore the error job
        'errntry': 3
    }
}

if summary_only:
    PyPPL(pyppl_config).start(pSummary).run()
elif biers_run:
    PyPPL(pyppl_config).start(pBiers).run()
else:
    PyPPL(pyppl_config).start(pBpRNA).run()

# endregion


# ----------------------------------
# region TEST

# Biers
# from py_scripts.rna_lib_pipeline.script.proc_biers_rna_structure import biers_rna_structure
# biers_rna_structure(rna_lib_item_object_json_list[1], cwd, sequence_start=14, max_bootstrap=20)

# Summary
# from py_scripts.rna_lib_pipeline.script.proc_rna_lib_structure_results import generate_rna_lib_structure_results
# generate_rna_lib_structure_results(
#     rna_library_object_json, editing_level_file_path, biers_results_folder, bprna_results_folder, structure_summary_file,
#     sequence_start, sequence_end)

# endregion



