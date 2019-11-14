#!/usr/bin/env python

# ================================================================
# RNA Lib Analysis Pipeline - Pairwise Analysis for RNA Lib Profiling Replicates
#
# It does "pairwise comparison" on the "profiling" results of two replicates to help exam the "reproducibility".
#
#

import os
import argparse

# Add "py_scripts" into module path, relative to "current" script
import sys
local_module_path = \
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.dirname(
                    os.path.realpath(__file__)))))
sys.path.append(local_module_path)

import logging
from py_scripts import setup_logging

import json

from pyppl import PyPPL, Proc

from neoRNA.library.rna_library import RnaLibrary
from neoRNA.util.file_utils import FileUtils
from neoRNA.util.json_serializable import as_python_object


# ----------------------------------
# region Parsing Argument
#
# :link: https://docs.python.org/dev/library/argparse.html
#

arguments_parser \
    = argparse.ArgumentParser(description='RNA Lib Analysis - Replicate Pairwise Analysis')

# Inputs
arguments_parser.add_argument('--replicate_1',
                              action='store',
                              metavar='rna_lib_profiling_file_replicate_1',
                              help='The file path to RNA Lib profiling file, for replicate #1.')
arguments_parser.add_argument('--replicate_2',
                              action='store',
                              metavar='rna_lib_profiling_file_replicate_2',
                              help='The file path to RNA Lib profiling file, for replicate #2.')

#
arguments_parser.add_argument('--name',
                              action='store', default='dataset',
                              metavar='dataset_name',
                              help='The "name" of the dataset')
arguments_parser.add_argument('--ac_only',
                              action="store_true", default=True,
                              help='If only consider "A"/"C" nt.')

# parse the arguments
args = arguments_parser.parse_args()

# Get the "absolute path" for file / folder
replicate_1_file_path = os.path.abspath(args.replicate_1)
replicate_2_file_path = os.path.abspath(args.replicate_2)

# Validation
if not os.path.exists(replicate_1_file_path):
    raise ValueError('RNA Lib result file (replicate #1) file does not exist.')
if not os.path.exists(replicate_2_file_path):
    raise ValueError('RNA Lib result file (replicate #2) file does not exist.')

#
dataset_name = args.name
ac_only = args.ac_only

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
# region Prep - RNA Lib Object

rna_library_1_str = FileUtils.load_file_as_str(replicate_1_file_path)
rna_library_1: RnaLibrary = json.loads(rna_library_1_str, object_hook=as_python_object)


rna_library_2_str = FileUtils.load_file_as_str(replicate_2_file_path)
rna_library_2: RnaLibrary = json.loads(rna_library_2_str, object_hook=as_python_object)

# endregion


# ----------------------------------
# region Pipeline - Data Prep

overall_reactivity_list_1 = []
overall_reactivity_list_2 = []

# Stats
total_nt = 0
# Total count for "load QC" nt.
total_nt_low_qc_1 = 0
total_nt_low_qc_2 = 0

for rna_item in rna_library_1.rna_items:
    #
    rna_item.calculate_reactivity_v1(nt_a_c_only=ac_only)
    overall_reactivity_list_1 += rna_item.flatten_reactivity_list()

    # Stats
    total_nt += rna_item.total_nt_with_condition(nt_a_c_only=ac_only)
    total_nt_low_qc_1 += len(rna_item.shape_profile_list_low_quality())

for rna_item in rna_library_2.rna_items:
    #
    rna_item.calculate_reactivity_v1(nt_a_c_only=ac_only)
    overall_reactivity_list_2 += rna_item.flatten_reactivity_list()

    # Stats
    total_nt_low_qc_2 += len(rna_item.shape_profile_list_low_quality())

# Discard "None" data point
overall_reactivity_list_filtered_1 = []
overall_reactivity_list_filtered_2 = []
for index, value in enumerate(overall_reactivity_list_1):
    if value != -999.0 and overall_reactivity_list_2[index] != -999.0:
        #
        overall_reactivity_list_filtered_1.append(value)
        overall_reactivity_list_filtered_2.append(overall_reactivity_list_2[index])

#
print('total: {} | replicate_1_low_qc: {} | replicate_2_low_qc: {}'
      .format(total_nt, total_nt_low_qc_1, total_nt_low_qc_2))

# Generate a temp file to capture all "pairwise data"
cwd = os.getcwd()
pairwise_data_file = os.path.join(cwd, 'pairwise_data.tsv')
with open(pairwise_data_file, 'w') as outfile:
    for index, value in enumerate(overall_reactivity_list_filtered_1):
        line = '\t'.join([str(value), str(overall_reactivity_list_filtered_2[index])])
        line += '\n'
        outfile.write(line)

# endregion

# ----------------------------------
# region Pipeline - Plot (by Python)

# import math
# import numpy
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
#
# from scipy.stats import linregress

# slope, intercept, r_value, p_value, std_err \
#     = linregress(overall_reactivity_list_filtered_1, overall_reactivity_list_filtered_2)
#
# data = {
#         'replicate_1': overall_reactivity_list_filtered_1,
#         'replicate_2': overall_reactivity_list_filtered_2,
# }
# data_frame = pd.DataFrame.from_dict(data)
# lm = sns.lmplot(x="replicate_1", y="replicate_2", data=data_frame)
#
# fig = lm.fig
# annotation = """
# R: {:.2f} | P: {:.2f} | SE {:.2f}
# """.format(r_value, p_value, std_err)
# fig.suptitle(annotation, fontsize=8)
# plt.xlabel('Neil1-TGR Replicate #1')
# plt.ylabel('Neil1-TGR Replicate #2')
# #
# plt.subplots_adjust(bottom=0.1)
#
# #
# cwd = os.getcwd()
# plt.savefig(os.path.join(cwd, '{}.png'.format('Neil1-TGR-pairwise-byShape2')))


#
# slope, intercept, r_value, p_value, std_err \
#     = linregress(overall_reactivity_list_by_own_1, overall_reactivity_list_by_own_2)
#
# data = {
#         'replicate_1': overall_reactivity_list_by_own_1,
#         'replicate_2': overall_reactivity_list_by_own_2,
# }
# data_frame = pd.DataFrame.from_dict(data)
# lm = sns.lmplot(x="replicate_1", y="replicate_2", data=data_frame)
#
# fig = lm.fig
# annotation = """
# Reactivity by Own Method
# R: {:.2f} | P: {:.2f} | SE {:.2f}
# """.format(r_value, p_value, std_err)
# fig.suptitle(annotation, fontsize=8)
# plt.xlabel('Neil1-TGR Replicate #1')
# plt.ylabel('Neil1-TGR Replicate #2')
# #
# plt.subplots_adjust(bottom=0.1)
#
# #
# cwd = os.getcwd()
# plt.savefig(os.path.join(cwd, '{}.png'.format('Neil1-TGR-pairwise-byOwn')))

# endregion


# ----------------------------------
# region Pipeline - Plot by R

pPairwise = Proc(desc='Pairwise comparison')
pPairwise.input = {'infile:file': pairwise_data_file}
pPairwise.output = 'outfile:file:{}-pairwise-byOwn-byR.png'.format(dataset_name)

#
pPairwise.args.name = dataset_name

#
# Use full path "/path/to/Rscript" if it's not in $PATH
# You can also use a shebang in script
# in this case: #!/usr/bin/env Rscript
pPairwise.lang = 'Rscript'
pPairwise.script = """
#!/usr/bin/env Rscript

library(ggplot2)

reactvitiy_replicates <- read.table("{{in.infile}}", sep="\t", header=FALSE)
names(reactvitiy_replicates) <- c("r1", "r2")

r1_r2 = lm(reactvitiy_replicates$r1 ~ reactvitiy_replicates$r2)
r1_r2
r1_r2_label = paste("R2=", format(summary(r1_r2)$adj.r.squared, digits=3))
r1_r2_label

png("{{out.outfile}}",
    width = 8, height = 8, units = "in",
    res = 300)
op <- par(mar = c(2,2,2,6), cex = 0.64)
plot1 <- ggplot(data = reactvitiy_replicates, aes(x = reactvitiy_replicates$r1, y = reactvitiy_replicates$r2), axis.ticks = 10) +
        geom_abline(slope = 1, intercept = 0) +
        geom_point(data = reactvitiy_replicates, aes(x = reactvitiy_replicates$r1, y = reactvitiy_replicates$r2), shape=21, fill='grey17', size=1.5) +
        # geom_point(data = sig1_2, aes(x = sig1_2$r1, y = sig1_2$r2), fill = 'grey67', shape=21, size=1.5) +
        theme_bw() +
        annotate('text', x = 0.05, y = 0.90, label = r1_r2_label, size = 3) +
        xlab('{{args.name}} #1') +
        ylab('{{args.name}} #2') +
        ggtitle("Pairwise Comparison") +
        theme(axis.text = element_text(size = '8'),
              panel.border = element_rect(colour='black'),
              axis.title.y = element_text(size = 8),
              axis.title.x = element_text(size = 8),
              panel.background = element_blank())

plot1
par(op) # At end of plotting, reset to previous settings
dev.off()

"""

PyPPL().start(pPairwise).run()

# endregion
