# -*- coding: utf-8 -*-

"""
Command Runner - ShapeMapper
================

CMD runner for "ShapeMapper".
"""

import json
import os
import subprocess

from typing import Dict, Any

from neoRNA.library.library_config import RnaLibConfig
from neoRNA.library.library_item import LibraryItem

from neoRNA.util.json_serializable import as_python_object
from neoRNA.util.file_utils import FileUtils


class ShapeMapperRunner(object):
    """
    CMD runner for "ShapeMapper".

    Currently "ShapeMapper" has two versions - 1.x and 2.x.
    """

    # ----------------------------------
    # region ShapeMapper v1.x

    @classmethod
    def shape_mapper_v1(cls, working_folder, config_file):
        r"""
        Run "ShapeMapper" pipeline.

        Parameters
        ----------
        working_folder: str
            Path to "working folder".
        config_file: str
            Path to "config file" of ShapeMapper.

        Returns
        -------

        """
        #
        command_shape_mapper = 'ShapeMapper.py ' + config_file

        #
        try:
            os.chdir(working_folder)
            subprocess.call([command_shape_mapper])
        except OSError:
            print ('ShapeMapper command does not exist!')

    # endregion

    # ----------------------------------
    # region ShapeMapper v2x

    @classmethod
    def shape_mapper_v2(cls, configs_json: str, rna_item_json: str):
        """
        Run "ShapeMapper 2.x" pipeline.

        It depends on the "RNA Lib" configs.

        ## Steps
        - Create "target" sequence file - used for the "CMD"
        - Prepare "CMD"
        - Run "CMD"

        Parameters
        ----------
        configs_json: str
            The "RNA Lib" configs.
        rna_item_json: str
            The RNA Lib Item object

        Returns
        -------

        """

        #
        configs: RnaLibConfig = json.loads(configs_json, object_hook=as_python_object)
        rna_item: LibraryItem = json.loads(rna_item_json, object_hook=as_python_object)

        rna_id = rna_item.rna_id
        barcode = rna_item.barcode
        sequence = rna_item.sequence

        # working folder
        working_folder = configs.working_folder
        if not working_folder:
            raise ValueError('"Working Folder" does not exist. ')

        # Create the folder for "shapemapper" results
        shapemapper_results_folder_path = os.path.join(working_folder, 'shapemapper_results')
        if not os.path.exists(shapemapper_results_folder_path):
            os.mkdir(shapemapper_results_folder_path)
        # Output folder for ShapeMapper results
        output_folder_path = \
            os.path.join(shapemapper_results_folder_path, '{}_{}_out'.format(rna_id, barcode.barcode))

        # Create the folder for "shapemapper" temp results
        shapemapper_temp_folder_path = os.path.join(working_folder, 'shapemapper_temp')
        if not os.path.exists(shapemapper_temp_folder_path):
            os.mkdir(shapemapper_temp_folder_path)
        # Output folder for ShapeMapper results
        temp_folder_path = \
            os.path.join(shapemapper_temp_folder_path, '{}_{}_temp'.format(rna_id, barcode.barcode))

        # Create "target" sequence file - "RNA_ID.fa"
        # Here, it needs to be "DNA" sequence
        rna_sequence_content = '>' + rna_id + '\n' + str(sequence.get_dna_sequence())
        rna_sequence_file_path = os.path.join(shapemapper_results_folder_path, rna_id + '.fa')
        FileUtils.save_file(rna_sequence_file_path, rna_sequence_content)

        # Data Source Path
        if not configs.use_barcode_as_folder_name():
            # Use "rna id" as the "folder name" of each rna item.
            data_folder_path = os.path.join(configs['data_source_folder_path'], rna_id)
        else:
            data_folder_path = os.path.join(configs['data_source_folder_path'], barcode.barcode)

        # Prepare "CMD"
        cmd = 'shapemapper --target {} '.format(rna_sequence_file_path)
        cmd += ' --out {} '.format(output_folder_path)
        cmd += ' --log {} '.format(os.path.join(output_folder_path, 'shapemapper_log.txt'))
        cmd += ' --temp {} '.format(temp_folder_path)
        if configs['shape_min_read_depth']:
            cmd += ' --min-depth {} '.format(configs['shape_min_read_depth'])

        if configs['shape_use_folder'] is True:
            # Use folder parameter
            cmd += ' --modified --folder {} '\
                .format(os.path.join(data_folder_path, configs['shape_modified_folder']))
            cmd += ' --untreated --folder {} ' \
                .format(os.path.join(data_folder_path, configs['shape_untreated_folder']))
            # Check if need to include `denatured` data
            if configs['shape_denatured_folder'] != '':
                cmd += ' --denatured --folder {} '\
                    .format(os.path.join(data_folder_path, configs['shape_denatured_folder']))
        else:
            cmd += ' --modified --R1 {} --R2 {}'\
                .format(os.path.join(data_folder_path, configs['shape_modified_r1']),
                        os.path.join(data_folder_path, configs['shape_modified_r2']))
            cmd += ' --untreated --R1 {} --R2 {}' \
                .format(os.path.join(data_folder_path, configs['shape_untreated_r1']),
                        os.path.join(data_folder_path, configs['shape_untreated_r2']))
            if configs['shape_denatured_r1'] != '':
                cmd += ' --denatured --R1 {} --R2 {}' \
                    .format(os.path.join(data_folder_path, configs['shape_denatured_r1']),
                            os.path.join(data_folder_path, configs['shape_denatured_r2']))

        # Check if there is extra flags
        if configs['shape_extra_flags'] != '':
            cmd += ' {} '.format(configs['shape_extra_flags'])

        # Run CMD
        try:
            print(cmd)
            os.chdir(data_folder_path)
            subprocess.call(cmd, shell=True)
        except OSError as error:
            print('ShapeMapper 2.x command running error - ', error.args)

    # endregion


