# -*- coding: utf-8 -*-

"""
VARNA Runners
================

The running wrapper for "VARNA.jar"
"""

import os
import subprocess

from typing import List


class VarnaRunner(object):
    r"""
    VARNA Runner

    This runner is based on the "jar" paackage of VERNA.

    Ref: http://varna.lri.fr/index.php?lang=en&page=command&css=varna
    """

    BASE_CMD_STR_TEMPLATE = 'java -cp {} fr.orsay.lri.varna.applications.VARNAcmd {}'

    # Default options used for "VARNA.jar"
    DEFAULT_OPTIONS = {
        'algorithm': "radiate",
        'flat': 'true',

        # Base Pair related
        'bpStyle': 'Line',
        'baseInner': '#FFFFFF',         # Remove fill color inside each "nt circle"
        'baseOutline': '#FFFFFF',       # Remove "nt circle"
        'bp': '#000000',                # nt color
        'spaceBetweenBases': '0.8',     # Space between each BP, in "ratio"
        'baseNum': '#FFFFFF',           # Hide "number" of nt counting

        #
        'resolution': '4.0'             # Size ratio
    }

    # Default highlight region options
    DEFAULT_HIGHLIGHT_REGION_OPTIONS = {
        'fill': "#C55337",
        'outline': '#FFFFFF',
    }

    # Default annotation options
    DEFAULT_ANNOTATION_OPTIONS = {
        'type': "B",
        'anchor': "66",
        'outline': '#000000',
        'size': '6',
    }

    # ----------------------------------
    # region Init

    def __init__(self, varna_location: str):
        r"""
        Init

        Parameters
        ----------
        varna_location: str
            The "path" to the VARNA jar package.

        """

        self.varna_location = varna_location
        self.varna_options = self.DEFAULT_OPTIONS

    # endregion

    # ----------------------------------
    # region Methods - RNA Structure Image

    def gen_struct_image(self, rna_sequence: str, rna_structure: str, out_file: str,
                         options: dict = None):
        r"""
        Generate RNA Structre image.

        Parameters
        ----------
        rna_sequence: str
            The RNA sequence
        rna_structure: str
            A well-parenthesized expression with dots whose size matches that of the input sequence.
        out_file: str
            An output file whose format is guessed from the extension.
        options: dict
            The CMD options for "VARNA.jar". Optional.

        Returns
        -------
        """

        # Update options if needed
        if options:
            self.varna_options.update(options)

        # Prep the "final" parameters
        params = {
            'sequenceDBN': rna_sequence,
            'structureDBN': rna_structure,
            'o': out_file,
        }
        params.update(self.varna_options)

        # Call
        cmd = self.__gen_cmd_str(params)
        try:
            print(cmd)
            subprocess.call(cmd, shell=True)
        except OSError as error:
            print('VARNA.jar command running error - ', error.args)

    # endregion

    # ----------------------------------
    # region Methods - Utils

    def gen_highlight_region_str(self, regions: list):
        r"""
        Generate "parameter string" for VARNA Highlight Region.

        Each element in "region" list is a "definition dict" with the following properties:
        - nt_range, like "47-47", "1-10"
        - fill: The filling color.
        - outline: The color of "outline".

        Parameters
        ----------
        regions: list
            A list of "region" definitions.

        Returns
        -------

        """

        region_str_list = []
        for region_options in regions:
            #
            options = self.DEFAULT_HIGHLIGHT_REGION_OPTIONS
            options.update(region_options)

            #
            if 'nt_range' not in options:
                continue
            nt_range_str = options['nt_range']
            options.pop('nt_range', None)  # Remove key 'nt_range'

            #
            region_str = '{}:{}'.format(nt_range_str, self.__gen_param_str(options,
                                                                           prefix='', equal_sign='=', separator=',', use_quote=False))
            region_str_list.append(region_str)

        #
        return ';'.join(region_str_list)

    def gen_annotation_str(self, annotations: list):
        r"""
        Generate "parameter string" for VARNA Annotations.

        Each element in "region" list is a "definition dict" with the following properties:
        - annotation_str, like "E-47", "Mutation"
        - type: Annotation Type, type=[P,B,H,L]
        - anchor: The anchor "nt".
        - color: The color of "outline".
        - size: The size of annotation string.

        Parameters
        ----------
        annotations: list
            A list of "annotations" definitions.

        Returns
        -------

        """

        annotation_str_list = []
        for annotation_options in annotations:
            #
            options = self.DEFAULT_ANNOTATION_OPTIONS
            options.update(annotation_options)

            #
            if 'annotation_str' not in options:
                continue
            annotation_str = options['annotation_str']
            options.pop('annotation_str', None)  # Remove key 'annotation_str'

            #
            annotation_str = '{}:{}'.format(annotation_str, self.__gen_param_str(options,
                                                                                 prefix='', equal_sign='=', separator=',', use_quote=False))
            annotation_str_list.append(annotation_str)

        #
        return ';'.join(annotation_str_list)

    # endregion

    # ----------------------------------
    # region Internal Methods

    def __gen_cmd_str(self, params):
        r"""
        Generate "CMD string" for VARNA package.

        Parameters
        ----------
        params: dict
            The parameter object.

        Returns
        -------
        cmd_str: str
            The CMD string.
        """

        #
        params_str = self.__gen_param_str(params)

        #
        return self.BASE_CMD_STR_TEMPLATE.format(self.varna_location, params_str)

    def __gen_param_str(self, params: dict,
                        prefix: str = '-', equal_sign: str = ' ', separator: str = ' ',
                        use_quote: bool = True) -> str:
        r"""
        Generate parameter string.

        Each element in the "parameter" dict will be string.

        Parameters
        ----------
        params: dict
            The parameter object.
        prefix: str
            The "prefix" used when generate "each" parameter pair.
        equal_sign: str
            The "equal sign" used between the "key" and "value".
        separator: str
            The "separator string" between each parameter item.
        use_quote: bool
            If use "quote" on the "value".

        Returns
        -------
        param_str: str
            The parameter string.
        """

        if not params:
            return ''

        param_str_items = []
        for param, value in params.items():
            #
            str_template = '{}{}{}"{}"' if use_quote else '{}{}{}{}'

            item_str = str_template.format(prefix, param, equal_sign, value)
            param_str_items.append(item_str)

        #
        return separator.join(param_str_items)

    # endregion
