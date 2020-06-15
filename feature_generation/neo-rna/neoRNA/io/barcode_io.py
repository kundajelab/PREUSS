# -*- coding: utf-8 -*-

from neoRNA.sequence.barcode import Barcode


class BarcodeIO(object):
    """
    Parse the "Barcode file" into a list of "barcodes".

    The format of "barcode file":

    ```
    # Comment-line, to be ignored when parsing
    # Each line uses `\t` to separate
    ID_1    ATCGCGC
    ID_@    AAATCCGC
    ```
    """

    # region Init
    # -----------------------------------------
    # Init
    # -----------------------------------------
    def __init__(self, barcode_file):
        """
        Parse the barcode file into a list of "barcodes".

        :param barcode_file: The barcode file.
        """
        self.__barcode_file = barcode_file

        # Parse the file
        self.__barcodes = self.__parse(self.__barcode_file)

    # endregion

    # region Methods
    # -----------------------------------------
    # Methods
    # -----------------------------------------
    def barcodes(self):
        """
        Return the parsed list of barcodes, in `Barcode` type.

        :rtype: Barcode
        :return: The list of barcodes, in `Barcode` type.
        """
        return self.__barcodes

    # endregion

    # region Private Methods
    # -----------------------------------------
    # Private Methods
    # -----------------------------------------
    def __parse(self, barcode_file):
        """
        Parse the barcode file into a list of "barcodes", in `Barcode` type.

        :rtype: Barcode
        :return: The list of barcodes.
        """
        barcodes = []

        with open(barcode_file) as infile:
            for line in infile:
                # Ignore line starting with "#"
                if line.startswith('#'):
                    continue

                # Parse the barcode and reference_id
                reference_id, barcode_string = line.strip().split('\t')
                barcode = Barcode(barcode_string, reference_id)
                barcodes.append(barcode)

        return barcodes

    # endregion
