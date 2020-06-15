# -*- coding: utf-8 -*-

"""
DNA/RNA Base Pair
================

"""


class BasePair(object):
    """
    An Object definition of DNA/RNA Base Pair.

    For each of the Base Pair, it include the following info:
    - The "left" nucleotide and its position
    - The "right" nucleotide and its position

    """

    # ----------------------------------
    # region Init

    def __init__(self, nt_pair, position_pair):
        """

        :param nt_pair: `tuple`, Nucleotide pair,
        :param position_pair: `tuple`, Position pair, need to be paired with "Nucleotide pair"
        """
        #
        self.__nt_pair = nt_pair
        self.__position_pair = position_pair

        # Validate the pair info
        self._validate()

    # endregion

    # ----------------------------------
    # region Properties

    @property
    def left_nt(self):
        return self.__nt_pair[0]

    @property
    def left_position(self):
        return self.__position_pair[0]

    @property
    def right_nt(self):
        return self.__nt_pair[1]

    @property
    def right_position(self):
        return self.__position_pair[1]

    # endregion

    # ----------------------------------
    # region Methods

    def get_nt_pair_string(self, format_string='{}:{}'):
        """
        Get "nt pair" string, based on the given string "format".

        :param format_string:
        :return:
        """

        return format_string.format(self.left_nt, self.right_nt)

    # endregion

    # ----------------------------------
    # region Internal Methods

    def _validate(self):
        """
        Validate the "pair info"

        - The pair needs to have "TWO" elements.

        :return:
        """

        if len(self.__nt_pair) != 2 or len(self.__position_pair) != 2:
            raise ValueError('The nucleotide pair is not valid. It must have 2 elements. ')

    # endregion

    # ----------------------------------
    # region Methods - Magic ones

    def __str__(self):
        """
        Return a string format of the base pair

        :return:
        """

        # Example: (47,30) (G, C)
        return '{} {}'.format(str(self.__position_pair), str(self.__nt_pair))

    # endregion
