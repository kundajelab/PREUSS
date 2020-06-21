# -*- coding: utf-8 -*-

"""
Base Pair Probability
================

It help defines the "Dot Bracket Notation" info for a RNA sequence
"""

import json
from neoRNA.sequence.sequence import Sequence


class BasePairProbability(object):
    """
    Base Pair Probability
    """

    # ----------------------------------
    # region Init

    def __init__(self, comment, bbp_str, sequence_str=None):
        """
        Init

        Parameters
        ----------
        comment: str
            The comment string. It may also include some "basic info" which needs further parsing.
        bbp_str: str
            The string that includes "bpp" info.
        sequence_str: str
            The sequence string.
        """

        #
        self.__comment = comment
        self.__bbp_str = bbp_str
        self.__sequence = Sequence(sequence_str) if sequence_str else None

        #
        self.__bpp_list = self._parse_bpp_list(self.__bbp_str)

    # endregion

    # ----------------------------------
    # region Property

    @property
    def comment(self):
        return self.__comment

    @property
    def bpp_list(self):
        return self.__bpp_list

    @property
    def sequence(self):
        return self.__sequence

    @sequence.setter
    def sequence(self, sequence_str):
        self.__sequence = Sequence(sequence_str) if sequence_str else None

    # endregion

    # ----------------------------------
    # region Private Methods

    def _parse_bpp_list(self, bpp_str):
        r"""
        Parse a bpp list from the "string" list.

        Parameters
        ----------
        bpp_str: str
            The string that includes "bpp" info.

        Returns
        -------
        bpp_list: list
            A list of bpp values, in float.
        """
        #
        if not bpp_str:
            return None

        import re
        return [float(item.strip()) for item in re.split("\s+", bpp_str)]

    # endregion

    # ----------------------------------
    # region JSON

    def to_json(self):
        r"""
        Create a "JSON compatible" object when using "JSON Dump" to export this class.

        Typical Usage
        -------
        >>> def dumper(obj):
        ... try:
        ...     return obj.toJSON()
        ... except:
        ...     return obj.__dict__
        ...
        ... json_data = json.dumps(json_dict, default=dumper, sort_keys=False, indent=4, ensure_ascii=True)


        Returns
        -------
        json_object: json
            The JSON object.

        """

        json_data = {
            "comment": self.__comment,
            "sequence": self.__sequence.sequence_str,
            "bpp": self.__bbp_str,
        }

        return json_data

        # return json.dumps(json_data,
        #                   sort_keys=False,
        #                   indent=4, ensure_ascii=True)
        # return json.dumps(self,
        #                   default=lambda o: o.__dict__,
        #                   sort_keys=True, indent=4)

    # endregion
