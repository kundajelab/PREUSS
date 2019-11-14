# -*- coding: utf-8 -*-

"""
Dot Bracket Notation
================

It help defines the "Dot Bracket Notation" info for a RNA sequence
"""

import json
from neoRNA.sequence.sequence import Sequence
from neoRNA.util.file_utils import FileUtils


class DotBracketNotation(object):
    """
    Dot Bracket Notation
    """

    # ----------------------------------
    # region Init

    def __init__(self, comment, sequence_str, notation_str):
        """
        Init

        Parameters
        ----------
        comment: str
            The comment string. It may also include some "basic info" which needs further parsing.
        sequence_str: str
            The sequence string.
        notation_str: str
            The "dot-bracket" string.
        """

        #
        self.__comment = comment
        self.__sequence = Sequence(sequence_str)
        self.__notation_str = notation_str

    # endregion

    # ----------------------------------
    # region Property

    @property
    def comment(self):
        return self.__comment

    @property
    def sequence_str(self):
        return self.__sequence.sequence_str

    @property
    def notation_str(self):
        return self.__notation_str

    # endregion

    # ----------------------------------
    # region Converter

    def to_ct(self):
        r"""
        Convert the "notation string" to "ct" data

        ref:
        - ct - http://projects.binf.ku.dk/pgardner/bralibase/RNAformats.html

        Returns
        -------
        ct_data: List[str]
            CT string list
        """

        ct_str_list = []

        # Helper stacks
        stack1 = []
        stack2 = []
        stack3 = []
        stack4 = []
        stack5 = []
        stack6 = []
        stack7 = []
        stack8 = []
        stack9 = []
        stack10 = []
        stack11 = []

        # Pair into
        pairs = {}

        #
        sequence = self.__sequence.sequence_str

        for i, c in enumerate(self.__notation_str):
            if c == '(':
                stack1.append(i + 1)
            elif c == '[':
                stack2.append(i + 1)
            elif c == '{':
                stack3.append(i + 1)
            elif c == '<':
                stack4.append(i + 1)
            elif c == 'A':
                stack5.append(i + 1)
            elif c == 'B':
                stack6.append(i + 1)
            elif c == 'C':
                stack7.append(i + 1)
            elif c == 'D':
                stack8.append(i + 1)
            elif c == 'E':
                stack9.append(i + 1)
            elif c == 'F':
                stack10.append(i + 1)
            elif c == 'G':
                stack11.append(i + 1)
            elif c == ')':
                pairs[i + 1] = stack1.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == ']':
                pairs[i + 1] = stack2.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == '}':
                pairs[i + 1] = stack3.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == '>':
                pairs[i + 1] = stack4.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'a':
                pairs[i + 1] = stack5.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'b':
                pairs[i + 1] = stack6.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'c':
                pairs[i + 1] = stack7.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'd':
                pairs[i + 1] = stack8.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'e':
                pairs[i + 1] = stack9.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'f':
                pairs[i + 1] = stack10.pop()
                pairs[pairs[i + 1]] = i + 1
            elif c == 'g':
                pairs[i + 1] = stack11.pop()
                pairs[pairs[i + 1]] = i + 1

        for i in range(1, len(self.__notation_str) + 1):
            ct_str_list.append(
                "%d%s%s%s%d%s%d%s%d%s%d" % (i, ' ', sequence[i - 1], ' ', i - 1, ' ', i + 1, ' ', pairs.get(i, 0), ' ', i))

        return ct_str_list

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
            "dot_notation": self.__notation_str,
        }

        return json_data

        # return json.dumps(json_data,
        #                   sort_keys=False,
        #                   indent=4, ensure_ascii=True)
        # return json.dumps(self,
        #                   default=lambda o: o.__dict__,
        #                   sort_keys=True, indent=4)

    # endregion

    # ----------------------------------
    # region Output

    def to_file(self, file_path: str):
        r"""
        Output the "Dot Bracket Notation" record as a single file.

        The file is in "dbn" format.

        Parameters
        ----------
        file_path: str
            The "target" file path of the output.

        Returns
        -------

        """

        content = '>{}\n{}\n{}'.format(self.__comment, self.__sequence.sequence_str, self.__notation_str)
        #
        FileUtils.save_file(file_path, content)

    # endregion
