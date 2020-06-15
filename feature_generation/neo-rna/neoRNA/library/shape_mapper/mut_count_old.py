# -*- coding: utf-8 -*-

"""Mutation Count record (Old Style)

Each record object represents **one** line from the file.

The line string includes **three** parts, separated by `\t`.

- start position
- end position
- mutation status - a string with each character representing one nt.

----

"""


class MutCountOld(object):
    """RNA library Item Object.

    Each item should include the following elements:

    - RNA ID: The unique ID (usually a number) throughout the library.
    - RNA Barcode: The unique "barcode" sequence that represents this RNA
    - RNA sequence: The actual target sequence of this RNA
    - Notes: A string to describe this RNA, supplementary info.
    """

    __NORMAL_MARKER = '|'
    __MUT_MARKERS = ['A','T','G','C','-']
    __IGNORE_MARKERS = ['s', '~']
    __NON_OVERLAP_MARKERS = ['s|', '~|']
    __IGNORE_NON_OVERLAP = True

    def __init__(self, start_position, end_position, mutation_string):
        self.__start = start_position
        self.__end = end_position
        self.__mutation_string = mutation_string

        # Check if it is a valid record
        # Be sure to check it "before" doing self adjustment
        self.__is_valid = self.__is_valid_record()

        # Adjust
        if self.__is_valid:
            self.__adjust()

    # region Methods

    @property
    def start_position(self):
        """The start position of the sequence.

        :note: It starts from "1", not "0" as the list.

        :return:
        """
        return self.__start

    @property
    def end_position(self):
        return self.__end

    @property
    def mutation_string(self):
        return self.__mutation_string

    @property
    def is_valid(self):
        return self.__is_valid

    def simple_mutation_string(self):
        """Generate "simple" version of mutation string.

         In the "simple" version:

         - All "mutation types" are marked as "1"
         - All others are marked as "0"

        :return: The mutation string (simple version).
        """

        import re

        replace_matchers = {}

        for marker in self.__MUT_MARKERS:
            replace_matchers[marker] = '1'

        for marker in self.__IGNORE_MARKERS:
            replace_matchers[marker] = '0'

        replace_matchers[self.__NORMAL_MARKER] = '0'

        return re.sub('|'.join(re.escape(key) for key in replace_matchers.keys()),
                      lambda k: replace_matchers[k.group(0)], self.__mutation_string)

    # endregion

    # region Private Methods

    def __is_valid_record(self):
        """Decide if this mutation count record is valid.

        If it is valid,

        - the `length` of `mutation string` should equal `start_position - end_position + 1`
        - it should not contain any "non-overlap" case (if "ignore non-overlap" flag is set to be ON)

        :return: True | False
        """

        __valid = True

        if not len(self.__mutation_string) == (int(self.__end) - int(self.__start) + 1):
            __valid = False

        if self.__IGNORE_NON_OVERLAP:
            #
            __mutation_string_trimmed = self.__mutation_string.strip(''.join(self.__IGNORE_MARKERS))
            for marker in self.__NON_OVERLAP_MARKERS:
                if __mutation_string_trimmed.find(marker) >= 0:
                    __valid = False
                    break

        return __valid

    def __adjust(self):
        """Adjust "mutation string" to trim "ignore markers" ('s', '~')

        :return:
        """
        # Trim from the beginning
        for i in range(len(self.__mutation_string)):
            if self.__mutation_string[i] in self.__IGNORE_MARKERS:
                self.__start += 1
            else:
                break

        # Trim from the end
        for i in reversed(range(len(self.__mutation_string))):
            if self.__mutation_string[i] in self.__IGNORE_MARKERS:
                self.__end -= 1
            else:
                break

        # Also update the new mutation string based on the new "start" and "end" position
        self.__mutation_string = self.__mutation_string[(self.__start - 1):self.__end]

        return

    # endregion