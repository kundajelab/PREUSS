# -*- coding: utf-8 -*-

"""
Comment Parser
================

Parser for the "comment" string.
"""


class CommentParser(object):
    """
    "Comment string" parser.

    A "comment string" is commonly used in various of NGS file format (like "fasta", "dot-bracket" files, etc.) which
    helps convey some "basic info".

    This parser is to help retrieve this info.

    """

    # ----------------------------------
    # region RSample

    @classmethod
    def parse_rsample_dot_bracket(cls, comment_str):
        """
        Parse the "comment string" for a "RSample dot-bracket" result file.

        An example for the "RSample comment string":
            >001 centroid 1 655

        Elements:
        - '001' - RNA ID
        - 'centroid' - clustering result type, centroid | bpp
        - '1' - clustering #
        - '655' - clustering size

        Parameters
        ----------
        comment_str: str
            The "comment string".

        Returns
        -------
        parsing_elements: dict
            The parsing results.

        """

        import re

        # Regex format to match the string
        regex = "^(?P<id>\w+)\s+(?P<type>\w+)\s+(?P<cluster_id>\w+)\s+(?P<cluster_count>\w+)$"
        matched = re.match(regex, comment_str.strip())

        if not matched:
            return None

        #
        return matched.groupdict()

    @classmethod
    def parse_rsample_bpp(cls, comment_str):
        """
        Parse the "comment string" for a "RSample bpp" result file.

        An example for the "RSample comment string":
            >81 001 bpp 1 655

        Elements:
        - '81' - sequence length
        - '001' - RNA ID
        - 'centroid' - clustering result type, centroid | bpp
        - '1' - clustering #
        - '655' - clustering size

        Parameters
        ----------
        comment_str: str
            The "comment string".

        Returns
        -------
        parsing_elements: dict
            The parsing results.

        """

        import re

        # Regex format to match the string
        regex = "^(?P<length>\w+)\s+(?P<id>\w+)\s+(?P<type>\w+)\s+(?P<cluster_id>\w+)\s+(?P<cluster_count>\w+)$"
        matched = re.match(regex, comment_str.strip())

        if not matched:
            return None

        #
        return matched.groupdict()

    # endregion

    # ----------------------------------
    # region bpRNA

    @classmethod
    def parse_bp_rna(cls, comment_str):
        """
        Parse the "comment string" for a "bpRNA" result file.

        An example for the "RSample comment string":
            001.centroid.1

        Elements:
        - '001' - RNA ID
        - 'centroid' - clustering result type, centroid | bpp
        - '1' - clustering #

        Parameters
        ----------
        comment_str: str
            The "comment string".

        Returns
        -------
        parsing_elements: dict
            The parsing results.

        """

        import re

        # Regex format to match the string
        # regex = "^(?P<id>\w+)\.(?P<type>\w+)\.(?P<cluster_id>\w+)$"
        regex = "^(?P<id>[^_]+)_.*$"
        matched = re.match(regex, comment_str.strip())

        if not matched:
            return None

        #
        return matched.groupdict()

    # endregion
