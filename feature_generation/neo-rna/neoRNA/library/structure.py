# -*- coding: utf-8 -*-

"""Final object to contain RNA Secondary Structure Info.

The info includes the following parts:

- The RNA Item, in `RNA library` type
- The Reference Structure, in `dot bracket` notation
- Reactivity data
- The "Predicted" Structure (with Reactivity data), in `dot bracket` notation

NOTE:
    The current structure prediction is done by "RNAStructure" (http://rna.urmc.rochester.edu/RNAstructure.html).

"""

from Bio.Seq import Seq
import numpy as numpy

from neoRNA.util.dot_bracket_notation import DotBracketNotation


class Structure(object):

    __STRUCTURE_LINE_NUMBER = 3

    def __init__(self, rna_item,
                 reference_structure, reference_free_energy,
                 reactivity, predicted_structure, free_energy,
                 non_mod_rates, mod_rates,
                 sequence_slice=None,
                 editing_level=None):
        """

        :param rna_item:
        :param reference_structure:
        :param reactivity:
        :type reactivity: class: `Reactivity`
        :param predicted_structure:
        :param non_mod_rates:
        :type non_mod_rates: class: `MutRateIO`
        :param mod_rates:
        :type mod_rates: class: `MutRateIO`
        :param sequence_slice: If need to do a slice on sequence
        :type sequence_slice: class: `slice`
        :param editing_level: The editing level score
        """

        self.__rna_item = rna_item
        self.__reference_structure = reference_structure
        self.__reactivity = reactivity
        self.__predicted_structure = predicted_structure
        self.__non_mod_rates = non_mod_rates
        self.__mod_rates = mod_rates

        # Free Energy
        self.__free_energy = free_energy
        self.__reference_free_energy = reference_free_energy

        if sequence_slice is not None:
            self.__rna_item.apply_sequence_slice(sequence_slice)
            self.__reference_structure = reference_structure
            self.__reactivity = reactivity[sequence_slice]
            self.__predicted_structure = predicted_structure
            self.__non_mod_rates.apply_sequence_slice(sequence_slice)
            self.__mod_rates.apply_sequence_slice(sequence_slice)

        # Process "dot-bracket" notation of reference structure
        self.__predicted_structure_dot_bracket_notation = DotBracketNotation(self.__predicted_structure, self.sequence_string)

        # Secondary Structure Info
        self.__secondary_structure = None

        self.__editing_level = editing_level

    # region Methods

    @property
    def rna_item(self):
        return self.__rna_item

    @property
    def sequence_string(self):
        return str(self.__rna_item.sequence)

    @property
    def rna_sequence_string(self):
        return str(self.__rna_item.get_rna_sequence())

    @property
    def dna_sequence_string(self):
        return str(self.__rna_item.get_dna_sequence())

    @property
    def rna_id(self):
        return self.__rna_item.id

    @property
    def barcode_string(self):
        return str(self.__rna_item.barcode)

    @property
    def reference_structure(self):
        return self.__reference_structure

    @property
    def predicted_structure(self):
        return self.__predicted_structure

    @property
    def reactivity(self):
        return self.__reactivity

    @property
    def mod_rates(self):
        return self.__mod_rates

    @property
    def non_mod_rates(self):
        return self.__non_mod_rates

    @property
    def free_energy(self):
        return self.__free_energy

    @free_energy.setter
    def free_energy(self, free_energy):
        self.__free_energy = free_energy

    @property
    def editing_level(self):
        return self.__editing_level

    @editing_level.setter
    def editing_level(self, editing_level):
        self.__editing_level = editing_level

    @property
    def secondary_structure(self):
        return self.__secondary_structure

    @secondary_structure.setter
    def secondary_structure(self, secondary_structure):
        self.__secondary_structure = secondary_structure

    @property
    def stem_length(self):
        return self.__predicted_structure_dot_bracket_notation.stem_length

    @property
    def hairpin_length(self):
        return self.__predicted_structure_dot_bracket_notation.hairpin_length

    @property
    def predicted_structure_stats(self):
        return self.__predicted_structure_dot_bracket_notation

    def get_dict(self):
        """Return a `dictionary` object including the following info:

        - RNA ID
        - RNA Barcode
        - RNA sequence
        - Mutation Syntax
        - Editing Level Score
        - Free Energy (Reference Structure)
        - RNA Reference Structure
        - RNA Reference Structure Stats Info, like Stem Length, Hairpin Length
        - Free Energy (Structure)
        - RNA Predicted Structure
        - RNA Reactivity
        - Mod Reads Total
        - Mod Read Depths
        - Non-mod Reads Total
        - Non-mod Read Depths

        :return: A `dictionary` object
        """

        _dict = {
            'rna_id': self.rna_id,
            'barcode_string': self.barcode_string,
            'sequence_string': self.rna_sequence_string,
            'mutation_syntax': self.__rna_item.mutation_syntax,
            'editing_level': self.__editing_level,
            'mod_reads': self.__rna_item.mod_reads,
            'non_mod_reads': self.__rna_item.non_mod_reads,
            'reference_free_energy': self.__reference_free_energy,
            'reference_structure': self.reference_structure,
            'free_energy': self.__free_energy,
            'predicted_structure': self.predicted_structure,
            'predicted_structure_stem_length': self.__predicted_structure_dot_bracket_notation.stem_length,
            'predicted_structure_hairpin_length': self.__predicted_structure_dot_bracket_notation.hairpin_length,
            'reactivity': self.reactivity.tolist(),
            'mod_depths': self.mod_rates.depths.tolist(),
            'non_mod_depths': self.non_mod_rates.depths.tolist()
        }

        return _dict

    # endregion


    # region Class Methods

    @classmethod
    def parse_dot_file(cls, dot_structure_file):
        """Parse structure file (dot format) to get the structure string as well as free energy.

        :param dot_structure_file: The filename to the "structure file", in `dot bracket` notation
        :rtype dot_structure_file: str
        :return: list of "structure string" and "free energy"
        """

        import re

        with open(dot_structure_file, 'rU') as infile:
            #
            free_energy = None

            # Get the "free energy" line
            line = infile.readline()
            matches = re.search(r'ENERGY =\s+(\S+)\s+', line)
            if matches:
                free_energy = matches.group(1)

            # Skip the second line
            infile.readline()

            # The third line
            structure = infile.readline().strip()

        return [structure, free_energy]

    # endregion


    # region Private Methods


    # endregion