#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main
from qiime.make_library_id_lists import (get_ids, get_first_id)
"""Tests of the make_library_id_lists.py file.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


class TopLevelTests(TestCase):

    """Top-level tests of functions in make_library_id_lists.py"""

    def test_get_ids(self):
        """get_ids should identify which ids are in which library"""
        lines = """>S74_1 E86FECS01CEVAV orig_bc=ACATGTCACGTG new_bc=ACATGTCACGTG bc_diffs=0
CTCCTC
>Unassigned_2 E86FECS01EKKMF orig_bc=AGCGCTGATGTA new_bc=None bc_diffs=1
GGTGCCTCCCTCGC
>S80_3 E86FECS01EKKMF orig_bc=AGCGCTGATGTA new_bc=None bc_diffs=1
GGTGCCTCCCTCGC
>S80_4 E86FECS01CW66X orig_bc=AGTCCATAGCTG new_bc=AGTCCATAGCTG bc_diffs=0
GTCCTGGCAG""".splitlines()
        result = get_ids(lines, 1)
        self.assertEqual(
            dict(result),
            {'S74': ['E86FECS01CEVAV'],
             'Unassigned': ['E86FECS01EKKMF'],
             'S80': ['E86FECS01EKKMF',
                     'E86FECS01CW66X']})

    def test_get_first_id(self):
        """get_first_id should identify first id in fasta file"""
        lines = """>S74_1 E86FECS01CEVAV orig_bc=ACATGTCACGTG new_bc=ACATGTCACGTG bc_diffs=0
CTCCTC
>Unassigned_2 E86FECS01EKKMF orig_bc=AGCGCTGATGTA new_bc=None bc_diffs=1
GGTGCCTCCCTCGC
>S80_3 E86FECS01EKKMF orig_bc=AGCGCTGATGTA new_bc=None bc_diffs=1
GGTGCCTCCCTCGC
>S80_4 E86FECS01CW66X orig_bc=AGTCCATAGCTG new_bc=AGTCCATAGCTG bc_diffs=0
GTCCTGGCAG""".splitlines()
        self.assertEqual(
            get_first_id(lines),
            set(['S74_1',
                 'Unassigned_2',
                 'S80_3',
                 'S80_4']))

if __name__ == '__main__':
    main()
