#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski", "Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main

import qiime.golay as golay
""" tests the golay DNA barcode decode/encode functionality"""

class GolayTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Set up shared variables"""
        pass

    def test_golay1(self):
        """switching the last base should recover the original barcode"""
        sent = golay.encode([0,0,0,0,0,0,0,0,0,1,0,0])
        rec = sent[:-1] + 'C' # possible error here
        decoded, errors = golay.decode(rec)
        self.assertEqual(decoded, sent)
        self.assertLessThan(errors, 3)
        rec = sent[:-1] + 'T' # possible error here
        decoded, errors = golay.decode(rec)
        self.assertEqual(decoded, sent)
        self.assertLessThan(errors, 3)

    def test_golay_matches_old_code(self):
        """ should behave as micah's code did, i.e., same golay encoding
        this requires 
        DEFAULT_NT_TO_BITS = { "A":"11",  "C":"00", "T":"10", "G":"01"}
        """
        original = 'GCATCGTCAACA'
        error = 'GCATCGTCCACA'
        self.assertEqual(golay.decode(error)[0], original)

if __name__ == '__main__':
    main()
