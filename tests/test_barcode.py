#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Justin Kuczynski", "Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from cogent.util.unit_test import TestCase, main

import qiime.barcode as barcode

class BarcodeTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """Set up shared variables"""

    def test_correct_barcode(self):
        """ correct_barcode should correctly assign barcode to 2nd possibility, 
        2 errors
        """
        original = 'ATTTTTTTTTCG'
        recieved = 'ATTTTTTTTTTT'
        possibilities = ['TGTATTCGTGTA','ATTTTTTTTTCG','TGTAGGCGTGTA',
            'TGTAGAAGTGTA','TGTAGGCGTATA','TGTAAAAAAAAA']
        decoded, num_errors = barcode.correct_barcode(recieved, possibilities)
        self.assertEqual(decoded, original)
        self.assertEqual(num_errors, 2)
    
    def test_correct_barcode_bitwise(self):
        """ correct_barcode should correctly assign barcode to 2nd possibility, 
        2 base errors,
        3 bit errors with this nt_to_bits
        """
        nt_to_bits = { "A":"11",  "C":"00", "T":"10", "G":"01"}

        original = 'ATTTTTTTTTCG'
        recieved = 'ATTTTTTTTTTT'
        possibilities = ['TGTATTCGTGTA','ATTTTTTTTTCG','TGTAGGCGTGTA',
            'TGTAGAAGTGTA','TGTAGGCGTATA','TGTAAAAAAAAA']
        decoded, num_errors = barcode.correct_barcode_bitwise(\
            recieved, possibilities, nt_to_bits)
        self.assertEqual(decoded, original)
        self.assertEqual(num_errors, 3)
    
    def test_correct_barcode_bitwise_tie(self):
        """ correct_barcode should not assign barcode to to a tie situation
        """
        nt_to_bits = { "A":"11",  "C":"00", "T":"10", "G":"01"}

        #~ original = 'ATTTTTTTTTCG' #doesn't matter, last entry is just
        # as close to recieved as original is
        recieved = 'ATTTTTTTTTTT'
        possibilities = ['TGTATTCGTGTA','ATTTTTTTTTCG','TGTAGGCGTGTA',
            'TGTAGAAGTGTA','TGTAGGCGTATA','TGTAAAAAAAAA', 'ATTTTTTTTAAA']
        decoded, num_errors = barcode.correct_barcode_bitwise(\
            recieved, possibilities, nt_to_bits)
        self.assertEqual(decoded, None)
        self.assertEqual(num_errors, 3)

if __name__ == '__main__':
    main()

