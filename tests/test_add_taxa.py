#!/usr/bin/env python

"""Tests of code for adding taxa to OTU table"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main
from qiime.add_taxa import fix_taxonomy_delimiters

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_fix_taxonomy_delimiters(self):
        """fix_taxonomy_delimiters should remove commas from RDP"""
        rdp = { '2483210':['Root,Bacteria','0.940'],
                '2381498':['Root,Bacteria,Firmicutes,"Clostridia",Clostridiales,"Lachnospiraceae"','1.000']}
        self.assertEqual(fix_taxonomy_delimiters(rdp),
            {'2483210':'Root;Bacteria',
             '2381498':'Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae',
             })

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
