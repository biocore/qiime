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
from qiime.summarize_taxa import parse_command_line_parameters #no functions yet

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_import(self):
        """empty test just verifies import of module"""
        pass

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
