#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["justin kuczynski", "Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""Contains tests for producing rarefied OTU tables."""

from qiime.rarefaction import RarefactionMaker
from cogent.util.unit_test import TestCase, main
import numpy

class FunctionTests(TestCase):
    def setUp(self):
        self.otu_table_transpose = numpy.array([
                                    [2,0,0,1],
                                    [1,1,1,1],
                                    [0,0,0,0]])
        self.sample_names = list('XYZ')
        self.taxon_names = list('abcd')
        self.otu_tuple = (self.sample_names, self.taxon_names, 
        self.otu_table_transpose.T,
        None)
    
    def test_rarefy_to_list(self):
        """should rarefy correctly, same names, and rm empty samples
        
        """
        maker = RarefactionMaker(self.otu_tuple, 0, 1, 1, 1)
        res = maker.rarefy_to_list(include_full=True)
        self.assertFloatEqual(res[-1][2], self.sample_names)
        self.assertFloatEqual(res[-1][3], self.taxon_names)
        self.assertFloatEqual(res[-1][4], self.otu_table_transpose.T)
        
        # each sample should have 1 seq, sample z should be removed
        self.assertFloatEqual((res[1][4]).sum(0),[1.0,1.0] )
        
#run tests if called from command line
if __name__ == '__main__':
    main()
