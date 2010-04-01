#!/usr/bin/env python
# File created on 1 Apr 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"

import qiime.biplots as bp
import numpy as np

from cogent.util.unit_test import TestCase, main


class BiplotTests(TestCase):
    
    def setUp(self):
        pass

    def tearDown(self):
        pass
        
    def test_get_taxa(self):
        pass
        
    def test_get_taxa_coords(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        sample_names = list('WXYZ')
        otu_names = list('abc')
    
        res = bp.get_taxa_coords(otu_table, [.4,.2,.1,.9])
        otu_coords= range(3)
        otu_coords[0] = .4*2/3 + .9*1/3
        otu_coords[1] = .4*1/4 + .2*1/4 + .1*1/4 + .9*1/4
        otu_coords[2] = .4*0/5 + .2*2/5 + .1*2/5 + .9*1/5
        self.assertFloatEqual(res, otu_coords)
    
    def test_get_taxa_prevalence(self):
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,0,0,0]],float)
        sample_weights = [3,1,1,2]
        res = bp.get_taxa_prevalence(otu_table)
        # print res
        # self.assertFloatEqual(res, np.array([(2/3) + 1/2, 1/3+1+1+1/2, 0])/4) 
        self.assertFloatEqual(res, np.array([(2/3) + 1/2, 1/3+1+1+1/2, 0])/4\
            * 4/(2.5+1/3))                    
        otu_table = np.array([  [2,0,0,1],
                                [1,1,1,1],
                                [0,2,2,1]],float)
        res = bp.get_taxa_prevalence(otu_table)
        # print res
        # self.assertFloatEqual(res, np.array([3,4,5])/12) # if no normalize
        self.assertFloatEqual(res, [0,.5,1])
        
    def test_remove_rare_taxa(self):
        pass
if __name__ == "__main__":
    main()