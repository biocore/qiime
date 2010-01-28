#!/usr/bin/env python
#file test_make_rarefaction_plots.py
#from __future__ import division
__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Meg Pirrung"] 
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Meg Pirrung"
__email__ = "meg.pirrung@colorado.edu"
__status__ = "Pre-release"

"""
Unit tests for make_rarefaction_plots.py
"""

from cogent.util.unit_test import TestCase, main
from qiime.make_rarefaction_plots import *
from qiime import parse

class makeRarefactionPlotsTests(TestCase):
    
    def setUp(self):
        self.mappingfile = ['#SampleID\tSex\tAge',
                            '123\tF\t32',
                            '234\tM\t30',
                            '345\tM\t32']
        self.p_mappingfile = parse.parse_map(self.mappingfile,return_header=True, strip_quotes=True)
        self.p_mappingfile[0][0] = [h.strip('#').strip(' ') for h in  self.p_mappingfile[0][0]]
                            
        self.rarefactionfile = ['\tsequences per sample\titeration\t123\t234\t345',
                                'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
                                'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
                                'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
                                'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
                                'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
                                'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642']
                                
        self.matrix, self.seqs_per_samp, self.sampleIDs = parse_rarefaction(self.rarefactionfile)
        
        self.ave_seqs_per_sample = {'123':[2.03172,9.4417849999999994,12.508435],
                                    '234':[0.42876999999999998,0.42876999999999998,0.42876999999999998],
                                    '345':[2.255255,9.625995,11.58785]}
        
        self.collapsed_ser_sex = {'M':[1.3420125000000001,5.0273824999999999,6.0083099999999998], 
                                    'F':[2.03172,9.4417849999999994,12.508435]}
        self.err_ser_sex = {'M':[0.91324250000000007,4.5986124999999998,5.5795399999999997],
                            'F':[0.0,0.0,0.0]}
        self.overall_averages = {'123':7.9939800000000005,
                                    '234':0.42876999999999993,
                                    '345':7.8230333333333322}
    
    def test_ave_seqs_per_sample(self):
        test = ave_seqs_per_sample(self.matrix,self.seqs_per_samp,self.sampleIDs)
        self.assertEqual(test, self.ave_seqs_per_sample)
    
    def test_is_max_category_ops_neg(self):
            test = is_max_category_ops(self.p_mappingfile, 'Sex')
            self.assertEqual(test[0], False)
            
    def test_is_max_category_ops_pos(self):
            test = is_max_category_ops(self.p_mappingfile, 'SampleID')
            self.assertEqual(test[0], True)
            
    def test_make_error_series(self):
        #make_error_series(rare_mat, sampleIDs, mapping, mapping_category)
        test = make_error_series(self.ave_seqs_per_sample,self.sampleIDs,self.p_mappingfile,'Sex')
        self.assertEqual(test[0], self.collapsed_ser_sex)
        self.assertEqual(test[1], self.err_ser_sex)
    
    def test_get_overall_averages(self):
        test = get_overall_averages(self.ave_seqs_per_sample,self.sampleIDs)
        self.assertEqual(test, self.overall_averages)

#run tests if called from command line
if __name__ == '__main__':
    main()