#!/usr/bin/env python
# File created on 27 Sep 2011
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from qiime.distance_matrix_from_mapping import distance_matrix
from numpy import array
from cogent.util.unit_test import TestCase, main
import StringIO
        

class FunctionTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        self.fasting_map = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Float_Col	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	.1	Control_mouse__I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	.2	Control_mouse__I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	.3	Control_mouse__I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	.4	Control_mouse__I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	.5	Control_mouse__I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	.6	Fasting_mouse__I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.7	Fasting_mouse__I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.8	Fasting_mouse__I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	.9	Fasting_mouse__I.D._636"""

    def test_distance_int(self):
        """ distance calculations on ints should throw no errors"""
        exp_out = "\tPC.481\tPC.607\tPC.634\tPC.635\tPC.593\tPC.636\tPC.355\tPC.354\tPC.356\nPC.481\t0.0\t798.0\t9802.0\t9802.0\t896.0\t9802.0\t9096.0\t9096.0\t9188.0\nPC.607\t798.0\t0.0\t9004.0\t9004.0\t98.0\t9004.0\t9894.0\t9894.0\t9986.0\nPC.634\t9802.0\t9004.0\t0.0\t0.0\t8906.0\t0.0\t18898.0\t18898.0\t18990.0\nPC.635\t9802.0\t9004.0\t0.0\t0.0\t8906.0\t0.0\t18898.0\t18898.0\t18990.0\nPC.593\t896.0\t98.0\t8906.0\t8906.0\t0.0\t8906.0\t9992.0\t9992.0\t10084.0\nPC.636\t9802.0\t9004.0\t0.0\t0.0\t8906.0\t0.0\t18898.0\t18898.0\t18990.0\nPC.355\t9096.0\t9894.0\t18898.0\t18898.0\t9992.0\t18898.0\t0.0\t0.0\t92.0\nPC.354\t9096.0\t9894.0\t18898.0\t18898.0\t9992.0\t18898.0\t0.0\t0.0\t92.0\nPC.356\t9188.0\t9986.0\t18990.0\t18990.0\t10084.0\t18990.0\t92.0\t92.0\t0.0"
        res_out = distance_matrix(StringIO.StringIO(self.fasting_map), "DOB")        
        self.assertEqual(exp_out, res_out)
        
    def test_distance_floats(self):
        """ distance calculations on floats should throw no errors"""
        # testing floats
        exp_out = "\tPC.481\tPC.607\tPC.634\tPC.635\tPC.593\tPC.636\tPC.355\tPC.354\tPC.356\nPC.481\t0.0\t0.2\t0.3\t0.4\t0.1\t0.5\t0.2\t0.3\t0.1\nPC.607\t0.2\t0.0\t0.1\t0.2\t0.1\t0.3\t0.4\t0.5\t0.3\nPC.634\t0.3\t0.1\t0.0\t0.1\t0.2\t0.2\t0.5\t0.6\t0.4\nPC.635\t0.4\t0.2\t0.1\t0.0\t0.3\t0.1\t0.6\t0.7\t0.5\nPC.593\t0.1\t0.1\t0.2\t0.3\t0.0\t0.4\t0.3\t0.4\t0.2\nPC.636\t0.5\t0.3\t0.2\t0.1\t0.4\t0.0\t0.7\t0.8\t0.6\nPC.355\t0.2\t0.4\t0.5\t0.6\t0.3\t0.7\t0.0\t0.1\t0.1\nPC.354\t0.3\t0.5\t0.6\t0.7\t0.4\t0.8\t0.1\t0.0\t0.2\nPC.356\t0.1\t0.3\t0.4\t0.5\t0.2\t0.6\t0.1\t0.2\t0.0"
        res_out = distance_matrix(StringIO.StringIO(self.fasting_map), "Float_Col")        
        self.assertEqual(exp_out, res_out)


#run tests if called from command line
if __name__ == '__main__':
    main()
    
    
    

