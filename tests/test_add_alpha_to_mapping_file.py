#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"

from numpy import array, median
from numpy.random import shuffle
from cogent.util.unit_test import TestCase, main
from qiime.add_alpha_to_mapping_file import (add_alpha_diversity_values_to_mapping_file,\
                                            _get_level, quantile, _quantile,\
                                            alpha_diversity_data_to_dict,\
                                            alpha_diversity_dict_to_data)

class TopLevelTests(TestCase):
    
    def setUp(self):
        self.metrics = ['chao1', 'PD_whole_tree']
        self.alpha_diversity_data =  array([[ 173., 6.39901], [332.5, 7.48089],\
            [189.9375, 5.5103],[223.58333333, 6.26648], [176.8, 5.40341],\
            [90., 4.84129], [127., 4.50866], [211., 7.3172], [146., 6.57543]])
        self.sample_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354',\
            'PC.593', 'PC.355', 'PC.607', 'PC.634']

        self.alpha_diversity_dict = {\
            'PD_whole_tree':{\
                'PC.355': 4.5086599999999999, 'PC.607': 7.3171999999999997,\
                'PC.634': 6.5754299999999999, 'PC.635': 7.4808899999999996,\
                'PC.593': 4.8412899999999999, 'PC.636': 6.3990099999999996,\
                'PC.481': 6.2664799999999996, 'PC.354': 5.40341,\
                'PC.356': 5.5103}, 
            'chao1': {\
                'PC.355': 127.0, 'PC.607': 211.0, 'PC.634': 146.0,\
                'PC.635': 332.5, 'PC.593': 90.0, 'PC.636': 173.0,\
                'PC.481': 223.58333332999999, 'PC.354': 176.80000000000001,\
                'PC.356': 189.9375}\
            }

        self.mapping_file_data = MAPPING_FILE_DATA
        self.mapping_file_headers = ['SampleID', 'BarcodeSequence',\
            'LinkerPrimerSequence', 'Treatment', 'DOB', 'Description']

    def test_add_alpha_diversity_values_to_mapping_file(self):
        """checks a mapping file is added with the proper fields """

        # regular case no special cases
        expected_mapping_file_data = MAPPING_FILE_DATA_WITH_ALPHA
        expected_mapping_file_headers = ['SampleID', 'BarcodeSequence',\
            'LinkerPrimerSequence', 'Treatment', 'DOB', 'Description',\
            'chao1_alpha', 'chao1_normalized_alpha', 'chao1_alpha_label',\
            'PD_whole_tree_alpha', 'PD_whole_tree_normalized_alpha',\
            'PD_whole_tree_alpha_label']

        out_mapping_file_data, out_mapping_file_headers =\
            add_alpha_diversity_values_to_mapping_file(self.metrics,\
                self.sample_ids, self.alpha_diversity_data,\
                self.mapping_file_headers, self.mapping_file_data, 4, 'N/A')

        self.assertEquals(out_mapping_file_data, expected_mapping_file_data)
        self.assertEquals(out_mapping_file_headers, expected_mapping_file_headers)


    def test__get_level(self):
        """ checks the level assignment is done correctly """

        # check regular case with and without prefix tags
        expected_output = 1
        output = _get_level(0.20, 4)
        self.assertEquals(output, expected_output)

        expected_output = 'level_bin_1_of_4'
        output = _get_level(0.20, 4, 'level_bin')
        self.assertEquals(output, expected_output)

        # edge cases with and without prefix tags
        expected_output = 2
        output = _get_level(0.25, 4)
        self.assertEquals(output, expected_output)

        expected_output = 4
        output = _get_level(1, 4)
        self.assertEquals(output, expected_output)

        expected_output = 'testing_bin_2_of_4'
        output = _get_level(0.25, 4, 'testing_bin')
        self.assertEquals(output, expected_output)

        expected_output = 'testing_bin_4_of_4'
        output = _get_level(1, 4, 'testing_bin')
        self.assertEquals(output, expected_output)

        # unwanted cases, greater than one and negative values
        with self.assertRaises(AssertionError):
            output = _get_level(1.3, 5)

        with self.assertRaises(AssertionError):
            output = _get_level(-1, 4)

        with self.assertRaises(AssertionError):
            output = _get_level(0.2, 4.4)

        with self.assertRaises(AssertionError):
            output = _get_level(0.2, -1)

    def test_quantile(self):
        """checks for correct quantile statistic values"""
        
        # suffle the data to be sure, it is getting sorted
        sample_data = array(range(1, 11))
        shuffle(sample_data)

        # regular cases
        expected_output = [1.9, 2.8, 3.25, 5.5, 7.75, 7.93]
        list_of_quantiles = [0.1, 0.2, 0.25, 0.5, 0.75, 0.77]
        output = quantile(sample_data, list_of_quantiles)
        self.assertFloatEqual(expected_output, output)

        sample_data = array([42, 32, 24, 57, 15, 34, 83, 24, 60, 67, 55, 17,
            83, 17, 80, 65, 14, 34, 39, 53])
        list_of_quantiles = [0.5]
        output = quantile(sample_data, list_of_quantiles)
        self.assertFloatEqual(output, median(sample_data))

        # quantiles must be between [0, 1]
        with self.assertRaises(AssertionError):
            output = quantile(sample_data, [0.1, 0.2, -0.1, 2, 0.3, 0.5])

        # quantiles must be a list or a numpy array
        with self.assertRaises(AssertionError):
            output = quantile(sample_data, 1)

        # the data must be a list or a numpy array
        with self.assertRaises(AssertionError):
            output = quantile(1, [0])

    def test__quantile(self):
        """checks for correct quantiles according to R. type 7 algorithm"""
        # regular cases
        sample_data = array(range(25, 42))
        self.assertFloatEqual(_quantile(sample_data, 0.5), median(sample_data))

        # sorted data is assumed for this function
        sample_data = array([ 0.17483293,  0.99891939,  0.81377467,  0.8137437 ,
            0.51990174, 0.35521497,  0.98751461])
        sample_data.sort()
        self.assertFloatEqual(_quantile(sample_data, 0.10), 0.283062154)

    def test_alpha_diversity_data_to_dict(self):
        """checks the conversion from data to dict is working correctly """

        # test a case with regular inputs and outputs
        output = alpha_diversity_data_to_dict(self.metrics, self.sample_ids,\
            self.alpha_diversity_data)
        self.assertEquals(output, self.alpha_diversity_dict)

        # test a case with fewer metrics than needed
        with self.assertRaises(AssertionError):
            output = alpha_diversity_data_to_dict([], self.sample_ids,\
            self.alpha_diversity_data)

        # test a case with fewer sample_ids than needed
        with self.assertRaises(AssertionError):
            output = alpha_diversity_data_to_dict(self.metrics, [],\
            self.alpha_diversity_data)

    def test_alpha_diversity_dict_to_data(self):
        """checks the conversion from dict to data is working correctly """

        # test the common case
        out_metrics, out_sample_ids, out_data = alpha_diversity_dict_to_data(\
            self.alpha_diversity_dict)

        # the order of matters for comparison, hence, check this repeated data
        expected_data = array([[ 6.26648, 223.58333333], [6.39901, 173.],\
            [4.50866, 127.], [5.40341, 176.8], [7.3172, 211.], [6.57543, 146.],\
            [7.48089, 332.5], [4.84129, 90.], [5.5103, 189.9375]])

        self.assertEqualItems(out_metrics, self.metrics)
        self.assertEqualItems(out_sample_ids, self.sample_ids)
        self.assertEquals(out_data, expected_data)

MAPPING_FILE_DATA = [\
    ['PC.354','AGCACGAGCCTA','YATGCTGCCTCCCGTAGGAGT','Control','20061218','Control_mouse_I.D._354'],\
    ['PC.355','AACTCGTCGATG','YATGCTGCCTCCCGTAGGAGT','Control','20061218','Control_mouse_I.D._355'],\
    ['PC.356','ACAGACCACTCA','YATGCTGCCTCCCGTAGGAGT','Control','20061126','Control_mouse_I.D._356'],\
    ['PC.481','ACCAGCGACTAG','YATGCTGCCTCCCGTAGGAGT','Control','20070314','Control_mouse_I.D._481'],\
    ['PC.593','AGCAGCACTTGT','YATGCTGCCTCCCGTAGGAGT','Control','20071210','Control_mouse_I.D._593'],\
    ['PC.607','AACTGTGCGTAC','YATGCTGCCTCCCGTAGGAGT','Fast','20071112','Fasting_mouse_I.D._607'],\
    ['PC.634','ACAGAGTCGGCT','YATGCTGCCTCCCGTAGGAGT','Fast','20080116','Fasting_mouse_I.D._634'],\
    ['PC.635','ACCGCAGAGTCA','YATGCTGCCTCCCGTAGGAGT','Fast','20080116','Fasting_mouse_I.D._635'],\
    ['PC.636','ACGGTGAGTGTC','YATGCTGCCTCCCGTAGGAGT','Fast','20080116','Fasting_mouse_I.D._636']\
]

MAPPING_FILE_DATA_WITH_ALPHA = [\
['PC.354', 'AGCACGAGCCTA', 'YATGCTGCCTCCCGTAGGAGT', 'Control', '20061218', 'Control_mouse_I.D._354', '176.8', '0.35793814433', 'bin_2_of_4', '5.40341', '0.301036595418', 'bin_2_of_4'],\
['PC.355', 'AACTCGTCGATG', 'YATGCTGCCTCCCGTAGGAGT', 'Control', '20061218', 'Control_mouse_I.D._355', '127.0', '0.152577319588', 'bin_1_of_4', '4.50866', '0.0', 'bin_1_of_4'],\
['PC.356', 'ACAGACCACTCA', 'YATGCTGCCTCCCGTAGGAGT', 'Control', '20061126', 'Control_mouse_I.D._356', '189.9375', '0.412113402062', 'bin_2_of_4', '5.5103', '0.336999491964', 'bin_2_of_4'],\
['PC.481', 'ACCAGCGACTAG', 'YATGCTGCCTCCCGTAGGAGT', 'Control', '20070314', 'Control_mouse_I.D._481', '223.58333333', '0.550859106515', 'bin_3_of_4', '6.26648', '0.59141452714', 'bin_3_of_4'],\
['PC.593', 'AGCAGCACTTGT', 'YATGCTGCCTCCCGTAGGAGT', 'Control', '20071210', 'Control_mouse_I.D._593', '90.0', '0.0', 'bin_1_of_4', '4.84129', '0.111912604341', 'bin_1_of_4'],\
['PC.607', 'AACTGTGCGTAC', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20071112', 'Fasting_mouse_I.D._607', '211.0', '0.498969072165', 'bin_2_of_4', '7.3172', '0.944926873089', 'bin_4_of_4'],\
['PC.634', 'ACAGAGTCGGCT', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._634', '146.0', '0.230927835052', 'bin_1_of_4', '6.57543', '0.695360049525', 'bin_3_of_4'],\
['PC.635', 'ACCGCAGAGTCA', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._635', '332.5', '1.0', 'bin_4_of_4', '7.48089', '1.0', 'bin_4_of_4'],\
['PC.636', 'ACGGTGAGTGTC', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._636', '173.0', '0.342268041237', 'bin_2_of_4', '6.39901', '0.636003943167', 'bin_3_of_4']]


if __name__ == "__main__":
    main()