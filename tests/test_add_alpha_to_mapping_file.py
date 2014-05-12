#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"

from numpy import array, median
from unittest import TestCase, main
from qiime.add_alpha_to_mapping_file import (
    add_alpha_diversity_values_to_mapping_file,
    _get_level, mean_alpha)


class TopLevelTests(TestCase):

    def setUp(self):
        self.metrics = ['chao1', 'PD_whole_tree']
        self.alpha_diversity_data = array([[173., 6.39901], [332.5, 7.48089],
                                           [189.9375, 5.5103], [223.58333333,
                                                                6.26648], [
                                               176.8, 5.40341],
                                           [90., 4.84129], [127., 4.50866], [211., 7.3172], [146., 6.57543]])
        self.sample_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354',
                           'PC.593', 'PC.355', 'PC.607', 'PC.634']

        self.collated_alpha_dict_a = COLLATED_ALPHA_DICT_A
        self.collated_alpha_dict_b = COLLATED_ALPHA_DICT_B

        self.mapping_file_data = MAPPING_FILE_DATA
        self.mapping_file_headers = ['SampleID', 'BarcodeSequence',
                                     'LinkerPrimerSequence', 'Treatment', 'DOB', 'Description']

    def test_add_alpha_diversity_values_to_mapping_file(self):
        """checks a mapping file is added with the proper fields """

        # regular case no special cases for avg method
        expected_mapping_file_data = MAPPING_FILE_DATA_WITH_ALPHA_A
        expected_mapping_file_headers = ['SampleID', 'BarcodeSequence',
                                         'LinkerPrimerSequence', 'Treatment', 'DOB', 'Description',
                                         'chao1_alpha', 'chao1_normalized_alpha', 'chao1_alpha_label',
                                         'PD_whole_tree_alpha', 'PD_whole_tree_normalized_alpha',
                                         'PD_whole_tree_alpha_label']

        out_mapping_file_data, out_mapping_file_headers =\
            add_alpha_diversity_values_to_mapping_file(self.metrics,
                                                       self.sample_ids, self.alpha_diversity_data,
                                                       self.mapping_file_headers, self.mapping_file_data, 4, 'equal')

        self.assertEquals(out_mapping_file_data, expected_mapping_file_data)
        self.assertEquals(
            out_mapping_file_headers,
            expected_mapping_file_headers)

        # regular case no special cases for quantile method
        expected_mapping_file_data = MAPPING_FILE_DATA_WITH_ALPHA_B
        out_mapping_file_data, out_mapping_file_headers =\
            add_alpha_diversity_values_to_mapping_file(self.metrics,
                                                       self.sample_ids, self.alpha_diversity_data,
                                                       self.mapping_file_headers, self.mapping_file_data, 4, 'quantile')

        self.assertEquals(out_mapping_file_data, expected_mapping_file_data)
        self.assertEquals(
            out_mapping_file_headers,
            expected_mapping_file_headers)

    def test__get_level(self):
        """ checks the level assignment is done correctly """

        # check regular case with and without prefix tags
        expected_output = 1
        output = _get_level(0.20, [0.25, 0.5, 0.75])
        self.assertEquals(output, expected_output)

        expected_output = 'level_bin_1_of_4'
        output = _get_level(0.20, [0.25, 0.5, 0.75], 'level_bin')
        self.assertEquals(output, expected_output)

        expected_output = 'level_bin_3_of_6'
        output = _get_level(0.20, [0.05, 0.15, 0.35, 0.8, 0.95], 'level_bin')
        self.assertEquals(output, expected_output)

        # edge cases with and without prefix tags
        expected_output = 2
        output = _get_level(0.25, [0.25, 0.5, 0.75])
        self.assertEquals(output, expected_output)

        expected_output = 4
        output = _get_level(1, [0.25, 0.5, 0.75])
        self.assertEquals(output, expected_output)

        expected_output = 'testing_bin_2_of_4'
        output = _get_level(0.25, [0.25, 0.5, 0.75], 'testing_bin')
        self.assertEquals(output, expected_output)

        expected_output = 'testing_bin_4_of_4'
        output = _get_level(1, [0.25, 0.5, 0.75], 'testing_bin')
        self.assertEquals(output, expected_output)

        # unwanted cases, greater than one and negative values
        with self.assertRaises(ValueError):
            output = _get_level(1.3, [0.5])

        with self.assertRaises(ValueError):
            output = _get_level(-1, [0.25, 0.5, 0.75])

    def test_mean_alpha(self):
        """checks data is being correctly averaged"""

        # regular use-cases for this function
        expected_data = [[9.441785, 82.93],
                         [0.42877, 5.2006], [9.625995, 8.18]]
        expected_metrics = ['PD_whole_tree_even_310', 'chao1_even_310']
        expected_sample_ids = ['s1', 's2', 's3']

        o_metrics, o_sample_ids, o_data = mean_alpha(
            self.collated_alpha_dict_a, 310)

        self.assertEquals(o_metrics, expected_metrics)
        self.assertEquals(o_sample_ids, expected_sample_ids)
        self.assertEquals(o_data, expected_data)

        expected_data = [[12.508435, 11.6105],
                         [0.42877, 8.42], [11.58785, 1.0]]
        expected_metrics = ['PD_whole_tree_even_610', 'chao1_even_610']

        o_metrics, o_sample_ids, o_data = mean_alpha(
            self.collated_alpha_dict_a, 610)

        self.assertEquals(o_metrics, expected_metrics)
        self.assertEquals(o_sample_ids, expected_sample_ids)
        self.assertEquals(o_data, expected_data)

        # should default to the highest depth
        o_metrics, o_sample_ids, o_data = mean_alpha(
            self.collated_alpha_dict_a,
            None)
        self.assertEquals(o_metrics, expected_metrics)
        self.assertEquals(o_sample_ids, expected_sample_ids)
        self.assertEquals(o_data, expected_data)

        # non-existant depth
        with self.assertRaises(ValueError):
            o_metrics, o_sample_ids, o_data = mean_alpha(
                self.collated_alpha_dict_b, 111111)

        # files with non-matching sample ids should raise an exception
        with self.assertRaises(ValueError):
            o_metrics, o_sample_ids, o_data = mean_alpha(
                self.collated_alpha_dict_b, 310)

        # input types that should not be processed
        with self.assertRaises(AssertionError):
            output = mean_alpha([1, 2, 3], 5)

        with self.assertRaises(AssertionError):
            output = mean_alpha({'a': 'b'}, -1.4)


MAPPING_FILE_DATA = [
    ['PC.354',
     'AGCACGAGCCTA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._354'],
    ['PC.355',
     'AACTCGTCGATG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._355'],
    ['PC.356',
     'ACAGACCACTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061126',
     'Control_mouse_I.D._356'],
    ['PC.481',
     'ACCAGCGACTAG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20070314',
     'Control_mouse_I.D._481'],
    ['PC.593',
     'AGCAGCACTTGT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20071210',
     'Control_mouse_I.D._593'],
    ['PC.607',
     'AACTGTGCGTAC',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20071112',
     'Fasting_mouse_I.D._607'],
    ['PC.634',
     'ACAGAGTCGGCT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._634'],
    ['PC.635',
     'ACCGCAGAGTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._635'],
    ['PC.636', 'ACGGTGAGTGTC', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._636']]

MAPPING_FILE_DATA_WITH_ALPHA_A = [
    ['PC.354',
     'AGCACGAGCCTA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._354',
     '176.8',
     '0.35793814433',
     'bin_2_of_4',
     '5.40341',
     '0.301036595418',
     'bin_2_of_4'],
    ['PC.355',
     'AACTCGTCGATG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._355',
     '127.0',
     '0.152577319588',
     'bin_1_of_4',
     '4.50866',
     '0.0',
     'bin_1_of_4'],
    ['PC.356',
     'ACAGACCACTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061126',
     'Control_mouse_I.D._356',
     '189.9375',
     '0.412113402062',
     'bin_2_of_4',
     '5.5103',
     '0.336999491964',
     'bin_2_of_4'],
    ['PC.481',
     'ACCAGCGACTAG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20070314',
     'Control_mouse_I.D._481',
     '223.58333333',
     '0.550859106515',
     'bin_3_of_4',
     '6.26648',
     '0.59141452714',
     'bin_3_of_4'],
    ['PC.593',
     'AGCAGCACTTGT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20071210',
     'Control_mouse_I.D._593',
     '90.0',
     '0.0',
     'bin_1_of_4',
     '4.84129',
     '0.111912604341',
     'bin_1_of_4'],
    ['PC.607',
     'AACTGTGCGTAC',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20071112',
     'Fasting_mouse_I.D._607',
     '211.0',
     '0.498969072165',
     'bin_2_of_4',
     '7.3172',
     '0.944926873089',
     'bin_4_of_4'],
    ['PC.634',
     'ACAGAGTCGGCT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._634',
     '146.0',
     '0.230927835052',
     'bin_1_of_4',
     '6.57543',
     '0.695360049525',
     'bin_3_of_4'],
    ['PC.635',
     'ACCGCAGAGTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._635',
     '332.5',
     '1.0',
     'bin_4_of_4',
     '7.48089',
     '1.0',
     'bin_4_of_4'],
    ['PC.636', 'ACGGTGAGTGTC', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._636', '173.0', '0.342268041237', 'bin_2_of_4', '6.39901', '0.636003943167', 'bin_3_of_4']]

MAPPING_FILE_DATA_WITH_ALPHA_B = [
    ['PC.354',
     'AGCACGAGCCTA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._354',
     '176.8',
     '0.35793814433',
     'bin_3_of_4',
     '5.40341',
     '0.301036595418',
     'bin_2_of_4'],
    ['PC.355',
     'AACTCGTCGATG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061218',
     'Control_mouse_I.D._355',
     '127.0',
     '0.152577319588',
     'bin_1_of_4',
     '4.50866',
     '0.0',
     'bin_1_of_4'],
    ['PC.356',
     'ACAGACCACTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20061126',
     'Control_mouse_I.D._356',
     '189.9375',
     '0.412113402062',
     'bin_3_of_4',
     '5.5103',
     '0.336999491964',
     'bin_2_of_4'],
    ['PC.481',
     'ACCAGCGACTAG',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20070314',
     'Control_mouse_I.D._481',
     '223.58333333',
     '0.550859106515',
     'bin_4_of_4',
     '6.26648',
     '0.59141452714',
     'bin_3_of_4'],
    ['PC.593',
     'AGCAGCACTTGT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Control',
     '20071210',
     'Control_mouse_I.D._593',
     '90.0',
     '0.0',
     'bin_1_of_4',
     '4.84129',
     '0.111912604341',
     'bin_1_of_4'],
    ['PC.607',
     'AACTGTGCGTAC',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20071112',
     'Fasting_mouse_I.D._607',
     '211.0',
     '0.498969072165',
     'bin_4_of_4',
     '7.3172',
     '0.944926873089',
     'bin_4_of_4'],
    ['PC.634',
     'ACAGAGTCGGCT',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._634',
     '146.0',
     '0.230927835052',
     'bin_2_of_4',
     '6.57543',
     '0.695360049525',
     'bin_4_of_4'],
    ['PC.635',
     'ACCGCAGAGTCA',
     'YATGCTGCCTCCCGTAGGAGT',
     'Fast',
     '20080116',
     'Fasting_mouse_I.D._635',
     '332.5',
     '1.0',
     'bin_4_of_4',
     '7.48089',
     '1.0',
     'bin_4_of_4'],
    ['PC.636', 'ACGGTGAGTGTC', 'YATGCTGCCTCCCGTAGGAGT', 'Fast', '20080116', 'Fasting_mouse_I.D._636', '173.0', '0.342268041237', 'bin_2_of_4', '6.39901', '0.636003943167', 'bin_3_of_4']]

COLLATED_ALPHA_DICT_A = {
    'PD_whole_tree': ['\tsequences per sample\titeration\ts1\ts2\ts3',
                      'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
                      'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
                      'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
                      'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
                      'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
                      'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642'],
    'chao1': ['\tsequences per sample\titeration\ts1\ts2\ts3',
              'rare10.txt\t10\t0\t4.2\t3.1415\t9.11',
              'rare10.txt\t10\t1\t5.6\t3.15\t9.62',
              'rare310.txt\t310\t0\t83.11\t5.2012\t8.12',
              'rare310.txt\t310\t1\t82.75\t5.2000\t8.24',
              'rare610.txt\t610\t0\t11.11\t8.42\t1',
              'rare610.txt\t610\t1\t12.111\t8.42\t1']
}

COLLATED_ALPHA_DICT_B = {
    'PD_whole_tree': ['\tsequences per sample\titeration\ts1\ts2\ts3',
                      'rare10.txt\t10\t0\t1.99181\t0.42877\t2.13996',
                      'rare10.txt\t10\t1\t2.07163\t0.42877\t2.37055',
                      'rare310.txt\t310\t0\t8.83115\t0.42877\t11.00725',
                      'rare310.txt\t310\t1\t10.05242\t0.42877\t8.24474',
                      'rare610.txt\t610\t0\t12.03067\t0.42877\t11.58928',
                      'rare610.txt\t610\t1\t12.9862\t0.42877\t11.58642'],
    'chao1': ['\tsequences per sample\titeration\ts511\ts512\ts3',
              'rare10.txt\t10\t0\t4.2\t3.1415\t9.11',
              'rare10.txt\t10\t1\t5.6\t3.15\t9.62',
              'rare310.txt\t310\t0\t83.11\t5.2012\t8.12',
              'rare310.txt\t310\t1\t82.75\t5.2000\t8.24',
              'rare610.txt\t610\t0\t11.11\t8.42\t1',
              'rare610.txt\t610\t1\t12.111\t8.42\t1']
}

if __name__ == "__main__":
    main()
