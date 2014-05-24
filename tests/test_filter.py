#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from StringIO import StringIO
from tempfile import mkstemp
from os import close

from numpy import inf
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from skbio.util.misc import remove_files
from biom.parse import Table
from qiime.parse import (parse_distmat, parse_mapping_file,
                         parse_metadata_state_descriptions)
from qiime.filter import (filter_fasta, filter_samples_from_otu_table,
                          filter_otus_from_otu_table, get_sample_ids,
                          sample_ids_from_category_state_coverage,
                          filter_samples_from_distance_matrix,
                          negate_tips_to_keep,
                          filter_mapping_file, filter_tree,
                          filter_fastq, filter_otus_from_otu_map,
                          filter_otu_table_to_n_samples,
                          filter_mapping_file_from_mapping_f,
                          filter_mapping_file_by_metadata_states,
                          get_otu_ids_from_taxonomy_f,
                          sample_ids_from_metadata_description)
from qiime.test import FakeFile
from qiime.util import load_qiime_config


class fake_output_f():

    def __init__(self):
        self.s = ""

    def write(self, s):
        self.s += s

    def close(self):
        pass


class FilterTests(TestCase):

    def setUp(self):
        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir']
        self.files_to_remove = []

        self.filter_fasta_expected1 = filter_fasta_expected1
        self.filter_fasta_expected2 = filter_fasta_expected2
        self.filter_fastq_expected1 = filter_fastq_expected1
        self.filter_fastq_expected2 = filter_fastq_expected2
        self.input_dm1 = input_dm1.split('\n')
        self.expected_dm1a = expected_dm1a.split('\n')
        self.expected_dm1b = expected_dm1b.split('\n')
        self.map_str1 = map_str1
        self.map_str2 = map_str2.split('\n')
        self.map_data, self.map_headers, self.map_comments =\
            parse_mapping_file(StringIO(self.map_str1))
        self.tree1 = DndParser(tree1)
        self.tree2 = DndParser(tree2)
        self.tutorial_mapping_f = FakeFile(tutorial_mapping_f)

        self.otu_table2 = Table.from_json(sparse_otu_table2)

        # For sample_ids_from_category_state_coverage() tests.
        self.exp_empty = (set([]), 0, set([]))
        self.exp_all = (set(['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593',
                             'PC.607', 'PC.634', 'PC.635', 'PC.636']), 6,
                        set(['Control', 'Fast']))

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_negate_tips_to_keep(self):
        """ negate_tips_to_keep functions as expected """
        t = DndParser("((S5:0.00014,S7:0.00015)0.752:0.45762,(S3:0.00014,"
                      "seq6:0.00014)0.180:0.00015,(Seq1:0.00014,s2:0.00014)0.528:1.0466);")

        tips_to_keep = ["S5", "Seq1", "s2"]
        expected = ["S7", "S3", "seq6"]
        self.assertItemsEqual(negate_tips_to_keep(tips_to_keep, t), expected)

        tips_to_keep = ["S5", "Seq1"]
        expected = ["S7", "S3", "seq6", "s2"]
        self.assertItemsEqual(negate_tips_to_keep(tips_to_keep, t), expected)

        tips_to_keep = []
        expected = ["S7", "S3", "seq6", "s2", "S5", "Seq1"]
        self.assertItemsEqual(negate_tips_to_keep(tips_to_keep, t), expected)

        tips_to_keep = ["S7", "S3", "seq6", "s2", "S5", "Seq1"]
        expected = []
        self.assertItemsEqual(negate_tips_to_keep(tips_to_keep, t), expected)

    def test_filter_mapping_file(self):
        """filter_mapping_file should filter map file according to sample ids"""
        self.assertEqual(filter_mapping_file(self.map_data, self.map_headers,
                                             ['a', 'b', 'c', 'd', 'e', 'f']), (self.map_headers, self.map_data))
        self.assertEqual(
            filter_mapping_file(self.map_data, self.map_headers, ['a']),
            (['SampleID', 'Description'], ['a\tx'.split('\t')]))

    def test_filter_mapping_file_from_mapping_f(self):
        """ filter_mapping_file_from_mapping_f functions as expected """
        actual = filter_mapping_file_from_mapping_f(
            self.tutorial_mapping_f, ["PC.354", "PC.355"])
        expected = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355"""
        self.assertEqual(actual, expected)

    def test_filter_mapping_file_from_mapping_f_negate(self):
        """ filter_mapping_file_from_mapping_f functions as expected when negate is True """
        actual = filter_mapping_file_from_mapping_f(self.tutorial_mapping_f,
                                                    ["PC.356",
                                                     "PC.481",
                                                     "PC.593",
                                                     "PC.607",
                                                     "PC.634",
                                                     "PC.635",
                                                     "PC.636"],
                                                    negate=True)
        expected = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355"""
        self.assertEqual(actual, expected)

    def test_filter_mapping_file_by_metadata_states(self):
        """ filter_mapping_file_by_metadata_states functions as expected """
        actual = filter_mapping_file_by_metadata_states(
            self.tutorial_mapping_f,
            "Treatment:Control")
        expected = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593"""
        self.assertEqual(actual, expected)

    def test_filter_fasta(self):
        """filter_fasta functions as expected"""
        input_seqs = [('Seq1 some comment', 'ACCTTGG'),
                      ('s2 some other comment', 'TTGG'),
                      ('S3', 'AAGGCCGG'),
                      ('S5 some comment', 'CGT'),
                      ('seq6 some other comment', 'AA'),
                      ('S7', 'T')]
        seqs_to_keep = {}.fromkeys(['Seq1',
                                    's2 some other comment',
                                    'S3 no comment'])

        actual = fake_output_f()
        filter_fasta(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=False)
        self.assertEqual(actual.s, self.filter_fasta_expected1)

        actual = fake_output_f()
        filter_fasta(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=True)
        self.assertEqual(actual.s, self.filter_fasta_expected2)

    def test_filter_fastq(self):
        """filter_fastq functions as expected"""
        input_seqs = [('Seq1 some comment', 'ACCTTGG', 'BBBBBBB'),
                      ('s2 some other comment', 'TTGG', 'BBBB'),
                      ('S3', 'AAGGCCGG', 'BBCtatcc'),
                      ('S5 some comment', 'CGT', 'BBB'),
                      ('seq6 some other comment', 'AA', 'BB'),
                      ('S7', 'T', 's')]
        seqs_to_keep = {}.fromkeys(['Seq1',
                                    's2 some other comment',
                                    'S3 no comment'])

        actual = fake_output_f()
        filter_fastq(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=False)
        self.assertEqual(actual.s, self.filter_fastq_expected1)

        actual = fake_output_f()
        filter_fastq(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=True)
        self.assertEqual(actual.s, self.filter_fastq_expected2)

    def test_filter_tree(self):
        """filter_tree functions as expected"""
        actual = [e.Name for e in filter_tree(
            self.tree1, ['bbb', 'ccc']).tips()]
        #(a_a:10,(b_b:2,c_c:4):5);
        expected = [
            e.Name for e in DndParser(
                "((bbb:2,ccc:4));",
                PhyloNode).tips(
            )]
        self.assertEqual(actual, expected)

        actual = [e.Name for e in filter_tree(
            self.tree2, ['bbb', 'ccc']).tips()]
        #(a_a:10,(b_b:2,c_c:4):5);
        expected = [
            e.Name for e in DndParser(
                "(('bbb':2,ccc:4));",
                PhyloNode).tips(
            )]
        self.assertEqual(actual, expected)

    def test_filter_samples_from_distance_matrix(self):
        """filter_samples_from_distance_matrix functions as expected """
        actual = filter_samples_from_distance_matrix(
            parse_distmat(self.input_dm1),
            ["GHI blah", "XYZ"])
        self.assertEqual(actual, expected_dm1a)
        actual = filter_samples_from_distance_matrix(
            parse_distmat(self.input_dm1),
            ["GHI", "DEF"])
        self.assertEqual(actual, expected_dm1b)

    def test_filter_samples_from_distance_matrix_file_input(self):
        """filter_samples_from_distance_matrix handles file input """
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                                     ["GHI blah", "XYZ"])
        self.assertEqual(actual, expected_dm1a)
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                                     ["GHI", "DEF"])
        self.assertEqual(actual, expected_dm1b)

    def test_filter_samples_from_distance_matrix_negate(self):
        """filter_samples_from_distance_matrix functions w negate """
        actual = filter_samples_from_distance_matrix(
            parse_distmat(self.input_dm1),
            ["ABC blah", "DEF"],
            negate=True)
        self.assertEqual(actual, expected_dm1a)
        actual = filter_samples_from_distance_matrix(
            parse_distmat(self.input_dm1),
            ["ABC", "XYZ"],
            negate=True)
        self.assertEqual(actual, expected_dm1b)

    def test_get_otu_ids_from_taxonomy_f(self):
        """get_otu_ids_from_taxonomy_f returns functions that work as expected"""
        # positive list only
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['a'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['b'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(['d'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['d', 'a'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['b', 'a'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        # only full-length match works
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(['c'])([], 42, {'taxonomy': ['a', 'b', 'cc']}))

        # negative list only
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(None, ['d'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(None, ['x'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(None, ['a'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(None, ['b'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(None, ['c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(get_otu_ids_from_taxonomy_f(
            None, ['x', 'c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(get_otu_ids_from_taxonomy_f(
            None, ['b', 'c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))

        # positive and negative list
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['a'], ['d'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(['a'], ['c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(get_otu_ids_from_taxonomy_f(
            ['a'], ['d', 'c'])([], 42, {'taxonomy': ['a', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(['x'], ['y'])([], 42, {'taxonomy': ['a', 'b', 'c']}))

        # invalid input
        self.assertRaises(
            ValueError,
            get_otu_ids_from_taxonomy_f,
            ['x'],
            ['x'])
        self.assertRaises(KeyError, get_otu_ids_from_taxonomy_f(
            ['x'], metadata_field='x'), [], 42, {'taxonomy': ['a', 'b', 'c']})

        # alt metadata field
        self.assertTrue(get_otu_ids_from_taxonomy_f(
            ['b', 'a'], metadata_field='x')([], 42, {'x': ['a', 'b', 'c']}))

        # case insensitive
        self.assertTrue(
            get_otu_ids_from_taxonomy_f(['a'])([], 42, {'taxonomy': ['A', 'b', 'c']}))
        self.assertFalse(
            get_otu_ids_from_taxonomy_f(None, ['B'])([], 42, {'taxonomy': ['a', 'b', 'c']}))

    def test_filter_otus_from_otu_table_ids(self):
        """filter_otus_from_otu_table functions with list of OTU ids"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)
        filtered_otu_table = filter_otus_from_otu_table(otu_table,
                                                        set(otu_table.observation_ids) - set(['34', '155', '152']), 0, inf, 0, inf)
        expected_otu_ids = set(otu_table.observation_ids) - \
            set(['34', '155', '152'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_otus_from_otu_table_ids_negate(self):
        """filter_otus_from_otu_table functions with list of OTU ids and negate option"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)
        filtered_otu_table = filter_otus_from_otu_table(otu_table,
                                                        set(otu_table.observation_ids) - set(['34', '155', '152']), 0, inf, 0, inf, negate_ids_to_keep=True)
        expected_otu_ids = set(['34', '155', '152'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_otus_from_otu_table_counts_dense(self):
        """filter_otus_from_otu_table functions with count-based filtering (dense OTU table)"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)

        # min and max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            20,
            25,
            0,
            inf)
        expected_otu_ids = set(['34', '155', '152'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            43,
            inf,
            0,
            inf)
        expected_otu_ids = set(['267', '154', '254', '17'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no min
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            1,
            0,
            inf)
        expected_otu_ids = set(
            ['0', '1', '10', '100', '102', '104', '105', '106', '107', '108',
             '11', '111', '112', '113', '114', '115', '116', '118', '119', '12', '121', '123', '124',
             '125', '127', '128', '129', '132', '133', '134', '135', '136', '137', '138', '139',
             '141', '142', '143', '144', '148', '149', '15', '150', '157', '160', '161', '163', '164',
             '166', '167', '168', '170', '171', '172', '173', '175', '176', '177', '179', '18', '180',
             '182', '183', '185', '186', '188', '189', '19', '190', '192', '193', '195', '197', '2',
             '20', '202', '205', '206', '207', '209', '210', '212', '214', '215', '216', '219', '221',
             '222', '224', '226', '230', '232', '233', '234', '237', '238', '239', '24', '240', '242',
             '243', '244', '246', '247', '249', '25', '252', '255', '256', '258', '259', '260', '261',
             '263', '264', '268', '269', '27', '270', '271', '272', '273', '274', '275', '276', '277',
             '278', '279', '28', '280', '281', '284', '285', '288', '291', '292', '293', '294', '296',
             '297', '298', '30', '300', '302', '303', '304', '305', '306', '307', '308', '309', '31',
             '310', '311', '312', '314', '316', '317', '318', '32', '320', '321', '322', '323', '324',
             '325', '327', '328', '33', '330', '331', '332', '334', '335', '336', '337', '338', '339',
             '342', '343', '344', '345', '346', '347', '348', '350', '354', '355', '356', '358', '359',
             '364', '366', '367', '368', '369', '37', '372', '374', '376', '377', '378', '379', '38',
             '380', '382', '384', '385', '386', '387', '388', '389', '39', '390', '391', '392', '393',
             '394', '397', '398', '4', '40', '400', '401', '402', '403', '404', '405', '406', '410',
             '411', '413', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '55', '56',
             '57', '59', '6', '60', '62', '64', '66', '67', '68', '69', '70', '71', '72', '74', '76',
             '77', '80', '81', '85', '86', '88', '89', '91', '92', '94', '97', '98', '99'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_otus_from_otu_table_counts_sparse(self):
        """filter_otus_from_otu_table functions with count-based filtering (sparse OTU table)"""
        otu_table = Table(sparse_otu_table1)

        # min and max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            20,
            25,
            0,
            inf)
        expected_otu_ids = set(['34', '155', '152'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            43,
            inf,
            0,
            inf)
        expected_otu_ids = set(['267', '154', '254', '17'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no min
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            1,
            0,
            inf)
        expected_otu_ids = set(
            ['0', '1', '10', '100', '102', '104', '105', '106', '107', '108',
             '11', '111', '112', '113', '114', '115', '116', '118', '119', '12', '121', '123', '124',
             '125', '127', '128', '129', '132', '133', '134', '135', '136', '137', '138', '139',
             '141', '142', '143', '144', '148', '149', '15', '150', '157', '160', '161', '163', '164',
             '166', '167', '168', '170', '171', '172', '173', '175', '176', '177', '179', '18', '180',
             '182', '183', '185', '186', '188', '189', '19', '190', '192', '193', '195', '197', '2',
             '20', '202', '205', '206', '207', '209', '210', '212', '214', '215', '216', '219', '221',
             '222', '224', '226', '230', '232', '233', '234', '237', '238', '239', '24', '240', '242',
             '243', '244', '246', '247', '249', '25', '252', '255', '256', '258', '259', '260', '261',
             '263', '264', '268', '269', '27', '270', '271', '272', '273', '274', '275', '276', '277',
             '278', '279', '28', '280', '281', '284', '285', '288', '291', '292', '293', '294', '296',
             '297', '298', '30', '300', '302', '303', '304', '305', '306', '307', '308', '309', '31',
             '310', '311', '312', '314', '316', '317', '318', '32', '320', '321', '322', '323', '324',
             '325', '327', '328', '33', '330', '331', '332', '334', '335', '336', '337', '338', '339',
             '342', '343', '344', '345', '346', '347', '348', '350', '354', '355', '356', '358', '359',
             '364', '366', '367', '368', '369', '37', '372', '374', '376', '377', '378', '379', '38',
             '380', '382', '384', '385', '386', '387', '388', '389', '39', '390', '391', '392', '393',
             '394', '397', '398', '4', '40', '400', '401', '402', '403', '404', '405', '406', '410',
             '411', '413', '42', '43', '44', '45', '46', '47', '48', '49', '5', '50', '51', '55', '56',
             '57', '59', '6', '60', '62', '64', '66', '67', '68', '69', '70', '71', '72', '74', '76',
             '77', '80', '81', '85', '86', '88', '89', '91', '92', '94', '97', '98', '99'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_otus_from_otu_table_samples_sparse(self):
        """filter_otus_from_otu_table functions with sample-based filtering (sparse OTU table)"""
        otu_table = Table(sparse_otu_table1)

        # min and max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            6,
            7)
        expected_otu_ids = set(
            ['153',
             '154',
             '203',
             '286',
             '353',
             '191',
             '227'])

        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            6,
            inf)
        expected_otu_ids = set(
            ['153',
             '154',
             '203',
             '286',
             '353',
             '191',
             '227',
             '120'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no min
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            0,
            1)
        expected_otu_ids = set(
            ['0', '1', '2', '4', '5', '6', '9', '10', '11', '12', '15', '18', '19',
             '20', '24', '25', '27', '28', '30', '31', '32', '33', '37', '38', '39',
             '40', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52',
             '54', '55', '56', '57', '59', '60', '62', '64', '66', '67', '68', '69',
             '70', '71', '72', '73', '74', '76', '77', '80', '81', '82', '85', '86',
             '88', '89', '91', '92', '94', '95', '97', '98', '99', '100', '101', '102',
             '104', '105', '106', '107', '108', '110', '111', '112', '113', '114', '115',
             '116', '118', '119', '121', '123', '124', '125', '127', '128', '129', '132',
             '133', '134', '135', '136', '137', '138', '139', '141', '142', '143', '144',
             '145', '148', '149', '150', '157', '160', '161', '163', '164', '166', '167',
             '168', '170', '171', '172', '173', '175', '176', '177', '178', '179', '180',
             '182', '183', '185', '186', '188', '189', '190', '192', '193', '194', '195',
             '197', '200', '202', '205', '206', '207', '209', '210', '212', '213', '214',
             '215', '216', '219', '221', '222', '224', '226', '230', '232', '233', '234',
             '236', '237', '238', '239', '240', '241', '242', '243', '244', '246', '247',
             '249', '252', '255', '256', '257', '258', '259', '260', '261', '263', '264',
             '265', '266', '268', '269', '270', '271', '272', '273', '274', '275', '276',
             '277', '278', '279', '280', '281', '282', '284', '285', '288', '289', '290',
             '291', '292', '293', '294', '296', '297', '298', '300', '301', '302', '303',
             '304', '305', '306', '307', '308', '309', '310', '311', '312', '314', '316',
             '317', '318', '320', '321', '322', '323', '324', '325', '327', '328', '330',
             '331', '332', '334', '335', '336', '337', '338', '339', '340', '342', '343',
             '344', '345', '346', '347', '348', '350', '352', '354', '355', '356', '357',
             '358', '359', '364', '365', '366', '367', '368', '369', '372', '374', '375',
             '376', '377', '378', '379', '380', '381', '382', '384', '385', '386', '387',
             '388', '389', '390', '391', '392', '393', '394', '396', '397', '398', '399',
             '400', '401', '402', '403', '404', '405', '406', '410', '411', '412', '413'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_otus_from_otu_table_samples_dense(self):
        """filter_otus_from_otu_table functions with sample-based filtering (dense OTU table)"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)

        # min and max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            6,
            7)
        expected_otu_ids = set(
            ['153',
             '154',
             '203',
             '286',
             '353',
             '191',
             '227'])

        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no max
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            6,
            inf)
        expected_otu_ids = set(
            ['153',
             '154',
             '203',
             '286',
             '353',
             '191',
             '227',
             '120'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)
        # no min
        filtered_otu_table = filter_otus_from_otu_table(
            otu_table,
            otu_table.observation_ids,
            0,
            inf,
            0,
            1)
        expected_otu_ids = set(
            ['0', '1', '2', '4', '5', '6', '9', '10', '11', '12', '15', '18', '19',
             '20', '24', '25', '27', '28', '30', '31', '32', '33', '37', '38', '39',
             '40', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52',
             '54', '55', '56', '57', '59', '60', '62', '64', '66', '67', '68', '69',
             '70', '71', '72', '73', '74', '76', '77', '80', '81', '82', '85', '86',
             '88', '89', '91', '92', '94', '95', '97', '98', '99', '100', '101', '102',
             '104', '105', '106', '107', '108', '110', '111', '112', '113', '114', '115',
             '116', '118', '119', '121', '123', '124', '125', '127', '128', '129', '132',
             '133', '134', '135', '136', '137', '138', '139', '141', '142', '143', '144',
             '145', '148', '149', '150', '157', '160', '161', '163', '164', '166', '167',
             '168', '170', '171', '172', '173', '175', '176', '177', '178', '179', '180',
             '182', '183', '185', '186', '188', '189', '190', '192', '193', '194', '195',
             '197', '200', '202', '205', '206', '207', '209', '210', '212', '213', '214',
             '215', '216', '219', '221', '222', '224', '226', '230', '232', '233', '234',
             '236', '237', '238', '239', '240', '241', '242', '243', '244', '246', '247',
             '249', '252', '255', '256', '257', '258', '259', '260', '261', '263', '264',
             '265', '266', '268', '269', '270', '271', '272', '273', '274', '275', '276',
             '277', '278', '279', '280', '281', '282', '284', '285', '288', '289', '290',
             '291', '292', '293', '294', '296', '297', '298', '300', '301', '302', '303',
             '304', '305', '306', '307', '308', '309', '310', '311', '312', '314', '316',
             '317', '318', '320', '321', '322', '323', '324', '325', '327', '328', '330',
             '331', '332', '334', '335', '336', '337', '338', '339', '340', '342', '343',
             '344', '345', '346', '347', '348', '350', '352', '354', '355', '356', '357',
             '358', '359', '364', '365', '366', '367', '368', '369', '372', '374', '375',
             '376', '377', '378', '379', '380', '381', '382', '384', '385', '386', '387',
             '388', '389', '390', '391', '392', '393', '394', '396', '397', '398', '399',
             '400', '401', '402', '403', '404', '405', '406', '410', '411', '412', '413'])
        self.assertEqual(
            set(filtered_otu_table.observation_ids),
            expected_otu_ids)

    def test_filter_samples_from_otu_table_counts_dense(self):
        """filter_samples_from_otu_table functions with count-based filtering (dense OTU table)"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)

        # min and max
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            148,
            149)
        expected_sample_ids = set(['PC.354', 'PC.635', 'PC.593', 'PC.607'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)
        # min only
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            148,
            inf)
        expected_sample_ids = set(
            ['PC.354',
             'PC.635',
             'PC.593',
             'PC.607',
             'PC.356',
             'PC.634'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)
        # max only
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            0,
            149)
        expected_sample_ids = set(
            ['PC.355',
             'PC.481',
             'PC.636',
             'PC.354',
             'PC.635',
             'PC.593',
             'PC.607'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

    def test_filter_samples_from_otu_table_sample_ids_dense(self):
        """filter_samples_from_otu_table functions with count-based filtering (dense OTU table)"""
        otu_table = Table(dense_otu_table1, input_is_dense=True)

        # keep two samples
        expected_sample_ids = set(['PC.593', 'PC.607'])
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            expected_sample_ids,
            0,
            inf)
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

        # keep some other samples
        expected_sample_ids = set(['PC.354', 'PC.635', 'PC.593'])
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            expected_sample_ids,
            0,
            inf)
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

    def test_filter_otu_table_to_n_samples(self):
        """filter_otu_table_to_n_samples returns randomly selected subset of samples
        """
        otu_table = Table(dense_otu_table1, input_is_dense=True)

        # keep two samples
        filtered_otu_table = filter_otu_table_to_n_samples(otu_table, 2)
        self.assertEqual(len(filtered_otu_table.sample_ids), 2)

        # keep three samples
        filtered_otu_table = filter_otu_table_to_n_samples(otu_table, 3)
        self.assertEqual(len(filtered_otu_table.sample_ids), 3)

        # ValueError on invalid n
        self.assertRaises(
            ValueError,
            filter_otu_table_to_n_samples,
            otu_table,
            -1)
        self.assertRaises(
            ValueError,
            filter_otu_table_to_n_samples,
            otu_table,
            10)

        # Multiple iterations yield different OTU tables - check that in 100 iterations
        # we get at least two different results to avoid random test failures.
        results = []
        for i in range(100):
            filtered_otu_table = filter_otu_table_to_n_samples(otu_table, 3)
            results.append(tuple(filtered_otu_table.sample_ids))
        self.assertTrue(len(set(results)) > 1)

    def test_filter_samples_from_otu_table_counts_sparse(self):
        """filter_samples_from_otu_table functions with count-based filtering (sparse OTU table)"""
        otu_table = Table(sparse_otu_table1)

        # min and max
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            148,
            149)
        expected_sample_ids = set(['PC.354', 'PC.635', 'PC.593', 'PC.607'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)
        # min only
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            148,
            inf)
        expected_sample_ids = set(
            ['PC.354',
             'PC.635',
             'PC.593',
             'PC.607',
             'PC.356',
             'PC.634'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)
        # max only
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            otu_table.sample_ids,
            0,
            149)
        expected_sample_ids = set(
            ['PC.355',
             'PC.481',
             'PC.636',
             'PC.354',
             'PC.635',
             'PC.593',
             'PC.607'])
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

    def test_filter_samples_from_otu_table_sample_ids_sparse(self):
        """filter_samples_from_otu_table functions with count-based filtering (sparse OTU table)"""
        otu_table = Table(sparse_otu_table1)

        # keep two samples
        expected_sample_ids = set(['PC.593', 'PC.607'])
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            expected_sample_ids,
            0,
            inf)
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

        # keep some other samples
        expected_sample_ids = set(['PC.354', 'PC.635', 'PC.593'])
        filtered_otu_table = filter_samples_from_otu_table(
            otu_table,
            expected_sample_ids,
            0,
            inf)
        self.assertEqual(
            set(filtered_otu_table.sample_ids),
            expected_sample_ids)

    def test_sample_ids_from_metadata_description(self):
        """Testing sample_ids_from_metadata_description fails on an empty set"""
        self.assertRaises(ValueError, sample_ids_from_metadata_description,
                          self.tutorial_mapping_f, "Treatment:Foo")
        self.assertRaises(ValueError, sample_ids_from_metadata_description,
                          self.tutorial_mapping_f, "DOB:!20061218,!20070314,!20071112,"
                          "!20080116")

    def test_get_sample_ids(self):
        """get_sample_ids should return sample ids matching criteria."""
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,
                                        parse_metadata_state_descriptions('Study:Twin')), [])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,
                                        parse_metadata_state_descriptions('Study:Dog')), ['a', 'b'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,
                                        parse_metadata_state_descriptions('Study:*,!Dog')), ['c', 'd', 'e'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,
                                        parse_metadata_state_descriptions('Study:*,!Dog;BodySite:Stool')), ['e'])
        self.assertEqual(get_sample_ids(self.map_data, self.map_headers,
                                        parse_metadata_state_descriptions('BodySite:Stool')), ['a', 'b', 'e'])

    def test_sample_ids_from_category_state_coverage_min_num_states(self):
        """Test returns samp IDs based on number of states that are covered."""
        # Filter out all samples.
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=2)
        self.assertEqual(obs, self.exp_empty)

        # Don't filter out any samples.
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1)
        self.assertEqual(obs, self.exp_all)

        # Filter out some samples.
        exp = (set(['d', 'e']), 1, set(['Palm', 'Stool']))
        obs = sample_ids_from_category_state_coverage(
            self.map_str1.split('\n'),
            'BodySite', 'Study', min_num_states=2)
        self.assertEqual(obs, exp)

        # Keep subject that has more than specified minimum number and has more
        # than one sample at a single coverage state.
        exp = (set(['a', 'c', 'd', 'e', 'f', 'g']), 2, set(['1', '3', '2']))
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', min_num_states=2)
        self.assertEqual(obs, exp)

    def test_sample_ids_from_category_state_coverage_min_num_states_w_considered_states(
            self):
        """Test returns samp IDs based on number of considered states that are covered."""

        # Filter out all samples
        # min_num_states too high
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=2,
                                                      considered_states=["Control", "Fast"])
        self.assertEqual(obs, self.exp_empty)
        # considered_states too restrictive
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      considered_states=[])
        self.assertEqual(obs, self.exp_empty)

        # Don't filter out any samples.
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      considered_states=["Control", "Fast"])
        self.assertEqual(obs, self.exp_all)

        # Some samples filtered when considered states is partially restrictive
        exp = (set(['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']), 4,
               set(['Control']))
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      considered_states=["Control"])
        self.assertEqual(obs, exp)

        exp = (set(['PC.607', 'PC.634', 'PC.635', 'PC.636']), 2, set(['Fast']))
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      considered_states=["Fast"])
        self.assertEqual(obs, exp)

    def test_sample_ids_from_category_state_coverage_required_states(self):
        """Test returns samp IDs based on specific category states covered."""
        # Filter out all samples.
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', required_states=['Control', 'Fast'])
        self.assertEqual(obs, self.exp_empty)

        # Don't filter out any samples.
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', required_states=[])
        self.assertEqual(obs, self.exp_all)

        # Filter out some samples.
        exp = (set(['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']), 4,
               set(['Control']))
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', required_states=['Control'])
        self.assertEqual(obs, exp)

        exp = (set(['d', 'e']), 1, set(['Palm', 'Stool']))
        obs = sample_ids_from_category_state_coverage(
            self.map_str1.split('\n'),
            'BodySite', 'Study', required_states=['Stool', 'Palm'])
        self.assertEqual(obs, exp)

        # Keep subject that has more than specified covered states and has more
        # than one sample at a single coverage state.
        exp = (set(['c', 'f', 'a', 'g']), 1, set(['1', '2', '3']))
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', required_states=['3', '2'])
        self.assertEqual(obs, exp)

        # Should convert states to strings.
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', required_states=[3, 2])
        self.assertEqual(obs, exp)

    def test_sample_ids_from_category_state_combined_filters(self):
        """Test returns samp IDs using both supported filters."""
        # Filter out all samples (fails both filters).
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=2,
                                                      required_states=['Control', 'Fast'])
        self.assertEqual(obs, self.exp_empty)

        # Filter out all samples (fails one filter).
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=2, required_states=['Control'])
        self.assertEqual(obs, self.exp_empty)

        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      required_states=['Control', 'Fast'])
        self.assertEqual(obs, self.exp_empty)

        # Don't filter out any samples (passes both filters).
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=0, required_states=[])
        self.assertEqual(obs, self.exp_all)

        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1, required_states=[])
        self.assertEqual(obs, self.exp_all)

        # Filter out some samples.
        exp = (set(['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']), 4,
               set(['Control']))
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1, required_states=['Control'])
        self.assertEqual(obs, exp)

        exp = (set(['PC.607', 'PC.634', 'PC.635', 'PC.636']), 2, set(['Fast']))
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1, required_states=['Fast'])
        self.assertEqual(obs, exp)

        exp = (set(['d', 'e']), 1, set(['Palm', 'Stool']))
        obs = sample_ids_from_category_state_coverage(
            self.map_str1.split('\n'),
            'BodySite', 'Study', required_states=['Stool', 'Palm'],
            min_num_states=1)
        self.assertEqual(obs, exp)

        # Keep subject that has more than specified covered states and has more
        # than one sample at a single coverage state (i.e. timepoint).
        exp = (set(['c', 'f', 'a', 'g']), 1, set(['1', '2', '3']))
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', min_num_states=3, required_states=['3', '2'])
        self.assertEqual(obs, exp)

        # Test filtering out the subject (from the above test) that has four
        # timepoints, but only three are unique.
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', min_num_states=4, required_states=['3', '2'])
        self.assertEqual(obs, self.exp_empty)

    def test_sample_ids_from_category_state_splitter_category(self):
        """Test correctly filters using a category to split sample IDs on."""
        # Split on category that has only one state (essentially not splitting
        # at all).
        exp = {'YATGCTGCCTCCCGTAGGAGT':
               (set(['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593']), 4,
                set(['Control']))}
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      required_states=[
                                                          'Control'],
                                                      splitter_category='LinkerPrimerSequence')
        self.assertEqual(obs, exp)

        # Split on category that has a unique state for each sample ID.
        exp = {'ACCGCAGAGTCA': (set(['PC.635']), 1, set(['Fast'])),
               'ACCAGCGACTAG': (set(['PC.481']), 1, set(['Control'])),
               'AACTGTGCGTAC': (set(['PC.607']), 1, set(['Fast'])),
               'AGCAGCACTTGT': (set(['PC.593']), 1, set(['Control'])),
               'ACAGAGTCGGCT': (set(['PC.634']), 1, set(['Fast'])),
               'AACTCGTCGATG': (set(['PC.355']), 1, set(['Control'])),
               'ACGGTGAGTGTC': (set(['PC.636']), 1, set(['Fast'])),
               'AGCACGAGCCTA': (set(['PC.354']), 1, set(['Control'])),
               'ACAGACCACTCA': (set(['PC.356']), 1, set(['Control']))}
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=1,
                                                      splitter_category='BarcodeSequence')
        self.assertEqual(obs, exp)

        # Filter all samples.
        exp = {'ACCGCAGAGTCA': self.exp_empty,
               'ACCAGCGACTAG': self.exp_empty,
               'AACTGTGCGTAC': self.exp_empty,
               'AGCAGCACTTGT': self.exp_empty,
               'ACAGAGTCGGCT': self.exp_empty,
               'AACTCGTCGATG': self.exp_empty,
               'ACGGTGAGTGTC': self.exp_empty,
               'AGCACGAGCCTA': self.exp_empty,
               'ACAGACCACTCA': self.exp_empty}
        obs = sample_ids_from_category_state_coverage(self.tutorial_mapping_f,
                                                      'Treatment', 'DOB', min_num_states=2,
                                                      splitter_category='BarcodeSequence')
        self.assertEqual(obs, exp)

        # Split by body site to filter samples/individuals out that would have
        # passed if the splitter category wasn't provided.
        exp = {'Palm': (set(['a', 'f']), 1, set(['3', '2'])),
               'Stool': (set(['c', 'g']), 1, set(['1', '2']))}
        obs = sample_ids_from_category_state_coverage(self.map_str2,
                                                      'Time', 'Individual', min_num_states=2,
                                                      splitter_category='BodySite')
        self.assertEqual(obs, exp)

    def test_sample_ids_from_category_state_coverage_invalid_input(self):
        """Test that errors are thrown on bad input."""
        # Using SampleID for coverage, subject, or splitter category.
        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'SampleID', 'DOB',
                          required_states=['Control', 'Fast'])

        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'SampleID',
                          required_states=['Control', 'Fast'])

        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'DOB',
                          required_states=['Control', 'Fast'], splitter_category='SampleID')

        # Nonexisting coverage, subject, or splitter category.
        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'foo', 'DOB',
                          required_states=['Control', 'Fast'])

        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'foo',
                          required_states=['Control', 'Fast'])

        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'DOB',
                          required_states=['Control', 'Fast'], splitter_category='foo')

        # Reusing same category for coverage, subject, or splitter category.
        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'Treatment',
                          required_states=['Control', 'Fast'])

        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'DOB',
                          required_states=['Control', 'Fast'], splitter_category='Treatment')

        # Nonexisting required coverage category state.
        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'DOB',
                          required_states=['Fast', 'foo'])

        # No filters are provided.
        self.assertRaises(ValueError, sample_ids_from_category_state_coverage,
                          self.tutorial_mapping_f, 'Treatment', 'DOB')

    def test_filter_otus_from_otu_map(self):
        """ filter_otus_from_otu_map functions as expected """
        otu_map_in = """o1 some comment	s1_1	s1_2
o2	s1_3	s1_4	s2_5
o3	s2_3
"""
        otu_map_no_single = """o1 some comment	s1_1	s1_2
o2	s1_3	s1_4	s2_5
"""
        otu_map_no_single_double = """o2	s1_3	s1_4	s2_5
"""
        otu_map_no_single_min_sample2 = """o2	s1_3	s1_4	s2_5
"""

        # write the test files
        fd, in_fp = mkstemp(dir=self.tmp_dir,
                        prefix='qiime_filter_test', suffix='.txt')
        close(fd)
        fasting_seqs_f = open(in_fp, 'w')
        fasting_seqs_f.write(otu_map_in)
        fasting_seqs_f.close()
        self.files_to_remove.append(in_fp)

        fd, actual_fp = mkstemp(dir=self.tmp_dir,
                            prefix='qiime_filter_test', suffix='.txt')
        close(fd)
        self.files_to_remove.append(actual_fp)

        retained_otus = filter_otus_from_otu_map(in_fp, actual_fp, 2)
        self.assertEqual(open(actual_fp).read(), otu_map_no_single)
        self.assertEqual(retained_otus, set(['o1 some comment', 'o2']))

        retained_otus = filter_otus_from_otu_map(in_fp, actual_fp, 3)
        self.assertEqual(open(actual_fp).read(), otu_map_no_single_double)
        self.assertEqual(retained_otus, set(['o2']))

        retained_otus = filter_otus_from_otu_map(in_fp, actual_fp, 2, 2)
        self.assertEqual(open(actual_fp).read(), otu_map_no_single_min_sample2)
        self.assertEqual(retained_otus, set(['o2']))


tree1 = "(aaa:10,(bbb:2,ccc:4):5);"
tree2 = "(aaa:10,('bbb':2,ccc:4):5);"

map_str1 = """#SampleID\tStudy\tBodySite\tDescription
a\tDog\tStool\tx
b\tDog\tStool\ty
c\tHand\tPalm\tz
d\tWholeBody\tPalm\ta
e\tWholeBody\tStool\tb"""

map_str2 = """#SampleID\tIndividual\tTime\tBodySite\tDescription
a\tI1\t2\tPalm\tx
b\tI2\t3\tStool\ty
c\tI1\t1\tStool\tz
d\tI3\t3\tStool\ta
e\tI3\t1\tPalm\tb
f\tI1\t3\tPalm\tc
g\tI1\t2\tStool\td"""

filter_fasta_expected1 = """>Seq1 some comment
ACCTTGG
>s2 some other comment
TTGG
>S3
AAGGCCGG
"""
filter_fasta_expected2 = """>S5 some comment
CGT
>seq6 some other comment
AA
>S7
T
"""

filter_fastq_expected1 = """@Seq1 some comment
ACCTTGG
+
BBBBBBB
@s2 some other comment
TTGG
+
BBBB
@S3
AAGGCCGG
+
BBCtatcc
"""
filter_fastq_expected2 = """@S5 some comment
CGT
+
BBB
@seq6 some other comment
AA
+
BB
@S7
T
+
s
"""

input_seqs_to_discard1 = """x
1 some comment
42 not a real otu id"""

input_dm1 = """\tABC\tDEF\tGHI\tXYZ
ABC\t0.0\t0.75\t0.00\t0.0063
DEF\t0.75\t0.0\t0.01\t0.65
GHI\t0.00\t0.01\t0.0\t1.0
XYZ\t0.0063\t0.065\t1.0\t0.0"""

expected_dm1a = """\tABC\tDEF
ABC\t0.0\t0.75
DEF\t0.75\t0.0"""

expected_dm1b = """\tABC\tXYZ
ABC\t0.0\t0.0063
XYZ\t0.0063\t0.0"""

tutorial_mapping_f = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636"""

expected_mapping_f1 = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636"""

sparse_otu_table1 = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "1", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "2", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "3", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "4", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "5", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "6", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "7", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "8", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "9", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "10", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "11", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "12", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "13", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "14", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "15", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "16", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "17", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "18", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "19", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "20", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "21", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "22", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "23", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "24", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "25", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "26", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "27", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "28", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "29", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "30", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "31", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "32", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "33", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "34", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "35", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "36", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "37", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "38", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "39", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "40", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "41", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "42", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "43", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "44", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "45", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Coprobacillus"]}}, {"id": "46", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "47", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "48", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "49", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "50", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "51", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "52", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "53", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Deltaproteobacteria"]}}, {"id": "54", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "55", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "56", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "57", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "58", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "59", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "60", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "61", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "62", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "63", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "64", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "65", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "66", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "67", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "68", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "69", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "70", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "71", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "72", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "73", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "74", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "75", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "76", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "77", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "78", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "79", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "80", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "81", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "82", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "83", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "84", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcus"]}}, {"id": "85", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "86", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "87", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "88", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "89", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "90", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Turicibacter"]}}, {"id": "91", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "92", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "93", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "94", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "95", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "96", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "97", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "98", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "99", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "100", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "101", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "102", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "103", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "104", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "105", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "106", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "107", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "108", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII", "Anaerovorax"]}}, {"id": "109", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "110", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Olsenella"]}}, {"id": "111", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "112", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "113", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "114", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "115", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "116", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "117", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "118", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "119", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "120", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "121", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "122", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "123", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "124", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "125", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "126", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "127", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "128", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "129", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "130", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Epsilonproteobacteria", "Campylobacterales", "Helicobacteraceae", "Helicobacter"]}}, {"id": "131", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "132", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "133", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "134", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "135", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "136", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "137", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "138", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "139", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "140", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "141", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "142", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "143", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "144", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "145", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "146", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "147", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "148", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "149", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "150", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "151", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "152", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "153", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "154", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "155", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "156", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "157", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "158", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "159", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "160", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "161", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "162", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "163", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "164", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "165", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "166", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "167", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "168", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "169", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "170", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "171", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "172", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "173", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "174", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Peptostreptococcaceae", "Peptostreptococcaceae Incertae Sedis"]}}, {"id": "175", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "176", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "177", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "178", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "179", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "180", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "181", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "182", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "183", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "184", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "185", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "186", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "187", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "188", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "189", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "190", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "191", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "192", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "193", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "194", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Acetanaerobacterium"]}}, {"id": "195", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "196", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "197", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "198", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "199", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "200", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "201", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "202", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "203", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "204", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "205", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "206", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "207", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "208", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "209", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "210", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "211", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "212", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "213", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "214", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "215", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "216", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "217", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "218", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "219", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "220", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "221", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "222", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "223", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "224", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "225", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "226", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "227", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "228", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "229", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII"]}}, {"id": "230", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "231", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "232", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "233", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "234", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "235", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "236", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "237", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "238", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "239", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "240", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "241", "metadata": {"taxonomy": ["Root", "Bacteria", "TM7", "TM7_genera_incertae_sedis"]}}, {"id": "242", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "243", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "244", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "245", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "246", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "247", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "248", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "249", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "250", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "251", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "252", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "253", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "254", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "255", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "256", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "257", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "258", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "259", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "260", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "261", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "262", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Bryantella"]}}, {"id": "263", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "264", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "265", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "266", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "267", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "268", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "269", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "270", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "271", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "272", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "273", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "274", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "275", "metadata": {"taxonomy": ["Root", "Bacteria", "Verrucomicrobia", "Verrucomicrobiae", "Verrucomicrobiales", "Verrucomicrobiaceae", "Akkermansia"]}}, {"id": "276", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "277", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "278", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "279", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "280", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "281", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "282", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "283", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "284", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "285", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "286", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "287", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "288", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "289", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "290", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "291", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "292", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "293", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "294", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "295", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "296", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "297", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "298", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "299", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "300", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "301", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "302", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "303", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "304", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "305", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "306", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "307", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "308", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcaceae Incertae Sedis"]}}, {"id": "309", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Denitrobacterium"]}}, {"id": "310", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "311", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "312", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "313", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "314", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "315", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "316", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "317", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "318", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria"]}}, {"id": "319", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "320", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "321", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "322", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "323", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "324", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "325", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "326", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "327", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "328", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "329", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "330", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "331", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "332", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "333", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "334", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "335", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "336", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "337", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "338", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "339", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "340", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "341", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "342", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "343", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "344", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "345", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "346", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "347", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "348", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "349", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "350", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "351", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "352", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "353", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "354", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "355", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "356", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "357", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "358", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "359", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "360", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "361", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "362", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "363", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae"]}}, {"id": "364", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "365", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "366", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Roseburia"]}}, {"id": "367", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "368", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "369", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "370", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "371", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "372", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "373", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Clostridiaceae", "Clostridiaceae 1", "Clostridium"]}}, {"id": "374", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "375", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "376", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "377", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "378", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "379", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "380", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "381", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "382", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "383", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "384", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "385", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Carnobacteriaceae", "Carnobacteriaceae 1"]}}, {"id": "386", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "387", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "388", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "389", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "390", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "391", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "392", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "393", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "394", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "395", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "396", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "397", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "398", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "399", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "400", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "401", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "402", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "403", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Prevotellaceae"]}}, {"id": "404", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "405", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "406", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "407", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "408", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "409", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "410", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "411", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "412", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "413", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "414", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "415", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "416", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 7, 1.0], [1, 5, 1.0], [2, 8, 1.0], [3, 0, 2.0], [3, 1, 1.0], [4, 0, 1.0], [5, 8, 1.0], [6, 7, 1.0], [7, 2, 2.0], [7, 8, 2.0], [8, 0, 1.0], [8, 1, 1.0], [8, 3, 2.0], [8, 4, 4.0], [9, 2, 2.0], [10, 1, 1.0], [11, 6, 1.0], [12, 6, 1.0], [13, 0, 1.0], [13, 3, 1.0], [13, 5, 1.0], [14, 2, 1.0], [14, 3, 1.0], [15, 4, 1.0], [16, 0, 1.0], [16, 2, 2.0], [17, 3, 1.0], [17, 6, 4.0], [17, 7, 10.0], [17, 8, 37.0], [18, 1, 1.0], [19, 8, 1.0], [20, 4, 1.0], [21, 6, 2.0], [21, 7, 3.0], [21, 8, 2.0], [22, 4, 2.0], [22, 6, 1.0], [23, 0, 14.0], [23, 1, 1.0], [23, 2, 14.0], [23, 3, 1.0], [24, 0, 1.0], [25, 3, 1.0], [26, 7, 1.0], [26, 8, 1.0], [27, 8, 1.0], [28, 1, 1.0], [29, 0, 6.0], [29, 2, 4.0], [29, 4, 2.0], [30, 5, 1.0], [31, 0, 1.0], [32, 4, 1.0], [33, 3, 1.0], [34, 6, 8.0], [34, 7, 10.0], [34, 8, 2.0], [35, 0, 1.0], [35, 2, 1.0], [36, 0, 1.0], [36, 2, 1.0], [36, 7, 1.0], [36, 8, 1.0], [37, 5, 1.0], [38, 2, 1.0], [39, 7, 1.0], [40, 2, 1.0], [41, 2, 1.0], [41, 7, 1.0], [42, 5, 1.0], [43, 5, 1.0], [44, 2, 1.0], [45, 0, 1.0], [46, 8, 1.0], [47, 3, 1.0], [48, 4, 1.0], [49, 3, 1.0], [50, 1, 1.0], [51, 1, 1.0], [52, 1, 2.0], [53, 6, 2.0], [53, 8, 1.0], [54, 6, 5.0], [55, 6, 1.0], [56, 5, 1.0], [57, 7, 1.0], [58, 0, 1.0], [58, 2, 1.0], [59, 8, 1.0], [60, 7, 1.0], [61, 2, 1.0], [61, 7, 1.0], [62, 2, 1.0], [63, 0, 1.0], [63, 2, 1.0], [64, 8, 1.0], [65, 3, 6.0], [65, 7, 1.0], [66, 2, 1.0], [67, 2, 1.0], [68, 0, 1.0], [69, 2, 1.0], [70, 5, 1.0], [71, 2, 1.0], [72, 5, 1.0], [73, 5, 5.0], [74, 3, 1.0], [75, 0, 1.0], [75, 2, 1.0], [76, 3, 1.0], [77, 3, 1.0], [78, 0, 1.0], [78, 2, 1.0], [78, 3, 1.0], [79, 0, 2.0], [79, 1, 3.0], [79, 2, 8.0], [79, 4, 1.0], [80, 8, 1.0], [81, 0, 1.0], [82, 5, 2.0], [83, 3, 1.0], [83, 7, 1.0], [84, 0, 1.0], [84, 7, 2.0], [85, 8, 1.0], [86, 7, 1.0], [87, 2, 1.0], [87, 5, 2.0], [87, 7, 1.0], [88, 7, 1.0], [89, 2, 1.0], [90, 3, 9.0], [90, 6, 3.0], [91, 3, 1.0], [92, 6, 1.0], [93, 6, 2.0], [93, 7, 1.0], [94, 7, 1.0], [95, 3, 2.0], [96, 3, 1.0], [96, 5, 1.0], [96, 7, 1.0], [96, 8, 1.0], [97, 5, 1.0], [98, 7, 1.0], [99, 3, 1.0], [100, 3, 1.0], [101, 3, 3.0], [102, 1, 1.0], [103, 1, 1.0], [103, 6, 1.0], [104, 5, 1.0], [105, 1, 1.0], [106, 5, 1.0], [107, 5, 1.0], [108, 6, 1.0], [109, 3, 1.0], [109, 6, 1.0], [109, 7, 5.0], [109, 8, 2.0], [110, 5, 2.0], [111, 6, 1.0], [112, 6, 1.0], [113, 5, 1.0], [114, 5, 1.0], [115, 5, 1.0], [116, 1, 1.0], [117, 0, 1.0], [117, 2, 2.0], [117, 5, 6.0], [118, 3, 1.0], [119, 7, 1.0], [120, 0, 1.0], [120, 1, 3.0], [120, 2, 1.0], [120, 3, 2.0], [120, 4, 1.0], [120, 5, 9.0], [120, 6, 2.0], [120, 7, 4.0], [120, 8, 5.0], [121, 7, 1.0], [122, 3, 1.0], [122, 5, 2.0], [123, 6, 1.0], [124, 6, 1.0], [125, 6, 1.0], [126, 2, 2.0], [126, 7, 1.0], [127, 5, 1.0], [128, 6, 1.0], [129, 3, 1.0], [130, 4, 5.0], [130, 5, 2.0], [131, 2, 1.0], [131, 3, 3.0], [132, 4, 1.0], [133, 2, 1.0], [134, 8, 1.0], [135, 2, 1.0], [136, 0, 1.0], [137, 7, 1.0], [138, 2, 1.0], [139, 0, 1.0], [140, 6, 1.0], [140, 7, 3.0], [141, 4, 1.0], [142, 4, 1.0], [143, 2, 1.0], [144, 5, 1.0], [145, 2, 2.0], [146, 0, 1.0], [146, 4, 2.0], [146, 6, 2.0], [146, 8, 3.0], [147, 1, 1.0], [147, 3, 1.0], [147, 4, 1.0], [147, 8, 3.0], [148, 5, 1.0], [149, 7, 1.0], [150, 4, 1.0], [151, 3, 1.0], [151, 7, 1.0], [152, 3, 1.0], [152, 6, 1.0], [152, 7, 2.0], [152, 8, 19.0], [153, 1, 2.0], [153, 2, 1.0], [153, 3, 2.0], [153, 6, 1.0], [153, 7, 1.0], [153, 8, 1.0], [154, 0, 2.0], [154, 1, 18.0], [154, 3, 1.0], [154, 6, 21.0], [154, 7, 4.0], [154, 8, 4.0], [155, 5, 5.0], [155, 6, 9.0], [155, 7, 5.0], [155, 8, 3.0], [156, 2, 1.0], [156, 7, 1.0], [157, 2, 1.0], [158, 0, 1.0], [158, 2, 1.0], [159, 7, 1.0], [159, 8, 1.0], [160, 6, 1.0], [161, 2, 1.0], [162, 5, 3.0], [162, 6, 5.0], [162, 7, 2.0], [162, 8, 6.0], [163, 8, 1.0], [164, 5, 1.0], [165, 0, 2.0], [165, 1, 1.0], [165, 2, 1.0], [166, 7, 1.0], [167, 0, 1.0], [168, 3, 1.0], [169, 1, 2.0], [169, 3, 7.0], [169, 7, 2.0], [170, 3, 1.0], [171, 3, 1.0], [172, 0, 1.0], [173, 5, 1.0], [174, 0, 1.0], [174, 4, 10.0], [175, 4, 1.0], [176, 5, 1.0], [177, 3, 1.0], [178, 3, 2.0], [179, 3, 1.0], [180, 4, 1.0], [181, 0, 1.0], [181, 1, 4.0], [181, 2, 2.0], [181, 3, 6.0], [182, 5, 1.0], [183, 6, 1.0], [184, 3, 1.0], [184, 6, 3.0], [184, 7, 1.0], [185, 8, 1.0], [186, 2, 1.0], [187, 1, 1.0], [187, 8, 1.0], [188, 7, 1.0], [189, 3, 1.0], [190, 7, 1.0], [191, 0, 2.0], [191, 1, 1.0], [191, 2, 10.0], [191, 3, 2.0], [191, 4, 24.0], [191, 7, 1.0], [191, 8, 1.0], [192, 5, 1.0], [193, 5, 1.0], [194, 2, 2.0], [195, 5, 1.0], [196, 5, 1.0], [196, 7, 1.0], [197, 1, 1.0], [198, 1, 2.0], [198, 5, 1.0], [199, 5, 1.0], [199, 6, 1.0], [200, 3, 2.0], [201, 3, 1.0], [201, 5, 1.0], [202, 6, 1.0], [203, 1, 2.0], [203, 2, 2.0], [203, 3, 4.0], [203, 5, 5.0], [203, 6, 1.0], [203, 7, 5.0], [204, 0, 1.0], [204, 1, 4.0], [204, 3, 1.0], [205, 7, 1.0], [206, 1, 1.0], [207, 7, 1.0], [208, 1, 2.0], [208, 3, 2.0], [208, 7, 1.0], [209, 2, 1.0], [210, 8, 1.0], [211, 0, 1.0], [211, 3, 1.0], [212, 8, 1.0], [213, 7, 2.0], [214, 7, 1.0], [215, 7, 1.0], [216, 7, 1.0], [217, 5, 2.0], [217, 7, 1.0], [218, 4, 9.0], [218, 5, 1.0], [219, 4, 1.0], [220, 0, 1.0], [220, 4, 1.0], [221, 7, 1.0], [222, 1, 1.0], [223, 7, 2.0], [223, 8, 2.0], [224, 3, 1.0], [225, 1, 2.0], [225, 2, 1.0], [226, 5, 1.0], [227, 1, 1.0], [227, 2, 2.0], [227, 4, 9.0], [227, 5, 1.0], [227, 6, 1.0], [227, 7, 1.0], [227, 8, 3.0], [228, 0, 16.0], [228, 4, 12.0], [229, 5, 1.0], [229, 6, 1.0], [230, 3, 1.0], [231, 1, 19.0], [231, 2, 2.0], [231, 4, 2.0], [231, 6, 3.0], [232, 6, 1.0], [233, 4, 1.0], [234, 4, 1.0], [235, 1, 1.0], [235, 2, 1.0], [235, 4, 1.0], [236, 5, 2.0], [237, 4, 1.0], [238, 7, 1.0], [239, 5, 1.0], [240, 5, 1.0], [241, 6, 2.0], [242, 6, 1.0], [243, 6, 1.0], [244, 8, 1.0], [245, 3, 1.0], [245, 7, 1.0], [246, 7, 1.0], [247, 2, 1.0], [248, 0, 1.0], [248, 3, 1.0], [249, 0, 1.0], [250, 0, 1.0], [250, 7, 1.0], [251, 3, 1.0], [251, 4, 4.0], [252, 3, 1.0], [253, 4, 2.0], [253, 7, 5.0], [254, 0, 11.0], [254, 1, 13.0], [254, 2, 6.0], [254, 3, 13.0], [254, 4, 2.0], [255, 5, 1.0], [256, 6, 1.0], [257, 6, 5.0], [258, 2, 1.0], [259, 7, 1.0], [260, 7, 1.0], [261, 7, 1.0], [262, 1, 1.0], [262, 8, 1.0], [263, 4, 1.0], [264, 5, 1.0], [265, 5, 2.0], [266, 3, 2.0], [267, 0, 1.0], [267, 3, 5.0], [267, 4, 17.0], [267, 5, 20.0], [268, 6, 1.0], [269, 3, 1.0], [270, 2, 1.0], [271, 8, 1.0], [272, 3, 1.0], [273, 6, 1.0], [274, 6, 1.0], [275, 6, 1.0], [276, 7, 1.0], [277, 0, 1.0], [278, 5, 1.0], [279, 5, 1.0], [280, 1, 1.0], [281, 0, 1.0], [282, 6, 2.0], [283, 6, 2.0], [283, 7, 1.0], [284, 3, 1.0], [285, 6, 1.0], [286, 1, 2.0], [286, 2, 3.0], [286, 3, 1.0], [286, 4, 4.0], [286, 6, 5.0], [286, 8, 4.0], [287, 6, 1.0], [287, 7, 1.0], [287, 8, 1.0], [288, 5, 1.0], [289, 4, 3.0], [290, 8, 2.0], [291, 4, 1.0], [292, 4, 1.0], [293, 5, 1.0], [294, 1, 1.0], [295, 0, 29.0], [295, 1, 1.0], [295, 2, 10.0], [296, 4, 1.0], [297, 3, 1.0], [298, 6, 1.0], [299, 6, 1.0], [299, 8, 1.0], [300, 5, 1.0], [301, 6, 2.0], [302, 5, 1.0], [303, 8, 1.0], [304, 7, 1.0], [305, 0, 1.0], [306, 8, 1.0], [307, 2, 1.0], [308, 1, 1.0], [309, 3, 1.0], [310, 2, 1.0], [311, 5, 1.0], [312, 2, 1.0], [313, 1, 1.0], [313, 8, 1.0], [314, 2, 1.0], [315, 0, 1.0], [315, 1, 3.0], [315, 2, 1.0], [316, 1, 1.0], [317, 6, 1.0], [318, 5, 1.0], [319, 1, 2.0], [319, 2, 1.0], [320, 3, 1.0], [321, 8, 1.0], [322, 3, 1.0], [323, 2, 1.0], [324, 2, 1.0], [325, 1, 1.0], [326, 4, 4.0], [326, 8, 2.0], [327, 7, 1.0], [328, 3, 1.0], [329, 0, 2.0], [329, 1, 2.0], [329, 3, 1.0], [330, 2, 1.0], [331, 4, 1.0], [332, 1, 1.0], [333, 5, 6.0], [333, 7, 3.0], [334, 0, 1.0], [335, 7, 1.0], [336, 2, 1.0], [337, 3, 1.0], [338, 7, 1.0], [339, 2, 1.0], [340, 2, 2.0], [341, 2, 1.0], [341, 7, 1.0], [342, 5, 1.0], [343, 8, 1.0], [344, 2, 1.0], [345, 0, 1.0], [346, 1, 1.0], [347, 3, 1.0], [348, 3, 1.0], [349, 6, 1.0], [349, 8, 1.0], [350, 0, 1.0], [351, 4, 2.0], [351, 5, 2.0], [351, 6, 1.0], [351, 7, 4.0], [351, 8, 1.0], [352, 0, 3.0], [353, 1, 4.0], [353, 2, 4.0], [353, 4, 1.0], [353, 5, 2.0], [353, 7, 2.0], [353, 8, 1.0], [354, 5, 1.0], [355, 7, 1.0], [356, 5, 1.0], [357, 3, 4.0], [358, 2, 1.0], [359, 2, 1.0], [360, 2, 1.0], [360, 7, 1.0], [360, 8, 1.0], [361, 0, 2.0], [361, 2, 2.0], [361, 3, 1.0], [362, 0, 1.0], [362, 3, 1.0], [363, 5, 1.0], [363, 7, 1.0], [364, 0, 1.0], [365, 5, 2.0], [366, 3, 1.0], [367, 4, 1.0], [368, 5, 1.0], [369, 5, 1.0], [370, 0, 2.0], [370, 1, 1.0], [370, 3, 5.0], [370, 5, 1.0], [371, 0, 1.0], [371, 1, 1.0], [372, 1, 1.0], [373, 1, 1.0], [373, 6, 3.0], [374, 6, 1.0], [375, 6, 4.0], [376, 8, 1.0], [377, 7, 1.0], [378, 8, 1.0], [379, 5, 1.0], [380, 8, 1.0], [381, 2, 2.0], [382, 7, 1.0], [383, 0, 4.0], [383, 1, 9.0], [383, 3, 2.0], [383, 7, 2.0], [384, 2, 1.0], [385, 8, 1.0], [386, 2, 1.0], [387, 2, 1.0], [388, 3, 1.0], [389, 1, 1.0], [390, 8, 1.0], [391, 8, 1.0], [392, 1, 1.0], [393, 5, 1.0], [394, 2, 1.0], [395, 0, 1.0], [395, 1, 1.0], [395, 2, 1.0], [396, 0, 2.0], [397, 7, 1.0], [398, 7, 1.0], [399, 6, 13.0], [400, 6, 1.0], [401, 1, 1.0], [402, 1, 1.0], [403, 7, 1.0], [404, 7, 1.0], [405, 7, 1.0], [406, 5, 1.0], [407, 0, 1.0], [407, 5, 4.0], [408, 0, 1.0], [408, 1, 5.0], [408, 2, 3.0], [408, 3, 2.0], [408, 8, 1.0], [409, 7, 1.0], [409, 8, 1.0], [410, 4, 1.0], [411, 3, 1.0], [412, 4, 2.0], [413, 7, 1.0], [414, 0, 1.0], [414, 2, 1.0], [415, 5, 7.0], [415, 7, 2.0], [415, 8, 2.0], [416, 1, 1.0], [416, 4, 1.0]], "columns": [{"id": "PC.354", "metadata": null}, {"id": "PC.355", "metadata": null}, {"id": "PC.356", "metadata": null}, {"id": "PC.481", "metadata": null}, {"id": "PC.593", "metadata": null}, {"id": "PC.607", "metadata": null}, {"id": "PC.634", "metadata": null}, {"id": "PC.635", "metadata": null}, {"id": "PC.636", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2544", "matrix_type": "sparse", "shape": [417, 9], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T12:53:24.597295", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

dense_otu_table1 = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "1", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "2", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "3", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "4", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "5", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "6", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "7", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "8", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "9", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "10", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "11", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "12", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "13", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "14", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "15", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "16", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "17", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "18", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "19", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "20", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "21", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "22", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "23", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "24", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "25", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "26", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "27", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "28", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "29", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "30", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "31", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "32", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "33", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "34", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "35", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "36", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "37", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "38", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "39", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "40", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "41", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "42", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "43", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "44", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "45", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Coprobacillus"]}}, {"id": "46", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "47", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "48", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "49", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "50", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "51", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "52", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "53", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Deltaproteobacteria"]}}, {"id": "54", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "55", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "56", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "57", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "58", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "59", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "60", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "61", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "62", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "63", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "64", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "65", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "66", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "67", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "68", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "69", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "70", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "71", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "72", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "73", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "74", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "75", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "76", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "77", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "78", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "79", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "80", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "81", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "82", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "83", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "84", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcus"]}}, {"id": "85", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "86", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "87", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "88", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "89", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "90", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Turicibacter"]}}, {"id": "91", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "92", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "93", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "94", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "95", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "96", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "97", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "98", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "99", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "100", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "101", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "102", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "103", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "104", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "105", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "106", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "107", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "108", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII", "Anaerovorax"]}}, {"id": "109", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "110", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Olsenella"]}}, {"id": "111", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "112", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "113", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "114", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "115", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "116", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "117", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "118", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "119", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "120", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "121", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "122", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "123", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "124", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "125", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "126", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "127", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "128", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "129", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "130", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria", "Epsilonproteobacteria", "Campylobacterales", "Helicobacteraceae", "Helicobacter"]}}, {"id": "131", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "132", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "133", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "134", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "135", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "136", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "137", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "138", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "139", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "140", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "141", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "142", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "143", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "144", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "145", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "146", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "147", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "148", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "149", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "150", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "151", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "152", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "153", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "154", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "155", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "156", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "157", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "158", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "159", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "160", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "161", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "162", "metadata": {"taxonomy": ["Root", "Bacteria", "Deferribacteres", "Deferribacteres", "Deferribacterales", "Deferribacteraceae", "Mucispirillum"]}}, {"id": "163", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "164", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "165", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "166", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "167", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "168", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "169", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "170", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "171", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "172", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "173", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "174", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Peptostreptococcaceae", "Peptostreptococcaceae Incertae Sedis"]}}, {"id": "175", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "176", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "177", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "178", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "179", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "180", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "181", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "182", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "183", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "184", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "185", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "186", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "187", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "188", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "189", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "190", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "191", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "192", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "193", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Butyrivibrio"]}}, {"id": "194", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Acetanaerobacterium"]}}, {"id": "195", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "196", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "197", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "198", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "199", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "200", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "201", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "202", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "203", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "204", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "205", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "206", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "207", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "208", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "209", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "210", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "211", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "212", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "213", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "214", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "215", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "216", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "217", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "218", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "219", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "220", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "221", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "222", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "223", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "224", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "225", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "226", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "227", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "228", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "229", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Incertae Sedis XIII"]}}, {"id": "230", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "231", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "232", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "233", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "234", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "235", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "236", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales"]}}, {"id": "237", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "238", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "239", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "240", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "241", "metadata": {"taxonomy": ["Root", "Bacteria", "TM7", "TM7_genera_incertae_sedis"]}}, {"id": "242", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "243", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "244", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "245", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "246", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "247", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "248", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"]}}, {"id": "249", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "250", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "251", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "252", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "253", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "254", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "255", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "256", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "257", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "258", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "259", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "260", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "261", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "262", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Bryantella"]}}, {"id": "263", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "264", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "265", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "266", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae", "Alistipes"]}}, {"id": "267", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "268", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "269", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "270", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "271", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "272", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "273", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "274", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "275", "metadata": {"taxonomy": ["Root", "Bacteria", "Verrucomicrobia", "Verrucomicrobiae", "Verrucomicrobiales", "Verrucomicrobiaceae", "Akkermansia"]}}, {"id": "276", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "277", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "278", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "279", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "280", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "281", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "282", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "283", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "284", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "285", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "286", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "287", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "288", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "289", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "290", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "291", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "292", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "293", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "294", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "295", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "296", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "297", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "298", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria"]}}, {"id": "299", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "300", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia"]}}, {"id": "301", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "302", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "303", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "304", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "305", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "306", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "307", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "308", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae", "Ruminococcaceae Incertae Sedis"]}}, {"id": "309", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae", "Denitrobacterium"]}}, {"id": "310", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "311", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "312", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "313", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Porphyromonadaceae", "Parabacteroides"]}}, {"id": "314", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "315", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "316", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "317", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "318", "metadata": {"taxonomy": ["Root", "Bacteria", "Proteobacteria"]}}, {"id": "319", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "320", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "321", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "322", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "323", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "324", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "325", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "326", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "327", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "328", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "329", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "330", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "331", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "332", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "333", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "334", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "335", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "336", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "337", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "338", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "339", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "340", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "341", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "342", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "343", "metadata": {"taxonomy": ["Root", "Bacteria", "Actinobacteria", "Actinobacteria", "Coriobacteridae", "Coriobacteriales", "Coriobacterineae", "Coriobacteriaceae"]}}, {"id": "344", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "345", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "346", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "347", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "348", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "349", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "350", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "351", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "352", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "353", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "354", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "355", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "356", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "357", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "358", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "359", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "360", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "361", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "362", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "363", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Rikenellaceae"]}}, {"id": "364", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "365", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "366", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Roseburia"]}}, {"id": "367", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "368", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "369", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "370", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "371", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "372", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "373", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Clostridiaceae", "Clostridiaceae 1", "Clostridium"]}}, {"id": "374", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "375", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Erysipelotrichi", "Erysipelotrichales", "Erysipelotrichaceae", "Erysipelotrichaceae Incertae Sedis"]}}, {"id": "376", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "377", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "378", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "379", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Ruminococcaceae"]}}, {"id": "380", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"]}}, {"id": "381", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "382", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "383", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "384", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "385", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Carnobacteriaceae", "Carnobacteriaceae 1"]}}, {"id": "386", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "387", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "388", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "389", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "390", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "391", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes"]}}, {"id": "392", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "393", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "394", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "395", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "396", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "397", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "398", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "399", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Bacteroidaceae", "Bacteroides"]}}, {"id": "400", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "401", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "402", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "403", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes", "Bacteroidetes", "Bacteroidales", "Prevotellaceae"]}}, {"id": "404", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae", "Lachnospiraceae Incertae Sedis"]}}, {"id": "405", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "406", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "407", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "408", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "409", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "410", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "411", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "412", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "413", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}, {"id": "414", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales", "Lachnospiraceae"]}}, {"id": "415", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "416", "metadata": {"taxonomy": ["Root", "Bacteria", "Firmicutes", "Clostridia", "Clostridiales"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 2, 0, 0, 0, 0, 0, 2], [1, 1, 0, 2, 4, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [1, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 1, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 4, 10, 37], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 3, 2], [0, 0, 0, 0, 2, 0, 1, 0, 0], [14, 1, 14, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [6, 0, 4, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 8, 10, 2], [1, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 1], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 6, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 5, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 1, 0, 0, 0, 0, 0], [2, 3, 8, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 9, 0, 0, 3, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 3, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 1, 5, 2], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 2, 0, 0, 6, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 3, 1, 2, 1, 9, 2, 4, 5], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 2, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 2, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 5, 2, 0, 0, 0], [0, 0, 1, 3, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 3, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 2, 0, 2, 0, 3], [0, 1, 0, 1, 1, 0, 0, 0, 3], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 1, 2, 19], [0, 2, 1, 2, 0, 0, 1, 1, 1], [2, 18, 0, 1, 0, 0, 21, 4, 4], [0, 0, 0, 0, 0, 5, 9, 5, 3], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 3, 5, 2, 6], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 0, 7, 0, 0, 0, 2, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 10, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 4, 2, 6, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 3, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [2, 1, 10, 2, 24, 0, 0, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 2, 4, 0, 5, 1, 5, 0], [1, 4, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 2, 0, 2, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 2, 0, 1, 0], [0, 0, 0, 0, 9, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 2, 2], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 2, 0, 9, 1, 1, 1, 3], [16, 0, 0, 0, 12, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 19, 2, 0, 2, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 1, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 4, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 5, 0], [11, 13, 6, 13, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 5, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 2, 0, 0, 0, 0, 0], [1, 0, 0, 5, 17, 20, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 2, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 2, 3, 1, 4, 0, 5, 0, 4], [0, 0, 0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 3, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [29, 1, 10, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 3, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 2, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 4, 0, 0, 0, 2], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [2, 2, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 6, 0, 3, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 2, 1, 4, 1], [3, 0, 0, 0, 0, 0, 0, 0, 0], [0, 4, 4, 0, 1, 2, 0, 2, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 4, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 1], [2, 0, 2, 1, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [2, 1, 0, 5, 0, 1, 0, 0, 0], [1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 3, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 4, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [4, 9, 0, 2, 0, 0, 0, 2, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0, 0, 0], [2, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 13, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 4, 0, 0, 0], [1, 5, 3, 2, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0], [1, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 7, 0, 2, 2], [0, 1, 0, 0, 1, 0, 0, 0, 0]], "columns": [{"id": "PC.354", "metadata": null}, {"id": "PC.355", "metadata": null}, {"id": "PC.356", "metadata": null}, {"id": "PC.481", "metadata": null}, {"id": "PC.593", "metadata": null}, {"id": "PC.607", "metadata": null}, {"id": "PC.634", "metadata": null}, {"id": "PC.635", "metadata": null}, {"id": "PC.636", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2544", "matrix_type": "dense", "shape": [417, 9], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T12:53:14.862973", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""

sparse_otu_table2 = """{"rows": [{"id": "o1", "metadata": {"taxonomy": "a;b;c"}}, {"id": "o2", "metadata": {"taxonomy": "a;b;d"}}, {"id": "o3", "metadata": {"taxonomy": "gg;b;e"}}, {"id": "o4", "metadata": {"taxonomy": "a;b;f"}}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 0, 66.0], [0, 1, 65.0], [0, 2, 2.0], [0, 3, 3.0], [0, 4, 6.0], [1, 0, 5.0], [1, 1, 6.0], [1, 3, 12.0], [1, 4, 3.0], [2, 0, 1.0], [2, 3, 85.0], [2, 4, 91.0], [3, 0, 28.0], [3, 1, 29.0], [3, 2, 98.0]], "columns": [{"id": "s1", "metadata": null}, {"id": "s2", "metadata": null}, {"id": "s3", "metadata": null}, {"id": "c1", "metadata": null}, {"id": "c2", "metadata": null}], "generated_by": "BIOM-Format 1.0.0c", "matrix_type": "sparse", "shape": [4, 5], "format_url": "http://biom-format.org", "date": "2012-12-03T12:42:38.166092", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
