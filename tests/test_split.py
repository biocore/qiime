#!/usr/bin/env python
# File created on 24 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import close
from tempfile import mkstemp
from unittest import TestCase, main

from biom.parse import parse_biom_table
from biom.table import Table
from skbio.sequence import DNA
from skbio.alignment import SequenceCollection
from skbio.parse.sequences import parse_fasta

from qiime.split import split_fasta
from qiime.util import get_qiime_temp_dir, remove_files

from itertools import product
from numpy import array, arange
from qiime.parse import parse_mapping_file
from qiime.split import (make_field_value_list, make_field_set_iterable,
                         make_non_empty_sample_lists, subset_mapping_data)

from numpy.testing import assert_array_equal


# create test table without metadata
_feature_ids = ['f%s' % i for i in range(10)]
_sample_ids = ['s%s' % i for i in range(5)]
_feature_vals = arange(50).reshape(10, 5)
TEST_TABLE1 = Table(_feature_vals, _feature_ids, _sample_ids)

# create test table with metadata 
_feature_ids = ['f%s' % i for i in range(3)]
_sample_ids = ['s%s' % i for i in range(5)]
_feature_vals = arange(15).reshape(3, 5)
observ_metadata = [{'taxonomy': ['Bacteria', 'A']},
                   {'taxonomy': ['Bacteria', 'B']},
                   {'taxonomy': ['Bacteria', 'C']}]
TEST_TABLE2 = Table(_feature_vals, _feature_ids, _sample_ids, observ_metadata)

# create a test mapping file
TEST_MF = '''#SampleID	color	temp	size
#USELESS COMMENTS
s0	blue	hot	13
s1	blue	cold	1
s2	green	cold	12
s3	cyan	hot	1
s4	blue	0	0'''


class SplitTests(TestCase):
    """Test biom table splitting functions."""

    def setUp(self):
        """Load data created on the fly with the biom.table.Table."""
        self.bt1 = TEST_TABLE1
        self.bt2 = TEST_TABLE2
        mdata, mheaders, _ = parse_mapping_file(TEST_MF.split('\n'))
        self.mdata = array(mdata)
        self.mheaders = mheaders

    def test_make_field_value_list(self):
        """Test that field values are correctly returned."""
        field = 'color'
        obs = make_field_value_list(self.mheaders, field, self.mdata)
        exp = ['blue', 'cyan', 'green']
        self.assertEqual(obs, exp)
        field = 'temp'
        obs = make_field_value_list(self.mheaders, field, self.mdata)
        exp = ['0', 'cold', 'hot']
        self.assertEqual(obs, exp)
        field = 'size'
        obs = make_field_value_list(self.mheaders, field, self.mdata)
        exp = ['0', '1', '12', '13']
        self.assertEqual(obs, exp)
        field = 'SampleID'
        obs = make_field_value_list(self.mheaders, field, self.mdata)
        exp = ['s0', 's1', 's2', 's3', 's4']
        self.assertEqual(obs, exp)

    def test_make_field_set_iterable(self):
        """Test that iteration order for field_setiterable is correct."""
        fields = ['color', 'temp']
        obs = list(make_field_set_iterable(self.mheaders, fields, self.mdata))
        exp = [('blue', '0'), ('blue', 'cold'), ('blue', 'hot'),
               ('cyan', '0'), ('cyan', 'cold'), ('cyan', 'hot'),
               ('green', '0'), ('green', 'cold'), ('green', 'hot')]
        self.assertEqual(obs, exp)
        fields = ['temp', 'color']
        obs = list(make_field_set_iterable(self.mheaders, fields, self.mdata))
        exp = [('0', 'blue'), ('0', 'cyan'), ('0', 'green'),
               ('cold', 'blue'), ('cold', 'cyan'), ('cold', 'green'),
               ('hot', 'blue'), ('hot', 'cyan'), ('hot', 'green')]
        self.assertEqual(obs, exp)
        fields = ['color', 'temp', 'SampleID']
        obs = list(make_field_set_iterable(self.mheaders, fields, self.mdata))
        exp = list(product(['blue', 'cyan', 'green'], ['0', 'cold', 'hot'],
                           ['s0', 's1', 's2', 's3', 's4']))
        self.assertEqual(obs, exp)

    def test_make_non_empty_sample_lists(self):
        """Test that sample lists are created correctly."""
        fields = ['color', 'temp']
        obs_sgs, obs_vgs = make_non_empty_sample_lists(fields, self.mheaders,
                                                       self.mdata)
        exp_sgs = [array(['s4'], dtype='|S5'),
                   array(['s1'], dtype='|S5'),
                   array(['s0'], dtype='|S5'),
                   array(['s3'], dtype='|S5'),
                   array(['s2'], dtype='|S5')]
        exp_vgs = [('blue', '0'),
                   ('blue', 'cold'),
                   ('blue', 'hot'),
                   ('cyan', 'hot'),
                   ('green', 'cold')]
        for i,j in zip(obs_sgs, exp_sgs):
            assert_array_equal(i, j)
        for i,j in zip(obs_vgs, exp_vgs):
            assert_array_equal(i, j)
        fields = ['color']
        obs_sgs, obs_vgs = make_non_empty_sample_lists(fields, self.mheaders,
                                                       self.mdata)
        exp_sgs = [array(['s0', 's1', 's4'], dtype='|S5'),
                   array(['s3'], dtype='|S5'),
                   array(['s2'], dtype='|S5')]
        exp_vgs = [('blue',), ('cyan',), ('green',)]
        for i,j in zip(obs_sgs, exp_sgs):
            assert_array_equal(i, j)
        for i,j in zip(obs_vgs, exp_vgs):
            assert_array_equal(i, j)

    def test_subset_mapping_data(self):
        """Test that mapping data is subset correctly."""
        samples_of_interest = ['s0', 's1', 's2']
        obs = subset_mapping_data(self.mdata, samples_of_interest)
        exp = array([['s0', 'blue', 'hot', '13'],
                     ['s1', 'blue', 'cold', '1'],
                     ['s2', 'green', 'cold', '12']],
                    dtype='|S5')
        assert_array_equal(obs, exp)
        samples_of_interest = ['s1', 's2', 's0'] #order change won't affect out
        obs = subset_mapping_data(self.mdata, samples_of_interest)
        assert_array_equal(obs, exp)
        samples_of_interest = ['s4', 's0'] #order change won't affect out
        obs = subset_mapping_data(self.mdata, samples_of_interest)
        exp = array([['s0', 'blue', 'hot', '13'],
                     ['s4', 'blue', '0', '0']],
                    dtype='|S5')
        assert_array_equal(obs, exp)


class SplitTestsFasta(TestCase):

    """ Tests of the fasta splitting capabilities of the split module """


    def test_split_fasta_equal_num_seqs_per_file(self):
        """split_fasta funcs as expected when equal num seqs go to each file
        """
        fd, filename_prefix = mkstemp(dir=get_qiime_temp_dir(),
                                     prefix='split_fasta_tests',
                                     suffix='')
        close(fd)
        infile = ['>seq1', 'AACCTTAA', '>seq2', 'TTAACC', 'AATTAA',
                  '>seq3', 'CCTT--AA']

        actual = split_fasta(infile, 1, filename_prefix)
        actual_seqs = []
        for fp in actual:
            actual_seqs += list(open(fp))
        remove_files(actual)

        expected = ['%s.%d.fasta' % (filename_prefix, i) for i in range(3)]

        self.assertEqual(actual, expected)
        self.assertEqual(
            SequenceCollection.from_fasta_records(parse_fasta(infile), DNA),
            SequenceCollection.from_fasta_records(parse_fasta(actual_seqs), DNA))

    def test_split_fasta_diff_num_seqs_per_file(self):
        """split_fasta funcs as expected when diff num seqs go to each file
        """
        fd, filename_prefix = mkstemp(dir=get_qiime_temp_dir(),
                                     prefix='split_fasta_tests',
                                     suffix='')
        close(fd)
        infile = ['>seq1', 'AACCTTAA', '>seq2', 'TTAACC', 'AATTAA',
                  '>seq3', 'CCTT--AA']

        actual = split_fasta(infile, 2, filename_prefix)

        actual_seqs = []
        for fp in actual:
            actual_seqs += list(open(fp))
        remove_files(actual)

        expected = ['%s.%d.fasta' % (filename_prefix, i) for i in range(2)]
        # list of file paths is as expected
        self.assertEqual(actual, expected)
        # building seq collections from infile and the split files result in
        # equivalent seq collections
        self.assertEqual(
            SequenceCollection.from_fasta_records(parse_fasta(infile), DNA),
            SequenceCollection.from_fasta_records(parse_fasta(actual_seqs), DNA))

    def test_split_fasta_diff_num_seqs_per_file_alt(self):
        """split_fasta funcs always catches all seqs
        """
        # start with 59 seqs (b/c it's prime, so should make more
        # confusing splits)
        in_seqs = SequenceCollection.from_fasta_records(
            [('seq%s' % k, 'AACCTTAA') for k in range(59)], DNA)
        infile = in_seqs.to_fasta().split('\n')

        # test seqs_per_file from 1 to 1000
        for i in range(1, 1000):
            fd, filename_prefix = mkstemp(dir=get_qiime_temp_dir(),
                                         prefix='split_fasta_tests',
                                         suffix='')
            close(fd)

            actual = split_fasta(infile, i, filename_prefix)

            actual_seqs = []
            for fp in actual:
                actual_seqs += list(open(fp))
            # remove the files now, so if the test fails they still get
            # cleaned up
            remove_files(actual)

            # building seq collections from infile and the split files result in
            # equivalent seq collections
            self.assertEqual(
                SequenceCollection.from_fasta_records(parse_fasta(infile), DNA),
                SequenceCollection.from_fasta_records(parse_fasta(actual_seqs), DNA))


if __name__ == "__main__":
    main()
