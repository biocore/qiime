#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Release"

"""Test suite for the compare_distance_matrices.py module."""

from string import digits
from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_distmat
from qiime.compare_distance_matrices import (run_mantel_correlogram,
                                             run_mantel_test)

class CompareDistanceMatricesTests(TestCase):
    """Tests for the compare_distance_matrices.py module.
    
    For these tests, we are interested in the structure of the strings that are
    constructed by the functions in compare_distance_matrices.py.  Thus, we
    simply remove all numbers from the strings so that we don't have to parse
    the strings for numbers, convert them, and then test them. The code that
    actually generates these numbers is already thoroughly tested in the stats
    module.
    """

    def remove_nums(self, text):
        """Removes all digits from the given string.

        Returns the string will all digits removed. Useful for testing strings
        for equality in unit tests where you don't care about numeric values,
        or if some values are random.

        This code was taken from http://bytes.com/topic/python/answers/
            850562-finding-all-numbers-string-replacing

        Arguments:
            text - the string to remove digits from
        """
        return text.translate(None, digits)
    

    def setUp(self):
        """Define some distance matrices that will be used by the tests."""
        self.dm1_str = ["\ts1\ts2\ts3", "s1\t0\t0.5\t0.2", "s2\t0.5\t0\t0.3",
                        "s3\t0.2\t0.3\t0"]
        self.dm1 = parse_distmat(self.dm1_str)
        self.dm2_str = ["\ts1\ts2\ts3", "s1\t0\t0.8\t0.25", "s2\t0.8\t0\t0.4",
                        "s3\t0.25\t0.4\t0"]
        self.dm2 = parse_distmat(self.dm2_str)
        self.dm3_str = ["\ts1\ts2\ts3", "s1\t0\t0.1\t0.2", "s2\t0.1\t0\t0.9",
                        "s3\t0.2\t0.9\t0"]
        self.dm3 = parse_distmat(self.dm3_str)
        self.dm4_str = ["\tz1\tz2\tz3", "z1\t0\t0.1\t0.2", "z2\t0.1\t0\t0.9",
                        "z3\t0.2\t0.9\t0"]
        self.dm4 = parse_distmat(self.dm4_str)
        self.distmats = [self.dm1, self.dm2, self.dm3]

        # Sample filepaths (these aren't created or modified, just used as
        # strings to be added to the results).
        self.fp1 = "foo.txt"
        self.fp2 = "bar.txt"
        self.fp3 = "baz.txt"
        self.fps = [self.fp1, self.fp2, self.fp3]

        # Some sample parameters to use for many of the tests.
        self.num_perms = 999
        self.comment = "# A sample comment.\n"
        self.alpha = 0.01
        self.tail_type = 'greater'
        self.sample_id_map = {'z1':'s1', 'z2':'s2', 'z3':'s3', 's1':'s1',
                              's2':'s2', 's3':'s3'}

    def test_run_mantel_test(self):
        """Test running mantel test on two distmats."""
        exp = '# A sample comment.\nDM\tDM\tNumber of entries\tMantel r ' + \
              'statistic\tp-value\tNumber of permutations\tTail type\n' + \
              'foo.txt\tbar.txt\t\t.\t.\t\tgreater\n'
        obs = run_mantel_test('mantel', [self.fp1, self.fp2],
                [self.dm1, self.dm2], self.num_perms, self.tail_type,
                self.comment, self.alpha)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_multiple(self):
        """Test running mantel test on three distmats."""
        exp = '# A sample comment.\nDM\tDM\tNumber of entries\tMantel r ' + \
              'statistic\tp-value\tNumber of permutations\tTail type\n' + \
              'foo.txt\tbar.txt\t\t.\t.\t\tgreater\nfoo.txt\tbaz.txt\t\t-.' + \
              '\t.\t\tgreater\nbar.txt\tbaz.txt\t\t-.\t.\t\tgreater\n'
        obs = run_mantel_test('mantel', [self.fp1, self.fp2, self.fp3],
                [self.dm1, self.dm2, self.dm3], self.num_perms, self.tail_type,
                self.comment, self.alpha)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_sample_id_map(self):
        """Test running mantel test on two distmats that need IDs mapped."""
        exp = '# A sample comment.\nDM\tDM\tNumber of entries\tMantel r ' + \
              'statistic\tp-value\tNumber of permutations\tTail type\n' + \
              'foo.txt\tbar.txt\t\t.\t.\t\tless\n'
        obs = run_mantel_test('mantel', [self.fp1, self.fp2],
                [self.dm3, self.dm4], self.num_perms, 'less',
                self.comment, 0.5, sample_id_map=self.sample_id_map)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_too_small(self):
        """Test running mantel test on two distmats that are incompatible."""
        exp = '# A sample comment.\nDM\tDM\tNumber of entries\tMantel r ' + \
        'statistic\tp-value\tNumber of permutations\tTail type\nfoo.txt\t' + \
        'bar.txt\t\tToo few samples\n'
        obs = run_mantel_test('mantel', [self.fp1, self.fp2],
                [self.dm3, self.dm4], self.num_perms, 'less',
                self.comment, 0.5)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_single_matrix(self):
        """Test running mantel test on a single dm."""
        exp = '# A sample comment.\nDM\tDM\tNumber of entries\tMantel r ' + \
              'statistic\tp-value\tNumber of permutations\tTail type\n'
        obs = run_mantel_test('mantel', [self.fp1], [self.dm3], self.num_perms,
                'less', self.comment, 0.5)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_invaid_input(self):
        """Test running mantel test on invalid input."""
        self.assertRaises(ValueError, run_mantel_test, 'mantel',
                [self.fp1, self.fp2], [self.dm1], self.num_perms, 'two sided',
                None, self.alpha)

    def test_run_mantel_test_bad_method(self):
        """Test running mantel test using a bad method."""
        self.assertRaises(ValueError, run_mantel_test, 'foo',
                [self.fp1, self.fp2], [self.dm1, self.dm2], self.num_perms,
                'two sided', None, self.alpha)

    def test_run_mantel_test_no_cdm(self):
        """Test running partial mantel test with no control dm."""
        self.assertRaises(ValueError, run_mantel_test, 'partial_mantel',
                [self.fp1, self.fp2], [self.dm1, self.dm2], self.num_perms,
                'two sided', self.comment, self.alpha)

    def test_run_mantel_test_partial_mantel(self):
        """Test running partial mantel test with two dms and a control dm."""
        exp = '# A sample comment.\nDM\tDM\tCDM\tNumber of entries\t' + \
              'Mantel r statistic\tp-value\tNumber of permutations\tTail ' + \
              'type\nfoo.txt\tbar.txt\tbaz.txt\t\t.\t.\t\tgreater\n'
        obs = run_mantel_test('partial_mantel', [self.fp1, self.fp2],
                [self.dm1, self.dm2], self.num_perms, self.tail_type,
                self.comment, self.fp3, self.dm3)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_partial_mantel_too_small(self):
        """Test running partial mantel test with incompatible dms."""
        exp = '# A sample comment.\nDM\tDM\tCDM\tNumber of entries\t' + \
              'Mantel r statistic\tp-value\tNumber of permutations\t' + \
              'Tail type\nfoo.txt\tbar.txt\tbaz.txt\t\tToo few samples\n'
        obs = run_mantel_test('partial_mantel', [self.fp1, self.fp2],
                [self.dm1, self.dm2], self.num_perms, self.tail_type,
                self.comment, self.fp3, self.dm4)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_test_partial_mantel_sample_id_map(self):
        """Test running partial mantel test with incompatible dms to map."""
        exp = '# A sample comment.\nDM\tDM\tCDM\tNumber of entries\tMantel' + \
              ' r statistic\tp-value\tNumber of permutations\tTail type\n' + \
              'foo.txt\tbar.txt\tbaz.txt\t\t.\t.\t\tgreater\n'
        obs = run_mantel_test('partial_mantel', [self.fp1, self.fp2],
                [self.dm1, self.dm2], self.num_perms, self.tail_type,
                self.comment, self.fp3, self.dm4, self.sample_id_map)
        self.assertEqual(self.remove_nums(obs), exp)

    def test_run_mantel_correlogram(self):
        """Test running mantel correlogram on two distmats."""
        exp = ('# A sample comment.\nDM\tDM\tNumber of entries\tNumber of '
        'permutations\tClass index\tNumber of distances\tMantel r statistic\t'
        'p-value\tp-value (Bonferroni corrected)\tTail type\nfoo.txt\tbar.txt'
        '\t\t\t.\t\t.\t.\t.\tgreater\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n\t'
        '\t\t\t.\t\tNone\tNone\tNone\tNone\n',
        ['foo.txt_AND_bar.txt_mantel_correlogram.'], 1)
        obs = run_mantel_correlogram([self.fp1, self.fp2],
                [self.dm1, self.dm2], self.num_perms, self.comment, self.alpha)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_multiple(self):
        """Test running mantel correlogram on three distmats."""
        exp = ('# A sample comment.\nDM\tDM\tNumber of entries\tNumber of '
        'permutations\tClass index\tNumber of distances\tMantel r statistic\t'
        'p-value\tp-value (Bonferroni corrected)\tTail type\nfoo.txt\tbar.txt'
        '\t\t\t.\t\t.\t.\t.\tgreater\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n\t'
        '\t\t\t.\t\tNone\tNone\tNone\tNone\nfoo.txt\tbaz.txt\t\t\t.\t\t-.\t.\t'
        '.\tless\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n\t\t\t\t.\t\tNone\tNone'
        '\tNone\tNone\nbar.txt\tbaz.txt\t\t\t.\t\t-.\t.\t.\tless\n\t\t\t\t.\t'
        '\tNone\tNone\tNone\tNone\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n',
        ['foo.txt_AND_bar.txt_mantel_correlogram.',
        'foo.txt_AND_baz.txt_mantel_correlogram.',
        'bar.txt_AND_baz.txt_mantel_correlogram.'], 3)
        obs = run_mantel_correlogram(self.fps, self.distmats, self.num_perms,
                                     self.comment, self.alpha)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_sample_id_map(self):
        """Test running mantel correlogram on two dms that need IDs mapped."""
        exp = ('# A sample comment.\nDM\tDM\tNumber of entries\tNumber of '
               'permutations\tClass index\tNumber of distances\tMantel r '
               'statistic\tp-value\tp-value (Bonferroni corrected)\tTail type'
               '\nfoo.txt\tbar.txt\t\t\t.\t\t.\t.\t.\tgreater\n\t\t\t\t.\t\t'
               'None\tNone\tNone\tNone\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n',
               ['foo.txt_AND_bar.txt_mantel_correlogram.'], 1)
        obs = run_mantel_correlogram([self.fp1, self.fp2],
                [self.dm3, self.dm4], self.num_perms, self.comment, self.alpha,
                sample_id_map=self.sample_id_map)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_too_small(self):
        """Test running mantel correlogram on two incompatible dms."""
        exp = ('# A sample comment.\nDM\tDM\tNumber of entries\tNumber of '
        'permutations\tClass index\tNumber of distances\tMantel r statistic\t'
        'p-value\tp-value (Bonferroni corrected)\tTail type\nfoo.txt\tbaz.txt'
        '\t\tToo few samples\n', [], 0)
        obs = run_mantel_correlogram([self.fp1, self.fp3],
                [self.dm2, self.dm4], self.num_perms, self.comment, self.alpha)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_single_matrix(self):
        """Test running mantel correlogram on one dm."""
        exp = ('# A sample comment.\nDM\tDM\tNumber of entries\tNumber of '
        'permutations\tClass index\tNumber of distances\tMantel r statistic\t'
        'p-value\tp-value (Bonferroni corrected)\tTail type\n', [], 0)
        obs = run_mantel_correlogram([self.fp1], [self.dm1], self.num_perms,
                self.comment, self.alpha, sample_id_map=self.sample_id_map)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_no_comment(self):
        """Test running mantel correlogram without supplying a comment."""
        exp = ('DM\tDM\tNumber of entries\tNumber of permutations\tClass index'
        '\tNumber of distances\tMantel r statistic\tp-value\tp-value '
        '(Bonferroni corrected)\tTail type\nfoo.txt\tbar.txt\t\t\t.\t\t-.\t.'
        '\t.\tless\n\t\t\t\t.\t\tNone\tNone\tNone\tNone\n\t\t\t\t.\t\tNone\t'
        'None\tNone\tNone\n', ['foo.txt_AND_bar.txt_mantel_correlogram.'], 1)
        obs = run_mantel_correlogram([self.fp1, self.fp2],
                [self.dm1, self.dm3], self.num_perms, None, self.alpha,
                sample_id_map=self.sample_id_map)
        self.assertEqual((self.remove_nums(obs[0]), obs[1], len(obs[2])), exp)

    def test_run_mantel_correlogram_invaid_input(self):
        """Test running mantel correlogram on invalid input."""
        self.assertRaises(ValueError, run_mantel_correlogram, [self.fp1],
                [self.dm1, self.dm3], self.num_perms, None, self.alpha)


if __name__ == "__main__":
    main()
