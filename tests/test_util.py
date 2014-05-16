#!/usr/bin/env python
from __future__ import division
# unit tests for util.py

from os import chdir, getcwd, mkdir, rmdir, remove, close
from os.path import split, abspath, dirname, exists, isdir, join
from glob import glob
from random import seed
from shutil import rmtree
from StringIO import StringIO
from tempfile import mkdtemp, mkstemp
from collections import defaultdict
from unittest import TestCase, main
import gzip

from biom.parse import parse_biom_table_str, parse_biom_table

from numpy.testing import assert_almost_equal

from skbio.core.sequence import DNASequence
from skbio.parse.sequences import parse_fasta
from skbio.util.misc import remove_files

from cogent.cluster.procrustes import procrustes

from brokit.formatdb import build_blast_db_from_fasta_file

from qiime.parse import (fields_to_dict, parse_distmat, parse_mapping_file,
                         parse_mapping_file_to_dict, parse_otu_table,
                         QiimeParseError)
from qiime.util import (make_safe_f, FunctionWithParams, qiime_blast_seqs,
                        extract_seqs_by_sample_id, get_qiime_project_dir,
                        get_qiime_scripts_dir, matrix_stats,
                        raise_error_on_parallel_unavailable,
                        convert_OTU_table_relative_abundance, create_dir,
                        summarize_pcoas, _compute_jn_pcoa_avg_ranges, _flip_vectors, IQR,
                        idealfourths, isarray, matrix_IQR, degap_fasta_aln,
                        write_degapped_fasta_to_file, compare_otu_maps, get_diff_for_otu_maps,
                        convert_otu_table_relative, write_seqs_to_fasta,
                        split_fasta_on_sample_ids, split_fasta_on_sample_ids_to_dict,
                        split_fasta_on_sample_ids_to_files, median_absolute_deviation,
                        guess_even_sampling_depth, compute_days_since_epoch,
                        get_interesting_mapping_fields, inflate_denoiser_output,
                        flowgram_id_to_seq_id_map, count_seqs, count_seqs_from_file,
                        count_seqs_in_filepaths,
                        iseq_to_qseq_fields,
                        make_compatible_distance_matrices, stderr, _chk_asarray, expand_otu_ids,
                        subsample_fasta, summarize_otu_sizes_from_otu_map,
                        load_qiime_config, MetadataMap,
                        RExecutor, duplicates_indices, trim_fasta, get_qiime_temp_dir,
                        qiime_blastx_seqs, add_filename_suffix, is_valid_git_refname,
                        is_valid_git_sha1, sync_biom_and_mf,
                        biom_taxonomy_formatter, invert_dict)

import numpy
from numpy import array, asarray


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso",
               "Justin Kuczynski", "Jens Reeder", "Catherine Lozupone",
               "Jai Ram Rideout", "Logan Knecht", "Michael Dwan",
               "Levi McCracken", "Damien Coy", "Yoshiki Vazquez Baeza",
               "Will Van Treuren", "Luke Ursell"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


class TopLevelTests(TestCase):

    """Tests of top-level module functions."""

    def setUp(self):
        self.otu_table_f1 = otu_table_fake1.split('\n')
        self.otu_table_f2 = otu_table_fake2.split('\n')
        self.otu_table_f3 = otu_table_fake3.split('\n')
        self.otu_table_f4 = otu_table_fake4.split('\n')
        self.otu_table_f1_no_tax = otu_table_fake1_no_tax.split('\n')
        self.otu_table_f2_no_tax = otu_table_fake2_no_tax.split('\n')
        self.otu_table_f3_no_tax = otu_table_fake3_no_tax.split('\n')
        self.fasta1 = fasta1.split('\n')
        self.fasta2 = fasta2.split('\n')
        self.mapping_f1 = mapping_f1.split('\n')
        self.mapping_f2 = mapping_f2.split('\n')
        self.mapping_f3 = mapping_f3.split('\n')
        self.dirs_to_remove = []
        self.files_to_remove = []
        self.centroid_seqs1 = centroid_seqs1.split('\n')
        self.singleton_seqs1 = singleton_seqs1.split('\n')
        self.denoiser_mapping1 = denoiser_mapping1.split('\n')
        self.raw_seqs1 = raw_seqs1.split('\n')
        self.fastq_barcodes = fastq_barcodes
        self.fasta_barcodes = fasta_barcodes

    def tearDown(self):
        remove_files(self.files_to_remove)
        for dir in self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)

    def test_expand_otu_ids(self):
        """expand otu ids functions as expected """
        otu_map = {'o1': ['s1', 's2'],
                   'o2': ['s3'],
                   '3': ['s4', 's5']}
        otus_to_expand = ['3', 'o1 test']
        actual = expand_otu_ids(otu_map, otus_to_expand)
        expected = ['s4', 's5', 's1', 's2']
        self.assertEqual(actual, expected)

        # ignore missing
        otus_to_expand = ['3', 'o1 test', '99']
        actual = expand_otu_ids(otu_map,
                                otus_to_expand,
                                ignore_missing=True)
        expected = ['s4', 's5', 's1', 's2']
        self.assertEqual(actual, expected)

    def test_median_absolute_deviation(self):
        """ median_absolute_deviation returns MAD and median """
        data = [0, 0, 0, 0, 0, 0]
        expected = (0, 0)
        self.assertEqual(median_absolute_deviation(data), expected)
        data = [1, 1, 1, 1, 1, 1]
        expected = (0, 1)
        self.assertEqual(median_absolute_deviation(data), expected)
        data = [6, 1, 2, 9, 4, 1, 2]
        expected = (1, 2)
        self.assertEqual(median_absolute_deviation(data), expected)
        data = [-6, -1, -2, -9, -4, -1, -2]
        expected = (1, -2)
        self.assertEqual(median_absolute_deviation(data), expected)

    def test_guess_even_sampling_depth(self):
        """ guess_even_sampling_depth functions as expected """
        data = [6, 1, 2, 9, 4, 1, 2]
        expected = 1  # MAD = 2.25; med - MAD = -0.25
        self.assertEqual(guess_even_sampling_depth(data), expected)

    def test_write_seqs_to_fasta(self):
        """ write_seqs_to_fasta functions as expected """
        fd, output_fp = mkstemp(
            prefix="qiime_util_write_seqs_to_fasta_test",
            suffix='.fasta')
        close(fd)
        self.files_to_remove.append(output_fp)
        seqs = [('s1', 'ACCGGTTGG'), ('s2', 'CCTTGG'),
                ('S4 some comment string', 'A')]
        exp = ">s1\nACCGGTTGG\n>s2\nCCTTGG\n>S4 some comment string\nA\n"
        # works in write mode
        write_seqs_to_fasta(output_fp, seqs, 'w')
        self.assertEqual(open(output_fp).read(), exp)
        # calling again in write mode overwrites original file
        write_seqs_to_fasta(output_fp, seqs, 'w')
        self.assertEqual(open(output_fp).read(), exp)
        # works in append mode
        exp2 = exp + exp
        write_seqs_to_fasta(output_fp, seqs, 'a')
        self.assertEqual(open(output_fp).read(), exp2)

    def test_summarize_otu_sizes_from_otu_map(self):
        """ summarize_otu_sizes_from_otu_map functions as expected """
        otu_map_f = """O1	seq1
o2	seq2	seq3	seq4	seq5
o3	seq5
o4	seq6	seq7""".split('\n')
        expected = [(1, 2), (2, 1), (4, 1)]
        self.assertEqual(summarize_otu_sizes_from_otu_map(otu_map_f), expected)

    def test_split_fasta_on_sample_ids(self):
        """ split_fasta_on_sample_ids functions as expected
        """
        actual = list(split_fasta_on_sample_ids(
                      parse_fasta(self.fasta1)))
        expected = [('Samp1', 'Samp1_42', 'ACCGGTT'),
                    ('s2_a', 's2_a_50', 'GGGCCC'),
                    ('Samp1', 'Samp1_43 some comme_nt', 'AACCG'),
                    ('s3', 's3_25', 'AAACCC')]
        self.assertEqual(actual, expected)

    def test_split_fasta_on_sample_ids_to_dict(self):
        """ split_fasta_on_sample_ids_to_dict functions as expected
        """
        actual = split_fasta_on_sample_ids_to_dict(
            parse_fasta(self.fasta1))
        expected = {'Samp1': [('Samp1_42', 'ACCGGTT'),
                              ('Samp1_43 some comme_nt', 'AACCG')],
                    's2_a': [('s2_a_50', 'GGGCCC')],
                    's3': [('s3_25', 'AAACCC')]}
        self.assertEqual(actual, expected)

    def test_split_fasta_on_sample_ids_to_files(self):
        """ split_fasta_on_sample_ids_to_files functions as expected
        """
        temp_output_dir = mkdtemp()
        self.dirs_to_remove.append(temp_output_dir)

        split_fasta_on_sample_ids_to_files(
            parse_fasta(self.fasta2),
            output_dir=temp_output_dir,
            per_sample_buffer_size=2)
        self.files_to_remove.extend(glob('%s/*fasta' % temp_output_dir))

        # confirm that all files are as expected
        self.assertEqual(open('%s/Samp1.fasta' % temp_output_dir).read(),
                         ">Samp1_42\nACCGGTT\n>Samp1_43 some comme_nt\nAACCG\n>Samp1_44\nA\n")
        self.assertEqual(open('%s/s2_a.fasta' % temp_output_dir).read(),
                         ">s2_a_50\nGGGCCC\n")
        self.assertEqual(open('%s/s3.fasta' % temp_output_dir).read(),
                         ">s3_25\nAAACCC\n")
        # confirm number of files is as expected
        self.assertEqual(len(glob('%s/*' % temp_output_dir)), 3)

    def test_convert_otu_table_relative(self):
        """should convert a parsed otu table into relative abundances"""
        otu_table = parse_otu_table(self.otu_table_f1)
        exp_counts = array([[1.0 / 6, 0],
                            [1.0 / 6, 0],
                            [4.0 / 6, 0]])
        rel_otu_table = convert_otu_table_relative(otu_table)
        self.assertEqual(rel_otu_table[0], otu_table[0])
        self.assertEqual(rel_otu_table[1], otu_table[1])
        assert_almost_equal(rel_otu_table[2], exp_counts)
        self.assertEqual(rel_otu_table[3], otu_table[3])

    def test_make_safe_f(self):
        """make_safe_f should return version of f that ignores extra kwargs."""
        def f(x, y):
            return x * y
        self.assertEqual(f(3, 4), 12)
        g = make_safe_f(f, ['x', 'y'])
        self.assertEqual(g(3, 4), 12)
        self.assertEqual(g(x=3, y=4, z=10, xxx=11), 12)

    def test_extract_seqs_by_sample_id(self):
        """extract_seqs_by_sample_id: functions as expected """

        seqs = [('Samp1_109', 'ACGG'),
                ('Samp1_110', 'CCGG'),
                ('samp1_109', 'GCGG'),
                ('S2', 'AA'),
                ('S3', 'CC'),
                ('S4', 'GG'),
                ('S44', 'TT'),
                ('S4', 'TAAT')]
        sample_ids = ['Samp1', 'S44']
        expected = [('Samp1_109', 'ACGG'),
                    ('Samp1_110', 'CCGG'),
                    ('S44', 'TT')]
        actual = list(extract_seqs_by_sample_id(seqs, sample_ids))
        self.assertEqual(actual, expected)

        # negated
        expected_neg = [('samp1_109', 'GCGG'),
                        ('S2', 'AA'),
                        ('S3', 'CC'),
                        ('S4', 'GG'),
                        ('S4', 'TAAT')]
        actual = list(extract_seqs_by_sample_id(seqs, sample_ids, negate=True))
        self.assertEqual(actual, expected_neg)

        # OK if user passes dict of sample ids
        sample_ids = {'samp1': 25}
        expected = [('samp1_109', 'GCGG')]
        actual = list(extract_seqs_by_sample_id(seqs, sample_ids))
        self.assertEqual(actual, expected)

    def test_get_qiime_project_dir(self):
        """getting the qiime project directory functions as expected """

        # Do an explicit check on whether the file system containing
        # the current file is case insensitive. This is in response
        # to SF bug #2945548, where this test would fail on certain
        # unusual circumstances on case-insensitive file systems
        # because the case of abspath(__file__) was inconsistent.
        # (If you don't believe this, set case_insensitive_filesystem
        # to False, and rename your top-level Qiime directory as
        # qiime on OS X. That sould cause this test to fail as
        # actual will be path/to/qiime and expected will be
        # path/to/Qiime.) Note that we don't need to change anything
        # in the get_qiime_project_dir() function as if the
        # file system is case insenstive, the case of the returned
        # string is irrelevant.
        case_insensitive_filesystem = \
            exists(__file__.upper()) and exists(__file__.lower())

        actual = get_qiime_project_dir()
        # I base the expected here off the imported location of
        # qiime/util.py here, to handle cases where either the user has
        # Qiime in their PYTHONPATH, or when they've installed it with
        # setup.py.
        # If util.py moves this test will fail -- that
        # is what we want in this case, as the get_qiime_project_dir()
        # function would need to be modified.
        import qiime.util
        util_py_filepath = abspath(abspath(qiime.util.__file__))
        expected = dirname(dirname(util_py_filepath))

        if case_insensitive_filesystem:
            # make both lowercase if the file system is case insensitive
            actual = actual.lower()
            expected = expected.lower()
        self.assertEqual(actual, expected)

    def test_get_qiime_scripts_dir(self):
        """Test inferring the QIIME scripts directory."""
        obs = get_qiime_scripts_dir()

        # We can't do much testing of the observed value, but let's at least
        # check that the directory exists.
        self.assertTrue(isdir(obs))

    def test_matrix_stats1(self):
        """ matrix_stats should match mean, median, stdev calc'd by hand"""
        headers_list = [['a', 'c', 'b'], ['a', 'c', 'b']]
        d1 = numpy.array([[0, .2, .9],
                          [.2, 0, .8],
                          [.9, .8, 0]], 'float')
        d2 = numpy.array([[0, .3, 1.1],
                          [.3, 0, .8],
                          [1.1, .8, 0]], 'float')
        distmats_list = [d1, d2]

        exp_mean = numpy.array([[0, .25, 1.0],
                                [.25, 0, .8],
                                [1.0, .8, 0]], 'float')
        exp_median = numpy.array([[0, .25, 1.0],
                                  [.25, 0, .8],
                                  [1.0, .8, 0]], 'float')
        exp_std = numpy.array([[0, .05, .1],
                               [.05, 0, 0],
                               [.1, 0, 0]], 'float')
        results = matrix_stats(headers_list, distmats_list)
        assert_almost_equal(results[1:], [exp_mean, exp_median, exp_std])
        self.assertEqual(results[0], ['a', 'c', 'b'])

    def test_matrix_stats2(self):
        """ matrix_stats should raise valerr if headers are in different orders
        """
        headers_list = [['a', 'c', 'b'], ['b', 'c', 'a']]
        d1 = numpy.array([[0, .2, .9],
                          [.2, 0, .8],
                          [.9, .8, 0]], 'float')
        d2 = numpy.array([[0, .3, 1.1],
                          [.3, 0, .8],
                          [1.1, .8, 0]], 'float')
        distmats_list = [d1, d2]

        exp_mean = numpy.array([[0, .25, 1.0],
                                [.25, 0, .8],
                                [1.0, .8, 0]], 'float')
        exp_median = numpy.array([[0, .25, 1.0],
                                  [.25, 0, .8],
                                  [1.0, .8, 0]], 'float')
        exp_std = numpy.array([[0, .05, .1],
                               [.05, 0, 0],
                               [.1, 0, 0]], 'float')
        self.assertRaises(
            ValueError,
            matrix_stats,
            headers_list,
            distmats_list)

    def test_raise_error_on_parallel_unavailable(self):
        """raise_error_on_parallel_unavailable functions as expected """
        self.assertRaises(
            RuntimeError,
            raise_error_on_parallel_unavailable,
            {})
        self.assertRaises(RuntimeError, raise_error_on_parallel_unavailable,
                          {'jobs_to_start': '1'})
        raise_error_on_parallel_unavailable({'jobs_to_start': '2'})
        raise_error_on_parallel_unavailable({'jobs_to_start': '24'})

    def test_support_files_available(self):
        """support_files are available """
        # check that the qiime/support_files directory exists
        support_files_dir = \
            join(get_qiime_project_dir(), 'qiime', 'support_files')
        self.assertTrue(exists(support_files_dir))

        # check that a file in qiime/support_files exists
        default_qiime_config_fp = join(support_files_dir, 'qiime_config')
        self.assertTrue(exists(default_qiime_config_fp))

    def test_compute_days_since_epoch(self):
        """compute_days_since_epoch functions as expected """
        self.assertEqual(compute_days_since_epoch(29, 10, 2002), 11989)
        # can pass keyword arguments
        self.assertEqual(
            compute_days_since_epoch(
                year=2002,
                month=10,
                day=29),
            11989)
        self.assertEqual(
            compute_days_since_epoch(
                day=27,
                month=3,
                year=2009),
            14330)
        self.assertEqual(
            compute_days_since_epoch(
                day="27",
                month="3",
                year="2009"),
            14330)

    def test_get_interesting_mapping_fields(self):
        """get_interesting_mapping_fields returns expected fields """
        # all columns are completely unique
        d = parse_mapping_file(self.mapping_f1)
        actual = get_interesting_mapping_fields(d[0], d[1])
        expected = []
        self.assertEqual(actual, expected)

        # all columns are completely identical
        d = parse_mapping_file(self.mapping_f2)
        actual = get_interesting_mapping_fields(d[0], d[1])
        expected = []
        self.assertEqual(actual, expected)

        # some columns retained
        d = parse_mapping_file(self.mapping_f3)
        actual = get_interesting_mapping_fields(d[0], d[1])
        expected = ['Something', 'days_since_epoch']
        self.assertEqual(actual, expected)

    def test_inflate_denoiser_output(self):
        """ inflate_denoiser_output expands denoiser results as expected """
        actual = list(inflate_denoiser_output(
            parse_fasta(self.centroid_seqs1),
            parse_fasta(self.singleton_seqs1),
            self.denoiser_mapping1,
            parse_fasta(self.raw_seqs1)))
        expected = [("S1_0 FXX111 some comments", "TTTT"),
                    ("S1_2 FXX113 some other comments", "TTTT"),
                    ("S3_7 FXX117", "TTTT"),
                    ("S2_1 FXX112 some comments", "TATT"),
                    ("S3_5 FXX114", "TATT"),
                    ("S3_6 FXX115", "TTGA"),
                    ("S3_6 FXX116", "TAGA")]
        self.assertEqual(actual, expected)

    def test_flowgram_id_to_seq_id_map(self):
        """ flowgram_id_to_seq_id_map functions as expected """
        actual = flowgram_id_to_seq_id_map(parse_fasta(self.raw_seqs1))
        expected = {'FXX111': 'S1_0 FXX111 some comments',
                    'FXX112': 'S2_1 FXX112 some comments',
                    'FXX113': 'S1_2 FXX113 some other comments',
                    'FXX114': 'S3_5 FXX114',
                    'FXX115': 'S3_6 FXX115',
                    'FXX116': 'S3_6 FXX116',
                    'FXX117': 'S3_7 FXX117'}
        self.assertEqual(actual, expected)

    def test_count_seqs_from_file(self):
        """ count_seqs: functions as expected with varied data
        """
        f1 = ['>seq1', 'AACCTT', 'ACTGGT',
              '>seq2', 'CCAATT',
              '>seq3', 'CCC---GG']
        f2 = ['> s42', 'ABCDEFG',
              '>s33', 'A',
              '> 4>', 'AA>',
              '>blah', 'AA']
        assert_almost_equal(
            count_seqs_from_file(f1),
            (3,
             8.666,
             2.4944), decimal=3)
        assert_almost_equal(
            count_seqs_from_file(f2),
            (4,
             3.25,
             2.2776), decimal=3)
        self.assertEqual(count_seqs_from_file([]), (0, None, None))

    def test_count_seqs(self):
        """ count_seqs functions as expected with fake seq_counter
        """
        def seq_counter(filepath, parser=None):
            # Fake sequence counter to test count_seqs without
            # having to write files to disk (note don't need to
            # test actual sequence counters here as they're tested
            # elsewhere)
            if filepath.startswith('fake'):
                raise IOError
            else:
                return len(filepath), 0, 0

        in_fps = ['1.fasta', 'fake1.fasta', 'fake.fasta', '2.fa']
        expected = [((7, 0, 0), '1.fasta'),
                    ((4, 0, 0), '2.fa')],\
            11, ['fake1.fasta', 'fake.fasta']
        self.assertEqual(count_seqs_in_filepaths(
            in_fps, seq_counter), expected)

        in_fps = ['fake1.fasta', 'fake.fasta']
        expected = [], 0, ['fake1.fasta', 'fake.fasta']
        self.assertEqual(count_seqs_in_filepaths(
            in_fps, seq_counter), expected)

        in_fps = ['1.fasta', '2.fa', '12.txt']
        expected = [((7, 0, 0), '1.fasta'),
                    ((4, 0, 0), '2.fa'),
                    ((6, 0, 0), '12.txt')], 17, []
        self.assertEqual(count_seqs_in_filepaths(
            in_fps, seq_counter), expected)


    def test_make_compatible_distance_matrices(self):
        """make_compatible_distance_matrices: functions as expected"""
        dm1 = (['A', 'B', 'C', 'D'],
               array([[0.0, 2.3, 3.3, 4.3],
                      [2.9, 0.0, 5.3, 6.3],
                      [3.9, 5.9, 0.0, 4.3],
                      [4.9, 8.9, 9.9, 0.0]]))

        dm2 = (['C', 'A', 'T'],
               array([[10.0, 12.3, 13.3],
                      [12.9, 10.0, 15.3],
                      [13.9, 15.9, 10.0]]))

        expected_dm1 = (['A', 'C'],
                        array([[0.0, 3.3],
                               [3.9, 0.0]]))

        expected_dm2 = (['A', 'C'],
                        array([[10.0, 12.9],
                               [12.3, 10.0]]))

        actual_dm1, actual_dm2 = make_compatible_distance_matrices(dm1, dm2)
        self.assertItemsEqual(actual_dm1[0], expected_dm1[0])
        assert_almost_equal(actual_dm1[1], expected_dm1[1])
        self.assertItemsEqual(actual_dm2[0], expected_dm2[0])
        assert_almost_equal(actual_dm2[1], expected_dm2[1])

    def test_make_compatible_distance_matrices_w_lookup(self):
        """make_compatible_distance_matrices: functions as expected with lookup"""
        dm1 = (['A', 'B', 'C', 'D'],
               array([[0.0, 2.3, 3.3, 4.3],
                      [2.9, 0.0, 5.3, 6.3],
                      [3.9, 5.9, 0.0, 4.3],
                      [4.9, 8.9, 9.9, 0.0]]))

        dm2 = (['C', 'A', 'T'],
               array([[10.0, 12.3, 13.3],
                      [12.9, 10.0, 15.3],
                      [13.9, 15.9, 10.0]]))

        lookup = {'C': 'C', 'A': 'A', 'B': 'B', 'T': 'B', 'D': 'D'}

        expected_dm1 = (['A', 'B', 'C'],
                        array([[0.0, 2.3, 3.3],
                               [2.9, 0.0, 5.3],
                               [3.9, 5.9, 0.0]]))

        expected_dm2 = (['A', 'B', 'C'],
                        array([[10.0, 15.3, 12.9],
                               [15.9, 10.0, 13.9],
                               [12.3, 13.3, 10.0]]))

        actual_dm1, actual_dm2 = make_compatible_distance_matrices(
            dm1, dm2, lookup)
        self.assertItemsEqual(actual_dm1[0], expected_dm1[0])
        assert_almost_equal(actual_dm1[1], expected_dm1[1])
        self.assertItemsEqual(actual_dm2[0], expected_dm2[0])
        assert_almost_equal(actual_dm2[1], expected_dm2[1])

        lookup = {'C': 'C', 'B': 'B', 'T': 'B', 'D': 'D'}
        self.assertRaises(KeyError,
                          make_compatible_distance_matrices, dm1, dm2, lookup)

    def test_trim_fasta(self):
        """ trim_fasta functions as expected """
        expected = [""">HWUSI-EAS552R_0357:8:1:10040:6364#0/1
GACGAG
""",
                    """>HWUSI-EAS552R_0357:8:1:10184:6365#0/1
GTCTGA
"""]

        self.assertEqual(list(trim_fasta(self.fasta_barcodes, 6)), expected)

    def test_duplicates_indices(self):
        """ Properly returns dict of duplicates and their indices """

        no_dups = ['1', '2', '3', '4']

        results = duplicates_indices(no_dups)

        expected_results = defaultdict(list)

        self.assertEqual(results, expected_results)

        dups = ['1', '2', '3', '4', '2']

        results = duplicates_indices(dups)

        expected_results = defaultdict(list)
        expected_results['2'] = [1, 4]

        self.assertEqual(results, expected_results)

    def test_add_filename_suffix(self):
        """Test adding a suffix to a filename works correctly."""
        self.assertEqual(add_filename_suffix('/foo/bar/baz.txt', 'z'),
                         'bazz.txt')
        self.assertEqual(add_filename_suffix('baz.txt', 'z'),
                         'bazz.txt')
        self.assertEqual(add_filename_suffix('/foo/bar/baz', 'z'),
                         'bazz')
        self.assertEqual(add_filename_suffix('baz', 'z'),
                         'bazz')
        self.assertEqual(add_filename_suffix('/baz.fasta.txt', 'z'),
                         'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('baz.fasta.txt', 'z'),
                         'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('/foo/', 'z'), 'z')

    def test_is_valid_git_refname(self):
        """Test correct validation of refnames"""
        # valid branchnames
        self.assertTrue(is_valid_git_refname('master'))
        self.assertTrue(is_valid_git_refname('debuggatron_2000'))
        self.assertTrue(is_valid_git_refname('refname/bar'))
        self.assertTrue(is_valid_git_refname('ref.nameslu/_eggs_/spam'))
        self.assertTrue(is_valid_git_refname('valid{0}char'.format(
            unichr(40))))
        self.assertTrue(is_valid_git_refname('master@head'))
        self.assertTrue(is_valid_git_refname('bar{thing}foo'))

        # case happening with git < 1.6.6
        self.assertFalse(is_valid_git_refname(
            '--abbrev-ref\nbaa350d7b7063d585ca293fc16ef15e0765dc9ee'))

        # different invalid refnames, for a description of each group see the
        # man page of git check-ref-format
        self.assertFalse(is_valid_git_refname('bar/.spam/eggs'))
        self.assertFalse(is_valid_git_refname('bar.lock/spam/eggs'))
        self.assertFalse(is_valid_git_refname('bar.lock'))
        self.assertFalse(is_valid_git_refname('.foobar'))

        self.assertFalse(is_valid_git_refname('ref..name'))

        self.assertFalse(is_valid_git_refname(u'invalid{0}char'.format(
            unichr(177))))
        self.assertFalse(is_valid_git_refname('invalid{0}char'.format(
            unichr(39))))
        self.assertFalse(is_valid_git_refname('ref~name/bar'))
        self.assertFalse(is_valid_git_refname('refname spam'))
        self.assertFalse(is_valid_git_refname('bar/foo/eggs~spam'))
        self.assertFalse(is_valid_git_refname('bar:_spam_'))
        self.assertFalse(is_valid_git_refname('eggtastic^2'))

        self.assertFalse(is_valid_git_refname('areyourandy?'))
        self.assertFalse(is_valid_git_refname('bar/*/spam'))
        self.assertFalse(is_valid_git_refname('bar[spam]/eggs'))

        self.assertFalse(is_valid_git_refname('/barfooeggs'))
        self.assertFalse(is_valid_git_refname('barfooeggs/'))
        self.assertFalse(is_valid_git_refname('bar/foo//////eggs'))

        self.assertFalse(is_valid_git_refname('dotEnding.'))

        self.assertFalse(is_valid_git_refname('@{branch'))

        self.assertFalse(is_valid_git_refname('contains\\slash'))

        self.assertFalse(is_valid_git_refname('$newbranch'))

    def test_is_valid_git_sha1(self):
        """ """

        # valid sha1 strings
        self.assertTrue(is_valid_git_sha1(
            '65a9ba2ef4b126fb5b054ea6b89b457463db4ec6'))
        self.assertTrue(is_valid_git_sha1(
            'a29a9911e41253405494c43889925a6d79ca26db'))
        self.assertTrue(is_valid_git_sha1(
            'e099cd5fdea89eba929d6051fbd26cc9e7a0c961'))
        self.assertTrue(is_valid_git_sha1(
            '44235d322c3386bd5ce872d9d7ea2e10d27c86cb'))
        self.assertTrue(is_valid_git_sha1(
            '7d2fc23E04540EE92c742948cca9ed5bc54d08d1'))
        self.assertTrue(is_valid_git_sha1(
            'fb5dc0285a8b11f199c4f3a7547a2da38138373f'))
        self.assertTrue(is_valid_git_sha1(
            '0b2abAEb195ba7ebc5cfdb53213a66fbaddefdb8'))

        # invalid length
        self.assertFalse(is_valid_git_sha1('cca9ed5bc54d08d1'))
        self.assertFalse(is_valid_git_sha1(''))

        # invalid characters
        self.assertFalse(is_valid_git_sha1(
            'fb5dy0f85a8b11f199c4f3a75474a2das8138373'))
        self.assertFalse(is_valid_git_sha1(
            '0x5dcc816fbc1c2e8eX087d7d2ed8d2950a7c16b'))

    def test_invert_dict(self):
        """invert_dict should invert keys and values, keeping all keys

        Ported from PyCogent's cogent.util.misc.InverseDictMulti unit tests.
        """
        self.assertEqual(invert_dict({}), {})
        self.assertEqual(invert_dict({'3':4}), {4:['3']})
        self.assertEqual(invert_dict(\
            {'a':'x','b':1,'c':None,'d':('a','b')}), \
            {'x':['a'],1:['b'],None:['c'],('a','b'):['d']})
        self.assertRaises(TypeError, invert_dict, {'a':['a','b','c']})
        d = invert_dict({'a':3, 'b':3, 'c':3, 'd':'3', 'e':'3'})
        self.assertEqual(len(d), 2)
        assert 3 in d
        d3_items = d[3][:]
        self.assertEqual(len(d3_items), 3)
        d3_items.sort()
        self.assertEqual(''.join(d3_items), 'abc')
        assert '3' in d
        d3_items = d['3'][:]
        self.assertEqual(len(d3_items), 2)
        d3_items.sort()
        self.assertEqual(''.join(d3_items), 'de')


raw_seqs1 = """>S1_0 FXX111 some comments
TTTT
>S2_1 FXX112 some comments
TATT
>S1_2 FXX113 some other comments
GGGG
>S3_5 FXX114
GGGA
>S3_6 FXX115
TTGA
>S3_6 FXX116
TAGA
>S3_7 FXX117
TAGT"""

centroid_seqs1 = """>FXX111 | cluster size: 3
TTTT
>FXX112 | cluster size: 2
TATT"""

singleton_seqs1 = """>FXX115
TTGA
>FXX116
TAGA"""

denoiser_mapping1 = """FXX111:\tFXX113\tFXX117
FXX115:
FXX112:\tFXX114
FXX116:"""


otu_table_fake1 = """#Full OTU Counts
#OTU ID	S1	S2	Consensus Lineage
0	1	0	Root;Bacteria
1	1	0	Root;Bacteria;Verrucomicrobia
2	4	0	Root;Bacteria"""

otu_table_fake2 = """#Full OTU Counts
#OTU ID	S3	S4	S5	Consensus Lineage
0	1	0	1	Root;Bacteria
3	2	0	1	Root;Bacteria;Acidobacteria
4	1	0	9	Root;Bacteria;Bacteroidetes
2	1	0	1	Root;Bacteria;Acidobacteria;Acidobacteria;Gp5
6	1	25	42	Root;Archaea"""

otu_table_fake3 = """#Full OTU Counts
#OTU ID	samp7	Consensus Lineage
6	1	Root;Archaea"""

otu_table_fake4 = """#Full OTU Counts
#OTU ID	S3	S4	S1	Consensus Lineage
0	1	0	1	Root;Bacteria
3	2	0	1	Root;Bacteria;Acidobacteria
4	1	0	9	Root;Bacteria;Bacteroidetes
2	1	0	1	Root;Bacteria;Acidobacteria;Acidobacteria;Gp5
6	1	25	42	Root;Archaea"""

otu_table_fake1_no_tax = """#Full OTU Counts
#OTU ID	S1	S2
0	1	0
1	1	0
2	4	0"""

otu_table_fake2_no_tax = """#Full OTU Counts
#OTU ID	S3	S4	S5
0	1	0	1
3	2	0	1
4	1	0	9
2	1	0	1
6	1	25	42"""

otu_table_fake3_no_tax = """#Full OTU Counts
#OTU ID	samp7
6	1"""


class FunctionWithParamsTests(TestCase):

    """Tests of the FunctionWithParams class.

    Note: FunctionWithParams itself is abstract, so need to subclass here."""

    def setUp(self):
        """Set up standard items for testing."""
        class FWP(FunctionWithParams):
            Name = 'FWP'
            t = True
            _tracked_properties = ['t']

            def getResult(self):
                return 3
        self.FWP = FWP

        self.files_to_remove = []

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_init(self):
        """FunctionWithParams __init__ should be successful."""
        x = self.FWP({'x': 3})
        self.assertEqual(x.Name, 'FWP')
        self.assertEqual(x.Params, {'x': 3})

    def test_str(self):
        """FunctionWithParams __str__ should produce expected string"""
        x = self.FWP({'x': 3})
        lines = ['FWP parameters:', 't:True', 'Application:None',
                 'Algorithm:None', 'Citation:None', 'x:3']
        self.assertEqual(str(x), '\n'.join(lines))

    def test_call(self):
        """FunctionWithParams __str__ should produce expected string"""
        x = self.FWP({'x': 3})
        self.assertEqual(x(), 3)

    def test_formatResult(self):
        """FunctionWithParams formatResult should produce expected format"""
        x = self.FWP({'x': 3})
        self.assertEqual(x.formatResult(3), '3')

    def test_getBiomData(self):
        """FunctionWithParams getBiomData should return biom"""
        bt_string = '''{
        "id":null,
        "format": "Biological Observation Matrix 0.9.1-dev",
        "format_url": "http://biom-format.org",
        "type": "OTU table",
        "generated_by": "QIIME revision XYZ",
        "date": "2011-12-19T19:00:00",
        "rows":[
                {"id":"GG_OTU_1", "metadata":null},
                {"id":"GG_OTU_2", "metadata":null},
                {"id":"GG_OTU_3", "metadata":null},
                {"id":"GG_OTU_4", "metadata":null},
                {"id":"GG_OTU_5", "metadata":null}
            ],
        "columns": [
                {"id":"Sample1", "metadata":null},
                {"id":"Sample2", "metadata":null},
                {"id":"Sample3", "metadata":null},
                {"id":"Sample4", "metadata":null},
                {"id":"Sample5", "metadata":null},
                {"id":"Sample6", "metadata":null}
            ],
        "matrix_type": "dense",
        "matrix_element_type": "int",
        "shape": [5,6],
        "data":  [[0,0,1,0,0,0],
                  [5,1,0,2,3,1],
                  [0,0,1,4,2,0],
                  [2,1,1,0,0,1],
                  [0,1,1,0,0,0]]
    }'''
        biom_data = parse_biom_table_str(bt_string)
        F = FunctionWithParams('')

        self.assertEqual(biom_data, F.getBiomData(biom_data))

        # write biom_data to temp location
        fd, bt_path = mkstemp(suffix='.biom')
        close(fd)
        biom_file = open(bt_path, 'w')
        biom_file.writelines(bt_string)
        biom_file.close()
        self.assertEqual(biom_data, F.getBiomData(bt_path))

        # cleanup
        remove(bt_path)

    def test_convert_OTU_table_relative_abundance(self):
        """convert_OTU_table_relative_abundance works
        """
        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""".split('\n')
        result = convert_OTU_table_relative_abundance(otu_table)
        self.assertEqual(
            result,
            ['#Full OTU Counts',
             '#OTU ID\tsample1\tsample2\tsample3',
             '0\t0.0\t0.666666666667\t0.0',
             '1\t0.5\t0.0\t0.0',
             '2\t0.5\t0.333333333333\t1.0'])

        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage
0\t0\t2\t0\tBacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475
1\t1\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369
2\t1\t1\t1\tBacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121""".split('\n')
        result = convert_OTU_table_relative_abundance(otu_table)
        self.assertEqual(
            result,
            ['#Full OTU Counts',
             '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage',
             '0\t0.0\t0.666666666667\t0.0\tBacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475',
             '1\t0.5\t0.0\t0.0\tBacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369',
             '2\t0.5\t0.333333333333\t1.0\tBacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121'])

    def test_flip_vectors(self):
        """_flip_vectors makes a new PCA matrix with correct signs"""
        m_matrix = array([[1.0, 0.0, 1.0], [2.0, 4.0, 4.0]])
        jn_matrix = array([[1.2, 0.1, -1.2], [2.5, 4.0, -4.5]])
        new_matrix = _flip_vectors(jn_matrix, m_matrix)
        assert_almost_equal(new_matrix, array([[1.2, 0.1, 1.2], [2.5, 4.0, 4.5]]))

    def test_compute_jn_pcoa_avg_ranges(self):
        """_compute_jn_pcoa_avg_ranges works
        """
        jn_flipped_matrices = [array([[2.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[3.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[4.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[5.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[6.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[7.0, 4.0, -4.5], [-1.2, -0.1, 1.2]]),
                               array([[1.0, 4.0, -4.5], [-1.2, -0.1, 1.2]])]
        avg_matrix, low_matrix, high_matrix = _compute_jn_pcoa_avg_ranges(
            jn_flipped_matrices, 'ideal_fourths')
        assert_almost_equal(avg_matrix[(0, 0)], 4.0)
        assert_almost_equal(avg_matrix[(0, 2)], -4.5)
        assert_almost_equal(low_matrix[(0, 0)], 2.16666667)
        assert_almost_equal(high_matrix[(0, 0)], 5.83333333)

        avg_matrix, low_matrix, high_matrix = _compute_jn_pcoa_avg_ranges(
            jn_flipped_matrices, 'sdev')
        x = array([m[0, 0] for m in jn_flipped_matrices])
        self.assertEqual(x.mean(), avg_matrix[0, 0])
        self.assertEqual(-x.std(ddof=1) / 2, low_matrix[0, 0])
        self.assertEqual(x.std(ddof=1) / 2, high_matrix[0, 0])

    def test_summarize_pcoas(self):
        """summarize_pcoas works
        """
        master_pcoa = [['1', '2', '3'],
                       array([[-1.0, 0.0, 1.0], [2.0, 4.0, -4.0]]),
                       array([.76, .24])]
        jn1 = [['1', '2', '3'],
               array([[1.2, 0.1, -1.2], [-2.5, -4.0, 4.5]]),
               array([0.80, .20])]
        jn2 = [['1', '2', '3'],
               array([[-1.4, 0.05, 1.3], [2.6, 4.1, -4.7]]),
               array([0.76, .24])]
        jn3 = [['1', '2', '3'],
               array([[-1.5, 0.05, 1.6], [2.4, 4.0, -4.8]]),
               array([0.84, .16])]
        jn4 = [['1', '2', '3'],
               array([[-1.5, 0.05, 1.6], [2.4, 4.0, -4.8]]),
               array([0.84, .16])]
        support_pcoas = [jn1, jn2, jn3, jn4]
        # test with the ideal_fourths option
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, 'ideal_fourths',
                            apply_procrustes=False)
        self.assertEqual(m_names, ['1', '2', '3'])
        assert_almost_equal(matrix_average[(0, 0)], -1.4)
        assert_almost_equal(matrix_average[(0, 1)], 0.0125)
        assert_almost_equal(matrix_low[(0, 0)], -1.5)
        assert_almost_equal(matrix_high[(0, 0)], -1.28333333)
        assert_almost_equal(matrix_low[(0, 1)], -0.0375)
        assert_almost_equal(matrix_high[(0, 1)], 0.05)
        assert_almost_equal(eigval_average[0], 0.81)
        assert_almost_equal(eigval_average[1], 0.19)
        # test with the IQR option
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, method='IQR',
                            apply_procrustes=False)
        assert_almost_equal(matrix_low[(0, 0)], -1.5)
        assert_almost_equal(matrix_high[(0, 0)], -1.3)

        # test with procrustes option followed by sdev
        m, m1, msq = procrustes(master_pcoa[1], jn1[1])
        m, m2, msq = procrustes(master_pcoa[1], jn2[1])
        m, m3, msq = procrustes(master_pcoa[1], jn3[1])
        m, m4, msq = procrustes(master_pcoa[1], jn4[1])
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, method='sdev',
                            apply_procrustes=True)

        x = array([m1[0, 0], m2[0, 0], m3[0, 0], m4[0, 0]])
        self.assertEqual(x.mean(), matrix_average[0, 0])
        self.assertEqual(-x.std(ddof=1) / 2, matrix_low[0, 0])
        self.assertEqual(x.std(ddof=1) / 2, matrix_high[0, 0])

    def test_IQR(self):
        "IQR returns the interquartile range for list x"
        # works for odd with odd split
        x = [2, 3, 4, 5, 6, 7, 1]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 6)
        # works for even with odd split
        x = [1, 2, 3, 4, 5, 6]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 5)
        # works for even with even split
        x = [1, 2, 3, 4, 5, 6, 7, 8]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2.5)
        self.assertEqual(maxv, 6.5)
        # works with array
        # works for odd with odd split
        x = array([2, 3, 4, 5, 6, 7, 1])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 6)
        # works for even with odd split
        x = array([1, 2, 3, 4, 5, 6])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 5)
        # works for even with even split
        x = array([1, 2, 3, 4, 5, 6, 7, 8])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2.5)
        self.assertEqual(maxv, 6.5)

    def test_matrix_IQR(self):
        """matrix_IQR calcs the IQR for each column in an array correctly
        """
        x = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
        min_vals, max_vals = matrix_IQR(x)
        assert_almost_equal(min_vals, array([2.5, 3.5, 4.5]))
        assert_almost_equal(max_vals, array([8.5, 9.5, 10.5]))

    def test_idealfourths(self):
        """idealfourths: tests the ideal-fourths function which was imported from scipy
        at the following location (http://projects.scipy.org/scipy/browser/trunk/scipy/stats/tests/test_mmorestats.py?rev=4154)
        """
        test = numpy.arange(100)
        self.assertEqual(idealfourths(test),
                         [24.416666666666668, 74.583333333333343])
        test_2D = test.repeat(3).reshape(-1, 3)

        assert_almost_equal(numpy.asarray(idealfourths(test_2D, axis=0)),
                                 numpy.array(
                                     [[24.41666667, 24.41666667, 24.41666667],
                                      [74.58333333, 74.58333333, 74.58333333]]))

        assert_almost_equal(idealfourths(test_2D, axis=1),
                         test.repeat(2).reshape(-1, 2))
        test = [0, 0]
        _result = idealfourths(test)
        assert_almost_equal(numpy.isnan(_result).all(), True)

    def test_isarray(self):
        "isarray: tests the isarray function"
        obs1 = asarray('Test')
        obs2 = 'Test'

        exp1 = True
        exp2 = False

        self.assertEqual(isarray(obs1), exp1)
        self.assertEqual(isarray(obs2), exp2)

    def test_degap_fasta_aln(self):
        """degap_fasta_align removes gps from fasta seqs."""

        test_aln = [("a", "AAAAAAAAAGGGG"),
                    ("b", "-A-A-G-G-A-G-C."),
                    ('c', "..-----------"),
                    ('d', "--------AAAAAAA"),
                    ('e', "")]

        expected_result = map(lambda a_b: DNASequence(a_b[1],
                                                      id=a_b[0]),
                              [("a", "AAAAAAAAAGGGG"),
                               ("b", "AAGGAGC"),
                               ('c', ""),
                               ('d', "AAAAAAA"),
                               ('e', "")])

        self.assertEqual(list(degap_fasta_aln(test_aln)), expected_result)

        self.assertEqual(list(degap_fasta_aln([])), [])

    def test_write_degapped_fasta_to_file(self):

        test_aln = [("a", "AAAAAAAAAGGGG"),
                    ("b", "-A-A-G-G-A-G-C"),
                    ('c', "---------------"),
                    ('d', "--------AAAAAAA"),
                    ('e', "")]

        expected_result = """>a
AAAAAAAAAGGGG
>b
AAGGAGC
>c

>d
AAAAAAA
>e

"""
        tmp_filename = write_degapped_fasta_to_file(test_aln)
        self.files_to_remove.append(tmp_filename)
        observed = "".join(list(open(tmp_filename, "U")))
        self.assertEqual(observed, expected_result)

    def test_get_diff_for_otu_maps(self):
        """get_diff_for_otu_map return correct set difference"""

        # compare to self
        self.assertEqual(get_diff_for_otu_maps(otu_map1, otu_map1),
                         (set([]), set([])))

        # compare to otu_map with one difference
        self.assertEqual(get_diff_for_otu_maps(otu_map1, otu_map2),
                         (set(['b']), set([])))

        # compare to empty
        self.assertEqual(get_diff_for_otu_maps(otu_map1, fields_to_dict("")),
                         (set(['a', 'b', 'c', 'd', 'e', 'f']), set([])))

    def test_compare_otu_maps(self):
        """compare_otu_maps computes correct values"""

        assert_almost_equal(compare_otu_maps(otu_map1, otu_map1), 0.0)
        assert_almost_equal(compare_otu_maps(otu_map1, otu_map3), 0.0)
        assert_almost_equal(
            compare_otu_maps(
                otu_map1,
                otu_map4),
            0.33333333333)
        assert_almost_equal(
            compare_otu_maps(
                otu_map3,
                otu_map4),
            0.33333333333)
        assert_almost_equal(compare_otu_maps(otu_map1, otu_map5), 1)

    def test_iseq_to_qseq_fields(self):
        """iseq_to_qseq_fields functions as expected"""
        i = "HWI-ST753_50:6:1101:15435:9071#0/1:ACCAGACGATGCTACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC:gggggggfggdegggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg"
        # barcode in sequence, barcode length = 12
        expected = (
            ("HWI-ST753", "50", "6", "1101", "15435", "9071", "0", "1"),
            "TACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC", "gggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg", "ACCAGACGATGC", "gggggggfggde")
        self.assertEqual(
            iseq_to_qseq_fields(i, barcode_in_header=False, barcode_length=12),
            expected)
        # barcode in sequence, barcode length = 6
        expected = (
            ("HWI-ST753", "50", "6", "1101", "15435", "9071", "0", "1"),
            "CGATGCTACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC", "gfggdegggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg", "ACCAGA", "gggggg")
        self.assertEqual(
            iseq_to_qseq_fields(i, barcode_in_header=False, barcode_length=6),
            expected)

        # barcode in header, barcode length = 6
        i = "HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB"
        expected = (("HWI-6X", "9267", "1", "1", "4", "1699", "ACCACCC", "1"),
                    "TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA", "abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB", "ACCACC", "bbbbbb")
        self.assertEqual(
            iseq_to_qseq_fields(i, barcode_in_header=True, barcode_length=6),
            expected)
        # barcode in header, barcode length = 3
        expected = (("HWI-6X", "9267", "1", "1", "4", "1699", "ACCACCC", "1"),
                    "TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA", "abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB", "ACC", "bbb")
        self.assertEqual(
            iseq_to_qseq_fields(i, barcode_in_header=True, barcode_length=3),
            expected)

    def test_stderr(self):
        """stderr computes standard error for an array"""

        exp = array([0.57735026918962584, 0.57735026918962584,
                     0.57735026918962584, 0.57735026918962584])
        obs = stderr([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3]])
        assert_almost_equal(obs, exp)

    def test__chk_asarray(self):
        """_chk_asarray converts list into a numpy array"""

        exp = (array([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3]]), 0)
        obs = _chk_asarray([[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3]], 0)
        assert_almost_equal(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])

class BlastSeqsTests(TestCase):

    """ Tests of the qiime_blast_seqs function (will move to PyCogent eventually)
    """

    def setUp(self):
        """
        """
        self.refseqs1 = refseqs1.split('\n')
        self.inseqs1 = inseqs1.split('\n')
        self.blast_db, db_files_to_remove =\
            build_blast_db_from_fasta_file(self.refseqs1,
                                           output_dir=get_qiime_temp_dir())
        self.files_to_remove = db_files_to_remove

        fd, self.refseqs1_fp = mkstemp(dir=get_qiime_temp_dir(),
                                      prefix="BLAST_temp_db_",
                                      suffix=".fasta")
        close(fd)
        fasta_f = open(self.refseqs1_fp, 'w')
        fasta_f.write(refseqs1)
        fasta_f.close()

        self.files_to_remove = db_files_to_remove + [self.refseqs1_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_w_refseqs_file(self):
        """qiime_blast_seqs functions with refseqs file
        """
        inseqs = parse_fasta(self.inseqs1)
        actual = qiime_blast_seqs(inseqs, refseqs=self.refseqs1)
        self.assertEqual(len(actual), 5)

        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'], 's2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'], 's1')

    def test_w_refseqs_fp(self):
        """qiime_blast_seqs functions refseqs_fp
        """
        inseqs = parse_fasta(self.inseqs1)
        actual = qiime_blast_seqs(inseqs, refseqs_fp=self.refseqs1_fp)
        self.assertEqual(len(actual), 5)

        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'], 's2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'], 's1')

    def test_w_preexising_blastdb(self):
        """qiime_blast_seqs functions with pre-existing blast_db
        """
        # pre-existing blast db
        inseqs = parse_fasta(self.inseqs1)
        actual = qiime_blast_seqs(inseqs, blast_db=self.blast_db)
        self.assertEqual(len(actual), 5)

        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'], 's2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'], 's1')

    def test_w_alt_seqs_per_blast_run(self):
        """qiime_blast_seqs: functions with alt seqs_per_blast_run
        """
        for i in range(1, 20):
            inseqs = parse_fasta(self.inseqs1)
            actual = qiime_blast_seqs(
                inseqs, blast_db=self.blast_db, seqs_per_blast_run=i)
            self.assertEqual(len(actual), 5)

            # couple of sanity checks against command line blast
            self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'], 's2')
            self.assertEqual(actual['s105'][0][2]['SUBJECT ID'], 's1')

    def test_alt_blast_param(self):
        """qiime_blast_seqs: alt blast params give alt results"""
        # Fewer blast hits with stricter e-value
        inseqs = parse_fasta(self.inseqs1)
        actual = qiime_blast_seqs(
            inseqs,
            blast_db=self.blast_db,
            params={'-e': '1e-4'})
        # NOTE: A BUG (?) IN THE BlastResult OR THE PARSER RESULTS IN AN EXTRA
        # BLAST RESULT IN THE DICT (actual HERE) WITH KEY ''

        self.assertTrue('s2_like_seq' in actual)
        self.assertFalse('s100' in actual)

    def test_error_on_bad_param_set(self):
        inseqs = parse_fasta(self.inseqs1)
        # no blastdb or refseqs
        self.assertRaises(AssertionError, qiime_blast_seqs, inseqs)


class BlastXSeqsTests(TestCase):

    """ Tests of the qiime_blastx_seqs function (will move to PyCogent eventually)
    """

    def setUp(self):
        """
        """
        self.nt_inseqs1 = nt_inseqs1.split('\n')

        fd, self.pr_refseqs1_fp = mkstemp(dir=get_qiime_temp_dir(),
                                        prefix="BLAST_temp_db_",
                                        suffix=".fasta")
        close(fd)
        fasta_f = open(self.pr_refseqs1_fp, 'w')
        fasta_f.write(pr_refseqs1)
        fasta_f.close()

        self.files_to_remove = [self.pr_refseqs1_fp]

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_w_refseqs_file(self):
        """qiime_blastx_seqs functions with refseqs file
        """
        inseqs = parse_fasta(self.nt_inseqs1)
        actual = qiime_blastx_seqs(inseqs, refseqs_fp=self.pr_refseqs1_fp)
        self.assertEqual(len(actual), 3)

        # couple of sanity checks against command line blast
        self.assertEqual(actual['eco:b0001'][0][0]['SUBJECT ID'], 'eco:b0001')
        self.assertEqual(actual['eco:b0122'][0][0]['SUBJECT ID'], 'eco:b0122')
        self.assertEqual(actual['eco:b0122'][0][1]['SUBJECT ID'], 'eco:b0015')


pr_refseqs1 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (A)
MKRISTTITTTITITTGNGAG
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (A)
MAKQDYYEILGVSKTAEEREIRKAYKRLAMKYHPDRNQGDKEAEAKFKEIKEAYEVLTDS
QKRAAYDQYGHAAFEQGGMGGGGFGGGADFSDIFGDVFGDIFGGGRGRQRAARGADLRYN
MELTLEEAVRGVTKEIRIPTLEECDVCHGSGAKPGTQPQTCPTCHGSGQVQMRQGFFAVQ
QTCPHCQGRGTLIKDPCNKCHGHGRVERSKTLSVKIPAGVDTGDRIRLAGEGEAGEHGAP
AGDLYVQVQVKQHPIFEREGNNLYCEVPINFAMAALGGEIEVPTLDGRVKLKVPGETQTG
KLFRMRGKGVKSVRGGAQGDLLCRVVVETPVGLNERQKQLLQELQESFGGPTGEHNSPRS
KSFFDGVKKFFDDLTR
>eco:b0122 yacC; conserved protein, PulS_OutS family (A)
MKTFFRTVLFGSLMAVCANSYALSESEAEDMADLTAVFVFLKNDCGYQNLPNGQIRRALV
FFAQQNQWDLSNYDTFDMKALGEDSYRDLSGIGIPVAKKCKALARDSLSLLAYVK"""

nt_inseqs1 = """>eco:b0001 thrL; thr operon leader peptide; K08278 thr operon leader peptide (N)
atgaaacgcattagcaccaccattaccaccaccatcaccattaccacaggtaacggtgcg
ggctga
>eco:b0015 dnaJ; chaperone Hsp40, co-chaperone with DnaK; K03686 molecular chaperone DnaJ (N)
atggctaagcaagattattacgagattttaggcgtttccaaaacagcggaagagcgtgaa
atcagaaaggcctacaaacgcctggccatgaaataccacccggaccgtaaccagggtgac
aaagaggccgaggcgaaatttaaagagatcaaggaagcttatgaagttctgaccgactcg
caaaaacgtgcggcatacgatcagtatggtcatgctgcgtttgagcaaggtggcatgggc
ggcggcggttttggcggcggcgcagacttcagcgatatttttggtgacgttttcggcgat
atttttggcggcggacgtggtcgtcaacgtgcggcgcgcggtgctgatttacgctataac
atggagctcaccctcgaagaagctgtacgtggcgtgaccaaagagatccgcattccgact
ctggaagagtgtgacgtttgccacggtagcggtgcaaaaccaggtacacagccgcagact
tgtccgacctgtcatggttctggtcaggtgcagatgcgccagggattcttcgctgtacag
cagacctgtccacactgtcagggccgcggtacgctgatcaaagatccgtgcaacaaatgt
catggtcatggtcgtgttgagcgcagcaaaacgctgtccgttaaaatcccggcaggggtg
gacactggagaccgcatccgtcttgcgggcgaaggtgaagcgggcgagcatggcgcaccg
gcaggcgatctgtacgttcaggttcaggttaaacagcacccgattttcgagcgtgaaggc
aacaacctgtattgcgaagtcccgatcaacttcgctatggcggcgctgggtggcgaaatc
gaagtaccgacccttgatggtcgcgtcaaactgaaagtgcctggcgaaacccagaccggt
aagctattccgtatgcgcggtaaaggcgtcaagtctgtccgcggtggcgcacagggtgat
ttgctgtgccgcgttgtcgtcgaaacaccggtaggcctgaacgaaaggcagaaacagctg
ctgcaagagctgcaagaaagcttcggtggcccaaccggcgagcacaacagcccgcgctca
aagagcttctttgatggtgtgaagaagttttttgacgacctgacccgctaa
>eco:b0122 yacC; conserved protein, PulS_OutS family (N)
atgaagacgtttttcagaacagtgttattcggcagcctgatggccgtctgcgcaaacagt
tacgcgctcagcgagtctgaagccgaagatatggccgatttaacggcagtttttgtcttt
ctgaagaacgattgtggttaccagaacttacctaacgggcaaattcgtcgcgcactggtc
tttttcgctcagcaaaaccagtgggacctcagtaattacgacaccttcgacatgaaagcc
ctcggtgaagacagctaccgcgatctcagcggcattggcattcccgtcgctaaaaaatgc
aaagccctggcccgcgattccttaagcctgcttgcctacgtcaaataa"""

otu_map1 = fields_to_dict("""1:\ta\tb\tc
2:\td
3:\te\tf
""".split("\n"))

# b missing
otu_map2 = fields_to_dict("""1:\ta\tc
3:\te\tf\td
""".split("\n"))

# several reads swapped
otu_map3 = fields_to_dict("""1:\tc\ta\tb
3:\te\tf
2:\td
""".split("\n"))

# several reads swapped
otu_map4 = fields_to_dict("""1:\tc\ta\tb\tf
3:\te\td
""".split("\n"))

# everything differs
otu_map5 = fields_to_dict("""4:\ta\tb\tc\td\te\tf""".split("\n"))


inseqs1 = """>s2_like_seq
TGCAGCTTGAGCACAGGTTAGAGCCTTC
>s100
TGCAGCTTGAGCACAGGTTAGCCTTC
>s101
TGCAGCTTGAGCACAGGTTTTTCAGAGCCTTC
>s104
TGCAGCTTGAGCACAGGTTAGCCTTC
>s105
TGCAGCTTGAGCACAGGTTAGATC"""

refseqs1 = """>s0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
>s1
TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC
>s2
TGCAGCTTGAGCCACAGGAGAGAGCCTTC
>s3
TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC
>s4
ACCGATGAGATATTAGCACAGGGGAATTAGAACCA
>s5
TGTCGAGAGTGAGATGAGATGAGAACA
>s6
ACGTATTTTAATTTGGCATGGT
>s7
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>s8
CCAGAGCGAGTGAGATAGACACCCAC
"""

mapping_f1 = """#SampleID\tSomething\tdays_since_epoch
Z1\t42\t23
Z2\thello\t10
A\t4\t400000
1\tr\t5.7
NotInOtuTable\tf\t0"""

fasta1 = """>Samp1_42
ACCGGTT
>s2_a_50
GGGCCC
>Samp1_43 some comme_nt
AACCG
>s3_25
AAACCC"""

fasta2 = """>Samp1_42
ACCGGTT
>s2_a_50
GGGCCC
>Samp1_43 some comme_nt
AACCG
>s3_25
AAACCC
>Samp1_44
A"""

mapping_f2 = """#SampleID\tSomething\tdays_since_epoch
Z1\thello\t23
Z2\thello\t23
A\thello\t23
1\thello\t23
NotInOtuTable\thello\t23"""

mapping_f3 = """#SampleID\tSomething\tdays_since_epoch
Z1\t42\t23
Z2\t42\t10
A\t4\t400000
1\t4\t5.7
NotInOtuTable\t9\t5.7"""

fastq_barcodes = ["@HWUSI-EAS552R_0357:8:1:10040:6364#0/1",
                  "GACGAGTCAGTCA",
                  "+HWUSI-EAS552R_0357:8:1:10040:6364#0/1",
                  "hhhhhhhhhhhhh",
                  "@HWUSI-EAS552R_0357:8:1:10184:6365#0/1",
                  "GTCTGACAGTTGA",
                  "+HWUSI-EAS552R_0357:8:1:10184:6365#0/1",
                  "hhhhhhhhhhhhh"]
fasta_barcodes = [">HWUSI-EAS552R_0357:8:1:10040:6364#0/1",
                  "GACGAGTCAGTCA",
                  ">HWUSI-EAS552R_0357:8:1:10184:6365#0/1",
                  "GTCTGACAGTTGA"
                  ]
fastq_seqs = ["@HWUSI-EAS552R_0357:8:1:10040:6364#0/2",
              "TACAGGGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTCTGCAGGTTGCTCACTAAGTCTTTTGTTAAATCTTCGGGCTTAACCCGAAACCTGCAAAAGAAACTAGTGCTCTCGAGTATGGTAGAGGTAAAGGGAATTTCCAG",
              "+HWUSI-EAS552R_0357:8:1:10040:6364#0/2",
              "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhfhhhhhhgghhhhhWhfcffehf]hhdhhhhhgcghhhchhhhfhcfhhgggdfhgdcffadccfdcccca]^b``ccfdd_caccWbb[b_dfdcdeaec`^`^_daba_b_WdY^`",
              "@HWUSI-EAS552R_0357:8:1:10184:6365#0/2",
              "TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGGCTCTTAAGTCGGAGGTGAAATCCCAAGGCTCAACCTTGGAACTGCCTTCGATACTGAGAGTCTTGAGTCCGGAAGAGGTAAGTGGAACTCCAAG",
              "+HWUSI-EAS552R_0357:8:1:10184:6365#0/2",
              "hfhhchhghhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhfghghhhggfhhfghfeghfggfdhdfdfffacbfcddgcfccddbcddbccada_aadWaaaacddccccdacaaa_acbc]c`aa[a\\a_a^V\T_^^^^X^R_BBBB"]

fastq_mapping_rev = [
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
    "sample1\tGACTGACTCGTC\tCCGGACTACHVGGGTWTCTAAT\tsample1",
    "sample2\tCAACTGTCAGAC\tCCGGACTACHVGGGTWTCTAAT\tsample2"]

fastq_mapping_fwd = [
    "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
    "sample1\tGACGAGTCAGTC\tCCGGACTACHVGGGTWTCTAAT\tsample1",
    "sample2\tGTCTGACAGTTG\tCCGGACTACHVGGGTWTCTAAT\tsample2"]


class SubSampleFastaTests(TestCase):

    """ """

    def setUp(self):
        """ """
        self.expected_lines_50_perc = expected_lines_50_perc
        self.expected_lines_20_perc = expected_lines_20_perc

        self.temp_dir = load_qiime_config()['temp_dir']

        self.fasta_lines = fasta_lines
        fd, self.fasta_filepath = mkstemp(
            prefix='subsample_test_', suffix='.fasta')
        close(fd)
        self.fasta_file = open(self.fasta_filepath, "w")
        self.fasta_file.write(self.fasta_lines)
        self.fasta_file.close()

        fd, self.output_filepath = mkstemp(prefix='subsample_output_',
                                          suffix='.fasta')
        close(fd)

        self._files_to_remove =\
            [self.fasta_filepath]

    def tearDown(self):
        remove_files(self._files_to_remove)

    def test_subsample_fasta_50(self):
        """ subsample_fasta correctly subsamples input fasta file """

        # fixed seed for consistent calls with random()
        seed(128)

        subsample_fasta(self.fasta_filepath, self.output_filepath,
                        percent_subsample=0.50)

        self._files_to_remove.append(self.output_filepath)

        actual_results =\
            [line.strip() for line in open(self.output_filepath, "U")]

        self.assertEqual(actual_results, self.expected_lines_50_perc)

    def test_subsample_fasta_20(self):
        """ subsample_fasta correctly subsamples input fasta file """

        seed(12210)

        subsample_fasta(self.fasta_filepath, self.output_filepath,
                        percent_subsample=0.20)

        self._files_to_remove.append(self.output_filepath)

        actual_results =\
            [line.strip() for line in open(self.output_filepath, "U")]

        self.assertEqual(actual_results, self.expected_lines_20_perc)


class MetadataMapTests(TestCase):

    """Tests for the MetadataMap class."""

    def setUp(self):
        """Create MetadataMap objects that will be used in the tests."""
        # Create a map using the overview tutorial mapping file.
        self.overview_map_str = [
            "#SampleID\tBarcodeSequence\tTreatment\tDOB\tDescription",
            "PC.354\tAGCACGAGCCTA\tControl\t20061218\t354",
            "PC.355\tAACTCGTCGATG\tControl\t20061218\t355",
            "PC.356\tACAGACCACTCA\tControl\t20061126\t356",
            "PC.481\tACCAGCGACTAG\tControl\t20070314\t481",
            "PC.593\tAGCAGCACTTGT\tControl\t20071210\t593",
            "PC.607\tAACTGTGCGTAC\tFast\t20071112\t607",
            "PC.634\tACAGAGTCGGCT\tFast\t20080116\t634",
            "PC.635\tACCGCAGAGTCA\tFast\t20080116\t635",
            "PC.636\tACGGTGAGTGTC\tFast\t20080116\t636"]
        self.overview_map = MetadataMap(
            *parse_mapping_file_to_dict(self.overview_map_str))

        # Create the same overview tutorial map, but this time with some
        # comments.
        self.comment = "# Some comments about this mapping file"
        self.map_with_comments_str = self.overview_map_str[:]
        self.map_with_comments_str.insert(1, self.comment)
        self.map_with_comments = MetadataMap(*parse_mapping_file_to_dict(
            self.map_with_comments_str))

        # Create a MetadataMap object that has no metadata (i.e. no sample IDs,
        # so no metadata about samples).
        self.empty_map = MetadataMap({}, [])

        # Create a MetadataMap object that has samples (i.e. sample IDs) but
        # not associated metadata (i.e. no columns other than SampleID).
        self.no_metadata_str = ["#SampleID",
                                "PC.354",
                                "PC.355",
                                "PC.356",
                                "PC.481",
                                "PC.593",
                                "PC.607",
                                "PC.634",
                                "PC.635",
                                "PC.636"]
        self.no_metadata = MetadataMap(*parse_mapping_file_to_dict(
            self.no_metadata_str))

        # Create a MetadataMap object that has a category with only one value
        # throughout the entire column.
        self.single_value_str = ["#SampleID\tFoo",
                                 "PC.354\tfoo",
                                 "PC.355\tfoo",
                                 "PC.356\tfoo",
                                 "PC.481\tfoo",
                                 "PC.593\tfoo",
                                 "PC.607\tfoo",
                                 "PC.634\tfoo",
                                 "PC.635\tfoo",
                                 "PC.636\tfoo"]
        self.single_value = MetadataMap(*parse_mapping_file_to_dict(
            self.single_value_str))

        self.m1 = StringIO(m1)
        self.m2 = StringIO(m2)
        self.m3 = StringIO(m3)

        self.m1_dup_bad = StringIO(m1_dup_bad)
        self.m1_dup_good = StringIO(m1_dup_good)


    def test_parseMetadataMap(self):
        """Test parsing a mapping file into a MetadataMap instance."""
        obs = MetadataMap.parseMetadataMap(self.overview_map_str)
        self.assertEqual(obs, self.overview_map)

    def test_parseMetadataMap_empty(self):
        """Test parsing empty mapping file contents."""
        self.assertRaises(QiimeParseError, MetadataMap.parseMetadataMap, [])

    def test_eq(self):
        """Test whether two MetadataMaps are equal."""
        self.assertTrue(self.empty_map == MetadataMap({}, []))
        self.assertTrue(self.overview_map == MetadataMap(
            self.overview_map._metadata, self.overview_map.Comments))

    def test_ne(self):
        """Test whether two MetadataMaps are not equal."""
        self.assertTrue(self.empty_map != MetadataMap({}, ["foo"]))
        self.assertTrue(self.overview_map != MetadataMap(
            self.overview_map._metadata, ["foo"]))
        self.assertTrue(self.overview_map != MetadataMap({},
                                                         self.overview_map.Comments))
        self.assertTrue(self.overview_map != self.empty_map)
        self.assertTrue(self.overview_map != self.map_with_comments)
        self.assertTrue(self.overview_map != self.no_metadata)

    def test_getSampleMetadata(self):
        """Test metadata by sample ID accessor with valid sample IDs."""
        exp = {'BarcodeSequence': 'AGCACGAGCCTA', 'Treatment': 'Control',
               'DOB': '20061218', 'Description': '354'}
        obs = self.overview_map.getSampleMetadata('PC.354')
        self.assertEqual(obs, exp)

        exp = {'BarcodeSequence': 'ACCAGCGACTAG', 'Treatment': 'Control',
               'DOB': '20070314', 'Description': '481'}
        obs = self.map_with_comments.getSampleMetadata('PC.481')
        self.assertEqual(obs, exp)

        exp = {'BarcodeSequence': 'ACGGTGAGTGTC', 'Treatment': 'Fast',
               'DOB': '20080116', 'Description': '636'}
        obs = self.map_with_comments.getSampleMetadata('PC.636')
        self.assertEqual(obs, exp)

        exp = {}
        obs = self.no_metadata.getSampleMetadata('PC.636')
        self.assertEqual(obs, exp)

    def test_getSampleMetadata_bad_sample_id(self):
        """Test metadata by sample ID accessor with invalid sample IDs."""
        # Nonexistent sample ID.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata,
                          'PC.000')
        self.assertRaises(KeyError, self.no_metadata.getSampleMetadata,
                          'PC.000')
        # Integer sample ID.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata, 42)
        # Sample ID of type None.
        self.assertRaises(KeyError, self.overview_map.getSampleMetadata, None)

        # Sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, 's1')
        # Integer sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, 1)
        # Sample ID of None on empty map.
        self.assertRaises(KeyError, self.empty_map.getSampleMetadata, None)

    def test_getCategoryValue(self):
        """Test category value by sample ID/category name accessor."""
        exp = "Fast"
        obs = self.overview_map.getCategoryValue('PC.634', 'Treatment')
        self.assertEqual(obs, exp)

        exp = "20070314"
        obs = self.overview_map.getCategoryValue('PC.481', 'DOB')
        self.assertEqual(obs, exp)

        exp = "ACGGTGAGTGTC"
        obs = self.map_with_comments.getCategoryValue(
            'PC.636', 'BarcodeSequence')
        self.assertEqual(obs, exp)

    def test_getCategoryValues(self):
        """Test category value list by sample ID/category name accessor."""
        smpl_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
                    'PC.634', 'PC.635', 'PC.636']

        exp = [
            'Control',
            'Control',
            'Control',
            'Control',
            'Control',
            'Fast',
            'Fast',
            'Fast',
            'Fast']
        obs = self.overview_map.getCategoryValues(smpl_ids, 'Treatment')
        self.assertEqual(obs, exp)

    def test_isNumericCategory(self):
        """Test checking if a category is numeric."""
        obs = self.overview_map.isNumericCategory('Treatment')
        self.assertEqual(obs, False)

        obs = self.overview_map.isNumericCategory('DOB')
        self.assertEqual(obs, True)

    def test_hasUniqueCategoryValues(self):
        """Test checking if a category has unique values."""
        obs = self.overview_map.hasUniqueCategoryValues('Treatment')
        self.assertEqual(obs, False)

        obs = self.overview_map.hasUniqueCategoryValues('DOB')
        self.assertEqual(obs, False)

        obs = self.overview_map.hasUniqueCategoryValues('Description')
        self.assertEqual(obs, True)

    def test_hasSingleCategoryValue(self):
        """Test checking if a category has only a single value."""
        obs = self.overview_map.hasSingleCategoryValue('Treatment')
        self.assertEqual(obs, False)

        obs = self.single_value.hasSingleCategoryValue('Foo')
        self.assertEqual(obs, True)

    def test_getCategoryValue_bad_sample_id(self):
        """Test category value by sample ID accessor with bad sample IDs."""
        # Nonexistent sample ID.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
                          'PC.000', 'Treatment')
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
                          'PC.000', 'Treatment')
        # Integer sample ID.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue, 42,
                          'DOB')
        # Sample ID of type None.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue, None,
                          'Treatment')

        # Sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, 's1',
                          'foo')
        # Integer sample ID on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, 1,
                          'bar')
        # Sample ID of None on empty map.
        self.assertRaises(KeyError, self.empty_map.getCategoryValue, None,
                          'baz')

    def test_getCategoryValue_bad_category(self):
        """Test category value by sample ID accessor with bad categories."""
        # Nonexistent category.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
                          'PC.354', 'foo')
        # Integer category.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
                          'PC.354', 42)
        # Category of type None.
        self.assertRaises(KeyError, self.overview_map.getCategoryValue,
                          'PC.354', None)

        # Category on map with no metadata, but that has sample IDs.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
                          'PC.354', 'Treatment')
        # Integer category on map with no metadata.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
                          'PC.354', 34)
        # Category of type None on map with no metadata.
        self.assertRaises(KeyError, self.no_metadata.getCategoryValue,
                          'PC.354', None)

    def test_SampleIds(self):
        """Test sample IDs accessor."""
        exp = ["PC.354", "PC.355", "PC.356", "PC.481", "PC.593", "PC.607",
               "PC.634", "PC.635", "PC.636"]
        obs = self.overview_map.SampleIds
        self.assertEqual(obs, exp)

        obs = self.no_metadata.SampleIds
        self.assertEqual(obs, exp)

        obs = self.empty_map.SampleIds
        self.assertEqual(obs, [])

    def test_CategoryNames(self):
        """Test category names accessor."""
        exp = ["BarcodeSequence", "DOB", "Description", "Treatment"]
        obs = self.overview_map.CategoryNames
        self.assertEqual(obs, exp)

        obs = self.no_metadata.CategoryNames
        self.assertEqual(obs, [])

        obs = self.empty_map.CategoryNames
        self.assertEqual(obs, [])

    def test_filterSamples(self):
        """Test filtering out samples from metadata map."""
        exp = ['PC.356', 'PC.593']
        self.overview_map.filterSamples(['PC.593', 'PC.356'])
        obs = self.overview_map.SampleIds
        self.assertEqual(obs, exp)

        self.overview_map.filterSamples([])
        self.assertEqual(self.overview_map.SampleIds, [])

    def test_filterSamples_strict(self):
        """Test strict checking of sample prescence when filtering."""
        with self.assertRaises(ValueError):
            self.overview_map.filterSamples(['PC.356', 'abc123'])

        with self.assertRaises(ValueError):
            self.empty_map.filterSamples(['foo'])

    def test_filterSamples_no_strict(self):
        """Test missing samples does not raise error."""
        self.overview_map.filterSamples(['PC.356', 'abc123'], strict=False)
        self.assertEqual(self.overview_map.SampleIds, ['PC.356'])

        self.empty_map.filterSamples(['foo'], strict=False)
        self.assertEqual(self.empty_map.SampleIds, [])

    def test_str(self):
        """Test conversion to string representation
        """
        map_obj = MetadataMap.parseMetadataMap(self.m1)
        string_rep = StringIO(map_obj)

        self.m1.seek(0)

        # The string representation of the map_obj is unpredictable, since
        # it is generated from a dict, so we have to get a little clever
        exp_headers = self.m1.readline().strip().split('\t')
        obs_headers = string_rep.readline().strip().split('\t')

        # make sure that they have the same columns
        self.assertEqual(set(exp_headers), set(obs_headers))

        # make sure they have the same values for the same columns
        for obs, exp in zip(string_rep, self.m1):
            obs_elements = obs.strip().split('\t')
            exp_elements = exp.strip().split('\t')

            for exp_i, exp_header in enumerate(exp_headers):
                obs_i = obs_headers.index(exp_header)

                self.assertEqual(obs_elements[exp_i],
                                 exp_elements[obs_i])

    def test_merge_mapping_file(self):
        """merge_mapping_file: functions with default parameters
        """
        observed = MetadataMap.mergeMappingFiles([self.m1,self.m3])
        self.assertEqual(observed, m1_m3_exp)

    def test_merge_mapping_file_different_no_data_value(self):
        """merge_mapping_file: functions with different no_data_value
        """
        observed = MetadataMap.mergeMappingFiles([self.m1,self.m2],
                                                  "TESTING_NA")
        self.assertEqual(observed, m1_m2_exp)

    def test_merge_mapping_file_three_mapping_files(self):
        """merge_mapping_file: 3 mapping files
        """
        observed = MetadataMap.mergeMappingFiles([self.m1,self.m2,self.m3])
        self.assertEqual(observed, m1_m2_m3_exp)

    def test_merge_mapping_file_bad_duplicates(self):
        """merge_mapping_file: error raised when merging mapping files where
        same sample ids has different values
        """
        self.assertRaises(
            ValueError,
            MetadataMap.mergeMappingFiles,
            [self.m1, self.m1_dup_bad])

    def test_merge_mapping_file_good_duplicates(self):
        """merge_mapping_file: same sample ids merged correctly when they have
        mergable data
        """
        observed = MetadataMap.mergeMappingFiles([self.m1,self.m1_dup_good])
        self.assertEqual(observed, m1_m1_dup_good_exp)


class SyncBiomTests(TestCase):

    """Tests of sync_biom_and_mf."""

    def setUp(self):
        """Define data needed by all tests."""
        self.BT_IN_1 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        self.MF_IN_1 = ['#SampleIttest_cat\ttest_corr',
                        'Sample1\tcat1\t1',
                        'Sample2\tcat1\t2',
                        'Sample3\tcat2\t3',
                        'Sample4\tcat2\t4',
                        'Sample5\tcat3\t5',
                        'Sample6\tcat3\t6',
                        'NotInOtuTable1\tcat5\t7',
                        'NotInOtuTable2\tcat5\t8']
        self.MF_OUT_1 = ['#SampleIttest_cat\ttest_corr',
                         'Sample1\tcat1\t1',
                         'Sample2\tcat1\t2',
                         'Sample3\tcat2\t3',
                         'Sample4\tcat2\t4',
                         'Sample5\tcat3\t5',
                         'Sample6\tcat3\t6']
        self.BT_OUT_2 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T16:00:53.390212","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 4],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null}]}'
        self.MF_IN_2 = ['#SampleIttest_cat\ttest_corr',
                        'Sample1\tcat1\t1',
                        'Sample2\tcat1\t2',
                        'Sample3\tcat2\t3',
                        'Sample4\tcat2\t4']
        self.BT_OUT_3 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T16:09:55.670254","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 5],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,16.0],[0,4,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,48.0],[1,4,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,64.0],[2,4,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,33.0],[3,4,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,15.0],[4,4,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,29.0],[5,4,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        self.MF_IN_3 = ['#SampleIttest_cat\ttest_corr',
                        'Sample1\tcat1\t1',
                        'Sample2\tcat1\t2',
                        'NotInOtuTable2\tcat5\t8',
                        'Sample3\tcat2\t3',
                        'Sample5\tcat3\t5',
                        'Sample6\tcat3\t6',
                        'NotInOtuTable1\tcat5\t7']
        self.MF_OUT_3 = ['#SampleIttest_cat\ttest_corr',
                         'Sample1\tcat1\t1',
                         'Sample2\tcat1\t2',
                         'Sample3\tcat2\t3',
                         'Sample5\tcat3\t5',
                         'Sample6\tcat3\t6']
        self.BT_4 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

    def test_mapping_file_removal(self):
        """Test syncing removes samples from mf only."""
        mf_exp, _ = parse_mapping_file_to_dict(self.MF_OUT_1)
        # bt doesn't change in this test
        bt_exp = parse_biom_table(self.BT_IN_1)
        mf_obs, bt_obs, nonshared_samples = \
            sync_biom_and_mf(parse_mapping_file_to_dict(self.MF_IN_1)[0],
                             parse_biom_table(self.BT_IN_1))
        self.assertEqual(mf_exp, mf_obs)
        self.assertEqual(bt_exp, bt_obs)
        self.assertEqual(nonshared_samples, set(['NotInOtuTable2',
                                                 'NotInOtuTable1']))

    def test_mapping_file_removal(self):
        """Test syncing removes samples bt only."""
        # doesn't change in test
        mf_exp, _ = parse_mapping_file_to_dict(self.MF_IN_2)
        bt_exp = parse_biom_table(self.BT_OUT_2)
        mf_obs, bt_obs, nonshared_samples = \
            sync_biom_and_mf(parse_mapping_file_to_dict(self.MF_IN_2)[0],
                             parse_biom_table(self.BT_IN_1))
        self.assertEqual(mf_exp, mf_obs)
        self.assertEqual(bt_exp, bt_obs)
        self.assertEqual(nonshared_samples, set(['Sample5', 'Sample6']))

    def test_mapping_file_removal(self):
        """Test samples removed correctly when nonoverlapping."""
        mf_exp, _ = parse_mapping_file_to_dict(self.MF_OUT_3)
        # bt doesn't change in this test
        bt_exp = parse_biom_table(self.BT_OUT_3)
        mf_obs, bt_obs, nonshared_samples = \
            sync_biom_and_mf(parse_mapping_file_to_dict(self.MF_IN_3)[0],
                             parse_biom_table(self.BT_IN_1))
        self.assertEqual(mf_exp, mf_obs)
        self.assertEqual(bt_exp, bt_obs)
        self.assertEqual(nonshared_samples, set(['Sample4', 'NotInOtuTable2',
                                                 'NotInOtuTable1']))

    def test_biom_taxonomy_formatter(self):
        """Test that different metadata types have taxonomy returned."""
        # biom tables with different types of taxonomy metadata
        bt_list_taxonomy = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One","p__testCode"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        bt_str_taxonomy = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        bt_dict_taxonomy = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": {"k":"k__One"}}},{"id": "OTU2", "metadata": {"taxonomy": {"k":"k__Two"}}},{"id": "OTU3", "metadata": {"taxonomy": {"k":"k__Three"}}},{"id": "OTU4", "metadata": {"taxonomy": {"k":"k__Four"}}},{"id": "OTU5", "metadata": {"taxonomy": {"k":"k__Five"}}},{"id": "OTU6", "metadata": {"taxonomy": {"kingdom":"k__Six", "phylum":"TM7"}}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        bt_no_metadata = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": null},{"id": "OTU2", "metadata": null},{"id": "OTU3", "metadata": null},{"id": "OTU4", "metadata": null},{"id": "OTU5", "metadata": null},{"id": "OTU6", "metadata": null}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
        # test with list taxonomy
        bt = parse_biom_table(bt_list_taxonomy)
        obs = biom_taxonomy_formatter(bt, 'taxonomy')
        exp = ['k__One;p__testCode',
               'k__Two',
               'k__Three',
               'k__Four',
               'k__Five',
               'k__Six']
        self.assertEqual(obs, exp)
        # test with string taxonomy
        bt = parse_biom_table(bt_str_taxonomy)
        obs = biom_taxonomy_formatter(bt, 'taxonomy')
        exp = ['k__One', 'k__Two', 'k__Three', 'k__Four', 'k__Five', 'k__Six']
        self.assertEqual(obs, exp)
        # tests with dict taxonomy
        bt = parse_biom_table(bt_dict_taxonomy)
        obs = biom_taxonomy_formatter(bt, 'taxonomy')
        exp = ['k_k__One',
               'k_k__Two',
               'k_k__Three',
               'k_k__Four',
               'k_k__Five',
               'kingdom_k__Six phylum_TM7']
        self.assertEqual(obs, exp)
        # test that returns none when the taxonomy key is incorrect
        self.assertEqual(None, biom_taxonomy_formatter(bt, 'Nonexistent_MD'))
        # test performs correctly when no metadata in the biom table
        bt = parse_biom_table(bt_no_metadata)
        obs = biom_taxonomy_formatter(bt, 'taxonomy')
        exp = None
        self.assertEqual(obs, exp)


class RExecutorTests(TestCase):

    """Tests of the RExecutor class."""

    def setUp(self):
        """Define some useful test objects."""
        # The unweighted unifrac distance matrix from the overview tutorial.
        self.overview_dm_str = ["\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\
                                \tPC.607\tPC.634\tPC.635\tPC.636",
                                "PC.354\t0.0\t0.595483768391\t0.618074717633\
                                \t0.582763100909\t0.566949022108\
                                \t0.714717232268\t0.772001731764\
                                \t0.690237118413\t0.740681707488",
                                "PC.355\t0.595483768391\t0.0\t0.581427669668\
                                \t0.613726772383\t0.65945132763\
                                \t0.745176523638\t0.733836123821\
                                \t0.720305073505\t0.680785600439",
                                "PC.356\t0.618074717633\t0.581427669668\t0.0\
                                \t0.672149021573\t0.699416863323\
                                \t0.71405573754\t0.759178215168\
                                \t0.689701276341\t0.725100672826",
                                "PC.481\t0.582763100909\t0.613726772383\
                                \t0.672149021573\t0.0\t0.64756120797\
                                \t0.666018240373\t0.66532968784\
                                \t0.650464714994\t0.632524644216",
                                "PC.593\t0.566949022108\t0.65945132763\
                                \t0.699416863323\t0.64756120797\t0.0\
                                \t0.703720200713\t0.748240937349\
                                \t0.73416971958\t0.727154987937",
                                "PC.607\t0.714717232268\t0.745176523638\
                                \t0.71405573754\t0.666018240373\
                                \t0.703720200713\t0.0\t0.707316869557\
                                \t0.636288883818\t0.699880573956",
                                "PC.634\t0.772001731764\t0.733836123821\
                                \t0.759178215168\t0.66532968784\
                                \t0.748240937349\t0.707316869557\t0.0\
                                \t0.565875193399\t0.560605525642",
                                "PC.635\t0.690237118413\t0.720305073505\
                                \t0.689701276341\t0.650464714994\
                                \t0.73416971958\t0.636288883818\
                                \t0.565875193399\t0.0\t0.575788039321",
                                "PC.636\t0.740681707488\t0.680785600439\
                                \t0.725100672826\t0.632524644216\
                                \t0.727154987937\t0.699880573956\
                                \t0.560605525642\t0.575788039321\t0.0"]

        # The overview tutorial's metadata mapping file.
        self.overview_map_str = ["#SampleID\tBarcodeSequence\tTreatment\tDOB",
                                 "PC.354\tAGCACGAGCCTA\tControl\t20061218",
                                 "PC.355\tAACTCGTCGATG\tControl\t20061218",
                                 "PC.356\tACAGACCACTCA\tControl\t20061126",
                                 "PC.481\tACCAGCGACTAG\tControl\t20070314",
                                 "PC.593\tAGCAGCACTTGT\tControl\t20071210",
                                 "PC.607\tAACTGTGCGTAC\tFast\t20071112",
                                 "PC.634\tACAGAGTCGGCT\tFast\t20080116",
                                 "PC.635\tACCGCAGAGTCA\tFast\t20080116",
                                 "PC.636\tACGGTGAGTGTC\tFast\t20080116"]

        # The prefix to use for temporary files/dirs. This prefix may be added
        # to, but all temp dirs and files created by the tests will have this
        # prefix at a minimum.
        self.prefix = 'qiime_RExecutor_tests'

        self.start_dir = getcwd()
        self.dirs_to_remove = []
        self.files_to_remove = []

        self.tmp_dir = get_qiime_temp_dir()

        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)

            # If test creates the temp dir, also remove it.
            self.dirs_to_remove.append(self.tmp_dir)

        # Create temporary input dir/files.
        self.input_dir = mkdtemp(dir=self.tmp_dir,
                                 prefix='%s_input_dir_' % self.prefix)
        self.dirs_to_remove.append(self.input_dir)

        self.dm_fp = join(self.input_dir, 'dm.txt')
        dm_f = open(self.dm_fp, 'w')
        for line in self.overview_dm_str:
            dm_f.write(line + "\n")
        dm_f.close()
        self.files_to_remove.append(self.dm_fp)

        self.map_fp = join(self.input_dir, 'map.txt')
        map_f = open(self.map_fp, 'w')
        for line in self.overview_map_str:
            map_f.write(line + "\n")
        map_f.close()
        self.files_to_remove.append(self.map_fp)

        # Create temporary output directory.
        self.output_dir = mkdtemp(dir=self.tmp_dir,
                                  prefix='%s_output_dir_' % self.prefix)
        self.dirs_to_remove.append(self.output_dir)

    def tearDown(self):
        """Remove temporary files/dirs created by tests."""
        # Change back to the start dir.
        chdir(self.start_dir)
        remove_files(self.files_to_remove)

        # Remove directories last, so we don't get errors trying to remove
        # files which may be in the directories.
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_call(self):
        """Test executing an arbitrary command."""
        args = ["-d " + self.dm_fp + " -m " + self.map_fp +
                " -c DOB -o " + self.output_dir]

        rex = RExecutor(TmpDir=self.tmp_dir)
        results = rex(args, "permdisp.r", self.output_dir)

        self.files_to_remove.append(join(self.tmp_dir, 'R.stdout'))
        self.files_to_remove.append(join(self.tmp_dir, 'R.stderr'))

        # Make sure an output file was created with something in it.
        results_fp = join(self.output_dir, 'permdisp_results.txt')
        results_f = open(results_fp, 'U')
        results = results_f.read()
        results_f.close()
        self.files_to_remove.append(results_fp)
        self.assertTrue(len(results) > 0)


# Long strings of test data go here
fasta_lines = """>seq1
ACCAGCGGAGAC
>seq2
ACAGAGAGACCC
>seq3
ATTACCAGATTAC
>seq4
ACAGGAGACCGAGAAGA
>seq5
ACCAGAGACCGAGA
"""

expected_lines_50_perc = """>seq2
ACAGAGAGACCC
>seq4
ACAGGAGACCGAGAAGA
>seq5
ACCAGAGACCGAGA""".split('\n')

expected_lines_20_perc = """>seq4
ACAGGAGACCGAGAAGA""".split('\n')

m1="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tfirst1111\tTHE FIRST"""

# when merged with m1, this should raise a ValueError
m1_dup_bad="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tdupe11111\tDUPE BAD"""

# when merged with m1, this should NOT raise a ValueError
m1_dup_good="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional9\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tsomthing\tTHE FIRST"""

m2="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional2\tDescription
222222222\tGGGGGGGGGGGGGGG\tCCCCCCCCCCCCCCCCCCCC\tsecond222\tTHE SECOND"""

m3="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
333333333\tFFFFFFFFFFFFFFF\tIIIIIIIIIIIIIIIIIIII\tthird3333\tTHE THIRD"""

# the dict representation of the object with no_data_value set to TESTING_NA
m1_m2_exp = MetadataMap({'111111111': {'optional1': 'first1111', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'Description': 'THE FIRST', 'optional2': None}, '222222222': {'optional1': None, 'LinkerPrimerSequence': 'CCCCCCCCCCCCCCCCCCCC', 'BarcodeSequence': 'GGGGGGGGGGGGGGG', 'Description': 'THE SECOND', 'optional2': 'second222'}}, [])
m1_m2_exp.no_data_value = 'TESTING_NA'

# dict representation, with no_data_value left as default ("no_data")
m1_m3_exp = MetadataMap({'111111111': {'optional1': 'first1111', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'Description': 'THE FIRST'}, '333333333': {'optional1': 'third3333', 'LinkerPrimerSequence': 'IIIIIIIIIIIIIIIIIIII', 'BarcodeSequence': 'FFFFFFFFFFFFFFF', 'Description': 'THE THIRD'}}, [])
m1_m3_exp.no_data_value = 'no_data'

# dict representation, with no_data_value left as default ("no_data")
m1_m2_m3_exp = MetadataMap({'111111111': {'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'optional1': 'first1111', 'Description': 'THE FIRST', 'optional2': None}, '333333333': {'BarcodeSequence': 'FFFFFFFFFFFFFFF', 'LinkerPrimerSequence': 'IIIIIIIIIIIIIIIIIIII', 'optional1': 'third3333', 'Description': 'THE THIRD', 'optional2': None}, '222222222': {'BarcodeSequence': 'GGGGGGGGGGGGGGG', 'LinkerPrimerSequence': 'CCCCCCCCCCCCCCCCCCCC', 'optional1': None, 'Description': 'THE SECOND', 'optional2': 'second222'}}, [])
m1_m2_m3_exp.no_data_value = 'no_data'

# dict representation, m1 and m1_dup_bad should be mergeable
m1_m1_dup_good_exp = MetadataMap({'111111111': {'optional9': 'somthing', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'optional1': 'first1111', 'Description': 'THE FIRST', 'BarcodeSequence': 'AAAAAAAAAAAAAAA'}}, [])
m1_m1_dup_good_exp.no_data_value = 'no_data'


# run unit tests if run from command-line
if __name__ == '__main__':
    main()
