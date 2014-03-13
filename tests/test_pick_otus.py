#!/usr/bin/env python

"""Tests of code for OTU picking"""

__author__ = "Kyle Bittinger, Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = [
    "Kyle Bittinger",
    "Greg Caporaso",
    "Rob Knight",
    "Jens Reeder",
    "William Walters",
    "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import remove
from os.path import abspath, join, exists
from shutil import rmtree

from cogent.util.misc import create_dir
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from cogent.util.misc import remove_files
from cogent import DNA
from cogent.app.formatdb import build_blast_db_from_fasta_path

from qiime.util import get_tmp_filename, load_qiime_config, create_dir
from qiime.pick_otus import (CdHitOtuPicker, OtuPicker,
                             MothurOtuPicker, PrefixSuffixOtuPicker, TrieOtuPicker, BlastOtuPicker,
                             expand_otu_map_seq_ids, map_otu_map_files, UclustOtuPicker,
                             UclustReferenceOtuPicker, expand_failures, UsearchOtuPicker,
                             UsearchReferenceOtuPicker, get_blast_hits, BlastxOtuPicker,
                             Usearch610DeNovoOtuPicker, Usearch61ReferenceOtuPicker)


class OtuPickerTests(TestCase):

    """Tests of the abstract OtuPicker class"""

    def test_init(self):
        """Abstract OtuPicker __init__ should store name, params"""
        p = OtuPicker({})
        self.assertEqual(p.Name, 'OtuPicker')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract OtuPicker __call__ should raise NotImplementedError"""
        p = OtuPicker({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs')

    def test_prefilter_exact_matches(self):
        """Abstract OtuPicker _prefilter_exact_matches functions as expected
        """
        seqs = [('s1 comment1', 'ACCTTGTTACTTT'),  # three copies
                ('s2 comment2', 'ACCTTGTTACTTTC'),  # one copy
                ('s3 comment3', 'ACCTTGTTACTTTCC'),  # two copies
                ('s4 comment4', 'ACCTTGTTACTTT'),
                ('s5 comment5', 'ACCTTGTTACTTTCC'),
                ('s6 comment6', 'ACCTTGTTACTTT')]
        expected0 = [('QiimeExactMatch.s1', 'ACCTTGTTACTTT'),
                     ('QiimeExactMatch.s2', 'ACCTTGTTACTTTC'),
                     ('QiimeExactMatch.s3', 'ACCTTGTTACTTTCC')]
        expected1 = {'QiimeExactMatch.s1': ['s1', 's4', 's6'],
                     'QiimeExactMatch.s2': ['s2'],
                     'QiimeExactMatch.s3': ['s3', 's5']}
        expected = (expected0, expected1)
        p = OtuPicker({})
        actual = p._prefilter_exact_matches(seqs)
        self.assertEqual(actual, expected)


class MothurOtuPickerTests(TestCase):

    def setUp(self):
        self.small_seq_path = get_tmp_filename(
            prefix='MothurOtuPickerTest_', suffix='.fasta')
        f = open(self.small_seq_path, 'w')
        f.write(
            '>aaaaaa\nTAGGCTCTGATATAATAGCTCTC---------\n'
            '>cccccc\n------------TGACTACGCAT---------\n'
            '>bbbbbb\n----TATCGCTTCGACGATTCTCTGATAGAGA\n'
        )
        f.close()

    def tearDown(self):
        remove(self.small_seq_path)

    def test_call(self):
        app = MothurOtuPicker({})
        observed_otus = app(self.small_seq_path)
        expected_otus = [['cccccc'], ['bbbbbb'], ['aaaaaa']]
        assert_almost_equal(observed_otus.keys(),
                              [0, 1, 2])
        assert_almost_equal(observed_otus.values(),
                              expected_otus)

    def test_call_low_similarity(self):
        app = MothurOtuPicker({'Similarity': 0.35})
        observed_otus = app(self.small_seq_path)
        expected_otus = [['bbbbbb', 'cccccc'], ['aaaaaa']]
        assert_almost_equal(observed_otus.keys(),
                              [0, 1])
        assert_almost_equal(observed_otus.values(),
                              expected_otus)

    def test_call_nearest_neighbor(self):
        app = MothurOtuPicker({'Algorithm': 'nearest', 'Similarity': 0.35})
        observed_otus = app(self.small_seq_path)
        expected_otus = [['bbbbbb', 'cccccc'], ['aaaaaa']]
        assert_almost_equal(observed_otus.keys(),
                              [0, 1])
        assert_almost_equal(observed_otus.values(),
                              expected_otus)


class BlastxOtuPickerTests(TestCase):

    """ Tests of the blastx-based otu picker """

    def setUp(self):
        """
        """
        self.otu_picker = BlastxOtuPicker({'max_e_value': 0.001})
        self.seqs = [
            ('s0  some description', 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'),
            ('s1', 'TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),
            ('s2', 'TGCAGCTTGAGCCACAGGAGAGAGCCTTC'),
            ('s3', 'TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),
            ('s4', 'ACCGATGAGATATTAGCACAGGGGAATTAGAACCA'),
            ('s5', 'TGTCGAGAGTGAGATGAGATGAGAACA'),
            ('s6', 'ACGTATTTTAATTTGGCATGGT'),
            ('s7', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
        ]

        self.ref_seqs_pr = [
            ('ref1', 'CSLSHRRERA'),
            ('ref2', 'TDEILAQGN'),
            ('ref3', 'CRE'),
            ('ref4', 'TYFNGAW'),
            ('ref5', 'RATGEREL'),
        ]

        self.seqs_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')
        self.reference_seqs_pr_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')

        f = open(self.seqs_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.seqs]))
        f.close()

        f = open(self.reference_seqs_pr_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.ref_seqs_pr]))
        f.close()

        self.blast_db_pr, self.pr_db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_pr_fp,
                                           is_protein=True)

        self._files_to_remove = self.pr_db_files_to_remove +\
            [self.seqs_fp,
             self.reference_seqs_pr_fp]

    def tearDown(self):
        """
        """
        remove_files(self._files_to_remove, error_on_missing=False)

    def test_get_blast_hits_blastx(self):
        """get_blast_hits functions as expected with blastx """

        actual = get_blast_hits(
            self.seqs,
            self.blast_db_pr,
            max_e_value=0.01,
            min_pct_identity=0.5,
            min_aligned_percent=0.5,
            blast_program='blastx')

        # couple of sanity checks against command line blast
        self.assertEqual(len(actual['s3']), 2)
        self.assertEqual(actual['s3'][0]['SUBJECT ID'], 'ref1')
        self.assertEqual(actual['s3'][1]['SUBJECT ID'], 'ref5')

        # increase stringency reduces number of blast hits
        actual = get_blast_hits(
            self.seqs,
            self.blast_db_pr,
            max_e_value=0.001,
            min_pct_identity=0.5,
            min_aligned_percent=0.5,
            blast_program='blastx')
        # couple of sanity checks against command line blast
        self.assertEqual(len(actual['s3']), 1)
        self.assertEqual(actual['s3'][0]['SUBJECT ID'], 'ref1')

    def test_call(self):
        """BLASTX OTU Picker functions as expected
        """

        expected = {'ref1': ['s3', 's2', 's1'],
                    'ref2': ['s4']}
        actual = self.otu_picker(self.seqs_fp,
                                 refseqs_fp=self.reference_seqs_pr_fp)
        self.assertEqual(actual, expected)


class BlastOtuPickerTests(TestCase):

    """ Tests of the blast-based otu picker """

    def setUp(self):
        """
        """
        self.otu_picker = BlastOtuPicker({'max_e_value': 1e-3})
        self.seqs = [
            ('s0  some description', 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'),
            ('s1', 'TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),
            ('s2', 'TGCAGCTTGAGCCACAGGAGAGAGCCTTC'),
            ('s3', 'TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),
            ('s4', 'ACCGATGAGATATTAGCACAGGGGAATTAGAACCA'),
            ('s5', 'TGTCGAGAGTGAGATGAGATGAGAACA'),
            ('s6', 'ACGTATTTTAATTTGGCATGGT'),
            ('s7', 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),
        ]

        self.ref_seqs = [
            ('ref1', 'TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),
            ('ref2', 'ACCGATGAGATATTAGCACAGGGGAATTAGAACCA'),
            ('ref3', 'TGTCGAGAGTGAGATGAGATGAGAACA'),
            ('ref4', 'ACGTATTTTAATGGGGCATGGT'),
            ('ref5', 'AGAGCCACAGGAGAGAGAGAGCTTC'),
        ]

        self.ref_seqs_rc = [
            ('ref1', DNA.rc('TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC')),
            ('ref2', DNA.rc('ACCGATGAGATATTAGCACAGGGGAATTAGAACCA')),
            ('ref3', DNA.rc('TGTCGAGAGTGAGATGAGATGAGAACA')),
            ('ref4', DNA.rc('ACGTATTTTAATGGGGCATGGT')),
        ]

        self.seqs_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')
        self.reference_seqs_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')
        self.reference_seqs_rc_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')

        f = open(self.seqs_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.seqs]))
        f.close()

        f = open(self.reference_seqs_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.ref_seqs]))
        f.close()

        f = open(self.reference_seqs_rc_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.ref_seqs_rc]))
        f.close()

        self.blast_db, self.db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_fp)

        self._files_to_remove = self.db_files_to_remove +\
            [self.seqs_fp,
             self.reference_seqs_fp,
             self.reference_seqs_rc_fp]

    def tearDown(self):
        """
        """
        remove_files(self._files_to_remove, error_on_missing=False)

    def test_blast_seqs(self):
        """ blast_seqs: functions as expected
        """
        blast_db, db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._files_to_remove += db_files_to_remove
        self.otu_picker.blast_db = blast_db

        actual_clusters, actual_failures =\
            self.otu_picker._blast_seqs(self.seqs)

        for v in actual_clusters.values():
            v.sort()
        actual_failures.sort()

        expected_clusters = {'ref1': ['s1', 's2', 's3'], 'ref2': ['s4'],
                             'ref3': ['s5']}
        expected_failures = ['s0', 's6', 's7']

        self.assertEqual(actual_clusters, expected_clusters)
        self.assertEqual(actual_failures, expected_failures)

    def test_update_cluster_map(self):
        """update_cluster_map: functions as expected
        """
        # nothing in original cm
        cm = {}
        new_cm = {'c1': ['1', '2', '5'], 'c2': ['4', '3']}
        expected = new_cm
        actual = self.otu_picker._update_cluster_map(cm, new_cm)
        self.assertEqual(actual, expected)

        # no new clusters
        cm = {'c1': ['1', '2', '5'], 'c2': ['4', '3']}
        new_cm = {}
        expected = cm
        actual = self.otu_picker._update_cluster_map(cm, new_cm)
        self.assertEqual(actual, expected)

        # overlapping clusters
        cm = {'c1': ['1', '2', '5'], 'c2': ['4', '3']}
        new_cm = {'c1': ['8'], 'c2': ['10', '14'], '3': ['42']}
        expected = {'c1': ['1', '2', '5', '8'], 'c2':
                    ['4', '3', '10', '14'], '3': ['42']}
        actual = self.otu_picker._update_cluster_map(cm, new_cm)
        self.assertEqual(actual, expected)

        # no duplicate seq_id checking
        cm = {'c1': ['1']}
        new_cm = cm
        expected = {'c1': ['1', '1']}
        actual = self.otu_picker._update_cluster_map(cm, new_cm)
        self.assertEqual(actual, expected)

        # no clusters at all
        actual = self.otu_picker._update_cluster_map({}, {})
        self.assertEqual(actual, {})

    def test_get_blast_hits_blastn(self):
        """get_blast_hits functions as expected with blastn """

        actual = get_blast_hits(
            self.seqs,
            self.blast_db,
            max_e_value=1e-10,
            min_pct_identity=0.5,
            min_aligned_percent=0.5)
        # couple of sanity checks against command line blast
        self.assertEqual(len(actual['s3']), 2)
        self.assertEqual(actual['s3'][0]['SUBJECT ID'], 'ref1')
        self.assertEqual(actual['s3'][1]['SUBJECT ID'], 'ref5')

        # increase stringency reduces number of blast hits
        actual = get_blast_hits(
            self.seqs,
            self.blast_db,
            max_e_value=1e-10,
            min_pct_identity=0.5,
            min_aligned_percent=0.8)
        # couple of sanity checks against command line blast
        self.assertEqual(len(actual['s3']), 1)
        self.assertEqual(actual['s3'][0]['SUBJECT ID'], 'ref1')

    def test_call(self):
        """BLAST OTU Picker functions as expected
        """

        expected = {'ref1': ['s3', 's2', 's1'],
                    'ref2': ['s4'],
                    'ref3': ['s5']}
        actual = self.otu_picker(self.seqs_fp,
                                 refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual, expected)

    def test_call_alt_min_aligned_length(self):
        """BLAST OTU picker handles alt min_aligned_percent values """
        # first 12 bases match perfect, and no alignment from there
        seqs = [('s1', 'TGCAGCTTGAGCGTTGTTACCGCTTT')]
        ref_seqs = [
            ('r1', 'TGCAGCTTGAGCCACGCCGAATAGCCGAGTTTGACCGGGCCCAGGAGGAGAGAGAGAGCTTC')]

        seqs_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')
        reference_seqs_fp = get_tmp_filename(
            prefix='BlastOtuPickerTest_', suffix='.fasta')

        f = open(seqs_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in seqs]))
        f.close()

        f = open(reference_seqs_fp, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in ref_seqs]))
        f.close()

        self._files_to_remove.append(seqs_fp)
        self._files_to_remove.append(reference_seqs_fp)

        # with low min_aligned_percent s1 matches r1
        otu_picker = BlastOtuPicker({'max_e_value': 1e-3,
                                     'min_aligned_percent': 0.10})
        expected = {'r1': ['s1']}
        actual = otu_picker(seqs_fp,
                            refseqs_fp=reference_seqs_fp)
        self.assertEqual(actual, expected)

        # with min_aligned_percent s1 doesn't match r1
        otu_picker = BlastOtuPicker({'max_e_value': 1e-3,
                                     'min_aligned_percent': 0.50})
        expected = {}
        actual = otu_picker(seqs_fp,
                            refseqs_fp=reference_seqs_fp)
        self.assertEqual(actual, expected)

    def test_call_rc(self):
        """BLAST OTU picker: RC seqs cluster to same OTU as forward orientation
        """

        expected = {'ref1': ['s3', 's2', 's1'],
                    'ref2': ['s4'],
                    'ref3': ['s5']}
        actual = self.otu_picker(self.seqs_fp,
                                 refseqs_fp=self.reference_seqs_rc_fp)
        self.assertEqual(actual, expected)

    def test_call_alt_params(self):
        """BLAST OTU Picker functions as expected with alt params
        """
        otu_picker = BlastOtuPicker({'max_e_value': 1e-30})
        expected = {}
        actual = otu_picker(self.seqs_fp,
                            refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual, expected)

        self.otu_picker = BlastOtuPicker(
            {'max_e_value': 1e-3, 'Similarity': 0.90})
        expected_90 = {'ref1': ['s3', 's2', 's1'],
                       'ref2': ['s4'],
                       'ref3': ['s5'],
                       'ref4': ['s6']}
        actual = self.otu_picker(self.seqs_fp,
                                 refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual, expected_90)

    def test_call_preexisting_blast_db(self):
        """BLAST OTU Picker functions w preexisting blast db
        """
        blast_db, db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._files_to_remove += db_files_to_remove
        expected = {'ref1': ['s3', 's2', 's1'],
                    'ref2': ['s4'],
                    'ref3': ['s5']}
        actual = self.otu_picker(self.seqs_fp, blast_db=blast_db)
        self.assertEqual(actual, expected)

    def test_call_multiple_blast_runs(self):
        """BLAST OTU Picker not affected by alt SeqsPerBlastRun
        """
        expected = {'ref1': ['s1', 's2', 's3'],
                    'ref2': ['s4'],
                    'ref3': ['s5']}
        for v in expected.values():
            v.sort()
        for SeqsPerBlastRun in [1, 2, 4, 6, 7, 8, 100]:
            self.otu_picker.Params['seqs_per_blast_run'] \
                = SeqsPerBlastRun
            actual = self.otu_picker(self.seqs_fp,
                                     refseqs_fp=self.reference_seqs_fp)
            for v in actual.values():
                v.sort()
            self.assertEqual(actual, expected)


class PrefixSuffixOtuPickerTests(TestCase):

    """ Tests of the prefix/suffix-based OTU picker """

    def setUp(self):
        """
        """
        self.otu_picker = PrefixSuffixOtuPicker({})
        self.seqs = [
            ('s1  some description', 'ACGTAATGGT'),
            ('s2', 'ATTTAATGGT'),
            ('s3', 'ACGTAATTTT'),
            ('s4', 'AAATAAAAA'),
            ('s5', 'ACGTTGGT'),
            ('s6', 'ACGTATTTTAATTTGGCATGGT'),
        ]

        self.small_seq_path = get_tmp_filename(
            prefix='PrefixSuffixOtuPickerTest_', suffix='.fasta')
        self._files_to_remove = [self.small_seq_path]
        f = open(self.small_seq_path, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in self.seqs]))
        f.close()

    def tearDown(self):
        """
        """
        remove_files(self._files_to_remove)

    def test_call(self):
        """Prefix/suffix OTU Picker functions as expected
        """
        expected = {3: ['s1', 's5', 's6'],
                    1: ['s2'],
                    0: ['s3'],
                    2: ['s4']}
        actual = self.otu_picker(self.small_seq_path,
                                 prefix_length=4, suffix_length=4)
        self.assertEqual(actual, expected)

    def test_call_extra_long_lengths(self):
        """Prefix/suffix OTU Picker functions as expected
        """
        seqs = [
            ('s1  some description',
             'ACGTAATGGTCCCCCCCCCGGGGGGGGCCCCCCGGG'),
            ('s2', 'ATTTAATGGT'),
            ('s3', 'ACGTAATTTT'),
            ('s4', 'AAATAAAAA'),
            ('s5', 'ACGTTGGT'),
            ('s6', 'ACGTATTTTAATTTGGCATGGT'),
            ('s7', 'ACGTATTTTAATTTGGCATGG'),
            ('s1_dup',
             'ACGTAATGGTCCCCCCCCCGGGGGGGGCCCCCCGGG'),
            ('s2_dup', 'ATTTAATGGT'),
        ]
        seq_path = get_tmp_filename(
            prefix='PrefixSuffixOtuPickerTest_', suffix='.fasta')
        self._files_to_remove.append(seq_path)
        f = open(seq_path, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in seqs]))
        f.close()
        expected = {1: ['s1', 's1_dup'],
                    2: ['s2', 's2_dup'],
                    3: ['s3'],
                    4: ['s4'],
                    5: ['s5'],
                    6: ['s6'],
                    7: ['s7']}

        # long prefix collapses identical sequences
        actual = self.otu_picker(seq_path,
                                 prefix_length=400, suffix_length=0)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        assert_almost_equal(actual_clusters, expected_clusters)

        # long suffixes collapses identical sequences
        actual = self.otu_picker(seq_path,
                                 prefix_length=0, suffix_length=400)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        assert_almost_equal(actual_clusters, expected_clusters)

        # long prefix and suffixes collapses identical sequences
        actual = self.otu_picker(seq_path,
                                 prefix_length=400, suffix_length=400)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        assert_almost_equal(actual_clusters, expected_clusters)

    def test_collapse_exact_matches_prefix_and_suffix(self):
        """Prefix/suffix: collapse_exact_matches fns with pref/suf len > 0
        """
        expected = [['s1', 's5', 's6'], ['s2'], ['s3'], ['s4']]
        actual = sorted(
            self.otu_picker._collapse_exact_matches(self.seqs, 4, 4))
        expected.sort()
        self.assertEqual(actual, expected)

        expected = [['s1', 's2', 's3', 's5', 's6'], ['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs, 1, 1)
        actual.sort()
        expected.sort()
        self.assertEqual(actual, expected)

    def test_collapse_exact_matches_prefix_zero(self):
        """Prefix/suffix: collapse_exact_matches fns with prefix len = 0
        """
        expected = [['s1', 's2', 's5', 's6'], ['s3'], ['s4']]
        actual = sorted(
            self.otu_picker._collapse_exact_matches(self.seqs, 0, 4))
        expected.sort()
        self.assertEqual(actual, expected)

        expected = [['s1', 's2', 's3', 's5', 's6'], ['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs, 0, 1)
        actual.sort()
        expected.sort()
        self.assertEqual(actual, expected)

    def test_collapse_exact_matches_suffix_zero(self):
        """Prefix/suffix: collapse_exact_matches fns with suffix len = 0
        """
        expected = [['s1', 's3', 's5', 's6'], ['s2'], ['s4']]
        actual = sorted(
            self.otu_picker._collapse_exact_matches(self.seqs, 4, 0))
        expected.sort()
        self.assertEqual(actual, expected)

        expected = [['s1', 's2', 's3', 's4', 's5', 's6']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs, 1, 0)
        actual.sort()
        expected.sort()
        self.assertEqual(actual, expected)

    def test_build_seq_hash(self):
        """ """
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 0, 0), '')

        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 2, 2), 'ATGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATTACGT', 2, 1), 'ATT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 1, 2), 'AGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 1, 1), 'AT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 4, 3), 'ATGTCGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 3, 4), 'ATGACGT')

        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 4, 4), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 5, 3), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 8, 0), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 3, 5), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 0, 8), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 4, 5), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 5, 4), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 300, 0), 'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(
            'ATGTACGT', 0, 300), 'ATGTACGT')


class TrieOtuPickerTests(TestCase):

    """ Tests of the Trie-based OTU picker """

    def setUp(self):
        """
        """
        self.otu_picker = TrieOtuPicker({})
        self.otu_picker_rev = TrieOtuPicker({'Reverse': True})
        seqs = [
            ('s1 some description', 'ACGTAATGGT'),
            ('s2', 'ACGTATTTTAATTTGGCATGGT'),
            ('s3', 'ACGTAAT'),
            ('s4', 'ACGTA'),
            ('s5', 'ATTTAATGGT'),
            ('s6', 'ATTTAAT'),
            ('s7', 'AAATAAAAA')
        ]
        seqs_rev = [
            ('s1 some description', 'TGGTAATGCA'),
            ('s2', 'TGGTACGGTTTAATTTTATGCA'),
            ('s3', 'TAATGCA'),
            ('s4', 'ATGCA'),
            ('s5', 'TGGTAATTTA'),
            ('s6', 'TAATTTA'),
            ('s7', 'AAAAATAAA')
        ]

        self.small_seq_path = get_tmp_filename(
            prefix='TrieOtuPickerTest_', suffix='.fasta')
        self._files_to_remove = [self.small_seq_path]
        f = open(self.small_seq_path, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in seqs]))
        f.close()

        self.small_seq_path_rev = get_tmp_filename(
            prefix='TrieOtuPickerTest_', suffix='.fasta')
        self._files_to_remove.append(self.small_seq_path_rev)
        f = open(self.small_seq_path_rev, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in seqs_rev]))
        f.close()

    def tearDown(self):
        """
        """
        remove_files(self._files_to_remove)

    def test_call(self):
        """Trie OTU Picker functions as expected
        """
        expected = {0: ['s2'],
                    1: ['s3', 's4', 's1'],
                    2: ['s7'],
                    3: ['s6', 's5']}
        actual = self.otu_picker(self.small_seq_path)
        self.assertEqual(actual, expected)

    def test_call_reverse(self):
        """Trie OTU Picker functions as expected with the 'Reverse' option
        """
        expected = {0: ['s2'],
                    1: ['s3', 's4', 's1'],
                    2: ['s7'],
                    3: ['s6', 's5']}
        actual = self.otu_picker_rev(self.small_seq_path_rev)
        self.assertEqual(actual, expected)


class Usearch610DeNovoOtuPickerTests(TestCase):

    """ Tests for usearch 6.1 de novo functionality """

    def setUp(self):
        # create the temporary input files

        self.output_dir = load_qiime_config()['temp_dir']

        self.dna_seqs_usearch_97perc_id = dna_seqs_usearch_97perc_id
        self.dna_seqs_usearch_97perc_id_rc = dna_seqs_usearch_97perc_id_rc
        self.dna_seqs_usearch_97perc_id_len_diff =\
            dna_seqs_usearch_97perc_id_len_diff
        self.dna_seqs_usearch_97perc_dups = dna_seqs_usearch_97perc_dups

        self.tmp_seq_filepath_97perc_id = get_tmp_filename(
            prefix='Usearch610DeNovoOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath_97perc_id, 'w')
        seq_file.write(self.dna_seqs_usearch_97perc_id)
        seq_file.close()

        self.tmp_seq_filepath_97perc_id_rc = get_tmp_filename(
            prefix='Usearch610DeNovoOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath_97perc_id_rc, 'w')
        seq_file.write(self.dna_seqs_usearch_97perc_id_rc)
        seq_file.close()

        self.tmp_seqs_usearch97perc_id_len_diff = get_tmp_filename(
            prefix="Usearch610DeNovoOtuPickerTest_",
            suffix=".fasta")
        seq_file = open(self.tmp_seqs_usearch97perc_id_len_diff, "w")
        seq_file.write(self.dna_seqs_usearch_97perc_id_len_diff)
        seq_file.close()

        self.tmp_seqs_usearch_97perc_dups = get_tmp_filename(
            prefix="Usearch610DeNovoOtuPickerTest_",
            suffix=".fasta")
        seq_file = open(self.tmp_seqs_usearch_97perc_dups, "w")
        seq_file.write(self.dna_seqs_usearch_97perc_dups)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath_97perc_id, self.tmp_seq_filepath_97perc_id_rc,
             self.tmp_seqs_usearch97perc_id_len_diff,
             self.tmp_seqs_usearch_97perc_dups]

        self._dirs_to_remove = []

    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)

    def test_call_default_params(self):
        """ clusters seqs within 97% identity with default parameters """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        # All seqs should fall into a single cluster
        expected_clusters = {'denovo0': ['usearch_ecoli_seq',
                                         'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq_1bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

    def test_call_default_params_and_lower_id(self):
        """ clusters seqs within 95% identity with default parameters """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'percent_id': 0.95
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        # All seqs should fall into a single cluster
        expected_clusters = {'denovo0': ['usearch_ecoli_seq',
                                         'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq_1bp_change']}

        assert_almost_equal(obs_clusters, expected_clusters)

    def test_call_default_params_and_higher_id(self):
        """ clusters seqs within 99% identity with default parameters """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'percent_id': 0.99
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        # Seqs should fall into separate clusters
        expected_clusters = {'denovo0': ['usearch_ecoli_seq'],
                             'denovo1': ['usearch_ecoli_seq_2bp_change'],
                             'denovo2': ['usearch_ecoli_seq_1bp_change']}

        # should be exactly 3 clusters
        self.assertEqual(len(obs_clusters), 3)
        assert_almost_equal(obs_clusters.keys(), expected_clusters.keys())
        assert_almost_equal(
            obs_clusters.values(),
            expected_clusters.values())

    def test_call_default_params_reversed_seq(self):
        """ Does not cluster reverse complemented sequence without --rev """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id_rc)

        # RC seq should fall into its own cluster
        expected_clusters = [set(['usearch_ecoli_seq',
                                  'usearch_ecoli_seq_1bp_change']),
                             set(['usearch_ecoli_seq_2bp_change_rc'])]

        self.assertEqual(len(obs_clusters), 2)
        for result in obs_clusters:
            self.assertTrue(set(obs_clusters[result]) in expected_clusters)

    def test_call_default_params_reversed_seq_w_rev(self):
        """ Does not cluster reverse complemented sequence without --rev """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'rev': True
                    })

        obs_clusters =\
            set(app(self.tmp_seq_filepath_97perc_id_rc).values()[0])

        # All seqs should fall into a single cluster
        expected_clusters = set(['usearch_ecoli_seq',
                                 'usearch_ecoli_seq_2bp_change_rc', 'usearch_ecoli_seq_1bp_change'])

        self.assertEqual(obs_clusters, expected_clusters)

    def test_call_default_params_save_intermediate_files(self):
        """ Preserves files if save_intermediate_files/logs is True """

        intermediate_files_dir = self.output_dir + "/test_usearch61/"
        create_dir(intermediate_files_dir)
        self._dirs_to_remove.append(intermediate_files_dir)

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': True,
                    'output_dir':
                    intermediate_files_dir,
                    'remove_usearch_logs': False
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        expected_intermediate_fps =\
            [intermediate_files_dir + "denovo_abundance_sorted.fna",
             intermediate_files_dir + "denovo_abundance_sorted.uc",
             intermediate_files_dir + "denovo_smallmem_clustered.uc",
             intermediate_files_dir + "abundance_sorted.log",
             intermediate_files_dir + "smallmem_clustered.log"]

        for curr_file in expected_intermediate_fps:
            self.assertTrue(exists(curr_file))

    def test_call_default_params_save_intermediate_files_fast_cluster(self):
        """ Preserves files if save_intermediate_files/logs is True """

        intermediate_files_dir = self.output_dir + "/test_usearch61_fast/"
        create_dir(intermediate_files_dir)
        self._dirs_to_remove.append(intermediate_files_dir)

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': True,
                    'output_dir':
                    intermediate_files_dir,
                    'remove_usearch_logs': False,
                    'usearch61_sort_method':
                    'length',
                    'usearch_fast_cluster': True
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        expected_intermediate_fps =\
            [intermediate_files_dir + "denovo_fast_clustered.uc",
             intermediate_files_dir + "fast_clustered.log"]

        for curr_file in expected_intermediate_fps:
            self.assertTrue(exists(curr_file))

        # All seqs should fall into a single cluster
        expected_clusters = {'denovo0': ['usearch_ecoli_seq',
                                         'usearch_ecoli_seq_1bp_change', 'usearch_ecoli_seq_2bp_change']}

        assert_almost_equal(obs_clusters, expected_clusters)

    def test_call_default_params_minlen(self):
        """ Discards reads that fall below minlen setting """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'minlen': 101
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id)

        # Should get no results
        expected_clusters = {}

        self.assertEqual(obs_clusters, expected_clusters)

    def test_usearch61_params(self):
        """ usearch61 handles changes to other parameters """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'wordlength': 25,
                    'usearch61_maxrejects': 200,
                    'usearch61_maxaccepts': 5
                    })

        obs_clusters = app(self.tmp_seq_filepath_97perc_id, otu_prefix="test")

        # All seqs should fall into a single cluster
        expected_clusters = {'test0': ['usearch_ecoli_seq',
                                       'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq_1bp_change']}

        assert_almost_equal(obs_clusters, expected_clusters)

    def test_usearch61_length_sorting(self):
        """ Sorting according to length, clusters seqs """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'usearch61_sort_method':
                    'length'
                    })

        obs_clusters = app(self.tmp_seqs_usearch97perc_id_len_diff)

        # All seqs should fall into a single cluster
        expected_clusters = {'denovo0': ['usearch_ecoli_seq_2bp_change',
                                         'usearch_ecoli_seq_1bp_change', 'usearch_ecoli_seq']}

        self.assertEqual(obs_clusters, expected_clusters)

    def test_usearch61_sizeorder(self):
        """ Handles sizeorder option, clusters seqs """

        app = Usearch610DeNovoOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir': self.output_dir,
                    'remove_usearch_logs': True,
                    'sizeorder': True
                    })

        obs_clusters = app(self.tmp_seqs_usearch_97perc_dups)

        # All seqs should fall into a single cluster
        expected_clusters = {'denovo0': ['usearch_ecoli_seq_1bp_change',
                                         'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq',
                                         'usearch_ecoli_seq_1bp_change_dup1',
                                         'usearch_ecoli_seq_1bp_change_dup2']}

        self.assertEqual(obs_clusters, expected_clusters)


class Usearch61ReferenceOtuPickerTests(TestCase):

    """ Tests for usearch 6.1 reference functionality """

    def setUp(self):
        # create the temporary input files

        self.output_dir = load_qiime_config()['temp_dir']

        self.dna_seqs_usearch_97perc_id = dna_seqs_usearch_97perc_id
        self.dna_seqs_usearch_97perc_id_rc = dna_seqs_usearch_97perc_id_rc
        self.dna_seqs_usearch_97perc_id_len_diff =\
            dna_seqs_usearch_97perc_id_len_diff
        self.dna_seqs_usearch_97perc_dups = dna_seqs_usearch_97perc_dups
        self.dna_seqs_rc_single_seq = dna_seqs_rc_single_seq

        self.tmp_seq_filepath_97perc_id = get_tmp_filename(
            prefix='Usearch610DeNovoOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath_97perc_id, 'w')
        seq_file.write(self.dna_seqs_usearch_97perc_id)
        seq_file.close()

        self.tmp_seq_filepath_97perc_id_rc = get_tmp_filename(
            prefix='Usearch610DeNovoOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath_97perc_id_rc, 'w')
        seq_file.write(self.dna_seqs_usearch_97perc_id_rc)
        seq_file.close()

        self.tmp_seqs_usearch97perc_id_len_diff = get_tmp_filename(
            prefix="Usearch610DeNovoOtuPickerTest_",
            suffix=".fasta")
        seq_file = open(self.tmp_seqs_usearch97perc_id_len_diff, "w")
        seq_file.write(self.dna_seqs_usearch_97perc_id_len_diff)
        seq_file.close()

        self.tmp_seqs_usearch_97perc_dups = get_tmp_filename(
            prefix="Usearch610DeNovoOtuPickerTest_",
            suffix=".fasta")
        seq_file = open(self.tmp_seqs_usearch_97perc_dups, "w")
        seq_file.write(self.dna_seqs_usearch_97perc_dups)
        seq_file.close()

        self.tmp_seqs_rc_single_seq = get_tmp_filename(
            prefix="Usearch610DeNovoOtuPickerTest_",
            suffix=".fasta")
        seq_file = open(self.tmp_seqs_rc_single_seq, "w")
        seq_file.write(self.dna_seqs_rc_single_seq)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath_97perc_id, self.tmp_seq_filepath_97perc_id_rc,
             self.tmp_seqs_usearch97perc_id_len_diff,
             self.tmp_seqs_usearch_97perc_dups, self.tmp_seqs_rc_single_seq]

        self._dirs_to_remove = []

    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)

    def test_call_default_params(self):
        """ clusters seqs within 97% identity with default parameters """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        # Randomly selected match is used for equivalent matches, so need to
        # test for results without order affecting output
        expected_clusters =\
            {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change',
                                              'usearch_ecoli_seq_2bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_default_params_and_lower_id(self):
        """ clusters seqs within 95% identity with default parameters """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'percent_id': 0.95
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        expected_clusters = {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change',
                                                              'usearch_ecoli_seq_2bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_default_params_and_higher_id(self):
        """ clusters seqs within 99% identity with default parameters """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'percent_id': 0.99
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        # Seqs should fall into separate clusters
        expected_clusters = {'denovo0': ['usearch_ecoli_seq_2bp_change'],
                             'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change']}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_default_params_reversed_seq(self):
        """ Does not cluster reverse complemented sequence without --rev """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'suppress_new_clusters':
                    True,
                    'rev': False
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seqs_rc_single_seq)

        # As seqs are not in same frame, should all fail.
        expected_clusters =\
            {}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = ['usearch_ecoli_seq',
                             'usearch_ecoli_seq_2bp_change',
                             'usearch_ecoli_seq_1bp_change']
        self.assertEqual(len(failures), 3)
        for curr_failure in failures:
            self.assertTrue(curr_failure in expected_failures)

    def test_call_default_params_reversed_seq_w_rev(self):
        """ Does not cluster reverse complemented sequence without --rev """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'rev': True,
                    'suppress_new_clusters': True
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        # All seqs should fall into a single cluster
        expected_clusters =\
            {'usearch_ecoli_seq_2bp_change_rc': ['usearch_ecoli_seq_2bp_change'],
             'usearch_ecoli_seq': ['usearch_ecoli_seq'],
             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change']}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_default_params_save_intermediate_files(self):
        """ Preserves files if save_intermediate_files/logs is True """

        intermediate_files_dir = self.output_dir + "/test_usearch61/"
        create_dir(intermediate_files_dir)
        self._dirs_to_remove.append(intermediate_files_dir)

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': True,
                    'output_dir':
                    intermediate_files_dir,
                    'remove_usearch_logs': False
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        expected_intermediate_fps =\
            [join(intermediate_files_dir, "abundance_sorted.fna"),
             join(intermediate_files_dir, "abundance_sorted.log"),
             join(intermediate_files_dir, "abundance_sorted.uc"),
             join(intermediate_files_dir, "ref_clustered.log"),
             join(intermediate_files_dir, "ref_clustered.uc")]

        for curr_file in expected_intermediate_fps:
            self.assertTrue(exists(curr_file))

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_default_params_save_intermediate_files_fast_cluster(self):
        """ Preserves files if save_intermediate_files/logs is True """

        intermediate_files_dir = self.output_dir + "/test_usearch61_fast_1160/"
        create_dir(intermediate_files_dir)
        self._dirs_to_remove.append(intermediate_files_dir)

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': True,
                    'output_dir':
                    intermediate_files_dir,
                    'remove_usearch_logs': False,
                    'usearch61_sort_method':
                    'length',
                    'usearch_fast_cluster': True
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        expected_intermediate_fps =\
            [join(intermediate_files_dir, "abundance_sorted.fna"),
             join(intermediate_files_dir, "abundance_sorted.log"),
             join(intermediate_files_dir, "abundance_sorted.uc"),
             join(intermediate_files_dir, "ref_clustered.log"),
             join(intermediate_files_dir, "ref_clustered.uc")]

        for curr_file in expected_intermediate_fps:
            self.assertTrue(exists(curr_file))

        expected_clusters = {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_2bp_change',
                                                              'usearch_ecoli_seq_1bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    # Don't have a good way to catch this error currently
    '''def test_call_default_params_minlen(self):
        """ Discards reads that fall below minlen setting """

        app = Usearch61ReferenceOtuPicker(params={'save_intermediate_files':False,
                                            'output_dir':self.output_dir,
                                            'remove_usearch_logs':True,
                                            'minlen':101
                                           })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
         refseqs_fp = self.tmp_seq_filepath_97perc_id_rc)

        # Should get no results
        expected_clusters = {}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = []
        self.assertEqual(failures, expected_failures)'''

    def test_usearch61_params(self):
        """ usearch61 handles changes to other parameters """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'wordlength': 25,
                    'usearch61_maxrejects': 200,
                    'usearch61_maxaccepts': 5
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     otu_prefix="test", refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        # won't get 2bp_change as reference, due to RC status
        expected_clusters = {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change',
                                                              'usearch_ecoli_seq_2bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_usearch61_length_sorting(self):
        """ Sorting according to length, clusters seqs """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'usearch61_sort_method':
                    'length'
                    })

        obs_clusters, failures = app(self.tmp_seqs_usearch97perc_id_len_diff,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        expected_clusters = {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change',
                                                              'usearch_ecoli_seq_2bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_usearch61_sizeorder(self):
        """ Handles sizeorder option, clusters seqs """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'sizeorder': True
                    })

        obs_clusters, failures = app(self.tmp_seqs_usearch_97perc_dups,
                                     refseqs_fp=self.tmp_seq_filepath_97perc_id_rc)

        # Should have ecoli match ecoli, and remaining seqs match 1bp change.
        expected_clusters = {'usearch_ecoli_seq': ['usearch_ecoli_seq'],
                             'usearch_ecoli_seq_1bp_change': ['usearch_ecoli_seq_1bp_change',
                                                              'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq_1bp_change_dup1',
                                                              'usearch_ecoli_seq_1bp_change_dup2']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_closed_reference_usearch61(self):
        """ usearch61 does closed reference OTU picking successfully """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'suppress_new_clusters': True
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seqs_rc_single_seq)

        # Randomly selected match is used for equivalent matches, so need to
        # test for results without order affecting output
        expected_clusters = {}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = ['usearch_ecoli_seq',
                             'usearch_ecoli_seq_2bp_change', 'usearch_ecoli_seq_1bp_change']
        assert_almost_equal(failures, expected_failures)

    def test_closed_reference_with_match_usearch61(self):
        """ usearch61 does closed reference OTU picking successfully """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'suppress_new_clusters': True
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id_rc,
                                     refseqs_fp=self.tmp_seqs_rc_single_seq)

        # Randomly selected match is used for equivalent matches, so need to
        # test for results without order affecting output
        expected_clusters = {'usearch_ecoli_seq_2bp_change_rc':
                             ['usearch_ecoli_seq_2bp_change_rc']}

        self.assertEqual(obs_clusters, expected_clusters)

        expected_failures = ['usearch_ecoli_seq',
                             'usearch_ecoli_seq_1bp_change']
        self.assertEqual(set(failures), set(expected_failures))

    def test_call_open_reference_usearch61(self):
        """ usearch61 does open reference OTU picking successfully """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'suppress_new_clusters':
                    False
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id,
                                     refseqs_fp=self.tmp_seqs_rc_single_seq)

        # Randomly selected match is used for equivalent matches, so need to
        # test for results without order affecting output
        expected_clusters = {'denovo0': ['usearch_ecoli_seq',
                                         'usearch_ecoli_seq_2bp_change',
                                         'usearch_ecoli_seq_1bp_change']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)

    def test_call_open_reference_with_match_usearch61(self):
        """ usearch61 does open reference OTU picking successfully """

        app = Usearch61ReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'output_dir':
                    self.output_dir,
                    'remove_usearch_logs': True,
                    'suppress_new_clusters':
                    False
                    })

        obs_clusters, failures = app(self.tmp_seq_filepath_97perc_id_rc,
                                     refseqs_fp=self.tmp_seqs_rc_single_seq)

        # Randomly selected match is used for equivalent matches, so need to
        # test for results without order affecting output
        expected_clusters = {'denovo0': ['usearch_ecoli_seq',
                                         'usearch_ecoli_seq_1bp_change'],
                             'usearch_ecoli_seq_2bp_change_rc':
                             ['usearch_ecoli_seq_2bp_change_rc']}

        for result in obs_clusters:
            for cluster in obs_clusters[result]:
                self.assertTrue(cluster in expected_clusters[result])

        expected_failures = []
        self.assertEqual(failures, expected_failures)


class UsearchOtuPickerTests(TestCase):

    """ Tests of the usearch-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.dna_seqs_3 = dna_seqs_3
        self.dna_seqs_3_derep = dna_seqs_3_derep
        self.dna_seqs_4 = dna_seqs_usearch
        self.ref_database = usearch_ref_seqs1

        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1, 'w')
        seq_file.write(self.dna_seqs_3)
        seq_file.close()

        self.tmp_seq_filepath1_derep = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1_derep, 'w')
        seq_file.write(self.dna_seqs_3_derep)
        seq_file.close()

        self.tmp_seq_filepath2 = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2, 'w')
        seq_file.write(self.dna_seqs_4)
        seq_file.close()

        self.tmp_ref_database = get_tmp_filename(
            prefix='UsearchRefDatabase_',
            suffix='.fasta')
        seq_file = open(self.tmp_ref_database, 'w')
        seq_file.write(self.ref_database)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath1, self.tmp_seq_filepath1_derep,
             self.tmp_seq_filepath2, self.tmp_ref_database]

        self._dirs_to_remove = []

    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)

    def seqs_to_temp_fasta(self, seqs):
        """ """
        fp = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(fp, 'w')
        self._files_to_remove.append(fp)
        for s in seqs:
            seq_file.write('>%s\n%s\n' % s)
        seq_file.close()
        return fp

    def test_call_default_params(self):
        """UsearchOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # All seqs should create own cluster

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': True,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1
                                       })

        obs = app(self.tmp_seq_filepath1)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_derep(self):
        """UsearchOtuPicker.__call__ returns expected clusters when using
        --derep_fullseq"""

        # adapted from test_call_default_params
        # Sequences 1 and 9 have exact replicates

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1',
                         'uclust_test_seqs_1rep',
                         'uclust_test_seqs_1rep2'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9',
                         'uclust_test_seqs_9rep']]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'derep_fullseq': True
                                       })

        obs = app(self.tmp_seq_filepath1_derep)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_default_no_reference(self):
        """UsearchOtuPicker.__call__ returns expected clusters no referencedb"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # All seqs should create own cluster

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1
                                       })

        obs = app(self.tmp_seq_filepath1)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_low_cluster_identity(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should only get 6 clusters
        exp_otu_ids = [str(x) for x in range(6)]

        exp_clusters =\
            [['uclust_test_seqs_0', 'uclust_test_seqs_6',
              'uclust_test_seqs_9'], ['uclust_test_seqs_1'], ['uclust_test_seqs_2'],
             ['uclust_test_seqs_3',
              'uclust_test_seqs_5',
              'uclust_test_seqs_8'],
                ['uclust_test_seqs_4'], ['uclust_test_seqs_7']]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': False,
                                       'de_novo_chimera_detection': False,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.80,
                                       'percent_id_err': 0.97
                                       })

        obs = app(self.tmp_seq_filepath1)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_detects_de_novo_chimeras(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should

        exp_otu_ids = ['2', '3', '4', '5']

        exp_clusters = [['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_7']
                        ]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': False,
                                       'de_novo_chimera_detection': True,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.97,
                                       'percent_id_err': 0.80,
                                       'abundance_skew': 2
                                       })

        obs = app(self.tmp_seq_filepath1)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_detects_reference_chimeras(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1']

        exp_clusters = [['Solemya', 'Solemya_seq2'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': True,
                                       'de_novo_chimera_detection': False,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.97,
                                       'percent_id_err': 0.97,
                                       'abundance_skew': 2
                                       })

        obs = app(self.tmp_seq_filepath2)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_usearch_handles_intersections(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1']

        exp_clusters = [['Solemya', 'Solemya_seq2'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': True,
                                       'de_novo_chimera_detection': True,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.97,
                                       'percent_id_err': 0.97,
                                       'abundance_skew': 2,
                                       'chimeras_retention': 'intersection'
                                       })

        obs = app(self.tmp_seq_filepath2)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_usearch_handles_unions(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1', '2']

        # will retain 'chimera' with union option.
        exp_clusters = [['Solemya', 'Solemya_seq2'], ['chimera'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': True,
                                       'de_novo_chimera_detection': True,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.97,
                                       'percent_id_err': 0.97,
                                       'abundance_skew': 2,
                                       'chimeras_retention': 'union'
                                       })

        obs = app(self.tmp_seq_filepath2)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_writes_output(self):
        """UsearchOtuPicker.__call__ writes expected output clusters file"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection, then write the OTU mapping file in
        # QIIME format.

        self.tmp_result_path = get_tmp_filename(
            prefix='UsearchOTUMapping_',
            suffix='.txt')
        f = open(self.tmp_result_path, "w")

        self.tmp_failures_path = get_tmp_filename(
            prefix='UsearchFailures_',
            suffix='.txt')
        f = open(self.tmp_failures_path, "w")

        self._files_to_remove.append(self.tmp_result_path)
        self._files_to_remove.append(self.tmp_failures_path)

        app = UsearchOtuPicker(params={'save_intermediate_files': False,
                                       'db_filepath': self.tmp_ref_database,
                                       'output_dir': self.temp_dir,
                                       'remove_usearch_logs': True,
                                       'reference_chimera_detection': True,
                                       'de_novo_chimera_detection': False,
                                       'cluster_size_filtering': False,
                                       'minlen': 12,
                                       'w': 12,
                                       'minsize': 1,
                                       'percent_id': 0.97,
                                       'percent_id_err': 0.97,
                                       'abundance_skew': 2
                                       })

        obs = app(self.tmp_seq_filepath2, result_path=self.tmp_result_path,
                  failure_path=self.tmp_failures_path)

        expected_otu_mapping =\
            ["1\tSolemya\tSolemya_seq2\n",
             "0\tusearch_ecoli_seq\tusearch_ecoli_seq2\n"""
             ]

        f = open(self.tmp_result_path, "U")

        actual_otu_mapping = f.readlines()

        self.assertEqual(actual_otu_mapping, expected_otu_mapping)

        expected_failures = ["chimera"]

        f = open(self.tmp_failures_path, "U")

        actual_failures = f.readlines()

        self.assertEqual(actual_failures, expected_failures)


class UsearchReferenceOtuPickerTests(TestCase):

    """ Tests of the usearch-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.dna_seqs_3 = dna_seqs_3
        self.dna_seqs_3_derep = dna_seqs_3_derep
        self.dna_seqs_4 = dna_seqs_usearch
        self.ref_database = usearch_ref_seqs1
        self.otu_ref_database = uclustref_query_seqs1

        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1, 'w')
        seq_file.write(self.dna_seqs_3)
        seq_file.close()

        self.tmp_seq_filepath1_derep = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1_derep, 'w')
        seq_file.write(self.dna_seqs_3_derep)
        seq_file.close()

        self.tmp_seq_filepath2 = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2, 'w')
        seq_file.write(self.dna_seqs_4)
        seq_file.close()

        self.tmp_ref_database = get_tmp_filename(
            prefix='UsearchRefDatabase_',
            suffix='.fasta')
        seq_file = open(self.tmp_ref_database, 'w')
        seq_file.write(self.ref_database)
        seq_file.close()

        self.tmp_otu_ref_database = get_tmp_filename(
            prefix='UsearchRefOtuDatabase_',
            suffix='.fasta')
        seq_file = open(self.tmp_otu_ref_database, 'w')
        seq_file.write(self.otu_ref_database)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath1, self.tmp_seq_filepath2,
             self.tmp_seq_filepath1_derep, self.tmp_ref_database,
             self.tmp_otu_ref_database]

        self._dirs_to_remove = []

    def tearDown(self):
        remove_files(self._files_to_remove)
        if self._dirs_to_remove:
            for curr_dir in self._dirs_to_remove:
                rmtree(curr_dir)

    def seqs_to_temp_fasta(self, seqs):
        """ """
        fp = get_tmp_filename(
            prefix='UsearchOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(fp, 'w')
        self._files_to_remove.append(fp)
        for s in seqs:
            seq_file.write('>%s\n%s\n' % s)
        seq_file.close()
        return fp

    def test_call_default_params(self):
        """UsearchOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # All seqs should create own cluster

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1
                    })

        obs = app(self.tmp_seq_filepath1, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_derep(self):
        """UsearchOtuPicker.__call__ returns expected clusters when using
        --derep_fullseq"""

        # adapted from test_call_default_params
        # Sequences 1 and 9 have exact replicates

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1',
                         'uclust_test_seqs_1rep',
                         'uclust_test_seqs_1rep2'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9',
                         'uclust_test_seqs_9rep']]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'derep_fullseq': True
                    })

        obs = app(self.tmp_seq_filepath1_derep, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_default_no_reference(self):
        """UsearchOtuPicker.__call__ returns expected clusters no referencedb"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # All seqs should create own cluster

        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1
                    })

        obs = app(self.tmp_seq_filepath1, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_low_cluster_identity(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should only get 6 clusters
        exp_otu_ids = [str(x) for x in range(10)]

        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    False,
                    'de_novo_chimera_detection':
                    False,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.80,
                    'percent_id_err': 0.97
                    })

        obs = app(self.tmp_seq_filepath1, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_detects_de_novo_chimeras(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should

        exp_otu_ids = ['0', '1', '2', '3']

        exp_clusters = [['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_7']
                        ]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    False,
                    'de_novo_chimera_detection':
                    True,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.97,
                    'percent_id_err': 0.80,
                    'abundance_skew': 2
                    })

        obs = app(self.tmp_seq_filepath1, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_detects_reference_chimeras(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1']

        exp_clusters = [['Solemya', 'Solemya_seq2'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'de_novo_chimera_detection':
                    False,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.97,
                    'percent_id_err': 0.97,
                    'abundance_skew': 2
                    })

        obs = app(self.tmp_seq_filepath2, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_usearch_handles_intersections(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1']

        exp_clusters = [['Solemya', 'Solemya_seq2'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'de_novo_chimera_detection':
                    True,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.97,
                    'percent_id_err': 0.97,
                    'abundance_skew': 2,
                    'chimeras_retention':
                    'intersection'
                    })

        obs = app(self.tmp_seq_filepath2, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_usearch_handles_unions(self):
        """UsearchOtuPicker.__call__ returns expected clusters"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection

        exp_otu_ids = ['0', '1', '2']

        # will retain 'chimera' with union option.
        exp_clusters = [['Solemya', 'Solemya_seq2'], ['chimera'],
                        ['usearch_ecoli_seq', 'usearch_ecoli_seq2']
                        ]

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'de_novo_chimera_detection':
                    True,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.97,
                    'percent_id_err': 0.97,
                    'abundance_skew': 2,
                    'chimeras_retention': 'union'
                    })

        obs = app(self.tmp_seq_filepath2, self.tmp_otu_ref_database)

        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_writes_output(self):
        """UsearchOtuPicker.__call__ writes expected output clusters file"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        # Should detect and remove chimeric sequence based
        # during ref based detection, then write the OTU mapping file in
        # QIIME format.

        self.tmp_result_path = get_tmp_filename(
            prefix='UsearchOTUMapping_',
            suffix='.txt')
        f = open(self.tmp_result_path, "w")

        self.tmp_failures_path = get_tmp_filename(
            prefix='UsearchFailures_',
            suffix='.txt')
        f = open(self.tmp_failures_path, "w")

        self._files_to_remove.append(self.tmp_result_path)
        self._files_to_remove.append(self.tmp_failures_path)

        app = UsearchReferenceOtuPicker(
            params={'save_intermediate_files': False,
                    'db_filepath':
                    self.tmp_ref_database,
                    'output_dir': self.temp_dir,
                    'remove_usearch_logs': True,
                    'reference_chimera_detection':
                    True,
                    'de_novo_chimera_detection':
                    False,
                    'cluster_size_filtering':
                    False,
                    'minlen': 12,
                    'w': 12,
                    'minsize': 1,
                    'percent_id': 0.97,
                    'percent_id_err': 0.97,
                    'abundance_skew': 2
                    })

        obs = app(self.tmp_seq_filepath2, self.tmp_otu_ref_database,
                  result_path=self.tmp_result_path,
                  failure_path=self.tmp_failures_path)

        expected_otu_mapping =\
            ["1\tSolemya\tSolemya_seq2\n",
             "0\tusearch_ecoli_seq\tusearch_ecoli_seq2\n"""
             ]

        f = open(self.tmp_result_path, "U")

        actual_otu_mapping = f.readlines()

        self.assertEqual(actual_otu_mapping, expected_otu_mapping)

        expected_failures = ["chimera"]

        f = open(self.tmp_failures_path, "U")

        actual_failures = f.readlines()

        self.assertEqual(actual_failures, expected_failures)


class UclustOtuPickerTests(TestCase):

    """ Tests of the uclust-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(
            prefix='UclustOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1, 'w')
        seq_file.write(dna_seqs_3)
        seq_file.close()

        self.tmp_seq_filepath2 = get_tmp_filename(
            prefix='UclustOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2, 'w')
        seq_file.write(dna_seqs_4)
        seq_file.close()

        self.tmp_seq_filepath3 = get_tmp_filename(
            prefix='UclustOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath3, 'w')
        seq_file.write(dna_seqs_5)
        seq_file.close()

        self.tmp_seq_filepath4 = get_tmp_filename(
            prefix='UclustOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath4, 'w')
        seq_file.write(dna_seqs_6)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath1, self.tmp_seq_filepath2,
             self.tmp_seq_filepath3, self.tmp_seq_filepath4]

    def tearDown(self):
        remove_files(self._files_to_remove)

    def seqs_to_temp_fasta(self, seqs):
        """ """
        fp = get_tmp_filename(
            prefix='UclustReferenceOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(fp, 'w')
        self._files_to_remove.append(fp)
        for s in seqs:
            seq_file.write('>%s\n%s\n' % s)
        seq_file.close()
        return fp

    def test_toggle_collapse_identical_sequences(self):
        """UclustOtuPicker: toggle prefilter identical seqs doesn't affect clusters
        """

        # generate result including prefilter
        app_w_collapse_identical =\
            UclustOtuPicker(params={'Similarity': 0.90,
                                    'save_uc_files': False,
                                    'prefilter_identical_sequences': True,
                                    'output_dir': self.temp_dir})
        result_w_collapse_identical = \
            sorted(app_w_collapse_identical(self.tmp_seq_filepath4).values())

        # generate result excluding prefilter
        app_wo_collapse_identical =\
            UclustOtuPicker(params={'Similarity': 0.90,
                                    'save_uc_files': False,
                                    'prefilter_identical_sequences': False,
                                    'output_dir': self.temp_dir})
        result_wo_collapse_identical = \
            sorted(app_wo_collapse_identical(self.tmp_seq_filepath4).values())

        self.assertEqual(result_w_collapse_identical,
                         result_wo_collapse_identical)

    def test_toggle_suppress_sort(self):
        """UclustOtuPicker: togging suppress sort functions as expected
        """
        seqs = [('s1', 'ACCTTGTTACTTT'),  # three copies
                ('s2', 'ACCTTGTTACTTTC'),  # one copy
                ('s3', 'ACCTTGTTACTTTCC'),  # two copies
                ('s4', 'ACCTTGTTACTTT'),
                ('s5', 'ACCTTGTTACTTTCC'),
                ('s6', 'ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)

        # no abundance sorting and uclust's sorting enabled
        # so length-based sorting
        app = UclustOtuPicker(params={'Similarity': 0.80,
                                      'enable_rev_strand_matching': False,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(seqs_fp)
        exp = {0: ['s3', 's5', 's2', 's1', 's4', 's6']}
        self.assertEqual(obs, exp)

        # no abundance sorting and uclust's sorting enabled
        # so no sorting at all
        app = UclustOtuPicker(params={'Similarity': 0.80,
                                      'enable_rev_strand_matching': False,
                                      'suppress_sort': True,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(seqs_fp)
        exp = {0: ['s1', 's4', 's6', 's2', 's3', 's5']}
        self.assertEqual(obs, exp)

    def test_abundance_sort(self):
        """UclustOtuPicker: abundance sort functions as expected
        """
        # enable abundance sorting with suppress sort = False (it gets
        # set to True internally, otherwise uclust's length sort would
        # override the abundance sorting)
        seqs = [('s1 comment1', 'ACCTTGTTACTTT'),  # three copies
                ('s2 comment2', 'ACCTTGTTACTTTC'),  # one copy
                ('s3 comment3', 'ACCTTGTTACTTTCC'),  # two copies
                ('s4 comment4', 'ACCTTGTTACTTT'),
                ('s5 comment5', 'ACCTTGTTACTTTCC'),
                ('s6 comment6', 'ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)

        # abundance sorting changes order
        app = UclustOtuPicker(params={'Similarity': 0.80,
                                      'enable_rev_strand_matching': False,
                                      'suppress_sort': False,
                                      'presort_by_abundance': True,
                                      'save_uc_files': False})
        obs = app(seqs_fp)
        exp = {0: ['s1', 's4', 's6', 's3', 's5', 's2']}
        self.assertEqual(obs, exp)

        # abundance sorting changes order -- same results with suppress_sort =
        # True b/c (it gets set to True to when presorting by abundance)
        app = UclustOtuPicker(params={'Similarity': 0.80,
                                      'enable_rev_strand_matching': False,
                                      'suppress_sort': True,
                                      'presort_by_abundance': True,
                                      'save_uc_files': False})
        obs = app(seqs_fp)
        exp = {0: ['s1', 's4', 's6', 's3', 's5', 's2']}
        self.assertEqual(obs, exp)

    def test_call_default_params(self):
        """UclustOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp_otu_ids = range(10)
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'save_uc_files': False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_default_params_suppress_sort(self):
        """UclustOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp_otu_ids = range(10)
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'save_uc_files': False,
                                      'suppress_sort': True})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_default_params_save_uc_file(self):
        """ returns expected clusters default params, writes correct .uc file"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        exp_otu_ids = range(10)
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_8'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'save_uc_files': True,
                                      'output_dir': self.temp_dir})
        obs = app(self.tmp_seq_filepath1)

        uc_fasta_fp = "_".join(self.tmp_seq_filepath1.split('_')[0:2])
        uc_output_fp = uc_fasta_fp.replace('.fasta', '_clusters.uc')

        uc_output_f = open(uc_output_fp, "U")
        self._files_to_remove.append(uc_output_fp)

        # Testing content of file minus header (tmp filename of sorted fasta
        # file difficult to access here).  Also not testing the version number
        # of uclust that could vary between systems but still function for the
        # purpose of generating appropriate clusters.

        uc_result = [line.strip() for line in uc_output_f][2:]

        self.assertEqual(uc_result, expected_uc_output)

        # Make sure other results are correct with uc file being saved.
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_alt_threshold(self):
        """UclustOtuPicker.__call__ returns expected clusters with alt threshold
        """
        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp_otu_ids = range(9)
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6', 'uclust_test_seqs_8'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_otu_id_prefix(self):
        """UclustOtuPicker.__call__ returns expected clusters with alt threshold
        """
        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp_otu_ids = ['my_otu_%d' % i for i in range(9)]
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6', 'uclust_test_seqs_8'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'new_cluster_identifier': 'my_otu_',
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_suppress_sort(self):
        """UclustOtuPicker.__call__ handles suppress sort
        """

        exp_otu_ids = range(3)
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2']]

        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'suppress_sort': True,
                                      'optimal': True,
                                      'enable_rev_strand_matching': True,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath2)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_call_rev_matching(self):
        """UclustOtuPicker.__call__ handles reverse strand matching
        """
        exp_otu_ids = range(2)
        exp_clusters = [['uclust_test_seqs_0'], ['uclust_test_seqs_0_rc']]
        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'enable_rev_strand_matching': False,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath3)
        obs_otu_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)

        exp = {0: ['uclust_test_seqs_0', 'uclust_test_seqs_0_rc']}
        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'enable_rev_strand_matching': True,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath3)
        self.assertEqual(obs, exp)

    def test_call_output_to_file(self):
        """UclustHitOtuPicker.__call__ output to file functions as expected
        """

        tmp_result_filepath = get_tmp_filename(
            prefix='UclustOtuPickerTest.test_call_output_to_file_',
            suffix='.txt')

        app = UclustOtuPicker(params={'Similarity': 0.90,
                                      'suppress_sort': False,
                                      'presort_by_abundance': False,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath1, result_path=tmp_result_filepath)

        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)

        exp_otu_ids = map(str, range(9))
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6', 'uclust_test_seqs_8'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_9']]
        obs_otu_ids = []
        obs_clusters = []
        for line in result_file_str.split('\n'):
            if line:
                fields = line.split('\t')
                obs_otu_ids.append(fields[0])
                obs_clusters.append(fields[1:])
        obs_otu_ids.sort()
        obs_clusters.sort()
        # The relation between otu ids and clusters is abitrary, and
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)
        # confirm that nothing is returned when result_path is specified
        self.assertEqual(obs, None)

    def test_call_log_file(self):
        """UclustOtuPicker.__call__ writes log when expected
        """

        tmp_log_filepath = get_tmp_filename(
            prefix='UclustOtuPickerTest.test_call_output_to_file_l_',
            suffix='.txt')
        tmp_result_filepath = get_tmp_filename(
            prefix='UclustOtuPickerTest.test_call_output_to_file_r_',
            suffix='.txt')

        app = UclustOtuPicker(params={'Similarity': 0.99,
                                      'save_uc_files': False})
        obs = app(self.tmp_seq_filepath1,
                  result_path=tmp_result_filepath, log_path=tmp_log_filepath)

        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        # remove the temp files before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)

        log_file_99_exp = ["UclustOtuPicker parameters:",
                           "Similarity:0.99", "Application:uclust",
                           "enable_rev_strand_matching:False",
                           "suppress_sort:True",
                           "optimal:False",
                           'max_accepts:20',
                           'max_rejects:500',
                           'stepwords:20',
                           'word_length:12',
                           "exact:False",
                           "Num OTUs:10",
                           "new_cluster_identifier:None",
                           "presort_by_abundance:True",
                           "stable_sort:True",
                           "output_dir:.",
                           "save_uc_files:False",
                           "prefilter_identical_sequences:True",
                           "Result path: %s" % tmp_result_filepath]
        # compare data in log file to fake expected log file
        # NOTE: Since app.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that
        # the equal unordered lists of lines is present in actual and expected
        assert_almost_equal(log_file_str.split('\n'), log_file_99_exp)

    def test_map_filtered_clusters_to_full_clusters(self):
        """UclustOtuPicker._map_filtered_clusters_to_full_clusters functions as expected
        """
        # original and mapped full clusters are the same
        app = UclustOtuPicker(params={})
        filter_map = {'s1': ['s1'], 's2': ['s2'],
                      's3': ['s3'], 's4': ['s4'],
                      's5': ['s5'], 's6': ['s6']}
        clusters = [['s1'], ['s2'], ['s3'], ['s4'], ['s5'], ['s6']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        expected = clusters
        self.assertEqual(actual, expected)

        # original and mapped full clusters are not the same
        filter_map = {'s1': ['s1', 's2', 's3', 's4'], 's5': ['s5', 's6']}
        clusters = [['s1', 's5']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        for e in actual:
            e.sort()
        expected = [['s1', 's2', 's3', 's4', 's5', 's6']]
        self.assertEqual(actual, expected)

        filter_map = {'s1': ['s1', 's2', 's6'],
                      's3': ['s3'], 's5': ['s4', 's5']}
        clusters = [['s1', 's3'], ['s5']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        for e in actual:
            e.sort()
        expected = [['s1', 's2', 's3', 's6'], ['s4', 's5']]
        self.assertEqual(actual, expected)


class UclustReferenceOtuPickerTests(TestCase):

    """ Tests of the uclust reference-based OTU picker """

    def setUp(self):
        """ """
        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(
            prefix='UclustReferenceOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1, 'w')
        seq_file.write(uclustref_query_seqs1)
        seq_file.close()

        self.temp_ref_filepath1 = get_tmp_filename(
            prefix='UclustReferenceOtuPickerTest_',
            suffix='.fasta')
        ref_file = open(self.temp_ref_filepath1, 'w')
        ref_file.write(uclustref_ref_seqs1)
        ref_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath1,
             self.temp_ref_filepath1]

    def tearDown(self):
        remove_files(self._files_to_remove)

    def seqs_to_temp_fasta(self, seqs):
        """ """
        fp = get_tmp_filename(
            prefix='UclustReferenceOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(fp, 'w')
        self._files_to_remove.append(fp)
        for s in seqs:
            seq_file.write('>%s\n%s\n' % s)
        seq_file.close()
        return fp

    def test_toggle_suppress_sort(self):
        """UclustReferenceOtuPicker: togging suppress sort functions as expected
        """
        seqs = [('s1 comment1', 'ACCTTGTTACTTT'),  # three copies
                ('s2 comment2', 'ACCTTGTTACTTTC'),  # one copy
                ('s3 comment3', 'ACCTTGTTACTTTCC'),  # two copies
                ('s4 comment4', 'ACCTTGTTACTTT'),
                ('s5 comment5', 'ACCTTGTTACTTTCC'),
                ('s6 comment6', 'ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        ref_seqs = [('r1 blah', 'ACCTTGTTACTTT')]
        ref_seqs_fp = self.seqs_to_temp_fasta(ref_seqs)

        # no abundance sorting and uclust's sorting enabled
        # so length-based sorting
        app = UclustReferenceOtuPicker(params={'Similarity': 0.80,
                                               'enable_rev_strand_matching':
                                               False,
                                               'suppress_sort': False,
                                               'presort_by_abundance': False,
                                               'save_uc_files': False})
        obs = app(seqs_fp, ref_seqs_fp)
        exp = {'r1': ['s3', 's5', 's2', 's1', 's4', 's6']}
        self.assertEqual(obs, exp)

        # no abundance sorting and uclust's sorting enabled
        # so no sorting at all
        app = UclustReferenceOtuPicker(params={'Similarity': 0.80,
                                               'enable_rev_strand_matching':
                                               False,
                                               'suppress_sort': True,
                                               'presort_by_abundance': False,
                                               'save_uc_files': False})
        obs = app(seqs_fp, ref_seqs_fp)
        exp = {'r1': ['s1', 's4', 's6', 's2', 's3', 's5']}
        self.assertEqual(obs, exp)

    def test_abundance_sort(self):
        """UclustReferenceOtuPicker: abundance sort functions as expected
        """
        # enable abundance sorting with suppress sort = False (it gets
        # set to True internally, otherwise uclust's length sort would
        # override the abundance sorting)
        seqs = [('s1 comment1', 'ACCTTGTTACTTT'),  # three copies
                ('s2 comment2', 'ACCTTGTTACTTTC'),  # one copy
                ('s3 comment3', 'ACCTTGTTACTTTCC'),  # two copies
                ('s4 comment4', 'ACCTTGTTACTTT'),
                ('s5 comment5', 'ACCTTGTTACTTTCC'),
                ('s6 comment6', 'ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        ref_seqs = [('r1 blah', 'ACCTTGTTACTTT')]
        ref_seqs_fp = self.seqs_to_temp_fasta(ref_seqs)

        # abundance sorting changes order
        app = UclustReferenceOtuPicker(params={'Similarity': 0.80,
                                               'enable_rev_strand_matching':
                                               False,
                                               'suppress_sort': True,
                                               'presort_by_abundance': True,
                                               'save_uc_files': False})
        obs = app(seqs_fp, ref_seqs_fp)
        exp = {'r1': ['s1', 's4', 's6', 's3', 's5', 's2']}
        self.assertEqual(obs, exp)

        # abundance sorting changes order -- same results with suppress_sort =
        # True b/c (it gets set to True to when presorting by abundance)
        app = UclustReferenceOtuPicker(params={'Similarity': 0.80,
                                               'enable_rev_strand_matching':
                                               False,
                                               'suppress_sort': True,
                                               'presort_by_abundance': True,
                                               'save_uc_files': False})
        obs = app(seqs_fp, ref_seqs_fp)
        exp = {'r1': ['s1', 's4', 's6', 's3', 's5', 's2']}
        self.assertEqual(obs, exp)

    def test_toggle_suppress_new_clusters(self):
        """UclustReferenceOtuPicker: toggle suppress new clusters
        """
        seqs = [('s1 a', 'ACCTTGTTACTTT'),
                ('s2 bb', 'ACCTAGTTACTTT'),
                ('s3 c  c', 'TTGCGTAACGTTTGAC')]
        ref_seqs = [
            ('r1 d', 'ACCTCGTTACTTT')]
        # these seqs should match at 0.90, but don't -- I can confirm this
        # running uclust directly, and have contacted Robert Edgar for
        # clarification
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': True,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        exp = {'r1': ['s1', 's2']}
        self.assertEqual(obs, exp)

        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': False,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        exp = {'r1': ['s1', 's2'], 'new_42': ['s3']}
        self.assertEqual(obs, exp)

    def test_toggle_collapse_identical_sequences_prefilter_w_new_clusters(
            self):
        """UclustReferenceOtuPicker: ident. seqs prefilter fns w new clusters
        """
        # s4 == s2 and s3 == s5
        seqs = [('s1 a', 'ACCTTGTTACTTT'),
                ('s2 bb', 'ACCTAGTTACTTT'),
                ('s4 bb', 'ACCTAGTTACTTT'),
                ('s3 c  c', 'TTGCGTAACGTTTGAC'),
                ('s5 c  c', 'TTGCGTAACGTTTGAC')]
        ref_seqs = [
            ('r1 d', 'ACCTCGTTACTTT')]
        exp = {'r1': ['s2', 's4', 's1'],
               'new_42': ['s3', 's5']}

        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': False,
                                       'save_uc_files': False,
                                       'prefilter_identical_sequences': False})
        obs_no_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                              self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        self.assertEqual(obs_no_prefilter, exp)

        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': False,
                                       'save_uc_files': False,
                                       'prefilter_identical_sequences': True})
        obs_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                           self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        self.assertEqual(obs_prefilter, exp)

        # a little paranoia never hurt anyone
        self.assertEqual(obs_prefilter, obs_no_prefilter)

    def test_toggle_collapse_identical_sequences_prefilter_wo_new_clusters(
            self):
        """UclustReferenceOtuPicker: ident. seqs prefilter fns wo new clusters
        """
        # s4 == s2 and s3 == s5
        seqs = [('s1 a', 'ACCTTGTTACTTT'),
                ('s2 bb', 'ACCTAGTTACTTT'),
                ('s4 bb', 'ACCTAGTTACTTT'),
                ('s3 c  c', 'TTGCGTAACGTTTGAC'),
                ('s5 c  c', 'TTGCGTAACGTTTGAC')]
        ref_seqs = [
            ('r1 d', 'ACCTCGTTACTTT')]

        exp = {'r1': ['s2', 's4', 's1']}

        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': True,
                                       'save_uc_files': False,
                                       'prefilter_identical_sequences': False})
        fail_path_no_prefilter = get_tmp_filename(
            prefix='UclustRefOtuPickerFailures', suffix='.txt')
        self._files_to_remove.append(fail_path_no_prefilter)
        obs_no_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                              self.seqs_to_temp_fasta(ref_seqs),
                              failure_path=fail_path_no_prefilter,
                              HALT_EXEC=False)
        self.assertEqual(obs_no_prefilter, exp)
        self.assertEqual(open(fail_path_no_prefilter).read(),
                         "s3\ns5")

        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': True,
                                       'save_uc_files': False,
                                       'prefilter_identical_sequences': True})
        fail_path_prefilter = get_tmp_filename(
            prefix='UclustRefOtuPickerFailures', suffix='.txt')
        self._files_to_remove.append(fail_path_prefilter)
        obs_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                           self.seqs_to_temp_fasta(ref_seqs),
                           failure_path=fail_path_prefilter,
                           HALT_EXEC=False)
        self.assertEqual(obs_prefilter, exp)
        self.assertEqual(open(fail_path_prefilter).read(),
                         "s3\ns5")

        # a little paranoia never hurt anyone
        self.assertEqual(obs_prefilter, obs_no_prefilter)
        self.assertEqual(open(fail_path_prefilter).read(),
                         open(fail_path_no_prefilter).read())

    def test_varied_similarity(self):
        """UclustReferenceOtuPicker: varying similarity affects clustering
        """
        seqs = [('s1', 'ACCTTGTTACTTT'),
                ('s2', 'ACCTAGTTACTTT')]
        ref_seqs = [
            ('r1', 'ACCTCGTTACTTT')]
        # these seqs should match at 0.90, but don't -- I can confirm this
        # running uclust directly, and have contacted Robert Edgar for
        # clarification
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'suppress_new_clusters': False,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        exp = {'r1': ['s1', 's2']}
        self.assertEqual(obs, exp)

        # set similarity to 100%
        uc = UclustReferenceOtuPicker({'Similarity': 1.0,
                                       'suppress_new_clusters': False,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of
        # iteration over a dict
        exp1 = {'QiimeOTU1': ['s1'], 'QiimeOTU2': ['s2']}
        exp2 = {'QiimeOTU2': ['s1'], 'QiimeOTU1': ['s2']}
        self.assertTrue(obs == exp1 or obs == exp2)

    def test_toggle_rev_strand_matching(self):
        """UclustReferenceOtuPicker: toggle rev strand matching
        """
        # s3 and s4 are rc of one another
        seqs = [('s1', 'ACCTTGTTACTTT'),
                ('s2', 'ACCTAGTTACTTT'),
                ('s3', 'TTGCGTAACGTTTGAC'),
                ('s4', 'GTCAAACGTTACGCAA')]
        ref_seqs = [
            ('r1', 'ACCTCGTTACTTT')]

        # rev strand matching disabled
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'enable_rev_strand_matching': False,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        exp = {'r1': ['s1', 's2'], 'new_42': ['s3'], 'new_43': ['s4']}
        self.assertEqual(obs, exp)

        # enable rev strand matching
        uc = UclustReferenceOtuPicker({'Similarity': 0.80,
                                       'new_cluster_identifier': 'new_',
                                       'next_new_cluster_number': 42,
                                       'enable_rev_strand_matching': True,
                                       'save_uc_files': False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs), HALT_EXEC=False)
        exp = {'r1': ['s1', 's2'], 'new_42': ['s3', 's4']}
        self.assertEqual(obs, exp)

    def test_call_log_file(self):
        """UclustReferenceOtuPicker.__call__ writes log when expected
        """
        tmp_log_filepath = get_tmp_filename(prefix='UclustReferenceOtuPicker',
                                            suffix='log')
        tmp_result_filepath = get_tmp_filename(
            prefix='UclustReferenceOtuPicker',
            suffix='txt')
        tmp_failure_filepath = get_tmp_filename(
            prefix='UclustReferenceOtuPicker',
            suffix='txt')
        seqs = [('s1', 'ACCTTGTTACTTT'),
                ('s2', 'ACCTAGTTACTTT'),
                ('s3', 'TTGCGTAACGTTTGAC'),
                ('s4', 'GTCAAACGTTACGCAA')]
        ref_seqs = [
            ('r1', 'ACCTCGTTACTTT')]

        # rev strand matching disabled
        uc = UclustReferenceOtuPicker({'Similarity': 0.8,
                                       'suppress_new_clusters': True,
                                       'save_uc_files': False,
                                       'suppress_sort': True})
        ref_seqs_fp = self.seqs_to_temp_fasta(ref_seqs)
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 ref_seqs_fp,
                 result_path=tmp_result_filepath,
                 log_path=tmp_log_filepath,
                 failure_path=tmp_failure_filepath)

        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        fail_file = open(tmp_failure_filepath)
        fail_file_str = fail_file.read()
        fail_file.close()
        # remove the temp files before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)
        remove(tmp_failure_filepath)

        log_file_99_exp = ["OtuPicker parameters:",
                           "Reference seqs:%s" % abspath(ref_seqs_fp),
                           "Similarity:0.8",
                           "Application:uclust",
                           "enable_rev_strand_matching:False",
                           "suppress_sort:True",
                           "suppress_new_clusters:True",
                           "optimal:False",
                           "exact:False",
                           "Num OTUs:1",
                           "Num new OTUs:0",
                           "Num failures:2",
                           'max_accepts:20',
                           'max_rejects:500',
                           'stepwords:20',
                           'word_length:12',
                           "stable_sort:True",
                           "new_cluster_identifier:QiimeOTU",
                           "next_new_cluster_number:1",
                           "presort_by_abundance:True",
                           'save_uc_files:False',
                           'output_dir:.',
                           'prefilter_identical_sequences:True',
                           "Result path: %s" % tmp_result_filepath]
        # compare data in log file to fake expected log file
        # NOTE: Since app.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that
        # the equal unordered lists of lines is present in actual and expected

        assert_almost_equal(log_file_str.split('\n'), log_file_99_exp)

        failures_file_99_exp = ["s3", "s4"]
        assert_almost_equal(fail_file_str.split('\n'), failures_file_99_exp)

    def test_default_parameters_new_clusters_allowed(self):
        """UclustReferenceOtuPicker: default parameters, new clusters allowed
        """
        uc = UclustReferenceOtuPicker({'save_uc_files': False})
        obs = uc(self.tmp_seq_filepath1, self.temp_ref_filepath1)
        exp = {'ref1': ['uclust_test_seqs_0'],
               'ref2': ['uclust_test_seqs_1'],
               'ref3': ['uclust_test_seqs_2'],
               'ref4': ['uclust_test_seqs_3'],
               'QiimeOTU1': ['uclust_test_seqs_4'],
               'QiimeOTU2': ['uclust_test_seqs_5'],
               'QiimeOTU3': ['uclust_test_seqs_6'],
               'QiimeOTU4': ['uclust_test_seqs_7'],
               'QiimeOTU5': ['uclust_test_seqs_8'],
               'QiimeOTU6': ['uclust_test_seqs_9']}

        # expected number of clusters observed
        self.assertEqual(len(obs), len(exp))

        expected_ref_hits = ['ref1', 'ref2', 'ref3', 'ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k], exp[k])

        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of
        # iteration over a dict
        exp_cluster_ids = sorted(exp.keys())
        exp_clusters = sorted(exp.values())

        obs_cluster_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())

        self.assertEqual(obs_cluster_ids, exp_cluster_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_default_parameters_new_clusters_allowed_save_uc_files(self):
        """UclustReferenceOtuPicker: default parameters, saves uc file
        """

        uc = UclustReferenceOtuPicker({'save_uc_files': True,
                                       'output_dir': self.temp_dir,
                                       'suppress_sort': True})

        obs = uc(self.tmp_seq_filepath1, self.temp_ref_filepath1)
        exp = {'ref1': ['uclust_test_seqs_0'],
               'ref2': ['uclust_test_seqs_1'],
               'ref3': ['uclust_test_seqs_2'],
               'ref4': ['uclust_test_seqs_3'],
               'QiimeOTU1': ['uclust_test_seqs_4'],
               'QiimeOTU2': ['uclust_test_seqs_5'],
               'QiimeOTU3': ['uclust_test_seqs_6'],
               'QiimeOTU4': ['uclust_test_seqs_7'],
               'QiimeOTU5': ['uclust_test_seqs_8'],
               'QiimeOTU6': ['uclust_test_seqs_9']}

        # expected number of clusters observed
        self.assertEqual(len(obs), len(exp))

        expected_ref_hits = ['ref1', 'ref2', 'ref3', 'ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k], exp[k])

        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of
        # iteration over a dict
        exp_cluster_ids = sorted(exp.keys())
        exp_clusters = sorted(exp.values())

        obs_cluster_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())

        self.assertEqual(obs_cluster_ids, exp_cluster_ids)
        self.assertEqual(obs_clusters, exp_clusters)

        uc_fasta_fp = "_".join(self.tmp_seq_filepath1.split('_')[0:2])
        uc_output_fp = uc_fasta_fp.replace('.fasta', '_clusters.uc')

        uc_output_f = open(uc_output_fp, "U")
        self._files_to_remove.append(uc_output_fp)

        # Testing content of file minus header (tmp filename of sorted fasta
        # file difficult to access here), and second line which could contain
        # slight variations in uclust versions but still function for the
        # purpose of generating correct clusters

        uc_result = [line.strip() for line in uc_output_f][2:]

        self.assertEqual(uc_result, expected_ref_uc_file)

    def test_alt_similarity_new_clusters_allowed(self):
        """UclustReferenceOtuPicker: alt parameters, new clusters allowed
        """
        uc = UclustReferenceOtuPicker({'Similarity': 0.90,
                                      'suppress_sort': False,
                                       'presort_by_abundance': False,
                                       'save_uc_files': False,
                                       'output_dir': self.temp_dir})
        obs = uc(self.tmp_seq_filepath1, self.temp_ref_filepath1)
        exp = {'ref1': ['uclust_test_seqs_0'],
               'ref2': ['uclust_test_seqs_1'],
               'ref3': ['uclust_test_seqs_2'],
               'ref4': ['uclust_test_seqs_3'],
               'QiimeOTU1': ['uclust_test_seqs_4'],
               'QiimeOTU2': ['uclust_test_seqs_5'],
               'QiimeOTU3': ['uclust_test_seqs_6', 'uclust_test_seqs_8'],
               'QiimeOTU4': ['uclust_test_seqs_7'],
               'QiimeOTU5': ['uclust_test_seqs_9']}

        # expected number of clusters observed
        self.assertEqual(len(obs), len(exp))

        expected_ref_hits = ['ref1', 'ref2', 'ref3', 'ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k], exp[k])

        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of
        # iteration over a dict
        exp_cluster_ids = sorted(exp.keys())
        exp_clusters = sorted(exp.values())

        obs_cluster_ids = sorted(obs.keys())
        obs_clusters = sorted(obs.values())

        self.assertEqual(obs_cluster_ids, exp_cluster_ids)
        self.assertEqual(obs_clusters, exp_clusters)

    def test_default_parameters_new_clusters_disallowed(self):
        """UclustReferenceOtuPicker: default params, new clusters not allowed
        """
        uc = UclustReferenceOtuPicker({'suppress_new_clusters': True,
                                       'save_uc_files': False})
        obs = uc(self.tmp_seq_filepath1, self.temp_ref_filepath1)
        exp = {'ref1': ['uclust_test_seqs_0'],
               'ref2': ['uclust_test_seqs_1'],
               'ref3': ['uclust_test_seqs_2'],
               'ref4': ['uclust_test_seqs_3']}

        # expected number of clusters observed
        self.assertEqual(obs, exp)

    def test_alt_parameters_new_clusters_disallowed(self):
        """UclustReferenceOtuPicker: alt params, new clusters not allowed
        """
        uc = UclustReferenceOtuPicker({'suppress_new_clusters': True,
                                       'Similarity': 1.0,
                                       'save_uc_files': False})
        obs = uc(self.tmp_seq_filepath1, self.temp_ref_filepath1)
        exp = {'ref3': ['uclust_test_seqs_2'], 'ref4': ['uclust_test_seqs_3']}

        # expected number of clusters observed
        self.assertEqual(obs, exp)


class CdHitOtuPickerTests(TestCase):

    """ Tests of the cd-hit-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.tmp_seq_filepath1 = get_tmp_filename(
            prefix='CdHitOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1, 'w')
        seq_file.write(dna_seqs_1)
        seq_file.close()

        self.tmp_seq_filepath2 = get_tmp_filename(
            prefix='CdHitOtuPickerTest_',
            suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2, 'w')
        seq_file.write(dna_seqs_2)
        seq_file.close()

        self._files_to_remove =\
            [self.tmp_seq_filepath1, self.tmp_seq_filepath2]

    def tearDown(self):
        remove_files(self._files_to_remove)

    def test_call_default_params(self):
        """CdHitOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp = {0: ['cdhit_test_seqs_0'],
               1: ['cdhit_test_seqs_1'],
               2: ['cdhit_test_seqs_2'],
               3: ['cdhit_test_seqs_3'],
               4: ['cdhit_test_seqs_4'],
               5: ['cdhit_test_seqs_5'],
               6: ['cdhit_test_seqs_6'],
               7: ['cdhit_test_seqs_7'],
               8: ['cdhit_test_seqs_8'],
               9: ['cdhit_test_seqs_9']}

        app = CdHitOtuPicker(params={})
        obs = app(self.tmp_seq_filepath1)
        self.assertEqual(obs, exp)

    def test_call_alt_threshold(self):
        """CdHitOtuPicker.__call__ returns expected clusters with alt threshold
        """
        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs

        exp = {0: ['cdhit_test_seqs_0'],
               1: ['cdhit_test_seqs_1'],
               2: ['cdhit_test_seqs_2'],
               3: ['cdhit_test_seqs_3'],
               4: ['cdhit_test_seqs_4'],
               5: ['cdhit_test_seqs_5'],
               6: ['cdhit_test_seqs_6', 'cdhit_test_seqs_8'],
               7: ['cdhit_test_seqs_7'],
               8: ['cdhit_test_seqs_9']}

        app = CdHitOtuPicker(params={'Similarity': 0.90})
        obs = app(self.tmp_seq_filepath1)
        self.assertEqual(obs, exp)

    def test_call_output_to_file(self):
        """CdHitOtuPicker.__call__ output to file functions as expected
        """

        tmp_result_filepath = get_tmp_filename(
            prefix='CdHitOtuPickerTest.test_call_output_to_file_',
            suffix='.txt')

        app = CdHitOtuPicker(params={'Similarity': 0.90})
        obs = app(self.tmp_seq_filepath1, result_path=tmp_result_filepath)

        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)

        # compare data in result file to fake expected file
        self.assertEqual(result_file_str, dna_seqs_result_file_90_exp)
        # confirm that nothing is returned when result_path is specified
        self.assertEqual(obs, None)

    def test_call_log_file(self):
        """CdHitOtuPicker.__call__ writes log when expected
        """

        tmp_log_filepath = get_tmp_filename(
            prefix='CdHitOtuPickerTest.test_call_output_to_file_l_',
            suffix='.txt')
        tmp_result_filepath = get_tmp_filename(
            prefix='CdHitOtuPickerTest.test_call_output_to_file_r_',
            suffix='.txt')

        app = CdHitOtuPicker(params={'Similarity': 0.99})
        obs = app(self.tmp_seq_filepath1,
                  result_path=tmp_result_filepath, log_path=tmp_log_filepath)

        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        # remove the temp files before running the test, so in
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)

        log_file_99_exp = ["CdHitOtuPicker parameters:",
                           "Similarity:0.99", "Application:cdhit",
                           'Algorithm:cdhit: "longest-sequence-first list removal algorithm"',
                           'No prefix-based prefiltering.',
                           "Result path: %s" % tmp_result_filepath]
        # compare data in log file to fake expected log file
        # NOTE: Since app.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that
        # the equal unordered lists of lines is present in actual and expected
        assert_almost_equal(log_file_str.split('\n'), log_file_99_exp)

    def test_prefilter_exact_prefixes_no_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when no seqs get filtered
        """
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTACAA'),
                ('s3', 'ACGTAG'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCAA'),
                ('s6', 'ACGTCCAAAAAAAAAAAA')]

        prefix_length = 6
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = seqs, {'s1': ['s1'], 's2': ['s2'],
                          's3': ['s3'], 's4': ['s4'],
                          's5': ['s5'], 's6': ['s6']}
        self.assertEqual(actual, expected)

        # same result if prefix_length is too long
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTACAA'),
                ('s3', 'ACGTAG'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCAA'),
                ('s6', 'ACGTCCAAAAAAAAAAAA')]
        prefix_length = 42
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = seqs, {'s1': ['s1'], 's2': ['s2'],
                          's3': ['s3'], 's4': ['s4'],
                          's5': ['s5'], 's6': ['s6']}
        self.assertEqual(actual, expected)

    def test_prefilter_exact_prefixes_all_to_one_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when all seqs map to one
        """
        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1 comment', 'ACGTAA'),
                ('s2', 'ACGTAC'),
                ('s3', 'ACGTAG'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCA'),
                ('s6', 'ACGTCC')]

        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = [('s1', 'ACGTAA')], {'s1':
                                        ['s1', 's2', 's3', 's4', 's5', 's6']}
        self.assertEqual(actual, expected)

        # maps to longest seq
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTACA'),
                ('s3', 'ACGTAG'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCA'),
                ('s6', 'ACGTCC')]

        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = [('s2', 'ACGTACA')], {'s2':
                                         ['s1', 's2', 's3', 's4', 's5', 's6']}
        self.assertEqual(actual, expected)

        # maps to longest seq
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTACA'),
                ('s3', 'ACGTAGAA'),
                ('s4', 'ACGTATAAA'),
                ('s5', 'ACGTCAAAAA'),
                ('s6', 'ACGTCCAAAAA')]

        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = [('s6', 'ACGTCCAAAAA')
                    ], {'s6': ['s1', 's2', 's3', 's4', 's5', 's6']}
        self.assertEqual(actual, expected)

    def test_prefilter_exact_prefixes_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when filtering occurs
        """
        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTAC'),
                ('s3', 'ACGTAG'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCA'),
                ('s6', 'ACGTCC')]

        prefix_length = 5
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = [('s1', 'ACGTAA'), ('s5', 'ACGTCA')], \
            {'s1': ['s1', 's2', 's3', 's4'], 's5': ['s5', 's6']}
        self.assertEqual(actual, expected)

        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1', 'ACGTAA'),
                ('s2', 'ACGTAC'),
                ('s3', 'ACGTAGAAAA'),
                ('s4', 'ACGTAT'),
                ('s5', 'ACGTCA'),
                ('s6', 'ACGTCC')]

        prefix_length = 5
        actual = app._prefilter_exact_prefixes(seqs, prefix_length)
        actual[0].sort()
        expected = [('s3', 'ACGTAGAAAA'), ('s5', 'ACGTCA')], \
            {'s3': ['s1', 's2', 's3', 's4'], 's5': ['s5', 's6']}
        self.assertEqual(actual, expected)

    def test_map_filtered_clusters_to_full_clusters(self):
        """CdHitOtuPicker._map_filtered_clusters_to_full_clusters functions as expected
        """
        # original and mapped full clusters are the same
        app = CdHitOtuPicker(params={})
        filter_map = {'s1': ['s1'], 's2': ['s2'],
                      's3': ['s3'], 's4': ['s4'],
                      's5': ['s5'], 's6': ['s6']}
        clusters = [['s1'], ['s2'], ['s3'], ['s4'], ['s5'], ['s6']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        expected = clusters
        self.assertEqual(actual, expected)

        # original and mapped full clusters are not the same
        filter_map = {'s1': ['s1', 's2', 's3', 's4'], 's5': ['s5', 's6']}
        clusters = [['s1', 's5']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        for e in actual:
            e.sort()
        expected = [['s1', 's2', 's3', 's4', 's5', 's6']]
        self.assertEqual(actual, expected)

        filter_map = {'s1': ['s1', 's2', 's6'],
                      's3': ['s3'], 's5': ['s4', 's5']}
        clusters = [['s1', 's3'], ['s5']]
        actual = app._map_filtered_clusters_to_full_clusters(
            clusters,
            filter_map)
        for e in actual:
            e.sort()
        expected = [['s1', 's2', 's3', 's6'], ['s4', 's5']]
        self.assertEqual(actual, expected)

    def test_call_prefilters_when_requested(self):
        """ CdHitOtuPicker.__call__ prefilters when requested
        """
        # no pre-filtering results in one cluster per sequence as they all
        # differ at their 3' ends
        app = CdHitOtuPicker(params={})
        app = CdHitOtuPicker(params={'Similarity': 0.99})
        self.assertEqual(app(self.tmp_seq_filepath2), dna_seqs_2_result)

        # no pre-filtering results in one cluster per sequence as they are all
        # the same at their 5' ends
        app = CdHitOtuPicker(params={})
        app = CdHitOtuPicker(params={'Similarity': 0.99},)
        self.assertEqual(
            app(self.tmp_seq_filepath2, prefix_prefilter_length=5),
            dna_seqs_2_result_prefilter)


class PickOtusStandaloneFunctions(TestCase):

    """ Tests of stand-alone functions in pick_otus.py """

    def setUp(self):
        """
        """
        self.otu_map1 = {'0': ['seq1', 'seq2', 'seq5'],
                         '1': ['seq3', 'seq4'],
                         '2': ['seq6', 'seq7', 'seq8']}
        self.otu_map2 = {'110': ['0', '2'],
                         '221': ['1']}
        self.otu_map3 = {'a': ['110', '221']}

        self.otu_map1_file = ['0\tseq1\tseq2\tseq5',
                              '1\tseq3\tseq4',
                              '2\tseq6\tseq7\tseq8']
        self.otu_map2_file = ['110\t0\t2',
                              '221\t1']
        self.otu_map3_file = ['a\t110\t221']

        self.failures1 = ['110']
        self.failures2 = ['110\n', '221']
        self.failures3 = ['a']

    def test_expand_failures_one_otu_map(self):
        """expanding failures generated by chained otu picking fns as expected
        """
        expected_f1 = ['0', '2']
        assert_almost_equal(expand_failures(self.failures1, self.otu_map2),
                              expected_f1)
        expected_f2 = ['0', '1', '2']
        assert_almost_equal(expand_failures(self.failures2, self.otu_map2),
                              expected_f2)

    def test_expand_failures_two_otu_maps(self):
        """expanding failures generated by chained otu picking fns as expected
        """
        expected_f1 = ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8']

        actual = expand_failures(self.failures1,
                                 expand_otu_map_seq_ids(self.otu_map2, self.otu_map1))
        assert_almost_equal(actual, expected_f1)

    def test_map_otu_map_files_failures_file_two_otu_maps1(self):
        """map_otu_map_files: correctly maps two otu files and failures
        """
        exp = ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8']
        actual = map_otu_map_files(
            [self.otu_map1_file, self.otu_map2_file],
            self.failures1)
        self.assertEqual(actual, exp)

    def test_map_otu_map_files_failures_file_two_otu_maps2(self):
        """map_otu_map_files: correctly maps two otu files and failures (alt failures)
        """
        exp = ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8', 'seq3', 'seq4']
        actual = map_otu_map_files(
            [self.otu_map1_file, self.otu_map2_file],
            self.failures2)
        self.assertEqual(actual, exp)

    def test_map_otu_map_files_failures_file_three_otu_maps(self):
        """map_otu_map_files: correctly maps three otu files and failures
        """
        exp = ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8', 'seq3', 'seq4']
        actual = map_otu_map_files(
            [self.otu_map1_file, self.otu_map2_file, self.otu_map3_file],
            self.failures3)
        self.assertEqual(actual, exp)

    def test_expand_otu_map_seq_ids_error(self):
        """expand_otu_map_seq_ids: error on missing seq_ids
        """
        self.assertRaises(KeyError, expand_otu_map_seq_ids,
                          self.otu_map3, self.otu_map1)

    def test_expand_otu_map_seq_ids_two(self):
        """expand_otu_map_seq_ids: correctly maps seq_ids from two otu maps
        """
        exp12 = {'110': ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8'],
                 '221': ['seq3', 'seq4']}
        actual12 = expand_otu_map_seq_ids(self.otu_map2, self.otu_map1)
        self.assertEqual(actual12, exp12)

    def test_expand_otu_map_seq_ids_three(self):
        """expand_otu_map_seq_ids: correctly maps seq_ids from three otu maps
        """
        exp123 = {'a': ['seq1', 'seq2', 'seq5', 'seq6',
                        'seq7', 'seq8', 'seq3', 'seq4']}
        actual123 = expand_otu_map_seq_ids(self.otu_map3,
                                           expand_otu_map_seq_ids(self.otu_map2, self.otu_map1))
        self.assertEqual(actual123, exp123)

    def test_map_otu_map_files_two(self):
        """map_otu_map_files: correctly maps two otu files
        """
        exp12 = {'110': ['seq1', 'seq2', 'seq5', 'seq6', 'seq7', 'seq8'],
                 '221': ['seq3', 'seq4']}
        actual12 = map_otu_map_files([self.otu_map1_file, self.otu_map2_file])
        self.assertEqual(exp12, actual12)

    def test_map_otu_map_files_three(self):
        """map_otu_map_files: correctly maps three otu files
        """
        exp123 = {'a': ['seq1', 'seq2', 'seq5', 'seq6',
                        'seq7', 'seq8', 'seq3', 'seq4']}
        actual123 = map_otu_map_files(
            [self.otu_map1_file, self.otu_map2_file, self.otu_map3_file])
        self.assertEqual(exp123, actual123)

        # third 'file' contains mixed tabs and spaces
        actual123 = map_otu_map_files(
            [self.otu_map1_file, self.otu_map2_file, ['a\t110 221']])
        self.assertEqual(exp123, actual123)

dna_seqs_1 = """>cdhit_test_seqs_0 comment fields, not part of sequence identifiers
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
> cdhit_test_seqs_1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>cdhit_test_seqs_2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>cdhit_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>cdhit_test_seqs_4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>cdhit_test_seqs_5
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>cdhit_test_seqs_6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>cdhit_test_seqs_7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>cdhit_test_seqs_8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>cdhit_test_seqs_9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA"""

dna_seqs_result_file_90_exp = """0\tcdhit_test_seqs_0
1\tcdhit_test_seqs_1
2\tcdhit_test_seqs_2
3\tcdhit_test_seqs_3
4\tcdhit_test_seqs_4
5\tcdhit_test_seqs_5
6\tcdhit_test_seqs_6\tcdhit_test_seqs_8
7\tcdhit_test_seqs_7
8\tcdhit_test_seqs_9
"""

dna_seqs_2 = """>cdhit_test_seqs_0 comment fields, not part of sequence identifiers
ACACCCCGGGGGTTTACATTTTTTTTTTTTTTTTTTTTTTTT
>cdhit_test_seqs_1
ACACCCCGGGGGTTTACACCAACATACACCGAGTTGGA
>cdhit_test_seqs_2
ACACCCCGGGGGTTTACGGGGGGGGGGGGGGGGGGGGGGGGGG"""

# results are in length order
dna_seqs_2_result = {0: ['cdhit_test_seqs_2'],
                     1: ['cdhit_test_seqs_0'],
                     2: ['cdhit_test_seqs_1']}

dna_seqs_2_result_prefilter =\
    {0: ['cdhit_test_seqs_0', 'cdhit_test_seqs_1', 'cdhit_test_seqs_2']}


dna_seqs_3 = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA"""

dna_seqs_3_derep = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_1rep some comment1rep
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_1rep2 some comment1rep2
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5 some comment4_again
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment7
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
>uclust_test_seqs_9rep some comment9rep
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
"""

uclustref_query_seqs1 = """>uclust_test_seqs_0 some comment aaa
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_1 some comment bbb
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_2 some comment vv
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment
ACCCACACGGTGGATGCAACAGATCCCATACACCGAGTTGGATGCTTAAGACGCATCGCGTGAGTTTTGCGTCAAGGCT
>uclust_test_seqs_5 some comment
CCGCGGTAGGTGCAACACGTCCCATACAACGGGTTGGAAGGTTAAGACACAACGCGTTAATTTTGTGTCAGGGCA
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA
"""

uclustref_ref_seqs1 = """>ref1 25 random bases appended to uclust_test_seqs_0 and one mismatch
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATATTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCTATAGCAGCCCCAGCGTTTACTTCTA
>ref2 15 random bases prepended to uclust_test_seqs_1 and one mismatch
GCTGCGGCGTCCTGCGCCACGGTGGGTACAACACGTCCACTACATCTGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>ref3 5 random bases prepended and 10 random bases appended to uclust_test_seqs_2
ATAGGCCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACTGCCTGATTCA
>ref4 exact match to uclust_test_seqs_3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
"""

dna_seqs_3_result_file_90_exp = """0\tuclust_test_seqs_0
1\tuclust_test_seqs_1
2\tuclust_test_seqs_2
3\tuclust_test_seqs_3
4\tuclust_test_seqs_4
5\tuclust_test_seqs_5
6\tuclust_test_seqs_6\tuclust_test_seqs_8
7\tuclust_test_seqs_7
8\tuclust_test_seqs_9
"""

dna_seqs_4 = """>uclust_test_seqs_0 comment fields, not part of sequence identifiers
ACACCCCGGGGGTTTACATTTTTTTTTTTTTTTTTTTTTTTT
>uclust_test_seqs_1 blah blah blah
ACACCCCGGGGGTTTACACCAACATACACCGAGTTGGA
>uclust_test_seqs_2 blah blah
ACACCCCGGGGGTTTACGGGGGGGGGGGGGGGGGGGGGGGGGG"""

# results are in length order
dna_seqs_4_result = {0: ['uclust_test_seqs_2'],
                     1: ['uclust_test_seqs_0'],
                     2: ['uclust_test_seqs_1']}

dna_seqs_5 = """>uclust_test_seqs_0 some comment
ACGGTGGCTACAAGACGTCCCATCCAACGGGTTGGATACTTAAGGCACATCACGTCAGTTTTGTGTCAGAGCT
>uclust_test_seqs_0_rc some other comment
AGCTCTGACACAAAACTGACGTGATGTGCCTTAAGTATCCAACCCGTTGGATGGGACGTCTTGTAGCCACCGT
"""

dna_seqs_6 = """>uclust_test_seqs_0 some comment0
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_1 some comment1
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_2 some comment2
CCCCCACGGTGGCAGCAACACGTCACATACAACGGGTTGGATTCTAAAGACAAACCGCGTCAAAGTTGTGTCAGAACT
>uclust_test_seqs_3 some comment3
CCCCACGGTAGCTGCAACACGTCCCATACCACGGGTAGGATGCTAAAGACACATCGGGTCTGTTTTGTGTCAGGGCT
>uclust_test_seqs_4 some comment4
GCCACGGTGGGTACAACACGTCCACTACATCGGCTTGGAAGGTAAAGACACGTCGCGTCAGTATTGCGTCAGGGCT
>uclust_test_seqs_5 some comment4_again
AACCCCCACGGTGGATGCCACACGCCCCATACAAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_6 some comment6
CGCGGTGGCTGCAAGACGTCCCATACAACGGGTTGGATGCTTAAGACACATCGCAACAGTTTTGAGTCAGGGCT
>uclust_test_seqs_7 some comment7
AACCCCCACGGTGGATGCCACACGCCCCATACCAAGGGTAGGATGCTTAAGACACATCGCGTCAGGTTTGTGTCAGGCCT
>uclust_test_seqs_8 some comment8
CGGTGGCTGCAACACGTGGCATACAACGGGTTGGATGCTTAAGACACATCGCCTCAGTTTTGTGTCAGGGCT
>uclust_test_seqs_9 some comment9
GGTGGCTGAAACACATCCCATACAACGGGTTGGATGCTTAAGACACATCGCATCAGTTTTATGTCAGGGGA"""

dna_seqs_4_result_prefilter =\
    {0: ['uclust_test_seqs_0', 'uclust_test_seqs_1', 'uclust_test_seqs_2']}

expected_uc_output =\
    ['# Tab-separated fields:',
     '# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel',
     '# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit',
     '# For C and D types, PctId is average id with seed.',
     '# QueryStart and SeedStart are zero-based relative to start of sequence.',
     '# If minus strand, SeedStart is relative to reverse-complemented seed.',
     'S\t0\t71\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_9\t*',
     'S\t1\t76\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_4\t*',
     'S\t2\t72\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_8\t*',
     'S\t3\t74\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_6\t*',
     'S\t4\t75\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_5\t*',
     'S\t5\t78\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_2\t*',
     'S\t6\t77\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_3\t*',
     'S\t7\t73\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_7\t*',
     'S\t8\t79\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_1\t*',
     'S\t9\t80\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_0\t*',
     'C\t0\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_9\t*',
     'C\t1\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_4\t*',
     'C\t2\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_8\t*',
     'C\t3\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_6\t*',
     'C\t4\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_5\t*',
     'C\t5\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_2\t*',
     'C\t6\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_3\t*',
     'C\t7\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_7\t*',
     'C\t8\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_1\t*',
     'C\t9\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_0\t*']

expected_ref_uc_file =\
    ['# Tab-separated fields:',
     '# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel',
     '# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit', '# For C and D types, PctId is average id with seed.',
     '# QueryStart and SeedStart are zero-based relative to start of sequence.',
     '# If minus strand, SeedStart is relative to reverse-complemented seed.',
     'S\t4\t71\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_9\t*',
     'L\t1\t91\t*\t*\t*\t*\t*\tref2 15 random bases prepended to uclust_test_seqs_1 and one mismatch\t*',
     'H\t1\t76\t98.7\t+\t0\t0\t15I76M\tQiimeExactMatch.uclust_test_seqs_1\tref2 15 random bases prepended to uclust_test_seqs_1 and one mismatch',
     'S\t5\t72\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_8\t*',
     'S\t6\t74\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_6\t*',
     'S\t7\t75\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_5\t*',
     'L\t2\t93\t*\t*\t*\t*\t*\tref3 5 random bases prepended and 10 random bases appended to uclust_test_seqs_2\t*',
     'H\t2\t78\t100.0\t+\t0\t0\t5I78M10I\tQiimeExactMatch.uclust_test_seqs_2\tref3 5 random bases prepended and 10 random bases appended to uclust_test_seqs_2',
     'L\t3\t77\t*\t*\t*\t*\t*\tref4 exact match to uclust_test_seqs_3\t*',
     'H\t3\t77\t100.0\t+\t0\t0\t77M\tQiimeExactMatch.uclust_test_seqs_3\tref4 exact match to uclust_test_seqs_3',
     'L\t0\t98\t*\t*\t*\t*\t*\tref1 25 random bases appended to uclust_test_seqs_0 and one mismatch\t*',
     'H\t0\t73\t98.6\t+\t0\t0\t73M25I\tQiimeExactMatch.uclust_test_seqs_0\tref1 25 random bases appended to uclust_test_seqs_0 and one mismatch',
     'S\t8\t79\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_4\t*',
     'S\t9\t80\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_7\t*',
     'D\t0\t2\t*\t*\t*\t*\t98.6\tref1 25 random bases appended to uclust_test_seqs_0 and one mismatch\t*',
     'D\t1\t2\t*\t*\t*\t*\t98.7\tref2 15 random bases prepended to uclust_test_seqs_1 and one mismatch\t*',
     'D\t2\t2\t*\t*\t*\t*\t100.0\tref3 5 random bases prepended and 10 random bases appended to uclust_test_seqs_2\t*',
     'D\t3\t2\t*\t*\t*\t*\t100.0\tref4 exact match to uclust_test_seqs_3\t*',
     'C\t4\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_9\t*',
     'C\t5\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_8\t*',
     'C\t6\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_6\t*',
     'C\t7\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_5\t*',
     'C\t8\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_4\t*',
     'C\t9\t1\t*\t*\t*\t*\t*\tQiimeExactMatch.uclust_test_seqs_7\t*']

usearch_ref_seqs1 = """>ref1 ecoli sequence
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCA
>EU199232 1 1236 Bacteria/Deltaproteobacteria/Desulfurella - Hippea/uncultured
TACGCGCGGAAATCGAGCGAGATTGGGAACGCAAGTTCCTGAGTATTGCGGCGAACGGGTGAGTAAGACGTGGGTGATCTACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC
>L07864 1 1200 Bacteria/Beta Gammaproteobacteria/Solemya symbiont
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTAATGCATGGGAATCTGCCATATAGTGGGGGACAACTGGGGAAACCCAGGCTAATACCGCATAATCTCTACGGAGGAAAGGCTTC
"""

dna_seqs_usearch = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGT
>Solemya seq
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTA
>usearch_ecoli_seq2
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTCCAT
>Solemya_seq2
GGCTCAGATTGAACGCTGGCGGCATGCCTAACACATGCAAGTCGAACGGTAACAGGCGGAGCTTGCTCTGCGCTGACGAGTGGCGGACGGGTGAGTATCAAG
>chimera
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACCCCTAGGGTGGGAATAACCCGGGGAAACCCGGGCTAATACCGAATAAGACCACAGGAGGCGACTCCAGAGGGTCAAAGGGAGCCTTGGCCTCCCCC
"""

dna_seqs_usearch_97perc_id = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_1bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_2bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGCGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
"""

dna_seqs_usearch_97perc_id_len_diff = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_1bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGA
>usearch_ecoli_seq_2bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGCGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGAAA
"""

dna_seqs_usearch_97perc_dups = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_1bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGA
>usearch_ecoli_seq_2bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGCGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGAAA
>usearch_ecoli_seq_1bp_change_dup1
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGA
>usearch_ecoli_seq_1bp_change_dup2
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGA
"""

dna_seqs_usearch_97perc_id_rc = """>usearch_ecoli_seq
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_1bp_change
CGCGTGTATGAAGAAGGCCTACGGGTTGTAAAGTACTTTCAGCGGGGAGGAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAG
>usearch_ecoli_seq_2bp_change_rc
CTTCTTCTGCGGGTAACGTCAATGAGCAAAGGTATTAACTTTACTCCCTCCGCCCCGCTGAAAGTACTTTACAACCCGTAGGCCTTCTTCATACACGCG
"""

dna_seqs_rc_single_seq = """>usearch_ecoli_seq_2bp_change_rc
CTTCTTCTGCGGGTAACGTCAATGAGCAAAGGTATTAACTTTACTCCCTCCGCCCCGCTGAAAGTACTTTACAACCCGTAGGCCTTCTTCATACACGCG
"""

# run unit tests if run from command-line
if __name__ == '__main__':
    main()
