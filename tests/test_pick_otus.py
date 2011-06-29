#!/usr/bin/env python

"""Tests of code for OTU picking"""

__author__ = "Kyle Bittinger, Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Kyle Bittinger", "Greg Caporaso", "Rob Knight", "Jens Reeder","William Walters"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from os import remove

from cogent.util.misc import create_dir
from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename, load_qiime_config
from cogent.util.misc import remove_files
from cogent import DNA
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.pick_otus import (CdHitOtuPicker, OtuPicker,
    MothurOtuPicker, PrefixSuffixOtuPicker, TrieOtuPicker, BlastOtuPicker,
    expand_otu_map_seq_ids, map_otu_map_files, UclustOtuPicker,
    UclustReferenceOtuPicker, expand_failures)


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
        seqs = [('s1 comment1','ACCTTGTTACTTT'),  # three copies
                ('s2 comment2','ACCTTGTTACTTTC'), # one copy
                ('s3 comment3','ACCTTGTTACTTTCC'),# two copies
                ('s4 comment4','ACCTTGTTACTTT'),
                ('s5 comment5','ACCTTGTTACTTTCC'),
                ('s6 comment6','ACCTTGTTACTTT')]
        expected0 = [('QiimeExactMatch.s1','ACCTTGTTACTTT'),
                     ('QiimeExactMatch.s2','ACCTTGTTACTTTC'),
                     ('QiimeExactMatch.s3','ACCTTGTTACTTTCC')]
        expected1 = {'QiimeExactMatch.s1':['s1','s4','s6'],
                     'QiimeExactMatch.s2':['s2'],
                     'QiimeExactMatch.s3':['s3','s5']}
        expected = (expected0,expected1)
        p = OtuPicker({})
        actual = p._prefilter_exact_matches(seqs)
        self.assertEqual(actual,expected)


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
        expected_otus = {0: ['cccccc'], 1: ['bbbbbb'], 2: ['aaaaaa']}
        self.assertEqual(observed_otus, expected_otus)

    def test_call_low_similarity(self):
        app = MothurOtuPicker({'Similarity': 0.35})
        observed_otus = app(self.small_seq_path)
        expected_otus = {0: ['bbbbbb', 'cccccc'], 1: ['aaaaaa']}
        self.assertEqual(observed_otus, expected_otus)

    def test_call_nearest_neighbor(self):
        app = MothurOtuPicker({'Algorithm': 'nearest', 'Similarity': 0.35})
        observed_otus = app(self.small_seq_path)
        expected_otus = {0: ['bbbbbb', 'cccccc'], 1: ['aaaaaa']}
        self.assertEqual(observed_otus, expected_otus)

class BlastOtuPickerTests(TestCase):
    """ Tests of the blast-based otu picker """
    
    def setUp(self):
        """
        """
        self.otu_picker = BlastOtuPicker({'max_e_value':1e-3})
        self.seqs = [\
         ('s0  some description','CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'),\
         ('s1','TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),\
         ('s2','TGCAGCTTGAGCCACAGGAGAGAGCCTTC'),\
         ('s3','TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),\
         ('s4','ACCGATGAGATATTAGCACAGGGGAATTAGAACCA'),\
         ('s5','TGTCGAGAGTGAGATGAGATGAGAACA'),\
         ('s6','ACGTATTTTAATTTGGCATGGT'),\
         ('s7','TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'),\
        ]
        
        self.ref_seqs = [\
         ('ref1','TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC'),\
         ('ref2','ACCGATGAGATATTAGCACAGGGGAATTAGAACCA'),\
         ('ref3','TGTCGAGAGTGAGATGAGATGAGAACA'),\
         ('ref4','ACGTATTTTAATGGGGCATGGT'),\
        ]
        
        self.ref_seqs_rc = [\
         ('ref1',DNA.rc('TGCAGCTTGAGCCACAGGAGAGAGAGAGCTTC')),\
         ('ref2',DNA.rc('ACCGATGAGATATTAGCACAGGGGAATTAGAACCA')),\
         ('ref3',DNA.rc('TGTCGAGAGTGAGATGAGATGAGAACA')),\
         ('ref4',DNA.rc('ACGTATTTTAATGGGGCATGGT')),\
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
        
        self._files_to_remove = \
         [self.seqs_fp,self.reference_seqs_fp,self.reference_seqs_rc_fp]
        
    def tearDown(self):
        """
        """
        remove_files(self._files_to_remove)
        
    def test_blast_seqs(self):
        """ blast_seqs: functions as expected
        """
        blast_db, db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._files_to_remove += db_files_to_remove
        self.otu_picker.blast_db = blast_db
        
        actual_clusters, actual_failures =\
         self.otu_picker._blast_seqs(self.seqs)
         
        for v in actual_clusters.values(): v.sort()
        actual_failures.sort()
        
        expected_clusters = {'ref1':['s1','s2','s3'],'ref2':['s4'],\
                    'ref3':['s5']}
        expected_failures = ['s0','s6','s7']
        
        self.assertEqual(actual_clusters,expected_clusters)
        self.assertEqual(actual_failures,expected_failures)
        
    def test_update_cluster_map(self):
        """update_cluster_map: functions as expected
        """
        # nothing in original cm
        cm = {}
        new_cm = {'c1':['1','2','5'],'c2':['4','3']}
        expected = new_cm
        actual = self.otu_picker._update_cluster_map(cm,new_cm)
        self.assertEqual(actual,expected)
        
        # no new clusters
        cm = {'c1':['1','2','5'],'c2':['4','3']}
        new_cm = {}
        expected = cm
        actual = self.otu_picker._update_cluster_map(cm,new_cm)
        self.assertEqual(actual,expected)
        
        # overlapping clusters
        cm = {'c1':['1','2','5'],'c2':['4','3']}
        new_cm = {'c1':['8'],'c2':['10','14'],'3':['42']}
        expected = {'c1':['1','2','5','8'],'c2':['4','3','10','14'],'3':['42']}
        actual = self.otu_picker._update_cluster_map(cm,new_cm)
        self.assertEqual(actual,expected)
        
        # no duplicate seq_id checking 
        cm = {'c1':['1']}
        new_cm = cm
        expected = {'c1':['1','1']}
        actual = self.otu_picker._update_cluster_map(cm,new_cm)
        self.assertEqual(actual,expected)
        
        # no clusters at all
        actual = self.otu_picker._update_cluster_map({},{})
        self.assertEqual(actual,{})
        
        
    def test_call(self):
        """BLAST OTU Picker functions as expected
        """
        
        expected = {'ref1':['s3','s2','s1'],\
                    'ref2':['s4'],\
                    'ref3':['s5']}
        actual = self.otu_picker(self.seqs_fp,\
            refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual,expected)
        
    def test_call_alt_min_aligned_length(self):
        """BLAST OTU picker handles alt min_aligned_percent values """
        # first 12 bases match perfect, and no alignment from there
        seqs = [('s1','TGCAGCTTGAGCGTTGTTACCGCTTT')]
        ref_seqs = [\
         ('r1','TGCAGCTTGAGCCACGCCGAATAGCCGAGTTTGACCGGGCCCAGGAGGAGAGAGAGAGCTTC')]
        
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
        otu_picker = BlastOtuPicker({'max_e_value':1e-3,
                                     'min_aligned_percent':0.10})
        expected = {'r1':['s1']}
        actual = otu_picker(seqs_fp,\
            refseqs_fp=reference_seqs_fp)
        self.assertEqual(actual,expected)
        
        # with min_aligned_percent s1 doesn't match r1
        otu_picker = BlastOtuPicker({'max_e_value':1e-3,
                                     'min_aligned_percent':0.50})
        expected = {}
        actual = otu_picker(seqs_fp,\
            refseqs_fp=reference_seqs_fp)
        self.assertEqual(actual,expected)
         
        
    def test_call_rc(self):
        """BLAST OTU picker: RC seqs cluster to same OTU as forward orientation
        """
        
        expected = {'ref1':['s3','s2','s1'],\
                    'ref2':['s4'],\
                    'ref3':['s5']}
        actual = self.otu_picker(self.seqs_fp,\
            refseqs_fp=self.reference_seqs_rc_fp)
        self.assertEqual(actual,expected)
        
        
    def test_call_alt_params(self):
        """BLAST OTU Picker functions as expected with alt params
        """
        otu_picker = BlastOtuPicker({'max_e_value':1e-30})
        expected = {}
        actual = otu_picker(self.seqs_fp,\
            refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual,expected)
        
        self.otu_picker = BlastOtuPicker({'max_e_value':1e-3,'Similarity':0.90})
        expected_90 = {'ref1':['s3','s2','s1'],\
                    'ref2':['s4'],\
                    'ref3':['s5'],\
                    'ref4':['s6']}
        actual = self.otu_picker(self.seqs_fp,\
            refseqs_fp=self.reference_seqs_fp)
        self.assertEqual(actual,expected_90)
        
    def test_call_preexisting_blast_db(self):
        """BLAST OTU Picker functions w preexisting blast db
        """
        blast_db, db_files_to_remove = \
            build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._files_to_remove += db_files_to_remove
        expected = {'ref1':['s3','s2','s1'],\
                    'ref2':['s4'],\
                    'ref3':['s5']}
        actual = self.otu_picker(self.seqs_fp,blast_db=blast_db)
        self.assertEqual(actual,expected)
        
    def test_call_multiple_blast_runs(self):
        """BLAST OTU Picker not affected by alt SeqsPerBlastRun
        """
        expected = {'ref1':['s1','s2','s3'],\
                    'ref2':['s4'],\
                    'ref3':['s5']}
        for v in expected.values():
            v.sort()
        for SeqsPerBlastRun in [1,2,4,6,7,8,100]:
            self.otu_picker.Params['seqs_per_blast_run'] \
             = SeqsPerBlastRun
            actual = self.otu_picker(self.seqs_fp,\
                refseqs_fp=self.reference_seqs_fp)
            for v in actual.values():
                v.sort()
            self.assertEqual(actual,expected)

class PrefixSuffixOtuPickerTests(TestCase):
    """ Tests of the prefix/suffix-based OTU picker """
    
    def setUp(self):
        """
        """
        self.otu_picker = PrefixSuffixOtuPicker({})
        self.seqs = [\
         ('s1  some description','ACGTAATGGT'),\
         ('s2','ATTTAATGGT'),\
         ('s3','ACGTAATTTT'),\
         ('s4','AAATAAAAA'),\
         ('s5','ACGTTGGT'),\
         ('s6','ACGTATTTTAATTTGGCATGGT'),\
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
        expected = {3:['s1','s5','s6'],\
                    1:['s2'],\
                    0:['s3'],\
                    2:['s4']}
        actual = self.otu_picker(self.small_seq_path,\
            prefix_length=4,suffix_length=4)
        self.assertEqual(actual,expected)
        
    def test_call_extra_long_lengths(self):
        """Prefix/suffix OTU Picker functions as expected
        """
        seqs = [\
         ('s1  some description',
          'ACGTAATGGTCCCCCCCCCGGGGGGGGCCCCCCGGG'),\
         ('s2','ATTTAATGGT'),\
         ('s3','ACGTAATTTT'),\
         ('s4','AAATAAAAA'),\
         ('s5','ACGTTGGT'),\
         ('s6','ACGTATTTTAATTTGGCATGGT'),\
         ('s7','ACGTATTTTAATTTGGCATGG'),\
         ('s1_dup',
          'ACGTAATGGTCCCCCCCCCGGGGGGGGCCCCCCGGG'),\
         ('s2_dup','ATTTAATGGT'),\
        ]
        seq_path = get_tmp_filename(
            prefix='PrefixSuffixOtuPickerTest_', suffix='.fasta')
        self._files_to_remove.append(seq_path)
        f = open(seq_path, 'w')
        f.write('\n'.join(['>%s\n%s' % s for s in seqs]))
        f.close()
        expected = {1:['s1','s1_dup'],\
                    2:['s2','s2_dup'],\
                    3:['s3'],\
                    4:['s4'],\
                    5:['s5'],\
                    6:['s6'],
                    7:['s7']}
        
        # long prefix collapses identical sequences
        actual = self.otu_picker(seq_path,\
            prefix_length=400,suffix_length=0)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        self.assertEqualItems(actual_clusters,expected_clusters)
        
        # long suffixes collapses identical sequences
        actual = self.otu_picker(seq_path,\
            prefix_length=0,suffix_length=400)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        self.assertEqualItems(actual_clusters,expected_clusters)
        
        # long prefix and suffixes collapses identical sequences
        actual = self.otu_picker(seq_path,\
            prefix_length=400,suffix_length=400)
        actual_clusters = actual.values()
        expected_clusters = expected.values()
        self.assertEqualItems(actual_clusters,expected_clusters)
        
        
    def test_collapse_exact_matches_prefix_and_suffix(self):
        """Prefix/suffix: collapse_exact_matches fns with pref/suf len > 0
        """
        expected = [['s1','s5','s6'],['s2'],['s3'],['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,4,4)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)
        
        expected = [['s1','s2','s3','s5','s6'],['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,1,1)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)   
             
    def test_collapse_exact_matches_prefix_zero(self):
        """Prefix/suffix: collapse_exact_matches fns with prefix len = 0
        """
        expected = [['s1','s2','s5','s6'],['s3'],['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,0,4)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)
        
        expected = [['s1','s2','s3','s5','s6'],['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,0,1)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)
             
    def test_collapse_exact_matches_suffix_zero(self):
        """Prefix/suffix: collapse_exact_matches fns with suffix len = 0
        """
        expected = [['s1','s3','s5','s6'],['s2'],['s4']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,4,0)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)
        
        expected = [['s1','s2','s3','s4','s5','s6']]
        actual = self.otu_picker._collapse_exact_matches(self.seqs,1,0)
        actual.sort()
        expected.sort()
        self.assertEqual(actual,expected)
    
    def test_build_seq_hash(self):
        """ """
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',0,0),'')
            
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',2,2),'ATGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATTACGT',2,1),'ATT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',1,2),'AGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',1,1),'AT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',4,3),'ATGTCGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',3,4),'ATGACGT')
            
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',4,4),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',5,3),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',8,0),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',3,5),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',0,8),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',4,5),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',5,4),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',300,0),'ATGTACGT')
        self.assertEqual(self.otu_picker._build_seq_hash(\
            'ATGTACGT',0,300),'ATGTACGT')

class TrieOtuPickerTests(TestCase):
    """ Tests of the Trie-based OTU picker """
    
    def setUp(self):
        """
        """
        self.otu_picker = TrieOtuPicker({})
        self.otu_picker_rev = TrieOtuPicker({'Reverse':True})
        seqs = [\
         ('s1 some description','ACGTAATGGT'),\
         ('s2','ACGTATTTTAATTTGGCATGGT'),\
         ('s3','ACGTAAT'),\
         ('s4','ACGTA'),\
         ('s5','ATTTAATGGT'),\
         ('s6','ATTTAAT'),\
         ('s7','AAATAAAAA')
        ]
        seqs_rev = [\
         ('s1 some description','TGGTAATGCA'),\
         ('s2','TGGTACGGTTTAATTTTATGCA'),\
         ('s3','TAATGCA'),\
         ('s4','ATGCA'),\
         ('s5','TGGTAATTTA'),\
         ('s6','TAATTTA'),\
         ('s7','AAAAATAAA')
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
        expected = {0:['s2'],\
                    1:['s3','s4','s1'],\
                    2:['s7'],\
                    3:['s6','s5']}
        actual = self.otu_picker(self.small_seq_path)
        self.assertEqual(actual,expected)
        
    def test_call_reverse(self):
        """Trie OTU Picker functions as expected with the 'Reverse' option
        """
        expected = {0:['s2'],\
                    1:['s3','s4','s1'],\
                    2:['s7'],\
                    3:['s6','s5']}
        actual = self.otu_picker_rev(self.small_seq_path_rev)
        self.assertEqual(actual,expected)
        
        

class UclustOtuPickerTests(TestCase):
    """ Tests of the uclust-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(\
         prefix='UclustOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1,'w')
        seq_file.write(dna_seqs_3)
        seq_file.close()        
        
        self.tmp_seq_filepath2 = get_tmp_filename(\
         prefix='UclustOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2,'w')
        seq_file.write(dna_seqs_4)
        seq_file.close()
        
        self.tmp_seq_filepath3 = get_tmp_filename(\
         prefix='UclustOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath3,'w')
        seq_file.write(dna_seqs_5)
        seq_file.close()
        
        self.tmp_seq_filepath4 = get_tmp_filename(\
         prefix='UclustOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath4,'w')
        seq_file.write(dna_seqs_6)
        seq_file.close()
        
        self._files_to_remove =\
         [self.tmp_seq_filepath1, self.tmp_seq_filepath2, 
          self.tmp_seq_filepath3, self.tmp_seq_filepath4]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        
    def seqs_to_temp_fasta(self,seqs):
        """ """
        fp = get_tmp_filename(
         prefix='UclustReferenceOtuPickerTest_',
         suffix='.fasta')
        seq_file = open(fp,'w')
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
          UclustOtuPicker(params={'Similarity':0.90,
                                  'save_uc_files':True,
                                  'prefilter_identical_sequences':True})
        result_w_collapse_identical = \
         app_w_collapse_identical(self.tmp_seq_filepath4).values()
        result_w_collapse_identical.sort()
        
        # generate result excluding prefilter
        app_wo_collapse_identical =\
          UclustOtuPicker(params={'Similarity':0.90,
                                  'save_uc_files':False,
                                  'prefilter_identical_sequences':False})
        result_wo_collapse_identical = \
         app_wo_collapse_identical(self.tmp_seq_filepath4).values()
        result_wo_collapse_identical.sort()
        
        self.assertEqual(result_w_collapse_identical,
                         result_wo_collapse_identical)
        
    def test_toggle_suppress_sort(self):
        """UclustOtuPicker: togging suppress sort functions as expected
        """
        seqs = [('s1','ACCTTGTTACTTT'),  # three copies
                ('s2','ACCTTGTTACTTTC'), # one copy
                ('s3','ACCTTGTTACTTTCC'),# two copies
                ('s4','ACCTTGTTACTTT'),
                ('s5','ACCTTGTTACTTTCC'),
                ('s6','ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        
        # no abundance sorting and uclust's sorting enabled 
        # so length-based sorting
        app = UclustOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(seqs_fp)
        exp = {0:['s3','s5','s2','s1','s4','s6']}
        self.assertEqual(obs,exp)
        
        # no abundance sorting and uclust's sorting enabled 
        # so no sorting at all
        app = UclustOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':True,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(seqs_fp)
        exp = {0:['s1','s4','s6','s2','s3','s5']}
        self.assertEqual(obs,exp)
        
    def test_abundance_sort(self):
        """UclustOtuPicker: abundance sort functions as expected
        """
        #enable abundance sorting with suppress sort = False (it gets
        # set to True internally, otherwise uclust's length sort would
        # override the abundance sorting)
        seqs = [('s1 comment1','ACCTTGTTACTTT'),  # three copies
                ('s2 comment2','ACCTTGTTACTTTC'), # one copy
                ('s3 comment3','ACCTTGTTACTTTCC'),# two copies
                ('s4 comment4','ACCTTGTTACTTT'),
                ('s5 comment5','ACCTTGTTACTTTCC'),
                ('s6 comment6','ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        
        # abundance sorting changes order 
        app = UclustOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':False,
                                      'presort_by_abundance':True,
                                      'save_uc_files':False})
        obs = app(seqs_fp)
        exp = {0:['s1','s4','s6','s3','s5','s2']}
        self.assertEqual(obs,exp)
        
        # abundance sorting changes order -- same results with suppress_sort =
        # True b/c (it gets set to True to when presorting by abundance)
        app = UclustOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':True,
                                      'presort_by_abundance':True,
                                      'save_uc_files':False})
        obs = app(seqs_fp)
        exp = {0:['s1','s4','s6','s3','s5','s2']}
        self.assertEqual(obs,exp)

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
        
        app = UclustOtuPicker(params={'save_uc_files':False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
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
        
        

        
        app = UclustOtuPicker(params={'save_uc_files':True,
                                      'output_dir':self.temp_dir})
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
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
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
                        ['uclust_test_seqs_6','uclust_test_seqs_8'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
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
                        ['uclust_test_seqs_6','uclust_test_seqs_8'],
                        ['uclust_test_seqs_7'],
                        ['uclust_test_seqs_9']]

        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'new_cluster_identifier':'my_otu_',
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath1)
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
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

        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'suppress_sort':True,
                                      'optimal':True,
                                      'enable_rev_strand_matching':True,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath2)
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
        # The relation between otu ids and clusters is abitrary, and 
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)
        
    def test_call_rev_matching(self):
        """UclustOtuPicker.__call__ handles reverse strand matching
        """
        exp_otu_ids = range(2)
        exp_clusters = [['uclust_test_seqs_0'],['uclust_test_seqs_0_rc']]
        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath3)
        obs_otu_ids = obs.keys()
        obs_otu_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
        # The relation between otu ids and clusters is abitrary, and 
        # is not stable due to use of dicts when parsing clusters -- therefore
        # just checks that we have the expected group of each
        self.assertEqual(obs_otu_ids, exp_otu_ids)
        self.assertEqual(obs_clusters, exp_clusters)
        
        exp = {0: ['uclust_test_seqs_0','uclust_test_seqs_0_rc']}
        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'enable_rev_strand_matching':True,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath3)
        self.assertEqual(obs, exp)
        
    def test_call_output_to_file(self):
        """UclustHitOtuPicker.__call__ output to file functions as expected
        """
        
        tmp_result_filepath = get_tmp_filename(\
         prefix='UclustOtuPickerTest.test_call_output_to_file_',\
         suffix='.txt')
        
        app = UclustOtuPicker(params={'Similarity':0.90,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath1,result_path=tmp_result_filepath)
        
        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in 
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)

        exp_otu_ids = map(str,range(9))
        exp_clusters = [['uclust_test_seqs_0'],
                        ['uclust_test_seqs_1'],
                        ['uclust_test_seqs_2'],
                        ['uclust_test_seqs_3'],
                        ['uclust_test_seqs_4'],
                        ['uclust_test_seqs_5'],
                        ['uclust_test_seqs_6','uclust_test_seqs_8'],
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
        self.assertEqual(obs,None)
        
    def test_call_log_file(self):
        """UclustOtuPicker.__call__ writes log when expected
        """
        
        tmp_log_filepath = get_tmp_filename(\
         prefix='UclustOtuPickerTest.test_call_output_to_file_l_',\
         suffix='.txt')
        tmp_result_filepath = get_tmp_filename(\
         prefix='UclustOtuPickerTest.test_call_output_to_file_r_',\
         suffix='.txt')
        
        app = UclustOtuPicker(params={'Similarity':0.99,
                                      'save_uc_files':False})
        obs = app(self.tmp_seq_filepath1,\
         result_path=tmp_result_filepath,log_path=tmp_log_filepath)
        
        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        # remove the temp files before running the test, so in 
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)
        
        log_file_99_exp = ["UclustOtuPicker parameters:",
         "Similarity:0.99","Application:uclust",
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
        self.assertEqualItems(log_file_str.split('\n'), log_file_99_exp)
        
        
    def test_map_filtered_clusters_to_full_clusters(self):
        """UclustOtuPicker._map_filtered_clusters_to_full_clusters functions as expected
        """
        # original and mapped full clusters are the same
        app = UclustOtuPicker(params={})
        filter_map = {'s1':['s1'], 's2':['s2'], \
                      's3':['s3'], 's4':['s4'], \
                      's5':['s5'], 's6':['s6']}
        clusters = [['s1'], ['s2'], ['s3'], ['s4'], ['s5'], ['s6']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        expected = clusters
        self.assertEqual(actual,expected)
        
        # original and mapped full clusters are not the same
        filter_map = {'s1':['s1','s2','s3','s4'],'s5':['s5','s6']}
        clusters = [['s1','s5']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        for e in actual: e.sort()
        expected = [['s1','s2','s3','s4','s5','s6']]
        self.assertEqual(actual,expected)
        
        filter_map = {'s1':['s1','s2','s6'],'s3':['s3'],'s5':['s4','s5']}
        clusters = [['s1','s3'],['s5']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        for e in actual: e.sort()
        expected = [['s1','s2','s3','s6'],['s4','s5']]
        self.assertEqual(actual,expected)
        
class UclustReferenceOtuPickerTests(TestCase):
    """ Tests of the uclust reference-based OTU picker """
    
    def setUp(self):
        """ """
        self.temp_dir = load_qiime_config()['temp_dir']
        self.tmp_seq_filepath1 = get_tmp_filename(
         prefix='UclustReferenceOtuPickerTest_',
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1,'w')
        seq_file.write(uclustref_query_seqs1)
        seq_file.close()        
        
        self.temp_ref_filepath1 = get_tmp_filename(
         prefix='UclustReferenceOtuPickerTest_',
         suffix='.fasta')
        ref_file = open(self.temp_ref_filepath1,'w')
        ref_file.write(uclustref_ref_seqs1)
        ref_file.close()
        
        self._files_to_remove =\
         [self.tmp_seq_filepath1,
          self.temp_ref_filepath1]
        
    def tearDown(self):
        remove_files(self._files_to_remove)
        
    def seqs_to_temp_fasta(self,seqs):
        """ """
        fp = get_tmp_filename(
         prefix='UclustReferenceOtuPickerTest_',
         suffix='.fasta')
        seq_file = open(fp,'w')
        self._files_to_remove.append(fp)
        for s in seqs:
            seq_file.write('>%s\n%s\n' % s)
        seq_file.close()
        return fp
        
    def test_toggle_suppress_sort(self):
        """UclustReferenceOtuPicker: togging suppress sort functions as expected
        """
        seqs = [('s1 comment1','ACCTTGTTACTTT'),  # three copies
                ('s2 comment2','ACCTTGTTACTTTC'), # one copy
                ('s3 comment3','ACCTTGTTACTTTCC'),# two copies
                ('s4 comment4','ACCTTGTTACTTT'),
                ('s5 comment5','ACCTTGTTACTTTCC'),
                ('s6 comment6','ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        ref_seqs = [('r1 blah','ACCTTGTTACTTT')]
        ref_seqs_fp = self.seqs_to_temp_fasta(ref_seqs)
        
        # no abundance sorting and uclust's sorting enabled 
        # so length-based sorting
        app = UclustReferenceOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(seqs_fp,ref_seqs_fp)
        exp = {'r1':['s3','s5','s2','s1','s4','s6']}
        self.assertEqual(obs,exp)
        
        # no abundance sorting and uclust's sorting enabled 
        # so no sorting at all
        app = UclustReferenceOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':True,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = app(seqs_fp,ref_seqs_fp)
        exp = {'r1':['s1','s4','s6','s2','s3','s5']}
        self.assertEqual(obs,exp)
        
    def test_abundance_sort(self):
        """UclustReferenceOtuPicker: abundance sort functions as expected
        """
        #enable abundance sorting with suppress sort = False (it gets
        # set to True internally, otherwise uclust's length sort would
        # override the abundance sorting)
        seqs = [('s1 comment1','ACCTTGTTACTTT'),  # three copies
                ('s2 comment2','ACCTTGTTACTTTC'), # one copy
                ('s3 comment3','ACCTTGTTACTTTCC'),# two copies
                ('s4 comment4','ACCTTGTTACTTT'),
                ('s5 comment5','ACCTTGTTACTTTCC'),
                ('s6 comment6','ACCTTGTTACTTT')]
        seqs_fp = self.seqs_to_temp_fasta(seqs)
        ref_seqs = [('r1 blah','ACCTTGTTACTTT')]
        ref_seqs_fp = self.seqs_to_temp_fasta(ref_seqs)
        
        # abundance sorting changes order 
        app = UclustReferenceOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':False,
                                      'presort_by_abundance':True,
                                      'save_uc_files':False})
        obs = app(seqs_fp,ref_seqs_fp)
        exp = {'r1':['s1','s4','s6','s3','s5','s2']}
        self.assertEqual(obs,exp)
        
        # abundance sorting changes order -- same results with suppress_sort =
        # True b/c (it gets set to True to when presorting by abundance)
        app = UclustReferenceOtuPicker(params={'Similarity':0.80,
                                      'enable_rev_strand_matching':False,
                                      'suppress_sort':True,
                                      'presort_by_abundance':True,
                                      'save_uc_files':False})
        obs = app(seqs_fp,ref_seqs_fp)
        exp = {'r1':['s1','s4','s6','s3','s5','s2']}
        self.assertEqual(obs,exp)
        
    def test_toggle_suppress_new_clusters(self):
        """UclustReferenceOtuPicker: toggle suppress new clusters 
        """
        seqs = [('s1 a','ACCTTGTTACTTT'),
                ('s2 bb','ACCTAGTTACTTT'),
                ('s3 c  c','TTGCGTAACGTTTGAC')]
        ref_seqs = [
                ('r1 d','ACCTCGTTACTTT')]
        # these seqs should match at 0.90, but don't -- I can confirm this
        # running uclust directly, and have contacted Robert Edgar for
        # clarification
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':True,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        exp = {'r1':['s1','s2']}
        self.assertEqual(obs,exp)
        
        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':False,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        exp = {'r1':['s1','s2'],'new_42':['s3']}
        self.assertEqual(obs,exp)
        
        
    def test_toggle_collapse_identical_sequences_prefilter_w_new_clusters(self):
        """UclustReferenceOtuPicker: ident. seqs prefilter fns w new clusters
        """
        # s4 == s2 and s3 == s5
        seqs = [('s1 a','ACCTTGTTACTTT'),
                ('s2 bb','ACCTAGTTACTTT'),
                ('s4 bb','ACCTAGTTACTTT'),
                ('s3 c  c','TTGCGTAACGTTTGAC'),
                ('s5 c  c','TTGCGTAACGTTTGAC')]
        ref_seqs = [
                ('r1 d','ACCTCGTTACTTT')]
        exp = {'r1':['s2','s4','s1'],
               'new_42':['s3','s5']}
               
        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':False,
                                       'save_uc_files':False,
                                       'prefilter_identical_sequences':False})
        obs_no_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        self.assertEqual(obs_no_prefilter,exp)
        
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':False,
                                       'save_uc_files':False,
                                       'prefilter_identical_sequences':True})
        obs_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        self.assertEqual(obs_prefilter,exp)
        
        # a little paranoia never hurt anyone
        self.assertEqual(obs_prefilter,obs_no_prefilter)

    def test_toggle_collapse_identical_sequences_prefilter_wo_new_clusters(self):
        """UclustReferenceOtuPicker: ident. seqs prefilter fns wo new clusters
        """
        # s4 == s2 and s3 == s5
        seqs = [('s1 a','ACCTTGTTACTTT'),
                ('s2 bb','ACCTAGTTACTTT'),
                ('s4 bb','ACCTAGTTACTTT'),
                ('s3 c  c','TTGCGTAACGTTTGAC'),
                ('s5 c  c','TTGCGTAACGTTTGAC')]
        ref_seqs = [
                ('r1 d','ACCTCGTTACTTT')]
        
        exp = {'r1':['s2','s4','s1']}
        
        # add seq that clusters independently
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':True,
                                       'save_uc_files':False,
                                       'prefilter_identical_sequences':False})
        fail_path_no_prefilter = get_tmp_filename(
         prefix='UclustRefOtuPickerFailures',suffix='.txt')
        self._files_to_remove.append(fail_path_no_prefilter)
        obs_no_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                              self.seqs_to_temp_fasta(ref_seqs),
                              failure_path=fail_path_no_prefilter,
                              HALT_EXEC=False)
        self.assertEqual(obs_no_prefilter,exp)
        self.assertEqual(open(fail_path_no_prefilter).read(),
                         "s3\ns5")
        
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':True,
                                       'save_uc_files':False,
                                       'prefilter_identical_sequences':True})
        fail_path_prefilter = get_tmp_filename(
         prefix='UclustRefOtuPickerFailures',suffix='.txt')
        self._files_to_remove.append(fail_path_prefilter)
        obs_prefilter = uc(self.seqs_to_temp_fasta(seqs),
                           self.seqs_to_temp_fasta(ref_seqs),
                           failure_path=fail_path_prefilter,
                           HALT_EXEC=False)
        self.assertEqual(obs_prefilter,exp)
        self.assertEqual(open(fail_path_prefilter).read(),
                         "s3\ns5")
        
        # a little paranoia never hurt anyone
        self.assertEqual(obs_prefilter,obs_no_prefilter)
        self.assertEqual(open(fail_path_prefilter).read(),
                         open(fail_path_no_prefilter).read())




    def test_varied_similarity(self):
        """UclustReferenceOtuPicker: varying similarity affects clustering
        """
        seqs = [('s1','ACCTTGTTACTTT'),
                ('s2','ACCTAGTTACTTT')]
        ref_seqs = [
                ('r1','ACCTCGTTACTTT')]
        # these seqs should match at 0.90, but don't -- I can confirm this
        # running uclust directly, and have contacted Robert Edgar for
        # clarification
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'suppress_new_clusters':False,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        exp = {'r1':['s1','s2']}
        self.assertEqual(obs,exp)
        
        # set similarity to 100%
        uc = UclustReferenceOtuPicker({'Similarity':1.0,
                                       'suppress_new_clusters':False,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of 
        # iteration over a dict
        exp1 = {'QiimeOTU1':['s1'],'QiimeOTU2':['s2']}
        exp2 = {'QiimeOTU2':['s1'],'QiimeOTU1':['s2']}
        self.assertTrue(obs == exp1 or obs == exp2)
 
    def test_toggle_rev_strand_matching(self):
        """UclustReferenceOtuPicker: toggle rev strand matching
        """
        # s3 and s4 are rc of one another
        seqs = [('s1','ACCTTGTTACTTT'),
                ('s2','ACCTAGTTACTTT'),
                ('s3','TTGCGTAACGTTTGAC'),
                ('s4','GTCAAACGTTACGCAA')]
        ref_seqs = [
                ('r1','ACCTCGTTACTTT')]
        
        # rev strand matching disabled
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'enable_rev_strand_matching':False,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        exp = {'r1':['s1','s2'],'new_42':['s3'],'new_43':['s4']}
        self.assertEqual(obs,exp)
        
        # enable rev strand matching
        uc = UclustReferenceOtuPicker({'Similarity':0.80,
                                       'new_cluster_identifier':'new_',
                                       'next_new_cluster_number':42,
                                       'enable_rev_strand_matching':True,
                                       'save_uc_files':False})
        obs = uc(self.seqs_to_temp_fasta(seqs),
                 self.seqs_to_temp_fasta(ref_seqs),HALT_EXEC=False)
        exp = {'r1':['s1','s2'],'new_42':['s3','s4']}
        self.assertEqual(obs,exp)
        
    def test_call_log_file(self):
        """UclustReferenceOtuPicker.__call__ writes log when expected
        """
        tmp_log_filepath = get_tmp_filename(prefix='UclustReferenceOtuPicker',
                                            suffix='log')
        tmp_result_filepath = get_tmp_filename(prefix='UclustReferenceOtuPicker',
                                            suffix='txt')
        tmp_failure_filepath = get_tmp_filename(prefix='UclustReferenceOtuPicker',
                                            suffix='txt')
        seqs = [('s1','ACCTTGTTACTTT'),
                ('s2','ACCTAGTTACTTT'),
                ('s3','TTGCGTAACGTTTGAC'),
                ('s4','GTCAAACGTTACGCAA')]
        ref_seqs = [
                ('r1','ACCTCGTTACTTT')]
        
        # rev strand matching disabled
        uc = UclustReferenceOtuPicker({'Similarity':0.8,
                                       'suppress_new_clusters':True,
                                       'save_uc_files':False})
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
         "Reference seqs:%s" % ref_seqs_fp,
         "Similarity:0.8",
         "Application:uclust",
         "enable_rev_strand_matching:True",
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

        self.assertEqualItems(log_file_str.split('\n'), log_file_99_exp)
        
        failures_file_99_exp = ["s3","s4"]
        self.assertEqualItems(fail_file_str.split('\n'), failures_file_99_exp)
        
        
    def test_default_parameters_new_clusters_allowed(self):
        """UclustReferenceOtuPicker: default parameters, new clusters allowed
        """
        uc = UclustReferenceOtuPicker({'save_uc_files':False})
        obs = uc(self.tmp_seq_filepath1,self.temp_ref_filepath1)
        exp = {'ref1':['uclust_test_seqs_0'],
               'ref2':['uclust_test_seqs_1'],
               'ref3':['uclust_test_seqs_2'],
               'ref4':['uclust_test_seqs_3'],
               'QiimeOTU1':['uclust_test_seqs_4'],
               'QiimeOTU2':['uclust_test_seqs_5'],
               'QiimeOTU3':['uclust_test_seqs_6'],
               'QiimeOTU4':['uclust_test_seqs_7'],
               'QiimeOTU5':['uclust_test_seqs_8'],
               'QiimeOTU6':['uclust_test_seqs_9']}
        
        # expected number of clusters observed
        self.assertEqual(len(obs),len(exp))
        
        expected_ref_hits = ['ref1','ref2','ref3','ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k],exp[k])
        
        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of 
        # iteration over a dict
        exp_cluster_ids = exp.keys()
        exp_cluster_ids.sort()
        exp_clusters = exp.values()
        exp_clusters.sort()
        
        obs_cluster_ids = obs.keys()
        obs_cluster_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
        
        self.assertEqual(obs_cluster_ids,exp_cluster_ids)
        self.assertEqual(obs_clusters,exp_clusters)
        
    def test_default_parameters_new_clusters_allowed_save_uc_files(self):
        """UclustReferenceOtuPicker: default parameters, saves uc file
        """
        uc = UclustReferenceOtuPicker({'save_uc_files':True,
                                       'output_dir':self.temp_dir})
                                       
        obs = uc(self.tmp_seq_filepath1,self.temp_ref_filepath1)
        exp = {'ref1':['uclust_test_seqs_0'],
               'ref2':['uclust_test_seqs_1'],
               'ref3':['uclust_test_seqs_2'],
               'ref4':['uclust_test_seqs_3'],
               'QiimeOTU1':['uclust_test_seqs_4'],
               'QiimeOTU2':['uclust_test_seqs_5'],
               'QiimeOTU3':['uclust_test_seqs_6'],
               'QiimeOTU4':['uclust_test_seqs_7'],
               'QiimeOTU5':['uclust_test_seqs_8'],
               'QiimeOTU6':['uclust_test_seqs_9']}
        
        # expected number of clusters observed
        self.assertEqual(len(obs),len(exp))
        
        expected_ref_hits = ['ref1','ref2','ref3','ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k],exp[k])
        
        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of 
        # iteration over a dict
        exp_cluster_ids = exp.keys()
        exp_cluster_ids.sort()
        exp_clusters = exp.values()
        exp_clusters.sort()
        
        obs_cluster_ids = obs.keys()
        obs_cluster_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
        
        self.assertEqual(obs_cluster_ids,exp_cluster_ids)
        self.assertEqual(obs_clusters,exp_clusters)
        
        
    
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
        uc = UclustReferenceOtuPicker({'Similarity':0.90,
                                      'suppress_sort':False,
                                      'presort_by_abundance':False,
                                      'save_uc_files':False})
        obs = uc(self.tmp_seq_filepath1,self.temp_ref_filepath1)
        exp = {'ref1':['uclust_test_seqs_0'],
               'ref2':['uclust_test_seqs_1'],
               'ref3':['uclust_test_seqs_2'],
               'ref4':['uclust_test_seqs_3'],
               'QiimeOTU1':['uclust_test_seqs_4'],
               'QiimeOTU2':['uclust_test_seqs_5'],
               'QiimeOTU3':['uclust_test_seqs_6','uclust_test_seqs_8'],
               'QiimeOTU4':['uclust_test_seqs_7'],
               'QiimeOTU5':['uclust_test_seqs_9']}
        
        # expected number of clusters observed
        self.assertEqual(len(obs),len(exp))
        
        expected_ref_hits = ['ref1','ref2','ref3','ref4']
        for k in expected_ref_hits:
            # seqs that hit refs should have same otu_id and cluster
            self.assertEqual(obs[k],exp[k])
        
        # testing is harder for new clusters, since the otu identifiers are
        # arbitrary, and otu identifier assignment is based on order of 
        # iteration over a dict
        exp_cluster_ids = exp.keys()
        exp_cluster_ids.sort()
        exp_clusters = exp.values()
        exp_clusters.sort()
        
        obs_cluster_ids = obs.keys()
        obs_cluster_ids.sort()
        obs_clusters = obs.values()
        obs_clusters.sort()
        
        self.assertEqual(obs_cluster_ids,exp_cluster_ids)
        self.assertEqual(obs_clusters,exp_clusters)
        
    def test_default_parameters_new_clusters_disallowed(self):
        """UclustReferenceOtuPicker: default params, new clusters not allowed
        """
        uc = UclustReferenceOtuPicker({'suppress_new_clusters':True,
                                       'save_uc_files':False})
        obs = uc(self.tmp_seq_filepath1,self.temp_ref_filepath1)
        exp = {'ref1':['uclust_test_seqs_0'],
               'ref2':['uclust_test_seqs_1'],
               'ref3':['uclust_test_seqs_2'],
               'ref4':['uclust_test_seqs_3']}
        
        # expected number of clusters observed
        self.assertEqual(obs,exp)
        
        
    def test_alt_parameters_new_clusters_disallowed(self):
        """UclustReferenceOtuPicker: alt params, new clusters not allowed
        """
        uc = UclustReferenceOtuPicker({'suppress_new_clusters':True,
                                       'Similarity':1.0,
                                       'save_uc_files':False})
        obs = uc(self.tmp_seq_filepath1,self.temp_ref_filepath1)
        exp = {'ref3':['uclust_test_seqs_2'],'ref4':['uclust_test_seqs_3']}
        
        # expected number of clusters observed
        self.assertEqual(obs,exp)
        

class CdHitOtuPickerTests(TestCase):
    """ Tests of the cd-hit-based OTU picker """

    def setUp(self):
        # create the temporary input files
        self.tmp_seq_filepath1 = get_tmp_filename(\
         prefix='CdHitOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath1,'w')
        seq_file.write(dna_seqs_1)
        seq_file.close()        
        
        self.tmp_seq_filepath2 = get_tmp_filename(\
         prefix='CdHitOtuPickerTest_',\
         suffix='.fasta')
        seq_file = open(self.tmp_seq_filepath2,'w')
        seq_file.write(dna_seqs_2)
        seq_file.close()
        
        self._files_to_remove =\
         [self.tmp_seq_filepath1, self.tmp_seq_filepath2]
        
    def tearDown(self):
        remove_files(self._files_to_remove)

    def test_call_default_params(self):
        """CdHitOtuPicker.__call__ returns expected clusters default params"""

        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        
        exp = {0:['cdhit_test_seqs_0'],\
               1:['cdhit_test_seqs_1'],\
               2:['cdhit_test_seqs_2'],\
               3:['cdhit_test_seqs_3'],\
               4:['cdhit_test_seqs_4'],\
               5:['cdhit_test_seqs_5'],\
               6:['cdhit_test_seqs_6'],\
               7:['cdhit_test_seqs_7'],\
               8:['cdhit_test_seqs_8'],\
               9:['cdhit_test_seqs_9']}
        
        app = CdHitOtuPicker(params={})
        obs = app(self.tmp_seq_filepath1)
        self.assertEqual(obs, exp)
        
    def test_call_alt_threshold(self):
        """CdHitOtuPicker.__call__ returns expected clusters with alt threshold
        """
        # adapted from test_app.test_cd_hit.test_cdhit_clusters_from_seqs
        
        exp = {0:['cdhit_test_seqs_0'],\
               1:['cdhit_test_seqs_1'],\
               2:['cdhit_test_seqs_2'],\
               3:['cdhit_test_seqs_3'],\
               4:['cdhit_test_seqs_4'],\
               5:['cdhit_test_seqs_5'],\
               6:['cdhit_test_seqs_6','cdhit_test_seqs_8'],\
               7:['cdhit_test_seqs_7'],\
               8:['cdhit_test_seqs_9']}

        app = CdHitOtuPicker(params={'Similarity':0.90})
        obs = app(self.tmp_seq_filepath1)
        self.assertEqual(obs, exp)
        
    def test_call_output_to_file(self):
        """CdHitOtuPicker.__call__ output to file functions as expected
        """
        
        tmp_result_filepath = get_tmp_filename(\
         prefix='CdHitOtuPickerTest.test_call_output_to_file_',\
         suffix='.txt')
        
        app = CdHitOtuPicker(params={'Similarity':0.90})
        obs = app(self.tmp_seq_filepath1,result_path=tmp_result_filepath)
        
        result_file = open(tmp_result_filepath)
        result_file_str = result_file.read()
        result_file.close()
        # remove the result file before running the test, so in 
        # case it fails the temp file is still cleaned up
        remove(tmp_result_filepath)
        
        # compare data in result file to fake expected file
        self.assertEqual(result_file_str, dna_seqs_result_file_90_exp)
        # confirm that nothing is returned when result_path is specified
        self.assertEqual(obs,None)
        
    def test_call_log_file(self):
        """CdHitOtuPicker.__call__ writes log when expected
        """
        
        tmp_log_filepath = get_tmp_filename(\
         prefix='CdHitOtuPickerTest.test_call_output_to_file_l_',\
         suffix='.txt')
        tmp_result_filepath = get_tmp_filename(\
         prefix='CdHitOtuPickerTest.test_call_output_to_file_r_',\
         suffix='.txt')
        
        app = CdHitOtuPicker(params={'Similarity':0.99})
        obs = app(self.tmp_seq_filepath1,\
         result_path=tmp_result_filepath,log_path=tmp_log_filepath)
        
        log_file = open(tmp_log_filepath)
        log_file_str = log_file.read()
        log_file.close()
        # remove the temp files before running the test, so in 
        # case it fails the temp file is still cleaned up
        remove(tmp_log_filepath)
        remove(tmp_result_filepath)
        
        log_file_99_exp = ["CdHitOtuPicker parameters:",\
         "Similarity:0.99","Application:cdhit",\
         'Algorithm:cdhit: "longest-sequence-first list removal algorithm"',\
         'No prefix-based prefiltering.',\
         "Result path: %s" % tmp_result_filepath]
        # compare data in log file to fake expected log file
        # NOTE: Since app.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that 
        # the equal unordered lists of lines is present in actual and expected
        self.assertEqualItems(log_file_str.split('\n'), log_file_99_exp)
        
    def test_prefilter_exact_prefixes_no_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when no seqs get filtered
        """
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTACAA'),\
                ('s3','ACGTAG'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCAA'),\
                ('s6','ACGTCCAAAAAAAAAAAA')]
        
        prefix_length = 6
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = seqs, {'s1':['s1'], 's2':['s2'], \
                            's3':['s3'], 's4':['s4'], \
                            's5':['s5'], 's6':['s6']}
        self.assertEqual(actual,expected)
        
        # same result if prefix_length is too long
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTACAA'),\
                ('s3','ACGTAG'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCAA'),\
                ('s6','ACGTCCAAAAAAAAAAAA')]
        prefix_length = 42
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = seqs, {'s1':['s1'], 's2':['s2'], \
                            's3':['s3'], 's4':['s4'], \
                            's5':['s5'], 's6':['s6']}
        self.assertEqual(actual,expected)
        
    def test_prefilter_exact_prefixes_all_to_one_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when all seqs map to one
        """
        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1 comment','ACGTAA'),\
                ('s2','ACGTAC'),\
                ('s3','ACGTAG'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCA'),\
                ('s6','ACGTCC')]
        
        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = [('s1','ACGTAA')], {'s1':['s1','s2','s3','s4','s5','s6']}
        self.assertEqual(actual,expected)
        
        # maps to longest seq
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTACA'),\
                ('s3','ACGTAG'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCA'),\
                ('s6','ACGTCC')]
        
        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = [('s2','ACGTACA')], {'s2':['s1','s2','s3','s4','s5','s6']}
        self.assertEqual(actual,expected)
        
        # maps to longest seq
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTACA'),\
                ('s3','ACGTAGAA'),\
                ('s4','ACGTATAAA'),\
                ('s5','ACGTCAAAAA'),\
                ('s6','ACGTCCAAAAA')]
        
        prefix_length = 4
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = [('s6','ACGTCCAAAAA')], {'s6':['s1','s2','s3','s4','s5','s6']}
        self.assertEqual(actual,expected)
        
    def test_prefilter_exact_prefixes_filtering(self):
        """ CdHitOtuPicker._prefilter_exact_prefixes fns as expected when filtering occurs
        """
        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTAC'),\
                ('s3','ACGTAG'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCA'),\
                ('s6','ACGTCC')]
        
        prefix_length = 5
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = [('s1','ACGTAA'),('s5','ACGTCA')], \
                   {'s1':['s1','s2','s3','s4'],'s5':['s5','s6']}
        self.assertEqual(actual,expected) 
        
        # maps to first when all are same length
        app = CdHitOtuPicker(params={})
        seqs = [('s1','ACGTAA'),\
                ('s2','ACGTAC'),\
                ('s3','ACGTAGAAAA'),\
                ('s4','ACGTAT'),\
                ('s5','ACGTCA'),\
                ('s6','ACGTCC')]
        
        prefix_length = 5
        actual = app._prefilter_exact_prefixes(seqs,prefix_length)
        actual[0].sort()
        expected = [('s3','ACGTAGAAAA'),('s5','ACGTCA')], \
                   {'s3':['s1','s2','s3','s4'],'s5':['s5','s6']}
        self.assertEqual(actual,expected)  
        
    def test_map_filtered_clusters_to_full_clusters(self):
        """CdHitOtuPicker._map_filtered_clusters_to_full_clusters functions as expected
        """
        # original and mapped full clusters are the same
        app = CdHitOtuPicker(params={})
        filter_map = {'s1':['s1'], 's2':['s2'], \
                      's3':['s3'], 's4':['s4'], \
                      's5':['s5'], 's6':['s6']}
        clusters = [['s1'], ['s2'], ['s3'], ['s4'], ['s5'], ['s6']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        expected = clusters
        self.assertEqual(actual,expected)
        
        # original and mapped full clusters are not the same
        filter_map = {'s1':['s1','s2','s3','s4'],'s5':['s5','s6']}
        clusters = [['s1','s5']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        for e in actual: e.sort()
        expected = [['s1','s2','s3','s4','s5','s6']]
        self.assertEqual(actual,expected)
        
        filter_map = {'s1':['s1','s2','s6'],'s3':['s3'],'s5':['s4','s5']}
        clusters = [['s1','s3'],['s5']]
        actual = app._map_filtered_clusters_to_full_clusters(clusters,filter_map)
        for e in actual: e.sort()
        expected = [['s1','s2','s3','s6'],['s4','s5']]
        self.assertEqual(actual,expected)
        
    def test_call_prefilters_when_requested(self):
        """ CdHitOtuPicker.__call__ prefilters when requested
        """
        # no pre-filtering results in one cluster per sequence as they all
        # differ at their 3' ends
        app = CdHitOtuPicker(params={})
        app = CdHitOtuPicker(params={'Similarity':0.99})
        self.assertEqual(app(self.tmp_seq_filepath2),dna_seqs_2_result)
        
        # no pre-filtering results in one cluster per sequence as they are all
        # the same at their 5' ends
        app = CdHitOtuPicker(params={})
        app = CdHitOtuPicker(params={'Similarity':0.99},)
        self.assertEqual(app(self.tmp_seq_filepath2,prefix_prefilter_length=5),\
                             dna_seqs_2_result_prefilter)


class PickOtusStandaloneFunctions(TestCase):
    """ Tests of stand-alone functions in pick_otus.py """
    
    def setUp(self):
        """
        """
        self.otu_map1 = {'0':['seq1','seq2','seq5'],\
                         '1':['seq3','seq4'],\
                         '2':['seq6','seq7','seq8']}
        self.otu_map2 = {'110':['0','2'],\
                         '221':['1']}
        self.otu_map3 = {'a':['110','221']}
        
        
        self.otu_map1_file = ['0\tseq1\tseq2\tseq5',\
                              '1\tseq3\tseq4',\
                              '2\tseq6\tseq7\tseq8']
        self.otu_map2_file = ['110\t0\t2',\
                              '221\t1']
        self.otu_map3_file = ['a\t110\t221']
        
        self.failures1 = ['110']
        self.failures2 = ['110\n','221']
        self.failures3 = ['a']
        
        
    def test_expand_failures_one_otu_map(self):
        """expanding failures generated by chained otu picking fns as expected
        """
        expected_f1 = ['0','2']
        self.assertEqualItems(expand_failures(self.failures1,self.otu_map2),
                         expected_f1)
        expected_f2 = ['0','1','2']
        self.assertEqualItems(expand_failures(self.failures2,self.otu_map2),
                         expected_f2)

    def test_expand_failures_two_otu_maps(self):
        """expanding failures generated by chained otu picking fns as expected
        """
        expected_f1 = ['seq1','seq2','seq5','seq6','seq7','seq8']
        
        actual   = expand_failures(self.failures1,\
         expand_otu_map_seq_ids(self.otu_map2,self.otu_map1))
        self.assertEqualItems(actual,expected_f1)
        
    def test_map_otu_map_files_failures_file_two_otu_maps1(self):
        """map_otu_map_files: correctly maps two otu files and failures
        """
        exp = ['seq1','seq2','seq5','seq6','seq7','seq8']
        actual = map_otu_map_files(\
         [self.otu_map1_file,self.otu_map2_file],
         self.failures1)
        self.assertEqual(actual,exp)
        
    def test_map_otu_map_files_failures_file_two_otu_maps2(self):
        """map_otu_map_files: correctly maps two otu files and failures (alt failures)
        """
        exp = ['seq1','seq2','seq5','seq6','seq7','seq8','seq3','seq4']
        actual = map_otu_map_files(\
         [self.otu_map1_file,self.otu_map2_file],
         self.failures2)
        self.assertEqual(actual,exp)
        
    def test_map_otu_map_files_failures_file_three_otu_maps(self):
        """map_otu_map_files: correctly maps three otu files and failures
        """
        exp = ['seq1','seq2','seq5','seq6','seq7','seq8','seq3','seq4']
        actual = map_otu_map_files(\
         [self.otu_map1_file,self.otu_map2_file,self.otu_map3_file],
         self.failures3)
        self.assertEqual(actual,exp)
        
        
    def test_expand_otu_map_seq_ids_error(self):
        """expand_otu_map_seq_ids: error on missing seq_ids
        """
        self.assertRaises(KeyError,expand_otu_map_seq_ids,\
         self.otu_map3,self.otu_map1)
        
    def test_expand_otu_map_seq_ids_two(self):
        """expand_otu_map_seq_ids: correctly maps seq_ids from two otu maps
        """
        exp12 = {'110':['seq1','seq2','seq5','seq6','seq7','seq8'],\
                 '221':['seq3','seq4']}
        actual12 = expand_otu_map_seq_ids(self.otu_map2,self.otu_map1)
        self.assertEqual(actual12,exp12)
                 
    def test_expand_otu_map_seq_ids_three(self):
        """expand_otu_map_seq_ids: correctly maps seq_ids from three otu maps
        """
        exp123 = {'a':['seq1','seq2','seq5','seq6',\
                       'seq7','seq8','seq3','seq4']}
        actual123 = expand_otu_map_seq_ids(self.otu_map3,\
         expand_otu_map_seq_ids(self.otu_map2,self.otu_map1))
        self.assertEqual(actual123,exp123)
        
    def test_map_otu_map_files_two(self):
        """map_otu_map_files: correctly maps two otu files
        """
        exp12 = {'110':['seq1','seq2','seq5','seq6','seq7','seq8'],\
                 '221':['seq3','seq4']}
        actual12 = map_otu_map_files([self.otu_map1_file,self.otu_map2_file])
        self.assertEqual(exp12,actual12)
        
    def test_map_otu_map_files_three(self):
        """map_otu_map_files: correctly maps three otu files
        """
        exp123 = {'a':['seq1','seq2','seq5','seq6',\
                       'seq7','seq8','seq3','seq4']}
        actual123 = map_otu_map_files(\
         [self.otu_map1_file,self.otu_map2_file,self.otu_map3_file])
        self.assertEqual(exp123,actual123)
        
        # third 'file' contains mixed tabs and spaces
        actual123 = map_otu_map_files(\
         [self.otu_map1_file,self.otu_map2_file,['a\t110 221']])
        self.assertEqual(exp123,actual123)
    
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
dna_seqs_2_result = {0: ['cdhit_test_seqs_2'],\
                     1: ['cdhit_test_seqs_0'],\
                     2: ['cdhit_test_seqs_1']}

dna_seqs_2_result_prefilter =\
 {0: ['cdhit_test_seqs_0','cdhit_test_seqs_1','cdhit_test_seqs_2']}


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
dna_seqs_4_result = {0: ['uclust_test_seqs_2'],\
                     1: ['uclust_test_seqs_0'],\
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
 {0: ['uclust_test_seqs_0','uclust_test_seqs_1','uclust_test_seqs_2']}
 
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

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
