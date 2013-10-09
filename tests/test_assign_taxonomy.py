#!/usr/bin/env python

"""Tests of code for assigning taxonomy"""

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Greg Caporaso", "Kyle Bittinger", "David Soergel",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from cStringIO import StringIO
from os import remove, system, path, getenv
from os.path import exists
from glob import glob
from tempfile import NamedTemporaryFile, mkdtemp
from shutil import copy as copy_file, rmtree

from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs
from cogent.app.util import ApplicationError
from cogent.app.formatdb import build_blast_db_from_fasta_path
from cogent.app.rdp_classifier import train_rdp_classifier
from cogent.util.misc import remove_files, create_dir
from cogent.parse.fasta import MinimalFastaParser

from qiime.util import get_tmp_filename, get_qiime_temp_dir
from qiime.test import initiate_timeout, disable_timeout

from qiime.assign_taxonomy import (
    TaxonAssigner, BlastTaxonAssigner, RdpTaxonAssigner, RtaxTaxonAssigner,
    RdpTrainingSet, RdpTree, _QIIME_RDP_TAXON_TAG, validate_rdp_version,
    MothurTaxonAssigner, UclustConsensusTaxonAssigner)
from sys import stderr


class TopLevelFunctionTests(TestCase):
    """ """
    def test_validate_rdp_version(self):
        """validate_rdp_version raises Exception if version is not valid"""
        rdp22_fp = "/Applications/rdp_classifier_2.2/rdp_classifier-2.2.jar"
        self.assertEqual(validate_rdp_version(rdp22_fp), 2.2)

        invalid_fp = "/Applications/rdp_classifier_2.2/rdp_classifier.jar"
        self.assertRaises(RuntimeError, validate_rdp_version, invalid_fp)

        rdp20_fp = "/Applications/rdp_classifier_2.2/rdp_classifier-2.0.jar"
        self.assertRaises(RuntimeError, validate_rdp_version, rdp20_fp)


class TaxonAssignerTests(TestCase):
    """Tests of the abstract TaxonAssigner class"""

    def test_init(self):
        """Abstract TaxonAssigner __init__ should store name, params"""
        p = TaxonAssigner({})
        self.assertEqual(p.Name, 'TaxonAssigner')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract TaxonAssigner __call__ should raise NotImplementedError"""
        p = TaxonAssigner({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs')

class UclustConsensusTaxonAssignerTests(TestCase):
    
    def setUp(self):
        
        self.files_to_remove = []
        self.dirs_to_remove = []
        
        # Create test output directory
        tmp_dir = get_qiime_temp_dir()
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_uclust_tax_tests',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)
        
        self.inseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='in',
                                            suffix='.fasta')
        self.refseqs1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                            prefix='in',
                                            suffix='.fasta')
        self.id_to_tax1_fp = get_tmp_filename(tmp_dir=self.test_out,
                                              prefix='id-to-tax',
                                              suffix='.txt')
        
        self.output_log_fp = get_tmp_filename(tmp_dir=self.test_out,
                                              prefix='out',
                                              suffix='.log')
        self.output_uc_fp = get_tmp_filename(tmp_dir=self.test_out,
                                              prefix='out',
                                              suffix='.uc')
        self.output_txt_fp = get_tmp_filename(tmp_dir=self.test_out,
                                              prefix='out',
                                              suffix='.txt')
        
        inseqs1_f = open(self.inseqs1_fp,'w')
        inseqs1_f.write(uclust_inseqs1)
        inseqs1_f.close()
        self.files_to_remove.append(self.inseqs1_fp)

        refseqs1_f = open(self.refseqs1_fp,'w')
        refseqs1_f.write(uclust_refseqs1)
        refseqs1_f.close()
        self.files_to_remove.append(self.refseqs1_fp)

        id_to_tax1_f = open(self.id_to_tax1_fp,'w')
        id_to_tax1_f.write(uclust_id_to_tax1)
        id_to_tax1_f.close()
        self.files_to_remove.append(self.id_to_tax1_fp)
        
        self.uc1_lines = uc1.split('\n')
        
        # Define number of seconds a test can run for before timing out 
        # and failing
        initiate_timeout(60)

    
    def tearDown(self):
        
        disable_timeout()
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)
    
    def test_uclust_assigner_write_to_file(self):
        """UclustConsensusTaxonAssigner returns without error, writing results
        """
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        
        t = UclustConsensusTaxonAssigner(params)
        result = t(seq_path=self.inseqs1_fp,
                   result_path=self.output_txt_fp,
                   uc_path=self.output_uc_fp,
                   log_path=self.output_log_fp)
        del t
        # result files exist after the UclustConsensusTaxonAssigner
        # no longer exists
        self.assertTrue(exists(self.output_txt_fp))
        self.assertTrue(exists(self.output_uc_fp))
        self.assertTrue(exists(self.output_log_fp))
        
        # check that result has the expected lines
        output_lines = list(open(self.output_txt_fp,'U'))
        self.assertTrue('q1\tA;F;G\t1.00\t1\n' in output_lines)
        self.assertTrue('q2\tA;H;I;J\t1.00\t1\n' in output_lines)

    def test_uclust_assigner(self):
        """UclustConsensusTaxonAssigner returns without error, returning dict
        """
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        
        t = UclustConsensusTaxonAssigner(params)
        result = t(seq_path=self.inseqs1_fp,
                   result_path=None,
                   uc_path=self.output_uc_fp,
                   log_path=self.output_log_fp)
                   
        self.assertEqual(result['q1'],(['A','F','G'],1.0,1))
        self.assertEqual(result['q2'],(['A','H','I','J'],1.0,1))

        # no result paths provided
        t = UclustConsensusTaxonAssigner(params)
        result = t(seq_path=self.inseqs1_fp,
                   result_path=None,
                   uc_path=None,
                   log_path=None)
                   
        self.assertEqual(result['q1'],(['A','F','G'],1.0,1))
        self.assertEqual(result['q2'],(['A','H','I','J'],1.0,1))
    
    def test_get_consensus_assignment(self):
        """_get_consensus_assignment fuctions as expected """
        in1 = [['Ab','Bc','De'],
               ['Ab','Bc','Fg','Hi'],
               ['Ab','Bc','Fg','Jk']]
        
        # defaults 
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        expected = (['Ab','Bc','Fg'],2./3.,3)
        t = UclustConsensusTaxonAssigner(params)
        self.assertEqual(t._get_consensus_assignment(in1),
                         expected)
        
        # increased confidence yields decreased specificity
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp,
                  'confidence':0.99}
        expected = (['Ab','Bc'],1.0,3)
        t = UclustConsensusTaxonAssigner(params)
        self.assertEqual(t._get_consensus_assignment(in1),
                         expected)
        
        ## if only a single input assignment, it is returned as consensus
        in2 = [['Ab','Bc','De']]
        
        # increased confidence
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp,
                  'confidence':1.0}
        expected = (['Ab','Bc','De'],1.0,1)
        t = UclustConsensusTaxonAssigner(params)
        self.assertEqual(t._get_consensus_assignment(in2),
                         expected)
                         
        # decreased confidence
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp,
                  'confidence':0.0}
        expected = (['Ab','Bc','De'],1.0,1)
        t = UclustConsensusTaxonAssigner(params)
        self.assertEqual(t._get_consensus_assignment(in2),
                         expected)
    
    def test_get_consensus_assignment_adjusts_resolution(self):
        """_get_consensus_assignment max result depth is that of shallowest assignment
        """
        in1 = [['Ab','Bc','Fg'],
               ['Ab','Bc','Fg','Hi'],
               ['Ab','Bc','Fg','Hi']]
        
        # defaults 
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        expected = (['Ab','Bc','Fg'],1.0,3)
        t = UclustConsensusTaxonAssigner(params)
        self.assertEqual(t._get_consensus_assignment(in1),
                         expected)
    
    def test_uc_to_assignments(self):
        """_uc_to_assignments functions as expected"""
        expected = {'q1':[['A','B','C','D'],
                          ['A','B','C','E']],
                    'q2':[['A','H','I','J'],
                          ['A','H','K','L','M'],
                          ['A','H','I','J']],
                    'q3':[[None]],
                    'q4':[[None]],
                    'q5':[[None]]
                    }
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        t = UclustConsensusTaxonAssigner(params)
        actual = t._uc_to_assignments(self.uc1_lines)
        self.assertEqual(actual,expected)

    def test_uc_to_assignment(self):
        """_uc_to_assignment functions as expected"""
        expected = {'q1':(['A','B','C'],1.0,2),
                    'q2':(['A','H','I','J'],2./3.,3),
                    'q3':(['Unassigned'],1.0,0),
                    'q4':(['Unassigned'],1.0,0),
                    'q5':(['Unassigned'],1.0,0)
                    }
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp}
        t = UclustConsensusTaxonAssigner(params)
        actual = t._uc_to_assignment(self.uc1_lines)
        self.assertEqual(actual,expected)
        
        # change label for unassignable
        expected = {'q1':(['A','B','C'],1.0,2),
                    'q2':(['A','H','I','J'],2./3.,3),
                    'q3':(['x'],1.0,0),
                    'q4':(['x'],1.0,0),
                    'q5':(['x'],1.0,0)
                    }
        params = {'id_to_taxonomy_fp':self.id_to_tax1_fp,
                  'reference_sequences_fp':self.refseqs1_fp,
                  'unassignable_label':'x'}
        t = UclustConsensusTaxonAssigner(params)
        actual = t._uc_to_assignment(self.uc1_lines)
        self.assertEqual(actual,expected)
        

uc1 = """# uclust --input /Users/caporaso/Dropbox/code/short-read-tax-assignment/data/qiime-mock-community/Broad-1/rep_set.fna --lib /Users/caporaso/data/gg_13_5_otus/rep_set/97_otus.fasta --uc /Users/caporaso/outbox/uclust_tax_parameter_sweep/Broad-1/gg_13_5_otus/uclust/id1.000000_ma3.uc --id 1.00 --maxaccepts 3 --libonly --allhits
# version=1.2.22
# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=QueryLabel, 10=TargetLabel
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NoHit
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.
N	*	195	*	*	*	*	*	q3	*
N	*	191	*	*	*	*	*	q4	*
N	*	192	*	*	*	*	*	q5	*
L	748	1374	*	*	*	*	*	1081058	*
H	r3	193	100.0	+	0	0	534I193M787I	q2	r3
H	r5	193	97.0	+	0	0	534I193M787I	q2	r5
H	r6	193	97.0	+	0	0	534I193M787I	q2	r6
L	92734	1541	*	*	*	*	*	4440404	*
H	r2	189	99.0	+	0	0	531I189M821I	q1	r2
H	r4	189	100.0	+	0	0	531I189M821I	q1	r4
"""

uclust_id_to_tax1 = """r1	A;F;G
r2	A;B;C;D
r3	A;H;I;J
r4	A;B;C;E
r5	A;H;K;L;M
r6	A;H;I;J
"""

uclust_inseqs1 = """>q1 equal to r1
GAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGCAGCATGACTTAGCTTGCTAAGTTGATGGCGAGTGGCGAACGGGTGAGTAACGCGTAGGAATATGCCTTAAAGAGGGGGACAACTTGGGGAAACTCAAGCTAATACCGCATAAACTCTTCGGAGAAAAGCTGGGGACTTTCGAGCCTGGCGCTTTAAGATTAGCCTGCGTCCGATTAGCTAGTTGGTAGGGTAAAGGCCTACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGGGCAACCCTGATCCAGCAATGCCGCGTGTGTGAAGAAGGCCTGAGGGTTGTAAAGCACTTTCAGTGGGGAGGAGGGTTTCCCGGTTAAGAGCTAGGGGCATTGGACGTTACCCACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCCGCGGTAATACGGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCCGTTAAAANGGTGCCTAAGGTGGTTTGGATNAGTTATGTGTTAAATTCCCTGGCGCCTCCACCCTGGNGCCAGGTCCATANTAAAAACTGTTAAACTCCGAAGTATGGGCACAAGGTAANTTGGAAANTTCCGGTGGTNANCCGNTGAAAATGCGCTTAGAGATNCGGGAAGGGACCACCCCAGTGGGGAAGGCGGCTACCTGGCCTAATAACTGACATTGAGGCACGAAAAGCGTGGGGAGCAACCAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTCAACTAGCTGTNGGTTATATGAATATAATTAGTGGCGAAGCTAACGCGATAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATNGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTACCCTTGACATACAGTAAATCTTTCAGAGATGAGAGAGTGCCTTCGGGAATACTGATACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTATCTCTAGTTGCCAGCGAGTAATGTCGGGAACTCTAAAGAGACTGCCGGTGACAAACCGGAGGAAGGCGGGGACGACGTCAAGTCATCATGGCCCTTACGGGTAGGGCTACACACGTGCTACAATGGCCGATACAGAGGGGCGCGAAGGAGCGATCTGGAGCAAATCTTATAAAGTCGGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCGAATCAGCATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGCTGCACCAGAAGTAGATAGTCTAACCGCAAGGGGGACGTTTACCACGGTGTGGTTCATGACTGGGGTGAAGTCGTAACAAGGTAGCCG
>q2 equal to r3
AGAGTTTGATCATGGCTCAGGATGAACGTTGGTAGTATGCCTAACACATGCAAGTCGAGCGGAAAGTAGTAGCAATATTACCTTTAGCGGCGAACGGGTGAGTAATACTTATCTTACCTGCCATTTAGTGGGGGATAAAATTCCGCATATTATGAGTAATCATGAAAGGGCTTTCGGGCTCGCTAATTGATGGGGATAAGTCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATGCGTAGCCAGTCTGAGAGGATGAATGGCCACAAAGGAACTGAGACACGGTCCTTACTCCTATTGAGGAGGCAGCAGTGGGGAATATTCTGCAATGGGGGAAACCCTGACAGAGCAATACTACGTGAAGGAGGAAGGTCTACGGATTGTAAACTTCTTTACTTAAATGTACTAACCGCAAGGTTTTGCATTTGAGCAAAAAGCGACGACTAACTATGTGCCAGCAGTCGCGGTAAGACATAGGTCGCGAACGTTATCCGGAATTATTGGGCGTAAAGGATGCGTAGATGGTTCAGTAAGTTACTGGTGGGAAATCGAGGCCTAACCTCGTGGAAGTCAGTAATACTGTTGAACTTGAGTGCAGGAGAGGTTAACGGAACTTCATGTGGAGCGGTAAAATGCGTAGATATATGAAAGAACATCAATATAGCGAAGGCAGTTAACTATTCTGCTACTGACATTGAGGCATGAAAGCGTGGGGAGCAAAACGGATTAGATACCCGTGTAGTCCACGCCCTAAACGATGAGTGCTAGATATTGGGAACTTGATTCTCAGTGTCGCTAGGTAACCCGGTAAGCACTCCGCCTGGGGAGTACGCTCGCAAGAGTGAAACTTAAAGGAATTGACGGGGACTCGCACAAGCAGTGGAGCATGTGGTTTAATCCGATACAACGCGTAGAACCTTACCAAGGCTTGACATGTAAGGTAGCAATACTAAATTAGGGAAACCTAGTTGAATTACACAGGTGCTGCATGGCCGTCGTCAGTTCGTGCCGTGAGGTGTATGGTTAAGTCCTATAACGAACGCAACCCTTGTCCTTAGTTGCCAGCATTCAGTTGGGGACTCTAAGGAGACTGCAAGGCATGTTTATGCATGTCTGAACTGAGTAATCAGGAGGAAGGAAAGGATGACGCCTGGTCGTCATGGCCCTTATGCCTTGGGCGACACACGTGCTACAATGGTTAGTACAAAGGGTCGCCAACCCGCGAGGGGGAGCTAATCTCAAAAAACTAACCGCAGTTCAGATTGGAGGCTGCAACTCGCCTCCATGAAGGTGGAATTGCTAGTAATCGTAAATCAGCCATGTTACGGTGAATACGTTCTCGAGTCTTGTACACACCGCCCGTCAAATCACGAAAGTCGGTAATGCCTAAAACCAGTTGATTAACCCGCAAGGGATACAGCTGACAAGGGTAGGATTGGCAATTGGGGTTAAGTCGTAACAAGGTAGCCGTAG
"""

uclust_refseqs1 = """>r1 229854
GAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCATGCTTAACACATGCAAGTCGAACGGCAGCATGACTTAGCTTGCTAAGTTGATGGCGAGTGGCGAACGGGTGAGTAACGCGTAGGAATATGCCTTAAAGAGGGGGACAACTTGGGGAAACTCAAGCTAATACCGCATAAACTCTTCGGAGAAAAGCTGGGGACTTTCGAGCCTGGCGCTTTAAGATTAGCCTGCGTCCGATTAGCTAGTTGGTAGGGTAAAGGCCTACCAAGGCGACGATCAGTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGGGCAACCCTGATCCAGCAATGCCGCGTGTGTGAAGAAGGCCTGAGGGTTGTAAAGCACTTTCAGTGGGGAGGAGGGTTTCCCGGTTAAGAGCTAGGGGCATTGGACGTTACCCACAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCCGCGGTAATACGGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCCGTTAAAANGGTGCCTAAGGTGGTTTGGATNAGTTATGTGTTAAATTCCCTGGCGCCTCCACCCTGGNGCCAGGTCCATANTAAAAACTGTTAAACTCCGAAGTATGGGCACAAGGTAANTTGGAAANTTCCGGTGGTNANCCGNTGAAAATGCGCTTAGAGATNCGGGAAGGGACCACCCCAGTGGGGAAGGCGGCTACCTGGCCTAATAACTGACATTGAGGCACGAAAAGCGTGGGGAGCAACCAGGATTAGATACCCTGGTAGTCCACGCTGTAAACGATGTCAACTAGCTGTNGGTTATATGAATATAATTAGTGGCGAAGCTAACGCGATAAGTTGACCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAAGGAATNGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTACCCTTGACATACAGTAAATCTTTCAGAGATGAGAGAGTGCCTTCGGGAATACTGATACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGTAACGAGCGCAACCCTTATCTCTAGTTGCCAGCGAGTAATGTCGGGAACTCTAAAGAGACTGCCGGTGACAAACCGGAGGAAGGCGGGGACGACGTCAAGTCATCATGGCCCTTACGGGTAGGGCTACACACGTGCTACAATGGCCGATACAGAGGGGCGCGAAGGAGCGATCTGGAGCAAATCTTATAAAGTCGGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGCGAATCAGCATGTCGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGCTGCACCAGAAGTAGATAGTCTAACCGCAAGGGGGACGTTTACCACGGTGTGGTTCATGACTGGGGTGAAGTCGTAACAAGGTAGCCG
>r2 107103
TTTTCTTGGATTTGATTCTGGTCCAGAGTAAACGCTTGAGATATGTTGATACATGTTAGTTAAACGTGAATATTTGGTTTTTATGCCAACTTTATTTAAGTAGCGTATAGGTGAGTAATATGCAAGAATCCTACCTTTTAGTTTATGTAGCTCGTAAATTTATAAAAGATTTTTTCGCTAAAAGATGGGCTTGCACAAGATTAGGTTTTTGGTTTGCTAAAAACGTTCCAAGCCTAAGATCTTTAGCCGGCTTTCGTGAGTGACCGGCCACATAGGGACTGAGACAATGCCCTAGCTCCTTTTCTGGAGGCATCAGTACAAAGCATTGGACAATGAACGAAAGTTTGATCCAGTAATATCTCGTGAATGATGAAGGGTTTTTGCTCGTAAATTTCTTTTAGTTGAAAGAAAAAAGATATATTTCAACAGAAAAAATCCTGGCAAATCCTCGTGCCAGCAGCCGCGGTAATACGAGAAGGGTTAGCGTTACTCGAAATTATTGGGCGTAAAGTGCGTGAACAGCTGCTTTTTAAGCTATAGGCAGAAAAATCAAGGGTTAATCTTGTTTTTGTCATAGTTCTGATAAGCTTGAGTTTGGAAGAAGATAATAGAACATTTTATGGAGCGATGAAATGCTATGATATAAAAGAGAATACCAAAAGCGAAGGCAGTTATCTAGTACAAAACTGACGCCTATACGCGAAGGCTTAGGTAGCAAAAAGGATTAGGGACCCTTGTAGTCTAAGCTGTCAACGATGAACACTCGTTTTTGGATCACTTTTTTTCAGAAACTAAGCTAACGCGTTAAGTGTTTCGCCTGGGTACTACGGTCGCAAGACTAAAACTTAAAGAAATTGGCGGGAGTAAAAACAAGCAGTGGAGCGTGTGGTTTAATTCGATAGTACACGCAAATCTTACCATTACTTGACTCAAACATTGAAATGCACTATGTTTATGGTGTTGTTTAAGTATTATTTTACTTATAGATGTGCAGGCGCTGCATGGTTGTCGTCAGTTCGTGTCGTGAGATGTTTGGTTAATTCCCTTAACGAACGTAACCCTCAAAGCATATTCAAAACATTTTGTTTTTTTGTTAAACAGTCGGGGAAACCTGAATGTAGAGGGGTAGACGTCTAAATCTTTATGGCCCTTATGTATTTGGGCTACTCATGCGCTACAATGGGTGTATTCTACAAAAAGACGCAAAAACTCTTCAGTTTGAGCAAAACTTGAAAAGCACCCTCTAGTTCGGATTGAACTCTGGAACTCGAGTTCATAAAGTTGGAATTGCTAGTAATCGTGAGTTAGCGTATCGCGGTGAATCGAAAATTTACTTTGTACATACCGCCCGTCAAGTACTGAAAATTTGTATTGCAAGAAATTTTTGGAGAATTTACTTAACTCTTTTTTTTTTTAAGTTGGCTGTATCAGTCTTTTAAAAACTTTGAGTTAGGTTTTAAGCATCCGAGGGTAAAAGCAACATTTTTTATTGGTATTAAGTCGTAACAAGGTAGCCCTACGGG
>r3 696036
AGAGTTTGATCATGGCTCAGGATGAACGTTGGTAGTATGCCTAACACATGCAAGTCGAGCGGAAAGTAGTAGCAATATTACCTTTAGCGGCGAACGGGTGAGTAATACTTATCTTACCTGCCATTTAGTGGGGGATAAAATTCCGCATATTATGAGTAATCATGAAAGGGCTTTCGGGCTCGCTAATTGATGGGGATAAGTCGTATTAGTTAGTTGGTGGGGTAATGGCCTACCAAGACAATGATGCGTAGCCAGTCTGAGAGGATGAATGGCCACAAAGGAACTGAGACACGGTCCTTACTCCTATTGAGGAGGCAGCAGTGGGGAATATTCTGCAATGGGGGAAACCCTGACAGAGCAATACTACGTGAAGGAGGAAGGTCTACGGATTGTAAACTTCTTTACTTAAATGTACTAACCGCAAGGTTTTGCATTTGAGCAAAAAGCGACGACTAACTATGTGCCAGCAGTCGCGGTAAGACATAGGTCGCGAACGTTATCCGGAATTATTGGGCGTAAAGGATGCGTAGATGGTTCAGTAAGTTACTGGTGGGAAATCGAGGCCTAACCTCGTGGAAGTCAGTAATACTGTTGAACTTGAGTGCAGGAGAGGTTAACGGAACTTCATGTGGAGCGGTAAAATGCGTAGATATATGAAAGAACATCAATATAGCGAAGGCAGTTAACTATTCTGCTACTGACATTGAGGCATGAAAGCGTGGGGAGCAAAACGGATTAGATACCCGTGTAGTCCACGCCCTAAACGATGAGTGCTAGATATTGGGAACTTGATTCTCAGTGTCGCTAGGTAACCCGGTAAGCACTCCGCCTGGGGAGTACGCTCGCAAGAGTGAAACTTAAAGGAATTGACGGGGACTCGCACAAGCAGTGGAGCATGTGGTTTAATCCGATACAACGCGTAGAACCTTACCAAGGCTTGACATGTAAGGTAGCAATACTAAATTAGGGAAACCTAGTTGAATTACACAGGTGCTGCATGGCCGTCGTCAGTTCGTGCCGTGAGGTGTATGGTTAAGTCCTATAACGAACGCAACCCTTGTCCTTAGTTGCCAGCATTCAGTTGGGGACTCTAAGGAGACTGCAAGGCATGTTTATGCATGTCTGAACTGAGTAATCAGGAGGAAGGAAAGGATGACGCCTGGTCGTCATGGCCCTTATGCCTTGGGCGACACACGTGCTACAATGGTTAGTACAAAGGGTCGCCAACCCGCGAGGGGGAGCTAATCTCAAAAAACTAACCGCAGTTCAGATTGGAGGCTGCAACTCGCCTCCATGAAGGTGGAATTGCTAGTAATCGTAAATCAGCCATGTTACGGTGAATACGTTCTCGAGTCTTGTACACACCGCCCGTCAAATCACGAAAGTCGGTAATGCCTAAAACCAGTTGATTAACCCGCAAGGGATACAGCTGACAAGGGTAGGATTGGCAATTGGGGTTAAGTCGTAACAAGGTAGCCGTAG
>r4 3190878
TTAATACATGCGAGTTGAACGTGAATTTTTTAATTAAAATGAAAGTAGCGTACTGGTGAGTAACACGTGAGAATCTACCTTTCAAATCAACATAAAATGTTGAATAAAAGCTTCTAAAGCTATAAAGATATGTTTTCGTTGAAAGATGAGCTTGCGCAAGATTAGGTAGTTGGTAAGGTAACGGCTTACCAAGCCAAAGATCTTTAGCTGGTTTGAGAAAATGATCAGCCACATTGGAACTGAAACACAGTCCAAACGTAATATAACGGCAGCAGTAGGGAATTTTGAACACTGAGCGAAAGCTTGATTCAGCCAAGTATCGTGGATGAAGAAGGCTGTCTTTTGGTCGTAAAATCCATTTATATAGTCACATGAAATGTGTCTTTTATTTCGATAAAAGGAAAGATTATGACTTTCTATTGAAAAGTCCCGGCTAATCTCGTGCCAGCAGCCGCGGTAATACGAGAGGGGCAAACGATGTTTAGCATGATTGGGCGTAAAGAGCTTGTAGATGGTTTCTTTTAATTTTATATAAAAGCTCTAAGCTTAACTTTGATTATATATAAAGGAAAGATAACTTGAGTTATGGAAAGGAAAGTAGAATTCTTGGAGGAGAGGTAGAATTTGGTGATATCAAGAGGAATTCCAAAAGCGAAGGCAGCTTTCTTGCCATATACTGACATTGAAGGGCGAAAGCGTGGGTAGCGACAGGGATTAGATACCCCATTAGTCCACGCCGTCAACGATGACCTTTATTTATTGGTTTCTCTTAAAATAAATAAATTATTTTTTAGTTTGATCAGTGAAACAGTTAACGCGTTAAAAGGTCCGCCTGAGGAGTACGATCGCAAGATTAAAACTCAAAAGAATAGACGGGAGCGTTCACAAGTGGTGGAGCATGAAGTTTAATGCGATACAACACGCAAAACCTTACCATTTTTTGATATTTTACTTATCAGTTATTTCTCATGAAATAATGTTTTTTACTAAAGTAAAAATTTGTTTGTATAACAGGCGTTGCATGGCTGTCGTAAGTTCGTACTGTGAAGTGTTGGATTAATTTCCTTAACGAACGTAACCCCTTGGTTTTGTTAAAACTAAAATCTACCGCTAGTCATAAACTAGAGGAAGGGAGGGATCACGTCAAGTCCTCATGACCCTTATAAAATGGGCTACGCTTTTCGTGCTACAATGATAAATACAATAAGAAGCAATAACGAAAGTTGGAGCAAATCTATAAAATTTATCTCAGTTCAGATTGTTCTCTGCAATTCGAGAACATGAAGATGGAATCACTAGTAATCGTAGATCAGCATGCTACGGTGAATATGTAATTACGCTCTGTACTCACAGCCCGTCACACAATGGAAGTAAAATGTATCGGAAATTTGTCAAATATTGTTAGATTTTCTTTTTTAAATTTATTGAATAAATTATTTTAATTAATATCTTTCAACTAAATGGGAACTGATGATATGTTTCATGACTGTTGTGAAGTCGTAACAAGGTAGCGCTAGCGGAAGCTGGTGCTGGAT
>r5 3761685
GATGAACGCTCGCGGCGTGCCTAAGGTATGCAAGTCGAACGCCGACCTTCGGGTCGTGCGTGGCGAACGGGTGAGTAACGCGTGGGAACGTACCGCCGAGATGGGGACAACTCCGCGAAAGCGGAGCTAATACCCAATGGTCCCGGGTCGTACCTCACTGCGTTCGGTACGGATACGATGGACATTCGTATTCCCGCGCCGAGCGTAGCGAGGTGCGACTTGGGTAAAGGAGCAATCCGCTCGGTGAGCGGCCCGCGTCCTATCAGCTTGTTGGTGGGGTAATGGCCTACCAAGGCGATGACGGGTAGGGGGTGTGAGAGCATGATCCCCCACAATGGCACTGAGACAAGGGCCATACACCTACGGGTGGCAGCAATCAGGAATCTTGCGCAATGGGCGAAAGCCTGACGCAGCGACGCCGCGTGGAGGATGAAGGTCGAGAGATTGTAAACTCCTAGCCTTGTGCGTAGCACAAGGGGTTCGTATATCGTTAGTCGGTTACGAAGCCCGTATCGTAAGGCGGTTACGGTTACGTGNAAAAACACGTAACCCAATCACGTGACCGAAACGGGCATCGTTACCGTAACCGGCAGACGTTTCACGAATCTCTCGTGTTACGCGCGAGGGGATAAAGTCCCGGCTAACTACGTGCCAGCAGCAGCGGTAAAACGTANGGGGCGAGCGCTATCCGGATTTATTGGGCGTAAAGAGTGNNCAGGCNGCCAGNNNNGTTTTTCGTTAAAGCCCCCCGCTCAACGGGGGAACTGCGGAGAATACTGTTTGGCTGGAGGAAGGTAGAGGCTGTCGGAACTCGGGGAGGAGGGGTGAAATCCTATGATCCCCCAAGGAACACCAAAGGCGAAGGCAGACAGCTGGGCCTTTCCTGACGCTGAGGCACGAAAGCTAGGGGAGCGAAACGGATTAGATACCCGTGTAGTCCTAGCCGTAAACGATCCCCGCTAGATTTTGGCATCTGTCAGAGTCGAAGCTAACGCATTAAGCGGGGCGCCTGGGAAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGGGACCCGCACAAGCGGTGGAGCGCCTGGTTTAATTCGATGATAACCGAAGAACCTTACCAGGGCTTGACATCCCGACAGAAATCCCGACGAAAGTTGGGCCCCCCGCAAGGGCTGTCGGGACAGGTGTTGCATGGCTGTCGTCAGCTCGTGCCGTGAGGTGTATGCTTAAGTGCCGAAACGAGCGCAACCCCTACCCCATGTTAGAAATGTCGTGGGGGACTGCCTCCGTAAGGGGGAGGAAGGTGGGGACGACGTCAAGTCAGTATGGCCCTTATGCCCTGGGCAACACAGACGCTACAATGGCGGCTACAACGGGTTGCAATGTCGTAAGGCGGAGCTAATCCTAAAAAGCCGTCTCAGTTCGGATTGGGGGCTGAAACCCGCCCCCATGAAGTTGGAATCGCTAGTAATGGCCGGTCAGCTATACGGCCGTGAATACGTTCTCGGGTCTTGTACACACCGCCCGTCA
>r6 4251079
AGTGGCGAACGGGTGCGTAACGCGTGGGAATCTGCCGAACAGTTCGGGCCAAATCCTGAATCAAGCTAAAAAGCGCTGTTTGATGAGCCTGCGTAGTATTAGGTAGTTGGTCAGGTAAAGGCTGACCAAGCCAATGATGCTTAGCTGGTCTTTTCGGATGATCAGCCACACTGGGACTGAGACACGGCCCGGACTCCCACGGGGGGCAGCAGTGGGGAATCTTGGACAATGGGCGAAAGCCCGATCCAGCAATATCGCGTGAGTGAAGAAGGGCAATGCCGCTTGTAAAGCTCTTTCGTCGAGTGCGCGATCATGACAGGACTCGAGGAAGAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAAGACGGGGGGGGCAAGTGTTCTTCGGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTTCAAGTCGCCAAAAACTGGCTCCATGCTCTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGGAGGGGTGAAATCCGGAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGGGGTGCGAAAGCATGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCATGCCGTAAACGATGAGTGTTCGCCCTTGGTCTACGCGGATCAGGGGCCCAGCTAACGCGTGAAACACTCCGCCTGGGGAGTACGGTCGCAAGACCGAAACTCAAAGGAATTGACGGGGGCCTGCACAAGCGGTGGAGCATGTGGTTTAATTCGATACAACGCGCAAAACCTTACCAGCCCTTGACATATGAACAAAAAAACCTGTCCTTAACGGGATGGTACTGACTTTCATACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTTGGTCAAGTCCTATAACGAGCGAAACCCTCGTTTTGTGTTGCTGAGACATGCGCCTAAGGAAAAAGTCTTTGCAACCGAAGTGAGCCGAGGAGCCGAGTGACGTGCCAGCGCTACTAATTGAGTGCCAGCACGTAGCTGTGCTGTCAGTAAGAAGGGAGCCGGCGCCTTTCGAAGCACTTTCTAGTGTGAACCGAACGTCCCGCGTTCCGGCTTGTTCGAATCGCGTCACGAGTCTACAAGCCCCACTGATACCTACCTATAGTGACGTCAAAGTACCAGTGACGGTGACTTGGTTGCGTAACGTAATATGGATTCAGTCAGCGAAACTCCCTCCAACTCAATCAATATCAACAACATGTCGTGACGAGCCTGACTACACTCCACCTGACGGTACTACACTATCTATCAGTAGAGCGCCTTGCGCGATACGGCTTTTTGGCCGTATCTTGCAGGTGACGACGACGTCGAGTTGGCGGCGGAGGAAGACTCGGCATTCAGGCGAGCCGCCCGGTGGTGTGGTACGTAGTGGTAATAGTACGCCCCGCCCGTCACACCCTTCAAGCTCCGAAACAAACAAAAAGGTGCGTGCCGCACTCACGAGGGACTTGCCAGTGATATACTGAAGGAAGGTGGGGATGACGTCAAGTCCGCATGGCCCTTATGGGCTGGGCCACACACGTGCTACAATGGCAATTACAATGGGAAGCAAGGCTGTAAGGCGGAGCGAATCCGGAAAGATTGCCTCAGTTCGGATTGTTCTCTGCAACTCGGGAACATGAAGTTGGAATCGCTAGTAATCGCGGATCAGCATGCCGCGGTGAATATGTACCCGGGCCCTGTACACACCGCCCGTCACACCCTGGGAATTGGTTTCGCCCGAAGCATCGGACCAATGATCACCCATTCCTTGTGTGTACCACTAGTGCCACAAAGGCTTTTGGTGGTCTTATTGGCGCATACCACGGTGGGGTCTTCGACTGGGGTGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGTG
"""


class BlastTaxonAssignerTests(TestCase):
    """Tests of the BlastTaxonAssigner class"""

    def setUp(self):
        self.id_to_taxonomy_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.txt')
        self.input_seqs_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self.reference_seqs_fp = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')

        self._paths_to_clean_up =\
         [self.id_to_taxonomy_fp,self.input_seqs_fp,self.reference_seqs_fp]

        open(self.id_to_taxonomy_fp,'w').write(id_to_taxonomy_string)
        open(self.input_seqs_fp,'w').write(test_seq_coll.toFasta())
        self.test_seqs = test_seq_coll.items()
        open(self.reference_seqs_fp,'w').write(test_refseq_coll.toFasta())

        self.expected1 = {
            's1': ('Archaea;Euryarchaeota;Halobacteriales;uncultured', 0.0, "AY800210"),
            's2': ('Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.', 0.0, "EU883771"),
            's3': ('Archaea;Crenarchaeota;uncultured;uncultured', 0.0, "EF503699"),
            's4': ('Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium', 0.0, "DQ260310"),
            's5': ('Archaea;Crenarchaeota;uncultured;uncultured', 0.0, "EF503697"),
            's6': ('No blast hit', None, None),
            }

    def tearDown(self):
        remove_files(set(self._paths_to_clean_up))

    def test_init(self):
        """BlastTaxonAssigner __init__ should store name, params"""
        p = BlastTaxonAssigner({})
        self.assertEqual(p.Name, 'BlastTaxonAssigner')
        # default parameters correctly initialized
        default_params = {'Min percent identity':90.0,\
         'Max E value':1e-30,\
         'Application':'blastn/megablast'}
        self.assertEqual(p.Params, default_params)

    def test_parse_id_to_taxonomy_file(self):
        """Parsing taxonomy files functions as expected
        """
        lines = id_to_taxonomy_string.splitlines()
        p = BlastTaxonAssigner({})
        expected = {\
         "AY800210":"Archaea;Euryarchaeota;Halobacteriales;uncultured",\
         "EU883771":"Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",\
         "EF503699":"Archaea;Crenarchaeota;uncultured;uncultured",\
         "DQ260310":"Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",\
         "EF503697":"Archaea;Crenarchaeota;uncultured;uncultured"}
        self.assertEqual(p._parse_id_to_taxonomy_file(lines),expected)

    def test_map_ids_to_taxonomy(self):
        """Mapping sequence ids to taxonomy functions as expected
        """
        p = BlastTaxonAssigner({})
        id_to_taxonomy_map = {
            "AY800210": "Archaea;Euryarchaeota;Halobacteriales;uncultured",
            "EU883771": "Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",
            "EF503699": "Archaea;Crenarchaeota;uncultured;uncultured",
            "DQ260310": "Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",
            "EF503697": "Archaea;Crenarchaeota;uncultured;uncultured",
            }
        hits = {
            's1': ("AY800210", 1e-99),
            's5': ("EU883771", 'weird confidence value'),
            's3': ("DQ260310", 42.),
            's4': None,
            }
        expected = {
            's1': ("Archaea;Euryarchaeota;Halobacteriales;uncultured", 1e-99, "AY800210"),
            's5': ('Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.',
                   'weird confidence value',"EU883771"),
            's3': ("Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium", 42.,"DQ260310"),
            's4': ('No blast hit', None, None),
            }
        actual = p._map_ids_to_taxonomy(hits,id_to_taxonomy_map)
        self.assertEqual(actual,expected)

    def test_get_first_blast_hit_per_seq(self):
        """Extracting the first blast hit for each seq functions as expected
        """
        p = BlastTaxonAssigner({})
        blast_hits = {'s1':[('blah',0.0)],\
                      's3':[('dsasd',1e-42),('rrr',1e-12),('qqq',0.001)],\
                      's2':[]}
        expected = {'s1':('blah',0.0),\
                      's3':('dsasd',1e-42),\
                      's2':None}
        actual = p._get_first_blast_hit_per_seq(blast_hits)
        self.assertEqual(actual,expected)

    def test_get_blast_hits(self):
        """BlastTaxonAssigner._get_blast_hits functions w existing db

        """
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove

        p = BlastTaxonAssigner({})
        seq_coll_blast_results = p._get_blast_hits(blast_db,self.test_seqs)
        # mapping from identifier in test_seq_coll to the id of the sequence
        # in the refseq collection (a silva derivative)
        expected_matches = {\
         's1':'AY800210',
         's2':'EU883771',\
         's3':'EF503699',\
         's4':'DQ260310',\
         's5':'EF503697'}

        # no results for s6 (which is a randomly-generated sequence)
        s6_blast_results = seq_coll_blast_results['s6']
        self.assertEqual(s6_blast_results,[])

        # expected results for all other query sequences
        for seq_id in expected_matches:
            blast_results = seq_coll_blast_results[seq_id]
            blast_results_d = dict(blast_results)
            # explicitly checks that the result is in the data before
            # pulling it out (this is redundant, but allows for a useful
            # error message if the data wasn't in there b/c e.g. there
            # were no blast results returned)
            self.assertTrue(expected_matches[seq_id] in blast_results_d)
            # now check that the perfect match got a 0.0 e-value as it should
            # on this data
            self.assertEqual(blast_results_d[expected_matches[seq_id]],0.0)

    def test_call_existing_blast_db(self):
        """BlastTaxonAssigner.__call__ functions w existing db
        """
        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove

        p = BlastTaxonAssigner({'blast_db':blast_db,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        actual = p(self.input_seqs_fp)

        self.assertEqual(actual,self.expected1)

    def test_call_alt_input_types(self):
        """BlastTaxonAssigner.__call__ functions w alt input types """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})

        # neither seqs or seq_fp passed results in AssertionError
        self.assertRaises(AssertionError,p)

        # Functions with a list of (seq_id, seq) pairs
        seqs = list(MinimalFastaParser(open(self.input_seqs_fp)))
        actual = p(seqs=seqs)
        self.assertEqual(actual,self.expected1)

        # Functions with input path
        actual = p(self.input_seqs_fp)
        self.assertEqual(actual,self.expected1)

        # same result when passing fp or seqs
        self.assertEqual(p(seqs=seqs),p(self.input_seqs_fp))

    def test_seqs_to_taxonomy(self):
        """BlastTaxonAssigner._seqs_to_taxonomy: functions as expected
        """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})

        # build the id_to_taxonomy_map as this test doesn't execute __call__
        id_to_taxonomy_map = {
            "AY800210": \
             "Archaea;Euryarchaeota;Halobacteriales;uncultured",
            "EU883771": \
             "Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.",
            "EF503699": \
             "Archaea;Crenarchaeota;uncultured;uncultured",
            "DQ260310": \
             "Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium",
            "EF503697": \
             "Archaea;Crenarchaeota;uncultured;uncultured",
            }

        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove

        # read the input file into (seq_id, seq) pairs
        seqs = list(MinimalFastaParser(open(self.input_seqs_fp)))

        actual = p._seqs_to_taxonomy(seqs,blast_db,id_to_taxonomy_map)
        self.assertEqual(actual,self.expected1)

        # passing empty list of seqs functions as expected
        actual = p._seqs_to_taxonomy([],blast_db,id_to_taxonomy_map)
        self.assertEqual(actual,{})


    def test_call_on_the_fly_blast_db(self):
        """BlastTaxonAssigner.__call__ functions w creating blast db
        """
        p = BlastTaxonAssigner({\
         'reference_seqs_filepath':self.reference_seqs_fp,\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp})
        actual = p(self.input_seqs_fp)

        self.assertEqual(actual,self.expected1)

    def test_call_output_to_file(self):
        """BlastTaxonAssigner.__call__ functions w output to file
        """
        result_path = get_tmp_filename(
            prefix='BlastTaxonAssignerTests_', suffix='.fasta')
        self._paths_to_clean_up.append(result_path)

        p = BlastTaxonAssigner({
            'reference_seqs_filepath': self.reference_seqs_fp,
            'id_to_taxonomy_filepath': self.id_to_taxonomy_fp,
            })
        actual = p(self.input_seqs_fp, result_path=result_path)

        expected_lines = set([
            's1\tArchaea;Euryarchaeota;Halobacteriales;uncultured\t0.0\tAY800210\n',
            's2\tArchaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.\t0.0\tEU883771\n',
            's3\tArchaea;Crenarchaeota;uncultured;uncultured\t0.0\tEF503699\n',
            's4\tArchaea;Euryarchaeota;Methanobacteriales;Methanobacterium\t0.0\tDQ260310\n',
            's5\tArchaea;Crenarchaeota;uncultured;uncultured\t0.0\tEF503697\n',
            's6\tNo blast hit\tNone\tNone\n',
            ])
        f = open(result_path)
        observed_lines = set(f.readlines())
        f.close()
        self.assertEqual(observed_lines, expected_lines)

        # Return value is None when result_path is provided (Not sure
        # if this is what we want yet, or if we would want both so
        # results could be logged to file...)
        self.assertEqual(actual, None)

    def test_call_logs_run(self):
        """BlastTaxonAssigner.__call__ logs the run when expected
        """
        log_path = get_tmp_filename(\
         prefix='BlastTaxonAssignerTests_',suffix='.fasta')
        self._paths_to_clean_up.append(log_path)

        # build the blast database and keep track of the files to clean up
        blast_db, files_to_remove = \
         build_blast_db_from_fasta_path(self.reference_seqs_fp)
        self._paths_to_clean_up += files_to_remove

        p = BlastTaxonAssigner({\
         'id_to_taxonomy_filepath':self.id_to_taxonomy_fp,\
         'blast_db':blast_db})
        actual = p(self.input_seqs_fp,log_path=log_path)

        log_file = open(log_path)
        log_file_str = log_file.read()
        log_file.close()

        log_file_exp = [
            "BlastTaxonAssigner parameters:",
            'Min percent identity:90.0',
            'Application:blastn/megablast',
            'Max E value:1e-30',
            'Result path: None, returned as dict.',
            'blast_db:%s' % str(self.reference_seqs_fp)[1:-1],
            'id_to_taxonomy_filepath:%s' % self.id_to_taxonomy_fp,
            'Number of sequences inspected: 6',
            'Number with no blast hits: 1',
            '',
         ]
        # compare data in log file to fake expected log file
        # NOTE: Since p.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that
        # the equal unordered lists of lines is present in actual and expected
        self.assertEqualItems(log_file_str.split('\n'), log_file_exp)



class RtaxTaxonAssignerTests(TestCase):
    """Tests for the RTAX taxonomy assigner."""

    def setUp(self):
      self.id_to_taxonomy_fp = get_tmp_filename(\
       prefix='RtaxTaxonAssignerTests_',suffix='.txt')
      self.input_seqs_fp = get_tmp_filename(\
       prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
      self.reference_seqs_fp = get_tmp_filename(\
       prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
      self.read_1_seqs_fp = get_tmp_filename(\
       prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
      self.read_2_seqs_fp = get_tmp_filename(\
       prefix='RtaxTaxonAssignerTests_',suffix='.fasta')

      self._paths_to_clean_up =\
       [self.id_to_taxonomy_fp,self.input_seqs_fp,self.reference_seqs_fp,self.read_1_seqs_fp,self.read_2_seqs_fp]

      open(self.id_to_taxonomy_fp,'w').write(rtax_reference_taxonomy)
      open(self.input_seqs_fp,'w').write(rtax_test_repset_fasta)
      open(self.reference_seqs_fp,'w').write(rtax_reference_fasta)
      open(self.read_1_seqs_fp,'w').write(rtax_test_read1_fasta)
      open(self.read_2_seqs_fp,'w').write(rtax_test_read2_fasta)


    def tearDown(self):
        remove_files(set(self._paths_to_clean_up),error_on_missing=False)

    def cleanAll(self, path):
        return [path, path + ".pos.db", path + ".pos.dir", path + ".pos.pag", path + ".lines.dir", path + ".lines.db", path + ".lines.pag"]



    def test_init(self):
        """RtaxTaxonAssigner.__init__ should set default attributes and params
        """
        p = RtaxTaxonAssigner({})
        self.assertEqual(p.Name, 'RtaxTaxonAssigner')
        # default parameters correctly initialized
        default_params = {
            'id_to_taxonomy_fp': None,
            'reference_sequences_fp': None,
            'header_id_regex' : "\\S+\\s+(\\S+?)\/",
            'read_id_regex' : "\\S+\\s+(\\S+)",
            'amplicon_id_regex' : "(\\S+)\\s+(\\S+?)\/",
            'read_1_seqs_fp' : None,
            'read_2_seqs_fp' : None,
            'single_ok' : False,
            'no_single_ok_generic' : False
            }
        self.assertEqual(p.Params, default_params)

    def test_call_requires_read_1_file(self):
        """RtaxTaxonAssigner.__call__ functions w alt input types """
        p = RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp})

        # no read_1_seqs_fp passed results in AssertionError
        self.assertRaises(AssertionError, p, self.input_seqs_fp)

    def test_call_single_result(self):
        p = RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp,
         'read_1_seqs_fp':self.read_1_seqs_fp})

        self._paths_to_clean_up += self.cleanAll(self.read_1_seqs_fp)

        actual = p(self.input_seqs_fp)
        self.assertEqual(actual,rtax_expected_result_single)

    def test_call_paired_with_fallback_result(self):
        p = RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp,
         'read_1_seqs_fp':self.read_1_seqs_fp,
         'read_2_seqs_fp':self.read_2_seqs_fp,
         'single_ok':True})

        self._paths_to_clean_up += self.cleanAll(self.read_1_seqs_fp)
        self._paths_to_clean_up += self.cleanAll(self.read_2_seqs_fp)

        actual = p(self.input_seqs_fp)
        self.assertEqual(actual,rtax_expected_result_paired_with_fallback)

    def test_call_paired_result(self):
        p = RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp,
         'read_1_seqs_fp':self.read_1_seqs_fp,
         'read_2_seqs_fp':self.read_2_seqs_fp})

        self._paths_to_clean_up += self.cleanAll(self.read_1_seqs_fp)
        self._paths_to_clean_up += self.cleanAll(self.read_2_seqs_fp)

        actual = p(self.input_seqs_fp)
        self.assertEqual(actual,rtax_expected_result_paired)

    def test_call_output_to_file(self):
        """RtaxTaxonAssigner.__call__ functions w output to file
        """
        result_path = get_tmp_filename(
            prefix='RtaxTaxonAssignerTests_', suffix='.fasta')
        self._paths_to_clean_up.append(result_path)

        p =  RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp,
         'read_1_seqs_fp':self.read_1_seqs_fp})

        self._paths_to_clean_up += self.cleanAll(self.read_1_seqs_fp)

        actual = p(self.input_seqs_fp, result_path=result_path)

        f = open(result_path)
        observed_lines = set(f.readlines())
        f.close()
        self.assertEqual(observed_lines, rtax_expected_result_single_lines)

        # Return value is None when result_path is provided (Not sure
        # if this is what we want yet, or if we would want both so
        # results could be logged to file...)
        self.assertEqual(actual, None)

    def test_call_logs_run(self):
        """RtaxTaxonAssigner.__call__ logs the run when expected
        """
        log_path = get_tmp_filename(
            prefix='RtaxTaxonAssignerTests_',suffix='.fasta')
        self._paths_to_clean_up.append(log_path)

        p =  RtaxTaxonAssigner({
         'reference_sequences_fp':self.reference_seqs_fp,
         'id_to_taxonomy_fp':self.id_to_taxonomy_fp,
         'read_1_seqs_fp':self.read_1_seqs_fp})

        self._paths_to_clean_up += self.cleanAll(self.read_1_seqs_fp)

        actual = p(self.input_seqs_fp,log_path=log_path)

        log_file = open(log_path)
        log_file_str = log_file.read()
        log_file.close()
        # stderr.write(log_file_str)

        log_file_exp = [
            "RtaxTaxonAssigner parameters:",
            "Application:RTAX classifier",
            "Citation:Soergel D.A.W., Dey N., Knight R., and Brenner S.E.  2012.  Selection of primers for optimal taxonomic classification of environmental 16S rRNA gene sequences.  ISME J (6), 1440-1444",
            "amplicon_id_regex:(\S+)\s+(\S+?)\/",
            "header_id_regex:\S+\s+(\S+?)\/",
            "id_to_taxonomy_fp:%s" % self.id_to_taxonomy_fp,
            "read_1_seqs_fp:%s" % self.read_1_seqs_fp,
            "read_2_seqs_fp:None",
            "read_id_regex:\S+\s+(\S+)",
            "reference_sequences_fp:%s" % self.reference_seqs_fp,
            "single_ok:False",
            "no_single_ok_generic:False"
         ]
        # compare data in log file to fake expected log file
        # NOTE: Since p.params is a dict, the order of lines is not
        # guaranteed, so testing is performed to make sure that
        # the equal unordered lists of lines is present in actual and expected
        self.assertEqualItems(log_file_str.split('\n')[0:12], log_file_exp)


class MothurTaxonAssignerTests(TestCase):
    """Tests for the Mothur taxon assigner.
    """
    def setUp(self):
        def tfp(suffix):
            return get_tmp_filename(
                prefix='MothurTaxonAssigner_',
                suffix=suffix,
                result_constructor=str,
                )
        
        tax_fp = tfp('.txt')
        f = open(tax_fp, "w")
        f.write(rdp_id_to_taxonomy)
        f.close()

        ref_fp = tfp('.fna')
        g = open(ref_fp, "w")
        g.write(rdp_reference_seqs)
        g.close()

        self.seq_fp1 = tfp(".fasta")
        h = open(self.seq_fp1, "w")
        h.write(rdp_test1_fasta)
        h.close()

        self.params = {
            'id_to_taxonomy_fp': tax_fp,
            'reference_sequences_fp': ref_fp,
            }
        
        self._paths_to_clean_up = [
            tax_fp, ref_fp, self.seq_fp1,
            ]

    def tearDown(self):
        remove_files(self._paths_to_clean_up)

    def test_assignment_minimum_conf(self):
        self.params["Confidence"] = 0.00
        assigner = MothurTaxonAssigner(self.params)
        result = assigner(self.seq_fp1)

        x_lineage, x_conf = result['X67228']
        self.assertEqual(x_lineage, [
            'Bacteria', 'Proteobacteria', 'Alphaproteobacteria',
            'Rhizobiales', 'Rhizobiaceae', 'Rhizobium',
            ])
        self.assertTrue(x_conf > 0.5)
        
        e_lineage, e_conf = result['EF503697']
        self.assertEqual(e_lineage, [
            'Bacteria', 'Proteobacteria', 'Gammaproteobacteria',
            'Vibrionales', 'Vibrionaceae', 'Vibrio'
            ])
        self.assertTrue(e_conf < 0.5)

    def test_assignment_maximum_conf(self):
        self.params["Confidence"] = 0.95
        assigner = MothurTaxonAssigner(self.params)
        result = assigner(self.seq_fp1)

        x_lineage, x_conf = result['X67228']
        self.assertEqual(x_lineage, [
            'Bacteria', 'Proteobacteria', 'Alphaproteobacteria',
            'Rhizobiales', 'Rhizobiaceae', 'Rhizobium',
            ])
        self.assertTrue(x_conf > 0.94)
        
        e_lineage, e_conf = result['EF503697']
        self.assertTrue(len(e_lineage) < 3)

    def test_unassignable(self):
        f = open(self.seq_fp1, "w")
        f.write(
            ">MostlyTs\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")
        f.close()
        
        assigner = MothurTaxonAssigner(self.params)
        result = assigner(self.seq_fp1)
        self.assertEqual(result, {"MostlyTs": (["Unknown"], 0.0)})


class RdpTaxonAssignerTests(TestCase):
    """Tests for the Rdp-based taxonomy assigner.

    Slow tests are a bug, and currently these tests take about 38s on
    a Dual 2.3GHz Mac.  Presumably most of the time is spent
    initializing the Java VM on each run.  If so, this problem should
    be fixed upstream in PyCogent.
    """

    def setUp(self):
        # Temporary input file
        self.tmp_seq_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTest_',
            suffix='.fasta'
            )
        seq_file = open(self.tmp_seq_filepath, 'w')
        seq_file.write(rdp_test1_fasta)
        seq_file.close()

        # Temporary results filename
        self.tmp_res_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTestResult_',
            suffix='.tsv',
            )
        # touch the file so we don't get an error trying to close it
        open(self.tmp_res_filepath,'w').close()

        # Temporary log filename
        self.tmp_log_filepath = get_tmp_filename(
            prefix='RdpTaxonAssignerTestLog_',
            suffix='.txt',
            )
        # touch the file so we don't get an error trying to close it
        open(self.tmp_log_filepath,'w').close()

        self._paths_to_clean_up = \
         [self.tmp_seq_filepath, self.tmp_res_filepath, self.tmp_log_filepath]

        self.id_to_taxonomy_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.txt')
        self.id_to_taxonomy_file.write(rdp_id_to_taxonomy)
        self.id_to_taxonomy_file.seek(0)

        self.reference_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.fasta')
        self.reference_seqs_file.write(rdp_reference_seqs)
        self.reference_seqs_file.seek(0)

        self.default_app = RdpTaxonAssigner({})

    def tearDown(self):
        remove_files(self._paths_to_clean_up)

    def test_init(self):
        """RdpTaxonAssigner.__init__ should set default attributes and params
        """
        self.assertEqual(self.default_app.Name, 'RdpTaxonAssigner')

    def test_train_on_the_fly(self):
        """Training on-the-fly classifies reference sequence correctly with 100% certainty
        """
        input_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.fasta')
        input_seqs_file.write(test_seq_coll.toFasta())
        input_seqs_file.seek(0)

        exp_assignments = rdp_trained_test1_expected_dict

        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': self.id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        obs_assignments = app(self.tmp_seq_filepath)

        key = 'X67228 some description'
        self.assertEqual(obs_assignments[key], exp_assignments[key])

    def test_taxa_with_special_characters(self):
        """Special characters in taxa do not cause RDP errors
        """
        taxonomy_fp = get_tmp_filename()
        f = open(taxonomy_fp, "w")
        f.write(rdp_id_to_taxonomy_special_chars)
        f.close()
        self._paths_to_clean_up.append(taxonomy_fp)

        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': taxonomy_fp,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        res = app(self.tmp_seq_filepath)
        obs_lineage, obs_confidence = res['X67228 some description']

        self.assertEqual(obs_lineage, (
            "Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales<What;"
            "Rhizobiaceae&Huh?;Rhizobium"))
        self.assertEqual(obs_confidence, 1.0)

    def test_train_on_the_fly_low_memory(self):
        """Training on-the-fly with lower heap size classifies reference sequence correctly with 100% certainty
        """
        input_seqs_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.fasta')
        input_seqs_file.write(test_seq_coll.toFasta())
        input_seqs_file.seek(0)

        exp_assignments = rdp_trained_test1_expected_dict

        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': self.id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                'max_memory': '75M'
                })
        obs_assignments = app(self.tmp_seq_filepath)

        key = 'X67228 some description'
        self.assertEqual(obs_assignments[key], exp_assignments[key])

    def test_generate_training_files(self):
        app = RdpTaxonAssigner({
                'id_to_taxonomy_fp': self.id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        actual_taxonomy_file, actual_training_seqs_file = \
            app._generate_training_files()

        # see note in test_build_tree()
        self.assertEqual(actual_taxonomy_file.read(), rdp_expected_taxonomy)

    def test_call_result_as_dict(self):
        """RdpTaxonAssigner should return correct taxonomic assignment

           This test may periodically fail, but should be rare.

        """
        exp_assignments = rdp_test1_expected_dict
        min_confidence = self.default_app.Params['Confidence']

        # Since there is some variation in the assignments, run
        # 10 trials and make sure we get the expected result at least once
        num_trials = 10
        unverified_seq_ids = set(exp_assignments.keys())
        for i in range(num_trials):
            obs_assignments = self.default_app(self.tmp_seq_filepath)
            for seq_id in list(unverified_seq_ids):
                obs_assignment, obs_confidence = obs_assignments[seq_id]
                exp_assignment, exp_confidence = exp_assignments[seq_id]
                self.assertTrue(obs_confidence >= min_confidence)
                if obs_assignment == exp_assignment:
                    unverified_seq_ids.remove(seq_id)
            if not unverified_seq_ids:
                break

        messages = []
        for seq_id in unverified_seq_ids:
            messages.append(
                "Unable to verify %s in %s trials" % (seq_id, num_trials))
            messages.append("  Expected: %s" % exp_assignments[seq_id][0])
            messages.append("  Observed: %s" % obs_assignments[seq_id][0])
            messages.append("  Confidence: %s" % obs_assignments[seq_id][1])

        # make sure all taxonomic results were correct at least once
        self.assertFalse(unverified_seq_ids, msg='\n'.join(messages))

    def test_call_with_properties_file(self):
        """RdpTaxonAssigner should return correct taxonomic assignment

           This test may periodically fail, but should be rare.

        """
        id_to_taxonomy_file = NamedTemporaryFile(
            prefix='RdpTaxonAssignerTest_', suffix='.txt')
        id_to_taxonomy_file.write(rdp_id_to_taxonomy2)
        id_to_taxonomy_file.seek(0)

        app1 = RdpTaxonAssigner({
                'id_to_taxonomy_fp': id_to_taxonomy_file.name,
                'reference_sequences_fp': self.reference_seqs_file.name,
                })
        taxonomy_file, training_seqs_file = app1._generate_training_files()

        training_dir = mkdtemp(prefix='RdpTrainer_')
        training_results = train_rdp_classifier(
            training_seqs_file, taxonomy_file, training_dir)

        training_data_fp = training_results['properties'].name
        min_confidence = 0.80
        app2 = RdpTaxonAssigner({
            'training_data_properties_fp': training_data_fp,
            'Confidence': min_confidence,
            })

        expected = rdp_trained_test2_expected_dict

        # Since there is some variation in the assignments, run
        # 10 trials and make sure we get the expected result at least once
        num_trials = 10
        num_seqs = len(expected)
        seq_ids = expected.keys()
        assignment_comp_results = [False] * num_seqs
        expected_assignment_comp_results = [True] * num_seqs

        for i in range(num_trials):
            actual = app2(self.tmp_seq_filepath)
            # seq ids are the same, and all input sequences get a result
            self.assertEqual(actual.keys(),expected.keys())
            for j,seq_id in enumerate(seq_ids):
                # confidence is above threshold
                self.assertTrue(actual[seq_id][1] >= min_confidence)
                # confidence roughly matches expected
                self.assertFloatEqual(\
                 actual[seq_id][1],expected[seq_id][1],0.1)
                # check if the assignment is correct -- this must happen
                # at least once per seq_id for the test to pass
                if actual[seq_id][0] == expected[seq_id][0]:
                    assignment_comp_results[j] = True
            if assignment_comp_results == expected_assignment_comp_results:
                # break once we've seen a correct assignment for each seq
                break

        self.assertEqual(\
         assignment_comp_results,\
         expected_assignment_comp_results,\
         "Taxonomic assignments never correct in %d trials." % num_trials)

    def test_call_result_to_file(self):
        """RdpTaxonAssigner should save results to file

           This test may periodically fail, but should be rare.

        """
        expected_lines = rdp_test1_expected_lines

        # Since there is some variation in the assignments, run
        # 10 trials and make sure we get the expected result at least once
        # for each sequence
        num_trials = 10
        num_seqs = len(expected_lines)
        assignment_comp_results = [False] * num_seqs
        expected_assignment_comp_results = [True] * num_seqs

        for i in range(num_trials):
            retval = self.default_app(
             seq_path=self.tmp_seq_filepath,
             result_path=self.tmp_res_filepath,
             log_path=None)
            actual = [l.strip() for l in open(self.tmp_res_filepath, 'r')]
            message = "Expected return value of None but observed %s" % retval
            self.assertTrue(retval is None, message)
            for j in range(num_seqs):
                a = actual[j]
                e = expected_lines[j]
                # note we're testing using startswith here to allow
                # for some variability in confidence
                if a.startswith(e):
                    assignment_comp_results[j] = True
            if assignment_comp_results == expected_assignment_comp_results:
                break

        self.assertEqual(\
         assignment_comp_results,\
         expected_assignment_comp_results,\
         "Taxonomic assignments never correct in %d trials." % num_trials)

    def test_log(self):
        """RdpTaxonAssigner should write correct message to log file"""
        # expected result when no result_path is provided
        self.default_app(
            seq_path=self.tmp_seq_filepath,
            result_path=None,
            log_path=self.tmp_log_filepath,
            )

        # open the actual log file and the expected file, and pass into lists
        obs = [l.strip() for l in list(open(self.tmp_log_filepath, 'r'))]
        exp = rdp_test1_log_file_contents.split('\n')
        # sort the lists as the entries are written from a dict,
        # so order may vary
        obs.sort()
        exp.sort()
        self.assertEqual(obs, exp)

class RdpTrainingSetTests(TestCase):
    def setUp(self):
        self.tagged_str = (
            '>a1\tA;B;C' + _QIIME_RDP_TAXON_TAG +
            '\tD_' + _QIIME_RDP_TAXON_TAG +
            ';E;F' + _QIIME_RDP_TAXON_TAG + '\n' +
            'GGGCCC\n'
            )
        self.untagged_str = '>a1\tA;B;C\tD_;E;F\nGGGCCC\n'

    def test_add_sequence(self):
        s = RdpTrainingSet()
        s.add_sequence('a1', 'GGCCTT')
        self.assertEqual(s.sequences['a1'], 'GGCCTT')

    def test_add_lineage(self):
        s = RdpTrainingSet()
        s.add_lineage('a1', 'A;B;C;D;E;F')
        n = s.sequence_nodes['a1']
        self.assertEqual(n.name, 'F')

        # Raise a ValueError if lineage does not have 6 ranks
        self.assertRaises(ValueError, s.add_lineage, 'a2', 'A;B;C')

    def test_get_training_seqs(self):
        s = RdpTrainingSet()
        s.add_sequence('a1', 'GGCCTT')
        s.add_sequence('b1', 'CCCCGG')
        s.add_lineage('a1', 'A;B;C;D;E;F')
        s.add_lineage('c1', 'G;H;I;J;K;L')

        # Neither b1 or c1 appear, since they do not have both
        # sequence and lineage.
        obs = s.get_training_seqs()
        self.assertEqual(list(obs), [("a1 Root;A;B;C;D;E;F", "GGCCTT")])

    def test_get_rdp_taxonomy(self):
        s = RdpTrainingSet()
        s.add_lineage('a1', 'A;B;C;D;E;F')
        s.add_lineage('c1', 'A;B;I;J;K;L')

        # All taxa appear, regardless of whether a sequence was
        # registered.
        expected = (
            '0*Root*0*0*norank\n'
            '1*A*0*1*a\n'
            '2*B*1*2*b\n'
            '3*C*2*3*c\n'
            '4*D*3*4*d\n'
            '5*E*4*5*e\n'
            '6*F*5*6*genus\n'
            '7*I*2*3*c\n'
            '8*J*7*4*d\n'
            '9*K*8*5*e\n'
            '10*L*9*6*genus\n'
            )
        self.assertEqual(s.get_rdp_taxonomy(), expected)

    def test_fix_output_file(self):
        fp = get_tmp_filename()
        open(fp, 'w').write(self.tagged_str)

        s = RdpTrainingSet()
        s.fix_output_file(fp)
        obs = open(fp).read()
        remove(fp)

        self.assertEqual(obs, self.untagged_str)

    def test_fix_results(self):
        tagged_lineage = self.tagged_str[4:]
        untagged_lineage = self.untagged_str[4:]
        
        s = RdpTrainingSet()
        obs = s.fix_results({'a1': (tagged_lineage, 1.00)})
        self.assertEqual(obs, {'a1': (untagged_lineage, 1.00)})


class RdpTreeTests(TestCase):
    def test_insert_lineage(self):
        t = RdpTree()
        t.insert_lineage(['a', 'b', 'c'])

        self.assertEqual(t.children.keys(), ['a'])
        self.assertEqual(t.children['a'].children['b'].children['c'].name, 'c')

    def test_get_lineage(self):
        t = RdpTree()
        cnode = t.insert_lineage(['a', 'b', 'c'])
        self.assertEqual(cnode.get_lineage(), ['Root', 'a', 'b', 'c'])

    def test_get_nodes(self):
        t = RdpTree()
        t.insert_lineage(['a', 'b', 'c'])
        obs_names = [n.name for n in t.get_nodes()]
        self.assertEqual(obs_names, ['Root', 'a', 'b', 'c'])

    def test_dereplicate_taxa(self):
        t = RdpTree()
        t.insert_lineage(['a', 'x'])
        t.insert_lineage(['b', 'x'])
        t.dereplicate_taxa()

        c1 = t.children['a'].children['x']
        c2 = t.children['b'].children['x']
        self.assertNotEqual(c1.name, c2.name)

        # Test case-insensitive dereplication.
        t = RdpTree()
        t.insert_lineage(['a', 'x'])
        t.insert_lineage(['b', 'X'])
        t.dereplicate_taxa()

        c1 = t.children['a'].children['x']
        c2 = t.children['b'].children['X']
        self.assertTrue((c1.name == 'x' and c2.name != 'X') or
                        (c1.name != 'x' and c2.name == 'X'))

    def test_get_rdp_taxonomy(self):
        t = RdpTree()
        t.insert_lineage(['A', 'B', 'C', 'D', 'E', 'F'])
        expected = (
            '0*Root*0*0*norank\n'
            '1*A*0*1*a\n'
            '2*B*1*2*b\n'
            '3*C*2*3*c\n'
            '4*D*3*4*d\n'
            '5*E*4*5*e\n'
            '6*F*5*6*genus\n'
            )
        self.assertEqual(t.get_rdp_taxonomy(), expected)


rdp_test1_fasta = \
""">X67228 some description
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>EF503697
TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC
"""

rdp_test1_log_file_contents = \
"""RdpTaxonAssigner parameters:
Application:RDP classfier
Citation:Wang, Q, G. M. Garrity, J. M. Tiedje, and J. R. Cole. 2007. Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl Environ Microbiol. 73(16):5261-7.
Taxonomy:RDP
Confidence:0.8
id_to_taxonomy_fp:None
reference_sequences_fp:None
training_data_properties_fp:None
max_memory:None"""

rdp_test1_expected_dict = {
    'X67228 some description': (
        'Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales', 0.95),
    'EF503697': (
        'Archaea;Crenarchaeota;Thermoprotei', 0.88),
    }

rdp_test1_expected_lines = [\
 "\t".join(["X67228 some description",\
  "Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales",\
  "0."]),
 "\t".join(['EF503697','Archaea;Crenarchaeota;Thermoprotei','0.'])]

rdp_trained_test1_expected_dict = {
    'X67228 some description': ('Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium', 1.0),
    'EF503697': ('Bacteria;Proteobacteria', 0.93000000000000005),
    }

rdp_trained_test2_expected_dict = {
    'X67228 some description': ('Bacteria;Proteobacteria;Alphaproteobacteria2;Rhizobiales;Rhizobiaceae;Rhizobium', 1.0),
    'EF503697': ('Bacteria;Proteobacteria', 0.93000000000000005),
    }

rdp_id_to_taxonomy = \
"""X67228	Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium
X73443	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium
AB004750	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
xxxxxx	Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
AB004748	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
AB000278	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium
AB000390	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio
"""

rdp_id_to_taxonomy2 = \
"""X67228	Bacteria;Proteobacteria;Alphaproteobacteria2;Rhizobiales;Rhizobiaceae;Rhizobium
X73443	Bacteria;Firmicutes;Clostridia2;Clostridiales;Clostridiaceae;Clostridium
AB004750	Bacteria;Proteobacteria;Gammaproteobacteria2;Enterobacteriales;Enterobacteriaceae;Enterobacter
xxxxxx	Bacteria;Proteobacteria;Gammaproteobacteria2;Pseudomonadales;Pseudomonadaceae;Pseudomonas
AB004748	Bacteria;Proteobacteria;Gammaproteobacteria2;Enterobacteriales;Enterobacteriaceae;Enterobacter
AB000278	Bacteria;Proteobacteria;Gammaproteobacteria2;Vibrionales;Vibrionaceae;Photobacterium
AB000390	Bacteria;Proteobacteria;Gammaproteobacteria2;Vibrionales;Vibrionaceae;Vibrio
"""

rdp_id_to_taxonomy_special_chars = \
"""X67228	Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales<What;Rhizobiaceae&Huh?;Rhizobium
X73443	Bacteria;Firmicutes;Clostridia2;Clostridiales;Clostridiaceae;Clostridium
AB004750	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
xxxxxx	Bacteria;Proteobacteria;Gammaproteobacteria2;Pseudomonadales;Pseudomonadaceae;Pseudomonas
AB004748	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
AB000278	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrio>nales;Vibrionaceae;Photobacterium
AB000390	Bacteria;Proteobacteria;Gammaproteobacteria;Vibrio>nales;Vibrionaceae;Vibrio
"""

rdp_reference_seqs = \
""">X67228
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>X73443
nnnnnnngagatttgatcctggctcaggatgaacgctggccggccgtgcttacacatgcagtcgaacgaagcgcttaaactggatttcttcggattgaagtttttgctgactgagtggcggacgggtgagtaacgcgtgggtaacctgcctcatacagggggataacagttagaaatgactgctaataccnnataagcgcacagtgctgcatggcacagtgtaaaaactccggtggtatgagatggacccgcgtctgattagctagttggtggggt
>AB004750
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgca
>xxxxxx
ttgaacgctggcggcaggcctaacacatgcaagtcgagcggcagcannnncttcgggaggctggcgagcggcggacgggtgagtaacgcatgggaacttacccagtagtgggggatagcccggggaaacccggattaataccgcatacgccctgagggggaaagcgggctccggtcgcgctattggatgggcccatgtcggattagttagttggtggggtaatggcctaccaaggcgacgatccgtagctggtctgagaggatgatcagccacaccgggactgagacacggcccggactcctacgggaggcagcagtggggaatattggacaatgggggcaaccctgatccagccatgccg
>AB004748
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggccttcgggttg
>AB000278
caggcctaacacatgcaagtcgaacggtaanagattgatagcttgctatcaatgctgacgancggcggacgggtgagtaatgcctgggaatataccctgatgtgggggataactattggaaacgatagctaataccgcataatctcttcggagcaaagagggggaccttcgggcctctcgcgtcaggattagcccaggtgggattagctagttggtggggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggggaaaccctgatgcagccatgccgcgtgta
>AB000390
tggctcagattgaacgctggcggcaggcctaacacatgcaagtcgagcggaaacgantnntntgaaccttcggggnacgatnacggcgtcgagcggcggacgggtgagtaatgcctgggaaattgccctgatgtgggggataactattggaaacgatagctaataccgcataatgtctacggaccaaagagggggaccttcgggcctctcgcttcaggatatgcccaggtgggattagctagttggtgaggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgag
"""


rdp_expected_taxonomy = """\
0*Root*0*0*norank
1*Bacteria*0*1*a
7*Firmicutes*1*2*b
8*Clostridia*7*3*c
9*Clostridiales*8*4*d
10*Clostridiaceae*9*5*e
11*Clostridium*10*6*genus
2*Proteobacteria*1*2*b
3*Alphaproteobacteria*2*3*c
4*Rhizobiales*3*4*d
5*Rhizobiaceae*4*5*e
6*Rhizobium*5*6*genus
12*Gammaproteobacteria*2*3*c
13*Enterobacteriales*12*4*d
14*Enterobacteriaceae*13*5*e
15*Enterobacter*14*6*genus
16*Pseudomonadales*12*4*d
17*Pseudomonadaceae*16*5*e
18*Pseudomonas*17*6*genus
19*Vibrionales*12*4*d
20*Vibrionaceae*19*5*e
21*Photobacterium*20*6*genus
22*Vibrio*20*6*genus
"""

# newline at the end makes a difference
rdp_expected_training_seqs = \
""">AB000278 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Photobacterium
caggcctaacacatgcaagtcgaacggtaanagattgatagcttgctatcaatgctgacgancggcggacgggtgagtaatgcctgggaatataccctgatgtgggggataactattggaaacgatagctaataccgcataatctcttcggagcaaagagggggaccttcgggcctctcgcgtcaggattagcccaggtgggattagctagttggtggggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggggaaaccctgatgcagccatgccgcgtgta
>AB000390 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Vibrionales;Vibrionaceae;Vibrio
tggctcagattgaacgctggcggcaggcctaacacatgcaagtcgagcggaaacgantnntntgaaccttcggggnacgatnacggcgtcgagcggcggacgggtgagtaatgcctgggaaattgccctgatgtgggggataactattggaaacgatagctaataccgcataatgtctacggaccaaagagggggaccttcgggcctctcgcttcaggatatgcccaggtgggattagctagttggtgaggtaatggctcaccaaggcgacgatccctagctggtctgagaggatgatcagccacactggaactgag
>AB004748 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgcacaatgggcgcaagcctgatgcagccatgccgcgtgtatgaagaaggccttcgggttg
>AB004750 Root;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Enterobacter
acgctggcggcaggcctaacacatgcaagtcgaacggtagcagaaagaagcttgcttctttgctgacgagtggcggacgggtgagtaatgtctgggaaactgcccgatggagggggataactactggaaacggtagctaataccgcataacgtcttcggaccaaagagggggaccttcgggcctcttgccatcggatgtgcccagatgggattagctagtaggtggggtaacggctcacctaggcgacgatccctagctggtctgagaggatgaccagccacactggaactgagacacggtccagactcctacgggaggcagcagtggggaatattgca
>X67228 Root;Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Rhizobium
aacgaacgctggcggcaggcttaacacatgcaagtcgaacgctccgcaaggagagtggcagacgggtgagtaacgcgtgggaatctacccaaccctgcggaatagctctgggaaactggaattaataccgcatacgccctacgggggaaagatttatcggggatggatgagcccgcgttggattagctagttggtggggtaaaggcctaccaaggcgacgatccatagctggtctgagaggatgatcagccacattgggactgagacacggcccaaa
>X73443 Root;Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium
nnnnnnngagatttgatcctggctcaggatgaacgctggccggccgtgcttacacatgcagtcgaacgaagcgcttaaactggatttcttcggattgaagtttttgctgactgagtggcggacgggtgagtaacgcgtgggtaacctgcctcatacagggggataacagttagaaatgactgctaataccnnataagcgcacagtgctgcatggcacagtgtaaaaactccggtggtatgagatggacccgcgtctgattagctagttggtggggt
>xxxxxx Root;Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas
ttgaacgctggcggcaggcctaacacatgcaagtcgagcggcagcannnncttcgggaggctggcgagcggcggacgggtgagtaacgcatgggaacttacccagtagtgggggatagcccggggaaacccggattaataccgcatacgccctgagggggaaagcgggctccggtcgcgctattggatgggcccatgtcggattagttagttggtggggtaatggcctaccaaggcgacgatccgtagctggtctgagaggatgatcagccacaccgggactgagacacggcccggactcctacgggaggcagcagtggggaatattggacaatgggggcaaccctgatccagccatgccg"""

id_to_taxonomy_string = \
"""AY800210\tArchaea;Euryarchaeota;Halobacteriales;uncultured
EU883771\tArchaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.
EF503699\tArchaea;Crenarchaeota;uncultured;uncultured
DQ260310\tArchaea;Euryarchaeota;Methanobacteriales;Methanobacterium
EF503697\tArchaea;Crenarchaeota;uncultured;uncultured"""

test_seq_coll = LoadSeqs(data=[\
 ('s1','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('s2','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('s3','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('s4','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('s5','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC'),\
 ('s6','ATAGTAGGTGATTGCGAAGACCGCGGAACCGGGACCTAGCACCCAGCCTGTACCGAGGGATGGGGAGCTGTGGCGGTCCACCGACGACCCTTTGTGACAGCCGATTCCTACAATCCCAGCAACTGCAATGATCCACTCTAGTCGGCATAACCGGGAATCGTTAACCTGGTAGGGTTCTCTACGTCTGAGTCTACAGCCCAGAGCAGTCAGGCTACTATACGGTTTGCTGCATTGCATAGGCATCGGTCGCGGGCACTCCTCGCGGTTTCAGCTAGGGTTTAAATGGAGGGTCGCTGCATGAGTATGCAAATAGTGCCACTGCTCTGATACAGAGAAGTGTTGATATGACACCTAAGACCTGGTCACAGTTTTAACCTGCCTACGCACACCAGTGTGCTATTGATTAACGATATCGGTAGACACGACCTTGGTAACCTGACTAACCTCATGGAAAGTGACTAGATAAATGGACCGGAGCCAACTTTCACCCGGAAAACGGACCGACGAATCGTCGTAGACTACCGATCTGACAAAATAAGCACGAGGGAGCATGTTTTGCGCAGGCTAGCCTATTCCCACCTCAAGCCTCGAGAACCAAGACGCCTGATCCGGTGCTGCACGAAGGGTCGCCTCTAGGTAAGGAGAGCTGGCATCTCCAGATCCGATATTTTACCCAACCTTTGCGCGCTCAGATTGTTATAGTGAAACGATTTAAGCCTGAACGGAGTTCCGCTCCATATGTGGGTTATATATGTGAGATGTATTAACTTCCGCAGTTGTCTCTTTCGGTGCAGTACGCTTGGTATGTGTCTCAAATAATCGGTATTATAGTGATCTGAGAGGTTTTAAG')],aligned=False)

test_refseq_coll = LoadSeqs(data=[\
 ('AY800210','TTCCGGTTGATCCTGCCGGACCCGACTGCTATCCGGATGCGACTAAGCCATGCTAGTCTAACGGATCTTCGGATCCGTGGCATACCGCTCTGTAACACGTAGATAACCTACCCTGAGGTCGGGGAAACTCCCGGGAAACTGGGCCTAATCCCCGATAGATAATTTGTACTGGAATGTCTTTTTATTGAAACCTCCGAGGCCTCAGGATGGGTCTGCGCCAGATTATGGTCGTAGGTGGGGTAACGGCCCACCTAGCCTTTGATCTGTACCGGACATGAGAGTGTGTGCCGGGAGATGGCCACTGAGACAAGGGGCCAGGCCCTACGGGGCGCAGCAGGCGCGAAAACTTCACAATGCCCGCAAGGGTGATGAGGGTATCCGAGTGCTACCTTAGCCGGTAGCTTTTATTCAGTGTAAATAGCTAGATGAATAAGGGGAGGGCAAGGCTGGTGCCAGCCGCCGCGGTAAAACCAGCTCCCGAGTGGTCGGGATTTTTATTGGGCCTAAAGCGTCCGTAGCCGGGCGTGCAAGTCATTGGTTAAATATCGGGTCTTAAGCCCGAACCTGCTAGTGATACTACACGCCTTGGGACCGGAAGAGGCAAATGGTACGTTGAGGGTAGGGGTGAAATCCTGTAATCCCCAACGGACCACCGGTGGCGAAGCTTGTTCAGTCATGAACAACTCTACACAAGGCGATTTGCTGGGACGGATCCGACGGTGAGGGACGAAACCCAGGGGAGCGAGCGGGATTAGATACCCCGGTAGTCCTGGGCGTAAACGATGCGAACTAGGTGTTGGCGGAGCCACGAGCTCTGTCGGTGCCGAAGCGAAGGCGTTAAGTTCGCCGCCAGGGGAGTACGGCCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCAC'),\
 ('EU883771','TGGCGTACGGCTCAGTAACACGTGGATAACTTACCCTTAGGACTGGGATAACTCTGGGAAACTGGGGATAATACTGGATATTAGGCTATGCCTGGAATGGTTTGCCTTTGAAATGTTTTTTTTCGCCTAAGGATAGGTCTGCGGCTGATTAGGTCGTTGGTGGGGTAATGGCCCACCAAGCCGATGATCGGTACGGGTTGTGAGAGCAAGGGCCCGGAGATGGAACCTGAGACAAGGTTCCAGACCCTACGGGGTGCAGCAGGCGCGAAACCTCCGCAATGTACGAAAGTGCGACGGGGGGATCCCAAGTGTTATGCTTTTTTGTATGACTTTTCATTAGTGTAAAAAGCTTTTAGAATAAGAGCTGGGCAAGACCGGTGCCAGCCGCCGCGGTAACACCGGCAGCTCGAGTGGTGACCACTTTTATTGGGCTTAAAGCGTTCGTAGCTTGATTTTTAAGTCTCTTGGGAAATCTCACGGCTTAACTGTGAGGCGTCTAAGAGATACTGGGAATCTAGGGACCGGGAGAGGTAAGAGGTACTTCAGGGGTAGAAGTGAAATTCTGTAATCCTTGAGGGACCACCGATGGCGAAGGCATCTTACCAGAACGGCTTCGACAGTGAGGAACGAAAGCTGGGGGAGCGAACGGGATTAGATACCCCGGTAGTCCCAGCCGTAAACTATGCGCGTTAGGTGTGCCTGTAACTACGAGTTACCGGGGTGCCGAAGTGAAAACGTGAAACGTGCCGCCTGGGAAGTACGGTCGCAAGGCTGAAACTTAAAGGAATTGGCGGGGGAGCACCACAACGGGTGGAGCCTGCGGTTTAATTGGACTCAACGCCGGGCAGCTCACCGGATAGGACAGCGGAATGATAGCCGGGCTGAAGACCTTGCTTGACCAGCTGAGA'),\
 ('EF503699','AAGAATGGGGATAGCATGCGAGTCACGCCGCAATGTGTGGCATACGGCTCAGTAACACGTAGTCAACATGCCCAGAGGACGTGGACACCTCGGGAAACTGAGGATAAACCGCGATAGGCCACTACTTCTGGAATGAGCCATGACCCAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGGAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCACGAAACCTCTGCAATAGGCGAAAGCTTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCCGCTTAACGGATGGGCTGCGGAGGATACTGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCTTTGATCTACTGAAGACCACCAGTGGTGAAGGCGGTTCGCCAGAACGCGCTCGAACGGTGAGGATGAAAGCTGGGGGAGCAAACCGGAATAGATACCCGAGTAATCCCAACTGTAAACGATGGCAACTCGGGGATGGGTTGGCCTCCAACCAACCCCATGGCCGCAGGGAAGCCGTTTAGCTCTCCCGCCTGGGGAATACGGTCCGCAGAATTGAACCTTAAAGGAATTTGGCGGGGAACCCCCACAAGGGGGAAAACCGTGCGGTTCAATTGGAATCCACCCCCCGGAAACTTTACCCGGGCGCG'),\
 ('DQ260310','GATACCCCCGGAAACTGGGGATTATACCGGATATGTGGGGCTGCCTGGAATGGTACCTCATTGAAATGCTCCCGCGCCTAAAGATGGATCTGCCGCAGAATAAGTAGTTTGCGGGGTAAATGGCCACCCAGCCAGTAATCCGTACCGGTTGTGAAAACCAGAACCCCGAGATGGAAACTGAAACAAAGGTTCAAGGCCTACCGGGCACAACAAGCGCCAAAACTCCGCCATGCGAGCCATCGCGACGGGGGAAAACCAAGTACCACTCCTAACGGGGTGGTTTTTCCGAAGTGGAAAAAGCCTCCAGGAATAAGAACCTGGGCCAGAACCGTGGCCAGCCGCCGCCGTTACACCCGCCAGCTCGAGTTGTTGGCCGGTTTTATTGGGGCCTAAAGCCGGTCCGTAGCCCGTTTTGATAAGGTCTCTCTGGTGAAATTCTACAGCTTAACCTGTGGGAATTGCTGGAGGATACTATTCAAGCTTGAAGCCGGGAGAAGCCTGGAAGTACTCCCGGGGGTAAGGGGTGAAATTCTATTATCCCCGGAAGACCAACTGGTGCCGAAGCGGTCCAGCCTGGAACCGAACTTGACCGTGAGTTACGAAAAGCCAAGGGGCGCGGACCGGAATAAAATAACCAGGGTAGTCCTGGCCGTAAACGATGTGAACTTGGTGGTGGGAATGGCTTCGAACTGCCCAATTGCCGAAAGGAAGCTGTAAATTCACCCGCCTTGGAAGTACGGTCGCAAGACTGGAACCTAAAAGGAATTGGCGGGGGGACACCACAACGCGTGGAGCCTGGCGGTTTTATTGGGATTCCACGCAGACATCTCACTCAGGGGCGACAGCAGAAATGATGGGCAGGTTGATGACCTTGCTTGACAAGCTGAAAAGGAGGTGCAT'),\
 ('EF503697','TAAAATGACTAGCCTGCGAGTCACGCCGTAAGGCGTGGCATACAGGCTCAGTAACACGTAGTCAACATGCCCAAAGGACGTGGATAACCTCGGGAAACTGAGGATAAACCGCGATAGGCCAAGGTTTCTGGAATGAGCTATGGCCGAAATCTATATGGCCTTTGGATTGGACTGCGGCCGATCAGGCTGTTGGTGAGGTAATGGCCCACCAAACCTGTAACCGGTACGGGCTTTGAGAGAAGTAGCCCGGAGATGGGCACTGAGACAAGGGCCCAGGCCCTATGGGGCGCAGCAGGCGCGAAACCTCTGCAATAGGCGAAAGCCTGACAGGGTTACTCTGAGTGATGCCCGCTAAGGGTATCTTTTGGCACCTCTAAAAATGGTGCAGAATAAGGGGTGGGCAAGTCTGGTGTCAGCCGCCGCGGTAATACCAGCACCCCGAGTTGTCGGGACGATTATTGGGCCTAAAGCATCCGTAGCCTGTTCTGCAAGTCCTCCGTTAAATCCACCTGCTCAACGGATGGGCTGCGGAGGATACCGCAGAGCTAGGAGGCGGGAGAGGCAAACGGTACTCAGTGGGTAGGGGTAAAATCCATTGATCTACTGAAGACCACCAGTGGCGAAGGCGGTTTGCCAGAACGCGCTCGACGGTGAGGGATGAAAGCTGGGGGAGCAAACCGGATTAGATACCCGGGGTAGTCCCAGCTGTAAACGGATGCAGACTCGGGTGATGGGGTTGGCTTCCGGCCCAACCCCAATTGCCCCCAGGCGAAGCCCGTTAAGATCTTGCCGCCCTGTCAGATGTCAGGGCCGCCAATACTCGAAACCTTAAAAGGAAATTGGGCGCGGGAAAAGTCACCAAAAGGGGGTTGAAACCCTGCGGGTTATATATTGTAAACC')],aligned=False)




# sample data copied from GreenGenes


rtax_reference_taxonomy = """508720	99.0	k__Bacteria	 p__Actinobacteria	 c__Actinobacteria	 o__Actinomycetales	 f__Propionibacteriaceae	 g__Propionibacterium	 s__Propionibacterium acnes
508050	99.0	k__Bacteria	 p__Proteobacteria	 c__Betaproteobacteria	 o__Burkholderiales	 f__Comamonadaceae	 g__Diaphorobacter	 s__
502492	99.0	k__Bacteria	 p__Proteobacteria	 c__Betaproteobacteria	 o__Burkholderiales	 f__	 g__Aquabacterium	 s__
"""

rtax_reference_fasta = """>508720
GACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGTGGGGTGCTCGAGTGGCGAACG
GGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGGCTAATACCGGATAGGAGCTCCT
GCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCAGCTTGTTGGTGGGGTAGTGGCTT
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
GGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCGTGCGGGATGACGGCCTTCGGGTTGT
AAACCGCTTTCGCCTGTGACGAAGCGTGAGTGACGGTAATGGGTAAAGAAGCACCGGCTAACTACGTGCCAGCAGCCGCG
GTGATACGTAGGGTGCGAGCGTTGTCCGGATTTATTGGGCGTAAAGGGCTCGTAGGTGGTTGATCGCGTCGGAAGTGTAA
TCTTGGGGCTTAACCCTGAGCGTGCTTTCGATACGGGTTGACTTGAGGAAGGTAGGGGAGAATGGAATTCCTGGTGGAGC
GGTGGAATGCGCAGATATCAGGAGGAACACCAGTGGCGAAGGCGGTTCTCTGGGCCTTTCCTGACGCTGAGGAGCGAAAG
CGTGGGGAGCGAACAGGCTTAGATACCCTGGTAGTCCACGCTGTAAACGGTGGGTACTAGGTGTGGGGTCCATTCCACGG
GTTCCGTGCCGTAGCTAACGCTTTAAGTACCCCGCCTGGGGAGTACGGCCGCAAGGCTAAAACTCAAAGGAATTGACGGG
GCCCCGCACAAGCGGCGGAGCATGCGGATTAATTCGATGCAACGCGTAGAACCTTACCTGGGTTTGACATGGATCGGGAG
TGCTCAGAGATGGGTGTGCCTCTTTTGGGGTCGGTTCACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTT
GGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
GTCAACTCGGAGGAAGGTGGGGATGACGTCAAGTCATCATGCCCCTTATGTCCAGGGCTTCACGCATGCTACAATGGCTG
GTACAGAGAGTGGCGAGCCTGTGAGGGTGAGCGAATCTCGGAAAGCCGGTCTCAGTTCGGATTGGGGTCTGCAACTCGAC
CTCATGAAGTCGGAGTCGCTAGTAATCGCAGATCAGCAACGCTGCGGTGAATACGTTCCCGGGGCT
>508050
ATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGTAACAGGTCTTCGGATGCTGACGAGTGGCGAACGGGTG
AGTAATACATCGGAACGTGCCCGATCGTGGGGGATAACGAGGCGAAAGCTTTGCTAATACCGCATACGATCTACGGATGA
AAGCGGGGGATCTTCGGACCTCGCGCGGACGGAGCGGCCGATGGCAGATTAGGTAGTTGGTGGGATAAAAGCTTACCAAG
CCGACGATCTGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
AGTGGGGAATTTTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGCAGGATGAAGGCCTTCGGGTTGTAAACTG
CTTTTGTACGGAACGAAAAGCCTCTTTCTAATAAAGAGGGGTCATGACGGTACCGTAAGAATAAGCACCGGCTAACTACG
TGCCAGCAGCCGCGGTAATACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTA
AGACAGAGGTGAAATCCCCGGGCTCAACCTGGGAACTGCCTTTGTGACTGCAAGGCTGGAGTGCGGCAGAGGGGGATGGA
ATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAATCCCCTGGGCCTGCACTGACG
CTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTGGTTGTTG
GGTCTTCACTGACTCAGTAACGAAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAAAG
GAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGTTTAATTCGATGCAACGCGAAAAACCTTACCCACCTTTGACA
TGGCAGGAAGTTTCCAGAGATGGATTCGTGCCCGAAAGGGAACCTGCACACAGGTGCTGCATGGCTGTCGTCAGCTCGTG
TCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
CCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAGTCCTCATGGCCCTTATAGGTGGGGCTACACACGTCATACAAT
GGCTGGTACAGAGGGTTGCCAACCCGCGAGGGGGAGCTAATCCCATAAAGCCAGTCGTAGTCCGGATCGCAGTCTGCAAC
TCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGGATCAGAATGTCGCGGTGAATACGTTCCCGGGTCT
>502492
ATTGAACGCTGGCGGCATGCCTTACACATGCAAGTCGAACGGTAACGGGTCCTTCGGGATGCCGACGAGTGGCGAACGGG
TGAGTAATATATCGGAACGTGCCCAGTAGTGGGGGATAACTGCTCGAAAGAGCAGCTAATACCGCATACGACCTGAGGGT
GAAAGGGGGGGATCGCAAGACCTCTCGCTATTGGAGCGGCCGATATCAGATTAGCTAGTTGGTGGGGTAAAGGCCTACCA
AGGCAACGATCTGTAGTTGGTCTGAGAGGACGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
GCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCAATGCCGCGTGCAGGAAGAAGGCCTTCGGGTTGTAAAC
TGCTTTTGTCAGGGAAGAAATCTTCTGGGCTAATACCCCGGGAGGATGACGGTACCTGAAGAATAAGCACCGGCTAACTA
CGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCTTTG
CAAGACAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTTGTGACTGCAAGGCTAGAGTACGGCAGAGGGGGATG
GAATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCAATGGCGAAGGCAATCCCCTGGGCCTGTACTGA
CGCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCCTAAACGATGTCAACTGGTTGT
TGGACGGCTTGCTGTTCAGTAACGAAGCTAACGCGTGAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTGAAACTCAA
AGGAATTGACGGGGACCCGCACAAGCGGTGGATGATGTGGTTTAATTCGATGCAACGCGAAAAACCTTACCTACCCTTGA
CATGTCAAGAATTCTGCAGAGATGTGGAAGTGCTCGAAAGAGAACTTGAACACAGGTGCTGCATGGCCGTCGTCAGCTCG
TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
TGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAGGTCCTCATGGCCCTTATGGGTAGGGCTACACACGTCATACA
ATGGCCGGTACAGAGGGCTGCCAACCCGCGAGGGGGAGCCAATCCCAGAAAACCGGTCGTAGTCCGGATCGTAGTCTGCA
ACTCGACTGCGTGAAGTCGGAATCGCTAGTAATCGCGGATCAGCTTGCCGCGGTGAATACGTTCCCGGGTCT
"""


rtax_test_repset_fasta = """>clusterIdA splitRead1IdA
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGTGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
>clusterIdB splitRead1IdB
CCGACGATCTGTAGCTGGTCTGAGAGGATGTTCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
>clusterIdC splitRead1IdC
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>clusterIdD splitRead1IdD
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACGGGGGGGGGGCCCAGACTCCTACGGGAGGCA
"""

# these reads are the 4th and 14th lines from the reference seqs

#rtax_test_read1_fasta = """>splitRead1IdA ampliconId_34563456/1
#ACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
#>splitRead1IdB ampliconId_
#CCGACGATCTGTAGCTGGTCTGAGAGGATGATCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
#>splitRead1IdC ampliconId_
#AGGCAACGATCTGTAGTTGGTCTGAGAGGACGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
#"""
#
#rtax_test_read2_fasta = """>splitRead2IdA ampliconId_34563456/3
#GGGTTAAGTCCCGCAACGAGCGCAACCCTTGTTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
#>splitRead2IdB ampliconId_
#TCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
#>splitRead2IdC ampliconId_
#TGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
#"""


# these reads are the 4th and 14th lines from the reference seqs, with one nucleotide changed each
# except D and E, which are based on the C reads, altered with a bunch of 'A's,
# and with IDs that are unique to one read or the other and so can't be paired.
# F and G are just decoys

rtax_test_read1_fasta = """>splitRead1IdA ampliconId_34563456/1
ACCAAGGCTTTGACGGGTAGCCGGCCTGAGTGGGTGACCGGCCACATTGGGACTGAGATACGGCCCAGACTCCTACGGGA
>splitRead1IdB ampliconId_12341234/1
CCGACGATCTGTAGCTGGTCTGAGAGGATGTTCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGC
>splitRead1IdC ampliconId_23452345/1
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAGCCACACTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>splitRead1IdD ampliconId_45674567/1
AGGCAACGATCTGTAGTTGGTCTGAGAGGAGGACCAAAAAAAAAAAGACTGAGACACGGCCCAGACTCCTACGGGAGGCA
>splitRead1IdF ampliconId_56785678/1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
"""

rtax_test_read2_fasta = """>splitRead2IdA ampliconId_34563456/3
GGGTTAAGTCCCGCAACGAGCGCAACCCTTATTCACTGTTGCCAGCACGTTATGGTGGGGACTCAGTGGAGACCGCCGGG
>splitRead2IdB ampliconId_12341234/3
TCGTGAGATGTTGGGTTAAGTCCCGCAACGTGCGCAACCCTTGCCATTAGTTGCTACGAAAGGGCACTCTAATGGGACTG
>splitRead2IdC ampliconId_23452345/3
TGTCGTGAGATGTTGGGTTAAGTCCCGCAAAGAGCGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
>splitRead2IdE ampliconId_67896789/3
TGTCGTGAGATGTTGGGTTAAAAAAAAAAAAAAACGCAACCCTTGTCATTAGTTGCTACGCAAGAGCACTCTAATGAGAC
>splitRead2IdG ampliconId_78907890/3
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
"""


rtax_expected_result_paired = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

rtax_expected_result_paired_with_fallback = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    'clusterIdD splitRead1IdD': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

rtax_expected_result_single = {
    'clusterIdA splitRead1IdA': ('k__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes', 1.0),
    'clusterIdB splitRead1IdB': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__', 1.0),
    'clusterIdC splitRead1IdC': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    'clusterIdD splitRead1IdD': ('k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__', 1.0),
    }

rtax_expected_result_single_lines = set([
    'clusterIdA splitRead1IdA\tk__Bacteria; p__Actinobacteria; c__Actinobacteria; o__Actinomycetales; f__Propionibacteriaceae; g__Propionibacterium; s__Propionibacterium acnes\t1.000\n',
    'clusterIdB splitRead1IdB\tk__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Diaphorobacter; s__\t1.000\n',
    'clusterIdC splitRead1IdC\tk__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__\t1.000\n',
    'clusterIdD splitRead1IdD\tk__Bacteria; p__Proteobacteria; c__Betaproteobacteria; o__Burkholderiales; f__; g__Aquabacterium; s__\t1.000\n',
    ])

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
