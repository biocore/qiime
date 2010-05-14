#!/usr/bin/env python
#unit tests for util.py

from os import rmdir
from os.path import split, abspath, dirname, exists, join
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.util import make_safe_f, FunctionWithParams, qiime_blast_seqs,\
    extract_seqs_by_sample_id, get_qiime_project_dir, matrix_stats,\
    raise_error_on_parallel_unavailable, merge_otu_tables,\
    convert_OTU_table_relative_abundance, create_dir, handle_error_codes,\
    summarize_pcoas, _compute_jn_pcoa_avg_ranges, _flip_vectors, IQR, \
    idealfourths, isarray, matrix_IQR, sort_fasta_by_abundance
from cogent.app.formatdb import build_blast_db_from_fasta_file
from cogent.util.misc import get_random_directory_name
import numpy
from numpy import array, asarray
from cogent.cluster.procrustes import procrustes

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Daniel McDonald","Greg Caporaso", 
               "Justin Kuczynski", "Jens Reeder", "Catherine Lozupone"] 
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

class TopLevelTests(TestCase):
    """Tests of top-level module functions."""
    
    def setUp(self):
        self.otu_table_f1 = otu_table_fake1.split('\n')
        self.otu_table_f2 = otu_table_fake2.split('\n')
        self.dirs_to_remove = []
    
    def tearDown(self):
        for dir in  self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)

    def test_make_safe_f(self):
        """make_safe_f should return version of f that ignores extra kwargs."""
        def f(x,y): return x*y
        self.assertEqual(f(3,4), 12)
        g = make_safe_f(f, ['x','y'])
        self.assertEqual(g(3,4), 12)
        self.assertEqual(g(x=3,y=4,z=10,xxx=11), 12)
        
    def test_extract_seqs_by_sample_id(self):
        """extract_seqs_by_sample_id: functions as expected """
        
        seqs = [('Samp1_109','ACGG'),\
                ('Samp1_110','CCGG'),\
                ('samp1_109','GCGG'),\
                ('S2','AA'),\
                ('S3','CC'),\
                ('S4','GG'),\
                ('S44','TT'),\
                ('S4','TAAT')]
        sample_ids = ['Samp1','S44']
        expected = [('Samp1_109','ACGG'),\
                    ('Samp1_110','CCGG'),\
                    ('S44','TT')]
        actual = list(extract_seqs_by_sample_id(seqs,sample_ids))
        self.assertEqual(actual,expected)
        
        #negated
        expected_neg = [('samp1_109','GCGG'),\
                ('S2','AA'),\
                ('S3','CC'),\
                ('S4','GG'),\
                ('S4','TAAT')]
        actual = list(extract_seqs_by_sample_id(seqs,sample_ids,negate=True))
        self.assertEqual(actual,expected_neg)
        
        # OK if user passes dict of sample ids
        sample_ids = {'samp1':25}
        expected = [('samp1_109','GCGG')]
        actual = list(extract_seqs_by_sample_id(seqs,sample_ids))
        self.assertEqual(actual,expected)
        
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
        self.assertEqual(actual,expected)
    
    def test_matrix_stats1(self):
        """ matrix_stats should match mean, median, stdev calc'd by hand"""
        headers_list = [['a','c','b'],['a','c','b']]
        d1 = numpy.array([  [ 0,.2,.9],
                            [.2,0,.8],
                            [.9,.8,0]],'float')
        d2 = numpy.array([  [ 0,.3,1.1],
                            [.3,0,.8],
                            [1.1,.8,0]],'float')
        distmats_list = [d1,d2]
        
        exp_mean = numpy.array([    [ 0,.25,1.0],
                                    [.25,0,.8],
                                    [1.0,.8,0]],'float')
        exp_median = numpy.array([  [ 0,.25,1.0],
                                    [.25,0,.8],
                                    [1.0,.8,0]],'float')
        exp_std = numpy.array([     [ 0,.05,.1],
                                    [.05,0,0],
                                    [.1,0,0]],'float')
        results = matrix_stats(headers_list, distmats_list)
        self.assertFloatEqual(results[1:], [exp_mean,exp_median,exp_std])
        self.assertEqual(results[0],['a','c','b'])
        
    def test_matrix_stats2(self):
        """ matrix_stats should raise valerr if headers are in different orders
        """
        headers_list = [['a','c','b'],['b','c','a']]
        d1 = numpy.array([  [ 0,.2,.9],
                            [.2,0,.8],
                            [.9,.8,0]],'float')
        d2 = numpy.array([  [ 0,.3,1.1],
                            [.3,0,.8],
                            [1.1,.8,0]],'float')
        distmats_list = [d1,d2]

        exp_mean = numpy.array([    [ 0,.25,1.0],
                                    [.25,0,.8],
                                    [1.0,.8,0]],'float')
        exp_median = numpy.array([  [ 0,.25,1.0],
                                    [.25,0,.8],
                                    [1.0,.8,0]],'float')
        exp_std = numpy.array([     [ 0,.05,.1],
                                    [.05,0,0],
                                    [.1,0,0]],'float')
        self.assertRaises(ValueError, matrix_stats, headers_list, distmats_list)
        
    def test_raise_error_on_parallel_unavailable(self):
        """raise_error_on_parallel_unavailable functions as expected """
        self.assertRaises(RuntimeError,raise_error_on_parallel_unavailable,{})
        self.assertRaises(RuntimeError,raise_error_on_parallel_unavailable,\
         {'jobs_to_start':'1'})
        raise_error_on_parallel_unavailable({'jobs_to_start':'2'})
        raise_error_on_parallel_unavailable({'jobs_to_start':'24'})
        
    def test_support_files_available(self):
        """support_files are available """
        # check that the qiime/support_files directory exists
        support_files_dir = \
         join(get_qiime_project_dir(),'qiime','support_files')
        self.assertTrue(exists(support_files_dir))
        
        # check that a file in qiime/support_files exists
        default_qiime_config_fp = join(support_files_dir,'qiime_config')
        self.assertTrue(exists(default_qiime_config_fp))
        
    def test_merge_otu_tables_error(self):
        """merge_otu_tables throws error on overlapping sample IDs"""
        # error on overlapping sample ids
        self.assertRaises(AssertionError,\
         merge_otu_tables,iter(self.otu_table_f1),iter(self.otu_table_f1))
    
    def test_merge_otu_tables(self):
        """merge_otu_tables functions as expected"""
        otu_table_f1 = iter(self.otu_table_f1)
        otu_table_f2 = iter(self.otu_table_f2)
        exp_sample_ids = ['S1','S2','S3','S4','S5']
        exp_otu_ids = ['0','1','2','3','4','6']
        exp_otu_table = array([[1,0,1,0,1],\
                           [1,0,0,0,0],\
                           [4,0,1,0,1],\
                           [0,0,2,0,1],\
                           [0,0,1,0,9],\
                           [0,0,1,25,42]])
        exp_lineages = [['Root','Bacteria'],\
                    ['Root','Bacteria','Verrucomicrobia'],
                    ['Root','Bacteria'],\
                    ['Root','Bacteria','Acidobacteria'],\
                    ['Root','Bacteria','Bacteroidetes'],\
                    ['Root','Archaea']]
        actual = merge_otu_tables(otu_table_f1,otu_table_f2)
        self.assertEqual(actual[0],exp_sample_ids)
        self.assertEqual(actual[1],exp_otu_ids)
        self.assertEqual(actual[2],exp_otu_table)
        self.assertEqual(actual[3],exp_lineages)

    def test_create_dir(self):
        """create_dir creates dir and fails meaningful."""

        tmp_dir_path = get_random_directory_name()
        tmp_dir_path2 = get_random_directory_name(suppress_mkdir=True)
        tmp_dir_path3 = get_random_directory_name(suppress_mkdir=True)

        self.dirs_to_remove.append(tmp_dir_path)
        self.dirs_to_remove.append(tmp_dir_path2)
        self.dirs_to_remove.append(tmp_dir_path3)

        # create on existing dir raises OSError if fail_on_exist=True
        self.assertRaises(OSError, create_dir, tmp_dir_path,
                          fail_on_exist=True)
        self.assertEquals(create_dir(tmp_dir_path,
                                     fail_on_exist=True,
                                     handle_errors_externally=True), 1)

        # return should be 1 if dir exist and fail_on_exist=False 
        self.assertEqual(create_dir(tmp_dir_path, fail_on_exist=False), 1)

        # if dir not there make it and return always 0
        self.assertEqual(create_dir(tmp_dir_path2), 0)
        self.assertEqual(create_dir(tmp_dir_path3, fail_on_exist=True), 0)

    def test_handle_error_codes(self):
        """handle_error_codes raises the right error."""

        self.assertRaises(OSError, handle_error_codes, "test", False,1)
        self.assertEqual(handle_error_codes("test", True, 1), 1)
        self.assertEqual(handle_error_codes("test", False, 0), 0)
        self.assertEqual(handle_error_codes("test"), 0)
        
    def test_sort_fasta_by_abundance(self):
        """sort_fasta_by_abundance functions as expected"""
        class FakeOutF(object):
            def __init__(self):
                self.s = ""
            def write(self,s):
                self.s += s
        
        actual = FakeOutF()
        expected = ""
        sort_fasta_by_abundance([],actual)
        self.assertEqual(actual.s,expected)
        
        # no sorting necessary
        actual = FakeOutF()
        expected1 = "\n".join(['>s1','ACCGT','>s2 comment','ATTGC',''])
        expected2 = "\n".join(['>s2 comment','ATTGC','>s1','ACCGT',''])
        sort_fasta_by_abundance(['>s1','ACCGT','>s2 comment','ATTGC'],actual)
        # order is unimportant here
        self.assertTrue(actual.s == expected1 or actual.s == expected2)
        
        # sorting necessary
        actual = FakeOutF()
        inseqs = ['>s1','ACCGT',
                   '>s2 comment','ATTGC',
                   '>s3 blah','ATTGC']
        expected = "\n".join(['>s2 comment','ATTGC',
                              '>s3 blah','ATTGC',
                              '>s1','ACCGT',''])
        sort_fasta_by_abundance(inseqs,actual)
        self.assertEqual(actual.s,expected)
        
        # sorting necessary, but skipped due to low sampling
        # actual = FakeOutF()
        # inseqs = ['>s1','ACCGT',
        #            '>s2 comment','ATTGC'
        #            '>s3 blah','ATTGC']
        # expected = "\n".join(['>s1','ACCGT',
        #                       '>s2 comment','ATTGC'
        #                       '>s3 blah','ATTGC',''])
        # sort_fasta_by_abundance(inseqs,actual,sample=2)
        # self.assertEqual(actual.s,expected)
        



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

    def test_init(self):
        """FunctionWithParams __init__ should be successful."""
        x = self.FWP({'x':3})
        self.assertEqual(x.Name, 'FWP')
        self.assertEqual(x.Params, {'x':3})

    def test_str(self):
        """FunctionWithParams __str__ should produce expected string"""
        x = self.FWP({'x':3})
        lines = ['FWP parameters:','t:True','Application:None',\
                'Algorithm:None','Citation:None','x:3']
        self.assertEqual(str(x), '\n'.join(lines))

    def test_call(self):
        """FunctionWithParams __str__ should produce expected string"""
        x = self.FWP({'x':3})
        self.assertEqual(x(), 3)

    def test_formatResult(self):
        """FunctionWithParams formatResult should produce expected format"""
        x = self.FWP({'x':3})
        self.assertEqual(x.formatResult(3), '3')

class BlastSeqsTests(TestCase):
    """ Tests of the qiime_blast_seqs function (will move to PyCogent eventually)
    """

    def setUp(self):
        """ 
        """
        self.refseqs1 = refseqs1.split('\n')
        self.inseqs1 = inseqs1.split('\n')
        self.blast_db, db_files_to_remove =\
          build_blast_db_from_fasta_file(self.refseqs1,output_dir='/tmp/')
        self.files_to_remove = db_files_to_remove
        
        self.refseqs1_fp = get_tmp_filename(\
         tmp_dir='/tmp/', prefix="BLAST_temp_db_", suffix=".fasta")    
        fasta_f = open(self.refseqs1_fp,'w')
        fasta_f.write(refseqs1)
        fasta_f.close()
        
        self.files_to_remove = db_files_to_remove + [self.refseqs1_fp]
          
    def tearDown(self):
        remove_files(self.files_to_remove)
        
    def test_w_refseqs_file(self):
        """qiime_blast_seqs functions with refseqs file 
        """
        inseqs = MinimalFastaParser(self.inseqs1)
        actual = qiime_blast_seqs(inseqs,refseqs=self.refseqs1)
        self.assertEqual(len(actual),5)
        
        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'],'s2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'],'s1')
        
    def test_w_refseqs_fp(self):
        """qiime_blast_seqs functions refseqs_fp
        """
        inseqs = MinimalFastaParser(self.inseqs1)
        actual = qiime_blast_seqs(inseqs,refseqs_fp=self.refseqs1_fp)
        self.assertEqual(len(actual),5)
        
        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'],'s2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'],'s1')
    
    def test_w_preexising_blastdb(self):
        """qiime_blast_seqs functions with pre-existing blast_db
        """        
        # pre-existing blast db
        inseqs = MinimalFastaParser(self.inseqs1)
        actual = qiime_blast_seqs(inseqs,blast_db=self.blast_db)
        self.assertEqual(len(actual),5)
        
        # couple of sanity checks against command line blast
        self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'],'s2')
        self.assertEqual(actual['s105'][0][2]['SUBJECT ID'],'s1')
    
    def test_w_alt_seqs_per_blast_run(self):
        """qiime_blast_seqs: functions with alt seqs_per_blast_run
        """
        for i in range(1,20):
            inseqs = MinimalFastaParser(self.inseqs1)
            actual = qiime_blast_seqs(\
             inseqs,blast_db=self.blast_db,seqs_per_blast_run=i)
            self.assertEqual(len(actual),5)
    
            # couple of sanity checks against command line blast
            self.assertEqual(actual['s2_like_seq'][0][0]['SUBJECT ID'],'s2')
            self.assertEqual(actual['s105'][0][2]['SUBJECT ID'],'s1')
            
    def test_alt_blast_param(self):
        """qiime_blast_seqs: alt blast params give alt results"""
        # Fewer blast hits with stricter e-value
        inseqs = MinimalFastaParser(self.inseqs1)
        actual = qiime_blast_seqs(inseqs,blast_db=self.blast_db,params={'-e':'1e-4'})
        # NOTE: A BUG (?) IN THE BlastResult OR THE PARSER RESULTS IN AN EXTRA
        # BLAST RESULT IN THE DICT (actual HERE) WITH KEY ''
        
        self.assertTrue('s2_like_seq' in actual)
        self.assertFalse('s100' in actual)
        
    def test_error_on_bad_param_set(self):
        inseqs = MinimalFastaParser(self.inseqs1)
        # no blastdb or refseqs
        self.assertRaises(AssertionError,qiime_blast_seqs,inseqs)
        
    def test_convert_OTU_table_relative_abundance(self):
        """convert_OTU_table_relative_abundance works
        """
        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""".split('\n')
        result = convert_OTU_table_relative_abundance(otu_table)
        self.assertEqual(result, ['#Full OTU Counts', '#OTU ID\tsample1\tsample2\tsample3', '0\t0.0\t0.666666666667\t0.0', '1\t0.5\t0.0\t0.0', '2\t0.5\t0.333333333333\t1.0'])

        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage
0\t0\t2\t0\tBacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475
1\t1\t0\t0\tBacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369
2\t1\t1\t1\tBacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121""".split('\n')
        result = convert_OTU_table_relative_abundance(otu_table)
        self.assertEqual(result, ['#Full OTU Counts', '#OTU ID\tsample1\tsample2\tsample3\tConsensus Lineage', '0\t0.0\t0.666666666667\t0.0\tBacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; Unclassified; otu_475', '1\t0.5\t0.0\t0.0\tBacteria; Bacteroidetes; Bacteroidales; adhufec77-25; Barnesiella; Barnesiella_viscericola; otu_369', '2\t0.5\t0.333333333333\t1.0\tBacteria; Firmicutes; Clostridia; Clostridiales; Faecalibacterium; Unclassified; otu_1121'])

    def test_flip_vectors(self):
        """_flip_vectors makes a new PCA matrix with correct signs"""
        m_matrix = array([[1.0, 0.0, 1.0], [2.0, 4.0, 4.0]])
        jn_matrix = array([[1.2, 0.1, -1.2], [2.5, 4.0, -4.5]])
        new_matrix = _flip_vectors(jn_matrix, m_matrix)
        self.assertEqual(new_matrix, array([[1.2, 0.1, 1.2], [2.5, 4.0, 4.5]]))

    def test_compute_jn_pcoa_avg_ranges(self):
        """_compute_jn_pcoa_avg_ranges works
        """
        jn_flipped_matrices = [array([[2.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[3.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[4.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[5.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[6.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[7.0,4.0, -4.5],[-1.2,-0.1,1.2]]),\
                array([[1.0,4.0, -4.5],[-1.2,-0.1,1.2]])]
        avg_matrix, low_matrix, high_matrix = _compute_jn_pcoa_avg_ranges(\
                jn_flipped_matrices, 'ideal_fourths')
        self.assertFloatEqual(avg_matrix[(0,0)], 4.0)
        self.assertFloatEqual(avg_matrix[(0,2)], -4.5)
        self.assertFloatEqual(low_matrix[(0,0)], 2.16666667)
        self.assertFloatEqual(high_matrix[(0,0)], 5.83333333)

        avg_matrix, low_matrix, high_matrix = _compute_jn_pcoa_avg_ranges(\
                jn_flipped_matrices, 'sdev')
        x = array([m[0,0] for m in jn_flipped_matrices])
        self.assertEqual(x.mean(),avg_matrix[0,0])
        self.assertEqual(-x.std(ddof=1)/2,low_matrix[0,0])
        self.assertEqual(x.std(ddof=1)/2,high_matrix[0,0])
        
    def test_summarize_pcoas(self):
        """summarize_pcoas works
        """
        master_pcoa = [['1', '2', '3'], \
            array([[-1.0, 0.0, 1.0], [2.0, 4.0, -4.0]]), \
            array([.76, .24])]
        jn1 = [['1', '2', '3'], \
            array([[1.2, 0.1, -1.2],[-2.5, -4.0, 4.5]]), \
            array([0.80, .20])]
        jn2 = [['1', '2', '3'], \
            array([[-1.4, 0.05, 1.3],[2.6, 4.1, -4.7]]), \
            array([0.76, .24])]
        jn3 = [['1', '2', '3'], \
            array([[-1.5, 0.05, 1.6],[2.4, 4.0, -4.8]]), \
            array([0.84, .16])]
        jn4 = [['1', '2', '3'], \
            array([[-1.5, 0.05, 1.6],[2.4, 4.0, -4.8]]), \
            array([0.84, .16])]
        support_pcoas = [jn1, jn2, jn3, jn4]
        #test with the ideal_fourths option
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, 'ideal_fourths',
                            apply_procrustes=False)
        self.assertEqual(m_names, ['1', '2', '3'])
        self.assertFloatEqual(matrix_average[(0,0)], -1.4)
        self.assertFloatEqual(matrix_average[(0,1)], 0.0125)
        self.assertFloatEqual(matrix_low[(0,0)], -1.5)
        self.assertFloatEqual(matrix_high[(0,0)], -1.28333333)
        self.assertFloatEqual(matrix_low[(0,1)], -0.0375)
        self.assertFloatEqual(matrix_high[(0,1)], 0.05)
        self.assertFloatEqual(eigval_average[0], 0.81)
        self.assertFloatEqual(eigval_average[1], 0.19)
        #test with the IQR option
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, method='IQR',
                            apply_procrustes=False)
        self.assertFloatEqual(matrix_low[(0,0)], -1.5)
        self.assertFloatEqual(matrix_high[(0,0)], -1.3)

        #test with procrustes option followed by sdev
        m, m1, msq = procrustes(master_pcoa[1],jn1[1])
        m, m2, msq = procrustes(master_pcoa[1],jn2[1])
        m, m3, msq = procrustes(master_pcoa[1],jn3[1])
        m, m4, msq = procrustes(master_pcoa[1],jn4[1])
        matrix_average, matrix_low, matrix_high, eigval_average, m_names = \
            summarize_pcoas(master_pcoa, support_pcoas, method='sdev',
                            apply_procrustes=True)

        x = array([m1[0,0],m2[0,0],m3[0,0],m4[0,0]])
        self.assertEqual(x.mean(),matrix_average[0,0])
        self.assertEqual(-x.std(ddof=1)/2,matrix_low[0,0])
        self.assertEqual(x.std(ddof=1)/2,matrix_high[0,0])

    def test_IQR(self):
        "IQR returns the interquartile range for list x"
        #works for odd with odd split
        x = [2,3,4,5,6,7,1]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 6)
        #works for even with odd split
        x = [1,2,3,4,5,6]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 5)
        #works for even with even split
        x = [1,2,3,4,5,6,7,8]
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2.5)
        self.assertEqual(maxv, 6.5)
        #works with array
        #works for odd with odd split
        x = array([2,3,4,5,6,7,1])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 6)
        #works for even with odd split
        x = array([1,2,3,4,5,6])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2)
        self.assertEqual(maxv, 5)
        #works for even with even split
        x = array([1,2,3,4,5,6,7,8])
        minv, maxv = IQR(x)
        self.assertEqual(minv, 2.5)
        self.assertEqual(maxv, 6.5)
        
    def test_matrix_IQR(self):
        """matrix_IQR calcs the IQR for each column in an array correctly
        """
        x = array([[1,2,3],[4,5,6],[7,8,9], [10,11,12]])
        min_vals, max_vals = matrix_IQR(x)
        self.assertEqual(min_vals, array([2.5,3.5,4.5]))
        self.assertEqual(max_vals, array([8.5,9.5,10.5]))

    def test_idealfourths(self):
        "idealfourths: tests the ideal-fourths function which was imported from scipy"
        test = numpy.arange(100)
        self.assertEqual(idealfourths(test),
                            [24.416666666666668, 74.583333333333343])
        test_2D = test.repeat(3).reshape(-1,3)
        
        self.assertFloatEqualRel(numpy.asarray(idealfourths(test_2D, axis=0)),\
                    numpy.array([[24.41666667, 24.41666667, 24.41666667], \
                                 [74.58333333, 74.58333333, 74.58333333]]))
        
        self.assertEqual(idealfourths(test_2D, axis=1),
                            test.repeat(2).reshape(-1,2))
        test = [0,0]
        _result = idealfourths(test)
        self.assertEqual(numpy.isnan(_result).all(),True)
        
    def test_isarray(self):
        "isarray: tests the isarray function"
        test1=asarray('Test')
        test2 = 'Test'
        
        exp1 = True
        exp2 = False
        
        self.assertEqual(isarray(test1),exp1)
        self.assertEqual(isarray(test2),exp2)
        
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


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
