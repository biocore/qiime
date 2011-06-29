#!/usr/bin/env python
#unit tests for util.py

from os import rmdir
from os.path import split, abspath, dirname, exists, join
from glob import glob
from cogent import Sequence
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files
from cogent.cluster.procrustes import procrustes
from cogent.app.formatdb import build_blast_db_from_fasta_file
from cogent.util.misc import get_random_directory_name, remove_files
from StringIO import StringIO
from qiime.parse import fields_to_dict, parse_otu_table, parse_mapping_file
from qiime.util import (make_safe_f, FunctionWithParams, qiime_blast_seqs,
    extract_seqs_by_sample_id, get_qiime_project_dir, matrix_stats,
    raise_error_on_parallel_unavailable, merge_otu_tables,
    convert_OTU_table_relative_abundance, create_dir, handle_error_codes,
    summarize_pcoas, _compute_jn_pcoa_avg_ranges, _flip_vectors, IQR,
    idealfourths, isarray, matrix_IQR, degap_fasta_aln,
    write_degapped_fasta_to_file, compare_otu_maps, get_diff_for_otu_maps,
    merge_n_otu_tables, convert_otu_table_relative, write_seqs_to_fasta,
    split_fasta_on_sample_ids, split_fasta_on_sample_ids_to_dict,
    split_fasta_on_sample_ids_to_files, median_absolute_deviation,
    guess_even_sampling_depth, compute_days_since_epoch,
    get_interesting_mapping_fields,inflate_denoiser_output,
    flowgram_id_to_seq_id_map, count_seqs, count_seqs_from_file,
    count_seqs_in_filepaths,get_split_libraries_fastq_params_and_file_types,
    iseq_to_qseq_fields,get_top_fastq_two_lines,make_compatible_distance_matrices)

import numpy
from numpy import array, asarray


__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Daniel McDonald","Greg Caporaso", 
               "Justin Kuczynski", "Jens Reeder", "Catherine Lozupone"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

class TopLevelTests(TestCase):
    """Tests of top-level module functions."""
    
    def setUp(self):
        self.otu_table_f1 = otu_table_fake1.split('\n')
        self.otu_table_f2 = otu_table_fake2.split('\n')
        self.otu_table_f3 = otu_table_fake3.split('\n')
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

    def tearDown(self):    
        remove_files(self.files_to_remove)
        for dir in  self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)
                
    def test_median_absolute_deviation(self):
        """ median_absolute_deviation returns MAD and median """
        data = [0,0,0,0,0,0]
        expected = (0,0)
        self.assertEqual(median_absolute_deviation(data),expected)
        data = [1,1,1,1,1,1]
        expected = (0,1)
        self.assertEqual(median_absolute_deviation(data),expected)
        data = [6,1,2,9,4,1,2]
        expected = (1,2)
        self.assertEqual(median_absolute_deviation(data),expected)
        data = [-6,-1,-2,-9,-4,-1,-2]
        expected = (1,-2)
        self.assertEqual(median_absolute_deviation(data),expected)
        
    def test_guess_even_sampling_depth(self):
        """ guess_even_sampling_depth functions as expected """
        data = [6,1,2,9,4,1,2]
        expected = 1 # MAD = 2.25; med - MAD = -0.25
        self.assertEqual(guess_even_sampling_depth(data),expected)
    
    def test_write_seqs_to_fasta(self):
        """ write_seqs_to_fasta functions as expected """
        output_fp = get_tmp_filename(
                     prefix="qiime_util_write_seqs_to_fasta_test",
                     suffix='.fasta')
        self.files_to_remove.append(output_fp)
        seqs = [('s1','ACCGGTTGG'),('s2','CCTTGG'),('S4 some comment string','A')]
        exp = ">s1\nACCGGTTGG\n>s2\nCCTTGG\n>S4 some comment string\nA\n"
        # works in write mode
        write_seqs_to_fasta(output_fp,seqs,'w')
        self.assertEqual(open(output_fp).read(),exp)
        # calling again in write mode overwrites original file
        write_seqs_to_fasta(output_fp,seqs,'w')
        self.assertEqual(open(output_fp).read(),exp)
        # works in append mode
        exp2 = exp + exp
        write_seqs_to_fasta(output_fp,seqs,'a')
        self.assertEqual(open(output_fp).read(),exp2)
        
    
    def test_split_fasta_on_sample_ids(self):
        """ split_fasta_on_sample_ids functions as expected 
        """
        actual = list(split_fasta_on_sample_ids(\
                      MinimalFastaParser(self.fasta1)))
        expected = [('Samp1','Samp1_42','ACCGGTT'),
                       ('s2_a','s2_a_50','GGGCCC'),
                       ('Samp1','Samp1_43 some comme_nt','AACCG'),
                       ('s3','s3_25','AAACCC')]
        self.assertEqual(actual,expected)
    
    def test_split_fasta_on_sample_ids_to_dict(self):
        """ split_fasta_on_sample_ids_to_dict functions as expected
        """
        actual = split_fasta_on_sample_ids_to_dict(\
                      MinimalFastaParser(self.fasta1))
        expected = {'Samp1':[('Samp1_42','ACCGGTT'),
                             ('Samp1_43 some comme_nt','AACCG')],
                    's2_a':[('s2_a_50','GGGCCC')],
                    's3':[('s3_25','AAACCC')]}
        self.assertEqual(actual,expected)

    def test_split_fasta_on_sample_ids_to_files(self):
        """ split_fasta_on_sample_ids_to_files functions as expected 
        """
        temp_output_dir = get_random_directory_name(output_dir='/tmp/')
        self.dirs_to_remove.append(temp_output_dir)
        
        split_fasta_on_sample_ids_to_files(
         MinimalFastaParser(self.fasta2),
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
        self.assertEqual(len(glob('%s/*' % temp_output_dir)),3)

                
    def test_convert_otu_table_relative(self):
        """should convert a parsed otu table into relative abundances"""
        otu_table = parse_otu_table(self.otu_table_f1)
        exp_counts = array([[1.0/6, 0],
                            [1.0/6, 0],
                            [4.0/6, 0]])
        rel_otu_table = convert_otu_table_relative(otu_table)
        self.assertEqual(rel_otu_table[0], otu_table[0])
        self.assertEqual(rel_otu_table[1], otu_table[1])
        self.assertEqual(rel_otu_table[2], exp_counts)
        self.assertEqual(rel_otu_table[3], otu_table[3])

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
        
    def test_merge_n_otu_tables_error(self):
        """merge_n_otu_tables throws error on overlapping sample IDs"""
        # error on overlapping sample ids
        self.assertRaises(AssertionError,merge_n_otu_tables,
                          [iter(self.otu_table_f1),
                           iter(self.otu_table_f2),
                           iter(self.otu_table_f3),
                           iter(self.otu_table_f1)])
    
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
    
    def test_merge_otu_tables_passed_in_tables(self):
        """merge_otu_tables functions with otu_table tuples"""
        otu_table1 = parse_otu_table(self.otu_table_f1)
        otu_table2 = parse_otu_table(self.otu_table_f2)
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
        exp = (exp_sample_ids,exp_otu_ids,exp_otu_table,exp_lineages)
        obs = merge_otu_tables(otu_table1,otu_table2)
        self.assertEqual(obs,exp)

    def test_merge_n_otu_tables(self):
        """merge_n_otu_tables functions as expected"""
        otu_table_f1 = iter(self.otu_table_f1)
        otu_table_f2 = iter(self.otu_table_f2)
        otu_table_f3 = iter(self.otu_table_f3)
        exp_sample_ids = ['S1','S2','S3','S4','S5','samp7']
        exp_otu_ids = ['0','1','2','3','4','6']
        exp_otu_table = array([[1,0,1,0,1,0],\
                           [1,0,0,0,0,0],\
                           [4,0,1,0,1,0],\
                           [0,0,2,0,1,0],\
                           [0,0,1,0,9,0],\
                           [0,0,1,25,42,1]])
        exp_lineages = [['Root','Bacteria'],\
                    ['Root','Bacteria','Verrucomicrobia'],
                    ['Root','Bacteria'],\
                    ['Root','Bacteria','Acidobacteria'],\
                    ['Root','Bacteria','Bacteroidetes'],\
                    ['Root','Archaea']]
        actual = merge_n_otu_tables([otu_table_f1,otu_table_f2,otu_table_f3])
        self.assertEqual(actual[0],exp_sample_ids)
        self.assertEqual(actual[1],exp_otu_ids)
        self.assertEqual(actual[2],exp_otu_table)
        self.assertEqual(actual[3],exp_lineages)

    def test_merge_n_otu_tables_no_tax(self):
        """merge_n_otu_tables functions as expected with no taxonomy"""
        otu_table_f1 = iter(self.otu_table_f1_no_tax)
        otu_table_f2 = iter(self.otu_table_f2_no_tax)
        otu_table_f3 = iter(self.otu_table_f3_no_tax)
        exp_sample_ids = ['S1','S2','S3','S4','S5','samp7']
        exp_otu_ids = ['0','1','2','3','4','6']
        exp_otu_table = array([[1,0,1,0,1,0],\
                           [1,0,0,0,0,0],\
                           [4,0,1,0,1,0],\
                           [0,0,2,0,1,0],\
                           [0,0,1,0,9,0],\
                           [0,0,1,25,42,1]])
        exp_lineages = None
        actual = merge_n_otu_tables([otu_table_f1,otu_table_f2,otu_table_f3])
        self.assertEqual(actual[0],exp_sample_ids)
        self.assertEqual(actual[1],exp_otu_ids)
        self.assertEqual(actual[2],exp_otu_table)
        self.assertEqual(actual[3],exp_lineages)
        
    def test_merge_n_otu_tables_error_on_mixed_tax(self):
        """merge_n_otu_tables fails with some tax, some no tax"""
        otu_table_f1 = iter(self.otu_table_f1)
        otu_table_f2 = iter(self.otu_table_f2_no_tax)
        otu_table_f3 = iter(self.otu_table_f3_no_tax)
        self.assertRaises(ValueError,
         merge_n_otu_tables,[otu_table_f1,otu_table_f2,otu_table_f3])
         
        otu_table_f1 = iter(self.otu_table_f1)
        otu_table_f2 = iter(self.otu_table_f2)
        otu_table_f3 = iter(self.otu_table_f3_no_tax)
        self.assertRaises(ValueError,
         merge_n_otu_tables,[otu_table_f1,otu_table_f2,otu_table_f3])
         
        otu_table_f1 = iter(self.otu_table_f1_no_tax)
        otu_table_f2 = iter(self.otu_table_f2)
        otu_table_f3 = iter(self.otu_table_f3)
        self.assertRaises(ValueError,
         merge_n_otu_tables,[otu_table_f1,otu_table_f2,otu_table_f3])
        
    def test_compute_days_since_epoch(self):
        """compute_days_since_epoch functions as expected """
        self.assertEqual(compute_days_since_epoch(29,10,2002),11989)
        # can pass keyword arguments
        self.assertEqual(compute_days_since_epoch(year=2002,month=10,day=29),11989)
        self.assertEqual(compute_days_since_epoch(day=27,month=3,year=2009),14330)
        self.assertEqual(compute_days_since_epoch(day="27",month="3",year="2009"),14330)
    
    def test_get_interesting_mapping_fields(self):
        """get_interesting_mapping_fields returns expected fields """
        # all columns are completely unique
        d = parse_mapping_file(self.mapping_f1)
        actual = get_interesting_mapping_fields(d[0],d[1])
        expected = []
        self.assertEqual(actual,expected)
        
        # all columns are completely identical
        d = parse_mapping_file(self.mapping_f2)
        actual = get_interesting_mapping_fields(d[0],d[1])
        expected = []
        self.assertEqual(actual,expected)
        
        # some columns retained
        d = parse_mapping_file(self.mapping_f3)
        actual = get_interesting_mapping_fields(d[0],d[1])
        expected = ['Something','days_since_epoch']
        self.assertEqual(actual,expected)
        
    def test_inflate_denoiser_output(self):
        """ inflate_denoiser_output expands denoiser results as expected """
        actual = list(inflate_denoiser_output(
         MinimalFastaParser(self.centroid_seqs1),
         MinimalFastaParser(self.singleton_seqs1),
         self.denoiser_mapping1,
         MinimalFastaParser(self.raw_seqs1)))
        expected = [("S1_0 FXX111 some comments","TTTT"),
                    ("S1_2 FXX113 some other comments","TTTT"),
                    ("S3_7 FXX117","TTTT"),
                    ("S2_1 FXX112 some comments","TATT"),
                    ("S3_5 FXX114","TATT"),
                    ("S3_6 FXX115","TTGA"),
                    ("S3_6 FXX116","TAGA")]
        self.assertEqual(actual,expected)
    
    def test_flowgram_id_to_seq_id_map(self):
        """ flowgram_id_to_seq_id_map functions as expected """
        actual = flowgram_id_to_seq_id_map(MinimalFastaParser(self.raw_seqs1))
        expected = {'FXX111':'S1_0 FXX111 some comments',
                    'FXX112':'S2_1 FXX112 some comments',
                    'FXX113':'S1_2 FXX113 some other comments',
                    'FXX114':'S3_5 FXX114',
                    'FXX115':'S3_6 FXX115',
                    'FXX116':'S3_6 FXX116',
                    'FXX117':'S3_7 FXX117'}
        self.assertEqual(actual,expected)


    def test_count_seqs_from_file(self):
        """ count_seqs: functions as expected with varied data
        """
        f1 = ['>seq1','AACCTT','ACTGGT',
              '>seq2','CCAATT',
              '>seq3','CCC---GG']
        f2 = ['> s42','ABCDEFG',
              '>s33','A',
              '> 4>','AA>',
              '>blah','AA']
        self.assertFloatEqual(count_seqs_from_file(f1),(3,8.666,2.4944),0.001)
        self.assertFloatEqual(count_seqs_from_file(f2),(4,3.25,2.2776),0.001)
        self.assertEqual(count_seqs_from_file([]),(0,None,None))
        
    def test_count_seqs(self):
        """ count_seqs functions as expected with fake seq_counter
        """
        def seq_counter(filepath):
            # Fake sequence counter to test count_seqs without
            # having to write files to disk (note don't need to
            # test actual sequence counters here as they're tested
            # elsewhere)
            if filepath.startswith('fake'):
                raise IOError
            else:
                return len(filepath), 0, 0
                
        in_fps = ['1.fasta','fake1.fasta','fake.fasta','2.fa']
        expected = [((7,0,0),'1.fasta'),
                    ((4,0,0),'2.fa')],\
                    11, ['fake1.fasta','fake.fasta']
        self.assertEqual(count_seqs_in_filepaths(\
         in_fps,seq_counter),expected)
        
        in_fps = ['fake1.fasta','fake.fasta']
        expected = [], 0, ['fake1.fasta','fake.fasta']
        self.assertEqual(count_seqs_in_filepaths(\
         in_fps,seq_counter),expected)
                
        in_fps = ['1.fasta','2.fa','12.txt']
        expected = [((7,0,0),'1.fasta'),
                    ((4,0,0),'2.fa'),
                    ((6,0,0),'12.txt')], 17, []
        self.assertEqual(count_seqs_in_filepaths(\
         in_fps,seq_counter),expected)

    def test_get_split_libraries_fastq_params_and_file_types_reverse(self):
        """get_split_libraries_fastq_params_and_file_types using reverse 
           barcodes computes correct values"""
    
        temp_output_dir = get_random_directory_name(output_dir='/tmp/')
        self.dirs_to_remove.append(temp_output_dir)

        #generate the fastq mapping file
        map_fpath=join(temp_output_dir,'map.txt')
        map_fopen=open(map_fpath,'w')
        map_fopen.write('\n'.join(fastq_mapping_rev))
        map_fopen.close()
        self.files_to_remove.append(map_fpath)

        fastq_files=[]
        #generate fastq seqs file
        seq_fpath=join(temp_output_dir,'seqs.fastq')
        seqs_fopen=open(seq_fpath,'w')
        seqs_fopen.write('\n'.join(fastq_seqs))
        seqs_fopen.close()

        fastq_files.append(seq_fpath)
        self.files_to_remove.append(seq_fpath)

        #generate fastq seqs file
        barcode_fpath=join(temp_output_dir,'barcodes.fastq')
        barcode_fopen=open(barcode_fpath,'w')
        barcode_fopen.write('\n'.join(fastq_barcodes))
        barcode_fopen.close()

        fastq_files.append(barcode_fpath)
        self.files_to_remove.append(barcode_fpath)

        exp='-i %s -b %s --rev_comp_barcode' % (seq_fpath,barcode_fpath)

        obs=get_split_libraries_fastq_params_and_file_types(fastq_files,
                                                           map_fpath)

        self.assertEqual(obs,exp)

    def test_get_split_libraries_fastq_params_and_file_types_forward(self):
        """get_split_libraries_fastq_params_and_file_types using forward
           barcodes computes correct values"""
    
        temp_output_dir = get_random_directory_name(output_dir='/tmp/')
        self.dirs_to_remove.append(temp_output_dir)
    
        #generate the fastq mapping file
        map_fpath=join(temp_output_dir,'map.txt')
        map_fopen=open(map_fpath,'w')
        map_fopen.write('\n'.join(fastq_mapping_fwd))
        map_fopen.close()
        self.files_to_remove.append(map_fpath)
        
        fastq_files=[]
        #generate fastq seqs file
        seq_fpath=join(temp_output_dir,'seqs.fastq')
        seqs_fopen=open(seq_fpath,'w')
        seqs_fopen.write('\n'.join(fastq_seqs))
        seqs_fopen.close()
        
        fastq_files.append(seq_fpath)
        self.files_to_remove.append(seq_fpath)
        
        #generate fastq seqs file
        barcode_fpath=join(temp_output_dir,'barcodes.fastq')
        barcode_fopen=open(barcode_fpath,'w')
        barcode_fopen.write('\n'.join(fastq_barcodes))
        barcode_fopen.close()
        
        fastq_files.append(barcode_fpath)
        self.files_to_remove.append(barcode_fpath)
        
        exp='-i %s -b %s ' % (seq_fpath,barcode_fpath)
    
        obs=get_split_libraries_fastq_params_and_file_types(fastq_files,
                                                            map_fpath)
    
        self.assertEqual(obs,exp)

    def test_get_top_fastq_two_lines(self):
        """ get_top_fastq_two_lines: this function gets the first 4 lines of 
            the open fastq file
        """
        
        temp_output_dir = get_random_directory_name(output_dir='/tmp/')
        self.dirs_to_remove.append(temp_output_dir)

        fastq_files=[]
        #generate fastq seqs file
        seq_fpath=join(temp_output_dir,'seqs.fastq')
        seqs_fopen=open(seq_fpath,'w')
        seqs_fopen.write('\n'.join(fastq_seqs))
        seqs_fopen.close()

        fastq_files.append(seq_fpath)
        self.files_to_remove.append(seq_fpath)

        #generate fastq seqs file
        barcode_fpath=join(temp_output_dir,'barcodes.fastq')
        barcode_fopen=open(barcode_fpath,'w')
        barcode_fopen.write('\n'.join(fastq_barcodes))
        barcode_fopen.close()

        fastq_files.append(barcode_fpath)
        self.files_to_remove.append(barcode_fpath)
        exp=[('@HWUSI-EAS552R_0357:8:1:10040:6364#0/1\n', 'GACGAGTCAGTC\n', 
              '+HWUSI-EAS552R_0357:8:1:10040:6364#0/1\n', 'hhhhhhhhhhhh\n'),
             ('@HWUSI-EAS552R_0357:8:1:10040:6364#0/2\n', 
              'TACAGGGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTCTGCAGGTTGCTCACTAAGTCTTTTGTTAAATCTTCGGGCTTAACCCGAAACCTGCAAAAGAAACTAGTGCTCTCGAGTATGGTAGAGGTAAAGGGAATTTCCAG\n', 
              '+HWUSI-EAS552R_0357:8:1:10040:6364#0/2\n', 
              'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhfhhhhhhgghhhhhWhfcffehf]hhdhhhhhgcghhhchhhhfhcfhhgggdfhgdcffadccfdcccca]^b``ccfdd_caccWbb[b_dfdcdeaec`^`^_daba_b_WdY^`\n')]
        
        #iterate of dict and make sure the top 4 lines of each file are in the
        #expected list
        for i in fastq_files:
            obs=get_top_fastq_two_lines(open(i))
            self.assertTrue(obs in exp)
            
    def test_make_compatible_distance_matrices(self):
        """make_compatible_distance_matrices: functions as expected"""
        dm1 = (['A','B','C','D'],
               array([[0.0,2.3,3.3,4.3],
                      [2.9,0.0,5.3,6.3],
                      [3.9,5.9,0.0,4.3],
                      [4.9,8.9,9.9,0.0]]))
        
        dm2 = (['C','A','T'],
               array([[10.0,12.3,13.3],
                      [12.9,10.0,15.3],
                      [13.9,15.9,10.0]]))
        
        expected_dm1 = (['A','C'],
               array([[0.0,3.3],
                      [3.9,0.0]]))
        
        expected_dm2 = (['A','C'],
               array([[10.0,12.9],
                      [12.3,10.0]]))
        
        actual_dm1, actual_dm2 = make_compatible_distance_matrices(dm1,dm2)
        self.assertEqual(actual_dm1,expected_dm1)
        self.assertEqual(actual_dm2,expected_dm2)

    def test_make_compatible_distance_matrices_w_lookup(self):
        """make_compatible_distance_matrices: functions as expected with lookup"""
        dm1 = (['A','B','C','D'],
               array([[0.0,2.3,3.3,4.3],
                      [2.9,0.0,5.3,6.3],
                      [3.9,5.9,0.0,4.3],
                      [4.9,8.9,9.9,0.0]]))
        
        dm2 = (['C','A','T'],
               array([[10.0,12.3,13.3], 
                      [12.9,10.0,15.3], 
                      [13.9,15.9,10.0]]))
                      
        lookup = {'C':'C','A':'A','B':'B','T':'B','D':'D'}
        
        expected_dm1 = (['A','B','C'],
               array([[0.0,2.3,3.3],
                      [2.9,0.0,5.3],
                      [3.9,5.9,0.0]]))
        
        expected_dm2 = (['A','B','C'],
               array([[10.0,15.3,12.9],
                      [15.9,10.0,13.9],
                      [12.3,13.3,10.0]]))
        
        actual_dm1, actual_dm2 = make_compatible_distance_matrices(dm1,dm2,lookup)
        self.assertEqual(actual_dm1,expected_dm1)
        self.assertEqual(actual_dm2,expected_dm2)
        
        lookup = {'C':'C','B':'B','T':'B','D':'D'}
        self.assertRaises(KeyError,
         make_compatible_distance_matrices,dm1,dm2,lookup)


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
        """idealfourths: tests the ideal-fourths function which was imported from scipy
        at the following location (http://projects.scipy.org/scipy/browser/trunk/scipy/stats/tests/test_mmorestats.py?rev=4154)
        """
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
        obs1=asarray('Test')
        obs2 = 'Test'
        
        exp1 = True
        exp2 = False
        
        self.assertEqual(isarray(obs1),exp1)
        self.assertEqual(isarray(obs2),exp2)

    def test_degap_fasta_aln(self):
        """degap_fasta_align removes gps from fasta seqs."""

        test_aln = [("a","AAAAAAAAAGGGG"),
                    ("b","-A-A-G-G-A-G-C."),
                    ('c',"..-----------"),
                    ('d',"--------AAAAAAA"),
                    ('e',"")]

        expected_result = map(lambda (a, b): Sequence(name=a, seq=b),
                              [("a","AAAAAAAAAGGGG"),
                               ("b","AAGGAGC"),
                               ('c',""),
                               ('d',"AAAAAAA"),
                               ('e',"")])
        
        self.assertEqual(list(degap_fasta_aln(test_aln)), expected_result)

        self.assertEqual(list(degap_fasta_aln([])),[])

    def test_write_degapped_fasta_to_file(self):
     
        test_aln = [("a","AAAAAAAAAGGGG"),
                    ("b","-A-A-G-G-A-G-C"),
                    ('c',"---------------"),
                    ('d',"--------AAAAAAA"),
                    ('e',"")]

        expected_result =""">a
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
        observed = "".join(list(open(tmp_filename,"U")))
        self.assertEqual(observed, expected_result)


    def test_get_diff_for_otu_maps(self):
        """get_diff_for_otu_map return correct set difference"""

        #compare to self
        self.assertEqual(get_diff_for_otu_maps(otu_map1, otu_map1),
                         (set([]), set([])) )
        
        #compare to otu_map with one difference
        self.assertEqual(get_diff_for_otu_maps(otu_map1, otu_map2),
                         (set(['b']), set([])))
        
        #compare to empty
        self.assertEqual(get_diff_for_otu_maps(otu_map1, fields_to_dict("")),
                         (set(['a','b','c','d','e','f']), set([])))

    def test_compare_otu_maps(self):
        """compare_otu_maps computes correct values"""

        self.assertFloatEqual(compare_otu_maps(otu_map1, otu_map1), 0.0)
        self.assertFloatEqual(compare_otu_maps(otu_map1, otu_map3), 0.0)
        self.assertFloatEqual(compare_otu_maps(otu_map1, otu_map4), 0.33333333333)
        self.assertFloatEqual(compare_otu_maps(otu_map3, otu_map4), 0.33333333333)
        self.assertFloatEqual(compare_otu_maps(otu_map1, otu_map5), 1)

    def test_iseq_to_qseq_fields(self):
        """iseq_to_qseq_fields functions as expected"""
        i = "HWI-ST753_50:6:1101:15435:9071#0/1:ACCAGACGATGCTACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC:gggggggfggdegggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg"
        # barcode in sequence, barcode length = 12
        expected = (("HWI-ST753","50","6","1101","15435","9071","0","1"),
                    "TACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC","gggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg","ACCAGACGATGC","gggggggfggde")
        self.assertEqual(iseq_to_qseq_fields(i,barcode_in_header=False,barcode_length=12),
                         expected)
        # barcode in sequence, barcode length = 6
        expected = (("HWI-ST753","50","6","1101","15435","9071","0","1"),
                    "CGATGCTACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC","gfggdegggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg","ACCAGA","gggggg")
        self.assertEqual(iseq_to_qseq_fields(i,barcode_in_header=False,barcode_length=6),
                         expected)
                         
        # barcode in header, barcode length = 6
        i = "HWI-6X_9267:1:1:4:1699#ACCACCC/1:TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA:abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB"
        expected = (("HWI-6X","9267","1","1","4","1699","ACCACCC", "1"),
         "TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA","abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB","ACCACC","bbbbbb")
        self.assertEqual(iseq_to_qseq_fields(i,barcode_in_header=True,barcode_length=6),
                         expected)
        # barcode in header, barcode length = 3
        expected = (("HWI-6X","9267","1","1","4","1699","ACCACCC", "1"),
         "TACGGAGGGTGCGAGCGTTAATCGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCGAAAAAAAAAAAAAAAAAAAAAAA","abbbbbbbbbb`_`bbbbbb`bb^aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaDaabbBBBBBBBBBBBBBBBBBBB","ACC","bbb")
        self.assertEqual(iseq_to_qseq_fields(i,barcode_in_header=True,barcode_length=3),
                         expected)

 
otu_map1 = fields_to_dict("""1:\ta\tb\tc
2:\td
3:\te\tf
""".split("\n"))

#b missing
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

fastq_barcodes=["@HWUSI-EAS552R_0357:8:1:10040:6364#0/1",
"GACGAGTCAGTC",
"+HWUSI-EAS552R_0357:8:1:10040:6364#0/1",
"hhhhhhhhhhhh",
"@HWUSI-EAS552R_0357:8:1:10184:6365#0/1",
"GTCTGACAGTTG",
"+HWUSI-EAS552R_0357:8:1:10184:6365#0/1",
"hhhhhhhhhhhh"]
fastq_seqs=["@HWUSI-EAS552R_0357:8:1:10040:6364#0/2",
"TACAGGGGATGCAAGTGTTATCCGGAATTATTGGGCGTAAAGCGTCTGCAGGTTGCTCACTAAGTCTTTTGTTAAATCTTCGGGCTTAACCCGAAACCTGCAAAAGAAACTAGTGCTCTCGAGTATGGTAGAGGTAAAGGGAATTTCCAG",
"+HWUSI-EAS552R_0357:8:1:10040:6364#0/2",
"hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhfhhhhhhgghhhhhWhfcffehf]hhdhhhhhgcghhhchhhhfhcfhhgggdfhgdcffadccfdcccca]^b``ccfdd_caccWbb[b_dfdcdeaec`^`^_daba_b_WdY^`",
"@HWUSI-EAS552R_0357:8:1:10184:6365#0/2",
"TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGGCTCTTAAGTCGGAGGTGAAATCCCAAGGCTCAACCTTGGAACTGCCTTCGATACTGAGAGTCTTGAGTCCGGAAGAGGTAAGTGGAACTCCAAG",
"+HWUSI-EAS552R_0357:8:1:10184:6365#0/2",
"hfhhchhghhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhfghghhhggfhhfghfeghfggfdhdfdfffacbfcddgcfccddbcddbccada_aadWaaaacddccccdacaaa_acbc]c`aa[a\\a_a^V\T_^^^^X^R_BBBB"]

fastq_mapping_rev=["#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
"sample1\tGACTGACTCGTC\tCCGGACTACHVGGGTWTCTAAT\tsample1",
"sample2\tCAACTGTCAGAC\tCCGGACTACHVGGGTWTCTAAT\tsample2"]

fastq_mapping_fwd=["#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription",
"sample1\tGACGAGTCAGTC\tCCGGACTACHVGGGTWTCTAAT\tsample1",
"sample2\tGTCTGACAGTTG\tCCGGACTACHVGGGTWTCTAAT\tsample2"]


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
