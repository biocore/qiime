#!/usr/bin/env python
#unit tests for format.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","Jeremy Widmann","Jens Reeder", "Daniel McDonald"] 
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import json
from os import remove
from numpy import array, nan
from cogent.util.misc import remove_files
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import get_tmp_filename, get_qiime_library_version
from qiime.parse import fields_to_dict, parse_mapping_file
from qiime.format import (format_distance_matrix, format_otu_table,
    format_coords, build_prefs_string, format_matrix, format_map_file,
    format_histograms, write_Fasta_from_name_seq_pairs, 
    format_unifrac_sample_mapping,format_otu_map,write_otu_map, 
    format_summarize_taxa, write_summarize_taxa, 
    format_add_taxa_summary_mapping, write_add_taxa_summary_mapping,
    format_qiime_parameters, format_p_value_for_num_iters,
    format_mapping_file,illumina_data_to_fastq, format_biom_table,
    format_mapping_html_data)
from biom.parse import parse_biom_table

class TopLevelTests(TestCase):
    """Tests of top-level module functions."""
    
    def setUp(self):
        self.otu_map1 = [('0',['seq1','seq2','seq5']),
                         ('1',['seq3','seq4']),
                         ('2',['seq6','seq7','seq8'])]
        self.tmp_fp1 = get_tmp_filename(prefix='FormatTests_',suffix='.txt')
        self.files_to_remove = []

        self.taxa_summary = [[('a','b','c'),0,1,2],
                             [('d','e','f'),3,4,5]]
        self.taxa_header = ['Taxon','foo','bar','foobar']
       
        self.add_taxa_summary = {'s1':[1,2],'s2':[3,4]}
        self.add_taxa_header = ['sample_id','foo','bar']
        self.add_taxa_order = [('a','b','c'),('d','e','f')]
        self.add_taxa_mapping = [['s1','something1','something2'],
                                 ['s2','something3','something4'],
                                 ['s3','something5','something6']]
        self.biom1 = parse_biom_table(biom1.split('\n'))
        
        self.expected_formatted_html_no_errors_warnings =\
         expected_formatted_html_no_errors_warnings
        self.expected_formatted_html_errors =\
         expected_formatted_html_errors
        self.expected_formatted_html_warnings =\
         expected_formatted_html_warnings
        self.expected_formatted_html_data_nonloc_error =\
         expected_formatted_html_data_nonloc_error
        
    def tearDown(self):
        remove_files(self.files_to_remove)
        
    
    def test_format_biom_table(self):
        """ Formatting of BIOM table correctly includes "generated-by" information 
        """
        generated_by = "QIIME " + get_qiime_library_version()
        self.assertTrue(generated_by in format_biom_table(self.biom1))

    def test_format_mapping_file(self):
        """ format_mapping file should match expected result"""
        headers = ['SampleID','col1','col0','Description']
        samples =\
         [['bsample','v1_3','v0_3','d1'],['asample','aval','another','d2']]
        comments = ['this goes after headers','this too']
        self.assertEqual(format_mapping_file(headers,samples,comments),
         example_mapping_file)
        # need file or stringIO for roundtrip test
        # roundtrip = parse_mapping_file(format_mapping_file(headers,samples,comments))
        # self.assertEqual(roundtrip, [headers,samples,comments])

    def test_format_p_value_for_num_iters(self):
        """ format_p_value_for_num_iters functions as expected """
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,100),"0.12")
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,250),"0.12")
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,1000),"0.119")
        # test num_iters too low still returns a string (this can
        # be the last step of a long process, so we don't want to fail)
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,9),
         "Too few iters to compute p-value (num_iters=9)")
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,1),
         "Too few iters to compute p-value (num_iters=1)")
        self.assertEqual(\
         format_p_value_for_num_iters(0.119123123123,0),
         "Too few iters to compute p-value (num_iters=0)")

    def test_format_summarize_taxa(self):
        """format_summarize_taxa functions as expected"""
        exp = '\n'.join(['Taxon\tfoo\tbar\tfoobar',
                         'a;b;c\t0\t1\t2',
                         'd;e;f\t3\t4\t5\n'])
        obs = ''.join(list(format_summarize_taxa(self.taxa_summary, \
                                                 self.taxa_header)))
        self.assertEqual(obs, exp)
        
    def test_write_summarize_taxa(self):
        """write_summarize_taxa functions as expected"""
        write_summarize_taxa(self.taxa_summary, self.taxa_header, self.tmp_fp1)
        obs = open(self.tmp_fp1).read()
        exp = '\n'.join(['Taxon\tfoo\tbar\tfoobar',
                         'a;b;c\t0\t1\t2',
                         'd;e;f\t3\t4\t5\n'])
        self.assertEqual(obs,exp)
        self.files_to_remove.append(self.tmp_fp1)
        
        
    def test_write_summarize_taxa_transposed_output(self):
        """write_summarize_taxa_transposed_output functions as expected"""
        write_summarize_taxa(self.taxa_summary, self.taxa_header, self.tmp_fp1, transposed_output=True)
        obs = open(self.tmp_fp1).read()
        exp = '\n'.join(['SampleID\ta;b;c\td;e;f',
                         'foo\t0\t3\nbar\t1\t4',
                         'foobar\t2\t5\n'])
        self.assertEqual(obs,exp)
        self.files_to_remove.append(self.tmp_fp1)
        

    def test_format_add_taxa_summary_mapping(self):
        """format_add_taxa_summary_mapping functions as expected"""
        exp = '\n'.join(['#sample_id\tfoo\tbar\ta;b;c\td;e;f',
                         's1\tsomething1\tsomething2\t1\t2',
                         's2\tsomething3\tsomething4\t3\t4\n'])
        tmp = format_add_taxa_summary_mapping(self.add_taxa_summary,\
                                              self.add_taxa_order,\
                                              self.add_taxa_mapping, \
                                              self.add_taxa_header)
        obs = ''.join(list(tmp))
        self.assertEqual(obs,exp)

    def test_format_qiime_parameters(self):
        """format_qiime_parameters: returns lines in qiime_parameters format"""
        params = {'pick_otus':
                    {'similarity':'0.94','otu_picking_method':'cdhit'},
                 'assign_taxonomy':
                    {'use_rdp':None}}
        obs = format_qiime_parameters(params)
        exp = ["#QIIME parameters",
               "assign_taxonomy:use_rdp\tTrue",
               "pick_otus:otu_picking_method\tcdhit",
               "pick_otus:similarity\t0.94"]
        self.assertEqual(obs, exp)

    def test_write_add_taxa_summary_mapping(self):
        """write_add_taxa_summary_mapping functions as expected"""
        write_add_taxa_summary_mapping(self.add_taxa_summary,\
                                       self.add_taxa_order,\
                                       self.add_taxa_mapping,\
                                       self.add_taxa_header,\
                                       self.tmp_fp1)
        obs = open(self.tmp_fp1).read()
        exp = '\n'.join(['#sample_id\tfoo\tbar\ta;b;c\td;e;f',
                         's1\tsomething1\tsomething2\t1\t2',
                         's2\tsomething3\tsomething4\t3\t4\n'])
        self.assertEqual(obs, exp)
        self.files_to_remove.append(self.tmp_fp1)

    def test_format_otu_map(self):
        """format_otu_map functions as expected """
        actual = list(format_otu_map(self.otu_map1,''))
        actual.sort()
        expected = ['0\tseq1\tseq2\tseq5\n',
                    '1\tseq3\tseq4\n',
                    '2\tseq6\tseq7\tseq8\n']
        expected.sort()
        self.assertEqual(actual,expected)
        
    def test_write_otu_map(self):
        """write_otu_map functions as expected """
        write_otu_map(self.otu_map1,self.tmp_fp1)
        actual = fields_to_dict(open(self.tmp_fp1))
        self.files_to_remove.append(self.tmp_fp1)
        self.assertEqual(actual,dict(self.otu_map1))
        
        
    def test_write_otu_map_prefix(self):
        """write_otu_map functions as expected w otu prefix """
        write_otu_map(self.otu_map1,self.tmp_fp1,'my.otu.')
        actual = fields_to_dict(open(self.tmp_fp1))
        self.files_to_remove.append(self.tmp_fp1)
        
        exp = {'my.otu.0':['seq1','seq2','seq5'],
               'my.otu.1':['seq3','seq4'],
               'my.otu.2':['seq6','seq7','seq8']}
        self.assertEqual(actual,exp)
        
    
    def test_format_otu_map_prefix(self):
        """format_otu_map functions as expected w prefix"""
        actual = list(format_otu_map(self.otu_map1,'my.otu.'))
        actual.sort()
        expected = ['my.otu.0\tseq1\tseq2\tseq5\n',
                    'my.otu.1\tseq3\tseq4\n',
                    'my.otu.2\tseq6\tseq7\tseq8\n']
        expected.sort()
        self.assertEqual(actual,expected)
    
    def test_format_otu_map_error_on_bad_prefix(self):
        """format_otu_map functions as expected with bad prefix char"""
        self.assertRaises(ValueError,list,
                          format_otu_map(self.otu_map1,'my_otu_'))

    def test_format_distance_matrix(self):
        """format_distance_matrix should return tab-delimited dist mat"""
        a = array([[1,2,3],[4,5,6],[7,8,9]])
        labels = [11,22,33]
        res = format_distance_matrix(labels, a)
        self.assertEqual(res, 
            '\t11\t22\t33\n11\t1\t2\t3\n22\t4\t5\t6\n33\t7\t8\t9')
        self.assertRaises(ValueError, format_distance_matrix, labels[:2], a)

    def test_format_matrix(self):
        """format_matrix should return tab-delimited mat"""
        a = [[1,2,3], [4,5,6], [7,8,9]]
        row_labels = ['a','b','c']
        col_labels = [11,22,33]
        res = format_matrix(a, row_labels, col_labels)
        
        #test as list
        self.assertEqual(res, 
            '\t11\t22\t33\na\t1\t2\t3\nb\t4\t5\t6\nc\t7\t8\t9')
        self.assertRaises(ValueError, format_matrix, a, row_labels[:2], col_labels)
        self.assertRaises(ValueError, format_matrix, None, row_labels, col_labels)

        #tes as array
        a = array(a)
        self.assertEqual(res, 
                         '\t11\t22\t33\na\t1\t2\t3\nb\t4\t5\t6\nc\t7\t8\t9')
        self.assertRaises(ValueError, format_matrix, a, row_labels[:2], col_labels)
        self.assertRaises(ValueError, format_matrix, None, row_labels, col_labels)

    def assertEqualOtuTable(self,obs,exp):
        """ """
        obs = json.loads(obs)
        exp = json.loads(exp)
        for e in ['generated_by','date']:
            del obs[e]
            del exp[e]
        self.assertEqual(obs,exp)

    def test_format_otu_table(self):
        """format_otu_table should return biom-formatted string"""
        a = array([[1,2,3],
                   [4,5,2718281828459045]])
        samples = ['a','b','c']
        otus = [1,2]
        taxa = ['Bacteria','Archaea']
        res = format_otu_table(samples, otus, a)
        # confirm that parsing the res gives us a valid biom file with 
        # expected observation and sample ids
        t = parse_biom_table(res.split('\n'))
        self.assertEqual(t.ObservationIds,['1','2'])
        self.assertEqual(t.SampleIds,['a','b','c'])

    def test_format_coords(self):
        """format_coords should return tab-delimited table of coords"""
        a = array([[1,2,3],[4,5,6],[7,8,9]])
        header = list('abc')
        eigvals = [2,4,6]
        pct_var = [3,2,1]
        res = format_coords(header, a, eigvals, pct_var)
        self.assertEqual(res, "pc vector number\t1\t2\t3\na\t1\t2\t3\nb\t4\t5\t6\nc\t7\t8\t9\n\n\neigvals\t2\t4\t6\n% variation explained\t3\t2\t1")
    
    def test_build_prefs_string(self):
        """build_prefs_string should return a properly formatted prefs string.
        """
        #Try with correctly formatted color_by_string
        mapping_headers_to_use='First,Second'
        background_color='black'
        monte_carlo_dist=10
        otu_ids=['Root;Bacteria']
        headers=['First','Second']
        ball_size=2.5
        arrow_head_color='red'
        arrow_line_color='white'
        exp_string = \
        """{\n'background_color':'black',\n\n'sample_coloring':\n\t{\n\t\t'First':\n\t\t{\n\t\t\t'column':'First',\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t\t},\n\t\t'Second':\n\t\t{\n\t\t\t'column':'Second',\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t\t}\n\t},\n'MONTE_CARLO_GROUP_DISTANCES':\n\t{\n\t\t'First': 10,\n\t\t'Second': 10\n\t},\n'FIELDS':\n\t[\n\t\t'Second',\n\t\t'First'\n\t],\n'taxonomy_coloring':\n\t{\n\t\t'Level_1':\n\t\t{\n\t\t\t'column':'1',\n\t\t\t'colors':\n\t\t\t{\n\t\t\t\t'Root;Bacteria':('red0',(0,100,100))\n\t\t\t}\n\t\t}\n\t},\n'ball_scale':'2.500000',\n'arrow_line_color':'white',\n'arrow_head_color':'red'\n}"""
        obs_string = build_prefs_string(mapping_headers_to_use, \
                                    background_color, monte_carlo_dist, headers,
                                    otu_ids, ball_size, arrow_line_color,
                                    arrow_head_color)
        
        self.assertEqual(obs_string,exp_string)
        
    def test_format_map_file(self):
        """format_map_file should produce correct result"""
        
        desc_key = "Description"
        sample_id = "SampleID"
        headers = ['SampleID', 'a', 'Description', 'b']
        id_map = {'x':{'a':3,'b':4}, 'y':{'a':5,'b':6}}
        desc_map = {'x':'sample x','y':'sample y'}
        run_desc = 'run desc'
        self.assertEqual(format_map_file(headers, id_map, desc_key, sample_id,\
         desc_map, run_desc),
"""#SampleID\ta\tb\tDescription
#run desc
x\t3\t4\tsample x
y\t5\t6\tsample y""")

    def test_format_histograms(self):
        """format_histograms should print histograms correctly"""
        self.assertEqual(format_histograms(array([0,1,0,2,2,3]),
         array([2,1,0,2,0,0]), array([0,0,0,2,0,1]),
         array([100,110,120,130,140,150,160])),
            """# bins raw sequence lengths, length of sequences that pass quality filters before processing, and lengths of sequences that pass quality filters post processing.\nLength\tRaw\tBefore\tAfter\n100\t0\t2\t0\n110\t1\t1\t0\n120\t0\t0\t0\n130\t2\t2\t2\n140\t2\t0\t0\n150\t3\t0\t1""")

    def test_write_Fasta_from_name_seqs_pairs(self):
        """write_Fasta_from_name_seqs_pairs write proper FASTA string."""
        
        seqs = [('1',"AAA"),('2',"CCCCC"),('3',"GGGG")]

        #None fh raises Error
        self.assertRaises(ValueError, write_Fasta_from_name_seq_pairs,seqs,None)

        tmp_filename = get_tmp_filename(prefix="test_write_Fasta", suffix=".fna")
        fh = open(tmp_filename,"w")
        write_Fasta_from_name_seq_pairs(seqs,fh)
        fh.close()
        actual_seqs = list(MinimalFastaParser(open(tmp_filename,"U")))
        remove(tmp_filename)
        
        self.assertEqual(actual_seqs, seqs)
        
    def test_format_unifrac_sample_mapping(self):
        """format sample mapping works
        """
        a = [[1,0,0], [0,2,4], [7,0,9.0]]
        otu_ids = ['OTUa','OTUb','OTUc']
        sample_ids = ['Sa','Sb','Sc']
        result = format_unifrac_sample_mapping(sample_ids, otu_ids, a)
        self.assertEqual(result, ['OTUa\tSa\t1', 'OTUb\tSb\t2', 'OTUb\tSc\t4', 'OTUc\tSa\t7', 'OTUc\tSc\t9.0'])
        
    def test_illumina_data_to_fastq(self):
        """illumina_data_to_fastq functions as expected """
        in1 = ("M10","68","1","1","28680","29475","0","1","AACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACG.","BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB","0")
        expected = """@M10_68:1:1:28680:29475#0/1\nAACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACGN\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB""", 0
        
        self.assertEqual(illumina_data_to_fastq(in1),expected)

        expected12 = """@M10_68:1:1:28680:29475#0/1\nAACGAAAGGCAG\n+\nBBBBBBBBBBBB""", 0
        self.assertEqual(illumina_data_to_fastq(in1,number_of_bases=12),expected12)
        
        # different value in the pass filter field
        in2 = ("M10","68","1","1","28680","29475","0","1","AACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACG.","BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB","1")
        expected = """@M10_68:1:1:28680:29475#0/1\nAACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACGN\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB""", 1
        
        self.assertEqual(illumina_data_to_fastq(in2),expected)
        
    def test_illumina_data_to_fastq_no_pass_filter_field(self):
        """illumina_data_to_fastq functions as expected with no pass filter field"""
        in1 = ("M10","68","1","1","28680","29475","0","1","AACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACG.","BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
        expected = """@M10_68:1:1:28680:29475#0/1\nAACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACGN\n+\nBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB""", 2
        
        self.assertEqual(illumina_data_to_fastq(in1),expected)
        
    def test_format_mapping_html_data(self):
        """ Properly formats html string for mapping file errors/warnings """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'Description']
        mapping_data = [['Sample1','AACCGGTT','ACATATT','Desc_1'],
                        ['Sample2','CCAATTGG','ACATATT','Desc_2']
                       ]
        errors = []
        warnings = []
        
        # no errors or warnings, shouldn't get any popup mouseover data
        
        actual_formatted_html_data = format_mapping_html_data(header,
         mapping_data, errors, warnings)
         
        self.assertEqual(actual_formatted_html_data,
         self.expected_formatted_html_no_errors_warnings)
         
    def test_format_mapping_html_data_errors(self):
        """ Properly formats html string for mapping file errors/warnings """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'Description']
        mapping_data = [['Sample1','AACCGGTT','ACATATT','Desc_1'],
                        ['Sample2','CCAATTGG','ACATATT','Desc_2']
                       ]
        errors = ['problem1\t1,2']
        warnings = []
        
        # Should create a an error popup in the right location
        
        actual_formatted_html_data = format_mapping_html_data(header,
         mapping_data, errors, warnings)
         
        self.assertEqual(actual_formatted_html_data,
         self.expected_formatted_html_errors)
         
    def test_format_mapping_html_data_warnings(self):
        """ Properly formats html string for mapping file errors/warnings """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'Description']
        mapping_data = [['Sample1','AACCGGTT','ACATATT','Desc_1'],
                        ['Sample2','CCAATTGG','ACATATT','Desc_2']
                       ]
        errors = []
        warnings = ['warning1\t2,2']
        
        # Should create a an warning popup in the right location
        
        actual_formatted_html_data = format_mapping_html_data(header,
         mapping_data, errors, warnings)
         
         
        self.assertEqual(actual_formatted_html_data,
         self.expected_formatted_html_warnings)

    def test_format_mapping_html_data_non_location_error(self):
        """ Properly formats html string for mapping file errors/warnings """
        
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
         'Description']
        mapping_data = [['Sample1','AACCGGTT','ACATATT','Desc_1'],
                        ['Sample2','CCAATTGG','ACATATT','Desc_2']
                       ]
        errors = ['error1\t-1,-1']
        warnings = []
        
        # Should list errors with location -1,-1 outside of table
        
        actual_formatted_html_data = format_mapping_html_data(header,
         mapping_data, errors, warnings)
         
        
        self.assertEqual(actual_formatted_html_data,
         self.expected_formatted_html_data_nonloc_error)
         
example_mapping_file = """#SampleID\tcol1\tcol0\tDescription
#this goes after headers
#this too
bsample\tv1_3\tv0_3\td1
asample\taval\tanother\td2"""

expected_formatted_html_no_errors_warnings = """<html>
<head>

<script type="text/javascript" src="./overlib.js"></script>
</head>
<body bgcolor="white"> <h1>No errors or warnings detected.<br></h1><h1>Mapping file error and warning details.</h1>
Notes for interpreting this report:
<ul>
    <li>Errors will be listed in red, warnings in yellow.  
    <li>Mouse over an error or warning in a cell for more details.
    <li>Errors in the header row may mask other errors, so these should be corrected first.
    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>check_id_map.py</tt> until no errors (nor warnings, ideally) are found.
</ul>
<p>
Some general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):
<ul> 
    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.
    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>
    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.
</ul>
General issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>
<table border="2" cellspacing="0" cellpadding="5">

<tr></tr>
<tr>
<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Description</th>
</tr>

<tr>
<tr><th><tt>Sample1</tt></th><th><tt>AACCGGTT</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_1</tt></th></tr><tr><th><tt>Sample2</tt></th><th><tt>CCAATTGG</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_2</tt></th></tr>
</tr>
</table>

</body>
</html>"""

expected_formatted_html_errors = """<html>
<head>

<script type="text/javascript" src="./overlib.js"></script>
</head>
<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>
Notes for interpreting this report:
<ul>
    <li>Errors will be listed in red, warnings in yellow.  
    <li>Mouse over an error or warning in a cell for more details.
    <li>Errors in the header row may mask other errors, so these should be corrected first.
    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>check_id_map.py</tt> until no errors (nor warnings, ideally) are found.
</ul>
<p>
Some general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):
<ul> 
    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.
    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>
    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.
</ul>
General issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>
<table border="2" cellspacing="0" cellpadding="5">

<tr></tr>
<tr>
<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Description</th>
</tr>

<tr>
<tr><th><tt>Sample1</tt></th><th><tt>AACCGGTT</tt></th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib('problem1<br>Location (SampleID,Header Field)<br>Sample1,LinkerPrimerSequence');" onmouseout="return nd();"><font color=white><tt>ACATATT</tt></a></th><th><tt>Desc_1</tt></th></tr><tr><th><tt>Sample2</tt></th><th><tt>CCAATTGG</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_2</tt></th></tr>
</tr>
</table>

</body>
</html>"""

expected_formatted_html_warnings = """<html>
<head>

<script type="text/javascript" src="./overlib.js"></script>
</head>
<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>
Notes for interpreting this report:
<ul>
    <li>Errors will be listed in red, warnings in yellow.  
    <li>Mouse over an error or warning in a cell for more details.
    <li>Errors in the header row may mask other errors, so these should be corrected first.
    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>check_id_map.py</tt> until no errors (nor warnings, ideally) are found.
</ul>
<p>
Some general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):
<ul> 
    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.
    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>
    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.
</ul>
General issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>
<table border="2" cellspacing="0" cellpadding="5">

<tr></tr>
<tr>
<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Description</th>
</tr>

<tr>
<tr><th><tt>Sample1</tt></th><th><tt>AACCGGTT</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_1</tt></th></tr><tr><th><tt>Sample2</tt></th><th><tt>CCAATTGG</tt></th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib('warning1<br>Location (SampleID,Header Field)<br>Sample2,LinkerPrimerSequence');" onmouseout="return nd();"><font color=black><tt>ACATATT</tt></a></th><th><tt>Desc_2</tt></th></tr>
</tr>
</table>

</body>
</html>"""

expected_formatted_html_data_nonloc_error = """<html>
<head>

<script type="text/javascript" src="./overlib.js"></script>
</head>
<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>
Notes for interpreting this report:
<ul>
    <li>Errors will be listed in red, warnings in yellow.  
    <li>Mouse over an error or warning in a cell for more details.
    <li>Errors in the header row may mask other errors, so these should be corrected first.
    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>check_id_map.py</tt> until no errors (nor warnings, ideally) are found.
</ul>
<p>
Some general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):
<ul> 
    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.
    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>
    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.
</ul>
General issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr><td bgcolor="red"><font color="white">error1<font color="black"></td></tr></table><br>
<table border="2" cellspacing="0" cellpadding="5">

<tr></tr>
<tr>
<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Description</th>
</tr>

<tr>
<tr><th><tt>Sample1</tt></th><th><tt>AACCGGTT</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_1</tt></th></tr><tr><th><tt>Sample2</tt></th><th><tt>CCAATTGG</tt></th><th><tt>ACATATT</tt></th><th><tt>Desc_2</tt></th></tr>
</tr>
</table>

</body>
</html>"""

biom1 = """{"rows": [{"id": "tax1", "metadata": {}}, {"id": "tax2", "metadata": {}}, {"id": "tax3", "metadata": {}}, {"id": "tax4", "metadata": {}}, {"id": "endbigtaxon", "metadata": {}}, {"id": "tax6", "metadata": {}}, {"id": "tax7", "metadata": {}}, {"id": "tax8", "metadata": {}}, {"id": "tax9", "metadata": {}}], "format": "Biological Observation Matrix 0.9.0-dev", "data": [[0, 0, 7.0], [0, 1, 4.0], [0, 2, 2.0], [0, 3, 1.0], [1, 0, 1.0], [1, 1, 2.0], [1, 2, 4.0], [1, 3, 7.0], [1, 4, 8.0], [1, 5, 7.0], [1, 6, 4.0], [1, 7, 2.0], [1, 8, 1.0], [2, 5, 1.0], [2, 6, 2.0], [2, 7, 4.0], [2, 8, 7.0], [2, 9, 8.0], [2, 10, 7.0], [2, 11, 4.0], [2, 12, 2.0], [2, 13, 1.0], [3, 10, 1.0], [3, 11, 2.0], [3, 12, 4.0], [3, 13, 7.0], [3, 14, 8.0], [3, 15, 7.0], [3, 16, 4.0], [3, 17, 2.0], [3, 18, 1.0], [4, 15, 1.0], [4, 16, 2.0], [4, 17, 4.0], [4, 18, 7.0], [5, 1, 1.0], [5, 2, 1.0], [6, 6, 2.0], [6, 7, 1.0], [7, 11, 3.0], [7, 12, 1.0], [8, 16, 4.0], [8, 17, 1.0]], "columns": [{"id": "sam1", "metadata": null}, {"id": "sam2", "metadata": null}, {"id": "sam3", "metadata": null}, {"id": "sam4", "metadata": null}, {"id": "sam5", "metadata": null}, {"id": "sam6", "metadata": null}, {"id": "sam7", "metadata": null}, {"id": "sam8", "metadata": null}, {"id": "sam9", "metadata": null}, {"id": "sam_middle", "metadata": null}, {"id": "sam11", "metadata": null}, {"id": "sam12", "metadata": null}, {"id": "sam13", "metadata": null}, {"id": "sam14", "metadata": null}, {"id": "sam15", "metadata": null}, {"id": "sam16", "metadata": null}, {"id": "sam17", "metadata": null}, {"id": "sam18", "metadata": null}, {"id": "sam19", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2520", "matrix_type": "sparse", "shape": [9, 19], "format_url": "http://biom-format.org", "date": "2011-12-20T19:03:28.130403", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

#run unit tests if run from command-line
if __name__ == '__main__':
    main()




