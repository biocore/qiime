#!/usr/bin/env python
#unit tests for format.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","Jeremy Widmann","Jens Reeder"] 
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

from os import remove
from cogent.util.unit_test import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from cogent.app.util import get_tmp_filename
from numpy import array, nan
from qiime.format import (format_distance_matrix, format_otu_table,
    format_coords, build_prefs_string, format_matrix, format_map_file,
    format_histograms, write_Fasta_from_name_seq_pairs, 
    format_unifrac_sample_mapping)

class TopLevelTests(TestCase):
    """Tests of top-level module functions."""

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

    def test_format_otu_table(self):
        """format_otu_table should return tab-delimited table"""
        a = array([[1,2,3],[4,5,2718281828459045]])
        samples = ['a','b','c']
        otus = [1,2]
        taxa = ['Bacteria','Archaea']
        res = format_otu_table(samples, otus, a)
        self.assertEqual(res,
            '#Full OTU Counts\n#OTU ID\ta\tb\tc\n1\t1\t2\t3\n2\t4\t5\t2718281828459045')
        res = format_otu_table(samples, otus, a, taxa)
        self.assertEqual(res,
            '#Full OTU Counts\n#OTU ID\ta\tb\tc\tConsensus Lineage\n1\t1\t2\t3\tBacteria\n2\t4\t5\t2718281828459045\tArchaea')
        self.assertRaises(ValueError, format_otu_table, samples, [1,2,3], a)

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
        exp_string = \
        """{\n'background_color':'black',\n\n'sample_coloring':\n\t{\n\t\t'First':\n\t\t{\n\t\t\t'column':'First',\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t\t},\n\t\t'Second':\n\t\t{\n\t\t\t'column':'Second',\n\t\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t\t}\n\t},\n'MONTE_CARLO_GROUP_DISTANCES':\n\t{\n\t\t'First': 10,\n\t\t'Second': 10\n\t},\n'FIELDS':\n\t[\n\t\t'Second',\n\t\t'First'\n\t],\n'taxonomy_coloring':\n\t{\n\t\t'Level_1':\n\t\t{\n\t\t\t'column':'1',\n\t\t\t'colors':\n\t\t\t{\n\t\t\t\t'Root;Bacteria':('red0',(0,100,100))\n\t\t\t}\n\t\t}\n\t}\n}"""
        obs_string = build_prefs_string(mapping_headers_to_use, \
                                    background_color, monte_carlo_dist, headers,
                                    otu_ids)
                                    
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
        self.assertEqual(format_histograms(array([2,1,0,2,0,0]),
            array([0,0,0,2,0,1]), array([100,110,120,130,140,150,160])),
            """Length\tBefore\tAfter
100\t2\t0
110\t1\t0
120\t0\t0
130\t2\t2
140\t0\t0
150\t0\t1""")

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

#run unit tests if run from command-line
if __name__ == '__main__':
    main()




