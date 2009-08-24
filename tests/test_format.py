#!/usr/bin/env python
#unit tests for format.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main
from numpy import array, nan
from qiime.format import (format_distance_matrix, format_otu_table,
    format_coords, format_array_row, format_rarefaction_table)

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

    def test_format_otu_table(self):
        """format_otu_table should return tab-delimited table"""
        a = array([[1,2,3],[4,5,6]])
        samples = ['a','b','c']
        otus = [1,2]
        taxa = ['Bacteria','Archaea']
        res = format_otu_table(samples, otus, a)
        self.assertEqual(res,
            '#Full OTU Counts\n#OTU ID\ta\tb\tc\n1\t1\t2\t3\n2\t4\t5\t6')
        res = format_otu_table(samples, otus, a, taxa)
        self.assertEqual(res,
            '#Full OTU Counts\n#OTU ID\ta\tb\tc\tConsensus Lineage\n1\t1\t2\t3\tBacteria\n2\t4\t5\t6\tArchaea')
        self.assertRaises(ValueError, format_otu_table, samples, [1,2,3], a)

    def test_format_coords(self):
        """format_coords should return tab-delimited table of coords"""
        a = array([[1,2,3],[4,5,6],[7,8,9]])
        header = list('abc')
        eigvals = [2,4,6]
        pct_var = [3,2,1]
        res = format_coords(header, a, eigvals, pct_var)
        self.assertEqual(res, "pc vector number\t1\t2\t3\na\t1\t2\t3\nb\t4\t5\t6\nc\t7\t8\t9\n\n\neigvals\t2\t4\t6\n% variation explained\t3\t2\t1")

    def test_format_array_row(self):
        """format_array_row should strip zeros and nans out of array"""
        a = array([1,0,50,nan,3], int)
        res = format_array_row(a)
        self.assertEqual(res, ['1','','50','','3'])

    def test_format_rarefaction_table(self):
        """format_rarefaction_table should strip bad vals, omit zeros"""
        header = list('abc')
        sizes = [100, 200, 300, 400]
        table = array([[50,60,70],[80,0,100],[110,120,130],[140,0,160]], int)
        res = format_rarefaction_table(header, sizes, 3, table)
        self.assertEqual(res, 'n\ta\tb\tc\n100\t50\t60\t70\n300\t110\t120\t130')
        res = format_rarefaction_table(header, sizes, 0, table)
        self.assertEqual(res, 'n\ta\tb\tc\n100\t50\t60\t70\n200\t80\t\t100\n300\t110\t120\t130\n400\t140\t\t160')

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
