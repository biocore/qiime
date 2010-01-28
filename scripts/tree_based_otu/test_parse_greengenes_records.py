#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from parse_greengenes_records import is_blank, split_delim_f, DemarkedParser, \
        taxonomy_delim_string

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

class ParseGreengenesRecordsTests(TestCase):
    def setUp(self):
        pass

    def test_is_blank(self):
        """Determine if a line is blank or not"""
        blank_line = ''
        with_spaces = '       '
        not_empty = '   asdasd'

        self.assertTrue(is_blank(blank_line))
        self.assertTrue(is_blank(with_spaces))
        self.assertFalse(is_blank(not_empty))

    def test_split_delim_f(self):
        """Splits a delimited line once"""
        ### this has a check for a specific field...
        f = split_delim_f('=')
        input1 = 'foo=bar'
        input2 = 'foo='
        exp1 = ('foo','bar')
        exp2 = ('foo', None)
        self.assertEqual(f(input1), exp1)
        self.assertEqual(f(input2), exp2)

    def test_DemarkedParser(self):
        """Returns records between starting and ending points"""
        lines = rec_data.splitlines()
        records = list(DemarkedParser(lines, '=', start='my_starting', \
                end='my_ending'))
        exp_r1 = {'a':'1','b':'2','c':'3','d':None,'e':'5'}
        exp_r2 = {'q':'asdasd','c':'taco'}
        self.assertEqual(records[0], exp_r1)
        self.assertEqual(records[1], exp_r2)

    def test_taxonomy_delim_string(self):
        """Returns a tab delim string of the taxonomy information"""
        self.assertRaises(KeyError, taxonomy_delim_string, {})
        exp = '\t'.join(map(str, range(9)))
        rec = dict([(k,i) for i,k in enumerate(taxonomy_fields)])
        rec['foo'] = 'bar'
        self.assertEqual(taxonomy_delim_string(rec), exp)

taxonomy_fields = ['prokMSA_id','ncbi_tax_id','ncbi_acc_w_ver',
                   'core_set_member','ncbi_tax_string_format_2',
                   'Hugenholtz_tax_string_format_2',
                   'Ludwig_tax_string_format_2',
                   'Pace_tax_string_format_2',
                   'RDP_tax_string_format_2']

rec_data = """my_starting
a=1
b=2
c=3
d=
e=5
my_ending

my_starting
q=asdasd
c=taco
my_ending
"""

if __name__ == '__main__':
    main()

