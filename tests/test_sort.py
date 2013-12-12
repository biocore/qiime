#!/usr/bin/env python
# File created on 15 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.parse import parse_mapping_file
from biom.parse import parse_biom_table_str
from qiime.sort import (sort_sample_ids_by_mapping_value,
                        sort_fasta_by_abundance, natsort,
                        natsort_case_insensitive, sort_otu_table,
                        sort_otu_table_by_mapping_field, signed_natsort)


class SortTests(TestCase):

    def setUp(self):
        self.mapping_f1 = mapping_f1.split('\n')
        self.otu_table1 = otu_table1
        self.mapping_f2 = mapping_f2.split('\n')
        self.age_sorted_otu_table1 = age_sorted_otu_table1
        self.name_sorted_otu_table1 = name_sorted_otu_table1
        self.nothing_sorted_otu_table1 = nothing_sorted_otu_table1
        self.otu_table1_bad_sampleID = otu_table1_bad_sampleID
        self.dirs_to_remove = []
        self.files_to_remove = []

    def tearDown(self):
        for dir in self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)
        remove_files(self.files_to_remove)

    def test_sort_sample_ids_by_mapping_value(self):
        """ sort_sample_ids_by_mapping_value functions as expected """
        actual = sort_sample_ids_by_mapping_value(mapping_file=self.mapping_f1,
                                                  field='days_since_epoch',
                                                  field_type_f=float)
        expected = zip(['NotInOtuTable', '1', 'Z2', 'Z1', 'A'],
                       [0.0, 5.7, 10, 23, 400000])
        self.assertEqual(actual, expected)

    def test_sort_sample_ids_by_mapping_value_error(self):
        """ sort_sample_ids_by_mapping_value handles errors """
        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='years_since_spoch',
                          field_type_f=float)

        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='Something',
                          field_type_f=float)

    def test_sort_fasta_by_abundance(self):
        """sort_fasta_by_abundance functions as expected"""
        class FakeOutF(object):

            def __init__(self):
                self.s = ""

            def write(self, s):
                self.s += s

        actual = FakeOutF()
        expected = ""
        sort_fasta_by_abundance([], actual)
        self.assertEqual(actual.s, expected)

        # no sorting necessary
        actual = FakeOutF()
        expected1 = "\n".join(['>s1', 'ACCGT', '>s2 comment', 'ATTGC', ''])
        expected2 = "\n".join(['>s2 comment', 'ATTGC', '>s1', 'ACCGT', ''])
        sort_fasta_by_abundance(
            ['>s1',
             'ACCGT',
             '>s2 comment',
             'ATTGC'],
            actual)
        # order is unimportant here
        self.assertTrue(actual.s == expected1 or actual.s == expected2)

        # sorting necessary
        actual = FakeOutF()
        inseqs = ['>s1', 'ACCGT',
                  '>s2 comment', 'ATTGC',
                  '>s3 blah', 'ATTGC']
        expected = "\n".join(['>s2 comment', 'ATTGC',
                              '>s3 blah', 'ATTGC',
                              '>s1', 'ACCGT', ''])
        sort_fasta_by_abundance(inseqs, actual)
        self.assertEqual(actual.s, expected)

    def test_natsort(self):
        """natsort should perform numeric comparisons on strings"""
        # string with alpha and numerics sort correctly
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s),
                         'sample1 sample2 sample11 sample12'.split())
        s.reverse()
        self.assertEqual(natsort(s),
                         'sample1 sample2 sample11 sample12'.split())
        self.assertEqual(natsort(list('cba321')), list('123abc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')), list('abcd'))

        # string of ints sort correctly
        self.assertEqual(natsort(['11', '2', '1', '0']),
                         ['0', '1', '2', '11'])

        # strings of floats sort correctly
        self.assertEqual(natsort(['1.11', '1.12', '1.00', '0.009']),
                         ['0.009', '1.00', '1.11', '1.12'])

        # string of ints sort correctly
        self.assertEqual(
            natsort([('11', 'A'), ('2', 'B'), ('1', 'C'), ('0', 'D')]),
            [('0', 'D'), ('1', 'C'), ('2', 'B'), ('11', 'A')])

    #
    def test_natsort_case_insensitive(self):
        """natsort should perform numeric comparisons on strings and is
           _not_ case-sensitive"""

        # string with alpha and numerics sort correctly
        s = [
            'sample1',
            'sample2',
            'sample11',
            'sample12',
            'SAmple1',
            'Sample2']

        # expected values
        exp_natsort = ['SAmple1', 'Sample2', 'sample1', 'sample2', 'sample11',
                       'sample12']
        exp_natsort_case_insensitive = ['sample1', 'SAmple1', 'sample2',
                                        'Sample2', 'sample11', 'sample12']

        # test natsort
        self.assertEqual(natsort(s), exp_natsort)
        # test natsort_case_insensitive
        self.assertEqual(natsort_case_insensitive(s),
                         exp_natsort_case_insensitive)

        s.reverse()
        # test natsort
        self.assertEqual(natsort(s), exp_natsort)
        # test natsort_case_insensitive
        self.assertEqual(natsort(list('cbaA321')), list('123Aabc'))

        # strings with alpha only sort correctly
        self.assertEqual(natsort_case_insensitive(list('cdBa')), list('aBcd'))

        # string of ints sort correctly
        self.assertEqual(natsort_case_insensitive(['11', '2', '1', '0']),
                         ['0', '1', '2', '11'])

        # strings of floats sort correctly
        self.assertEqual(natsort_case_insensitive(['1.11', '1.12', '1.00',
                                                  '0.009']), ['0.009', '1.00',
                                                              '1.11', '1.12'])

        # string of ints sort correctly
        self.assertEqual(natsort_case_insensitive([('11', 'A'), ('2', 'B'),
                                                  ('1', 'C'), ('0', 'D')]),
                         [('0', 'D'), ('1', 'C'),
                          ('2', 'B'), ('11', 'A')])

    def test_sort_otu_table_by_mapping_field_all_values_differ(self):
        """ sort_otu_table_by_mapping_field fns when all values differ"""

        actual = sort_otu_table_by_mapping_field(
            parse_biom_table_str(self.otu_table1),
            parse_mapping_file(
                self.mapping_f2),
            sort_field="Age")
        expected = parse_biom_table_str(self.age_sorted_otu_table1)
        self.assertEqual(actual, expected)

    def test_sort_otu_table(self):
        """ sort_otu_table fns as expected """

        actual = sort_otu_table(parse_biom_table_str(self.otu_table1),
                                ['NA', 'Key', 'Fing'])
        expected = parse_biom_table_str(self.age_sorted_otu_table1)
        self.assertEqual(actual, expected)

    def test_sort_otu_table_error(self):
        """ sort_otu_table handles errors """

        self.assertRaises(ValueError, sort_otu_table,
                          parse_biom_table_str(self.otu_table1), ['NA', 'Key', 'Fing', 'Key'])
        self.assertRaises(KeyError, sort_otu_table,
                          parse_biom_table_str(self.otu_table1), ['NA', 'Key'])

    def test_sort_otu_table_by_mapping_field_some_values_differ(self):
        """ sort_otu_table fns when some values differ"""

        actual = sort_otu_table_by_mapping_field(
            parse_biom_table_str(self.otu_table1),
            parse_mapping_file(
                self.mapping_f2),
            sort_field="Nothing")
        expected = parse_biom_table_str(self.nothing_sorted_otu_table1)
        self.assertEqual(actual, expected)

    def test_sort_otu_table_by_mapping_field_some_values_same(self):
        """ sort_otu_table_by_mapping_field fns when all values are the same"""

        actual = sort_otu_table_by_mapping_field(
            parse_biom_table_str(self.otu_table1),
            parse_mapping_file(
                self.mapping_f2),
            sort_field="Name")
        expected = parse_biom_table_str(self.name_sorted_otu_table1)
        self.assertEqual(actual, expected)

    def test_sort_otu_table_by_mapping_field_error(self):
        """ sort_otu_table_by_mapping_field fails on samples in otu table but not mapping"""

        self.assertRaises(KeyError, sort_otu_table_by_mapping_field,
                          parse_biom_table_str(
                              self.otu_table1_bad_sampleID),
                          parse_mapping_file(self.mapping_f2),
                          sort_field="Age")

    def test_signed_sort(self):
        """Test correct sorting of different data types"""

        # an empty list must be returned when an empty list needs to be sorted
        self.assertEqual(signed_natsort([]), [])

        # tuples that can be sorted by type-casting the first element
        test_list = [('9', 'SampleA'), ('-1', 'SampleD'), ('7', 'SampleC'),
                     ('-2', 'SampleE'), ('-0.11',
                                         'SampleF'), ('17.11', 'SampleB'),
                     ('100', 'SampleG'), ('13', 'SampleH')]
        expected_result = [('-2', 'SampleE'), ('-1', 'SampleD'),
                           ('-0.11', 'SampleF'), ('7',
                                                  'SampleC'), ('9', 'SampleA'),
                           ('13', 'SampleH'), ('17.11', 'SampleB'), ('100', 'SampleG')]

        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # tuples that must be sorted alphabetically
        test_list = [('Cygnus', 'SampleA'), ('Cepheus', 'SampleD'),
                     ('Auriga', 'SampleC'), ('Grus',
                                             'SampleE'), ('Hydra', 'SampleF'),
                     ('Carina', 'SampleB'), ('Orion', 'SampleG'), ('Lynx', 'SampleH')]
        expected_result = [('Auriga', 'SampleC'), ('Carina', 'SampleB'),
                           ('Cepheus', 'SampleD'), ('Cygnus',
                                                    'SampleA'), ('Grus', 'SampleE'),
                           ('Hydra', 'SampleF'), ('Lynx', 'SampleH'), ('Orion', 'SampleG')]

        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed case, tuples will be sorted alpha-numerically
        test_list = [('Cygnus', 'SampleA'), ('Cepheus', 'SampleD'),
                     ('Auriga', 'SampleC'), ('Grus',
                                             'SampleE'), ('-0.11', 'SampleF'),
                     ('17.11', 'SampleB'), ('100', 'SampleG'), ('Lynx', 'SampleH')]
        expected_result = [('17.11', 'SampleB'), ('100', 'SampleG'),
                           ('-0.11', 'SampleF'), ('Auriga',
                                                  'SampleC'), ('Cepheus', 'SampleD'),
                           ('Cygnus', 'SampleA'), ('Grus', 'SampleE'), ('Lynx', 'SampleH')]

        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed case just a list
        test_list = ['foo', 'bar', '-100', '12', 'spam', '4', '-1']
        expected_result = ['4', '12', '-1', '-100', 'bar', 'foo', 'spam']

        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # list of elements that can be type-casted
        test_list = ['0', '1', '14', '12', '-15', '4', '-1']
        expected_result = ['-15', '-1', '0', '1', '4', '12', '14']

        output = signed_natsort(test_list)
        self.assertEquals(output, expected_result)

        # mixed dict case
        test_dict = {
            'foo': 'a', 'bar': 'b', '-100': '1', '12': '11', 'spam': 'q',
            '4': '11', '-1': 'e'}
        expected_result = ['4', '12', '-1', '-100', 'bar', 'foo', 'spam']

        output = signed_natsort(test_dict)
        self.assertEquals(output, expected_result)

        # dict where the keys can be type-casted
        test_dict = {
            '0': 'foo', '1': 'bar', '14': 'stand', '12': 'eggs', '-15': 'q',
            '4': 'b', '-1': 'h'}
        expected_result = ['-15', '-1', '0', '1', '4', '12', '14']

        output = signed_natsort(test_dict)
        self.assertEquals(output, expected_result)


mapping_f1 = """#SampleID\tSomething\tdays_since_epoch
Z1\t42\t23
Z2\thello\t10
A\t4\t400000
1\tr\t5.7
NotInOtuTable\tf\t0"""

otu_table1 = '{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Propionibacterineae", "Propionibacterium"]}}, {"id": "1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Lactobacillales", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "2", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Gordoniaceae", "Corynebacteriaceae"]}}, {"id": "3", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Staphylococcaceae"]}}, {"id": "4", "metadata": {"taxonomy": ["Bacteria", "Cyanobacteria", "Chloroplasts", "vectors"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 19111.0], [0, 1, 44536.0], [0, 2, 42.0], [1, 0, 1216.0], [1, 1, 3500.0], [1, 2, 6.0], [2, 0, 1803.0], [2, 1, 1184.0], [2, 2, 2.0], [3, 0, 1722.0], [3, 1, 4903.0], [3, 2, 17.0], [4, 0, 589.0], [4, 1, 2074.0], [4, 2, 34.0]], "columns": [{"id": "Fing", "metadata": null}, {"id": "Key", "metadata": null}, {"id": "NA", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2572", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T20:26:01.408138", "type": "OTU table", "id": null, "matrix_element_type": "float"}'

otu_table1_bad_sampleID = '{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Propionibacterineae", "Propionibacterium"]}}, {"id": "1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Lactobacillales", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "2", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Gordoniaceae", "Corynebacteriaceae"]}}, {"id": "3", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Staphylococcaceae"]}}, {"id": "4", "metadata": {"taxonomy": ["Bacteria", "Cyanobacteria", "Chloroplasts", "vectors"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 19111.0], [0, 1, 44536.0], [0, 2, 42.0], [1, 0, 1216.0], [1, 1, 3500.0], [1, 2, 6.0], [2, 0, 1803.0], [2, 1, 1184.0], [2, 2, 2.0], [3, 0, 1722.0], [3, 1, 4903.0], [3, 2, 17.0], [4, 0, 589.0], [4, 1, 2074.0], [4, 2, 34.0]], "columns": [{"id": "Fing", "metadata": null}, {"id": "Key", "metadata": null}, {"id": "NotInMapping", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2572", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T20:19:02.153603", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
age_sorted_otu_table1 = '{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Propionibacterineae", "Propionibacterium"]}}, {"id": "1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Lactobacillales", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "2", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Gordoniaceae", "Corynebacteriaceae"]}}, {"id": "3", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Staphylococcaceae"]}}, {"id": "4", "metadata": {"taxonomy": ["Bacteria", "Cyanobacteria", "Chloroplasts", "vectors"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 42.0], [0, 1, 44536.0], [0, 2, 19111.0], [1, 0, 6.0], [1, 1, 3500.0], [1, 2, 1216.0], [2, 0, 2.0], [2, 1, 1184.0], [2, 2, 1803.0], [3, 0, 17.0], [3, 1, 4903.0], [3, 2, 1722.0], [4, 0, 34.0], [4, 1, 2074.0], [4, 2, 589.0]], "columns": [{"id": "NA", "metadata": null}, {"id": "Key", "metadata": null}, {"id": "Fing", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2572", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T20:19:06.444921", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
nothing_sorted_otu_table1 = '{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Propionibacterineae", "Propionibacterium"]}}, {"id": "1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Lactobacillales", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "2", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Gordoniaceae", "Corynebacteriaceae"]}}, {"id": "3", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Staphylococcaceae"]}}, {"id": "4", "metadata": {"taxonomy": ["Bacteria", "Cyanobacteria", "Chloroplasts", "vectors"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 19111.0], [0, 1, 42.0], [0, 2, 44536.0], [1, 0, 1216.0], [1, 1, 6.0], [1, 2, 3500.0], [2, 0, 1803.0], [2, 1, 2.0], [2, 2, 1184.0], [3, 0, 1722.0], [3, 1, 17.0], [3, 2, 4903.0], [4, 0, 589.0], [4, 1, 34.0], [4, 2, 2074.0]], "columns": [{"id": "Fing", "metadata": null}, {"id": "NA", "metadata": null}, {"id": "Key", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2572", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T20:19:10.397288", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
name_sorted_otu_table1 = '{"rows": [{"id": "0", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Propionibacterineae", "Propionibacterium"]}}, {"id": "1", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Lactobacillales", "Lactobacillales", "Streptococcaceae", "Streptococcus"]}}, {"id": "2", "metadata": {"taxonomy": ["Bacteria", "Actinobacteria", "Actinobacteridae", "Gordoniaceae", "Corynebacteriaceae"]}}, {"id": "3", "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Alicyclobacillaceae", "Bacilli", "Staphylococcaceae"]}}, {"id": "4", "metadata": {"taxonomy": ["Bacteria", "Cyanobacteria", "Chloroplasts", "vectors"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 19111.0], [0, 1, 44536.0], [0, 2, 42.0], [1, 0, 1216.0], [1, 1, 3500.0], [1, 2, 6.0], [2, 0, 1803.0], [2, 1, 1184.0], [2, 2, 2.0], [3, 0, 1722.0], [3, 1, 4903.0], [3, 2, 17.0], [4, 0, 589.0], [4, 1, 2074.0], [4, 2, 34.0]], "columns": [{"id": "Fing", "metadata": null}, {"id": "Key", "metadata": null}, {"id": "NA", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2572", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T20:19:15.978029", "type": "OTU table", "id": null, "matrix_element_type": "float"}'
# values in 'Age' column sort differently incorrectly if sorted as strings
mapping_f2 = """
#SampleID	Name	Age	Nothing
Fing	Blah	11	4
Key	Blah	2	5
NA	Blah	1	4
IgnoredSample	Blah	1	4"""


if __name__ == "__main__":
    main()
