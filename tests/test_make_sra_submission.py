#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main
from qiime.make_sra_submission import (md5_path, safe_for_xml, 
    read_tabular_data, 
    rows_data_as_dicts,
    make_study_links, twocol_data_to_dict, make_study, make_submission,
    make_sample)
"""Tests of the make_study_and_experiment.py file.

Note: these tests assume you are running from within the directory that also
has the templates and sample tabular data.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

class TopLevelTests(TestCase):
    """Top-level tests of functions in make_study_and_experiment"""

    def test_md5_path(self):
        """md5_path should match hand-calculated value"""
        result = md5_path('sra_xml_templates/study_template.xml')
        self.assertEqual(result, 'bcd7ec3afb9fe75ea09f5ac1cfeeb450')

    def test_safe_for_xml(self):
        """safe_for_xml should return string without 'bad' xml entities."""
        s = '''ab&c'd"e<f>'''
        t = safe_for_xml(s)
        self.assertEqual(t, 'ab&amp;c&apos;d&quot;e&lt;f&gt;')

    def test_read_tabular_data(self):
        """read_tabular_data should read simple table"""
        data = """#a\tb\tc d
x\ty\tz

x x\ty y\tc c 

"""
        self.assertEqual(read_tabular_data(data.splitlines()), 
            ([['#a', 'b', 'c d']], [['x','y','z'],['x x', 'y y', 'c c ']]))

    def test_make_study_links(self):
        """make_study_links should return correct study links from pmid."""
        self.assertEqual(make_study_links(19004758), """    <STUDY_LINKS>
      <STUDY_LINK>
        <ENTREZ_LINK>
         <DB>pubmed</DB>
         <ID>19004758</ID>
        </ENTREZ_LINK>
      </STUDY_LINK>
    </STUDY_LINKS>""")

    def test_twocol_data_to_dict(self):
        """twocol_data_to_dict should produce expected result"""
        data = """#a\tb\tc d
x\ty\tz

x x\ty y\tc c 

"""
        header, body = read_tabular_data(data.splitlines())
        self.assertEqual(twocol_data_to_dict(body), {'x':'y', 'x x':'y y'})
        #test the case where data are multiple
        data = """#a\tb\tc d
x\ty\tz

x x\ty y\tc c 
x\ty\tz
x\tc\tv

"""
        header, body = read_tabular_data(data.splitlines())
        self.assertEqual(twocol_data_to_dict(body, is_multiple=True),
            {'x':['y','y','c'], 'x x':['y y']})

    def test_rows_data_as_dicts(self):
        """rows_data_as_dicts should iterate over trimmed row data"""
        data = """#a\tb\tc d
x\t\tz

x x\ty y\t 
aa\tbb\tcc
"""
        header, body = read_tabular_data(data.splitlines())
        result = list(rows_data_as_dicts(header[0], body))
        self.assertEqual(result, [{'a':'x','c d':'z'}, {'a':'x x','b':'y y'},
            {'a':'aa', 'b':'bb','c d':'cc'}])

    def test_make_study(self):
        """make_study should produce expected results given info/template"""
        study_sample_data = open('sra_test_files/study.txt', 'U')
        study_template = open(
            'sra_xml_templates/study_template.xml', 'U').read()
        result = make_study(study_sample_data, study_template)
        expected = open('sra_test_files/example_study.xml', 'U').read()
        self.assertEqual(result, expected)
        #check that it works if pmid field empty
        study_sample_data = open('sra_test_files/study_empty_pmid.txt', 'U')
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, open(
            'sra_test_files/example_study_no_pmid.xml', 'U').read())
        #check that it works if pmid field not supplied
        study_sample_data = open('sra_test_files/study_missing_pmid.txt', 'U')
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, open(
            'sra_test_files/example_study_no_pmid.xml', 'U').read())

    def test_make_submission(self):
        """make_submission should produce expected results given info/template"""
        submission_sample_data = open('sra_test_files/submission.txt', 'U')
        submission_template = open('sra_xml_templates/submission_template.xml',
            'U').read()
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
        expected = open(
            'sra_test_files/example_submission_study_sample_only.xml').read()
        self.assertEqual(result, expected)
        #test that it works with a file
        submission_sample_data = open(
            'sra_test_files/submission_with_file.txt', 'U')
        submission_template = open(
            'sra_xml_templates/submission_template.xml', 'U').read()
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
        expected = open(
            'sra_test_files/example_submission_study_sample_only_with_file.xml').read()
        self.assertEqual(result, expected)

    def test_make_sample(self):
        """make_sample should produce expected reuslts given info/template"""
        sample_data = open('sra_test_files/sample_small.txt', 'U')
        sample_template = open(
            'sra_xml_templates/sample_template.xml', 'U').read()
        result = make_sample(sample_data, sample_template)
        expected = open('sra_test_files/example_sample.xml', 'U').read()
        self.assertEqual(result, expected)

if __name__ == '__main__':
    main()
