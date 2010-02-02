#!/usr/bin/env python
import os
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.make_sra_submission import (
    md5_path, safe_for_xml, read_tabular_data, rows_data_as_dicts,
    make_study_links, twocol_data_to_dict, make_study, make_submission,
    make_sample, trim_quotes, defaultdict, group_lines_by_field,
    parse_command_line_parameters, write_xml_generic,
    make_run_and_experiment)
from qiime.util import get_qiime_project_dir

"""Tests of the make_study_and_experiment.py file.

Note: these tests assume you are running from within the directory that also
has the templates and sample tabular data.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

class TopLevelTests(TestCase):
    """Top-level tests of functions in make_study_and_experiment"""

    def setUp(self):
        """ """
        qiime_dir = get_qiime_project_dir()
        self.sra_xml_templates_dir = os.path.join(
            qiime_dir, 'tests', 'sra_xml_templates')
        self.sra_test_files_dir = os.path.join(
            qiime_dir, 'tests', 'sra_test_files')

        self.submission_with_file_fp = \
         get_tmp_filename(prefix='make_sra_submission_tests')
        open(self.submission_with_file_fp,'w').write(\
         submission_with_file_text %  self.sra_test_files_dir)
        self.files_to_remove = [self.submission_with_file_fp]
        
    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_md5_path(self):
        """md5_path should match hand-calculated value"""
        result = md5_path('%s/study_template.xml' % self.sra_xml_templates_dir)
        self.assertEqual(result, 'bcd7ec3afb9fe75ea09f5ac1cfeeb450')

    def test_safe_for_xml(self):
        """safe_for_xml should return string without 'bad' xml entities."""
        s = 'ab&c\'d"e<f>'
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
        study_sample_data = \
         open('%s/study.txt' % self.sra_test_files_dir, 'U')
        study_template = open(
            '%s/study_template.xml' % self.sra_xml_templates_dir, 'U').read()
        result = make_study(study_sample_data, study_template)
        expected = \
         open('%s/example_study.xml' % self.sra_test_files_dir, 'U').read()
        self.assertEqual(result, expected)
        #check that it works if pmid field empty
        study_sample_data = \
         open('%s/study_empty_pmid.txt' % self.sra_test_files_dir, 'U')
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, open('%s/example_study_no_pmid.xml' %\
            self.sra_test_files_dir, 'U').read())
        #check that it works if pmid field not supplied
        study_sample_data = \
         open('%s/study_missing_pmid.txt' % self.sra_test_files_dir, 'U')
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, open('%s/example_study_no_pmid.xml' %\
         self.sra_test_files_dir, 'U').read())

    def test_make_submission(self):
        """make_submission should produce expected results given info/template"""
        submission_sample_data = \
         open('%s/submission.txt' % self.sra_test_files_dir, 'U')
        submission_template = \
         open('%s/submission_template.xml' % self.sra_xml_templates_dir,'U').read()
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
        expected = open('%s/example_submission_study_sample_only.xml' \
            % self.sra_test_files_dir).read()
        self.assertEqual(result, expected)
        #test that it works with a file
        submission_sample_data = open(
            self.submission_with_file_fp, 'U')
        submission_template = open('%s/submission_template.xml' \
         % self.sra_xml_templates_dir, 'U').read()
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
            
        expected = submission_with_file_expected % self.sra_test_files_dir
        self.assertEqual(result, expected)

    def test_make_sample(self):
        """make_sample should produce expected reuslts given info/template"""
        sample_data = open('%s/sample_small.txt' \
         % self.sra_test_files_dir, 'U')
        sample_template = open('%s/sample_template.xml' %\
         self.sra_xml_templates_dir, 'U').read()
        result = make_sample(sample_data, sample_template)
        expected = open('%s/example_sample.xml' %\
         self.sra_test_files_dir, 'U').read()
        self.assertEqual(result, expected)

    def test_trim_quotes(self):
        self.assertEqual(trim_quotes('"abcd"'), 'abcd')

    def test_group_lines_by_field(self):
        lines = [
            ['x',   'why', 'b'],
            ['x x', 'y y', 'a'],
            ['wx ', 'y y', 'c']
            ]
        observed = group_lines_by_field(lines, 1).items()
        observed.sort()
        expected = [
            ('why', [['x',   'why', 'b']]),
            ('y y', [['x x', 'y y', 'a'], ['wx ', 'y y', 'c']]),
            ]
        self.assertEqual(observed, expected)

    def test_parse_command_line_parameters(self):
        argv = [
            '-a', 'sample.tsv',
            '-t', 'study.tsv',
            '-u', 'submission.tsv',
            '-e', 'experiment.tsv',
            '-s', 'sffs/',
            ]
        opts, args = parse_command_line_parameters(argv)
        self.assertEqual(opts.input_sample_fp, 'sample.tsv')
        self.assertEqual(opts.input_study_fp, 'study.tsv')
        self.assertEqual(opts.input_submission_fp, 'submission.tsv')
        self.assertEqual(opts.input_experiment_fp, 'experiment.tsv')
        self.assertEqual(opts.sff_dir, 'sffs/')

    def test_write_xml_generic(self):
        input_file = tempfile.NamedTemporaryFile()
        input_file.write('abc')
        input_file.seek(0)
        
        template_file = tempfile.NamedTemporaryFile()
        template_file.write('<xml>%s</xml>')
        template_file.seek(0)

        def simple_xml_f(lines, template):
            return template % ''.join(lines)
        
        observed_fp = write_xml_generic(
            input_file.name, template_file.name, simple_xml_f)
        observed = open(observed_fp).read()
        self.assertEqual(observed, '<xml>abc</xml>')

        self.files_to_remove = [observed_fp]

    def test_make_run_and_experiment(self):
        expt_lines = ['#STUDY_REF\tEXPERIMENT_ALIAS\n']
        observed_experiment_xml, observed_run_xml = make_run_and_experiment(
            expt_lines, '/tmp')
        expected_experiment_xml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
            '\n'
            '</EXPERIMENT_SET>'
            )
        self.assertEqual(observed_experiment_xml, expected_experiment_xml)
        expected_run_xml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
            '\n'
            '</RUN_SET>'
            )
        self.assertEqual(observed_run_xml, expected_run_xml)

submission_with_file_text = '''#Field	Value	Example	Comments
accession	SRA003492	SRA003492	"leave blank if not assigned yet, e.g. if new submission"
submission_id	fierer_hand_study	fierer_hand_study	internally unique id for the submission
center_name	CCME	CCME	name of the center preparing the submission
submission_comment	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Free-text comments regarding submission
lab_name	Knight	Knight	"name of lab preparing submission, can differ from center (usually refers to the PI's info, not the sequencing center's)"
submission_date	2009-10-22T01:23:00-05:00	2009-10-22T01:23:00-05:00	timestamp of submission
CONTACT	Rob Knight;Rob.Knight@Colorado.edu	Rob Knight;Rob.Knight@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
CONTACT	Noah Fierer;Noah.Fierer@Colorado.edu	Noah Fierer;Noah.Fierer@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
study	study.xml	fierer_hand_study.study.xml	"leave blank if not submitting study, put in filename otherwise"
sample	sample.xml	fierer_hand_study.sample.xml	"leave blank if not submitting sample, put in filename otherwise"
experiment		fierer_hand_study.experiment.xml	"leave blank if not submitting experiment, put in filename otherwise"
run		fierer_hand_study.run.xml	"leave blank if not submitting run, put in filename otherwise"
file	%s/fake_hand_data_for_sra.tgz	fierer_hand_study.seqs.tgz	"leave blank if not submitting sequence data, put in filename otherwise"'''

submission_with_file_expected = '''<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
 accession="SRA003492"
 submission_id="fierer_hand_study"
 center_name="CCME"
 submission_comment="Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"
 lab_name="Knight"
 submission_date="2009-10-22T01:23:00-05:00"
>
 <CONTACTS>
    <CONTACT name="Rob Knight" inform_on_status="Rob.Knight@Colorado.edu" inform_on_error="Rob.Knight@Colorado.edu"/>
    <CONTACT name="Noah Fierer" inform_on_status="Noah.Fierer@Colorado.edu" inform_on_error="Noah.Fierer@Colorado.edu"/>
 </CONTACTS>
 <ACTIONS>
   <ACTION><ADD source="sample.xml" schema="sample" notes="sample metadata"/></ACTION>
   <ACTION><ADD source="study.xml" schema="study" notes="study metadata"/></ACTION>
   <ACTION><RELEASE/></ACTION>
 </ACTIONS>
 <FILES>
 <FILE filename="%s/fake_hand_data_for_sra.tgz" checksum_method="MD5" checksum="d41d8cd98f00b204e9800998ecf8427e"/>
 </FILES>
</SUBMISSION>
'''

if __name__ == '__main__':
    main()

