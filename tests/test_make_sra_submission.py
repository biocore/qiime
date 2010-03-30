#!/usr/bin/env python
from cStringIO import StringIO
import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.make_sra_submission import (
    md5_path, safe_for_xml, read_tabular_data, rows_data_as_dicts,
    make_study_links, twocol_data_to_dict, make_study, make_submission,
    make_sample, trim_quotes, defaultdict, group_lines_by_field,
    write_xml_generic, make_run_and_experiment)
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
__version__ = "0.92-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

class TopLevelTests(TestCase):
    """Top-level tests of functions in make_study_and_experiment"""
    def setUp(self):
        """ """
        self.files_to_remove = []

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_md5_path(self):
        """md5_path should match hand-calculated value"""
        template_fp = tempfile.mktemp(suffix='.xml')
        open(template_fp, 'wb').write(study_template)
        result = md5_path(template_fp)
        self.assertEqual(result, 'bcd7ec3afb9fe75ea09f5ac1cfeeb450')

        self.files_to_remove.append(template_fp)

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
        study_sample_data = StringIO(study_txt)
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, study_xml)

        # Test when pmid field is empty
        study_sample_data = StringIO(study_pmid_empty_txt)
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, study_pmid_missing_xml)

        # Test when pmid field is missing
        study_sample_data = StringIO(study_pmid_missing_txt)
        result = make_study(study_sample_data, study_template)
        self.assertEqual(result, study_pmid_missing_xml)

    def test_make_submission(self):
        """make_submission should produce expected results given info/template"""
        submission_sample_data = StringIO(submission_txt)
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
        self.assertEqual(result, submission_xml)

        # TODO Rewrite using only temp files.
        fake_tgz_file = tempfile.NamedTemporaryFile(suffix='.tgz')
        submission_sample_data = StringIO(
            submission_with_file_txt % fake_tgz_file.name)
        result = make_submission(submission_sample_data, submission_template,
            {'study':'study.xml', 'sample':'sample.xml'})
        self.assertEqual(result, submission_with_file_xml % fake_tgz_file.name)

    def test_make_sample(self):
        """make_sample should produce expected reuslts given info/template"""
        sample_data = StringIO(sample_txt)
        result = make_sample(sample_data, sample_template)
        self.assertEqual(result, sample_xml)

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
            '</RUN_SET>\n'
            )
        self.assertEqual(observed_run_xml, expected_run_xml)

submission_with_file_txt = '''#Field	Value	Example	Comments
accession	SRA003492	SRA003492	"leave blank if not assigned yet, e.g. if new submission"
submission_id	fierer_hand_study	fierer_hand_study	internally unique id for the submission
center_name	CCME	CCME	name of the center preparing the submission
submission_comment	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Free-text comments regarding submission
lab_name	Knight	Knight	"name of lab preparing submission, can differ from center (usually refers to the PI\'s info, not the sequencing center\'s)"
submission_date	2009-10-22T01:23:00-05:00	2009-10-22T01:23:00-05:00	timestamp of submission
CONTACT	Rob Knight;Rob.Knight@Colorado.edu	Rob Knight;Rob.Knight@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
CONTACT	Noah Fierer;Noah.Fierer@Colorado.edu	Noah Fierer;Noah.Fierer@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
study	study.xml	fierer_hand_study.study.xml	"leave blank if not submitting study, put in filename otherwise"
sample	sample.xml	fierer_hand_study.sample.xml	"leave blank if not submitting sample, put in filename otherwise"
experiment		fierer_hand_study.experiment.xml	"leave blank if not submitting experiment, put in filename otherwise"
run		fierer_hand_study.run.xml	"leave blank if not submitting run, put in filename otherwise"
file	%s	fierer_hand_study.seqs.tgz	"leave blank if not submitting sequence data, put in filename otherwise"'''

submission_with_file_xml = '''<?xml version="1.0" encoding="UTF-8"?>
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
 <FILE filename="%s" checksum_method="MD5" checksum="d41d8cd98f00b204e9800998ecf8427e"/>
 </FILES>
</SUBMISSION>
'''

study_txt = '''#Field	Value	Example	Comments
STUDY_alias	fierer_hand_study	fierer_handstudy	One study per publication: this is used as an id to link files
STUDY_TITLE	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	Expected (or actual) title of the paper
STUDY_TYPE	Metagenomics	Metagenomics	"Should be ""Metagenomics"" for 16S surveys"
STUDY_ABSTRACT	"Short \'abstract\' with special characters <10%."	"Abstract, e.g. of the publication"
STUDY_DESCRIPTION	Targeted Gene Survey from Human Skin	Targeted Gene Survey from Human Skin	"Use ""Targeted Gene Survey"" for 16S or other target gene studies"
CENTER_NAME	CCME	CCME	"NCBI-approved name of sequencing center, e.g. WUGSC"
CENTER_PROJECT_NAME	NULL	NULL	"Name of project as used by the sequencing center, NULL if none."
PROJECT_ID	34527	34527	"Project ID, assigned by SRA, leave blank if not yet assigned."
PMID	19004758	19004758	"PubMed ID of paper describing project, if supplied will write out STUDY_LINK block, can be multiple (comma-delimited)"
'''


study_template = '''<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <STUDY alias="%(STUDY_alias)s">
    <DESCRIPTOR>
        <STUDY_TITLE>%(STUDY_TITLE)s</STUDY_TITLE>
        <STUDY_TYPE existing_study_type="%(STUDY_TYPE)s"/>
        <STUDY_ABSTRACT>%(STUDY_ABSTRACT)s</STUDY_ABSTRACT>
        <STUDY_DESCRIPTION>%(STUDY_DESCRIPTION)s</STUDY_DESCRIPTION>
        <CENTER_NAME>%(CENTER_NAME)s</CENTER_NAME>
        <CENTER_PROJECT_NAME>%(CENTER_PROJECT_NAME)s</CENTER_PROJECT_NAME>
        <PROJECT_ID>%(PROJECT_ID)s</PROJECT_ID>
    </DESCRIPTOR>%(XML_STUDY_LINKS_BLOCK)s
  </STUDY>
</STUDY_SET>
'''

study_xml = '''<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <STUDY alias="fierer_hand_study">
    <DESCRIPTOR>
        <STUDY_TITLE>The influence of sex, handedness, and washing on the diversity of hand surface bacteria</STUDY_TITLE>
        <STUDY_TYPE existing_study_type="Metagenomics"/>
        <STUDY_ABSTRACT>Short &apos;abstract&apos; with special characters &lt;10%.</STUDY_ABSTRACT>
        <STUDY_DESCRIPTION>Targeted Gene Survey from Human Skin</STUDY_DESCRIPTION>
        <CENTER_NAME>CCME</CENTER_NAME>
        <CENTER_PROJECT_NAME>NULL</CENTER_PROJECT_NAME>
        <PROJECT_ID>34527</PROJECT_ID>
    </DESCRIPTOR>
    <STUDY_LINKS>
      <STUDY_LINK>
        <ENTREZ_LINK>
         <DB>pubmed</DB>
         <ID>19004758</ID>
        </ENTREZ_LINK>
      </STUDY_LINK>
    </STUDY_LINKS>
  </STUDY>
</STUDY_SET>
'''

study_pmid_empty_txt = '''#Field	Value	Example	Comments
STUDY_alias	fierer_hand_study	fierer_handstudy	One study per publication: this is used as an id to link files
STUDY_TITLE	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	Expected (or actual) title of the paper
STUDY_TYPE	Metagenomics	Metagenomics	"Should be ""Metagenomics"" for 16S surveys"
STUDY_ABSTRACT	"Short \'abstract\' with special characters <10%."	"Abstract, e.g. of the publication"
STUDY_DESCRIPTION	Targeted Gene Survey from Human Skin	Targeted Gene Survey from Human Skin	"Use ""Targeted Gene Survey"" for 16S or other target gene studies"
CENTER_NAME	CCME	CCME	"NCBI-approved name of sequencing center, e.g. WUGSC"
CENTER_PROJECT_NAME	NULL	NULL	"Name of project as used by the sequencing center, NULL if none."
PROJECT_ID	34527	34527	"Project ID, assigned by SRA, leave blank if not yet assigned."
PMID		19004758	"PubMed ID of paper describing project, if supplied will write out STUDY_LINK block, can be multiple (comma-delimited)"
'''


study_pmid_missing_xml = '''<?xml version="1.0" encoding="UTF-8"?>
<STUDY_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <STUDY alias="fierer_hand_study">
    <DESCRIPTOR>
        <STUDY_TITLE>The influence of sex, handedness, and washing on the diversity of hand surface bacteria</STUDY_TITLE>
        <STUDY_TYPE existing_study_type="Metagenomics"/>
        <STUDY_ABSTRACT>Short &apos;abstract&apos; with special characters &lt;10%.</STUDY_ABSTRACT>
        <STUDY_DESCRIPTION>Targeted Gene Survey from Human Skin</STUDY_DESCRIPTION>
        <CENTER_NAME>CCME</CENTER_NAME>
        <CENTER_PROJECT_NAME>NULL</CENTER_PROJECT_NAME>
        <PROJECT_ID>34527</PROJECT_ID>
    </DESCRIPTOR>
  </STUDY>
</STUDY_SET>
'''

study_pmid_missing_txt = '''#Field	Value	Example	Comments
STUDY_alias	fierer_hand_study	fierer_handstudy	One study per publication: this is used as an id to link files
STUDY_TITLE	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	Expected (or actual) title of the paper
STUDY_TYPE	Metagenomics	Metagenomics	"Should be ""Metagenomics"" for 16S surveys"
STUDY_ABSTRACT	"Short \'abstract\' with special characters <10%."	"Abstract, e.g. of the publication"
STUDY_DESCRIPTION	Targeted Gene Survey from Human Skin	Targeted Gene Survey from Human Skin	"Use ""Targeted Gene Survey"" for 16S or other target gene studies"
CENTER_NAME	CCME	CCME	"NCBI-approved name of sequencing center, e.g. WUGSC"
CENTER_PROJECT_NAME	NULL	NULL	"Name of project as used by the sequencing center, NULL if none."
PROJECT_ID	34527	34527	"Project ID, assigned by SRA, leave blank if not yet assigned."
'''

submission_txt = '''#Field	Value	Example	Comments
accession	SRA003492	SRA003492	"leave blank if not assigned yet, e.g. if new submission"
submission_id	fierer_hand_study	fierer_hand_study	internally unique id for the submission
center_name	CCME	CCME	name of the center preparing the submission
submission_comment	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Free-text comments regarding submission
lab_name	Knight	Knight	"name of lab preparing submission, can differ from center (usually refers to the PI\'s info, not the sequencing center\'s)"
submission_date	2009-10-22T01:23:00-05:00	2009-10-22T01:23:00-05:00	timestamp of submission
CONTACT	Rob Knight;Rob.Knight@Colorado.edu	Rob Knight;Rob.Knight@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
CONTACT	Noah Fierer;Noah.Fierer@Colorado.edu	Noah Fierer;Noah.Fierer@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
'''

submission_template = '''<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"%(ACCESSION_STRING)s
 submission_id="%(submission_id)s"
 center_name="%(center_name)s"
 submission_comment="%(submission_comment)s"
 lab_name="%(lab_name)s"
 submission_date="%(submission_date)s"
>%(XML_CONTACT_BLOCK)s%(XML_ACTION_BLOCK)s%(XML_FILE_BLOCK)s
</SUBMISSION>
'''

submission_xml = '''<?xml version="1.0" encoding="UTF-8"?>
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
</SUBMISSION>
'''

sample_txt = '''#SAMPLE_ALIAS	TITLE	TAXON_ID	COMMON_NAME	ANONYMIZED_NAME	DESCRIPTION	host_taxon_id	subject	sex	hand	age	palm size	dominant hand	hours since wash
fierer_hand_study_default	human hand microbiome	539655	human skin metagenome		"Human palm microbiome, default sample for unclassified reads"								
S1	human hand microbiome	539655	human skin metagenome	subject 1	female right palm	9606	1	female	right	18	9.5	right	less than 2
S2	human hand microbiome	539655	human skin metagenome	subject 1	female left palm	9606	1	female	left	18	9.5	right	less than 2
'''

sample_template = '''<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
%(XML_SAMPLE_BLOCK)s
</SAMPLE_SET>
'''

sample_xml = '''<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <SAMPLE alias="fierer_hand_study_default">
    <TITLE>human hand microbiome</TITLE>
    <SAMPLE_NAME>
      <COMMON_NAME>human skin metagenome</COMMON_NAME>
      <TAXON_ID>539655</TAXON_ID>
    </SAMPLE_NAME>
    <DESCRIPTION>Human palm microbiome, default sample for unclassified reads</DESCRIPTION>
  </SAMPLE>
  <SAMPLE alias="S1">
    <TITLE>human hand microbiome</TITLE>
    <SAMPLE_NAME>
      <ANONYMIZED_NAME>subject 1</ANONYMIZED_NAME>
      <COMMON_NAME>human skin metagenome</COMMON_NAME>
      <TAXON_ID>539655</TAXON_ID>
    </SAMPLE_NAME>
    <DESCRIPTION>female right palm</DESCRIPTION>
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE> <TAG>age</TAG> <VALUE>18</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>dominant hand</TAG> <VALUE>right</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>hand</TAG> <VALUE>right</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>host_taxon_id</TAG> <VALUE>9606</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>hours since wash</TAG> <VALUE>less than 2</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>palm size</TAG> <VALUE>9.5</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>sex</TAG> <VALUE>female</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>subject</TAG> <VALUE>1</VALUE> </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
  <SAMPLE alias="S2">
    <TITLE>human hand microbiome</TITLE>
    <SAMPLE_NAME>
      <ANONYMIZED_NAME>subject 1</ANONYMIZED_NAME>
      <COMMON_NAME>human skin metagenome</COMMON_NAME>
      <TAXON_ID>539655</TAXON_ID>
    </SAMPLE_NAME>
    <DESCRIPTION>female left palm</DESCRIPTION>
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE> <TAG>age</TAG> <VALUE>18</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>dominant hand</TAG> <VALUE>right</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>hand</TAG> <VALUE>left</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>host_taxon_id</TAG> <VALUE>9606</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>hours since wash</TAG> <VALUE>less than 2</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>palm size</TAG> <VALUE>9.5</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>sex</TAG> <VALUE>female</VALUE> </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE> <TAG>subject</TAG> <VALUE>1</VALUE> </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</SAMPLE_SET>
'''


if __name__ == '__main__':
    main()

