#!/usr/bin/env python
from cStringIO import StringIO
import os
import re
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files
from qiime.make_sra_submission import (
    detect_missing_experiment_fields, detect_missing_study_fields,
    detect_missing_submission_fields, detect_missing_sample_fields,
    md5_path, safe_for_xml, read_tabular_data,
    make_study, make_submission,
    make_sample, write_xml_generic,
    make_run, make_experiment, threecol_data_to_dict,
    _experiment_link_xml, _experiment_attribute_xml, pretty_xml,
    generate_output_fp, #update_entry_with_derived_fields,
    SraEntity, SraRun, SraDataBlock, SraFile, SraPoolMember,
    SraSpotDescriptor, SraSampleDescriptor, SraExperimentSet, SraExperiment,
    SraStudyRef, SraPlatform, #update_entry_with_deprecated_fields,
    SraSample, SraSampleName, SraSampleAttributes, SraSampleSet, SraStudy,
    SraContacts, SraActions, SraSubmissionFiles, SraSubmission,
    SraSubmissionTable, SraInputTable, SraExperimentTable,
    )
from qiime.util import get_qiime_project_dir
import xml.etree.ElementTree as ET
from cStringIO import StringIO

"""Tests of the make_study_and_experiment.py file.

Note: these tests assume you are running from within the directory that also
has the templates and sample tabular data.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Rohini Sinha"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Development"


def standardize_xml(xml_str):
    """Standardize an XML string with the ElementTree package."""
    return ET.tostring(ET.XML(xml_str))

class TopLevelTests(TestCase):
    """Top-level tests of functions in make_study_and_experiment"""
    def setUp(self):
        """ """
        self.files_to_remove = []

        # Cannot use get_qiime_project_dir() due to test errors in virtual box
        test_dir = os.path.dirname(os.path.abspath(__file__))
        self.sff_dir = os.path.join(test_dir, 'sra_test_files', 'F6AVWTA')

    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_detect_missing_experiment_fields(self):
        """detect_missing_experiment_fields should return a list of all required fields not found in the input file header."""
        input_file = StringIO('#NONSENSE_FIELD\nnonsense_value\n')
        observed = detect_missing_experiment_fields(input_file)
        observed.sort()
        expected = [
            'BARCODE',
            'EXPERIMENT_CENTER',
            'EXPERIMENT_DESIGN_DESCRIPTION',
            'EXPERIMENT_TITLE',
            'LIBRARY_CONSTRUCTION_PROTOCOL',
            'POOL_PROPORTION',
            'RUN_PREFIX',
            'SAMPLE_ALIAS',
            'STUDY_REF',
            ]
        self.assertEqual(observed, expected)
        
        # EXPERIMENT_TITLE is present, so doesn't show up in missing fields
        input_file = StringIO('#NONSENSE_FIELD\tEXPERIMENT_TITLE\nnonsense_value\tfake_title\n')
        observed = detect_missing_experiment_fields(input_file)
        observed.sort()
        expected = [
            'BARCODE',
            'EXPERIMENT_CENTER',
            'EXPERIMENT_DESIGN_DESCRIPTION',
            'LIBRARY_CONSTRUCTION_PROTOCOL',
            'POOL_PROPORTION',
            'RUN_PREFIX',
            'SAMPLE_ALIAS',
            'STUDY_REF',
            ]
        self.assertEqual(observed, expected)
        
        # case of EXPERIMENT_TITLE in headers doesn't matter
        input_file = StringIO('#NONSENSE_FIELD\tExPeRiMeNt_TiTlE\nnonsense_value\tfake_title\n')
        observed = detect_missing_experiment_fields(input_file)
        observed.sort()
        expected = [
            'BARCODE',
            'EXPERIMENT_CENTER',
            'EXPERIMENT_DESIGN_DESCRIPTION',
            'LIBRARY_CONSTRUCTION_PROTOCOL',
            'POOL_PROPORTION',
            'RUN_PREFIX',
            'SAMPLE_ALIAS',
            'STUDY_REF',
            ]
        self.assertEqual(observed, expected)

    def test_detect_missing_study_fields(self):
        """detect_missing_study_fields should return a list of all required fields not found in the input file header."""
        input_file = StringIO('#NONSENSE_FIELD\nnonsense_value\n')
        observed = detect_missing_study_fields(input_file)
        observed.sort()
        expected = [
            'CENTER_NAME',
            'CENTER_PROJECT_NAME',
            'STUDY_ABSTRACT',
            'STUDY_ALIAS',
            'STUDY_DESCRIPTION',
            'STUDY_TITLE',
            'STUDY_TYPE',
            ]
        self.assertEqual(observed, expected)

    def test_detect_missing_submission_fields(self):
        """detect_missing_submission_fields should return a list of all required fields not found in the input file header."""
        input_file = StringIO('#NONSENSE_FIELD\nnonsense_value\n')
        observed = detect_missing_submission_fields(input_file)
        observed.sort()
        expected = [
            'CENTER_NAME',
            'CONTACT',
            'LAB_NAME',
            'SUBMISSION_DATE',
            'SUBMISSION_ID',
            ]
        self.assertEqual(observed, expected)

    def test_detect_missing_sample_fields(self):
        """detect_missing_sample_fields should return a list of all required fields not found in the input file header."""
        input_file = StringIO('#NONSENSE_FIELD\nnonsense_value\n')
        observed = detect_missing_sample_fields(input_file)
        observed.sort()
        expected = [
            'SAMPLE_ALIAS',
            'TAXON_ID',
            'TITLE',
            ]
        self.assertEqual(observed, expected)

    def test_generate_output_fp(self):
        input_fp = os.path.join('tmp', 'hello.txt')
        self.assertEqual(
            generate_output_fp(input_fp, '.xml'),
            os.path.join('tmp', 'hello.xml'))

        output_dir = os.path.join('home', 'user1')
        self.assertEqual(
            generate_output_fp(input_fp, '.xml', output_dir),
            os.path.join('home', 'user1', 'hello.xml'))

    def test_md5_path(self):
        """md5_path should match hand-calculated value"""
        template_fp = tempfile.mktemp(suffix='.xml')
        open(template_fp, 'wb').write('ABCDEFGHIJKLMNOP')
        result = md5_path(template_fp)
        self.assertEqual(result, '19fc8eff82037f1fc0d8ea1d32b5e339')

        self.files_to_remove.append(template_fp)

    def test_safe_for_xml(self):
        """safe_for_xml should return string without 'bad' xml entities."""
        s = 'ab&c\'d"e<f>'
        t = safe_for_xml(s)
        self.assertEqual(t, 'ab&amp;c&apos;d&quot;e&lt;f&gt;')

    def test_read_tabular_data(self):
        """read_tabular_data should read simple table"""
        data = [
            'a\tb\tc d',
            'x\ty\tz',
            '',
            'x x\ty y\tc c ',
            ]
        expected = (
            ['a', 'b', 'c d'],
            [['x','y','z'], ['x x', 'y y', 'c c ']],
            )
        self.assertEqual(read_tabular_data(data), expected)

        # Test with pound sign in header
        data[0] = '#' + data[0]
        self.assertEqual(read_tabular_data(data), expected)

    def test_make_study(self):
        """make_study should produce expected results given info/template"""
        study_sample_data = StringIO(study_twocol_txt)
        result = make_study(study_sample_data)
        self.assertEqual(standardize_xml(result), standardize_xml(study_xml))

        # Test when pmid field is empty
        study_sample_data = StringIO(study_pmid_empty_txt)
        result = make_study(study_sample_data)
        self.assertEqual(standardize_xml(result), standardize_xml(study_pmid_missing_xml))

        # Test when pmid field is missing
        study_sample_data = StringIO(study_pmid_missing_txt)
        result = make_study(study_sample_data)
        self.assertEqual(standardize_xml(result), standardize_xml(study_pmid_missing_xml))

        # Test for row-based data
        study_sample_data = StringIO(study_manycol_txt)
        result = make_study(study_sample_data, twocol_input_format=False)
        self.assertEqual(standardize_xml(result), standardize_xml(study_xml))        

    def test_make_submission(self):
        """make_submission should produce expected results given info/template"""
        submission_sample_data = StringIO(submission_twocol_txt)
        my_docnames = {'study':'study.xml', 'sample':'sample.xml'}
        result = make_submission(submission_sample_data, docnames=my_docnames)
        self.assertEqual(standardize_xml(result), standardize_xml(submission_xml))

        # Test for row-based data
        submission_sample_data = StringIO(submission_manycol_txt)
        result = make_submission(submission_sample_data, docnames=my_docnames,
                                 twocol_input_format=False)
        self.assertEqual(standardize_xml(result), standardize_xml(submission_xml))

        # Test when data file is provided
        fake_tgz_file = tempfile.NamedTemporaryFile(suffix='.tgz')
        submission_sample_data = StringIO(
            submission_with_file_txt % fake_tgz_file.name)
        result = make_submission(submission_sample_data, docnames=my_docnames)
        self.assertEqual(standardize_xml(result), standardize_xml(
            submission_with_file_xml  % fake_tgz_file.name))

    def test_make_sample(self):
        """make_sample should produce expected reuslts given info/template"""
        sample_data = StringIO(sample_txt)
        result = make_sample(sample_data)
        self.assertEqual(standardize_xml(result), standardize_xml(sample_xml))

        # test with lowercase field names
        sample_data = StringIO(sample_with_lowercase_txt)
        result = make_sample(sample_data)
        self.assertEqual(standardize_xml(result), standardize_xml(sample_xml))

    def test_write_xml_generic(self):
        input_file = tempfile.NamedTemporaryFile()
        input_file.write('abc')
        input_file.seek(0)
        
        def simple_xml_f(lines, template=None):
            return '<xml>%s</xml>' % ''.join(lines)
        
        observed_fp = write_xml_generic(
            input_file.name, simple_xml_f)
        observed = open(observed_fp).read()
        self.assertEqual(
            standardize_xml(observed),
            standardize_xml('<xml>abc</xml>')
            )
        self.files_to_remove.append(observed_fp)

        observed_fp2 = write_xml_generic(
            input_file.name, simple_xml_f, output_suffix='testsuffix')
        self.assertTrue(observed_fp2.endswith('testsuffix.xml'))
        observed = open(observed_fp2).read()
        self.assertEqual(
            standardize_xml(observed),
            standardize_xml('<xml>abc</xml>')
            )
        self.files_to_remove.append(observed_fp2)



    def test_make_run(self):
        observed_run_xml = make_run(StringIO(experiment), self.sff_dir)
        self.assertEqual(
            standardize_xml(observed_run_xml),
            standardize_xml(run_xml_str),
            )

    def test_make_experiment(self):
        observed_exp_xml = make_experiment(
            StringIO(experiment),
            attribute_file=StringIO(attrs),
            link_file=StringIO(links),
            )
        self.assertEqual(
            standardize_xml(observed_exp_xml),
            standardize_xml(experiment_xml_str),
            )

    def test_experiment_with_accession_xml(self):
        """make_run_and_experiment should return correct XML for experiment with accession numbers."""
        observed_exp_xml = make_experiment(
            StringIO(experiment_with_accessions),
            attribute_file=StringIO(attrs),
            link_file=StringIO(links),
            )
        self.assertEqual(
            standardize_xml(observed_exp_xml),
            standardize_xml(experiment_with_accessions_xml_str),
            )

    def test_experiment_with_default_sample_xml(self):
        """make_experiment should return correct XML for experiment with user-specified default sample."""
        observed_exp_xml = make_experiment(
            StringIO(experiment_with_default_sample),
            attribute_file=StringIO(attrs),
            link_file=StringIO(links),
            )
        self.assertEqual(
            standardize_xml(observed_exp_xml),
            standardize_xml(experiment_with_default_sample_xml_str),
            )

    def test_metagenomic_experiment_xml(self):
        """make_experiment should return correct XML for metagenomic experiment."""
        observed_exp_xml = make_experiment(
             StringIO(metagenomic_experiment))
        self.assertEqual(
            standardize_xml(observed_exp_xml),
            standardize_xml(metagenomic_experiment_xml_str),
            )

    def test_validate_experiment_xml(self):
        """make_experiment should produce valid XML according to official schema"""
        try:
            import lxml.etree
        except:
            return None

        observed_exp_xml = make_experiment(
            StringIO(experiment),
            attribute_file=StringIO(attrs),
            link_file=StringIO(links),
            )
        experiment_schema = lxml.etree.XMLSchema(
            lxml.etree.fromstring(experiment_xsd))
        self.assertTrue(experiment_schema.validate(
            lxml.etree.fromstring(observed_exp_xml)))

    def test_validate_run_xml(self):
        """make_run should produce valid XML according to official schema"""
        try:
            import lxml.etree
        except:
            return None

        observed_run_xml = make_run(
            StringIO(experiment), self.sff_dir)
        run_schema = lxml.etree.XMLSchema(lxml.etree.fromstring(run_xsd))
        self.assertTrue(run_schema.validate(
            lxml.etree.fromstring(observed_run_xml)))

    def test_threecol_data_to_dict(self):
        """threecol_data_to_dict should parse experiment links table correctly"""
        header, body = read_tabular_data(StringIO(links))
        observed = threecol_data_to_dict(body)
        expected = {
            'bodysites_F6AVWTA02': [(
                'bodysites Library Construction Protocol',
                'http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf'
                )],
            'bodysites_F6AVWTA01': [(
                'bodysites Library Construction Protocol',
                'http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf'
                )],
            }
        self.assertEqual(observed.items(), expected.items())

    def test_experiment_link_xml(self):
        """experiment_link_xml should generate correct XML fragment for list of links."""
        links = [('link1', 'http://google.com'),
                 ('links2', 'http://www.ncbi.nlm.nih.gov')]
        expected = '''
        <EXPERIMENT_LINKS>
          <EXPERIMENT_LINK>
            <URL_LINK>
              <LABEL>link1</LABEL>
              <URL>http://google.com</URL>
            </URL_LINK>
          </EXPERIMENT_LINK>
          <EXPERIMENT_LINK>
            <URL_LINK>
              <LABEL>links2</LABEL>
              <URL>http://www.ncbi.nlm.nih.gov</URL>
            </URL_LINK>
          </EXPERIMENT_LINK>
        </EXPERIMENT_LINKS>'''
        observed = _experiment_link_xml(links)
        self.assertEqual(pretty_xml(observed, 4), expected)

    def test_experiment_attribute_xml(self):
        """experiment_attribute_xml should generate correct XML fragment for list of attributes."""
        attrs = [('a1', 'val1'),
                 ('attr2', 'something else')]
        expected = '''
        <EXPERIMENT_ATTRIBUTES>
          <EXPERIMENT_ATTRIBUTE>
            <TAG>a1</TAG>
            <VALUE>val1</VALUE>
          </EXPERIMENT_ATTRIBUTE>
          <EXPERIMENT_ATTRIBUTE>
            <TAG>attr2</TAG>
            <VALUE>something else</VALUE>
          </EXPERIMENT_ATTRIBUTE>
        </EXPERIMENT_ATTRIBUTES>'''
        observed = _experiment_attribute_xml(attrs)
        self.assertEqual(pretty_xml(observed, 4), expected)

        observed = _experiment_attribute_xml([])
        self.assertEqual(pretty_xml(observed, 4), '')

    def test_pretty_xml(self):
        """pretty_xml should modify XML ElementTree to include correct indentation."""
        a = ET.Element('a')
        b = ET.SubElement(a, 'b')
        b.text = 'hello'
        observed = pretty_xml(a)
        expected = '\n<a>\n  <b>hello</b>\n</a>'
        self.assertEqual(observed, expected)

        c = ET.Element('c')
        c.set('myattr', '12345')
        observed = pretty_xml(c)
        expected = '\n<c myattr="12345" />'
        self.assertEqual(observed, expected)

        c = ET.Element('c')
        c.text = '\n'
        observed = pretty_xml(c, 1)
        expected = '\n  <c>\n</c>'
        self.assertEqual(observed, expected)

    ## def test_update_entry_with_deprecated_fields(self):
    ##     """update_entry_with_deprecated_fields should move deprecated field values to new fields."""
    ##     a = {'SAMPLE_ACCESSION': '12345'}
    ##     expected = {'SAMPLE_ACCESSION': '12345',
    ##                 'DEFAULT_SAMPLE_ACCESSION': '12345'}
    ##     update_entry_with_deprecated_fields(a, warn=False)
    ##     self.assertEqual(a, expected)

    ## def test_update_entry_with_derived_fields(self):
    ##     """update_entry_with_derived_fields should derive correct values for optional fields."""
    ##     a = {
    ##         'STUDY_REF': 'mystudy',
    ##         'RUN_PREFIX': 'RUN01',
    ##         'SAMPLE_ALIAS': 'mysample',
    ##         'BARCODE': 'AAAAGGGG',
    ##         'PRIMER': 'ATTACCGCGGCTGCTGGC',
    ##         'EXPERIMENT_CENTER': 'NCBI',
    ##         }
    ##     update_entry_with_derived_fields(a)
    ##     self.assertEqual(a['EXPERIMENT_ALIAS'], 'mystudy_RUN01')
    ##     self.assertEqual(a['RUN_ALIAS'], 'mystudy_mysample_RUN01')
    ##     self.assertEqual(a['BARCODE_READ_GROUP_TAG'], 'RUN01_AAAAGGGG')
    ##     self.assertEqual(a['PRIMER_READ_GROUP_TAG'], 'V1-V3')
    ##     self.assertEqual(a['POOL_MEMBER_NAME'], 'RUN01_mysample_V1-V3')
    ##     self.assertEqual(a['POOL_MEMBER_FILENAME'], 'RUN01_mysample_V1-V3.sff')
    ##     self.assertEqual(a['DEFAULT_SAMPLE_CENTER'], 'NCBI')
    ##     self.assertEqual(a['DEFAULT_SAMPLE_NAME'], 'mystudy_default')
    ##     self.assertEqual(a['DEFAULT_SAMPLE_FILENAME'], 'mystudy_default_RUN01.sff')
    ##     self.assertEqual(a['DEFAULT_RUN_ALIAS'], 'mystudy_default_RUN01')
    ##     self.assertEqual(a['LIBRARY_STRATEGY'], 'AMPLICON')
    ##     self.assertEqual(a['LIBRARY_SOURCE'], 'GENOMIC')
    ##     self.assertEqual(a['LIBRARY_SELECTION'], 'PCR')

class SraInputTableTests(TestCase):
    def setUp(self):
        self.input_file = StringIO(submission_manycol_txt)
        self.twocol_input_file = StringIO(submission_twocol_txt)
        self.expected_header = [
            'ACCESSION',
            'SUBMISSION_ID',
            'CENTER_NAME',
            'SUBMISSION_COMMENT',
            'LAB_NAME',
            'SUBMISSION_DATE',
            'CONTACT',
            ]
        self.expected_rows = [[ 
            'SRA003492',
            'fierer_hand_study',
            'CCME',
            'Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov',
            'Knight',
            '2009-10-22T01:23:00-05:00',
            'Rob Knight;Rob.Knight@Colorado.edu,Noah Fierer;Noah.Fierer@Colorado.edu',
            ]]
    
    def test_entries(self):
        t = SraInputTable(['Col1', 'Col2'], [
            ['a', 'b'],
            ['c', 'd'],
            ])
        observed = list(t.entries)
        expected = [
            {'Col1': 'a', 'Col2': 'b'},
            {'Col1': 'c', 'Col2': 'd'},
            ]
        self.assertEqual(observed, expected)

    def test_parse(self):
        t = SraInputTable.parse(self.input_file)
        self.assertEqual(t.header, self.expected_header)
        self.assertEqual(t.rows, self.expected_rows)

    def test_parse_twocol_format(self):
        t = SraInputTable.parse_twocol_format(self.twocol_input_file)
        self.assertEqual(t.header, self.expected_header)
        self.assertEqual(t.rows, self.expected_rows)
        
        multiple_data = [
            '#ignored\n',
            'G\ty\tz\n',
            'G x\ty y\tc c\n',
            'G\t1\tzzz\n',
            'G\tc\tcomment\n',
            ]
        t2 = SraInputTable.parse_twocol_format(multiple_data)
        self.assertEqual(t2.header, ['G', 'G X'])
        self.assertEqual(t2.rows, [['y,1,c', 'y y']])

    def test_get_entry(self):
        """_get_entry should return dict from list of row values"""
        t = SraInputTable(['a', 'b', 'c d'], None)
        self.assertEqual(
            t._get_entry(['aa', 'bb', 'cc']),
            {'a': 'aa', 'b': 'bb', 'c d': 'cc'},
            )
        # missing second field
        self.assertEqual(
            t._get_entry(['x', '', 'z']),
            {'a': 'x', 'c d': 'z'},
            )
        # missing final field
        self.assertEqual(
            t._get_entry(['x x', 'y y']),
            {'a': 'x x', 'b': 'y y'},
            )


class SraSubmissionTableTests(TestCase):
    def test_derive_optional_fields(self):
        t = SraSubmissionTable([], [[]])
        t.derive_optional_fields()
        entry = t.entries.next()
        observed_datetime = entry.get('SUBMISSION_DATE')
        # Regular expression constructed from w3 XML documentation
        datetime_regex = re.compile(
            '-?\d\d\d\d-\d\d-\d\dT\d\d:\d\d:\d\d'
            '(\.\d+)?'
            '(Z|([+-]\d\d:\d\d))?'
            '$'
            )
        self.assertTrue(datetime_regex.match(observed_datetime))

        # TODO: Test with file, check that correct checksum is derived

 

class SraEntityTests(TestCase):
    def test_check_attributes(self):
        """SraEntity.__init__ should raise ValueError if unknown attributes are provided."""
        self.assertRaises(ValueError, SraEntity, {'a': 'b'})

    def test_attribute_names(self):
        """SraEntity.attribute_names should list all known attributes."""
        self.assertEqual(SraEntity({}).attribute_names, [])

    def test_from_entry(self):
        """SraEntity.from_entry should return new instance."""
        a = SraEntity.from_entry({})
        self.assertEqual(a.__class__, SraEntity)

    def test_set_xml_attributes(self):
        """SraEntity.set_xml_attributes should transfer instance attributes to an ETree XML element."""
        element = ET.Element('a')
        entity = SraEntity({})
        entity.name = 'John'
        entity.set_xml_attributes(element, ['name'])
        self.assertEqual(element.get('name'), 'John')


class SraSubmissionTests(TestCase):
    def test_with_date(self):
        class mocktable:
            first_entry = {
                'SUBMISSION_ID': 'fierer_hand_study',
                'CENTER_NAME': 'CCME',
                'LAB_NAME': 'Knight',
                'SUBMISSION_DATE': '2009-10-22T01:23:00-05:00',
                'CONTACT': (
                    'Rob Knight;Rob.Knight@Colorado.edu,'
                    'Noah Fierer;Noah.Fierer@Colorado.edu'),
                }
        s = SraSubmission.from_table(mocktable)
        expected = '''<?xml version="1.0" encoding="UTF-8"?>
        <SUBMISSION center_name="CCME" lab_name="Knight" submission_date="2009-10-22T01:23:00-05:00" submission_id="fierer_hand_study" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
          <CONTACTS>
            <CONTACT inform_on_error="Rob.Knight@Colorado.edu" inform_on_status="Rob.Knight@Colorado.edu" name="Rob Knight" />
            <CONTACT inform_on_error="Noah.Fierer@Colorado.edu" inform_on_status="Noah.Fierer@Colorado.edu" name="Noah Fierer" />
          </CONTACTS>
        </SUBMISSION>'''
        self.assertEqual(s.to_xml_string(level=4), expected)

    def test_without_date(self):
        t = SraSubmissionTable(
            ['SUBMISSION_ID', 'CENTER_NAME', 'LAB_NAME', 'CONTACT'],
            [['test', 'ABCD', 'Lab1', 'Bob;bob@abcd.com']],
            )
        t.derive_optional_fields()
        date_str = t.entries.next()['SUBMISSION_DATE']
        expected = '''
        <SUBMISSION center_name="ABCD" lab_name="Lab1" submission_date="%s" submission_id="test" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
          <CONTACTS>
            <CONTACT inform_on_error="bob@abcd.com" inform_on_status="bob@abcd.com" name="Bob" />
          </CONTACTS>
        </SUBMISSION>'''
        submission = SraSubmission.from_table(t)
        self.assertEqual(pretty_xml(submission.to_xml(), 4), expected % date_str)


class SraSubmissionFilesTests(TestCase):
    def test_to_xml(self):
        entry = {
            'FILE': 'my_submission.tgz',
            'SUBMISSION_FILE_CHECKSUM': 'abcdefg12345',
            }
        f = SraSubmissionFiles.from_entry(entry)
        expected = '''
        <FILES>
          <FILE checksum="abcdefg12345" checksum_method="MD5" filename="my_submission.tgz" />
        </FILES>'''
        self.assertEqual(pretty_xml(f.to_xml(), 4), expected)


class SraActionsTests(TestCase):
    def test_to_xml(self):
        documents = {'study':'study.xml', 'sample':'sample.xml'}
        a = SraActions.from_entry(documents)
        expected = '''
        <ACTIONS>
          <ACTION>
            <ADD notes="study metadata" schema="study" source="study.xml" />
          </ACTION>
          <ACTION>
            <ADD notes="sample metadata" schema="sample" source="sample.xml" />
          </ACTION>
          <ACTION>
            <RELEASE />
          </ACTION>
        </ACTIONS>'''
        self.assertEqual(pretty_xml(a.to_xml(), 4), expected)


class SraContactsTests(TestCase):
    def test_parse_contacts(self):
        observed = SraContacts._parse_contacts(
            'jane;smith@college.edu, john; doe@company.com ')
        expected = [
            ('jane', 'smith@college.edu'),
            ('john', 'doe@company.com'),
            ]
        self.assertEqual(observed, expected)

    def test_to_xml(self):
        entry = {
            'CONTACT': 'Rob Knight;Rob.Knight@Colorado.edu, Noah Fierer;Noah.Fierer@Colorado.edu'
            }
        c = SraContacts.from_entry(entry)
        expected = '''
        <CONTACTS>
          <CONTACT inform_on_error="Rob.Knight@Colorado.edu" inform_on_status="Rob.Knight@Colorado.edu" name="Rob Knight" />
          <CONTACT inform_on_error="Noah.Fierer@Colorado.edu" inform_on_status="Noah.Fierer@Colorado.edu" name="Noah Fierer" />
        </CONTACTS>'''
        self.assertEqual(pretty_xml(c.to_xml(), 4), expected)


class SraStudyTests(TestCase):
    def test_to_xml(self):
        """ """
        entry = {
            'STUDY_ALIAS': 'fierer_hand_study',
            'STUDY_TITLE': (
                'The influence of sex, handedness, and washing on '
                'the diversity of hand surface bacteria'),
            'STUDY_ABSTRACT': (
                'Short &apos;abstract&apos; with special characters &lt;10%.'),
            'STUDY_DESCRIPTION': 'Targeted Gene Survey from Human Skin',
            'STUDY_TYPE': 'Metagenomics',
            'CENTER_NAME': 'CCME',
            'CENTER_PROJECT_NAME': 'NULL',
            'PMID': '19004758',
            }
        s = SraStudy.from_entry(entry)
        expected = '''
        <STUDY alias="fierer_hand_study">
          <DESCRIPTOR>
            <STUDY_TITLE>The influence of sex, handedness, and washing on the diversity of hand surface bacteria</STUDY_TITLE>
            <STUDY_TYPE existing_study_type="Metagenomics" />
            <STUDY_ABSTRACT>Short &apos;abstract&apos; with special characters &lt;10%.</STUDY_ABSTRACT>
            <STUDY_DESCRIPTION>Targeted Gene Survey from Human Skin</STUDY_DESCRIPTION>
            <CENTER_NAME>CCME</CENTER_NAME>
            <CENTER_PROJECT_NAME>NULL</CENTER_PROJECT_NAME>
          </DESCRIPTOR>
          <STUDY_LINKS>
            <STUDY_LINK>
              <ENTREZ_LINK>
                <DB>pubmed</DB>
                <ID>19004758</ID>
              </ENTREZ_LINK>
            </STUDY_LINK>
          </STUDY_LINKS>
        </STUDY>'''
        self.assertEqual(pretty_xml(s.to_xml(), 4), expected)

class SraSampleSetTests(TestCase):
    def test_to_xml(self):
        """ """
        ss = SraSampleSet()
        expected = (
            '\n<SAMPLE_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '\n</SAMPLE_SET>')
        self.assertEqual(pretty_xml(ss.to_xml()), expected)


class SraSampleTests(TestCase):
    def test_to_xml(self):
        """ """
        entry = {
            'SAMPLE_ALIAS': 'fierer_hand_study_default',
            'DESCRIPTION': 'default sample for unclassified reads',
            'TITLE': 'human hand microbiome'
            }
        s = SraSample.from_entry(entry)
        expected =  '''
        <SAMPLE alias="fierer_hand_study_default">
          <TITLE>human hand microbiome</TITLE>
          <DESCRIPTION>default sample for unclassified reads</DESCRIPTION>
        </SAMPLE>'''
        self.assertEqual(pretty_xml(s.to_xml(), 4), expected)


class SraSampleNameTests(TestCase):
    def test_to_xml(self):
        """SraSamplename.to_xml should produce valid XML."""
        entry = {
            'TAXON_ID': '66545',
            'COMMON_NAME': 'aname',
            'ANONYMIZED_NAME': 'hhyot',
            }
        s = SraSampleName.from_entry(entry)
        expected = '''
        <SAMPLE_NAME>
          <ANONYMIZED_NAME>hhyot</ANONYMIZED_NAME>
          <COMMON_NAME>aname</COMMON_NAME>
          <TAXON_ID>66545</TAXON_ID>
        </SAMPLE_NAME>'''
        self.assertEqual(pretty_xml(s.to_xml(), 4), expected)


class SraSampleAttributesTests(TestCase):
    def test_to_xml(self):
        """SraSampleAttributes.to_xml should produce valid XML."""
        entry = {
            'TAXON_ID': '66545',
            'COMMON_NAME': 'aname',
            'ANONYMIZED_NAME': 'hhyot',
            'hand': 'left',
            'age': '45',
            }
        s = SraSampleAttributes.from_entry(entry)
        expected = '''
        <SAMPLE_ATTRIBUTES>
          <SAMPLE_ATTRIBUTE>
            <TAG>age</TAG>
            <VALUE>45</VALUE>
          </SAMPLE_ATTRIBUTE>
          <SAMPLE_ATTRIBUTE>
            <TAG>hand</TAG>
            <VALUE>left</VALUE>
          </SAMPLE_ATTRIBUTE>
        </SAMPLE_ATTRIBUTES>'''

        self.assertEqual(pretty_xml(s.to_xml(), 4), expected)


class SraExperimentSetTests(TestCase):
    def test_register(self):
        """SraExperimentSet.register should preserve order of experiment entries."""
        s = SraExperimentSet()
        common_fields = {
            'EXPERIMENT_CENTER': 'UPENNBL',
            'EXPERIMENT_TITLE': 'my_experiment',
            'EXPERIMENT_DESIGN_DESCRIPTION': 'A b c',
            'LIBRARY_STRATEGY': 'D e f',
            'LIBRARY_SOURCE': 'G h i',
            'LIBRARY_SELECTION': 'J k l',
            'LIBRARY_CONSTRUCTION_PROTOCOL': 'M n o',
            'PLATFORM': 'FLX',
            'KEY_SEQ': 'TCAG',
            'BARCODE': 'AAAAGGGG',
            'BARCODE_READ_GROUP_TAG': 'my_AAAAGGGG',
            'SAMPLE_ALIAS': 'P q r',
            'SAMPLE_CENTER': 'NCBI',
            'POOL_MEMBER_NAME': 'S t u',
            'POOL_PROPORTION': '0.1',
            }
        e1 = {'EXPERIMENT_ALIAS': 'experiment1'}
        e1.update(common_fields)
        s.register(e1)
        e2 = {'EXPERIMENT_ALIAS': 'experiment2'}
        e2.update(common_fields)
        s.register(e2)
        e3 = {'EXPERIMENT_ALIAS': 'experiment3'}
        e3.update(common_fields)
        s.register(e3)
        self.assertEqual([a for a, b in s.experiments],
                         ['experiment1', 'experiment2', 'experiment3'])

    def test_to_xml(self):
        """SraExperimentSet.to_xml should return blank experiment set."""
        s = SraExperimentSet()
        expected = (
            '\n<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '\n</EXPERIMENT_SET>'
            )
        self.assertEqual(pretty_xml(s.to_xml()), expected)

class SraExperimentTests(TestCase):
    def setUp(self):
        self.entry = {
            'EXPERIMENT_ALIAS': 'my_experiment',
            'EXPERIMENT_CENTER': 'UPENNBL',
            'EXPERIMENT_TITLE': 'My first experiment',
            'EXPERIMENT_DESIGN_DESCRIPTION': 'A b c',
            'LIBRARY_STRATEGY': 'D e f',
            'LIBRARY_SOURCE': 'G h i',
            'LIBRARY_SELECTION': 'J k l',
            'LIBRARY_CONSTRUCTION_PROTOCOL': 'M n o',
            'PLATFORM': 'FLX',
            'KEY_SEQ': 'TCAG',
            'BARCODE': 'AAAAGGGG',
            'BARCODE_READ_GROUP_TAG': 'my_AAAAGGGG',
            'SAMPLE_ALIAS': 'P q r',
            'SAMPLE_CENTER': 'NCBI',
            'POOL_MEMBER_NAME': 'S t u',
            'POOL_PROPORTION': '0.1',
            }
        self.experiment = SraExperiment.from_entry(self.entry)
    
    def test_from_entry(self):
        """SraExperiment.from_entry should return new instance."""
        self.assertEqual(self.experiment.__class__, SraExperiment)

        del self.entry['EXPERIMENT_ALIAS']
        self.assertRaises(KeyError, SraExperiment.from_entry, self.entry)

    def test_register(self):
        """SraExperiment.register should not raise an error for valid entry."""
        self.entry['BARCODE'] = 'CCCCTTTT'
        self.entry['BARCODE_READ_GROUP_TAG'] = 'my_CCCCTTTT'
        self.assertEqual(self.experiment.register(self.entry), None)

    def test_to_xml(self):
        """SraExperiment.to_xml should return valid XML for experiment with one entry."""
        observed = self.experiment.to_xml()
        expected = '''
        <EXPERIMENT alias="my_experiment" center_name="UPENNBL">
          <TITLE>My first experiment</TITLE>
          <STUDY_REF refcenter="NCBI" />
          <DESIGN>
            <DESIGN_DESCRIPTION>A b c</DESIGN_DESCRIPTION>
            <SAMPLE_DESCRIPTOR />
            <LIBRARY_DESCRIPTOR>
              <LIBRARY_NAME>my_experiment</LIBRARY_NAME>
              <LIBRARY_STRATEGY>D e f</LIBRARY_STRATEGY>
              <LIBRARY_SOURCE>G h i</LIBRARY_SOURCE>
              <LIBRARY_SELECTION>J k l</LIBRARY_SELECTION>
              <LIBRARY_LAYOUT>
                <SINGLE />
              </LIBRARY_LAYOUT>
              <LIBRARY_CONSTRUCTION_PROTOCOL>M n o</LIBRARY_CONSTRUCTION_PROTOCOL>
            </LIBRARY_DESCRIPTOR>
            <SPOT_DESCRIPTOR>
              <SPOT_DECODE_SPEC>
                <READ_SPEC>
                  <READ_INDEX>0</READ_INDEX>
                  <READ_CLASS>Technical Read</READ_CLASS>
                  <READ_TYPE>Adapter</READ_TYPE>
                  <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
                </READ_SPEC>
                <READ_SPEC>
                  <READ_INDEX>1</READ_INDEX>
                  <READ_CLASS>Application Read</READ_CLASS>
                  <READ_TYPE>Forward</READ_TYPE>
                  <RELATIVE_ORDER follows_read_index="0" />
                </READ_SPEC>
              </SPOT_DECODE_SPEC>
            </SPOT_DESCRIPTOR>
          </DESIGN>
          <PLATFORM>
            <LS454>
              <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
              <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
              <FLOW_COUNT>400</FLOW_COUNT>
            </LS454>
          </PLATFORM>
          <PROCESSING>
            <BASE_CALLS>
              <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
              <BASE_CALLER>454 BaseCaller</BASE_CALLER>
            </BASE_CALLS>
            <QUALITY_SCORES qtype="phred">
              <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
              <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
              <MULTIPLIER>1.0</MULTIPLIER>
            </QUALITY_SCORES>
          </PROCESSING>
        </EXPERIMENT>'''
        self.assertEqual(pretty_xml(observed, 4), expected)

class SraStudyRefTests(TestCase):
    def test_to_xml(self):
        """SraStudyRef.to_xml should return valid XML fragment."""
        a = SraStudyRef.from_entry({
            'SAMPLE_CENTER': 'mycenter',
            'STUDY_REF': 'mystudy',
            })
        expected = '''
        <STUDY_REF refcenter="mycenter" refname="mystudy" />'''
        self.assertEqual(pretty_xml(a.to_xml(), 4), expected)

        b = SraStudyRef.from_entry({
            'SAMPLE_CENTER': 'mycenter',
            'STUDY_ACCESSION': '12345'
            })
        expected = '''
        <STUDY_REF accession="12345" refcenter="mycenter" />'''
        self.assertEqual(pretty_xml(b.to_xml(), 4), expected)

class SraPlatformTests(TestCase):
    def test_from_entry(self):
        """SraPlatform.from_entry should raise ValueError if PLATFORM is unrecognized."""
        self.assertRaises(ValueError, SraPlatform.from_entry, {'PLATFORM': 'unknown'})

    def test_to_xml(self):
        """SraPlatform.to_xml should produce valid XML fragment."""
        a = SraPlatform.from_entry({'PLATFORM': 'FLX'})
        expected = '''
        <PLATFORM>
          <LS454>
            <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
            <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
            <FLOW_COUNT>400</FLOW_COUNT>
          </LS454>
        </PLATFORM>'''
        self.assertEqual(pretty_xml(a.to_xml(), 4), expected)

class SraSpotDescriptorTests(TestCase):
    def setUp(self):
        entry = {
            'KEY_SEQ': 'TACG',
            'BARCODE_READ_GROUP_TAG': 'Sample1',
            'BARCODE': 'TTAACCGG',
            'PRIMER_READ_GROUP_TAG': 'V1-V2',
            'PRIMER': 'CTGCTGCCTYCCGTA',
            }
        self.spot_descriptor = SraSpotDescriptor.from_entry(entry)
        self.spot_descriptor.register(entry)
    
    def test_register(self):
        """SraSpotDescriptor.register should record primers and barcodes."""
        self.assertEqual(self.spot_descriptor.adapter, 'TACG')
        self.assertEqual(self.spot_descriptor.linker, None)
        self.assertEqual(self.spot_descriptor.barcodes, set(['TTAACCGG']))

    def test_to_xml(self):
        """SraSpotDescriptor.to_xml should return valid XML fragment."""
        expected = '''
        <SPOT_DESCRIPTOR>
          <SPOT_DECODE_SPEC>
            <READ_SPEC>
              <READ_INDEX>0</READ_INDEX>
              <READ_CLASS>Technical Read</READ_CLASS>
              <READ_TYPE>Adapter</READ_TYPE>
              <EXPECTED_BASECALL>TACG</EXPECTED_BASECALL>
            </READ_SPEC>
            <READ_SPEC>
              <READ_INDEX>1</READ_INDEX>
              <READ_LABEL>barcode</READ_LABEL>
              <READ_CLASS>Technical Read</READ_CLASS>
              <READ_TYPE>BarCode</READ_TYPE>
              <EXPECTED_BASECALL_TABLE>
                <BASECALL match_edge="full" max_mismatch="0" min_match="8" read_group_tag="Sample1">TTAACCGG</BASECALL>
              </EXPECTED_BASECALL_TABLE>
            </READ_SPEC>
            <READ_SPEC>
              <READ_INDEX>2</READ_INDEX>
              <READ_LABEL>rRNA_primer</READ_LABEL>
              <READ_CLASS>Technical Read</READ_CLASS>
              <READ_TYPE>Primer</READ_TYPE>
              <EXPECTED_BASECALL_TABLE>
                <BASECALL match_edge="full" max_mismatch="0" min_match="15" read_group_tag="V1-V2">CTGCTGCCTYCCGTA</BASECALL>
              </EXPECTED_BASECALL_TABLE>
            </READ_SPEC>
            <READ_SPEC>
              <READ_INDEX>3</READ_INDEX>
              <READ_CLASS>Application Read</READ_CLASS>
              <READ_TYPE>Forward</READ_TYPE>
              <RELATIVE_ORDER follows_read_index="2" />
            </READ_SPEC>
          </SPOT_DECODE_SPEC>
        </SPOT_DESCRIPTOR>'''
        self.assertEqual(pretty_xml(self.spot_descriptor.to_xml(), 4), expected)


class SraSampleDescriptorTests(TestCase):
    def setUp(self):
        self.attrs1 = {
            'DEFAULT_SAMPLE_CENTER': 'NCBI',
            'DEFAULT_SAMPLE_NAME': 'protected_unassigned',
            'DEFAULT_SAMPLE_ACCESSION': '12345',
            'SAMPLE_ALIAS': 'mysample',
            'SAMPLE_CENTER': 'UPENNBL',
            'POOL_MEMBER_NAME': 'pm',
            'POOL_PROPORTION': '1.0',
            'BARCODE': 'AAAATTTT',
            'BARCODE_READ_GROUP_TAG': 'run1_AAAATTTT',
            }
    
    def test_register(self):
        """SraSampleDescriptor.register should raise ValueError for inconsistent entries."""
        a = SraSampleDescriptor.from_entry(self.attrs1)
        a.register(self.attrs1)
        self.assertEquals(len(a.pool_members), 1)

        incompatible_attrs = {
            'DEFAULT_SAMPLE_CENTER': 'notNCBI',
            'DEFAULT_SAMPLE_NAME': 'protected_unassigned',
            'DEFAULT_SAMPLE_ACCESSION': '12345',
            'SAMPLE_ALIAS': 'mysample_2',
            'SAMPLE_CENTER': 'UPENNBL',
            'POOL_MEMBER_NAME': 'pm2',
            'POOL_PROPORTION': '1.0',
            'BARCODE': 'AAAATTTC',
            'BARCODE_READ_GROUP_TAG': 'run1_AAAATTTC',
            }
        self.assertRaises(ValueError, a.register, incompatible_attrs)

    def test_to_xml(self):
        """SraSampleDescriptor.to_xml should return valid XML fragment."""
        a = SraSampleDescriptor.from_entry(self.attrs1)
        a.register(self.attrs1)
        expected = '''
        <SAMPLE_DESCRIPTOR accession="12345" refcenter="NCBI" refname="protected_unassigned">
          <POOL>
            <MEMBER member_name="pm" proportion="1.0" refcenter="UPENNBL" refname="mysample">
              <READ_LABEL read_group_tag="run1_AAAATTTT">barcode</READ_LABEL>
            </MEMBER>
          </POOL>
        </SAMPLE_DESCRIPTOR>'''
        self.assertEquals(pretty_xml(a.to_xml(), 4), expected)

class SraPoolMemberTests(TestCase):
    def test_from_entry(self):
        """SraPoolMember.from_entry should raise ValueError if neither barcode_tag or primer_tag are provided."""
        vals = {
            'SAMPLE_ALIAS': 'x',
            'SAMPLE_CENTER': 'NCBI',
            'POOL_MEMBER_NAME': 'y',
            'POOL_PROPORTION': '0.1',
            'POOL_MEMBER_ACCESSION': '12345',
            }
        self.assertRaises(ValueError, SraPoolMember.from_entry, vals)
    
    def test_to_xml(self):
        attrs = {
            'refname': 'x',
            'refcenter': 'NCBI',
            'member_name': 'y',
            'proportion': '0.1',
            'barcode_tag': 'bc1',
            }
        a = SraPoolMember(attrs)
        expected = '''
        <MEMBER member_name="y" proportion="0.1" refcenter="NCBI" refname="x">
          <READ_LABEL read_group_tag="bc1">barcode</READ_LABEL>
        </MEMBER>'''
        self.assertEqual(pretty_xml(a.to_xml(), 4), expected)

        attrs['primer_tag'] = 'V1-V3'
        b = SraPoolMember(attrs)
        expected_with_primer = '''
        <MEMBER member_name="y" proportion="0.1" refcenter="NCBI" refname="x">
          <READ_LABEL read_group_tag="bc1">barcode</READ_LABEL>
          <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
        </MEMBER>'''
        self.assertEqual(
            pretty_xml(b.to_xml(), 4), expected_with_primer)

class SraRunTests(TestCase):
    def test_to_xml(self):
        """SraRun.to_xml should produce valid XML fragment."""
        class MockSraDataBlock:
            def to_xml(self):
                return ET.fromstring(
                    '<DATA_BLOCK serial="1" name="F6AVWTA02" region="0" '
                    'member_name=""><FILES>'
                    '<FILE filename="bodysites_study_default_F6AVWTA02.sff" '
                    'filetype="sff" checksum_method="MD5" '
                    'checksum="d41d8cd98f00b204e9800998ecf8427e"/>'
                    '</FILES></DATA_BLOCK>'
                    )
        run = SraRun.from_entry({
            'RUN_ALIAS': 'GHNNEQQ',
            'RUN_CENTER': 'UPENNBL',
            'EXPERIMENT_ALIAS': 'my_expt_alias',
            'STUDY_CENTER': 'NCBI',
            })
        run.data_blocks.append(MockSraDataBlock())
        observed = run.to_xml()
        expected = '''
        <RUN alias="GHNNEQQ" center_name="UPENNBL" run_center="UPENNBL">
          <EXPERIMENT_REF refcenter="NCBI" refname="my_expt_alias" />
          <DATA_BLOCK member_name="" name="F6AVWTA02" region="0" serial="1">
            <FILES>
              <FILE checksum="d41d8cd98f00b204e9800998ecf8427e" checksum_method="MD5" filename="bodysites_study_default_F6AVWTA02.sff" filetype="sff" />
            </FILES>
          </DATA_BLOCK>
        </RUN>'''
        self.assertEqual(pretty_xml(observed, 4), expected)

class SraDataBlockTests(TestCase):
    def test_from_entry(self):
        """SraDataBlock.from_entry should preserve difference between blank and null member_name."""
        vals = {
            'RUN_PREFIX': 'F6AVWTA02',
            'REGION': '2',
            }
        a = SraDataBlock.from_entry(vals)
        self.assertEqual(a.member_name, None)

        vals['POOL_MEMBER_NAME'] = ''
        b = SraDataBlock.from_entry(vals)
        self.assertEqual(b.member_name, '')
        self.assertEqual(b.serial, '1')
        self.assertEqual(b.name, 'F6AVWTA02')
        self.assertEqual(b.region, '2')
    
    def test_to_xml(self):
        """SraDataBlock.to_xml should return valid XML fragment."""
        class MockSraFile:
            def to_xml(self):
                return ET.fromstring(
                    '<FILE checksum="d41d8cd98f00b204e9800998ecf8427e" '
                    'checksum_method="MD5" '
                    'filename="bodysites_study_default_F6AVWTA02.sff" '
                    'filetype="sff" />'
                    )
        attrs = {
            'member_name': '',
            'name': 'F6AVWTA02',
            'region': '0',
            }
        block = SraDataBlock(attrs)
        block.files.append(MockSraFile())
        observed = block.to_xml()
        expected = '''
        <DATA_BLOCK member_name="" name="F6AVWTA02" region="0" serial="1">
          <FILES>
            <FILE checksum="d41d8cd98f00b204e9800998ecf8427e" checksum_method="MD5" filename="bodysites_study_default_F6AVWTA02.sff" filetype="sff" />
          </FILES>
        </DATA_BLOCK>'''
        self.assertEqual(pretty_xml(observed, 4), expected)   

class SraFileTests(TestCase):
    def test_from_entry(self):
        """SraFile.from_entry should return new instance, using existing CHECKSUM value."""
        test_dir = os.path.dirname(os.path.abspath(__file__))

        f = SraFile.from_entry({
            'CHECKSUM': 'f05c1bb96759a0c3c9bcc0196ceac3bb',
            'POOL_MEMBER_FILENAME': 'B-2004-03-S1.sff',
            })
        self.assertEqual(f.filename, 'B-2004-03-S1.sff')
        self.assertEqual(f.checksum, 'f05c1bb96759a0c3c9bcc0196ceac3bb')

    def test_to_xml(self):
        """SraFile.to_xml should return valid XML fragment."""
        f = SraFile({
            'checksum': '12345',
            'filename': 'abc.sff',
            })
        observed = f.to_xml()
        expected = (
            '<FILE checksum="12345" checksum_method="MD5" filename="abc.sff" '
            'filetype="sff" />')
        self.assertEqual(ET.tostring(observed), expected)


experiment = '''\
#EXPERIMENT_ALIAS	EXPERIMENT_CENTER	EXPERIMENT_TITLE	STUDY_REF	STUDY_CENTER	EXPERIMENT_DESIGN_DESCRIPTION	LIBRARY_CONSTRUCTION_PROTOCOL	SAMPLE_ALIAS	SAMPLE_CENTER	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	POOL_PROPORTION	BARCODE_READ_GROUP_TAG	BARCODE	LINKER	PRIMER_READ_GROUP_TAG	KEY_SEQ	PRIMER	RUN_PREFIX	REGION	PLATFORM	RUN_ALIAS	RUN_CENTER	RUN_DATE	INSTRUMENT_NAME
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015438	NCBI	F6AVWTA01_2878_700015438_V1-V3	B-2004-03-S1.sff	0.014492754	F6AVWTA01_ATGTTCGATG	ATGTTCGATG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2878_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015438	NCBI	F6AVWTA02_2878_700015438_V1-V3	B-2008-05-S1.sff	0.014492754	F6AVWTA02_ATGTTCTAGT	ATGTTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2878_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015470	NCBI	F6AVWTA01_2866_700015470_V1-V3	B-2004-04-S1.sff	0.014492754	F6AVWTA01_GCTCTACGTC	GCTCTACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2866_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015470	NCBI	F6AVWTA02_2866_700015470_V1-V3	B-2008-08-S1.sff	0.014492754	F6AVWTA02_GCTCTGTACT	GCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2866_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015766	NCBI	F6AVWTA01_2898_700015766_V1-V3	B-2004-08-S1.sff	0.014492754	F6AVWTA01_CATGAGCGTC	CATGAGCGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2898_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015766	NCBI	F6AVWTA02_2898_700015766_V1-V3	B-2009-06-S1.sff	0.014492754	F6AVWTA02_CATGAGCGTG	CATGAGCGTG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2898_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015468	NCBI	F6AVWTA01_2865_700015468_V1-V3	B-2005-06-S1.sff	0.014492754	F6AVWTA01_AGTACGTACT	AGTACGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2865_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700015468	NCBI	F6AVWTA02_2865_700015468_V1-V3	B-2011-01-S1.sff	0.014492754	F6AVWTA02_AGTACACGTC	AGTACACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2865_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700016371	NCBI	F6AVWTA01_2907_700016371_V1-V3	B-2006-03-S1.sff	0.014492754	F6AVWTA01_TCTCTCTAGT	TCTCTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2907_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	700016371	NCBI	F6AVWTA02_2907_700016371_V1-V3	B-2011-02-S1.sff	0.014492754	F6AVWTA02_TCTCTGTACT	TCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2907_F6AVWTA02	JCVI	NULL	NULL
'''

metagenomic_experiment = '''\
#EXPERIMENT_ALIAS	EXPERIMENT_ACCESSION	EXPERIMENT_CENTER	EXPERIMENT_TITLE	STUDY_REF	STUDY_CENTER	EXPERIMENT_DESIGN_DESCRIPTION	LIBRARY_CONSTRUCTION_PROTOCOL	LIBRARY_STRATEGY	LIBRARY_SELECTION	SAMPLE_ALIAS	SAMPLE_CENTER	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	POOL_PROPORTION	BARCODE_READ_GROUP_TAG	BARCODE	LINKER	PRIMER_READ_GROUP_TAG	KEY_SEQ	PRIMER	RUN_PREFIX	REGION	PLATFORM	RUN_ALIAS	RUN_CENTER	RUN_DATE	INSTRUMENT_NAME
bodysites_F6AVWTA01	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015438	NCBI	F6AVWTA01_2878_700015438_V1-V3	B-2004-03-S1.sff	0.014492754	F6AVWTA01_ATGTTCGATG	ATGTTCGATG			TCAG		F6AVWTA01	0	FLX	bodysites_lib2878_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015438	NCBI	F6AVWTA02_2878_700015438_V1-V3	B-2008-05-S1.sff	0.014492754	F6AVWTA02_ATGTTCTAGT	ATGTTCTAGT			TCAG		F6AVWTA02	0	FLX	bodysites_lib2878_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015470	NCBI	F6AVWTA01_2866_700015470_V1-V3	B-2004-04-S1.sff	0.014492754	F6AVWTA01_GCTCTACGTC	GCTCTACGTC			TCAG		F6AVWTA01	0	FLX	bodysites_lib2866_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015470	NCBI	F6AVWTA02_2866_700015470_V1-V3	B-2008-08-S1.sff	0.014492754	F6AVWTA02_GCTCTGTACT	GCTCTGTACT			TCAG		F6AVWTA02	0	FLX	bodysites_lib2866_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015766	NCBI	F6AVWTA01_2898_700015766_V1-V3	B-2004-08-S1.sff	0.014492754	F6AVWTA01_CATGAGCGTC	CATGAGCGTC			TCAG		F6AVWTA01	0	FLX	bodysites_lib2898_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015766	NCBI	F6AVWTA02_2898_700015766_V1-V3	B-2009-06-S1.sff	0.014492754	F6AVWTA02_CATGAGCGTG	CATGAGCGTG			TCAG		F6AVWTA02	0	FLX	bodysites_lib2898_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015468	NCBI	F6AVWTA01_2865_700015468_V1-V3	B-2005-06-S1.sff	0.014492754	F6AVWTA01_AGTACGTACT	AGTACGTACT			TCAG		F6AVWTA01	0	FLX	bodysites_lib2865_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700015468	NCBI	F6AVWTA02_2865_700015468_V1-V3	B-2011-01-S1.sff	0.014492754	F6AVWTA02_AGTACACGTC	AGTACACGTC			TCAG		F6AVWTA02	0	FLX	bodysites_lib2865_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700016371	NCBI	F6AVWTA01_2907_700016371_V1-V3	B-2006-03-S1.sff	0.014492754	F6AVWTA01_TCTCTCTAGT	TCTCTCTAGT			TCAG		F6AVWTA01	0	FLX	bodysites_lib2907_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	SRX01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	WGS	RANDOM	700016371	NCBI	F6AVWTA02_2907_700016371_V1-V3	B-2011-02-S1.sff	0.014492754	F6AVWTA02_TCTCTGTACT	TCTCTGTACT			TCAG		F6AVWTA02	0	FLX	bodysites_lib2907_F6AVWTA02	JCVI	NULL	NULL
'''

experiment_with_accessions = '''\
#EXPERIMENT_ALIAS	EXPERIMENT_CENTER	EXPERIMENT_TITLE	STUDY_REF	STUDY_CENTER	STUDY_ACCESSION	EXPERIMENT_DESIGN_DESCRIPTION	LIBRARY_CONSTRUCTION_PROTOCOL	SAMPLE_ALIAS	DEFAULT_SAMPLE_ACCESSION	SAMPLE_CENTER	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	POOL_MEMBER_ACCESSION	POOL_PROPORTION	BARCODE_READ_GROUP_TAG	BARCODE	LINKER	PRIMER_READ_GROUP_TAG	KEY_SEQ	PRIMER	RUN_PREFIX	REGION	PLATFORM	RUN_ALIAS	RUN_CENTER	RUN_DATE	INSTRUMENT_NAME
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015438	SRS077	NCBI	F6AVWTA01_2878_700015438_V1-V3	B-2004-03-S1.sff	SRS001	0.014492754	F6AVWTA01_ATGTTCGATG	ATGTTCGATG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2878_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015438	SRS077	NCBI	F6AVWTA02_2878_700015438_V1-V3	B-2008-05-S1.sff	SRS002	0.014492754	F6AVWTA02_ATGTTCTAGT	ATGTTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2878_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015470	SRS077	NCBI	F6AVWTA01_2866_700015470_V1-V3	B-2004-04-S1.sff	SRS003	0.014492754	F6AVWTA01_GCTCTACGTC	GCTCTACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2866_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015470	SRS077	NCBI	F6AVWTA02_2866_700015470_V1-V3	B-2008-08-S1.sff	SRS004	0.014492754	F6AVWTA02_GCTCTGTACT	GCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2866_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015766	SRS077	NCBI	F6AVWTA01_2898_700015766_V1-V3	B-2004-08-S1.sff	SRS005	0.014492754	F6AVWTA01_CATGAGCGTC	CATGAGCGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2898_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015766	SRS077	NCBI	F6AVWTA02_2898_700015766_V1-V3	B-2009-06-S1.sff	SRS006	0.014492754	F6AVWTA02_CATGAGCGTG	CATGAGCGTG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2898_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015468	SRS077	NCBI	F6AVWTA01_2865_700015468_V1-V3	B-2005-06-S1.sff	SRS007	0.014492754	F6AVWTA01_AGTACGTACT	AGTACGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2865_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700015468	SRS077	NCBI	F6AVWTA02_2865_700015468_V1-V3	B-2011-01-S1.sff	SRS008	0.014492754	F6AVWTA02_AGTACACGTC	AGTACACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2865_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700016371	SRS077	NCBI	F6AVWTA01_2907_700016371_V1-V3	B-2006-03-S1.sff	SRS009	0.014492754	F6AVWTA01_TCTCTCTAGT	TCTCTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2907_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	SRP001	Pool of samples from different individual subjects	Dummy Protocol	700016371	SRS077	NCBI	F6AVWTA02_2907_700016371_V1-V3	B-2011-02-S1.sff	SRS010	0.014492754	F6AVWTA02_TCTCTGTACT	TCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2907_F6AVWTA02	JCVI	NULL	NULL
'''

experiment_with_default_sample = '''\
#EXPERIMENT_ALIAS	EXPERIMENT_CENTER	EXPERIMENT_TITLE	STUDY_REF	STUDY_CENTER	EXPERIMENT_DESIGN_DESCRIPTION	LIBRARY_CONSTRUCTION_PROTOCOL	DEFAULT_SAMPLE_NAME	DEFAULT_SAMPLE_ACCESSION	DEFAULT_SAMPLE_CENTER	SAMPLE_ALIAS	SAMPLE_CENTER	POOL_MEMBER_NAME	POOL_MEMBER_FILENAME	POOL_PROPORTION	BARCODE_READ_GROUP_TAG	BARCODE	LINKER	PRIMER_READ_GROUP_TAG	KEY_SEQ	PRIMER	RUN_PREFIX	REGION	PLATFORM	RUN_ALIAS	RUN_CENTER	RUN_DATE	INSTRUMENT_NAME
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015438	NCBI	F6AVWTA01_2878_700015438_V1-V3	B-2004-03-S1.sff	0.014492754	F6AVWTA01_ATGTTCGATG	ATGTTCGATG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2878_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015438	NCBI	F6AVWTA02_2878_700015438_V1-V3	B-2008-05-S1.sff	0.014492754	F6AVWTA02_ATGTTCTAGT	ATGTTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2878_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015470	NCBI	F6AVWTA01_2866_700015470_V1-V3	B-2004-04-S1.sff	0.014492754	F6AVWTA01_GCTCTACGTC	GCTCTACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2866_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015470	NCBI	F6AVWTA02_2866_700015470_V1-V3	B-2008-08-S1.sff	0.014492754	F6AVWTA02_GCTCTGTACT	GCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2866_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015766	NCBI	F6AVWTA01_2898_700015766_V1-V3	B-2004-08-S1.sff	0.014492754	F6AVWTA01_CATGAGCGTC	CATGAGCGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2898_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015766	NCBI	F6AVWTA02_2898_700015766_V1-V3	B-2009-06-S1.sff	0.014492754	F6AVWTA02_CATGAGCGTG	CATGAGCGTG		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2898_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015468	NCBI	F6AVWTA01_2865_700015468_V1-V3	B-2005-06-S1.sff	0.014492754	F6AVWTA01_AGTACGTACT	AGTACGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2865_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700015468	NCBI	F6AVWTA02_2865_700015468_V1-V3	B-2011-01-S1.sff	0.014492754	F6AVWTA02_AGTACACGTC	AGTACACGTC		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2865_F6AVWTA02	JCVI	NULL	NULL
bodysites_F6AVWTA01	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700016371	NCBI	F6AVWTA01_2907_700016371_V1-V3	B-2006-03-S1.sff	0.014492754	F6AVWTA01_TCTCTCTAGT	TCTCTCTAGT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA01	0	FLX	bodysites_lib2907_F6AVWTA01	JCVI	NULL	NULL
bodysites_F6AVWTA02	JCVI	Survey of multiple body sites	bodysites_study	bodysites	Pool of samples from different individual subjects	Dummy Protocol	unidentified-protected	SRS026543	NCBI	700016371	NCBI	F6AVWTA02_2907_700016371_V1-V3	B-2011-02-S1.sff	0.014492754	F6AVWTA02_TCTCTGTACT	TCTCTGTACT		V1-V3	TCAG	TAATCCGCGGCTGCTGG	F6AVWTA02	0	FLX	bodysites_lib2907_F6AVWTA02	JCVI	NULL	NULL
'''

attrs = '''\
#Experiment	Attribute	Value
bodysites_F6AVWTA01	library_strategy	targeted-locus
bodysites_F6AVWTA01	gene	16S rRNA V1-V3 region
bodysites_F6AVWTA02	library_strategy	targeted-locus
bodysites_F6AVWTA02	gene	16S rRNA V1-V3 region'''

links = '''\
#Experiment	Link Name	Link URL
bodysites_F6AVWTA01	bodysites Library Construction Protocol	http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf
bodysites_F6AVWTA02	bodysites Library Construction Protocol	http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf
'''

experiment_with_default_sample_xml_str = '''<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EXPERIMENT alias="bodysites_F6AVWTA01" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="unidentified-protected" refcenter="NCBI" accession="SRS026543">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA01_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_ATGTTCGATG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA01_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_GCTCTACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA01_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_CATGAGCGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA01_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_AGTACGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA01_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_TCTCTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA01</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL match_edge="full" max_mismatch="0" read_group_tag="F6AVWTA01_ATGTTCGATG" min_match="10">ATGTTCGATG</BASECALL>
              <BASECALL match_edge="full" max_mismatch="0" read_group_tag="F6AVWTA01_GCTCTACGTC" min_match="10">GCTCTACGTC</BASECALL>
              <BASECALL match_edge="full" max_mismatch="0" read_group_tag="F6AVWTA01_CATGAGCGTC" min_match="10">CATGAGCGTC</BASECALL>
              <BASECALL match_edge="full" max_mismatch="0" read_group_tag="F6AVWTA01_AGTACGTACT" min_match="10">AGTACGTACT</BASECALL>
              <BASECALL match_edge="full" max_mismatch="0" read_group_tag="F6AVWTA01_TCTCTCTAGT" min_match="10">TCTCTCTAGT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL match_edge="full" max_mismatch="0" min_match="17" read_group_tag="V1-V3">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
  <EXPERIMENT alias="bodysites_F6AVWTA02" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="unidentified-protected" refcenter="NCBI" accession="SRS026543">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA02_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_ATGTTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA02_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_GCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA02_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_CATGAGCGTG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA02_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_AGTACACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA02_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_TCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA02</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA02_ATGTTCTAGT" min_match="10" max_mismatch="0" match_edge="full">ATGTTCTAGT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_GCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">GCTCTGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_CATGAGCGTG" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_AGTACACGTC" min_match="10" max_mismatch="0" match_edge="full">AGTACACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_TCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">TCTCTGTACT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="V1-V3" min_match="17" max_mismatch="0" match_edge="full">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
</EXPERIMENT_SET>'''

experiment_with_accessions_xml_str = '''<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EXPERIMENT alias="bodysites_F6AVWTA01" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI" accession="SRP001"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refcenter="NCBI" accession="SRS077">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA01_2878_700015438_V1-V3" proportion="0.014492754" accession="SRS001">
            <READ_LABEL read_group_tag="F6AVWTA01_ATGTTCGATG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA01_2866_700015470_V1-V3" proportion="0.014492754" accession="SRS003">
            <READ_LABEL read_group_tag="F6AVWTA01_GCTCTACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA01_2898_700015766_V1-V3" proportion="0.014492754" accession="SRS005">
            <READ_LABEL read_group_tag="F6AVWTA01_CATGAGCGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA01_2865_700015468_V1-V3" proportion="0.014492754" accession="SRS007">
            <READ_LABEL read_group_tag="F6AVWTA01_AGTACGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA01_2907_700016371_V1-V3" proportion="0.014492754" accession="SRS009">
            <READ_LABEL read_group_tag="F6AVWTA01_TCTCTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA01</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA01_ATGTTCGATG" min_match="10" max_mismatch="0" match_edge="full">ATGTTCGATG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_GCTCTACGTC" min_match="10" max_mismatch="0" match_edge="full">GCTCTACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_CATGAGCGTC" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_AGTACGTACT" min_match="10" max_mismatch="0" match_edge="full">AGTACGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_TCTCTCTAGT" min_match="10" max_mismatch="0" match_edge="full">TCTCTCTAGT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="V1-V3" min_match="17" max_mismatch="0" match_edge="full">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
  <EXPERIMENT alias="bodysites_F6AVWTA02" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI" accession="SRP001"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refcenter="NCBI" accession="SRS077">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA02_2878_700015438_V1-V3" proportion="0.014492754" accession="SRS002">
            <READ_LABEL read_group_tag="F6AVWTA02_ATGTTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA02_2866_700015470_V1-V3" proportion="0.014492754" accession="SRS004">
            <READ_LABEL read_group_tag="F6AVWTA02_GCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA02_2898_700015766_V1-V3" proportion="0.014492754" accession="SRS006">
            <READ_LABEL read_group_tag="F6AVWTA02_CATGAGCGTG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA02_2865_700015468_V1-V3" proportion="0.014492754" accession="SRS008">
            <READ_LABEL read_group_tag="F6AVWTA02_AGTACACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA02_2907_700016371_V1-V3" proportion="0.014492754" accession="SRS010">
            <READ_LABEL read_group_tag="F6AVWTA02_TCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA02</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA02_ATGTTCTAGT" min_match="10" max_mismatch="0" match_edge="full">ATGTTCTAGT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_GCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">GCTCTGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_CATGAGCGTG" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_AGTACACGTC" min_match="10" max_mismatch="0" match_edge="full">AGTACACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_TCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">TCTCTGTACT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="V1-V3" min_match="17" max_mismatch="0" match_edge="full">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
</EXPERIMENT_SET>'''

metagenomic_experiment_xml_str = '''<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EXPERIMENT alias="bodysites_F6AVWTA01" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="bodysites_study_default" refcenter="NCBI">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA01_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_ATGTTCGATG">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA01_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_GCTCTACGTC">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA01_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_CATGAGCGTC">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA01_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_AGTACGTACT">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA01_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_TCTCTCTAGT">barcode</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA01</LIBRARY_NAME>
        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA01_ATGTTCGATG" min_match="10" max_mismatch="0" match_edge="full">ATGTTCGATG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_GCTCTACGTC" min_match="10" max_mismatch="0" match_edge="full">GCTCTACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_CATGAGCGTC" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_AGTACGTACT" min_match="10" max_mismatch="0" match_edge="full">AGTACGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_TCTCTCTAGT" min_match="10" max_mismatch="0" match_edge="full">TCTCTCTAGT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="1"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
  </EXPERIMENT>
  <EXPERIMENT alias="bodysites_F6AVWTA02" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="bodysites_study_default" refcenter="NCBI">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA02_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_ATGTTCTAGT">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA02_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_GCTCTGTACT">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA02_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_CATGAGCGTG">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA02_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_AGTACACGTC">barcode</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA02_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_TCTCTGTACT">barcode</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA02</LIBRARY_NAME>
        <LIBRARY_STRATEGY>WGS</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>RANDOM</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA02_ATGTTCTAGT" min_match="10" max_mismatch="0" match_edge="full">ATGTTCTAGT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_GCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">GCTCTGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_CATGAGCGTG" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_AGTACACGTC" min_match="10" max_mismatch="0" match_edge="full">AGTACACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_TCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">TCTCTGTACT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="1"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
  </EXPERIMENT>
</EXPERIMENT_SET>'''

experiment_xml_str = '''<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <EXPERIMENT alias="bodysites_F6AVWTA01" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="bodysites_study_default" refcenter="NCBI">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA01_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_ATGTTCGATG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA01_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_GCTCTACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA01_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_CATGAGCGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA01_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_AGTACGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA01_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA01_TCTCTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA01</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA01_ATGTTCGATG" min_match="10" max_mismatch="0" match_edge="full">ATGTTCGATG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_GCTCTACGTC" min_match="10" max_mismatch="0" match_edge="full">GCTCTACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_CATGAGCGTC" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_AGTACGTACT" min_match="10" max_mismatch="0" match_edge="full">AGTACGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA01_TCTCTCTAGT" min_match="10" max_mismatch="0" match_edge="full">TCTCTCTAGT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="V1-V3" min_match="17" max_mismatch="0" match_edge="full">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
  <EXPERIMENT alias="bodysites_F6AVWTA02" center_name="JCVI">
    <TITLE>Survey of multiple body sites</TITLE>
    <STUDY_REF refname="bodysites_study" refcenter="NCBI"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>Pool of samples from different individual subjects</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="bodysites_study_default" refcenter="NCBI">
        <POOL>
          <MEMBER refname="700015438" refcenter="NCBI" member_name="F6AVWTA02_2878_700015438_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_ATGTTCTAGT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015470" refcenter="NCBI" member_name="F6AVWTA02_2866_700015470_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_GCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015766" refcenter="NCBI" member_name="F6AVWTA02_2898_700015766_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_CATGAGCGTG">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700015468" refcenter="NCBI" member_name="F6AVWTA02_2865_700015468_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_AGTACACGTC">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
          <MEMBER refname="700016371" refcenter="NCBI" member_name="F6AVWTA02_2907_700016371_V1-V3" proportion="0.014492754">
            <READ_LABEL read_group_tag="F6AVWTA02_TCTCTGTACT">barcode</READ_LABEL>
            <READ_LABEL read_group_tag="V1-V3">rRNA_primer</READ_LABEL>
          </MEMBER>
        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>bodysites_F6AVWTA02</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>Dummy Protocol</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>
      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
            <EXPECTED_BASECALL>TCAG</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="F6AVWTA02_ATGTTCTAGT" min_match="10" max_mismatch="0" match_edge="full">ATGTTCTAGT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_GCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">GCTCTGTACT</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_CATGAGCGTG" min_match="10" max_mismatch="0" match_edge="full">CATGAGCGTG</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_AGTACACGTC" min_match="10" max_mismatch="0" match_edge="full">AGTACACGTC</BASECALL>
              <BASECALL read_group_tag="F6AVWTA02_TCTCTGTACT" min_match="10" max_mismatch="0" match_edge="full">TCTCTGTACT</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>
              <BASECALL read_group_tag="V1-V3" min_match="17" max_mismatch="0" match_edge="full">TAATCCGCGGCTGCTGG</BASECALL>
            </EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
            <RELATIVE_ORDER follows_read_index="2"/>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <LS454>
        <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
        <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
        <FLOW_COUNT>400</FLOW_COUNT>
      </LS454>
    </PLATFORM>
    <PROCESSING>
      <BASE_CALLS>
        <SEQUENCE_SPACE>Base Space</SEQUENCE_SPACE>
        <BASE_CALLER>454 BaseCaller</BASE_CALLER>
      </BASE_CALLS>
      <QUALITY_SCORES qtype="phred">
        <QUALITY_SCORER>454 BaseCaller</QUALITY_SCORER>
        <NUMBER_OF_LEVELS>40</NUMBER_OF_LEVELS>
        <MULTIPLIER>1.0</MULTIPLIER>
      </QUALITY_SCORES>
    </PROCESSING>
    <EXPERIMENT_LINKS>
      <EXPERIMENT_LINK>
        <URL_LINK>
          <LABEL>bodysites Library Construction Protocol</LABEL>
          <URL>http://hmpdacc.org/doc/HMP_MDG_454_16S_Protocol_V4_2_102109.pdf</URL>
        </URL_LINK>
      </EXPERIMENT_LINK>
    </EXPERIMENT_LINKS>
    <EXPERIMENT_ATTRIBUTES>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>library_strategy</TAG>
        <VALUE>targeted-locus</VALUE>
      </EXPERIMENT_ATTRIBUTE>
      <EXPERIMENT_ATTRIBUTE>
        <TAG>gene</TAG>
        <VALUE>16S rRNA V1-V3 region</VALUE>
      </EXPERIMENT_ATTRIBUTE>
    </EXPERIMENT_ATTRIBUTES>
  </EXPERIMENT>
</EXPERIMENT_SET>
'''

run_xml_str = '''<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <RUN alias="bodysites_study_default_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="">
      <FILES>
        <FILE filename="bodysites_study_default_F6AVWTA01.sff" filetype="sff" checksum_method="MD5" checksum="d41d8cd98f00b204e9800998ecf8427e"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2878_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="F6AVWTA01_2878_700015438_V1-V3">
      <FILES>
        <FILE filename="B-2004-03-S1.sff" filetype="sff" checksum_method="MD5" checksum="f05c1bb96759a0c3c9bcc0196ceac3bb"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_study_default_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="">
      <FILES>
        <FILE filename="bodysites_study_default_F6AVWTA02.sff" filetype="sff" checksum_method="MD5" checksum="d41d8cd98f00b204e9800998ecf8427e"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2878_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="F6AVWTA02_2878_700015438_V1-V3">
      <FILES>
        <FILE filename="B-2008-05-S1.sff" filetype="sff" checksum_method="MD5" checksum="9fcf84e2a06a1175124e15064e9b63a1"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2866_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="F6AVWTA01_2866_700015470_V1-V3">
      <FILES>
        <FILE filename="B-2004-04-S1.sff" filetype="sff" checksum_method="MD5" checksum="c693c0f79bfd05d7fb348f0182fbf808"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2866_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="F6AVWTA02_2866_700015470_V1-V3">
      <FILES>
        <FILE filename="B-2008-08-S1.sff" filetype="sff" checksum_method="MD5" checksum="34f4185163c4ce10610a2427ba554ba3"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2898_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="F6AVWTA01_2898_700015766_V1-V3">
      <FILES>
        <FILE filename="B-2004-08-S1.sff" filetype="sff" checksum_method="MD5" checksum="b6f607366b77c2d5b58837870d6b07c1"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2898_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="F6AVWTA02_2898_700015766_V1-V3">
      <FILES>
        <FILE filename="B-2009-06-S1.sff" filetype="sff" checksum_method="MD5" checksum="02d0751ce913ca796a5916803c574636"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2865_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="F6AVWTA01_2865_700015468_V1-V3">
      <FILES>
        <FILE filename="B-2005-06-S1.sff" filetype="sff" checksum_method="MD5" checksum="8d2ef95655a4f0cf4b6e1e2966c6ac30"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2865_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="F6AVWTA02_2865_700015468_V1-V3">
      <FILES>
        <FILE filename="B-2011-01-S1.sff" filetype="sff" checksum_method="MD5" checksum="7852588b980ba08c1ff0d0ca7b686b16"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2907_F6AVWTA01" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA01" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA01" region="0" member_name="F6AVWTA01_2907_700016371_V1-V3">
      <FILES>
        <FILE filename="B-2006-03-S1.sff" filetype="sff" checksum_method="MD5" checksum="a20fa67736ffc9f966827275036954b5"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
  <RUN alias="bodysites_lib2907_F6AVWTA02" center_name="JCVI" run_center="JCVI">
    <EXPERIMENT_REF refname="bodysites_F6AVWTA02" refcenter="bodysites"/>
    <DATA_BLOCK serial="1" name="F6AVWTA02" region="0" member_name="F6AVWTA02_2907_700016371_V1-V3">
      <FILES>
        <FILE filename="B-2011-02-S1.sff" filetype="sff" checksum_method="MD5" checksum="6b2c7045be67a4cf4958d22c5b6ab790"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
</RUN_SET>'''

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
    <ACTION>
      <ADD source="study.xml" schema="study" notes="study metadata"/>
    </ACTION>
    <ACTION>
      <ADD source="sample.xml" schema="sample" notes="sample metadata"/>
    </ACTION>
    <ACTION>
      <RELEASE/>
    </ACTION>
  </ACTIONS>
  <FILES>
    <FILE filename="%s" checksum_method="MD5" checksum="d41d8cd98f00b204e9800998ecf8427e"/>
  </FILES>
</SUBMISSION>
'''

study_twocol_txt = '''#Field	Value	Example	Comments
STUDY_alias	fierer_hand_study	fierer_handstudy	One study per publication: this is used as an id to link files
STUDY_TITLE	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	Expected (or actual) title of the paper
STUDY_TYPE	Metagenomics	Metagenomics	"Should be ""Metagenomics"" for 16S surveys"
STUDY_ABSTRACT	"Short \'abstract\' with special characters <10%."	"Abstract, e.g. of the publication"
STUDY_DESCRIPTION	Targeted Gene Survey from Human Skin	Targeted Gene Survey from Human Skin	"Use ""Targeted Gene Survey"" for 16S or other target gene studies"
CENTER_NAME	CCME	CCME	"NCBI-approved name of sequencing center, e.g. WUGSC"
CENTER_PROJECT_NAME	NULL	NULL	"Name of project as used by the sequencing center, NULL if none."
PMID	19004758	19004758	"PubMed ID of paper describing project, if supplied will write out STUDY_LINK block, can be multiple (comma-delimited)"
'''

study_manycol_txt = '''
#STUDY_alias	STUDY_TITLE	STUDY_TYPE	STUDY_ABSTRACT	STUDY_DESCRIPTION	CENTER_NAME	CENTER_PROJECT_NAME	PMID
fierer_hand_study	"The influence of sex, handedness, and washing on the diversity of hand surface bacteria"	Metagenomics	"Short \'abstract\' with special characters <10%."	Targeted Gene Survey from Human Skin	CCME	NULL	19004758
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

submission_twocol_txt = '''#Field	Value	Example	Comments
accession	SRA003492	SRA003492	"leave blank if not assigned yet, e.g. if new submission"
submission_id	fierer_hand_study	fierer_hand_study	internally unique id for the submission
center_name	CCME	CCME	name of the center preparing the submission
submission_comment	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Free-text comments regarding submission
lab_name	Knight	Knight	"name of lab preparing submission, can differ from center (usually refers to the PI\'s info, not the sequencing center\'s)"
submission_date	2009-10-22T01:23:00-05:00	2009-10-22T01:23:00-05:00	timestamp of submission
CONTACT	Rob Knight;Rob.Knight@Colorado.edu	Rob Knight;Rob.Knight@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
CONTACT	Noah Fierer;Noah.Fierer@Colorado.edu	Noah Fierer;Noah.Fierer@Colorado.edu	"Use semicolon to separate email address from name, can be multiple contacts."
'''

submission_manycol_txt = '''\
#accession	submission_id	center_name	submission_comment	lab_name	submission_date	CONTACT
SRA003492	fierer_hand_study	CCME	"Barcode submission prepared by osulliva@ncbi.nlm.nih.gov, shumwaym@ncbi.nlm.nih.gov"	Knight	2009-10-22T01:23:00-05:00	Rob Knight;Rob.Knight@Colorado.edu,Noah Fierer;Noah.Fierer@Colorado.edu
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
    <ACTION>
      <ADD source="study.xml" schema="study" notes="study metadata" />
    </ACTION>
    <ACTION>
      <ADD source="sample.xml" schema="sample" notes="sample metadata" />
    </ACTION>
    <ACTION>
      <RELEASE />
    </ACTION>
  </ACTIONS>
</SUBMISSION>
'''

minimal_submission_txt = '''\
#SUBMISSION_ID	CENTER_NAME	LAB_NAME	SUBMISSION_DATE	CONTACT
fierer_hand_study	CCME	Knight	2009-10-22T01:23:00-05:00	Rob Knight;Rob.Knight@Colorado.edu,Noah Fierer;Noah.Fierer@Colorado.edu
'''

sample_txt = '''#SAMPLE_ALIAS	TITLE	TAXON_ID	COMMON_NAME	ANONYMIZED_NAME	DESCRIPTION	host_taxon_id	subject	sex	hand	age	palm size	dominant hand	hours since wash
fierer_hand_study_default	human hand microbiome	539655	human skin metagenome		"Human palm microbiome, default sample for unclassified reads"								
S1	human hand microbiome	539655	human skin metagenome	subject 1	female right palm	9606	1	female	right	18	9.5	right	less than 2
S2	human hand microbiome	539655	human skin metagenome	subject 1	female left palm	9606	1	female	left	18	9.5	right	less than 2
'''

sample_with_lowercase_txt = '''#sample_alias	title	taxon_id	COMMON_NAME	anonymized_name	description	host_taxon_id	subject	sex	hand	age	palm size	dominant hand	hours since wash
fierer_hand_study_default	human hand microbiome	539655	human skin metagenome		"Human palm microbiome, default sample for unclassified reads"								
S1	human hand microbiome	539655	human skin metagenome	subject 1	female right palm	9606	1	female	right	18	9.5	right	less than 2
S2	human hand microbiome	539655	human skin metagenome	subject 1	female left palm	9606	1	female	left	18	9.5	right	less than 2
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
      <SAMPLE_ATTRIBUTE>
        <TAG>age</TAG>
        <VALUE>18</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>dominant hand</TAG>
        <VALUE>right</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>hand</TAG>
        <VALUE>right</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>host_taxon_id</TAG>
        <VALUE>9606</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>hours since wash</TAG>
        <VALUE>less than 2</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>palm size</TAG>
        <VALUE>9.5</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>sex</TAG>
        <VALUE>female</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>subject</TAG>
        <VALUE>1</VALUE>
      </SAMPLE_ATTRIBUTE>
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
      <SAMPLE_ATTRIBUTE>
        <TAG>age</TAG>
        <VALUE>18</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>dominant hand</TAG>
        <VALUE>right</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>hand</TAG>
        <VALUE>left</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>host_taxon_id</TAG>
        <VALUE>9606</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>hours since wash</TAG>
        <VALUE>less than 2</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>palm size</TAG>
        <VALUE>9.5</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>sex</TAG>
        <VALUE>female</VALUE>
      </SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE>
        <TAG>subject</TAG>
        <VALUE>1</VALUE>
      </SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
</SAMPLE_SET>
'''

# SRA Production release 1.1
experiment_xsd = '''\
<?xml version="1.0" encoding="UTF-8"?>
<!-- INSDC Short Read Archive resource Experiment (SRX) object XML specification -->
<!-- $Id$ -->
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <!-- BEGIN COMMON BLOCK            TO BE NORMALIZED IN THE FUTURE -->
    <xs:attributeGroup name="NameGroup">
        <xs:attribute name="alias" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Submitter designated name of the SRA document of this type.  At minimum alias should
                    be unique throughout the submission of this document type.  If center_name is specified, the name should
                    be unique in all submissions from that center of this document type.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="center_name" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Owner authority of this document and namespace for submitter\'s name of this document. 
                    If not provided, then the submitter is regarded as "Individual" and document resolution
                    can only happen within the submission.
                </xs:documentation>
            </xs:annotation>  
        </xs:attribute>
        <xs:attribute name="broker_name" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Broker authority of this document.  If not provided, then the broker is considered "direct".
                </xs:documentation>
            </xs:annotation>  
        </xs:attribute>
        <xs:attribute name="accession" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    The document\'s accession as assigned by the Home Archive.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
    </xs:attributeGroup>
    <xs:attributeGroup name="RefNameGroup">
        <xs:attribute name="refname" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Identifies a record by name that is known within the namespace defined by attribute "refcenter"
                    Use this field when referencing an object for which an accession has not yet been issued.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="refcenter" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    The namespace of the attribute "refname". When absent, the namespace is assumed to be the current submission.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="accession" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Identifies a record by its accession.  The scope of resolution is the entire Archive.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
    </xs:attributeGroup>
    <xs:complexType name="AttributeType">
        <xs:annotation>
            <xs:documentation>
                Reusable attributes to encode tag-value pairs with optional units.
            </xs:documentation>
        </xs:annotation>
        <xs:all>
            <xs:element name="TAG" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                <xs:annotation>
                    <xs:documentation>
                        Name of the attribute.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
            <xs:element name="VALUE" type="xs:string" minOccurs="1" maxOccurs="1">
                <xs:annotation>
                    <xs:documentation>
                        Value of the attribute.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
            <xs:element name="UNITS" type="xs:string" minOccurs="0" maxOccurs="1">          
                <xs:annotation>
                    <xs:documentation>
                        Optional scientific units.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
        </xs:all>
    </xs:complexType>
    <xs:complexType name="SraLinkType">
        <xs:annotation>
            <xs:documentation>
                The SraLinkType mechanism encodes local references between SRA objects.  These 
                references are local to the Home Archive and during archival are resolved to
                exportable references suitable for mirroring between Archives.  SraLinks can be used
                as temporary place holders for pre-published or post-suppressed relationships that
                will not be mirrored between Archives.
            </xs:documentation>
        </xs:annotation>
        <xs:attribute name="sra_object_type" use="optional">
            <xs:annotation>
                <xs:documentation>
                    SRA link type.
                </xs:documentation>
            </xs:annotation>
            <xs:simpleType>
                <xs:restriction base="xs:string">
                    <xs:enumeration value="STUDY"/>
                    <xs:enumeration value="SAMPLE"/>
                    <xs:enumeration value="ANALYSIS"/>
                    <xs:enumeration value="EXPERIMENT"/>
                    <xs:enumeration value="RUN"/>            
                </xs:restriction>
            </xs:simpleType>
        </xs:attribute>
        <xs:attributeGroup ref="RefNameGroup"/>
        
    </xs:complexType>
    
    
    
    
    <xs:complexType name="LinkType">
        <xs:annotation>
            <xs:documentation>
                Reusable external links type to encode URL links, Entrez links, and db_xref links.
            </xs:documentation>
        </xs:annotation>
        <xs:choice>                   
            <xs:element name="SRA_LINK" type="SraLinkType"/>
            <xs:element name="URL_LINK">
                <xs:complexType>
                    <xs:all>
                        <xs:element name="LABEL" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    Text label to display for the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="URL" minOccurs="1" maxOccurs="1" type="xs:anyURI">
                            <xs:annotation>
                                <xs:documentation>
                                    The internet service link (file:, http:, ftp:, etc).
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:all>
                </xs:complexType>
            </xs:element>
            <xs:element name="XREF_LINK">
                <xs:complexType>
                    <xs:all>
                        <xs:element name="DB" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    INSDC controlled vocabulary of permitted cross references.  Please see http://www.insdc.org/page.php?page=db_xref .
                                    For example, FLYBASE.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="ID" minOccurs="1" maxOccurs="1" type="xs:string">
                            <xs:annotation>
                                <xs:documentation>
                                    Accession in the referenced database.    For example,  FBtr0080008 (in FLYBASE).
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="LABEL" type="xs:string" minOccurs="0" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    Text label to display for the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:all>
                </xs:complexType>
            </xs:element>
            <xs:element name="ENTREZ_LINK">
                <xs:complexType>
                    <xs:sequence>
                        <xs:element name="DB" type="xs:string" minOccurs="1" maxOccurs="1">
                            <xs:annotation>
                                <xs:documentation>
                                    NCBI controlled vocabulary of permitted cross references.  Please see http://www.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi? .
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:choice>
                            <xs:element name="ID" type="xs:nonNegativeInteger" minOccurs="1" maxOccurs="1">
                                <xs:annotation>
                                    <xs:documentation>
                                        Numeric record id meaningful to the NCBI Entrez system.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:element>
                            <xs:element name="QUERY" type="xs:string" minOccurs="1" maxOccurs="1">
                                <xs:annotation>
                                    <xs:documentation>
                                        Accession string meaningful to the NCBI Entrez system.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:element>
                        </xs:choice>
                        <xs:element name="LABEL" type="xs:string" minOccurs="0" maxOccurs="1">
                            <xs:annotation>
                                <xs:documentation>
                                    How to label the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:sequence>
                </xs:complexType>
            </xs:element>
        </xs:choice>
    </xs:complexType>
    
    <!--  END COMMON BLOCK -->
 
   
    
    <xs:complexType name="ExperimentType">

              <xs:annotation>
                <xs:documentation>
                  An Experiment specifies of what will be sequenced and how the sequencing will be performed.  
                  It does not contain results.  
                  An Experiment is composed of a design, a platform selection, and processing parameters.
                </xs:documentation>
              </xs:annotation>  
        
                <xs:sequence >
                  <xs:element name="TITLE" type="xs:string" minOccurs="0" maxOccurs="1">
                    <xs:annotation>
                      <xs:documentation>
                        Short text that can be used to call out experiment records in searches or in displays.
                        This element is technically optional but should be used for all new records.
                      </xs:documentation>
                    </xs:annotation>
                  </xs:element>
                  <xs:element name="STUDY_REF" minOccurs="1" maxOccurs="1">
                    <xs:annotation>
                      <xs:documentation>
                        The STUDY_REF descriptor establishes the relationship of the experiment to the parent
                        study.  This can either be the accession of an existing archived study record, or
                        a reference to a new study record in the same submission (which does not yet have an
                        accession).
                      </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                       <xs:attributeGroup ref="RefNameGroup"/>
                    
                    </xs:complexType>
                  </xs:element>
                  <xs:element name="DESIGN" maxOccurs="1" minOccurs="1">
                    <xs:complexType>
                      <xs:sequence>
                        <xs:element name="DESIGN_DESCRIPTION" type="xs:string">
                          <xs:annotation>
                          <xs:documentation>
                              More details about the setup and goals of the experiment as supplied by the Investigator.
                          </xs:documentation>
                          </xs:annotation>                                              
                        </xs:element>

 
                         
                        <xs:element name="SAMPLE_DESCRIPTOR" minOccurs="1" maxOccurs="1">
                          <xs:annotation>
                            <xs:documentation>
                              Pick a sample to associate this experiment with.  
                              The sample may be an individual or a pool, depending on how it is specified.
                            </xs:documentation>
                          </xs:annotation>  
                          <xs:complexType>
                              <xs:sequence>
                              <xs:element name="POOL" minOccurs="0" maxOccurs="1"> 
                                  <xs:annotation>
                                      <xs:documentation>
                                          Identifies a list of group/pool/multiplex sample members.  This implies that
                                          this sample record is a group, pool, or multiplex, but is continues to receive
                                          its own accession and can be referenced by an experiment.  By default if
                                          no match to any of the listed members can be determined, then the default
                                          sampel reference is used.
                                      </xs:documentation>
                                  </xs:annotation>
                                  <xs:complexType>
                                      <xs:sequence>
                                          <xs:element name="MEMBER" minOccurs="1" maxOccurs="unbounded">
                                              <xs:complexType>
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Impementation of lookup table between Sample Pool member and identified read_group_tags for a given READ_LABEL
                                                      </xs:documentation>
                                                  </xs:annotation>
                                                  <xs:sequence minOccurs="1" maxOccurs="unbounded">
                                                      <xs:element name=\'READ_LABEL\'>
                                                          <xs:complexType>
                                                              <xs:simpleContent>
                                                                  <xs:extension base="xs:string">
                                                                      <xs:attribute name="read_group_tag" type="xs:string">
                                                                          <xs:annotation>
                                                                              <xs:documentation>
                                                                                  Assignment of read_group_tag to decoded read
                                                                              </xs:documentation>
                                                                          </xs:annotation>
                                                                      </xs:attribute>
                                                                  </xs:extension>
                                                              </xs:simpleContent>
                                                          </xs:complexType>
                                                              
                                                      </xs:element>
                                                      
                                                  </xs:sequence>
                                                  <xs:attributeGroup ref="RefNameGroup">
                                                      <xs:annotation>
                                                          <xs:documentation>
                                                              Reference to SRA sample
                                                          </xs:documentation>
                                                      </xs:annotation>
                                                  </xs:attributeGroup>
                                                 
                                                  <xs:attribute name="member_name" type="xs:string" use="optional">
                                                      <xs:annotation>
                                                          <xs:documentation>
                                                             Label a sample within a scope of the pool 
                                                          </xs:documentation>
                                                      </xs:annotation>
                                                  </xs:attribute>
                                                  <xs:attribute name="proportion" type="xs:float" use="optional">
                                                      <xs:annotation>
                                                          <xs:documentation>
                                                              Proportion of this sample (in percent) that was included in sample pool.
                                                          </xs:documentation>
                                                      </xs:annotation>
                                                  </xs:attribute>
                                              </xs:complexType>
                                              
                                          </xs:element>
                                      </xs:sequence>                              
                                  </xs:complexType>
                              </xs:element>                      
                              </xs:sequence>
                              <xs:attributeGroup ref="RefNameGroup"/>                  
                          </xs:complexType>
                        </xs:element>

                        <xs:element name="LIBRARY_DESCRIPTOR">
                          <xs:annotation>
                              <xs:documentation>
                                  The LIBRARY_DESCRIPTOR specifies the origin of the material being sequenced and any treatments that the 
                                  material might have undergone that affect the sequencing result.  This specification is needed even if the platform
                                  does not require a library construction step per se.
                              </xs:documentation>
                          </xs:annotation>  
                          <xs:complexType>
                              <xs:sequence>
                                  <xs:element name="LIBRARY_NAME" type="xs:string">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The submitter\'s name for this library.
                                          </xs:documentation>
                                      </xs:annotation>  
                                  </xs:element>
                                  <xs:element name="LIBRARY_STRATEGY" maxOccurs="1" minOccurs="0">
                                      <xs:annotation>
                                          <xs:documentation>
                                              Sequencing technique intended for this library.
                                          </xs:documentation>
                                      </xs:annotation> 
                                      <xs:simpleType>
                                          <xs:restriction base="xs:string">
                                              <xs:enumeration value="WGS">                                                                    
                                              <xs:annotation>
                                                  <xs:documentation>
                                                      Whole genome shotgun.
                                                  </xs:documentation>
                                              </xs:annotation> 
                                              </xs:enumeration>    
                                              <xs:enumeration value="WCS" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Whole chromosome (or other replicon) shotgun.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="CLONE" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Genomic clone based (hierarchical) sequencing.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="POOLCLONE" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Shotgun of pooled clones (usually BACs and Fosmids).
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="AMPLICON" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                         Sequencing of overlapping or distinct PCR or RT-PCR products.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="BARCODE" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          DEPRECATED.  Sequencing of overlapping or distinct  products that have been
                                                          tagged with a short identifying sequence (barcode).  Each sequence read can
                                                          therefore be assigned to an individual product.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="CLONEEND" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Clone end (5', 3', or both) sequencing.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="FINISHING" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Sequencing intended to finish (close) gaps in existing coverage.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="ChIP-Seq">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Direct sequencing of chromatin immunoprecipitates.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                 
                                              </xs:enumeration>
                                              <xs:enumeration value="MNase-Seq">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Direct sequencing following MNase digestion.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration>
                                              <xs:enumeration value="DNase-Hypersensitivity">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Sequencing of hypersensitive sites, or segments of open chromatin that are more readily cleaved by DNaseI.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration>
                                              <xs:enumeration value="Bisulfite-Seq">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Sequencing following treatment of DNA with bisulfite to convert cytosine residues to uracil depending on methylation status.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration>
                                              <xs:enumeration value="EST">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Single pass sequencing of cDNA templates
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration> 
                                              <xs:enumeration value="FL-cDNA">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Full-length sequencing of cDNA templates
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration>
                                              <xs:enumeration value="CTS">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Concatenated Tag Sequencing
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                               
                                              </xs:enumeration>  
                                              <xs:enumeration value="OTHER">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Library strategy not listed.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                                                                                
                                              </xs:enumeration>
                                          </xs:restriction>
                                      </xs:simpleType>
                                  </xs:element>
                                  <xs:element name="LIBRARY_SOURCE" >
                                      <xs:annotation>
                                          <xs:documentation>
                                              The LIBRARY_SOURCE specifies the type of source material that is being sequenced.
                                          </xs:documentation>
                                      </xs:annotation>
                                      <xs:simpleType>
                                          <xs:restriction base="xs:string">
                                              <xs:enumeration value="GENOMIC"/>
                                              <xs:enumeration value="NON GENOMIC" />
                                              <xs:enumeration value="SYNTHETIC"/>
                                              <xs:enumeration value="VIRAL RNA"/>                                                               
                                              <xs:enumeration value="OTHER"/>
                                          </xs:restriction>
                                      </xs:simpleType>
                                  </xs:element>
 
                                  <xs:element name="LIBRARY_SELECTION">
                                      <xs:annotation>
                                          <xs:documentation>
                                              Whether any method was used to select and/or enrich the material being sequenced.     
                                          </xs:documentation>
                                      </xs:annotation> 
                                      <xs:simpleType>
                                          <xs:restriction base="xs:string">
                                              <xs:enumeration value="RANDOM">
                                              <xs:annotation>
                                                  <xs:documentation>Random shearing only.</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="PCR">
                                              <xs:annotation>
                                                  <xs:documentation>Source material was selected by designed primers.</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="RANDOM PCR">
                                              <xs:annotation>
                                                  <xs:documentation>Source material was selected by randomly generated primers.</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="RT-PCR">                    
                                              <xs:annotation>
                                                  <xs:documentation>Source material was selected by reverse transcription PCR</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="HMPR">
                                              <xs:annotation>
                                                  <xs:documentation>Hypo-methylated partial restriction digest</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="MF" >
                                              <xs:annotation>
                                                  <xs:documentation>Methyl Filtrated</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="CF-S">
                                              <xs:annotation>
                                                  <xs:documentation>Cot-filtered single/low-copy genomic DNA</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="CF-M">
                                              <xs:annotation>
                                                  <xs:documentation>Cot-filtered moderately repetitive genomic DNA</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="CF-H">
                                              <xs:annotation>
                                                  <xs:documentation>Cot-filtered highly repetitive genomic DNA</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="CF-T">
                                              <xs:annotation>
                                                  <xs:documentation>Cot-filtered theoretical single-copy genomic DNA</xs:documentation>
                                              </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="MSLL">
                                                  <xs:annotation>
                                                      <xs:documentation>Methylation Spanning Linking Library</xs:documentation>
                                                  </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="cDNA">
                                                  <xs:annotation>
                                                      <xs:documentation>complementary DNA</xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>  
                                              <xs:enumeration value="ChIP">
                                                  <xs:annotation>
                                                      <xs:documentation>Chromatin immunoprecipitation</xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>
                                              <xs:enumeration value="MNase">
                                                  <xs:annotation>
                                                      <xs:documentation>Micrococcal Nuclease (MNase) digestion</xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>                                                               
                                              <xs:enumeration value="DNAse">
                                                  <xs:annotation>
                                                      <xs:documentation>Deoxyribonuclease (MNase) digestion</xs:documentation>
                                                  </xs:annotation>                                                                   
                                              </xs:enumeration>                                                               
                                              <xs:enumeration value="Hybrid Selection">
                                                  <xs:annotation>
                                                      <xs:documentation>Selection by hybridization in array or solution.</xs:documentation>
                                                  </xs:annotation>
                                              </xs:enumeration>
                                              <xs:enumeration value="Reduced Representation">
                                                  <xs:annotation>
                                                      <xs:documentation>Reproducible genomic subsets, often generated by restriction fragment size selection, 
                                                          containing a manageable number of loci to facilitate re-sampling.
                                                      </xs:documentation>
                                                  </xs:annotation>
                                              </xs:enumeration>                                           
                                              <xs:enumeration value="other"/>                                                                
                                              <xs:enumeration value="unspecified"/>
                                          </xs:restriction>
                                      </xs:simpleType>
                                  </xs:element>
                                  <xs:element name="LIBRARY_LAYOUT">
                                      <xs:annotation>
                                          <xs:documentation>
                                              LIBRARY_LAYOUT specifies whether to expect single, paired, or other configuration of reads.  
                                              In the case of paired reads, information about the relative distance and orientation is specified.
                                          </xs:documentation>
                                      </xs:annotation>                                                        
                                      <xs:complexType>
                                          <xs:choice>
                                              <xs:element name="SINGLE">
                                                <xs:complexType>
                                                 <xs:annotation>
                                                  <xs:documentation>
                                                    Reads are unpaired (usual case).
                                                  </xs:documentation>
                                                 </xs:annotation>                                                        
                                                </xs:complexType>
                                              </xs:element>
                                              <xs:element name="PAIRED">
                                                  <xs:complexType>
                                                    <xs:attribute name="ORIENTATION" type="xs:string">
                                                      <xs:annotation>
                                                        <xs:documentation>
                                                        </xs:documentation>
                                                      </xs:annotation>                                                        
                                                    </xs:attribute>
                                                    <xs:attribute name="NOMINAL_LENGTH" type="xs:nonNegativeInteger">
                                                      <xs:annotation>
                                                        <xs:documentation>
                                                        </xs:documentation>
                                                      </xs:annotation>                                                        
                                                    </xs:attribute>
                                                    <xs:attribute name="NOMINAL_SDEV" type="xs:double">
                                                      <xs:annotation>
                                                        <xs:documentation>
                                                        </xs:documentation>
                                                      </xs:annotation>                                                        
                                                    </xs:attribute>
                                                  </xs:complexType>
                                              </xs:element>
                                          </xs:choice>
                                      </xs:complexType>
                                  </xs:element>
                                  
                                  <xs:element name="LIBRARY_CONSTRUCTION_PROTOCOL" type="xs:string" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              Free form text describing the protocol by which the sequencing library was constructed.                             
                                          </xs:documentation>
                                      </xs:annotation>                                                        
                                  </xs:element>
                              </xs:sequence>
                          </xs:complexType>
                      </xs:element>
                      <xs:element name="SPOT_DESCRIPTOR">
                          <xs:annotation>
                              <xs:documentation>
                                  The SPOT_DESCRIPTOR specifies how to decode the individual reads of interest from the 
                                  monolithic spot sequence.  The spot descriptor contains aspects of the experimental design, 
                                  platform, and processing information.  There will be two methods of specification: one 
                                  will be an index into a table of typical decodings, the other being an exact specification.                                      
                              </xs:documentation>
                          </xs:annotation>
                          <xs:complexType>
                              <xs:choice>
                                  <xs:element name="SPOT_DECODE_METHOD" type="xs:unsignedInt" >
                                      <xs:annotation>
                                          <xs:documentation>
                                             DEPRECATED.                                    
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="SPOT_DECODE_SPEC" >
                                      <xs:complexType   >
                                          <xs:sequence>
                                              <xs:element name="NUMBER_OF_READS_PER_SPOT" type="xs:unsignedInt" minOccurs="0" maxOccurs="1" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          DEPRECATED.  Number of tags (reads) per spot.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                   
                                              </xs:element>
                                              <xs:element name="SPOT_LENGTH" type="xs:unsignedInt" minOccurs="0" maxOccurs="1" >
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Expected number of base calls or cycles per spot.
                                                      </xs:documentation>
                                                  </xs:annotation>                                                    
                                              </xs:element>
                                              <xs:element name="ADAPTER_SPEC" maxOccurs="1" minOccurs="0" nillable="false" type="xs:string">
                                                  <xs:annotation>
                                                      <xs:documentation>
                                                          Some technologies will require knowledge of the sequencing adapter or the last base of the adapter in order to decode the spot.
                                                      </xs:documentation>
                                                  </xs:annotation>   
                                              </xs:element>
                                              <xs:element name="READ_SPEC" minOccurs="1" maxOccurs="unbounded" >
                                                  <xs:complexType>                              
                                                      <xs:sequence>
                                                          <xs:element name="READ_INDEX" type="xs:nonNegativeInteger" nillable="false">
                                                              <xs:annotation>
                                                                  <xs:documentation>READ_INDEX starts at 0 and is incrementally increased for each sequential READ_SPEC within a SPOT_DECODE_SPEC</xs:documentation>
                                                              </xs:annotation> 
                                                          </xs:element>
                                                          <xs:element name="READ_LABEL" type="xs:string" minOccurs="0" maxOccurs="1">
                                                              <xs:annotation>
                                                                  <xs:documentation>READ_LABEL is a name for this tag, and can be used to on output to determine read name, for example F or R.</xs:documentation>
                                                              </xs:annotation> 
                                                          </xs:element>                                  
                                                          <xs:element name="READ_CLASS" >
                                                              <xs:simpleType>
                                                                  <xs:restriction base = "xs:string">
                                                                      <xs:enumeration value = "Application Read"/>
                                                                      <xs:enumeration value = "Technical Read"/>
                                                                  </xs:restriction>                                   
                                                              </xs:simpleType>   
                                                          </xs:element>
                                                          <xs:element name="READ_TYPE"  default="Forward">
                                                              <xs:simpleType>
                                                                  <xs:restriction base = "xs:string">
                                                                      <xs:enumeration value = "Forward"/>
                                                                      <xs:enumeration value = "Reverse"/>
                                                                      <xs:enumeration value = "Adapter"/>                                   
                                                                      <xs:enumeration value = "Primer"/>
                                                                      <xs:enumeration value = "Linker"/>
                                                                      <xs:enumeration value = "BarCode"/>                                   
                                                                      <xs:enumeration value = "Other"/>
                                                                  </xs:restriction>                                   
                                                              </xs:simpleType>   
                                                          </xs:element>
                                                          <xs:choice>
                                                              <xs:annotation>
                                                                  <xs:documentation>
                                                                      There are various methods to ordering the reads on the spot.
                                                                  </xs:documentation>
                                                              </xs:annotation>   
                                                              <xs:element name="RELATIVE_ORDER" >
                                                                  <xs:annotation>
                                                                      <xs:documentation>
                                                                          The read is located beginning at the offset or cycle relative to another read.  
                                                                          This choice is appropriate for example when specifying a read
                                                                          that follows a variable length expected sequence(s).
                                                                      </xs:documentation>
                                                                  </xs:annotation>  
                                                                  <xs:complexType>
                                                                      <xs:attribute name="follows_read_index" type="xs:nonNegativeInteger" use="optional">
                                                                          <xs:annotation>
                                                                              <xs:documentation>
                                                                                  Specify the read index that precedes this read.
                                                                              </xs:documentation>
                                                                          </xs:annotation>  
                                                                      </xs:attribute>
                                                                      <xs:attribute name="precedes_read_index" type="xs:nonNegativeInteger" use="optional">
                                                                          <xs:annotation>
                                                                              <xs:documentation>
                                                                                  Specify the read index that follows this read.
                                                                              </xs:documentation>
                                                                          </xs:annotation>  
                                                                      </xs:attribute>
                                                                  </xs:complexType>
                                                              </xs:element>
                                                              <xs:element name="BASE_COORD" type="xs:integer" >
                                                              <xs:annotation>
                                                                  <xs:documentation>
                                                                      The location of the read start in terms of base count (1 is beginning of spot).
                                                                  </xs:documentation>
                                                              </xs:annotation>   
                                                              </xs:element>
                                                              <xs:element name="CYCLE_COORD" type="xs:integer" >
                                                              <xs:annotation>
                                                                  <xs:documentation>
                                                                      The location of the read start in terms of cycle count (1 is beginning of spot).
                                                                  </xs:documentation>
                                                              </xs:annotation>
                                                              </xs:element>
                                                              
                                                              <xs:element name="EXPECTED_BASECALL" type="xs:string" >
                                                              <xs:annotation>
                                                                  <xs:documentation>
                                                                      An expected basecall for a current read. Read will be zero-length if basecall is not present.
                                                                  </xs:documentation>
                                                              </xs:annotation> 
                                                                  
                                                              </xs:element>
                                                                  <xs:element name="EXPECTED_BASECALL_TABLE">
                                                                      <xs:annotation>
                                                                          <xs:documentation>
                                                                              A set of choices of expected basecalls for a current read. Read will be zero-length if none is found.
                                                                          </xs:documentation>
                                                                      </xs:annotation> 
                                                                      <xs:complexType>
                                                                          <xs:sequence minOccurs="1" maxOccurs="unbounded">
                                                                              <xs:element name="BASECALL">
                                                                                  <xs:annotation>
                                                                                      <xs:documentation>
                                                                                          Element\'s body contains a basecall, attribute provide description of this read meaning as well as matching rules.
                                                                                      </xs:documentation>
                                                                                  </xs:annotation> 
                                                                                  <xs:complexType>
                                                                                      <xs:simpleContent>
                                                                                          <xs:extension base  ="xs:string">
                                                                                              <xs:attribute name="read_group_tag" type="xs:string" use="optional">
                                                                                                  <xs:annotation>
                                                                                                      <xs:documentation>
                                                                                                          When match occurs, the read will be tagged with this group membership
                                                                                                      </xs:documentation>
                                                                                                  </xs:annotation>                                                                                          
                                                                                              </xs:attribute>
                                                                                              <xs:attribute name="min_match" type="xs:nonNegativeInteger" use="optional">
                                                                                                  <xs:annotation>
                                                                                                      <xs:documentation>
                                                                                                          Minimum number of matches to trigger identification.
                                                                                                      </xs:documentation>
                                                                                                  </xs:annotation>                                                                                          
                                                                                              </xs:attribute>
                                                                                              <xs:attribute name="max_mismatch" type="xs:nonNegativeInteger" use="optional">
                                                                                                  <xs:annotation>
                                                                                                      <xs:documentation>
                                                                                                          Maximum number of mismatches 
                                                                                                      </xs:documentation>
                                                                                                  </xs:annotation>                                                                                          
                                                                                              </xs:attribute>                                                                                              
                                                                                              <xs:attribute name="match_edge">
                                                                                                  <xs:annotation>
                                                                                                      <xs:documentation>
                                                                                                          Where the match should occur. Changes the rules on how min_match and max_mismatch are counted.                                                                                                          
                                                                                                      </xs:documentation>
                                                                                                  </xs:annotation>                              
                                                                                                  <xs:simpleType>
                                                                                                      <xs:restriction base="xs:string">
                                                                                                          <xs:enumeration value="full">
                                                                                                              <xs:annotation>
                                                                                                                  <xs:documentation>
                                                                                                                      Only @max_mismatch influences matching process                                                                                                          
                                                                                                                  </xs:documentation>
                                                                                                              </xs:annotation>
                                                                                                          </xs:enumeration>
                                                                                                          <xs:enumeration value="start">
                                                                                                              <xs:annotation>
                                                                                                                  <xs:documentation>
                                                                                                                      Both matches and mismatches are counted. 
                                                                                                                      When @max_mismatch is exceeded - it is not a match.
                                                                                                                      When @min_match is reached - match is declared.                                                                                                                                                                                                                           
                                                                                                                  </xs:documentation>
                                                                                                              </xs:annotation>                                                                                                              
                                                                                                          </xs:enumeration>
                                                                                                          <xs:enumeration value="end">
                                                                                                              <xs:annotation>
                                                                                                                  <xs:documentation>
                                                                                                                      Both matches and mismatches are counted. 
                                                                                                                      When @max_mismatch is exceeded - it is not a match.
                                                                                                                      When @min_match is reached - match is declared.                                                                                                                                                                                                                           
                                                                                                                  </xs:documentation>
                                                                                                              </xs:annotation>                                                                                                                      
                                                                                                          </xs:enumeration>
                                                                                                      </xs:restriction>
                                                                                                  </xs:simpleType>
                                                                                              </xs:attribute>                                                                                                                                                                                            
                                                                                          </xs:extension>                                                                                          
                                                                                      </xs:simpleContent>                                                                                                                                                                                                                                                                                                                                                   
                                                                                  </xs:complexType>                                                                                                                                                                   
                                                                              </xs:element>
                                                                          </xs:sequence>
                                                                      </xs:complexType>                                                                      
                                                                  </xs:element>
                                                                  </xs:choice>

                                                          
                                                      </xs:sequence>
                                                  </xs:complexType>
                                              </xs:element>

                      </xs:sequence>
                                      </xs:complexType>
                                  </xs:element>
                              </xs:choice>

                          </xs:complexType>
                      </xs:element>
                      </xs:sequence>
                        </xs:complexType>
                  </xs:element>
      

              <xs:element name="PLATFORM" maxOccurs="1" minOccurs="1">
                  <xs:annotation>
                    <xs:documentation>
                      The PLATFORM record selects which sequencing platform and platform-specific runtime parameters.  
                      This will be determined by the Center.
                    </xs:documentation>
                  </xs:annotation>  
                  <xs:complexType>
                    <xs:choice>
                      <xs:element name="LS454">
                          <xs:annotation>
                              <xs:documentation>
                                  454 technology use 1-color sequential flows 
                              </xs:documentation>
                          </xs:annotation>
                          <xs:complexType>
                              <xs:all>
                                  <xs:element name="INSTRUMENT_MODEL" maxOccurs="1" minOccurs="1">
                                      <xs:simpleType>
                                          <xs:restriction base = "xs:string">
                                              <xs:enumeration value="454 GS"/>
                                              <xs:enumeration value="454 GS 20"/>
                                              <xs:enumeration value="454 GS FLX"/>
                                              <xs:enumeration value="454 Titanium"/>
                                              <xs:enumeration value = "GS 20">
                                                <xs:annotation><xs:documentation>DEPRECATED</xs:documentation></xs:annotation></xs:enumeration>
                                              <xs:enumeration value = "GS FLX">
                                                <xs:annotation><xs:documentation>DEPRECATED</xs:documentation></xs:annotation></xs:enumeration>
                                              <xs:enumeration value = "unspecified"/>                   
                                          </xs:restriction>                                   
                                      </xs:simpleType>  
                                  </xs:element>
                                  <xs:element name="KEY_SEQUENCE" type="xs:string" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The first bases that are expected to be produced by the challenge bases.  
                                              This is optional in the schema now but will be required by business rules and future schema versions.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="FLOW_SEQUENCE" type="xs:string"  minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The fixed sequence of challenge bases that flow across the picotiter plate.  
                                              This is optional in the schema now but will be required by business rules and future schema versions.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="FLOW_COUNT" type="xs:positiveInteger" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The number of flows of challenge bases.  This is a constraint on maximum read length, but not equivalent.
                                              This is optional in the schema now but will be required by business rules and future schema versions.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                              </xs:all>
                          </xs:complexType>
                      </xs:element>
                      <xs:element name="ILLUMINA">
                          <xs:annotation>
                              <xs:documentation>
                                  Illumina is 4-channel flowgram with 1-to-1 mapping between basecalls and flows
                              </xs:documentation>
                          </xs:annotation>
                          <xs:complexType>
                              <xs:sequence>
                                  <xs:element name="INSTRUMENT_MODEL" maxOccurs="1" minOccurs="1">
                                      <xs:simpleType>
                                          <xs:restriction base = "xs:string">
                                            <xs:enumeration value="Solexa 1G Genome Analyzer">
                                              <xs:annotation><xs:documentation>DEPRECATED</xs:documentation></xs:annotation></xs:enumeration>
                                            <xs:enumeration value="Illumina Genome Analyzer"/>
                                            <xs:enumeration value="Illumina Genome Analyzer II"/>
                                            <xs:enumeration value = "unspecified"/>                   
                                          </xs:restriction>                                   
                                      </xs:simpleType>  
                                  </xs:element>
                                  <xs:element name="CYCLE_SEQUENCE" type="xs:string" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              DEPRECATED.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="CYCLE_COUNT" type="xs:positiveInteger" minOccurs="1" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              DEPRECATED, use SEQUENCE_LENGTH instead.  The fixed number of bases  in each raw sequence, including both mate pairs and any technical reads.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="SEQUENCE_LENGTH" type="xs:positiveInteger" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The fixed number of bases expected in each raw sequence, including both mate pairs and any technical reads.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                              </xs:sequence>
                          </xs:complexType>
                      </xs:element>
                      <xs:element name="HELICOS">
                          <xs:annotation>
                              <xs:documentation>
                                  Helicos is similar to 454 technology - uses 1-color sequential flows   
                              </xs:documentation>
                          </xs:annotation>
                          <xs:complexType>
                              <xs:all>
                                  <xs:element name="INSTRUMENT_MODEL" maxOccurs="1" minOccurs="1">
                                      <xs:simpleType>
                                          <xs:restriction base = "xs:string">
                                            <xs:enumeration value="Helicos HeliScope"/>
                                            <xs:enumeration value = "unspecified"/>                   
                                          </xs:restriction>                                   
                                      </xs:simpleType>  
                                  </xs:element>
                                  <xs:element name="FLOW_SEQUENCE" type="xs:string" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The fixed sequence of challenge bases that flow across the flowcell. 
                                              This is optional in the schema now but will be required by business rules and future schema versions.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                  <xs:element name="FLOW_COUNT" type="xs:positiveInteger" minOccurs="0" maxOccurs="1">
                                      <xs:annotation>
                                          <xs:documentation>
                                              The number of flows of challenge bases.  This is a constraint on maximum read length, but not equivalent. 
                                              This is optional in the schema now but will be required by business rules and future schema versions.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                              </xs:all>
                          </xs:complexType>
                      </xs:element>
                      <xs:element name="ABI_SOLID">
                          <xs:annotation>
                              <xs:documentation>
                                  ABI is 4-channel flowgram with 1-to-1 mapping between basecalls and flows
                              </xs:documentation>
                          </xs:annotation>
                          <xs:complexType>
                              <xs:sequence>
                                  <xs:element name="INSTRUMENT_MODEL" maxOccurs="1" minOccurs="1">
                                      <xs:simpleType>
                                          <xs:restriction base = "xs:string">
                                            <xs:enumeration value="AB SOLiD System"/>
                                            <xs:enumeration value="AB SOLiD System 2.0"/>
                                            <xs:enumeration value="AB SOLiD System 3.0"/>
                                            <xs:enumeration value = "unspecified"/>                   
                                          </xs:restriction>                                   
                                      </xs:simpleType>  
                                  </xs:element>
                                  <xs:element name="COLOR_MATRIX"  minOccurs="0"  maxOccurs="1">
                                      <xs:complexType>
                                          <xs:sequence >
                                              <xs:element name="COLOR"  minOccurs="1" maxOccurs="unbounded">
                                                  <xs:complexType mixed="true">
                                                      <xs:simpleContent>
                                                          <xs:extension base="xs:string">
                                                              <xs:attribute name="dibase" type="xs:string" />
                                                          </xs:extension>
                                                      </xs:simpleContent>
                                                  </xs:complexType>
                                              </xs:element>
                                          </xs:sequence>
                                      </xs:complexType>
                                  </xs:element>
                                  <xs:element name="COLOR_MATRIX_CODE" type="xs:string" minOccurs="0" maxOccurs="1">
                                  </xs:element>
                                  <xs:element name="CYCLE_COUNT" type="xs:positiveInteger">
                                      <xs:annotation>
                                          <xs:documentation>
                                              DEPRECATED.  Use SEQUENCE_LENGTH instead.
                                          </xs:documentation>
                                      </xs:annotation>
                                  </xs:element>
                                      <xs:element name="SEQUENCE_LENGTH" type="xs:positiveInteger" minOccurs="0" maxOccurs="1">
                                          <xs:annotation>
                                              <xs:documentation>
                                                  The fixed number of bases expected in each raw sequence, including both mate pairs and any technical reads.
                                                  This is optional in the schema now but will be required by business rules and future schema versions.
                                              </xs:documentation>
                                          </xs:annotation>
                                      </xs:element>
                              </xs:sequence>
                          </xs:complexType>
                      </xs:element>
                      <xs:element name="COMPLETE_GENOMICS">
                          <xs:annotation>
                              <xs:documentation>
                                  Placeholder for CompleteGenomics platform type.   
                              </xs:documentation>
                          </xs:annotation>
                      </xs:element>
                        <xs:element name="PACBIO_SMRT">
                            <xs:annotation>
                                <xs:documentation>
                                    Placeholder for PacificBiosciences platform type.   
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:choice>
                  </xs:complexType>
              </xs:element>
  
                <xs:element name="PROCESSING" maxOccurs="1" minOccurs="1">
                  <xs:annotation>
                    <xs:documentation> 
                      The PROCESSING block specifies how sequencing is processed from image files.  
                      This information should be  vendor determined based on the PLATFORM selection.
                      </xs:documentation>
                    </xs:annotation>  
                  <xs:complexType>
                    <xs:sequence>
                      <xs:element name="BASE_CALLS" maxOccurs="1" minOccurs="1">
                          <xs:complexType>
                              <xs:all>
                                  <xs:element name="SEQUENCE_SPACE">
                                                        <xs:simpleType>
                                                            <xs:restriction base = "xs:string">
                                                                <xs:enumeration value = "Base Space"/>
                                                                <xs:enumeration value = "Color Space"/>
                                                            </xs:restriction>                                   
                                                        </xs:simpleType>
                                                        </xs:element>
                                                        <xs:element name="BASE_CALLER" maxOccurs="1" minOccurs="1" type="xs:string">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Name and version of the base or color calling software.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:element>
                                                </xs:all>
                                            </xs:complexType>
                                        </xs:element>
                                        <xs:element name="QUALITY_SCORES" maxOccurs="unbounded" minOccurs="1">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    In the future there will be only one instance allowed of the QUALITY_SCORES spec.
                                                </xs:documentation>
                                            </xs:annotation>                                           
                                            <xs:complexType>
                                                <xs:all>
                                                    <xs:element name="QUALITY_SCORER" maxOccurs="1" minOccurs="1" type="xs:string">              
                                                    <xs:annotation>
                                                        <xs:documentation>
                                                            Name and version of the quality scoring software.
                                                        </xs:documentation>
                                                    </xs:annotation>
                                                    </xs:element>
                                                    <xs:element name="NUMBER_OF_LEVELS" maxOccurs="1" minOccurs="1" type="xs:int">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                DEPRECATED.  Number of distinct values possible with this scoring system.
                                                            </xs:documentation>
                                                        </xs:annotation>     
                                                    </xs:element>
                                                    <xs:element name="MULTIPLIER" maxOccurs="1" minOccurs="1" type="xs:double">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                DEPRECATED.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:element>
                                                </xs:all>
                                                <xs:attribute name="qtype" use="optional">
                                                   <xs:simpleType>
                                                     <xs:restriction base = "xs:string">                             
                                                             <xs:enumeration value = "phred">
                                                                 <xs:annotation>
                                                                     <xs:documentation>
                                                                         The quality score is expressed as a probability of error in log form:
                                                                         -10 log(1/p) where p is the probability of error, with value range 0..63,
                                                                         0 meaning no base call.
                                                                     </xs:documentation>
                                                                 </xs:annotation>
                                                             </xs:enumeration>
                                                             <xs:enumeration value = "other">
                                                                 <xs:annotation>
                                                                     <xs:documentation>
                                                                         Another quality scoring system is used.  Please contact NCBI with details so that correct conversion can take place.
                                                                     </xs:documentation>
                                                                 </xs:annotation>
                                                             </xs:enumeration>
                                                         </xs:restriction>
                                                   </xs:simpleType>
                                                </xs:attribute>
                                            </xs:complexType>
                                        </xs:element>
 
                </xs:sequence>
              </xs:complexType>
            </xs:element>
	    <xs:element name="EXPERIMENT_LINKS" minOccurs="0" maxOccurs="1">
		  <xs:annotation>
		      <xs:documentation>
			  Links to resources related to this experiment or experiment set (publication, datasets, online databases).
		      </xs:documentation>
		  </xs:annotation>
		  <xs:complexType>
		    <xs:sequence minOccurs="1" maxOccurs="unbounded">
		      <xs:element name="EXPERIMENT_LINK" type="LinkType"/>
		    </xs:sequence>
		  </xs:complexType>
	        </xs:element>
                <xs:element name="EXPERIMENT_ATTRIBUTES" minOccurs="0" maxOccurs="1">
               	  <xs:annotation>
                    <xs:documentation>
                       Properties and attributes of the experiment.  These can be entered as free-form 
                       tag-value pairs. 
                    </xs:documentation>
                  </xs:annotation>
                  <xs:complexType>
                    <xs:sequence maxOccurs="unbounded" minOccurs="1">
                      <xs:element name="EXPERIMENT_ATTRIBUTE" type="AttributeType"/>
                    </xs:sequence>
                  </xs:complexType>
                </xs:element>
                </xs:sequence>
        <xs:attributeGroup ref="NameGroup"/>
             
     
                  
              <xs:attribute name="expected_number_runs" use="optional" type="xs:positiveInteger">
                <xs:annotation>
                  <xs:documentation>
                      DEPRECATED.  Number of runs expected to be submitted  for this experiment. 
                  </xs:documentation>
                </xs:annotation>
              </xs:attribute>
              <xs:attribute name="expected_number_spots" use="optional" type="xs:positiveInteger">
                  <xs:annotation>
                     <xs:documentation>
                       DEPRECATED. Number of spots expected to be submitted  for this experiment. 
                     </xs:documentation>
                  </xs:annotation>
              </xs:attribute>
              <xs:attribute name="expected_number_reads" use="optional" type="xs:positiveInteger">
                  <xs:annotation>
                    <xs:documentation>
                      DEPRECATED. Number of reads expected to be submitted  for this experiment. 
                    </xs:documentation>
                  </xs:annotation>
              </xs:attribute>

</xs:complexType>
    
    
<xs:element name="EXPERIMENT_SET" >
  <xs:annotation>
    <xs:documentation>
      An EXPERMENT_SET is a container for a set of experiments and a common namespace.
    </xs:documentation>
  </xs:annotation>  
  <xs:complexType>
    <xs:sequence minOccurs="1" maxOccurs="unbounded">
        <xs:element name="EXPERIMENT" type="ExperimentType"/>        
    </xs:sequence>
  </xs:complexType>
</xs:element> 
 
<xs:element name="EXPERIMENT" type="ExperimentType"/>        

</xs:schema>
'''

# SRA Production release 1.1
run_xsd = '''<?xml version="1.0" encoding="UTF-8"?>
<!-- NCBI Short Read Archive resource Run (SRR) object XML specification -->
<!-- $Id$ -->
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <!-- BEGIN COMMON BLOCK            TO BE NORMALIZED IN THE FUTURE -->
    <xs:attributeGroup name="NameGroup">
        <xs:attribute name="alias" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Submitter designated name of the SRA document of this type.  At minimum alias should
                    be unique throughout the submission of this document type.  If center_name is specified, the name should
                    be unique in all submissions from that center of this document type.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="center_name" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Owner authority of this document and namespace for submitter\'s name of this document. 
                    If not provided, then the submitter is regarded as "Individual" and document resolution
                    can only happen within the submission.
                </xs:documentation>
            </xs:annotation>  
        </xs:attribute>
        <xs:attribute name="broker_name" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Broker authority of this document.  If not provided, then the broker is considered "direct".
                </xs:documentation>
            </xs:annotation>  
        </xs:attribute>
        <xs:attribute name="accession" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    The document\'s accession as assigned by the Home Archive.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
    </xs:attributeGroup>
    <xs:attributeGroup name="RefNameGroup">
        <xs:attribute name="refname" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Identifies a record by name that is known within the namespace defined by attribute "refcenter"
                    Use this field when referencing an object for which an accession has not yet been issued.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="refcenter" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    The namespace of the attribute "refname". When absent, the namespace is assumed to be the current submission.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
        <xs:attribute name="accession" type="xs:string" use="optional">
            <xs:annotation>
                <xs:documentation>
                    Identifies a record by its accession.  The scope of resolution is the entire Archive.
                </xs:documentation>
            </xs:annotation>
        </xs:attribute>
    </xs:attributeGroup>
    <xs:complexType name="AttributeType">
        <xs:annotation>
            <xs:documentation>
                Reusable attributes to encode tag-value pairs with optional units.
            </xs:documentation>
        </xs:annotation>
        <xs:all>
            <xs:element name="TAG" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                <xs:annotation>
                    <xs:documentation>
                        Name of the attribute.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
            <xs:element name="VALUE" type="xs:string" minOccurs="1" maxOccurs="1">
                <xs:annotation>
                    <xs:documentation>
                        Value of the attribute.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
            <xs:element name="UNITS" type="xs:string" minOccurs="0" maxOccurs="1">          
                <xs:annotation>
                    <xs:documentation>
                        Optional scientific units.
                    </xs:documentation>
                </xs:annotation>
            </xs:element>
        </xs:all>
    </xs:complexType>
    <xs:complexType name="SraLinkType">
        <xs:annotation>
            <xs:documentation>
                The SraLinkType mechanism encodes local references between SRA objects.  These 
                references are local to the Home Archive and during archival are resolved to
                exportable references suitable for mirroring between Archives.  SraLinks can be used
                as temporary place holders for pre-published or post-suppressed relationships that
                will not be mirrored between Archives.
            </xs:documentation>
        </xs:annotation>
        <xs:attribute name="sra_object_type" use="optional">
            <xs:annotation>
                <xs:documentation>
                    SRA link type.
                </xs:documentation>
            </xs:annotation>
            <xs:simpleType>
                <xs:restriction base="xs:string">
                    <xs:enumeration value="STUDY"/>
                    <xs:enumeration value="SAMPLE"/>
                    <xs:enumeration value="ANALYSIS"/>
                    <xs:enumeration value="EXPERIMENT"/>
                    <xs:enumeration value="RUN"/>            
                </xs:restriction>
            </xs:simpleType>
        </xs:attribute>
        <xs:attributeGroup ref="RefNameGroup"/>
        
    </xs:complexType>
    
    
    
    
    <xs:complexType name="LinkType">
        <xs:annotation>
            <xs:documentation>
                Reusable external links type to encode URL links, Entrez links, and db_xref links.
            </xs:documentation>
        </xs:annotation>
        <xs:choice>                   
            <xs:element name="SRA_LINK" type="SraLinkType"/>
            <xs:element name="URL_LINK">
                <xs:complexType>
                    <xs:all>
                        <xs:element name="LABEL" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    Text label to display for the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="URL" minOccurs="1" maxOccurs="1" type="xs:anyURI">
                            <xs:annotation>
                                <xs:documentation>
                                    The internet service link (file:, http:, ftp:, etc).
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:all>
                </xs:complexType>
            </xs:element>
            <xs:element name="XREF_LINK">
                <xs:complexType>
                    <xs:all>
                        <xs:element name="DB" type="xs:string" minOccurs="1" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    INSDC controlled vocabulary of permitted cross references.  Please see http://www.insdc.org/page.php?page=db_xref .
                                    For example, FLYBASE.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="ID" minOccurs="1" maxOccurs="1" type="xs:string">
                            <xs:annotation>
                                <xs:documentation>
                                    Accession in the referenced database.    For example,  FBtr0080008 (in FLYBASE).
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:element name="LABEL" type="xs:string" minOccurs="0" maxOccurs="1">                                                    
                            <xs:annotation>
                                <xs:documentation>
                                    Text label to display for the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:all>
                </xs:complexType>
            </xs:element>
            <xs:element name="ENTREZ_LINK">
                <xs:complexType>
                    <xs:sequence>
                        <xs:element name="DB" type="xs:string" minOccurs="1" maxOccurs="1">
                            <xs:annotation>
                                <xs:documentation>
                                    NCBI controlled vocabulary of permitted cross references.  Please see http://www.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi? .
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                        <xs:choice>
                            <xs:element name="ID" type="xs:nonNegativeInteger" minOccurs="1" maxOccurs="1">
                                <xs:annotation>
                                    <xs:documentation>
                                        Numeric record id meaningful to the NCBI Entrez system.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:element>
                            <xs:element name="QUERY" type="xs:string" minOccurs="1" maxOccurs="1">
                                <xs:annotation>
                                    <xs:documentation>
                                        Accession string meaningful to the NCBI Entrez system.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:element>
                        </xs:choice>
                        <xs:element name="LABEL" type="xs:string" minOccurs="0" maxOccurs="1">
                            <xs:annotation>
                                <xs:documentation>
                                    How to label the link.
                                </xs:documentation>
                            </xs:annotation>
                        </xs:element>
                    </xs:sequence>
                </xs:complexType>
            </xs:element>
        </xs:choice>
    </xs:complexType>
    
    <!--  END COMMON BLOCK -->
    
    <xs:complexType name="RunType">
        <xs:annotation>
            <xs:documentation>
                A Run contains the sequencing results from a particular run on a sequencing instrument.  The Run was done in fulfillment 
                of an Experiment that was designed into the Study.  
                One SRA run may consist of only a portion of the plate or slide if the instrument run was partitioned among distinct experiments 
                or between distinct RUNs.  By policy the Center should separate these data into distinct SRA Run objects.
            </xs:documentation>
        </xs:annotation>        
            <xs:sequence>
                <xs:element name="EXPERIMENT_REF" nillable="false"  maxOccurs="1" minOccurs="1">
                    <xs:annotation>
                        <xs:documentation> The EXPERIMENT_REF descriptor identifies the parent experiment
                            to which this run pertains.
                            The Experiment object contains all the mapping information needed to decode each spot and map application reads
                            to RUN objects.
                        </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                        <xs:attributeGroup ref="RefNameGroup"/>                                                                    
                    </xs:complexType>
                </xs:element> 
                <xs:sequence maxOccurs="unbounded" minOccurs="0" >
                    <xs:element name="DATA_BLOCK" nillable="false">           
                        <xs:annotation>                           
                             <xs:documentation>
                                 Convenience partition for processing large datasets.         
                             </xs:documentation>
                        </xs:annotation>
                        <xs:complexType>                        

                            <xs:sequence>
                               <xs:element name="FILES">
                                   <xs:annotation>
                                       <xs:documentation> Actual run data are contained in one of the files listed in the submission manifest. 
                                           Each data block is represented by one SRF file, one SFF file, one compressed fastq file, 
                                           or one compressed tar archive file.
                                       </xs:documentation>
                                   </xs:annotation>
                                  <xs:complexType>
                                      <xs:sequence maxOccurs="unbounded" minOccurs="1">
                                     <xs:element name="FILE">                                        
                                    <xs:complexType>
                                        <xs:sequence>
                                            <xs:element name="READ_LABEL"  type="xs:string" minOccurs="0" maxOccurs="unbounded">
                                                <xs:annotation>
                                                    <xs:documentation>
                                                        The READ_LABEL can associate a certain file to a certain read_label defined in the SPOT_DESCRIPTOR.
                                                        For example, the file "slide1_F3.csfasta" can be associated with read labeled F3 (the first forward read in a mate pair).
                                                        The FILE may contain data from multiple READ_LABELs.
                                                    </xs:documentation>
                                                </xs:annotation>
                                            </xs:element>
                                            
                                            <xs:element name="DATA_SERIES_LABEL" minOccurs="0" maxOccurs="unbounded">
                                                <xs:annotation>
                                                    <xs:documentation>
                                                        The DATA_SERIES_LABEL can associate a certain file to a certain data series (column) defined in the SPOT_DESCRIPTOR.
                                                        The FILE may contain data from multiple DATA_SERIES_LABELs.
                                                    </xs:documentation>
                                                </xs:annotation>
                                                <xs:simpleType>
                                                    <xs:restriction base = "xs:string">
                                                        <xs:enumeration value = "INSDC:read">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Base/color calls data series, one value per call.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:read_filter">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Read filter flags, one value per read.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:quality">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Quality scores, one score per base/color call.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:intensity">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Spot image intensity measurements, 4 values per call.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:signal">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Flow measurements, one value per flow.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:noise">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Spot image noise measurements, 4 values per call.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:position">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Position calls, one value per call.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:clip_quality_left">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Clip quality left, 1-based coordinate, inclusive.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:clip_quality_right">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Clip quality right, 1-based coordinate, inclusive.
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:readname">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Name of the individual read (spot sequence).
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                        <xs:enumeration value = "INSDC:read_seg">
                                                            <xs:annotation>
                                                                <xs:documentation>
                                                                    Vector of coordinates describing the location of a certain tag (spot sub-sequence).
                                                                </xs:documentation>
                                                            </xs:annotation>
                                                        </xs:enumeration>
                                                    </xs:restriction>
                                                </xs:simpleType>
                                            </xs:element>
                                        </xs:sequence>
                                        <xs:attribute name="filename" type="xs:string" use="required">
                                            <xs:annotation>
                                                <xs:documentation>The name or relative pathname of a run data file.</xs:documentation>
                                            </xs:annotation>
                                        </xs:attribute>
                                        <xs:attribute name="filetype" use="required">
                                            <xs:annotation>
                                                <xs:documentation> The run data file model.</xs:documentation>
                                            </xs:annotation>
                                            <xs:simpleType>
                                                <xs:restriction base="xs:string">
                                                    <xs:enumeration value="sra">                                                                 
                                                        <xs:annotation>
                                                            <xs:documentation>Sequence Read Archives native format in serialized (single file) form.</xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="srf">                                                                 
                                                      <xs:annotation>
                                                        <xs:documentation>Standard Short Read Format file (.srf), all platforms</xs:documentation>
                                                      </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="sff">
                                                      <xs:annotation>
                                                        <xs:documentation>454 Standard Flowgram Format file (.sff)</xs:documentation>
                                                      </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="fastq">               
                                                      <xs:annotation>
                                                          <xs:documentation>
                                                              Combined nucleotide/qualities sequence file in .fastq form.
                                                              Please see SRA File Formats Guide for definitions of the definition and restrictions on this form.
                                                          </xs:documentation>
                                                      </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="tab">               
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                Tab delimited text file used to deliver certain auxiliary data along with sequencing submissions (only needed for certain
                                                                use cases).   The first line is devoted to column headers.  Each column is dedicated to an INDSC
                                                                data series type.
                                                                Please see SRA File Formats Guide for definitions of the definition and restrictions on this form.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="_seq.txt, _prb.txt, _sig2.txt, _qhg.txt">               
                                                        <xs:annotation>
                                                            <xs:documentation>DEPRECATED</xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="454_native">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                A combination of 454 primary analysis output files, including 
                                                                seq
                                                                qual
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="454_native_seq">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                454 base calls (for example  .seq or .fna).
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="454_native_qual">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                454 quality scores  (for example  .qual).
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Helicos_native">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                A kind of fastq format specific to the Helicos platform.
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                               Please see SRA File Formats Guide for definitions of these file formats, 
                                                               and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                               Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_seq">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_prb">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_int">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_qseq">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_fastq">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>                                                      
                                                    </xs:enumeration>
                                                    <xs:enumeration value="Illumina_native_scarf">
                                                        <xs:annotation>
                                                            <xs:documentation> 
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                        </xs:enumeration>
                                                    <xs:enumeration value="SOLiD_native">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                A combination of SOLiD  primary analysis output files, including:
                                                                csfasta
                                                                _QV.qual
                                                                _intensity.ScaledCY3.fasta
                                                                _intensity.ScaledCY5.fasta
                                                                _intensity.ScaledFTC.fasta
                                                                _intensity.ScaledTXR.fasta
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="SOLiD_native_csfasta">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                Colorspace calls (for example .csfasta)
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value="SOLiD_native_qual">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                Colorspace quality scores (for example .qual)
                                                                Please see SRA File Formats Guide for definitions of these file formats, 
                                                                and the SRA Submission Guidelines document for data series that are appropriate for your study.
                                                                Sequence and qualities are minimally required.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                </xs:restriction>
                                            </xs:simpleType>
                                        </xs:attribute>

                                        <xs:attribute name="quality_scoring_system" use="optional">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    How the input data are scored for quality.  
                                                </xs:documentation>
                                            </xs:annotation>
                                            <xs:simpleType>
                                                <xs:restriction base = "xs:string">                             
                                                    <xs:enumeration value = "phred">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                The quality score is expressed as a probability of error in log form:
                                                                -10 log(1/p) where p is the probability of error, with value range 0..63,
                                                                0 meaning no base call.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value = "log-odds">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                The quality score is expressed as the ratio of error to non-error in log form:
                                                                -10 log(p/(1-p)) where p is the probability of error, with value range -40..40.
                                                                The SRA will convert these into phred scale during loadtime.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                </xs:restriction>
                                            </xs:simpleType>
                                        </xs:attribute>
                                        <xs:attribute name="quality_encoding"  use="optional">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    Character used in representing the minimum quality value.  Helps specify how to decode text rendering of quality data.
                                                </xs:documentation>
                                            </xs:annotation>
                                            <xs:simpleType>
                                                <xs:restriction base = "xs:string">
                                                    <xs:enumeration value = "ascii">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                ASCII character based encoding.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value = "decimal">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                               Single decimal value per quality score.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration> 
                                                    <xs:enumeration value = "hexadecimal">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                Single hexadecimal value per quality score.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>                                           
                                                </xs:restriction>
                                            </xs:simpleType>
                                        </xs:attribute>
                                        <xs:attribute name="ascii_offset"  use="optional">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    Character used in representing the minimum quality value.  Helps specify how to decode text rendering of quality data.
                                                </xs:documentation>
                                            </xs:annotation>
                                            <xs:simpleType>
                                                <xs:restriction base = "xs:string">
                                                    <xs:enumeration value = "!">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                ASCII value 33.  Typically used for range 0..63.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                    <xs:enumeration value = "@">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                ASCII value 64.  Typically used for range 0..60.
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
 
                                                </xs:restriction>
                                            </xs:simpleType>
                                        </xs:attribute>
                                        <xs:attribute name="checksum_method"  use="optional">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    Checksum method used.
                                                </xs:documentation>
                                            </xs:annotation>
                                            <xs:simpleType>
                                                <xs:restriction base = "xs:string">
                                                    <xs:enumeration value = "MD5">
                                                        <xs:annotation>
                                                            <xs:documentation>
                                                                Checksum generated by the MD5 method (md5sum in unix). 
                                                            </xs:documentation>
                                                        </xs:annotation>
                                                    </xs:enumeration>
                                                </xs:restriction>
                                            </xs:simpleType>
                                        </xs:attribute>
                                        <xs:attribute name="checksum"  type="xs:string" use="optional">
                                            <xs:annotation>
                                                <xs:documentation>
                                                    Checksum of uncompressed file.
                                                </xs:documentation>
                                            </xs:annotation>
                                        </xs:attribute>
                                    </xs:complexType>
                                     </xs:element>
                                      </xs:sequence>                                      
                                    </xs:complexType>
                                </xs:element>
                            </xs:sequence>

                            <xs:attribute name="name" type="xs:string" use="optional">
                                <xs:annotation>
                                    <xs:documentation>The plate/slide/flowcell name for this data block.
                                                                         (454) plate name
                                                                         (Illumina) flowcell name
                                                                         (SOLiD) slide name
                                                                         (Helicos) flowcell
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="sector" type="xs:nonNegativeInteger" use="optional">
                                <xs:annotation>
                                    <xs:documentation>Higher level partition of run data to which this data block pertains.
                                                                         (454) not used
                                                                         (Illumina) Lane number
                                                                         (SOLiD) slide
                                                                         (Helicos) channel
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="region" type="xs:nonNegativeInteger" use="optional">
                                <xs:annotation>
                                    <xs:documentation>Lower level partition of run data to which this data block pertains,
                                        typically the field of view for the imaging camera.
                                        (454) 0 if whole plate is used, 1..16 for gasket partition.
                                        (Illumina) Tile number (1..200+), or use 0 if file contains all tiles. 
                                        (SOLiD) Panel number (1..4096), or use 0 if file contains all panels. 
                                        (Helicos)  Field
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="total_spots" type="xs:nonNegativeInteger" use="optional">
                                <xs:annotation>
                                    <xs:documentation>DEPRECATED. The number of spots in this data block.</xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="total_reads" type="xs:nonNegativeInteger" use="optional">
                                <xs:annotation>
                                    <xs:documentation>DEPRECATED. The number of read records in the data block.</xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="number_channels" type="xs:nonNegativeInteger" use="optional">
                                <xs:annotation>
                                    <xs:documentation>DEPRECATED. The number of channels for which instrument data are provided. 
                                        This may be different from the data the instrument actually produced during the run.
                                        0 - no instrumentation data provided, bases only are provided.
                                        1 - one channel of instrumentation data per spot.
                                        4 - four channels of instrumentation data per spot.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="format_code" type="xs:int" use="optional">
                                <xs:annotation>
                                    <xs:documentation>
                                        DEPRECATED.  The format code determines alternate interpretations of the data.
                                        1 - Use the default for each platform.
                                        (454) Equivalent to published SFF flowgram format code #1.  
                                        Bases are interpreted as is.  Peaks are delta encoded for indexing into the flow sequence.
                                        The peak code is used to lookup the signal intensity (amplitude) for that flow to verify the decision
                                        of the base caller in calling base runs.  The intensity measurement is expressed as a two
                                        decimal place fixed point  integer value with maximum measured value 655.35.
                                        Qualities are estimated phred quality scores for each
                                        base call.  The spot name is set to the read name of the input data.  The cartesian coordinates
                                        of the spot are derived from the read name, which has base 36 encoding of x and y in the latter
                                        characters.  Flow sequence is factory default (TACG)*42 and the key sequence is factory default
                                        TCAG, unless these values are overridden in the experiment descriptor.
                                        (Illumina) 
                                        Bases are presented masked by binary peak values.   Qualities are delivered in 4-channel integer values which
                                        indicate probability of substitution.  Intensities are delivered in 4-channel decimal values scaled
                                        such that the minimum base-called channel intensity is still non-negative.  Spot coordinates pertain
                                        to x-y sites on the tile.  
                                        (SOLiD)
                                        Bases are presented in color space with the last base of the adapter sequence provided, which
                                        long with the appropriate color space conversion matrix can be used to decode the sequence.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="serial" type="xs:int" use="optional">
                                <xs:annotation>
                                    <xs:documentation>
                                        Integer value used to order loading of DATA_BLOCKs.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                            <xs:attribute name="member_name" type="xs:string" use="optional">
                                <xs:annotation>
                                    <xs:documentation>
                                        Allow for an individual DATA_BLOCK to be associated with a member of a sample pool.
                                    </xs:documentation>
                                </xs:annotation>
                            </xs:attribute>
                        </xs:complexType>
                    </xs:element>   
                </xs:sequence>
                <xs:element name="RUN_LINKS" minOccurs="0" maxOccurs="1">
                    <xs:annotation>
                        <xs:documentation>
                            Links to resources related to this RUN or RUN set (publication, datasets, online databases).
                        </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                        <xs:sequence minOccurs="1" maxOccurs="unbounded">
                            <xs:element name="RUN_LINK" type="LinkType"/>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element>
                
                <xs:element name="RUN_ATTRIBUTES" minOccurs="0" maxOccurs="1">
                    <xs:annotation>
                        <xs:documentation>
                            Properties and attributes of a RUN.  These can be entered as free-form 
                            tag-value pairs. For certain studies, submitters may be asked to follow a
                            community established ontology when describing the work.
                        </xs:documentation>
                    </xs:annotation>
                    <xs:complexType>
                        <xs:sequence maxOccurs="unbounded" minOccurs="1">
                            <xs:element name="RUN_ATTRIBUTE" type="AttributeType"/>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element> 
            </xs:sequence>
        <xs:attributeGroup ref="NameGroup"/>                       
            <xs:attribute name="instrument_model" use="optional">
                    <xs:annotation>
                        <xs:documentation>
                            DEPRECATED (use EXPERIMENT.PLATFORM..Instrument_model).  Instrument model actually used in the run.  Normally this is inherited from
                            the Experiment but may be overidden with the actual model.
                        </xs:documentation>
                    </xs:annotation>
                    <xs:simpleType>
                        <xs:restriction base="xs:string">
                            <xs:enumeration value="454 GS"/>
                            <xs:enumeration value="454 GS 20"/>
                            <xs:enumeration value="454 GS FLX"/>
                            <xs:enumeration value="Solexa 1G Genome Analyzer"><xs:annotation><xs:documentation>DEPRECATED</xs:documentation></xs:annotation></xs:enumeration>
                            <xs:enumeration value="Illumina Genome Analyzer"/>
                            <xs:enumeration value="Illumina Genome Analyzer II"/>
                            <xs:enumeration value="AB SOLiD System"/>
                            <xs:enumeration value="AB SOLiD System 2.0"/>
                            <xs:enumeration value="Other"/>
                        </xs:restriction>
                    </xs:simpleType>
            </xs:attribute>
            <xs:attribute name="instrument_name" type="xs:string" use="optional">
                <xs:annotation >
                    <xs:documentation >
                        Center-assigned name or id of the instrument used in the run.
                    </xs:documentation>
                </xs:annotation>
            </xs:attribute>
            <xs:attribute name="run_date" use="optional" type="xs:dateTime">
                <xs:annotation >
                    <xs:documentation >
                        ISO date when the run took place.  
                    </xs:documentation>
                </xs:annotation>
            </xs:attribute> 
            <xs:attribute name="run_file" use="optional" type="xs:string">
            <xs:annotation >
                <xs:documentation >
                    DEPRECATED.  Name of the submission file containing the run.  This file may have
                    included multiple runs.
                </xs:documentation>
            </xs:annotation>
            </xs:attribute>            
            <xs:attribute name="run_center" use="optional" type="xs:string">
                <xs:annotation >
                    <xs:documentation >
                        If applicable, the name of the contract sequencing center that executed the run.
                        Example: 454MSC.
                    </xs:documentation>
                </xs:annotation>
            </xs:attribute>
            <xs:attribute name="total_data_blocks" type="xs:positiveInteger" use="optional">
                <xs:annotation >
                    <xs:documentation >
                        DEPRECATED.  The number of data blocks to expect in this run.
                    </xs:documentation>
                </xs:annotation>
            </xs:attribute>
        </xs:complexType>

    <xs:element name="RUN_SET">
        <xs:annotation>
            <xs:documentation>
                RUN_SET serves as a container for a set of runs and a name space
                for establishing referential integrity between them. 
            </xs:documentation>
        </xs:annotation>  
        <xs:complexType >
            <xs:sequence minOccurs="1" maxOccurs="unbounded">
                <xs:element name="RUN" type="RunType"/>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
    
    <xs:element name="RUN" type="RunType"/>
</xs:schema>
'''

if __name__ == '__main__':
    main()

