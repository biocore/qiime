#!/usr/bin/env python
import copy
import datetime
import os
from string import strip
from collections import defaultdict
from cStringIO import StringIO
from hashlib import md5
from os.path import splitext, join
from sys import stderr
import xml.etree.ElementTree as ET

"""This module makes the submission xml files for SRA (study, experiment, etc.).

Assumes simple tab-delimited text input (allowing examples/comments; produces
xml output.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Kyle Bittinger", "Rohini Sinha"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"


def detect_missing_experiment_fields(input_file):
    """Return a list of required fields missing from an experiment input file."""
    return SraExperimentTable.parse(input_file).detect_missing_fields()


def detect_missing_study_fields(input_file):
    """Return a list of required fields missing from a study input file."""
    return SraStudyTable.parse(input_file).detect_missing_fields()


def detect_missing_submission_fields(input_file):
    """Return a list of required fields missing from a submission input file."""
    return SraSubmissionTable.parse(input_file).detect_missing_fields()


def detect_missing_sample_fields(input_file):
    """Return a list of required fields missing from a sample input file."""
    return SraSampleTable.parse(input_file).detect_missing_fields()


def md5_path(filename, block_size=8192):
    """Returns md5 hash from fileame without reading whole thing into memory"""
    m = md5('')
    infile = open(filename)
    curr = True
    while curr:
        curr = infile.read(block_size)
        m.update(curr)
    infile.close()
    return m.hexdigest()


def safe_for_xml(s):
    """Makes string s safe for xml by replacing entities."""
    return s.replace('&', '&amp;').replace('"','&quot;').replace("'",'&apos;').replace('<','&lt;').replace('>','&gt;')


def parse_tsv_with_header(tsv_file, data_fcn=None, header_fcn=None):
    """Parser for TSV files with commented headers

    File format:
      header line with field names
      optionally other comment lines starting with #
      tab-delimited fields

    Accepts keyword arguments data_fcn and header_fcn, which specify
    string-processing functions for each token found in the body and
    header of the tsv file.

    Returns a triple of the data, header, and comment lines
    """
    if data_fcn is None:
        data_fcn = lambda x: x
    if header_fcn is None:
        header_fcn = lambda x: x

    data = []
    header = []
    comments = []

    for line in tsv_file:
        # Skip lines containing only whitespace
        if not line.strip():
            continue
        line = line.rstrip('\n\r')
        if not header:
            if line.startswith('#'):
                line = line[1:]
            header = map(header_fcn, line.split('\t'))
            continue
        if line.startswith('#'):
            line = line[1:]
            comments.append(line)
        else:
            items = map(data_fcn, line.split('\t'))
            data.append(items)

    return data, header, comments


def read_tabular_data(tabular_file):
    """Reads tabular data from lines, skipping blanks"""
    def f(s):
        return safe_for_xml(s.strip('"'))
    def g(s):
        return s.strip()
    body, header, comments = parse_tsv_with_header(
        tabular_file, data_fcn=f, header_fcn=g)
    return header, body


class SraDerivedFieldError(Exception): pass

class SraInputFormatError(Exception): pass

class SraInputTable(object):
    """Class representing an input table for SRA submissions.
    """
    required_fields = []
    deprecated_fields = []
    
    def __init__(self, header, rows):
        """Create a new table from a header and a sequence of rows.
        """
        self.header = header
        self.rows = rows

    def __eq__(self, other):
        return (self.header == other.header) and (self.rows == other.rows)

    @property
    def entries(self):
        """Iterate over each row in the table as a dict.
        """
        for row in self.rows:
            yield self._get_entry(row)

    @property
    def first_entry(self):
        """Retrieve the first entry in the table.
        """
        return self._get_entry(self.rows[0])

    def _get_entry(self, row):
        """Transform a row of the table into a dict.
        """
        return dict([(k, v) for k, v in zip(self.header, row) if v])

    @classmethod
    def parse(cls, input_file):
        """Parse an input file to create a new Input Table.
        """
        header, rows = read_tabular_data(input_file)
        header = map(cls.canonicalize_field_name, header)
        return cls(header, rows)

    @classmethod
    def parse_twocol_format(cls, input_file):
        """Parse a legacy-format file to create a new Input Table.
        """
        _, records = read_tabular_data(input_file)
        table = cls([], [[]])
        for rec in records:
            field_name = cls.canonicalize_field_name(rec[0])
            table._ensure_field_exists(field_name)
            try:
                val = rec[1].strip()
            except IndexError:
                raise SraInputFormatError(
                    'Only one column found in two-column input file: %s' % rec)
            table.append_with_value(field_name, val)
        return table

    @classmethod
    def canonicalize_field_name(cls, field_name):
        return field_name.upper()

    def to_tsv(self):
        """Format an input table as a string of tab separated values.
        """
        lines = []
        lines.append('#' + '\t'.join(self.header))
        for row in self.rows:
            lines.append('\t'.join(row))
        return '\n'.join(lines)

    def derive_with_function(self, field, fcn):
        """Update a column of the table using a function.

        The funcion should take a single argument, the current entry
        in the table.
        """
        self._ensure_field_exists(field)

        field_idx = self.header.index(field)
        for row_num, row in enumerate(self.rows):
            if not row[field_idx]:
                entry = self._get_entry(row)
                try:
                    val = fcn(entry)
                except:
                    raise SraDerivedFieldError(
                        'Derivation of field %s failed in row %s, '
                        'parsed as %s.' % (field, row_num, entry))
                if val:
                    row[field_idx] = val

    def append_with_value(self, field, val, sep=','):
        self._ensure_field_exists(field)

        field_idx = self.header.index(field)
        for row in self.rows:
            if row[field_idx]:
                row[field_idx] = row[field_idx] + sep + val
            else:
                row[field_idx] = val

    def derive_with_value(self, field, val):
        self.derive_with_function(field, lambda x: val)

    def derive_with_format(self, field, format_string):
        """Update a column of the table using a format string.

        For each row, the format string will be evaluated using the
        current entry in the table.
        """
        try:
            self.derive_with_function(field, lambda x: format_string % x)
        except SraDerivedFieldError as e:
            raise SraDerivedFieldError(
                '%s Format string: %s' % (e, format_string))

    def _ensure_field_exists(self, fieldname):
        """Add a new field to the table.

        The value in each row will be set to an empty string.
        """
        if fieldname not in self.header:
            self.header.append(fieldname)
            for row in self.rows:
                row.append('')

    def derive_optional_fields(self):
        """Derive default values for optional fields.
        """
        pass

    def fix_deprecated_fields(self):
        """Move values from deprecated fields to new locations.
        """
        for old_field, new_field in self.deprecated_fields:
            if old_field in self.header:
                self._ensure_field_exists(new_field)
                stderr.write(
                    'Warning: The %s field has been deprecated. Please rename '
                    'the field to %s.\n' % (old_field, new_field))
                self.derive_with_format(new_field, '%(' + old_field + ')s')

    def detect_missing_fields(self):
        missing_fields = []
        for field in self.required_fields:
            if field not in self.header:
                missing_fields.append(field)
        return missing_fields


class SraStudyTable(SraInputTable):
    required_fields = [
        'CENTER_NAME',
        'CENTER_PROJECT_NAME',
        'STUDY_ABSTRACT',
        'STUDY_ALIAS',
        'STUDY_DESCRIPTION',
        'STUDY_TITLE',
        'STUDY_TYPE',
        ]


def make_study(input_file, twocol_input_format=True):
    """Returns string for study xml."""
    if twocol_input_format:
        table = SraStudyTable.parse_twocol_format(input_file)
    else:
        table = SraStudyTable.parse(input_file)
    study = SraStudySet.from_table(table)
    return study.to_xml_string()


class SraSubmissionTable(SraInputTable):
    """Class representing the input table for an SRA Submission
    """
    required_fields = [
        'CENTER_NAME',
        'CONTACT',
        'LAB_NAME',
        'SUBMISSION_DATE',
        'SUBMISSION_ID',
        ]

    def derive_optional_fields(self, data_dir=None):
        """Derive default values for optional fields.
        """
        date_obj = datetime.datetime.now().replace(microsecond=0)
        self.derive_with_format(
            'SUBMISSION_DATE', date_obj.isoformat())

        def derive_checksum(entry):
            filename = entry.get('FILE')
            if filename:
                if data_dir:
                    filepath = os.path.join(data_dir, filename)
                else:
                    filepath = filename
                return md5_path(filepath)
            return None
        
        self.derive_with_function(
            'SUBMISSION_FILE_CHECKSUM', derive_checksum)


def make_submission(submission_file, docnames=None, submission_dir=None,
                    twocol_input_format=True):
    """Returns string for submission xml.

    The docnames keyword specifies a dictionary of document filenames,
    which will be added to the submission.  The submission_dir keyword
    specifies the directory where these documents are to be found.  If
    the submission_dir is None, the current working directory is used.

    Finally, legacy two-column formatted files amy be used as input.
    The twocol_input_format keyword specifies if they are to be parsed
    using this format.
    """
    if twocol_input_format:
        table = SraSubmissionTable.parse_twocol_format(submission_file)
    else:
        table = SraSubmissionTable.parse(submission_file)
    table.derive_optional_fields(submission_dir)
    submission = SraSubmission.from_table(table)
    if docnames:
        submission.register_documents(docnames)
    return submission.to_xml_string()


class SraSampleTable(SraInputTable):
    """Class representing the input table for an SRA Sample
    """
    required_fields = [
        'SAMPLE_ALIAS',
        'TAXON_ID',
        'TITLE',
        ]

    @classmethod
    def canonicalize_field_name(cls, field_name):
        cn = super(SraSampleTable, cls).canonicalize_field_name(field_name)
        if cn in SraSampleAttributes.standard_fields():
            return cn
        return field_name


def make_sample(input_file):
    """Returns string for sample xml."""
    table = SraSampleTable.parse(input_file)
    sample_set = SraSampleSet.from_table(table)
    return sample_set.to_xml_string()


class SraExperimentTable(SraInputTable):
    required_fields = [
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
    deprecated_fields = [
        ('SAMPLE_ACCESSION', 'DEFAULT_SAMPLE_ACCESSION'),
        ]
    
    def derive_optional_fields(self):
        """Derive default values for optional fields.
        The optional field CHECKSUM is not handled by this function,
        but is left for the SraRunSet class to derive.
        """        
        f = self.derive_with_format
        g = self.derive_with_function

        f('REGION', '0')
        f('EXPERIMENT_ALIAS', '%(STUDY_REF)s_%(RUN_PREFIX)s')
        f('RUN_ALIAS', '%(STUDY_REF)s_%(SAMPLE_ALIAS)s_%(RUN_PREFIX)s')
        f('BARCODE_READ_GROUP_TAG', '%(RUN_PREFIX)s_%(BARCODE)s')
        g('PRIMER_READ_GROUP_TAG', self.__derive_primer_tag)
        g('POOL_MEMBER_NAME', self.__derive_pool_member_name)
        f('POOL_MEMBER_FILENAME', '%(POOL_MEMBER_NAME)s.sff')
        f('RUN_CENTER', '%(EXPERIMENT_CENTER)s')
        f('STUDY_CENTER', '%(EXPERIMENT_CENTER)s')
        f('SAMPLE_CENTER', '%(EXPERIMENT_CENTER)s')
        f('DEFAULT_SAMPLE_CENTER', '%(SAMPLE_CENTER)s')
        g('DEFAULT_SAMPLE_NAME', self.__derive_default_sample_name)
        f('DEFAULT_SAMPLE_FILENAME', '%(STUDY_REF)s_default_%(RUN_PREFIX)s.sff')
        f('DEFAULT_RUN_ALIAS', '%(STUDY_REF)s_default_%(RUN_PREFIX)s')
        f('PLATFORM', 'Titanium')
        f('KEY_SEQ', 'TCAG')
        f('LIBRARY_STRATEGY', 'AMPLICON')
        f('LIBRARY_SOURCE', 'GENOMIC')
        f('LIBRARY_SELECTION', 'PCR')

    def __derive_primer_tag(self, entry):
        primer = entry.get('PRIMER')
        if primer:
            try:
                return self.PRIMER_READ_GROUP_TAGS[primer]
            except KeyError:
                raise KeyError(
                    'No PRIMER_READ_GROUP_TAG has been provided, and '
                    'primer %s is not found in table of standard '
                    'primers. Please provide a value for the '
                    'PRIMER_READ_GROUP_TAG.' % primer)

    PRIMER_READ_GROUP_TAGS = {
        # Standard primers
        # BSR357
        'CTGCTGCCTYCCGTA': 'V1-V2',
        # BSR534
        'ATTACCGCGGCTGCTGGC': 'V1-V3',
        # BSR926
        'CCGTCAATTYYTTTRAGTTT': 'V3-V5',
        # BSR1492
        'GYTACCTTGTTAYGACTT': 'V6-V9',
        # Broad Institute Primers
        # 534r
        'ATTACCGCGGCTGCTGG': 'V1-V3b',
        # 926r
        'CCGTCAATTCMTTTRAGT': 'V3-V5b',
        # 1492r
        'TACGGYTACCTTGTTAYGACTT': 'V6-V9b',
        }

    def __derive_pool_member_name(self, entry):
        primer_tag = entry.get('PRIMER_READ_GROUP_TAG')
        if primer_tag:
            s = '%(RUN_PREFIX)s_%(SAMPLE_ALIAS)s_%(PRIMER_READ_GROUP_TAG)s'
        else:
            s = '%(RUN_PREFIX)s_%(SAMPLE_ALIAS)s'
        return s % entry

    def __derive_default_sample_name(self, entry):
        if not entry.get('DEFAULT_SAMPLE_ACCESSION'):
            return '%(STUDY_REF)s_default' % entry


def make_run(input_file, sff_dir=None):
    table = SraExperimentTable.parse(input_file)
    table.derive_optional_fields()
    table.fix_deprecated_fields()
    run_set = SraRunSet(sff_dir)
    for entry in table.entries:
        run_set.register(entry)
    return run_set.to_xml_string()


def make_experiment(input_file, attribute_file=None, link_file=None):
    table = SraExperimentTable.parse(input_file)
    table.derive_optional_fields()
    table.fix_deprecated_fields()
    experiment_set = SraExperimentSet()

    for entry in table.entries:
        experiment_set.register(entry)

    experiment_attributes = defaultdict(list)
    if attribute_file is not None:
        _, attr_body = read_tabular_data(attribute_file)
        experiment_attributes = threecol_data_to_dict(attr_body)

    experiment_links = defaultdict(list)
    if link_file is not None:
        _, link_body = read_tabular_data(link_file)
        experiment_links = threecol_data_to_dict(link_body)

    for alias, expt in experiment_set.experiments:
        expt.experiment_attributes = experiment_attributes.get(alias)
        expt.experiment_links = experiment_links.get(alias)
        
    return experiment_set.to_xml_string()


class SraEntity(object):
    """Base class for elements defined in the SRA's XML schema.

    The SraEntity class forms the foundation for generating valid XML
    from input files.  The overall strategy is to create a simplified
    representation of the XML document with a set of SraEntity
    classes.  If this is implemented properly, users can (1)
    instantiate an SraEntity object representing the top-level XML
    document (2) register each line of the input file with the
    top-level object's register() method, and (3) generate an
    ElementTree XML representation of the document from the top-level
    object's to_xml() method.  Each step is described in detail below.

    (1) An instance of the SraEntity class corresponding to the
    top-level document is created.  Ususally, this constructor has no
    arguments.

    (2) Each row of an input file is parsed as a dictionary and passed
    to the top-level SraEntity object, using the register() method.
    Inside the register() method, the field values of the new entry
    represented at this level of the document are checked for
    consistency with previous entries.  Where a new SraEntity
    sub-element is needed in the document, it is created from the new
    entry using the from_entry() factory method.  Where an SraEntity
    sub-element already exists, the new entry is passed to it via the
    register() method.

    The process of generating a valid document continues as each new
    entry is registered with the top-level SraEntity.  New elements of
    the document are generated as the new entry propogates down the
    document tree, and each entry is checked for consistency with
    existing field values.  Once all entries are registered with the
    top-level SraEntity, this object should contain a valid
    representation of the XML document as a tree of SraEntity objects.

    (3) The XML output is generated via the to_xml() method of the
    top-level entity.  This method calls the to_xml() method on all
    sub-elements, and returns an ElementTree representation of the XML
    document.
    """

    """Maps object attributes to required fields in a table entry."""
    required_fields = {}
    """Maps object attributes to optional fields in a table entry."""
    optional_fields = {}

    def __init__(self, attrs):
        """Create a new SRA entity from a dict of attributes.

        Sets instance attributes to the values provided in attrs.

        Using a class-level list of attribute names, it ensures that
        all necessary instance attributes are set, even if they are
        not provided in attrs (None is the default value).  Provided
        attributes must be found in the class-level list of attribute
        names, or a ValueError is raised.

        This method does NOT check to make sure that all required
        attributes are present in the attrs dictionary.  This
        responsibility falls on the register() and from_entry()
        methods.
        """
        self._check_attrs(attrs)
        for attribute_name in self.attribute_names:
            val = attrs.get(attribute_name)
            setattr(self, attribute_name, val)

    @property
    def attribute_names(self):
        return self.required_fields.keys() + self.optional_fields.keys()

    def _check_attrs(self, attrs):
        """Ensure that provided attribute keys are valid.

        To be valid, a provided key must be found in the class-level
        list of attribute names.
        """
        registered_attrs = set(self.attribute_names)
        provided_attrs = set(attrs.keys())
        unrecognized_attrs = provided_attrs.difference(registered_attrs)
        if unrecognized_attrs:
            raise ValueError(
                'Unrecognized keys provided in attribute dictionary: %s' % \
                ', '.join(unrecognized_attrs))

    def register(self, entry):
        """Register a new table entry with an existing SRA element.

        In the default case, this method checks that the relevant
        attributes of the new entry are consistent with the existing
        attributes.  This method raises a ValueError if an
        inconsistency is found.
        """
        for attr_name, field_name in self.required_fields.items():
            observed_value = entry[field_name]
            expected_value = getattr(self, attr_name)
            if observed_value != expected_value:
                raise ValueError(
                    'Observed value "%s" for required field %s is not '
                    'consistent with existing value "%s".' % (
                        observed_value, field_name, expected_value))
        for attr_name, field_name in self.optional_fields.items():
            observed_value = entry.get(field_name)
            expected_value = getattr(self, attr_name)
            if observed_value != expected_value:
                raise ValueError(
                    'Observed value "%s" for optional field %s is not '
                    'consistent with existing value "%s".' % (
                        observed_value, field_name, expected_value))

    @classmethod
    def from_entry(cls, entry):
        """Factory method to create a new SraEntity from a table entry."""
        attrs = cls.gather_instance_attributes(entry)
        return cls(attrs)

    def to_xml(self):
        """Generate an ElementTree XML representation of the object."""
        raise NotImplementedError(
            'This method is not implemented in the base class.')

    def to_xml_string(self, level=0):
        return pretty_xml(self.to_xml(), level=level, encoding='UTF-8')

    @classmethod
    def gather_instance_attributes(cls, entry):
        """Gather object attribute values from a table entry.

        The mappings between object attributes and field names of the
        table entry are specified in two class-level dictionaries:
        required_fields and optional_fields.

        Raises a KeyError if a required field is not present.
        """
        attrs = {}
        for attr_name, field_name in cls.required_fields.items():
            try:
                attrs[attr_name] = entry[field_name]
            except KeyError:
                raise KeyError(
                    'Required field %s was not found for the following line '
                    'item in the input file: %s' % (field_name, entry))
        for attr_name, field_name in cls.optional_fields.items():
            attrs[attr_name] = entry.get(field_name)
        return attrs
    
    def set_xml_attributes(self, xml_element, instance_attributes):
        """Transfer attributes from an instance to an XML element.

        Transfers instance attributes from the list that evaluate as
        non-false.  Blank attributes are not transferred.

        This method mutates the XML element as a side effect.  To make
        this clear, we do not provide a return value.
        """
        for attr in instance_attributes:
            val = getattr(self, attr)
            if val:
                xml_element.set(attr, val)


class SraSubmission(SraEntity):
    """Class representing an SRA Submission entity."""

    required_fields = {
        'alias': 'SUBMISSION_ID',
        'center': 'CENTER_NAME',
        'lab': 'LAB_NAME',
        'date': 'SUBMISSION_DATE',
        }

    optional_fields = {
        'comment': 'SUBMISSION_COMMENT',
        'accession': 'ACCESSION',
        }
    
    def __init__(self, attrs, contacts, files=None):
        """Create a new SRA Submission Entity."""
        super(SraSubmission, self).__init__(attrs)
        self.contacts = contacts
        self.files = files
        self.actions = None

    @classmethod
    def from_table(cls, table):
        """Create a new SRA Submission from an input table.
        """
        entry = table.first_entry
        attrs = cls.gather_instance_attributes(entry)
        contacts = SraContacts.from_entry(entry)

        if entry.get('SUBMISSION_FILE_CHECKSUM'):
            files = SraSubmissionFiles.from_entry(entry)
        else:
            files = None
        return cls(attrs, contacts, files)

    def register_documents(self, document_dict):
        """Register a set of documents with the submission."""
        self.actions = SraActions.from_entry(document_dict)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SUBMISSION entity."""
        root = ET.Element('SUBMISSION')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        if self.accession:
            root.set('accession', self.accession)
        root.set('submission_id', self.alias)
        root.set('center_name', self.center)
        if self.comment:
            root.set('submission_comment', self.comment)
        root.set('lab_name', self.lab)
        root.set('submission_date', self.date)
        root.append(self.contacts.to_xml())
        if self.actions:
            root.append(self.actions.to_xml())
        if self.files:
            root.append(self.files.to_xml())
        return root


class SraSubmissionFiles(SraEntity):
    """Class representing an SRA Files entity."""
    
    required_fields = {
        'filename': 'FILE',
        'checksum': 'SUBMISSION_FILE_CHECKSUM',
        }

    def to_xml(self):
        """Create an ElementTree XML object for the SRA FILES entity."""
        root = ET.Element('FILES')
        ET.SubElement(root, 'FILE', filename=self.filename, checksum_method='MD5',
                      checksum=self.checksum)
        return root


class SraActions(SraEntity):
    """Class representing an SRA Actions entity."""
    
    schema = ['study', 'sample', 'experiment', 'run']
    optional_fields = dict([(x, x) for x in schema])

    def to_xml(self):
        """Create an ElementTree XML object for the SRA ACTIONS entity."""
        root = ET.Element('ACTIONS')
        for sch in self.schema:
            source = getattr(self, sch)
            if source:
                action_elem = ET.SubElement(root, 'ACTION')
                add_elem = ET.SubElement(
                    action_elem, 'ADD', schema=sch, source=source)
                add_elem.set('notes', '%s metadata' % sch)
        action_elem = ET.SubElement(root, 'ACTION')
        ET.SubElement(action_elem, 'RELEASE')
        return root


class SraContacts(SraEntity):
    """Class representing an SRA Contacts Entity."""
    
    required_fields = {
        'contact_string': 'CONTACT'
        }
    
    def __init__(self, attrs):
        """Create a new set of SRA Contacts."""
        super(SraContacts, self).__init__(attrs)
        self.contacts = self._parse_contacts(self.contact_string)

    @staticmethod
    def _parse_contacts(contact_string):
        """Parse a list of contacts (name, email pairs) from a string."""
        contacts = []
        for c in contact_string.split(','):
            toks = map(strip, c.split(';'))
            if len(toks) != 2:
                raise ValueError(
                    'Improperly formatted contacts string: "%s" (name and '
                    'email should be separated by a semicolon in substring '
                    '"%s")' % (self.contact_string, c))
            name = toks[0]
            email = toks[1]
            contacts.append((name, email))
        return contacts

    def to_xml(self):
        """Create an ElementTree XML object for the SRA CONTACTS entity."""
        root = ET.Element('CONTACTS')
        for name, email in self.contacts:
            ET.SubElement(
                root, 'CONTACT', name=name, inform_on_status=email,
                inform_on_error=email)
        return root


class SraStudySet(SraEntity):
    """Class repersenting a set of SRA Studies."""

    def __init__(self):
        """Create a new set of SRA Studies"""
        self.studies = []

    @classmethod
    def from_table(cls, table):
        instance = cls()
        for entry in table.entries:
            instance.register(entry)
        return instance

    def register(self, entry):
        """Register a new table entry with the study set."""
        new_study = SraStudy.from_entry(entry)
        self.studies.append(new_study)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA STUDY_SET entity."""
        root = ET.Element('STUDY_SET')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        for s in self.studies:
            root.append(s.to_xml())
        if not self.studies:
            root.text = '\n'
        return root


class SraStudy(SraEntity):
    """Class representing an SRA Study entity."""

    required_fields = {
        'alias': 'STUDY_ALIAS',
        'title': 'STUDY_TITLE',
        'type': 'STUDY_TYPE',
        'abstract': 'STUDY_ABSTRACT',
        'description': 'STUDY_DESCRIPTION',
        'center': 'CENTER_NAME',
        'project': 'CENTER_PROJECT_NAME',
        }
    optional_fields = {
        'pubmed_id': 'PMID',
        }

    def to_xml(self):
        """Create an ElementTree XML object for the SRA STUDY entity."""
        root = ET.Element('STUDY')
        self.set_xml_attributes(root, ['alias'])
        desc = ET.SubElement(root, 'DESCRIPTOR')
        ET.SubElement(desc, 'STUDY_TITLE').text = self.title
        ET.SubElement(desc, 'STUDY_TYPE', existing_study_type=self.type)
        ET.SubElement(desc, 'STUDY_ABSTRACT').text = self.abstract
        ET.SubElement(desc, 'STUDY_DESCRIPTION').text = self.description
        ET.SubElement(desc, 'CENTER_NAME').text = self.center
        ET.SubElement(desc, 'CENTER_PROJECT_NAME').text = self.project
        if self.pubmed_id:
            link_elem = ET.SubElement(ET.SubElement(ET.SubElement(
                root, 'STUDY_LINKS'), 'STUDY_LINK'), 'ENTREZ_LINK')
            ET.SubElement(link_elem, 'DB').text = 'pubmed'
            ET.SubElement(link_elem, 'ID').text = self.pubmed_id
        return root

class SraSampleSet(SraEntity):
    """Class repersenting a set of SRA samples."""

    def __init__(self):
        """Create a new set of SRA Samples"""
        self.samples = []

    @classmethod
    def from_table(cls, table):
        instance = cls()
        for entry in table.entries:
            instance.register(entry)
        return instance

    def register(self, entry):
        new_sample = SraSample.from_entry(entry)
        self.samples.append(new_sample)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SAMPLE_SET entity."""
        root = ET.Element('SAMPLE_SET')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        for s in self.samples:
            root.append(s.to_xml())
        if not self.samples:
            root.text = '\n'
        return root


class SraSample(SraEntity):
    """Class representing an SRA Sample entity."""
    
    required_fields = {
        'alias': 'SAMPLE_ALIAS',
        'title': 'TITLE',
        }
    optional_fields = {
        'description': 'DESCRIPTION',
        }

    def __init__(self, attrs, sample_name, sample_attributes):
        """Create a new SRA sample."""
        super(SraSample, self).__init__(attrs)
        self.sample_name = sample_name
        self.sample_attributes = sample_attributes

    @classmethod
    def from_entry(cls, entry):
        """Factory method to creat a new sample from a table entry."""
        attrs = cls.gather_instance_attributes(entry)
        sample_name = SraSampleName.from_entry(entry)
        sample_attributes = SraSampleAttributes.from_entry(entry)
        return cls(attrs, sample_name, sample_attributes)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SAMPLE entity."""
        root = ET.Element('SAMPLE')
        self.set_xml_attributes(root, ['alias'])
        ET.SubElement(root, 'TITLE').text = self.title
        if self.sample_name:
            root.append(self.sample_name.to_xml())
        ET.SubElement(root, 'DESCRIPTION').text = self.description
        if self.sample_attributes:
            root.append(self.sample_attributes.to_xml())
        return root


class SraSampleName(SraEntity):
    """Class representing an SRA Sample name."""

    optional_fields = {
        'taxon_id': 'TAXON_ID',
        'common_name': 'COMMON_NAME',
        'anonymized_name': 'ANONYMIZED_NAME',
        }

    def __nonzero__(self):
        return any([getattr(self, k) for k in self.optional_fields.keys()])

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SAMPLE_NAME entity."""
        root = ET.Element('SAMPLE_NAME')
        subelements = [(y, x) for (x, y) in self.optional_fields.items()]
        if not subelements:
            return None
        # sort by element name
        subelements.sort()
        for elem_name, attr in subelements:
            val = getattr(self, attr)
            if val:
                x = ET.SubElement(root, elem_name)
                x.text = val
        return root


class SraSampleAttributes(SraEntity):
    """Class representing a set of SRA Sample Attribute entities."""
    
    def __init__(self, attrs):
        """Create a new SRA sample attribute set.

        Attribute tags may be any string, so we place no restrictions
        on them.
        """
        self.sample_attributes = attrs

    def __nonzero__(self):
        return bool(self.sample_attributes)

    @classmethod
    def from_entry(cls, entry):
        """Factory method to create new sample attributes from a table entry.
        """
        attrs = {}
        for k in cls.attribute_fields(entry):
            attrs[k] = entry[k]
        return cls(attrs)

    @classmethod
    def attribute_fields(cls, entry):
        """Determine the field names that correspond to SRA Attribute tags.

        This is computed by comparing the set of field names in the
        entry to the set of standard field names for an SRA Sample.
        If there is no match, the field is assumed to be an SRA Sample
        Atribute.
        """
        entry_keys = set(entry.keys())
        return entry_keys.difference(cls.standard_fields())

    @staticmethod
    def standard_fields():
        """Find the list of all standard field names for an SRA Sample.
        """
        fs = set()
        sample_classes = [SraSampleName, SraSample]
        for sample_class in sample_classes:
            fs.update(sample_class.required_fields.values())
            fs.update(sample_class.optional_fields.values())
        return fs

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SAMPLE_ATTRIBUTES entity."""
        if not self.sample_attributes:
            return None
        root = ET.Element('SAMPLE_ATTRIBUTES')
        attributes = self.sample_attributes.items()
        attributes.sort()
        for tag, value in attributes:
            attr_elem = ET.SubElement(root, 'SAMPLE_ATTRIBUTE')
            ET.SubElement(attr_elem, 'TAG').text = tag
            ET.SubElement(attr_elem, 'VALUE').text = value
        return root
        

class SraExperimentSet(SraEntity):
    """Class representing a set of SRA Experiment entities."""

    def __init__(self):
        """Create a new SRA experiment set.
        """
        self.experiments = []

    def register(self, entry):
        """Register a new entry with the experiment set.

        New entries are assigned to experiments based on the value of
        the EXPERIMENT_ALIAS field.
        """
        is_registered = False
        # Avoid raising a KeyError if EXPERIMENT_ALIAS is missing.
        # The error will be raised with a better message in
        # SraExperiment.from_entry()
        entry_alias = entry.get('EXPERIMENT_ALIAS')
        # Use list of tuples instead of dict to preserve ordering
        for alias, experiment in self.experiments:
            if alias == entry_alias:
                experiment.register(entry)
                is_registered = True
                continue
        if not is_registered:
            new_experiment = SraExperiment.from_entry(entry)
            new_experiment.register(entry)
            self.experiments.append((entry_alias, new_experiment))

    def to_xml(self):
        """Create an ElementTree XML object for the SRA EXPERIMENT_SET entity."""
        root = ET.Element('EXPERIMENT_SET')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        for alias, expt in self.experiments:
            root.append(expt.to_xml())
        if not self.experiments:
            root.text = '\n'
        return root

class SraExperiment(SraEntity):
    """Class representing an SRA Experiment entity."""
    
    required_fields = {
        'alias': 'EXPERIMENT_ALIAS',
        'center_name': 'EXPERIMENT_CENTER',
        'title': 'EXPERIMENT_TITLE',
        'design_description': 'EXPERIMENT_DESIGN_DESCRIPTION',
        'library_strategy': 'LIBRARY_STRATEGY',
        'library_source': 'LIBRARY_SOURCE',
        'library_selection': 'LIBRARY_SELECTION',
        'library_construction_protocol': 'LIBRARY_CONSTRUCTION_PROTOCOL',
        }

    def __init__(self, attrs, study_ref, platform, sample_descriptor,
                 spot_descriptor, experiment_links=None,
                 experiment_attributes=None):
        """Create a new SRA Experiment.

        A valid experiment requires the following child entities,
        which are provided as positional arguments: SraStudyRef,
        SraPlatform, SraSampleDescriptor, and SraSpotDescriptor.
        """
        super(SraExperiment, self).__init__(attrs)
        self.study_ref = study_ref
        self.platform = platform
        self.sample_descriptor = sample_descriptor
        self.spot_descriptor = spot_descriptor
        self.experiment_links = experiment_links
        self.experiment_attributes = experiment_attributes

    @classmethod
    def from_entry(cls, entry):
        """Factory method to create a new experiment from a table entry."""
        attrs = cls.gather_instance_attributes(entry)
        platform = SraPlatform.from_entry(entry)
        study_ref = SraStudyRef.from_entry(entry)
        sample_descriptor = SraSampleDescriptor.from_entry(entry)
        spot_descriptor = SraSpotDescriptor.from_entry(entry)
        return cls(attrs, study_ref, platform, sample_descriptor, spot_descriptor)

    def register(self, entry):
        """Register a new table entry with the experiment."""
        super(SraExperiment, self).register(entry)
        self.platform.register(entry)
        self.study_ref.register(entry)
        self.sample_descriptor.register(entry)
        self.spot_descriptor.register(entry)

    def __library_descriptor_xml(self):
        ld_elem = ET.Element('LIBRARY_DESCRIPTOR')
        ET.SubElement(ld_elem, 'LIBRARY_NAME').text = self.alias
        ET.SubElement(ld_elem, 'LIBRARY_STRATEGY').text = self.library_strategy
        ET.SubElement(ld_elem, 'LIBRARY_SOURCE').text = self.library_source
        ET.SubElement(ld_elem, 'LIBRARY_SELECTION').text = self.library_selection
        ll_elem = ET.SubElement(ld_elem, 'LIBRARY_LAYOUT')
        ET.SubElement(ll_elem, 'SINGLE')
        lcp_elem = ET.SubElement(ld_elem, 'LIBRARY_CONSTRUCTION_PROTOCOL')
        lcp_elem.text = self.library_construction_protocol
        return ld_elem

    @staticmethod
    def __processing_xml():
        proc_elem = ET.Element('PROCESSING')
        bc_elem = ET.SubElement(proc_elem, 'BASE_CALLS')
        ET.SubElement(bc_elem, 'SEQUENCE_SPACE').text = 'Base Space'
        ET.SubElement(bc_elem, 'BASE_CALLER').text = '454 BaseCaller'
        qs_elem = ET.SubElement(proc_elem, 'QUALITY_SCORES', qtype='phred')
        ET.SubElement(qs_elem, 'QUALITY_SCORER').text = '454 BaseCaller'
        ET.SubElement(qs_elem, 'NUMBER_OF_LEVELS').text = '40'
        ET.SubElement(qs_elem, 'MULTIPLIER').text = '1.0'
        return proc_elem

    def to_xml(self):
        """Create an ElementTree XML object for the SRA EXPERIMENT entity."""
        root = ET.Element('EXPERIMENT')
        self.set_xml_attributes(root, ['alias', 'center_name'])
        ET.SubElement(root, 'TITLE').text = self.title
        root.append(self.study_ref.to_xml())
        design_elem = ET.SubElement(root, 'DESIGN')
        dd_elem = ET.SubElement(design_elem, 'DESIGN_DESCRIPTION')
        dd_elem.text = self.design_description

        design_elem.append(self.sample_descriptor.to_xml())
        design_elem.append(self.__library_descriptor_xml())
        design_elem.append(self.spot_descriptor.to_xml())
        root.append(self.platform.to_xml())
        root.append(self.__processing_xml())
        
        if self.experiment_links:
            root.append(_experiment_link_xml(self.experiment_links))
        if self.experiment_attributes:
            root.append(_experiment_attribute_xml(self.experiment_attributes))
        
        return root


class SraStudyRef(SraEntity):
    """Class representing an SRA Study Ref entity."""

    optional_fields = {
        'refcenter': 'SAMPLE_CENTER',
        'refname': 'STUDY_REF',
        'accession': 'STUDY_ACCESSION',
        }

    def to_xml(self):
        """Create an ElementTree XML object for the SRA STUDY_REF entity."""
        root = ET.Element('STUDY_REF')
        self.set_xml_attributes(root, ['refname', 'refcenter', 'accession'])
        return root


class SraPlatform(SraEntity):
    """Class representing an SRA Platform entity."""
    
    required_fields = {
        'platform': 'PLATFORM',
        }
    standard_platforms = {
        'FLX': {
            'instrument_model': '454 GS FLX',
            'flow_sequence': 'TACG',
            'flow_count': '400',
            },
        'Titanium': {
            'instrument_model': '454 Titanium',
            'flow_sequence': 'TACG',
            'flow_count': '800',            
            },
        }

    @classmethod
    def from_entry(cls, entry):
        """Factory method to create a new platform from a table entry."""
        platform = entry['PLATFORM']
        if platform not in cls.standard_platforms:
            raise ValueError(
                'Platform %s not recognized (supported options are %s)' % \
                (platform, cls.standard_platforms.keys()))
        return super(SraPlatform, cls).from_entry(entry)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA PLATFORM entity."""
        platform_info = self.standard_platforms[self.platform]
        root = ET.Element('PLATFORM')
        ls454_elem = ET.SubElement(root, 'LS454')
        ET.SubElement(ls454_elem, 'INSTRUMENT_MODEL').text = \
            platform_info['instrument_model']
        ET.SubElement(ls454_elem, 'FLOW_SEQUENCE').text = \
            platform_info['flow_sequence']
        ET.SubElement(ls454_elem, 'FLOW_COUNT').text = \
            platform_info['flow_count']
        return root


class SraSpotDescriptor(SraEntity):
    """Class representing an SRA Spot Descriptor entity."""

    required_fields = {
        'adapter': 'KEY_SEQ',
        }
    optional_fields = {
        'linker': 'LINKER',
        }
    
    def __init__(self, attrs):
        """Create a new SRA Spot Descriptor."""
        super(SraSpotDescriptor, self).__init__(attrs)
        self.barcode_basecalls = []
        self.primer_basecalls = []
        self.barcodes = set()
        self.primers = set()

    def register(self, entry):
        """Register a new table entry with the spot descriptor."""
        super(SraSpotDescriptor, self).register(entry)
        
        barcode = entry.get('BARCODE')
        if barcode and barcode not in self.barcodes:
            self.barcode_basecalls.append(SraBarcodeBasecall.from_entry(entry))
            self.barcodes.add(barcode)

        primer = entry.get('PRIMER')
        if primer and primer not in self.primers:
            self.primer_basecalls.append(SraPrimerBasecall.from_entry(entry))
            self.primers.add(primer)

    def __adapter_spec_xml(self, read_index):
        root = ET.Element('READ_SPEC')
        ET.SubElement(root, 'READ_INDEX').text = str(read_index)
        ET.SubElement(root, 'READ_CLASS').text = 'Technical Read'
        ET.SubElement(root, 'READ_TYPE').text = 'Adapter'
        ET.SubElement(root, 'EXPECTED_BASECALL').text = self.adapter
        return root

    def __barcode_spec_xml(self, read_index):
        root = ET.Element('READ_SPEC')
        ET.SubElement(root, 'READ_INDEX').text = str(read_index)
        # According to SRA.experiment.xsd v1.1, the READ_LABEL element
        # must appear immediately following the READ_INDEX
        ET.SubElement(root, 'READ_LABEL').text = 'barcode'
        ET.SubElement(root, 'READ_CLASS').text = 'Technical Read'
        ET.SubElement(root, 'READ_TYPE').text = 'BarCode'
        table_elem = ET.SubElement(root, 'EXPECTED_BASECALL_TABLE')        
        for b in self.barcode_basecalls:
            table_elem.append(b.to_xml())
        return root

    def __linker_spec_xml(self, read_index):
        root = ET.Element('READ_SPEC')
        ET.SubElement(root, 'READ_INDEX').text = str(read_index)
        ET.SubElement(root, 'READ_LABEL').text = 'linker'
        ET.SubElement(root, 'READ_CLASS').text = 'Technical Read'
        ET.SubElement(root, 'READ_TYPE').text = 'Linker'
        ET.SubElement(root, 'EXPECTED_BASECALL').text = self.linker
        return root

    def __primer_spec_xml(self, read_index):
        root = ET.Element('READ_SPEC')
        ET.SubElement(root, 'READ_INDEX').text = str(read_index)
        ET.SubElement(root, 'READ_LABEL').text = 'rRNA_primer'
        ET.SubElement(root, 'READ_CLASS').text = 'Technical Read'
        ET.SubElement(root, 'READ_TYPE').text = 'Primer'
        table_elem = ET.SubElement(root, 'EXPECTED_BASECALL_TABLE')        
        for b in self.primer_basecalls:
            table_elem.append(b.to_xml())
        return root

    def __application_spec_xml(self, read_index):
        root = ET.Element('READ_SPEC')
        ET.SubElement(root, 'READ_INDEX').text = str(read_index)
        ET.SubElement(root, 'READ_CLASS').text = 'Application Read'
        ET.SubElement(root, 'READ_TYPE').text = 'Forward'
        prev_index = int(read_index) - 1
        order_elem = ET.SubElement(root, 'RELATIVE_ORDER')
        order_elem.set('follows_read_index', str(prev_index))
        return root

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SPOT_DESCRIPTOR entity."""
        root = ET.Element('SPOT_DESCRIPTOR')
        decode_elem = ET.SubElement(root, 'SPOT_DECODE_SPEC')
        i = 0
        decode_elem.append(self.__adapter_spec_xml(i))
        i += 1
        if self.barcode_basecalls:
            decode_elem.append(self.__barcode_spec_xml(i))
            i += 1
        if self.linker:
            decode_elem.append(self.__linker_spec_xml(i))
            i += 1
        if self.primers:
            decode_elem.append(self.__primer_spec_xml(i))
            i += 1
        decode_elem.append(self.__application_spec_xml(i))
        return root


class SraBarcodeBasecall(SraEntity):
    """Class representing an SRA Basecall entity for a barcode."""

    required_fields = {
        'basecall': 'BARCODE',
        'read_group_tag': 'BARCODE_READ_GROUP_TAG',
        }

    def to_xml(self):
        """Create an ElementTree XML object for the SRA BASECALL entity."""
        root = ET.Element('BASECALL')
        root.set('read_group_tag', self.read_group_tag)
        root.set('min_match', str(len(self.basecall)))
        root.set('max_mismatch', '0')
        root.set('match_edge', 'full')
        root.text = self.basecall
        return root


class SraPrimerBasecall(SraEntity):
    """Class representing an SRA Basecall entity for a primer."""

    required_fields = {
        'basecall': 'PRIMER',
        'read_group_tag': 'PRIMER_READ_GROUP_TAG',
        }

    def to_xml(self):
        """Create an ElementTree XML object for the SRA BASECALL entity."""
        root = ET.Element('BASECALL')
        root.set('read_group_tag', self.read_group_tag)
        root.set('min_match', str(len(self.basecall)))
        root.set('max_mismatch', '0')
        root.set('match_edge', 'full')
        root.text = self.basecall
        return root


class SraSampleDescriptor(SraEntity):
    """Class representing an SRA Sample Descriptor entity."""

    optional_fields = {
        'refcenter': 'DEFAULT_SAMPLE_CENTER',
        'refname': 'DEFAULT_SAMPLE_NAME',
        'accession': 'DEFAULT_SAMPLE_ACCESSION',
        }

    def __init__(self, attrs):
        """Create a new SRA dample descriptor."""
        super(SraSampleDescriptor, self).__init__(attrs)
        self.pool_members = []

    def register(self, entry):
        """Register a new table entry with the sample descriptor."""
        super(SraSampleDescriptor, self).register(entry)
        self.pool_members.append(SraPoolMember.from_entry(entry))

    def to_xml(self):
        """Create an ElementTree XML object for the SRA SAMPLE_DESCRIPTOR entity."""
        root = ET.Element('SAMPLE_DESCRIPTOR')
        self.set_xml_attributes(root, ['refname', 'refcenter', 'accession'])
        if self.pool_members:
            pool_elem = ET.SubElement(root, 'POOL')
            for m in self.pool_members:
                pool_elem.append(m.to_xml())
        return root


class SraPoolMember(SraEntity):
    """Class representing an SRA Pool Member entity."""

    required_fields = {
        'refname': 'SAMPLE_ALIAS',
        'refcenter': 'SAMPLE_CENTER',
        'member_name': 'POOL_MEMBER_NAME',
        'proportion': 'POOL_PROPORTION',
        }
    optional_fields = {
        'accession': 'POOL_MEMBER_ACCESSION',
        'barcode_tag': 'BARCODE_READ_GROUP_TAG',
        'primer_tag': 'PRIMER_READ_GROUP_TAG',
        }

    def __init__(self, attrs, barcode_tag=None, primer_tag=None):
        """Create a new SRA Pool Member.

        According to the SRA schema, all MEMBER attributes are
        optional, including member_name.  We require a minimum of
        refname, refcenter, member_name, and proportion.

        According to the SRA schema, a MEMBER element must contain at
        least one READ_LABEL.  Therefore, this method raises a
        ValueError if the barcode read group tag and the primer read
        group tag are both None.
        """
        super(SraPoolMember, self).__init__(attrs)
        if (self.barcode_tag is None) and (self.primer_tag is None):
            raise ValueError(
                'Must provide either barcode read group tag or primer read '
                'group tag to create a pool member.')

    def to_xml(self):
        """Create an ElementTree XML object for the SRA MEMBER entity."""
        root = ET.Element('MEMBER')
        self.set_xml_attributes(
            root, ['refname', 'refcenter', 'member_name', 'proportion', 'accession'])
        if self.barcode_tag:
            barcode_elem = ET.SubElement(
                root, 'READ_LABEL', read_group_tag=self.barcode_tag)
            barcode_elem.text = 'barcode'
        if self.primer_tag:
            primer_elem = ET.SubElement(
                root, 'READ_LABEL', read_group_tag=self.primer_tag)
            primer_elem.text = 'rRNA_primer'
        return root


class SraRunSet(SraEntity):
    """Class representing a set of SRA Run entities.

    SraRunSet manages the list of runs, as well as default files for
    each run.
    """

    def __init__(self, data_dir):
        """Create a new set of SRA Run entities."""
        self.data_dir = data_dir
        self.runs = []
        self.default_entries = set()

    def register(self, entry):
        """Register a new table entry with the SRA Run Set.

        Searches for each file in the data directory under a folder
        named RUN_PREFIX.  If the file is found, this method computes
        the checksum and returns the SraFile object.  If the file is
        not found, the method prints a warning to stderr and returns
        None.
        """
        default_entry_identifier = (entry['STUDY_REF'], entry['EXPERIMENT_ALIAS'])
        if not default_entry_identifier in self.default_entries:
            default_entry = self._create_default_entry(entry)
            self._update_entry_with_checksum(default_entry)
            if default_entry.get('CHECKSUM'):
                default_run = SraRun.from_entry(default_entry)
                default_run.register(default_entry)
                self.runs.append(default_run)
                self.default_entries.add(default_entry_identifier)

        self._update_entry_with_checksum(entry)
        if entry.get('CHECKSUM'):
            run = SraRun.from_entry(entry)
            run.register(entry)
            self.runs.append(run)

    def _create_default_entry(self, entry):
        default_entry = copy.deepcopy(entry)
        default_entry['RUN_ALIAS'] = default_entry['DEFAULT_RUN_ALIAS']
        default_entry['POOL_MEMBER_FILENAME'] = \
            default_entry['DEFAULT_SAMPLE_FILENAME']
        default_entry['POOL_MEMBER_NAME'] = ''
        return default_entry

    def _update_entry_with_checksum(self, entry, warn=True):
        """Find the data file and update the entry dictionary with its checksum.

        If the data file is not found, the CHECKSUM field is not set,
        and a warning is printed to stderr if the warning keyword
        argument is set to True.
        """
        if not entry.get('CHECKSUM'):
            filepath = os.path.join(
                self.data_dir, entry['RUN_PREFIX'], entry['POOL_MEMBER_FILENAME'])
            if os.path.exists(filepath):
                entry['CHECKSUM'] = md5_path(filepath)
            else:
                if warn:
                    stderr.write(
                        'SFF file %s was not found (maybe because no '
                        'sequences were recovered).\n' % filepath)
        return entry.get('CHECKSUM')

    def to_xml(self):
        """Create an ElementTree XML object for the SRA RUN_SET entity."""
        root = ET.Element('RUN_SET')
        root.set('xmlns:xsi', 'http://www.w3.org/2001/XMLSchema-instance')
        for run in self.runs:
            root.append(run.to_xml())
        if not self.runs:
            root.text = '\n'
        return root


class SraRun(SraEntity):
    """Class representing an SRA Run entity."""

    required_fields = {
        'alias': 'RUN_ALIAS',
        'center_name': 'RUN_CENTER',
        'run_center': 'RUN_CENTER',
        'refname': 'EXPERIMENT_ALIAS',
        'refcenter': 'STUDY_CENTER',
        }

    def __init__(self, attrs):
        """Create a new SRA Run entity."""
        super(SraRun, self).__init__(attrs)
        self.data_blocks = []

    def register(self, entry):
        """Register a new table entry with the run."""
        data_block = SraDataBlock.from_entry(entry)
        data_block.register(entry)
        data_block.serial = str(len(self.data_blocks) + 1)
        self.data_blocks.append(data_block)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA RUN entity."""
        root = ET.Element('RUN')
        self.set_xml_attributes(root, ['alias', 'center_name', 'run_center'])
        expt_elem = ET.SubElement(root, 'EXPERIMENT_REF')
        self.set_xml_attributes(expt_elem, ['refname', 'refcenter'])
        if self.data_blocks:
            for data_block in self.data_blocks:
                root.append(data_block.to_xml())
        return root


class SraDataBlock(SraEntity):
    """Class representing an SRA Data Block entity."""

    optional_fields = {
        'member_name': 'POOL_MEMBER_NAME',
        'name': 'RUN_PREFIX',
        'region': 'REGION',
        }
    
    def __init__(self, attrs):
        """Create a new SRA Data Block entity.

        The member_name attribute specifies the pool member to which
        this data block applies.  If an empty string is passed, this
        value is kept in the XML output.  An empty member_name
        attribute is used by the SRA to refer to the default sample.
        If member_name is None, the attribute is not set in the XML.

        All attributes of the DATA_BLOCK, including member_name, are
        optional according to the SRA schema.  As a convenience, the
        default value of the 'serial' attribute is set to 1.
        """
        super(SraDataBlock, self).__init__(attrs)
        self.serial = '1'
        self.files = []

    def register(self, entry):
        """Register a new table entry with the data block."""
        super(SraDataBlock, self).register(entry)
        self.files.append(SraFile.from_entry(entry))

    def to_xml(self):
        """Create an ElementTree XML object for the SRA DATA_BLOCK entity."""
        root = ET.Element('DATA_BLOCK')
        if self.member_name is not None:
            root.set('member_name', self.member_name)
        self.set_xml_attributes(root, ['name', 'serial', 'region'])
        if self.files:
            files_elem = ET.SubElement(root, 'FILES')
            for f in self.files:
                files_elem.append(f.to_xml())
        return root


class SraFile(SraEntity):
    """Class representing an SRA File entity."""

    required_fields = {
        'filename': 'POOL_MEMBER_FILENAME',
        'checksum': 'CHECKSUM',
        }

    @classmethod
    def from_entry(cls, entry):
        """Factory method to create a new file from an experiment entry."""
        return super(SraFile, cls).from_entry(entry)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA FILE element."""
        root = ET.Element('FILE')
        self.set_xml_attributes(root, ['filename', 'checksum'])
        root.set('filetype', 'sff')
        root.set('checksum_method', 'MD5')
        return root


def threecol_data_to_dict(body, warn=False):
    """Converts three-col data to dict of key: (value1, value2) , ignoring other cols"""
    result = defaultdict(list)
    for rec in body:
        try:
            key = rec[0]
            v1 = rec[1]
            v2 = rec[2]
            result[key].append((v1, v2))
        except IndexError:
            if warn:
                stderr.write(
                    'Less than 3 fields found in three-column input: %s' % rec)
    return result


def _experiment_link_xml(links):
    """Creates the EXPERIMENT_LINKS subtree for SRA Experiment XML."""
    if not links:
        return None
    root = ET.Element('EXPERIMENT_LINKS')
    for label, url in links:
        link_elem = ET.SubElement(root, 'EXPERIMENT_LINK')
        url_link_elem = ET.SubElement(link_elem, 'URL_LINK')
        label_elem = ET.SubElement(url_link_elem, 'LABEL')
        label_elem.text = label
        url_elem = ET.SubElement(url_link_elem, 'URL')
        url_elem.text = url
    return root


def _experiment_attribute_xml(attrs):
    """Creates the EXPERIMENT_ATTRIBUTES subtree for SRA Experiment XML."""
    if not attrs:
        return None
    root = ET.Element('EXPERIMENT_ATTRIBUTES')
    for tag, val in attrs:
        attr_elem = ET.SubElement(root, 'EXPERIMENT_ATTRIBUTE')
        tag_elem = ET.SubElement(attr_elem, 'TAG')
        tag_elem.text = tag
        val_elem = ET.SubElement(attr_elem, 'VALUE')
        val_elem.text = val
    return root


def pretty_xml(element, level=0, encoding=None):
    """Formats XML tree as a string with proper indentation.

    The level kwarg specifies the indentation level for the root
    element.  Child elements are indented with two additional spaces
    per level.

    If the encoding keyword argument is provided, an XML encoding
    declaration is prepended to the resultant string.  The value of
    the encoding argument is used as the encoding attribute, inside
    the declaration.
    """
    if element is None:
        return ''
    element = copy.deepcopy(element)
    __pretty_xml_helper(element, level)

    # The lxml example function places a newline and spaces at the end
    # of the parent element. We'd like these spaces to appear at the
    # beginning of the parent element, because this helps with
    # inclusion of pretty-printed XML in pre-formatted strings.  Since
    # the ElementTree library does not support pre-element text, we
    # store the whitespace from the tree and prepend it to the string
    # output.
    if element.tail is None:
        indentation = '\n'
    else:
        indentation = element.tail
        element.tail = ''

    # HACK: The ElementTree library escapes ampersand characters, but
    # not other characters that we want to escape in the XML (like >,
    # <, ", or ').  Support for this is improved in python 2.7.  For
    # now, we handle the escaping BEFORE stringifying the XML.  As a
    # consequence, we must go back and un-escape the ampersand
    # characters now.
    xml_string = indentation + ET.tostring(element).replace('&amp;', '&')

    if encoding is None:
        encoding_string = ''
    else:
        encoding_string = '<?xml version="1.0" encoding="%s"?>' % encoding
    return encoding_string + xml_string


def __pretty_xml_helper(element, level=0):
    """Private helper function for pretty_xml()

    Implements the recursive part of the pretty_xml function,
    indenting the XML tree in a mutable fashion.  It has no return
    value; the tree is indented as a side effect.  Although we copied
    it from a vendor's website, it is clear that this function is a
    shameless hack and should never be exposed outside the current
    module.

    Adapted from the lxml documentation at
    http://effbot.org/zone/element-lib.htm

    The original source code for this function is listed as public
    domain. See http://effbot.org/zone/copyright.htm for details.
    """
    i = "\n" + (level * "  ")
    if len(element):
        if not element.text or not element.text.strip():
            element.text = i + "  "
        if not element.tail or not element.tail.strip():
            element.tail = i
        for element in element:
            __pretty_xml_helper(element, level + 1)
        if not element.tail or not element.tail.strip():
            element.tail = i
    else:
        if level and (not element.tail or not element.tail.strip()):
            element.tail = i


def generate_output_fp(input_fp, ext, output_dir=None):
    """Generate new filepath by replacing the file's extension."""
    input_dir, input_filename = os.path.split(input_fp)
    basename, _ = os.path.splitext(input_filename)
    output_filename = basename + ext
    if output_dir is None:
        output_dir = input_dir
    return os.path.join(output_dir, output_filename)


def write_xml_generic(infile_path, xml_f, xml_kwargs=None, output_dir=None, output_suffix=''):
    """Writes generic xml based on contents of infilepath, returns filename."""
    if xml_kwargs is None:
        xml_kwargs = {}
    output_filepath_ext = output_suffix + '.xml'
    outfile_path = generate_output_fp(infile_path, output_filepath_ext, output_dir=output_dir)
    infile = open(infile_path, 'U')
    open(outfile_path, 'w').write(xml_f(infile, **xml_kwargs))
    return outfile_path

