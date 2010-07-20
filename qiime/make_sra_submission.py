#!/usr/bin/env python
import copy
import os
from string import strip
from collections import defaultdict
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
__version__ = "1.1.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Development"

study_links_wrapper = """    <STUDY_LINKS>%s
    </STUDY_LINKS>"""

study_link_wrapper = """
      <STUDY_LINK>
        <ENTREZ_LINK>
         <DB>pubmed</DB>
         <ID>%s</ID>
        </ENTREZ_LINK>
      </STUDY_LINK>"""

sample_attribute_wrapper = """      <SAMPLE_ATTRIBUTE> <TAG>%s</TAG> <VALUE>%s</VALUE> </SAMPLE_ATTRIBUTE>"""

sample_wrapper = """  <SAMPLE alias="%(SAMPLE_ALIAS)s">
    <TITLE>%(TITLE)s</TITLE>
    <SAMPLE_NAME>%(OPTIONAL_TITLE_FIELDS)s
    </SAMPLE_NAME>
    <DESCRIPTION>%(DESCRIPTION)s</DESCRIPTION>%(SAMPLE_ATTRIBUTES_XML)s
  </SAMPLE>"""
opt_sample_field_wrapper =  '      <%s>%s</%s>'

contact_wrapper = """    <CONTACT name="%s" inform_on_status="%s" inform_on_error="%s"/>"""

contacts_wrapper = """ <CONTACTS>
%s
 </CONTACTS>"""

action_wrapper = """   <ACTION><ADD source="%s" schema="%s" notes="%s metadata"/></ACTION>"""

actions_wrapper = """ <ACTIONS>\n%s\n   <ACTION><RELEASE/></ACTION>\n </ACTIONS>"""

file_wrapper = """ <FILES>
 <FILE filename="%s" checksum_method="MD5" checksum="%s"/>
 </FILES>"""



def generate_output_fp(input_fp, ext, output_dir=None):
    """Generate new filepath by replacing the file's extension."""
    input_dir, input_filename = os.path.split(input_fp)
    basename, _ = os.path.splitext(input_filename)
    output_filename = basename + ext
    if output_dir is None:
        output_dir = input_dir
    return os.path.join(output_dir, output_filename)

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
        line = line.rstrip('\n\r')
        # Skip lines containing only whitespace
        if not line.strip():
            continue
        if line.startswith('#'):
            line = line[1:]
            if not header:
                header = map(header_fcn, line.split('\t'))
            else:
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

def make_study_links(pmid):
    """Makes study links comment block given pmids"""
    return study_links_wrapper % (study_link_wrapper % pmid)

def twocol_data_to_dict(body, is_multiple=False):
    """Converts two-col data to dict of key-value pairs, ignoring other cols"""
    if is_multiple:
        result = defaultdict(list)
    else:
        result = {}
    for rec in body:
        try:
            if is_multiple:
                result[rec[0].strip()].append(rec[1].strip())
            else:
                result[rec[0].strip()] = rec[1].strip()
        except IndexError:
            print rec
    return result

def threecol_data_to_dict(body):
    result = defaultdict(list)
    for rec in body:
        try:
            key = rec[0]
            v1 = rec[1]
            v2 = rec[2]
            result[key].append((v1, v2))
        except IndexError:
            print rec
    return result

def row_data_to_dict(header, row):
    return dict([(k, v) for k, v in zip(header, row) if v])

def rows_data_as_dicts(header, body):
    """Iterates over rows as dicts where header has keys, each row has vals.

    Omits key-value pairs where the values are empty.
    Assumes that header will be passed as a single row.
    """
    for row in body:
        yield row_data_to_dict(header, row)

def make_study(study_lines, study_template, twocol_input_format=True):
    """Returns string for study xml."""
    header, rows = read_tabular_data(study_lines)
    if twocol_input_format:
        info = twocol_data_to_dict(rows)
    else:
        info_generator = rows_data_as_dicts(header, rows)
        info = info_generator.next()
    pmid = info.get('PMID', '').strip()
    if pmid:
        study_links_block = '\n'+make_study_links(pmid)
    else:
        study_links_block = ''
    info['XML_STUDY_LINKS_BLOCK'] = study_links_block
    return study_template % info

def make_submission(submission_lines, submission_template, docnames=None,
    submission_dir=None, twocol_input_format=True):
    """Returns string for submission xml."""
    header, rows = read_tabular_data(submission_lines)
    if twocol_input_format:
        info = twocol_data_to_dict(rows, True)
    else:
        info_generator = rows_data_as_dicts(header, rows)
        info = info_generator.next()
        if 'CONTACT' in info:
            info['CONTACT'] = info['CONTACT'].split(',')
    docnames = docnames or {}
    #build up contacts strings
    contacts = []
    if 'CONTACT' in info:
        for c in info['CONTACT']:
            name, address = map(strip, c.split(';'))
            contacts.append(contact_wrapper % (name, address, address))
        contacts_str = '\n' + (contacts_wrapper % '\n'.join(contacts)) + '\n'
    else:
        contacts_str = '\n'
    info['XML_CONTACT_BLOCK'] = contacts_str
    
    #convert info vals to scalar at this point since none multiple
    new_info = {}
    for k, v in info.items():
        if isinstance(v, list):
            try:
                new_info[k] = v[0]
            except IndexError:
                pass  # Empty list: do not make entry in new_info
        else:
            new_info[k] = v
    info = new_info
    accession = info.get('accession', '')
    if accession:
        info['ACCESSION_STRING'] = '\n accession="%s"' % accession
    else:
        info['ACCESSION_STRING'] = ''

    actions=[]
    for k, v in docnames.items():
        if v:
            actions.append(action_wrapper % (v, k, k))
    actions_str = actions_wrapper % '\n'.join(actions)
    info['XML_ACTION_BLOCK'] = actions_str

    filename = info.get('file', '')
    if filename:
        if submission_dir:
            checksum = md5_path(os.path.join(submission_dir, filename))
        else:
            checksum = md5_path(filename)
        info['XML_FILE_BLOCK'] = '\n' + file_wrapper % (filename, checksum)
    else:
        info['XML_FILE_BLOCK'] = ''
    return submission_template % info

def make_sample(sample_lines, sample_template):
    """Returns string for sample xml."""
    title_fields = ['SAMPLE_ALIAS', 'TITLE', 'TAXON_ID', 'COMMON_NAME', 
        'ANONYMIZED_NAME', 'DESCRIPTION'] #these go in the title, not the record
    optional_title_fields = ['TAXON_ID', 'COMMON_NAME', 'ANONYMIZED_NAME']
    header, body = read_tabular_data(sample_lines)
    samples = []
    for d in rows_data_as_dicts(header, body):
        attrs = [sample_attribute_wrapper % (k,v) for k, v in 
            sorted(d.items()) if not k in title_fields]
        title_attrs = [opt_sample_field_wrapper % (k, v, k) for k, v in
            sorted(d.items()) if k in optional_title_fields]
        if attrs:
            attrs = ['\n    <SAMPLE_ATTRIBUTES>'] + attrs + \
                ['    </SAMPLE_ATTRIBUTES>']
        d['SAMPLE_ATTRIBUTES_XML'] = '\n'.join(attrs)
        d['OPTIONAL_TITLE_FIELDS'] = '\n' + '\n'.join(title_attrs)
        samples.append(sample_wrapper % d)
    return sample_template % {'XML_SAMPLE_BLOCK':'\n'.join(samples)}

def group_lines_by_field(lines, field):
    """Returns dict of {state:[lines]} by state of field."""
    result = defaultdict(list)
    for line in lines:
        result[line[field]].append(line)
    if '' in result:
        del result['']
    return result

def make_run_and_experiment(experiment_lines, sff_dir, attribute_file=None,
                            link_file=None):
    """Returns strings for experiment and run xml."""
    columns, body = read_tabular_data(experiment_lines)

    run_set = SraRunSet(sff_dir)
    experiment_set = SraExperimentSet()

    for line in body:
        field_dict = dict(zip(columns, line))
        update_entry_with_deprecated_fields(field_dict)
        update_entry_with_derived_fields(field_dict)
        experiment_set.register(field_dict)
        run_set.register(field_dict)

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
        
    run_set_string = '<?xml version="1.0" encoding="UTF-8"?>' + \
                     pretty_xml(run_set.to_xml())
    experiment_set_string = '<?xml version="1.0" encoding="UTF-8"?>' + \
                            pretty_xml(experiment_set.to_xml())
    return experiment_set_string, run_set_string


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


class SraExperimentSet(SraEntity):
    """Class representing a set of SRA Experiment entities.
    """

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
        optional, including member_name.  This constructor is true to
        the SRA specification, although we can't think of a case where
        it would be useful to create a MEMBER element with no
        attributes.

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

    SraRunSet manages the list of runs, as well as default files for each run.
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
            checksum = self._update_entry_with_checksum(default_entry)
            if checksum:
                default_run = SraRun.from_entry(default_entry)
                default_run.register(default_entry)
                self.runs.append(default_run)
                self.default_entries.add(default_entry_identifier)

        checksum = self._update_entry_with_checksum(entry)
        if checksum:
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

        If the data file is not found, the CHECKSUM field is not set.

        In addition to the side-effect of updating the CHECKSUM field of the entry, the resulting value of the
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
        """Factory method to register a new file from an experiment entry."""
        return super(SraFile, cls).from_entry(entry)

    def to_xml(self):
        """Create an ElementTree XML object for the SRA FILE element."""
        root = ET.Element('FILE')
        self.set_xml_attributes(root, ['filename', 'checksum'])
        root.set('filetype', 'sff')
        root.set('checksum_method', 'MD5')
        return root


def update_entry_with_deprecated_fields(entry, warn=True):
    """Move values from deprecated fields to valid fields.

    Deprecated fields:
      - SAMPLE_ACCESSION (has become DEFAULT_SAMPLE_ACCESSION)
    """
    deprecated_fields = [
        ('SAMPLE_ACCESSION', 'DEFAULT_SAMPLE_ACCESSION'),
        ]
    for old_field, new_field in deprecated_fields:
        if old_field in entry:
            if warn:
                stderr.write(
                    'Warning: The %s field has been deprecated. Please rename '
                    'the field to %s.\n' % (old_field, new_field))
            if not entry.get(new_field):
                entry[new_field] = entry[old_field]
    return entry

def update_entry_with_derived_fields(entry):
    """Derive default values for blank/missing fields in input file.

    Derives the following fields:
      - EXPERIMENT_ALIAS (<STUDY_REF>_<RUN_PREFIX>)
      - RUN_ALIAS (<STUDY_REF>_<SAMPLE_ALIAS>_<RUN_PREFIX>)
      - BARCODE_READ_GROUP_TAG (<RUN_PREFIX>_<BARCODE>)
      - PRIMER_READ_GROUP_TAG (derived from table of standard primer
        read group tags if a primer is found, raises KeyError if
        primer is specified in entry but not found in the table)
      - POOL_MEMBER_NAME (if a primer read group tag is found,
        <RUN_PREFIX>_<SAMPLE_ALIAS>_<PRIMER_READ_GROUP_TAG>;
        otherwise, <RUN_PREFIX>_<SAMPLE_ALIAS>)
      - POOL_MEMBER_FILENAME (<POOL_MEMBER_NAME>.sff)
      - DEFAULT_SAMPLE_CENTER (<SAMPLE_CENTER>)
      - DEFAULT_SAMPLE_NAME (if default sample accession number is not
        found, <STUDY_REF>_default)
      - DEFAULT_SAMPLE_FILENAME (<STUDY_REF>_default_<RUN_PREFIX>.sff)
      - DEFAULT_RUN_ALIAS (<STUDY_REF>_default_<RUN_PREFIX>)
      - LIBRARY_STRATEGY (AMPLICON)
      - LIBRARY_SOURCE (GENOMIC)
      - LIBRARY_SELECTION (PCR)

    The optional field CHECKSUM is not handled by this function, but
    is left for the SraRunSet class to derive.
    """
    # Values of optional fields may depend on values of other optional
    # fields, so order is important.  Probably want the order to be
    # consistent with documentation.
    default_format_strings = [
        ('EXPERIMENT_ALIAS', '%(STUDY_REF)s_%(RUN_PREFIX)s'),
        ('RUN_ALIAS', '%(STUDY_REF)s_%(SAMPLE_ALIAS)s_%(RUN_PREFIX)s'),
        ('BARCODE_READ_GROUP_TAG', '%(RUN_PREFIX)s_%(BARCODE)s'),
        ]

    # Derive PRIMER_READ_GROUP_TAG, if necessary
    primer = entry.get('PRIMER')
    if primer:
        if not entry.get('PRIMER_READ_GROUP_TAG'):
            try:
                default_format_strings.append(
                    ('PRIMER_READ_GROUP_TAG', PRIMER_READ_GROUP_TAGS[primer]))
            except KeyError:
                raise KeyError(
                    'No PRIMER_READ_GROUP_TAG has been provided, and primer %s '
                    'is not found in table of standard primers. Please provide '
                    'a value for the PRIMER_READ_GROUP_TAG.' % primer)
        default_format_strings.append(
            ('POOL_MEMBER_NAME',
             '%(RUN_PREFIX)s_%(SAMPLE_ALIAS)s_%(PRIMER_READ_GROUP_TAG)s'))
    else:
        default_format_strings.append(
            ('POOL_MEMBER_NAME', '%(RUN_PREFIX)s_%(SAMPLE_ALIAS)s'))

    default_format_strings.extend([
        ('POOL_MEMBER_FILENAME', '%(POOL_MEMBER_NAME)s.sff'),
        ('DEFAULT_SAMPLE_CENTER', '%(SAMPLE_CENTER)s'),
        ])

    # Derive DEFAULT_SAMPLE_NAME
    if not entry.get('DEFAULT_SAMPLE_ACCESSION'):
        default_format_strings.append(
            ('DEFAULT_SAMPLE_NAME', '%(STUDY_REF)s_default'))

    default_format_strings.extend([
        ('DEFAULT_SAMPLE_FILENAME', '%(STUDY_REF)s_default_%(RUN_PREFIX)s.sff'),
        ('DEFAULT_RUN_ALIAS', '%(STUDY_REF)s_default_%(RUN_PREFIX)s'),
        ('LIBRARY_STRATEGY', 'AMPLICON'),
        ('LIBRARY_SOURCE', 'GENOMIC'),
        ('LIBRARY_SELECTION', 'PCR')
        ])

    for field_name, format_string in default_format_strings:
        if not entry.get(field_name):
            entry[field_name] = format_string % entry
    return entry

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
    'ATTACCGCGGCTGCTGG': 'V1-V3',
    # 926r
    'CCGTCAATTCMTTTRAGT': 'V3-V5',
    # 1492r
    'TACGGYTACCTTGTTAYGACTT': 'V6-V9',
    }


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


def pretty_xml(element, level=0):
    """Formats XML tree as a string with proper indentation.

    The level kwarg specifies the indentation level for the root
    element.  Child elements are indented with two additional spaces
    per level.
    """
    if element is None:
        return ''
    element = copy.deepcopy(element)
    __pretty_xml_helper(element, level)
    if element.tail is None:
        return '\n' + ET.tostring(element)
    # The lxml example function places a newline and spaces at the end
    # of the parent element. We'd like these spaces to appear at the
    # beginning of the parent element, because this helps with
    # inclusion of pretty-printed XML in pre-formatted strings.  Since
    # the ElementTree library does not support pre-element text, we
    # store the whitespace from the tree and prepend it to the string
    # output.
    indentation = element.tail
    element.tail = ''
    return indentation + ET.tostring(element)

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

def write_xml_generic(infile_path, template_path, xml_f, xml_kwargs=None):
    """Writes generic xml based on contents of infilepath, returns filename."""
    if xml_kwargs is None:
        xml_kwargs = {}
    template = open(template_path, 'U').read()
    base_path, ext = splitext(infile_path)
    outfile_path = base_path + '.xml'
    outfile = open(outfile_path, 'w')
    result = xml_f(open(infile_path, 'U'), template, **xml_kwargs)
    outfile.write(result)
    outfile.close()
    return outfile_path

