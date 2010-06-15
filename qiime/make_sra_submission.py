#!/usr/bin/env python
import copy
import os
from string import strip
from collections import defaultdict
from hashlib import md5
from os.path import splitext, join
from sys import stderr
import xml.etree.ElementTree as ET

"""This script makes the submission xml files for SRA (study, experiment, etc.).

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

run_set_wrapper = '''\
<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">%s</RUN_SET>'''

run_wrapper = '''\
  <RUN alias="%(RUN_ALIAS)s" center_name="%(RUN_CENTER)s" run_center="%(RUN_CENTER)s">
    <EXPERIMENT_REF refname="%(EXPERIMENT_ALIAS)s" refcenter="%(STUDY_CENTER)s"/>%(DATA_BLOCK_XML)s
  </RUN>'''

experiment_set_wrapper = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">%s</EXPERIMENT_SET>"""

#note: the experiment wrapper is attributing default reads at the study level, not
#at the experiment level. we might want to revisit this design decision later.
experiment_wrapper = """\
  <EXPERIMENT alias="%(EXPERIMENT_ALIAS)s" center_name="%(EXPERIMENT_CENTER)s">
    <TITLE>%(EXPERIMENT_TITLE)s</TITLE>
    <STUDY_REF refname="%(STUDY_REF)s" refcenter="%(SAMPLE_CENTER)s"%(STUDY_ACCESSION_ATTRIBUTE)s/>
    <DESIGN>
      <DESIGN_DESCRIPTION>%(EXPERIMENT_DESIGN_DESCRIPTION)s</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR%(DEFAULT_SAMPLE_NAME_ATTRIBUTE)s refcenter="%(DEFAULT_SAMPLE_CENTER)s"%(DEFAULT_SAMPLE_ACCESSION_ATTRIBUTE)s>
        <POOL>%(POOL_MEMBERS_XML)s        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>%(EXPERIMENT_ALIAS)s</LIBRARY_NAME>
        <LIBRARY_STRATEGY>%(LIBRARY_STRATEGY)s</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>%(LIBRARY_SOURCE)s</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>%(LIBRARY_SELECTION)s</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>%(LIBRARY_CONSTRUCTION_PROTOCOL)s</LIBRARY_CONSTRUCTION_PROTOCOL>
      </LIBRARY_DESCRIPTOR>%(SPOT_DESCRIPTORS_XML)s
      </DESIGN>
      %(PLATFORM_XML)s
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
      </PROCESSING>%(LINK_XML)s%(ATTRIBUTE_XML)s
  </EXPERIMENT>"""

platform_blocks = { 'Titanium':
"""    <PLATFORM>
        <LS454>
            <INSTRUMENT_MODEL>454 Titanium</INSTRUMENT_MODEL>
            <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
            <FLOW_COUNT>800</FLOW_COUNT>
        </LS454>
    </PLATFORM>""",
                    'FLX':
"""    <PLATFORM>
        <LS454>
            <INSTRUMENT_MODEL>454 GS FLX</INSTRUMENT_MODEL>
            <FLOW_SEQUENCE>TACG</FLOW_SEQUENCE>
            <FLOW_COUNT>400</FLOW_COUNT>
        </LS454>
    </PLATFORM>"""}


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

def read_tabular_data(lines):
    """Reads tabular data from lines, skipping blanks"""
    header, body = [], []
    for line in lines:
        line = line.rstrip('\n')
        if not line:
            continue
        fields = [safe_for_xml(i.strip('"')) for i in line.split('\t')]
        if line.startswith('#'):
            header.append(fields)
        else:
            body.append(fields)
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

def trim_quotes(s):
    """Trims quotes and whitespace from string."""
    return s.strip('"').strip()

def rows_data_as_dicts(header, body):
    """Iterates over rows as dicts where header has keys, each row has vals.

    Omits key-value pairs where the values are empty.
    Assumes that header will be passed as a single row.
    """
    header = [h.strip('#').strip('"') for h in header] #get rid of junk
    for row in body:
        yield dict([(k, v) for k, v in zip(header, map(trim_quotes, row)) if v])

def make_study(study_lines, study_template):
    """Returns string for study xml."""
    info = twocol_data_to_dict(read_tabular_data(study_lines)[1])
    pmid = info.get('PMID', '').strip()
    if pmid:
        study_links_block = '\n'+make_study_links(pmid)
    else:
        study_links_block = ''
    info['XML_STUDY_LINKS_BLOCK'] = study_links_block
    return study_template % info

def make_submission(submission_lines, submission_template, docnames=None,
    submission_dir=None):
    """Returns string for submission xml."""
    docnames = docnames or {}
    info = twocol_data_to_dict(read_tabular_data(submission_lines)[1], True)
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
    for d in rows_data_as_dicts(header[0], body):
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
    header, body = read_tabular_data(experiment_lines)
    columns = [i.lstrip('#').strip() for i in header[0]]
    study_ref_index = columns.index('STUDY_REF')
    experiment_ref_index = columns.index('EXPERIMENT_ALIAS')
    studies = group_lines_by_field(body, study_ref_index)

    # Process attribute file -- to be refactored.
    experiment_attributes = defaultdict(list)
    if attribute_file is not None:
        _, attr_body = read_tabular_data(attribute_file)
        for elems in attr_body:
            try:
                expt_name = elems[0]
                attribute_name = elems[1]
                attribute_val = elems[2]
                experiment_attributes[expt_name].append(
                    (attribute_name, attribute_val))
            except IndexError:
                print ("Not enough items in attribute line (%s), must define "
                       "experiment name, attribute name, and attribute value. "
                       "Skipping..." % elems)

    # Process links file -- to be refactored.
    experiment_links = defaultdict(list)
    if link_file is not None:
        _, attr_body = read_tabular_data(link_file)
        for elems in attr_body:
            try:
                expt_name = elems[0]
                link_name = elems[1]
                link_url = elems[2]
                experiment_links[expt_name].append(
                    (link_name, link_url))
            except IndexError:
                print ("Not enough items in link line (%s), must define "
                       "experiment name, link name, and link url. "
                       "Skipping..." % elems)

    experiments = []
    runs = []
    for study_id, study_lines in studies.items():
        experiment_groups = group_lines_by_field(study_lines, experiment_ref_index)
        field_dict = {} # to keep the last one in scope for outer block
        for experiment_id, experiment_lines in experiment_groups.items():
            #collect unique pool members
            pool_member_dict = defaultdict(list)
            for line in experiment_lines:
                field_dict = dict(zip(columns, line))
                pool_member_dict[field_dict['POOL_MEMBER_NAME']].append(field_dict)

            # Set up default sample using optional fields
            if not field_dict.get('DEFAULT_SAMPLE_CENTER'):
                field_dict['DEFAULT_SAMPLE_CENTER'] = field_dict['SAMPLE_CENTER']

            # Still use SAMPLE_ACCESSION field, but announce
            # deprecation in favor of DEFAULT_SAMPLE_ACCESSION
            if 'SAMPLE_ACCESSION' in field_dict:
                stderr.write(
                    'Warning: The SAMPLE_ACCESSION field has been deprecated. '
                    'Please rename the field to DEFAULT_SAMPLE_ACCESSION.\n')
                sample_acc_DEPRECATED = field_dict.get('SAMPLE_ACCESSION')
                if sample_acc_DEPRECATED and (not field_dict.get('DEFAULT_SAMPLE_ACCESSION')):
                    field_dict['DEFAULT_SAMPLE_ACCESSION'] = sample_acc_DEPRECATED
            default_acc = field_dict.get('DEFAULT_SAMPLE_ACCESSION')
            field_dict['DEFAULT_SAMPLE_ACCESSION_ATTRIBUTE'] = (
                ' accession="%s"' % default_acc if default_acc else '')

            # If necessary, derive a default sample name
            default_name = field_dict.get('DEFAULT_SAMPLE_NAME')
            # If a default accession number has been provided, we
            # should not derive the default sample name automatically.
            if (not default_name) and (not default_acc):
                default_name = field_dict['STUDY_REF'] + '_default'
            field_dict['DEFAULT_SAMPLE_NAME_ATTRIBUTE'] = (
                ' refname="%s"' % default_name if default_name else '')

            #make default pool member dict
            default_field_dict = dict(zip(columns, experiment_lines[0]))
            default_pool_member_id = field_dict['STUDY_REF'] + '_default_' + field_dict['RUN_PREFIX']
            default_field_dict['RUN_ALIAS'] = default_pool_member_id
            default_field_dict['POOL_MEMBER_FILENAME'] = default_pool_member_id + '.sff'
            default_field_dict['POOL_MEMBER_NAME'] = ''

            barcode_basecall_table = []
            primer_basecall_table = []
            barcodes = set()
            primers = set()
            pool_members = []
            data_blocks = []
            MEMBER_ORDER = 1

            pool_member_list = [('',[default_field_dict])] + list(sorted(pool_member_dict.items()))
            for pool_name, pool_field_dicts in pool_member_list:
                # Assume fields not related to data blocks are
                # identical, read from first entry
                field_dict = pool_field_dicts[0]

                key_seq = field_dict['KEY_SEQ']
                barcode = field_dict['BARCODE']
                if barcode not in barcodes:
                    barcodes.add(barcode)
                    barcode_basecall_table.append(
                        (field_dict['BARCODE_READ_GROUP_TAG'], barcode))
                primer = field_dict['PRIMER']
                if primer and primer not in primers:
                    primers.add(primer)
                    primer_basecall_table.append(
                        (field_dict['PRIMER_READ_GROUP_TAG'], primer))
                linker = field_dict['LINKER']

                if pool_name:
                    pool_member_xml = '          ' + ET.tostring(_pool_member_xml(
                        refname=field_dict.get('SAMPLE_ALIAS'),
                        refcenter=field_dict.get('SAMPLE_CENTER'),
                        member_name=field_dict.get('POOL_MEMBER_NAME'),
                        proportion=field_dict.get('POOL_PROPORTION'),
                        barcode_tag=field_dict.get('BARCODE_READ_GROUP_TAG'),
                        primer_tag=field_dict.get('PRIMER_READ_GROUP_TAG'),
                        accession=field_dict.get('POOL_MEMBER_ACCESSION'),
                        ))
                    pool_members.append(pool_member_xml)

                ####################################################
                # RUN XML
                ####################################################
                for pool_field_dict in pool_field_dicts:
                    data_block = _make_data_block(pool_field_dict, sff_dir)
                    if data_block is not None:
                        field_dict['DATA_BLOCK_XML'] = pretty_xml(data_block, 2)
                        runs.append(run_wrapper % field_dict)

            ########################
            # EXPERIMENT XML
            ########################
            field_dict['POOL_MEMBERS_XML'] = '\n' + '\n'.join(pool_members) + '\n'

            spot_descriptor = pretty_xml(_spot_descriptor_xml(
                key_seq, barcode_basecall_table, linker, primer_basecall_table), 3)
            field_dict['SPOT_DESCRIPTORS_XML'] = spot_descriptor

            field_dict['PLATFORM_XML'] = platform_blocks[field_dict['PLATFORM']]

            field_dict['ATTRIBUTE_XML'] = pretty_xml(
                _experiment_attribute_xml(experiment_attributes[experiment_id]), 3)
            field_dict['LINK_XML'] = pretty_xml(
                _experiment_link_xml(experiment_links[experiment_id]), 3)

            # Insert study accession attribute, if present.
            study_acc = field_dict.get('STUDY_ACCESSION')
            field_dict['STUDY_ACCESSION_ATTRIBUTE'] = (
                ' accession="%s"' % study_acc if study_acc else '')

            # Insert sample accession attribute, if present.
            sample_acc = field_dict.get('SAMPLE_ACCESSION')
            field_dict['SAMPLE_ACCESSION_ATTRIBUTE'] = (
                ' accession="%s"' % sample_acc if sample_acc else '')

            # Utilize optional fields for library descriptor block
            if 'LIBRARY_SELECTION' not in field_dict:
                field_dict['LIBRARY_SELECTION'] = 'PCR'
            if 'LIBRARY_STRATEGY' not in field_dict:
                field_dict['LIBRARY_STRATEGY'] = 'AMPLICON'
            if 'LIBRARY_SOURCE' not in field_dict:
                field_dict['LIBRARY_SOURCE'] = 'GENOMIC'

            experiments.append(experiment_wrapper % field_dict)
    return experiment_set_wrapper % ('\n'+'\n'.join(experiments)+'\n'), run_set_wrapper % ('\n'+'\n'.join(runs)+'\n')

def derive_default_run_fields(experiment_entry):
    """Derive default values for missing or blank fields in Run XML.

    Derives the following fields:
    POOL_MEMBER_FILENAME (<POOL_MEMBER_NAME>.sff)
    """
    vals = copy.deepcopy(experiment_entry)
    if not vals.get('POOL_MEMBER_FILENAME'):
        vals['POOL_MEMBER_FILENAME'] = vals['POOL_MEMBER_NAME'] + '.sff'
    return vals

def _make_data_block(experiment_entry, sff_dir):
    """Make a data block from a single line in the experiment table.

    Returns an ElementTree XML object representing the RUN_BLOCK, or
    None if the SFF files were not found.
    """
    vals = derive_default_run_fields(experiment_entry)
    relative_sff_path = join(
        sff_dir, vals['RUN_PREFIX'], vals['POOL_MEMBER_FILENAME'])
    if os.path.exists(relative_sff_path):
        md5sum = md5_path(relative_sff_path)
        data_block_file = (vals['POOL_MEMBER_FILENAME'], md5sum)
        return _data_block_xml(
            vals.get('POOL_MEMBER_NAME'),
            name=vals.get('RUN_PREFIX'),
            region=vals.get('REGION'),
            files=[data_block_file],
            )
    else:
        stderr.write(
            'Pool member file %s not found, maybe because no sequences '
            'were recovered.\n' % relative_sff_path)
        return None

def _data_block_xml(
    member_name, serial='1', name=None, region=None, files=[]):
    """Create a DATA_BLOCK subtree for SRA Run XML.

    The member_name argument specifies the pool member to which this
    data block applies.  If an empty string is passed, this value is
    kept in the XML output.  An empty member_name attribute is used by
    the SRA to refer to the default sample.  If member_name is None,
    the attribute is not set.

    The remaining attributes of the DATA_BLOCK element are provided as
    optional keyword arguments.  All attributes, including
    member_name, are optional according to the SRA schema.  As a
    convenience, the default value of the 'serial' attribute is set to
    1.

    FILE elements are created from the list passed via the 'files'
    keyword -- files should be provided as a list of (filename,
    md5sum) tuples.
    """
    root = ET.Element('DATA_BLOCK')
    if member_name is not None:
        root.set('member_name', member_name)
    block_attributes = {'serial': serial, 'name': name, 'region': region}
    for attr, val in block_attributes.items():
        if val:
            root.set(attr, val)
    if files:
        files_elem = ET.SubElement(root, 'FILES')
        for filename, md5sum in files:
            ET.SubElement(files_elem, 'FILE', filename=filename, filetype='sff',
                          checksum_method='MD5', checksum=md5sum)
    return root

def _pool_member_xml(
    member_name=None, proportion=None, refcenter=None, refname=None,
    barcode_tag=None, primer_tag=None, accession=None):
    """Creates a MEMBER subtree for SRA Experiment XML.

    The attributes of the MEMBER element are provided as optional
    keyword arguments.  According to the SRA schema, absolutely all
    MEMBER attributes are optional, including member_name.  This
    function implementation is true to the SRA specification, although
    we can't think of a case where it would be useful to create a
    MEMBER element with no attributes.

    READ_LABELs are created using the keyword args barcode_tag and
    primer_tag.  According to the SRA schema, a MEMBER element must
    contain at least one READ_LABEL.  Therefore, this function raises
    a ValueError if barcode_tag and primer_tag are both None.
    """
    if not (barcode_tag or primer_tag):
        raise ValueError('Must provide either barcode_tag or primer_tag.')
    member_attributes = {
        'member_name': member_name,
        'refcenter': refcenter,
        'refname': refname,
        'proportion': proportion,
        'accession': accession,
        }
    root = ET.Element('MEMBER')
    for attr, val in member_attributes.items():
        if val:
            root.set(attr, val)
    if barcode_tag:
        barcode_elem = ET.SubElement(
            root, 'READ_LABEL', read_group_tag=barcode_tag)
        barcode_elem.text = 'barcode'
    if primer_tag:
        primer_elem = ET.SubElement(
            root, 'READ_LABEL', read_group_tag=primer_tag)
        primer_elem.text = 'rRNA_primer'
    return root

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

def _spot_descriptor_xml(adapter, barcodes, linker, primers):
    """Creates the SPOT_DESCRIPTOR subtree for SRA Experiment XML.

    The adaptor and linker are to be passed as strings.  Barcodes and
    primers should be passed as a list of tuples, mapping read group
    tags to basecalls.
    """
    root = ET.Element('SPOT_DESCRIPTOR')
    decode_elem = ET.SubElement(root, 'SPOT_DECODE_SPEC')
    read_index = 0
    decode_elem.append(_read_spec_xml(
        read_index, 'Technical Read', 'Adapter', expected_basecall=adapter))
    read_index += 1
    if barcodes:
        decode_elem.append(_read_spec_xml(
            read_index, 'Technical Read', 'BarCode', read_label='barcode',
            expected_basecall_table=barcodes))
        read_index += 1
    if linker:
        decode_elem.append(_read_spec_xml(
            read_index, 'Technical Read', 'Linker', read_label='linker',
            expected_basecall=linker))
        read_index += 1
    if primers:
        decode_elem.append(_read_spec_xml(
            read_index, 'Technical Read', 'Primer', read_label='rRNA_primer',
            expected_basecall_table=primers))
        read_index += 1
    decode_elem.append(_read_spec_xml(
        read_index, 'Application Read', 'Forward',
        in_relative_order=True))
    return root

def _read_spec_xml(
    read_index, read_class, read_type, read_label=None,
    expected_basecall=None, expected_basecall_table=None,
    in_relative_order=False):
    """Returns XML for the READ_SPEC element of an SRA experiment.

    The read_index, read_class, and read_type arguments specify the
    text of the corresponding XML subelements of READ_SPEC.

    The read_label element is used for barcodes and primers, but not
    the technical or application reads.

    Expected basecalls may be provided as a single string
    (expected_basecall) or as a sequence of tuples mapping read group
    tags to basecalls (expected_basecall_table).  If both are
    provided, the expected_basecall argument is used and the
    expected_basecall_table is discarded.

    For the application read spec, it is customary to provide a
    relative_order.  Passing a True value to the in_relative_order
    keyword argument will generate this element with the correct
    value.
    """
    root = ET.Element('READ_SPEC')
    ET.SubElement(root, 'READ_INDEX').text = str(read_index)
    # According to SRA.experiment.xsd v1.1, the READ_LABEL element
    # must appear immediately following the READ_INDEX
    if read_label:
        ET.SubElement(root, 'READ_LABEL').text = read_label
    ET.SubElement(root, 'READ_CLASS').text = read_class
    ET.SubElement(root, 'READ_TYPE').text = read_type
    if expected_basecall:
        ET.SubElement(root, 'EXPECTED_BASECALL').text = expected_basecall
    elif expected_basecall_table:
        table_elem = ET.SubElement(root, 'EXPECTED_BASECALL_TABLE')
        for read_group_tag, basecall in expected_basecall_table:
            min_match = len(basecall)
            basecall_elem = ET.SubElement(table_elem, 'BASECALL')
            basecall_elem.set('read_group_tag', read_group_tag)
            basecall_elem.set('min_match', str(min_match))
            basecall_elem.set('max_mismatch', '0')
            basecall_elem.set('match_edge', 'full')
            basecall_elem.text = basecall
    if in_relative_order:
        prevoius_index = int(read_index) - 1
        order_elem = ET.SubElement(root, 'RELATIVE_ORDER')
        order_elem.set('follows_read_index', str(prevoius_index))
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

def write_xml_generic(infile_path, template_path, xml_f):
    """Writes generic xml based on contents of infilepath, returns filename."""
    template = open(template_path, 'U').read()
    base_path, ext = splitext(infile_path)
    outfile_path = base_path + '.xml'
    outfile = open(outfile_path, 'w')
    result = xml_f(open(infile_path, 'U'), template)
    outfile.write(result)
    outfile.close()
    return outfile_path

