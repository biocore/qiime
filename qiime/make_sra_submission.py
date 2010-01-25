#!/usr/bin/env python
from string import strip
from collections import defaultdict
from md5 import md5
from optparse import OptionParser
from os.path import splitext, join
from sys import stderr

"""This script makes the submission xml files for study and experiment.

Assumes simple tab-delimited text input (allowing examples/comments; produces
xml output.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

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

run_set_wrapper = """<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">%s</RUN_SET>"""

run_wrapper = """  <RUN
    alias = "%(EXPERIMENT_ALIAS)s"
    center_name = "%(EXPERIMENT_CENTER)s" 
    run_date = "%(RUN_DATE)s"    
    run_center = "%(RUN_CENTER)s"    
    instrument_name = "%(INSTRUMENT_NAME)s"
  >
    <EXPERIMENT_REF refname="%(EXPERIMENT_ALIAS)s" refcenter="%(EXPERIMENT_CENTER)s" />%(DATA_BLOCK_XML)s    <RUN_ATTRIBUTES>
      <RUN_ATTRIBUTE>
        <TAG>notes</TAG>
        <VALUE>Submitter demultiplexed reads.  Each read was assigned to a sample pool member for those samples that yielded data. </VALUE>
      </RUN_ATTRIBUTE>
    </RUN_ATTRIBUTES>
  </RUN>"""

data_block_wrapper = """    <DATA_BLOCK
      serial = "%(MEMBER_ORDER)s"
      name = "%(RUN_PREFIX)s"
      region = "%(REGION)s"
      member_name = "%(POOL_MEMBER_NAME)s"
    >
      <FILES>
        <FILE filename="%(POOL_MEMBER_FILENAME)s" filetype="sff" checksum_method="MD5" checksum="%(CHECKSUM)s"  />
      </FILES>
    </DATA_BLOCK>"""
 

experiment_set_wrapper = """<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">%s</EXPERIMENT_SET>"""

pool_member_wrapper = """            <MEMBER refname="%(SAMPLE_ALIAS)s" refcenter="%(SAMPLE_CENTER)s" member_name="%(POOL_MEMBER_NAME)s" proportion="%(POOL_PROPORTION)s"><READ_LABEL read_group_tag="%(POOL_MEMBER_NAME)s">barcode</READ_LABEL><READ_LABEL read_group_tag="%(PRIMER_READ_GROUP_TAG)s">rRNA_primer</READ_LABEL></MEMBER>"""

basecall_wrapper = """               <BASECALL read_group_tag="%(READ_GROUP)s" min_match="%(MATCH_LEN)s" max_mismatch="%(NUM_MISMATCHES)s" match_edge="full">%(MATCH_SEQ)s</BASECALL>"""

spot_descriptor_with_linker_wrapper = """      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
          <EXPECTED_BASECALL>%(KEY_SEQ)s</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>
            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>%(BARCODE_TABLE_XML)s</EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>linker</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Linker</READ_TYPE>
            <EXPECTED_BASECALL>%(LINKER)s</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>%(PRIMER_TABLE_XML)s</EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>4</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>"""

spot_descriptor_without_linker_wrapper = """      <SPOT_DESCRIPTOR>
        <SPOT_DECODE_SPEC>
          <READ_SPEC>
            <READ_INDEX>0</READ_INDEX>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Adapter</READ_TYPE>
          <EXPECTED_BASECALL>%(KEY_SEQ)s</EXPECTED_BASECALL>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>1</READ_INDEX>            <READ_LABEL>barcode</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>BarCode</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>%(BARCODE_TABLE_XML)s</EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>2</READ_INDEX>
            <READ_LABEL>rRNA_primer</READ_LABEL>
            <READ_CLASS>Technical Read</READ_CLASS>
            <READ_TYPE>Primer</READ_TYPE>
            <EXPECTED_BASECALL_TABLE>%(PRIMER_TABLE_XML)s</EXPECTED_BASECALL_TABLE>
          </READ_SPEC>
          <READ_SPEC>
            <READ_INDEX>3</READ_INDEX>
            <READ_CLASS>Application Read</READ_CLASS>
            <READ_TYPE>Forward</READ_TYPE>
          </READ_SPEC>
        </SPOT_DECODE_SPEC>
      </SPOT_DESCRIPTOR>"""

#note: the experiment wrapper is attributing default reads at the study level, not
#at the experiment level. we might want to revisit this design decision later.
experiment_wrapper = """  <EXPERIMENT
    alias="%(EXPERIMENT_ALIAS)s"
    center_name="%(EXPERIMENT_CENTER)s"
  >
    <TITLE>%(EXPERIMENT_TITLE)s</TITLE>
    <STUDY_REF refname="%(STUDY_REF)s" refcenter="%(STUDY_CENTER)s"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>%(EXPERIMENT_DESIGN_DESCRIPTION)s</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR refname="%(STUDY_REF)s_default" refcenter="%(EXPERIMENT_CENTER)s">
        <POOL>%(POOL_MEMBERS_XML)s        </POOL>
      </SAMPLE_DESCRIPTOR>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>%(EXPERIMENT_ALIAS)s</LIBRARY_NAME>
        <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>GENOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <SINGLE></SINGLE>
        </LIBRARY_LAYOUT>
        <LIBRARY_CONSTRUCTION_PROTOCOL>
          %(LIBRARY_CONSTRUCTION_PROTOCOL)s
        </LIBRARY_CONSTRUCTION_PROTOCOL>
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
      </PROCESSING>
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



def md5_path(filename, block_size=8192):
    """"Returns md5 hash from fileame without reading whole thing into memory"""
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

def make_submission(submission_lines, submission_template, docnames=None):
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
            new_info[k] = v[0]
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
        checksum = md5_path(filename)
        info['XML_FILE_BLOCK'] = '\n'+file_wrapper % (filename, checksum)
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

def make_run_and_experiment(experiment_lines, sff_dir):
    """Returns strings for experiment and run xml."""
    header, body = read_tabular_data(experiment_lines)
    columns = [i.lstrip('#').strip() for i in header[0]]
    study_ref_index = columns.index('STUDY_REF')
    experiment_ref_index = columns.index('EXPERIMENT_ALIAS')
    studies = group_lines_by_field(body, study_ref_index)
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
            #make default sample
            default_field_dict = dict(zip(columns, experiment_lines[0]))
            default_pool_member_filename = default_field_dict['STUDY_REF'] + '_default_' + field_dict['RUN_PREFIX'] + '.sff'
            default_field_dict['POOL_MEMBER_NAME'] = ''
            default_field_dict['POOL_MEMBER_FILENAME'] = default_pool_member_filename
            barcodes = set()
            primers = set()
            linkers = set()
            pool_members = []
            primer_basecalls = []
            barcode_basecalls = []
            data_blocks = []
            MEMBER_ORDER = 1
            pool_member_list = [('',[default_field_dict])] + list(sorted(pool_member_dict.items()))
            for pool_name, pool_field_dicts in pool_member_list:
                field_dict = pool_field_dicts[0]    #assume fields not related to data blocks are identical, read from first entry
                key_seq = field_dict['KEY_SEQ']
                barcode = field_dict['BARCODE']
                if barcode not in barcodes:
                    barcodes.add(barcode)
                    barcode_basecalls.append(basecall_wrapper % {
                        'READ_GROUP':field_dict['BARCODE_READ_GROUP_TAG'],
                        'MATCH_LEN':len(barcode),
                        'NUM_MISMATCHES':0,
                        'MATCH_SEQ':barcode})
                primer = field_dict['PRIMER']
                if primer not in primers:
                    primers.add(primer)
                    primer_basecalls.append(basecall_wrapper % {
                        'READ_GROUP':field_dict['PRIMER_READ_GROUP_TAG'],
                        'MATCH_LEN':len(primer),
                        'NUM_MISMATCHES':0,
                        'MATCH_SEQ':primer})
                linker = field_dict['LINKER']
                linkers.add(linker)
                #create and append the pool member
                field_dict['MEMBER_ORDER'] = MEMBER_ORDER
                if pool_name:
                    pool_members.append(pool_member_wrapper % field_dict)
                #create and append the data blocks
                for f in pool_field_dicts:
                    f['MEMBER_ORDER'] = MEMBER_ORDER
                    if not f.get('POOL_MEMBER_FILENAME',''):
                        f['POOL_MEMBER_FILENAME'] = f['POOL_MEMBER_NAME'] + '.sff'
                    try:
                        f['CHECKSUM'] = md5_path(join(sff_dir,f['RUN_PREFIX'],f['POOL_MEMBER_FILENAME']))
                        data_blocks.append(data_block_wrapper % f)
                        MEMBER_ORDER += 1   #skip members where we couldn't find the file
                    except IOError: #file missing, probably because no seqs were recovered
                        stderr.write("File failed with IOError:\n%s\n" % f['POOL_MEMBER_FILENAME'])
                        pass
                                    
            field_dict['BARCODE_TABLE_XML'] = '\n' + '\n'.join(barcode_basecalls) + '\n'
            field_dict['PRIMER_TABLE_XML'] = '\n' + '\n'.join(primer_basecalls) + '\n'
            field_dict['POOL_MEMBERS_XML'] = '\n' + '\n'.join(pool_members) + '\n'
            field_dict['DATA_BLOCK_XML'] = '\n' + '\n'.join(data_blocks) + '\n'
            runs.append(run_wrapper % field_dict)
            if linkers == (''): #not using linkers
               spot_descriptor_wrapper = spot_descriptor_without_linker_wrapper
            else:
                spot_descriptor_wrapper = spot_descriptor_with_linker_wrapper
            field_dict['TOTAL_TECHNICAL_READ_LENGTH'] = len(key_seq) + len(primer) + len(barcode) + len(linker) + 1 #note that SRA uses 1-indexed lengths
            spot_descriptor = spot_descriptor_wrapper % field_dict
            field_dict['SPOT_DESCRIPTORS_XML'] = spot_descriptor
            field_dict['PLATFORM_XML'] = platform_blocks[field_dict['PLATFORM']]
            experiments.append(experiment_wrapper % field_dict)
    return experiment_set_wrapper % ('\n'+'\n'.join(experiments)+'\n'), run_set_wrapper % ('\n'+'\n'.join(runs)+'\n')

usage_str = """usage: %prog [options]

[] indicates optional input (order unimportant)

Example usage:
 Read the sample data from sample.txt, the study data from study.txt,
 and the submission data from submission.txt. Write out the corresponding
 XML files.
 
 python make_study_and_sample_xml.py -a sample.txt -A sample_template.xml -t study.txt -T study_template.xml -u submission.txt -U submission_template.xml

 Produces files study.xml, submission.xml, sample.xml (based on filenames of
 the .txt files).
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-a','--input_sample_fp',\
        help='the tab-delimited text file with info about samples')

    parser.add_option('-A','--template_sample_fp',\
        help='the template file for samples [default: %default]')

    parser.add_option('-t','--input_study_fp',\
        help='the tab-delimited text file with info about the study')

    parser.add_option('-T','--template_study_fp',\
        help='the template file for the study [default: %default]')

    parser.add_option('-u','--input_submission_fp',\
        help='the tab-delimited text file with info about the submission')

    parser.add_option('-U','--template_submission_fp',\
        help='the template file for the submission [default: %default]')

    parser.add_option('-e', '--input_experiment_fp', \
        help ='the tab-delimited text file with info about the experiment')

    parser.add_option('-s', '--sff_dir', \
        help = 'the directory containing the demultiplexed sff files: 1 dir per run')

    # Set default values here if they should be other than None
    parser.set_defaults(template_sample_fp='sample_template.xml', 
        template_study_fp='study_template.xml', 
        template_submission_fp='submission_template.xml')

    opts,args = parser.parse_args()

    return opts,args

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


if __name__ == '__main__':
    #if run from the command-line, produce the appropriate xml files
    opts, args = parse_command_line_parameters()
    docnames = {}
    if opts.input_study_fp:
        docnames['study'] = write_xml_generic(opts.input_study_fp,
            opts.template_study_fp, make_study)

    if opts.input_sample_fp:
        docnames['sample'] = write_xml_generic(opts.input_sample_fp,
            opts.template_sample_fp, make_sample)

    if opts.input_experiment_fp:
        #in this case, we need to need to also get the sff dir
        if not opts.sff_dir:
            raise IOError, "Must specify an sff dir if making an experiment."
        base_name, ext = splitext(opts.input_experiment_fp)
        base_name = base_name.split('experiment')[-1]
        if base_name:
            base_name += '_'
        run_path = base_name + 'run.xml'
        run_file = open(run_path, 'w')
        experiment_path = base_name + 'experiment.xml'
        experiment_file = open(experiment_path, 'w')
        experiment_xml, run_xml = make_run_and_experiment(
            open(opts.input_experiment_fp, 'U'), opts.sff_dir)
        run_file.write(run_xml)
        experiment_file.write(experiment_xml)
        run_file.close()
        experiment_file.close()
        docnames['run'] = run_path
        docnames['experiment'] = experiment_path

    if opts.input_submission_fp:
        submission_template = open(opts.template_submission_fp, 'U').read()
        base_name, ext = splitext(opts.input_submission_fp)
        outfilename = base_name + '.xml'
        outfile = open(outfilename, 'w')
        outfile.write(make_submission(open(opts.input_submission_fp, 'U'), 
            submission_template, docnames))
        outfile.close()


