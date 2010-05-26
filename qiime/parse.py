#!/usr/bin/env python
#file parse.py: parsers for map file, distance matrix file, env file

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso",
    "Justin Kuczynski", "Cathy Lozupone", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

from string import strip
from collections import defaultdict
from cogent.util.misc import unzip
from cogent.maths.stats.rarefaction import subsample
from numpy import array, concatenate, repeat, zeros, nan
from numpy.random import permutation
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.fasta import FastaFinder
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from copy import deepcopy
import os
from cogent.util.misc import revComp


class QiimeParseError(Exception):
    pass

class IlluminaParseError(QiimeParseError):
    pass
    
def parse_newick(lines, constructor=PhyloNode):
    """Return PhyloNode from newick file handle stripping quotes from tip names
    
        This function wraps cogent.parse.tree.DndParser stripping 
         matched leading/trailing single quotes from tip names, and returning
         a PhyloNode object by default (alternate constructor can be passed 
         with constructor=).
         
        Sripping of quotes is essential for many applications in Qiime, as 
         the tip names are frequently matched to OTU ids, and if the tip name
         is read in with leading/trailing quotes, node.Name won't match to the
         corresponding OTU identifier. Disaster follows.
        
    """
    return DndParser(lines, constructor=constructor, unescape_name=True) 

def parse_mapping_file(lines, strip_quotes=True, suppress_stripping=False):
    """Parser for map file that relates samples to metadata.
    
    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    Result: list of lists of fields, incl. headers.
    """
    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"','')
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"','').strip()
    else:
        if suppress_stripping:
            # don't remove quotes or spaces
            strip_f = lambda x: x
        else:
            # remove spaces but not quotes
            strip_f = lambda x: x.strip()
    
    # Create lists to store the results
    mapping_data = []
    header = []
    comments = []
    
    # Begin iterating over lines
    for line in lines:
        line = strip_f(line)
        if not line or (suppress_stripping and not line.strip()):
            # skip blank lines when not stripping lines
            continue
        
        if line.startswith('#'):
            line = line[1:]
            if not header:
                header = line.strip().split('\t')
            else:
                comments.append(line)
        else:
            mapping_data.append(map(strip_f, line.split('\t')))

    return mapping_data, header, comments


def group_by_field(table, name):
    """Returns dict of field_state:[row_headers] from table.

    Use to extract info from table based on a single field.
    """
    try:
        col_index = table[0].index(name)
    except ValueError, e:
        raise ValueError, "Couldn't find name %s in headers: %s" % \
            (name, table[0])
    result = defaultdict(list)
    for row in table[1:]:
        header, state = row[0], row[col_index]
        result[state].append(header)
    return result

def group_by_fields(table, names):
    """Returns dict of (field_states):[row_headers] from table.

    Use to extract info from table based on combinations of fields.
    """
    col_indices = map(table[0].index, names)
    result = defaultdict(list)
    for row in table[1:]:
        header = row[0]
        states = tuple([row[i] for i in col_indices])
        result[states].append(header)
    return result

def parse_distmat(lines):
    """Parser for distance matrix file (e.g. UniFrac dist matrix).

    The examples I have of this file are just sample x sample tab-delimited
    text, so easiest way to handle is just to convert into a numpy array
    plus a list of field names.
    """
    header = None
    result = []
    for line in lines:
        if line[0] == '\t': #is header
            header = map(strip, line.split('\t')[1:])
        else:
            result.append(map(float, line.split('\t')[1:]))
    return header, array(result)

def parse_matrix(lines):
    """Parser for a matrix file Tab delimited. skips first lines if led
    by '#', assumes column headers line starts with a tab
    """
    col_headers = None
    result = []
    row_headers = []
    for line in lines:
        if line[0] == '#': continue
        if line[0] == '\t': #is header
            col_headers = map(strip, line.split('\t')[1:])
        else:
            entries = line.split('\t')
            result.append(map(float, entries[1:]))
            row_headers.append(entries[0])
    return col_headers, row_headers, array(result)

def parse_distmat_to_dict(table):
    """Parse a dist matrix into an 2d dict indexed by sample ids.
    
    table: table as lines
    """
    
    col_headers, row_headers, data = parse_matrix(table)
    
    assert(col_headers==row_headers)
    
    result = defaultdict(dict)
    for (sample_id_x, row) in zip (col_headers,data):
        for (sample_id_y, value) in zip(row_headers, row):
            result[sample_id_x][sample_id_y] = value
    return result
    
def parse_bootstrap_support(lines):
    """Parser for a bootstrap/jackknife support in tab delimited text
    """
    bootstraps = {}
    for line in lines:
        if line[0] == '#': continue
        wordlist = line.strip().split()
        bootstraps[wordlist[0]] = float(wordlist[1])
        
    return bootstraps

def parse_rarefaction_data(lines):
    data = {}
    data['headers'] = []
    data['options'] = []
    data['xaxis'] = []
    data['series'] = {}
    data['error'] = {}
    data['color'] = {}
    for l in lines:
        if l.startswith('#'):
            data['headers'].append(l.strip('#').strip())
            continue
        if l.startswith('xaxis'):
            data['xaxis'] = [float(v) for v in l[6:].strip().split('\t')]
            continue
        if l.startswith('>>'):
            data['options'].append(l.strip('>').strip())
            continue
        if l.startswith('series'):
            data['series'][data['options'][len(data['options'])-1]] = \
            [float(v) for v in l[7:].strip().split('\t')]
            continue
        if l.startswith('error'):
            data['error'][data['options'][len(data['options'])-1]] = \
            [float(v) for v in l[6:].strip().split('\t')]
        if l.startswith('color'):
            data['color'][data['options'][len(data['options'])-1]] = \
            str(l[6:].strip())
            if(len(str(l[6:].strip())) < 1):
                print data['options'][len(data['options'])-1]
    return data

def parse_rarefaction_record(line):
    """ Return (rarefaction_fn, [data])"""       
    
    def float_or_nan(v):
        try:
            return float(v)
        except ValueError:
            return nan
            
    entries = line.split('\t')
    return entries[0], map(float_or_nan, entries[1:])

def parse_rarefaction(lines):
    """Function for parsing rarefaction files specifically for use in
    make_rarefaction_plots.py"""
    col_headers = []
    comments = []
    rarefaction_data = []
    rarefaction_fns = []
    for line in lines:
        if line[0] == '#':
            # is comment
            comments.append(line)
        elif line[0] == '\t': 
            # is header
            col_headers = map(strip, line.split('\t'))
        else:
            # is rarefaction record
            rarefaction_fn, data = parse_rarefaction_record(line)
            rarefaction_fns.append(rarefaction_fn)
            rarefaction_data.append(data)

    return col_headers, comments, rarefaction_fns, rarefaction_data

def parse_coords(lines):
    """Parse unifrac coord file into coords, labels, eigvals, pct_explained.

    Returns:
    - list of sample labels in order
    - array of coords (rows = samples, cols = axes in descending order)
    - list of eigenvalues
    - list of percent variance explained

    File format is tab-delimited with following contents:
    - header line (starts 'pc vector number')
    - one-per-line per-sample coords
    - two blank lines
    - eigvals
    - % variation explained

    Strategy: just read the file into memory, find the lines we want
    """
    lines = list(lines)
    lines = map(strip, lines[1:])   #discard first line, which is a label
    lines = filter(None, lines) #remove any blank lines
    
    #now last 2 lines are eigvals and % variation, so read them
    eigvals = array(map(float, lines[-2].split('\t')[1:]))
    pct_var = array(map(float, lines[-1].split('\t')[1:]))
    
    #finally, dump the rest of the lines into a table
    header, result = [], []
    for line in lines[:-2]:
        fields = map(strip, line.split('\t'))
        header.append(fields[0])
        result.append(map(float, fields[1:]))

    return header, array(result), eigvals, pct_var

def parse_rarefaction_fname(name_string):
    """returns base, seqs/sam, iteration, extension.  seqs, iters as ints
    
    all as strings, some may be empty strings ('')"""

    root, ext = os.path.splitext(name_string)
    root_list = root.split("_")
    iters = int(root_list.pop())
    seqs_per_sam = int(root_list.pop())
    base_name = "_".join(root_list)
    return base_name, seqs_per_sam, iters, ext

def parse_taxonomy(infile):
    """parse a taxonomy file.


    Typically the lines in these files look like: 
      3 SAM1_32 \t Root;Bacteria;Fi... \t 0.9
      
     where the first field is the sequence identifier, the second field is the 
      taxonomy assignment separated by ; characters, and the third field is a
      quality score (e.g., confidence from the RDP classifier)
     
     when using the BLAST taxonomy assigner, an additional field is included, 
      containing the sequence identifier of the best blast hit or each input 
      sequence. these lines might look like:
      3 SAM1_32 \t Root;Bacteria;Fi... \t 1e-42 \t A1237756

    Returns: dict of otu id to taxonomy name.
    ignores other parts of the otu file, such as confidence and seq id (otu id
    only)
    """

    res = {}
    for line in infile:
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        otu = fields[0].split(' ')[0]
        res[otu] = fields[1]

    return res

def parse_otu_table(lines,count_map_f=int):
    """parses otu file

    Returns tuple: sample_ids, otu_ids, matrix of OTUs(rows) x samples(cols),
    and lineages from infile."""
    otu_table = []
    otu_ids = []
    lineages = []
    for i, line in enumerate(lines):
        if i == 0: continue # skip header line
        
        elif i == 1: # parse sample id line
            sample_ids = line.strip().split('\t')[1:]
            if len(sample_ids) == 0:
                    raise RuntimeError('no samples found in otu table')
            if sample_ids[-1] == 'Consensus Lineage':
                has_consensus = True
                sample_ids = sample_ids[:-1]
            else:
                has_consensus = False
                
        else: # parse each otu line
            line = line.strip()
            if line:
                fields = line.split('\t')
                #first and last col are otu id and consensus lineage respectively
                if has_consensus:
                    otu_table.append(array(map(count_map_f, fields[1:-1])))
                else:
                    otu_table.append(array(map(count_map_f, fields[1:])))
                otu_id = fields[0].strip()
                otu_ids.append(otu_id)
                if has_consensus:
                    lineages.append(map(strip, fields[-1].split(';')))
    return sample_ids, otu_ids, array(otu_table), lineages

def filter_otus_by_lineage(sample_ids, otu_ids, otu_table, lineages, \
    wanted_lineage, max_seqs_per_sample, min_seqs_per_sample):
    """Filter OTU table to keep only desired lineages and sample sizes."""
    #first step: figure out which OTUs we want to keep
    if wanted_lineage is not None:  #None = keep all
        if '&&' in wanted_lineage:
            wanted_lineage = set(wanted_lineage.split('&&'))
        else:
            wanted_lineage = set([wanted_lineage])
        good_indices = []
        for i,l in enumerate(lineages):
            if set(l).intersection(wanted_lineage):
                good_indices.append(i)
        otu_table = otu_table[good_indices]
        otu_ids = map(otu_ids.__getitem__, good_indices)
        lineages = map(lineages.__getitem__, good_indices)
    #now have reduced collection of OTUs filtered by lineage.
    #figure out which samples will be dropped because too small
    big_enough_samples = (otu_table.sum(0)>=min_seqs_per_sample).nonzero()
    otu_table = otu_table[:,big_enough_samples[0]]
    sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    #figure out which samples will be reduced because too big
    too_big_samples = (otu_table.sum(0)>max_seqs_per_sample).nonzero()[0]
    if too_big_samples.shape[0]:    #means that there were some
        for i in too_big_samples:
            otu_table[:,i] = subsample(otu_table[:,i].ravel(), \
                max_seqs_per_sample)
    return sample_ids, otu_ids, otu_table, lineages

def make_envs_dict(abund_mtx, sample_names, taxon_names):
    """ makes an envs dict suitable for unifrac from an abundance matrix

    abund_mtx is samples (rows) by seqs (colunmns) numpy 2d array
    sample_names is a list, length = num rows
    taxon_names is a list, length = num columns
    """
    num_samples, num_seqs = abund_mtx.shape
    if (num_samples, num_seqs) != (len(sample_names), len(taxon_names)):
        raise ValueError, \
            "Shape of matrix %s doesn't match # samples and # taxa (%s and %s)"%\
            (abund_mtx.shape, num_samples, num_seqs)
    envs_dict = {}
    sample_names=array(sample_names)
    for i, taxon in enumerate(abund_mtx.T):
        
        nonzeros=taxon.nonzero() # this removes zero values to reduce memory
        envs_dict[taxon_names[i]] = dict(zip(sample_names[nonzeros], \
                                             taxon[nonzeros]))
    return envs_dict

def fields_to_dict(lines, delim='\t', strip_f=strip):
    """makes a dict where first field is key, rest are vals."""
    result = {}
    for line in lines:
        #skip empty lines
        if strip_f:
            fields = map(strip_f, line.split(delim))
        else:
            fields = line.split(delim)
        if not fields[0]:   #empty string in first field implies problem
            continue
        result[fields[0]] = fields[1:]
    return result

def parse_qiime_parameters(lines):
    """ Return 2D dict of params (and values, if applicable) which should be on
    """
    # The qiime_config object is a default dict: if keys are not
    # present, {} is returned
    def return_empty_dict():
        return dict()
    result = defaultdict(return_empty_dict)
    
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            fields = line.split()
            script_id, parameter_id = fields[0].split(':')
            try:
                value = fields[1]
            except IndexError:
                continue
                
            if value.upper() == 'FALSE' or value.upper() == 'NONE':
                continue
            elif value.upper() == 'TRUE':
                value = None
            else:
                pass
            
            try:
                result[script_id][parameter_id] = value
            except KeyError:
                result[script_id] = {parameter_id:value}
    return result

def sample_mapping_to_otu_table(lines):
    """Converts the UniFrac sample mapping file to an OTU table
    
    The sample mapping file is a required input for the UniFrac web interface.
    """
    out = ["#Full OTU Counts"]
    header = ["#OTU ID"]

    OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
    all_sample_names = list(all_sample_names)
    all_sample_names.sort()
    header.extend(all_sample_names)
    out.append('\t'.join(header))
    for OTU in OTU_sample_info:
        new_line = []
        new_line.append(OTU)
        for sample in all_sample_names:
            new_line.append(OTU_sample_info[OTU][sample])
        out.append('\t'.join(new_line))
    return out

def parse_sample_mapping(lines):
    """Parses the UniFrac sample mapping file (environment file)

    The sample mapping file is a required input for the UniFrac web interface.
    Returns a dict of OTU names mapped to sample:count dictionaries.
    This code is used to convert this file to an OTU table for QIIME
    """
    #add the count of 1 if count info is not supplied
    new_lines = []
    for line in lines:
        line = line.strip().split('\t')
        if len(line) == 2:
            line.append('1')
        new_lines.append(line)

    all_sample_names = [line[1] for line in new_lines]
    all_sample_names = set(all_sample_names)
    #create a dict of dicts with the OTU name mapped to a dictionary of
    #sample names with counts
    OTU_sample_info = {}
    for line in new_lines:
        OTU_name = line[0]
        if OTU_name not in OTU_sample_info:
            sample_info = dict([(i,'0') for i in all_sample_names])
            OTU_sample_info[OTU_name] = deepcopy(sample_info)
        sample_name = line[1]
        count = line[2]
        OTU_sample_info[OTU_name][sample_name] = count
    return OTU_sample_info, all_sample_names

def parse_qiime_config_file(qiime_config_file):
    """ Parse lines in a qiime_config file
    """
    result = {}
    for line in qiime_config_file:
        line = line.strip()
        # ignore blank lines or lines beginning with '#'
        if not line or line.startswith('#'): continue
        fields = line.split()
        param_id = fields[0]
        param_value = ' '.join(fields[1:]) or None
        result[param_id] = param_value
    return result
    
def parse_qiime_config_files(qiime_config_files):
    """ Parse files in (ordered!) list of qiime_config_files
    
        The order of files must be least important to most important.
         Values defined in earlier files will be overwritten if the same 
         values are defined in later files.
    """
    # The qiime_config object is a default dict: if keys are not
    # present, none is returned
    def return_none():
        return None
    results = defaultdict(return_none)
    
    for qiime_config_file in qiime_config_files:
        try:
            results.update(parse_qiime_config_file(qiime_config_file))
        except IOError:
            pass
    
    return results

def parse_tmp_to_final_filepath_map_file(lines):
    """Parses poller maps of tmp -> final file names 
     
       For example, lines:
        tmpA1.txt tmpA2.txt tmpA3.txt A.txt
        B1.txt B2.txt B3.txt B.txt
        
       Would result in:
        ([[tmpA1.txt,tmpA2.txt,tmpA3.txt], [B1.txt,B2.txt,B3.txt]],
         [A.txt,B.txt])
    
    """
    infiles_lists = []
    out_filepaths = []
    for line in lines:
        fields = line.split()
        infiles_lists.append(fields[:-1])
        out_filepaths.append(fields[-1])
    return infiles_lists, out_filepaths
    
def parse_metadata_state_descriptions(state_string):
    """From string in format 'col1:good1,good2;col2:good1' return dict."""
    result = {}
    state_string = state_string.strip()
    if state_string:
        cols = map(strip, state_string.split(';'))
        for c in cols:
            colname, vals = map(strip, c.split(':'))
            vals = map(strip, vals.split(','))
            result[colname] = set(vals)
    return result

def parse_illumina_line(l,barcode_length,rev_comp_barcode):
    """Parses a single line of Illumina data
    """
    fields = l.strip().split(':')
    
    y_position_subfields = fields[4].split('#')
    y_position = int(y_position_subfields[0])
    barcode = y_position_subfields[1][:barcode_length]
    if rev_comp_barcode:
        barcode = revComp(barcode)
    
    result = {\
     'Full description':':'.join(fields[:5]),\
     'Machine Name':fields[0],\
     'Channel Number':int(fields[1]),\
     'Tile Number':int(fields[2]),\
     'X Position':int(fields[3]),\
     'Y Position':y_position,\
     'Barcode':barcode,\
     'Full Y Position Field':fields[4],\
     'Sequence':fields[5],\
     'Quality Score':fields[6]}
     
    return result

def parse_qual_score(infile):
    """Load quality scores into dict."""
    id_to_qual = {}
    for rec in FastaFinder(infile):
        curr_id = rec[0][1:]
        curr_qual = ' '.join(rec[1:])
        try:
            parts = array(map(int, curr_qual.split()))
        except ValueError:
            raise QiimeParseError,"Invalid qual file. Check the format of the qual files." 
        curr_pid = curr_id.split()[0]
        id_to_qual[curr_pid] = parts
    return id_to_qual

def parse_qual_scores(qual_files):
    """Load qual scores into dict of {id:qual_scores}.
    
    No filtering is performed at this step.
    """
    qual_mappings = {}
    for qual_file in qual_files:
        qual_mappings.update(parse_qual_score(qual_file))
    return qual_mappings
