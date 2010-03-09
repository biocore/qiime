#!/usr/bin/env python
#file parse.py: parsers for map file, distance matrix file, env file

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso",
    "Justin Kuczynski", "Cathy Lozupone", "Jens Reeder"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

from string import strip
from collections import defaultdict
from cogent.util.misc import unzip
from cogent.maths.stats.rarefaction import subsample
from numpy import array, concatenate, repeat, zeros
from numpy.random import permutation
from cogent.parse.record_finder import LabeledRecordFinder
from copy import deepcopy
import os

"""Parsers for internally used file formats from SIBS, OTUPicker, etc.

Note: this code initially copied over from MicrobePlots."""

def parse_map(lines, return_header=False, strip_quotes=True, \
 suppress_stripping=False):
    """Parser for map file that relates samples to metadata.
    
    Format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    Result: list of lists of fields, incl. headers.
    """
    if strip_quotes:
        filter_f = lambda x: x.strip().replace('"','').strip()
    else:
        filter_f = strip
    result = []
    header = []
    for line in lines:
        # Added for unadultered parsing, need to know if white space present
        if suppress_stripping:
            test_for_empty_line = line.strip()
            if not test_for_empty_line:
                continue
        else:
            line = line.strip()
            if not line:
                continue
        
        if line.startswith('"'): #quoted line from Excel?
            fields = line.split('\t')
            fields[0] = fields[0].strip('"')
            line = '\t'.join(fields)
        if line.startswith('#') and (header or result):
            header.append(line[1:].strip())
        else:
            result.append(map(filter_f, line.split('\t')))

    if return_header:
        return result, header
    return result

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

def parse_minimal_distmat(lines):
    """Parser for raw fast_unifrac output."""
    lines = list(lines)
    header = eval(lines[0])
    matrix = eval(''.join(lines[1:]))
    return header, matrix

def is_rarefaction_label_line(line):
    """Returns True if line is rarefaction label"""
    return line.startswith('#HEADER')

rrf = rarefaction_record_finder = LabeledRecordFinder(is_rarefaction_label_line)

def extract_otu_header_fields(line):
    """Extracts fields from header line."""
    fields = map(strip, line.split('\t'))
    pct_sim = float(fields[1])
    sample_id = fields[2]
    num_iters = int(fields[3])
    return pct_sim, sample_id, num_iters

def extract_index_fields(line):
    """Extracts index fields from a given line.

    Assumes first field is name (with #), then mean, lower, upper.
    """
    fields = map(strip, line.split('\t'))
    index_name = fields[0][1:]
    return index_name, map(float, fields[1:])

def extract_rarefaction_table(lines):
    """Extracts rarefaction table from the rest of the lines"""
    return array([map(float, l.split('\t')) for l in lines])

def parse_rarefaction_rec(rec):
    """Parses single rarefaction record"""
    result = {}
    for i, line in enumerate(rec):
        line = line.strip()
        if not line:
            continue
        if is_rarefaction_label_line(line):
            result['pct_sim'], result['sample_id'], result['num_iters'] = \
                extract_otu_header_fields(line)
        elif line.startswith('#n'): #rest of rec is table
            result['rarefaction_data'] = extract_rarefaction_table(rec[i+1:])
            break
        else:
            name, fields = extract_index_fields(line)
            result[name] = fields
    return result

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

def parse_rarefaction(lines):
    """Parser for rarefaction data file.

    The examples I have of this file share the following format:
    - header line starting with #HEADER, with %id, sample_id, num_reps(?)
    - lines for chao1, ace and shannon with mean, LC, UC
    - line for simpson with actual value of index
    - table header with n, rare, rare_lci, rare_uci
    - multiline 4-col output of this table.

    WARNING: vals in 1st col of table not same for all samples either because
    of limited # seqs or because sampling performed at equal # intervals not
    equal # seqs.

    For PD rarefaction, is the same but without the other stats calculated,
    i.e. just has per-sample data for the rarefaction, and the header.

    Desired result:
    dict keyed by sample id, vals are dicts containing all the info (e.g.
    rarefaction_data as 4-col array, other indices as either 1-col or 3-col
    lists of numbers, metadata about pct_sim and num_iters
    """
    result = {}
    for rec in rrf(lines):
        curr = parse_rarefaction_rec(rec)
        result[curr['sample_id']] = curr
    return result

def parse_rarefaction_fname(name_string):
    """returns base, seqs/sam, iteration, extension.  seqs, iters as ints
    
    all as strings, some may be empty strings ('')"""

    root, ext = os.path.splitext(name_string)
    root_list = root.split("_")
    iters = int(root_list.pop())
    seqs_per_sam = int(root_list.pop())
    base_name = "_".join(root_list)
    return base_name, seqs_per_sam, iters, ext

def otu_file_to_lineages(infile):
    """Returns lineage assignments for each otu in infile.
    
    Result is in format {otu_id:[taxonomy, support]}
    """
    result = {}
    for line in infile:
        if line.startswith('# OTU'): #is OTU line
            curr_otu = line.split()[2]
            next_line = infile.readline()
            support, tax = next_line.split(' from ', 1)
            result[curr_otu] = [map(strip, tax.split(';')),float(support[2:-1])]
    return result

def parse_taxonomy(infile):
    """parse a taxonomy file.

    Returns: dict of otu id to taxonomy name.
    ignores other parts of the otu file, such as confidence and seq id (otu id
    only)
    """

    res = {}
    for line in infile:
        fields = line.split('\t')
        # typically this looks like: 3 SAM1_32 \t Root,Bacteria,Fi... \t 0.9
        # implying otu 3; sample 1, seq 32 (the representative of otu 3);
        # followed by the taxonomy and confidence
        if not len(fields) == 3:
            continue
        otu = fields[0].split(' ')[0]
        res[otu] = fields[1]

    return res

def parse_otus(lines):
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
            fields = line.split('\t')
            #first and last col are otu id and consensus lineage respectively
            if has_consensus:
                otu_table.append(array(map(int, fields[1:-1])))
            else:
                otu_table.append(array(map(int, fields[1:])))
            otu_id = fields[0].strip()
            otu_ids.append(otu_id)
            if has_consensus:
                lineages.append(map(strip, fields[-1].split(';')))
    return sample_ids, otu_ids, array(otu_table), lineages

def otu_table_to_envs(sample_ids, otu_ids, otu_table):
    """Convert otu matrix to envs table.

    result is dict of {otu_name:{env_count}}
    """
    result = {}
    for otu, counts in zip(otu_ids, otu_table):
        result[otu] = dict([i for i in zip(sample_ids, counts) if i[1].any()])
    return result

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

def parse_sequences_by_otu(infile):
    """Parse sequences_by_otu file into two dicts: {seq:OTU} and {OTU:seqs}."""
    seq_to_otu = {}
    otu_to_seqs = defaultdict(list)
    curr_otu = None
    for line in infile:
        if line.startswith('# OTU '):
            start, rest = line.split('# OTU ')
            curr_otu = int(rest.split()[0])
        elif line.startswith('>'):
            label = line.split()[0][1:]
            orig_seqs = label.split('@@')
            otu_to_seqs[curr_otu].extend(orig_seqs)
            for s in orig_seqs:
                seq_to_otu[s]=curr_otu
    return otu_to_seqs, seq_to_otu

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
    
def envs_to_otu_counts(lines):
    """Reads envs lines into OTU counts {(sampleid,otu_id):count}."""
    result = defaultdict(int)
    for line in lines:
        fields = map(strip, line.split('\t'))
        if len(fields) != 3:
            continue
        result[(fields[1], fields[0])] += int(fields[2])
    return result
    
def otu_counts_to_matrix(otu_counts):
    """Build otu matrix from dict of {(sampleid, otu_id):count}.

    Adapted from Daniel McDonald's script.
    """
    all_sampleids, all_otus = unzip(otu_counts.keys())
    all_sampleids = sorted(set(all_sampleids))
    all_otus = sorted(set(all_otus))
    matrix = zeros((len(all_otus), len(all_sampleids)), int)
    for row, otu in enumerate(all_otus):
        for col, sampleid in enumerate(all_sampleids):
            matrix[row, col] += otu_counts.get((sampleid, otu), 0)
    return matrix, all_otus, all_sampleids

def envs_to_matrix(lines):
    """Reads envs lines into matrix of OTU counts, plus row/col"""
    return otu_counts_to_matrix(envs_to_otu_counts(lines))

# Start functions for handling qiime_parameters file
def parse_qiime_parameters(lines):
    """ Return 2D dict of params (and values, if applicable) which should be on
    """
    result = {}
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
    lines = [line.strip().split('\t') for line in lines]

    all_sample_names = [line[1] for line in lines]
    all_sample_names = set(all_sample_names)
    #create a dict of dicts with the OTU name mapped to a dictionary of
    #sample names with counts
    OTU_sample_info = {}
    for line in lines:
        OTU_name = line[0]
        if OTU_name not in OTU_sample_info:
            sample_info = dict([(i,'0') for i in all_sample_names])
            OTU_sample_info[OTU_name] = deepcopy(sample_info)
        sample_name = line[1]
        count = line[2]
        OTU_sample_info[OTU_name][sample_name] = count
    return OTU_sample_info, all_sample_names

# End functions for handling qiime_parameters file

