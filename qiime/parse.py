#!/usr/bin/env python
# file parse.py: parsers for map file, distance matrix file, env file

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald", "Greg Caporaso",
               "Justin Kuczynski", "Cathy Lozupone", "Jens Reeder",
               "Antonio Gonzalez Pena", "Jai Ram Rideout", "Will Van Treuren",
               "Yoshiki Vazquez-Baeza", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from string import strip
from collections import defaultdict
import os
from os.path import expandvars
import re
from cogent.util.dict2d import Dict2D
from cogent.util.misc import unzip
from cogent.maths.stats.rarefaction import subsample
from numpy import concatenate, repeat, zeros, nan, asarray
from numpy.random import permutation
from cogent.parse.record_finder import LabeledRecordFinder
from cogent.parse.fasta import FastaFinder
from cogent.parse.tree import DndParser
from bipy.parse.fastq import MinimalFastqParser as MinimalFastqParserBipy
from cogent.core.tree import PhyloNode
from cogent import DNA
from qiime.quality import ascii_to_phred33, ascii_to_phred64
from types import GeneratorType


def MinimalFastqParser(data, strict=False):
    return MinimalFastqParserBipy(data, strict=strict)

# this has to be here to avoid circular import


def is_casava_v180_or_later(header_line):
    """ True if this file is generated by Illumina software post-casava 1.8 """
    assert header_line.startswith('@'),\
        "Non-header line passed as input. Header must start with '@'."
    fields = header_line.split(':')
    if len(fields) == 10 and fields[7] in 'YN':
        return True
    return False


def MinimalSamParser(data):
    for line in data:
        line = line.strip()
        if not line or line.startswith('@'):
            continue
        else:
            yield line.strip().split('\t')


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
    if hasattr(lines, "upper"):
        # Try opening if a string was passed
        try:
            lines = open(lines, 'U')
        except IOError:
            raise QiimeParseError("A string was passed that doesn't refer "
                                  "to an accessible filepath.")

    if strip_quotes:
        if suppress_stripping:
            # remove quotes but not spaces
            strip_f = lambda x: x.replace('"', '')
        else:
            # remove quotes and spaces
            strip_f = lambda x: x.replace('"', '').strip()
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
            # Will add empty string to empty fields
            tmp_line = map(strip_f, line.split('\t'))
            if len(tmp_line) < len(header):
                tmp_line.extend([''] * (len(header) - len(tmp_line)))
            mapping_data.append(tmp_line)
    if not header:
        raise QiimeParseError("No header line was found in mapping file.")
    if not mapping_data:
        raise QiimeParseError("No data found in mapping file.")

    return mapping_data, header, comments


def parse_mapping_file_to_dict(*args, **kwargs):
    """Parser for map file that relates samples to metadata.

    input format: header line with fields
            optionally other comment lines starting with #
            tab-delimited fields

    calls parse_mapping_file, then processes the result into a 2d dict, assuming
    the first field is the sample id
    e.g.: {'sample1':{'age':'3','sex':'male'},'sample2':...

    returns the dict, and a list of comment lines
"""
    mapping_data, header, comments = parse_mapping_file(*args, **kwargs)
    return mapping_file_to_dict(mapping_data, header), comments


def mapping_file_to_dict(mapping_data, header):
    """processes mapping data in list of lists format into a 2 deep dict"""
    map_dict = {}
    for i in range(len(mapping_data)):
        sam = mapping_data[i]
        map_dict[sam[0]] = {}
        for j in range(len(header)):
            if j == 0:
                continue  # sampleID field
            map_dict[sam[0]][header[j]] = sam[j]
    return Dict2D(map_dict)


def parse_prefs_file(prefs_string):
    """Returns prefs dict evaluated from prefs_string.

        prefs_string: read buffer from prefs file or string containing prefs
            dict.  Must be able to evauluated as a dict using eval.
    """
    try:
        prefs = dict(eval(prefs_string))
    except TypeError:
        raise QiimeParseError(
            "Invalid prefs file. Prefs file must contain a valid prefs dictionary.")
    return prefs


def group_by_field(table, name):
    """Returns dict of field_state:[row_headers] from table.

    Use to extract info from table based on a single field.
    """
    try:
        col_index = table[0].index(name)
    except ValueError as e:
        raise ValueError("Couldn't find name %s in headers: %s" %
                         (name, table[0]))
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
        if line[0] == '\t':  # is header
            header = map(strip, line.split('\t')[1:])
        else:
            result.append(map(float, line.split('\t')[1:]))
    return header, asarray(result)


def parse_matrix(lines):
    """Parser for a matrix file Tab delimited. skips first lines if led
    by '#', assumes column headers line starts with a tab
    """
    col_headers = None
    result = []
    row_headers = []
    for line in lines:
        if line[0] == '#':
            continue
        if line[0] == '\t':  # is header
            col_headers = map(strip, line.split('\t')[1:])
        else:
            entries = line.split('\t')
            result.append(map(float, entries[1:]))
            row_headers.append(entries[0])
    return col_headers, row_headers, asarray(result)


def parse_distmat_to_dict(table):
    """Parse a dist matrix into an 2d dict indexed by sample ids.

    table: table as lines
    """

    col_headers, row_headers, data = parse_matrix(table)
    assert(col_headers == row_headers)

    result = defaultdict(dict)
    for (sample_id_x, row) in zip(col_headers, data):
        for (sample_id_y, value) in zip(row_headers, row):
            result[sample_id_x][sample_id_y] = value
    return result


def parse_bootstrap_support(lines):
    """Parser for a bootstrap/jackknife support in tab delimited text
    """
    bootstraps = {}
    for line in lines:
        if line[0] == '#':
            continue
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
            data['series'][data['options'][len(data['options']) - 1]] = \
                [float(v) for v in l[7:].strip().split('\t')]
            continue
        if l.startswith('error'):
            data['error'][data['options'][len(data['options']) - 1]] = \
                [float(v) for v in l[6:].strip().split('\t')]
        if l.startswith('color'):
            data['color'][data['options'][len(data['options']) - 1]] = \
                str(l[6:].strip())
            if(len(str(l[6:].strip())) < 1):
                print data['options'][len(data['options']) - 1]
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
    if hasattr(lines, 'next'):
        magic_check = lines.next().strip().split('\t')
    else:
        magic_check = lines[0].strip().split('\t')
        lines = lines[1:]

    if magic_check[0] != 'pc vector number':
        raise QiimeParseError("The line with the vector number was not "
                              "found, this information is required in "
                              "coordinates files")

    eigvals = None
    pct_var = None
    sample_ids = []
    result = []  # could determine n_samples, and preallocate...
    for line in lines:
        line = line.strip()
        if not line:
            continue

        fields = line.split('\t')
        values = asarray(fields[1:], dtype=float)

        if fields[0] == 'eigvals':
            eigvals = values
        elif fields[0] == '% variation explained':
            pct_var = values
        else:
            sample_ids.append(fields[0])
            result.append(values)

    # check on this information post removal of blank lines
    if eigvals is None:
        raise QiimeParseError("The line containing the eigenvalues was not "
                              "found, this information is required in coordinates files")
    if pct_var is None:
        raise QiimeParseError("The line with the percent of variation explained"
                              " was not found, this information is required in coordinates files")

    return sample_ids, asarray(result), eigvals, pct_var


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
        line = line.rstrip("\n")
        fields = line.split('\t')
        otu = fields[0].split(' ')[0]
        res[otu] = taxa_split(fields[1])

    return res
parse_observation_metadata = parse_taxonomy


def taxa_split(taxa_string):
    return [t.strip() for t in taxa_string.split(';')]


def parse_taxonomy_to_otu_metadata(
        lines, labels=['taxonomy', 'score'], process_fs=[taxa_split, float]):
    """ Return a dict mapping otu identifier to dict of otu metadata

         lines: file handle or list of lines - format should be:
          otu_id <tab> metadata entry 1 <tab> metadata entry 2 <tab> ...
         labels: list of lables for metadata entrys to be used in the
          internal dicts. each internal dict will have only as many entries
          as there are labels (extra metadata entries in the input file
          will be ignored)
         process_fs: functions which are applied to each metadata entry -
          if there are more process_fs than labels, the additional ones
          will be ignored
    """
    result = {}

    for line in lines:
        line = line.strip()
        fields = line.split('\t')
        id_ = fields[0].split()[0]
        result[id_] = {}
        for i, field in enumerate(fields[1:]):
            try:
                label = labels[i]
            except IndexError:
                continue
            try:
                value = process_fs[i](field)
            except IndexError:
                raise ValueError(
                    "Too few process functions provided (n=%d)." %
                    len(process_fs))
            result[id_][label] = value
    return result


def process_otu_table_sample_ids(sample_id_fields):
    """ process the sample IDs line of an OTU table """
    if len(sample_id_fields) == 0:
        raise ValueError('Error parsing sample ID line in OTU table. Fields are %s'
                         % ' '.join(sample_id_fields))

    # Detect if a metadata column is included as the last column. This
    # field will be named either 'Consensus Lineage' or 'OTU Metadata',
    # but we don't care about case or spaces.
    last_column_header = sample_id_fields[-1].strip().replace(' ', '').lower()
    if last_column_header in ['consensuslineage', 'otumetadata', 'taxonomy']:
        has_metadata = True
        sample_ids = sample_id_fields[:-1]
    else:
        has_metadata = False
        sample_ids = sample_id_fields

    # Return the list of sample IDs and boolean indicating if a metadata
    # column is included.
    return sample_ids, has_metadata


def parse_classic_otu_table(lines, count_map_f=int, remove_empty_rows=False):
    """parses a classic otu table (sample ID x OTU ID map)

    Returns tuple: sample_ids, otu_ids, matrix of OTUs(rows) x samples(cols),
    and lineages from infile.
    """
    otu_table = []
    otu_ids = []
    metadata = []
    sample_ids = []
    # iterate over lines in the OTU table -- keep track of line number
    # to support legacy (Qiime 1.2.0 and earlier) OTU tables
    for i, line in enumerate(lines):
        line = line.strip()
        if line:
            if (i == 1 or i == 0) and line.startswith('#OTU ID') and not sample_ids:
                # we've got a legacy OTU table
                try:
                    sample_ids, has_metadata = process_otu_table_sample_ids(
                        line.strip().split('\t')[1:])
                except ValueError:
                    raise ValueError("Error parsing sample IDs in OTU table. Appears to be a" +
                                     " legacy OTU table. Sample ID line:\n %s" % line)
            elif not line.startswith('#'):
                if not sample_ids:
                    # current line is the first non-space, non-comment line
                    # in OTU table, so contains the sample IDs
                    try:
                        sample_ids, has_metadata = process_otu_table_sample_ids(
                            line.strip().split('\t')[1:])
                    except ValueError:
                        raise ValueError("Error parsing sample IDs in OTU table." +
                                         " Sample ID line:\n %s" % line)
                else:
                    # current line is OTU line in OTU table
                    fields = line.split('\t')

                    if has_metadata:
                        # if there is OTU metadata the last column gets appended
                        # to the metadata list
                        # added in a try/except to handle OTU tables containing
                        # floating numbers
                        try:
                            valid_fields = asarray(
                                fields[1:-1],
                                dtype=count_map_f)
                        except ValueError:
                            valid_fields = asarray(fields[1:-1], dtype=float)
                        # validate that there are no empty rows
                        if remove_empty_rows and (valid_fields >= 0).all() and \
                           sum(valid_fields) == 0.0:
                            continue
                        metadata.append(map(strip, fields[-1].split(';')))
                    else:
                        # otherwise all columns are appended to otu_table
                        # added in a try/except to handle OTU tables containing
                        # floating numbers
                        try:
                            valid_fields = asarray(
                                fields[1:],
                                dtype=count_map_f)
                        except ValueError:
                            valid_fields = asarray(fields[1:], dtype=float)
                        # validate that there are no empty rows
                        if remove_empty_rows and (valid_fields >= 0.0).all() and \
                           sum(valid_fields) == 0.0:
                            continue
                    otu_table.append(valid_fields)
                    # grab the OTU ID
                    otu_id = fields[0].strip()
                    otu_ids.append(otu_id)

    return sample_ids, otu_ids, asarray(otu_table), metadata
parse_otu_table = parse_classic_otu_table


def parse_taxa_summary_table(lines):
    result = parse_classic_otu_table(lines, count_map_f=float)
    return result[0], result[1], result[2]


def filter_otus_by_lineage(sample_ids, otu_ids, otu_table, lineages,
                           wanted_lineage, max_seqs_per_sample, min_seqs_per_sample):
    """Filter OTU table to keep only desired lineages and sample sizes."""
    # first step: figure out which OTUs we want to keep
    if wanted_lineage is not None:  # None = keep all
        if '&&' in wanted_lineage:
            wanted_lineage = set(wanted_lineage.split('&&'))
        else:
            wanted_lineage = set([wanted_lineage])
        good_indices = []
        for i, l in enumerate(lineages):
            if set(l).intersection(wanted_lineage):
                good_indices.append(i)
        otu_table = otu_table[good_indices]
        otu_ids = map(otu_ids.__getitem__, good_indices)
        lineages = map(lineages.__getitem__, good_indices)
    # now have reduced collection of OTUs filtered by lineage.
    # figure out which samples will be dropped because too small
    big_enough_samples = (otu_table.sum(0) >= min_seqs_per_sample).nonzero()
    otu_table = otu_table[:, big_enough_samples[0]]
    sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    # figure out which samples will be reduced because too big
    too_big_samples = (otu_table.sum(0) > max_seqs_per_sample).nonzero()[0]
    if too_big_samples.shape[0]:  # means that there were some
        for i in too_big_samples:
            otu_table[:, i] = subsample(otu_table[:, i].ravel(),
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
        raise ValueError(
            "Shape of matrix %s doesn't match # samples and # taxa (%s and %s)" %
            (abund_mtx.shape, num_samples, num_seqs))
    envs_dict = {}
    sample_names = asarray(sample_names)
    for i, taxon in enumerate(abund_mtx.T):

        nonzeros = taxon.nonzero()  # this removes zero values to reduce memory
        envs_dict[taxon_names[i]] = dict(zip(sample_names[nonzeros],
                                             taxon[nonzeros]))
    return envs_dict


def fields_to_dict(lines, delim='\t', strip_f=strip):
    """makes a dict where first field is key, rest are vals."""
    result = {}
    for line in lines:
        # skip empty lines
        if strip_f:
            fields = map(strip_f, line.split(delim))
        else:
            fields = line.split(delim)
        if not fields[0]:  # empty string in first field implies problem
            continue
        result[fields[0]] = fields[1:]
    return result


def parse_qiime_parameters(lines):
    """ Return 2D dict of params (and values, if applicable) which should be on
    """
    # The result object is a default dict: if keys are not
    # present, {} is returned
    result = defaultdict(dict)

    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            pound_pos = line.find('#')

            # A pound sign only starts an inline comment if it is preceded by
            # whitespace.
            if pound_pos > 0 and line[pound_pos - 1].isspace():
                line = line[:pound_pos].rstrip()

            fields = line.split(None, 1)
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

            result[script_id][parameter_id] = value
    return result


def parse_qiime_config_file(qiime_config_file):
    """ Parse lines in a qiime_config file
    """
    result = {}
    for line in qiime_config_file:
        line = line.strip()
        # ignore blank lines or lines beginning with '#'
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        param_id = fields[0]
        param_value = expandvars(' '.join(fields[1:])) or None
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
            # split on the first colon to account for category names with
            # colons
            colname, vals = map(strip, c.split(':', 1))
            vals = map(strip, vals.split(','))
            result[colname] = set(vals)
    return result


def parse_illumina_line(l, barcode_length, rev_comp_barcode,
                        barcode_in_sequence=False):
    """Parses a single line of Illumina data
    """
    fields = l.strip().split(':')

    y_position_subfields = fields[4].split('#')
    y_position = int(y_position_subfields[0])
    sequence = fields[5]
    qual_string = fields[6]

    if barcode_in_sequence:
        barcode = sequence[:barcode_length]
        sequence = sequence[barcode_length:]
        qual_string = qual_string[barcode_length:]
    else:
        barcode = y_position_subfields[1][:barcode_length]

    if rev_comp_barcode:
        barcode = DNA.rc(barcode)

    result = {
        'Full description': ':'.join(fields[:5]),
        'Machine Name': fields[0],
        'Channel Number': int(fields[1]),
        'Tile Number': int(fields[2]),
        'X Position': int(fields[3]),
        'Y Position': y_position,
        'Barcode': barcode,
        'Full Y Position Field': fields[4],
        'Sequence': sequence,
        'Quality Score': qual_string}

    return result


def parse_qual_score(infile, value_cast_f=int):
    """Load quality scores into dict."""
    id_to_qual = dict([rec for rec in MinimalQualParser(infile, value_cast_f)])
    return id_to_qual


def parse_fastq_qual_score(fastq_lines):
    results = {}
    first_header = fastq_lines.readline()
    fastq_lines.seek(0)

    if is_casava_v180_or_later(first_header):
        ascii_to_phred_f = ascii_to_phred33
    else:
        ascii_to_phred_f = ascii_to_phred64

    for header, seq, qual in MinimalFastqParser(fastq_lines):
        results[header] = asarray(qual, dtype=ascii_to_phred_f)
    return results


def MinimalQualParser(infile, value_cast_f=int, full_header=False):
    """Yield quality scores"""
    for rec in FastaFinder(infile):
        curr_id = rec[0][1:]
        curr_qual = ' '.join(rec[1:])
        try:
            parts = asarray(curr_qual.split(), dtype=value_cast_f)
        except ValueError:
            raise QiimeParseError(
                "Invalid qual file. Check the format of the qual files.")
        if full_header:
            curr_pid = curr_id
        else:
            curr_pid = curr_id.split()[0]
        yield (curr_pid, parts)


def parse_qual_scores(qual_files):
    """Load qual scores into dict of {id:qual_scores}.

    No filtering is performed at this step.
    """
    qual_mappings = {}
    for qual_file in qual_files:
        qual_mappings.update(parse_qual_score(qual_file))
    return qual_mappings


def parse_trflp(lines):
    """Load a trflp file and returns a header and data lists"""

    sample_ids = []
    otu_ids = []
    data = []
    non_alphanum_mask = re.compile('[^\w|^\t]')
    # not sure why the above regex doesn't cover the following regex...
    dash_space_mask = re.compile('[_ -]')

    for i, line in enumerate(lines):
        elements = line.strip('\n').split('\t')

        # special handling for the first line only
        if i == 0:
            # validating if the file has a header
            if elements[0] == '':
                for otu_id in elements[1:]:
                    otu_ids.append(non_alphanum_mask.sub('_', otu_id))
                continue
            else:
                for j, otu_id in enumerate(elements[1:]):
                    otu_ids.append(non_alphanum_mask.sub('_', 'Bin%3d' % j))

        # handling of all other lines
        current_row = []

        # converting each value in the row to int
        for count in elements[1:]:
            try:
                current_row.append(int(round(float(count), 0)))
            except ValueError:
                current_row.append(0)

        # if the sum of all the values is equial to 0 ignore line
        if sum(current_row) == 0:
            continue

        # adding sample header to list
        sample_ids.append(non_alphanum_mask.sub('.',
                          dash_space_mask.sub('.', elements[0])))

        # validating the size of the headers to add missing columns
        # this is only valid when there is no header
        if len(current_row) > len(otu_ids):
            # modify header data
            extra_cols = []
            for j in range(len(otu_ids), len(current_row)):
                extra_cols.append(non_alphanum_mask.sub('_', 'Bin%3d' % j))
            # modify data
            for j in range(len(data)):
                data[j].extend([0] * (len(current_row) - len(otu_ids)))

            otu_ids.extend(extra_cols)
        elif len(current_row) < len(otu_ids):
            # modify data
            current_row.extend([0] * (len(otu_ids) - len(current_row)))

        data.append(current_row)

    return sample_ids, otu_ids, asarray(data).transpose()


def parse_denoiser_mapping(denoiser_map):
    """ read a denoiser mapping file into a dictionary """
    result = {}
    for line in denoiser_map:
        line = line.strip().split('\t')
        denoised_id = line[0].rstrip(':')
        original_ids = [denoised_id] + line[1:]
        if denoised_id in result:
            # just a healthy dose of paranoia
            raise ValueError("Duplicated identifiers in denoiser mapping file: "
                             "are you sure you merged the correct files?")
        else:
            result[denoised_id] = original_ids
    return result


def parse_otu_map(otu_map_f, otu_ids_to_exclude=None, delim='_'):
    """ parse otu map file into a sparse dict {(otu_idx,sample_idx):count}

        This function is much more memory efficent than fields_to_dict and
         and the result dict is of the correct format to be passed to
         table_factory for creating OtuTable objects.

    """
    if otu_ids_to_exclude is None:
        otu_ids_to_exclude = {}

    result = defaultdict(int)
    sample_ids = []
    sample_id_idx = {}
    otu_ids = []
    otu_count = 0
    sample_count = 0
    for line in otu_map_f:
        fields = line.strip().split('\t')
        otu_id = fields[0]
        if otu_id in otu_ids_to_exclude:
            continue
        for seq_id in fields[1:]:
            sample_id = seq_id.split(delim)[0]
            try:
                sample_index = sample_id_idx[sample_id]
            except KeyError:
                sample_index = sample_count
                sample_id_idx[sample_id] = sample_index
                sample_count += 1
                sample_ids.append(sample_id)
            # {(row,col):val}
            result[(otu_count, sample_index)] += 1
        otu_count += 1
        otu_ids.append(otu_id)
    return result, sample_ids, otu_ids


def parse_sample_id_map(sample_id_map_f):
    """Parses the lines of a sample ID map file into a dictionary.

    Returns a dictionary with original sample IDs as the keys and new sample
    IDs as the values.

    This function only allows a sample ID map to perform one-to-one mappings
    between sample IDs (e.g. S1 and T1 point to new ID 'a', but a third
    original ID, such as S2, cannot also point to 'a').

    Arguments:
        sample_id_map_f - the lines of a sample ID map file to parse. Each line
            should contain two sample IDs separated by a tab. Each value in the
            first column must be unique, since the returned data structure is a
            dictionary using those values as keys
    """
    result = {}
    new_samp_id_counts = defaultdict(int)
    for line in sample_id_map_f:
        # Only try to parse lines that aren't just whitespace.
        line = line.strip()
        if line:
            samp_id, mapped_id = line.split('\t')
            if samp_id in result:
                raise ValueError("The first column of the sample ID map must "
                                 "contain unique sample IDs ('%s' is "
                                 "repeated). The second column, however, may "
                                 "contain repeats." % samp_id)
            elif new_samp_id_counts[mapped_id] >= 2:
                raise ValueError("Only two original sample IDs may map to the "
                                 "same new sample ID. The new sample ID '%s' "
                                 "has more than two sample IDs mapping to it."
                                 % mapped_id)
            else:
                result[samp_id] = mapped_id
                new_samp_id_counts[mapped_id] += 1
    return result
