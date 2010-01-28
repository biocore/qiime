#!/usr/bin/env python

"""Parse fna and qual files into tab delimited text

Script for converting fna and qual files into tab delimited text for easy bulk
import into databases. This script uses the fasta parser code taken from 
PyCogent (http://pycogent.sourceforge.net), and is written to be run 
indepenedent of the PyCogent project for ease of portability.

output format:
seqid\tseq_data\tcomments

if has_sampleid is set:
seqid\tsampleid\tseq\tcomments
"""

from sys import argv
from string import strip
from optparse import make_option, OptionParser
import re, os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

#TODO: add useful help messages
options = [make_option('--input',dest='input',default=None),
           make_option('--output',dest='output',default=None),
           make_option('--has_sampleid',dest='has_sampleid',\
                   action='store_true', default=False)]

def is_fasta_label(x):
    """Checks if x looks like a FASTA label line."""
    return x.startswith('>')

def is_empty(line):
    """Returns True empty lines and lines consisting only of whitespace."""
    return (not line) or line.isspace()

def LabeledRecordFinder(is_label_line, constructor=strip, ignore=is_empty):
    """Returns function that returns successive labeled records from file.

    Includes label line in return value. Returns list of relevant lines.
    
    Default constructor is string.strip, but can supply another constructor
    to transform lines and/or coerce into correct type. If constructor is None,
    passes along the lines without alteration.

    Skips over any lines for which ignore(line) evaluates True (default is
    to skip empty lines).
    
    NOTE: Does _not_ raise an exception if the last line is a label line: for
    some formats, this is acceptable. It is the responsibility of whatever is
    parsing the sets of lines returned into records to complain if a record
    is incomplete.
    """
    def parser(lines):
        curr = []
        for l in lines:
            if constructor:
                line = constructor(l)
            else:
                line = l
            if ignore(line):
                continue
            #if we find the label, return the previous record
            if is_label_line(line):
                if curr:
                    yield curr
                    curr = []
            curr.append(line)
        #don't forget to return the last record in the file
        if curr:
            yield curr
    return parser

def is_blank_or_comment(x):
    """Checks if x is blank or a FASTA comment line."""
    return (not x) or x.startswith('#') or x.isspace()

qual_constructor = lambda s: strip(s, '\n')
QualFinder = LabeledRecordFinder(is_fasta_label, constructor=qual_constructor, \
        ignore=is_blank_or_comment)
FastaFinder = LabeledRecordFinder(is_fasta_label, ignore=is_blank_or_comment)

def MinimalFastaParser(infile, strict=True, \
    label_to_name=str, finder=FastaFinder, \
    is_label=None, label_characters='>'):
    """Yields successive sequences from infile as (label, seq) tuples.

    If strict is True (default), raises RecordError when label or seq missing.
    """
    if is_label:
        deprecated('argument', 'is_label', 'label_characters', 1.4)
    else:
        # use an re to search for line starting
        label_pattern = re.compile('^[%s]' % label_characters)
        is_label = label_pattern.search

    for rec in finder(infile):
        #first line must be a label line
        if not is_label(rec[0]):
            if strict:
                raise RecordError, "Found Fasta record without label line: %s"%\
                    rec
            else:
                continue
        #record must have at least one sequence
        if len(rec) < 2:
            if strict:
                raise RecordError, "Found label line without sequences: %s" % \
                    rec
            else:
                continue

        label = rec[0][1:].strip()
        label = label_to_name(label)

        if finder is QualFinder:
            clean_recs = []
            for r in rec[1:]:
                if not r.endswith(' '):
                    r = r + ' '
                clean_recs.append(r)
            rec[1:] = clean_recs
        seq = ''.join(rec[1:])

        yield label, seq

def load_input(input_file, is_qual=False):
    """Takes an open file and yields (id_with_comment, sequence)"""
    if is_qual:
        finder = QualFinder
    else:
        finder = FastaFinder

    for seqname, seq in MinimalFastaParser(input_file, finder=finder):
        yield (seqname, seq)

def delim_data(data_gen, has_sampleid, delim='\t'):
    """Takes (seqname, seq) from data_gen and returns delimited text
    
    if has_sampleid, returns seq_id, sampleid, seq, comment
    """
    for seqname, seq in data_gen:
        seqname_fields = seqname.split()
        seqid = seqname_fields[0]
        comment = ' '.join(seqname_fields[1:])

        if has_sampleid:
            sample_id = seqid.split('_',1)[0]
            yield delim.join([seqid, sample_id, seq, comment])
        else:
            yield delim.join([seqid, seq, comment])

def main():
    parser = OptionParser(option_list=options)
    opts, args = parser.parse_args()

    if opts.input is None:
        parser.error("Must specify either an input file or directory")
    if opts.output is None:
        opts.output = os.getcwd()
    
    if os.path.isdir(opts.input) and not os.path.exists(opts.output):
        try:
            os.mkdir(opts.output)
        except:
            parser.error("Unable to create directory %s" % opts.output)
    
    if os.path.isdir(opts.input):
        for f in os.listdir(opts.input):
            rel_in_f = os.path.join(opts.input, f)
            rel_out_f = os.path.join(opts.output, f + '.delim')

            if f.endswith('.qual'):
                data = load_input(open(rel_in_f), is_qual=True)
            else:
                data = load_input(open(rel_in_f), is_qual=False)

            data_gen = delim_data(data, opts.has_sampleid)

            output = open(rel_out_f, 'w')
            output.write('\n'.join(data_gen))
            output.close()
    else:
        if opts.input.endswith('.qual'):
            data = load_input(open(opts.input), is_qual=True)
        else:
            data = load_input(open(opts.input), is_qual=False)

        rel_out_f = opts.output + '.delim'
        data_gen = delim_data(data, opts.has_sampleid)
        output = open(rel_out_f, 'w')
        output.write('\n'.join(data_gen))
        output.close()

if __name__ == '__main__':
    main()
