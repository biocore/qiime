#!/usr/bin/env python
#file make_library_id_lists.py: make id list for each lib from fasta file
from optparse import OptionParser
from string import strip
from os.path import exists, join
from os import makedirs
from collections import defaultdict

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def get_ids(lines, field, bad_ids=None, debug=False):
    """Make dict of lib:ids"""
    result = defaultdict(list)
    for line in lines:
        if line.startswith('>'):
            fields = map(strip, line[1:].split())
            label = fields[0]
            if not '_' in label:   #no lib specified
                continue
            lib, id_ = label.rsplit('_', 1)
            if bad_ids and label in bad_ids:
                if debug:
                    print "Excluded bad id: %s" % label
            else:
                result[lib].append(fields[field])
    return result

def get_first_id(lines):
    """Gets first fasta id from each line in lines"""
    result = set()
    for line in lines:
        if line.startswith('>'):
            result.add(line[1:].split()[0])
    return result

def make_option_parser():
    """Generate a parser for command-line options"""
    
    usage = """\n\t python make_library_id_lists.py {-i input_fasta_fp -o
    output_dir} [options]

    [] indicates optional input (order unimportant)
    {} indicates required input (order unimportant)
    \nExample:\n\tpython make_library_id_lists.py -i
    seqs.fna -o /Users/rob/results
    
    (For help run: \n\tpython make_library_id_lists.py --help)
    """

    parser=OptionParser(usage=usage)
    parser.add_option("-i","--input_fasta",dest='in_fasta',default = None,\
        help="The path to a FASTA file containing input sequences [REQUIRED]")
    parser.add_option("-s", "--screened_rep_seqs",dest="screened_rep_seqs",
        default=None,
        help="The path to a FASTA file containing screened representative seqs" +
        "[DEFAULT: %default]")
    parser.add_option("-u", "--otus",dest="otus",
        default=None,
        help="The path to an OTU file mapping OTUs onto rep seqs" +
        "[DEFAULT: %default]")
    parser.add_option("-o","--outdir",dest='outdir',\
        default = '.',\
        help=""" The base directory to save results (one file per library).""") 
    parser.add_option("-f", "--field",dest="field", type=int,\
        default = 1,\
        help="Index of space-delimited field to read id from [DEFAULT: %default]")
    parser.add_option("--debug", dest="debug", action="store_true",
    default=False, help="Show debug output.")
    return parser
