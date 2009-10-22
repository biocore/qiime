#!/usr/bin/env python
#file make_library_id_lists.py: make id list for each lib from fasta file
from optparse import OptionParser
from string import strip
from os.path import exists, join
from os import makedirs
from collections import defaultdict

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

def get_ids(lines, field):
    """Made dict of lib:ids"""
    result = defaultdict(list)
    for line in lines:
        if line.startswith('>'):
            fields = map(strip, line[1:].split())
            label = fields[0]
            if not '_' in label:   #no lib specified
                continue
            lib, id_ = label.rsplit('_', 1)
            result[lib].append(fields[field])
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
    parser.add_option("-o","--outdir",dest='outdir',\
        default = '.',\
        help=""" The base directory to save results (one file per library).""") 
    parser.add_option("-f", "--field",dest="field", type=int,\
        default = 1,\
        help="Index of space-delimited field to read id from [DEFAULT: %default]")
    return parser

if __name__ == '__main__':
    option_parser = make_option_parser()
    options, args = option_parser.parse_args()
    ids = get_ids(open(options.in_fasta, 'U'), options.field)
    if not exists(options.outdir):
        makedirs(options.outdir)
    for k, idlist in ids.items():
        outfile = open(join(options.outdir, k + '.txt'), 'w')
        outfile.write('\n'.join(sorted(idlist)))
        outfile.close()
   
