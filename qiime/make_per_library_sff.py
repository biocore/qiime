#!/usr/bin/env python
#file make_library_id_lists.py: make id list for each lib from fasta file
from optparse import OptionParser
from string import strip
from os import walk, system
from os.path import splitext, join

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

cmd = "sfffile -i lib_ids/S100.txt -o S100.sff E86FECS01.sff E86FECS02.sff"
cmd = "%s -i %s -o %s %s"

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
    
    usage = """\n\t python make_per_library_sff.py {-i input_sff -l
    lib_dir} [options]

    [] indicates optional input (order unimportant)
    {} indicates required input (order unimportant)
    \nExample:\n\tpython make_per_library_sff.py -i
    input_sff -l /Users/rob/libs
    
    (For help run: \n\tpython make_library_id_lists.py --help)
    """

    parser=OptionParser(usage=usage)
    parser.add_option("-i","--input_sff",dest='in_sff',\
        help="The path to an input sff file (or files: separate w/ comma, no spaces) [REQUIRED]")
    parser.add_option("-l","--libdir",dest='libdir',\
        help=""" The directory containing per-library id files [REQUIRED]""")
    parser.add_option("-p","--sfffile_path",dest='sfffile_path',\
        help=""" Path to sfffile binary""", default='sfffile')

        

    parser.add_option("-f", "--field",dest="field", type=int,\
        default = 1,\
        help="Index of space-delimited field to read id from [DEFAULT: %default]")
    return parser

if __name__ == '__main__':
    option_parser = make_option_parser()
    options, args = option_parser.parse_args()
    input_sff_names = options.in_sff.replace(',',' ')
    for dirpath, dirnames, fnames in walk(options.libdir):
        for fname in fnames:
            if fname.startswith('.'):
                continue
            basename, ext = splitext(fname)
            system(cmd % (options.sfffile_path, join(dirpath,fname), 
                join(dirpath, basename+'.sff'), input_sff_names))
