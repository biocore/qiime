#!/usr/bin/env python
#file trim_sff_primers.py: resets trim values in sff file based on primers.
from optparse import OptionParser
from string import strip
from os import walk, popen, system
from os.path import splitext, join
from qiime.parse import parse_map

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"
"""Finds the technical read regions for each library, and resets the left trim.

Replaces the sff files with the trimmed versions.
"""
sfffile_cmd = "%s -t %s -o %s %s"

sffinfo_cmd = '%s %s %s'

def get_technical_lengths(input_map, debug=False):
    """Returns per-sample info on technical lengths.
    
    Note: KEY_SEQ, BARCODE and PRIMER fields are required. LINKER optional.
    """
    if debug:
        print "Making debug output"
    map_lines = parse_map(input_map)
    header, body = map_lines[0], map_lines[1:]
    if debug:
        print "HEADER:", header
    key_index = header.index('KEY_SEQ')
    bc_index = header.index('BARCODE')
    if 'LINKER' in header:
        linker_index = header.index('LINKER')
    else:
        linker_index = None
    primer_index = header.index('PRIMER')
    technical_lengths = {}
    for fields in body:
        curr_tech_len = len(fields[key_index]) + len(fields[bc_index]) + \
            len(fields[primer_index])
        if linker_index is not None:
            curr_tech_len += len(fields[linker_index]) 
        technical_lengths[fields[0]] = curr_tech_len
    if debug:
        print "Technical lengths:"
        print technical_lengths
    return technical_lengths


def make_option_parser():
    """Generate a parser for command-line options"""
    
    usage = """\n\t python trim_sff_primers.py {-m map_file -l
    lib_dir} [options]

    [] indicates optional input (order unimportant)
    {} indicates required input (order unimportant)
    
    (For help run: \n\tpython trim_sff_primers.py --help)
    """

    parser=OptionParser(usage=usage)
    parser.add_option("-l","--libdir",dest='libdir',\
        help=""" The directory containing per-library sff files [REQUIRED]""")
    parser.add_option("-m","--input_map",dest='input_map',\
        help=""" The input map describing the libraries [REQUIRED]""")
    parser.add_option("-p","--sfffile_path",dest='sfffile_path',\
        help=""" Path to sfffile binary""", default='sfffile')
    parser.add_option("-q","--sffinfo_path",dest='sffinfo_path',\
        help=""" Path to sffinfo binary""", default='sffinfo')
    parser.add_option('--debug', dest='debug', default=False, 
        action='store_true',
        help="Print command-line for debugging")
    return parser

if __name__ == '__main__':
    option_parser = make_option_parser()
    options, args = option_parser.parse_args()
    technical_lengths = get_technical_lengths(open(options.input_map, 'U'),
        options.debug)

    for dirpath, dirnames, fnames in walk(options.libdir):
        for fname in fnames:
            if fname.endswith('.sff'):
                sff_path = join(dirpath, fname)
                lib_id = fname.rsplit('.',1)[0]
                try:
                    readlength = technical_lengths[lib_id]
                except KeyError:
                    continue
                sffinfo_cmd_to_run = sffinfo_cmd % (options.sffinfo_path,'-s',
                    sff_path)
                if options.debug:
                    print "Running sffinfo command to get ids and lengths:", \
                        sffinfo_cmd_to_run
                lines = popen(sffinfo_cmd_to_run)
                seqlengths = {}
                for line in lines:
                    if line.startswith('>'):
                        fields = line[1:].split()
                        seqlengths[fields[0]] = fields[1].split('=')[1]

                outfile_path = sff_path + '.trim'
                outfile = open(outfile_path, 'w')
                for id_, length in seqlengths.items():
                    outfile.write("%s\t%s\t%s\n" %(id_,readlength + 1,  #need +1 for 1-based index 
                        seqlengths[id_]))
                outfile.close()

                sfffile_cmd_to_run = sfffile_cmd % (options.sfffile_path, 
                    outfile_path, sff_path+'.trimmed', sff_path)
                if options.debug:
                    print "Running sfffile command:", sfffile_cmd_to_run
                system(sfffile_cmd_to_run)    
                system('mv %s.trimmed %s' % (sff_path, sff_path))
