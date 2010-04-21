#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"
 
from qiime.util import parse_command_line_parameters
from qiime.trim_sff_primers import (sfffile_cmd, sffinfo_cmd, 
    get_technical_lengths)
from sys import stderr
from os.path import join
from os import walk, popen, system
from optparse import make_option

script_info={}
script_info['brief_description']="""Trim sff primers"""
script_info['script_description']="""Finds the technical read regions for each library, and resets the left trim."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""Trim a directory of per-sff files in sff_dir (-l sff_dir/) using an input map (-m input_map.txt). This script uses the sff utility binaries which must be in your path.""","""trim_sff_primers.py -l sff_dir/ -m input_map.txt"""))
script_info['output_description']="""This script replaces the original sff files with the trimmed versions."""

script_info['required_options'] = [
    make_option("-l", "--libdir", dest='libdir',
        help="The directory containing per-library sff files"),
    make_option("-m", "--input_map", dest='input_map',
        help="Path to the input mapping file describing the libraries"),
    ]

script_info['optional_options']=[
    make_option("-p", "--sfffile_path", dest='sfffile_path', default='sfffile', 
        help="Path to sfffile binary [default: %default]"),
    make_option("-q", "--sffinfo_path", dest='sffinfo_path', default='sffinfo',
        help="Path to sffinfo binary [default: %default]"),
    make_option('--debug', dest='debug', default=False, action='store_true',
        help="Print command-line output for debugging [default: %default]"),
    ]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    technical_lengths = get_technical_lengths(open(opts.input_map, 'U'),
        opts.debug)

    for dirpath, dirnames, fnames in walk(opts.libdir):
        for fname in fnames:
            if fname.endswith('.sff'):
                sff_path = join(dirpath, fname)
                lib_id = fname.rsplit('.',1)[0]
                try:
                    readlength = technical_lengths[lib_id]
                except KeyError:
                    continue
                sffinfo_cmd_to_run = sffinfo_cmd % (opts.sffinfo_path,'-s',
                    sff_path)
                if opts.debug:
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
                    curr_length = int(seqlengths[id_])
                    left_trim = readlength + 1
                    if curr_length > left_trim:
                        outfile.write("%s\t%s\t%s\n" %(id_,left_trim,  
                            #need +1 for 1-based index 
                            curr_length))
                    else:
                        stderr.write('Rejected read %s with trim points %s and %s (orig length %s)' % (id_, left_trim, curr_length, length))
                outfile.close()

                sfffile_cmd_to_run = sfffile_cmd % (opts.sfffile_path,
                    outfile_path, sff_path+'.trimmed', sff_path)
                if opts.debug:
                    print "Running sfffile command:", sfffile_cmd_to_run
                system(sfffile_cmd_to_run)
                system('mv %s.trimmed %s' % (sff_path, sff_path))


if __name__ == "__main__":
    main()
