#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 
from qiime.util import parse_command_line_parameters
from qiime.trim_sff_primers import (sfffile_cmd, sffinfo_cmd, 
    get_technical_lengths)
from os.path import join
from os import walk, popen, system
from optparse import make_option

script_description = """Finds the technical read regions for each library, and resets the left trim.

Replaces the sff files with the trimmed versions.
 """

script_usage = """Trim a directory of per-sff files at sff_dir (-l sff_dir)
using input map input_map.txt (-m input_map.txt) using sff utils binaries that are on your path already:

trim_sff_primers -l sff_dir -m input_map.txt
"""

required_options = [\
    make_option("-l","--libdir",dest='libdir',\
        help=""" The directory containing per-library sff files"""),
    make_option("-m","--input_map",dest='input_map',\
        help=""" The input map describing the libraries""")
]

optional_options = [\
    make_option("-p","--sfffile_path",dest='sfffile_path',\
        help=""" Path to sfffile binary [default: %default]""", 
        default='sfffile'),
    make_option("-q","--sffinfo_path",dest='sffinfo_path',\
        help=""" Path to sffinfo binary [default: %default]""", 
        default='sffinfo'),
    make_option('--debug', dest='debug', default=False,
        action='store_true',
        help="Print command-line for debugging [default: %default]")
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

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
                    outfile.write("%s\t%s\t%s\n" %(id_,readlength + 1,  
                        #need +1 for 1-based index 
                        seqlengths[id_]))
                outfile.close()

                sfffile_cmd_to_run = sfffile_cmd % (opts.sfffile_path,
                    outfile_path, sff_path+'.trimmed', sff_path)
                if opts.debug:
                    print "Running sfffile command:", sfffile_cmd_to_run
                system(sfffile_cmd_to_run)
                system('mv %s.trimmed %s' % (sff_path, sff_path))


if __name__ == "__main__":
    main()
