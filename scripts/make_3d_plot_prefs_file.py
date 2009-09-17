#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Greg Caporaso", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

from optparse import OptionParser

def parse_command_line_parameters():
    """ Parses command line arguments """
    
    usage = """usage: %prog [options] {-b COLOR_BY_STRING -p OUTPUT_PREFS_FILE}

This is a quick-and-dirty script to write prefs files to be passed via -p to 
gen_3d_plots.py. The prefs file allow for gradient coloring of continuous 
values in the 3D plots. The -b value passed in is the same as that passed in
via -b to gen_3D_plots.py: the command delimited list of fields that data 
should be included for. Currently there is only one color gradient: red to
blue, because, as mentioned, this is a quick-and-dirty script. If we decide to
stick with the pref file method for defining color gradients, we'll update
this script at that time.

"""
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    # An example string option
    parser.add_option('-b','--color_by',action='store',\
          type='string',dest='color_by',help='mapping fields to color by '+\
          '[default: %default]')
    parser.add_option('-p','--output_prefs_fp',action='store',\
          type='string',dest='output_prefs_fp',\
          help='path to store output file '+\
          '[default: %default]')
    opts,args = parser.parse_args()
    return opts,args

def build_prefs_string(color_by_string):
    fields = color_by_string.split(',')
    l = ['{']
    first = True
    entry_string = \
     "\t'%s':\n\t{\n\t\t'column':'%s',\n\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t}"
    for field in fields:
        if first:
            first=False
            l.append('\n')
        else:
            l.append(',\n')
        l.append(entry_string % (field, field))
    l.append('\n}')
    return ''.join(l)

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    out = build_prefs_string(opts.color_by)
    f = open(opts.output_prefs_fp,'w')
    f.write(out)
    f.close()