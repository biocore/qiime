#!/usr/bin/env python
# File created on 08 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the Qiime project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


from optparse import OptionParser
from qiime.pick_otus import expand_otu_map_seq_ids, map_otu_map_files,\
    write_otu_map

usage_str = """usage: %prog [options] {-i COMMA_SEPARATED_OTU_FPS -o OUTPUT_OTU_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

If the seq_ids in otu_table2.txt are otu_ids in otu_table1.txt, expand 
 the seq_ids in otu_table2.txt to be the full list of associated seq_ids from
 otu_table1.txt. Write the resulting otu table to otu_table.txt (-o).

python Qiime/scripts/merge_otu_maps.py -i otu_table1.txt,otu_table2.txt -o otu_table.txt

 For example, if otu_table1.txt contains:
    0	seq1	seq2	seq5
    1	seq3	seq4
    2	seq6	seq7	seq8
    
 and otu_table2.txt contains:
    110	0	2
    221	1
    
 The resulting OTU table will be:
    110	seq1	seq2	seq5	seq6	seq7	seq8
    221	seq3	seq4
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--otu_map_fps',\
         help='the otu map filepaths, comma-separated and '+\
          'ordered as the OTU pickers were run [REQUIRED]')
    parser.add_option('-o','--output_fp',\
         help='path to write output OTU map [REQUIRED]')
         
    parser.set_defaults(verbose=False)

    opts,args = parser.parse_args()
    required_options = ['otu_map_fps']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
    
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    otu_files = map(open,opts.otu_map_fps.split(','))
    
    try:
        otu_map = map_otu_map_files(otu_files)
    except KeyError,e:
        print 'Some keys do not map ('+ str(e) +') -- is the order of'+\
        ' your OTU maps equivalent to the order in which the OTU pickers'+\
        ' were run?'
        exit(1)
        
    write_otu_map(otu_map,opts.output_fp)
    
    