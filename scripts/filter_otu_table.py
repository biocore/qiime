#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight, Jesse Stombaugh, and Dan Knights"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from qiime.filter_otu_table import _filter_table, split_tax
from optparse import make_option
from string import strip

script_description = """Filters OTU table according to minimum OTU count and number of samples.
If OTU table has taxonomy assigned, can also use taxonomy to filter."""

script_usage = """Creates a new filtered OTU with the original filename + "_filtered.txt".
To retain only OTUs with at least this many sequences:
filter_otu_table -i otu_table.txt -c MIN_COUNT 

To retain only OTUs found in at least this many samples:
filter_otu_table -i otu_table.txt -s MIN_SAMPLES

To retain only OTUs containing particular taxonomic terms in their lineage:
filter_otu_table -i otu_table.txt -t "Actinobacteria,Bacteroidetes"

To exclude OTUs containing particular taxonomic terms in their lineage:
filter_otu_table -i otu_table.txt -e "Actinobacteria,Bacteroidetes"
"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
    make_option('-i', '--otu_filename', dest='otu_fname',
                help='otu file name')\
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
    make_option('-c', '--min_count', default=1, type=int,
        help='retain OTUs with at least this many sequences [default=%default]'),\
    make_option('-s', '--min_samples', default=2, type=int,
        help='retain OTUs found in at least this many samples [default=%default]'),\
    make_option('-t', '--include_taxonomy', default='',
        help='list of taxonomy terms to include [default=%default]'),\
    make_option('-e', '--exclude_taxonomy', default='', 
        help='list of taxonomy terms to exclude [default=%default]'),\
    make_option('-o', '--dir_path', default='./',
        help='directory prefix for all analyses [default=%default]')
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    # process options (was originally process_options())
    filepath=opts.otu_fname
    filename=filepath.strip().split('/')[-1]
    filename=filename.split('.')[0]

    params={}
    params['otu_file'] = opts.otu_fname
    params['min_otu_count'] = opts.min_count
    params['min_otu_samples'] = opts.min_samples

    if opts.include_taxonomy:
        included_taxa = set(map(strip, split_tax(opts.include_taxonomy)))
    else:
        included_taxa = set()

    if opts.exclude_taxonomy:
        excluded_taxa = set(map(strip, split_tax(opts.exclude_taxonomy)))
    else:
        excluded_taxa=set()

    params['included_taxa']=included_taxa
    params['excluded_taxa']=excluded_taxa
    params['dir_path']=opts.dir_path

    _filter_table(params)

if __name__ == "__main__":
    main()
