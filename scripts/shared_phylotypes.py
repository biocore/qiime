#!/usr/bin/env python
# File created on 12 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"
__status__ = "Release"
 
from glob import glob
from os.path import exists, split, splitext, isdir
from qiime.util import make_option

from qiime.util import parse_command_line_parameters,\
    get_options_lookup, create_dir
from qiime.shared_phylotypes import calc_shared_phylotypes


options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Compute shared OTUs between all pairs of samples"
script_info['script_description'] = "This script computes from an OTU table a matrix with the number of shared phylotypes between all pairs of samples."
script_info['script_usage'] = [
    ("Single example","Compute shared OTUs on one OTU table for all samples","%prog -i otu_table.biom -o shared_otus.txt"),
    ("Reference sample example","Compute shared OTUs with respect to a reference sample. Computes shared OTUs between all pairs of samples and the reference sample. E.g. in a transplant study this can be used to establish a base line count of shared OTUs with the Donor sample before and after the transplant.","%prog -i otu_table.biom -o shared_otus_PC.636.txt -r PC.636"),

    ("Batch mode example","Compute shared OTUs for a set of OTU tables, e.g. from running multiple_rarefactions.py, with an even number of sequences per sample. The resulting directory can be fed to dissimilarity_mtx_stats.py, which computes mean, median and the standard deviation on the provided tables.", "%prog -i rarefied_otu_tables/ -o shared_otus/")
]

script_info['output_description']= ""
script_info['required_options'] = [\
     make_option('-i','--otu_table_fp',type='existing_path',\
                     help='path to the input OTU table in biom format or a directory containing OTU tables'),
    options_lookup['output_fp']
]

script_info['optional_options'] = [\
 # Example optional option
 make_option('-r','--reference_sample',type='string',
             help='Name of reference sample to which all pairs of samples should be compared '
             + '[default: %default]', default=None),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if isdir(opts.otu_table_fp):
        ret_code = create_dir(opts.output_fp, fail_on_exist=False)
        #run on each file in dir
        for fp in glob(opts.otu_table_fp +'/*biom'):
            parent_dir_name, file_name = split(fp)
            basename, extension = splitext(file_name)
            out_fp = opts.output_fp +"/"+basename+"_shared_OTUs.txt"
            out_fh = open(out_fp, "w")
            
            out_fh.write(calc_shared_phylotypes(open(fp, "U"),
                                        opts.reference_sample))
            out_fh.close()
    else:
        #run in single file mode
        try:
            out_fh = open(opts.output_fp,"w")
        except IOError, message:
            exit(("Can't open output file %s for writing. Check the permissions "\
                     +"or existing directory with identical name.\n%s")
                 % (opts.output_fp,message))
        out_fh.write(calc_shared_phylotypes(open(opts.otu_table_fp, "U"),
                                            opts.reference_sample))
                
if __name__ == "__main__":
    main()
